
""" Multipole Math for the Fast Multipole Method 

In a 2D plane gravitational force takes the form 1 / (x - y), which follows
from a log(x - y) potential function where x and y are complex.

Upward pass components: P2M, M2M
Downwd pass components: M2L, L2L, NNC   """


include("Quadtree.jl")

"""

Wrapper for all components


"""

function FMM!(quadtree::Quadtree, points, masses, ω_p, binomial_table)
  # upward pass
  P2M!(quadtree, points, masses)
  M2M!(quadtree, binomial_table)
  # downward pass
  @btime M2L!(quadtree, $binomial_table)
  L2L!(quadtree)
  L2P!(quadtree, points, ω_p)
  NNC!(quadtree, points, masses, ω_p)
end

"""

Upward pass functions

"""


""" P2M : Particle to Multipole """

# TESTING PARAM:  P = 33, N = 1000
# UNOPTIMIZED:    516.538 μs (2240 allocations: 737.73 KiB)
# OPTIMIZED:      77.243 μs (128 allocations: 42.16 KiB)
function P2M!(quadtree::Quadtree, points, masses::Array{Float64, 1})
  leaf_offset::Int = getOffsetOfDepth(quadtree, quadtree.tree_depth)
  # determine multipole expansion at all leaf boxes
  for global_idx in leaf_offset+1:length(quadtree.tree)
    box::Box = quadtree.tree[global_idx]
    # builtin sum uses cumsum so using the builtin for better error
    if (boxHasPoints(box))
      relevant_points = @view points[box.start_idx:box.final_idx]
      relevant_masses = @view masses[box.start_idx:box.final_idx]
      box.a[1] = sum(relevant_masses)
      
      # OPTIMIZED array power operation
      diff = relevant_points .- box.center 
      tmp = diff.*relevant_masses
      for k in 1:P
        i = k + 1
        box.a[i] = -1/k * sum(tmp)
        tmp .*= diff
      end
      # UNOPTIMIZED
      # for i in 2:P+1
      #   # subtract box center to form multipole expansion about box center
      #   k = i - 1
      #   box.a[i] = -1/k*sum(((relevant_points .- box.center).^k).*relevant_masses)
      #   #@printf "%f + %fi, " real(box.a[i]) imag(box.a[i])
      # end
    else 
      # PROFILE, can track if already 0 so don't need to do this perhaps, but that would add branching
      box.a = zeros(P+1)
    end
  end
end


""" M2M : Multipole to multipole"""

# TESTING PARAM:  P = 33, N = 1000
# UNOPTIMIZED:    2.692 ms (1 allocation: 672 bytes)
# OPTIMIZED:      111.956 μs (1 allocation: 672 bytes)
function M2M!(quadtree::Quadtree, binomial_table::Array{Int64, 2})
  # propogate up multipole expansions from the leaves
  # to the highest boxes at depth 2
  depth_offsets::Array{Int, 1} = getDepthOffsets(quadtree)

  # PREALLOCATE all tmps used in multipole computation
  p = length(quadtree.tree[1].b)
  powers = Array{ComplexF64, 1}(undef, p+2)

  for depth in quadtree.tree_depth-1:-1:1
    for idx in 1:4^depth
      global_idx::Int = depth_offsets[depth] + idx
      parent_box::Box = quadtree.tree[global_idx]
      # Zero out as will be adding on contributions from children and not previous timesteps
      # QUESTION: when to zero
      parent_box.a .= zero(parent_box.a[1])
      # QUESTION: Faster to evaluate all children at once as longer array?
      
      # OPTIMIZED power operations and binomial lookup
      for child_idx in parent_box.children_idxs
        child_global_idx::Int = depth_offsets[depth+1] + child_idx
        child_box::Box = quadtree.tree[child_global_idx]
        parent_box.a[1] += child_box.a[1]

        diff = child_box.center - parent_box.center
        powers[1] = 1.0
        powers[2] = diff

        for l in 1:P
          i = l + 1
          parent_box.a[i] -= 1/l*child_box.a[1]*powers[i]
          powers[i+1] = powers[i]*diff
          for k in 1:l
            j = k + 1
            # PROFILED: binomial is expensive! Use a lookup table as only using small k (30 or 50 and under)!
            parent_box.a[i] += binomial_table[k, l]*child_box.a[j]*powers[l-k+1]
          end
        end
      end

      # UNOPTIMIZED
      #for child_idx in parent_box.children_idxs
      #  child_global_idx::Int = depth_offsets[depth+1] + child_idx
      #  child_box::Box = quadtree.tree[child_global_idx]
      #  parent_box.a[1] += child_box.a[1]
      #  for l in 1:P
      #    i = l + 1
      #    parent_box.a[i] -= 1/l*child_box.a[1]*(child_box.center - parent_box.center).^l
      #    for k in 1:l
      #      j = k + 1
      #      parent_box.a[i] += binomial(l-1, k-1)*child_box.a[j]*(child_box.center - parent_box.center).^(l-k)
      #    end
      #  end
      #end
    end
  end
end



""" Downward pass multipole stages """


""" M2L : Multipole to Local """

# TESTING PARAM:  P = 33, N = 1000
# UNOPTIMIZED:    164.646 ms (44521 allocations: 26.49 MiB)
# OPTIMIZED:      
function M2L!(quadtree::Quadtree, binomial_table::Array{Int64, 2})
  # Add contribution of boxes in interaction list to expansion of potential of each box
  depth_offsets::Array{Int, 1} = getDepthOffsets(quadtree)

  # PREALLOCATE all tmps used in multipole computation
  #Ps = [x for x in 1:P]
  #Ps_idxs = 2:P+1
  csa::Array{ComplexF64, 1} = Array{ComplexF64, 1}(undef, P)
  powers = Array{ComplexF64, 1}(undef, P+1)
  
  # propogate all the way to the leaf level (inclusive)
  for depth in 2:quadtree.tree_depth
    for global_idx in depth_offsets[depth]+1:depth_offsets[depth]+4^depth
      box::Box = quadtree.tree[global_idx]
      box.b .= zero(box.b[1])
      for interaction_idx in box.interaction_idxs
        # interacting box
        inter_box::Box = quadtree.tree[depth_offsets[depth] + interaction_idx]
        # common sub array 
        # skip 0th term as computing force, not potential
        
        # UNOPTIMIZED
        #csa = (-1).^Ps.*inter_box.a[Ps_idxs]./((inter_box.center - box.center).^Ps)
        #for l in 1:P
        #  # first term
        #  box.b[l] += -inter_box.a[1]/(l*(inter_box.center - box.center)^l)
        #  # summation term
        #  box.b[l] += 1/(inter_box.center - box.center)^l * sum((binomial.(Ps.+l.-1, Ps.-1).*csa))
        #end
        
        # OPTIMIZED
        diff = inter_box.center - box.center
        powers[1] = 1/diff
        sign = -1
        for i in 1:P
          csa[i] = sign*inter_box.a[i+1] * powers[i]
          powers[i+1] = powers[i]/diff
          sign *= -1
        end
        @inbounds for l in 1:P
          # first term
          box.b[l] += -1/l*inter_box.a[1]*powers[l]
          # summation term
          for k in 1:P
            box.b[l] += powers[l] * (binomial_table[k, k+l] * csa[k])
          end
        end
      end
    end
  end
end


""" L2L : Local to Local """


function L2L!(quadtree::Quadtree)
  # propogate down information to children
  # do for every level above leaf level
  depth_offsets::Array{Int, 1} = getDepthOffsets(quadtree)
  for depth in 2:quadtree.tree_depth-1
    for global_idx in depth_offsets[depth]+1:depth_offsets[depth]+4^depth
      parent_box::Box = quadtree.tree[global_idx]
      #println()
      for child_idx in parent_box.children_idxs
        child_global_idx::Int = depth_offsets[depth+1] + child_idx
        child_box::Box = quadtree.tree[child_global_idx]
        # NOTE, may want to vectorize this double loop if possible
        for l in 1:P
          for k in l:P
            child_box.b[l] += parent_box.b[k]*binomial(k,l)*(child_box.center - parent_box.center)^(k-l);
          end
        end
        # DEBUG
        #if depth == 2
        #  @printf "center: %f + %fi, " real(child_box.center) imag(child_box.center)
        #  println(child_box.b)
        #end
      end
    end
  end
end


""" L2P : Local to Particle """


function L2P!(quadtree::Quadtree, points, ω_p::Array{ComplexF64, 1})
  # Pass to particles the local information
  # multiply out all of the coefficients for the expansion
  leaf_offset::Int = getOffsetOfDepth(quadtree, quadtree.tree_depth)
  for global_idx in leaf_offset+1:length(quadtree.tree)
    box::Box = quadtree.tree[global_idx]
    if (boxHasPoints(box))
      # can potentially reorder looping when thinking about vectorization
      # this should be good vectorization but double check
      #relevant_potentials = @view potentials[box.start_idx:box.final_idx]
      #relevant_points = @view points[box.start_idx:box.final_idx]
      #for k in 0:P
      #  # likely want to perform an outer product here and sum along one dimension
      #  #relevant_potentials += sum(((relevant_points .- box.center).^k .* box.a))
      #    for i = Box(ibox_global).particlelist
      for idx in box.start_idx:box.final_idx
        # for forces don't use the first term
        ω_p[idx] = 0
        for k = 1:P
          # multiply by k because taking derivative
          ω_p[idx] += k*box.b[k]*(points[idx] - box.center)^(k-1);
        end
      end        
    end
  end
  #println(potentials)
end


""" NNC : Near Neighbor Contribution """


function NNC!(quadtree::Quadtree, points, masses::Array{Float64, 1}, ω_p::Array{ComplexF64, 1})
  # Final step in computation
  # The multipole expansions have been passed up the tree from sources
  # The low rank procedure applied to well separated 
  # Now must do O(N^2) computation for close points but since uniform distribution is assumed
  # then this is a constant computation
  
  # debugging copy
  #o = copy(potentials)
  #println("before")
  #println(potentials)

  leaf_offset::Int = getOffsetOfDepth(quadtree, quadtree.tree_depth)
  for global_idx in leaf_offset+1:length(quadtree.tree)
    box::Box = quadtree.tree[global_idx]
    if (boxHasPoints(box))
      # do not want to construct large matrices out of memory concerns
      # though if the matrix construction is not a huge penalty may want to look into because of vectorization
      # remebering julia is column-major
      relevant_points = @view points[box.start_idx:box.final_idx]
      relevant_ω_p = @view ω_p[box.start_idx:box.final_idx]
      relevant_masses = @view masses[box.start_idx:box.final_idx]

      
      # contribution from within same box
      # zero out computation with of body with itself as log(0) = -Inf
      # this is also 0 if only one point in box
      #kernel_mtx::Array{ComplexF64, 2} = log.(transpose(relevant_points) .- relevant_points)
      kernel_mtx::Array{ComplexF64, 2} = 1 ./ (transpose(relevant_points) .- relevant_points .+ 1e-32)
      foreach(i -> kernel_mtx[i, i] = zero(kernel_mtx[1, 1]), 1:length(relevant_points))
      # can't do simple dot product unfortunately
      relevant_ω_p .+= vec(sum(kernel_mtx.*relevant_masses, dims=1))
      
      # contribution from neighbor boxes
      for neighbor_idx in box.neighbor_idxs
        neighbor_box::Box = quadtree.tree[leaf_offset + neighbor_idx]
         # evaluate whether branching might be a bigger hit to performance
         # when assuming uniform distrubition should almost always have points
         # but also want to avoid zero dimensional arrays
         if boxHasPoints(neighbor_box)
           # construct matrix from set of points in box 
           # CHECK 
           neighbor_points = @view points[neighbor_box.start_idx:neighbor_box.final_idx]
           neighbor_masses = @view masses[neighbor_box.start_idx:neighbor_box.final_idx]
           #relevant_ω .+= vec(sum(neighbor_masses.*log.(transpose(relevant_points) .- neighbor_points), dims=1))
           relevant_ω_p .+= vec(sum(neighbor_masses .* (1 ./ (transpose(relevant_points) .- neighbor_points)), dims=1))
         end
      end

    end
  end 
  #println("NNC contribution")
  #println(potentials - o)
  #println("after")
  #println(potentials)
end

