
""" Multipole Math for the Fast Multipole Method 

In a 2D plane gravitational force takes the form 1 / (x - y), which follows
from a log(x - y) potential function where x and y are complex.

Upward pass components: P2M, M2M
Downwd pass components: M2L, L2L, NNC   """


using LoopVectorization


include("Quadtree.jl")

"""

Wrapper for all components


"""

function FMM!(quadtree::Quadtree, points, masses, ω_p, binomial_table, binomial_table_t)
  # upward pass
  P2M!(quadtree, points, masses)
  M2M!(quadtree, binomial_table)
  # downward pass
  M2L!(quadtree, binomial_table)
  # pass in both binomial and transpose as haven't decided entirely which to use
  # may need the tranpose if want to use BigInt
  L2L!(quadtree, binomial_table, binomial_table_t)
  L2P!(quadtree, points, ω_p)
  #@btime NNC!(quadtree, $points, $masses, $ω_p)
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
          # TRY @simd and @inbounds, especially here
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
# OPTIMIZED:      1.565 ms (2 allocations: 1.22 KiB)
function M2L!(quadtree::Quadtree, binomial_table::Array{Int64, 2})
  # Add contribution of boxes in interaction list to expansion of potential of each box
  depth_offsets::Array{Int, 1} = getDepthOffsets(quadtree)

  # PREALLOCATE all tmps used in multipole computation
  csa::Array{ComplexF64, 1} = Array{ComplexF64, 1}(undef, P)
  powers = Array{ComplexF64, 1}(undef, P+1)
  
  # propogate all the way to the leaf level (inclusive)
  for depth in 2:quadtree.tree_depth
    @inbounds for global_idx in depth_offsets[depth]+1:depth_offsets[depth]+4^depth
      box::Box = quadtree.tree[global_idx]
      box.b .= zero(box.b[1])
      @inbounds for interaction_idx in box.interaction_idxs
        # interacting box
        inter_box::Box = quadtree.tree[depth_offsets[depth] + interaction_idx]
        # common sub array 
        # skip 0th term as computing force, not potential
        
        # UNOPTIMIZED (main portion, also was not preallocating and using @inbounds previously)
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
          box.b[l] -= 1/l*inter_box.a[1]*powers[l]
          # summation term
          tmp = 0.0 + 0.0im
          # CURIOUS, fastest without @simd and @inbounds, @avx was not working at the time
          for k in 1:P
             tmp += binomial_table[k, k+l] * csa[k]
          end
          box.b[l] += powers[l]*tmp
        end
      end
    end
  end
end


""" L2L : Local to Local """


# TESTING PARAM:  P = 32, N = 1000
# UNOPTIMIZED:    2.074 ms (0 allocations: 0 bytes)  
# OPTIMIZED:      47.533 μs (1 allocation: 624 bytes)
function L2L!(quadtree::Quadtree, binomial_table::Array{Int64, 2}, binomial_table_t::Array{Int64, 2})
  # propogate down information to children
  # do for every level above leaf level
  depth_offsets::Array{Int, 1} = getDepthOffsets(quadtree)
  # PREALLOCATE
  powers = Array{ComplexF64, 1}(undef, P)
  powers[1] = 1.0

  for depth in 2:quadtree.tree_depth-1
    for global_idx in depth_offsets[depth]+1:depth_offsets[depth]+4^depth
      parent_box::Box = quadtree.tree[global_idx]
      for child_idx in parent_box.children_idxs
        child_global_idx::Int = depth_offsets[depth+1] + child_idx
        child_box::Box = quadtree.tree[child_global_idx]

        # OPTIMIZED
        diff = child_box.center - parent_box.center
        powers[2] = diff
        @inbounds @simd for i in 3:P
          powers[i] = powers[i-1]*diff
        end
        for l in 1:P
          @inbounds @simd for k in l:P
            # transposed binomial table for different access pattern
            # May be useful if start using BigInts
            # Normal access pattern:
            # child_box.b[l] += parent_box.b[k]*binomial_table[l+1,k+1]*powers[k-l+1];
            child_box.b[l] += parent_box.b[k]*binomial_table_t[k,l]*powers[k-l+1];
          end
        end

        #UNOPTIMIZED
        #for l in 1:P
        #  for k in l:P
        #    child_box.b[l] += parent_box.b[k]*binomial(k,l)*(child_box.center - parent_box.center)^(k-l);
        #  end
        #end

      end
    end
  end
end


""" L2P : Local to Particle """

# TESTING PARAM:  P = 32, N = 1000
# UNOPTIMIZED:    306.131 μs (0 allocations: 0 bytes) 
# OPTIMIZED:      86.411 μs (0 allocations: 0 bytes)
function L2P!(quadtree::Quadtree, points, ω_p::Array{ComplexF64, 1})
  # Pass to particles the local information
  # multiply out all of the coefficients for the expansion
  # for forces don't use the first term
  # multiply by k because taking derivative of potential expansion
  # NOTE: @inbounds and @simd didn't make any performance difference here
  
  leaf_offset::Int = getOffsetOfDepth(quadtree, quadtree.tree_depth)
  for global_idx in leaf_offset+1:length(quadtree.tree)
    box::Box = quadtree.tree[global_idx]
    if (boxHasPoints(box))
      # Not using outerproducts or creating vectors as too memory intensive
      for idx in box.start_idx:box.final_idx
        ω_p[idx] = 0
        diff = points[idx] - box.center
        tmp = 1.0 
        for k = 1:P
          ω_p[idx] += k*box.b[k]*tmp
          tmp *= diff
        end
      end        

      # VECTOR FORMAT (slower due to allocation, preallocating requires searching through boxes for max number of points in a box)
      # idxs = box.start_idx:box.final_idx
      # ω_p[idxs] .= zero(ComplexF64)
      # diff = points[idxs] .- box.center
      # tmp = ones(ComplexF64, length(diff)) 
      # @inbounds @simd for k = 1:P
      #   ω_p[idxs] .+= k*box.b[k].*tmp
      #   tmp .*= diff
      # end

      # UNOPTIMIZED
      #for idx in box.start_idx:box.final_idx
      #  ω_p[idx] = 0
      #  for k = 1:P
      #    ω_p[idx] += k*box.b[k]*(points[idx] - box.center)^(k-1);
      #  end
      #end        
    end
  end
end


""" NNC : Near Neighbor Contribution """

# TESTING PARAM:  P = 32, N = 1000
# UNOPTIMIZED:    2.159 ms (2000 allocations: 2.31 MiB)
# OPTIMIZED:      

function NNC!(quadtree::Quadtree, points, masses::Array{Float64, 1}, ω_p::Array{ComplexF64, 1})
  # Final step in computation
  # The multipole expansions have been passed up the tree from sources.
  # The low rank procedure applied to well separated regions.
  # Now must do O(N^2) computation for close points but since uniform distribution is assumed,
  # then this is designed to be an O(1) computation.
  # Zero out computation with of body with itself as 1/(x-x) = Inf.

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

      # OPTIMIZED
      # Contribution from within same box
      kernel_mtx::Array{ComplexF64, 2} = 1 ./ (transpose(relevant_points) .- relevant_points .+ 1e-32)
      foreach(i -> kernel_mtx[i, i] = zero(kernel_mtx[1, 1]), 1:length(relevant_points))
      # can't do simple dot product unfortunately
      relevant_ω_p .+= vec(sum(kernel_mtx.*relevant_masses, dims=1))
      # Contribution from neighbor boxes
      for neighbor_idx in box.neighbor_idxs
        neighbor_box::Box = quadtree.tree[leaf_offset + neighbor_idx]
         # Evaluate whether branching might be a bigger hit to performance
         # when assuming uniform distrubition should almost always have points
         # but also want to avoid zero dimensional arrays
         if boxHasPoints(neighbor_box)
           # construct matrix from set of points in box 
           neighbor_points = @view points[neighbor_box.start_idx:neighbor_box.final_idx]
           neighbor_masses = @view masses[neighbor_box.start_idx:neighbor_box.final_idx]
           relevant_ω_p .+= vec(sum(neighbor_masses .* (1 ./ (transpose(relevant_points) .- neighbor_points)), dims=1))
         end
      end

      # UNOPTIMIZED
      # Contribution from within same box
      # kernel_mtx::Array{ComplexF64, 2} = 1 ./ (transpose(relevant_points) .- relevant_points .+ 1e-32)
      # foreach(i -> kernel_mtx[i, i] = zero(kernel_mtx[1, 1]), 1:length(relevant_points))
      # # Can't do simple dot product unfortunately
      # relevant_ω_p .+= vec(sum(kernel_mtx.*relevant_masses, dims=1))
      # # Contribution from neighbor boxes
      # for neighbor_idx in box.neighbor_idxs
      #   neighbor_box::Box = quadtree.tree[leaf_offset + neighbor_idx]
      #    if boxHasPoints(neighbor_box)
      #      # construct matrix from set of points in box 
      #      neighbor_points = @view points[neighbor_box.start_idx:neighbor_box.final_idx]
      #      neighbor_masses = @view masses[neighbor_box.start_idx:neighbor_box.final_idx]
      #      relevant_ω_p .+= vec(sum(neighbor_masses .* (1 ./ (transpose(relevant_points) .- neighbor_points)), dims=1))
      #    end
      # end

    end
  end 
end

