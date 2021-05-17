
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

# TESTING PARAM:    P = 33, N = 1000
# CURRENT RUNTIME:  3.783 ms (5490 allocations: 484.67 KiB)
function FMM!(quadtree::Quadtree, points, masses, ω_p, binomial_table, binomial_table_t, large_binomial_table_t, preallocated_mtx)
  # upward pass
  P2M!(quadtree, points, masses)
  M2M!(quadtree, binomial_table)
  # downward pass
  M2L!(quadtree, binomial_table, large_binomial_table_t)
  # pass in both binomial and transpose as haven't decided entirely which to use
  # may need the tranpose if want to use BigInt
  L2L!(quadtree, binomial_table, binomial_table_t)
  L2P!(quadtree, points, ω_p)
  NNC!(quadtree, points, masses, ω_p, preallocated_mtx)
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
  p = length(quadtree.tree[1].b)
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
      for k in 1:p
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
      box.a = zeros(p+1)
    end
  end
end


""" M2M : Multipole to multipole"""

# TESTING PARAM:  P = 33, N = 1000
# UNOPTIMIZED:    2.692 ms (1 allocation: 672 bytes)
# OPTIMIZED:      61.929 μs (1 allocation: 624 bytes)#
function M2M!(quadtree::Quadtree, binomial_table)
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

        for l in 1:p
          i = l + 1
          parent_box.a[i] -= 1/l*child_box.a[1]*powers[i]
          powers[i+1] = powers[i]*diff
          @inbounds for k in 1:l
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
# UNFORTUNATELY, using BigInt binomial table for high powers of P leads to a jump in runtime for 
# P = 34
# 17.214 ms (2 allocations: 1.27 KiB)
# after introducing another large transposed binomial table and switching loop order obtain:
# 12.449 ms (2 allocations: 1.27 KiB)
function M2L!(quadtree::Quadtree, binomial_table, large_binomial_table_t)
  # Add contribution of boxes in interaction list to expansion of potential of each box
  depth_offsets::Array{Int, 1} = getDepthOffsets(quadtree)
  p = length(quadtree.tree[1].b)

  # PREALLOCATE all tmps used in multipole computation
  csa::Array{ComplexF64, 1} = Array{ComplexF64, 1}(undef, p)
  powers = Array{ComplexF64, 1}(undef, p+1)
  #tmp = Array{ComplexF64, 1}(undef, p)
  
  if p > 33
    # propogate all the way to the leaf level (inclusive)
    for depth in 2:quadtree.tree_depth
      @inbounds for global_idx in depth_offsets[depth]+1:depth_offsets[depth]+4^depth
        box::Box = quadtree.tree[global_idx]
        box.b .= zero(box.b[1])
        @inbounds for interaction_idx in box.interaction_idxs
          inter_box::Box = quadtree.tree[depth_offsets[depth] + interaction_idx]
          diff = inter_box.center - box.center
          powers[1] = 1/diff
          sign = -1
          for i in 1:p
            csa[i] = sign*inter_box.a[i+1] * powers[i]
            powers[i+1] = powers[i]/diff
            sign *= -1
          end
          # FURTHER OPTIMIZED for large P
          @inbounds @simd for l in 1:p
            box.b[l] -= 1/l*inter_box.a[1]*powers[l]
          end
          #box.b .= inter_box.a[1] ./ tmp .* powers[1:p]         
          #tmp .= powers[1:p] .* csa[1:p]
          @inbounds for k in 1:p
            @inbounds for l in 1:p
               #box.b[l] += powers[l] * binomial_table[k, k+l] * csa[k]
               box.b[l] += powers[l] * large_binomial_table_t[k+l, k] * csa[k]
            end
            # Trying to vectorize like this allocates too much memory
            #box.b[1:p] .+= tmp .* large_binomial_table_t[k:k+p-1, k]
          end
        end
      end
    end
  else
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
          for i in 1:p
            csa[i] = sign*inter_box.a[i+1] * powers[i]
            powers[i+1] = powers[i]/diff
            sign *= -1
          end

          @inbounds for l in 1:p
            # first term
            box.b[l] -= 1/l*inter_box.a[1]*powers[l]
            # summation term
            tmp = 0.0 + 0.0im
            # CURIOUS, fastest without @simd and @inbounds, @avx was not working at the time
            for k in 1:p
               tmp += binomial_table[k, k+l] * csa[k]
            end
            box.b[l] += powers[l]*tmp
          end
        end
      end
    end
  end
end


""" L2L : Local to Local """


# TESTING PARAM:  P = 32, N = 1000
# UNOPTIMIZED:    2.074 ms (0 allocations: 0 bytes)  
# OPTIMIZED:      47.533 μs (1 allocation: 624 bytes)
function L2L!(quadtree::Quadtree, binomial_table, binomial_table_t::Array{Int64, 2})
  # propogate down information to children
  # do for every level above leaf level
  depth_offsets::Array{Int, 1} = getDepthOffsets(quadtree)
  p = length(quadtree.tree[1].b)

  # PREALLOCATE
  powers = Array{ComplexF64, 1}(undef, p)
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
        @inbounds @simd for i in 3:p
          powers[i] = powers[i-1]*diff
        end
        for l in 1:p
          @inbounds @simd for k in l:p
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
  p = length(quadtree.tree[1].b)

  for global_idx in leaf_offset+1:length(quadtree.tree)
    box::Box = quadtree.tree[global_idx]
    if (boxHasPoints(box))
      # Not using outerproducts or creating vectors as too memory intensive
      for idx in box.start_idx:box.final_idx
        ω_p[idx] = 0
        diff = points[idx] - box.center
        tmp = 1.0 
        for k = 1:p
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
# OPTIMIZED:      1.859 ms (5358 allocations: 440.08 KiB)
# Large gains in memory using preallocated matrix, allocations are from views
function NNC!(quadtree::Quadtree, points, masses::Array{Float64, 1}, ω_p::Array{ComplexF64, 1}, preallocated_mtx)
  # Final step in computation
  # The multipole expansions have been passed up the tree from sources.
  # The low rank procedure applied to well separated regions.
  # Now must do O(N^2) computation for close points but since uniform distribution is assumed,
  # then this is designed to be an O(1) computation.
  # Zero out computation with of body with itself as 1/(x-x) = Inf.

  # PREFETCH full preallocated mtx
  mtx = @view preallocated_mtx[1:size(preallocated_mtx)[1], 1:size(preallocated_mtx)[2]]

  leaf_offset::Int = getOffsetOfDepth(quadtree, quadtree.tree_depth)
  @inbounds for global_idx in leaf_offset+1:length(quadtree.tree)
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
      space_needed = box.final_idx - box.start_idx + 1
      mtx = @view preallocated_mtx[1:space_needed, 1:space_needed]
      mtx .= one(ComplexF64) ./ (transpose(relevant_points) .- relevant_points .+ S) .* relevant_masses
      #mtx .=  transpose(relevant_points) .- relevant_points .+ S
      #mtx .= one(ComplexF64) ./ mtx
      foreach(i -> mtx[i, i] = zero(mtx[1, 1]), 1:length(relevant_points))
      # # Can't do simple dot product unfortunately
      #mtx .*= relevant_masses
      relevant_ω_p .+= vec(sum(mtx, dims=1))
      # # Contribution from neighbor boxes
      for neighbor_idx in box.neighbor_idxs
        neighbor_box::Box = quadtree.tree[leaf_offset + neighbor_idx]
         if boxHasPoints(neighbor_box)
           # construct matrix from set of points in box 
           neighbor_points = @view points[neighbor_box.start_idx:neighbor_box.final_idx]
           neighbor_masses = @view masses[neighbor_box.start_idx:neighbor_box.final_idx]
           neighbor_space_needed = neighbor_box.final_idx - neighbor_box.start_idx + 1
           #println(@allocated mtx = @view preallocated_mtx[1:neighbor_space_needed, 1:space_needed])
           mtx = @view preallocated_mtx[1:neighbor_space_needed, 1:space_needed]
           mtx .= one(ComplexF64) ./ (transpose(relevant_points) .- neighbor_points) .* neighbor_masses
           #mtx .= one(ComplexF64) ./ (transpose(relevant_points) .- neighbor_points)
           #mtx .*= neighbor_masses
           relevant_ω_p .+= vec(sum(mtx, dims=1))
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

      # SLOW, turns out it is faster to construct matrices
      # for i in box.start_idx:box.final_idx
      #   for j in box.start_idx:box.final_idx
      #     if i == j
      #       continue
      #     else
      #       ω_p[i] += masses[j] / (points[i] - points[j] + 1e-32)
      #     end
      #   end
      # end
      # # Contribution from neighbor boxes
      # for neighbor_idx in box.neighbor_idxs
      #   neighbor_box::Box = quadtree.tree[leaf_offset + neighbor_idx]
      #    if boxHasPoints(neighbor_box)
      #      for i in box.start_idx:box.final_idx
      #        for j in neighbor_box.start_idx:2:neighbor_box.final_idx
      #          ω_p[i] += masses[j]* 1 / (points[i] - points[j] + 1e-32)
      #        end
      #      end
      #    end
      # end
       
    end
  end 
end

