"""

Multipole math
Assuming a log(x - y) potential function where x and y are complex (and multiplied by masses as well)

Upward: P2M, M2M
Downward: M2L, L2L, ??, neighbor contribution

"""

include("ArrayQuadtreeV2.jl")


# helper constants
const Ps = [x for x in 1:P]
const Ps_idxs = [x for x in 2:P+1]


"""

Upward pass functions

"""


""" P2M : Particle to Multipole """

function P2M(quadtree::Array{Box, 1}, points::Array{ComplexF64, 1}, masses::Array{Float64, 1}, tree_depth::Int)
  leaf_offset::Int = getOffsetOfDepth(tree_depth)
  # determine multipole expansion at all leaf boxes
  for global_idx in leaf_offset+1:length(quadtree)
    box::Box = quadtree[global_idx]
    # builtin sum uses cumsum so using the builtin for better error
    if (boxHasPoints(box))
      relevant_points = @view points[box.start_idx:box.final_idx]
      relevant_masses = @view masses[box.start_idx:box.final_idx]
      box.a[1] = sum(relevant_masses)
      #println(global_idx - leaf_offset)
      #@printf "center: %f + %fi\n" real(box.center) imag(box.center)
      #@printf "%f + %fi, " real(box.a[1]) imag(box.a[1])
      for i in 2:P+1
        # subtract box center to form multipole expansion about box center
        k = i - 1
        box.a[i] = -1/k*sum(((relevant_points .- box.center).^k).*relevant_masses)
        #@printf "%f + %fi, " real(box.a[i]) imag(box.a[i])
      end
      #@printf "\n"
    else 
      box.a = zeros(P+1)
    end
  end
end


""" M2M : Multipole to multipole"""

function M2M(quadtree::Array{Box, 1}, tree_depth::Int)
  # propogate up multipole expansions from the leaves
  # to the highest boxes at depth 2
  depth_offsets::Array{Int, 1} = getDepthOffsets(tree_depth)
  for depth in tree_depth-1:-1:1
    for idx in 1:4^depth
      global_idx::Int = depth_offsets[depth] + idx
      parent_box::Box = quadtree[global_idx]
      # zero out as will be adding on contributions from children
      # do not want data from previous timestep
      # should not need an extra temporary array for new contents
      # of a as processing one level at a time and always looking at one level deeper
      parent_box.a .= zero(parent_box.a[1])
      for child_idx in parent_box.children_idxs
        child_global_idx::Int = depth_offsets[depth+1] + child_idx
        child_box::Box = quadtree[child_global_idx]
        parent_box.a[1] += child_box.a[1]
        for l in 1:P
          i = l + 1
          parent_box.a[i] -= 1/l*child_box.a[1]*(child_box.center - parent_box.center).^l
          for k in 1:l
            j = k + 1
            parent_box.a[i] += binomial(l-1, k-1)*child_box.a[j]*(child_box.center - parent_box.center).^(l-k)
          end
        end
      end
      # DEBUG
      #if depth == 2
      #  @printf "center: %f + %fi\n" real(parent_box.center) imag(parent_box.center)
      #  println(parent_box.a)
      #end
    end
  end
end



""" Downward pass multipole stages """


""" M2L : Multipole to Local """


function M2L(quadtree::Array{Box, 1}, tree_depth::Int, start_depth::Int, end_depth::Int)
  # Add contribution of boxes in interaction list to multipole expansion of each box
  # Need to use separate array b as cannot update a mid computation as that would affect
  # later box interaction computations
  depth_offsets::Array{Int, 1} = getDepthOffsets(tree_depth)
  Ps = [x for x in 1:P]
  Ps_idxs = 2:P+1
  # propogate all the way to the leaf level (inclusive)
  for depth in start_depth:end_depth
    for global_idx in depth_offsets[depth]+1:depth_offsets[depth]+4^depth
      box::Box = quadtree[global_idx]
      box.b .= zero(box.b[1])
      for interaction_idx in box.interaction_idxs
        # interacting box
        inter_box::Box = quadtree[depth_offsets[depth] + interaction_idx]
        # b0 (1-indexing)
        box.b[1] += log(box.center - inter_box.center)*inter_box.a[1]
        # common sub array 
        csa = (-1).^Ps.*inter_box.a[Ps_idxs]./((inter_box.center - box.center).^Ps)
        box.b[1] += sum(csa)
        for l in 1:P
          j = l + 1
          # first term
          box.b[j] += -inter_box.a[1]/(l*(inter_box.center - box.center)^l)
          # summation term
          box.b[j] += 1/(inter_box.center - box.center)^l * sum((binomial.(Ps.+l.-1, Ps.-1).*csa))
        end
      end
      # DEBUG
      #if depth == 2
      #  @printf "center: %f + %fi\n" real(box.center) imag(box.center)
      #  println(box.b)
      #end
    end
  end
end


""" L2L : Local to Local """


function L2L(quadtree::Array{Box, 1}, tree_depth::Int)
  # propogate down information to children
  # do for every level above leaf level
  depth_offsets::Array{Int, 1} = getDepthOffsets(tree_depth)
  for depth in 2:tree_depth-1
    for global_idx in depth_offsets[depth]+1:depth_offsets[depth]+4^depth
      parent_box::Box = quadtree[global_idx]
      for child_idx in parent_box.children_idxs
        child_global_idx::Int = depth_offsets[depth+1] + child_idx
        child_box::Box = quadtree[child_global_idx]
        # can reuse a again as have completed addition of contributions from
        # interacting boxes
        # NOTE, may want to vectorize this double loop if possible
        for l in 0:P
          i = l+1
          for k in l:P
            j = k+1
            child_box.a[i] += parent_box.b[j]*binomial(k,l)*(child_box.center - parent_box.center)^(k-l);
          end
        end
      end
    end
  end
end


""" L2P : Local to Particle """


function L2P(quadtree::Array{Box, 1}, points::Array{ComplexF64, 1}, potentials::Array{ComplexF64, 1}, tree_depth::Int)
  # Pass to particles the local information
  # multiply out all of the coefficients for the expansion
  leaf_offset::Int = getOffsetOfDepth(tree_depth)
  for global_idx in leaf_offset+1:length(quadtree)
    box::Box = quadtree[global_idx]
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
        for k = 0:P
          j = k + 1
          potentials[idx] += box.a[j]*(points[idx] - box.center)^k;
        end
      end        
    end
  end
end


""" NNC : Near Neighbor Contribution """


function NNC(quadtree::Array{Box, 1}, points::Array{ComplexF64, 1}, masses::Array{Float64, 1}, potentials::Array{ComplexF64, 1}, tree_depth::Int)
  # Final step in computation
  # The multipole expansions have been passed up the tree from sources
  # The low rank procedure applied to well separated 
  # Now must do O(N^2) computation for close points but since uniform distribution is assumed
  # then this is a constant computation
  leaf_offset::Int = getOffsetOfDepth(tree_depth)
  for global_idx in leaf_offset+1:length(quadtree)
    box::Box = quadtree[global_idx]
    if (boxHasPoints(box))
      # do not want to construct large matrices out of memory concerns
      # though if the matrix construction is not a huge penalty may want to look into because of vectorization
      # remebering julia is column-major
      relevant_points = @view points[box.start_idx:box.final_idx]
      relevant_potentials = @view potentials[box.start_idx:box.final_idx]
      relevant_masses = @view masses[box.start_idx:box.final_idx]
      
      # contribution from within same box
      # zero out computation with of body with itself as log(0) = -Inf
      # this is also 0 if only one point in box
      kernel_mtx::Array{ComplexF64, 2} = log.(relevant_points .- relevant_points')
      foreach(i -> kernel_mtx[i, i] = zero(kernel_mtx[1, 1]), 1:length(relevant_points))
      # can't do simple dot product unfortunately
      relevant_potentials .+= vec(sum(kernel_mtx.*relevant_masses, dims=2))
      
      # contribution from neighbor boxes
      for neighbor_idx in box.neighbor_idxs
        neighbor_box::Box = quadtree[leaf_offset + neighbor_idx]
         # evaluate whether branching might be a bigger hit to performance
         # when assuming uniform distrubition should almost always have points
         # but also want to avoid zero dimensional arrays
         if boxHasPoints(neighbor_box)
           # construct matrix from set of points in box 
           # CHECK 
           neighbor_points = @view points[neighbor_box.start_idx:neighbor_box.final_idx]
           neighbor_masses = @view masses[neighbor_box.start_idx:neighbor_box.final_idx]
           relevant_potentials .+= vec(sum(log.(relevant_points' .- neighbor_points).*neighbor_masses, dims=2))
         end
      end

    end
  end 
end

