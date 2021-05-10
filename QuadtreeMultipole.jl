"""

Multipole math
Assuming a log(x - y) potential function where x and y are complex (and multiplied by masses as well)

Upward: P2M, M2M
Downward: M2L, L2L, ??, neighbor contribution

"""

include("ArrayQuadtree.jl")

# number of additional coefficients of multipole expansion to keep
# total number of coefficients will be P + 1
const P = 50
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
      for i in 2:P+1
        # subtract box center to form multipole expansion about box center
        k = i - 1
        box.a[i] = -1/k*sum(((relevant_points .- box.center).^k).*relevant_masses)
      end
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
          i = l + i
          parent_box.a[i] -= 1/l*child_box.a[1]*(child_box.center - parent_box.center).^l
          for k in 1:l
            j = k + 1
            parent_box.a[i] += binomial(l, k)*child_box.a[j]*(child_box.center - parent_box.center).^(l-k)
          end
        end
      end
    end
  end
end



""" Downward pass multipole stages """


""" M2L : Multipole to Local """


function M2L(quadtree::Array{Box, 1}, tree_depth::Int)
  # Add contribution of boxes in interaction list to multipole expansion of each box
  # Need to use separate array b as cannot update a mid computation as that would affect
  # later box interaction computations
  depth_offsets::Array{Int, 1} = getDepthOffsets(tree_depth)
  Ps = [x for x in 1:P]
  Ps_idxs = 2:P+1
  # propogate all the way to the leaf level (inclusive)
  for depth in 2:tree_depth
    for global_idx in depth_offsets[depth]+1:depth_offsets[depth+1]
      box::Box = quadtree[global_idx]
      box.b .= zero(box.b[1])
      for interaction_idx in box.interaction_idxs
        # interacting box
        inter_box::Box = depth_offsets[depth] + interaction_idx
        # b0 (1-indexing)
        box.b[1] += log(box.center - inter_box.center)*inter_box.a[1]
        # common sub array 
        csa = (-1).^Ps.*inter_box.a[Ps_idxs]./((inter_box.center - box.center).^Ps)
        box.b[1] += sum(cse)
        for l in 1:P
          j = l + 1
          # first term
          box.b[j] += -1/l*inter_box.a[1]*(inter_box.center - box.center)^l 
          # summation term
          box.b[j] += sum((binomial.(Ps+l, l+1)./((inter_box.center - box.center).^Ps).*cse))
        end
      end
    end
  end
end


""" L2L : Local to Local """


function L2L(quadtree::Array{Box, 1}, tree_depth::Int)
  # propogate down information to children
  # do for every level above leaf level
  depth_offsets::Array{Int, 1} = getDepthOffsets(tree_depth)
  for depth in 2:tree_depth-1
    for global_idx in depth_offsets[depth]+1:depth_offsets[depth+1]
      parent_box::Box = quadtree[global_idx]
      for child_idx in box.children_idxs
        child_global_idx::Int = depth_offsets[depth+1] + child_idx
        child_box::Box = quadtree[child_global_idx]
        # can reuse a again as have completed addition of contributions from
        # interacting boxes
        # NOTE, may want to vectorize this double loop if possible
        for l in 0:P
          i = l+1
          for k in l:P
            j = k+1                         # binomial(k+1, l+1) ??
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
      relevant_potentials = @view potentials[box.start_idx:box.final_idx]
      relevant_points = @view points[box.start_idx:box.final_idx]
      # can potentially reorder looping when thinking about vectorization
      # this should be good vectorization but double check
      for k in 0:P
        relevant_potentials .+= ((relevant_points .- box.center).^k .* box.a)
      end
    end
  end
end


""" NNC : Near Neighbor Contribution """


function NNC(quadtree::Array{Box, 1}, points::Array{ComplexF64, 1}, potentials::Array{ComplexF64, 1}, tree_depth::Int)
  # Final step in computation
  # The multipole expansions have been passed up the tree from sources
  # The low rank procedure applied to well separated 
  # Now must do O(N^2) computation for close points but since uniform distribution is assumed
  # then this is a constant computation
  leaf_offset::Int = getOffsetOfDepth(tree_depth)
  for global_idx in leaf_offset+1:length(quadtree)
    box::Box = quadtree[global_idx]
    # do not want to construct large matrices out of memory concerns
    for neighbor_idx in box.neighbor_idxs
       neighbor_box::Box = leaf_offset + neighbor_idx
    end
  end 
end

