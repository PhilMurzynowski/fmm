
""" Definitions and helper functions for a Quadtree """


include("FourColorSort.jl")


""" Basic unit of a Quadtree """

mutable struct Box
  # multipole expansion information
  # outer expansion and inner expansions
  outer_exp::Array{ComplexF64, 1}
  inner_exp::Array{ComplexF64, 1}
  # indices into point and mass arrays
  # cannot use this if keeping quadtree immutable and do not want to rebuild each timestep
  start_idx::Int
  final_idx::Int
  # center of box
  center::ComplexF64
  # all organizational info is kept in the footer
  depth::Int
  idx::Int
  children_idxs::Array{Int, 1}    # switch to static arrays, faster, # do not store for leaf nodes
  neighbor_idxs::Array{Int, 1}    # perhaps this is only stored for leaf nodes
  parent_idx::Int
  interaction_idxs::Array{Int, 1} # aka the interaction list
  # initialize start, final as 0, -1. -1 chosen so final < start
  Box(depth::Int, idx::Int, center::ComplexF64, p::Int) = new(zeros(ComplexF64, p+1), zeros(ComplexF64, p), 0, -1, center, depth, idx, findChildrenIdxs(idx, depth), findNeighborIdxs(idx, depth), findParentIdx(idx, depth), findInteractionIdxs(idx, depth)) 
end

struct Quadtree
  # header information
  tree_depth::Int
  depth_offsets::Array{Int, 1}
  # actual tree represented as an array
  tree::Array{Box, 1}
end


""" Function to build a quadtree of specific depth and size """


function buildQuadtree(tree_depth::Int, p::Int, side_length=1.0)
  @assert tree_depth >= 3
  # build as one array to potentially increase cache locality
  num_boxes_in_tree::Int = sum([4^depth for depth in 1:tree_depth])
  tree::Array{Box, 1} = Array{Box, 1}(undef, num_boxes_in_tree)
  depth_offsets = vcat([0], [sum([4^i for i in 1:depth]) for depth in 1:tree_depth-1])
  for depth in 1:tree_depth 
    for idx in 1:4^depth
      box_center::ComplexF64 = getBoxCenter(depth, idx, side_length)
      tree[depth_offsets[depth] + idx] = Box(depth, idx, box_center, p)
    end
  end
  return Quadtree(tree_depth, depth_offsets, tree)
end


""" Functions to update associated points and masses for boxes in quadtree """


# NOTE: doing this as BFS could potentially have better locality, unless switch quadtree to bfs ordering
# currently jumping and jumping to higher depth offsets
# RUNTIME: this ends up being O(Nlog(N))
function colorSortQuadtreePointMasses(box::Box, quadtree::Quadtree, points, masses::Array{Float64, 1}, prev_points, depth::Int, first::Int, last::Int)

  if depth == quadtree.tree_depth
    box.start_idx = first
    box.final_idx = last
    return
  end

  a::Int, c::Int, d::Int = fourColorSort!(points, masses, prev_points, box.center, first, last)
  depth_offset::Int = getOffsetOfDepth(quadtree, depth+1)
  tl_child::Box = quadtree.tree[depth_offset + box.children_idxs[1]]
  bl_child::Box = quadtree.tree[depth_offset + box.children_idxs[2]]
  tr_child::Box = quadtree.tree[depth_offset + box.children_idxs[3]]
  br_child::Box = quadtree.tree[depth_offset + box.children_idxs[4]]
  colorSortQuadtreePointMasses(tl_child, quadtree, points, masses, prev_points, depth+1, first, a-1)
  colorSortQuadtreePointMasses(bl_child, quadtree, points, masses, prev_points, depth+1, a, c)
  colorSortQuadtreePointMasses(tr_child, quadtree, points, masses, prev_points, depth+1, c+1, d)
  colorSortQuadtreePointMasses(br_child, quadtree, points, masses, prev_points, depth+1, d+1, last)

end

function updateQuadtreePointMasses(quadtree::Quadtree, points, masses::Array{Float64, 1}, prev_points, side_length::Float64=1.0)
  tree_center::ComplexF64 = side_length/2 + side_length/2*1im
  first::Int = 1
  last::Int = length(points)
  a::Int, c::Int, d::Int = fourColorSort!(points, masses, prev_points, tree_center, first, last)
  tl_child::Box = quadtree.tree[1]
  bl_child::Box = quadtree.tree[2]
  tr_child::Box = quadtree.tree[3]
  br_child::Box = quadtree.tree[4]
  colorSortQuadtreePointMasses(tl_child, quadtree, points, masses, prev_points, 1, first, a-1)
  colorSortQuadtreePointMasses(bl_child, quadtree, points, masses, prev_points, 1, a, c)
  colorSortQuadtreePointMasses(tr_child, quadtree, points, masses, prev_points, 1, c+1, d)
  colorSortQuadtreePointMasses(br_child, quadtree, points, masses, prev_points, 1, d+1, last)
end


""" Functions used to navigate quadtree """


function getBoxCenter(depth, idx, side_length)
  # reminder that using column-major orering for idxs
  box_side_length = side_length / 2^depth
  half_box_side_length = box_side_length / 2
  x::Float64 = floor((idx - 1) / 2^depth) * box_side_length + half_box_side_length 
  y::Float64 = (2^depth - mod1(idx, 2^depth)) * box_side_length + half_box_side_length
  return x + y*1im
end

function findChildrenIdxs(idx::Int, depth::Int)
  # find index for top left child
  # other children are simple offsets from first child
  y_contribution::Int = mod(idx - 1, 2^depth) * 2 + 1 
  x_contribution::Int = floor((idx - 1) / 2^depth) * 2*2^(depth+1)
  tl_child::Int = x_contribution + y_contribution
  bl_child::Int = tl_child + 1
  tr_child::Int = tl_child + 2^(depth+1)
  br_child::Int = tr_child + 1
  return [tl_child, bl_child, tr_child, br_child]
end

function findNeighborIdxs(idx::Int, depth::Int)
  column_size::Int = 2^depth
  num_boxes::Int = 4^depth
  neighbors::Array{Int, 1} = Array{Int, 1}[]
  # 2 middle side neighbors
  # check if not in left column
  if (idx > column_size)
    push!(neighbors, idx - column_size)
  end
  # check if not in right column
  if (idx <= num_boxes - column_size)
    push!(neighbors, idx + column_size)
  end
  # top 3 neighbors
  # check if not in top row
  if (mod1(idx, column_size) != 1)
    push!(neighbors, idx - 1)
    # check if not in left column
    if (idx > column_size)
      push!(neighbors, idx - column_size - 1)
    end
    # check if not in right column
    if (idx <= num_boxes - column_size)
      push!(neighbors, idx + column_size - 1)
    end
  end
  # bottom three neighbors
  # check if not in bottom row
  if (mod1(idx, column_size) != column_size)
    push!(neighbors, idx + 1)
    # check if not in left column
    if (idx > column_size)
      push!(neighbors, idx - column_size + 1)
    end
    # check if not in right column
    if (idx <= num_boxes - column_size)
      push!(neighbors, idx + column_size + 1)
    end
  end
  return neighbors
end

function boxHasPoints(box::Box)
  return box.final_idx >= box.start_idx
end

function getDepthOffsets(quadtree::Quadtree)
  return quadtree.depth_offsets
end

function getOffsetOfDepth(quadtree::Quadtree, depth::Int)
  return quadtree.depth_offsets[depth]
end

function findParentIdx(idx::Int, depth::Int)
  y::Int = ceil(mod1(idx, 2^depth) / 2)
  x::Int = floor((idx-1) / 2^(depth+1))
  return x*2^(depth-1) + y
end

function findInteractionIdxs(idx::Int, depth::Int)
  if (depth == 1)
    return Array{Int, 1}[]
  end
  parent_idx::Int = findParentIdx(idx, depth)
  parent_neighbor_idxs::Array{Int, 1} = findNeighborIdxs(parent_idx, depth-1)
  interacting_idxs::Array{Int, 1} = vcat(findChildrenIdxs.(parent_neighbor_idxs, depth-1)...)
  return setdiff(interacting_idxs, findNeighborIdxs(idx, depth))
end
