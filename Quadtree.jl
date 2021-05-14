

""" Basic unit of a Quadtree """


mutable struct Box
  depth::Int
  idx::Int
  # center of box
  center::ComplexF64
  parent_idx::Int
  children_idxs::Array{Int, 1}    # switch to static arrays, faster, # do not store for leaf nodes
  neighbor_idxs::Array{Int, 1}    # perhaps this is only stored for leaf nodes
  interaction_idxs::Array{Int, 1} # aka the interaction list
  # indices into point and mass arrays
  # cannot use this if keeping quadtree immutable and do not want to rebuild each timestep
  start_idx::Int
  final_idx::Int
  # multipole expansion information
  a::Array{ComplexF64, 1}
  b::Array{ComplexF64, 1}
  # initialize start, final as 0, -1. -1 chosen so final < start
  Box(depth::Int, idx::Int, center::ComplexF64, p::Int) = new(depth, idx, center, findParentIdx(depth, idx), findChildrenIdxs(depth, idx), findNeighborIdxs(depth, idx), findInteractionIdxs(depth, idx), 0, -1, zeros(ComplexF64, p+1), zeros(ComplexF64, p+1)) 
end


""" Function to build a quadtree of specific depth and size """


function buildQuadtree(max_depth::Int, side_length=1.0)
  @assert max_depth >= 1
  # build as one array to potentially increase cache locality
  num_boxes_in_tree::Int = sum([4^depth for depth in 1:max_depth])
  tree::Array{Box, 1} = Array{Box, 1}(undef, num_boxes_in_tree)
  depth_offsets = getDepthOffsets(max_depth)
  for depth in 1:max_depth 
    for idx in 1:4^depth
      box_center::ComplexF64 = getBoxCenter(depth, idx, side_length)
      tree[depth_offsets[depth] + idx] = Box(depth, idx, box_center)
    end
  end
  return tree
end


""" Functions to update associated points and masses for boxes in quadtree """


# NOTE: doing this as BFS could potentially have better locality, unless switch quadtree to bfs ordering
# currently jumping and jumping to higher depth offsets
function colorSortQuadtreePointMasses(box::Box, quadtree::Array{Box, 1}, points::Array{ComplexF64, 1}, masses::Array{Float64, 1}, max_depth::Int, depth::Int, first::Int, last::Int)

  if depth == max_depth
    box.start_idx = first
    box.final_idx = last
    return
  end

  a::Int, c::Int, d::Int = fourColorSort!(points, masses, box.center, first, last)
  depth_offset::Int = getOffsetOfDepth(depth+1)
  tl_child::Box = quadtree[depth_offset + box.children_idxs[1]]
  bl_child::Box = quadtree[depth_offset + box.children_idxs[2]]
  tr_child::Box = quadtree[depth_offset + box.children_idxs[3]]
  br_child::Box = quadtree[depth_offset + box.children_idxs[4]]
  colorSortQuadtreePointMasses(tl_child, quadtree, points, masses, max_depth, depth+1, first, a-1)
  colorSortQuadtreePointMasses(bl_child, quadtree, points, masses, max_depth, depth+1, a, c)
  colorSortQuadtreePointMasses(tr_child, quadtree, points, masses, max_depth, depth+1, c+1, d)
  colorSortQuadtreePointMasses(br_child, quadtree, points, masses, max_depth, depth+1, d+1, last)

end

function updateQuadtreePointMasses(quadtree::Array{Box, 1}, points::Array{ComplexF64, 1}, masses::Array{Float64, 1}, tree_depth::Int, side_length::Float64=1.0)
  tree_center::ComplexF64 = side_length/2 + side_length/2*1im
  first::Int = 1
  last::Int = length(points)
  a::Int, c::Int, d::Int = fourColorSort!(points, masses, tree_center, first, last)
  tl_child::Box = quadtree[1]
  bl_child::Box = quadtree[2]
  tr_child::Box = quadtree[3]
  br_child::Box = quadtree[4]
  colorSortQuadtreePointMasses(tl_child, quadtree, points, masses, tree_depth, 1, first, a-1)
  colorSortQuadtreePointMasses(bl_child, quadtree, points, masses, tree_depth, 1, a, c)
  colorSortQuadtreePointMasses(tr_child, quadtree, points, masses, tree_depth, 1, c+1, d)
  colorSortQuadtreePointMasses(br_child, quadtree, points, masses, tree_depth, 1, d+1, last)
end


""" Functions used to navigate quadtree """


function findParentIdx(depth::Int, idx::Int)
  y::Int = ceil(mod1(idx, 2^depth) / 2)
  x::Int = floor((idx-1) / 2^(depth+1))
  return x*2^(depth-1) + y
end

function findChildrenIdxs(depth::Int, idx::Int)
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

function findNeighborIdxs(depth::Int, idx::Int)
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

function findInteractionIdxs(depth::Int, idx::Int)

  if (depth == 1)
    return Array{Int, 1}[]
  end

  parent_idx::Int = findParentIdx(depth, idx)
  parent_neighbor_idxs::Array{Int, 1} = findNeighborIdxs(depth-1, parent_idx)
  interacting_idxs::Array{Int, 1} = vcat(findChildrenIdxs.(depth-1, parent_neighbor_idxs)...)

  return setdiff(interacting_idxs, findNeighborIdxs(depth, idx))

end

function getBoxCenter(depth, idx, side_length)
  # reminder that using column-major orering for idxs
  box_side_length = side_length / 2^depth
  half_box_side_length = box_side_length / 2
  x::Float64 = floor((idx - 1) / 2^depth) * box_side_length + half_box_side_length 
  y::Float64 = (2^depth - mod1(idx, 2^depth)) * box_side_length + half_box_side_length
  return x + y*1im
end

function boxHasPoints(box::Box)
  return box.final_idx >= box.start_idx
end

function getDepthOffsets(max_depth::Int)
  return vcat([0], [sum([4^i for i in 1:depth]) for depth in 1:max_depth-1])
end

function getOffsetOfDepth(depth::Int)
  if depth == 1
    return 0
  else 
    return sum([4^i for i in 1:depth-1])
  end
end

