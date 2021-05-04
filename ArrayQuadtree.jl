"""

Instead of mainly relying on recursive data structure,
store each level of quadtree in its own array.

Note: for simplicity not using any notion of depth 0 (one big node), lowest depth is 1 (space is split into 4 boxes)

"""

using Plots


struct Box

  depth::Int
  idx::Int

  # center of box
  center::ComplexF64

  parent_idx::Int
  children_idxs::Array{Int, 1}    # switch to static arrays, faster
                                  # do not store for leaf nodes
  neighbor_idxs::Array{Int, 1}    # perhaps this is only stored for leaf nodes
  interaction_idxs::Array{Int, 1} # aka the interaction list

  # indices into point and mass arrays
  start_idx::Int
  final_idx::Int

  Box(depth::Int, idx::Int, center::ComplexF64) = new(depth, idx, center, findParentIdx(depth, idx), findChildrenIdxs(depth, idx), findNeighborIdxs(depth, idx), findInteractionIdxs(depth, idx), 0, 0) 
end


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

function buildQuadtree(max_depth::Int, side_length=1.0)
  @assert max_depth >= 1
  num_boxes_in_tree::Int = sum([4^depth for depth in 1:max_depth])
  tree::Array{Box, 1} = Array{Box, 1}(undef, num_boxes_in_tree)
  depth_offsets = vcat([0], [sum([4^i for i in 1:depth]) for depth in 1:max_depth-1])
  for depth in 1:max_depth 
    for idx in 1:4^depth
      box_center::ComplexF64 = getBoxCenter(depth, idx, side_length)
      tree[depth_offsets[depth] + idx] = Box(depth, idx, box_center)
    end
  end

  return tree

end

function displayQuadtree(quadtree::Array{Box, 1}, max_depth::Int, side_length=1.0)
  gr(size=(1000, 1000))
  p = plot()
  bsl = side_length / 2^max_depth # box side length
  hbsl = bsl / 2

  # only plot at leaf level
  leaf_offset::Int = length(quadtree) - 4^max_depth + 1

  for global_idx in leaf_offset:length(quadtree) 
    center::ComplexF64 = quadtree[global_idx].center
    plot!(p, Shape(real(center) .+ [-hbsl, hbsl, hbsl, -hbsl], imag(center) .+ [-hbsl, -hbsl, hbsl, hbsl]), opacity=0.5, leg = false)
  end
  gui()

end
