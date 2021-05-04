"""

Instead of mainly relying on recursive data structure,
store each level of quadtree in its own array.

Note: for simplicity not using any notion of depth 0 (one big node), lowest depth is 1 (space is split into 4 boxes)

"""

struct Box

  depth::Int
  idx::Int
  parent_idx::Int
  children_idxs::Array{Int, 1}    # switch to static arrays, faster
  neighbor_idxs::Array{Int, 1}    # perhaps this is only stored for leaf nodes
  interaction_idxs::Array{Int, 1} # aka the interaction list

  Box(depth::Int, idx::Int) = new(depth, idx, findParentIdx(depth, idx), findChildrenIdxs(depth, idx), findNeighborIdxs(depth, idx), findInteractionIdxs(depth, idx)) 

end


function findParentIdx(depth::Int, idx::Int)
  y::Int = ceil(mod1(idx, 2^depth) / 2)
  x::Int = floor((idx-1) / 2^(depth+1))
  return x*2^(depth-1) + y
end

function findChildrenIdxs(depth::Int, idx::Int)
  # find index for top left child
  # other children are simple offsets from first child
  y_contribution::Int = mod(idx - 1, 2^depth) * 2
  x_contribution::Int = floor(idx / 2^depth) * 2^(depth + 1)
  tl_child::Int = x_contribution + y_contribution
  bl_child::Int = tl_child + 1
  tr_child::Int = tl_child + 2^(depth+1)
  br_child::Int = tr_child + 1
  return [tl_child, bl_child, tr_child, br_child]
end

function findNeighborIdxs(depth::Int, idx::Int)
  column_size::Int = 2^depth
  num_boxes::Int = 4^depth
  #neighbors::Array{Int, 1} = [idx - column_size - 1, idx - column_size, idx - column_size + 1, idx - 1, idx + 1, idx + column_size - 1, idx + column_size, idx + column_size + 1]
  neigbors::Array{Int, 1} = Array{Int, 1}[]

  # 2 middle side neighbors
  # check if not in left column
  if (idx > column_size)
    push!(neighbors, idx - column_size)
  end
  # check if not in right column
  if (idx <= num_boxes - column_size)
    push!(neighbors, idx + column_size1)
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
  interacting_idxs::Array{Int, 1} = vcat(findChildrenIdxs.(parent_neighbor_idxs))

  return setdiff(interacting_idxs, getNeighborIdxs(depth, int))

end
