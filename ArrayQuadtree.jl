"""

Instead of mainly relying on recursive data structure,
store each level of quadtree in its own array.

Note: for simplicity not using any notion of depth 0 (one big node), lowest depth is 1 (space is split into 4 boxes)
      minimum depth required is 2

"""

include("FourColorSort.jl")


using Plots

# NOTE:
# may want to try making Box immutable and keeping all mutable portions
# in some separate location
# could potentially also store in a DFS ordering!! may be much better for locality
# DFS ordering seems like a very good idea considering the four color sort!

mutable struct Box

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
  # cannot use this if keeping quadtree immutable and do not want to rebuild each timestep
  start_idx::Int
  final_idx::Int

  # source (notation of f is used for source)
  f::Float64
  # potentual (notation u is used for potential)
  u::Float64

  Box(depth::Int, idx::Int, center::ComplexF64) = new(depth, idx, center, findParentIdx(depth, idx), findChildrenIdxs(depth, idx), findNeighborIdxs(depth, idx), findInteractionIdxs(depth, idx), 0, 0, 0, 0) 

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

function displayQuadtreeBoxes(quadtree::Array{Box, 1}, max_depth::Int, side_length=1.0)
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

# shows leaf boxes and points with mass
function displayQuadtree(plot, tree::Array{Box, 1}, points::Array{ComplexF64, 1}, masses::Array{Float64, 1}, tree_depth::Int, side_length=1.0)
  bsl = side_length / 2^tree_depth # box side length
  hbsl = bsl / 2

  # only plot at leaf level
  leaf_offset::Int = length(tree) - 4^tree_depth + 1

  for global_idx in leaf_offset:length(tree) 
    box::Box = tree[global_idx]
    center::ComplexF64 = box.center
    plot!(plot, Shape(real(center) .+ [-hbsl, hbsl, hbsl, -hbsl], imag(center) .+ [-hbsl, -hbsl, hbsl, hbsl]), opacity=0.5, leg = false)
    # NOTE: try to optimize out this check
    if box.final_idx >= box.start_idx # if there is more than one particle 
      box_points = @view points[box.start_idx:box.final_idx]
      x = real.(box_points)
      y = imag.(box_points)
      box_masses = @view masses[box.start_idx:box.final_idx]
      scatter!(plot, (x, y), markersize=10*box_masses, legend=false, label="")
    end
  end
end

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

# NOTE: can optimize slighlty
# do not need to propagate up past depth 2 I believe
function propagateFUp(quadtree::Array{Box, 1}, masses::Array{Float64, 1}, max_depth::Int)
  depth_offsets::Array{Int, 1} = getDepthOffsets(max_depth)
  leaf_offset::Int = last(depth_offsets) + 1
  # sum all masses at leaf boxes
  for global_idx in leaf_offset:length(quadtree)
    fsum::Float64 = 0
    @simd for i in quadtree[global_idx].start_idx:quadtree[global_idx].final_idx
      @inbounds fsum += masses[i]
    end
    quadtree[global_idx].f = fsum
  end
  # propogate up the quadtree from smallest boxes to largest
  for depth in max_depth-1:-1:1
    for idx in 1:4^depth
      global_idx::Int = depth_offsets[depth] + idx
      fsum::Float64 = 0
      for child_idx in quadtree[global_idx].children_idxs
        child_global_idx::Int = depth_offsets[depth+1] + child_idx
        fsum += quadtree[child_global_idx].f
      end
      quadtree[global_idx].f = fsum
    end
  end

end

# NOTE: may make sense to be check ordering of interaction list for good locality
# Can also potentially take advantage of symmetry?
function computeBoxPotentials(quadtree::Array{Box, 1}, tree_depth::Int, kernelFunction::Function)
  depth_offsets::Array{Int, 1} = getDepthOffsets(max_depth)
  for depth in 2:tree_depth
    for global_idx in depth_offsets[depth]+1:depth_offsets[depth+1]
      box_b::Box = quadtree[global_idx]
      for interacting_idx in box_a.interaction_idxs
        box_a = quadtree[depth_offsets[depth] + interacting_idx]
        box_a.u += kernelFunction(box_a.center, box_b.center)*box_b.f
      end
    end
  end
end

function propagateDownPotential()

end

function computeNeighborPotentialContribution()

end
