"""
Construct a quadtree using the four color sort.
For clarity keeping as mutable struct (may have performance repurcussions)
In order to have a small constant number of particles in each leaf box will
assume near-uniform distribution then calculate necessary depth to achieve
expected constant number of particles in leaves.
"""

include("FourColorSort.jl")


# realize this is obscure, first draft notation
# red     =>  first:a-1
# orange  =>  a:c 
# green   =>  c+1:d
# blue    =>  d+1:last
mutable struct Box
  depth::Int                  # useful for verification
  potential::Float64
  data::Array{ComplexF64, 1}
  parent::Union{Box, Nothing}
  #top_left_child::Box         # potentially convert children to array
  #top_right_child::Box
  #bot_left_child::Box
  #bot_right_child::Box
  children::Union{Array{Box, 1}, Nothing}
  Box(depth::Int, potential::Float64, data::Array{ComplexF64, 1}) = new(depth, potential, data)
  Box(depth::Int, potential::Float64, data::Array{ComplexF64, 1}, parent::Box) = new(depth, potential, data, parent)
end


""" Construct the quadtree given an array of complex numbers """
function buildQuadtree!(array, parent, center, dim, first, last, depth)

  if depth == 0 return Nothing end

  a, c, d = fourColorSort!(array, center)
  top_left = @view array[first:a-1]
  top_right = @view array[a:c]
  bot_left = @view array[c+1:d]
  bot_right = @view array[d+1:last]
  half_dim = dim/2
  box = Box(depth - 1, 0.0, parent)
  box.children[1] = buildQuadtree!(top_left, center - half_dim + half_dim*1im, half_dim, first, a-1, depth-1)
  box.children[2] = buildQuadtree!(top_left, center + half_dim + half_dim*1im, half_dim, a, c, depth-1)
  box.children[3] = buildQuadtree!(top_left, center - half_dim - half_dim*1im, half_dim, c+1, d, depth-1)
  box.children[4] = buildQuadtree!(top_left, center + half_dim - half_dim*1im, half_dim, d+1, last, depth-1)

  return box

end
