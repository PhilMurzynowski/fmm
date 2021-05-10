"""

Test multipole code

"""


#include("ArrayQuadtreeV2.jl") #included by below file
include("QuadtreeMultipole.jl")


function testFullMultipoleComputation()
  # an end to end test for one time step

  # specify tree parameters and construct tree
  tree_depth::Int = 4
  side_length::Float64 = 1.0
  quadtree::Array{Box, 1} = buildQuadtree(tree_depth, side_length)
  # generate points and masses
  num_bodies::Int = 10
  points::Array{ComplexF64, 1} = rand(ComplexF64, num_bodies) * side_length
  masses::Array{Float64, 1} = rand(Float64, num_bodies)
  # associate quadtree with point masses
  updateQuadtreePointMasses(quadtree, points, masses, tree_depth, side_length)
  # upward pass
  P2M(quadtree, points, masses, tree_depth)
  M2M(quadtree, tree_depth)
  # downward pass
  M2L(quadtree, tree_depth)
  L2L(quadtree, tree_depth)
  # array to hold computed potentials
  potentials::Array{ComplexF64, 1} = Array{ComplexF64, 1}(undef, num_bodies)
  potentials .= zero(potentials[1])
  # remainder of downward pass
  L2P(quadtree, points, potentials, tree_depth)
  NNC(quadtree, points, masses, potentials, tree_depth)

  # actual results (construct full matrix, O(N^2) approach)
  # NOTE: points and masses will have been reordered by the color sort, but this is acceptable as the particles themselves have not been changed
  correct_potentials = similar(potentials)
  correct_potentials .= zero(correct_potentials[1])
  kernel_mtx::Array{ComplexF64, 2} = log.(points .- points')
  foreach(i -> kernel_mtx[i, i] = zero(kernel_mtx[1, 1]), 1:length(points))
  # can't do simple dot product unfortunately
  correct_potentials = vec(sum(kernel_mtx.*masses, dims=2))

  println(correct_potentials)
  println(potentials)
  @assert correct_potentials ≈ potentials

end


function runTests()
  testFullMultipoleComputation()
end



runTests()
