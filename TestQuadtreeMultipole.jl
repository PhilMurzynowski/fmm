"""

Test multipole code

"""


#include("ArrayQuadtreeV2.jl") #included by below file
include("QuadtreeMultipole.jl")


using Random


function testFullMultipoleComputation()
  # an end to end test for one time step

  # specify tree parameters and construct tree
  tree_depth::Int = 5
  side_length::Float64 = 1.0
  quadtree::Array{Box, 1} = buildQuadtree(tree_depth, side_length)
  # generate points and masses
  num_bodies::Int = 5000
  #Random.seed!(1)
  points::Array{ComplexF64, 1} = rand(ComplexF64, num_bodies) * side_length
  #println(points)
  masses::Array{Float64, 1} = rand(Float64, num_bodies)
  #masses::Array{Float64, 1} = ones(Float64, num_bodies)
  # associate quadtree with point masses
  updateQuadtreePointMasses(quadtree, points, masses, tree_depth, side_length)

  # visual
  gr(size = (3000, 3000))
  p = plot()
  displayQuadtree(p, quadtree, points, masses, tree_depth)
  xlimits = (-0.1, 1.1)
  ylimits = xlimits
  xlims!(xlimits)
  ylims!(ylimits)
  gui()

  # upward pass
  P2M(quadtree, points, masses, tree_depth)
  M2M(quadtree, tree_depth)
  # downward pass
  M2L(quadtree, tree_depth)
  L2L(quadtree, tree_depth)
  # run M2L on finest level to factor in contributions from L2L
  # array to hold computed potentials
  potentials::Array{ComplexF64, 1} = Array{ComplexF64, 1}(undef, num_bodies)
  potentials .= zero(potentials[1])
  # remainder of downward pass
  L2P(quadtree, points, potentials, tree_depth)
  NNC(quadtree, points, masses, potentials, tree_depth)

  # potential is actually only the real component
  potentials = real(potentials)

  # actual results (construct full matrix, O(N^2) approach)
  # NOTE: points and masses will have been reordered by the color sort, but this is acceptable as the particles themselves have not been changed
  correct_potentials = similar(potentials)
  correct_potentials .= zero(correct_potentials[1])
  kernel_mtx::Array{ComplexF64, 2} = log.(transpose(points) .- points)
  foreach(i -> kernel_mtx[i, i] = zero(kernel_mtx[1, 1]), 1:length(points))
  #show(stdout, "text/plain", kernel_mtx)
  # can't do simple dot product unfortunately
  correct_potentials = real(vec(sum(kernel_mtx.*masses, dims=1)))

  @assert correct_potentials ≈ potentials

end


function runTests()
  testFullMultipoleComputation()
  #testStatic()
end

function testStatic()
  tree_depth::Int = 3
  quadtree::Array{Box, 1} = buildQuadtree(tree_depth)
  # set specific masses in specific locations for testing
  # assuming default side_length = 1.0
  num_bodies::Int = 7
  points::Array{ComplexF64, 1} = Array{ComplexF64, 1}(undef, num_bodies)
  points[1] = 0.2 + 0.9im # box 9 in depth 3
  points[2] = 0.1 + 0.6im # box 4 in depth 3
  points[3] = 0.3 + 0.8im # box 18 in depth 3
  points[4] = 0.3 + 0.7im # box 19 in depth 3
  points[5] = 0.2 + 0.4im # box 13 in depth 3
  points[6] = 0.4 + 0.1im # box 32 in depth 3
  points[7] = 0.1 + 0.8im # box 2 in depth 3
  # masses are not relevant for this test
  masses::Array{Float64, 1} = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0]
  @assert length(points) == num_bodies
  @assert length(masses) == num_bodies
  updateQuadtreePointMasses(quadtree, points, masses, tree_depth)

  # array to hold computed potentials
  potentials::Array{ComplexF64, 1} = Array{ComplexF64, 1}(undef, num_bodies)
  potentials .= zero(potentials[1])

  P2M(quadtree, points, masses, tree_depth)
  M2M(quadtree, tree_depth)
  M2L(quadtree, tree_depth)
  L2L(quadtree, tree_depth)
  L2P(quadtree, points, potentials, tree_depth);
  NNC(quadtree, points, masses, potentials, tree_depth)

  # actual results (construct full matrix, O(N^2) approach)
  # NOTE: points and masses will have been reordered by the color sort, but this is acceptable as the particles themselves have not been changed
  correct_potentials = similar(potentials)
  correct_potentials .= zero(correct_potentials[1])
  kernel_mtx::Array{ComplexF64, 2} = log.(transpose(points) .- points)
  foreach(i -> kernel_mtx[i, i] = zero(kernel_mtx[1, 1]), 1:length(points))
  # can't do simple dot product unfortunately
  correct_potentials = vec(sum(kernel_mtx.*masses, dims=1))

  println(correct_potentials)
  println(potentials)
  @assert real(correct_potentials) ≈ real(potentials)

end


runTests()
