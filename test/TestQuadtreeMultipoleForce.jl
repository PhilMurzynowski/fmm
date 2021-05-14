"""

Test multipole code to calculate forces

"""


#include("ArrayQuadtreeV2.jl") #included by below file
include("QuadtreeMultipoleForce.jl")


using Random


function testFullForceComputation()
  # an end to end test of forces computation for one time step

  # specify tree parameters and construct tree
  tree_depth::Int = 3
  side_length::Float64 = 1.0
  quadtree::Array{Box, 1} = buildQuadtree(tree_depth, side_length)
  # generate points and masses
  num_bodies::Int = 100
  #Random.seed!(1)
  points::Array{complexf64, 1} = rand(complexf64, num_bodies) * side_length
  #println(points)
  masses::Array{Float64, 1} = rand(Float64, num_bodies)
  #masses::Array{Float64, 1} = ones(Float64, num_bodies)
  # associate quadtree with point masses
  updateQuadtreePointMasses(quadtree, points, masses, tree_depth, side_length)

  # visual
  """
  gr(size = (3000, 3000))
  p = plot()
  displayQuadtree(p, quadtree, points, masses, tree_depth)
  xlimits = (-0.1, 1.1)
  ylimits = xlimits
  xlims!(xlimits)
  ylims!(ylimits)
  gui()
  """

  # upward pass
  P2M_F(quadtree, points, masses, tree_depth)
  M2M_F(quadtree, tree_depth)
  # downward pass
  M2L_F(quadtree, tree_depth)
  L2L_F(quadtree, tree_depth)

  # using notation from greengard, ω
  # potential is Re(ω)
  # ω' = ω_p (derivative of ω, available analytically)
  # force is (Re(ω'), -Im(ω'))
  ω_p::Array{ComplexF64, 1} = Array{ComplexF64, 1}(undef, num_bodies)
  ω_p .= zero(ω_p[1])
  # remainder of downward pass
  L2P_F(quadtree, points, ω_p, tree_depth)
  NNC_F(quadtree, points, masses, ω_p, tree_depth)

  forces = [real(ω_p), -imag(ω_p)]

  # actual results (construct full matrix, O(N^2) approach)
  # NOTE: points and masses will have been reordered by the color sort, but this is acceptable as the particles themselves have not been changed
  correct_ω_p = similar(ω_p)
  correct_ω_p .= zero(correct_ω_p[1])
  kernel_mtx::Array{ComplexF64, 2} = 1 ./ (transpose(points) .- points)
  foreach(i -> kernel_mtx[i, i] = zero(kernel_mtx[1, 1]), 1:length(points))
  #show(stdout, "text/plain", kernel_mtx)
  # can't do simple dot product unfortunately
  correct_ω_p = vec(sum(kernel_mtx.*masses, dims=1))
  correct_forces = [real(correct_ω_p), -imag(correct_ω_p)]

  #println(correct_forces)
  #println(forces)
  @assert correct_forces ≈ forces 

end


function runTests()
  testFullForceComputation()
end

runTests()
