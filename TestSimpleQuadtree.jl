""" Test Quadtree """

include("SimpleQuadtree.jl")


function visualTestBoxesSpanFullSpace()
  array = rand(ComplexF64, 100) .- (0.5 + 0.5im)
  center = 0.0 + 0.0im
  dim = 1.0
  depth = 3
  displayQuadtreeBoxesWrapper(array, center, dim, depth)
end

function visualTestAllPointsInQuadtree()
  array = rand(ComplexF64, 100) .- (0.5 + 0.5im)
  center = 0.0 + 0.0im
  dim = 0.5
  depth = 4
  gr(size = (3000, 1500))
  xlimits = (-0.5, 0,5)
  ylimits = xlimits
  l = @layout [a  b]
  x = real.(array)
  y = imag.(array)
  p1 = scatter((x, y), leg=false, label="")
  xlims!(xlimits)
  ylims!(ylimits)
  p2 = displayQuadtreeBoxesAndParticlesWrapper(array, center, dim, depth)
  xlims!(xlimits)
  ylims!(ylimits)
  plot(p1, p2, layout=l)
  xlims!(xlimits)
  ylims!(ylimits)
  gui()
end


function testNumberOfPointsInQuadtree()
  # more specifically check if number of all points in quadtree
  # is the same as number of points in original array
  num_particles = 201
  array = rand(ComplexF64, num_particles) .- (0.5 + 0.5im)
  center = 0.0 + 0.0im
  dim = 0.5
  depth = 4
  tree = buildQuadtree!(array, nothing, center, dim, 1, length(array), depth)
  @assert countNumberOfPointsInQuadtree(tree) == num_particles
end

function testNumberOfPointsWithMass()
  num_particles = 201
  array = Array{Tuple{ComplexF64, Float64}, 1}(undef, num_particles)
  map(x -> (rand(ComplexF64) .- (0.5 + 0.5im), rand(Float64)), array)
  center = 0.0 + 0.0im
  dim = 0.5
  depth = 4
  tree = buildQuadtree!(array, nothing, center, dim, 1, length(array), depth)
  @assert countNumberOfPointsInQuadtree(tree) == num_particles
end

function visualTestAllPointsWithMass()
  num_particles = 201
  array = Array{Tuple{ComplexF64, Float64}, 1}(undef, num_particles)
  array = [(rand(ComplexF64) - (0.5 + 0.5im), rand(Float64)) for x in array]
  center = 0.0 + 0.0im
  dim = 0.5
  depth = 4
  gr(size = (3000, 1500))
  xlimits = (-0.5, 0,5)
  ylimits = xlimits
  l = @layout [a  b]
  points = first.(array)
  masses = getfield.(array, 2)
  x = real.(points)
  y = imag.(points)
  p1 = scatter((x, y), leg=false, label="", markersize=10*masses)
  xlims!(xlimits)
  ylims!(ylimits)
  p2 = displayQuadtreeBoxesParticlesMassWrapper(array, center, dim, depth)
  xlims!(xlimits)
  ylims!(ylimits)
  plot(p1, p2, layout=l)
  xlims!(xlimits)
  ylims!(ylimits)
  gui()
end

function testMass(tree, array)
  if depth > 0 
    for subtree in tree.children
      testMass(subtree, array)
    end
  end
  @assert tree.massExpansion = sum(getfield.(array[tree.first, tree.last], 2))
end

function testPropagateMassUp()
  num_particles = 201
  array = Array{Tuple{ComplexF64, Float64}, 1}(undef, num_particles)
  array = [(rand(ComplexF64) - (0.5 + 0.5im), rand(Float64)) for x in array]
  center = 0.0 + 0.0im
  dim = 0.5
  depth = 4

  tree = buildQuadtree!(array, nothing, center, dim, 1, length(array), depth)
  propagateMassUp!(tree, array)
  return
end


visualTestBoxesSpanFullSpace()
visualTestAllPointsInQuadtree()
testNumberOfPointsInQuadtree()
testNumberOfPointsWithMass()
visualTestAllPointsWithMass()
testPropagateMassUp()
