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
  array = rand(complexf64, 100) .- (0.5 + 0.5im)
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
  p2 = displayQuadtreeBoxesAndParticlesWrapper(array, center, dim, depth)
  plot(p1, p2, layout=l)
  xlims!(xlimits)
  ylims!(ylimits)
  gui()
end


function testNumberOfPointsInQuadtree()
  # more specifically check if number of all points in quadtree
  # is the same as number of points in original array
  array = rand(ComplexF64, 100) .- (0.5 + 0.5im)
  center = 0.0 + 0.0im
  dim = 0.5
  depth = 4
  tree = buildQuadtree!(array, nothing, center, dim, 1, length(array), depth)
  # TODO
end


visualTestBoxesSpanFullSpace()
visualTestAllPointsInQuadtree()
testNumberOfPointsInQuadtree()
