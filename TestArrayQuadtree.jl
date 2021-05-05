include("ArrayQuadtree.jl")

function testFindParentIdx()

  # depth 2
  @assert findParentIdx(2, 1) ==  1
  @assert findParentIdx(2, 2) ==  1
  @assert findParentIdx(2, 3) ==  2
  @assert findParentIdx(2, 4) ==  2
  @assert findParentIdx(2, 5) ==  1
  @assert findParentIdx(2, 6) ==  1
  @assert findParentIdx(2, 7) ==  2
  @assert findParentIdx(2, 8) ==  2
  @assert findParentIdx(2, 9) ==  3
  @assert findParentIdx(2, 10) == 3
  @assert findParentIdx(2, 11) == 4
  @assert findParentIdx(2, 12) == 4
  @assert findParentIdx(2, 13) == 3
  @assert findParentIdx(2, 14) == 3
  @assert findParentIdx(2, 15) == 4
  @assert findParentIdx(2, 16) == 4

  # depth 3
  @assert findParentIdx(3, 1) ==  1
  @assert findParentIdx(3, 2) ==  1
  @assert findParentIdx(3, 3) ==  2
  @assert findParentIdx(3, 4) ==  2
  @assert findParentIdx(3, 5) ==  3
  @assert findParentIdx(3, 6) ==  3
  @assert findParentIdx(3, 7) ==  4
  @assert findParentIdx(3, 8) ==  4
  @assert findParentIdx(3, 9) ==  1
  @assert findParentIdx(3, 10) == 1
  @assert findParentIdx(3, 11) == 2
  @assert findParentIdx(3, 12) == 2
  @assert findParentIdx(3, 13) == 3
  @assert findParentIdx(3, 14) == 3
  @assert findParentIdx(3, 15) == 4
  @assert findParentIdx(3, 16) == 4
  @assert findParentIdx(3, 17) == 5
  @assert findParentIdx(3, 18) == 5
  @assert findParentIdx(3, 19) == 6
  @assert findParentIdx(3, 20) == 6
  @assert findParentIdx(3, 21) == 7
  @assert findParentIdx(3, 22) == 7
  @assert findParentIdx(3, 23) == 8
  @assert findParentIdx(3, 24) == 8
  @assert findParentIdx(3, 25) == 5
  @assert findParentIdx(3, 26) == 5
  @assert findParentIdx(3, 27) == 6
  @assert findParentIdx(3, 28) == 6
  @assert findParentIdx(3, 29) == 7
  @assert findParentIdx(3, 30) == 7
  @assert findParentIdx(3, 31) == 8
  @assert findParentIdx(3, 32) == 8
  @assert findParentIdx(3, 33) == 9
  @assert findParentIdx(3, 34) == 9
  @assert findParentIdx(3, 35) == 10
  @assert findParentIdx(3, 36) == 10
  @assert findParentIdx(3, 37) == 11
  @assert findParentIdx(3, 38) == 11
  @assert findParentIdx(3, 39) == 12
  @assert findParentIdx(3, 40) == 12
  @assert findParentIdx(3, 41) == 9
  @assert findParentIdx(3, 42) == 9
  @assert findParentIdx(3, 43) == 10
  @assert findParentIdx(3, 44) == 10
  @assert findParentIdx(3, 45) == 11
  @assert findParentIdx(3, 46) == 11
  @assert findParentIdx(3, 47) == 12
  @assert findParentIdx(3, 48) == 12
  @assert findParentIdx(3, 49) == 13
  @assert findParentIdx(3, 50) == 13
  @assert findParentIdx(3, 51) == 14
  @assert findParentIdx(3, 52) == 14
  @assert findParentIdx(3, 53) == 15
  @assert findParentIdx(3, 54) == 15
  @assert findParentIdx(3, 55) == 16
  @assert findParentIdx(3, 56) == 16
  @assert findParentIdx(3, 57) == 13
  @assert findParentIdx(3, 58) == 13
  @assert findParentIdx(3, 59) == 14
  @assert findParentIdx(3, 60) == 14
  @assert findParentIdx(3, 61) == 15
  @assert findParentIdx(3, 62) == 15
  @assert findParentIdx(3, 63) == 16
  @assert findParentIdx(3, 64) == 16

  # a few quick depth 4 checks
  @assert findParentIdx(4, 33) == 9
  @assert findParentIdx(4, 65) == 17
  @assert findParentIdx(4, 84) == 18

end


function testFindChildrenIdxs()

  # depth 1
  @assert issetequal(findChildrenIdxs(1, 1), [1, 2, 5, 6])
  @assert issetequal(findChildrenIdxs(1, 2), [3, 4, 7, 8])
  @assert issetequal(findChildrenIdxs(1, 3), [9, 10, 13, 14])
  @assert issetequal(findChildrenIdxs(1, 4), [11, 12, 15, 16])

  # some depth 2 tests
  @assert issetequal(findChildrenIdxs(2, 1), [1, 2, 9, 10])
  @assert issetequal(findChildrenIdxs(2, 2), [3, 4, 11, 12])
  @assert issetequal(findChildrenIdxs(2, 6), [19, 20, 27, 28])

end


function testFindNeighborIdxs()
  
  # depth 1
  @assert issetequal(findNeighborIdxs(1, 1), [2, 3, 4])
  @assert issetequal(findNeighborIdxs(1, 2), [1, 3, 4])
  @assert issetequal(findNeighborIdxs(1, 3), [2, 1, 4])
  @assert issetequal(findNeighborIdxs(1, 4), [2, 3, 1])

  # depth 2 tests
  @assert issetequal(findNeighborIdxs(2, 1), [2, 5, 6])
  @assert issetequal(findNeighborIdxs(2, 2), [1, 3, 5, 6, 7])
  @assert issetequal(findNeighborIdxs(2, 3), [2, 4, 6, 7, 8])
  @assert issetequal(findNeighborIdxs(2, 4), [3, 7, 8])
  @assert issetequal(findNeighborIdxs(2, 5), [1, 9, 2, 6, 10])
  @assert issetequal(findNeighborIdxs(2, 6), [1, 2, 3, 5, 7, 9, 10, 11])
  @assert issetequal(findNeighborIdxs(2, 7), [2, 3, 4, 6, 8, 10, 11, 12])
  @assert issetequal(findNeighborIdxs(2, 8), [3, 4, 7, 11, 12])
  @assert issetequal(findNeighborIdxs(2, 9), [5, 6, 10, 13, 14])
  @assert issetequal(findNeighborIdxs(2, 10), [5, 6, 7, 9, 11, 13, 14, 15])
  @assert issetequal(findNeighborIdxs(2, 11), [6, 7, 8, 10, 12, 14, 15, 16])
  @assert issetequal(findNeighborIdxs(2, 12), [7, 8, 11, 15, 16])
  @assert issetequal(findNeighborIdxs(2, 13), [9, 10, 14])
  @assert issetequal(findNeighborIdxs(2, 14), [9, 10, 11, 13, 15])
  @assert issetequal(findNeighborIdxs(2, 15), [10, 11, 12, 14, 16])
  @assert issetequal(findNeighborIdxs(2, 16), [11, 12, 15])

  # some depth 3 tests
  @assert issetequal(findNeighborIdxs(3, 1), [2, 9, 10])
  @assert issetequal(findNeighborIdxs(3, 2), [1, 3, 9, 10, 11])
  @assert issetequal(findNeighborIdxs(3, 8), [7, 15, 16])
  @assert issetequal(findNeighborIdxs(3, 21), [12, 13, 14, 20, 22, 28, 29, 30])
  @assert issetequal(findNeighborIdxs(3, 32), [24, 40, 23, 31, 39])

end


function testFindInteractionIdxs()

  # depth 1
  @assert issetequal(findInteractionIdxs(1, 1), [])
  @assert issetequal(findInteractionIdxs(1, 2), [])
  @assert issetequal(findInteractionIdxs(1, 3), [])
  @assert issetequal(findInteractionIdxs(1, 4), [])

  # depth 2
  all_depth2_boxes::Array{Int, 1} = [x for x in 1:16]
  for i in all_depth2_boxes
    @assert issetequal(findInteractionIdxs(2, i), setdiff(all_depth2_boxes, vcat([i], findNeighborIdxs(2, i))))
  end

end

function testGetBoxCenter()

  side_length::Float64 = 1.0
  #depth 1
  @assert getBoxCenter(1, 1, side_length) ≈ 0.25 + 0.75im
  @assert getBoxCenter(1, 2, side_length) ≈ 0.25 + 0.25im
  @assert getBoxCenter(1, 3, side_length) ≈ 0.75 + 0.75im
  @assert getBoxCenter(1, 4, side_length) ≈ 0.75 + 0.25im

  # some depth 2 tests
  @assert getBoxCenter(2, 1, side_length) ≈ 0.125 + 0.875im
  @assert getBoxCenter(2, 4, side_length) ≈ 0.125 + 0.125im
  @assert getBoxCenter(2, 5, side_length) ≈ 0.375 + 0.875im

end

function testGetOffsetOfDepth()
  @assert getOffsetOfDepth(1) == 0
  @assert getOffsetOfDepth(2) == 4
  @assert getOffsetOfDepth(3) == 20
  @assert getOffsetOfDepth(4) == 84
end

function testPropagateFUp()
  depth::Int = 3
  tree::Array{Box, 1} = buildQuadtree(depth)
  num_particles::Int = 100
  masses::Array{Float64, 1} = rand(Float64, 100)
  propagateFUp(tree, masses, depth)
end

function testVisualBuildQuadtree()
  
  depth::Int = 4
  tree::Array{Box, 1} = buildQuadtree(depth)
  displayQuadtreeBoxes(tree, depth)

end

function testVisualQuadtreeColorSortWithMass()
  tree_depth::Int = 4
  side_length::Float64 = 1.0
  quadtree::Array{Box, 1} = buildQuadtree(tree_depth, side_length)
  # there is no depth 0 so must intially color sort
  # and past appropriate first and last to depth 1
  tree_center::ComplexF64 = side_length/2 + side_length/2*1im
  # generate points and masses
  num_bodies::Int = 10
  points::Array{ComplexF64, 1} = rand(ComplexF64, num_bodies) * side_length
  masses::Array{Float64, 1} = rand(Float64, num_bodies)

  #copies for testing
  original_points = copy(points)
  original_masses = copy(masses)
  
  
  a::Int, c::Int, d::Int = fourColorSort!(points, masses, tree_center, 1, length(points))
  tl_child::Box = quadtree[1]
  bl_child::Box = quadtree[2]
  tr_child::Box = quadtree[3]
  br_child::Box = quadtree[4]
  colorSortQuadtreePointMasses(tl_child, quadtree, points, masses, tree_depth, 1, first, a-1)
  colorSortQuadtreePointMasses(bl_child, quadtree, points, masses, tree_depth, 1, a, c)
  colorSortQuadtreePointMasses(tr_child, quadtree, points, masses, tree_depth, 1, c+1, d)
  colorSortQuadtreePointMasses(br_child, quadtree, points, masses, tree_depth, 1, d+1, last)

  # plot
  gr(size = (3000, 1500))
  xlimits = (-0.1, 1.1)
  ylimits = xlimits
  l = @layout [a  b]
  # plot original 
  x = real.(original_points)
  y = imag.(original_points)
  p1 = scatter((x, y), leg=false, label="", markersize=10*masses)
  xlims!(xlimits)
  ylims!(ylimits)
  
  # plot quadtree
  p2 = plot()
  displayQuadtree(p2, quadtree, points, masses, tree_depth, side_length)
  
  # show plot
  plot(p1, p2, layout=l)
  xlims!(xlimits)
  ylims!(ylimits)
  gui()
  
end

function runTests()
  testFindParentIdx()
  testFindChildrenIdxs()
  testFindNeighborIdxs()
  testFindInteractionIdxs()
  testGetBoxCenter()
  testGetOffsetOfDepth()
  testPropagateFUp()
end

function runVisualTests()
  testVisualBuildQuadtree()
  testVisualQuadtreeColorSortWithMass()
end


runTests()
runVisualTests()
