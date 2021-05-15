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

  # NOTE: may need to test deeper if optimize out propagating to depth 1
  # depth 2
  tree_depth::Int = 2
  tree::Array{Box, 1} = buildQuadtree(tree_depth)
  # set specific masses in specific locations for testing
  # assuming default side_length = 1.0
  num_bodies::Int = 7
  points::Array{ComplexF64, 1} = Array{ComplexF64, 1}(undef, num_bodies)
  points[1] = 0.2 + 0.9im # box 1 in depth 2
  points[2] = 0.1 + 0.6im # box 2 in depth 2
  points[3] = 0.3 + 0.8im # box 5 in depth 2
  points[4] = 0.3 + 0.7im # box 6 in depth 2
  points[5] = 0.2 + 0.4im # box 3 in depth 2
  points[6] = 0.4 + 0.1im # box 8 in depth 2
  points[7] = 0.1 + 0.8im # box 1 in depth 2
  masses::Array{Float64, 1} = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0]
  @assert length(points) == num_bodies
  @assert length(masses) == num_bodies

  updateQuadtreePointMasses(tree, points, masses, tree_depth)
  propagateFUp(tree, masses, tree_depth)
  
  depth_2_offset::Int = 4
  # leaves
  @assert tree[depth_2_offset + 1].f == 65.0 
  @assert tree[depth_2_offset + 2].f == 2.0 
  @assert tree[depth_2_offset + 5].f == 4.0 
  @assert tree[depth_2_offset + 6].f == 8.0 
  @assert tree[depth_2_offset + 3].f == 16.0 
  @assert tree[depth_2_offset + 8].f == 32.0 
  @assert tree[depth_2_offset + 4].f == 0.0 
  @assert tree[depth_2_offset + 7].f == 0.0 
  @assert tree[depth_2_offset + 9].f == 0.0 
  @assert tree[depth_2_offset + 10].f == 0.0 
  @assert tree[depth_2_offset + 11].f == 0.0 
  @assert tree[depth_2_offset + 12].f == 0.0 
  @assert tree[depth_2_offset + 13].f == 0.0 
  @assert tree[depth_2_offset + 14].f == 0.0 
  @assert tree[depth_2_offset + 15].f == 0.0 
  @assert tree[depth_2_offset + 16].f == 0.0 
  # depth 1
  @assert tree[1].f == 79.0 
  @assert tree[2].f == 48.0 
  @assert tree[3].f == 0.0 
  @assert tree[4].f == 0.0 

end

function testPropagateDownPotential()
  # test for depth 3 tree
  # propagate down from depth 2 to depth 3 boxes
  
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

  # array for potential of each point
  us::Array{Float64, 1} = zeros(Float64, num_bodies)

  # artificially set potentials in boxs depth 2 and higher
  # depth 3
  for i in 1:4^3
    global_idx::Int = 20+i
    quadtree[global_idx].u = 2^(i-1)
  end
  # depth 2
  for i in 1:4^2
    global_idx::Int = 4+i
    quadtree[global_idx].u = 2^(i+4^3-1)
  end
  propagateDownPotential(quadtree, us, tree_depth)

  # depth 2 should be unaffected
  for i in 1:4^2
    global_idx::Int = 4+i
    @assert quadtree[global_idx].u == 2^(i+4^3-1)
  end
  # depth 3 should have contribution from depth 2
  for i in 1:4^3
    global_idx::Int = 20+i
    box = quadtree[global_idx]
    @assert box.u == 2^(i-1) + quadtree[4+box.parent_idx].u
  end
  
  # check right boxes have or don't have points
  for i in 1:4^3
    global_idx::Int = 20+i
    box = quadtree[global_idx]
    #@printf "i: %d\n" i
    #@printf "box.start_idx: %d, box.final_idx: %d\n" box.start_idx box.final_idx
    if (i == 9 ||
        i == 4 ||
        i == 18 ||
        i == 19 ||
        i == 13 ||
        i == 32 ||
        i == 2)
      @assert boxHasPoints(box)
    else 
      @assert !boxHasPoints(box)
    end
  end

  # individual potential tests
  for i in 1:4^3
    global_idx::Int = 20+i
    box = quadtree[global_idx]
    #@printf "i: %d\n" i
    #@printf "box.start_idx: %d, box.final_idx: %d\n" box.start_idx box.final_idx
    if (i == 9 ||
        i == 4 ||
        i == 18 ||
        i == 19 ||
        i == 13 ||
        i == 32 ||
        i == 2)
      for idx in box.start_idx:box.final_idx
        @assert us[idx] == box.u
      end
    end
  end


end

function testComputeNeighborPotentialContribution()

end

function testFullTimestepComputation()

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
  # generate points and masses
  num_bodies::Int = 1000
  points::Array{ComplexF64, 1} = rand(ComplexF64, num_bodies) * side_length
  masses::Array{Float64, 1} = rand(Float64, num_bodies)

  #copies for testing
  original_points = copy(points)
  original_masses = copy(masses)

  updateQuadtreePointMasses(quadtree, points, masses, tree_depth, side_length)
  
  # plot
  gr(size = (3000, 1500))
  xlimits = (-0.1, 1.1)
  ylimits = xlimits
  l = @layout [a  b]
  # plot original 
  x = real.(original_points)
  y = imag.(original_points)
  p1 = scatter((x, y), leg=false, label="", markersize=10*original_masses)
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
  testPropagateDownPotential()
  testComputeNeighborPotentialContribution()
  testFullTimestepComputation()
end

function runVisualTests()
  # testVisualBuildQuadtree()
  testVisualQuadtreeColorSortWithMass()
end


runTests()
runVisualTests()
