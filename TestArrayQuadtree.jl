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





testFindParentIdx()
