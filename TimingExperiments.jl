
""" Exploring optimizations in Julia"""

using BenchmarkTools

function arrayPowerOperation!(ans::Array{ComplexF64, 1}, arr::Array{ComplexF64, 1}, n::Int)
  for k in 1:n 
    ans += arr.^k
  end
  return 
end

# OPTIMIZED:  array_size   n    speedup
#             100          50   ~6x
#             10           30   ~3.2
function fastArrayPowerOperation!(ans::Array{ComplexF64, 1}, arr::Array{ComplexF64, 1}, n::Int)
  tmp = copy(arr)
  for k in 1:n 
    ans += tmp
    tmp .*= arr
  end
  return 
end

# LOOK INTO DOT PRODUCTS, complex with complex, etc.
# to replace sum(x .* y)


const array_size = 10
const n = 30
arr = rand(ComplexF64, array_size)
arr2 = rand(ComplexF64, array_size)
ans = similar(arr)
ans .= zero(ans[1])
ans2 = similar(arr2)
ans2 .= zero(ans2[1])
@btime arrayPowerOperation!(ans, arr, n)
@btime fastArrayPowerOperation!(ans2, arr2, n)
println()
