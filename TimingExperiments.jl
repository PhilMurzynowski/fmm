
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

#const array_size = 10
#const n = 30
#arr = rand(ComplexF64, array_size)
#arr2 = rand(ComplexF64, array_size)
#ans = similar(arr)
#ans .= zero(ans[1])
#ans2 = similar(arr2)
#ans2 .= zero(ans2[1])
#@btime arrayPowerOperation!(ans, arr, n)
#@btime fastArrayPowerOperation!(ans2, arr2, n)
#println()

# UNOPTIMIZED n = 30
#  304.392 ns (0 allocations: 0 bytes)
function complexPowerOperation(in::ComplexF64, n::Int)
  out = zero(in)
  for k in 1:n 
    out += in^k
  end
  return out
end

# OPTIMIZED n = 30
#  60.182 ns (0 allocations: 0 bytes)
#  ~5x speedup
function fastComplexPowerOperation(in::ComplexF64, n::Int)
  out = zero(in)
  tmp = in
  for k in 1:n 
    out += tmp
    tmp *= in
  end
  return out
end


const n = 30
in1 = rand(ComplexF64)
in2 = rand(ComplexF64)
@btime complexPowerOperation(in1, n)
@btime fastComplexPowerOperation(in2, n)
@btime complexPowerOperation(in1, n)

# LOOK INTO DOT PRODUCTS, complex with complex, etc.
# to replace sum(x .* y)

# OPTIMIZE BINOMIAL using PASCAL's triangle, as incremeting one at a time, profile
# DONE: using LUT (look up tables)

