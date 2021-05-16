
""" Standard construction of a table holding binomial coefficients """

# Results for different k values and a given n are stored in a column
# as that is the most frequent access pattern.
# Not optimized but only called once in setup.
# For binomial(n, k), use table[k+1, n+1]
function binomialTable(n::Int)
  table = Array{Int64, 2}(undef, n+1, 2*n+1)
  for i in 0:2*n
    for j in 0:n
      table[j+1, i+1] = binomial(i, j)
    end
  end
  return table
end
