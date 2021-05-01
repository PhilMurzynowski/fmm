"""
Four color version of the dutch national flag problem
"""

# swap two elements in arra
function swap_elements(array, i, j)
  temp = array[j]
  array[j] = array[i]
  array[i] = temp
end

# for complex numbers
  # red     =>  1:a-1
  # orange  =>  a:c 
  # green   =>  c+1:d
  # blue    =>  d+1:end
function fourColorSort!(array::Array{ComplexF64, 1}, center_value, first, last)
  center_x = real(center_value)
  center_y = imag(center_value)
  len = length(array)
  # four pointers to the four colors
  a = first
  b = first
  c = last
  d = last
  while (b <= c)
    if (imag(array[b]) >= center_y)
      if (real(array[b]) <= center_x)
        swap_elements(array, a, b)
        a += 1
        b += 1
      else
        # b : top right corner
        b += 1
      end
    else 
      if (real(array[b]) <= center_x)
        # c : bottom left corner
        swap_elements(array, b, c)
        c -= 1
      else
        # d : bottom right corner
        swap_elements(array, b, c)
        swap_elements(array, c, d)
        c -= 1
        d -= 1
      end
    end
  end
  @assert b == c + 1
  return (a, c, d)
end

# now four color sort based off of the complex numbers but have format Array{Tuple{ComplexF64, Float64}, 1}
# as have some extra paramter such as mass to keep with the 
function fourColorSort!(array::Array{Tuple{ComplexF64, Float64}, 1}, center_value, first, last)
  center_x = real(center_value)
  center_y = imag(center_value)
  len = length(array)
  # four pointers to the four colors
  a = first
  b = first
  c = last
  d = last
  while (b <= c)
    if (imag(array[b][1]) >= center_y)
      if (real(array[b][1]) <= center_x)
        swap_elements(array, a, b)
        a += 1
        b += 1
      else
        # b : top right corner
        b += 1
      end
    else 
      if (real(array[b][1]) <= center_x)
        # c : bottom left corner
        swap_elements(array, b, c)
        c -= 1
      else
        # d : bottom right corner
        swap_elements(array, b, c)
        swap_elements(array, c, d)
        c -= 1
        d -= 1
      end
    end
  end
  @assert b == c + 1
  return (a, c, d)
end
