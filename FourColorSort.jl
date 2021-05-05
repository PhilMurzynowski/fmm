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

# four color sort with points and associated masses
# in separate ararys, sorting by points
# NOTE: modified for column-major order
# first color is tl, second is bl, third tr, fourth br
function fourColorSort!(points::Array{ComplexF64, 1}, masses::Array{Float64, 1}, center_value, first, last)
  center_x = real(center_value)
  center_y = imag(center_value)
  len = length(points)
  # four pointers to the four colors
  a = first
  b = first
  c = last
  d = last
  while (b <= c)
    if (real(points[b]) <= center_x)
      if (imag(points[b]) >= center_y)
        # top left corner
        swap_elements(points, a, b)
        swap_elements(masses, a, b)
        a += 1
        b += 1
      else
        # b : bottom left corner
        b += 1
      end
    else 
      if (imag(points[b]) >= center_y)
        # c : top right corner 
        swap_elements(points, b, c)
        swap_elements(masses, b, c)
        c -= 1
      else
        # d : bottom right corner
        swap_elements(points, b, c)
        swap_elements(points, c, d)
        swap_elements(masses, b, c)
        swap_elements(masses, c, d)
        c -= 1
        d -= 1
      end
    end
  end
  @assert b == c + 1
  return (a, c, d)
end
