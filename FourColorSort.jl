
""" Four color version of the dutch national flag problem """

# swap two elements in array
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
#
# four color sort with points and associated masses
# in separate ararys, sorting by points
# NOTE: modified for column-major order
# first color is tl, second is bl, third tr, fourth br
function fourColorSort!(points, masses::Array{Float64, 1}, center_value, first, last)
  center_x = real(center_value)
  center_y = imag(center_value)
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
  # @assert b == c + 1
  return (a, c, d)
end

# same concept, generalized to three arrays
function fourColorSort!(first, second, third, center_value, start, last) center_x = real(center_value)
  center_x = real(center_value)
  center_y = imag(center_value)
  # four pointers to the four colors
  a = start 
  b = start 
  c = last
  d = last
  while (b <= c)
    if (real(first[b]) <= center_x)
      if (imag(first[b]) >= center_y)
        # top left corner
        swap_elements(first, a, b)
        swap_elements(second, a, b)
        swap_elements(third, a, b)
        a += 1
        b += 1
      else
        # b : bottom left corner
        b += 1
      end
    else 
      if (imag(first[b]) >= center_y)
        # c : top right corner 
        swap_elements(first, b, c)
        swap_elements(second, b, c)
        swap_elements(third, b, c)
        c -= 1
      else
        # d : bottom right corner
        swap_elements(first, b, c)
        swap_elements(first, c, d)
        swap_elements(second, b, c)
        swap_elements(second, c, d)
        swap_elements(third, b, c)
        swap_elements(third, c, d)
        c -= 1
        d -= 1
      end
    end
  end
  # @assert b == c + 1
  return (a, c, d)
end
