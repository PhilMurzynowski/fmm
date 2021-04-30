""" Test four color sorting """


using Plots


include("FourColorSort.jl")



""" test functions """

function smallTest()
  # bottom right, bottom left, top right, top left
  array = [1 + -1im, -1 + -1im, 1 + 1im, -1 + 1im]
  center = 0 + 0im
  fourColorSort!(array, center)
  @assert array == [-1 + 1im, 1 + 1im, -1 - 1im, 1 - 1im]
end

function smallVisualTest()
  gr(size = (3000, 1500))
  l = @layout [a  b]
  array = rand(ComplexF64, 100)
  xlim = (-0.1, 1.1)
  ylim = xlim
  markersize = 20
  p1 = scatter(real.(array), imag.(array), xlim=xlim, ylim=ylim, markersize=markersize, legend=false)
  center = 0.5 + 0.5im
  a, c, d = fourColorSort!(array, center)
  x = real.(array[1:a-1])
  y = imag.(array[1:a-1])
  p2 = scatter((x, y), xlim = xlim, ylim = ylim, markersize=markersize, legend=false, color=:red, label="")
  x = real.(array[a:c])
  y = imag.(array[a:c])
  scatter!(p2, (x, y), xlim = xlim, ylim = ylim, markersize=markersize, legend=false, color=:orange, label="")
  x = real.(array[c+1:d])
  y = imag.(array[c+1:d])
  scatter!(p2, (x, y), xlim = xlim, ylim = ylim, markersize=markersize, legend=false, color=:green, label="")
  x = real.(array[d+1:end])
  y = imag.(array[d+1:end])
  scatter!(p2, (x, y), xlim = xlim, ylim = ylim, markersize=markersize, legend=false, color=:blue, label="")
  plot(p1, p2)
  gui()
end


""" tests run """

smallTest()
smallVisualTest()
