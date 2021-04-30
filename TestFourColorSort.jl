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

function largeVisualTest()
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
  plot(p1, p2, layout = l)
  gui()
end

function edgecaseVisualTest()
  gr(size = (3000, 1500))
  l = @layout [a b; c d]
  xlim = (-0.1, 1.1)
  ylim = xlim
  markersize = 20
  center = 0.5 + 0.5im

  array = [0.75 + 0.25im]
  a, c, d = fourColorSort!(array, center)
  x = real.(array[1:a-1])
  y = imag.(array[1:a-1])
  p1 = scatter((x, y), xlim = xlim, ylim = ylim, markersize=markersize, legend=false, color=:red, label="", aspectratio=1)
  x = real.(array[a:c])
  y = imag.(array[a:c])
  scatter!(p1, (x, y), xlim = xlim, ylim = ylim, markersize=markersize, legend=false, color=:orange, label="")
  x = real.(array[c+1:d])
  y = imag.(array[c+1:d])
  scatter!(p1, (x, y), xlim = xlim, ylim = ylim, markersize=markersize, legend=false, color=:green, label="")
  x = real.(array[d+1:end])
  y = imag.(array[d+1:end])
  scatter!(p1, (x, y), xlim = xlim, ylim = ylim, markersize=markersize, legend=false, color=:blue, label="")

  array = [0.25 + 0.25im]
  a, c, d = fourColorSort!(array, center)
  x = real.(array[1:a-1])
  y = imag.(array[1:a-1])
  p2 = scatter((x, y), xlim = xlim, ylim = ylim, markersize=markersize, legend=false, color=:red, label="", aspectratio=1)
  x = real.(array[a:c])
  y = imag.(array[a:c])
  scatter!(p2, (x, y), xlim = xlim, ylim = ylim, markersize=markersize, legend=false, color=:orange, label="")
  x = real.(array[c+1:d])
  y = imag.(array[c+1:d])
  scatter!(p2, (x, y), xlim = xlim, ylim = ylim, markersize=markersize, legend=false, color=:green, label="")
  x = real.(array[d+1:end])
  y = imag.(array[d+1:end])
  scatter!(p2, (x, y), xlim = xlim, ylim = ylim, markersize=markersize, legend=false, color=:blue, label="")

  array = [0.75 + 0.75im]
  a, c, d = fourColorSort!(array, center)
  x = real.(array[1:a-1])
  y = imag.(array[1:a-1])
  p3 = scatter((x, y), xlim = xlim, ylim = ylim, markersize=markersize, legend=false, color=:red, label="", aspectratio=1)
  x = real.(array[a:c])
  y = imag.(array[a:c])
  scatter!(p3, (x, y), xlim = xlim, ylim = ylim, markersize=markersize, legend=false, color=:orange, label="")
  x = real.(array[c+1:d])
  y = imag.(array[c+1:d])
  scatter!(p3, (x, y), xlim = xlim, ylim = ylim, markersize=markersize, legend=false, color=:green, label="")
  x = real.(array[d+1:end])
  y = imag.(array[d+1:end])
  scatter!(p3, (x, y), xlim = xlim, ylim = ylim, markersize=markersize, legend=false, color=:blue, label="")

  array = [0.25 + 0.75im]
  a, c, d = fourColorSort!(array, center)
  x = real.(array[1:a-1])
  y = imag.(array[1:a-1])
  p4 = scatter((x, y), xlim = xlim, ylim = ylim, markersize=markersize, legend=false, color=:red, label="", aspectratio=1)
  x = real.(array[a:c])
  y = imag.(array[a:c])
  scatter!(p4, (x, y), xlim = xlim, ylim = ylim, markersize=markersize, legend=false, color=:orange, label="")
  x = real.(array[c+1:d])
  y = imag.(array[c+1:d])
  scatter!(p4, (x, y), xlim = xlim, ylim = ylim, markersize=markersize, legend=false, color=:green, label="")
  x = real.(array[d+1:end])
  y = imag.(array[d+1:end])
  scatter!(p4, (x, y), xlim = xlim, ylim = ylim, markersize=markersize, legend=false, color=:blue, label="")

  plot(p1, p2, p3, p4, layout = l)
  gui()
end


""" tests run """

smallTest()
largeVisualTest()
edgecaseVisualTest()
