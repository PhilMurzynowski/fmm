"""
Construct a quadtree using the four color sort.
For clarity keeping as mutable struct (may have performance repurcussions)
In order to have a small constant number of particles in each leaf box will
assume near-uniform distribution then calculate necessary depth to achieve
expected constant number of particles in leaves.
"""

include("FourColorSort.jl")


using Plots


# realize this is obscure, first draft notation
# red     =>  first:a-1
# orange  =>  a:c 
# green   =>  c+1:d
# blue    =>  d+1:last
mutable struct Box
  depth::Int                  # useful for verification
  massExpansion::Float64
  first::Int
  last::Int
  parent::Union{Box, Nothing}
  #top_left_child::Box         # massExpansionly convert children to array
  #top_right_child::Box
  #bot_left_child::Box
  #bot_right_child::Box
  #children::Union{Array{Box, 1}, Nothing}
  children::Array{Box, 1}
  Box(depth::Int, massExpansion::Float64, first::Int, last::Int) = new(depth, massExpansion, first, last, nothing)
  Box(depth::Int, massExpansion::Float64, first::Int, last::Int, parent::Box) = new(depth, massExpansion, first, last, parent)
  Box(depth::Int, massExpansion::Float64, first::Int, last::Int, parent::Nothing) = new(depth, massExpansion, first, last, nothing)
end


""" Construct the quadtree given an array of complex numbers """
function buildQuadtree!(array, parent, center, dim, first, last, depth)

  if depth == 0 return Box(0, 0.0, first, last, parent) end
  a, c, d = fourColorSort!(array, center, first, last)
  #println("depth: $depth")
  #println("first: $first, a: $a, c:$c, d:$d, last:$last")
  half_dim = dim/2
  box = Box(depth, 0.0, first, last, parent)
  top_left_center = center - half_dim + half_dim*1im
  top_right_center = center + half_dim + half_dim*1im
  bot_left_center = center - half_dim - half_dim*1im 
  bot_right_center = center + half_dim - half_dim*1im

  box.children = Array{Box, 1}()
  push!(box.children, buildQuadtree!(array, box, top_left_center, half_dim, first, a-1, depth-1))
  push!(box.children, buildQuadtree!(array, box, top_right_center, half_dim, a, c, depth-1))
  push!(box.children, buildQuadtree!(array, box, bot_left_center, half_dim, c+1, d, depth-1))
  push!(box.children, buildQuadtree!(array, box, bot_right_center, half_dim, d+1, last, depth-1))

  return box

end


""" Display quadtree boxes helper"""
function displayQuadtreeBoxes(plot, tree, center, dim)
  half_dim = dim/2
  if tree.depth == 0
    plot!(plot, Shape(real(center) .+ [-dim, dim, dim, -dim], imag(center) .+ [-dim, -dim, dim, dim]), opacity=0.5, leg = false)
  else 
    top_left_center = center - half_dim + half_dim*1im
    top_right_center = center + half_dim + half_dim*1im
    bot_left_center = center - half_dim - half_dim*1im 
    bot_right_center = center + half_dim - half_dim*1im
    displayQuadtreeBoxes(plot, tree.children[1], top_left_center, half_dim)
    displayQuadtreeBoxes(plot, tree.children[2], top_right_center, half_dim)
    displayQuadtreeBoxes(plot, tree.children[3], bot_left_center, half_dim)
    displayQuadtreeBoxes(plot, tree.children[4], bot_right_center, half_dim)
  end
end

function displayQuadtreeBoxesWrapper(array, center, dim, depth)
  gr(size=(1000, 1000))
  tree = buildQuadtree!(array, nothing, center, dim, 1, length(array), depth)
  p = plot(-1:1, -1:1)
  displayQuadtreeBoxes(p, tree, center, dim)
end

""" Display quadtree boxes & particles helper"""
function displayQuadtreeBoxesAndParticles(plot, array, tree, center, dim)
  half_dim = dim/2
  if tree.depth == 0
    plot!(plot, Shape(real(center) .+ [-dim, dim, dim, -dim], imag(center) .+ [-dim, -dim, dim, dim]), opacity=0.5, leg = false)
    if tree.last >= tree.first # if there is more than one particle 
      points = array[tree.first:tree.last]
      x = real.(points)
      y = imag.(points)
      scatter!(plot, (x, y), legend=false, label="")
    end
  else 
    top_left_center = center - half_dim + half_dim*1im
    top_right_center = center + half_dim + half_dim*1im
    bot_left_center = center - half_dim - half_dim*1im 
    bot_right_center = center + half_dim - half_dim*1im
    displayQuadtreeBoxesAndParticles(plot, array, tree.children[1], top_left_center, half_dim)
    displayQuadtreeBoxesAndParticles(plot, array, tree.children[2], top_right_center, half_dim)
    displayQuadtreeBoxesAndParticles(plot, array, tree.children[3], bot_left_center, half_dim)
    displayQuadtreeBoxesAndParticles(plot, array, tree.children[4], bot_right_center, half_dim)
  end
end

function displayQuadtreeBoxesAndParticlesWrapper(array, center, dim, depth)
  gr(size=(1000, 1000))
  tree = buildQuadtree!(array, nothing, center, dim, 1, length(array), depth)
  p = plot()
  displayQuadtreeBoxesAndParticles(p, array, tree, center, dim)
  return p
end

""" Display quadtree boxes, particles, mass helper"""
function displayQuadtreeBoxesParticlesMass(plot, array, tree, center, dim)
  half_dim = dim/2
  if tree.depth == 0
    plot!(plot, Shape(real(center) .+ [-dim, dim, dim, -dim], imag(center) .+ [-dim, -dim, dim, dim]), opacity=0.5, leg = false)
    if tree.last > tree.first # if there is more than one particle 
      points = first.(array[tree.first:tree.last])
      x = real.(points)
      y = imag.(points)
      masses = getfield.(array[tree.first:tree.last], 2)
      scatter!(plot, (x, y), markersize=10*masses, legend=false, label="")
    elseif tree.last == tree.first # one particle
      point = first(array[tree.first])
      x = real(point)
      y = imag(point)
      mass = getfield(array[tree.first], 2)
      scatter!(plot, (x, y), markersize=10*mass, legend=false, label="")
    end
  else 
    top_left_center = center - half_dim + half_dim*1im
    top_right_center = center + half_dim + half_dim*1im
    bot_left_center = center - half_dim - half_dim*1im 
    bot_right_center = center + half_dim - half_dim*1im
    displayQuadtreeBoxesParticlesMass(plot, array, tree.children[1], top_left_center, half_dim)
    displayQuadtreeBoxesParticlesMass(plot, array, tree.children[2], top_right_center, half_dim)
    displayQuadtreeBoxesParticlesMass(plot, array, tree.children[3], bot_left_center, half_dim)
    displayQuadtreeBoxesParticlesMass(plot, array, tree.children[4], bot_right_center, half_dim)
  end
end

function displayQuadtreeBoxesParticlesMassWrapper(array, center, dim, depth)
  tree = buildQuadtree!(array, nothing, center, dim, 1, length(array), depth)
  p = plot()
  displayQuadtreeBoxesParticlesMass(p, array, tree, center, dim)
  return p
end
function countNumberOfPointsInQuadtree(tree) if tree.depth == 0
    if tree.last >= tree.first
      return tree.last - tree.first + 1
    else
      return 0
    end
  else
    return sum(map(subtree -> countNumberOfPointsInQuadtree(subtree), tree.children))
  end
end

function propagateMassUp!(tree, array)
  println("In propagate function")
  if tree.depth == 0
    println("base case")
    masses = getfield.(array[tree.first:tree.last], 2)
    println(masses)
    tree.massExpansion = sum(masses)
    return tree.massExpansion
  else
    println("non base case")
    #tree.massExpansion = mapreduce((subtree) -> propagateMassUp!(subtree, array), +, tree.children) 
    tree.massExpansion = sum([propagateMassUp!(subtree, array) for subtree in tree.children]) 
    return tree.massExpansion
  end
end
