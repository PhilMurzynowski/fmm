
"""
    Modified Julia code for the Barnes-Hut algorithm for comparison against the FMM
    Credit to: https://github.com/alexhad6/ParallelBarnesHut.jl

    Messily, oh very messily, ported to serial, two dimensional code for purposes of flopcount comparison to FMM using GFlops.jl.
    Converting to serial may not be necessary, but easier to debug and want to guarantee nothing is affecting
    flop count.

    This code is actually not that optimized turns out, despite being parallelized, flop count will be much better measurement of performance than runtime.
"""

# Imports
using Plots
using Random

""" NEW, added these constants """
# softening paramter to avoid a -> inf when particles close
const S = 1e-32

# Exports

"""
export Particle
export simulate!, grav_acc
export show_particles, show_quadview, show_boxes, animate_frames
export rand_particles
"""


# Structs
#
""" NEW 
Scrapped, am keeping points represented as arrays as don't feel like
writing broadcasting atm
struct Point
  x::Float64
  y::Float64
end
"""

"""
Represents a particle in a Barnes Hut simulation. For a galaxy simulation, a
particle is a star or a group of stars.
MODIFIED:
  to hold previous position instead of velocity
  also in 2D not 3d
"""
mutable struct Particle
    "Position coordinates of a particle (in meters)."
    pos::Array{Float64, 1}
    prev_pos::Array{Float64, 1}

    "Mass of a particle (in kilograms)."
    mass::Float64
end

"""
Used in the octree for the Barnes Hut algorithm. Represents a cubical octant of
space and stores the mass and center of mass of the contained particles.
"""



"""
MODIFIED to use Point
"""
struct Node
    "Side length of the cubical region of space represented by a node."
    size::Float64
    
    "Total mass of the particles contained in the region represented by a node."
    total_mass::Float64
    
    "Center of mass of the particles contained in the region represented by a node."
    center_of_mass::Array{Float64, 1}
    
    "Children of a node in the octree."
    children::Array{Union{Node, Nothing},1}
end


# Barnes Hut Tree Generation

"""
    generate_tree(particles; return_boxes = false)
Generate the Barnes Hut quadtree for a given array of particles. Optionally
return coordinates of the boxes represented by the nodes in the octree.

ELIMINATED return boxes for simplicity
"""
#function generate_tree(particles::Array{Particle,1}, return_boxes::Bool = false)
function generate_tree(particles::Array{Particle,1})
    # Calculate bounds of a box enclosing all particles
    # SKIP as for fair comparison with FMM will do box with corners (0, 0) & (1, 1)
    # As FMM already assumes a space of that size
    #minx, maxx = minimum(p -> p.pos[1], particles), maximum(p -> p.pos[1], particles)
    #miny, maxy = minimum(p -> p.pos[2], particles), maximum(p -> p.pos[2], particles)
    # cut out minz, maxz
    
    # Calculate corner and size of a cube enclosing all particles
    #size = max(maxx - minx, maxy - miny, maxz - minz)
    size = 1.0
    #corner = [minx, miny, minz]
    corner = [0.0, 0.0]
    
    #if return_boxes
    #    # Recursively compute tree of particles and return coordinates of boxes
    #    boxes = Tuple{Float64,Array{Float64,1}}[]
    #    # SERIAL only
    #    #boxes_threaded = [Tuple{Float64,Array{Float64,1}}[] for i = 1:Threads.nthreads()]
    #    #generate_tree_helper(particles, size, corner, return_boxes, boxes_threaded)
    #    generate_tree_helper(particles, size, corner, return_boxes, boxes)
    #    #for thread = eachindex(boxes_threaded)
    #    #   append!(boxes, boxes_threaded[thread])
    #    #end
    #    boxes
    #else
    #    # Recursively compute tree of particles and return that tree
    #    generate_tree_helper(particles, size, corner, false, Array{Tuple{Float64,Array{Float64,1}},1}[])
    #end
    generate_tree_helper(particles, size, corner)
end

"""
    generate_tree_helper(particles, size, corner, save_boxes, boxes)
Generate subtree of Barnes Hut quadtree given a list of particles and the corner
size of a square region of space. Optionally save coordinates of this region in
the array `boxes`. Helper function for `generate_tree`.

ELIMINATED save_boxes functionality for simplicity along with boxes as that was the only usecase
"""
function generate_tree_helper(particles::Array{Particle,1}, size::Float64, corner::Array{Float64,1})
    # If saving boxes, save current box coordinates
    #if save_boxes
    #    push!(boxes_threaded[Threads.threadid()], (size, corner))
    #end
    #println(size)
    
    if length(particles) == 0
        # Base case 1: if no particles in quadrant, return empty leaf node
        Node(size/2, 0, corner .+ size/2, [nothing for i = 1:4])
    elseif length(particles) == 1
        # Base case 2: if one particle in octant, return leaf node with that particle
        Node(size/2, particles[1].mass, particles[1].pos, [nothing for i = 1:4])
    else
        # Recursive case
        
        # Calculate total mass and center of mass
        total_mass = sum(particle.mass for particle in particles)
        center_of_mass = sum(particle.mass .* particle.pos for particle in particles) / total_mass
        
        # Sort particles into 4 quadrants (not using threading)
        quadrants = [Particle[] for quad = 1:4]
        #octants_threaded = [[Particle[] for oct = 1:8] for thread = 1:Threads.nthreads()]
        new_corners = vec([corner + [i, j] for i = (0, size/2), j = (0, size/2)])
        #@sync for particle in particles
        for particle in particles
            #Threads.@spawn begin
            begin
                quad_num = 0
                for quad = 1:4
                    if all(new_corners[quad] .<= particle.pos .<= new_corners[quad] .+ size/2)
                        quad_num = quad
                    end
                end
                #if oct_num != 0
                #    push!(octants_threaded[Threads.threadid()][oct_num], particle)
                #end
                # SERIAL EQUIVALENT to commented out sections above and below
                push!(quadrants[quad_num], particle)
            end
        end
        #for thread = eachindex(octants_threaded)
        #    for oct = 1:8
        #       append!(octants[oct], octants_threaded[thread][oct])
        #    end
        #end
        
        # Calculate child nodes recursively (NOT using threading)
        children::Array{Union{Node, Nothing},1} = fill(nothing, 4)
        #@sync for oct = 1:8
        for quad = 1:4
            #Threads.@spawn children[oct] = generate_tree_helper(octants[oct], size/2, new_corners[oct],
            #                                                    save_boxes, boxes_threaded)
            children[quad] = generate_tree_helper(quadrants[quad], size/2, new_corners[quad])
        end
        
        # Return node
        Node(size, total_mass, center_of_mass, children)
    end
end

"""
    net_acc(particle, node, θ, acc_func)
Recursively calculate the net acceleration on a particle given the Barnes Hut
quadtree, a threshold `θ`, and the acceleration function `acc_func`.
"""
function net_acc(particle::Particle, node::Union{Node, Nothing}, θ_squared::Float64, acc_func)
    # Recusively calculates net acceleration on a particle given a particle tree
    if node == nothing
        # Base case 1: if tree is empty, return 0 acceleration
        #zeros(Float64, 3)
        zeros(Float64, 2)
    else
        # Recursive case, Barnes Hut method
        s = node.size                            # Width of the square that this node represents
        r = node.center_of_mass .- particle.pos  # Vector from particle to node center of mass
        # couldn't resist fixing this
        #if s/√sum(r.^2) < θ
        if s*s/sum(r.^2) < θ_squared
            # Base case 2: if node size divided by distance is less than threshold θ,
            #              calculate gravitational acceleration
            acc_func(node.total_mass, r)
        else
            # Recursive case
            sum(net_acc(particle, child, θ_squared, acc_func) for child in node.children)
        end
    end
end

"""
Calculate the gravitation acceleration on a particle given the mass of another
object and the distance vector from particle to object. The usual gravitational
force is softened by a factor `ϵ` to prevent it from blowing up when the
distance is very small.

MODIFIED: using global constants, G, and S for gravitational constant and softening parameter respectively
"""
function grav_acc(mass::Float64, r::Array{Float64,1}; ϵ::Float64 = 0.02)
    # Calculate gravitational acceleration, using node center of mass
    # and softening factor ϵ (to prevent blow up at d2 = 0)
    #G::Float64 = 6.67430e-11
    #((G * mass) / ((sum(r.^2)+ϵ^2)^(3/2))) .* r
    # don't multiply by G here for consistency with FMM implementation
    #mass ./ (r .+ S) # this is not right, was only correct with complex numbers
    mass ./ (sum(r.^2) .+ S) .* r
end

""" New function """
function barnesHutUpdate!(acc, particles, θ_squared)
  tree = generate_tree(particles)
  for i in 1:length(particles)
    particle = particles[i]
    a = net_acc(particle, tree, θ_squared, grav_acc)
    acc[:, i] .= a
  end
end

"""
    simulate!(particles, steps; Δt = 0.1, θ = 0.5, acc_func = grav_acc)
Approximately solve N-body problem for an array of particles for number of time
steps into the future using the Barnes Hut method. Return an array of particle
arrays, where each particle array is the state of all particles at each time
step.
# Arguments
- `particles::Array{Particle,1}`: an array of particles.
- `steps::Int64`: the number of time steps to run the simulation.
- `Δt::Float64 = 0.1`: length of time step in seconds (keyword argument).
- `θ::Float64 = 0.5`: node size to particle distance threshold for Barnes Hut
  simulation (keyword argument).
- `acc_func(mass::Float64, r::Array{Float64,1}) = grav_acc`: calculates particle
  acceleration from mass and distance vector.
function simulate!(particles::Array{Particle,1}, steps::Int64;
                   Δt::Float64 = 0.1, θ::Float64 = 0.5, acc_func = grav_acc)
    # Initialize array of snapshots of particles at each time step
    #frames = Array{Particle,1}[]
    
    # Store initial conditions as first frame
    #push!(frames, deepcopy(particles))
    
    # Advance velocity by half time step to do leapfrog integration method
    #vel_step!(particles, generate_tree(particles), Δt/2, θ, acc_func)
    
    # Simulate steps and save each frame
    for i = 1:steps
        step!(particles, generate_tree(particles), Δt, θ, acc_func)
        #push!(frames, deepcopy(particles))
    end
    
    # Return frames
    #frames
end
"""

"""
    step!(particles, tree, Δt, θ, acc_func)
Calculate the acceleration of each particle using `net_acc` and approximate new
velocity and position using `Δt`.
function step!(particles::Array{Particle,1}, tree::Node, Δt::Float64, θ::Float64, acc_func)
    # Calculate acceleration, and approximately advance velocity and position (using threading)
    #@sync for particle in particles
    #    Threads.@spawn begin
    for particle in particles
      tmp = copy(particle.pos)
      particle.vel += net_acc(particle, tree, θ, acc_func) * Δt
      particle.pos += particle.vel * Δt
    end
end
"""
