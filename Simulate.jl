
""" Simulate the N-body problem using the Fast Multipole Method """


using Plots
using Random


include("Quadtree.jl")
include("Multipole.jl")


# scaled gravitational constant
# normal value is 6.67408e-11
const G = 1e-2
# softening paramter to avoid a -> inf when particles close
const S = 1e-10
# timestep
const Δt = 1e-1
# num timesteps
const TIMESTEPS = 20
# number of bodies
const N = 3
# number of past positions saved
const NUM_PAST_POSITIONS = 5
# depth of tree to construct
const TREE_DEPTH = 3
# number of terms to keep in multipole expansion
const P = 33

function runSimulation(quadtree, pos_memory, masses, ω_p, timesteps=TIMESTEPS, num_past_positions=NUM_PAST_POSITIONS)
  gr(reuse=true, size = (1000, 1000))
  lim = (-0.05, 1.05)
  # index to circular buffer
  next_idx = 3
  curr_idx = 2
  prev_idx = 1

  for i ∈ 1:timesteps
    # NOTE: make sure updated 
    #println(pos_memory)
    curr_points = @view pos_memory[:, curr_idx]
    prev_points = @view pos_memory[:, prev_idx]
   
    # OPTIMIZE, can zero within Multipole function
    ω_p .= zero(ω_p[1])

    # have to four color sort both current and previous points according
    # to current so that have proper alignment
    # another reason why position based verlet integration is good, as otherwise would
    # have to sort an additional array, velocity, as well
    # FMM
    updateQuadtreePointMasses(quadtree, curr_points, masses)
    # upward pass
    P2M(quadtree, curr_points, masses)
    M2M(quadtree)
    # downward pass
    M2L(quadtree)
    L2L(quadtree)
    L2P(quadtree, curr_points, ω_p)
    NNC(quadtree, curr_points, masses, ω_p)

    #VERIFY: verlet integration
    #OPTIMIZE : column layout, can reinterpte complex as two reals

    #println(ω_p)
    
    pos_memory[:, next_idx] .= 2*curr_points - prev_points + G * (-real(ω_p) + imag(ω_p)*1im) * Δt^2
    #pos_memory[:, next_idx] .= 2*curr_points - prev_points - G * imag(ω_p) * Δt^2
    

    # threshold so don't go out of bounds, can do this faster if reinterpret
    flattened = reinterpret(Float64, pos_memory[:, next_idx])
    flattened .= min.(flattened, 1.0) 
    pos_memory[:, next_idx] .= reinterpret(ComplexF64, max.(flattened, 0.0))

    if (i > 15)
      println("timestep")
      println(G * (-real(ω_p) + imag(ω_p)*1im))
      println(pos_memory[:, next_idx])
    end

    prev_idx = curr_idx
    curr_idx = next_idx
    next_idx = mod1(next_idx+1, num_past_positions)

    @assert curr_points != prev_points
    x = real(curr_points)
    y = imag(curr_points)
    #println(masses)
    # OPTIMIZE with reinterpreting / reshaping or storing differently
    # Can no longer add trailing affect because of sorting in place
    #plot((real(transpose(pos_memory)), imag(transpose(pos_memory))), xlim = lim, ylim = lim, color=:black, label="", legend=false)
    scatter((x, y), xlim = lim, ylim = lim, legend=false, markerstrokewidth=0, markersize=7*masses, color=:black, label="")
    gui() 
  end

end


""" Initialize Simulation """

quadtree = buildQuadtree(TREE_DEPTH, P)
position_memory = Array{ComplexF64}(undef, N, NUM_PAST_POSITIONS)
# initialize position at t=-dt and t=0 for verlet integration

Random.seed!(1)

position_memory[:, 1] = [0.1+0.7im, 0.3+0.8im, 0.5+0.9im]
position_memory[:, 2] = [0.11+0.71im, 0.31+0.81im, 0.51+0.91im]

#position_memory[:, 1]  = rand(ComplexF64, N)
#position_memory[:, 2]  = position_memory[:, 1] + 1e-10*rand(ComplexF64, N)
masses         = 0.5*rand(Float64, N) .+ 1
ω_p            = Array{ComplexF64, 1}(undef, N)

println("initial position")
show(stdout, "text/plain", position_memory[:, 1])
println()

runSimulation(quadtree, position_memory, masses, ω_p)

