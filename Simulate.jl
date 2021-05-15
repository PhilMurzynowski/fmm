
""" Simulate the N-body problem using the Fast Multipole Method """


using Plots


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
const TIMESTEPS = 1000
# number of bodies
const N = 10
# number of past positions saved
const NUM_PAST_POSITIONS = 5
# depth of tree to construct
const TREE_DEPTH = 3
# number of terms to keep in multipole expansion
const P = 30


function runSimulation(quadtree, pos_memory, masses, ω_p, timesteps=TIMESTEPS, num_past_positions=NUM_PAST_POSITIONS)
  gr(reuse=true, size = (1000, 1000))
  # index to circular buffer
  next_idx = 3
  curr_idx = 2
  prev_idx = 1

  for i ∈ 1:timesteps
    # NOTE: make sure updated 
    curr_points = @view pos_memory[:, curr_idx]
    prev_points = @view pos_memory[:, prev_idx]
   
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
    #pos_memory[next_idx, :, :] .= 2*curr_pos - prev_pos + acc * Δt^2
    pos_memory[:, next_idx] .= 2*curr_points - prev_points + real(ω_p) * Δt^2
    pos_memory[:, next_idx] .= 2*curr_points - prev_points - imag(ω_p) * Δt^2

    prev_idx = curr_idx
    curr_idx = next_idx
    next_idx = mod1(next_idx+1, num_past_positions)

    #@assert curr_points != prev_points
    plot_points = @view pos_memory[:, prev_idx]
    x = real(plot_points)
    y = imag(plot_points)
    # OPTIMIZE with reinterpreting / reshaping or storing differently
    #plot((real(pos_memory), imag(pos_memory)), xlim = (-2, 2), ylim = (-2, 2), color=:black, label="", legend=false)
    scatter!((x, y), xlim = (-2, 2), ylim = (-2, 2), legend=false, markerstrokewidth=0, markersize=7*masses, color=:black, label="")
    gui() 
  end

end


""" Initialize Simulation """

quadtree = buildQuadtree(TREE_DEPTH, P)
position_memory = Array{ComplexF64}(undef, N, NUM_PAST_POSITIONS)
# initialize position at t=-dt and t=0 for verlet integration
position_memory[:, 1]  = rand(ComplexF64, N)
position_memory[:, 2]  = rand(ComplexF64, N)
masses         = rand(Float64, N)
ω_p            = Array{ComplexF64, 1}(undef, N)

runSimulation(quadtree, position_memory, masses, ω_p)

