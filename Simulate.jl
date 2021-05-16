
""" Simulate the N-body problem using the Fast Multipole Method """


using Plots
using Random
using BenchmarkTools


include("Quadtree.jl")
include("Multipole.jl")
include("BinomialTable.jl")

# scaled gravitational constant
# normal value is 6.67408e-11
const G = 1e-2
# softening paramter to avoid a -> inf when particles close
const S = 1e-32
# timestep
const Δt = 1e-2
# num timesteps
const TIMESTEPS = 1
# number of bodies
const N = 1000
# number of past positions saved
const NUM_PAST_POSITIONS = 3 # do not support plotting history yet, set to 3
# depth of tree to construct
const TREE_DEPTH = 3
# number of terms to keep in multipole expansion
const P = 32 # keeping multiple of 4 for easier vectorization, actually more complicated as often need P+1

function runSimulation(quadtree, pos_memory, masses, ω_p, timesteps=TIMESTEPS, num_past_positions=NUM_PAST_POSITIONS)
  gr(reuse=true, size = (1000, 1000))
  lim = (-0.05, 1.05)
  # index to circular buffer
  next_idx = 3
  curr_idx = 2
  prev_idx = 1

  binomial_table = binomialTable(P)
  binomial_table_t = binomialTableTransposed(P)
  preallocated_size = floor(Int, N/2)
  preallocated_mtx::Array{ComplexF64, 2} = Array{ComplexF64, 2}(undef, preallocated_size, preallocated_size)

  for i ∈ 1:timesteps
    curr_points = @view pos_memory[:, curr_idx]
    prev_points = @view pos_memory[:, prev_idx]
   
    # TUNE, is it faster to zero within Multipole function, and which
    # ω_p .= zero(ω_p[1])

    # have to four color sort both current and previous points according
    # to current so that have proper alignment
    # another reason why position based verlet integration is good, as otherwise would
    # have to sort an additional array, velocity, as well
    updateQuadtreePointMasses(quadtree, curr_points, masses, prev_points)
    # FMM
    FMM!(quadtree, curr_points, masses, ω_p, binomial_table, binomial_table_t, preallocated_mtx)
    #@btime FMM!(quadtree, $curr_points, $masses, $ω_p, $binomial_table, $binomial_table_t, $preallocated_mtx)

    # TEST
    """
    println("Testing")
    forces = [real(ω_p), -imag(ω_p)]
    correct_ω_p = similar(ω_p)
    correct_ω_p .= zero(correct_ω_p[1])
    kernel_mtx::Array{ComplexF64, 2} = 1 ./ (transpose(curr_points) .- curr_points)
    foreach(i -> kernel_mtx[i, i] = zero(kernel_mtx[1, 1]), 1:length(curr_points))
    #show(stdout, "text/plain", kernel_mtx)
    # can't do simple dot product unfortunately
    correct_ω_p = vec(sum(kernel_mtx.*masses, dims=1))
    correct_forces = [real(correct_ω_p), -imag(correct_ω_p)]
    #println(correct_forces)
    #println(forces)
    @assert correct_forces ≈ forces 
    """

    #VERIFY: verlet integration
    #OPTIMIZE : column layout, can reinterpte complex as two reals

    pos_memory[:, next_idx] .= 2*curr_points - prev_points + G * (-real(ω_p) + imag(ω_p)*1im) * Δt^2
    

    # Threshold so don't go out of bounds.
    # OPTIMIEZE perhaps can do this faster if reinterpret.
    flattened = reinterpret(Float64, pos_memory[:, next_idx])
    flattened .= min.(flattened, 1.0) 
    pos_memory[:, next_idx] .= reinterpret(ComplexF64, max.(flattened, 0.0))

    #if (i == TIMESTEPS)
    #  println("timestep")
    #  println(G * (-real(ω_p) + imag(ω_p)*1im))
    #  println(pos_memory[:, next_idx])
    #  println(pos_memory[:, curr_idx])
    #  println(pos_memory[:, prev_idx])
    #end

    prev_idx = curr_idx
    curr_idx = next_idx
    next_idx = mod1(next_idx+1, num_past_positions)

    x = real(curr_points)
    y = imag(curr_points)
    # OPTIMIZE with reinterpreting / reshaping or storing differently
    
    # GRAPHICS: Can no longer add trailing affect because of sorting in place.
    # Perhaps in future update.
    # plot((real(transpose(pos_memory)), imag(transpose(pos_memory))), xlim = lim, ylim = lim, color=:black, label="", legend=false)
    
    scatter((x, y), xlim = lim, ylim = lim, legend=false, markerstrokewidth=0, markersize=7*masses, color=:black, label="")
    #scatter((x, y), xlim = lim, ylim = lim, legend=false, markerstrokewidth=0, markersize=7*masses, label="")
    gui() 
  end

end


""" Initialize Simulation """

quadtree = buildQuadtree(TREE_DEPTH, P)
position_memory = Array{ComplexF64}(undef, N, NUM_PAST_POSITIONS)
# initialize position at t=-dt and t=0 for verlet integration

Random.seed!(2)
#Random.seed!(3)

#position_memory[:, 1] = [0.1+0.7im, 0.3+0.8im, 0.5+0.9im]
#position_memory[:, 2] = [0.11+0.71im, 0.31+0.81im, 0.51+0.91im]

position_memory[:, 1]  = rand(ComplexF64, N)
# Give some small random initial velocity
# position_memory[:, 2]  = position_memory[:, 1] + 1e-3*rand(ComplexF64, N)
# Give some tangent velocity
# NOTE: this is assuming 1.0 x 1.0 space
tangent_velocity = imag(position_memory[:, 1] .- (0.5 + 0.5im)) .- 1im*real(position_memory[:, 1] .- (0.5 + 0.5im))
position_memory[:, 2]  = position_memory[:, 1] .+ 1e-2*tangent_velocity

masses         = 0.9*rand(Float64, N) .+ 1
ω_p            = Array{ComplexF64, 1}(undef, N)

#println("initial position")
#show(stdout, "text/plain", position_memory[:, 1])
#println()

runSimulation(quadtree, position_memory, masses, ω_p)

