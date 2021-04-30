"""

Naive O(n^2) force calculation.
Stores in arrays and uses builtin cumsum.

"""

using Plots
using ProfileView
# can compare Makie for speed


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
const N = 40
# number of past positions saved
const NUM_PAST_POSITIONS = 5


function squared_distance(p, q, soft=S)
  Δx = (p[1] - q[1])
  Δy = (p[2] - q[2])
  return Δx^2 + Δy^2 + soft
end

# update acceleration
function updateAcc!(acc, pos, mass, grav_constant=G)
  # clear previously computed forces
  fill!(acc, zero(acc[1]))
  # calculated force contribution on each body from each other body
  for i ∈  1:size(pos)[1]
    for j ∈  1:size(pos)[1] 
      p = @view pos[i, :]
      q = @view pos[j, :]
      Δ = q .- p
      # TODO: use cumsum
      acc[i, :] += grav_constant * Δ * mass[j] / squared_distance(p, q)
    end
  end
end


function simulate(pos_memory, mass, acc, timesteps=TIMESTEPS, num_past_positions=NUM_PAST_POSITIONS)
  gr(reuse=true, size = (1000, 1000))
  # index to circular buffer
  next_idx = 3
  curr_idx = 2
  prev_idx = 1
  for i ∈ 1:timesteps
    curr_pos = @view pos_memory[curr_idx, :, :]
    prev_pos = @view pos_memory[prev_idx, :, :]
    updateAcc!(acc, curr_pos, mass)
    # CHECK/VERIFY: verlet integration
    pos_memory[next_idx, :, :] .= 2*curr_pos - prev_pos + acc * Δt^2

    prev_idx = curr_idx
    curr_idx = next_idx
    next_idx = mod1(next_idx+1, num_past_positions)

    @assert curr_pos != prev_pos
    x = @view pos_memory[curr_idx, :, 1]
    y = @view pos_memory[curr_idx, :, 2]
    plot((pos_memory[:, :, 1], pos_memory[:, :, 2]), xlim = (-2, 2), ylim = (-2, 2), color=:black, label="", legend=false)
    scatter!((x, y), xlim = (-2, 2), ylim = (-2, 2), legend=false, markerstrokewidth=0, markersize=7*mass, color=:black, label="", )
    gui() 
  end
end

# use a circular buffer to hold NUM_PAST_POSITIONS
# initalize random position and random previous position
pos_memory = Array{Float64}(undef, NUM_PAST_POSITIONS, N, 2)
# initialize position at t=-dt and t=0 for verlet integration
pos_memory[1, :, :] = rand(Float64, (N, 2))
pos_memory[2, :, :] = pos_memory[1, :, :] + 1e-5*rand(Float64, (N, 2))
mass = 0.5*rand(Float64, N) .+ 1 # masses normalized with respect to gravitational constant
mass[1] *= 10; # make one mass large so keeps things mostly in place
acc = Array{Float64}(undef, N, 2)

# simulation loop
simulate(pos_memory, mass, acc)

# profile
#@profview simulate(curr_pos, prev_pos, tmp_pos, mass, acc, 1) # run once to compile
#@profview simulate(curr_pos, prev_pos, tmp_pos, mass, acc, TIMESTEPS)
