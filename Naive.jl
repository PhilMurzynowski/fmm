"""

Naive O(n^2) force calculation.
Stores in arrays and uses builtin cumsum.

"""

using Plots
# can compare Makie for speed


# scaled gravitational constant
# normal value is 6.67408e-11
const G = 1e-2
# softening paramter to avoid a -> inf when particles close
const S = 1e-10
# timestep
const Δt = 1e-1
# num timesteps
const TIMESTEPS = 100
# number of bodies
const N = 10


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

function simulate(curr_pos, prev_pos, tmp_pos, mass, acc)
  """
  # using Plots + GR to create a gif
  gr() # GR backend
  anim = Animation()
  for i ∈ 1:TIMESTEPS
    updateAcc!(acc, curr_pos, mass)
    # CHECK/VERIFY: verlet integration
    tmp_pos = curr_pos
    curr_pos .= 2*curr_pos - prev_pos + acc * Δt^2
    prev_pos = tmp_pos
    assert(curr_pos != prev_pos)
    x = @view curr_pos[:, 1]
    y = @view curr_pos[:, 2]
    scatter(x, y, xlim = (-10, 10), ylim = (-10, 10))
    frame(anim)
  end
  gif(anim, "fmm/n_body.gif", fps = 5)
  """

  # using Plots + GR realtime display
  #gr(show=true)
  gr(reuse=true, size = (1000, 1000))
  x = @view curr_pos[:, 1]
  y = @view curr_pos[:, 2]
  plot = scatter(x, y, xlim = (-2, 10), ylim = (-2, 2), legend=false)
  for i ∈ 1:TIMESTEPS
    updateAcc!(acc, curr_pos, mass)
    # CHECK/VERIFY: verlet integration
    tmp_pos = copy(curr_pos)
    curr_pos .= 2*curr_pos - prev_pos + acc * Δt^2
    prev_pos = tmp_pos
    @assert curr_pos != prev_pos
    x = @view curr_pos[:, 1]
    y = @view curr_pos[:, 2]
    scatter!(plot, x, y, xlim = (-2, 2), ylim = (-2, 2), legend=false)
    gui() 
    sleep(.1)
  end
end

# initalize random position and random previous position
prev_pos = rand(Float64, (N, 2))
curr_pos = prev_pos + 0.001*rand(Float64, (N, 2))
tmp_pos = similar(curr_pos) # preallocate
mass = ones(N) # masses normalized with respect to gravitational constant
acc = Array{Float64}(undef, N, 2)

# simulation loop
simulate(curr_pos, prev_pos, tmp_pos, mass, acc)
