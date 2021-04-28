"""

Naive O(n^2) force calculation.
Stores in arrays and uses builtin cumsum.

"""

using Plots
gr() # GR backend
# can compare Makie for speed


# gravitational constant
const G = 6.67408e-11
# softening paramter to avoid a -> inf when particles close
const S = 1e-10
# timestep
const Δt = 1e-3
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
      Δ = p .- q
      # TODO: use cumsum
      acc[i, :] += grav_constant * Δ * mass[j] / squared_distance(p, q)
    end
  end
end


# initalize random position and random previous position
prev_pos = rand(Float64, (N, 2))
curr_pos = rand(Float64, (N, 2))
mass = rand(Float64, (N, 2))
acc = Array{Float64}(undef, N, 2)

# simulation loop
anim = Animation()
for i ∈ 1:TIMESTEPS
  updateAcc!(acc, curr_pos, mass)
  # CHECK/VERIFY: verlet integration
  curr_pos .= 2*curr_pos - prev_pos + acc * Δt^2
  x = @view curr_pos[:, 1]
  y = @view curr_pos[:, 2]
  scatter!(x, y, xlim = (-10, 10), ylim = (-10, 10))
  frame(anim)
end

fig(anim, "n_body.gif", fps = 15)
