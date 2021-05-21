
""" Plot accuracy of Barnes Hut """

const G = 1e-2
# softening paramter to avoid a -> inf when particles close
const S = 1e-32
# timestep
const Δt = 1e-2
# num timesteps
const TIMESTEPS = 1
# number of bodies
const N = 10


include("BarnesHut.jl")


abs_error = Array{Float64, 1}(undef, length(Ps))
rel_error = Array{Float64, 1}(undef, length(Ps))


prev_points = rand(Float64, 2, N)
tangent_velocity = similar(prev_points)
tangent_velocity[1, :] = prev_points[2, :] .- 0.5
tangent_velocity[2, :] = -prev_points[1, :] .- 0.5
curr_points = prev_points .+ 1e-2 .* tangent_velocity
masses = 0.9*rand(Float64, N) .+ 1

particles = [Particle(curr_points[:, i], prev_points[:, i], masses[i]) for i in 1:N]

tree = generate_tree(particles)
for particle in particles
  tmp = copy(particle.pos)
  a = net_acc(particle, tree, θ_sq, acc_func)
  particle.pos .= 2*particle.pos - particle.prev_pos + G * a * Δt^2
  particle.prev_pos .= tmp
end
