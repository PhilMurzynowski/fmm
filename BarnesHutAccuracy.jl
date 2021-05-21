
""" Plot accuracy of Barnes Hut """

const G = 1e-2
# softening paramter to avoid a -> inf when particles close
const S = 1e-32
# timestep
const Δt = 1e-2
# num timesteps
const TIMESTEPS = 1
# number of bodies
const N = 1000


include("BarnesHut.jl")


abs_error = Array{Float64, 1}(undef, length(Ps))
rel_error = Array{Float64, 1}(undef, length(Ps))



  tree = generate_tree(particles)
  for particle in particles
    tmp = copy(particle.pos)
    a = net_acc(particle, tree, θ, acc_func)
    particle.pos .= 2*particle.pos - particle.prev_pos + G * a * Δt^2
    particle.prev_pos .= tmp
  end
