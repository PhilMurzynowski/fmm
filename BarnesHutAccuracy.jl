
""" Plot accuracy of Barnes Hut """

const G = 1e-2
# softening paramter to avoid a -> inf when particles close
const S = 1e-32
# timestep
const Δt = 1e-2
# num timesteps
const TIMESTEPS = 1
# number of bodies
const N = 2


include("BarnesHut.jl")


num_tests = 5


abs_error = Array{Float64, 1}(undef, length(num_tests))
rel_error = Array{Float64, 1}(undef, length(num_tests))


prev_points = rand(Float64, 2, N)
tangent_velocity = similar(prev_points)
tangent_velocity[1, :] = prev_points[2, :] .- 0.5
tangent_velocity[2, :] = -prev_points[1, :] .- 0.5
curr_points = prev_points .+ 1e-2 .* tangent_velocity
prev_points .= max.(min.(prev_points, 1.0), 0)
curr_points .= max.(min.(curr_points, 1.0), 0)
masses = 0.9*rand(Float64, N) .+ 1
particles = [Particle(curr_points[:, i], prev_points[:, i], masses[i]) for i in 1:N]

θ_sq = 0.25

gr(reuse=true, size = (1000, 1000))
println(curr_points)
scatter(curr_points[1, :], curr_points[2, :], xlim=lim, ylim=lim, legend=false, markerstrokewidth=0, markersize=7*masses, color=:black, label="")
lim = (-0.05, 1.05)

for i ∈ 1:TIMESTEPS
  tree = generate_tree(particles)
  for particle in particles
    tmp = copy(particle.pos)
    a = net_acc(particle, tree, θ_sq, grav_acc)
    particle.pos .= 2*particle.pos - particle.prev_pos + G * a * Δt^2
    particle.prev_pos .= tmp
  end
end


gui() 
