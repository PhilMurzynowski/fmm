
""" Plot accuracy of Barnes Hut 
    Inefficient Plotting pipeline but thats ok, just testing """

const G = 1e-2
# softening paramter to avoid a -> inf when particles close
const S = 1e-32
# timestep
const Δt = 1e-2
# num timesteps
const TIMESTEPS = 2
# number of bodies
const N = 2


include("BarnesHut.jl")


prev_points = rand(Float64, 2, N)
tangent_velocity = similar(prev_points)
tangent_velocity[1, :] = prev_points[2, :] .- 0.5
tangent_velocity[2, :] = -(prev_points[1, :] .- 0.5)
#tangent_velocity .= zero(Float64)
curr_points = prev_points .+ 1e-2 .* tangent_velocity
prev_points .= max.(min.(prev_points, 1.0), 0)
curr_points .= max.(min.(curr_points, 1.0), 0)
masses = 0.9*rand(Float64, N) .+ 1
acc = similar(curr_points)
particles = [Particle(curr_points[:, i], prev_points[:, i], masses[i]) for i in 1:N]

θ_sq = 0.25

gr(reuse=true, size = (1000, 1000))
#println(curr_points)
lim = (-0.05, 1.05)
scatter(curr_points[1, :], curr_points[2, :], xlim=lim, ylim=lim, legend=false, markerstrokewidth=0, markersize=7*masses, color=:black, label="")

for i ∈ 1:TIMESTEPS
  time = @elapsed barnesHutUpdate!(acc, particles, θ_sq)
  println(time)
  for i in 1:length(particles)
    particle = particles[i]
    tmp = copy(particle.pos)
    # verlet integration with boundary thresholding
    particle.pos .= max.(min.(2*particle.pos - particle.prev_pos + G * acc[:, i] * Δt^2, 1.0), 0.0)
    particle.prev_pos .= tmp
  end
  # put back into arrays for plotting, expensive but ok
  for i in 1:length(particles)
    curr_points[:, i] .= particles[i].pos
    masses[i] = particles[i].mass
  end

  #println("curr_points")
  #println(curr_points)
  #println("acc")
  #println(acc)

  # TEST
  """
  # hacky way, conversions here and back, keeping this way as have verified this test / good to have same test for everything
  complex_points = reinterpret(ComplexF64, curr_points)
  kernel_mtx::Array{ComplexF64, 2} = 1 ./ (transpose(complex_points) .- complex_points)
  foreach(i -> kernel_mtx[i, i] = zero(kernel_mtx[1, 1]), 1:length(complex_points))
  #show(stdout, "text/plain", kernel_mtx)
  # can't do simple dot product unfortunately
  # not really force, missing a mass term so actually acceleration
  correct_ω_p = vec(sum(kernel_mtx.*masses, dims=1))
  correct_acc = [real(correct_ω_p), -imag(correct_ω_p)]
  println(correct_acc)
  actual_acc = [acc[1, :], acc[2, :]]
  println(actual_acc)
  # Don't believe this will ever be close enough to not throw an assertion
  @assert correct_acc ≈ actual_acc 
  """
  
  scatter(curr_points[1, :], curr_points[2, :], xlim=lim, ylim=lim, legend=false, markerstrokewidth=0, markersize=7*masses, color=:black, label="")
  gui()
  #sleep(0.01)
end
