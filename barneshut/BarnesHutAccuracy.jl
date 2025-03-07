
""" Plot accuracy of Barnes Hut 
    Inefficient Plotting pipeline but thats ok, just testing """


using LinearAlgebra


const G = 1e-2
# softening paramter to avoid a -> inf when particles close
const S = 1e-32
# timestep
const Δt = 1e-2
# num timesteps
const TIMESTEPS = 1
# number of bodies
const N = 5000


include("BarnesHut.jl")


#θs = [x for x in 0.0:0.:5.0]
θs = exp10.(-3:0.5:0.5)
θs_sq = [x^2 for x in θs]
#Ns = exp10.(2:0.5:3)


abs_error = Array{Float64, 1}(undef, length(θs))
rel_error = Array{Float64, 1}(undef, length(θs))


for i in 1:length(θs)
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

  θ_sq = θs_sq[i]

  barnesHutUpdate!(acc, particles, θ_sq)

  # TEST
  # hacky way, conversions here and back, keeping this way as have verified this test / good to have same test for everything
  complex_points = reinterpret(ComplexF64, curr_points)
  kernel_mtx::Array{ComplexF64, 2} = 1 ./ (transpose(complex_points) .- complex_points)
  foreach(i -> kernel_mtx[i, i] = zero(kernel_mtx[1, 1]), 1:length(complex_points))
  #show(stdout, "text/plain", kernel_mtx)
  # can't do simple dot product unfortunately
  # not really force, missing a mass term so actually acceleration
  correct_ω_p = vec(sum(kernel_mtx.*masses, dims=1))
  correct_acc = [real(correct_ω_p), -imag(correct_ω_p)]
  actual_acc = [acc[1, :], acc[2, :]]
  #println(correct_acc)
  #println(actual_acc)
  abs_error[i] = norm(correct_acc .- actual_acc)
  rel_error[i] = norm(correct_acc .- actual_acc) / norm(correct_acc)
end

#print("absolute error: ")
#println(abs_error)
#print("relative error: ")
#println(rel_error)
gr(size = (1000, 1000))
Plots.resetfontsizes();
Plots.scalefontsizes(2.0);
p1 = plot(θs, abs_error, label="absolute error")
p2 = plot!(θs, rel_error, label="relative error") 
title!(".\nError (L2 norm) vs. Barnes Hut MAC (θ)\nNumber of bodies: $N")
p3 = plot!([eps(Float32)], seriestype="hline", label="eps(Float32)")
p4 = plot!([eps(Float64)], seriestype="hline", label="eps(Float64)")
yaxis!(:log)
#ylims!(0.001, 2.0)
xlabel!("Barnes Hut Multipole Acceptance Criterion (θ)")
ylabel!("Error (L2 norm)")
gui() 
