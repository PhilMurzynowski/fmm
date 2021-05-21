
""" Plot accuracy of Fast Multipole Method """

# Have some warnings due to ambiguous scope to clean up


using LinearAlgebra
using Plots


include("Quadtree.jl")
include("Multipole.jl")
include("BinomialTable.jl")


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

# Changing number of terms to keep in multipole expansion.
# System was designed around the assumption of constant P (number of multipole terms),
# so will have to do some inefficient rebuilds, but this inefficiency
# is currently acceptable for testing purposes.
Ps = 2:50
abs_error = Array{Float64, 1}(undef, length(Ps))
rel_error = Array{Float64, 1}(undef, length(Ps))

preallocated_size = floor(Int, N/2)
preallocated_mtx = Array{ComplexF64, 2}(undef, preallocated_size, preallocated_size)

for p in Ps
  quadtree = buildQuadtree(TREE_DEPTH, p)
  pos_memory = Array{ComplexF64}(undef, N, NUM_PAST_POSITIONS)
  # Not seeding for accuracy testing, random
  # Random.seed!(2)
  pos_memory[:, 1]  = rand(ComplexF64, N)
  tangent_velocity = imag(pos_memory[:, 1] .- (0.5 + 0.5im)) .- 1im*real(pos_memory[:, 1] .- (0.5 + 0.5im))
  pos_memory[:, 2]  = pos_memory[:, 1] .+ 1e-2*tangent_velocity
  masses         = 0.9*rand(Float64, N) .+ 1
  ω_p            = Array{ComplexF64, 1}(undef, N)
  binomial_table = largeBinomialTable(p)
  binomial_table_t = binomialTableTransposed(p)
  large_binomial_table_t = largeBinomialTableTransposed(p)
  next_idx = 3
  curr_idx = 2
  prev_idx = 1
  curr_points = @view pos_memory[:, curr_idx]
  prev_points = @view pos_memory[:, prev_idx]
  updateQuadtreePointMasses(quadtree, curr_points, masses, prev_points)
  FMM!(quadtree, curr_points, masses, ω_p, binomial_table, binomial_table_t, large_binomial_table_t, preallocated_mtx)

  # TEST
  # As a matter of fact optimize out the last mass multiplication so just calculating acceleration without G instead of
  # calculating force and then diving out mass again.
  forces = [real(ω_p), -imag(ω_p)]
  correct_ω_p = similar(ω_p)
  correct_ω_p .= zero(correct_ω_p[1])
  kernel_mtx::Array{ComplexF64, 2} = 1 ./ (transpose(curr_points) .- curr_points)
  foreach(i -> kernel_mtx[i, i] = zero(kernel_mtx[1, 1]), 1:length(curr_points))
  correct_ω_p = vec(sum(kernel_mtx.*masses, dims=1))
  correct_forces = [real(correct_ω_p), -imag(correct_ω_p)]
  abs_error[p - Ps[1] + 1] = norm(correct_forces .- forces)
  rel_error[p - Ps[1] + 1] = norm(correct_forces .- forces) / norm(correct_forces)
end

#print("absolute error: ")
#println(abs_error)
#print("relative error: ")
#println(rel_error)
gr(size = (1000, 1000))
Plots.resetfontsizes();
Plots.scalefontsizes(1.5);
p1 = plot(Ps, abs_error, label="absolute error")
p2 = plot!(Ps, rel_error, label="relative error") 
title!(".\nError (L2 norm) vs. Number of expansion Terms (P)\nNumber of bodies: $N, Tree depth: $TREE_DEPTH")
p3 = plot!([eps(Float32)], seriestype="hline", label="eps(Float32)")
p4 = plot!([eps(Float64)], seriestype="hline", label="eps(Float64)")
yaxis!(:log)
xlabel!("Number of expansion terms (P)")
ylabel!("Error (L2 norm)")
gui() 


