
""" Plot accuracy of Fast Multipole Method """

# BROKEN!!!! Currently not functional
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
const N = 30
# number of past positions saved
const NUM_PAST_POSITIONS = 3 # do not support plotting history yet, set to 3
# depth of tree to construct
const TREE_DEPTH = 3

# Changing number of terms to keep in multipole expansion.
# System was designed around the assumption of constant P (number of multipole terms),
# so will have to do some inefficient rebuilds, but this inefficiency
# is currently acceptable for testing purposes.
Ps = 2:60
abs_error = Array{Float64, 1}(undef, length(Ps))
rel_error = Array{Float64, 1}(undef, length(Ps))

# setup once so don't have to recompute the expensive O(N^2) calculation
pointsminusΔt  = rand(ComplexF64, N)
tan_vel        = imag(pointsminusΔt .- (0.5 + 0.5im)) .- 1im*real(pointsminusΔt .- (0.5 + 0.5im))
points0Δt      = pointsminusΔt .+ 1e-2*tan_vel
mass           = 0.9*rand(Float64, N) .+ 1

# PRE COLOR SORT ONCE so aligned for testing
quadtree = buildQuadtree(TREE_DEPTH, 4)
println(pointsminusΔt)
updateQuadtreePointMasses(quadtree, points0Δt, mass, pointsminusΔt)
println(pointsminusΔt)

# correct forces
correct_ω_p = zeros(ComplexF64, N)
#kernel_mtx::Array{ComplexF64, 2} = 1 ./ (transpose(curr_points) .- curr_points)
kernel_mtx = 1 ./ (transpose(points0Δt) .- points0Δt)
foreach(i -> kernel_mtx[i, i] = zero(kernel_mtx[1, 1]), 1:length(points0Δt))
correct_ω_p = vec(sum(kernel_mtx.*mass, dims=1))
correct_forces = [real(correct_ω_p), -imag(correct_ω_p)]

for p in Ps
  quadtree = buildQuadtree(TREE_DEPTH, p)
  pos_memory = Array{ComplexF64}(undef, N, NUM_PAST_POSITIONS)
  # Not seeding for accuracy testing, random
  # Random.seed!(2)
  pos_memory[:, 1]  = pointsminusΔt
  pos_memory[:, 2]  = points0Δt
  masses         = mass
  ω_p            = Array{ComplexF64, 1}(undef, N)
  binomial_table = largeBinomialTable(p)
  binomial_table_t = binomialTableTransposed(p)
  preallocated_size = floor(Int, N/2)
  preallocated_mtx::Array{ComplexF64, 2} = Array{ComplexF64, 2}(undef, preallocated_size, preallocated_size)
  next_idx = 3
  curr_idx = 2
  prev_idx = 1
  curr_points = @view pos_memory[:, curr_idx]
  prev_points = @view pos_memory[:, prev_idx]
  updateQuadtreePointMasses(quadtree, curr_points, masses, prev_points)
  FMM!(quadtree, curr_points, masses, ω_p, binomial_table, binomial_table_t, preallocated_mtx)

  # TEST
  forces = [real(ω_p), -imag(ω_p)]
  abs_error[p - Ps[1] + 1] = norm(correct_forces .- forces)
  println(abs_error[p - Ps[1] + 1])
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
title!(".\nNumber of expansion Terms (P) vs. Error (L2 norm)\nNumber of bodies: $N, Tree depth: $TREE_DEPTH")
p3 = plot!([eps(Float32)], seriestype="hline", label="eps(Float32)")
p4 = plot!([eps(Float64)], seriestype="hline", label="eps(Float64)")
yaxis!(:log)
xlabel!("Number of expansion terms (P)")
ylabel!("Error (L2 norm)")
gui() 


