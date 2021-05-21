
""" Plot flopcount of Fast Multipole Method In Comparison to Barnes-Hut"""

# Have some warnings due to ambiguous scope to clean up


using LinearAlgebra
using Plots
using GFlops


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
const N = 5000
# number of past positions saved
const NUM_PAST_POSITIONS = 3 # do not support plotting history yet, set to 3

# Changing number of terms to keep in multipole expansion.
# System was designed around the assumption of constant P (number of multipole terms),
# so will have to do some inefficient rebuilds, but this inefficiency
# is currently acceptable for testing purposes.
start_tree_depth = 3
end_tree_depth = 5
tree_depths = start_tree_depth:end_tree_depth
Ps = 2:10
flopcounts = Array{Float64, 2}(undef, length(Ps), length(tree_depths))

preallocated_size = floor(Int, N/2)
preallocated_mtx = Array{ComplexF64, 2}(undef, preallocated_size, preallocated_size)

# RUN ONCE FOR FAIR BENCHMARKING
p = 4
quadtree = buildQuadtree(start_tree_depth, p)
# Not seeding for accuracy testing, random
# Random.seed!(2)
prev_points      = rand(ComplexF64, N)
tangent_velocity = imag(prev_points .- (0.5 + 0.5im)) .- 1im*real(prev_points .- (0.5 + 0.5im))
curr_points      = prev_points .+ 1e-2*tangent_velocity
masses           = 0.9*rand(Float64, N) .+ 1
ω_p              = Array{ComplexF64, 1}(undef, N)
binomial_table   = binomialTable(p)
binomial_table_t = binomialTableTransposed(p)
large_binomial_table_t = largeBinomialTableTransposed(p)
updateQuadtreePointMasses(quadtree, curr_points, masses, prev_points)
FMM!(quadtree, curr_points, masses, ω_p, binomial_table, binomial_table_t, large_binomial_table_t, preallocated_mtx)

for tree_depth in tree_depths
  for p in Ps
    quadtree = buildQuadtree(tree_depth, p)
    # Not seeding for accuracy testing, random
    # Random.seed!(2)
    prev_points      = rand(ComplexF64, N)
    tangent_velocity = imag(prev_points .- (0.5 + 0.5im)) .- 1im*real(prev_points .- (0.5 + 0.5im))
    curr_points      = prev_points .+ 1e-2*tangent_velocity
    masses           = 0.9*rand(Float64, N) .+ 1
    ω_p              = Array{ComplexF64, 1}(undef, N)
    binomial_table   = largeBinomialTable(p)
    binomial_table_t = binomialTableTransposed(p)
    large_binomial_table_t = largeBinomialTableTransposed(p)
    # FLOPCOUNT THIS UPDATE AS WELL
    updateQuadtreePointMasses(quadtree, curr_points, masses, prev_points)
    # Benchmark tools takes too long atm
    flopcount = @count_ops FMM!(quadtree, curr_points, masses, ω_p, binomial_table, binomial_table_t, large_binomial_table_t, preallocated_mtx)
    println(flopcount)
    #flopcounts[p-Ps[1]+1, tree_depth-start_tree_depth+1] = flopcount
  end
end

gr(size = (1000, 1000))
Plots.resetfontsizes();
Plots.scalefontsizes(1.5);
p = scatter(Ps, flopcounts[:, 1], label="tree depth: $start_tree_depth")
for i in start_tree_depth+1:end_tree_depth
  depth = tree_depths[i-start_tree_depth+1]
  scatter!(Ps, flopcounts[:, i-start_tree_depth+1], label="tree depth: $depth")
end
title!(".\nFMM Flopcount vs. Number of expansion Terms (P)\nNumber of bodies: $N")
#yaxis!(:log)
ylabel!("FMM Flopcount")
xlabel!("Number of expansion terms (P)")
# Sometimes may want to set this
#ylims!((0, 0.1))
gui() 


