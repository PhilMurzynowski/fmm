
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
# const N = 100
# number of past positions saved
const NUM_PAST_POSITIONS = 3 # do not support plotting history yet, set to 3


function sumAllFlops(flopCounter::GFlops.Counter)
	FieldsInStruct=fieldnames(typeof(flopCounter));
  total_flops = 0
	for i=1:length(FieldsInStruct)
		total_flops += getfield(flopCounter, FieldsInStruct[i])
	end
  return total_flops
end


# Changing number of terms to keep in multipole expansion.
# System was designed around the assumption of constant P (number of multipole terms),
# so will have to do some inefficient rebuilds, but this inefficiency
# is currently acceptable for testing purposes.
start_tree_depth = 3
end_tree_depth = 4
tree_depths = start_tree_depth:end_tree_depth
Ps = [2, 10, 33]
#Ns = exp10.(1:0.25:2)
Ns = [10, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000]

flopcounts = Array{Float64, 3}(undef, length(Ns), length(Ps), length(tree_depths))

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

for i in 1:length(tree_depths)
  tree_depth = tree_depths[i]
  for j in 1:length(Ps)
    p = Ps[j]
    binomial_table   = largeBinomialTable(p)
    binomial_table_t = binomialTableTransposed(p)
    large_binomial_table_t = largeBinomialTableTransposed(p)
    for k in 1:length(Ns)
      N = trunc(Int, Ns[k])
      # not included in flop counts as this only need to be built once per simulation
      # COULD this be hoisted one loop?
      quadtree = buildQuadtree(tree_depth, p)
      # Not seeding for accuracy testing, random
      # Random.seed!(2)
      prev_points      = rand(ComplexF64, N)
      tangent_velocity = imag(prev_points .- (0.5 + 0.5im)) .- 1im*real(prev_points .- (0.5 + 0.5im))
      curr_points      = prev_points .+ 1e-2*tangent_velocity
      masses           = 0.9*rand(Float64, N) .+ 1
      ω_p              = Array{ComplexF64, 1}(undef, N)
      #println(tree_depth)
      #println(p)
      #println(N)
      # FLOPCOUNT THIS UPDATE AS WELL
      flopCounter1 = @count_ops updateQuadtreePointMasses(quadtree, curr_points, masses, prev_points)
      flops = sumAllFlops(flopCounter1)
      # Benchmark tools takes too long atm
      flopCounter2 = @count_ops FMM!(quadtree, curr_points, masses, ω_p, binomial_table, binomial_table_t, large_binomial_table_t, preallocated_mtx)
      flops += sumAllFlops(flopCounter2)
      println(flops)
      flopcounts[k, j, i] = flops
    end
  end
end

#yaxis!(:log)
gr(size = (1000, 1000))
Plots.resetfontsizes();
Plots.scalefontsizes(2.0);
p = plot()
for i in 1:length(tree_depths)
  tree_depth = tree_depths[i]
  for j in 1:length(Ps)
    p = Ps[j]
    scatter!(Ns, flopcounts[:, j, i], label="p: $p, tree depth: $tree_depth", markerstrokewidth=0, markersize=10, legend=:topleft)
  end
end
title!(".\nFMM FlopCount vs. Number of bodies (N)\n")
ylabel!("FMM Flopcount")
xlabel!("Number of bodies (N)")
# Sometimes may want to set this
#ylims!((0, 0.1))
gui()

