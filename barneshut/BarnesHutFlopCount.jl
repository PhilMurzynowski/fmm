
""" Plot accuracy of Barnes Hut 
    Inefficient Plotting pipeline but thats ok, just testing """


using GFlops


const G = 1e-2
# softening paramter to avoid a -> inf when particles close
const S = 1e-32
# timestep
const Δt = 1e-2
# num timesteps
const TIMESTEPS = 1
# number of bodies
# const N = 1000


include("BarnesHut.jl")


function sumAllFlops(flopCounter::GFlops.Counter)
	FieldsInStruct=fieldnames(typeof(flopCounter));
  total_flops = 0
	for i=1:length(FieldsInStruct)
		total_flops += getfield(flopCounter, FieldsInStruct[i])
	end
  return total_flops
end


# these thetas are arrange backwards for convenience
#θs = exp10.(-3:0.5:0.5)
#θs = exp10.(0.5:-2.5:-3.0)
θs = [0.5, 1.0]
θs_sq = [x^2 for x in θs]
start_theta = θs[1]

Ns = exp10.(1:0.25:2)

flopcounts = Array{Float64, 2}(undef, length(Ns), length(θs))

# warning, accidentally swapped i and j
for j in 1:length(θs)
  for i in 1:length(Ns)
    N = trunc(Int, Ns[i])
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

    θ_sq = θs_sq[j]

    barnesHutUpdate!(acc, particles, θ_sq)
    flopCounter = @count_ops barnesHutUpdate!($acc, $particles, $θ_sq)
    flops = sumAllFlops(flopCounter)
    println(flops)
    flopcounts[i, j] = flops
  end
end

title!(".\nBarnes Hut FlopCount vs. Number of bodies (N)\n")
#yaxis!(:log)
gr(size = (1000, 1000))
Plots.resetfontsizes();
Plots.scalefontsizes(1.75);
p = scatter(Ns, flopcounts[:, 1], label="θ: $start_theta", markerstrokewidth=0, markersize=10)
for i in 2:length(θs)
  θ = θs[i]
  scatter!(Ns, flopcounts[:, i], label="θ: $θ", markerstrokewidth=0, markersize=10)
end
ylabel!("Barnes Hut Flopcount")
xlabel!("Number of bodies (N)")
# Sometimes may want to set this
#ylims!((0, 0.1))
gui()
