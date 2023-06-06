
# ╔═╡ 9a26ccb8-ecb6-11ec-387d-2d2f24a4f5eb
using Agents, Random, Distributions, Statistics, StatsBase

# ╔═╡ 02d3a34f-7aa1-4a73-b675-606969f1f282
"""
Cluster agent class
"""
mutable struct Cluster <: AbstractAgent
	id::Int64
	size::Int64
	s::Float64 #sharing norm of cluster
	total_fitness::Float64
	average_fitness::Float64
	ϕ::Float64 #probability of receiving new member
	τ::Float64 #probability of losing a member (if whole pop is in the clusters)
	#fission_prob::Float64
end

# ╔═╡ f99fd667-3721-4580-9ac1-dd758de66493
function sharing_groups_model(;
	N::Int64 = 20, #number of sharing clusters
	max_N::Int64 = 200, #max number of sharing clusters
	n::Int64 = 200, #population size
	r::Int64 = 1, #number of agents that join/abandon clusters
	bc_ratio::Float64 = 200.0,
	T::Int64 = 100, #number of sub-periods in productive period
	init_size::Int64 = 10, #initial max size of clusters (binomial distribution)
	init_size_prob::Float64 = 0.5, #probability of cluster size
	B::Float64 = 20.0,
	w0::Float64 = 0.5,
	u::Float64 = 0.5, #average risk
	δ::Float64 = exp(1) - 1, #strength of selection
	#γ::Float64 = 10.0, #steepness of sigmoid for fission probability
	σ_small::Float64 = 0.05, #variance of inherited sharing norm after fission
	σ_large::Float64 = 0.05, #prob of large mutation in sharing norm after fission
	rep::Int64 = 1, #replicate number (paramscan only)
	size_bias::Bool = true, #whether there is a size bias in recruitment,
	seed::Int64 = 80085,
	true_random::Bool = false,
	total_ticks::Int64 = 5000,
	)

	model = ABM(
		Cluster,
		nothing,
		properties = Dict(
			:N => N,
			:max_N => max_N,
			:r => r,
			:T => T,
			:B => B,
			:C => B/bc_ratio,
			:w0 => w0,
			:size_dist => Binomial(init_size, init_size_prob),
			:share_dist => Beta(1, 1),
			:u => u,
			:δ => δ,
			#:γ => γ,
			:σ_small => σ_small,
			:σ_large => σ_large,
			:n => n,
			:current_n => 0,
			:current_N => N,
			:max_cluster_size => 0,
			:mean_sharing => 0.0,
			:median_sharing => 0.0,
			:mean_cluster_size => 0.0,
			:median_cluster_size => 0.0,
			:mean_sharing_vector => Float64[], #containers for pop-level data
			:median_sharing_vector => Float64[],
			:num_clusters_vector => Int64[],
			:mean_cluster_size_vector => Float64[],
			:median_cluster_size_vector => Float64[],
			:loner_fitness => (1+δ)^( T*( u*log(B) + (1-u)*log(1.0) ) ),
			:pop_fitness => 0.0,
			:tick => 0,
			:rep => rep,
			:size_bias => size_bias,
			:total_ticks => total_ticks,
		),
		rng = true_random ? RandomDevice() : MersenneTwister(seed),
	)

	for a in 1:N #add initial agents
		add_agent!(
			model,
			rand(model.rng, model.size_dist),
			rand(model.rng, model.share_dist),
			0.0,
			0.0,
			0.0,
			0.0,
		)
	end

	model.current_n = sum([a.size for a in allagents(model)])

	return model

end

# ╔═╡ b8e082ae-c735-4449-87c9-28ec88191f76
function generate_fitness!(model)

	model.pop_fitness = 0.0
	u = model.u
	B = model.B
	C = model.C
	w0 = model.w0

	for a in allagents(model) #iterate through all agents

		log_payoffs = zeros(a.size) #calculate payoffs on log scale for convenience

		for t in 1:model.T #iterate through every sub-period of productive period

			netcost = (a.size - 1)*C
			pool = 0
			members = Float64[]

			#calculate individual and pooled payoff for sub-period:
			for member in 1:(a.size) #iterate through every member of pool

				full_payoff = ( rand(model.rng) < u ) ?  B : 1.0
				pool += full_payoff > 1.0 ? a.s*B : 0
				member_payoff = full_payoff > 1.0 ? (1 - a.s) * B - netcost : 1.0 - netcost
				push!(members, member_payoff)

			end

			share = pool / a.size #give a share of resource to all members
			full_member_payoffs = clamp.(share .+ members, w0, Inf)

			#calculate the average log growth rate
			log_payoffs += log.( full_member_payoffs )

		end

		#calculate fitness
		a.total_fitness = sum( (1 + model.δ).^(log_payoffs) )

		a.average_fitness = a.total_fitness / a.size
		model.pop_fitness += a.total_fitness

	end

end

# ╔═╡ aed2940f-7fc1-4e1d-b916-82f4bbad555c
function growth!(model)

	#function for (more) stochastic version of fission NOT USED
	#fissfunc(fitness) = 1 / ( 1 + exp( -model.γ * (fitness - model.loner_fitness) ) )

	#calculate the total fitness of the population
	#pop_fitness = sum( [a.total_fitness for a in allagents(model)] )
	#pop_fission = sum( [fissfunc(a.average_fitness) for a in allagents(model)] )

	for a in allagents(model) #calculate the relative fitness of each cluster

		if model.size_bias #if we activate size bias, otherwise normal payoff bias
			a.ϕ = a.total_fitness / model.pop_fitness
		else
			a.ϕ = a.average_fitness / sum([a.average_fitness for a in allagents(model)])
		end

		a.τ = -(a.ϕ - 1.0) #also calculate the probability of losing a member
						   #this assumes a frequency-dependent advantage for clusters
		#a.fission_prob = fissfunc(a.average_fitness) / pop_fission
	end

	recruiter = sample( #sample a recruiter cluster that will grow
		model.rng,
		allagents(model) |> collect,
		Weights( [a.ϕ for a in allagents(model)] )
	)

	if model.current_n < model.n #if pop is not maxed out, add one and that's it

		recruiter.size += 1
		model.current_n += 1
		#rand(model.rng, Binomial(model.n, model.n_prob))

	else

		loser = sample( #if pop is maxed out, choose a loser
			model.rng,
			allagents(model) |> collect,
			Weights( [a.τ for a in allagents(model)] )
		)

		loser.size -= 1 #l0ser
		recruiter.size += 1 #winrar

		if loser.size < 1 #if the cluster is empty, end its misery

			kill_agent!(loser, model) #rip
			model.current_N -= 1

		end

	end

end


# ╔═╡ 4f39ca71-ae53-41f9-b067-40542b1467cd
function death_and_fission!(model)

	fission_candidates = [] #list of candidates for fission
	for a in allagents(model)
		if (a.average_fitness ≤ model.loner_fitness) & (a.size > 1)
			push!(fission_candidates, a)
		end
	end

	if length(fission_candidates) > 0

		fissioned = rand(model.rng, fission_candidates) #sample a cluster

		inh_sharing = clamp( rand(model.rng, Normal(fissioned.s, model.σ_small)), 0.0, 1.0 ) #inherited sharing norm

		mutated_sharing = rand(model.rng) #mutated sharing norm

		#get the size of the new offspring sharing cluster
		new_size = floor( rand(model.rng) * fissioned.size )

		if new_size > 0

			fissioned.size -= new_size #members leave from fissioned cluster

			add_agent!(
				model,
				new_size,
				rand(model.rng) < model.σ_large ? mutated_sharing : inh_sharing,
				0.0,
				0.0,
				0.0,
				0.0,
			)

			model.current_N += 1 #add a new cluster to our counter

			if model.current_N > model.max_N #if cluster pop is maxed out, time to die

				dead = sample(
					model.rng,
					allagents(model) |> collect,
					Weights( [-(a.ϕ - 1.0) for a in allagents(model)] )
				)

				model.current_n -= dead.size #remove the little ones

				model.current_N -= 1 #remove the big one

				kill_agent!(dead, model) #F

			end

		end

	end

end


# ╔═╡ 93a4cbd4-8154-4c25-9b69-82252acbdbcf
function sharing_step!(model)
	generate_fitness!(model)
	growth!(model)
	death_and_fission!(model)

	sharings = [
		a.size > 1 ? a.s : 0.0
		for a in allagents(model)
		]
	model.mean_sharing = mean(sharings)
	model.median_sharing = median(sharings)

	sizes = [a.size for a in allagents(model)]
	model.mean_cluster_size = mean(sizes)
	model.median_cluster_size = median(sizes)
	model.max_cluster_size = maximum(sizes)

	push!(model.mean_sharing_vector, model.mean_sharing)
	push!(model.median_sharing_vector, model.median_sharing)
	push!(model.mean_cluster_size_vector, model.mean_cluster_size)
	push!(model.median_cluster_size_vector, model.median_cluster_size)
	push!(model.num_clusters_vector, model.current_N)

	model.tick += 1

	if model.tick == model.total_ticks
		loners = filter(a -> a.size == 1, allagents(model)|>collect)
		for a in loners
			a.s = 0.0
		end
	end
end


function abm_plots(model; legend=true)

	function kmax(u; β=1.0, B=6.0, C=6.0/200, δ_base=0.2, s=0.7)
		δ_0 = δ_base + log(B - 1)

		r1 = ( u*(β - u^(δ_0-1))*s*B / C )
		r2 = ( log(1-u) * u * (1-u)^(r1) * s*β*B / C )

		if r2 > -1/exp(1)
			kbar = ( r1 - ( lambertw(r2) / log(1-u) ) ) / β
			kbar = kbar > 0 ? kbar : 0
		else
			kbar = 0
		end

		return kbar
	end

	optshare(b, k; u=0.5, c=0.1, β=1) = ( ( (β*k*c - 1) / ( 1 - (1-u)^(β*k) )*β ) + ( (1-u)*(b - β*k*c) / ( 1 - u*β*( 1 - (1-u)^(β*k) ) ) ) ) / b

	opt_s = optshare(model.B, model.median_cluster_size, u=model.u, c=model.C)
	s_time_plot = plot(
		model.mean_sharing_vector,
		xlabel="",
		ylab="sharing norm",
		label="mean",
		lw=1.5,
		legend=legend,
	)
	plot!(
		model.median_sharing_vector,
		xlabel="",
		label="median",
		lw=1.5
	)
	hline!(
		[clamp(opt_s, 0, 1)],
		label="s*",
		lw=1.5
	)

	t_calc = 2500:5000
	time_mean_size = sum(model.mean_cluster_size_vector[t_calc]) / length(t_calc)

	kbar_time_plot = plot(
		model.mean_cluster_size_vector,
		color="blue",
		xlabel="time",
		ylabel="cluster size",
		label="mean",
		lw=2,
		legend=legend
	)
	plot!(
		model.median_cluster_size_vector,
		color="dark red",
		xlabel="time",
		label="median",
		#legend=false
	)
	hline!(
		[1 .+ ( 0.5 .* kmax(model.u, B=model.B, C=model.C, s=model.median_sharing) )],
		lw=2,
		label="(k̄*/2) + 1"
	)
	#hline!([time_mean_size], lw=2.5, label = "time-averaged mean")

	abm_time_plots = plot(
		s_time_plot,
		kbar_time_plot,
		layout = (2,1)
	)

	return abm_time_plots

end
#end # module sharing_clusters_ABM
