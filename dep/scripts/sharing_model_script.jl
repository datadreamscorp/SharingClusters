### A Pluto.jl notebook ###
# v0.19.4

using Markdown
using InteractiveUtils

# ╔═╡ 289610f8-d533-4a19-a840-dcff5c2b8c18
using DrWatson

# ╔═╡ 8e07a598-1e5c-11ec-3ead-a79dca36b027
begin
	@quickactivate "Sharing Groups Model"
	using Agents, LightGraphs, Distributions, Pipe, Random
end

# ╔═╡ fcf3f7cd-2828-482b-a9a3-1ff9933b2aa5
using StatsBase

# ╔═╡ 205362df-1a2f-4ca9-acf0-e2bad4aee761
using Gadfly

# ╔═╡ 405db3b4-475d-4ebf-961b-906d55f15bca
"""
worker agent class
"""
mutable struct Worker <: AbstractAgent
	id::Int64 #agent id
	pos::Int64 #position in network
	niche::Int64 #chosen niche
	group::Int64 #agent group id
	age::Int64 #agent age in terms of ticks
	strategy::Bool #agent strategy - 1 if sharer, 0 if loner
	p::Float64 #probability of forming new connections
	raw_payoff::Float64 #agent raw payoff
	payoff::Float64 #agent final payoff
	fitness::Float64 #agent fitness
	occupancy::Int64 #agent niche occupancy
	success::Bool #success at current work attempt
end

# ╔═╡ b672bd0b-0552-4d98-aecf-e140ac88d5d6
"""
this function randomly accomodates unsuccessful agents into successful peers' niches, as long as there is enough niche capacity

payoffs are calculated for newly accomodated agents and niches are filled up as agents come in
"""
function accomodate_losers!(model, sharers)
	successful = filter(a -> a.success, sharers)
	failed = setdiff(sharers, successful)
	
	for worker in shuffle(model.rng, failed)
		neighbors = filter( 
					a -> a ≠ worker.id, 
					neighborhood(model.network, worker.id, 1)
				)
		if (neighbors|>length) > 0
			successful_neighbors = filter(a -> model[a].success, neighbors)
			
			if (successful_neighbors|>length) > 0
				for neighbor in shuffle(model.rng, successful_neighbors)
					if model[neighbor].occupancy < model.M
						model[neighbor].occupancy += 1
						payoff = rand(model.rng) < model.λ ? model.b/model.seasons : 0
						worker.raw_payoff += payoff
						break
					end
				end
			end
			
		end
	end
	
end

# ╔═╡ bda7b12b-8024-4cc7-bc38-a5fb41f07541
"""
calculates fitness of all model agents

stores fitness in each agent, overwriting previous fitness value 
"""
function fitness!(model, sharers)
	for a in setdiff(allagents(model)|>collect, sharers)
		a.fitness += ( (1 + model.δ)^model.u ) * (a.payoff) / model.seasons
		#a.fitness = exp(a.payoff)
		#a.fitness = model.δ*log(a.payoff)
	end
	for a in sharers
		if model.s > 0
			a.fitness += ( (1 + model.δ) ) * (a.payoff) / model.seasons
		else
			a.fitness += ( (1 + model.δ)^(model.u + model.u*model.λ) ) * (a.payoff) / model.seasons
		end
	end
end

# ╔═╡ 4c66326b-d7de-4e2b-a91e-b089cb4bcf8e
"""
calculates payoffs of all model agents

stores payoffs in each agent, overwriting previous payoff value 
"""
function payoffs!(model)
	
	for worker in allagents(model) #reset payoffs
		worker.payoff = 0
	end
	
	for season in 1:(model.seasons)
		
		for worker in allagents(model)
			worker.raw_payoff = 0
			if rand(model.rng) < model.u
				worker.raw_payoff += model.B / model.seasons
				(worker.strategy == 0) && (worker.payoff += worker.raw_payoff)
				worker.success = true
			else
				worker.success = false
			end
		end

		sharers = filter( a -> a.strategy, allagents(model)|>collect )

		if model.M > 1
			accomodate_losers!(model, sharers)
		end
		
		for worker in sharers
			friends = filter( 
					a -> a ≠ worker.id, neighborhood(model.network, worker.id, 1)
				)
			if friends |> length ≠ 0
				
				for friend in friends 
					friend_degree = degree(model.network, friend)
					shared_benefit = model[friend].raw_payoff / friend_degree
					cost = model.C / model.seasons
					worker.payoff += (model.s * shared_benefit) - cost
				end
				
				worker.payoff += (1 - model.s) * worker.raw_payoff
				
			else
				worker.payoff += worker.raw_payoff
			end
		end
		
		fitness!(model, sharers)
		
	end
	
end

# ╔═╡ b6150c76-7aa5-4690-a720-3ffb1586fee5
"""
sharing groups model initialization function, with the form

`sharing_groups(; kwargs...)`

see function comments and/or specification and cheat sheet documents for details on keyword argument parameters

returns an Agents.jl `Model` object that can be stepped
"""
function sharing_groups(;
		N::Int64 = 201, #number of agents
		Sh::Int64 = 30, #initial size of sharing cluster 
		M::Int64 = 25, #niche capacity
		B::Float64 = 10.0, #benefits gained from working at niche
		b::Float64 = 10.0, #decreased benefits
		C::Float64 = 0.01, #cost of connections
		u::Float64 = 0.1, #security, probability of niche success
		s::Float64 = 1.0, #sharing norm
		λ::Float64 = 0.1, #rate of benefit extraction from shared niche
		δ::Float64 = 1.0, #strength of selection
		p_init::Float64 = 1.0, #initial connection propensity of sharers
		seasons::Int64 = 1 #number of seasons in which working period is divided
	)
	
	model = ABM(
	        Worker,
			nothing;
	        properties = Dict(
				:N => N,
				:M => M,
				:B => B,
				:b => b,
				:C => C,
				:u => u,
				:s => s,
				:λ => λ,
				:δ => δ,
				:network => SimpleGraph(N),
				:seasons => seasons,
				:Sh => Sh,
				:p_init => p_init,
	        ),
			rng = RandomDevice()
			)
	
	#set up agents; N agents, N_l loners, 
	#a sharing cluster of size Sh = (N - N_l) and
	#mean degree k_bar ≃ p_init * (Sh - 1)
	#agents are set up with a strategy and their network
	#connections are saved in a Lightgraphs.jl graph
	#the graph is initialized first and then the
	#loner nodes are added
	
	for id in 1:N 
		add_agent!(		
			id, 
			model, 
			id, 
			id ≤ Sh ? 1 : 0,
			0,
			id ≤ Sh ? true : false,
			p_init,
			0,
			0,
			0,
			1,
			false,
		)	
	end
	
	#add the edges iteratively
	
	for id in 1:(Sh - 1)
		for node in (id+1):Sh
			if rand(model.rng) < p_init
				add_edge!( model.network, id, node )
			end
		end
	end
	
	payoffs!(model)
	
	return model
	
end

# ╔═╡ c464dff3-fe73-43a2-b271-1a45e69e773d
model = sharing_groups()

# ╔═╡ 5f142f87-e9d4-43e5-a0b8-4891b80060a3
sharers = filter( x -> x.strategy == 1, allagents(model)|>collect )

# ╔═╡ b56ee6fe-7dee-4194-8b30-471edd277a0c
non_sharers = filter( x -> x.strategy == 0, allagents(model)|>collect )

# ╔═╡ a669ac61-d780-4074-9746-0a7603b3b399
begin
	sharer_payoffs = []
	sharer_fitness = []
	for sh in sharers
		push!(sharer_payoffs, sh.payoff)
		push!(sharer_fitness, sh.fitness)
	end
end

# ╔═╡ 0a939cc0-9242-466d-80d9-87f92c5e0e22
plot(
	x = sharer_fitness,
	Geom.histogram
	)

# ╔═╡ 4887aea1-4ab2-481b-baee-04181c96fd2b
mean(sharer_payoffs)

# ╔═╡ a2bf7304-8b58-43ea-87ac-4abe88a04791
mean(sharer_fitness)

# ╔═╡ a84f2e50-f56d-4cd2-a60f-04f963114366
median(sharer_fitness)

# ╔═╡ d67f498b-f266-427c-8abd-1dceda5ebe31
maximum(sharer_fitness)

# ╔═╡ fa1d84f4-ef65-4038-bb80-adb307fdeb64
begin
	loner_payoffs = []
	loner_fitness = []
	for l in non_sharers
		push!(loner_payoffs, l.payoff)
		push!(loner_fitness, l.fitness)
	end
end

# ╔═╡ dae02780-c63f-4bb8-a11b-80fb074740a9
plot(
	x = loner_fitness,
	Geom.histogram()
	)

# ╔═╡ 61f258b1-a676-4780-a0b3-b335e17ba74c
mean(loner_payoffs)

# ╔═╡ edb244a1-0125-45ca-95f6-8c532502d3bc
mean(loner_fitness)

# ╔═╡ dfc70a5a-1582-447d-9f31-22dfce6d7115
median(loner_fitness)

# ╔═╡ 0fcc24e7-4c56-4e87-ad38-84fb08c64dc5
maximum(loner_fitness)

# ╔═╡ 1706eb61-69d9-4884-81b0-388a3672345b
begin
	share_mean_payoff = []
	share_median_payoff = []
	lone_mean_payoff = []
	
	share_mean_fitness = []
	share_median_fitness = []
	lone_mean_fitness = []
	
	for n in 1:2000
		model = sharing_groups()
		sharers = filter( x -> x.strategy == 1, allagents(model)|>collect )
		non_sharers = filter( x -> x.strategy == 0, allagents(model)|>collect )
	
		sharer_payoffs = []
		sharer_fitness = []
		for sh in sharers
			push!(sharer_payoffs, sh.payoff)
			push!(sharer_fitness, sh.fitness)
		end
		
		loner_payoffs = []
		loner_fitness = []
		for l in non_sharers
			push!(loner_payoffs, l.payoff)
			push!(loner_fitness, l.fitness)
		end
		
		push!(share_median_payoff, median(sharer_payoffs))
		push!(share_mean_payoff, mean(sharer_payoffs))
		push!(lone_mean_payoff, mean(loner_payoffs))
		
		push!(share_median_fitness, median(sharer_fitness))
		push!(share_mean_fitness, mean(sharer_fitness))
		push!(lone_mean_fitness, mean(loner_fitness))
	end
end

# ╔═╡ 1087b77a-8331-4c28-a1b3-454e302aed32
plot(
	x=share_median_fitness,
	Geom.histogram()
	)

# ╔═╡ 725a053a-7b5c-4d66-aaa3-1a3db8dd8f0f
mean(share_median_fitness)

# ╔═╡ a14f5b18-0027-4d3e-94c5-7e823059b899
plot(
	x=share_mean_fitness,
	Geom.histogram()
	)

# ╔═╡ 927c307e-ca56-4c3c-985d-42449c78c914
mean(share_mean_fitness)

# ╔═╡ 46dab85a-dce5-4b9e-8c15-c7c28af6f146
mean(share_mean_payoff)

# ╔═╡ 5d019f9e-fbad-48d4-aad9-91332b82c7d7
plot(
	x=lone_mean_fitness,
	Geom.histogram()
	)

# ╔═╡ 8e6119cc-eacb-4eef-b347-9979a2352952
mean(lone_mean_fitness)

# ╔═╡ 17284ac5-f02d-4fb2-bd9f-1f2da738092b
mean_deg = model.Sh * model.p_init

# ╔═╡ a292bc04-32c2-4f1d-ac37-e99aab47ad26
mean_deg_calc = mean( degree( model.network, [a.id for a in sharers] ) )

# ╔═╡ c4856a51-7f6b-4079-a117-40491b777b96
γ = (1 - (1 - model.u)^(mean_deg))

# ╔═╡ a343ea1d-72d9-43f7-a837-2d2f48c25b50
β = ( 1 - (1 - (1 / ( (1-model.u)*mean_deg) ) )^(model.u*mean_deg*(model.M - 1)) )

# ╔═╡ 8eae190b-52ac-4eea-b532-1d5631bee8eb
σ = γ * β

# ╔═╡ 33471e34-6e55-47be-9760-6654fd294059
predicted_mean_benefit = model.u*model.B + (1 - model.u)*σ*model.λ*model.b - mean_deg*model.C

# ╔═╡ 843acdbc-a3a1-49b9-9e58-6cf0b0376aac
predicted_mean_fitness = (1 + model.δ)*predicted_mean_benefit

# ╔═╡ 1b9e4d4e-7ca2-486a-9913-47963cbcc8b1
(model.u + (1 - model.u)*σ*model.λ)*model.B/mean_deg - model.C

# ╔═╡ 29723814-3826-4fca-89c5-dcce53ca95bd
model.u*model.B + (1 - model.u)*γ*model.λ*model.b - mean_deg*model.C

# ╔═╡ c714677e-89ed-49f3-a1dd-e42f60b573d4
model.u*model.B + (1 - model.u)*β*model.λ*model.b - mean_deg*model.C

# ╔═╡ Cell order:
# ╠═289610f8-d533-4a19-a840-dcff5c2b8c18
# ╠═8e07a598-1e5c-11ec-3ead-a79dca36b027
# ╠═405db3b4-475d-4ebf-961b-906d55f15bca
# ╠═b6150c76-7aa5-4690-a720-3ffb1586fee5
# ╠═4c66326b-d7de-4e2b-a91e-b089cb4bcf8e
# ╠═b672bd0b-0552-4d98-aecf-e140ac88d5d6
# ╠═bda7b12b-8024-4cc7-bc38-a5fb41f07541
# ╠═c464dff3-fe73-43a2-b271-1a45e69e773d
# ╠═5f142f87-e9d4-43e5-a0b8-4891b80060a3
# ╠═b56ee6fe-7dee-4194-8b30-471edd277a0c
# ╠═0a939cc0-9242-466d-80d9-87f92c5e0e22
# ╠═a669ac61-d780-4074-9746-0a7603b3b399
# ╠═fcf3f7cd-2828-482b-a9a3-1ff9933b2aa5
# ╠═4887aea1-4ab2-481b-baee-04181c96fd2b
# ╠═a2bf7304-8b58-43ea-87ac-4abe88a04791
# ╠═a84f2e50-f56d-4cd2-a60f-04f963114366
# ╠═d67f498b-f266-427c-8abd-1dceda5ebe31
# ╠═dae02780-c63f-4bb8-a11b-80fb074740a9
# ╠═fa1d84f4-ef65-4038-bb80-adb307fdeb64
# ╠═61f258b1-a676-4780-a0b3-b335e17ba74c
# ╠═edb244a1-0125-45ca-95f6-8c532502d3bc
# ╠═dfc70a5a-1582-447d-9f31-22dfce6d7115
# ╠═0fcc24e7-4c56-4e87-ad38-84fb08c64dc5
# ╠═205362df-1a2f-4ca9-acf0-e2bad4aee761
# ╠═1706eb61-69d9-4884-81b0-388a3672345b
# ╠═1087b77a-8331-4c28-a1b3-454e302aed32
# ╠═725a053a-7b5c-4d66-aaa3-1a3db8dd8f0f
# ╠═a14f5b18-0027-4d3e-94c5-7e823059b899
# ╠═927c307e-ca56-4c3c-985d-42449c78c914
# ╠═46dab85a-dce5-4b9e-8c15-c7c28af6f146
# ╠═5d019f9e-fbad-48d4-aad9-91332b82c7d7
# ╠═8e6119cc-eacb-4eef-b347-9979a2352952
# ╠═17284ac5-f02d-4fb2-bd9f-1f2da738092b
# ╠═a292bc04-32c2-4f1d-ac37-e99aab47ad26
# ╠═c4856a51-7f6b-4079-a117-40491b777b96
# ╠═a343ea1d-72d9-43f7-a837-2d2f48c25b50
# ╠═8eae190b-52ac-4eea-b532-1d5631bee8eb
# ╠═33471e34-6e55-47be-9760-6654fd294059
# ╠═843acdbc-a3a1-49b9-9e58-6cf0b0376aac
# ╠═1b9e4d4e-7ca2-486a-9913-47963cbcc8b1
# ╠═29723814-3826-4fca-89c5-dcce53ca95bd
# ╠═c714677e-89ed-49f3-a1dd-e42f60b573d4
