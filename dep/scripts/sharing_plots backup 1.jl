### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# ╔═╡ 284c4c16-bcf1-11ec-36d9-f724bda4cf9b
using Plots, Roots, Distributions, Colors, Random, Statistics, DataFrames, CSV

# ╔═╡ dc062df1-f6e5-4046-977e-bdb09258a3c6
using LambertW, ProgressLogging

# ╔═╡ 2958dbaf-e162-4f49-8d83-5a2ef0d7c824
function zero_difs(;B=5, C=0.05, β=1, s=0.75, δ_0=1.1, δ_1=2.1, δ_2=4, δ_inf=100)

	kbarV_inf = Vector{Float64}()
	kbarV_delta_0 = Vector{Float64}()
	kbarV_delta_1 = Vector{Float64}()
	kbarV_delta_2 = Vector{Float64}()
	kbarW = Vector{Float64}()
	
	for u in (0.0:0.01:1)

		W_dif_iter(kbar; u=u, B=B, C=C, β=β, s=s) = u*log( (1-s)*B + s*u*β*B - β*kbar*C ) + (1-u)*log( 1 + s*u*β*B - β*kbar*C ) - u*log(B)

		zero_V_inf = s*u*( (β - u^(δ_inf-1))/β )*B/C
		push!(kbarV_inf, zero_V_inf)

		zero_V_delta_0 = s*u*( (β - u^(δ_0-1))/β )*B/C
		push!(kbarV_delta_0, zero_V_delta_0)
		
		zero_V_delta_1 = s*u*( (β - u^(δ_1-1))/β )*B/C
		push!(kbarV_delta_1, zero_V_delta_1)

		zero_V_delta_2 = s*u*( (β - u^(δ_2-1))/β )*B/C
		push!(kbarV_delta_2, zero_V_delta_2)
		
		zero_W = find_zeros(W_dif_iter, 0, zero_V_inf)
		push!(kbarW, zero_W[1])
		
	end

	return [kbarV_inf, kbarV_delta_0, kbarV_delta_1, kbarV_delta_2, kbarW]
	
end

# ╔═╡ 6fdeae6a-1e1e-4301-baa7-91e47b072bcd
function zero_difs_raw(;B=10, C=0.25, β=1, s=9/10, δ_0 = 0.1 + (B - 1)^0.38, δ_1=B^(1/3))

	kbarV_delta_0 = Vector{Float64}()
	kbarV_delta_1 = Vector{Float64}()
	kbarW = Vector{Float64}()
	koptV = Vector{Float64}()
	kbar_vec_0 = Vector{Float64}()
	kbar_vec_1 = Vector{Float64}()
	
	for u in (0.001:0.01:1)

		W_dif_iter(kbar; u=u, B=B, C=C, β=β, s=s) = u*log( (1-s)*B + s*u*β*B - β*kbar*C ) + (1-u)*log( 1 + s*u*β*B - β*kbar*C ) - u*log(B)

		zero_V_delta_0 = s*u*( (β - u^(δ_0-1))/β )*B/C
		push!(kbarV_delta_0, zero_V_delta_0)

		zero_V_delta_1 = s*u*( (β - u^(δ_1-1))/β )*B/C
		push!(kbarV_delta_1, zero_V_delta_1)
		
		zero_V_delta_100 = s*u*( (β - u^(100-1))/β )*B/C
		
		zero_W = find_zeros(W_dif_iter, 0, zero_V_delta_100)
		push!(kbarW, zero_W[1])

		kopt = log( -C/( u*s*β*B*log(1-u) ) ) / (β*log(1-u))
		push!(koptV, kopt > 0 ? kopt : 0)

		r1 = ( u*(β - u^(δ_0-1))*s*B / C )
		r2 = ( log(1-u) * u * (1-u)^(r1) * s*β*B / C )
		
		if r2 > -1/exp(1)
			kbar = ( r1 - ( lambertw(r2) / log(1-u) ) ) / β
			kbar = kbar > 0 ? kbar : 0
		else
			kbar = 0
		end
		
		push!(kbar_vec_0, kbar)

		
		r1_1 = ( u*(β - u^(δ_1-1))*s*B / C )
		r2_1 = ( log(1-u) * u * (1-u)^(r1_1) * s*β*B / C )
		
		if r2_1 > -1/exp(1)
			kbar = ( r1_1 - ( lambertw(r2_1) / log(1-u) ) ) / β
			kbar = kbar > 0 ? kbar : 0
		else
			kbar = 0
		end
		
		push!(kbar_vec_1, kbar)
		
	end

	return [kbarW, kbarV_delta_0, kbarV_delta_1, koptV, kbar_vec_0, kbar_vec_1]
	
end

# ╔═╡ c7571899-7965-44d6-b8e6-1990d2ae08be
function zero_plots(kbar; leg=true, lab=true, notes=true)
	plot(
		collect(0:0.01:1),
		kbar[5],
		label="k̄* (fitness)",
		color=RGBA(0,0.3,1,0.7),
		lw=4,
		legend= leg ? :topleft : false,
		xlabel= lab ? "u (success rate)" : "",
		ylabel= lab ? "k̄*" : "",
	)
	
	plot!(
		collect(0:0.01:1),
		kbar[2],
		label="k̄* (utility) at δ = 1.1",
		color=RGBA(0,1,0,0.5),
		lw=2
	)

	plot!(
		collect(0:0.01:1),
		kbar[3],
		label="k̄* (utility) at δ = 2.1",
		color=RGBA(0,0.7,0,0.7),
		lw=2
	)
	
	plot!(
		collect(0:0.01:1),
		kbar[4],
		label="k̄* (utility) at δ = 4",
		color=RGBA(0,0.4,0,0.7),
		lw=2
	)

	plot!(
		collect(0:0.01:1), kbar[1],
		label="k̄* (utility) at δ = 100",
		color=RGBA(0,0.1,0,0.7),
		lw=2,
	)

end

# ╔═╡ 2572a88f-53e7-4958-b627-da52cd0b8145
function zero_plots_raw(kbar; leg=true, lab=true, notes=true, ylim=33, yax=true)
	
	plot(
		collect(0.01:0.01:1),
		kbar[5],
		label="",
		#line=:dash,
		color=RGBA(0.5,0.5,0.9,0.8),
		lw=2.5,
		yaxis = yax,
		legend = leg ? :topleft : false,
		xlabel = lab ? "u (success rate)" : "",
		ylabel = lab ? "k̄*" : "",
		ylim=[0.0, ylim],
		foreground_color_legend = nothing
	)
	"""
	plot!(
		collect(0.01:0.01:1),
		kbar[6],
		label="k̄* (utility) with δ = ³√B",
		line=:dot,
		color=RGBA(0,0.5,0.5,0.5),
		lw=2
	)
	"""
	plot!(
		collect(0.01:0.01:1),
		kbar[4],
		label="",
		#line=:dash,
		color=RGBA(0.9,0.5,0.5,0.9),
		lw=2.5
	)
end

# ╔═╡ 686c81e2-44af-48e2-9491-a4f770c336b8
function loner(u, B, T, v_0)

	v = Vector{Float64}()
	push!(v, v_0)
	
	for t in 1:T
		if rand(Float64) < u
			v_t = v[t]*B
		else
			v_t = v[t]
		end
		push!(v, v_t)
	end

	return v
	
end

# ╔═╡ fc1b104d-74f0-43ab-b009-8f948c24db08
function sharer(u, s, B, C, β, k_bar; T=200, v_0=1, w_0=0.5)

	v = Vector{Float64}()
	v_t = 0
	push!(v, v_0)

	for t in 1:T
		beta_kbar = floor(β*k_bar)
		k = sum( [rand() < u ? 1 : 0 for i in 1:beta_kbar] )
		lim1 = ( (1-s)*B + s*k*B/k_bar - beta_kbar*C > w_0 )
		lim2 = ( 1 + s*k*B/k_bar - beta_kbar*C > w_0 )
		
		if rand(Float64) < u
			v_t += log( lim1 ? (1-s)*B + s*k*B/k_bar - beta_kbar*C : w_0 )
		else
			v_t += log( lim2 ? 1 + s*k*B/k_bar - beta_kbar*C : w_0 )
		end
		
		push!(v, v_t)
		
	end

	return v

end

# ╔═╡ f5d88364-78d0-44f0-bd6c-14328cdb804a
function geometric_mean_fitness(u, s, B, C, β, k_bar; T=200, v_0=1, w_0=0.5, sharing=true, reps=1000)

	total_v = Vector{}()
	
	for i in 1:reps

		if sharing
			v = sharer(u, s, B, C, β, k_bar; T=T, v_0=v_0, w_0=w_0)
			push!( total_v, last(v)/(v_0*T) )
		else
			v = loner(u, B, T, v_0)
			push!( total_v, last(v)/(v_0*T) )
		end
		
	end

	return total_v
	
end

# ╔═╡ 8566ff2a-bf3a-4dcd-be5e-568bfccd973f
function mean_degree_fitness(u, s, B, C, β; T=500, v_0=1, w_0=0.0001, sharing=true, reps=500, max_kbar=35)

	w_mean = Vector{Float64}(undef, max_kbar)
	w_std = Vector{Float64}(undef, max_kbar)
	
	Threads.@threads for k_bar in 1:max_kbar
		fitness = geometric_mean_fitness(u, s, B, C, β, k_bar; T=T, v_0=v_0, w_0=w_0, sharing=sharing, reps=reps)
		w_mean[k_bar] = mean(fitness)
		w_std[k_bar] = std(fitness)
	end

	return w_mean, w_std

end

# ╔═╡ 0d0a0bcb-0853-4ad0-9ac3-27578252d3d4
optshare(b,k;u=0.5,c=0.1,β=1) = ( ( (β*k*c - 1) / ( 1 - (1-u)^(β*k) )*β ) + ( (1-u)*(b - β*k*c) / ( 1 - u*β*( 1 - (1-u)^(β*k) ) ) ) ) / b

# ╔═╡ 50fe2f05-abbc-44ae-a4ae-6db3e38def4d
function beta_fitness(u, s, B, C, k_bar, T, v_0; sharing=true, reps=5000, beta_max = 10)

	w_mean = Vector{}()
	w_std = Vector{}()
	
	for β in 0.1:0.1:beta_max
		fitness = geometric_mean_fitness(u, s, B, C, β, k_bar, T, v_0; sharing=sharing, reps=reps)
		push!(w_mean, mean(fitness))
		push!(w_std, std(fitness))
	end

	return w_mean, w_std

end

# ╔═╡ 58707353-5fc3-4359-b5d2-11951680391b
function beta_trajectories(u_max, s, B, C, k_bar, T, v_0; sharing=true, reps=5000, beta_max = 10)

	beta_vec = Vector{}()
	
	for u in 0.1:0.05:(u_max)
		beta_mean, beta_std = beta_fitness(u, s, B, C, k_bar,
		T, v_0; sharing=sharing, reps=reps, beta_max = beta_max)
		push!(beta_vec, beta_mean)
	end

	return beta_vec

end

# ╔═╡ 473b097b-6d2d-451b-b6ef-90dddbca3095
function beta_plots(beta_vec; 
	beta_max=10, 
	annotate1=true, 
	annotate2=false, 
	loc=(2.5,3.3), 
	s="s = s*",
)

	x_axis = 0.1:0.1:beta_max
	beta_plot = plot(x_axis, beta_vec[! , 1], color=RGBA(0,0,0,0.7),
		legend=false,
		lw=2,
	)
	
	for i in 1:( ncol(beta_vec) - 1 )
		dc = i*(1/ncol(beta_vec))
		plot!(x_axis, beta_vec[! , i+1], color=RGBA(dc,0,0,0.7), lw=2)
	end

	vline!([1], color=RGBA(0,0,0,0.3))
	ylims!(0, 3.5)

	if annotate1
		annotate!( loc[1], loc[2], text(s, 10) )
	end

	if annotate2
		annotate!([ (7.5, 2.5, text("u = 0.1", 8)), (7.5, 3, text("u = 0.9", 8)) ])
		points = 2.5:0.01:3
		for point in points
			scatter!([6], [point], 
				color=RGBA( (point-2.5)/0.5,0,0,0.1),
				legend=false,
				markerstrokewidth=0,
				#markersize=3
			)
		end
		annotate!(2, 2.7, text("k̄ = 10\nB = 10\nC = 0.1", 8))
		
	end
	
	return beta_plot
end

# ╔═╡ f89f30dd-7a04-4b6d-b0ae-a73edf7a80cf
function calculate_points(s, b, c, β; 
	T=500, 
	v_0=1, 
	w_0=1,
	max_kbar=35,
	reps=100,
	)

	axis = 0:0.01:1|>collect
	k_hat_subvector = Vector{Float64}()
	k_star_subvector = Vector{Float64}()
		
		@progress for i in 1:length(axis)
			
			u = axis[i]
			lone = u*log(b)
			
			mean_vec_u, std_vec_u = mean_degree_fitness(u, s, b, c, β, max_kbar=max_kbar, reps=reps, v_0=v_0, w_0=w_0)

			pos = filter( x -> x > 0, (mean_vec_u .- lone) )
			if length(pos) > 0
				mean_vec_exp = exp.(mean_vec_u)
				max_arg = argmax( mean_vec_exp )
				mean_vec_exp[1:(max_arg - 1)] .= mean_vec_exp[max_arg] + 1
				mean_vec_final = Vector{Float64}()
				for x in mean_vec_exp
					if x - exp(lone) ≥ 0
						push!(mean_vec_final, x)
					else
						break
					end
				end
				fission_point = argmin( exp.(mean_vec_final) )
			else
				max_arg = 0
				fission_point = 0
			end

			push!(k_hat_subvector, fission_point)
			push!(k_star_subvector, max_arg)
			
		end
	
	return [k_hat_subvector, k_star_subvector]
end

# ╔═╡ 855eb815-628c-4d8e-9993-e12e3d69f96b
begin
	u_rec = 0.5
	B_rec = 10
	C_rec = 0.1
	β_rec = 1.0
	s_rec = 0.7
	rng = MersenneTwister()
end

# ╔═╡ 62a59cec-4f54-491c-a743-e096f4433e95
W(kbar; u=u_rec, B=B_rec, C=C_rec, β=1, s=0.75) = kbar ≥ 1 ? u*log( (1-s)*B + (1 - (1 - u)^(β*kbar))*s*u*β*B - β*kbar*C ) + (1-u)*log( 1 + (1 - (1 - u)^(β*kbar))*s*u*β*B - β*kbar*C ) : u*log(B)

# ╔═╡ 011320ee-9010-4dfb-b13e-ad92e5e40e26
W4(kbar; u=u_rec, B=B_rec, C=C_rec, β=1, s=0.5) = kbar ≥ 1 ? ( 1 + (1 - (1 - u)^(β*kbar))*s*u*β*B - β*kbar*C > 0 ?  u*log( (1-s)*B + (1 - (1 - u)^(β*kbar))*s*u*β*B - β*kbar*C ) + (1-u)*log( 1 + (1 - (1 - u)^(β*kbar))*s*u*β*B - β*kbar*C ) : 0 ) : u*log(B)

# ╔═╡ 72ef239b-49fd-4257-bdf4-103aa178e854
x_axis = 1:25

# ╔═╡ 21a8107f-90b6-4334-a0f2-88a457e45970
begin
	u_test2 = 0.2
	b_test2 = 40
	c_test2 = b_test2/100
	s_test2 = 0.5
end

# ╔═╡ 8270e08e-e70c-40e7-90e9-35199867a695
m_test, std_test = mean_degree_fitness(u_test2, s_test2, b_test2, c_test2, 1, max_kbar=40, w_0=0.1)

# ╔═╡ 466d6f8b-3ab9-4bf2-90ce-6c904a209bb5
begin
	u_test = 0.2
	b_test = 30
	c_test = b_test / 100
	s_test = 0.7
end

# ╔═╡ 37367483-c9ce-44c1-a668-03376ab77bbd
khat_test, kstar_test = calculate_points(s_test, b_test, c_test, 1, w_0=0.1, max_kbar=35, reps=1000)

# ╔═╡ 0dd0c2c3-9261-463d-93c6-ddcffdcd51a9
begin
	khat_test2, kstar_test2 = calculate_points(s_test, 5, 0.05, 1, w_0=0.1, max_kbar=35, reps=1000)
	khat_test3, kstar_test3 = calculate_points(s_test, 15, 0.15, 1, w_0=0.1, max_kbar=35, reps=1000)
end

# ╔═╡ 08a8fa0b-42fe-4a9e-ab7c-8d6d808a5d2c
begin
	
	x = 0.01:0.001:1
	f(x) = log( -0.25/( log(1-x) ) ) / ( log(1-x) )

	cb_ratios = sort( collect(0.002:0.01:0.1), rev=true )
	y_step = 0.2
	
	kopt_plot = plot( x, f.(x), 
		ylim=(0,20), 
		legend=false, 
		color=RGBA(0,0.1,0.1,0.0),
		xlab="probability of success (u)",
		ylab="optimal mean degree (k̂) at β = 1"
	)

	count = 1
	for cb in cb_ratios
		yh = 11 + count*y_step
		g(x) = log( -cb/( x * log(1-x) ) ) / ( log(1-x) )
		col = RGBA((1 - cb/maximum(cb_ratios))*1,0.1,0.01,0.9)
		plot!( x, g.(x), 
			color=col,
			lw=2
		)
		scatter!([0.75], [yh], 
			color=col,
			markersize=5,
			markerstrokewidth=0
		)
		global count += 1
	end

	hline!( [1], color="dark blue", lw=1.5 )
	annotate!( [
		(0.78, 11, text("0.1", 6)),
		(0.785, 13.5, text("0.002", 6)),
		(0.75, 15.2, text("connection cost /\nexpected shared benefit\nratio", 8))
	]
	)
	
	kopt_plot
end

# ╔═╡ 50931fc3-7bae-4b46-8711-13ba4bf2a864
begin

	mksize = 3
	ylimit = 31
	
	zeros1 = zero_plots_raw( 
		zero_difs_raw(B=5, C=0.05, s=s_test), 
		ylim=ylimit, 
		lab=false 
	)
	scatter!(
		0:0.01:1, 
		khat_test2, 
		alpha=0.3, 
		markersize=mksize, 
		label="", 
		color=RGBA(0.5,0.5,0.9,0.8)
	)
	plot!(
		1:0, 
		label="Maximum mean degree (k̄*) \n with δ = 0.1 + (B - 1)^0.38",
		#line=:dash,
		color=RGBA(0.5,0.5,0.9,0.8),
		lw=2.5
	)
	scatter!(
		0:0.01:1, 
		kstar_test2, 
		alpha=0.3, 
		markersize=mksize, 
		label="", 
		color=RGBA(0.9,0.5,0.5,0.9)
	)
	plot!(
		1:0, 
		label="Optimal mean degree (k̂)",
		#line=:dash,
		color=RGBA(0.9,0.5,0.5,0.9),
		lw=2.5
	)
	
	#lel(x) = ( 1 - x^(0.1 + (5 - 1)^0.38) ) * (s_test * x * 5) / 0.05
	#plot!(
	#	collect(0:0.01:1),
	#	lel.(collect(0:0.01:1))
	#)


	
	zeros2 = zero_plots_raw( 
		zero_difs_raw(B=15, C=0.15, s=s_test),  
		lab=false, 
		leg=false,
		ylim=ylimit,
		yax=nothing
	)
	scatter!(
		0:0.01:1, 
		khat_test3, 
		alpha=0.3, 
		markersize=mksize, 
		label="k̄*", 
		color=RGBA(0.5,0.5,0.9,0.8)
	)
	scatter!(
		0:0.01:1, 
		kstar_test3, 
		alpha=0.3, 
		markersize=mksize, 
		label="k̃", 
		color=RGBA(0.9,0.5,0.5,0.9)
	)
	#annotate!([(0.6, 32, text("s = 0.7 \n β = 1", 7))])


	
	zeros3 = zero_plots_raw( 
		zero_difs_raw(B=b_test, C=c_test, s=s_test), 
		lab=false, 
		leg=false,
		ylim=ylimit,
		yax=nothing
	)
	scatter!(
		0:0.01:1, 
		khat_test, 
		alpha=0.3, 
		markersize=mksize, 
		label="k̄*", 
		color=RGBA(0.5,0.5,0.9,0.8)
	)
	scatter!(
		0:0.01:1, 
		kstar_test, 
		alpha=0.3, 
		markersize=mksize, 
		label="k̃", 
		color=RGBA(0.9,0.5,0.5,0.9)
	)


	zeros_plot = plot(
		zeros1, zeros2, zeros3, 
		layout=(1,3), 
		size=(700,400), 
		link=:all,
		title = ["Surplus (B) = 5 \n Per-connection cost (C) = 0.05" "Surplus (B) = 15 \n Per-connection cost (C) = 0.15" "Surplus (B) = 30 \n Per-connection cost (C) = 0.3"],
		titlefontsize=9
	)
	
end

# ╔═╡ c43f86cf-899b-4538-8685-b36a7e0eabe5
savefig(zeros_plot, "utility_vs_sim.pdf")

# ╔═╡ 1feca48b-c134-4a1e-9a29-af2adf60bd3d
begin
	mean_vec, std_vec = mean_degree_fitness(u_rec, 0.75, B_rec, C_rec, β_rec)
	#mean_vec2, std_vec2 = mean_degree_fitness(u_rec, s2, B_rec, C_rec, β_rec)
	#mean_vec3, std_vec3 = mean_degree_fitness(u_rec, s3, B_rec, C_rec, β_rec)
	#mean_vec4, std_vec4 = mean_degree_fitness(u_rec, s4, B_rec, C_rec, β_rec)
	mean_vec5, std_vec5 = mean_degree_fitness(u_rec, 0.5, B_rec, C_rec, β_rec)
end

# ╔═╡ afd92507-eeba-4f0f-b6fc-034fd46edf6f
begin

	ylim=[1.15, 1.6]
	
	k_plot_sim = plot(x_axis, mean_vec[x_axis],
		xlabel="mean degree (k̄), simulated", 
		ylabel="geometric mean fitness (log scale)",
		label="s = 0.9",
		legend=false,
		lw=1.5,
		color=RGBA(0,0.7,1,0.7),
		size=(550,450),
		ylim=ylim
	)
	plot!(x_axis, mean_vec5[x_axis],
		label="s = 0.5",
		lw=1.5,
		color=RGBA(0,0.1,0.1,0.9)
	)
	hline!([u_rec*log(B_rec)],
		label="loner",
		lw=1.5,
		color=RGBA(0.7,0,0,0.6)
	)

	k_plot_math = plot(x_axis, W,
		xlabel="approximation", 
		#ylabel="geometric mean fitness (log scale)",
		label="s = 0.9",
		legend=:topright,
		foreground_color_legend = nothing,
		lw=1.5,
		color=RGBA(0,0.6,1,0.7),
		size=(550,450),
		ylim=ylim
	)
	
	plot!(x_axis, W4,
		label="s = 0.5",
		lw=1.5,
		color=RGBA(0,0.1,0.1,0.9)
	)
	
	hline!([u_rec*log(B_rec)],
		label="loner",
		lw=1.5,
		color=RGBA(0.7,0,0,0.6)
	)
	
	annotate!(20, 1.45, text("B = 10\nC = 0.1\nu = 0.5\nβ = 1", :black, :right, 8))
	
	"paired plots"
end

# ╔═╡ 88abbec4-8cc2-46ae-870c-ccea1a579e2c
begin
	
	plot(x_axis, mean_vec[x_axis],
		xlabel="mean degree of sharing network (k̄)", 
		ylabel="geometric mean fitness (log scale)",
		label="s = 0.9",
		legend=:topright,
		lw=1.5,
		color=RGBA(0,0.7,1,0.7),
		size=(550,450),
		ylim=ylim
	)
	plot!(x_axis, mean_vec5[x_axis],
		label="s = 0.5",
		lw=1.5,
		color=RGBA(0,0.1,0.1,0.9)
	)
	hline!([u_rec*log(B_rec)],
		label="loner",
		lw=1.5,
		color=RGBA(0.7,0,0,0.6)
	)

	annotate!(23, 1.4, text("B = 10\nC = 0.1\nu = 0.5\nβ = 1", :black, :right, 8))
	
end

# ╔═╡ 8dc1936e-7ab7-418e-ba00-e31a6308ea3b
paired_plots = plot(k_plot_sim, k_plot_math, size=(600,400), link=:all)

# ╔═╡ 91035665-c615-4374-ad3c-e09f9a241e52
begin
	X = 1:0.01:20
	Y = 1:0.01:20

	Z1 = [clamp(optshare(x, y, β=0.75, u=0.8), 0, 1) for x in X, y in Y]
	pz1 = plot(
		X, Y, Z1, st=:contourf,
		#title = "β = 0.75",
		#xlabel="β = 0.75",
		ylabel="u = 0.8",
		xaxis=nothing,
		#yaxis=nothing,
		clims=(0,1),
		legend=false
	)

	Z2 = [clamp(optshare(x, y, u=0.8), 0, 1) for x in X, y in Y]
	pz2 = plot(
		X, Y, Z2, st=:contourf,
		#title = "β = 1",
		#xlabel="β = 1",
		#ylabel="k̄",
		xaxis=nothing,
		yaxis=nothing,
		clims=(0,1),
		legend=false
	)

	Z3 = [clamp(optshare(x, y, β=1.25, u=0.8), 0, 1) for x in X, y in Y]
	pz3 = plot(
		X, Y, Z3, st=:contourf,
		#xlabel="β = 1.25",
		#title = "β = 1.25",
		#yguide="u = 0.8",
		ymirror=true,
		xaxis=nothing,
		yaxis=nothing,
		clims=(0,1),
		legend=false
	)
	
	Z4 = [clamp(optshare(x, y, β=0.75), 0, 1) for x in X, y in Y]
	pz4 = plot(
		X, Y, Z4, st=:contourf,
		showaxis=false,
		#xlabel="β = 0.75",
		ylabel="u = 0.5",
		xtickfontsize=8,
		xaxis=nothing,
		#yaxis=nothing,
		clims=(0,1),
		legend=false
	)

	Z5 = [clamp(optshare(x, y), 0, 1) for x in X, y in Y]
	pz5 = plot(
		X, Y, Z5, st=:contourf,
		#xlabel="β = 1",
		#ylabel="k̄",
		xaxis=nothing,
		yaxis=nothing,
		clims=(0,1),
		legend=false
	)
	
	Z6 = [clamp(optshare(x, y, β=1.25), 0, 1) for x in X, y in Y]
	pz6 = plot(
		X, Y, Z6, st=:contourf,
		#xlabel="β = 1.25",
		#ylabel="u = 0.5",
		ymirror=true,
		xaxis=nothing,
		yaxis=nothing,
		clims=(0,1),
		legend=false
	)

	Z7 = [clamp(optshare(x, y, β=0.75, u=0.2), 0, 1) for x in X, y in Y]
	pz7 = plot(
		X, Y, Z7, st=:contourf,
		xlabel="β = 0.75",
		ylabel="u = 0.2",
		#yaxis=nothing,
		clims=(0,1),
		legend=false
	)

	Z8 = [clamp(optshare(x, y, u=0.2), 0, 1) for x in X, y in Y]
	pz8 = plot(
		X, Y, Z8, st=:contourf,
		xlabel="β = 1",
		#ylabel="k̄",
		yaxis=nothing,
		clims=(0,1),
		legend=false
	)

	Z9 = [clamp(optshare(x, y, β=1.25, u=0.2), 0, 1) for x in X, y in Y]
	pz9 = plot(
		X, Y, Z9, st=:contourf,
		xlabel="β = 1.25",
		#ylabel="u = 0.2",
		ymirror=true,
		yaxis=nothing,
		clims=(0,1),
		legend=false
	)

	

	Z = [clamp(optshare(x, y, β=0.75), 0, 1) for x in X, y in Y]
	pz = plot(
		X, Y, Z, st=:contourf,
		#xlabel="B",
		#ylabel="k̄",
		clims=(0,1),
	)
	
end

# ╔═╡ 91f11d58-acb1-432c-a0e7-c64124044bbf
s_opt_contour = plot( 
	pz1,pz2,pz3, 
	pz4,pz5,pz6, 
	pz7,pz8,pz9, 
	layout=(3,3), 
	main="optimal sharing norm (s*)", 
	size=(650,650), 
	linewidth=0, 
	xtickfontsize=9, 
	ytickfontsize=9, 
	xguidefontsize=12, 
	yguidefontsize=12, 
	titlefontsize=12,
	levels=0:0.05:1
)

# ╔═╡ a791f73d-4caf-4cda-828b-c1867e1e4d3b
begin
	Zk(b) = [clamp(optshare(b, y, β=1, u=0.2), 0, 1) for y in 1:15]
	
	zk_plot = plot(1:15, Zk(2), legend=:bottomright, label="", xguide="mean degree of sharing network (k̄)", yguide="optimal sharing norm (s*)", color=RGBA(0.2,0.5,0.9,0.7))
	scatter!(1:15, Zk(2), label = "B = 2", alpha=0.7, color=RGBA(0.2,0.5,0.9,0.7))

	plot!(1:15, Zk(4), label="", color=RGBA(0.5,0.2,0.7,0.7))
	scatter!(1:15, Zk(4), label = "B = 4", alpha=0.7)
	
	plot!(1:15, Zk(7), label="", color=RGBA(0.8,0.2,0.2,0.7))
	scatter!(1:15, Zk(7), label = "B = 7", alpha=0.7, color=RGBA(0.8,0.2,0.2,0.7))

	plot!(1:15, Zk(10), label="", color=RGBA(0.7,0.5,0.2,0.7))
	scatter!(1:15, Zk(10), label = "B = 10", alpha=0.7, color=RGBA(0.7,0.5,0.2,0.7))

	plot!(1:15, Zk(15), label="", color=RGBA(0.6,0.6,0.1,0.3))
	scatter!(1:15, Zk(15), label = "B = 15", alpha=0.5, color=RGBA(0.6,0.6,0.1,0.3))
end

# ╔═╡ a8cc023a-cc71-4d49-b715-4d455828cf45
begin
	Zb(k) = [clamp(optshare(y, k, β=1, u=0.2), 0, 1) for y in 1:15]
		
	zb_plot = plot(1:15, Zb(2), legend=:bottomright, label="", xguide="surplus (B)", yguide="optimal sharing norm (s*)", color=RGBA(0.2,0.5,0.9,0.7))
	scatter!(1:15, Zb(2), label = "k̄ = 2", alpha=0.7, color=RGBA(0.2,0.5,0.9,0.7))
	
	plot!(1:15, Zb(4), label="", color=RGBA(0.5,0.2,0.7,0.7))
	scatter!(1:15, Zb(4), label = "k̄ = 4", alpha=0.7)
		
	plot!(1:15, Zb(7), label="", color=RGBA(0.8,0.2,0.2,0.7))
	scatter!(1:15, Zb(7), label = "k̄ = 7", alpha=0.7, color=RGBA(0.8,0.2,0.2,0.7))
	
	plot!(1:15, Zb(10), label="", color=RGBA(0.7,0.5,0.2,0.7))
	scatter!(1:15, Zb(10), label = "k̄ = 10", alpha=0.7, color=RGBA(0.7,0.5,0.2,0.7))
	
	plot!(1:15, Zb(15), label="", color=RGBA(0.6,0.6,0.1,0.3))
	scatter!(1:15, Zb(15), label = "k̄ = 15", alpha=0.5, color=RGBA(0.6,0.6,0.1,0.3))
end

# ╔═╡ b37a1ea6-1760-4795-bd32-c99aa3c1a2db
function productlog_u(
	;B=10,
	C=0.25,
	β=1, 
	s=7/10, 
	δ_0=0.1 + (B-1)^(0.38), δ_1=B^(1/3),
)

	kbar_vec = Vector{Float64}()
	plogs = Vector{Float64}()
	kbarsicles = Vector{Float64}()
	
	for u in (0.001:0.01:1)

		r1 = ( u*(β - u^(δ_0-1))*s*B / C )
		r2 = ( log(1-u) * u * (1-u)^(r1) * s*β*B / C )
		
		if r2 > -1/exp(1)
			plog = lambertw(r2)
			kbar = ( r1 - ( plog / log(1-u) ) ) / β
			kbar = kbar > 0 ? kbar : 0
			kbarsicle = r1 / β
		else
			plog = 0
			kbar = 0
			kbarsicle = r1 / β
		end
		
		push!(kbar_vec, kbar)
		push!(plogs, plog)
		push!(kbarsicles, kbarsicle)
		
	end

	return (kbar_vec, plogs, kbarsicles)
	
end

# ╔═╡ 5e2124e2-c59f-4ec6-a1fc-17211b30c759
begin
	kbarcito, plog, kbarsicle = productlog_u(B=5, C=0.05)
	full_plog = plog ./ log.( 1 .- collect(0.001:0.01:1) )
end

# ╔═╡ 4de61246-6ca7-4f12-8d15-cd7630826a3e
begin
	kbarcito_plot = plot(
		0.001:0.01:1, kbarcito, 
		lw=2, alpha=0.5,
		ylab = "maximum mean degree (k̄*)",
		xlab = "probability of success (u)",
		label = "with boundaries"
	)
	plot!(
		0.001:0.01:1, kbarsicle,
		lw=2, alpha=0.5,
		label = "no boundaries"
	)
end

# ╔═╡ bd239232-4c9f-41a1-8cb7-e738f04e98cb
scatter(0.001:0.01:1, kbarsicle)

# ╔═╡ 8f9ff7b7-049a-4c50-8afa-0c5c4ddde9c4
"""
begin
	kstar1 = zero_difs(B=5, C=0.05)
	kstar2 = zero_difs(B=10, C=0.1)
	kstar3 = zero_difs(B=30, C=0.3)
	kstar4 = zero_difs(B=60, C=0.6)
end
"""

# ╔═╡ aec7dffe-1605-4aef-9a12-4f6a633f2bcf
"""
begin
	p1 = zero_plots(kstar1, lab=false, )
	p2 = zero_plots(kstar2, leg=false, lab=false )
	p3 = zero_plots(kstar3, leg=false, lab=false, )
	p4 = zero_plots(kstar4, leg=false, lab=false, )
	plot( 
		p1,p2,p3,p4,
		layout=(2,2), 
		size=(800,600),
		xlabel = ["Surplus (B) = 5; Cost (C) = 0.05" "Surplus (B) = 10; Cost (C) = 0.1" "Surplus (B) = 30; Cost (C) = 0.3" "Surplus (B) = 60; Cost (C) = 0.6"],
		plot_title = "Network carrying capacity (k̄*) as a function of security (u)",
	)
end
"""

# ╔═╡ b84c0c50-06b0-4955-9c8f-ff6db05f1177
"""
begin
	k_vector = Vector()
	benefits_costs = [
		[4,0.04],[7,0.07],
		[10,0.1],[15,0.15],
		[20,0.2],[25,0.25],
		[35,0.35],[55,0.55],
	]
	u_axis = 0:0.01:1|>collect
	
	for bc in benefits_costs
		k_hat_subvector = Vector{Float64}(undef, u_axis|>length)
		k_star_subvector = Vector{Float64}(undef, u_axis|>length)
		
		@progress for i in 1:length(u_axis)
			
			u = u_axis[i]
			lone = u*log(bc[1])
			
			mean_vec_u, std_vec_u = mean_degree_fitness(u, s_rec, bc[1], bc[2], 1, max_kbar=37)

			pos = filter( x -> x > 0, (mean_vec_u .- lone) )
			if length(pos) > 0
				mean_vec_exp = exp.(mean_vec_u)
				max_arg = argmax( mean_vec_exp )
				mean_vec_exp[1:(max_arg - 1)] .= mean_vec_exp[max_arg] + 1
				mean_vec_final = Vector{Float64}()
				for x in mean_vec_exp
					if x - exp(lone) ≥ 0
						push!(mean_vec_final, x)
					else
						break
					end
				end
				fission_point = argmin( exp.(mean_vec_final) )
			else
				max_arg = 0
				fission_point = 0
			end

			k_hat_subvector[i] = fission_point
			k_star_subvector[i] = max_arg
			
		end
		push!(k_vector, [k_star_subvector, k_hat_subvector])
	end
end
"""

# ╔═╡ 2e7deb04-48fb-4925-a89f-aadcb9954118
"""
begin
	kstar1_raw = zero_difs_raw(B=benefits_costs[1][1], C=benefits_costs[1][2], s=0.7)
	kstar2_raw = zero_difs_raw(B=benefits_costs[2][1], C=benefits_costs[2][2], s=0.7)
	kstar3_raw = zero_difs_raw(B=benefits_costs[3][1], C=benefits_costs[3][2], s=0.7)
	kstar4_raw = zero_difs_raw(B=benefits_costs[4][1], C=benefits_costs[4][2], s=0.7)
	kstar5_raw = zero_difs_raw(B=benefits_costs[5][1], C=benefits_costs[5][2], s=0.7)
	kstar6_raw = zero_difs_raw(B=benefits_costs[6][1], C=benefits_costs[6][2], s=0.7)
	kstar7_raw = zero_difs_raw(B=benefits_costs[7][1], C=benefits_costs[7][2], s=0.7)
	kstar8_raw = zero_difs_raw(B=benefits_costs[8][1], C=benefits_costs[8][2], s=0.7)
end
"""

# ╔═╡ e0ef219e-c837-4257-9dc7-73ed225306f6
"""
begin
	
	p1_raw = zero_plots_raw(kstar1_raw, lab=false, )
	scatter!(0:0.01:1, k_vector[1][2], alpha=0.15, markersize=2, lab="k̄* (simulated)", color=RGBA(0,0.5,0.5,0.5))
	scatter!(0:0.01:1, k_vector[1][1], alpha=0.15, markersize=2, lab="k̃ (simulated)", color=RGBA(1,0,0,0.4))
	
	p2_raw = zero_plots_raw(kstar2_raw, leg=false, lab=false )
	scatter!(0:0.01:1, k_vector[2][2], alpha=0.15, markersize=2, color=RGBA(0,0.5,0.5,0.5))
	scatter!(0:0.01:1, k_vector[2][1], alpha=0.15, markersize=2, color=RGBA(1,0,0,0.4))
	
	p3_raw = zero_plots_raw(kstar3_raw, leg=false, lab=false, )
	scatter!(0:0.01:1, k_vector[3][2], alpha=0.15, markersize=2, color=RGBA(0,0.5,0.5,0.5))
	scatter!(0:0.01:1, k_vector[3][1], alpha=0.15, markersize=2, color=RGBA(1,0,0,0.4))
	
	p4_raw = zero_plots_raw(kstar4_raw, leg=false, lab=false, )
	scatter!(0:0.01:1, k_vector[4][2], alpha=0.15, markersize=2, color=RGBA(0,0.5,0.5,0.5))
	scatter!(0:0.01:1, k_vector[4][1], alpha=0.15, markersize=2, color=RGBA(1,0,0,0.4))

	p5_raw = zero_plots_raw(kstar5_raw, leg=false, lab=false, )
	scatter!(0:0.01:1, k_vector[5][2], alpha=0.15, markersize=2, color=RGBA(0,0.5,0.5,0.5))
	scatter!(0:0.01:1, k_vector[5][1], alpha=0.15, markersize=2, color=RGBA(1,0,0,0.4))

	p6_raw = zero_plots_raw(kstar6_raw, leg=false, lab=false, )
	scatter!(0:0.01:1, k_vector[6][2], alpha=0.15, markersize=2, color=RGBA(0,0.5,0.5,0.5))
	scatter!(0:0.01:1, k_vector[6][1], alpha=0.15, markersize=2, color=RGBA(1,0,0,0.4))

	p7_raw = zero_plots_raw(kstar7_raw, leg=false, lab=false, )
	scatter!(0:0.01:1, k_vector[7][2], alpha=0.15, markersize=2, color=RGBA(0,0.5,0.5,0.5))
	scatter!(0:0.01:1, k_vector[7][1], alpha=0.15, markersize=2, color=RGBA(1,0,0,0.4))

	p8_raw = zero_plots_raw(kstar8_raw, leg=false, lab=false, )
	scatter!(0:0.01:1, k_vector[8][2], alpha=0.15, markersize=2, color=RGBA(0,0.5,0.5,0.5))
	scatter!(0:0.01:1, k_vector[8][1], alpha=0.15, markersize=2, color=RGBA(1,0,0,0.4))
	
	plot( 
		p1_raw, p2_raw, p3_raw, p4_raw, p5_raw, p6_raw, p7_raw, p8_raw,
		layout=(2,4), 
		size=(800,600),
		xlabel = ["B = 4; C = 0.04" "B = 7; C = 0.07" "B = 10; C = 0.1" "B = 15; C = 0.15" "B = 20; C = 0.2" "B = 25; C = 0.25" "B = 35; C = 0.35" "B = 55; C = 0.55"],
		#plot_title = "k̄* and k̄̂ as functions of u",
		link=:all
	)
end
"""

# ╔═╡ 9b02f8d7-93c4-4ecb-a549-5fd70daa16bb
"""
begin
	center_traj = plot(x_axis, mean_vec5[x_axis],
		xlabel="mean degree of sharing network (k̄)", 
		#ylabel="geometric mean fitness (log scale)",
		label="G(Sharer) at s = 1",
		legend=false,
		lw=2.5,
		color=RGBA(0,0.3,0.1,0.3),
		size=(600,500)
	)
	plot!(x_axis, mean_vec[x_axis],
		label="G(Sharer) at s = s*",
		lw=2.5,
		color=RGBA(0,0.3,1,0.7)
	)
	plot!(x_axis, mean_vec2[x_axis],
		label="G(Sharer) at s = 1/2",
		lw=2.5,
		color=RGBA(0,0.2,0.5,0.7)
	)
	plot!(x_axis, mean_vec3[x_axis],
		label="G(Sharer) at s = 1/3",
		lw=2.5,
		color=RGBA(0,0.1,0.1,0.7)
	)
	plot!(x_axis, mean_vec4[x_axis],
		label="G(Sharer) at s = 1/4",
		lw=2.5,
		color=RGBA(0,0.1,0.1,0.9)
	)
	hline!([u_rec*log(B_rec)],
		label="G(Loner)",
		lw=2,
		color=RGBA(0.7,0,0,0.6)
	)
	
	annotate!(5, 0.5, text("B = 10\nC = 0.15\nu = 0.5\nβ = 1", :black, :right, 8))
end
"""

# ╔═╡ dbf98495-6992-4b7c-a030-694c88c1a7fb
"""
begin
	β_rec_low = 0.7

	mean_vec5_low, std_vec5_low = mean_degree_fitness(u_rec, 1.0, B_rec, C_rec, β_rec_low)
	mean_vec_low, std_vec_low = mean_degree_fitness(u_rec, (B_rec - 1)/B_rec, B_rec, C_rec, β_rec_low)
	mean_vec2_low, std_vec2_low = mean_degree_fitness(u_rec, 1/2, B_rec, C_rec, β_rec_low)
	mean_vec3_low, std_vec3_low = mean_degree_fitness(u_rec, 1/3, B_rec, C_rec, β_rec_low)
	mean_vec4_low, std_vec4_low = mean_degree_fitness(u_rec, 1/4, B_rec, C_rec, β_rec_low)
end
"""

# ╔═╡ 6f4774dd-c729-427c-95de-1c2ce15425c2
"""
begin
	center_traj_low = plot(x_axis, mean_vec5_low[x_axis],
		xlabel="mean degree of sharing network (k̄)", 
		#ylabel="geometric mean fitness (log scale)",
		label="G(Sharer) at s = 1",
		legend=false,
		lw=2.5,
		color=RGBA(0,0.3,0.1,0.3),
		size=(600,500)
	)
	plot!(x_axis, mean_vec_low[x_axis],
		label="G(Sharer) at s = s*",
		lw=2.5,
		color=RGBA(0,0.3,1,0.7)
	)
	plot!(x_axis, mean_vec2_low[x_axis],
		label="G(Sharer) at s = 1/2",
		lw=2.5,
		color=RGBA(0,0.2,0.5,0.7)
	)
	plot!(x_axis, mean_vec3_low[x_axis],
		label="G(Sharer) at s = 1/3",
		lw=2.5,
		color=RGBA(0,0.1,0.1,0.7)
	)
	plot!(x_axis, mean_vec4_low[x_axis],
		label="G(Sharer) at s = 1/4",
		lw=2.5,
		color=RGBA(0,0.1,0.1,0.9)
	)
	hline!([u_rec*log(B_rec)],
		label="G(Loner)",
		lw=2,
		color=RGBA(0.7,0,0,0.6)
	)
	
	annotate!(5, 0.5, text("β = 0.7", :black, :right, 8))
end
"""

# ╔═╡ c0072be7-60bc-4764-9026-51a7d46c64b1
"""
begin
	β_rec_high = 1.3

	mean_vec5_high, std_vec5_high = mean_degree_fitness(u_rec, 1.0, B_rec, C_rec, β_rec_high)
	mean_vec_high, std_vec_high = mean_degree_fitness(u_rec, (B_rec - 1)/B_rec, B_rec, C_rec, β_rec_high)
	mean_vec2_high, std_vec2_high = mean_degree_fitness(u_rec, 1/2, B_rec, C_rec, β_rec_high)
	mean_vec3_high, std_vec3_high = mean_degree_fitness(u_rec, 1/3, B_rec, C_rec, β_rec_high)
	mean_vec4_high, std_vec4_high = mean_degree_fitness(u_rec, 1/4, B_rec, C_rec, β_rec_high)
end
"""

# ╔═╡ 68c159d6-75b8-4c9b-99c5-0206fe958a15
"""
begin
	center_traj_high = plot(x_axis, mean_vec5_high[x_axis],
		xlabel="mean degree of sharing network (k̄)", 
		#ylabel="geometric mean fitness (log scale)",
		label="G(Sharer) at s = 1",
		legend=:topright,
		lw=2.5,
		color=RGBA(0,0.3,0.1,0.3),
		size=(600,500)
	)
	plot!(x_axis, mean_vec_high[x_axis],
		label="G(Sharer) at s = s*",
		lw=2.5,
		color=RGBA(0,0.3,1,0.7)
	)
	plot!(x_axis, mean_vec2_high[x_axis],
		label="G(Sharer) at s = 1/2",
		lw=2.5,
		color=RGBA(0,0.2,0.5,0.7)
	)
	plot!(x_axis, mean_vec3_high[x_axis],
		label="G(Sharer) at s = 1/3",
		lw=2.5,
		color=RGBA(0,0.1,0.1,0.7)
	)
	plot!(x_axis, mean_vec4_high[x_axis],
		label="G(Sharer) at s = 1/4",
		lw=2.5,
		color=RGBA(0,0.1,0.1,0.9)
	)
	hline!([u_rec*log(B_rec)],
		label="G(Loner)",
		lw=2,
		color=RGBA(0.7,0,0,0.6)
	)
	
	annotate!(5, 0.5, text("β = 1.3", :black, :right, 8))
end
"""

# ╔═╡ 824763d1-9ea4-465e-8ebc-7345c5ffe579
"""
plot(
	center_traj_low, 
	center_traj, 
	center_traj_high,
	link=:all,
	layout=(1,3),
	size=(1500,600)
)
"""

# ╔═╡ bd084197-8c32-4380-b403-7314a603ad09
begin
	"""
	begin
		beta_vec_altF = beta_trajectories(0.9, 9/10, 10, 0.1, 10, 200, 1)
		beta_vec_altF2 = beta_trajectories(0.9, 1/2, 10, 0.1, 10, 200, 1)
		beta_vec_altF3 = beta_trajectories(0.9, 1/3, 10, 0.1, 10, 200, 1)
		beta_vec_altF4 = beta_trajectories(0.9, 1/4, 10, 0.1, 10, 200, 1)
		beta_vec_altF5 = beta_trajectories(0.9, 1/5, 10, 0.1, 10, 200, 1)
		beta_vec_altF6 = beta_trajectories(0.9, 1/8, 10, 0.1, 10, 200, 1)
	
		beta_df1 = DataFrame(beta_vec_altF, string.(1:17))
		beta_df2 = DataFrame(beta_vec_altF2, string.(1:17))
		beta_df3 = DataFrame(beta_vec_altF3, string.(1:17))
		beta_df4 = DataFrame(beta_vec_altF4, string.(1:17))
		beta_df5 = DataFrame(beta_vec_altF5, string.(1:17))
		beta_df6 = DataFrame(beta_vec_altF6, string.(1:17))
	
		CSV.write("data/beta_df1.csv", beta_df1)
		CSV.write("data/beta_df2.csv", beta_df2)
		CSV.write("data/beta_df3.csv", beta_df3)
		CSV.write("data/beta_df4.csv", beta_df4)
		CSV.write("data/beta_df5.csv", beta_df5)
		CSV.write("data/beta_df6.csv", beta_df6)
	end
	"""
	
	md"""
	###### uncomment and run to simulate and save data for "Fitness as a function of β"
	"""
end

# ╔═╡ 07332e59-4298-4dd4-86c0-7a4959b31b4d
"""
begin
	beta_df1 = DataFrame(CSV.File("data/beta_df1.csv"))
	beta_df2 = DataFrame(CSV.File("data/beta_df2.csv"))
	beta_df3 = DataFrame(CSV.File("data/beta_df3.csv"))
	beta_df4 = DataFrame(CSV.File("data/beta_df4.csv"))
	beta_df5 = DataFrame(CSV.File("data/beta_df5.csv"))
	beta_df6 = DataFrame(CSV.File("data/beta_df6.csv"))
	
	beta_p1F = beta_plots(beta_df1)
	beta_p2F = beta_plots(beta_df2, s = "s = 1/2")
	beta_p3F = beta_plots(beta_df3, s = "s = 1/3")
	beta_p4F = beta_plots(beta_df4, s = "s = 1/4")
	beta_p5F = beta_plots(beta_df5, s = "s = 1/5")
	beta_p6F = beta_plots(beta_df6, s = "s = 1/8", annotate2=true)

	beta_plotF = plot(
		beta_p1F,
		beta_p2F,
		beta_p3F,
		beta_p4F,
		beta_p5F,
		beta_p6F,
		layout=grid(2,3),
		link=:all,
		size=(900,700),
		plot_title="(log) fitness as a function of position with respect to mean degree (β)"
	)
end
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
LambertW = "984bce1d-4616-540c-a9ee-88d1112d94c9"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
ProgressLogging = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Roots = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
CSV = "~0.10.4"
Colors = "~0.12.8"
DataFrames = "~1.3.4"
Distributions = "~0.25.68"
LambertW = "~0.4.5"
Plots = "~1.27.5"
ProgressLogging = "~0.1.4"
Roots = "~2.0.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings"]
git-tree-sha1 = "873fb188a4b9d76549b81465b1f75c82aaf59238"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.4"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "9950387274246d08af38f6eef8cb5480862a435f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.14.0"

[[ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "12fc73e5e0af68ad3137b886e3f7c1eacfca2640"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.17.1"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[CommonSolve]]
git-tree-sha1 = "68a0743f578349ada8bc911a5cbd5a2ef6ed6d1f"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.0"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "b153278a25dd42c65abbf4e62344f9d22e59191b"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.43.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "daa21eb85147f72e41f6352a57fccea377e310a9"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.3.4"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "334a5896c1534bb1aa7aa2a642d30ba7707357ef"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.68"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "129b104185df66e408edd6625d480b7f9e9823a0"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.18"

[[FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "3399bbad4c9e9a2fd372a54d7b67b3c7121b6402"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.3"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "51d2dfe8e590fbd74e7a842cf6d13d8a2f45dc01"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.6+0"

[[GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "af237c08bda486b74318c8070adb96efa6952530"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.64.2"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "cd6efcf9dc746b06709df14e462f0a3fe0786b1e"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.64.2+0"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "83ea630384a13fc4f002b77690bc0afeb4255ac9"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.2"

[[Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions", "Test"]
git-tree-sha1 = "709d864e3ed6e3545230601f94e11ebc65994641"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.11"

[[IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "61feba885fac3a407465726d0c330b3055df897f"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.1.2"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "91b5dcf362c5add98049e6c29ee756910b03051d"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.3"

[[InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[LambertW]]
git-tree-sha1 = "2d9f4009c486ef676646bca06419ac02061c088e"
uuid = "984bce1d-4616-540c-a9ee-88d1112d94c9"
version = "0.4.5"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "6f14549f7760d84b2db7a9b10b88cd3cc3025730"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.14"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "c9551dd26e31ab17b86cbd00c2ede019c08758eb"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+1"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a970d55c2ad8084ca317a4658ba6ce99b7523571"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.12"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[NaNMath]]
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ab05aa4cc89736e95915b01e7279e61b1bfe33b8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.14+0"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "cf494dca75a69712a72b80bc48f59dcf3dea63ec"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.16"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "621f4f3b4977325b9128d5fae7a8b4829a0c2222"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.4"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "8162b2f8547bc23876edd0c5181b27702ae58dce"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.0.0"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "bb16469fd5224100e422f0b027d26c5a25de1200"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.2.0"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "88ee01b02fba3c771ac4dce0dfc4ecf0cb6fb772"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.27.5"

[[PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "d3538e7f8a790dc8903519090857ef8e1283eecd"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.5"

[[PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[ProgressLogging]]
deps = ["Logging", "SHA", "UUIDs"]
git-tree-sha1 = "80d919dee55b9c50e8d9e2da5eeafff3fe58b539"
uuid = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
version = "0.1.4"

[[Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

[[QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "dc1e451e15d90347a7decc4221842a022b011714"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.5.2"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "cdbd3b1338c72ce29d9584fdbe9e9b70eeb5adca"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.1.3"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[Roots]]
deps = ["CommonSolve", "Printf", "Setfield"]
git-tree-sha1 = "e382260f6482c27b5062eba923e36fde2f5ab0b9"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.0.0"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "6a2f7d70512d205ca8c7ee31bfa9f142fe74310c"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.12"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "38d88503f695eb0301479bc9b0d4320b378bafe5"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.8.2"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "4f6ec5d99a28e1a749559ef7dd518663c5eca3d5"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.3"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "8d7530a38dbd2c397be7ddd01a424e4f411dcc41"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.2"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8977b17906b0a1cc74ab2e3a05faa16cf08a8291"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.16"

[[StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "5783b877201a82fc0014cbf381e7e6eb130473a4"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.0.1"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "57617b34fa34f91d536eb265df67c2d4519b8b98"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.5"

[[SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

[[Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╠═284c4c16-bcf1-11ec-36d9-f724bda4cf9b
# ╠═dc062df1-f6e5-4046-977e-bdb09258a3c6
# ╟─4f044a8d-102d-4e79-9c80-fadab002f328
# ╟─ca41784c-dd23-4e57-ae2f-4ff7593f6966
# ╟─6f83a1ef-3e9a-4281-8ccc-2b0ac0dfec15
# ╟─8b62102b-9c47-46f7-bcdb-871bfc9bd84c
# ╟─11d19f76-fe57-4785-8262-563f5efc7547
# ╠═62a59cec-4f54-491c-a743-e096f4433e95
# ╟─a1c65637-1df2-47a0-80a6-480765670dc0
# ╟─cc9ad52e-b5b0-42cf-b084-8bfbab00d643
# ╟─841d32fd-7b33-453e-a425-2e50d4ab15f6
# ╟─011320ee-9010-4dfb-b13e-ad92e5e40e26
# ╠═2958dbaf-e162-4f49-8d83-5a2ef0d7c824
# ╟─6fdeae6a-1e1e-4301-baa7-91e47b072bcd
# ╟─c7571899-7965-44d6-b8e6-1990d2ae08be
# ╟─2572a88f-53e7-4958-b627-da52cd0b8145
# ╟─686c81e2-44af-48e2-9491-a4f770c336b8
# ╟─fc1b104d-74f0-43ab-b009-8f948c24db08
# ╟─f5d88364-78d0-44f0-bd6c-14328cdb804a
# ╟─8566ff2a-bf3a-4dcd-be5e-568bfccd973f
# ╠═0d0a0bcb-0853-4ad0-9ac3-27578252d3d4
# ╟─dd981f2d-687e-4fb4-8d9d-40d7dd8c0f15
# ╟─e49dbddf-2093-4014-b287-daec3e42216b
# ╟─50fe2f05-abbc-44ae-a4ae-6db3e38def4d
# ╟─58707353-5fc3-4359-b5d2-11951680391b
# ╟─473b097b-6d2d-451b-b6ef-90dddbca3095
# ╠═f89f30dd-7a04-4b6d-b0ae-a73edf7a80cf
# ╟─855eb815-628c-4d8e-9993-e12e3d69f96b
# ╠═80d92a5d-9969-4a5f-b4e3-c29e32dd4cd8
# ╠═72ef239b-49fd-4257-bdf4-103aa178e854
# ╠═21a8107f-90b6-4334-a0f2-88a457e45970
# ╠═8270e08e-e70c-40e7-90e9-35199867a695
# ╠═81784085-0f4c-45aa-9d5c-6eb0a38e2893
# ╠═466d6f8b-3ab9-4bf2-90ce-6c904a209bb5
# ╠═37367483-c9ce-44c1-a668-03376ab77bbd
# ╠═0dd0c2c3-9261-463d-93c6-ddcffdcd51a9
# ╠═08a8fa0b-42fe-4a9e-ab7c-8d6d808a5d2c
# ╟─50931fc3-7bae-4b46-8711-13ba4bf2a864
# ╠═c43f86cf-899b-4538-8685-b36a7e0eabe5
# ╟─1feca48b-c134-4a1e-9a29-af2adf60bd3d
# ╟─88abbec4-8cc2-46ae-870c-ccea1a579e2c
# ╟─62246fc0-9f7c-4c91-98c6-f73a246bf9d6
# ╟─84967d9c-3b16-4eb1-848a-fbb828d79fb8
# ╟─afd92507-eeba-4f0f-b6fc-034fd46edf6f
# ╠═8dc1936e-7ab7-418e-ba00-e31a6308ea3b
# ╟─91035665-c615-4374-ad3c-e09f9a241e52
# ╟─91f11d58-acb1-432c-a0e7-c64124044bbf
# ╟─a791f73d-4caf-4cda-828b-c1867e1e4d3b
# ╟─a8cc023a-cc71-4d49-b715-4d455828cf45
# ╠═b37a1ea6-1760-4795-bd32-c99aa3c1a2db
# ╠═5e2124e2-c59f-4ec6-a1fc-17211b30c759
# ╠═4de61246-6ca7-4f12-8d15-cd7630826a3e
# ╠═bd239232-4c9f-41a1-8cb7-e738f04e98cb
# ╟─69a41c27-6f35-4f77-9035-459a87fb1f82
# ╟─8f9ff7b7-049a-4c50-8afa-0c5c4ddde9c4
# ╟─aec7dffe-1605-4aef-9a12-4f6a633f2bcf
# ╟─b84c0c50-06b0-4955-9c8f-ff6db05f1177
# ╟─2e7deb04-48fb-4925-a89f-aadcb9954118
# ╟─e0ef219e-c837-4257-9dc7-73ed225306f6
# ╟─9b02f8d7-93c4-4ecb-a549-5fd70daa16bb
# ╟─dbf98495-6992-4b7c-a030-694c88c1a7fb
# ╟─6f4774dd-c729-427c-95de-1c2ce15425c2
# ╟─c0072be7-60bc-4764-9026-51a7d46c64b1
# ╟─68c159d6-75b8-4c9b-99c5-0206fe958a15
# ╟─824763d1-9ea4-465e-8ebc-7345c5ffe579
# ╟─96487e3c-b8bb-492c-8628-810026b6bb66
# ╟─bd084197-8c32-4380-b403-7314a603ad09
# ╟─07332e59-4298-4dd4-86c0-7a4959b31b4d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
