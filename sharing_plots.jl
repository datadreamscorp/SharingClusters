### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# ╔═╡ cb223404-03aa-11ee-0e4f-037151aa72ac
begin
	using Pkg
	Pkg.activate(".")
	using Revise
end

# ╔═╡ 860a03ef-32ff-4563-9e25-380359069f49
begin
	using Agents, Distributions, Random, Statistics
	using LambertW, Roots
	using Plots, DataFrames, CSV, ProgressLogging
end

# ╔═╡ f9c52e88-ef18-4bb5-976c-0eb909155570
import SharingClusters as sh

# ╔═╡ d6e849ab-77ff-45d6-8c50-b0117f95968e
md"
# The evolution of risk-pooling with reciprocity norms and networks under costly connections
"

# ╔═╡ dcb1a731-4658-4699-bedf-6174ca8c76dd
md"
### Alejandro Pérez Velilla & Paul Smaldino
"

# ╔═╡ e7cc19cb-14c9-48be-bb5c-a6399f555320
begin
	u = 0.5
	B = 10
	C = 0.1
	β = 1.0
	s = 0.75
end

# ╔═╡ c454acc0-4edf-4529-a655-1be1b6084346
md"
### Figure 1 - numeric simulations vs. mathematically-approximated growth rates
"

# ╔═╡ 1a86569a-f14f-4147-a489-cd5d4f24bd9f
begin
	ylim=[1.10, 1.6]
	max_kbar = 35
	x_axis = 1:(max_kbar)
	x_axis2 = 1:(max_kbar-10)
	
	mean_vec, std_vec = sh.mean_degree_fitness(u, 0.75, B, C, β)
	mean_vec2, std_vec2 = sh.mean_degree_fitness(u, 0.5, B, C, β)
	
	k_plot_sim = plot(x_axis2, mean_vec[x_axis2],
		xlabel="mean degree (k̄), simulated", 
		ylabel="geometric mean fitness (log scale)",
		label="s = 0.9",
		legend=false,
		lw=1.5,
		color=RGBA(0,0.7,1,0.7),
		size=(550,450),
		ylim=ylim
	)
	plot!(x_axis2, mean_vec2[x_axis2],
		label="s = 0.5",
		lw=1.5,
		color=RGBA(0,0.1,0.1,0.9)
	)
	hline!([u*log(B)],
		label="loner",
		lw=1.5,
		color=RGBA(0.7,0,0,0.6)
	)

	k_plot_math = plot(x_axis2, sh.W_s075,
		xlabel="approximation", 
		#ylabel="geometric mean fitness (log scale)",
		label="s = 0.75",
		legend=:topright,
		foreground_color_legend = nothing,
		lw=1.5,
		color=RGBA(0,0.6,1,0.7),
		size=(550,450),
		ylim=ylim
	)
	plot!(x_axis2, sh.W_s05,
		label="s = 0.5",
		lw=1.5,
		color=RGBA(0,0.1,0.1,0.9)
	)
	hline!([u*log(B)],
		label="loner",
		lw=1.5,
		color=RGBA(0.7,0,0,0.6)
	)
	annotate!(23, 1.4, text("B = 10\nC = 0.1\nu = 0.5\nβ = 1", :black, :right, 8))
	
	paired_plots = plot(k_plot_sim, k_plot_math, size=(700,400), link=:all)
end

# ╔═╡ 0c426787-3f6f-4b43-8140-f45752a49fe4
savefig(paired_plots, "./images/figure1.png")

# ╔═╡ 875b643e-b5e7-4497-9b10-5945ae490fd1
md"
### Figure 2 - approximated optimal sharing norms
"

# ╔═╡ 9e405011-abe4-4dd6-8ec7-bbb99e13cd3f
begin
	Zb(k) = [clamp(sh.optshare(y, k, β=1, u=0.2), 0, 1) for y in 1:15]
		
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

# ╔═╡ 2e28d384-2535-472f-bace-2d32f14aca81
savefig(zb_plot, "./images/figure2.png")

# ╔═╡ b694f876-d1f5-4977-b327-938ba31df886
md"
### Figure 3 - approximated optimal sharing norms (full view)
"

# ╔═╡ 3cb7426e-cb2b-472a-a598-edbfa58bdb63
begin
	X = 1:0.01:20
	Y = 1:0.01:20

	Z1 = [clamp(sh.optshare(x, y, β=0.75, u=0.8), 0, 1) for x in X, y in Y]
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

	Z2 = [clamp(sh.optshare(x, y, u=0.8), 0, 1) for x in X, y in Y]
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

	Z3 = [clamp(sh.optshare(x, y, β=1.25, u=0.8), 0, 1) for x in X, y in Y]
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
	
	Z4 = [clamp(sh.optshare(x, y, β=0.75), 0, 1) for x in X, y in Y]
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

	Z5 = [clamp(sh.optshare(x, y), 0, 1) for x in X, y in Y]
	pz5 = plot(
		X, Y, Z5, st=:contourf,
		#xlabel="β = 1",
		#ylabel="k̄",
		xaxis=nothing,
		yaxis=nothing,
		clims=(0,1),
		#legend=false
	)
	
	Z6 = [clamp(sh.optshare(x, y, β=1.25), 0, 1) for x in X, y in Y]
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

	Z7 = [clamp(sh.optshare(x, y, β=0.75, u=0.2), 0, 1) for x in X, y in Y]
	pz7 = plot(
		X, Y, Z7, st=:contourf,
		xlabel="β = 0.75",
		ylabel="u = 0.2",
		#yaxis=nothing,
		clims=(0,1),
		legend=false
	)

	Z8 = [clamp(sh.optshare(x, y, u=0.2), 0, 1) for x in X, y in Y]
	pz8 = plot(
		X, Y, Z8, st=:contourf,
		xlabel="β = 1",
		#ylabel="k̄",
		yaxis=nothing,
		clims=(0,1),
		legend=false
	)

	Z9 = [clamp(sh.optshare(x, y, β=1.25, u=0.2), 0, 1) for x in X, y in Y]
	pz9 = plot(
		X, Y, Z9, st=:contourf,
		xlabel="β = 1.25",
		#ylabel="u = 0.2",
		ymirror=true,
		yaxis=nothing,
		clims=(0,1),
		legend=false
	)

	

	Z = [clamp(sh.optshare(x, y, β=0.75), 0, 1) for x in X, y in Y]
	pz = plot(
		X, Y, Z, st=:contourf,
		#xlabel="B",
		#ylabel="k̄",
		clims=(0,1),
	)


	s_opt_contour = plot( 
	pz1,pz2,pz3, 
	pz4,pz5,pz6, 
	pz7,pz8,pz9, 
	layout=(3,3), 
	main="optimal sharing norm (s*)", 
	size=(650,650), 
	linewidth=1, 
	xtickfontsize=9, 
	ytickfontsize=9, 
	xguidefontsize=12, 
	yguidefontsize=12, 
	titlefontsize=12,
	levels=0:0.05:1,
	#legend=true
	)
end

# ╔═╡ eef9cd91-f1dd-40f2-aab6-954754d24f2d
savefig(s_opt_contour, "./images/figure3_unlabeled.png")

# ╔═╡ 005a77d0-4b11-418f-9e18-5aaaaa061632
md"
### Figure 4 - approximated optimal mean degree 
"

# ╔═╡ 1ca4b6c2-68c8-492f-88d4-e334af61d170
begin
	
	x = 0.01:0.001:1
	f(x) = log( -0.25/( log(1-x) ) ) / ( log(1-x) )

	bc_ratios = sort( collect(50:20:500), rev=true )
	
	cb_ratios = sort( collect(0.002:0.005:0.1), rev=true )
	y_step = 0.1
	
	kopt_plot = plot( x, f.(x), 
		ylim=(0,18), 
		legend=false, 
		color=RGBA(0,0.1,0.1,0.0),
		xlab="probability of success (u)",
		ylab="optimal mean degree (k̂) at β = 1"
	)

	count = 1
	for cb in bc_ratios
		yh = 13 - count*y_step
		#g(x) = log( -cb/( x * log(1-x) ) ) / ( log(1-x) )
		g(x) = log( -(1/cb)/( x * log(1-x) ) ) / ( log(1-x) )
		col = RGBA(
			1 - findall(x->x==cb, bc_ratios)[1]*(1/length(bc_ratios)),
			0.1,
			0.01,
			0.9
		)
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
		(0.78, 10.4, text("50", 6)),
		(0.785, 13.3, text("500", 6)),
		(0.75, 15.2, text("expected shared benefit /\nconnection cost\nratio", 8))
	]
	)
	
	kopt_plot
end

# ╔═╡ 13b589ad-40bd-4200-b0fb-e425bfeac925
savefig(kopt_plot, "./images/figure4.png")

# ╔═╡ c5db3f92-b233-4f26-8b47-2af8d56ac132
md"
### Figure 5 - approximated maximum mean degree
"

# ╔═╡ 24156d86-873e-476e-9b87-55664bbfcfa8
begin
	kbar_bound, plog, kbar_no_bound = sh.maxshare_Du(B=5, C=0.05)
	full_plog = plog ./ log.( 1 .- collect(0.001:0.01:1) )
	
	kbarcito_plot = plot(
		0.001:0.01:1, kbar_bound, 
		lw=3, alpha=0.5,
		ylab = "maximum mean degree (k̄*)",
		xlab = "probability of success (u)",
		label = "with boundaries"
	)
	plot!(
		0.001:0.01:1, kbar_no_bound,
		lw=3, alpha=0.5,
		label = "no boundaries"
	)
	annotate!(0.85, 11, text("B = 5\nC = 0.05\nβ = 1", 8))
end

# ╔═╡ 6e1100d8-3fd8-4453-9039-77054526452b
savefig(kbarcito_plot, "./images/figure5.png")

# ╔═╡ 7b18239f-a61a-4f48-80fb-3bb502a9a55b
begin
	"""
	begin
		khat_b5, kstar_b5 = sh.calculate_points_Sim(s, 5, 0.05, 1, w_0=0.1, max_kbar=35, reps=1000)
		khat_b15, kstar_b15 = sh.calculate_points_Sim(s, 15, 0.15, 1, w_0=0.1, max_kbar=35, reps=1000)
		khat_b30, kstar_b30 = sh.calculate_points_Sim(s, 30, 0.3, 1, w_0=0.1, max_kbar=35, reps=1000)
	
		points = [
				(khat_b5, kstar_b5),
				(khat_b15, kstar_b15),
				(khat_b30, kstar_b30),
			]
	
		Bs = [5, 15, 30]
		
		df_b = [
			DataFrame(
				khat = points[i][1],
				kstar = points[i][2],
				u = 0:0.01:1,
				B = Bs[i],
				C = Bs[i]/100,
				s = s
			)
			for i in 1:length(points)
			]
	
		df_b = vcat(df_b...)
	
		CSV.write("./data/fig6_data.csv", df_b)
	end
	"""
	
md"
##### The code in this cell can be uncommented and run to simulate the data for figure 6 and save it as a dataframe in the data folder.
	
> Be advised that these simulations might take **a while** to run to completion.
"
end

# ╔═╡ f4ea2fd8-0196-47e6-bd4e-6f155f2a55d0
md"
### Figure 6 - simulations vs. approximations for mean degree
"

# ╔═╡ 202fab42-ddf1-4f7d-bb43-ff4ea9b7b9e0
begin
	b_df_saved = DataFrame(CSV.File("./data/fig6_data.csv"))
	b_df5 = b_df_saved[b_df_saved.B .== 5, :]
	b_df15 = b_df_saved[b_df_saved.B .== 15, :]
	b_df30 = b_df_saved[b_df_saved.B .== 30, :]
	
	mksize = 3
	ylimit = 35
	
	zeros1 = sh.zero_plots( 
		sh.calculate_points_Du(B=5, C=0.05, s=s), 
		ylim=ylimit, 
		lab=false 
	)
	scatter!(
		b_df5.u, 
		b_df5.khat, 
		alpha=0.3, 
		markersize=mksize, 
		label="", 
		color=RGBA(0.5,0.5,0.9,0.8)
	)
	plot!(
		1:0, 
		label="Maximum mean degree (k̄*) \n with δ = 0.2 + log(B - 1)",
		#line=:dash,
		color=RGBA(0.5,0.5,0.9,0.8),
		lw=2.5
	)
	scatter!(
		b_df5.u, 
		b_df5.kstar, 
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
	annotate!(0.25, 23, text("s = $s\nβ = 1", 8))
	
	
	
	zeros2 = sh.zero_plots( 
		sh.calculate_points_Du(B=15, C=0.15, s=s),  
		lab=false, 
		leg=false,
		ylim=ylimit,
		yax=nothing
	)
	scatter!(
		b_df15.u, 
		b_df15.khat, 
		alpha=0.3, 
		markersize=mksize, 
		label="k̄*", 
		color=RGBA(0.5,0.5,0.9,0.8)
	)
	scatter!(
		b_df15.u, 
		b_df15.kstar, 
		alpha=0.3, 
		markersize=mksize, 
		label="k̃", 
		color=RGBA(0.9,0.5,0.5,0.9)
	)


	
	zeros3 = sh.zero_plots( 
		sh.calculate_points_Du(B=30, C=0.3, s=s), 
		lab=false, 
		leg=false,
		ylim=ylimit,
		yax=nothing
	)
	scatter!(
		b_df30.u, 
		b_df30.khat, 
		alpha=0.3, 
		markersize=mksize, 
		label="k̄*", 
		color=RGBA(0.5,0.5,0.9,0.8)
	)
	scatter!(
		b_df30.u, 
		b_df30.kstar, 
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

# ╔═╡ 4b82fe05-dbff-4ff5-9f8f-41a96df688c9
savefig(zeros_plot, "./images/figure6_unlabeled.png")

# ╔═╡ d6e31887-51e2-4a9c-92e5-05143a1871a3
md"
### Agent-based model of fission-fusion cluster dynamics
"

# ╔═╡ ced5cd17-939b-456b-a990-4eebdab1f5c5
model = sh.sharing_groups_model(;
	max_N=200,
	n=200,
	B=4.0,
	bc_ratio=100.0,
	u=u,
	T=100,
	w0=0.5,
	σ_small=0.01,
	δ = 0.01,
	size_bias=true,
	seed=28814091
)

# ╔═╡ 1be4fb02-a98f-4500-b9cc-230ecde4e061
begin
	s_hist_init = histogram( 
			[a.s for a in allagents(model)],
			bins = 0:0.01:1, 
			ylabel = "t = 0", 
			legend = false 
		)
	cluster_hist_init = histogram( 
			[a.size for a in allagents(model)],
			xlabel = "", 
			legend = false,
			color = "dark red"  
		)
	l1 = @layout [
				grid(2,2)
		     	]
end

# ╔═╡ 6b10a2ec-d006-4e47-9098-8f72eed0795e
step!(model, dummystep, sh.sharing_step!, 5000)

# ╔═╡ 0551d28b-b569-4983-b564-37adb611ac6c
md"
### Figure 7 - Prior vs. posterior distributions for sharing norms and cluster sizes
"

# ╔═╡ a2f3a45d-36d5-4add-b6db-60497f04de59
begin
	sizes = [a.size for a in allagents(model)]
	
	sharings = [a.s for a in allagents(model)]
	s_hist_post = histogram( 
		sharings, 
		bins = 0:0.01:1, 
		xlabel="sharing norm",
		ylabel="t = 5000",
		legend=false 
	)
	annotate!(
		0.25, 
		5, 
		text("B = $(model.B) \nC = $(model.C) \nu = $(model.u) \nβ = 1", 8)
	)
	
	cluster_hist_post = histogram( 
		sizes, 
		bins=0:1:maximum(sizes)+1, 
		xlabel="cluster size", 
		legend=false,
		color="dark red"
	)

	s_hist_complete = plot(
		s_hist_init,
		s_hist_post,
		layout=(2,1),
		link=:all,
	)

	cluster_hist_complete  = plot(
		cluster_hist_init,
		cluster_hist_post,
		layout=(2,1),
		link=:all,
	)
	
	abm_hist_grid = plot(
			s_hist_complete, 
			cluster_hist_complete,
			layout=(1,2)
		)
end

# ╔═╡ 363296c2-4a56-46d2-959d-39aff603e1be
savefig(abm_hist_grid, "./images/figure7.png")

# ╔═╡ 4efa0a3a-b935-4787-9186-4c0cb42dc060
md"
### Figure 8 - Time trajectories for sharing norms and cluster sizes
"

# ╔═╡ adb1cbd7-87c4-4755-919a-f3905cae2c5b
abm_time_plots = sh.abm_plots(model)

# ╔═╡ 7673ee2a-59fd-4aa3-937e-064e4cdbbd9c
savefig(abm_time_plots, "./images/figure8.png")

# ╔═╡ f7e2e305-d635-430c-b88b-059e9ba43ca7
md"
### Figure 9 - Expected cluster size estimate from approximation vs. ABM simulations
"

# ╔═╡ d71b0faf-0a78-47c1-9908-d7db9c1c5066
begin
	df = DataFrame( CSV.File("./data/sharing_clusters_data_Final.csv") )

	b = [6.0, 20.0, 40.0]
	kbar_sweep_plots = []
	idx = repeat(1:length(b), 2)
	counter = 1
	
	for i in idx
		df2 = df[df.step .== 5000, :]
		df2 = df2[df2.B .== b[i], :]
		df2 = counter <= length(idx)/2 ? df2[df2.δ .== 0.01, :] : df2[df2.δ .== exp(1)-1, :]
		df2 = groupby(df2, :u)
		df2 = combine(df2, valuecols(df2) .=> mean )

		if counter == 1
			y_label = "ζ = 0.01"
		elseif counter == length(idx)/2 + 1
			y_label = "ζ = e - 1"
		else
			y_label = ""
		end
		
		pi = plot(
			df2.u, 
			df2.mean_cluster_size_mean,
			legend=:top,
			label = counter == length(idx)/2 + 1 ? "mean" : "",
			lw=2,
			ylab=y_label
		)
		plot!(
			df2.u, 
			df2.median_cluster_size_mean,
			legend=:top,
			label = counter == length(idx)/2 + 1 ? "median" : "",
			lw=2,
			ylab=y_label
		)
		plot!(
			0:0.01:1.0, 
			1 .+ ( 0.5 .* sh.kmax.(0:0.01:1.0, B=b[i], C=b[i]/100) ),
			label = counter == length(idx)/2 + 1 ? "(k̄*/2) + 1" : "",
			title = counter <= length(idx)/2 ? "B = $(b[i])" : "",
			titlefontsize=9
		)
		#annotate!()
		
		push!(kbar_sweep_plots, pi)

		global counter += 1
	end

	counter = 1
	
	kbar_sweep_plot = plot(
		kbar_sweep_plots..., 
		layout=(2,length(b)),
		size=(600,400),
		link=:all
	)
	
end

# ╔═╡ afb1034b-db38-44fc-bfcd-5934afc6f51c
savefig(kbar_sweep_plot, "./images/figure9_unlabeled.png")

# ╔═╡ 3b1b541b-8027-449d-9177-b452d2c181da
md"
### Figure 10 - Optimal sharing norm estimate ($s^*$) vs. ABM simulations
"

# ╔═╡ aaf7d958-225f-43f2-b7ad-15fca6c58bbe
begin
	U =[0.2, 0.5, 0.8]

	sstar_sweep_plots = []
	idx2 = repeat(1:length(U), 2)
	counter2 = 1
	
	for i in idx2
	
		df3 = df[df.step .== 5000, :]
		df3 = df3[df3.u .== U[i], :]
		df3 = counter2 <= length(idx2)/2 ? df3[df3.δ .== 0.01, :] : df3[df3.δ .== exp(1)-1, :]
		df3 = groupby(df3, :B)
		df3 = combine(df3, valuecols(df3) .=> mean )
	
		opt_curve = clamp.( 
			sh.optshare2.(
				df3.B, 
				df3.median_cluster_size_mean,
				U[i],
				df3.C_mean), 
			0, 
			1
		)

		if counter2 == 1
			y_label = "ζ = 0.01"
		elseif counter2 == length(idx)/2 + 1
			y_label = "ζ = e - 1"
		else
			y_label = ""
		end
		
		pi = plot(
			df3.B,
			df3.mean_sharing_mean,
			lw=2,
			label = counter2 == length(idx)/2 ? "mean" : "",
			#xlab="surplus (B)",
			ylab=y_label,
			legend=:bottomright
		)
		plot!(
			df3.B, 
			df3.median_sharing_mean,
			lw=2,
			label = counter2 == length(idx)/2 ? "median" : ""
		)
		plot!(
			df3.B, 
			opt_curve,  
			lw=1.5,
			label = counter2 == length(idx)/2 ? "s*" : "",
			title = counter2 <= length(idx)/2 ? "u = $(U[i])" : "",
			titlefontsize=9
		)

		push!(sstar_sweep_plots, pi)
		global counter2 += 1
	end

	sstar_sweep_plot = plot(
		sstar_sweep_plots..., 
		layout=(2,length(b)),
		size=(600,400),
		link=:all
	)
end

# ╔═╡ 5343ae11-baa6-4c6c-a84a-1ab1be027958
savefig(sstar_sweep_plot, "./images/figure10_unlabeled.png")

# ╔═╡ 13da6b96-b1b0-426e-a4ec-0203f32834c3
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	fitnesses = [log(a.average_fitness)/model.T for a in allagents(model)]
	histogram( 
		fitnesses, 
		#bins=0:1:maximum(fitnesses)+1, 
		xlabel="average growth rate", 
		legend=false,
		color="dark green"
	)
	vline!([log(model.loner_fitness)/model.T])
	"deprecated"
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─cb223404-03aa-11ee-0e4f-037151aa72ac
# ╠═f9c52e88-ef18-4bb5-976c-0eb909155570
# ╟─860a03ef-32ff-4563-9e25-380359069f49
# ╟─d6e849ab-77ff-45d6-8c50-b0117f95968e
# ╟─dcb1a731-4658-4699-bedf-6174ca8c76dd
# ╠═e7cc19cb-14c9-48be-bb5c-a6399f555320
# ╟─c454acc0-4edf-4529-a655-1be1b6084346
# ╟─1a86569a-f14f-4147-a489-cd5d4f24bd9f
# ╟─0c426787-3f6f-4b43-8140-f45752a49fe4
# ╟─875b643e-b5e7-4497-9b10-5945ae490fd1
# ╟─9e405011-abe4-4dd6-8ec7-bbb99e13cd3f
# ╟─2e28d384-2535-472f-bace-2d32f14aca81
# ╟─b694f876-d1f5-4977-b327-938ba31df886
# ╟─3cb7426e-cb2b-472a-a598-edbfa58bdb63
# ╟─eef9cd91-f1dd-40f2-aab6-954754d24f2d
# ╟─005a77d0-4b11-418f-9e18-5aaaaa061632
# ╟─1ca4b6c2-68c8-492f-88d4-e334af61d170
# ╟─13b589ad-40bd-4200-b0fb-e425bfeac925
# ╟─c5db3f92-b233-4f26-8b47-2af8d56ac132
# ╟─24156d86-873e-476e-9b87-55664bbfcfa8
# ╟─6e1100d8-3fd8-4453-9039-77054526452b
# ╟─7b18239f-a61a-4f48-80fb-3bb502a9a55b
# ╟─f4ea2fd8-0196-47e6-bd4e-6f155f2a55d0
# ╟─202fab42-ddf1-4f7d-bb43-ff4ea9b7b9e0
# ╟─4b82fe05-dbff-4ff5-9f8f-41a96df688c9
# ╟─d6e31887-51e2-4a9c-92e5-05143a1871a3
# ╠═ced5cd17-939b-456b-a990-4eebdab1f5c5
# ╟─1be4fb02-a98f-4500-b9cc-230ecde4e061
# ╠═6b10a2ec-d006-4e47-9098-8f72eed0795e
# ╟─0551d28b-b569-4983-b564-37adb611ac6c
# ╟─a2f3a45d-36d5-4add-b6db-60497f04de59
# ╟─363296c2-4a56-46d2-959d-39aff603e1be
# ╟─4efa0a3a-b935-4787-9186-4c0cb42dc060
# ╟─adb1cbd7-87c4-4755-919a-f3905cae2c5b
# ╟─7673ee2a-59fd-4aa3-937e-064e4cdbbd9c
# ╟─f7e2e305-d635-430c-b88b-059e9ba43ca7
# ╟─d71b0faf-0a78-47c1-9908-d7db9c1c5066
# ╟─afb1034b-db38-44fc-bfcd-5934afc6f51c
# ╟─3b1b541b-8027-449d-9177-b452d2c181da
# ╟─aaf7d958-225f-43f2-b7ad-15fca6c58bbe
# ╟─5343ae11-baa6-4c6c-a84a-1ab1be027958
# ╟─13da6b96-b1b0-426e-a4ec-0203f32834c3
