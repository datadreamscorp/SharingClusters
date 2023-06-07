### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# ╔═╡ cb223404-03aa-11ee-0e4f-037151aa72ac
begin
	using Pkg
	Pkg.activate("..")
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
#### Year 2023 AD
"

# ╔═╡ e7cc19cb-14c9-48be-bb5c-a6399f555320
md"
We start with a population of individuals who procure resources during their lifetimes. Lifetimes are separated into discrete steps, each representing a resource acquisition event (i.e. a foraging trip), with an initial resource level $w_0$ and individual payoff representing productivity across the lifetime. At every step, an individual attempts to get resource, sometimes getting a base rate of resource $b$ with probability $1-u$ (risk) and a higher, surplus rate $B$ with probability $u$ (success rate). The payoff's growth is governed by multiplicative effects: every payoff gained at a particular step is multiplied by the payoff accumulated until that point. This assumes the payoff-accumulation process is not memory-less: my increase in payoff at a particular moment is not independent of the value of my payoffs up until that moment. This can be contrasted with an additive process, in which the increase in payoff at a time period is completely independent of the agent's payoff. We can set $b = 1$ for convenience, which assumes that, as far as the modeled individuals are concerned, the baseline is the same as getting nothing of value. The risk $1-u$ then represents the risk of getting nothing out of your productive effort.

There are two strategies that we care about. The first is the \textbf{loner}, who individually procures resource without interacting with peers. The second one is the \textbf{sharer} strategy, which characterizes individuals that are part of a group employing a sharing norm. This means that individuals within the group agree on a proportion $s$ of surplus $B$ that should be shared with other group individuals that they have a reciprocity connection with. Essentially, whenever a sharer is successful in obtaining surplus $B$, the break off $s$ of it and divide the resulting piece to give away to their connections within the group. In doing so, they can decrease the variance of their per-time-period resource acquisition rates. Furthermore, reciprocity connections are assumed to be costly, with sharers taking a per-connection loss $C$ in units of resource in order to maintain each reciprocity connection.

Sharers establish connections with others, so they can be characterized by a number of connections $k_i$ (where $i$ is the sharer's index). A sharer's payoff will depend not only on their number of connections, but also on the number of connections that each of their connections has. This is because each one of a sharer's network peers is dividing their resource among their own peers. A sharer that is connected to individuals who have a higher number of connections can be seen as investing more into their network than what they are obtaining from it, while a an individual who is connected to others who have less connections enjoys a privileged position in the network. Therefore, an individual's position within their sharing network depends on how their number of connections $k_i$ compares to the average number of connections of their network peers ($\bar{k}_i$). In fully regular networks, like lattices and cliques, $k_i$ will equal $\bar{k}_i$, and sharing network hierarchy is not a thing. We call the parameter $\frac{k_i}{\bar{k}_i} = \beta_i$ an individual's \textit{relative degree}, which tells us about their network position.

We simulate the production and sharing dynamic described above and construct mathematical approximations of it in order to get useful equations to describe it (see \textbf{Figure 1}). The simulations show us that:

- The payoff of sharers exceeds that of loners already in small networks. This resonates with Winterhalder's classic risk-sharing model, where the increase in payoff due to variance reduction diminishes with every additional connection.

- The payoff of sharers peaks ($\hat{k}$) with respect to their ego network's mean degree (fixing all other parameters). This means there is an optimal network density (network size, if the network is a clique). After this point, additional connection costs overtake the diminishing benefits of payoff variance reduction.

- After the optimal point, sharers in networks can still enjoy higher payoffs than loners until they reach an overcrowding limit ($\bar{k}^*$), the point at which the network is so dense/large that the payoff of sharers intersects the payoff of loners. This is the maximum degree of an agent's sharing network.

- Optimal sharing norms increase with higher surpluses, denser/larger networks and better position within the network hierarchy, while effect of risk is mediated by position within the hierarchy. If networks are large enough and networks are regular (so that hierarchy effects are muted), optimal sharing norms tend to full sharing ($s = 1$) as surplus size increases towards infinity.
"

# ╔═╡ c454acc0-4edf-4529-a655-1be1b6084346
md"
### Figure 1 - numeric simulations vs. approximated growth rates
"

# ╔═╡ 1a86569a-f14f-4147-a489-cd5d4f24bd9f
begin

	u = 0.5; B = 10; C = 0.1; β = 1.0; s = 0.75
	
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
		lw=2,
		color=RGBA(0,0.7,1,0.7),
		size=(550,450),
		ylim=ylim
	)
	plot!(x_axis2, mean_vec2[x_axis2],
		label="s = 0.5",
		lw=2,
		color=RGBA(0,0.1,0.1,0.9)
	)
	hline!([u*log(B)],
		label="loner",
		lw=2,
		color=RGBA(0.7,0,0,0.6)
	)

	k_plot_math = plot(x_axis2, sh.W_s075,
		xlabel="approximation", 
		#ylabel="geometric mean fitness (log scale)",
		label="s = 0.75",
		legend=:topright,
		foreground_color_legend = nothing,
		lw=2,
		color=RGBA(0,0.6,1,0.7),
		size=(550,450),
		ylim=ylim
	)
	plot!(x_axis2, sh.W_s05,
		label="s = 0.5",
		lw=2,
		color=RGBA(0,0.1,0.1,0.9)
	)
	hline!([u*log(B)],
		label="loner",
		lw=2,
		color=RGBA(0.7,0,0,0.6)
	)
	annotate!(23, 1.4, text("B = 10\nC = 0.1\nu = 0.5\nβ = 1", :black, :right, 8))
	
	paired_plots = plot(k_plot_sim, k_plot_math, size=(700,400), link=:all)
end

# ╔═╡ 0c426787-3f6f-4b43-8140-f45752a49fe4
begin
	#savefig(paired_plots, "../images/figure1.png")
md"
I propose using a **long-term geometric growth rate approximation** (essentially a geometric mean fitness, where sharing is a form of diversified bet-hedging), with a correction term penalizing the higher risk that smaller networks are exposed to (the smaller a network is, the higher the probability that no network peers are successful during a time period, yielding nothing but costs for the focal agent, and in a multiplicative dynamic the cases of heavy losses are very important).

To construct the approximation, we can consider the different components of a sharer's payoff. When successful during a time period, a sharer obtains


- A **piece of surplus** $(1 - s)B$ which is not shared among peers. The rest of this $sB$ goes to the sharer's peers, so it is not present in the payoff.

- The **shared surplus** from successful peers. Every peer of $i$ gets surplus $B$ with probability $u$, so on average there are $u k_i$ successful peers per time period. Peers, on average, have $\bar{k}_i$ connections, and they share $sB$ of their surplus among them. The average shared surplus received by a sharer $i$ is then $u k_i \frac{sB}{\bar{k_i}} = u \beta_i sB$. I include a correction factor $\rho_i$ into this average, meant to account for cases like correlated failures (see below). The whole thing is then $\rho_i u \beta_i sB$.

- The **cumulative network costs** of connection maintenance. If every connection costs $C$ to maintain, then the total cost of network maintenance can be written as $k_i C = \beta_i \bar{k}_i C$.

When a sharer fails, they get only get the shared surplus and pay the costs (remember, we are assuming here $b = 1$, which means the baseline is nothing). Given this, we can construct the approximate payoff. Or rather, its long-term average geometric growth rate, which is what we really care about:

$G_i (\textrm{Sharer}) = u \log{ \left\{ (1 - s)B + S_i \right\} } + (1 - u) \log{ \left\{ 1 + S_i \right\} }$

where the term

$S_i = \beta_i \left( \rho_i u s B - \bar{k}_i C \right)$

captures the network sharing components of the payoff. Here I set $\rho_i = 1 - (1-u)^{k_i}$ as the above-mentioned correction factor, which is nothing more than the probability that at least one of the focal agent's network peers is successful and has something to share. As $i$'s number of connections grows, this probability tends to 1. But having too little connections can put you in the situation of having a high chance of only experiencing costs in every time period, which is worse than just getting the baseline (worse than nothing, in this case). This happens because the chance of all of your peers failing is higher the less peers you have. In a multiplicative scenario like this one, this case is important, and this factor allows us to capture a bit of this importance, at least qualitatively. Here we assume that individuals' risks are uncorrelated, but we can also put the effect of risk correlation into this correction factor as well if we want to look at that.

By the same logic, the (geometric growth rate of the) payoff of a loner can be written as

$G(\textrm{Loner}) = u \log{B} + (1-u) \log{1} = u \log{B}$
"
end

# ╔═╡ 342da70f-d2e4-43c5-b8b9-f97a24ca2d23
html"<hr>"

# ╔═╡ 875b643e-b5e7-4497-9b10-5945ae490fd1
md"
### Figure 2 - approximated optimal sharing norms
"

# ╔═╡ 9e405011-abe4-4dd6-8ec7-bbb99e13cd3f
begin
	Zb(k) = [clamp(sh.optshare(y, k, β=1, u=0.2, c=k/100), 0, 1) for y in 1:15]
		
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

	annotate!(8, 0.2, text("B/C = 100 \nu = 0.2 \nβ = 1", 8))
end

# ╔═╡ 2e28d384-2535-472f-bace-2d32f14aca81
begin
	#savefig(zb_plot, "../images/figure2.png")
md"
Optimal sharing norms from the approximate growth rate are given by

$s^*_i = \frac{
   \frac{ (1 - u) (B - \beta_i \bar{k} C) }{ \beta_i u \left[(1-u)^{\beta_i  k}-1\right] + 1 } + 
   \frac{ 1-\beta_i \bar{k} C }{ \beta_i \left[(1-u)^{\beta_i k}-1\right] }
}{B}$

with the limiting cases, first for large-enough networks

$\lim_{\rho \rightarrow 1} s^*_i = \frac{\beta_i \left\{(\beta_i - 1) \bar{k} C - (1-u)B - u\right\} + 1}{ \beta_i  (\beta_i  u - 1) B}$

and, additionally, regular networks

$\lim_{\rho \rightarrow 1} s^*_i \left.\right|_{\beta_i=1} = \frac{B-1}{B}$
"
end

# ╔═╡ 5d1c8435-7c7c-445b-a5f8-e6cef400bf70
html"<hr>"

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
begin
	#savefig(s_opt_contour, "../images/figure3_unlabeled.png")
md"
Visualizing the full relationship requires contour plots. On the x-axis we have surplus $B$, and on the y-axis we have network mean degree $\bar{k}$. Here the effects of network inequality (through relative degree $\beta$) and risk (through $u$) can be appreciated. A disadvantage in relative degree (\beta < 1) will lead to lower optimal sharing norms with respect to individuals with higher advantage (\beta \geq 1). The effect of risk changes depending on the level of advantage. For advantaged individuals ($\beta > 1$) lower risk (higher $u$) leads to higher optimal sharing norms. For disadvantaged individuals ($\beta < 1$), higher risk does the opposite. This tells us that in non-regular networks there will be a conflict of interest between advantaged and disadvantaged individuals with respect to which sharing norms to use, which might make it more difficult for such networks to persist. For individuals sitting at their network mean degree (or, equivalently, individuals in regular networks) a change in $u$ merely changes the curvature of the optimal sharing norm function, such that higher risk means surplus $B$ must increase more before plateauing. In low risk environments, $s^*$ plateaus quickly with respect to $B$.
"
end

# ╔═╡ c508a66a-dda0-4be0-835b-8df1fecfa24c
html"<hr>"

# ╔═╡ 005a77d0-4b11-418f-9e18-5aaaaa061632
md"
### Figure 4 - approximated optimal mean degree 
"

# ╔═╡ 1ca4b6c2-68c8-492f-88d4-e334af61d170
begin
	
	x = 0.01:0.001:1
	f(x) = log( -0.25/( log(1-x) ) ) / ( log(1-x) )

	bc_ratios = sort( collect(50:25:500), rev=true )
	
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
		(0.78, 10.7, text("50", 6)),
		(0.785, 13.3, text("500", 6)),
		(0.75, 15.2, text("expected shared benefit /\nconnection cost\nratio (sB/C)", 8))
	]
	)
	
	kopt_plot
end

# ╔═╡ 13b589ad-40bd-4200-b0fb-e425bfeac925
begin
	#savefig(kopt_plot, "../images/figure4.png")
md"
For network density approximations, I propose a **rank-dependent utility model** which can be used to approximate geometric growth dynamics (see Price and Jones, 2020). This leads to payoffs that look like this:

$V_i(\textrm{Sharer}) = u^\delta \left\{ (1 - s)B + S_i \right\} + (1 - u^\delta) S_i$

$V(\textrm{Loner}) = u^\delta B$

We can think of this payoff structure as representing not an objective growth rate, but an individual decision-maker's perception of the payoffs in question. It includes an additional parameter $\delta$ that can be interpreted as an individual's attitude towards uncertainty. When $0 < \delta < 1$, individuals over-weight the success rate of acquiring surplus, leading to optimistic probability weighting. When $\delta = 1$, individuals take a standard arithmetic mean of the possible outcomes with respect to risk and success rates (neutral probability weighting). When $\delta > 1$, individuals increasingly over-weight risks, leading to pessimistic probability weighting. This allows for a tractable approximation of optimal network density, which we arrive by optimizing for $\bar{k}_i$, yielding the following equation:

$\hat{k}_i = \frac{ \log \left(\frac{-1}{\beta_i u \log (1-u)} \right) - \log \left(\frac{s B}{C}\right)}{\beta_i \log (1-u)}$

By fixing $\beta = 1$, we can see how the expect shared benefit to cost ratio ($\frac{s B}{C}$) mediates the effect of risk ($1 - u$). Importantly, at high enough values of $\frac{s B}{C}$, optimal cluster size becomes significantly higher at high risk versus low risk. At both ends of the spectrum, optimal cluster size collapses.
"
end

# ╔═╡ 6060b2c8-87de-4ba7-b34b-c2f2c58305a8
html"<hr>"

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
begin
	#savefig(kbarcito_plot, "../images/figure5.png")
md"
By asking when $V_i(\textrm{Sharer}) > V(\textrm{Loner})$ and solving for $\bar{k_i}$, the rank-dependent model spits out this approximation for the maximum density:

$\bar{k}^*_i = \frac{ u \left(\beta_i - u^{\delta-1}\right)    \frac{ s B }{C} - 
    \frac{W \left( \Sigma_i \right)}{\log(1-u)}
}{\beta_i}$

with

$\Sigma_i = (1-u)^{ u ( \beta_i - u^{\delta-1} ) \frac{s B}{C} } \cdot \beta_i u \log (1-u) \cdot \frac{s B}{C}$

the $\frac{W(\Sigma)}{\log(1-u)}$ term delimits the regions of parameter space for which sharing is possible as a strategy at all ($W$ is the Lambert W function, also known as the product-logarithm). If we assume that we are within the region where sharing is possible, then this term vanishes and we have

$\bar{k}^*_i \approx u \left(\beta_i - u^{\delta-1}\right) \frac{s B}{\beta_i C}$

which, for the case of regular networks, becomes

$\bar{k}^*_i \left.\right|_{\beta_i = 1} \approx u \left(1 - u^{\delta-1}\right) \frac{s B}{C}$

which tells us that we cannot ignore individuals' attitudes towards uncertainty when using this approximation. In particular, individuals who employ pessimistic probability weighting ($\delta > 1$) are the only ones that will choose to be sharers (in a regular network, this may change in an irregular one). As with the optimal density, the ratio of shared benefits to connection costs makes an appearance, but in this case its effects are linear, making the equation easier to read. The effect of risk at the extremes of the spectrum is easier to see as well. The term involving risk looks like the variance of a Bernoulli random variable, but incorporating the effect of probability weighting. But at any level of probability weighting, sharing collapses when risk is too high or too low.
"
end

# ╔═╡ d1e915da-a3f3-4d92-bf0f-80ede0485dc0
html"<hr>"

# ╔═╡ f4ea2fd8-0196-47e6-bd4e-6f155f2a55d0
md"
### Figure 6 - simulations vs. approximations for mean degree
"

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
>The code in this cell can be uncommented and run to simulate the data for figure 6 and save it as a dataframe in the data folder.
> Be advised that these simulations might take **a while** to run to completion.
"
end

# ╔═╡ 202fab42-ddf1-4f7d-bb43-ff4ea9b7b9e0
begin
	b_df_saved = DataFrame(CSV.File("../data/fig6_data.csv"))
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
	annotate!(0.5, 23, text("s = $s\nβ = 1", 8))
	
	
	
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
begin
	#savefig(zeros_plot, "../images/figure6_unlabeled.png")
md"
Optimal cluster size is extracted from simulations (red dots) and estimated using the rank-dependent approximation (red line). The same procedure is done with maximum cluster size (violet dots vs. violet line), and repeated for repeated surplus ($B$) values, mwhile keeping benefit-to-cost ratio constant at $\frac{B}{C} = 100$.

Optimal cluster size does not change in the approximation, and under increasing $B$ undergoes a transition from a lower regime, to a higher regime, where it settles. This is due to the effects of absorbing barriers in the numerical simulations: a sharer can end up losing a proportion of at most $w_0 = 0.5$ of their payoff in the worst possible time period. This generates a valley in the payoff function at higher values of risk, and a dual-peaked fitness function. As $B$ increases, the second peak overtakes the first and the optimal cluster size stabilizes.

Maximum cluster size increases even with a fixed benefit-cost ratio, due to the multiplicative nature of returns. This is well approximated by a probability weighting exponent that depends on the difference between baseline and surplus, $\delta = \alpha + \log(B - 1)$. This implies increasing pessimism (at a diminishing rate) for increasing surplus size (the higher the surplus size, the more advantage an agent loses when they don't obtain $B$ in a particular time period with respect to agents who do obtain $B$ in that time period). For the current plots, we choose $\alpha = 0.2$ for visualization purposes, but this value is otherwise arbitrary.
"
end

# ╔═╡ 261432c5-3193-4acb-8988-f6de4f130a5c
html"<hr>"

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
	u=0.5,
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
begin
	step!(model, dummystep, sh.sharing_step!, model.total_ticks)

md"
In this ABM, agents are risk-pooling clusters that play the sharing game $T$ times per time period. At the end of the time period, clusters attract an individual with probabilities proportional to their per-season total payoff. If the population of individuals or clusters has reached their respective maximum sizes, an individual and/or cluster will be chosen at random and removed from the simulation (a multi-level birth-death process). As they grow, the average within-cluster payoff increases, peaks and then decreases due to overcrowding. When the average cluster payoff intersects the average long-term payoff of a loner a cluster will fission, expelling a random proportion $p \sim U(0,1)$ of individuals.

The new cluster inherits the sharing norm of the parent cluster it fissioned from with some variation, such that $s_\text{inherited} \sim \text{Normal}(s_\text{parent}, \sigma_\text{small})$. With probability $\sigma_\text{large}$ it chooses a random new sharing norm $s_\text{new} \sim U(0,1)$, deviating from the parent cluster's norm.

>Use the code in this cell to control the stepping procedure for the ABM.
"
end

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
begin
	#savefig(abm_hist_grid, "../images/figure7.png")
md"
Distribution of sharing norms and cluster sizes at the first step versus the last step of the simulation (weak selection, $\zeta = 0.01$). Sharign norms evolve to optimal values and the population keeps a reserve of loners. The distribution of cluster size widens due to the growth and fission process.
"
end

# ╔═╡ 532908ba-ed0d-4015-ac90-a8fd0325c7cc
html"<hr>"

# ╔═╡ 4efa0a3a-b935-4787-9186-4c0cb42dc060
md"
### Figure 8 - Time trajectories for sharing norms and cluster sizes
"

# ╔═╡ adb1cbd7-87c4-4755-919a-f3905cae2c5b
abm_time_plots = sh.abm_plots(model)

# ╔═╡ 953d65a7-f847-4001-8572-9df62002a947
begin
	#savefig(abm_time_plots, "../images/figure8.png")

md"
Time trajectories of mean/median sharing norms and cluster size, alongside the estimates produced by the math from the approximate models (green horizontal lines).

Cluster size for a process where clusters are growing (assume constant growth rates for simplicity) and fissioning when their size gets to $\bar{k}^* + 1$ will have the distribution

$K \sim U(1, \bar{k}^* + 1)$

where $U$ is the uniform distribution. Therefore, the expected value of cluster size is

$\mathbf{E}(K) = \frac{1 + (\bar{k}^* + 1)}{2} = \frac{\bar{k}^*}{2} + 1$

Results are shown for $\zeta = 0.01$ (weak selection).
"
end

# ╔═╡ 8c3c739d-58e3-44dd-8bc1-67cc78a78aea
html"<hr>"

# ╔═╡ f7e2e305-d635-430c-b88b-059e9ba43ca7
md"
### Figure 9 - Expected cluster size estimate from approximation vs. ABM simulations
"

# ╔═╡ 2ce5b10b-5522-4dec-99a2-6b8295077aa3
md"
The plots are generated using the file 'sharing\_clusters\_data\_Final.csv' in the data folder.
To run the simulations that generated this dataset, open a julia session in the paramscan folder and run

>include(\"paramscan\_sharing\_clusters.jl\")

Take note that the parameter scanning is very resource intensive and will take several days when run in a multi-core computational cluster, so we do not recommend attempting to run it to completion on a personal-use machine.
"

# ╔═╡ d71b0faf-0a78-47c1-9908-d7db9c1c5066
begin
	df = DataFrame( CSV.File("../data/sharing_clusters_data_Final.csv") )

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
begin
	#savefig(kbar_sweep_plot, "../images/figure9_unlabeled.png")
md"
Differences between weak and strong selection in the evolution of cluster size distribution. Under strong selection groups grow and fission at higher and more homogenous rates, which leads to good fits from the estimated expected cluster size. Under weak selection, larger clusters can exist for longer without fissioning at increasing surplus values, overshooting the estimate.
"
end

# ╔═╡ 340a6464-f123-4f37-826f-e7bd1a7c6435
html"<hr>"

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
begin
	#savefig(sstar_sweep_plot, "../images/figure10_unlabeled.png")
md"
Differences between weak and strong selection cases for the evolution of sharing norms. Under strong selection, cluster size is the main driver of the evolutionary process, while sharing norms are allowed to drift over a larger area (while still remaining over 0.5 on average). Curiously, high risk promotes the evolution of high sharing norms even in the strong selection case.
"
end

# ╔═╡ 6eff6aa8-b06f-40eb-8797-23b235b14946
html"<hr>"

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
md"
## C'est fini!
"
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─cb223404-03aa-11ee-0e4f-037151aa72ac
# ╟─f9c52e88-ef18-4bb5-976c-0eb909155570
# ╟─860a03ef-32ff-4563-9e25-380359069f49
# ╟─d6e849ab-77ff-45d6-8c50-b0117f95968e
# ╟─dcb1a731-4658-4699-bedf-6174ca8c76dd
# ╟─e7cc19cb-14c9-48be-bb5c-a6399f555320
# ╟─c454acc0-4edf-4529-a655-1be1b6084346
# ╟─1a86569a-f14f-4147-a489-cd5d4f24bd9f
# ╟─0c426787-3f6f-4b43-8140-f45752a49fe4
# ╟─342da70f-d2e4-43c5-b8b9-f97a24ca2d23
# ╟─875b643e-b5e7-4497-9b10-5945ae490fd1
# ╟─9e405011-abe4-4dd6-8ec7-bbb99e13cd3f
# ╟─2e28d384-2535-472f-bace-2d32f14aca81
# ╟─5d1c8435-7c7c-445b-a5f8-e6cef400bf70
# ╟─b694f876-d1f5-4977-b327-938ba31df886
# ╟─3cb7426e-cb2b-472a-a598-edbfa58bdb63
# ╟─eef9cd91-f1dd-40f2-aab6-954754d24f2d
# ╟─c508a66a-dda0-4be0-835b-8df1fecfa24c
# ╟─005a77d0-4b11-418f-9e18-5aaaaa061632
# ╟─1ca4b6c2-68c8-492f-88d4-e334af61d170
# ╟─13b589ad-40bd-4200-b0fb-e425bfeac925
# ╟─6060b2c8-87de-4ba7-b34b-c2f2c58305a8
# ╟─c5db3f92-b233-4f26-8b47-2af8d56ac132
# ╟─24156d86-873e-476e-9b87-55664bbfcfa8
# ╟─6e1100d8-3fd8-4453-9039-77054526452b
# ╟─d1e915da-a3f3-4d92-bf0f-80ede0485dc0
# ╟─f4ea2fd8-0196-47e6-bd4e-6f155f2a55d0
# ╟─7b18239f-a61a-4f48-80fb-3bb502a9a55b
# ╟─202fab42-ddf1-4f7d-bb43-ff4ea9b7b9e0
# ╟─4b82fe05-dbff-4ff5-9f8f-41a96df688c9
# ╟─261432c5-3193-4acb-8988-f6de4f130a5c
# ╟─d6e31887-51e2-4a9c-92e5-05143a1871a3
# ╟─ced5cd17-939b-456b-a990-4eebdab1f5c5
# ╟─1be4fb02-a98f-4500-b9cc-230ecde4e061
# ╟─6b10a2ec-d006-4e47-9098-8f72eed0795e
# ╟─0551d28b-b569-4983-b564-37adb611ac6c
# ╟─a2f3a45d-36d5-4add-b6db-60497f04de59
# ╟─363296c2-4a56-46d2-959d-39aff603e1be
# ╟─532908ba-ed0d-4015-ac90-a8fd0325c7cc
# ╟─4efa0a3a-b935-4787-9186-4c0cb42dc060
# ╟─adb1cbd7-87c4-4755-919a-f3905cae2c5b
# ╟─953d65a7-f847-4001-8572-9df62002a947
# ╟─8c3c739d-58e3-44dd-8bc1-67cc78a78aea
# ╟─f7e2e305-d635-430c-b88b-059e9ba43ca7
# ╟─2ce5b10b-5522-4dec-99a2-6b8295077aa3
# ╟─d71b0faf-0a78-47c1-9908-d7db9c1c5066
# ╟─afb1034b-db38-44fc-bfcd-5934afc6f51c
# ╟─340a6464-f123-4f37-826f-e7bd1a7c6435
# ╟─3b1b541b-8027-449d-9177-b452d2c181da
# ╟─aaf7d958-225f-43f2-b7ad-15fca6c58bbe
# ╟─5343ae11-baa6-4c6c-a84a-1ab1be027958
# ╟─6eff6aa8-b06f-40eb-8797-23b235b14946
# ╟─13da6b96-b1b0-426e-a4ec-0203f32834c3
