using Plots, Roots, Distributions, Colors, Random, Statistics, DataFrames, CSV

using LambertW, ProgressLogging

W_s075(kbar; u=0.5, B=10, C=0.1, β=1, s=0.75) = kbar ≥ 1 ? u*log( (1-s)*B + (1 - (1 - u)^(β*kbar))*s*u*β*B - β*kbar*C ) + (1-u)*log( 1 + (1 - (1 - u)^(β*kbar))*s*u*β*B - β*kbar*C ) : u*log(B)

W_s05(kbar; u=0.5, B=10, C=0.1, β=1, s=0.5) = kbar ≥ 1 ? ( 1 + (1 - (1 - u)^(β*kbar))*s*u*β*B - β*kbar*C > 0 ?  u*log( (1-s)*B + (1 - (1 - u)^(β*kbar))*s*u*β*B - β*kbar*C ) + (1-u)*log( 1 + (1 - (1 - u)^(β*kbar))*s*u*β*B - β*kbar*C ) : 0 ) : u*log(B)

optshare(b, k; u=0.5, c=0.1, β=1) = ( ( (β*k*c - 1) / ( 1 - (1-u)^(β*k) )*β ) + ( (1-u)*(b - β*k*c) / ( 1 - u*β*( 1 - (1-u)^(β*k) ) ) ) ) / b

optshare2(b, k, u, c; β=1) = ( ( (β*k*c - 1) / ( 1 - (1-u)^(β*k) )*β ) + ( (1-u)*(b - β*k*c) / ( 1 - u*β*( 1 - (1-u)^(β*k) ) ) ) ) / b

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

function maxshare_Du(
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

function calculate_points_Du(;B=10, C=0.25, β=1, s=9/10, δ_0 = 0.2 + log(B - 1))

	kbarW = Vector{Float64}()
	koptV = Vector{Float64}()
	kbarV = Vector{Float64}()
	
	for u in (0.001:0.01:1)

		W_dif_iter(kbar; u=u, B=B, C=C, β=β, s=s) = u*log( (1-s)*B + s*u*β*B - β*kbar*C ) + (1-u)*log( 1 + s*u*β*B - β*kbar*C ) - u*log(B)
		
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
		
		push!(kbarV, kbar)
		
	end

	return [kbarW, kbarV, koptV]
	
end

function zero_plots(kbar; leg=true, lab=true, notes=true, ylim=35, yax=true)
	
	plot(
		collect(0.01:0.01:1),
		kbar[2],
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
		kbar[3],
		label="",
		#line=:dash,
		color=RGBA(0.9,0.5,0.5,0.9),
		lw=2.5
	)
end

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

function calculate_points_Sim(s, b, c, β; 
	T=500, 
	v_0=1, 
	w_0=1,
	max_kbar=35,
	reps=100,
	)

	axis = 0:0.01:1|>collect
	k_hat = Vector{Float64}()
	k_star = Vector{Float64}()
		
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

			push!(k_hat, fission_point)
			push!(k_star, max_arg)
			
		end
	
	return [k_hat, k_star]
end
