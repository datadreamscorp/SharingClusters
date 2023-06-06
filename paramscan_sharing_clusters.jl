using Pkg
Pkg.activate(".")

@everywhere begin #INCLUDE MODEL CODE AND NECESSARY LIBRARIES

	import SharingClusters as sh
	using Agents, Random, Distributions, Statistics, StatsBase
	
	B = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 20.0, 40.0]

	parameters = Dict( #ALTER THIS DICTIONARY TO DEFINE PARAMETER DISTRIBUTIONS
	    :B => B,
		:bc_ratio => 100.0,
		:n => [200],
		:u => collect(0.01:0.01:0.99),
		:rep => collect(1:100),
		:size_bias => true,
		:w0 => 0.5,
		:δ => [exp(1) - 1, 0.01]
	)


	mdata = [
		:current_n,
		:current_N,
		:max_cluster_size,
		:mean_sharing,
		:median_sharing,
		:mean_cluster_size,
		:median_cluster_size,
		:C,
	]

end

#USE THIS LINE AFTER DEFINITIONS TO BEGIN PARAMETER SCANNING
_, mdf = paramscan(
            parameters, sh.sharing_groups_model;
            mdata=mdata,
            agent_step! = dummystep,
        	model_step! = sh.sharing_step!,
            n = 5000,
			parallel=true,
			when_model = collect(0:100:5000),
			showprogress = true
	)

using CSV
CSV.write("./data/sharing_clusters_data_Replicate.csv", mdf)
