#!/usr/bin/env julia
using DataFrames, FlashWeave
using Pipe: @pipe
using LinearAlgebra, Distances, Clustering, Random

# This function uses spectral clustering to group communities identified through clique-percolation
# Code modified from https://rpubs.com/gargeejagtap/SpectralClustering

function sp_cluster(cliper::Vector{BitSet}, cliq_idx::Vector{Int64}, otudata::DataFrame, nc::Int64, sd::Int64)
    
    # prepare a table where rows are samples, columns are each clique (community), values are number of taxa
    tab=DataFrame(Label=String[], sum_reads=Int64[], n_taxa=Int64[], N_clq=Int64[])
    for (i,clq) in enumerate(collect.(cliper[cliq_idx]))
        splist = collect(skipmissing(tax_at_rank.(names(netw)[clq], "species")))
        sub = @pipe otudata |> filter(:tax_name => n-> n in splist, _) |>
                    combine(groupby(_,:Label), :N_reads=>sum=>:sum_reads, nrow=>:n_taxa) |>
                    insertcols(_,:N_clq=>i)
        append!(tab, sub)
    end
    wide = unstack(tab, :N_clq, :Label, :sum_reads, fill=0)

    # Create Similarity matrix of Euclidean distance between communities
    S = pairwise(euclidean,transpose(Matrix{Float64}(wide[!,2:end])))
    n = size(S)[1]
    
    # Create Degree matrix
    D = Matrix{Int64}(zeros(n,n))

    for i in 1:n
      # Find top 10 nearest neighbors using Euclidean distance
      index = sortperm(S[i,:])[2:11]

      # Assign value to neighbors
      D[i,index] .= 1
    end

    # find mutual neighbors
    D = D + transpose(D)
    D= (D.>0).*1

    degrees = eachcol(D) |> sum
    ## Compute Laplacian matrix
    # Since k > 2 clusters, we normalize the Laplacian matrix:
    laplacian = (Diagonal(ones(n)) - Diagonal(degrees.^(-1/2)) * D * Diagonal(degrees.^(-1/2)) )

    ## Compute eigenvectors
    eigenvectors = eigen(laplacian).vectors[:,2:(nc+1)] #2:(nc+1)

    ## Run Kmeans on eigenvectors
    Random.seed!(sd)
    sc = kmeans(transpose(eigenvectors), nc)

    #nclusters = 3:15; data = transpose(eigenvectors)
    #clusterings = kmeans.(Ref(data), nclusters)
    #test=clustering_quality.(Ref(data), clusterings, quality_index = :dunn)

    #return(test)
    return(sc.assignments)
end

