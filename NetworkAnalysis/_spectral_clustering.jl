#!/usr/bin/env julia
using DataFrames, FlashWeave
using Pipe: @pipe
using LinearAlgebra, Distances, Clustering, Random

# This function uses spectral clustering to group communities identified through clique-percolation
# Code modified from https://rpubs.com/gargeejagtap/SpectralClustering

function sp_cluster(cliper::Vector{BitSet}, cliq_idx::Vector{Int64}, otudata::DataFrame, nc::Int64)
    
    # prepare a table where rows are samples, columns are each clique (community), values are number of taxa
    tab=DataFrame(Label=String[], sum_reads=Int64[], n_taxa=Int64[], N_clq=Int64[])
    for (i,clq) in enumerate(collect.(cliper[cliq_idx]))
        splist = collect(skipmissing(tax_at_rank.(names(netw)[clq], "species")))
        sub = @pipe otudata |> filter(:tax_name => n-> n in splist, _) |>
                    combine(groupby(_,:Label), :N_reads=>sum=>:sum_reads, nrow=>:n_taxa) |>
                    insertcols(_,:N_clq=>i)
        append!(tab, sub)
    end
    wide = unstack(tab, :N_clq, :Label, :n_taxa, fill=0)

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
    # Since nc > 2 clusters, we normalize the Laplacian matrix:
    laplacian = (Diagonal(ones(n)) - Diagonal(degrees.^(-1/2)) * D * Diagonal(degrees.^(-1/2)) )

    ## Compute eigenvectors
    eigenvectors = eigen(laplacian).vectors[:,2:(nc+1)] #2:(nc+1)

    ## Run k-means on eigenvectors and repeat runs=1000 times to address randomness. 
    ## Find a consensus by applying a threshold in the frequencies of two cliq-communities being connected.
    # Initialize a co-occurrence matrix to zero (n x n)
    co_occurrence_matrix = zeros(Int, n, n)

    Random.seed!(1234)
    runs=1000
    for _ in 1:runs
        # Run Kmeans on eigenvectors
        result = kmeans(transpose(eigenvectors), nc; maxiter=100)  
        labels = result.assignments  # cluster labels for each vector

        # Update the co-occurrence matrix
        for i in 1:n
            for j in i+1:n
                if labels[i] == labels[j]
                    co_occurrence_matrix[i, j] += 1
                    co_occurrence_matrix[j, i] += 1  # ensure symmetry
                end
            end
        end
    end

    # Normalize co-occurrence matrix by the number of runs
    co_occurrence_matrix = co_occurrence_matrix / runs

    threshold = 0.9  

    # Thresholded adjacency matrix
    adjacency_matrix = co_occurrence_matrix .>= threshold

    # Find clusters as connected components
    clusters = conn_components(adjacency_matrix) # definition see below

    return(clusters)
end


function conn_components(adjacency_matrix)
    n = size(adjacency_matrix, 1)
    visited = falses(n)
    components = []
    
    function dfs(v, component)
        visited[v] = true
        push!(component, v)
        for u in 1:n
            if adjacency_matrix[v, u] && !visited[u]
                dfs(u, component)
            end
        end
    end
    
    for i in 1:n
        if !visited[i]
            component = []
            dfs(i, component)
            push!(components, component)
        end
    end
    
    return components
end

