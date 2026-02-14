#!/usr/bin/env julia
using DataFrames
using JLD2

pip="amaw"
add="_ANI92";
alg=".lca"

# Placeholder for the final result
summary_df = DataFrame()

# Loop through each taxonomic rank
#for alg in [".local",".lca"]
for rk in ["phylum","family","genus","species"]    
    # Load the dataset
    dt = load_object("./InitialExploration/data/$pip$alg$add.$rk.jld2")

    # Classify each row by domain based on `tax_path`
    dt.Domain = ifelse.(occursin.("Bacteria", dt.tax_path), "Bacteria",
                  ifelse.(occursin.("Archaea", dt.tax_path), "Archaea", "Other"))

    # Filter out "Other" if you want only Bacteria/Archaea
    filter!(row -> row.Domain != "Other", dt)

    # Group by Label and Domain, then sum reads
    grouped = combine(groupby(dt, [:Label, :Domain]), :N_reads => sum => :Total_Reads)

    # Pivot the domain to columns (e.g., "phylum_Bacteria", "phylum_Archaea")
    pivoted = unstack(grouped, :Domain, :Total_Reads)
    rename!(pivoted, Dict("Bacteria" => "$(rk)_Bacteria", "Archaea" => "$(rk)_Archaea"))

    # Merge into the summary DataFrame
    if isempty(summary_df)
        summary_df = pivoted
    else
        summary_df = outerjoin(summary_df, pivoted, on = :Label)
    end
end


# Rename :Label to :Sample if you prefer that terminology
rename!(summary_df, :Label => :Sample)

# Show result
println(summary_df)




