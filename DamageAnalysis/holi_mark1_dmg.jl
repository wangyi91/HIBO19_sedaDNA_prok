using StatsBase
using DataFrames
using JLD2
using Pipe

# This script calculates the mean and std of damage for abundant plant taxa, and use these stats as filtering criteria for all euk taxa

function calculate_threshold_stats(data::DataFrame, steps::StepRange{Int64, Int64})
       # "data", of one library only, contains a column of damage, one or multiple column(s) of read counts for filtering and weighting
       # steps is a range of filtering threshold
       tab = DataFrame(Label=Vector{Any}(missing, size(steps,1)),
                       min_Nreads=Vector{Int64}(undef, size(steps,1)),
                       Mean=Vector{Float64}(undef, size(steps,1)),
                       Median=Vector{Float64}(undef, size(steps,1)),
                       Std=Vector{Float64}(undef, size(steps,1)),
                       N_taxa = Vector{Int64}(undef, size(steps,1)))
       for (i,n) in enumerate(steps)
           sub = filter(:N_reads => x -> x > n, data)
           if size(sub,1) <= 3 break end
           println("Calculating for $(data.Label[1]), $n reads thres, $(size(sub,1)) taxa")
           #weights = Weights(Float64.(sqrt.(sub.N_reads)))
           tab.Mean[i] = mean(Float64.(sub.damage)) # if weighted: mean(Float64.(sub.damage), weights)
           tab.Median[i] = median(Float64.(sub.damage))
           tab.Std[i] = std(Float64.(sub.damage))
           tab.min_Nreads[i] = n
           tab.N_taxa[i] = size(sub,1)
           tab.Label[i] = data.Label[1]
       end    
       return(dropmissing(tab))
end


# Step 1: Plot guide data statistics for each sample with incrementing min read counts
libs = load("./metadata/HIBO_library_metadata_dated.csv") |> DataFrame |> dropmissing
tax = "Viridiplantae"; cl=:olivedrab; pip = "holi"; frank="genus";
dt = filter(row->any(occursin.(tax,row.tax_path)), x)
steps = 10:50:2010

using StatsPlots
Threads.@threads for (i,lib) in collect(enumerate(libs.Label))
    guide_data = filter(:Label => l -> l==lib, dt)
    if lib in ["HB_12","HB_18"]  # remove obvious contaminants in HB_12 and HB_18
        filter!(:damage => d -> d > 0.05, guide_data)
    end
    output = calculate_threshold_stats(guide_data, steps)
    gr(margins = 1.5Plots.cm, size=(500, 625))
    StatsPlots.plot(output.min_Nreads, output.N_taxa, seriestype=:scatter, marker = (stroke(0)), 
         label="genus", leg=:topleft, yminorticks=5, ylabel = "number of taxa", xlabel = "minimum N_reads")
    @df output StatsPlots.plot!(twinx(), :min_Nreads, :Mean, ribbon = :Std, fillalpha = 0.3,
                        ylims=(0,0.4), label="mean damage ± std", color=cl,
                        ylabel = "damage", title="$lib $tax $frank")
    savefig("./output_indiv/thresholds_$(lib)_$tax.pdf");
end



# Step 2: Decide for each sample a threshold, Mean±Std, using training data.

## guide data. n is minimum read count in guide data, chosen based on plots output from Step 1
tax = "Viridiplantae"; n = 600; pip="holi";

## first find b and c for each sample
stats = DataFrame(Label=Vector{Any}(missing, size(libs,1)),
                  mean_dmg_guide=Vector{Float64}(undef, size(libs,1)),
                  median_dmg_guide=Vector{Float64}(undef, size(libs,1)),
                  std_dmg_guide=Vector{Float64}(undef, size(libs,1)))

Threads.@threads for (i,lib) in collect(enumerate(libs.Label))
    guide_data = filter(:Label => l -> l==lib, dt)
    if lib in ["HB_12","HB_18"]  && tax == "Viridiplantae" # remove obvious contaminants in HB_12 and HB_18
        filter!(:damage => d -> d > 0.05, guide_data)
    end
    
    thres = calculate_threshold_stats(guide_data, n:1:n)
    stats.Label[i] = lib
    stats.mean_dmg_guide[i] = thres[1,:Mean]
    stats.median_dmg_guide[i] = thres[1, :Median]
    stats.std_dmg_guide[i] = thres[1, :Std]
end

transform!(stats, [:mean_dmg_guide,:std_dmg_guide] => ByRow((m,s) -> m-1.5s) => :lower)
transform!(stats, [:mean_dmg_guide,:std_dmg_guide] => ByRow((m,s) -> m+3s) => :higher)

# when needed: adjust negative lower boundaries
global_lower = minimum(l for l in stats.lower[stats.lower.>0])
transform!(stats, :lower => ByRow(l -> max(l, global_lower)) => :lower_adjusted)

save_object("./DmgMixtureModel/output/dmg_filter_thres_holi_euk.jld2", stats)



## Section 3: For each sample, label taxa as "damaged" or "not damaged"
stats = load_object("./DmgMixtureModel/output/dmg_filter_thres_holi_euk.jld2")
pip="holi";tax="";target="";alg = ".lca";add="";
tag = "$pip$tax$target$alg$add"
frank="genus";#frank="family";
x = load_object("./InitialExploration/data/$tag.$frank.jld2")

fdt = @pipe transform(x, [:Label,:damage]=>ByRow((l,d)->stats[stats.Label.==l,:lower][1]<d<stats[stats.Label.==l,:higher][1])=>:damaged) 
            #|> filter(row->any(occursin.(ftax,row.tax_path)), x) 
            
#summ = @pipe combine(groupby(fdt,:tax_name),nrow=>:noccur,:damaged=>sum=>:nvalid,:N_reads=>sum=>:nreads)|>
#        transform(_,[:noccur,:nvalid]=>ByRow((o,v)->v/o)=>:pvalid) |>
#        transform(_, :pvalid=>ByRow(p->p>0.4)=>:valid_taxa)

#@pipe filter(:pvalid=>p->p<0.4,summ)|> sort(_,:noccur) |>show(_,allrows=true)

#transform!(fdt, :tax_name=>ByRow(n->n in summ[summ.valid_taxa.==true,:tax_name])=>:valid)
save_object("./DmgMixtureModel/output/$tag.$frank.marked.dmg.jld2", fdt)
