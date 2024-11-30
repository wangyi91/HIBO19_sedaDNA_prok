using DataFrames
using JLD2
using Pipe


## Section 0: plot damage statistics for taxa with a series of N_reads, check if stable.
function calculate_threshold_stats(data::DataFrame, steps::StepRange{Int64, Int64})
       # data, of one library only, contains a column of damage, one or multiple column(s) of read counts for filtering and weighting
       # steps is a range of filtering threshold
       tab = DataFrame(Label=Vector{Any}(missing, size(steps,1)),
                       min_Nreads=Vector{Int64}(undef, size(steps,1)),
                       Mean=Vector{Float64}(undef, size(steps,1)),
                       Median=Vector{Float64}(undef, size(steps,1)),
                       Std=Vector{Float64}(undef, size(steps,1)),
                       N_taxa = Vector{Int64}(undef, size(steps,1)))

       for (i,n) in enumerate(steps)
           sub = filter(:N_reads => x -> x > n, data)
           if size(sub,1) <= 3
               break
           end
           println("Calculating for $(data.Label[1]), $n reads thres, $(size(sub,1)) taxa")
           weights = Weights(Float64.(sqrt.(sub.N_reads)))
           tab.Mean[i] = mean(Float64.(sub.damage), weights)
           tab.Median[i] = median(Float64.(sub.damage), weights)
           tab.Std[i] = std(Float64.(sub.damage), weights)
           tab.min_Nreads[i] = n
           tab.N_taxa[i] = size(sub,1)
           tab.Label[i] = data.Label[1]
       end    
       return(dropmissing(tab))
end


libs = load("./metadata/HIBO_library_metadata_dated.csv") |> DataFrame |> dropmissing
trtax = "Bacteria"; cl=:dodgerblue; pip = "amaw"; frank="species";
trtax = "Archaea"; cl=:goldenrod2; frank="genus";
steps = 10:50:2010
outdir = "/projects/wintherpedersen/people/tvg137/HIBO_shotgun/analysis/output_indiv"

Threads.@threads for (i,lib) in collect(enumerate(libs.Label))
    df = load_object("./DmgMixtureModel/data/obs_$(lib)_$(pip)_$(frank).jld2")
    
    training_data = filter(:tax_path => t -> occursin("$trtax", t), df)
    if lib in ["HB_12","HB_18"]  && trtax == "Viridiplantae" # remove obvious contaminants in HB_12 and HB_18
        training_data = @pipe training_data |> filter(:damage => d -> d > 0.05, _)
    elseif libs.yrBP[i] > 7500 && trtax == "Bacteria" # remove modern contaminants
        training_data = @pipe training_data |> filter(:damage => d -> d > 0.05, _)
    end

    output = calculate_threshold_stats(training_data, steps)

    gr(margins = 1.5Plots.cm, size=(500, 625))
    plot(output.min_Nreads, output.N_taxa, seriestype=:scatter, marker = (stroke(0)), 
         label="genus", leg=:topleft, yminorticks=5, ylabel = "number of taxa", xlabel = "minimum N_reads")
    @df output plot!(twinx(), :min_Nreads, :Mean, ribbon = :Std, fillalpha = 0.3,
                        ylims=(0,0.4), label="mean damage ± std", color=cl,
                        ylabel = "damage",
                        title="$(:Label[1]) $trtax")
                        savefig(outdir*"/thresholds_$(output.Label[1])_$(trtax)_ani90.png");
end

## Section 1: plot damage histogram of taxa N_reads > 1500  for each sample
        # use plot_dmg_by_sample.jl

## Section 2: Distribution is usually bimodal. Manually decide filtering thresholds for each sample, assisted by correlation with depth.
        # use threshold_interpolation.jl
        
## Section 3: Filter data for each sample

ftax = "Bacteria"; pip = "amaw"; frank="species";
ftax = "Archaea"; pip = "amaw"; frank="species"; frank="genus";

add=""; add2="";
add="_ANI95"; add2="_ANI95_strict";

stats = load("./DmgMixtureModel/output$(add)/$(pip)_$(ftax)_$(frank)_filter_thres_after$(add2).csv") |> DataFrame # this are manually + regression set thresholds

libs = load("./metadata/HIBO_library_metadata_dated.csv") |> DataFrame |> dropmissing

# option1: filter for ancient taxa
for (i,lib) in enumerate(stats.Label)
    df = load_object("./DmgMixtureModel/data$(add)/obs_$(lib)_$(pip)_$(frank).jld2")
    filtered_dt = @pipe df |> filter(:damage => d -> d > stats.lower[i], _) |>
                              filter(:tax_path => t -> occursin("$ftax", t), _)
    save_object("./DmgMixtureModel/output$(add)/obs_$(lib)_$(pip)_$(ftax)_$(frank)$(add2)_filtered.jld2", filtered_dt)
end

# option2: filter for non-ancient taxa
for (i,lib) in enumerate(stats.Label)
    df = load_object("./DmgMixtureModel/data$(add)/obs_$(lib)_$(pip)_$(frank).jld2")
    filtered_dt = @pipe df |> filter(:damage => d -> d <= stats.lower[i], _) |>
                              filter(:tax_path => t -> occursin("$ftax", t), _)
    save_object("./DmgMixtureModel/output$(add)/obs_$(lib)_$(pip)_$(ftax)_$(frank)$(add2)_dropped.jld2", filtered_dt)
end

# option3: omit damage filtering, just save each sample for $ftax
for lib in libs.Label
    df = load_object("./DmgMixtureModel/data$(add)/obs_$(lib)_$(pip)_$(frank).jld2")
    filtered_dt = @pipe df |> filter(:tax_path => t -> occursin("$ftax", t), _)
    save_object("./DmgMixtureModel/output$(add)/obs_$(lib)_$(pip)_$(ftax)_$(frank).jld2", filtered_dt)
end







