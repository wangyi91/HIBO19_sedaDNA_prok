#!/usr/bin/env julia

using JLD2, DataFrames
using CSVFiles, CSV
using CairoMakie, StatsPlots
using Pipe: @pipe

pip="amaw"; alg=".lca";add="_ANI90";
tag=pip*alg*add

x=load_object("./InitialExploration/data/$tag.jld2")


# plot ANI density by sample
fig=Figure()
Axis(fig[1,1], xlabel="Mean read ANI (%) of bacterial and archaeal species",xgridvisible=false,ygridvisible=false)
for lib in unique(x.Label)
    @pipe x |> filter(:Label=>l->l==lib,_) |>
    CairoMakie.density!(Vector{Float64}(_.read_ani_mean),color=(:white,0), strokecolor=ifelse(lib=="HB_22",:red,(:dodgerblue4,0.5)), strokewidth=1); # bins=0:0.005:0.45
end
save("MicrobeProfiling/output/density_ani_all_beta.pdf",fig);

# below are not included:





# plot breath vs. ANI by sample
for lib in libs.Label
    gr(margins = 1.5Plots.cm, size=(500, 625))
    @pipe tb |> filter(:label => l -> l==lib,_) |>
    filter(:n_reads => n->n>100,_) |>
    plot(_.read_ani_mean, _.breadth, seriestype=:scatter, mode="markers",
    marker=(log1p.(_.n_reads), 0.3, stroke(0)), marker_sizeref=0.3, label=lib); # bins=0:0.005:0.45
    xlabel!("mean read ANI of taxa, n_reads>100");
    ylabel!("breadth");
    savefig("MicrobeProfiling/output/ani_vs_breadth_$(lib).png");
end


# plot bratio vs. ANI by sample
for lib in libs.Label
    gr(margins = 1.5Plots.cm, size=(500, 625))
    @pipe tb |> filter(:label => l -> l==lib,_) |>
    filter(:n_reads => n->n>100,_) |>
    plot(_.read_ani_mean, _.breadth_exp_ratio, seriestype=:scatter, mode="markers",
    marker=(log1p.(_.n_reads), 0.3, stroke(0)), marker_sizeref=0.3, label=lib); # bins=0:0.005:0.45
    xlabel!("mean read ANI of taxa, n_reads>100");
    ylabel!("breadth_exp_ratio");
    savefig("MicrobeProfiling/output/ani_vs_bratio_$(lib).png");
end












# filtering criteria for sources
nreads = 200
ftax = "Bacteria"
frank = "species"

# filtered source data
src = @pipe tb |> filter(:tax_path => t -> occursin("$ftax", t), _) |>
          filter([:N_reads, :tax_rank] => (x,y) -> x > nreads && y==("$frank"), _) |>
          dropmissing(_, :damage) |>
          rename(_, :sample => :Label)

# load cleaned sink data
cd("/projects/wintherpedersen/people/tvg137/HIBO_shotgun/analysis") 
## then run ../InitialExploration/vcat_samples.jl and get object: data_filtered
cd("/projects/wintherpedersen/people/tvg137/HIBO_shotgun/SourceTracker2/prepare_st2_input")

# further filter sample (sink) data by read count
filter!(:N_reads => n -> n>200, data_filtered)



# cat source and sink data
otu_long = vcat(src, data_filtered, cols=:intersect)

# tranform to matrix as required by SourceTracker2
otu = unstack(otu_long, :tax_name, :Label, :N_reads, fill=0)

# output txt file
save(File(format"CSV", "../data/otu.txt"), otu; delim='\t',quotechar=nothing, header=true)
# then run under ../data the following command to convert to biom format:
conda activate st2
biom convert -i ../data/otu.txt -o ../data/otu_hdf5.biom --table-type="OTU table" --to-hdf5




