#!/usr/bin/env julia
using Diversity, Phylo, DataFramesMeta
using CSV, JLD2
include("./MicrobeProfiling/_get_annotated_df.jl")

pip="amaw"; alg=".lca";add="_ANI92";
tag=pip*alg*add
otudata = load_object("./InitialExploration/data/$tag.jld2")
ftax = "Bacteria"; ftax = "Archaea";

########### Plot taxonomy tree using GTDB taxonomy ##########

# output species names (as input for phyloT) and read count
using DelimitedFiles
# write taxa list for each sample
for lib in unique(otudata.Label)
    tmp = @pipe filter([:tax_path, :Label]=> (p,l)-> occursin(ftax,p) && l==lib, otudata)|> filter(:N_reads=>n->n>100,_)
    writedlm("./MicrobeProfiling/output/taxlist_$(ftax)_$(lib).csv",
             String.(tmp.tax_name))
    writedlm("./MicrobeProfiling/output/tax_Nread_$(ftax)_$(lib).csv",
             [String.(tmp.tax_name) tmp.N_reads], ',')
end

# taxa from all samples in one list
tmp = filter(:tax_path=>p->occursin(ftax,p),otudata);
writedlm("./MicrobeProfiling/output/taxlist_$(ftax).csv", String.(unique(tmp.tax_name)))

# send taxlist_xxx.csv file to local
rsync -avP tvg137@dandycomp01fl:/projects/wintherpedersen/people/tvg137/HIBO_shotgun/analysis/MicrobeProfiling/output/taxlist*.csv /Users/yiwang/PhD_Research/AG_Epp/Projects/Long_core/11_data_analysis/MicrobeProfiling/output/

# use phyloT to make trees, save in newick format and move to
# /Users/yiwang/PhD_Research/AG_Epp/Projects/Long_core/11_data_analysis/R.pd/data/Bacteria or ./Archaea

# use R to remove bootstrap values 
# /Users/yiwang/PhD_Research/AG_Epp/Projects/Long_core/11_data_analysis/R.pd/tree.R

# send output trees back to julia
rsync -avP /Users/yiwang/PhD_Research/AG_Epp/Projects/Long_core/11_data_analysis/R.pd/data/Bacteria/Bacteria_HB_{0..9}{0..9}.modified.newick tvg137@dandycomp01fl:/projects/wintherpedersen/people/tvg137/HIBO_shotgun/analysis/MicrobeProfiling/output



# calculate PD by sample
Libs=unique(otudata.Label)
pd = DataFrame(Label=Libs, pd_Bacteria=fill(0.0,34), pd_Archaea=fill(0.0,34))
Threads.@threads for (i,lib) in collect(enumerate(Libs))
    file="/maps/projects/wintherpedersen/people/tvg137/HIBO_shotgun/analysis/MicrobeProfiling/output/$(ftax)_$(lib).modified.newick"
    if !isfile(file) continue end
    tree = open(parsenewick, Phylo.path(file))
    leafnames = getleafnames(tree)
    ph = PhyloBranches(tree)
    abundance = CSV.File("/maps/projects/wintherpedersen/people/tvg137/HIBO_shotgun/analysis/MicrobeProfiling/output/tax_Nread_$(ftax)_$(lib).csv", header=0, delim=",") |> DataFrame
    rename!(abundance, [:tax, :N_read])
    abundance.tax = "s__".*replace.(abundance.tax, ' '=>'_')
    ab_sort = @rorderby abundance findfirst(==(:tax), leafnames)
    metaphylo = Metacommunity(ab_sort.N_read, ph)
    pd[i,"pd_$ftax"] = meta_gamma(metaphylo, 1).diversity[1]
    println(lib*" finished")
end
save_object("./MicrobeProfiling/output/pd.jld2",pd)


# calculate PD for subtrees
ranknames_file=CSV.File("/maps/projects/wintherpedersen/people/tvg137/HIBO_shotgun/analysis/MicrobeProfiling/output/ranknames_$(ftax)$(add3).csv", header=1, delim=",") |> DataFrame; # rankname file is generated with R.pd/tree.R
ranknames=unique(ranknames_file[!,"."])

for rank in ranknames, lib in Libs
    file="/maps/projects/wintherpedersen/people/tvg137/HIBO_shotgun/analysis/MicrobeProfiling/output/$(ftax)$(add3)_$(lib)_$(N).$(rank).subtree.newick"
    if !isfile(file) continue end
    tree = open(parsenewick, Phylo.path(file))
    leafnames = getleafnames(tree)
    ph = PhyloBranches(tree)
    abundance = CSV.File("/maps/projects/wintherpedersen/people/tvg137/HIBO_shotgun/analysis/MicrobeProfiling/output/tax_Nread_$(ftax)$(add3)_$(lib)_$(N).csv", header=0, delim=",") |> DataFrame
    rename!(abundance, [:tax, :N_read])

    abundance.tax = replace.(abundance.tax, ' '=>'_')
    ab_sub = filter(:tax => t -> t in leafnames, abundance)
    
    ab_sort = @rorderby ab_sub findfirst(==(:tax), leafnames)
    metaphylo = Metacommunity(ab_sort.N_read, ph)
    print(join([lib, rank, meta_gamma(metaphylo, 1).diversity[1], '\n'], ','))
end

# finally, save printout as subtree_PD_$(ftax)$(add3).csv to local, and plot using subtree_PD.R
