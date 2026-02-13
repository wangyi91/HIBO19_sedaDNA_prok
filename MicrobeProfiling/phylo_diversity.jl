#!/usr/bin/env julia
using Diversity, Phylo, DataFramesMeta
using CSV, JLD2

include("./MicrobeProfiling/_get_annotated_df.jl")
include("env.jl") # WORK_DIR
pip="amaw"; alg=".lca";add="_ANI92";rank="species";

tag=pip*alg*add
otudata = load_object("./deContamination§/data/$tag.$rank.jld2")
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

# copy taxlist_xxx.csv file as input for phyloT, e.g.:
`rsync -avP <path_in_HPC>/MicrobeProfiling/output/taxlist*.csv <path_local_PC>`

# use phyloT to make trees, with taxlist_xxx.csv as input. Save in newick format and move to R.pd/data/Bacteria or ./Archaea

# use R.pd/tree.R to remove bootstrap values 

# copy modified trees as input for phylo_distance.jl
`rsync -avP <path_local_PC>/R.pd/data/Bacteria/Bacteria_HB_{0..9}{0..9}.modified.newick <path_in_HPC>/MicrobeProfiling/output`




# --- Calculate per-sample Phylogenetic Diversity and merge ---
Libs=unique(otudata.Label)
pd = DataFrame(Label=Libs, pd_Bacteria=fill(0.0,34), pd_Archaea=fill(0.0,34))
Threads.@threads for (i,lib) in collect(enumerate(Libs))
    fil e= WORK_DIR*"HIBO_shotgun/analysis/MicrobeProfiling/output/$(ftax)_$(lib).modified.newick"
    if !isfile(file) continue end
    tree = open(parsenewick, Phylo.path(file))
    leafnames = getleafnames(tree)
    ph = PhyloBranches(tree)
    abundance = CSV.File(WORK_DIR*"HIBO_shotgun/analysis/MicrobeProfiling/output/tax_Nread_$(ftax)_$(lib).csv", header=0, delim=",") |> DataFrame
    rename!(abundance, [:tax, :N_read])
    abundance.tax = "s__".*replace.(abundance.tax, ' '=>'_')
    ab_sort = @rorderby abundance findfirst(==(:tax), leafnames)
    metaphylo = Metacommunity(ab_sort.N_read, ph)
    pd[i,"pd_$ftax"] = meta_gamma(metaphylo, 1).diversity[1]
    println(lib*" finished")
end
save_object("./MicrobeProfiling/output/pd.jld2",pd)

