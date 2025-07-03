#!/usr/bin/env julia

using Arrow
using DataFrames
using GZip
using CSVFiles
#using Statistics
using Pipe: @pipe
#using StatsPlots

include("./InitialExploration/_initial_load.jl")
include("./MicrobeProfiling/_get_annotated_df.jl") # to use the tax_at_rank function

pip="amaw"
pip="holi"

target=".cpnr"
target=".norwaycpnr"
target=".HB_31.dedup.bam"
target=""

alg = ".local"
alg = ".lca";

add="_ANI95";
add="_ANI92";
add=""; 

tag = "$pip$target$alg$add"
tag="holi.lca";

# ------------------ YOU ONLY NEED TO RUN THIS ONCE ----------------------- #
write_arrow(tag)
#---------------------------------------------------------------------------#


# join metadata; filter based on N_reads and tax_rank
nreads = 50
ftax = "Tracheophyta"
ftax = "Embryophyta"
ftax = "Mammalia"
ftax = "Insecta"
ftax = "Chordata"
ftax = "Aves"
ftax = "Actinopterygii"
ftax = "Virus"
ftax = "Bacillariophyceae"
ftax = "Viridiplantae"
ftax = "Cyanobacteria"
ftax = ["Bacteria","Archaea"]
frank = "species" 

nreads=10;ftax="";frank="genus";

## data of samples
(x,y) = load_arrow_as_df(tag, ftax, frank, nreads)


# load gtdb tax_path metadata
path207=@pipe DataFrame(load(File(format"TSV", "/projects/wintherpedersen/people/tvg137/gtdb/bac120_metadata_r207_path.txt"))) |>
          append!(_, DataFrame(load(File(format"TSV", "/projects/wintherpedersen/people/tvg137/gtdb/ar53_metadata_r207_path.txt"))));
path207.accession = chop.(path207.accession, head=3, tail=0);

# load mapping results
tb = CSV.read("/projects/wintherpedersen/people/tvg137/HIBO_shotgun/EnvironmenTracker/archive_results/results_latest_amaw$add/taxonomic-profiling/tp-mapping-filtered.summary.tsv.gz", DataFrame, delim='\t');
# for microbe analysis: remove plastid and mito data to view only prokaryotes
filter!(:reference => r-> !contains(r, "_mito"),tb);
filter!(:reference => r-> !contains(r, "_plas"),tb);

# read acc2taxid file
using CSVFiles
acc2taxid = DataFrame(load(File(format"TSV", "/projects/wintherpedersen/data/vanilla-organelles-virus_7nov2022/pkg/taxonomy/acc2taxid.map.gz")))

if startswith(tag,"amaw.local")
    # if e.g. in updated metaDMG local output, tax_path etc is not available, run the following to keep only bacteria and archaea:
    filter!(:accession => r-> !contains(r, "_mito"),x);
    filter!(:accession => r-> !contains(r, "_plas"),x);
    filter!(:accession => r->length(r)==15,x);
    # add path to table
    leftjoin!(x, path207, on=:accession)
    rename!(x, :gtdb_taxonomy=>:tax_path)

    leftjoin!(x, tb, on=[:accession=>:reference,:Label=>:label])
    x.tax_name.=tax_at_rank.(x.tax_path, "species")
elseif startswith(tag,"amaw.lca")
    x.tax_name = [n[4:end] for n in x.tax_name]
    tb1 = @pipe tb |> leftjoin(_, acc2taxid[!,[:taxid, :accession]], on=:reference=>:accession)
    leftjoin!(x,tb1, on=[:tax_id=>:taxid, :Label=>:label])
end

#filter!(:read_ani_mean=>a->a>93.8,x)
using JLD2
save_object("./InitialExploration/data/$tag.$frank.jld2",x)




# if you are looking for some (not-so-interesting) plotting scripts, go to: obsolete/some_basic_plotting.jl 
