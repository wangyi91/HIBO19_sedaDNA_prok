#  Analyses of prokaryote sedimentary aDNA in Lake Constance
## Background

This repository contains scripts for data analyses and visualisation for bacteria and archaea recorded in the sedimentary ancient DNA of Lake Constance.  

The manuscript of this study is currently under review in _The ISME Journal_.

## Prerequisite steps

Before using scripts in this repository, raw sequence data (fastq files) should be first processed with scripts in the repository [shotgun-data-processing](https://github.com/wangyi91/shotgun-data-processing.git). 

Then, fastq files are mapped against reference databases using the [aMAW-pipeline] for taxonomic profiling using the [GTDB database](https://gtdb.ecogenomic.org) (this study used the [r207 release](https://data.gtdb.ecogenomic.org/releases/release207/)), and DNA damage is subsequently estimated with [metaDMG](https://github.com/metaDMG-dev/metaDMG-core).

## Input files
The input files for analyses in this repository are: 

* The output file `results/taxonomic-profiling-dmg/tp-mdmg.lca.weight-1.csv.gz` of `metaDMG`. 
* Metadata are located in [metadata](https://github.com/wangyi91/sedaDNA-bacteria-archaea/tree/main/metadata)

## Steps of sequence analyses
### [InitialExploration](https://github.com/wangyi91/sedaDNA-bacteria-archaea/tree/main/InitialExploration)
This step is necessary for all other analyses. 

`load_data.jl` reads input data and provide options for filtering based on multiple criteria such as taxonomic groups, minimum DNA reads per taxon, etc. 

`_initial_load.jl` contains functions that are needed.

### [DamageAnalysis](https://github.com/wangyi91/sedaDNA-bacteria-archaea/tree/main/DamageAnalysis)
`plot_dmg_density_on_single_plot.jl` makes plots of DNA damage. `living_species.jl` defines living microbes.


The following analyses can be done in parallel to each other:
### [NetworkAnalysis](https://github.com/wangyi91/sedaDNA-bacteria-archaea/tree/main/NetworkAnalysis)

### [MicrobeProfilling](https://github.com/wangyi91/sedaDNA-bacteria-archaea/tree/main/MicrobeProfilling)

## Age-depth model
* [OxCal_model](https://github.com/wangyi91/sedaDNA-bacteria-archaea/tree/main/OxCal_model)
* [AgeDepthModel](https://github.com/wangyi91/sedaDNA-bacteria-archaea/tree/main/AgeDepthModel)
