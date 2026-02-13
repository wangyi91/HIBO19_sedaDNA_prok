#  Analyses of prokaryote sedimentary aDNA in Lake Constance
## Background

This repository contains scripts for data processing, analyses and visualisation for bacteria and archaea recorded in the sedimentary ancient DNA of Lake Constance.  

The manuscript of this study is provisionaly accepted in _The ISME Journal_.

## Prerequisite steps

Raw sequence data (fastq files) were processed with scripts in the directory `shotgun_data_processing`.

Then, fastq files were mapped against reference databases using the `aMAW-pipeline` (Fernandez-Guerra, in prep) for taxonomic profiling using the [GTDB database](https://gtdb.ecogenomic.org) (this study used the [r207 release](https://data.gtdb.ecogenomic.org/releases/release207/)), and DNA damage is subsequently estimated with `metaDMG v0.38.0`.

## Input files
The input files for analyses are: 
* File `tp-mdmg.amaw.lca_ANI92.csv.gz` output from `metaDMG`, saved under `./InitialExploration/data`
* Metadata are located in [metadata](https://github.com/wangyi91/sedaDNA-bacteria-archaea/tree/main/metadata)

## Analysis modules
### `InitialExploration`

`load_data.jl` reads input data and provide options for filtering based on multiple criteria such as taxonomic groups, minimum DNA reads per taxon, etc. 

`_initial_load.jl` contains functions that are needed.
### `deContamination` and `R.decontam`

### `LibrarySumary`

### `DamageAnalysis`
`plot_dmg_density.jl` makes plots of DNA damage. 

### `NetworkAnalysis`

### `MicrobeProfilling`
`R.pd` is for processing taxonomic trees and calculating phylogenetic diversity.

### `R.Bacon` for age-depth modeling
