#!/usr/bin/env julia
using Arrow
using DataFrames
using CSVFiles
using Pipe: @pipe

# This function loads output data from metaDMG (stored as arrow file for easy julia access) 
# and does basic filtering on tax, taxonomic rank and minimum reads

function write_arrow(tag::String)
    # load the data and save in Arrow format
    df = load(File(format"CSV", "./InitialExploration/data/tp-mdmg.$tag.csv.gz"));

    Arrow.write("./InitialExploration/data/$tag.arrow", df; compress = :lz4)
    println("Arrow file written as ./InitialExploration/data/$tag.arrow")
end

function load_arrow_as_df(tag::String, ftax, frank::String, nreads::Int)
    # read Arrow files
    tb = Arrow.Table("./InitialExploration/data/$tag.arrow");

    # create a df with only data that need to be checked
    ## here I rename "sample" to "Label"
    ## if using holi output, I relabel libraries
    if startswith(tag,"amaw.local")
        df = @pipe DataFrame([tb.sample tb.damage tb.significance tb.tax_id tb.N_reads tb.rho_Ac tb.phi tb.N_sum_total tb.N_min], [:Label, :damage, :significance, :accession, :N_reads, :rho_Ac, :phi, :N_sum_total, :N_min]) |> DataFrames.transform(_, :Label => ByRow(x -> "HB_"*last(x,2)) => :Label) 
    elseif startswith(tag,"amaw.lca")
        df = @pipe DataFrame([tb.sample tb.damage tb.significance tb.tax_id tb.tax_name tb.tax_rank tb.tax_path tb.N_reads tb.rho_Ac tb.phi tb.N_sum_total tb.N_min], [:Label, :damage, :significance, :tax_id, :tax_name, :tax_rank, :tax_path, :N_reads, :rho_Ac, :phi, :N_sum_total, :N_min]) |> DataFrames.transform(_, :Label => ByRow(x -> "HB_"*last(x,2)) => :Label)
    elseif startswith(tag,"holi")
        df = @pipe DataFrame([tb.sample tb.D_max tb.lambda_LR tb.tax_id tb.tax_name tb.tax_rank tb.tax_path tb.N_reads tb.rho_Ac tb.phi tb.N_sum_total tb.N_min], [:Label, :damage, :significance, :tax_id, :tax_name, :tax_rank,:tax_path, :N_reads, :rho_Ac, :phi, :N_sum_total, :N_min]) |> DataFrames.transform(_, :Label => ByRow(x -> "HB_"*last(x,2)) => :Label) |> dropmissing
    end
    
    ## data of samples
    mt = load("./metadata/HIBO_library_metadata.csv") |> DataFrame # library metadata
    x = @pipe leftjoin(df, select(mt,[:Label, :Middle_depth, :yrBP, :N_filtered_reads]), on = :Label) |>
        filter(:rho_Ac => r -> abs(r) <= 0.6, _) |>
        filter(:significance => s -> s>=0.3, _) |>
        filter(:N_reads => n -> n > nreads, _) |>
        sort(_, :Middle_depth) |>
        dropmissing(_, [:Middle_depth, :damage])
    if alg==".lca" 
        filter!(row->any(occursin.(ftax,row.tax_path)), x) 
        filter!(:tax_rank => r -> r==frank, x)
    end
    
    ## data of controls y = @pipe leftjoin(df, select(mt,[:Label, :Middle_depth, :yrBP, :sigma, :N_filtered_reads]), on = :Label) |>
    y = @pipe leftjoin(df, select(mt,[:Label, :Middle_depth, :yrBP, :N_filtered_reads]), on = :Label) |>
        filter(:significance => s -> s>=0.001, _) |>
        filter(:N_reads => n -> n > 10, _)|>
        filter(:Middle_depth=>d->ismissing(d),_) 
    if alg==".lca"
        filter!(:tax_rank => r -> r==frank, y)
    end
    
    ## remove contaminant taxa from sample data
    #if pip=="amaw" filter!(:tax_name=>n->!(n in y.tax_name), x) end
    return x, y
end


