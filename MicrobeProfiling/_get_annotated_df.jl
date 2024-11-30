#!/usr/bin/env julia
using JLD2
using DataFrames
using Pipe
using CSV
using DataFramesMeta

function tax_at_rank(path::String, rank::String)
    if occursin("__",path)
        regex= Regex("(?<=$(rank[1]*"__"))(.*?)(?=;|:|\$)") # match between e.g. "s__" and ; or : or end-of-line
        m1 = match(regex, path)
        if isnothing(m1) return(missing)
        elseif rank=="species" return(m1.match)
        else return(split(m1.match,"_")[1]) end
    elseif occursin("$rank",path)
        str=@pipe split(path,":\"$rank")[1] |> split(_,":")[end]
        return(str)
    else
        return(missing)
    end
end
#=
function tax_at_rank(path::String, rank::String)
    if occursin(";",path) # tax_path directly from gtdb taxonomy
        regex= Regex("(?<=$(rank[1]*"__"))(.*?)(?=;|:|\$)") # match between "s__" and ; or : or end-of-line
    else    
        regex = Regex("(?<=:)([^\"]+)(?=(:\"$rank))")
        m1 = match(regex, path)
        if isnothing(m1) return(missing)
        elseif occursin("__", m1.match)
            regex = Regex("(?<=__)([^\"]+)(?=(:\"$rank))")
            return(match(regex, path).match)
        else return(match(regex, path).match)
        end
    end
end
=#
function clean_tax(tax::String)
    if occursin("__", tax)
        out=match(r"(?<=__)(.*)", tax).match
    else out=tax
    end
    return(out)
end

function join_contig_taxassign_to_df(lib::String, pip::String, frank::String)
    sp = @pipe CSV.File("../MAG/diamond/output/$(lib).gtdb_$(first(frank,2)).txt", header=0, delim='\t') |> DataFrame |>
               rename(_, [:contigname, :taxpath]) |>
               DataFrames.transform(_, :taxpath => ByRow(t -> split.(t, ';')) => [:db, :domain, :phylum, :class,:order,:family,:genus,:species,:i]) |>
               select(_, Not([:db, :i, :taxpath]), :taxpath);

    df = load_object("./DmgMixtureModel/data/obs_$(lib)_$(pip)_$(frank).jld2");

    if frank=="species"
        sp = @pipe sp |> 
                   DataFrames.transform(_, :taxpath => ByRow(t -> split.(t, ';')) => [:db, :domain, :phylum, :class,:order,:family,:genus,:species,:i]) |>
                   select(_, Not([:db, :i, :taxpath]), :taxpath);
    elseif frank=="genus"
        sp = @pipe sp |>
                   DataFrames.transform(_, :taxpath => ByRow(t -> split.(t, ';')) => [:db, :domain, :phylum, :class,:order,:family,:genus,:i]) |>
                   select(_, Not([:db, :i, :taxpath]), :taxpath);
    end

    jnt = @pipe df |> filter(:tax_path => p -> occursin(ftax,p), _) |> 
        DataFrames.transform(_, :tax_path => ByRow(t -> tax_at_rank.(t, frank)) => Symbol(frank)) |>
        select(_, Not([:A, :phi, :N_sum_total, :N_min, :N_x1_forward, :N_x1_reverse])) |> 
        leftjoin(_, sp, on = Symbol(frank))

    jnt = @pipe jnt |> DataFrames.transform(_, :contigname => ByRow(n->ifelse(ismissing(n), 0,1) ) => :i) |> 
    combine(groupby(_, [:Label,:damage,:significance,Symbol(frank),:tax_path,:N_reads,:Middle_depth,:yrBP, :sigma]), :i .=> sum .=> "N_contigs")
    return(jnt)
end


function get_annotated_df(lib::String, pip::String, frank::String, DB::String="")
    sp = CSV.File("../MAG/diamond/output/$(lib).gtdb_$(first(frank,2)).txt", header=0, delim='\t') |> DataFrame;
    anno_name = CSV.File("../MAG/diamond/output/$(lib).anno.name.txt", header=0, delim='\t') |> DataFrame;
    anno_path = CSV.File("../MAG/diamond/output/$(lib).anno.path.txt", header=0, delim='\t') |> DataFrame;

    df = load_object("./DmgMixtureModel/data/obs_$(lib)_$(pip)_$(frank).jld2");

    # cleanup and reformat annotations
    apath = @pipe anno_path |> rename(_, [:contigname, :annopath]) |>
            DataFrames.transform(_, :annopath => ByRow(c -> split.(c, ';'; limit=3)) => [:annodb,:a1,:a2]);
    
    if frank=="species"
        aname = @pipe anno_name |> leftjoin(_, sp, on = :Column1, makeunique=true) |> dropmissing |>
                filter(:Column2 => c -> !occursin("Not ", c), _) |>
                rename(_, [:contigname, :annoname, :taxpath]) |>
                DataFrames.transform(_, :taxpath => ByRow(t -> split.(t, ';')) => [:db, :domain, :phylum, :class,:order,:family,:genus,:species,:i]) |> 
                select(_, Not([:db, :i, :taxpath]), :taxpath);

    elseif frank=="genus"
        aname = @pipe anno_name |> leftjoin(_, sp, on = :Column1, makeunique=true) |> dropmissing |>
            filter(:Column2 => c -> !occursin("Not ", c), _) |>
            rename(_, [:contigname, :annoname, :taxpath]) |>
            DataFrames.transform(_, :taxpath => ByRow(t -> split.(t, ';')) => [:db, :domain, :phylum, :class,:order,:family,:genus,:i]) |> 
            select(_, Not([:db, :i, :taxpath]), :taxpath);
    end

    # join annotation names
    jnt = @pipe df |> filter(:tax_path => p -> occursin(ftax,p), _) |> 
            DataFrames.transform(_, :tax_path => ByRow(t -> tax_at_rank.(t, frank)) => Symbol(frank)) |>
            select(_, Not([:A, :phi, :N_sum_total, :N_min, :N_x1_forward, :N_x1_reverse])) |> 
            leftjoin(_, aname, on = Symbol(frank)) |> dropmissing |>
            sort(_, :tax_id);

    # join annotation paths        
    jnt = @pipe jnt |> leftjoin(_, apath, on = :contigname) |> filter([:annoname,:annopath] => (n,p) -> occursin(n, p), _) |> 
            unique(_) |> sort(_,[:tax_name, :annopath]) |> 
            select(_,:Label, :tax_name,:N_reads, :annoname,:annopath,:);

    # filter to use only KEGG results
    if DB!=""
        dt = filter!(:annodb => d -> occursin(DB, d), jnt);
        DataFrames.transform!(dt, :a2 => ByRow(t -> split.(t, ';')) => [:a2,:a3,:a4,:i]);
        dt = @pipe select(dt, :a2, :a3, :a4,:) |> sort(_, [:a2,:a3,:a4]);
        return(dt)
    else
        return(jnt)
    end
end

