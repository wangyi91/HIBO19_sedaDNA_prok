#!/usr/bin/env julia

# This script generates table and heatmap of selected groups
# It requires as objects in key_communities.jl as input, therefore should only be run after it.
function output_tb_hmp_of_groups()

    for i in 1:length(Is)
    
    clqs = collect.(cliper[cliq_idx])[Is[i]]
    clq = reduce(vcat, clqs) |> sort |> unique
    splist = collect(skipmissing(tax_at_rank.(names(netw)[clq], "species")))
    sub = @pipe otudata |> filter(:tax_name => n-> n in splist, _) |> 
                rightjoin(_,DataFrame(Label=unique(otudata.Label), Middle_depth=unique(otudata.Middle_depth)),on=[:Label,:Middle_depth])
    sub = @pipe coalesce.(sub, 0)|>sort(_,:Middle_depth)
    wide = unstack(sub, :Label, :tax_name, :N_reads, fill=0)
    wide_dmg = unstack(sub, :Label, :tax_name, :damage, fill=0)
    x=1:size(wide)[2]
    y=35:-1:1
    fig=Figure(size=(1.5*size(wide)[2],250))
    ax=Axis(fig[1,1], title=titles[i]*"\n"*string(nsps[i])*" taxa")
    hidedecorations!(ax)
    CairoMakie.heatmap!(ax, x, y, transpose(Matrix(log1p.(wide[:,2:end]))))
    save("./NetworkAnalysis/output/heatmap_group$(i)of$(length(Is)).pdf", fig)
    # output table:
    CSV.write("./NetworkAnalysis/output/table_nreads_group$(i)of$(length(Is)).csv",wide)
    CSV.write("./NetworkAnalysis/output/table_dmg_group$(i)of$(length(Is)).csv",wide)
    end

end
