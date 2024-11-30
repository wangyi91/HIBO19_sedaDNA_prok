#!/usr/bin/env julia

# This script generates table and heatmap of selected groups
# It requires as objects in key_communities.jl as input, therefore should only be run after it.
function output_tb_hmp_of_groups()
for (ai,i) in enumerate(display_order)
    clqs = collect.(cliper[cliq_idx])[Is[i]]
    clq = reduce(vcat, clqs) |> sort |> unique
    splist = collect(skipmissing(tax_at_rank.(names(netw)[clq], "species")))
    sub = @pipe otudata |> filter(:tax_name => n-> n in splist, _) |> rightjoin(_,DataFrame(Label=unique(otudata.Label)),on=:Label)
    sub=coalesce.(sub, 0)
    wide = unstack(sub, :Label, :tax_name, :N_reads, fill=0)
    x=1:size(wide)[2]
    y=35:-1:1
    fig=Figure(size=(600,450))
    ax=Axis(fig[1,1], title=titles[i]*"\n"*string(nsps[i])*" taxa")
    hidedecorations!(ax)
    CairoMakie.heatmap!(ax, x, y, transpose(Matrix(log1p.(wide[:,2:end]))))
    save("./NetworkAnalysis/output/heatmap_group$(ai)of$(ngroups).pdf", fig)
    # output table:
    CSV.write("./NetworkAnalysis/output/table_group$(ai)of$(ngroups).csv",wide)
end

end
