#!/usr/bin/env julia

# This script generates table and heatmap of selected groups
# It requires objects in key_communities.jl as input
function output_hmp_tb_of_groups()
    upper_range=log1p(maximum(otudata.N_reads))
    for i in 1:length(Is)
    
    clqs = collect.(cliper[cliq_idx])[Is[i]]
    clq = reduce(vcat, clqs) |> sort |> unique
    splist = collect(skipmissing(tax_at_rank.(names(netw)[clq], "species")))
    sub = @pipe otudata |> filter(:tax_name => n-> n in splist, _) |> 
                rightjoin(_,DataFrame(Label=unique(otudata.Label), Middle_depth=unique(otudata.Middle_depth)),on=[:Label,:Middle_depth])
    sub0 = @pipe coalesce.(sub, 0)|>sort(_,:Middle_depth)
    wide = unstack(sub0, :Label, :tax_name, :N_reads, fill=0)
    wide_dmg = unstack(sub0, :Label, :tax_name, :damage, fill=0)
    tax_tb = sub[:,[:tax_name,:tax_path]] |> unique
    x=1:size(wide)[2]
    y=35:-1:1
    fig=Figure(size=(1.5*size(wide)[2],250))
    ax=Axis(fig[1,1], title=titles[i]*"\n"*string(nsps[i])*" taxa")
    hidedecorations!(ax)
    CairoMakie.heatmap!(ax, x, y, transpose(Matrix(log1p.(wide[:,2:end]))), colorrange=(0,upper_range))
    save("./NetworkAnalysis/output/heatmap_group$(i)of$(length(Is)).pdf", fig)
    # output table:
    CSV.write("./NetworkAnalysis/output/table_nreads_group$(i)of$(length(Is)).csv", wide)
    CSV.write("./NetworkAnalysis/output/table_dmg_group$(i)of$(length(Is)).csv", wide_dmg)
    CSV.write("./NetworkAnalysis/output/table_tax_group$(i)of$(length(Is)).csv", tax_tb)
    end

    # Create a separate figure for the colorbar
    colorbar_fig = Figure(size=(200, 100))  # Adjust size as needed

    # Add colorbar with the same color range as the heatmap
    Colorbar(colorbar_fig[1,1], colorrange=(0, upper_range), vertical=true, label=rich("log", subscript("e"),"(sequence count)"),height = 85, labelrotation=0)

    # Save the colorbar figure
    save("./NetworkAnalysis/output/legend_of_group_heatmaps.pdf", colorbar_fig)

end
