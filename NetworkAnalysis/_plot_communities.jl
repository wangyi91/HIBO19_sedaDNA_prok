function generate_plot()
npanel = length(Is)
depths=otudata.Middle_depth |> sort |> unique;

fig = Figure(size=(220*(npanel+0.8)+100,1000));
axs = [Axis(fig[2,i], limits=(0,0.33, -depths[end]-10, -depths[1]+200), yticks=(-depths, string.(round.(Int, depths))), xticklabelsize=13, yticklabelsize=13,
            xlabel="DNA damage", ylabel="Depth (cm)",
            xlabelpadding=0, ylabelpadding=0) for i in 1:npanel]; # on axs dmg-depth is plotted, plus α and phylogenetic diversity
hidedecorations!.(axs, label=false, ticklabels = false, grid=false);
hideydecorations!.(axs[2:end], grid = false);
#hidexdecorations!.(axs[2:end], label=false, ticklabels=false, grid=false);
linkyaxes!(axs...);

bxs = [Axis(fig[1,i], limits=(0,4,nothing,nothing), xlabel="Phylogenetic distance", ylabel="No. connections", 
            xticklabelsize=13, yticklabelsize=13, xlabelpadding=0, ylabelpadding=0, 
            title=titles[i]*"\n"*string(nsps[i])*" taxa", titlefont = :regular, titlesize = 15) for i in 1:length(Is)]; # on bxs hist of phylo distance is plotted
hidedecorations!.(bxs, label=false, ticklabels = false);
rowsize!(fig.layout, 2, Relative(4/5));
colgap!(fig.layout, 10); rowgap!(fig.layout, 10);

cxs = [Axis(fig[2,i], limits=(-100,230,-depths[end]-10,-depths[1]+200), yticks=-depths,
            xlabel="Species richness", ylabel="", xaxisposition=:top, xlabelsize=12,
            xlabelpadding=0, ylabelpadding=0) for i in 1:npanel]; # on cxs species richness are plotted
hideydecorations!.(cxs); hidexdecorations!.(cxs, ticklabels=true);

for i in 1:length(Is)
    clqs = collect.(cliper[cliq_idx])[Is[i]]
    clq = reduce(vcat, clqs) |> sort |> unique
    splist = collect(skipmissing(tax_at_rank.(names(netw)[clq], "species")))
    sub = @pipe otudata |> filter(:tax_name => n-> n in splist, _)
    # plot dmg histogram by depth
    for dep in depths
        dt = @pipe filter(:Middle_depth => d->d==dep, sub); if nrow(dt)==0 continue end
        #w = Float64.(sqrt.(dt.N_reads))*0.6
        CairoMakie.hist!(axs[i], Float64.(dt.damage); bins=0:0.005:0.4, normalization=:none, 
                         offset=-dep,weights=fill(5,size(dt)[1]), color=(:teal,0.5))
    end
    # plot species richness
    summ = combine(groupby(sub,:Middle_depth),nrow=>"n_species")
    deps = @pipe filter(d->d>=minimum(sub.Middle_depth), depths) |> filter(d->d<=maximum(sub.Middle_depth), _)
    app = DataFrame(Middle_depth=filter(d->!(d in sub.Middle_depth),deps)); app.n_species.=0;
    @pipe append!(summ,app) |> sort!(_,:Middle_depth)
    CairoMakie.lines!(cxs[i], summ.n_species, summ.Middle_depth*(-1); linewidth=3, color=(:tan1,0.6))
    # plot phylo dist
    deleteat!(clq, clq .> ntaxa);# remove MVs
    wsub = Graphs.weights(g)[clq,clq]; w=vec(wsub);
    psub = phylodist[clq,clq]
    a=vec(psub); a=a[findall(>(0),w)];a=a[findall(in(1e-10..5),a)]
    CairoMakie.hist!(bxs[i], a; bins=0:0.1:5, normalization=:none, color=(:firebrick, 0.6))
end
connect_color = [PolyElement(color = (clr, 0.4), strokecolor = :transparent) for clr in [colorant"firebrick", colorant"royalblue3"]];
bar_color=PolyElement(color = (:teal, 0.5), strokecolor = :transparent, points = Point2f[(0.25,0),(0.55,0),(0.55,1),(0.25,1)]);
save("./NetworkAnalysis/output/hist_dmg_groups.pdf", fig);
end

