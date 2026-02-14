function generate_plot()
npanel = length(Is)
depths=otudata.Middle_depth |> sort |> unique;

fig = Figure(size=(220*(npanel+0.8)+100,1000));
axs = [Axis(fig[2,i], limits=(0,0.33, -depths[end]-10, -depths[1]+200), yticks=(-depths, string.(round.(Int, depths))), xticklabelsize=13, yticklabelsize=13,
            xlabel="DNA damage", ylabel="Depth (cm)", yticklabelcolor=:white,
            xlabelpadding=0, ylabelpadding=0) for i in 1:npanel]; # on axs dmg-depth is plotted, plus α and phylogenetic diversity
hidedecorations!.(axs, label=false, ticklabels = false, grid=false);
hideydecorations!.(axs[2:end], ticklabels = false,grid = false);
linkyaxes!(axs...);

bxs = [Axis(fig[1,i], limits=(0,4,nothing,nothing), xlabel="Phylogenetic distance", ylabel="Network connections", 
            xticklabelsize=13, yticklabelsize=13, xlabelpadding=0, ylabelpadding=0, 
            title=titles[i]*"\n"*string(nsps[i])*" taxa", titlefont = :regular, titlesize = 15) for i in 1:length(Is)]; # on bxs hist of phylo distance is plotted
hidedecorations!.(bxs, label=false, ticklabels = false);
hideydecorations!.(bxs[2:end], ticklabels = false);

cxs = [Axis(fig[2,i], limits=(-100,230,-depths[end]-10,-depths[1]+200), yticks=-depths,
            xlabel="Species richness", ylabel="", xaxisposition=:top, xlabelsize=12,
            xlabelpadding=0, ylabelpadding=0) for i in 1:npanel]; # on cxs species richness are plotted
hideydecorations!.(cxs); hidexdecorations!.(cxs, ticklabels=true);

rowsize!(fig.layout, 2, Relative(4/5));
colgap!(fig.layout, 10); rowgap!(fig.layout, 10);
for i in 1:length(Is)
    colsize!(fig.layout, i, Relative(1/7))
end

for i in 1:length(Is)
    clqs = collect.(cliper[cliq_idx])[Is[i]]
    clq = reduce(vcat, clqs) |> sort |> unique
    splist = collect(skipmissing(tax_at_rank.(names(netw)[clq], "species")))
    sub = @pipe otudata |> filter(:tax_name => n-> n in splist, _)

    # plot dmg histogram by depth
    for dep in depths
        dt = @pipe filter(:Middle_depth => d->d==dep, sub); if nrow(dt)==0 continue end
        CairoMakie.hist!(axs[i], Float64.(dt.damage); bins=0:0.005:0.4, normalization=:none, 
                         offset=-dep, weights=fill(5,size(dt)[1]), color=(:deepskyblue4, 0.5))
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
    CairoMakie.hist!(bxs[i], a; bins=0:0.1:5, normalization=:none, color=(:palevioletred3, 1.0))
end

# --- scale bar ---
scale_ntaxa = 10
scale_height = 5*scale_ntaxa

# position (top-right corner, tweak if needed)
xpos = 0.22
y0   = -depths[1] + 130   # small offset from top

lines!(axs[7],
       [xpos, xpos],
       [y0, y0 + scale_height],
       color = (:deepskyblue4, 0.5),
       linewidth = 5)

text!(axs[7],
      xpos + 0.01,
      y0 + scale_height / 2,
      text = "$scale_ntaxa taxa",
      align = (:left, :center),
      fontsize = 11)


save("./NetworkAnalysis/output/hist_dmg_groups.pdf", fig);
end

