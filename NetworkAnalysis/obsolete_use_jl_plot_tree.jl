# use after tree_communities.jl
# alternatively use julia to create a less nice tree
using Plots
colors = Dict(titles[1]=>cgrad(:tab10,categorical=true)[1], # blue
              titles[2]=>cgrad(:tab10,categorical=true)[3], # green
              titles[3]=>cgrad(:tab10,categorical=true)[5], # purple
              titles[4]=>cgrad(:tab10,categorical=true)[2], # orange
              titles[5]=>cgrad(:tab10,categorical=true)[4], # red
              titles[6]=>cgrad(:tab10,categorical=true)[6], # brown
              ""=>:white, "rest"=>:white)
alphas= Dict([titles[i]=>1 for i in 1:length(Is)]); alphas[""]=0;alphas["rest"]=0;
getindex.(Ref(colors), df.nodegroup)



# plot tree
gr(margins = 4Plots.cm, size=(2000, 20*length(leafnames)+200));
Plots.plot(tree, treetype=:fan, showtips=false, legend=false,
           line_z = trait, linecolor=cgrad(clrscheme,rev=true), linewidth=0.5*ifelse(ftax=="Bacteria",1,2),
           markercolor=getindex.(Ref(colors), df.nodegroup), markeralpha=getindex.(Ref(alphas), df.nodegroup), markershape=:circle, markersize=2, markerstrokewidth=0);
savefig("./MicrobeProfiling/output/tree_comm_$(ftax).pdf");

# plot tree legends
using CairoMakie
fig=Figure(resolution=(1500,1000));
hidedecorations!(Axis(fig[1,1]));
hidedecorations!(Axis(fig[2,1]))
group_color = [MarkerElement(marker=:circle, color = clr, strokecolor = :transparent) for clr in collect(values(colors))];
Legend(fig[1,1], group_color, collect(keys(colors)), "Groups (tree leaves)", tellheight = false,
              patchsize=(13,13), labelsize=13)

phyla = filter(!ismissing,df.phylum);
phylaclrs = palette(clrscheme, 1:length(phyla), rev=true)
phyla_color = [PolyElement(color = clr, strokecolor = :transparent, points = Point2f[(0,0.2),(1,0.2),(1,0.6),(0,0.6)]) for clr in phylaclrs];
Legend(fig[2,1], phyla_color, phyla, ftax*rank*" (tree branches)", nbanks=2, patchsize=(13,13), labelsize=13)
save("./MicrobeProfiling/output/tree_comm_$(ftax)_legends.pdf",fig)

