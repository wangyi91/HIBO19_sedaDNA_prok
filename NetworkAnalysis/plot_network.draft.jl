
# load through flashweave.jl and key_communities.jl

taxmember = @pipe DataFrame(tax_path=names(otu)[2:end]) |>
        DataFrames.transform(_, :tax_path => ByRow(t -> tax_at_rank.(t, "kingdom")) => :kingdom) |>
        DataFrames.transform(_, :kingdom => ByRow(t -> ifelse(t=="Bacteria",1,ifelse(t=="Archaea",2,3))) => :idx)

member0 = vcat(taxmember.idx,fill(4,sum(netw.meta_variable_mask)));
member=copy(member0)
nodecolor = [colorant"turquoise",colorant"orange",colorant"lightgrey",colorant"cornsilk4",colorant"palevioletred"]; # bacteria, archaea, others, MVs, highlighted
    #nodecolor = [nothing,nothing,colorant"lightgrey",colorant"cornsilk4",colorant"palevioletred"]; # selective coloring
ndfill = nodecolor[member];
ndsize = [1,1,1,2,1][member];
ndlabel = hide_labels_w_str(names(netw), ["__",":\"","period"]);

strokecolor = [colorant"cornflowerblue",colorant"lightgrey"];
    #strokecolor = [nothing,colorant"lightgrey"];
wts = [Graphs.weights(g)[src.(edges(g))[i], dst.(edges(g))[i]] for i in 1:length(edges(g))];
strcolor = strokecolor[(wts.>0).+1];
    # nodes size proportional to their degree
    #nodesize = [Graphs.outdegree(g, v) for v in Graphs.vertices(g)]
    #alphas = nodesize/maximum(nodesize)
    #nodefills= [RGBA(nodecolor[member[i]],1) for i in 1:length(alphas)]

gp = gplot(g, nodefillc=ndfill, nodesize=ndsize, NODESIZE=0.01,
           nodelabel=ndlabel, NODELABELSIZE=2, nodelabeldist=2.5, nodelabelangleoffset=π/4,
           edgestrokec=strcolor, layout=spring_layout);
# no label
gp = gplot(g, nodefillc=ndfill, nodesize=ndsize, NODESIZE=0.01,
           edgestrokec=strcolor, layout=spring_layout);

draw(PDF("./NetworkAnalysis/output/graph$tag.pdf", 16cm, 16cm), gp)

# find cliques and color the nodes differently
cliper=clique_percolation(g, k=3)
save_object("./NetworkAnalysis/output/cliper_k3$(add3).jld2",cliper);

cliq_idx = findall(>(6), length.(clique_percolation(g, k=3)));
#cliq_idx = findall(in(6..10), length.(clique_percolation(g, k=4)));
node_hl = vcat(collect.(clique_percolation(g, k=3)[cliq_idx])...) |> unique |> sort;
member = member0; member[node_hl].=5; # hightlight these nodes with #5 color
ndfill = nodecolor[member]; # update ndfill

