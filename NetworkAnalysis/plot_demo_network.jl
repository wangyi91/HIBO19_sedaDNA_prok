using Graphs
using GraphPlot
using Colors
using StatsBase
using Compose#import Cairo, Fontconfig

# Step 1: Create a general graph
n_nodes = 100
g = Graph(n_nodes)  # Create a general graph with 50 nodes

# Initialize node degrees to zero
node_degrees = zeros(Int, n_nodes);

# Number of edges to add
n_edges = 650


# Add edges with preferential attachment
for _ in 1:n_edges
    # Choose a node based on degree (nodes with higher degrees have a higher chance of being selected)
    total_degree = sum(node_degrees)
    if total_degree == 0  # Handle the case where no edges are added yet
        u = rand(1:n_nodes)
    else
        #u = rand(1:n_nodes, p=node_degrees ./ total_degree)  # Preferential attachment
        u = sample(1:n_nodes, Weights(node_degrees ./ total_degree), 1)[1]
        #u = sample(1:n_nodes, Weights(2 .^ node_degrees), 1)[1]
    end

    # Randomly select another node, ensuring it's different from the first
    v = rand(1:n_nodes)
    while v == u
        v = rand(1:n_nodes)
    end

    # Add the edge
    add_edge!(g, u, v)

    # Update the node degrees
    node_degrees[u] += 1
    node_degrees[v] += 1
end


communities = clique_percolation(g, k=4) 

# Step 3: Map nodes to their community and assign colors
node_colors = fill(colorant"lightgrey", n_nodes);  # Default color for unassigned nodes

for (i, community) in enumerate(communities)
    color=sample([colorant"green3",colorant"magentar3",colorant"blue3"],1)[1]
    for node in community
        node_colors[node] = color
    end
end

# Step 4: Visualize the graph
gp = gplot(g, nodefillc=[color for color in node_colors]);
draw(PDF("./NetworkAnalysis/output/demo_graph.pdf", 16cm, 16cm), gp)

