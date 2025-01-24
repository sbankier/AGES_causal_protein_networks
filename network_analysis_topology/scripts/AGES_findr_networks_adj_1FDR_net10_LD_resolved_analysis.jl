using DrWatson
@quickactivate "AGES"

using DataFrames
using CSV
using Gadfly
using Cairo
using LinearAlgebra
using SparseArrays
using Graphs
using Test

using BioFindr

import ColorSchemes.seaborn_colorblind as cs

include(srcdir("motifs.jl"))

# network name
network_name = "AGES_findr_networks_adj_1FDR_net10_LD_resolved"

# Load the network file
network_file = datadir("networks", network_name * ".tsv")
net = CSV.read(network_file, DataFrame)

# theme for the plots
my_theme = Theme(
    major_label_font_size=14pt,
    minor_label_font_size=12pt,
    key_title_font_size=14pt,
    key_label_font_size=12pt,
    point_size=4pt);
Gadfly.push_theme(my_theme);

####################################
### Figures for complete network ###
####################################

# compute the outdegree of each A protein using split-apply-combine
outdegree = combine(groupby(net, :Protein_name_A), nrow)


# write outdegree to file
CSV.write(datadir("networks", network_name * "_nodes_outdegree_all.tsv"), outdegree, delim="\t")


# create heights between 0 and 40 for label annotation
outdegree.height = zeros(nrow(outdegree));
outdegree.height[1:10] = [14, 37, 33, 27, 23, 0, 0, 0, 0, 0]; 

# plot the outdegree distribution
pout = plot(outdegree, x=:nrow, Geom.histogram(bincount=10), 
    Scale.x_log10, color=[cs[1]],
    Coord.Cartesian(xmin=1, xmax=3.6),  
    Guide.xlabel("Outdegree"), 
    Guide.ylabel("Number of A proteins"))

# add the label annotation
push!(pout, layer(outdegree[1:5,:], x=:nrow, y=:height, label=:Protein_name_A, order=1,
Geom.hair, Geom.label(position=:left), color=[colorant"black"]))

# save the plot
draw(SVG(plotsdir(network_name * "_outdegree.svg"), 12cm, 9cm), pout)

# compute the indegree of each B protein using split-apply-combine and sort in reverse order
indegree = combine(groupby(net, :Protein_name_B), nrow)
sort!(indegree, :nrow, rev=true)

indegree.height = zeros(nrow(indegree));
indegree.height[1:11] = [3750, 3250, 2750, 2100, 1900, 1500, 1100, 900, 700, 500, 300]; 

pin = plot(indegree, x=:nrow, Geom.histogram(bincount=10),
    Coord.Cartesian(xmin=1, xmax=45, ymax=4100), color=[cs[1]],
    Guide.xlabel("Indegree"), 
    Guide.ylabel("Number of B proteins"))

# add the label annotation
push!(pin, layer(indegree[1:11,:], x=:nrow, y=:height, 
    label=:Protein_name_B, order=1,
    Geom.hair, Geom.label(position=:left), color=[colorant"black"]))

# save the plot
draw(SVG(plotsdir(network_name * "_indegree.svg"), 12cm, 9cm), pin)

##################################
### Figures for A -> A network ###
##################################

# optionally remove all B proteins that are not in the A protein list
u_Protein_name_A = unique(string.(net.Protein_name_A))
tf = [x in u_Protein_name_A for x in string.(net.Protein_name_B)]
netA = net[tf, :]


# compute the outdegree of each A protein using split-apply-combine
outdegreeA = combine(groupby(netA, :Protein_name_A), nrow)
sort!(outdegreeA, :nrow, rev=true)

# write outdegree to file
CSV.write(datadir("networks", network_name * "_nodes_outdegree_Aonly.tsv"), outdegreeA, delim="\t")

# create heights between 0 and 40 for label annotation
outdegreeA.height = zeros(nrow(outdegreeA));
outdegreeA.height[1:10] = [6, 18, 16, 14, 12, 0, 0, 0, 0, 0]; 

# plot the outdegree distribution
poutA = plot(outdegreeA, x=:nrow, Geom.histogram(bincount=10), 
    Scale.x_log10, color=[cs[2]],
    Coord.Cartesian(xmin=1, xmax=2.1, ymax=21),
    Guide.xticks(ticks=[1,2]),
    Guide.yticks(ticks=[0,5,10,15,20]),
    Guide.xlabel("Outdegree"), 
    Guide.ylabel("Number of A proteins"))

# add the label annotation
push!(poutA, layer(outdegreeA[1:5,:], x=:nrow, y=:height, label=:Protein_name_A, order=1,
Geom.hair, Geom.label(position=:left), color=[colorant"black"]))

# save the plot
draw(SVG(plotsdir(network_name * "_outdegree_Aonly.svg"), 12cm, 9cm), poutA)

# compute the indegree of each B protein using split-apply-combine and sort in reverse order
indegreeA = combine(groupby(netA, :Protein_name_B), nrow)
sort!(indegreeA, :nrow, rev=true)

indegreeA.height = zeros(nrow(indegreeA));
indegreeA.height[1:11] = [76, 69, 62, 57, 50, 43, 37, 31, 27, 23, 0]; 

pinA = plot(indegreeA, x=:nrow, Geom.histogram(bincount=10),
    Coord.Cartesian(xmin=1, xmax=45, ymax=80), 
    color=[cs[2]], 
    Guide.xlabel("Indegree"), 
    Guide.ylabel("Number of B proteins"))

# add the label annotation
push!(pinA, layer(indegreeA[1:10,:], x=:nrow, y=:height, label=:Protein_name_B,
Geom.hair, Geom.label(position=:left), color=[colorant"black"]))

# save the plot
draw(SVG(plotsdir(network_name * "_indegree_Aonly.svg"), 12cm, 9cm), pinA)

#######################################
### Graph analysis complete network ###
#######################################

# create list of unique proteins and put them  in a dataframe
nodes = union(string.(net.Protein_name_A), string.(net.Protein_name_B))
nodes = DataFrame(Protein_name=nodes)

# outerjoin nodes with outdegree 
leftjoin!(nodes, outdegree, on=:Protein_name=>:Protein_name_A)

# remove height column and set missing values to 0 and sort in reverse order
select!(nodes, Not(:height))
nodes.nrow[ismissing.(nodes.nrow)] .= 0;
sort!(nodes, :nrow, rev=true)

# we can now remove the number of targets column
select!(nodes, Not(:nrow))

# add a column ID with numerical node ids
nodes.ID = 1:nrow(nodes)

# select the protein name columns from net
net_edges = select(net, [:Protein_name_A, :Protein_name_B, :p2p5])

# in the net_edges dataframe, add columns with the ID of the A and B proteins
leftjoin!(net_edges, nodes, on=:Protein_name_A=>:Protein_name)
rename!(net_edges, :ID => :ID_A)
leftjoin!(net_edges, nodes, on=:Protein_name_B=>:Protein_name)
rename!(net_edges, :ID => :ID_B)

# create a sparse adjacency matrix from the ID columns
A = sparse(net_edges.ID_A, net_edges.ID_B, 1, nrow(nodes), nrow(nodes))

# create a directed graph from the adjacency matrix
G = DiGraph(A)

# double check that the node ids are still sorted by outdegree
@test issorted(Graphs.outdegree(G), rev=true)

# get the size of the largest weakly connected component
cc = weakly_connected_components(G)
resilience_all = DataFrame(hubs_removed=0, max_component=length(cc[1]))

# remove the top k nodes with the highest outdegree
K = 185
for k in 1:K
    # remove vertex k
    rem_vertex!(G, k)
    # get the size of the largest weakly connected component
    cc = weakly_connected_components(G)
    push!(resilience_all, (hubs_removed=k, max_component=length(cc[1])))
end

# add a column with the relative component size to the total number of nodes
resilience_all.max_component_perc = resilience_all.max_component ./ nrow(nodes)

####################################
### Graph analysis Aonly network ###
####################################

# create nodes and edges with only A proteins
tf = [x in u_Protein_name_A for x in nodes.Protein_name]
nodesA = nodes[tf, :]
tf = [x in u_Protein_name_A for x in string.(net_edges.Protein_name_B)]
net_edgesA = net_edges[tf, 1:3]


# outerjoin nodes with outdegree 
leftjoin!(nodesA, outdegreeA[:,1:2], on=:Protein_name=>:Protein_name_A)

# remove height column and set missing values to 0 and sort in reverse order
#select!(nodes, Not(:height))
nodesA.nrow[ismissing.(nodesA.nrow)] .= 0;
sort!(nodesA, :nrow, rev=true)

# we can now remove the number of targets column
select!(nodesA, Not(:nrow))

# add a column ID with numerical node ids
nodesA.ID = 1:nrow(nodesA)

# in the net_edges dataframe, add columns with the ID of the A and B proteins
leftjoin!(net_edgesA, nodesA, on=:Protein_name_A=>:Protein_name)
rename!(net_edgesA, :ID => :ID_A)
leftjoin!(net_edgesA, nodesA, on=:Protein_name_B=>:Protein_name)
rename!(net_edgesA, :ID => :ID_B)

# repeat resilience analysis for Aonly network
A = sparse(net_edgesA.ID_A, net_edgesA.ID_B, 1, nrow(nodesA), nrow(nodesA))
G = DiGraph(A)

# double check that the node ids are still sorted by outdegree
@test issorted(Graphs.outdegree(G), rev=true)

# get the size of the largest weakly connected component
cc = weakly_connected_components(G)
resilience_Aonly = DataFrame(hubs_removed=0, max_component=length(cc[1]))

# remove the top k nodes with the highest outdegree
K = 185
for k in 1:K
    # remove vertex k
    rem_vertex!(G, k)
    # get the size of the largest weakly connected component
    cc = weakly_connected_components(G)
    push!(resilience_Aonly, (hubs_removed=k, max_component=length(cc[1])))
end

# add a column with the relative component size
resilience_Aonly.max_component_perc = resilience_Aonly.max_component ./ nrow(nodesA)

######
### Resilience plot 

# merge the resilience results in a single dataframel
resilience_all.Network = fill("All", nrow(resilience_all))
resilience_Aonly.Network = fill("A only", nrow(resilience_Aonly))
resilience = vcat(resilience_all, resilience_Aonly)

# plot the resilience of the network

pcomp = plot(resilience[resilience.hubs_removed .<= 10, :], 
    x=:hubs_removed, y=:max_component_perc, color=:Network,
    Geom.line, Geom.point,
    Scale.color_discrete_manual(cs[1], cs[2]),
    #Geom.bar(position=:dodge),
    Guide.xlabel("Top hubs removed"), 
    Guide.ylabel("Max. component size (%)"))

# save the plot
draw(SVG(plotsdir(network_name * "_max_weak_component.svg"), 12cm, 9cm), pcomp)

# reset the network for the level analysis
A = sparse(net_edgesA.ID_A, net_edgesA.ID_B, 1, nrow(nodesA), nrow(nodesA))
G = DiGraph(A)

# get level information for the nodes
root_nodes = findall(Graphs.indegree(G) .== 0)
shp = dijkstra_shortest_paths(G, root_nodes)
# add level column to nodesA
nodesA.Level_noDAG = zeros(Int,nrow(nodesA))
for k in nodesA.ID
    nodesA.Level_noDAG[k] = shp.dists[k] 
end

# add outdegree column to nodesA
nodesA.outdegree = Graphs.outdegree(G)

# save nodesA to a file
CSV.write(datadir("networks", network_name * "_Aonly_nodes_all.tsv"), nodesA, delim="\t")

# add level difference column to net_edgesA
net_edgesA.Level_diff_noDAG = [nodesA.Level_noDAG[v] - nodesA.Level_noDAG[u] for (u,v) in zip(net_edgesA.ID_A, net_edgesA.ID_B)]

# find which edges are removed to make a DAG
rename!(net_edgesA, :Protein_name_A=>:Source, :Protein_name_B=>:Target, :p2p5 => :Probability)
sort!(net_edgesA, :Probability, rev=true)
BioFindr.globalfdr!(net_edgesA)
BioFindr.dagfindr!(net_edgesA)

# save net_edgesA to a file
CSV.write(datadir("networks", network_name * "_Aonly_all.tsv"), net_edgesA, delim="\t")

####################
### Motif counts ###
####################

A = sparse(net_edges.ID_A, net_edges.ID_B, 1, nrow(nodes), nrow(nodes))

# double check that A has no self-loops
@test sum(diag(A))==0

# Two-component feedback loops
fbl2 = tr(A*A)/2

# Three-component feedback and feedforward loops
fbl3 = round(Int,tr(A*A*A)/3)
ffl3 = tr(A*A*A')

A2 = sparse(net_edgesA.ID_A, net_edgesA.ID_B, 1, nrow(nodesA), nrow(nodesA))

# double check that A has no self-loops
@test sum(diag(A2))==0

# Two-component feedback loops
fbl2A = tr(A2*A2)/2

# Three-component feedback and feedforward loops
fbl3A = round(Int,tr(A2*A2*A2)/3)
ffl3A = tr(A2*A2*A2')

####################
### DAG analysis ###
####################

# create dataframe for the DAG analysis
netA_biofindr = select(netA, [:Protein_name_A, :Protein_name_B, :p2p5])
rename!(netA_biofindr, :Protein_name_A=>:Source, :Protein_name_B=>:Target, :p2p5 => :Probability)
sort!(netA_biofindr, :Probability, rev=true)

# compute q-values
BioFindr.globalfdr!(netA_biofindr)

# run dagfindr method
G, node2idx = BioFindr.dagfindr!(netA_biofindr)

# convert node2idx to a dataframe
nodesA_biofindr = DataFrame(Protein_name=[], idx=[])
for (k,v) in node2idx
    push!(nodesA_biofindr, (k,v))
end

# count motifs
tf = netA_biofindr.inDAG_greedy_edges
A3 = sparse(netA_biofindr.Source_idx[tf], netA_biofindr.Target_idx[tf], 1, nrow(nodesA_biofindr), nrow(nodesA_biofindr))

tr(A3*A3)/2
round(Int,tr(A3*A3*A3)/3)
tr(A3*A3*A3')

df_motif = createtensor3(A3,A3,A3')
df_motif_name = DataFrame(
    :I => string.(nodesA_biofindr[df_motif.I,1]),
    :J => string.(nodesA_biofindr[df_motif.J,1]),
    :K => string.(nodesA_biofindr[df_motif.K,1]),
)

# write to file
CSV.write(datadir("networks", network_name * "_Aonly_DAG_FFLs.tsv"), df_motif_name, delim="\t");

# add columns for the source and target level, where level is defined as the shortest path length to a root node
root_nodes = findall(Graphs.indegree(G) .== 0)
shp = dijkstra_shortest_paths(G, root_nodes)

# add the level of the source and target nodes
netA_biofindr.Source_level = [shp.dists[v] for v in netA_biofindr.Source_idx]
netA_biofindr.Target_level = [shp.dists[v] for v in netA_biofindr.Target_idx]
nodesA_biofindr.Level = [shp.dists[v] for v in nodesA_biofindr.idx]

# add a column with the level difference
netA_biofindr.Level_diff = netA_biofindr.Target_level .- netA_biofindr.Source_level

# save netA_biofindr  and nodesA_biofindr to a file
CSV.write(datadir("networks", network_name * "_Aonly_with_DAG.tsv"), netA_biofindr, delim="\t")
CSV.write(datadir("networks", network_name * "_Aonly_nodes.tsv"), nodesA_biofindr, delim="\t")

############################################
### Analysis of transitively reduced DAG ###
############################################

Gred = transitivereduction(G)

# add presence or absence of edges in the reduced graph to netA_biofindr
netA_biofindr.Edge_reduced = [has_edge(Gred, u, v) for (u,v) in zip(netA_biofindr.Source_idx, netA_biofindr.Target_idx)]