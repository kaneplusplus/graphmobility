library(Matrix)

# Get the mobility graphs 1 for each hour.
m = readRDS("ms.rds")

# Get the mobility graph correpsonding to 12 noon.
m1 = m[[12]]

# How many towers?
dim(m1)

# What does it look like?
m1[1:10, 1:10]

# Louvain clustering... symmetrize first.
mm = (m1 + t(m1)) / 2

# Create the graph.graph
library(igraph)
g = graph_from_adjacency_matrix(mm, diag = TRUE, weighted = TRUE, mode = "undirected")

# perform clustering.
lc = cluster_louvain(g)

# Get the number of vertices in each cluster. Show the first 50.
lc$membership |> table() |> sort(decreasing = TRUE) |> head(50)

# Louvain tries to optimize on modularity.
# names(lc)


lc$membership |> table() |> sort(decreasing = TRUE) |> head(50)
