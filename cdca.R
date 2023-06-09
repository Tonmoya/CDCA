# compile this file after compiling homogeneity.R

eigen_cen_community <- function(data, dataname){

library(e1071)
library(dplyr)
library(igraph)
library(apeglm)
library(ggplot2)
library(tibble)
library(stringr)
library(ggnetwork)
library(tidyverse)

setwd("CDCA/")  # set the working directing where all the files of the program will be saved
outfolder = "memberships/"    # specify a location where output files will be saved


#preparing the data--------
geneNames=data[,1]
rawData<-data[-1] 

rownames(rawData) <- data[,1]

corr_D1 <- as.matrix(rawData) 
corr_D1 <- abs(corr_D1)      # positive weights matrix
#corr_D1[lower.tri(corr_D1)] <- 0 
corr_D1[is.na(corr_D1)] <- 0
diag(corr_D1) <- 0


# Correlation threshold - find minimum value of correlation
min_corr <- data.frame(matrix(nrow = dim(corr_D1)[1], ncol = 1))
for (m in 1:dim(corr_D1)[1]) {
  min_corr[m,1] <- min(corr_D1[,m][which(corr_D1[,m]>0)])
}
min_corr <- unlist(min_corr)
min_corr[!is.finite(min_corr)] <- 0
max_corr_final <- max(min_corr)   # get maximum of the minimum values

corr_D1[corr_D1<=max_corr_final] <- 0  # threshold max_corr_final


# Create graph from correlation matrix
corr_graph_D <- graph_from_adjacency_matrix(corr_D1, weighted = TRUE, mode = "undirected")
corr_graph_D <- simplify(corr_graph_D)

comp_corr_graph_D <- components(corr_graph_D)  
mod1_D <- modularity(corr_graph_D, membership(comp_corr_graph_D))

corr_D1 <- as.data.frame(corr_D1)

adj_matrix <- as_adj(corr_graph_D, type = "both")
adj_matrix <- as.matrix(adj_matrix)


##---------------- CDCA ALGO ---------


#all_del_index <- 1:length(geneNames)

new_graph <- make_empty_graph(n=dim(corr_D1)[1], directed = FALSE)
comp_new_graph <- components(new_graph)  # membership, csize(cluster sizes), no(no. of clusters)
#mod1_new_graph <- modularity(new_graph, membership(comp_new_graph))
#mod1_temp = mod1_D


# STEP 1 - initialize each node with strength of the node 

deg_D <- degree(corr_graph_D) # DEGREES of all the nodes
strength_D <- strength(corr_graph_D) # weighted degree of all nodes
strength_D_vector <- c(strength_D)


# STEP 2 - Calculate eigen vector centrality of the graph

eigen_cen <- eigen_centrality(corr_graph_D, directed = FALSE)  # EIGEN VECTOR CENTRALITY
centrality <- eigen_cen$vector

centrality_nonzero <- subset(centrality, centrality>0) # remove outliers with 0 centrality
centrality_sort <- as.data.frame(sort(centrality_nonzero, decreasing = TRUE))

all_scores <- cbind(centrality,deg_D,strength_D,visited=0)

# Top N highest values as cluster centers
cluster_centers <- order(centrality, decreasing=TRUE)[1:5]

all_scores[cluster_centers,"visited"] = 1  # all clusters centers are made visited


# STEP 3 - Find SUB-COMMUNITY --------

# check correlations of cluster centers with all other nodes -
# connect node with cluster center with which it has highest correlation -
# initial number of clusters will be number of cluster centers found
# visited nodes to be mark visited - loop until no more unvisited nodes left

for (other_nodes in 1:length(centrality)) {
  
  get_nodes <- data.frame(matrix(ncol = 3)) # empty matrix to save nodes
  #get_nodes <- data.frame(matrix(ncol = 3, nrow =length(cluster_centers)))
  #print(paste0(all_scores[other_nodes,1]))
  
  if(all_scores[other_nodes,"visited"]==0 && all_scores[other_nodes,"deg_D"]!=0){
    all_scores[other_nodes,"visited"] = 1 
   
    for (k in 1:length(cluster_centers)) {
      # get correlation between two vertices
      get_nodes[k,] = c(other_nodes, cluster_centers[k], corr_graph_D[other_nodes,cluster_centers[k]])
      
    }
    
    max_value <- max(get_nodes[,3])
    if(max_value!=0){
      max_node <- get_nodes[which(get_nodes[,3]==max(get_nodes[,3])),]  # get pair of nodes with max correlation
      
      node1 = unlist(max_node[1])
      node2 = unlist(max_node[2])
      
      #print(paste0(node1,"----",node2))
      new_graph[node1,node2] = corr_graph_D[node1,node2]  # SUB-GRAPH
      
    }
    
    
  }
}

#Get modularity of graph after initial sub-communities 
comp_new_graph <- components(new_graph)  
mod1_new_graph <- modularity(new_graph, membership(comp_new_graph))



# STEP 4 - MERGE SUB-COMMUNITIES -------------------

# initialize change in modularity to a max limit
diff_mod = mod1_new_graph - mod1_D
#
org_graph = new_graph
comp_org_graph = comp_new_graph
merge_community = data.frame(matrix(ncol = 3)) # which communities to be merged


for (m in 1:((comp_new_graph$no)-1)) {
  print(paste0("community1 -- ",m))
  merging_nodes = data.frame(matrix(ncol = 3))  # matrix to take maximum merging factor
  
  inter_sub_correlation = data.frame(matrix(ncol = 3)) # matrix to save inter community correlation
  
  com_nodes1 <- communities(comp_new_graph)[m] #get nodes of a community
  sub_graph1 <- induced_subgraph(corr_graph_D, unlist(com_nodes1)) # get the first individual community 
  
  edge_sum1 = sum(E(sub_graph1)$weight)  #sum of correlation edges of own community - first
  avg_sum1 = edge_sum1/length(com_nodes1[[1]])  # average of correlation within community
  
  for (n in (m+1):(comp_new_graph$no)) {
    print(paste0("community2 -- ",n))
    sum_edge_cut = 0
    com_nodes2 <- communities(comp_new_graph)[n]  #get nodes of second community
    sub_graph2 <- induced_subgraph(corr_graph_D, unlist(com_nodes2)) # get the second individual community 
    
    edge_sum2 = sum(E(sub_graph2)$weight)  #sum of edges of own community - second
    avg_sum1 = edge_sum1/length(com_nodes2[[1]])  # average of correlation within community - second
    #merge_cen2 = intersect(unlist(com_nodes2),cluster_centers) # get the center node 1
    
    # sum of cprrelation edges shared by the two clusters
    for (p in 1: length(com_nodes1[[1]])) {
      for (q in 1:length(com_nodes2[[1]])) {
        sum_edge_cut = sum_edge_cut + corr_graph_D[com_nodes1[[1]][p],com_nodes2[[1]][q]]
        new_entry <- c(com_nodes1[[1]][p],com_nodes2[[1]][q],corr_graph_D[com_nodes1[[1]][p],com_nodes2[[1]][q]]) 
        inter_sub_correlation = rbind(inter_sub_correlation, new_entry) 
      }
    }
    
    
    # STEP 5 - Calculate MERGING FACTOR --------------------
    
    merging_factor = (sum_edge_cut/(gsize(sub_graph1)+gsize(sub_graph2)))/((edge_sum1/gsize(sub_graph1)) + (edge_sum2/gsize(sub_graph2)))
    merging_nodes = rbind(merging_nodes, c(m,n,merging_factor)) # community1, community2, value
    
  }
  
  max_merge = max(merging_nodes[,3],na.rm = TRUE)  # maximum merging factor among two communities
  merge_com1 = m
  merge_community = rbind(merge_community,merging_nodes[which(merging_nodes[,3]==max_merge),])  # get the two communities
  
}

sort_merge_community <- merge_community[order(merge_community[,3], decreasing = TRUE),]
sort_merge_community <- cbind(sort_merge_community, visited=0)
sort_merge_community <- na.omit(sort_merge_community)
sort_merge_community <- subset(sort_merge_community, sort_merge_community[,3]!=Inf)


org_graph = new_graph
comp_org_graph = comp_new_graph



for (r in 1:dim(sort_merge_community)[1]) { 
  
  temp_graph = new_graph
  comp_temp_graph = comp_new_graph
  diff_mod1 = diff_mod
  mod_temp = mod1_new_graph
  
  com_nodes1 <- communities(comp_org_graph)[sort_merge_community[r,1]] #get nodes of a community from original
  sub_graph1 <- induced_subgraph(corr_graph_D, unlist(com_nodes1)) # get the individual community 1
  
  com_nodes2 <- communities(comp_org_graph)[sort_merge_community[r,2]] #get nodes of a community from original
  sub_graph2 <- induced_subgraph(corr_graph_D, unlist(com_nodes2)) # get the individual community 2 
  
  
  for (e1 in 1:length(unlist(com_nodes1))) {
    merge_node1 = unlist(com_nodes1)[e1]
    for (e2 in 1:length(unlist(com_nodes2))) {
      merge_node2 = unlist(com_nodes2)[e2]
      
      new_graph[merge_node1, merge_node2] = corr_graph_D[merge_node1, merge_node2]    # CREATING FINAL GRAPH
    }
    
    sort_merge_community[r,"visited"]=1
  }
  
  comp_new_graph <- components(new_graph)  
  mod1_new_graph <- modularity(new_graph, membership(comp_new_graph))
  diff_mod = mod1_new_graph - mod_temp                                   # check difference in  modularity of current and previous graph
  print(paste0("diff_mod---  ",diff_mod,"----mod1_new_graph----",mod1_new_graph))
  
  if(diff_mod < diff_mod1){
    print(paste0(diff_mod," < ",diff_mod1))
    new_graph = temp_graph
    comp_new_graph = comp_temp_graph
    sort_merge_community[r,"visited"]=0   # if difference in modularity is less, merge of community is not accepted
    
    #comp_new_graph <- components(new_graph)  
    mod1_new_graph <- modularity(new_graph, membership(comp_new_graph))
    #diff_mod = mod1_new_graph - mod_temp # check difference in  modularity of current and previous graph
    print(paste0("diff_mod---  ",diff_mod,"----mod1_new_graph----",mod1_new_graph))
  }
}


comp_new_graph$csize


# GET GENENAMES OF COMMUNITIES-----------------

gene_name_list <- data.frame(matrix(ncol = length(communities(comp_new_graph)))) # communities with gene names
  
for (m in 1:length(communities(comp_new_graph))) {
  gene_num_list <- as.data.frame(communities(comp_new_graph)[m])
  for (n in 1:comp_new_graph$csize[m]) {
    gene_num <- gene_num_list[n,1]
    gene_name <- geneNames[gene_num]
    gene_name_list[n,m] = gene_name # add genenames to list
  }
}

member_genes <- data.frame(matrix(ncol = 2)) 
for (p in 1:length(membership(comp_new_graph))) {
  member_genes[p,1] = geneNames[p]
  member_genes[p,2] = membership(comp_new_graph)[p]
}


write.csv(gene_name_list, file = paste0(outfolder,"community_list_",dataname,".csv"), row.names = FALSE)
write.csv(member_genes, file = paste0(outfolder,"membership_",dataname,".csv"), row.names = FALSE)


homogen_result <- homogeneity_D(comp_new_graph,corr_graph_D,corr_D1)  


result_list <- list(comp_new_graph,modularity = mod1_new_graph,homogeneity = homogen_result)


return(result_list)

}
