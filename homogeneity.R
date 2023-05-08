# Function to calculate homogeneity of the resultant community structure of the graph. Run this file before running the main CDCA file.

homogeneity_D <- function(comm,corr_graph_D,corr_D1){

  total_sim = 0
  
  for (j in 1:length(communities(comm))) {
    comm_graph <- induced_subgraph(corr_graph_D,unlist(communities(comm)[j]))
    sim_clus = 0
    clus_node <- unlist(communities(comm)[j])
    str_comm <- as.data.frame(strength(comm_graph))
    centroid = rownames(str_comm)[which(str_comm == max(str_comm[,1]))]
    
    if(length(clus_node)>1){
      for(i in 1:length(clus_node)){
        sim_val <- corr_graph_D[centroid,clus_node[i]]
        sim_clus = sim_clus + sim_val
      }
      total_sim <- total_sim + (1/length(clus_node))*sim_clus
    }
  }
  homogen <- (1/length(communities(comm)))*total_sim
  return(homogen)
}
