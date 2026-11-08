# Author: Natasha Puzovic (puzovic@evolbio.mpg.de)
# Date created: 15 Jan 2023
# Last modified: 04 Feb 2023

### This is a script to generate network topologies.
#############  Erdos-Renyi networks

rm(list=ls())

library(igraph)

setwd("../sim_data/")
if(!dir.exists("generated_net_topos")){
  dir.create("generated_net_topos")
}
setwd("generated_net_topos")


################### params

# number of genes in the network
num_nodes <- 40

# number of network topologies to create
num_networks <- 2000

# names of nodes
node_names <- paste("G", 1:num_nodes, sep = "")

# tolerable level of density distance from specified, for generation
tol_dens_deviation <- 0.005

# all densities
ALLnetwork_densities <- c(0.05, seq(from = 0.1, to = 1, by = 0.1))
#ALLnetwork_densities <- c(0.05, 0.1, 0.3, 0.6, 1)
ALLprob_edge <- ALLnetwork_densities - 0.02

# this is for density = 0.1
param_ind = 1
#for(param_ind in 1:(length(ALLnetwork_densities)-1)) {
  
  network_density <- ALLnetwork_densities[param_ind]
  prob_edge <- ALLprob_edge[param_ind]
  print(paste(network_density))
  print(paste(prob_edge))

  if(!dir.exists(paste0("20230212_ER_", num_nodes, "node_", network_density, "dens"))){
    dir.create(paste0("20230212_ER_", num_nodes, "node_", network_density, "dens"))
  }
  setwd(paste0("20230212_ER_", num_nodes, "node_", network_density, "dens"))
  
  
  for(netNum in 1:num_networks)
    {
    # generate a random directed network with one component (no unconnected genes) 
    g <- erdos.renyi.game(n = num_nodes, p.or.m = prob_edge, directed = TRUE, loops = TRUE)
    while(components(g)$no != 1 | 
          ((graph.density(g) > network_density + tol_dens_deviation) | 
           (graph.density(g) < network_density - tol_dens_deviation))) {
        g <- erdos.renyi.game(n = num_nodes, p.or.m = prob_edge, directed = TRUE, loops = TRUE)
    }
    
    #plot(g)
    dens <- graph.density(g)
    print(dens)
    
    # convert to adj matrix & add weights
    rmtrx <- t(as_adjacency_matrix(g, sparse = FALSE))
    #rvals <- rbinom(n = length(rmtrx[rmtrx != 0]), size = 1, prob = 0.5)
    #rvals[rvals == 0] <- -0.5
    #rvals[rvals == 1] <- 0.5
    #rmtrx[rmtrx != 0] <- rvals 
    rmtrx[rmtrx != 0] <- runif(n = length(which(rmtrx != 0) == TRUE), min = -4, max = 4)
    
    # write out
    write.table(rmtrx, 
                file = sprintf("rmtrxFoundingNet_%d.txt", netNum), 
                sep = " ", row.names = FALSE, col.names = FALSE)
    
  }
  #setwd("..")
#}


# network density 1
#network_density <- 1
#prob_edge <- 1
#print(paste(network_density))
#print(paste(prob_edge))
#  
#if(!dir.exists(paste0("2211_ER_", num_nodes, "node_", network_density, "dens"))){
#  dir.create(paste0("2211_ER_", num_nodes, "node_", network_density, "dens"))
#}
#setwd(paste0("2211_ER_", num_nodes, "node_", network_density, "dens"))


#for(netNum in 1:num_networks){
    
  # generate a random directed network with one component (no unconnected genes) 
  #g <- erdos.renyi.game(n = num_nodes, p.or.m = prob_edge, directed = TRUE)
  #while(components(g)$no != 1 | graph.density(g) != network_density) {
  #  g <- erdos.renyi.game(n = num_nodes, p.or.m = prob_edge, directed = TRUE)
  #}
  
  # convert to adj matrix & add weights
  #rmtrx <- t(as_adjacency_matrix(g, sparse = FALSE))
  #rvals <- rbinom(n = length(rmtrx[rmtrx != 0]), size = 1, prob = 0.5)
  #rvals[rvals == 0] <- -0.5
  #rvals[rvals == 1] <- 0.5
  #rmtrx[rmtrx != 0] <- rvals 
  
  # write out
  #write.table(rmtrx, 
  #            file = sprintf("rmtrxFoundingNet_%d.txt", netNum), 
  #            sep = " ", row.names = FALSE, col.names = FALSE)
  
#}


##############################################  end   ############################################