##########################################################################
# Intervs                                                           #
# This is an R script designed to provide social network analyses for    #
# network-based interventions. As currently written it imports a         #
# a network stored in CSV format and computes network measures and       #
# analyses detailed in the book "Accelerating change:                    #
# Behavior Change Interventions." The code makes extensive use of the    #
# "igraph" package: Csardi (2010).
#                                                                        # 
# Users may want to set the working directory (setwd) path to a location #
# appropriate for thier computing environment                            #
##########################################################################

library(igraph)
library(foreign)

# input Stata file
setwd("c:/misc/accel/") 
att_data <- read.csv("know16_a.csv")
demo3 <- read.csv("know16.csv")
net_data <- demo3[, 1:8]
head(net_data)
dyad3 <-reshape(net_data, direction="long",  varying = 2:8, sep="",
                idvar = "net_data$id", timevar="rank")
dyad2  <- dyad3[!is.na(dyad3$nom),]  
dyad1  <- as.data.frame(dyad2)
vars   <- c("id", "nom")
dyad   <- as.matrix(dyad1[vars])
net_data <- graph_from_data_frame(dyad, vertices = att_data)
net_data  <- simplify(net_data)
adj_mat   <- as.matrix(as_adjacency_matrix(net_data))
net_size  <- gorder(net_data)
setwd("c:/misc/accel/graphs") 

################
# Individuals  #
################
# Leaders/centrals 

# Individual Metrics #
outdeg      <- degree(net_data, mode = "out")       # out degree
indeg       <- degree(net_data, mode = "in")        # in degree
between     <- betweenness(net_data)                # betweeness
outcls      <- closeness(net_data, mode = "out")    # out closeness
outcls      <- outcls*100
incls       <- closeness(net_data, mode = "in")     # in closeness
incls       <- incls*100
keyplay     <- rep.int(0, net_size)  #Key Players calculated in Key Player/UCINET 
inds        <- c(4, 16, 17, 18, 22)
keyplay[inds] <-1
constraintv <- (constraint(net_data))                        # constraint
constraintv[is.nan(constraintv)] = 1
inv_constr  <- (1/(constraintv))  
sym_mat_w <- sna::symmetrize(adj_mat, rule="weak")
sym_mat_s <- sna::symmetrize(adj_mat, rule="strong")
net_data_sym_w <- graph_from_adjacency_matrix(sym_mat_w, mode="undirected")
net_data_sym_s <- graph_from_adjacency_matrix(sym_mat_s, mode="undirected")
bridge_w       <- influenceR::bridging(net_data_sym_w)       # bridging - weak
bridge_w <- bridge_w*100
bridge_s       <- influenceR::bridging(net_data_sym_s)       # bridging - strong
bridge_s <- bridge_s*100
#Brokerage Run everett-valente-brokerage_v2 by Ardestani
broker <- (calculate.EV.brokerage(net_data) + .001)
marginals <-  as.integer(indeg <= 1)
persdens<- lapply(make_ego_graph(net_data, mode="out"), graph.density)
# reciprocity rate is the outdegree of net symmetrized on strong
recip        <- degree(net_data_sym_s, mode = "out") 
ids <-V(net_data)
measures<-cbind(ids,  outdeg, indeg, between, outcls, incls, inv_constr, bridge_w,
            bridge_s, persdens, recip)
write.table(measures, file="measures.csv", append=T, quote=FALSE, sep=",", row.names=F, col.names=F)

#Calculate how many leaders
numol <- as.integer(net_size/7)  
# Select those with highest degree
ids <- get.vertex.attribute(net_data, "id1")
ones<- rep.int(1, net_size)
V(net_data)$indeg    <- indeg
V(net_data)$between  <- between
V(net_data)$incls    <- incls
V(net_data)$keyplay  <- keyplay

layout = layout_nicely(net_data)
for (i in c(2, 4, 6, 8)) {
leaders <- cbind(ids, indeg, ones, between, ones, incls, ones, keyplay, ones) 
leaders <- leaders[rev(order(leaders[,(i)])), ] 
leaders[1:numol, (i+1)] <- 2
leaders <- leaders[sort.list(leaders[,1]), ]

#Plot the network
png(paste("Leader", i, ".png"))
plot(net_data,
     layout = layout,
     edge.arrow.size   =.1, 
     labels.scale      = 1,
     vertex.color      = leaders[,(i+1)],
     vertex.label.font =2, 
     edge.curved=TRUE,
     vertex.label.color = "black"
)
dev.off()
}

#Bridging
V(net_data)$invcons   <- inv_constr
V(net_data)$bridges   <- bridge_s
V(net_data)$broker    <- broker
V(net_data)$margins    <- marginals
for (i in c(2, 4, 6, 8)) {
  bridgers <- cbind(ids, inv_constr, ones, bridge_s, ones, broker, ones, marginals, ones) 
  bridgers <- bridgers[rev(order(bridgers[,(i)])), ] 
  bridgers[1:numol, (i+1)] <- 2
  bridgers <- bridgers[sort.list(bridgers[,1]), ]
  
#Plots for bridgers
png(paste("bridgers", i, ".png"))
plot(net_data,
       layout = layout,
       edge.arrow.size   =.1, 
       labels.scale      = 1,
       vertex.color      = bridgers[,(i+1)],
       vertex.label.font =2, 
       edge.curved=TRUE,
       vertex.label.color = "black"
  )
dev.off()
}

# Calculate exposure & thresholds
# net_data3 <- as.matrix.network(net_data)
# Exposure <-   ((net_data3 %*% node$adopt) / rowSums(net_data3))
# round(Exposure, digits = 4)
# T <- max(node$toa)
# Adopt_mat <- matrix(data = 0, nrow = net_size, ncol = T, byrow = TRUE)
# for(i in 1:net_size) {
#   Adopt_mat[i, node$toa[i]:T] <- rep.int(1, (T-node$toa[i]+1))
# }  
# Expos_mat <- matrix(data = 0, nrow = net_size, ncol = T, byrow = TRUE)
# # DIST <- geodist(net_data, inf.replace=net_size-1, count.paths=TRUE)
# for(j in 1:T) {
#   Expos_mat[,j] <- ((net_data3 %*% Adopt_mat[,j]) / (rowSums(net_data3)+.001))
# }
# round(Expos_mat, digits = 4)
# Thresh <- rep.int(0, net_size)
# for(i in 1:net_size) {
#   Thresh[i] <- Expos_mat[i, node$toa[i]]
# }
# Thresh
# round(Thresh, 3)

library(netdiffuseR)
toa <- att_data$toa
net_data_diff <-as_diffnet(adj_mat, toa)
net_data_exp  <- exposure(net_data_diff)
net_data_thr  <- threshold(net_data_diff)
V(net_data)$thresh <- net_data_thr  
V(net_data)$x <- toa
V(net_data)$y <- net_data_thr
# Plotting threshold adoptions 
png(paste("Thresholds", ".png"))
plot_threshold(net_data_diff, vertex.label = net_data$id1)
dev.off()

################
# Segmentation #
################
# Groups #
girnew  <- cluster_edge_betweenness(net_data_sym_s) # Girve-Newman only for symmetric graphs
optcom  <- optimal.community(net_data)
# devtools::install_github("aslez/concoR")
# concor_net <- cor(adj_mat)
# blks <- concor_hca(list(adj_mat), p = 2)
# From UCINET
concorU  <- t(cbind(1, 5, 1, 5, 7, 3, 3, 5, 4, 5, 5, 3, 6, 4, 4, 1, 8, 6, 6, 6, 7, 4, 5, 2, 1, 2, 3, 1, 3, 8, 7, 8, 3, 5, 1))

#Groups Plots
V(net_data)$girnew   <- girnew$membership
V(net_data)$optcom   <- optcom$membership
V(net_data)$concor   <- concorU
for (i in c(2:4)) {
  groups <- cbind(ids, girnew$membership, optcom$membership, concorU) 
  groups <- groups[sort.list(groups[,1]), ]
  
  #Plots for Groups
  png(paste("Groups", i, ".png"))
  plot(net_data,
       layout = layout,
       edge.arrow.size   =.3, 
       labels.scale      = 1,
       vertex.color      = groups[,(i)],
       vertex.label.font =2, 
       edge.curved=TRUE,
       vertex.label.color = "black"
  )
dev.off()
}

#Create Hiearchal mapping with In-degree as Y-variable 
#Plot with x-axis as TOA and y-axis as Indegree
coords <- as.matrix(data.frame(c(toa), c(indeg)))
png(paste("Hierachy Based on Indegree",  ".png"))
plot(net_data,
     layout = coords,
     edge.arrow.size   =.3, 
     labels.scale      = 1,
     vertex.color      = groups[,(i)],
     vertex.label.font =2, 
     edge.curved=TRUE,
     vertex.label.color = "black"
)
dev.off()
#Create blockmodel on Department in UCINET 
# blocks <- blockmodel(net_data, net_data$toa)
# summary(blocks)

################
#  Induction   #
################
# WOM, RDS & Snowball : Create a sub matrix of randomly selected cases and who they are linked to
# Select vertices that adopt early, grab their neigh edges, color them
ids_early <- which(toa < 3)
WOMlist <- incident_edges(net_data, ids_early, mode = "out")
ecolor <- c("gray", "tomato")[E(net_data) %in% unlist(WOMlist) + 1]
png(paste("Word of Mouth", ".png"))
plot(net_data,
     layout = layout,
     edge.arrow.size   =.1, 
     labels.scale      = 1,
     vertex.color      = "tan",
     vertex.label.font =2, 
     edge.curved=TRUE,
     vertex.label.color = "black",
     edge.color = ecolor,
     vertex.size = 1
)
dev.off()

#Network Outreach - generate a star
for(i in 1:numol) {
  stars <- make_star(6, mode = "out", center = 1)
  png(paste("Network Outreach", ".png"))
    plot(stars, 
       edge.arrow.size   =.3, 
       labels.scale      = 1,
       vertex.color      = "blue",
       vertex.size       = 25,
       vertex.label.font =2, 
       edge.curved=FALSE,
       vertex.label.color = "black"
  )
}
dev.off()
#Leader -learner matching Georges algorithm
ans <-netdiffuseR::mentor_matching(adj_mat, n=numol)
match_ll <- cbind(ans$name, ans$match)
match_net<- graph_from_edgelist(match_ll)
match_net<- simplify(match_net)
png(paste("Mentor Matching - Georges Algorithm", ".png"))
plot(match_net, 
     edge.arrow.size   =.3, 
     labels.scale      = 1,
     vertex.color      = "tan",
     vertex.size       = 12,
     vertex.label.font =2, 
     edge.curved=FALSE,
     vertex.label.color = "black"
)
dev.off()
# Leader - learner matching Valente code
# Select those with highest degree
in_deg2 <- cbind(ids, indeg) 
in_deg2 <- in_deg2[rev(order(in_deg2[,2])), ] 
# Create a sub matrix of the leaders and who they are linked to
lead_nodes <-in_deg2[1:numol, 1]                               # ID Leaders
lead_geos  <-distances(net_data)[ , lead_nodes]      # Get sub-matrix of distance
group_size <- ((net_size - numol)/ numol )                     # Calculate group size
match_assigns <- matrix(0, (group_size), numol)                # Shell for ID numbers
match_counter <- matrix(0, 1 , numol)                          # Shell to count groups
k<- 1                                                          # Row indicator for match_assigns
l<- 1                                                          # col indicator for match_assigns
earlyassigns<-0
assign_set  <- lead_nodes
for (j in 1:numol) {
  # calculate number that will be assigned
  coltotal <- min((sum(lead_geos[-assign_set,j]==1) + earlyassigns), group_size)
  for (i in 1:net_size) {
   if (i %in% assign_set) {next}
   if (lead_geos[i,j] == 1) {(match_assigns[k,j] <- ids[i])}  
   if (lead_geos[i,j] == 1) {assign_set <- c(assign_set, i)}
   if (lead_geos[i,j] == 1) {(match_counter[j] <- (match_counter[j])+1)}
   if (lead_geos[i,j] == 1) {k=k+1}
    if ((k==group_size+1) | (sum(match_assigns[,j] != 0)==coltotal)) {k=1}      # want to assign that k
    if ((k==group_size+1) | (sum(match_assigns[,j] != 0)==coltotal)) {earlyassigns<-0}
    if ((k==group_size+1) | (sum(match_assigns[,j] != 0)==coltotal)) {j=j+1}    # its size +1 because we still
    if (j==(numol+1))   {break}
    if ((lead_geos[i,j] == 1) & (l < j)) {earlyassigns=earlyassigns+1}
    if ((k==group_size+1) | ((sum(match_assigns[,(j-1)] != 0)==coltotal))) {coltotal <- min((sum(lead_geos[-assign_set,j]==1)), group_size)}
    if ((k==group_size+1) | (sum(match_assigns[,j] != 0)==coltotal)) {i=1}      # want to assign that k
  } # for j
}   # for i
match_assigns


match_assigns2<-match_assigns
lead_geos2<- lead_geos[-assign_set,]
lead_geos2<- as.matrix(lead_geos2)
assign_set2 <- assign_set
ids2 <- ids[-assign_set2]

for (k in 2:3) {
  for (j in 1:numol) {
    if (is.na(lead_geos2[(lead_geos2[, j] == k)][1]))  {next}
    for (i in 1:group_size) { 
      new_assign <- ids2[which(lead_geos2[, j] == k)][[1]]
      if (match_assigns2[i,j] == 0) {ids2 <- ids2[-which(ids2==new_assign)]}                # ditto
      if (match_assigns2[i,j] == 0) 
         {lead_geos2<- lead_geos2[(rownames(lead_geos2)!=new_assign),]}
      if (match_assigns2[i,j] == 0) {match_assigns2[i,j] <- new_assign}
      if (is.na(lead_geos2[(lead_geos2[, j] == k)][1]))  {break}
      #indexes2<- c(j, i, k, new_assign)
      #print(indexes2)
    } # for j
  } # for i
} # for k  

# This step assigns those not yet assiged at D<4

match_assigns3<-match_assigns2
lead_geos3<- lead_geos2
lead_geos3<- lead_geos2
assign_set3 <- assign_set
ids3 <- ids2
newone <- 0

for (j in 1:numol) {
  for (i in 1:group_size) { 
    #new_assign <- ids2[which(lead_geos2[, j] == k)][[1]]
    if (match_assigns3[i,j] == 0) {newone <- 1}
    if (match_assigns3[i,j] == 0) {match_assigns3[i,j] <- ids3[1]}    
    if (newone == 1) {ids3 <-ids3[-1]}
    newone <- 0
    if (length(ids3) == 0) {break}
  } # for j
} # for i
match_assigns <- match_assigns3

# Graph outcome
nodelist2 <- rbind(lead_nodes, match_assigns) 
nodelist   <-t(nodelist2)
nodelist   <- as.data.frame(nodelist)
colnames(nodelist)<- c("id", "nom1", "nom2", "nom3", "nom4", "nom5", "nom6")
dyad3 <-reshape(nodelist, direction="long",  varying = 2:(group_size+1), sep="",
                idvar = "id", timevar="rank")
dyad2  <- dyad3[!is.na(dyad3$nom),]  
dyad1  <- as.data.frame(dyad2)
vars   <- c("id", "nom")
dyad   <- as.matrix(dyad1[vars])
match_net_data <- graph_from_edgelist(dyad)
match_net_data  <- simplify(match_net_data)
match_adj_mat   <- as.matrix(as_adjacency_matrix(match_net_data))
match_net_size  <- gorder(match_net_data)
#Plot the network
png(paste("Optimal Mentor Matching Valente Algorithm", ".png"))
plot(match_net_data, 
     edge.arrow.size   =.3, 
     labels.scale      = 1,
     vertex.color      = "tan",
     vertex.size       = 12,
     vertex.label.font =2, 
     edge.curved=FALSE,
     vertex.label.color = "black"
)
dev.off()
write.table(match_assigns, file="match_assigns.csv", append=T, quote=FALSE, sep=",", row.names=F, col.names=F)
# End Matching #

################
#  Alteration  #
################
# Adding & Deleting Nodes Substituting Diameter+1 for infinite distances 
# Deleting Nodes: Size nodes based on those which reduce cohesion the most when deleted (vitality)
induc_scores <- cbind(ids, 0) 
for(i in 1:net_size) {
  net_data1 <- net_data
  apl <-mean_distance(net_data1)
  net_data_new <- delete.vertices(net_data1, i)
  apl_new <-mean_distance(net_data_new)
  apl_diff <-  apl_new - apl  
  induc_scores[i,2] <- apl_diff
  induc_scores
}
colnames(induc_scores)<- c("id", "vitality")
#We subtract original APL from the new one so large numbers indicate critical nodes 
#In order to sort that way need to reverse the sort largest to smallest


#Graph vitality 
V(net_data)$vitality <- induc_scores[,2]
ones<- rep.int(1, net_size)
vital <- cbind(ids, induc_scores[,2], ones) 
vital <- vital[rev(order(vital[,2])),] 
vital[1:numol, 3] <- 2
vital <- vital[sort.list(vital[,1]), ]
png(paste("Vitality Node Deletion", ".png"))
plot(net_data, 
     edge.arrow.size   =.3, 
     layout            = layout,
     labels.scale      = 1,
     vertex.label.font =2, 
     vertex.color      = vital[,3],
     vertex.label.color= "black",
     vertex.size       = 10,
     edge.curved       =TRUE
)
title(main=paste("Alteration - Node Deletion (Vitality)"))
dev.off()
# Adding nodes 
net_data_bak <- net_data
#pulls the ID numbers from the last no_adds node numbers from the sorted degree list
edge_list_to <- lead_nodes 
net_data_add <- add_vertices(net_data, numol)
edge_list_from <- c((net_size+1):(net_size+numol))
#Generate multiple links to and from new nodes
no_edges   <- numol*(numol-1)
edge_list_from <- sample(edge_list_from, no_edges, replace = TRUE, prob = NULL)
edge_list_to   <- sample(edge_list_to, no_edges, replace = TRUE, prob = NULL)
# Adding vertices does not update the ID number
#new_id <- c(1:(net_size+(length(numol))))
#V(net_data_bak)$new_id <- new_id              # here
# Add edges twice to make symmetric
net_data_add <- add.edges(net_data_add,  c(rbind(edge_list_from, edge_list_to)))

ones<- c(rep.int(1, net_size), rep.int(2, numol))

# Adding & deleting links
# Disconnected = ???
# in order to substitue distance+1 for unreachable nodes must use sna #
library(intergraph)
net_data_s <- intergraph::asNetwork(net_data)
geo_mat<- sna::geodist(net_data_s, inf.replace=(net_size-1))$gdist
#distance_plus1 <-max(geo_mat)+1
#geo_mat<- sna::geodist(net_data_s, inf.replace=distance_plus1)$gdist
apl <- mean(geo_mat)
i<-1
j<-1
Change  <- matrix(data = 0, nrow = net_size, ncol = net_size)
change_dyad_all <- rep.int(0, 3)
for (i in 1:net_size) {
  for (j in 1:net_size) {
    net_data2 <- net_data
    net_data2_adj <- as_adjacency_matrix(net_data2)
    net_data2_adj[i,j]<- abs(net_data2_adj[i,j]-1)
    net_data2 <- graph_from_adjacency_matrix(net_data2_adj)
    net_data2_s <- intergraph::asNetwork(net_data2)
    geo_mat<- sna::geodist(net_data2_s, inf.replace=(net_size-1))$gdist
    #distance_plus1 <-max(geo_mat)+1
    #geo_mat<- sna::geodist(net_data2_s, inf.replace=distance_plus1)$gdist
    apl_new <- mean(geo_mat)
    Change[i,j] <- (apl - apl_new)            # negative values are link deletions
    test <- c(i, j, apl, apl_new, Change[i,j])
    #if (distance_plus1==10) print(test)
    # print(test)
    change_dyad <- cbind(i, j, Change[i,j])
    change_dyad_all <- rbind(change_dyad_all, change_dyad)
  }
}
Change <-round(Change, digits = 4)
round(change_dyad_all, digits = 4)
# Top 10 Links to add 
change_dyad_all<- change_dyad_all[rev(order(change_dyad_all[,3])), ]
change_dyad_all[1:10, ]
# Top 10 Links to delete 
change_dyad_all<- change_dyad_all[order(change_dyad_all[,3]), ]
change_dyad_all[1:10, ]

################### Here
# Showing Links to be Removed
# Still to do: Graph Links to be deleted and added 
###############################
###############################

del_links <- change_dyad_all[1:10,]
#ecolor <- c("gray", "tomato")[E(net_data) %in% unlist(WOMlist) + 1]
#del_links[,3]<- NULL
# ecolor <- c("gray", "tomato")[E(net_data) %in% del_links[,1:2]] #runs but only 2 links
# ecolor <- c("gray", "tomato")[E(net_data), P= del_links[,1:2]]  # error incorrect number of dimensions
# ecolor <- c("gray", "tomato")[E(net_data), P= c(4, 11, 8, 4)]  # error incorrect number of dimensions
# ecolor <- c("gray", "tomato")[E(net_data), %in% (4 11, 8 4)]  # error incorrect number of dimensions
# ecolor <- c("tomato")[E(net_data), P= del_links[,1:2]]          # ditto
# ecolor <- c("gray", "tomato")[E(net_data) %in% del_links[,1:2]] #runs but only 2 links
png(paste("Edges to be Deleted", ".png"))
plot(net_data,
     layout = layout,
     edge.arrow.size   =.1, 
     labels.scale      = 1,
     vertex.color      = "tan",
     vertex.label.font =2, 
     edge.curved=TRUE,
     vertex.label.color = "black",
     edge.color = ecolor,
     vertex.size = 12
)
title(main=paste("Links to Deleted"))
dev.off()

# Graph Links to Add
change_dyad_all<- change_dyad_all[rev(order(change_dyad_all[,3])), ]
add_links <- change_dyad_all[1:10,.]
add_links[,3]<- NULL



title(main=paste("Links to be Added"))

# Rewire
# Use rewire function to make network a small world network
#net_data <- net_data_bak
# Rewire Networks
# Make random
rand_net<- rewire(net_data, keeping_degseq())
png(paste("Random Rewiring", ".png"))
plot(rand_net,
     layout = layout,
     edge.arrow.size   =.1, 
     labels.scale      = 1,
     vertex.color      = "tan",
     vertex.label.font =2, 
     edge.curved=TRUE,
     vertex.label.color = "black",
     vertex.size = 12
)
title(main=paste("Random Network-Same Layout"))
dev.off()

# Make Small World Network
library(netdiffuseR)
small_world <- rgraph_ws(35, k = 2.7, p = .25)
small_world <- graph_from_adjacency_matrix(small_world)
png(paste("Small World Network", ".png"))
plot(small_world,
     layout = layout,
     edge.arrow.size   =.1, 
     labels.scale      = 1,
     vertex.color      = "tan",
     vertex.label.font =2, 
     edge.curved=TRUE,
     vertex.label.color = "black",
     vertex.size = 12
)
title(main=paste("Small World Network-Same Layout"))
dev.off()

# Make Scale Free Network
scale_free <- rgraph_ba(m0=1, t=34, m = 3, self=FALSE)
scale_free <- graph_from_adjacency_matrix(scale_free)
png(paste("Scale Free Network", ".png"))
plot(scale_free,
     layout = layout,
     edge.arrow.size   =.1, 
     labels.scale      = 1,
     vertex.color      = "tan",
     vertex.label.font =2, 
     edge.curved=TRUE,
     vertex.label.color = "black",
     vertex.size = 12
)
title(main=paste("Scale Free Network-Same Layout"))
dev.off()

# Rewire network so that users are paired with non-users closest to them
V(net_data)$adopters<-as.integer(V(net_data)$toa<5)  # Here were setting adopters to TOA<5
adopters<-as.integer(V(net_data)$toa<5)
#blocks <- sna::blockmodel(net_data, V(net_data)$adopters) Old command

# create difference matrix then geodesic setting non-reachables to 1
adopt_mat_dif <- matrix(0, nr=net_size, nc=net_size)
i<-1 
for (i in 1:net_size) {
  j<-1
  for (j in 1:net_size) {
    #adopt_mat_dif[i, j] <- abs(adopters[i] - adopters[j])
    if ((adopters[i] - adopters[j])>0) {adopt_mat_dif[i, j] <- (adopters[i] - adopters[j])}
    #print(sum(adopt_mat_dif[,j])<=2)
  }
  #if (sum(adopt_mat_dif[i,])>=2) {next}
}

gdist <-sna::geodist(adj_mat, inf.replace=1)
gdist_mat <- gdist$gdist
#dist_lt3 <- gdist_mat < 3
dist_lt2 <- gdist_mat < 2
dif_match <- (adopt_mat_dif==1) & (dist_lt2==TRUE)  # Mmmmmm
#dif_match <- (adopt_mat_dif==1) # Mmmmmm
dif_match_net <- graph_from_adjacency_matrix(dif_match)
match_color <- adopters+1
png(paste("Matching Users with Non-users", ".png"))
plot(dif_match_net,
     layout = layout,
     edge.arrow.size   =.1, 
     labels.scale      = 1,
     vertex.color      = match_color,
     vertex.label.font =2, 
     edge.curved=TRUE,
     vertex.label.color = "black",
     vertex.size = 12
)
title(main=paste("Pairing Users with Non-users"))
dev.off()

############################################################
#                  The End The End                         #
############################################################
