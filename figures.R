library(foreign)
library(netplot)
library(grid)
library(igraph)


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
# broker <- (calculate.EV.brokerage(net_data) + .001)
marginals <-  as.integer(indeg <= 1)
persdens<- lapply(make_ego_graph(net_data, mode="out"), graph.density)
# reciprocity rate is the outdegree of net symmetrized on strong
recip        <- degree(net_data_sym_s, mode = "out") 
ids <-V(net_data)
measures<-cbind(ids,  outdeg, indeg, between, outcls, incls, inv_constr, bridge_w,
                bridge_s, persdens, recip)
write.table(measures, file="graphs/measures.csv", append=T, quote=FALSE, sep=",", row.names=F, col.names=F)

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
graphics.off()
for (i in c(2, 4, 6, 8)) {
  
  leaders <- cbind(ids, indeg, ones, between, ones, incls, ones, keyplay, ones) 
  leaders <- leaders[rev(order(leaders[,(i)])), ] 
  leaders[1:numol, (i+1)] <- 2
  leaders <- leaders[sort.list(leaders[,1]), ]
  
  #Plot the network
  
  lcolors <- grDevices::hcl.colors(
    2, palette="RdYlBu", alpha = .8, rev = TRUE
    )[leaders[,(i+1)]]
  
  png(paste("graphs/Leader", i, ".png"))
  p <- nplot(net_data,
       layout = layout,
       vertex.label.fontsize = ifelse(leaders[,(i+1)] == 1, 10, 20),
       vertex.label.fontfamily = "Times",
       vertex.label.show     = 1,
       vertex.size.range     = c(.02, .04),
       vertex.color          = lcolors,
       vertex.label.color    = "black", 
       bg.col                = "transparent",
       edge.curvature        = pi/4,
       edge.color            = ~ ego(mix=0, alpha=.1) + alter(mix=1, alpha=.9),
       edge.line.lwd         = .2,
       vertex.label.fontface = ifelse(leaders[,(i+1)] == 1, "plain", "bold")
  )
  
  # Netplot networks need to be printed!
  print(p)
  
  # p <- set_vertex_gpar(p, "core", fill = radialGradient(c("white", "black"),
  #                                                       cx1=.8, cy1=.8, r1=0))
  
  dev.off()
}
