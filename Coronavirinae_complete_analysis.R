
#Loading packages
library(MASS) #sammon  
library(SyNet) #reponderar las distancias
library(dummies) #dumificar
library(igraph)

#Read R objects (matrices of distances and data.frame with attributes)
e <- readRDS("Data/Protein_distances/e", refhook = NULL)
m <- readRDS("Data/Protein_distances/m", refhook = NULL)
n <- readRDS("Data/Protein_distances/n", refhook = NULL)
s <- readRDS("Data/Protein_distances/s", refhook = NULL)
atributo <- readRDS("Data/Virus_data/atributo", refhook = NULL) 

#Locating some relevant isolates
grep("Wuhan", atributo[,18])
grep("RaTG13", atributo[,18])
grep("Pangolin", atributo[,18])


#Calculate a single matrix of distances between virions
distE <- 0.5*(sweep(e, 2, apply(e, 2, max), "/") + sweep(e, 1, apply(e, 1, max), "/"))  
distM <- 0.5*(sweep(m, 2, apply(m, 2, max), "/") + sweep(m, 1, apply(m, 1, max), "/"))  
distN <- 0.5*(sweep(n, 2, apply(n, 2, max), "/") + sweep(n, 1, apply(n, 1, max), "/"))  
distS <- 0.5*(sweep(s, 2, apply(s, 2, max), "/") + sweep(s, 1, apply(s, 1, max), "/"))  
distvirus <- (distE + distM + distN + distS)/4
colnames(distvirus) <- rownames(distvirus) <- atributo[,1] #label by genome accession number


#############################################################
#Chunk of code leading to the basics of Figure 1 of the paper
#############################################################

#Preliminaries to obtain the proximity network (minimum spanning tree)
redv <- graph_from_adjacency_matrix(distvirus, mode = "undirected", diag = FALSE, weighted = T) -> redv
#Since edges with weight equals zero are overriden, we include them through the next lines 
olde <- E(redv) #saving edges
redv[V(redv), V(redv)] <- TRUE #adding all possible edges
E(redv)$weight <- 0 #all weights is 0 
E(redv)[olde]$weight <- olde$weight
redv <- simplify(redv) 
#Finally, obtain the proximity network
mstv <- igraph::mst(redv)
#Colouring of nodes based on the taxonomic genus of the sample
V(mstv)$color <-  c("gray", "cyan", "magenta", "yellow", "green")[as.integer(factor(atributo[,4]))]
#Genera are: "Alphacoronavirus" "Betacoronavirus"  "Deltacoronavirus" "Gammacoronavirus"
#There are two items of unkonwon genus. We predict they are betacoronavirus 
plot(mstv) #Draw the network. Try tkplot(mstv) to see an interactive network

#############################################################
#End of lines for producing Figure 1
#############################################################


#Loading phylogeny of hosts 
library(ape)
filo <- read.nexus("Data/Host_phylogeny/phyhost.nex")
filobrl <- compute.brlen(filo)
distancias <- cophenetic(filobrl) #Calculate the matrix of patristic distances
par(mar = rep(0, 4))
plot(filobrl) #Draw the cladogram

host <- unlist(lapply(strsplit(atributo$Host, " "), paste, collapse = "_"))
cual <- match(host, unique(c(filo$node.label, filo$tip.label))) #If cual is 1 the host is unknown
#Some manipulations are necessary to deal with uncertainties and 
#different degrees of taxonomic resolution 
listadesc <- list() #list of descendants
for(i in host){
  if(i == "") {listadesc <- c(listadesc, list(NA)); next}
  id <- match(i, c(filo$tip.label, filo$node.label))
  listadesc <- c(listadesc, phangorn::Descendants(filo, id, "tip"))
}

#Calculate the phylogenetic distance between hosts of sampled viruses
filodist <- matrix(0, length(host), length(host))
for(i in 1:(length(host) - 1))
      for(j in (i + 1):(length(host))){
         filodist[i, j] <- filodist[j, i] <- mean(distancias[listadesc[[i]], listadesc[[j]]])      
}


####################################################
#Chunk of code leading to the Figure 2 of the paper
####################################################
aristas <- ends(mstv, E(mstv), F) #Endpoints of edges
coordx1 <- unlist(lapply(listadesc[aristas[,1]], mean, na.rm = T)) + runif(nrow(aristas)) - 0.5
coordx2 <- unlist(lapply(listadesc[aristas[,2]], mean, na.rm = T)) + runif(nrow(aristas)) - 0.5
distvirus[aristas] -> distedges
parabolico <- function(x1, x2, altura, x) {
     a <- (-1*altura)/(x1 - 0.5*(x1 + x2))^2
     out <- a*(x - 0.5*(x1 + x2))^2 + altura
     out
}
layout(matrix(c(1,1,2,1,1,2), nrow = 2, byrow = T))
plot(filobrl, mar = rep(1, 4))
plot(seq(0, max(distedges), length.out = length(filo$tip.label)), 
     1:length(filo$tip.label),  type = "n", mar = rep(1, 4), xlab = "Viral distance")
for(i in 1:nrow(aristas)){
 if(is.na(coordx1[i])) next
 if(is.na(coordx2[i])) next
 nx <- seq(coordx1[i], coordx2[i], length.out = 100)
 ny <- parabolico(coordx1[i], coordx2[i], distedges[i], nx)
 lines(ny, nx)
}

####################################################
#End of lines for producing Figure 2
####################################################



####################################################
#Chunk of code leading to the Figure 3 of the paper
####################################################

par(mfrow = c(1, 1)) #Set the graphical device
aristas <- ends(mstv, E(mstv), F) #Endpoints of edges, two-columns matrix
#Observed sample quantiles
quantilex <- quantile(distvirus[aristas], p = seq(0, 1, 0.1))
quantiley <- quantile(filodist[aristas], p = seq(0, 1, 0.1), na.rm = TRUE)
out <- c()
#Random allocation of hosts using 10,000 runs
for(aux in 1:10000){
   azar <- sample(1:173) #simulate a random shuffling of indices to denote hosts
   quantfilo <- quantile(filodist[cbind(azar[aristas[,1]], azar[aristas[,2]])], p = seq(0, 1, 0.1), na.rm = TRUE)
   out <- c(out, quantfilo)
}

plot(rep(quantilex, 10000), out)
meanquantfil <- tapply(out, rep(1:length(quantilex), 10000), mean)
lowerbound <- tapply(out, rep(1:length(quantilex), 10000), quantile, p = 0.025)
upperbound <- tapply(out, rep(1:length(quantilex), 10000), quantile, p = 0.975)


plot(quantilex, quantiley, type = "l", lty = "dotted", 
     xlab = "Sample quantiles of viral distance", ylab = "Sample quantiles of phylogenetic distance")
lines(quantilex, meanquantfil, lwd = 1.5)
lines(quantilex, lowerbound, lty = "dotted")
lines(quantilex, upperbound, lty = "dotted")
polygon(c(quantilex, rev(quantilex)), c(lowerbound, rev(upperbound)), col = rgb(0, 0, 0, 0.14), border = NA)
points(quantilex, quantiley, pch = 21, bg = rev(rainbow(length(quantilex), end = 0.7)), cex = 1.5)
points(quantilex, meanquantfil, pch = 22, bg =  rev(rainbow(length(quantilex), end = 0.7)), cex = 1.5)
#Reference for the probabilities (edited then in the final generated figure)
points(rep(0.5, length(quantilex)) , seq(0.02, 1.2, length.out = length(quantilex)), pch = 21, bg = rev(rainbow(length(quantilex), end = 0.7)), cex = 1.5)
text(rep(0.5, length(quantilex)) , seq(0.02, 1.2, length.out = length(quantilex)), seq(0, 1, 0.1), pos = 4)

####################################################
#End of lines for producing Figure 3
####################################################


####################################################
#Chunk of code leading to the Figure 4C of the paper
####################################################
x <- c()
distfilo <- c()
distquimeral <- c()
for(xx in 1:100000){
       rnd <- sample(1:173, 4)
       heterotopicdisaff <- mean(c(distE[rnd[1], rnd[2:4]], distM[rnd[2], rnd[c(1, 3, 4)]], 
                                   distN[rnd[3], rnd[c(1, 2, 4)]], distS[rnd[4], rnd[1:3]]))
       x <- c(x, mean(distvirus[rnd, rnd], na.rm = T))
       distfilo <- c(distfilo, mean(filodist[rnd, rnd], na.rm = T))
       distquimeral <- c(distquimeral, heterotopicdisaff)
}
#Create the respective raster
library(raster)
library(RColorBrewer)
pts <- data.frame(lon = distfilo, lat = distquimeral, vals = rep(1, length(distfilo)))
coordinates(pts) <- ~lon+lat
rast <- raster(ncol = 50, nrow = 50)
extent(rast) <- extent(pts)
out <- rasterize(pts, rast, pts$vals, fun = sum)
plot(out)
out[] <- log(out[], 2)
cuts=seq(0, 11) #set breaks
pal <- colorRampPalette(c("blue","yellow", "red"))
par(mar = rep(4, 4))

plot(out, breaks=cuts, col = pal(12), lab.breaks = 2^(0:11), asp = 1.3, 
     xlab = "Phylogenetic distance", ylab = "Degree of chimerality") #plot with defined breaks
                                                                     #the color scale bar means frequency
####################################################
#End of lines for producing Figure 4C
####################################################


