library(phyloregion)



n20<-optimal_phyloregion(pb[[1]])
test<-select_linkage(pb[[3]])
#A good clustering algorithm for the distance matrix is:
#  UPGMA with cophenetic correlation = 0.6214848 

#select_linkage(pb[[2]])
#A good clustering algorithm for the distance matrix is:
#  UPGMC with cophenetic correlation = 0.881803 
n20_2<-optimal_phyloregion(pb[[2]], method="centroid")

#A good clustering algorithm for the distance matrix is:
#  UPGMA with cophenetic correlation = 0.9390076 
n20_3<-optimal_phyloregion(pb[[3]])
n50_3<-optimal_phyloregion(pb[[3]], k =50)


#Phyloregions show evolutionary affinities among disjunct assemblages (function plot.phyloregion). 
pb <- phylobeta(matPA, t)


phyreg1 <- phyloregion(pb[[1]], shp=spdf, k=10)
phyreg2 <- phyloregion(pb[[2]], shp=spdf, k=10)
phyreg3 <- phyloregion(pb[[3]], shp=spdf, k=10)


#Ordination of phyloregions in NMDS space shows that different phyloregions differ strongly in evolutionary uniqueness (function plot_NMDS). 
pal <- c(
  
  "#003946",
  
  "#057985",
  
  "#94d2bdff",
  "#C1E5D9",
  
  "#e9d8a6",
  "#ecba53",
  "#ee9b00ff", 
  "#ca6702",
  "#bb3e03ff",
  
  "#9b2226")



p1<-plot.phyloregion(phyreg1, col=pal)
p2<-plot.phyloregion(phyreg2,  col=pal)
p3<-plot.phyloregion(phyreg3,  col=pal)

par(mar=rep(4,4))
plot_NMDS(phyreg1, cex=3, palette(pal))
text_NMDS(phyreg1)


par(mar=rep(4,4))
plot_NMDS(phyreg2, cex=3, col=pal)
text_NMDS(phyreg2)


par(mar=rep(4,4))
plot_NMDS(phyreg2, cex=3, col=pal)
text_NMDS(phyreg2)

save.image(file="ucandoit.RData")
