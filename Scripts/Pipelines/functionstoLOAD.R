

### Function name=comdist.Swenson (= Dpw, MPFD).
### This function computes the mean pairwise distance among species between two communities. This metric is also called Dpw or comdist. This function use a fast code based on apply function developed by Nathan Swenson in Swenson 2014 "Functional and Phylogenetic Ecology in R" Springer p99-100
### This function can also be used to compute the Mean pairwise Functional Distance (MPFD) based on any functional distance matrix (e.g. Gower distance matrix).
### Input:
# comm is a community data matrix (all the communities have to contain at least 1 species).
# dis should be a matrix containing the distances among species (e.g. coming from cophenetic() function for instance). If the dis object is a "dist" object
# the function atomatically convert it to a "matrix" object

## The conspecific distances among communities are included in the computation (as for the comdist function in the R package picante).


comdist.Swenson=function(comm, dis){
  list.of.names=apply(comm, 1, function(x) names(x[x>0])) ### list of names species present in each community.
  
  if(class(list.of.names) == "matrix"){
    list.of.names = lapply(seq(1, dim(list.of.names)[2]), function(i) list.of.names[,i])
  }
  diss=as.matrix(dis)
  
  return(as.dist(do.call(cbind,lapply(list.of.names, function(x) lapply(list.of.names, function(z) mean(diss[x,z], na.rm=T))) )))
}


### Function name=comdistnt.Swenson (=Dnn, MNND, Gamma+)
### This function computes the mean nearest distance among species between two communities. This metric is also called Dnn in Swenson which is equivalent of the Gamma + metric/100 defined by Clarke and Warwick 2006 or the MNND (mean nearest neighbour distance).
### this is the phylo beta version of the MNTD, comparing the closest relative of each species between communities.
### This function can also be used to estimate functional beta diversity using the mean nearest neighbour distance (MNND) approach, using any type of functional distance matrix (e.g. Gower distance matrix).
### This metric is influcence by the terminal branches of the phylogenetic tree (or close functional distance).
### Input:
# comm is a community data matrix (all the communities have to contain at least 1 species).
# dis should be a matrix containing the distances among species (e.g. coming from cophenetic() function for instance). If the dis object is a "dist" object
# the function atomatically convert it to a "matrix" object.

## The conspecific distances among communities are included in the computation (as for the comdist function in the R package picante).

comdistnt.Swenson=function(comm, dis){
  list.of.names=apply(comm, 1, function(x) names(x[x>0])) ### list of names species present in each community.
  
  if(class(list.of.names) == "matrix"){
    list.of.names = lapply(seq(1, dim(list.of.names)[2]), function(i) list.of.names[,i])
  }
  diss=as.matrix(dis)
  matres=as.dist(do.call(cbind,lapply(list.of.names, function(x) lapply(list.of.names, function(z) {
    if(length(x)==1 | length(z)==1){
      mean(c(min(diss[x,z], na.rm=T), diss[x,z]), na.rm=T)
    } else {
      mean(c(apply(diss[x,z], MARGIN=1, min, na.rm=T), apply(diss[x,z], MARGIN=2, min, na.rm=T)), na.rm=T)
    }
  }
  ))))
  return(matres)
}

####################################
get_beta_val <- function(Data=Beta1,nb_rd=nb_rd,sites= sites_list,names="jac_diss"){
  
  DFtemp <- array(data=NA,c(nrow(as.matrix(Data[[1]])),ncol(as.matrix(Data[[1]])),nb_rd)); i=1
  for ( i in 1:nb_rd) DFtemp[,,i] <-  as.matrix(Data[[i]])
  
  res=apply(DFtemp,c(1,2),mean)  
  res_SD=apply(DFtemp,c(1,2),sd)  
  
  rownames(res) <- colnames(res) <- sites[-length(sites)]
  res[upper.tri(res,diag=T)] <- NA
  res_SD[upper.tri(res_SD,diag=T)] <- NA
  res_tot <- na.omit(cbind(expand.grid(dimnames(res)), Mean = as.vector(res),SD=as.vector(res_SD)))
  colnames(res_tot)[c(3,4)] <- paste(names,colnames(res_tot)[c(3,4)],sep="_")
  res_tot
}

########################################@
Spidercam <- function (df, axistype = 0, seg = 4, pty = 16, pcol = 1:8, plty = 1:6,
                       plwd = 1, pdensity = NULL, pangle = 45, pfcol = NA, cglty = 3,
                       cglwd = 1, cglcol = "navy", axislabcol = "blue", title = "",
                       maxmin = TRUE, na.itp = TRUE, centerzero = FALSE, vlabels = NULL,
                       vlcex = NULL, caxislabels = NULL, calcex = NULL, paxislabels = NULL,
                       palcex = NULL, ...)
{
  if (!is.data.frame(df)) {
    cat("The data must be given as dataframe.\n")
    return()
  }
  if ((n <- length(df)) < 3) {
    cat("The number of variables must be 3 or more.\n")
    return()
  }
  if (maxmin == FALSE) {
    dfmax <- apply(df, 2, max)
    dfmin <- apply(df, 2, min)
    df <- rbind(dfmax, dfmin, df)
  }
  plot(c(-1.2, 1.2), c(-1.2, 1.2), type = "n", frame.plot = FALSE,
       axes = FALSE, xlab = "", ylab = "", main = title, asp = 1,
       ...)
  theta <- seq(90, 450, length = n + 1) * pi/180
  theta <- theta[1:n]
  xx <- cos(theta)
  yy <- sin(theta)
  CGap <- ifelse(centerzero, 0, 1)
  for (i in 0:seg) {
    polygon(xx * (i + CGap)/(seg + CGap), yy * (i + CGap)/(seg + CGap), lty = cglty, lwd = cglwd, border = cglcol)
    if (axistype == 1 | axistype == 3)
      CAXISLABELS <- paste(i/seg * 100, "(%)")
    if (axistype == 4 | axistype == 5)
      CAXISLABELS <- sprintf("%3.2f", i/seg)
    if (!is.null(caxislabels) & (i < length(caxislabels)))
      CAXISLABELS <- caxislabels[i + 1]
    if (axistype == 1 | axistype == 3 | axistype == 4 | axistype ==
        5) {
      if (is.null(calcex))
        text(-0.05, (i + CGap)/(seg + CGap), CAXISLABELS,col = axislabcol)
      else text(-0.05, (i + CGap)/(seg + CGap), CAXISLABELS,
                col = axislabcol, cex = calcex)
    }
  }
  if (centerzero) {
    arrows(0, 0, xx * 1, yy * 1, lwd = cglwd, lty = cglty,length = 0, col = cglcol)
  }
  else {
    arrows(xx/(seg + CGap), yy/(seg + CGap), xx * 1, yy * 1, lwd = cglwd, lty = cglty, length = 0, col = cglcol)
  }
  PAXISLABELS <- df[1, 1:n]
  if (!is.null(paxislabels))
    PAXISLABELS <- paxislabels
  if (axistype == 2 | axistype == 3 | axistype == 5) {
    if (is.null(palcex))
      text(xx[1:n], yy[1:n], PAXISLABELS, col = axislabcol)
    else text(xx[1:n], yy[1:n], PAXISLABELS, col = axislabcol,
              cex = palcex)
  }
  VLABELS <- colnames(df)
  if (!is.null(vlabels))
    VLABELS <- vlabels
  if (is.null(vlcex))
    text(xx * 1.2, yy * 1.2, VLABELS)
  else text(xx * 1.2, yy * 1.2, VLABELS, cex = vlcex)
  series <- length(df[[1]])
  SX <- series - 2
  if (length(pty) < SX) {
    ptys <- rep(pty, SX)
  }
  else {
    ptys <- pty
  }
  if (length(pcol) < SX) {
    pcols <- rep(pcol, SX)
  }
  else {
    pcols <- pcol
  }
  if (length(plty) < SX) {
    pltys <- rep(plty, SX)
  }
  else {
    pltys <- plty
  }
  if (length(plwd) < SX) {
    plwds <- rep(plwd, SX)
  }
  else {
    plwds <- plwd
  }
  if (length(pdensity) < SX) {
    pdensities <- rep(pdensity, SX)
  }
  else {
    pdensities <- pdensity
  }
  if (length(pangle) < SX) {
    pangles <- rep(pangle, SX)
  }
  else {
    pangles <- pangle
  }
  if (length(pfcol) < SX) {
    pfcols <- rep(pfcol, SX)
  }
  else {
    pfcols <- pfcol
  }
  for (i in 3:series) {
    xxs <- xx
    yys <- yy
    scale <- CGap/(seg + CGap) + (df[i, ] - df[2, ])/(df[1,
    ] - df[2, ]) * seg/(seg + CGap)
    if (sum(!is.na(df[i, ])) < 3) {
      cat(sprintf("[DATA NOT ENOUGH] at %d\n%g\n", i, df[i,
      ]))
    }
    else {
      for (j in 1:n) {
        if (is.na(df[i, j])) {
          if (na.itp) {
            left <- ifelse(j > 1, j - 1, n)
            while (is.na(df[i, left])) {
              left <- ifelse(left > 1, left - 1, n)
            }
            right <- ifelse(j < n, j + 1, 1)
            while (is.na(df[i, right])) {
              right <- ifelse(right < n, right + 1, 1)
            }
            xxleft <- xx[left] * CGap/(seg + CGap) +
              xx[left] * (df[i, left] - df[2, left])/(df[1,left] - df[2, left]) * seg/(seg + CGap)
            yyleft <- yy[left] * CGap/(seg + CGap) +
              yy[left] * (df[i, left] - df[2, left])/(df[1,
                                                         left] - df[2, left]) * seg/(seg + CGap)
            xxright <- xx[right] * CGap/(seg + CGap) +
              xx[right] * (df[i, right] - df[2, right])/(df[1,
                                                            right] - df[2, right]) * seg/(seg + CGap)
            yyright <- yy[right] * CGap/(seg + CGap) +
              yy[right] * (df[i, right] - df[2, right])/(df[1,
                                                            right] - df[2, right]) * seg/(seg + CGap)
            if (xxleft > xxright) {
              xxtmp <- xxleft
              yytmp <- yyleft
              xxleft <- xxright
              yyleft <- yyright
              xxright <- xxtmp
              yyright <- yytmp
            }
            xxs[j] <- xx[j] * (yyleft * xxright - yyright * xxleft)/(yy[j] * (xxright - xxleft) - xx[j] * (yyright - yyleft))
            yys[j] <- (yy[j]/xx[j]) * xxs[j]
          }
          else {
            xxs[j] <- 0
            yys[j] <- 0
          }
        }
        else {
          xxs[j] <- xx[j]*CGap/(seg+CGap)+xx[j]*(df[i,j]-df[2,j])/(df[1,j]-df[2,j])*seg/(seg+CGap)
          yys[j] <- yy[j]*CGap/(seg+CGap)+yy[j]*(df[i,j]-df[2,j])/(df[1,j]-df[2,j])*seg/(seg+CGap)
        }
      }
      if (is.null(pdensities)) {
        polygon(xxs, yys, lty = pltys[i - 2], lwd = plwds[i - 2], border = pcols[i - 2], col = pfcols[i - 2])
      }
      else {
        polygon(xxs,yys,lty=pltys[i-2],lwd=plwds[i-2],border=pcols[i-2],density=pdensities[i-2],
                angle=pangles[i-2],col=pfcols[i-2])
      }
      points(xx * scale, yy * scale, pch =21, col = "grey30",bg=pcols[i - 2],lwd=0.3,cex=0.75)
    }
  }
}

###########################################################@
get_FAC <- function ( PCOA, samp, nb_rep=100){
  a <- seq(1,nrow(samp),1)
  
  FAC <- lapply(1:nb_rep, function(j){
    Indic <- vector()
    cat("j",j,"\n")
    for (i in 1: length(a)){
      comm <- samp[sample(a,i),]
      if(is.null(dim(comm))==T){
        comm <- comm[-which(comm==0)]
        comm <- t(as.matrix(comm)); rownames(comm) <- "test"
      } else {
        comm <- apply(comm,2,sum)
        if(length(which(comm==0))!=0){
          comm <- comm[-which(comm==0)]
        }
        comm <- t(as.matrix(comm)); rownames(comm) <- "test"
      }
      Indic[i] <- mFD::alpha.fd.multidim(as.matrix(PCOA[,c(1:5)]), comm,
                                         ind_vect = "fric",scaling=F,check_input=T,details_returned=F)$functional_diversity_indices[,"fric"]
    }
    Indic
  })
  Out <- do.call(rbind,FAC)
  list(Mean = apply(Out,2,mean), Sd = apply(Out,2,sd))
}

###########################################################@
get_SAC <- function (samp, nb_rep=100){
  a <- seq(1,nrow(samp),1)
  
  SAC <- lapply(1:nb_rep, function(j){
    Indic <- vector()
    cat("j",j,"\n")
    for (i in 1:length(a)){
      comm <- samp[sample(a,i),]
      if(is.null(dim(comm))==T){
        comm <- comm[-which(comm==0)]
        comm <- t(as.matrix(comm))
      } else {
        comm <- apply(comm,2,sum)
        if(length(which(comm==0))!=0){
          comm <- comm[-which(comm==0)]
        }
        comm <- t(as.matrix(comm))
      }
      Indic[i] <- length(comm)
    }
    Indic
  })
  
  Out <- do.call(rbind,SAC)
  list(Mean = apply(Out,2,mean), Sd = apply(Out,2,sd))
}

###########################################################@
get_PAC <- function (tree, samp=asb_sp_w, nb_rep=100){
  a <- seq(1,nrow(samp),1)
  
  PAC <- lapply(1:nb_rep, function(j){
    Indic <- vector()
    cat("j",j,"\n")
    for (i in 1: length(a)){
      comm <- samp[sample(a,i),]
      if(is.null(dim(comm))==T){
        comm <- comm[-which(comm==0)]
        comm <- t(as.matrix(comm)); rownames(comm) <- "test"
      } else {
        comm <- apply(comm,2,sum)
        if(length(which(comm==0))!=0){
          comm <- comm[-which(comm==0)]
        }
        comm <- t(as.matrix(comm)); rownames(comm) <- "test"
      }
      Indic[i] <- pd.query(tree=tree,matrix = comm)
    }
    Indic
  })
  Out <- do.call(rbind,PAC)
  list(Mean = apply(Out,2,mean), Sd = apply(Out,2,sd))
}

###########################################################@
get_SAC_Div <- function (samp,n=19, nb_rep=100){
  a <- seq(1,n,1)
  
  SAC <- lapply(1:nb_rep, function(j){
    Indic <- vector()
    cat("j",j,"\n")
    for (i in 1:length(a)){
      comm <- samp[sample(a,i),]
      if(is.null(dim(comm))==T){
        comm <- comm[-which(comm==0)]
        comm <- t(as.matrix(comm))
      } else {
        comm <- apply(comm,2,sum)
        if(length(which(comm==0))!=0){
          comm <- comm[-which(comm==0)]
        }
        comm <- t(as.matrix(comm))
      }
      Indic[i] <- length(comm)
    }
    Indic
  })
  
  Out <- do.call(rbind,SAC)
  list(Mean = apply(Out,2,mean), Sd = apply(Out,2,sd))
}

# function to extract the Average and sd and 0.025 and 97.5 quantiles of different statistics present in a tables across multiple rds tables 
# stored in a folder.

# path = "Results/PhyloAlpha/100trees/"

Summary.Stat = function(path = NULL, prefix = NULL){
  
  if(is.null(prefix)){
    NewInd_100tr = list.files(path)
  } else {
    NewInd_100tr = list.files(path)
    NewInd_100tr = NewInd_100tr[grep(prefix, NewInd_100tr, fixed = TRUE)]
  }
  
  # Convert all the rds tables into a single list.
  NewInd_list = list()
  i = 1
  for(i in 1:length(NewInd_100tr)){
    cat("i = ", i, "\n")
    NewInd_list[[i]] = readRDS(paste(path, NewInd_100tr[i], sep = ""))
  }
  
  nb.met = dim(NewInd_list[[1]])[2]
  
  # Compute the mean
  Met.Mean = as.data.frame(do.call(cbind, lapply(seq(1, nb.met), function(y){
    a.100 = do.call(cbind, lapply(NewInd_list, function(x) {x[,y]}))
    apply(a.100, 1, function(x) mean(x, na.rm = TRUE))
  })))
  names(Met.Mean) = paste(names(NewInd_list[[1]]), ".mean", sep = "")
  row.names(Met.Mean) = row.names(NewInd_list[[1]])
  
  # compute the sd.
  Met.Sd = as.data.frame(do.call(cbind, lapply(seq(1, nb.met), function(y){
    a.100 = do.call(cbind, lapply(NewInd_list, function(x) {x[,y]}))
    apply(a.100, 1, function(x) sd(x, na.rm = TRUE))
  })))
  names(Met.Sd) = paste(names(NewInd_list[[1]]), ".sd", sep = "")
  row.names(Met.Sd) = row.names(NewInd_list[[1]])
  
  # 0.025 quantiles
  Met.Quant0.025 = as.data.frame(do.call(cbind, lapply(seq(1, nb.met), function(y){
    a.100 = do.call(cbind, lapply(NewInd_list, function(x) {x[,y]}))
    apply(a.100, 1, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
  })))
  names(Met.Quant0.025) = paste(names(NewInd_list[[1]]), ".Quant0.025", sep = "")
  row.names(Met.Quant0.025) = row.names(NewInd_list[[1]])
  
  # 0.975 quantiles
  Met.Quant0.975 = as.data.frame(do.call(cbind, lapply(seq(1, nb.met), function(y){
    a.100 = do.call(cbind, lapply(NewInd_list, function(x) {x[,y]}))
    apply(a.100, 1, function(x) quantile(x, probs = 0.975, na.rm = TRUE))
  })))
  names(Met.Quant0.975) = paste(names(NewInd_list[[1]]), ".Quant0.975", sep = "")
  row.names(Met.Quant0.975) = row.names(NewInd_list[[1]])
  
  return(data.frame(Met.Mean, Met.Sd, Met.Quant0.025, Met.Quant0.975))
}
# Compute SR, MPD, MNTD, VPD and  VNTD metrics (See Tucker et al. 2017 for the definition of the metrics)
# Author David Eme, 29/10/2019.

# tree =  a phylogenetic tree, a phylo object,
# matrix =  a data.frame with only one row (one community),
# return a data frame with SR (species richness), MPD (Mean pairwise distance), MNTD (Mean nearest taxonomic distance), 
# VPD (Variance of the pairwise distance) and VNTD (Variance of the nearest taxonomic distance).

require(ape)
require(PhyloMeasures)

PD.NRI.NTI.MPD.MNTD.VPD.VNTD = function(tree = NULL, matrix = NULL, null.model = NULL, reps = 999){
  SR = sum(matrix)
  NTI<-mntd.query(tree=tree, matrix=matrix , standardize = TRUE, null.model=null.model, reps=reps)
  NRI<-mpd.query(tree=tree, matrix=matrix, standardize = TRUE, null.model=null.model, reps=reps)
  PD<-pd.query(tree=tree, matrix= matrix, null.model=null.model, reps=reps)
  
  sampleTaxa = colnames(matrix)[which(matrix > 0)]
  
  SpToRemove = setdiff(tree$tip.label, sampleTaxa) # species list to remove from the tree.
  tmp.tree=drop.tip(tree, SpToRemove) # prune tree
  
  if(SR > 1) {
    com.distD = cophenetic(tmp.tree)
    com.dist=as.dist(com.distD)
    MPD=mean(com.dist) ### MPD=AvTD=Delta+. ="Dispersion index"
    VPD=((sum(com.dist^2))/((SR*(SR-1))/2))-(mean(com.dist)^2) ### VPD=Var.MPD (Tucker et al. 2017) give the same results that the Lambda (VarTD: Clarke & Warwick 2001) computed by taxondive function from vegan. ="Regularity index"
    diag(com.distD)=NA
    MNTD=mean(apply(com.distD, MARGIN=1, min, na.rm=T)) ### Mean nearest taxon distance (mntd) See Swenson 2014. ="Dispersion index"
    VNTD=sum((apply(com.distD, MARGIN=1, min, na.rm=T)-mean(apply(com.distD, MARGIN=1, min, na.rm=T)))^2)/SR ### Variance of MNTD (Tucker et al. 2017). ="Regularity index"
  } else {
    MPD = NA
    VPD = NA
    MNTD = NA
    VNTD = NA
  }
  data.frame(SR = SR, PD=PD, NTI=NTI,NRI=NRI, MPD = MPD, MNTD = MNTD, VPD = VPD, VNTD = VNTD)
}

# Compute the function on multiple trees and large community data matrix, it runs in parallel.

# trees = a "multiPhylo" object,
# Com.Mat =  the community data matrix with the species names matching perfectly
# the tip label of the trees.
# save.path =  name of the path storing the .rds files (e.g. "Results/Indices_100trees/")
# nthreads  = number of threads that shoukd be used to run the analyses in parallel.

require(parallel)

PD.NRI.NTI.MPD.MNTD.VPD.VNTD.multi.trees = function(trees = NULL, Com.Mat = NULL, save.path = NULL, nthreads = NULL, null.model = NULL, reps = 999){ 
  
  if(class(trees) == "phylo" ){
    trees = list(trees)
  }
  
  lapply(1:length(trees), function(y){
    
    Index = lapply(1:nrow(Com.Mat), function(x){
      cat("x=",x,"\n")
      PD.NRI.NTI.MPD.MNTD.VPD.VNTD (tree = trees[[y]], matrix = Com.Mat[x,], null.model = null.model, reps = reps)
    })
    
    Index_results <- do.call(rbind,Index)
    row.names(Index_results) = row.names(Com.Mat)
    
    saveRDS(Index_results, paste0(save.path, "Index_results_Tree_",y,".rds"))
  })
}
