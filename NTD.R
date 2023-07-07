##Required Inputs
##phy = a phylogenetic tree in the form of a phylo object
##OTU = OTU table (ASV: rows, communities: columns, with raw sequence counts or cell counts)

library(phyloseq)
pseq_rarefy <- readRDS("pseq_rarefy.RDS")
phy = phy_tree(pseq_rarefy)
OTU = as.data.frame(otu_table(pseq_rarefy))

detach("package:phyloseq", unload = TRUE)

library(parallel)
library(picante)
library(data.table)
library(dplyr)

ncores = round(detectCores() / 2) ## the number of desired cores to use for parallelization.

timestamp()
print("Loading inputs")
OTU <- apply(OTU,2,FUN=function(x) x/sum(x)) # Relative abundances
species <- rownames(OTU)
phydist <- cophenetic(phy)
pdist <- phydist[species,species]

##Checking if the tree has 0-length branches
# c = pdist
# diag(c)<-NA
# min(c, na.rm=TRUE) #[1] 0
# rm(c)

##Optional fixing of the pdist object in cases where the tree has 0-length branches. We replace the 0s in pdist with the minimum phylogenetic distance across all the tree.
pdist[which(pdist==0)] <- NA
pdist[which(is.na(pdist))] <- head(sort(pdist))[1]
rm(phydist)
# taxonomy <- as.data.table(taxonomy)
# names(taxonomy) <- c("species", "taxonomy")

##Setting rows (species) without any absences to have one absence in the sample (column) where they have the lowest abundance. This allows the inclusion of all species in the analyses
for (i in 1:nrow(OTU)) {
  if (min(OTU[i,]>0)) {
    OTU[i, which(OTU[i,]==min(OTU[i,]))] <- 0
  }
}          

##Make data.table ix with rows like the upper diagonal of your (community*community) matrix
timestamp()
print("Making initial community-community data.table")
ix <- data.table('i'=rep(1:ncol(OTU),times=ncol(OTU):1-1))
ix[,j:=(i+1):ncol(OTU),by=i]

ix[,i:=colnames(OTU)[i]] #give names to community i in the table
ix[,j:=colnames(OTU)[j]] #give names to community j in the table

##Load sub-functions

NTD_calc <- function(i,j,spp,y,pdist,OTU){
  
  if (y[i]>0 & y[j]==0){
    other.species <- setdiff(colnames(pdist),spp)
    other.species <- intersect(other.species,rownames(OTU)[which(OTU[,j]>0)])
    min.dist <- min(pdist[spp,other.species])
    a <- y[i]
  } else if (y[j]>0 & y[i]==0){
    other.species <- setdiff(colnames(pdist),spp)
    other.species <- intersect(other.species,rownames(OTU)[which(OTU[,i]>0)])
    min.dist <- min(pdist[spp,other.species])
    a=y[j]
  }
  names(a) <- NULL
  return(c('abundance'=a,'distance'=min.dist))
}

getNTDvecs <- function(spp,OTU,ix,pdist.=pdist){
  y <- OTU[spp,]
  #limit community-community comparisons to where species is present in one, absent in another
  signed_abunds <- sign(as.numeric(y>0)-0.5)
  names(signed_abunds) <- names(y)
  contributing_pairs <- which(signed_abunds[ix$i]*signed_abunds[ix$j]<0)
  ix <- ix[contributing_pairs,]
  
  output <- matrix(NA,nrow=nrow(ix),ncol=2)
  for (ii in 1:nrow(ix)){
    output[ii,]<-unlist(NTD_calc(ix$i[ii],ix$j[ii],spp,y,pdist,OTU))
  }
  colnames(output) <- c('abundances','distances')
  output <- as.data.table(output)
  output <- cbind(output,ix)
  output[,pres_in_i:=y[i]>0]
  return(output)
}


cl <- makeCluster(ncores)
clusterEvalQ(cl,library(data.table))
clusterExport(cl,varlist=c('pdist','getNTDvecs','NTD_calc','ix'))
timestamp()
print("Making initial species table")
species_NTD_effects <- parLapply(cl,species,getNTDvecs,OTU=OTU,ix=ix)
stopCluster(cl)
rm('cl')

names(species_NTD_effects) <- species
NTD_data <- species_NTD_effects %>%  rbindlist(use.names = T,idcol = T)
names(NTD_data)[1] <- 'species'


NTD_data[,community_index:=paste(i,j,sep='_')]
NTD_data[,community_N:=.N,by=community_index]

setkey(NTD_data,species,community_index)

ntd_spp <- NTD_data[,list(mean_ntd=mean(distances)), by='species']

saveRDS(ntd_spp,file="ntd_spp.RDS")

