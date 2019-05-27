#############################################################################
################################################################################
# load packages
library(indicspecies)
library(vegan)
library(stringr)

#############################################################################
################################################################################
# previous steps: 
# 1. extract community data (abundances per species) from the North Amercan breeding bird survey for the years 
# 2004 to 2007, and selec only the first point of each route, i.e. the only one for each coordinates are available
# 2. from those data, suppress the points belonging to routes that will then be used for the statico analysis to limit circularity issues (1282 routes)
# 3. then select species with less than 10 points of occurence (i.e. mean abundance computed over the 4 years under 0.5)
# 4. For each count point, extract the local land use (Globcover data 2006), and the climate Koppen class (Peel, Finlayson, & McMahon, 2007)

#############################################################################
################################################################################
# load data
tab <- read.table("data_for_mrt.txt")

# Description of the variables of the table
# $rte_ID: id of the route to which belongs the count point (factor)
# X and Y: Latitude and Longitude (numeric)
# world_kopp: the Koppen class (factor)
# loc: the land cover/land use class (Globcover) (factor)

#############################################################################
################################################################################
# data preparation

table(tab$loc)
# 14  20  30  50  60  70  90 100 110 120 130 140 150 160 170 190 200 210 220 230 
# 85  53 224 247  26 320 105 133 207  90 170 314  68   2   9  10   2  46   1   2 

# deletion of loc 160, 170, 200, 220, 230, 210 for which the number of points is too low
tab <- tab[!is.element(tab$loc, c(160, 170, 200, 220, 230, 210, 190)),]  

tab$loc[tab$loc==60] <- 50 # grouping of classes 50 and 60 (closed and open broadleaved forests)

tab$loc[tab$loc==90] <- 70 # grouping of classes 90 and 70 (open and closed needleaved forests)

tab$loc[tab$loc==120] <- 140  # grouping of classes 120 and 140 (Mosaic grassland and close to open hebaceous vegetation)

tab$loc[tab$loc==20] <- 14 # grouping of classes 14 and 20 (rainfed croplands and mosaic croplands)

# suppression of the classes in which the dominant land cover is not fixed (see methods)
tab <- subset(tab, tab$loc!=30)   
tab <- subset(tab, tab$loc!=110)

table(tab$loc)
#  14  50  70 100 130 140 150 
# 138 273 425 133 170 404  68


########################
# simplification koppen climate classification

table(tab$world_kopp)
#  4   5   6   7   8   9  14  15  17  18  19  22  23  25  26  27 
# 11  15   9 103  11  33 148   4   2  33   6   2   1  76 266  60 

# suppression of the classes with less than 10 points
tab <- subset(tab, !tab$world_kopp==-9999)
tab <- subset(tab, !tab$world_kopp==2)
tab <- subset(tab, !tab$world_kopp==3)
tab <- subset(tab, !tab$world_kopp==22)
tab <- subset(tab, !tab$world_kopp==23)
tab <- subset(tab, !tab$world_kopp==17)
tab <- subset(tab, !tab$world_kopp==15)
tab <- subset(tab, !tab$world_kopp==19)

# 4   5   6   7   8   9  14  18  25  26  27 
# 17  35  17 216  24  67 275  58 153 530 138


############
# dividing the data in 2 parts: one for fitting the multivariate regression tree (and define the habitat classes), 
# the second one to compute the species specialization

id <- as.data.frame(unique(tab$rte_ID))
tree <- as.data.frame(sample(id[,1], 650)) 
for_ssi <- as.data.frame(id[!is.element(id[,1],tree[,1]),1]) 

tree <- tab[is.element(tab$rte_ID,tree[,1]),] # data for the multivariate regression tree
for_ssi <- tab[is.element(tab$rte_ID,for_ssi[,1]),] # data for computing species specialization

# write.table(tree, "data_for_tree.txt")
# write.table(for_ssi, "data_for_ssi.txt")

#############################################################################
################################################################################
# Defining the habitat classes using a multivariate regression tree approach (mrt) (library mvpart)

# mvpart is not longer active on CRAN, but can be dowloaded from the archives
# You may need to use a former version of r to run this part of the code

####
# load package
library(mvpart) # no longer available on cran, to be downloaded from the archives

# suppress species with less than 5 occurence points
sp2 <- tree[,10:215]

somme <- t(as.data.frame(colSums(sp2)))
sub1 <- subset(sp2, selec=somme>0.5)
sub1[sub1>0.5] <- 1
s1 <- colSums(sub1)  

s2 <- s1[s1>5]  # 153 sp left
sp_selec <- subset(sp2, select= names(s2))

tab <- cbind(sp_selec, tree[,1:9], tree[,216:220])

##############
# scaling the data

colMu <- colMeans(tab[,1:153]) 
spp2 <- sweep(tab[,1:153], 2, colMu, "-") 
rowMu <- rowMeans(spp2) 
spp2 <- sweep(spp2, 1, rowMu, "-") 

###############################
# run the mrt using presence/ absence

spp_pa <- tab[,1:153]
spp_pa[spp_pa>0] <- 1

esb2 <- mvpart(data.matrix(spp_pa)~ as.factor(tab$world_kopp) +as.factor(tab$loc), data=tab,xval=5, 
             xvmul=500, xv="pick",cp= 0.002, minbucket=30, method="mrt", rsq=T)  # 13 classes, good

###############################
## extracting the classes and the correspondance between class id and the Koppen and land cover classes

tab$cl <- esb2$where

corresp <- unique(tab[,c(160:161,168)])

#write.table(corresp,"correspondance.txt")
#write.table(tab, "class_for_ssi.txt")


#############################################################################
################################################################################
# computing species specialization (ssi) using indicator A of Legendre and Legendre (1997), see methods - package indicspecies

cl <- read.table("class_for_ssi.txt", header=T)
tab <- read.table("data_for_ssi.txt", header=T)
corresp <- read.table("correspondance.txt", header=T)
corresp$comb <- paste(corresp$world_kopp, corresp$loc,sep="_")
tab$comb <- paste(tab$world_kopp, tab$loc,sep="_")

# suppress species with less than 5 occurence points (ie same procedure than on the data used for fitting the mrt)
sp2 <- tab[,10:215] 

somme <- t(as.data.frame(colSums(sp2)))
sub1 <- subset(sp2, selec=somme>0.5)
sub1[sub1>0.5] <- 1
s1 <- colSums(sub1)  

s2 <- s1[s1>5]  # 162 sp left
sp_selec <- subset(sp2, select= names(s2))

tab <- cbind(sp_selec, tab[,1:9], tab[,216:221])

##################################
tab <- merge(tab, corresp, by="comb")
species <- tab[,c(2:163,180)]
comm <- species[,1:163]
class <- levels(as.factor(tab$cl))
sp <- colnames(species[,1:163])

indval <- multipatt(comm, species$cl, duleg=T)
summary(indval, indvalcomp=TRUE)

# species specialization (ssi) is the maximum of the Indicator A values accross the different habitat classes

ssi <- as.data.frame(cbind(row.names(indval$A), apply(indval$A, FUN=max, 1))) 
colnames(ssi) <- c("sp", "ssi")

#############################################################################
################################################################################
# computing community specialization (csi) and species richness at the route scale

####
# load and prepare data
ssi <- read.table("ssi.txt", header=T)
sp <- as.data.frame(ssi$sp)

bbs2 <- read.table("data_rtes.txt", header=T, sep=" ") 
routes <- bbs2

routes[is.na(routes)] <- 0
routes <- as.data.frame(routes)
ab <- routes[,2:733]

# fixing occurrence threashold at 1000 per route
ab[ab>999] <- 1000 

tab2 <- as.data.frame(t(ab))
tab2$sp <- as.factor(rownames(tab2))

tab3 <- merge(tab2, ssi, by.x="sp", by.y="sp")

selec <- t(tab3[,2:115962])

####
# computing community specialization
tab_ssi1 <- as.numeric(tab3[,115963]) * tab3[,2:115962]
ab_ssi1 <- t(tab_ssi1)
colnames(ab_ssi1) <- tab3$sp

tab_csi <- as.data.frame(colSums(tab_ssi1)/colSums(tab3[,2:115962]))

####
# computing richness
tab_csi$ric <- specnumber(selec)

####
# retrieve route id
tab_csi$id_rte_year <- routes[,1]
temp <- do.call("rbind", strsplit(as.character(tab_csi$id_rte_year), "_"))
tab_csi$year <- temp[,3]
tab_csi$route <- temp[,2]
tab_csi$statenum <- temp[,1]
tab_csi$rte_id <- str_sub(as.character(routes$Group.1), 1, nchar(as.character(routes$Group.1))-5)

####
# export database: communitiy and specialzation at the route scale 
colnames(tab_csi) <- c("csi", "ric", "id_rte_year","year","route", "statenum", "rte_id")
write.table(tab_csi, "Community_data.txt")


