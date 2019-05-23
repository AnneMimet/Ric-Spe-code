###############################################
# statico: analysis of coupled k-tables 

#############################################################################
################################################################################
# load packages
library(ade4)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(data.table)

#############################################################################
################################################################################
# load data

env_tab <- read.table("env.txt", header = T)
biod_tab <- read.table("Community_data_decades.txt", header = T)

#############################################################################
################################################################################
# statical analysis


##########################################
# intrastructure: common part to the three decades (biodiversity metrics)
wit1 <- withinpca(biod_tab[,2:3], as.factor(biod_tab$year), scaling="total")
2  

# structure environment
div_pca <- dudi.pca(data.frame(env_tab[,2:27]), scale = T) 
6  

# compromise environment
wit2 <- wca(div_pca, env_tab$per) #$ratio:  0.993373
6  

kta1 <- ktab.within(wit1) # first k-table: temporal change in biodiversity metrics
kta2 <- ktab.within(wit2) # first k-table: temporal change in environmental variables

# statico is aPartial triadic analysis on the series of covariance tables of each couple of tables
# i.e. pta on coinertie of the coupled tables
statico1 <- statico(kta1, kta2)
2  

######################################
# monte-carlo coinertia tests
statico.krandtest(kta1, kta2)

#############################################################################
################################################################################
# figures 

#######################################################
# figure showing the variables in the compromise 

# prepare the data
statico1$c1  # extract the correlated variables (normed scores)

var <- c("crop", "m_ppt_breed", "alt_mean", "range", "alt_rge",
       "r_ppt_wint", "m_ppt_wint", "impact", "tr_min_breed",
       "m_min_wint", "irrig", "m_mean_wint", "velo_rte", "npp0_rte")
var <- as.data.frame(var)

d <- statico1$co
d$var <- rownames(d)

sel <- merge(var, d, by="var")
sel <- sel[match(var$var,sel$var),]


sel$lab2 <- c("LU1", "P2", "To1", "LU3", "To2",
            "P3", "P1", "LU4", "T3",
            "T1", "LU2", "T2", "LT1", "LT2")

sel$lab3 <- c("Croplands", "Mean Precipitation (breeding season)", "Mean altitude", "Range", "Altitude range",
            "Spatial range of mean precipitation (winter)", "Mean Precipitation (winter)", "Naturalness", "Temporal range of minimum temperature (winter)",
            "Minimum Temperature (winter)", "Irrigation (croplands)", "Mean Temperature (winter)", "Climate stability", "NPP0")

sel$Variables <- factor(c("Land-use", "Precipitation", "Topography", "Land-use", "Topography",
                        "Precipitation", "Precipitation", "Land-use", "Temperature",
                        "Temperature", "Land-use", "Temperature", "Long-term", "Long-term"))


pal <- brewer.pal(5, "Set1")

paltopo <- c("#E41A1C", "#EE7621")

palLU1 <- brewer.pal(9, "Greys")  
palLU <- palLU1[6:9]  #  "#737373" "#525252" "#252525" "#000000"

palP1 <- brewer.pal(9,"Blues")  
palP <- c(palP1[5] , palP1[7], palP1[9])   #   "#6BAED6"  "#2171B5"  "#08306B" 

PalT1 <- brewer.pal(9,"Purples")  
PalTB <- PalT1[7:9]   #  "#6A51A3" "#54278F" "#3F007D"

pallong <- c("#00FFFF", "#00CDCD")

pal2 <- c("#E41A1C", "#EE7621","#737373", "#525252", "#252525", "#000000", "#6BAED6","#2171B5","#08306B",
        "#6A51A3", "#54278F", "#3F007D","#00FFFF", "#00CDCD")

barplot(rep(1,14), col=pal2)


sel2 <-as.data.frame(setorder(sel, cols=lab2))
rownames(sel2) <- 1:14

sel2$col1 <- c("#00FFFF", "#00CDCD", "#737373", "#525252", "#252525", "#000000", "#6BAED6",  "#2171B5", "#08306B",
             "#6A51A3", "#54278F", "#3F007D","#E41A1C", "#EE7621")

####
# Figure

vari <- statico1$li
rownames(vari) <- c("Specialisation", "Richness")
Community <- c("var", "var")

# create point cloud
p <- ggplot(sel2, aes(Comp1, Comp2)) +
  geom_point(color = 'black') +
  geom_abline(intercept=0,slope=0) +
  geom_vline(xintercept=0) +
  theme_classic(base_size = 16) +
  labs(x="Axis 1", y= "Axis 2") +
  theme(axis.title = element_text(color="black", face="bold", size=16),
        axis.title.x =element_text(margin=margin(10,10,10,10)),
        axis.title.y =element_text(margin=margin(10,10,10,10)),
        axis.text  = element_text(vjust=0.5, size=16))
p

p1 <- p + geom_label_repel(aes(label = sel2$lab2, fill = sel2$lab2), color = "white", size = 5.5, label.size=1.6, segment.colour="grey3") +
  scale_fill_manual(values = sel2$col1, name="Environmental Variables", levels(sel2$lab2)) 
p1


p2 <- p1+ geom_label_repel(data= vari, aes(x=Axis1, y=Axis2, label = rownames(vari),  fill = Community), color = "grey3", size = 5.5,label.size=1.6, segment.colour="grey3") +
  scale_fill_manual(values = c(sel2$col1, "white"))+
  geom_point(data= vari, aes(Axis1, Axis2), color = 'black')
theme(legend.position = "bottom") 
p2


#######################################################
# figure of the temporal changes of the variables in the compromise

var <- c("crop","impact", "m_ppt_breed", "m_ppt_wint","r_ppt_wint","tr_mean_breed","tr_min_breed")

var <- as.data.frame(var)

nom <- rownames(statico1$co)
d <- statico1$Tco
d$var <- nom
d$Periode <- factor(c(rep("1973-1983", 26), rep("1985-1995", 26), rep("2001-2011", 26)))


sel <- merge(var, d, by="var")


sel$lab <- c("Crop","Crop","Crop","Naturalness","Naturalness","Naturalness",
           "Prec.breed","Prec.breed","Prec.breed","Prec.wint","Prec.wint","Prec.wint",
           "Prec.var.wint","Prec.var.wint","Prec.var.wint",
           "Tempo.var.temp","Tempo.var.temp","Tempo.var.temp",
           "Tempo.var.temp.min","Tempo.var.temp.min","Tempo.var.temp.min" )

sel$categ <- c("Land-use","Land-use","Land-use","Land-use","Land-use","Land-use",
             "Precipitation","Precipitation","Precipitation",
             "Precipitation","Precipitation","Precipitation",
             "Precipitation","Precipitation","Precipitation",
             "Temperature","Temperature","Temperature",
             "Temperature","Temperature", "Temperature")

sel$lab2 <- c("LU1","LU1","LU1","LU4","LU4","LU4",
            "P2","P2","P2","P1","P1","P1",
            "P3", "P3","P3",
            "T3","T3","T3", 
            "T4", "T4", "T4")


sel$lab3 <- c("Croplands","Croplands","Croplands", "Naturalness","Naturalness","Naturalness",
            "Mean Precipitation (breeding season)","Mean Precipitation (breeding season)","Mean Precipitation (breeding season)",
            "Mean Precipitation (winter)", "Mean Precipitation (winter)","Mean Precipitation (winter)",
            "Spatial range of mean precipitation (winter)","Spatial range of mean precipitation (winter)","Spatial range of mean precipitation (winter)",
            "Temporal range of mean temperature (breeding season)", "Temporal range of mean temperature (breeding season)", "Temporal range of mean temperature (breeding season)",
            "Temporal range of minimum temperature (breeding season)","Temporal range of minimum temperature (breeding season)", "Temporal range of minimum temperature (breeding season)")

sel$col <- c( "#737373","#737373","#737373", "#000000",  "#000000",  "#000000", 
            "#2171B5","#2171B5","#2171B5", 
            "#6BAED6","#6BAED6","#6BAED6",
            "#08306B",  "#08306B", "#08306B",
            "#3F007D", "#3F007D","#3F007D",
            "#BCBDDC", "#BCBDDC","#BCBDDC")

sel2 <- sel[order(sel$lab2, sel$Periode),]
rownames(sel2) <- 1:21

sel2$lab4 <- c("LU1_P1","LU1_P2","LU1_P3","LU4_P1","LU4_P2","LU4_P3",
              "P2_P1","P2_P2","P2_P3","P1_P1","P1_P2","P1_P3",
              "P3_P1", "P3_P2","P3_P3",
              "T3_P1","T3_P2","T3_P3", 
              "T4_P1", "T4_P2", "T4_P3")

sel3 <- sel2[order(sel2$lab4),]


# Figure

sel3$class <- as.character(sel3$Periode)
sel3$class[sel3$class=="1973-1983"] <- "1"
sel3$class[sel3$class=="1985-1995"] <- "2"
sel3$class[sel3$class=="2001-2011"] <- "3"

vari<-statico1$Tli
vari$com <- c("Specialisation", "Richness", "Specialisation", "Richness", "Specialisation", "Richness")
vari$Periode <- c("1973-1983", "1973-1983","1985-1995", "1985-1995", "2001-2011","2001-2011")

vari3 <- vari[vari$Periode=="2001-2011",]

pt3 <- sel2[sel$Periode=="2001-2011",]
colnames(pt3) <- c("var","v1","v2","Periode","Variables", "var.names", "lab2","nom", "col", "per")

p1 <- ggplot(sel3, aes(RS1, RS2)) +
  geom_point(aes(group = lab4, size = Periode, color=lab4)) +
  scale_color_manual(values = sel3$col, guide="none") +
  geom_abline(intercept=0,slope=0) +
  geom_vline(xintercept=0) +
  geom_path(aes(group = var, color = lab4), size = 0.5) +
  theme_classic(base_size = 18) +
  labs(x="Axis 1", y= "Axis 2") +
  theme(axis.title = element_text(color="black", face="bold", size=18),
        axis.title.x =element_text(margin=margin(10,10,10,10)),
        axis.title.y =element_text(margin=margin(10,10,10,10)),
        axis.text  = element_text(vjust=0.5, size=18)) +
  theme(legend.position="none")
p1

p2 <- p1 + geom_label_repel(data=pt3, aes(x=v1, y=v2, label = lab2, fill = lab2), 
                           color = "white", size = 5.5, label.size=1.6, 
                           segment.colour="black",
                           nudge_y = -0.04, nudge_x = -0.2)+
  theme(legend.position = "bottom") +
  scale_fill_manual(values = pt3$col, name="", labels=pt3$nom)+
  theme(legend.position="none")
p2

p3 <- p2+ geom_point(data = vari, aes(CS1, CS2, size = Periode)) +
  geom_path(data = vari,aes(CS1, CS2,group = com), size = 0.5) +
  geom_label_repel(data=vari3, aes(x=CS1, y= CS2, label = com),size = 5.5, label.size=1.4, segment.colour="grey3")
p3


#######################################################
# Figures showing the variations of important environmental variables along the statico axes


################################
# include routes coordinates 
sper1 <- statico1$supIY[1:1282,]
sper2 <- statico1$supIY[1283:2564,]
sper3 <- statico1$supIY[2565:3846,]

sper <- cbind(sper1, sper2, sper3)
sper$rte_id <- biod_tab[1:1282,1]

# var environmentales influentes dsle tps

temp_tr_mean_breed <- as.data.frame((env3$tr_mean_breed - env1$tr_mean_breed)*100/env1$tr_mean_breed)
temp_tr_min_breed <- as.data.frame((env3$tr_min_breed - env1$tr_min_breed)*100/env1$tr_min_breed)
temp_m_ppt_breed <- as.data.frame((env3$m_ppt_breed - env1$m_ppt_breed)*100/env1$m_ppt_breed)
temp_crop <- as.data.frame((env3$crop - env1$crop)*100/env1$crop)
temp_naturalness <- as.data.frame((env3$impact - env1$impact)*100/env1$impact)
temp_m_ppt_wint <- as.data.frame((env3$m_ppt_wint - env1$m_ppt_wint)*100/env1$m_ppt_wint)
temp_r_ppt_wint <- as.data.frame((env3$r_ppt_wint - env1$r_ppt_wint)*100/env1$r_ppt_wint)

colnames(temp_tr_mean_breed) <- colnames(temp_tr_min_breed)<-colnames(temp_m_ppt_breed)<-"nom"
colnames(temp_crop) <- colnames(temp_naturalness)<-colnames(temp_m_ppt_wint)<-colnames(temp_r_ppt_wint)<-"nom"

tab_per1 <- cbind(env1$tr_mean_breed, env1$tr_min_breed,env1$m_ppt_breed,env1$crop, 
                env1$impact,env1$m_ppt_wint,env1$r_ppt_wint,sper1) 

tab_per3 <- cbind(env3$tr_mean_breed, env3$tr_min_breed,env3$m_ppt_breed,env3$crop, 
                env3$impact,env3$m_ppt_wint,env3$r_ppt_wint,sper3) 

tab_diff <- cbind(temp_tr_mean_breed, temp_tr_min_breed,temp_m_ppt_breed,temp_crop, 
                temp_naturalness,temp_m_ppt_wint,temp_r_ppt_wint,sper1)

colnames(tab_diff) <- colnames(tab_per1)<-colnames(tab_per3)<-c("tr_mean_breed","tr_min_breed","m_ppt_breed","crop",
                                                              "naturalness","m_ppt_wint","r_ppt_wint","sco1", "sco2")

tab_per1$per <- "per1"
tab_per3$per <- "per3"
tab_diff$per <- "diff"

tab <- rbind(tab_per1, tab_per3)
tab$rte_id <- env1$rte_id
tab_diff$rte_id <- env1$rte_id

## red: increase, blue: decrease
tab_diff$tr_mean_breed_sign <- ifelse(tab_diff$tr_mean_breed > 0, "#FF4040", "#6495ED")
tab_diff$tr_min_breed_sign <- ifelse(tab_diff$tr_min_breed > 0, "#FF4040", "#6495ED")
tab_diff$m_ppt_breed_sign <- ifelse(tab_diff$m_ppt_breed > 0, "#FF4040", "#6495ED")
tab_diff$crop_sign <- ifelse(tab_diff$crop > 0, "#FF4040", "#6495ED")
tab_diff$naturalness_sign <- ifelse(tab_diff$naturalness > 0, "#FF4040", "#6495ED")
tab_diff$m_ppt_wint_sign <- ifelse(tab_diff$m_ppt_wint > 0, "#FF4040", "#6495ED")
tab_diff$r_ppt_wint_sign <- ifelse(tab_diff$r_ppt_wint > 0, "#FF4040", "#6495ED")

tab_diff <- rbind(tab_diff, tab_diff)

###############################
# Figures

#Axis 1
# mean precipitation during the breeding season
p <- ggplot(tab, aes(sco1, m_ppt_breed)) +
  geom_abline(intercept=0,slope=0, col="grey") +
  geom_vline(xintercept=0, col="grey") +
  geom_path(aes(group = tab$rte_id), size = 0.5, col=tab_diff$m_ppt_breed_sign) +
  xlab("Axis 1")+
  ylab("Mean Precipitation (breeding season)")+
  theme(axis.title = element_text(color="black", face="bold", size=12),
        axis.text  = element_text(vjust=0.5, size=12))
p

# temporal variation of the mean temperature between breeding seasons
p<- ggplot(tab, aes(sco1, tr_mean_breed)) +
  geom_abline(intercept=0,slope=0, col="grey") +
  geom_vline(xintercept=0, col="grey") +
  geom_path(aes(group = tab$rte_id), size = 0.5, col=tab_diff$tr_mean_breed_sign)+
  xlab("Axis 1")+
  ylab("Mean Precipitation (breeding season)")+
  theme(axis.title = element_text(color="black", face="bold", size=12),
        axis.text  = element_text(vjust=0.5, size=12))
p
 
# temporal variation of the minimum temperature between breeding seasons
p<- ggplot(tab, aes(sco1, tr_min_breed)) +
  geom_abline(intercept=0,slope=0, col="grey") +
  geom_vline(xintercept=0, col="grey") +
  geom_path(aes(group = tab$rte_id), size = 0.5, col=tab_diff$tr_min_breed_sign) +
  xlab("Axis 1")+
  ylab("Mean Precipitation (breeding season)")+
  theme(axis.title = element_text(color="black", face="bold", size=12),
        axis.text  = element_text(vjust=0.5, size=12))
p
 
####################
## axis2

# naturalness
p<- ggplot(tab, aes(sco2, naturalness)) +
  geom_abline(intercept=0,slope=0, col="grey") +
  geom_vline(xintercept=0, col="grey") +
  geom_path(aes(group = tab$rte_id), size = 0.5, col=tab_diff$naturalness_sign) +
  xlab("Axis 2")+
  theme(axis.title = element_text(color="black", face="bold", size=12),
        axis.text  = element_text(vjust=0.5, size=12))
p

# cropping
p<- ggplot(tab, aes(sco2, crop)) +
  geom_abline(intercept=0,slope=0, col="grey") +
  geom_vline(xintercept=0, col="grey") +
  geom_path(aes(group = tab$rte_id), size = 0.5, col=tab_diff$crop_sign) +
  xlab("Axis 2")+
  theme(axis.title = element_text(color="black", face="bold", size=12),
        axis.text  = element_text(vjust=0.5, size=12))
p

# mean precipitatiopn in winter
p<- ggplot(tab, aes(sco2, m_ppt_wint)) +
  geom_abline(intercept=0,slope=0, col="grey") +
  geom_vline(xintercept=0, col="grey") +
  geom_path(aes(group = tab$rte_id), size = 0.5, col=tab_diff$m_ppt_wint_sign) +
  xlab("Axis 2")+
  theme(axis.title = element_text(color="black", face="bold", size=12),
        axis.text  = element_text(vjust=0.5, size=12))

p

# spatial range of precipitation in winter
p<- ggplot(tab, aes( sco2, r_ppt_wint)) +
  geom_abline(intercept=0,slope=0, col="grey") +
  geom_vline(xintercept=0, col="grey") +
  geom_path(aes(group = tab$rte_id), size = 0.5, col=tab_diff$r_ppt_wint_sign) +
  xlab("Axis 2")+
  theme(axis.title = element_text(color="black", face="bold", size=12),
        axis.text  = element_text(vjust=0.5, size=12))
p


