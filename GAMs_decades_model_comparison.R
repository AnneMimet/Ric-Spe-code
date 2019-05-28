#############################################################################
################################################################################
# load packages

library(visreg)
library(lme4)
library(nlme)
library(MASS)
library(ncf)
library(mgcv)
library(stringr)
library(snowfall)
library(ggplot2)
library(ggpubr)

#############################################################################
################################################################################
# functions

####
# # 1 breakpoint piecewise model: 
# Iteratively test the best combination for the breakpoint locations in a one breakpoint piecewise gam model (spline on the spatial coordinates)
# useful link for coding piecewise regressions: https://www.r-bloggers.com/r-for-ecologists-putting-together-a-piecewise-regression/

simple.bk <- function(tab, breaks, year)
{
  mse <- as.data.frame(matrix(nrow=0, ncol=2))
  colnames(mse) <- c("BP1", "dev.expl")
  temp <- tab[tab$year==year,]
  for (i in 1:length(breaks))
  {
    temp$csi.low <- temp$csi.up <- 0
    temp$above <- temp$below <- 0
    
    temp$csi.low[temp$csi < breaks[i]] <- temp$csi[temp$csi < breaks[i]]
    temp$csi.up[temp$csi > breaks[i]] <- temp$csi[temp$csi > breaks[i]]
    
    temp$above <- temp$csi > breaks[i]
    temp$below <- temp$csi < breaks[i]
    
    fit <- gam(ric ~ below * csi.low + above * csi.up + s(Lati, Longi, k=40), 
               data=temp)
    
    # selection of the best model based on explained deviance
    val <- cbind(breaks[i], summary(fit)[14])
    colnames(val) <- c("BP1", "dev.expl")
    mse <- rbind(mse, val)
    mse$dev.expl <- as.numeric(mse$dev.expl) 
  }
  return(unlist((mse$BP1[which(mse$dev.expl == max(mse$dev.expl))])))
}



####
# 2 breakpoints piecewise model: 
# Iteratively test the best combination of breakpoints locations in a two breakpoint piecewise gam model  (spline on the spatial coordinates)

double.bk <- function(tab, breaks, year)
{
  mse <- as.data.frame(matrix(nrow = 0, ncol = 4))
  colnames(mse) <- c("BP1", "BP2", "id", "dev.expl")
  temp <- tab[tab$year == year,]
  for (i in 1 : length(breaks))
  {
    #B2<-seq(breaks[i],0.37, by=0.001)  
    for(j in (i) : length(breaks))
    {
      temp$csi.mid <- temp$csi.low <- temp$csi.up <- 0
      temp$middle <- temp$above <- temp$below <- 0
      
      temp$csi.low[temp$csi < breaks[i]] <- temp$csi[temp$csi < breaks[i]]
      temp$csi.mid[temp$csi > breaks[i] & temp$csi<breaks[j]] <- temp$csi[temp$csi > breaks[i] & temp$csi < breaks[j]]
      temp$csi.up[temp$csi > breaks[j]] <- temp$csi[temp$csi > breaks[j]]
      
      temp$above <- temp$csi > breaks[j]
      temp$middle <- temp$csi > breaks[i] & temp$csi < breaks[j]
      temp$below <- temp$csi < breaks[i]
      
      fit <- gam(ric ~ below * csi.low + middle * csi.mid + above * csi.up + s(Lati, Longi, k = 40), 
                 data = temp)
      
      # selection of the best model based on explained deviance
      val <- cbind(breaks[i], breaks[j], paste(breaks[i], breaks[j], sep="_"), summary(fit)[14])
      colnames(val) <- c("BP1", "BP2", "id", "dev.expl")
      mse <- rbind(mse, val)
    }}
  mse$dev.expl <- as.numeric(mse$dev.expl)  
  return(unlist((mse$id[which(mse$dev.expl == max(mse$dev.expl))])))
}



#############################################################################
################################################################################
# load data
biod_tt <- read.table("Community_data.txt", header=T)
env <- read.table("env.txt", header=T)

selec <- unique(env$rte_id)

#############################################################################
################################################################################
# data preparation
# prepare richness and specialization variables at the decadal scale 

b1 <- biod_tt[biod_tt$year > 1972 & biod_tt$year < 1984,]
b2 <- biod_tt[biod_tt$year > 1984 & biod_tt$year < 1996,]
b3 <- biod_tt[biod_tt$year > 2000 & biod_tt$year < 2012,]

m1 <- aggregate(b1[,c(2:3, 8:9)], by = list(b1$rte_id), FUN = mean)
m2 <- aggregate(b2[,c(2:3, 8:9)], by = list(b2$rte_id), FUN = mean)
m3 <- aggregate(b3[,c(2:3, 8:9)], by = list(b3$rte_id), FUN = mean)

biod1 <- m1[is.element(m1$Group.1,selec),]
biod2 <- m2[is.element(m2$Group.1,selec),]
biod3 <- m3[is.element(m3$Group.1,selec),]

colnames(biod1) <- colnames(biod2) <- colnames(biod3) <- c("rte_id", "csi", "ric", "Lati", "Longi")
biod1$year <- "1"
biod2$year <- "2"
biod3$year <- "3"

biod_tab <- rbind(biod1, biod2, biod3)

#write.table(biod_tab, "Community_data_decades.txt")

#############################################################################
################################################################################
# data analyses
id <- unique(biod_tab$year)

setwd("E:/postdoc _aarhus/North_america/article/paper_trends/GEB/final_submission/data")

biod_tab <- read.table("Community_data_decades.txt")

####################################
# linear models

res_linear <- as.data.frame(matrix(nrow = 0, ncol = 5))
colnames(res_linear) <- c("year","estimate_lin", "se", "p", "aic_lin")

for (i in id) {
  temp <- biod_tab[biod_tab$year == i,]
  mod1 <- gam(ric ~ csi+ s(Lati, Longi, k = 40) , data = temp, family = gaussian)

  sum1 <- summary(mod1)

  vec <- as.data.frame(cbind(as.character(i), sum1$p.coeff[2], sum1$se[2], sum1$p.pv[2], mod1$aic))
  
  colnames(vec) <- c("year","estimate_lin", "se", "p", "aic_lin")
  res_linear <- rbind(res_linear, vec)
}


####################################
# one-breakpoint piecewise models

breaks <- seq(0.29, 0.39, by = 0.001)
bp_1BP <- sapply(id, function(x) simple.bk(biod_tab, breaks, x))   # 0.332 0.301 0.323 

####
# computing the slope of the best one-breakpoint gam piecewise models  

bp <- as.data.frame(bp_1BP)
bp$year <- id

resBP <- as.data.frame(matrix(nrow = 0, ncol = 10))
colnames(resBP) <- c("interc1", "interc2", "sl1", "sl2", "sd_interc1", "sd_interc2", "sd_sl1", "sd_sl2", "AIC", "BIC")
val_plot <- as.data.frame(matrix(nrow = 0, ncol = 5))
colnames(val_plot) <- c("csi_reg","ric_pred", "upr","lwr", "year")

for (i in bp$year){
  temp <- biod_tab[biod_tab$year == i,]
  loc <- bp[bp$year == i,]
  BP1 <- loc[,1]
  
  fit <- gam(ric ~ csi * I(csi > BP1) + s(Lati, Longi, k=40), data = temp)
  sfit <- summary(fit)
  
  tt <- visreg(fit, plot = F)
  prediction <- as.data.frame(cbind(tt[[1]]$fit[,1], tt[[1]]$fit[,7:9]))
  colnames(prediction) <- c("csi_reg", "ric_pred", "upr", "lwr")
  prediction$year <- i
  
  val_plot <- rbind(val_plot, prediction)
  
  interc2 <- fit$coefficients[1] + fit$coefficients[4]
  interc1 <- fit$coefficients[1] + fit$coefficients[2]
  sl2 <- fit$coefficients[5]
  sl1 <- fit$coefficients[3] 
  
  sd_interc2 <- sfit$se[1] + sfit$se[4]
  sd_interc1 <- sfit$se[1] + sfit$se[2]
  sd_sl2 <- sfit$se[5]
  sd_sl1 <- sfit$se[3] 
  
  var <- as.data.frame(cbind(interc1, interc2, sl1, sl2, sd_interc1, 
                           sd_interc2, sd_sl1, sd_sl2, AIC(fit), BIC(fit)))
  rownames(var) <- i
  resBP <- rbind(resBP, var)
}


####################################
# two-breakpoint piecewise models

breaks <- seq(0.29, 0.39, by = 0.001)
bp_double2 <- sapply(id, function(x)double.bk(biod_tab, breaks, x))  # 

# "0.337_0.338" "0.337_0.35"  "0.299_0.3"  

BP1 <- as.data.frame(c(0.337, 0.337, 0.299))
BP2 <- as.data.frame(c(0.338, 0.35, 0.3))
test <- cbind(BP1,BP2)
colnames(test) <- c("BP1", "BP2")
test$period <- c(1, 2, 3)

id <- c(1, 2, 3)

####
# compute the three slopes of the three parts of the piecewise regression with 2 breakpoints for the best model

res_dbBP <- as.data.frame(matrix(nrow = 0, ncol = 14))
colnames(res_dbBP) <- c("interc1", "interc2", "interc3", "sl1", "sl2", "sl3", "sd_interc1",
                      "sd_interc2", "sd_interc3", "sd_sl1", "sd_sl2", "sd_sl3", "AIC", "BIC")

for (i in id)
{
  temp <- biod_tab[biod_tab$year == i,]
  bp1 <- BP1[i,]
  bp2 <- BP2[i,]
  
  temp$csi.mid <- temp$csi.low <- temp$csi.up <- 0
  temp$middle <- temp$above <- temp$below <- 0
  
  temp$csi.low[temp$csi < bp1] <- temp$csi[temp$csi < bp1]
  temp$csi.mid[temp$csi > bp1 & temp$csi < bp2] <- temp$csi[temp$csi > bp1 & temp$csi < bp2]
  temp$csi.up[temp$csi > bp2] <- temp$csi[temp$csi > bp2]
  
  temp$above <- temp$csi > bp2
  temp$middle <- temp$csi > bp1 & temp$csi < bp2
  temp$below <- temp$csi < bp1
  
  fit <- gam(ric ~ below * csi.low + middle * csi.mid + above * csi.up + s(Lati, Longi, k = 40), data = temp)
  sfit <- summary(fit)
  
  interc3 <- fit$coefficients[1] + fit$coefficients[6]
  interc2 <- fit$coefficients[1] + fit$coefficients[4]
  interc1 <- fit$coefficients[1] + fit$coefficients[2]
  sl1 <- fit$coefficients[8] + fit$coefficients[3] 
  sl2 <- fit$coefficients[9] + fit$coefficients[5]  
  sl3 <- fit$coefficients[7] + fit$coefficients[10] 
  
  sd_interc3 <- sfit$se[1] + sfit$se[6]
  sd_interc2 <- sfit$se[1] + sfit$se[4]
  sd_interc1 <- sfit$se[1] + sfit$se[2]
  sd_sl1 <- sfit$se[8] + sfit$se[3] 
  sd_sl2 <- sfit$se[9] + sfit$se[5] 
  sd_sl3 <- sfit$se[7] + sfit$se[10] 
  
  var <- as.data.frame(cbind(interc1, interc2, interc3, sl1, sl2, sl3, sd_interc1, 
                       sd_interc2, sd_interc3, sd_sl1, sd_sl2, sd_sl3, AIC(fit), BIC(fit)))
  rownames(var) <- i
  res_dbBP <- rbind(res_dbBP,var)
}

res <- res_dbBP[order(rownames(res_dbBP)),]
res$year <- as.factor(rownames(res))
res <- merge(res, bp, by = "year")


####
# comparing the AIC from the 3 models (linear models and one-breakpoint and two-breakpoint piecewise models)

comp_aic <- as.data.frame(cbind(as.vector(res_linear$aic_lin), resBP$V9, res_dbBP$V13))
colnames(comp_aic) <- c('AIC_linear', 'AIC_1BP', 'AIC_2BP')

# best models:
# decade 1 : 2 BP
# decade 2 : 2 BP
# decade 3 : 2 BP

#############################################################################
################################################################################
# figure
# producing the plot for the 2 breakpoints piecewise models at the decade time scale 

xlabel <- "Specialization"
ylabel <- "Richness"

val_plot$Period <- as.character(val_plot$year)
val_plot$Period[val_plot$Period == "1"] <- "1973 - 1983"
val_plot$Period[val_plot$Period == "2"] <- "1985 - 1995"
val_plot$Period[val_plot$Period == "3"] <- "2001 - 2011"

tiff("GAMs_decades.tiff", units = "cm", width = 15, height = 12, res = 600)
p <- ggplot(data = val_plot, aes(x = csi_reg, y = ric_pred, ymin = 35, ymax = 55, fill = Period)) + 
  geom_line(aes(x = csi_reg, y = ric_pred, color = Period)) +
  xlab(xlabel) + 
  ylab(ylabel) +
  theme(axis.title = element_text(color="black", face="bold", size = 12),
        axis.title.x = element_text(margin = margin(10, 10, 10, 10)),
        axis.title.y = element_text(margin = margin(10, 10, 10, 10)),
        axis.text  = element_text(vjust = 0.5, size = 12)) +
  scale_fill_discrete(name = "Period") +
  geom_ribbon(data = val_plot, aes(ymin = lwr, ymax = upr, fill = year), alpha = 0.2)
p

dev.off()

