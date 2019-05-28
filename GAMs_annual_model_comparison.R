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
               data=temp, family=poisson)
    
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
  mse <- as.data.frame(matrix(nrow=0, ncol=4))
  colnames(mse) <- c("BP1", "BP2", "id", "dev.expl")
  temp <- tab[tab$year==year,]
  for (i in 1:length(breaks))
  {
    #B2<-seq(breaks[i],0.37, by=0.001)  
    for(j in (i):length(breaks))
    {
      temp$csi.mid <- temp$csi.low <- temp$csi.up <- 0
      temp$middle <- temp$above <- temp$below <- 0
      
      temp$csi.low[temp$csi < breaks[i]] <- temp$csi[temp$csi < breaks[i]]
      temp$csi.mid[temp$csi > breaks[i] & temp$csi<breaks[j]] <- temp$csi[temp$csi > breaks[i] & temp$csi<breaks[j]]
      temp$csi.up[temp$csi > breaks[j]] <- temp$csi[temp$csi > breaks[j]]
      
      temp$above <- temp$csi > breaks[j]
      temp$middle <- temp$csi > breaks[i] & temp$csi < breaks[j]
      temp$below <- temp$csi < breaks[i]
      
      fit <- gam(ric ~ below*csi.low + middle* csi.mid + above*csi.up + s(Lati, Longi, k=40), 
                 data=temp, family=poisson)
      
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
d_comm <- read.table("S:/2bp/Community_data.txt", header=T)
id <- unique(d_comm$year)

####################################
# linear models

res_linear <- as.data.frame(matrix(nrow=0, ncol=6))
colnames(res_linear) <- c("year","estimate_lin", "se", "p", "aic_lin", "bic_lin")

for (i in id) {
  temp <- d_comm[d_comm$year == i,]
  mod1 <- gam(ric ~ csi + s(Lati, Longi, k=40) , data=temp, family=poisson)

  sum1 <- summary(mod1)

  vec <- as.data.frame(cbind(as.character(i), sum1$p.coeff[2], sum1$se[2], sum1$p.pv[2], mod1$aic, BIC(mod1)))
  
  colnames(vec) <- c("year","estimate_lin", "se", "p", "aic_lin", "bic_lin")
  res_linear <- rbind(res_linear, vec)
}


####################################
# one-breakpoint piecewise models

breaks <- seq(0.29,0.39, by=0.001)
bp_1BP <- sapply(id, function(x) simple.bk(d_comm, breaks,x))  

####
# computing the slope of the best one-breakpoint gam piecewise models  

bp <- as.data.frame(bp_1BP)
bp$year <- id

resBP <- as.data.frame(matrix(nrow=0, ncol=10))
colnames(resBP) <- c("interc1","interc2","sl1","sl2","sd_interc1", "sd_interc2","sd_sl1","sd_sl2", "AIC", "BIC")
val_plot <- as.data.frame(matrix(nrow=0, ncol=5))
colnames(val_plot) <- c("csi_reg","ric_pred", "upr","lwr", "year")

for (i in bp$year){
  temp <- d_comm[d_comm$year==i,]
  loc <- bp[bp$year==i,]
  BP1 <- loc[,1]
  
  fit <- gam(ric ~ csi * I(csi > BP1) + s(Lati, Longi, k=40), data=temp, family=poisson)
  
  sfit <- summary(fit)
  
  tt <- visreg(fit, plot=F)
  prediction <- as.data.frame(cbind(tt[[1]]$fit[,1], tt[[1]]$fit[,7:9] ))
  colnames(prediction) <- c("csi_reg", "ric_pred", "upr","lwr")
  prediction$year <- i
  
  val_plot <- rbind(val_plot,prediction)
  
  interc2 <- fit$coefficients[1] + fit$coefficients[3]
  interc1 <- fit$coefficients[1]
  sl1 <- fit$coefficients[2] 
  sl2 <- fit$coefficients[2] + fit$coefficients[4] 
  
  sd_interc2 <- sfit$se[1] + sfit$se[3]
  sd_interc1 <- sfit$se[1]
  sd_sl1 <- sfit$se[2] 
  sd_sl2 <- sfit$se[2] + sfit$se[4] 
  
  var <- as.data.frame(cbind(interc1, interc2, sl1, sl2, sd_interc1, 
                           sd_interc2, sd_sl1, sd_sl2, AIC(fit), BIC(fit)))
  rownames(var) <- i
  resBP <- rbind(resBP, var)
}

write.table(resBP,"S:/2bp/slopes_1BP.txt")


####################################
# two-breakpoint piecewise models

breaks <- seq(0.29,0.39, by=0.001)

####
# parallelisation

sfInit(cpus=25, parallel=T)
sfExport(list=c("d_comm", "breaks", "double.bk", "id"))
sfLibrary(mgcv)

t1 <- Sys.time()
bk.pts <- sfSapply(id, function(x) double.bk(d_comm, breaks, x))
t2<-Sys.time()
sfStop()

bp <- as.data.frame(cbind(id, bk.pts))
colnames(bp) <- c("year", 'bk.pts')

bp <- bp[order(bp$year),]
l <- strsplit(as.character(bp$bk.pts), split="_")
l2 <- as.data.frame(do.call(rbind, l))
bp <- cbind(bp, l2)
bp$V1 <- as.numeric(as.character(bp$V1))
bp$V2 <- as.numeric(as.character(bp$V2))
bp$plat <- bp$V2 - bp$V1
bp$year <- as.factor(bp$year)

####
# compute the three slopes of the three parts of the piecewise regression with 2 breakpoints for the best model

res_dbBP <- as.data.frame(matrix(nrow=0, ncol=14))
colnames(res_dbBP) <- c("interc1","interc2","interc3","sl1","sl2","sl3","sd_interc1",
                      "sd_interc2","sd_interc3","sd_sl1","sd_sl2","sd_sl3", "AIC", "BIC")

for (i in id)
{
  temp <- d_comm[d_comm$year==i,]
  loc <- bp[bp$year==i,]
  BP1 <- loc$V1
  BP2 <- loc$V2
  
  temp$csi.mid <- temp$csi.low<-temp$csi.up <- 0
  temp$middle <- temp$above <- temp$below <- 0
  
  temp$csi.low[temp$csi < BP1] <- temp$csi[temp$csi < BP1]
  temp$csi.mid[temp$csi > BP1 & temp$csi < BP2] <- temp$csi[temp$csi > BP1 & temp$csi < BP2]
  temp$csi.up[temp$csi > BP2] <- temp$csi[temp$csi > BP2]
  
  temp$above <- temp$csi > BP2
  temp$middle <- temp$csi > BP1 & temp$csi < BP2
  temp$below <- temp$csi < BP1
  
  fit <- gam(ric ~ below * csi.low + middle * csi.mid + above * csi.up + s(Lati, Longi, k=40), 
             data=temp, family=poisson)
  
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
  res_dbBP <- rbind(res_dbBP, var)
}

res <- res_dbBP[order(rownames(res_dbBP)),]
res$year <- as.factor(rownames(res))
res <- merge(res, bp, by="year")

write.table(res, "S:/2bp/slopes_2BP.txt")

####
# comparing the AIC from the 3 models (linear models and one-breakpoint and two-breakpoint piecewise models)

comp_aic <- as.data.frame(cbind(as.vector(res_linear$aic_lin), resBP$V9, res_dbBP$V13))
colnames(comp_aic) <- c('AIC_linear', 'AIC_1BP', 'AIC_2BP')

comp_aic$best_mod <- colnames(comp_aic)[apply(comp_aic, 1, which.min)] # model with 1 breakpoints (for all years)
comp_aic$year <- id
comp_aic <- comp_aic[order(comp_aic$year),]

# write.table(comp_aic, "comp_mod_annual.txt")

#############################################################################
################################################################################
# producing plots for the one breakpoint piecewise models 

res_plot <- res_dbBP 
res_plot$year <- as.numeric(rownames(res_plot))
res_plot <- merge(res_plot, bp, by="year")

####
# Change in the location of the first breakpoint

tiff("loc_bp1.tiff", units="cm", width=15, height=12, res=600)
p<-ggplot(data=res_plot, aes(x=as.numeric(res_plot$year), y=as.numeric(res_plot$V1))) +
  labs(x="Year", y= "First breakpoint location") +
  xlim(1965,2015)+
  #ylim(0.28,0.38)+
  theme(axis.title = element_text(color="black", face="bold", size=14),
        axis.title.x = element_text(margin=margin(10,10,10,10)),
        axis.title.y = element_text(margin=margin(10,10,10,10)),
        axis.text  = element_text(vjust=0.5, size=12))


fit<-lm(res_plot$V1 ~ res_plot$year)

p <- p + geom_point(size=2) +
  geom_abline(intercept = fit$coefficients[1], slope = fit$coefficients[2]) 

p<- p +  stat_cor(method = "pearson", label.x = 1995, label.y = 0.37, size=5)  # Add correlation coefficient
print(p)

dev.off()

####
# Change in the location of the second breakpoint

tiff("loc_bp2.tiff", units="cm", width=15, height=12, res=600)
p<-ggplot(data=res_plot, aes(x=as.numeric(res_plot$year), y=as.numeric(res_plot$V2))) +
  labs(x="Year", y= "First breakpoint location") +
  xlim(1965,2015)+
  #ylim(0.28,0.38)+
  theme(axis.title = element_text(color="black", face="bold", size=14),
        axis.title.x = element_text(margin=margin(10,10,10,10)),
        axis.title.y = element_text(margin=margin(10,10,10,10)),
        axis.text  = element_text(vjust=0.5, size=12))


fit<-lm(res_plot$V2 ~ res_plot$year)

p <- p + geom_point(size=2) +
  geom_abline(intercept = fit$coefficients[1], slope = fit$coefficients[2]) 

p<- p +  stat_cor(method = "pearson", label.x = 1995, label.y = 0.37, size=5)  # Add correlation coefficient
print(p)

dev.off()


####
# change in the first slope 

tiff("slope1.tiff", units="cm", width=15, height=12, res=600)
p<-ggplot(data=res_plot, aes(x=as.numeric(res_plot$year), y=as.numeric(res_plot$sl1))) +
  labs(x="Year", y= "Slope part 1") +
  xlim(1965,2015)+
  ylim(0,10)+
  
  theme(axis.title = element_text(color="black", face="bold", size=14),
        axis.title.x =element_text(margin=margin(10,10,10,10)),
        axis.title.y =element_text(margin=margin(10,10,10,10)),
        axis.text  = element_text(vjust=0.5, size=12))

fit<-lm(res_plot$sl1~ res_plot$year, weights=1/res_plot$sd_sl1^2)

p <- p + geom_point(size=2) +  geom_errorbar(aes(ymin=sl1-sd_sl1, ymax=sl1+sd_sl1), width=.1) +
  geom_abline(intercept = fit$coefficients[1], slope = fit$coefficients[2]) 
p<- p +  stat_cor(method = "pearson", label.x = 1996, label.y = 9, size=5)  # Add correlation coefficient
print(p)

dev.off()


####
#slope 2

tiff("slope2.tiff", units="cm", width=15, height=12, res=600)
p<-ggplot(data=res_plot, aes(x=as.numeric(res_plot$year), y=as.numeric(res_plot$sl2))) +
  labs(x="Year", y= "Slope part 3") +
  ylim(-1000,1000) +
  xlim(1965,2015)+
  #geom_abline(intercept = 0, slope = 0) +
  theme(axis.title = element_text(color = "black", face = "bold", size = 14),
        axis.title.x = element_text(margin = margin(10,10,10,10)),
        axis.title.y = element_text(margin = margin(10,10,10,10)),
        axis.text  = element_text(vjust = 0.5, size=12))

fit<-lm(res_plot$sl2 ~ res_plot$year, weights= 1 / res_plot$sd_sl2^2)

p <- p + geom_point(size=2) +  geom_errorbar(aes(ymin = sl2 - sd_sl2, ymax = sl2 + sd_sl2), width =.1) +
  geom_abline(intercept = fit$coefficients[1], slope = fit$coefficients[2]) 
p<- p +  stat_cor(method = "pearson", label.x = 1996, label.y = 800, size = 5)  # Add correlation coefficient
print(p)

ggsave("S:/fig_1bp/slope2.png")

dev.off()

####
#slope 3

tiff("slope3.tiff", units="cm", width=15, height=12, res=600)
p<-ggplot(data=res_plot, aes(x=as.numeric(res_plot$year), y=as.numeric(res_plot$sl3))) +
  labs(x="Year", y= "Slope part 3") +
  #ylim(-1000,1000) +
  xlim(1965,2015)+
  #geom_abline(intercept = 0, slope = 0) +
  theme(axis.title = element_text(color = "black", face = "bold", size = 14),
        axis.title.x = element_text(margin = margin(10,10,10,10)),
        axis.title.y = element_text(margin = margin(10,10,10,10)),
        axis.text  = element_text(vjust = 0.5, size=12))

fit<-lm(res_plot$sl3 ~ res_plot$year, weights= 1 / res_plot$sd_sl3^2)

p <- p + geom_point(size=2) +  geom_errorbar(aes(ymin = sl3 - sd_sl3, ymax = sl3 + sd_sl3), width =.1) +
  geom_abline(intercept = fit$coefficients[1], slope = fit$coefficients[2]) 
p<- p +  stat_cor(method = "pearson", label.x = 1996, label.y = 2, size = 5)  # Add correlation coefficient
print(p)

dev.off()


####
# 
# graph pour 1971 and 2014 (values closest to the regression for the significant breakpoint and slopes)
library(visreg)
library(mgcv)

#1971
temp<-tab_ssi3[tab_ssi3$year==1971,]
loc<-bp[bp$year==1971,]
BP1<-loc$V1
BP2<-loc$V2

fit <-gam(ric~csi*I(csi<BP1)+ csi*I(csi>BP1& csi<BP2)+ s(Lati,Longi,k=40), data=temp, family=poisson)

tt<-visreg(fit, plot=F, scale='response')
prediction<- as.data.frame(cbind(tt[[1]]$fit[,1],tt[[1]]$fit[,9:11] ))
colnames(prediction)<-c("csi_reg","ric_pred", "upr","lwr")

library(ggplot2)
xlabel <- "Specialization"
ylabel <- "Richness"


#2014
temp2<-tab_ssi3[tab_ssi3$year==2014,]
loc<-bp[bp$year==2014,]
BP1<-loc$V1
BP2<-loc$V2

fit <-gam(ric~csi*I(csi<BP1)+ csi*I(csi>BP1& csi<BP2)+ s(Lati,Longi,k=40), data=temp2, family=poisson)

tt2<-visreg(fit, plot=F, scale='response', by=)
prediction2<- as.data.frame(cbind(tt2[[1]]$fit[,1],tt2[[1]]$fit[,9:11] ))
colnames(prediction2)<-c("csi_reg","ric_pred", "upr","lwr")

##############
library(ggplot2)
xlabel <- "Specialization"
ylabel <- "Richness"

prediction$Year<-"1971"
prediction2$Year<-"2014"

pred<-rbind(prediction, prediction2)
pred75_1<-pred[is.element(pred$Year,1971) & pred$csi_reg < 0.334,]
pred75_2<-pred[is.element(pred$Year,1971) & pred$csi_reg >= 0.334 & pred$csi_reg <=0.349 ,]
pred75_3<-pred[is.element(pred$Year,1971) & pred$csi_reg > 0.349 ,]

pred14_1<-pred[is.element(pred$Year,2014) & pred$csi_reg < 0.296,]
pred14_2<-pred[is.element(pred$Year,2014) & pred$csi_reg >= 0.296 & pred$csi_reg <=0.360 ,]
pred14_3<-pred[is.element(pred$Year,2014) & pred$csi_reg > 0.36 ,]


p <- ggplot(data=pred, aes(x=csi_reg, y=ric_pred, ymin=35, ymax=60, fill = Year)) + 
  geom_line(data=pred75_1,aes(x=csi_reg, y=ric_pred, color = Year)) +
  geom_line(data=pred75_2,aes(x=csi_reg, y=ric_pred, color = Year)) +
  geom_line(data=pred75_3,aes(x=csi_reg, y=ric_pred, color = Year)) +
  geom_line(data=pred14_1,aes(x=csi_reg, y=ric_pred, color = Year)) +
  geom_line(data=pred14_2,aes(x=csi_reg, y=ric_pred, color = Year)) +
  geom_line(data=pred14_3,aes(x=csi_reg, y=ric_pred, color = Year)) +  
  
  scale_fill_manual(values=c("#99CC66", "#FF6666"))+
  scale_color_manual(values=c("#99CC66", "#FF6666"), guide=F)+
  xlab(xlabel) + 
  ylab(ylabel) +
  geom_ribbon(data=pred,aes(ymin=lwr,ymax=upr, fill=Year),alpha=0.2)+
  theme(axis.title = element_text(color="black", face="bold", size=14),
        axis.title.x =element_text(margin=margin(10,10,10,10)),
        axis.title.y =element_text(margin=margin(10,10,10,10)),
        axis.text  = element_text(vjust=0.5, size=12),legend.position="bottom",
        legend.title=element_text(size=14),
        legend.text=element_text(size=14))
p 



ggsave("S:/relationship_1970_2014_ok.png")

