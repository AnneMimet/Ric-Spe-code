#############################################################################
################################################################################
# load packages

library(lme4)
library(nlme)
library(MASS)
library(ncf)
library(mgcv)
library(ggplot2)
library(car)
library(doParallel)
library(foreach)
library(MASS)
library(doSNOW)
library(Rmisc)

##########################################################
#########################################################
# GLMM with glmmPQL on the richness and specialization data
# load the data
d_comm <- read.table("Community_data.txt", header=T)

memory.size(20000000)

# glmmPQL mixed model

RS_mod <- glmmPQL(ric ~ csi + as.factor(year) + csi*as.factor(year) , random= ~1 | rte_id,
                 data=d_comm, family=poisson, niter=10)

# first look at the residuals
AT <- acf(resid(RS_mod, type="normalized")) 
plot(RS_mod)
plot(resid(RS_mod))

####
# testing residual spatial autocorrelation spline.correlog function, for each year, and print the correponding figure per year
# create a data frame with the needed variables
data <- cbind(d_comm[,8:9], RS_mod$fitted[,2], d_comm$year, resid(RS_mod, type="normalized"))
colnames(data) <- c("Lati", "Longi", "fitted", "year", "resid_norm")
time <- as.factor(sort(unique(as.numeric(d_comm$year))))

res <- data.frame(matrix(nrow=50, ncol=2))
colnames(res) <- c("year", "moy_sig")
res$year <- time

# build a pdf with the correlogram for each year + extract 
pdf("moran_glmm.pdf")
for(i in time)
{
  d <- data[data$year==i,]
  cor <- spline.correlog(x=d$Longi, y=d$Lati, z=d$resid_norm, latlon=T,resamp=2, xmax=400)
  title <- as.character(paste(i, min(cor$p), sep="_"))
    res[res$year==i,] <- cbind(i,mean(cor$real$predicted$y),sd(cor$real$predicted$y))
  
  # DF with correlation value, upper and lower CI (2.5 and 97.5% quantile) 
  # of the resampling distribution, and x values)
  df1 <- data.frame(lower_ci = cor$boot$boot.summary$predicted$y[2,],
                    upper_ci = cor$boot$boot.summary$predicted$y[10,],
                    pred_y = cor$real$predicted$y,
                    lags = cor$real$predicted$x)
  
  # Significativity of the correlation (logical, TRUE if significant)
  df1$signif <- sign(df1$lower_ci) == sign(df1$upper_ci)
  
  # Maximal significant correlation
  res[res$year==i,] <- cbind(i,sum(df1$pred_y[df1$signif == TRUE])/300)
  
  #plot residual autocorrelation
  plot(cor, main=titre, ylim=c(-1,1))
}
dev.off()


sum_glmm<-summary(RS_mod)
p_val<-sum_glmm$tTable[,"p-value"]
results<-sum_glmm$tTable

coef_temp<-as.data.frame(RS_mod$coefficients$fixed[52:100])
coef_temp$year<-c(1967:2015)

#####
# extracting coefficients for the final figure 
# function recomputing the coefficients of the specialization relatively to zero instead of the first year of the study

marginal.inter.pql <- function(mod, var)
{
  require(car)
  coef.sar <- mod$coefficients$fixed
  cnames <- names(coef.sar)
  
  names(coef.sar) <- paste("b", c(0:(length(coef.sar)-1)), sep="")
  
  #Extract the specific parameters of interest
  vars <- grep(var,cnames,value=T)
  
  names(coef.sar) <- paste("b", c(0:(length(coef.sar)-1)), sep="")
  
  #--Create Data Frame to store Main Effect
  var1 <- names(coef.sar)[which(cnames==vars[1])]
  
  interac <- deltaMethod(coef.sar,g=var1,vcov=vcov(mod)[cnames, cnames])
  
  for (i in c(2:length(vars)))
  {
    interac <- rbind(interac, 
                     deltaMethod(coef.sar,g=paste(var1, 
                    names(coef.sar)[which(cnames==vars[i])], sep="+"),vcov=vcov(mod)[cnames, cnames]))
  }
  
  rownames(interac) <- vars
  interac$tval <- interac$Estimate/interac$SE
  
  interac$pval <- 2*pnorm(abs(interac$tval), lower.tail=F)
  return(interac)
}

# computing the coeficients with reference to 0 for RS_mod
mod <- RS_mod
var <- "csi"

ref0 <- marginal.inter.pql(RS_mod, "csi")

coef_temp <- ref0[1:49,1]
coef_temp <- as.data.frame(coef_temp)
coef_temp$year <- c(1967:2015)
colnames(coef_temp) <- c("coef", "year")
coef_temp$cond <- "1"
coef_temp$se <- ref0[1:49,2]


##########################################################
#########################################################
# Building null models on randomized data to compare the slope  of the observed richness -specialization
# relationship with random slopes

# randomize the species specialization (SSI) 200 times and 
# compute 200 random community specialization (CSI) based on  randomized SSIs
mat_csi <- list()
mat_csi <- lapply (1:200, function(i) sample(as.matrix(ssi[,2])))

random_csi <- as.data.frame(matrix(nrow=nrow(selec), ncol=200))
colnames(random_csi) <- seq(1:200)
rownames(random_csi) <- selec$id_rte_year

for (j in seq(1:200))
{
  mult <- sweep(selec[,1:162], MARGIN=2,as.matrix(mat_csi[[j]]),`*`)
  random_csi[,j] <- rowSums(mult)/rowSums(selec[,1:162])
}

random_csi$ric <- selec$ric


# loop running a model per random distribution and extracting the coefficients 
random_csi[,1:200] <- random_csi[,1:200] / 100000

cl <- makeCluster(20)
registerDoSNOW(cl)

mods_random <- foreach (i=seq(200), .packages = "MASS", .export = 'random_csi')%dopar% {
  print(i)
  
  tab <- random_csi[,c(i, 1000:1002)]
  colnames(tab) <- c("csi", "ric", "year", "rte_id")
  
  AT_random <- glmmPQL(ric ~ csi + as.factor(year) + csi *as.factor(year) , random= ~1 | rte_id,
                       data=tab, family=poisson, niter=10)
  
}

stopCluster(cl)


###################################
# get all random values of the richness - specialization relationship per year and compute the confidence interval 
# over the 200 values


#########
# extract specialization coefficients of the random models relatively to 0

coef <- list()

for (i in seq(1:200)){
    mod <- mods_random[[i]]
  var <- "csi"
  ref0<-marginal.inter.pql(mod, "csi")
  
  coef_temp<-ref0[1:49,1]
  coef_temp<-as.data.frame(coef_temp)
  coef_temp$year<-c(1967:2015)
  colnames(coef_temp)<-c("coef", "year")
  coef_temp$cond<-"1"
  coef_temp$se<-ref0[1:49,2]
  coef[[i]] <- coef_temp
}

#########
# extract the mean and CI of the 200 coefficients from the 200 random models

list_coef <- lapply(coef, function(x) x[[1]])
tab_coef <- t(do.call(rbind, list_coef))
tab_coef <- as.data.frame(tab_coef)
tab_coef$year <- coef[[1]][,2]

ci_coef <- t(apply(tab_coef[,1:200], 1, CI))
ci_coef <- as.data.frame(ci_coef)
ci_coef$year <- tab_coef$year

write.table(ci_coef,"ci_coef_null.txt")


##########################################################
#########################################################
# final figure: glmmPQL coefficients with reference to 0 for each year and for the mean and CI from the random models

#############################
########################################
# load the file with the results of the 200 random models (coef already computed relatively to 0)
nulm <- read.table("ci_coef_null.txt")

# plot

p<-ggplot(data=coef_temp, aes(x=coef_temp$year, y=coef_temp$coef)) +
  xlab("Year")+
  ylab(expression(atop("Yearly slope of the", paste("richness-specialization relationship")))) +
  scale_y_continuous(limits=c(-0.1,1.1), breaks=seq(-0.1, 1.1, by = 0.1))+
  theme(axis.title = element_text(color="black", face="bold", size=13),
        axis.title.x =element_text(margin=margin(10,10,10,10)),
        axis.title.y =element_text(margin=margin(10,10,10,10), face="bold"),
        axis.text  = element_text(vjust=0.5, size=11))

p <- p + stat_smooth(color="cyan4", se=T)

limits <- aes(ymax = coef + se, ymin=coef - se)
p <- p + geom_point(size=2, color="cyan4")+ geom_errorbar(limits, width=0, color="cyan4")

p <- p + stat_smooth(data= nulm, aes(x=coef_temp$year, y=nulm$mean), color="black", se=F)
p <- p + geom_ribbon(aes(ymin=nulm$lower, ymax=nulm$upper), alpha=0.15)

print(p)  

ggsave("final_glmm_figure.png", width= 12, height= 12,unit="cm",dpi=600)
