rm(list=ls())
require(dagitty)
require(dplyr)
require(ggplot2)
require(ggdist)
require(tidyr)
require(forcats)

# Data coding issues ------------------------------------------------------
#' This section re-creates figure 1 from the commentary, and questions or
#' concerns about this section can be directed to adlightner at cas.au.dk
load('data/MSP_TableData.Rdata')
d <- AggrDat

rbb_densityPlot <- 
  d %>%
  mutate(MSP_sum=MSP_sum/max(d$MSP_sum, na.rm=TRUE)) %>%   # scaling additive values between 0-1
  rename(Multiplicative=MSP, Additive=MSP_sum) %>% 
  pivot_longer(c(Multiplicative, Additive), names_to='MSP_type', values_to='MSP') %>% 
  ggplot(aes(x=Time, y=MSP)) +
  theme_classic(base_size=18) +
  stat_density_2d(aes(fill = stat(density)), geom='raster', contour=FALSE) +
  geom_density_2d(colour='white', alpha=0.25) +
  scale_fill_viridis_c(name = "density") +
  geom_point(alpha=0.25, colour='white', size=3) +
  facet_wrap(~fct_rev(MSP_type)) +
  xlim(c(-5500,2000)) +
  labs(x='\nYear BCE/CE', y='MSP\n') +
  theme(legend.position = 'none')

d$infA2 <- NA
d$infA2[d$infA==1] <- 'Inferred'
d$infA2[d$infA==0] <- 'None inferred'

rbb_mspAbsences <- 
  d %>% 
  filter(!is.na(minMSP)) %>% 
  mutate(minMSP = case_when(
    minMSP == 1 ~ "MSP coded\nas present",
    minMSP == 0 ~ "MSP coded\nas absent"
  )) %>% 
  ggplot(aes(x=Time, y=factor(minMSP), fill=infA2, colour=infA2)) +
  theme_classic(base_size=18) +
  #ggdist::stat_histinterval(position=position_dodge(width=0.1)) +
  ggdist::stat_dots(position=position_dodge(width=0.65), size=4, justification=0.1) +  # , alpha=0.5
  scale_fill_manual(values=c(viridis::magma(11)[4],
                             viridis::magma(11)[8])) +
  scale_colour_manual(values=c(viridis::magma(11)[4],
                               viridis::magma(11)[8])) +
  labs(x='\nYear BCE/CE', y='', fill='Absences', colour='Absences') +
  theme(legend.position=c(0.125,0.85),
        legend.background=element_rect(fill='white', colour='black')) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

msp <- round(d$MSP[!is.na(d$MSP)], digits=4)   # taking all unique MSP values directly from Seshat dataset
msp <- sort(unique(msp))                       # for the upper left panel in figure 1

# Figure 1 panels ---------------------------------------------------------

#' ADL note: this was an alternative visualization that was cut for space, which
#' shows how the multiplicative scale (blue) compares to the additive scale (black)
par(mar=c(5,5,3,5))
plot(1:length(msp), msp, type='l',   # Upper left panel in figure 1
     xlab='Number of MSP features present', 
     ylab='MSP variable',
     lwd=2, col='blue', 
     cex=2.5, cex.lab=2, cex.axis=1.5)
lines(1:length(msp), seq(0,1,length.out=length(msp)), col='black', lwd=2)

# Saves the pdf's used in the commentary
ggsave('MSP-dot-histograms.pdf', plot=rbb_mspAbsences, device='pdf',
       width=10, height=6, units='in')   # Upper right panel in figure 1

ggsave('MSP-density-add-vs-mult.pdf', plot=rbb_densityPlot, device='pdf',
       width=10, height=6, units='in')   # Lower panel in figure 1


# Data analysis issues ----------------------------------------------------

#########################################################
### *Very* simple re-analysis of key results          ###
### in "Explaining the Rise of Moralizing Religions"  ###
###  by Turchin et al. (RBB, Target Article)          ###
#########################################################

##### Theiss Bendixen, Aaron D. Lightner, & Benjamin Grant Purzycki
##### Contact: tb@cas.au.dk
##### Last updated: 10th January 2022

##### Some notes:
# The following is a mix of the authors' script and your's truly's. 
# We try to delineate throughout; e.g., "TB:" denotes comments/script by us.
# Author's code and data can be retrieved from: https://osf.io/9w2t8/

###############
### Outline ###
###############

# 1) Inspect model comparison results and highlight why the authors' approach and conclusions are invalid
# 2) Reproduce Table 3 results for collinearity/masking effects
# 3) Reproduce Table 4 results for collinearity/masking effects

############################################################
### TB: Model comparison issues: MSP as outcome variable ###
############################################################

##### TB: Prepare data (the following code is from Turchin et al., file: "2MSP_SPC1.R")

#####  Dynamic Regressions with MSP as response variable
{                          ### Construct NGARegrDat
  load("data/MSP_TableData.Rdata")
  AggrDat <- AggrDat[AggrDat$MSP_approv == 1,]  ### Include only approved values **** Option
  AggrDat$Agri.sq <- AggrDat$Agri^2
  
  TableDat <- subset(AggrDat, select = c(NGA, PolID, Time, MSP, SPC1, MilTech, Cavalry, Agri, Agri.sq, Pastor, EnvPC1, EnvPC2))
  TableDat$MSP <- log10(TableDat$MSP) 
  
  # for(i in 1:ncol(TableDat)) { print(c(colnames(TableDat)[i], sum(is.na(TableDat[,i])))) }
  dpar <- 1000*1   ### d determines how rapidly geographic influence declines with distance
  source("data/fRegrDat.R")  
  NGARegrDat <- RegrDat
  rm( list = setdiff(ls(),c("NGARegrDat") ) )
}

##### TB: The "output" produced by the following code block (by the authors) is a data frame that contains the
##### various model permutations along with model comparison metrics (R^2, relative AIC scores, 
##### as well as number of parameters, p)

{                         ### Model Selection
  RegrDat <- subset(NGARegrDat, select = -c(NGA, PolID, Time) )
  summary(fit <- glm(RegrDat))
  
  RegrDat <- subset(NGARegrDat, select = -c(NGA, PolID, Time, Space, T, Phylogeny, Lag2) ) ## Drop non-sign or weak autocorr terms
  summary(fit <- lm(RegrDat))
  
  ####  Exhaustive regressions with these terms
  print(paste("Response variable =",colnames(RegrDat)[1]))
  Predictors <- 3:ncol(RegrDat)
  output <- data.frame()
  for (nPred in 1:length(Predictors)){ print(nPred)
    Preds<- combn(Predictors, nPred)
    for(j in 1:length(Preds[1,])){
      fit <- lm(RegrDat[, c(1:2, Preds[,j])])
      Pval <- summary(fit)$coefficients[,4]
      tval <- summary(fit)$coefficients[,3]
      out <- vector("numeric",length = length(RegrDat))
      out[c(1:2,Preds[,j])] <- tval
      out <- c(out,summary(fit)$r.sq)
      fit <- glm(RegrDat[, c(1:2, Preds[,j])])
      AIC <- summary(fit)$aic
      n <- length(fit$residuals)
      p <- length(fit$coefficients)  
      out <- c(out,AIC,n,p)
      output <- rbind(output,out)
    }
  }
  colnames(output) <- c(colnames(RegrDat),"R-sq","delAIC", "n","p")
  output$delAIC <- output$delAIC - min(output$delAIC)
  write.csv(output, file="output.csv",  row.names=FALSE)
  out <- output[order(output$delAIC),]
  out <- out[out$delAIC < 3,]
  # write.table(out, "clipboard", sep="\t", row.names=FALSE)
}

##### TB: Now, view "output". 
View(output)

##### TB: 
##### There are several things to note here: First, the favored models (high R^2, low delAIC) are also the most complex 
##### (highest number of parameters, p). This is to be expected, since adding variables usually increases predictive 
##### accuracy. However, these scores have no bearing on *causal* relationships, as detailed in our commentary.
##### Second, the models where SPC1 has a positive association with the outcome, MSP, are the models that
##### exclude MilTech and Cavalry (as well as Agri), supporting the suspicion that these variables share 
##### much information and might mask each others' associations (see below). Third, these models have lower AIC scores,
##### not because SPC1 is *included* (the impression that you might get from the target article; e.g., p. 17), but
##### because they *exclude* relevant explanatory variables, which naturally decrease predictive accuracy. 

#############################################################
### TB: Model comparison issues: SPC1 as outcome variable ###
#############################################################

##### TB: Prepare data (the following code is from Turchin et al., file: "2SPC_MSP.R")

#####  Dynamic Regressions with SPC1 as response variable
{                          ### Construct NGARegrDat
  load("data/MSP_TableData.Rdata")
  AggrDat <- AggrDat[is.na(AggrDat$MSP)==F,]
  AggrDat <- AggrDat[AggrDat$MSP_approv == 1,]  ### Include only approved values **** Option
  #AggrDat$MSP <- AggrDat$minMSP ### Option: use minMSP as the potential predictor; also turn off log-transform
  #AggrDat$MSP <- AggrDat$MSP_this  ### Option
  #AggrDat$MSP <- AggrDat$MSP_after ### Option
  #AggrDat$MSP <- AggrDat$MSP_agen  ### Option
  
  TableDat <- subset(AggrDat, select = c(NGA, PolID, Time, SPC1, MilTech, Cavalry, Agri, MSP)) 
  TableDat$MSP <- log10(TableDat$MSP)  ### Option: log-transform MSP
  
  dpar <- 1000*1   ### d determines how rapidly geographic influence declines with distance
  source("data/fRegrDat.R")  
  NGARegrDat <- RegrDat
  rm( list = setdiff(ls(),c("NGARegrDat") ) )
}

##### TB: The "output" produced by the following code block (by the authors) is a data frame that contains the
##### various model permutations along with model comparison metrics (R^2, relative AIC scores, 
##### as well as number of parameters, p)

{                         ### Model Selection
  RegrDat <- subset(NGARegrDat, select = -c(NGA, PolID, Time) )
  summary(fit <- glm(RegrDat))
  
  RegrDat <- subset(NGARegrDat, select = -c(NGA, PolID, Time, Space, Phylogeny, T) ) ## Drop non-sign or weak autocorr terms
  RegrDat <- RegrDat[is.na(RegrDat$Lag2)==FALSE,]  ### Since retaining Lag2 term, drop missing observations for this var
  summary(fit <- lm(RegrDat))
  out <- summary(fit <- glm(RegrDat))
  # write.table(out$coefficients, "clipboard", sep="\t")
  # clipr::write_clip(out$coefficients, "table", "\t")  # unix-friendly alternative
  
  ####  Exhaustive regressions with these terms
  print(paste("Response variable =",colnames(RegrDat)[1]))
  Predictors <- 3:ncol(RegrDat)
  output <- data.frame()
  for (nPred in 1:length(Predictors)){ print(nPred)
    Preds<- combn(Predictors, nPred)
    for(j in 1:length(Preds[1,])){
      fit <- lm(RegrDat[, c(1:2, Preds[,j])])
      Pval <- summary(fit)$coefficients[,4]
      tval <- summary(fit)$coefficients[,3]
      out <- vector("numeric",length = length(RegrDat))
      out[c(1:2,Preds[,j])] <- tval
      out <- c(out,summary(fit)$r.sq)
      fit <- glm(RegrDat[, c(1:2, Preds[,j])])
      AIC <- summary(fit)$aic
      n <- length(fit$residuals)
      p <- length(fit$coefficients)  
      out <- c(out,AIC,n,p)
      output <- rbind(output,out)
    }
  }
  colnames(output) <- c(colnames(RegrDat),"R-sq","delAIC", "n","p")
  output$delAIC <- output$delAIC - min(output$delAIC)
  write.csv(output, file="output.csv",  row.names=FALSE)
  out <- output[order(output$delAIC),]
  out <- out[out$delAIC < 3,]
  # write.table(out, "clipboard", sep="\t", row.names=FALSE)
}

##### TB: Now, view "output". 
View(output)

##### TB: 
##### There are several things to note here, overlapping strongly with the discussion above:
##### First, the favored models (high R^2, low delAIC) are again also the most complex models
##### (highest number of parameters, p).
##### Second, the models where MSP has a positive association with the outcome, SPC1, are the models that
##### exclude MilTech and Cavalry, again supporting the suspicion that these variables share 
##### much information and might mask each others' associations (see below). Third, these models again have lower AIC scores,
##### not because MSP is *included*, but because they *exclude* relevant explanatory variables, 
##### which naturally decrease predictive accuracy.

############################################
### TB: Table 3: MSP as outcome variable ###
############################################

##### TB: Prepare data again (the following code is from Turchin et al., file: "2MSP_SPC1.R")

#####  Dynamic Regressions with MSP as response variable
{                          ### Construct NGARegrDat
  load("data/MSP_TableData.Rdata")
  AggrDat <- AggrDat[AggrDat$MSP_approv == 1,]  ### Include only approved values **** Option
  AggrDat$Agri.sq <- AggrDat$Agri^2
  
  TableDat <- subset(AggrDat, select = c(NGA, PolID, Time, MSP, SPC1, MilTech, Cavalry, Agri, Agri.sq, Pastor, EnvPC1, EnvPC2))
  TableDat$MSP <- log10(TableDat$MSP) 
  
  # for(i in 1:ncol(TableDat)) { print(c(colnames(TableDat)[i], sum(is.na(TableDat[,i])))) }
  dpar <- 1000*1   ### d determines how rapidly geographic influence declines with distance
  source("data/fRegrDat.R")  
  NGARegrDat <- RegrDat
  rm( list = setdiff(ls(),c("NGARegrDat") ) )
}

# TB: Table 3 model: this by and large recovers the point estimates reported in Table 3
RegrDat <- NGARegrDat

tab3 <- glm(MSP ~ MilTech + Cavalry + Agri + Agri.sq + Pastor + EnvPC1 + EnvPC2, data = RegrDat)
summary(tab3)

# TB: What happens when we exclude MilTech and Cavalry -- which we suspect soak up much of the variation
# that is also in SPC1 -- and include SPC1 instead? SPC1 is now highly "significant"
tab3_new <- glm(MSP ~ SPC1 + Agri + Agri.sq + Pastor + EnvPC1 + EnvPC2, data = RegrDat)
summary(tab3_new)

#############################################
### TB: Table 4: SPC1 as outcome variable ###
#############################################

##### TB: Prepare data again (the following code is from Turchin et al., file: "2SPC_MSP.R")

#####  Dynamic Regressions with SPC1 as response variable
{                          ### Construct NGARegrDat
  load("data/MSP_TableData.Rdata")
  AggrDat <- AggrDat[is.na(AggrDat$MSP)==F,]
  AggrDat <- AggrDat[AggrDat$MSP_approv == 1,]  ### Include only approved values **** Option
  #AggrDat$MSP <- AggrDat$minMSP ### Option: use minMSP as the potential predictor; also turn off log-transform
  #AggrDat$MSP <- AggrDat$MSP_this  ### Option
  #AggrDat$MSP <- AggrDat$MSP_after ### Option
  #AggrDat$MSP <- AggrDat$MSP_agen  ### Option
  
  TableDat <- subset(AggrDat, select = c(NGA, PolID, Time, SPC1, MilTech, Cavalry, Agri, MSP)) 
  TableDat$MSP <- log10(TableDat$MSP)  ### Option: log-transform MSP
  
  dpar <- 1000*1   ### d determines how rapidly geographic influence declines with distance
  source("data/fRegrDat.R")  
  NGARegrDat <- RegrDat
  rm( list = setdiff(ls(),c("NGARegrDat") ) )
}

# TB: Recovering the authors' result reported in Table 4, using their code:
RegrDat <- subset(NGARegrDat, select = -c(NGA, PolID, Time) )
summary(fit <- glm(RegrDat))

# TB: Removing Miltech and Cavalry as well as the autoregression terms, 
# which are generally difficult to interpret in causal terms,
# reveals a highly "significant" association between SPC1 and MSP. 
# Here's a modification to the authors' own code: 
tab4_new_d <- subset(RegrDat, select = -c(SPC1, SPC1.sq, MilTech, Cavalry, Lag2))
summary(fit <- lm(tab4_new_d))
(tab4_new <- summary(fit <- glm(tab4_new_d)))

###################
### TB: Summary ###
###################

##### To sum up: The authors' analytic approach has no bearing on causal inference. R^2 and AIC
##### are (imperfect!) measures of predictive accuracy, *not* causal relations and directions.
##### Further, key variables share a substantial amount of information and hence mask each others' associations
##### in multiple regression models. Given that the target article's explicit goal is causal inference and that
##### the analytical approach is not aligned with this goal, the target article's conclusions are not even wrong,
##### but simply invalid. 

# Causal inference issues -------------------------------------------------

############################################################

############################################################

## Simulation of Turchin et al. data generation 

#####################
## Set-up 
#####################

rm(list = ls())

mycol1 <- rgb(255, 255, 255, max = 255, alpha = 100, names = "white")
mycol2 <- rgb(75, 75, 75, max = 255, alpha = 100, names = "lightgray") 
mycol3 <- rgb(0, 0, 0, max = 200, alpha = 125, names = "darkgray")

#####################
# Model 1
#####################
fd1 <- function(n, a) {
  MIL <- rnorm(n, 0, 1) # military
  e_u <- rnorm(n, 0, 1) # error in gods (unobserved)
  e_f <- rnorm(n, 0, 1) # error in military  
  e_s <- rnorm(n, 0, 1) # error in social complexity (unobserved)
  e_w <- rnorm(n, 0, 1) # error in writing
  e_m <- rnorm(n, 0, 1) # error in missingness
  e_q <- rnorm(n, 0, 1) # error in gods (observed)
  e_so <- rnorm(n, 0, 1) # error in social complexity (observed)
  MGU <- MIL * a + e_u # moralistic gods (unobserved)
  SCU <- MGU * a + MIL * a + e_s # social complexity (unobserved)
  WRI <- SCU * a + e_w # writing
  MIS <- WRI * a + e_m # missingness
  SCO <- SCU * a + MIS * a + e_so # social complexity (observed)
  MGO <- MGU * a + MIS * a + e_q # moralistic gods (observed)
  df <- data.frame(MGU, MIL, WRI, MIS, MGO, SCO) # data frame
  open <- coef(lm(SCO ~ MGO, data = df))[2] # social complexity ~ moralistic gods
  controlled <- coef(lm(SCO ~ MGO + MIL, data = df))[2] # ... + military
  return(c(open, controlled))
}

## Start here. Alter beta to examine shifts in estimates and magnitude of difference between models.

beta <- 0.5

trisim1 <- data.frame(t(replicate(1000, fd1(100, beta))))
names(trisim1) <- c("open", "controlled")

densop1 <- density(trisim1$open) # distribution of effect of moralistic gods on social complexity with no other variables
densco1 <- density(trisim1$controlled) # distributon effect when MIL is held constant

par(mfrow = c(2, 1), mar = c(2, 1, 1, 1)) # bottom, left, top, right

plot(dagitty('dag {
bb="0,0,1,1"
             "Moral Gods (Obs.)" [exposure,pos="0.267,0.135"]
             "Moral Gods (Unobs.)" [latent,pos="0.226,0.070"]
             "Social Complexity (Obs.)" [outcome,pos="0.180,0.136"]
             "Social Complexity (Unobs.)" [latent,pos="0.123,0.070"]
             Military [adjusted,pos="0.176,0.025"]
             Missingness [latent,pos="0.228,0.222"]
             Recording [latent,pos="0.123,0.222"]
             "Moral Gods (Unobs.)" -> "Moral Gods (Obs.)"
             "Moral Gods (Unobs.)" -> "Social Complexity (Unobs.)"
             "Social Complexity (Unobs.)" -> "Social Complexity (Obs.)"
             "Social Complexity (Unobs.)" -> Recording
             Military -> "Moral Gods (Unobs.)"
             Military -> "Social Complexity (Unobs.)"
             Missingness -> "Moral Gods (Obs.)"
             Missingness -> "Social Complexity (Obs.)"
             Recording -> Missingness
             }
             '))

plot(NA, xlab = NA, ylab = "", 
     xlim = c(0, 1), 
     ylim = c(0, 5), 
     cex.lab = 1.3, yaxt = 'n')
polygon(densop1, col = mycol1) # SCO ~ MGO
polygon(densco1, col = mycol3) # SCO ~ MGO + MIL
abline(v = .5, lty = 2)
legend("topright", legend = c("~ Moralistic Gods", 
                             "~ Moralistic Gods + Military"), 
       fill = c(mycol1, mycol3), cex = .85, horiz = F, bty = T, inset = c(0.05, .15))
