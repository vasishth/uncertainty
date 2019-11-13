## ----setup,include=FALSE,cache=FALSE,echo=FALSE--------------------------
library(MASS)
library(knitr)
library(xtable)
library(papaja)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
#library(sjPlot)
library(rstan)
options(mc.cores = parallel::detectCores())
library(brms)
library(gridExtra)
library(bayesplot)
library(ggridges)
library(lme4)
library(reshape2)

theme_set(theme_apa())

# set global chunk options, put figures into folder
options(warn=-1, replace.assign=TRUE)
opts_chunk$set(fig.path='figures/figure-', fig.align='center', fig.show='hold')
options(replace.assign=TRUE,width=75)
opts_chunk$set(dev='postscript')
opts_chunk$set(echo = TRUE)


## ----loadfunctions, include=TRUE, echo=FALSE, warning=FALSE, error=TRUE, message=TRUE----
source("../R/createStanDat.R")
source("../R/createStanDatAcc.R")
source("../R/magnifytext.R")
source("../R/multiplot.R")
source("../R/plotresults.R")
source("../R/plotpredictions.R")
source("../R/stan_results.R")


## ----powerplots,echo=FALSE,eval=FALSE------------------------------------
## nsim<-10000
## pow30<-pow40<-pow50<-pow60<-rep(NA,nsim)
## stddev <- 150
## d_mean <-30
## d_sd <- 10
## 
## for(i in 1:nsim){
##   d<-rnorm(1,mean=d_mean,sd=d_sd)
##   pow30[i]<-power.t.test(delta=d,n=30,sd=stddev,type="one.sample",strict=TRUE)$power
##   pow40[i]<-power.t.test(delta=d,n=40,sd=stddev,type="one.sample",strict=TRUE)$power
##   pow50[i]<-power.t.test(delta=d,n=50,sd=stddev,type="one.sample",strict=TRUE)$power
##   pow60[i]<-power.t.test(delta=d,n=60,sd=stddev,type="one.sample",strict=TRUE)$power
## }
## 
## powerdistrn<-data.frame(nsubj=factor(rep(c(30,40,50,60),each=nsim)),power=c(pow30,pow40,pow50,pow60))
## 
## sd150<-ggplot(powerdistrn,
##        aes(x = power, y = nsubj, fill = nsubj, height = ..density..)) +
##   geom_density_ridges(scale = 4, stat = "density") +
##   scale_y_discrete(expand = c(0.01, 0)) +
##   #scale_x_continuous(expand = c(0.05, 0)) +
##   xlim(0.05,1)+
##   scale_fill_brewer(palette = "PuBuGn") +
##   theme_ridges() + theme(legend.position = "none")+
##   xlab("prospective statistical power")+
##   ggtitle("Power estimates \n (sd=150)")+
##   ylab("")+
##   magnifytext(sze=16)
## 
## ## sd 200
## pow30<-pow40<-pow50<-pow60<-rep(NA,nsim)
## stddev <- 200
## 
## for(i in 1:nsim){
##   d<-rnorm(1,mean=d_mean,sd=d_sd)
##   pow30[i]<-power.t.test(delta=d,n=30,sd=stddev,type="one.sample",strict=TRUE)$power
##   pow40[i]<-power.t.test(delta=d,n=40,sd=stddev,type="one.sample",strict=TRUE)$power
##   pow50[i]<-power.t.test(delta=d,n=50,sd=stddev,type="one.sample",strict=TRUE)$power
##   pow60[i]<-power.t.test(delta=d,n=60,sd=stddev,type="one.sample",strict=TRUE)$power
## }
## 
## powerdistrn<-data.frame(nsubj=factor(rep(c(30,40,50,60),each=nsim)),power=c(pow30,pow40,pow50,pow60))
## 
## sd200<-ggplot(powerdistrn,
##        aes(x = power, y = nsubj, fill = nsubj, height = ..density..)) +
##   geom_density_ridges(scale = 4, stat = "density") +
##   scale_y_discrete(expand = c(0.01, 0)) +
##   #scale_x_continuous(expand = c(0.05, 0)) +
##   xlim(0.05,1)+
##   scale_fill_brewer(palette = "PuBuGn") +
##   theme_ridges() + theme(legend.position = "none")+
##   xlab("prospective statistical power")+
##   ggtitle("Power estimates \n (sd=200)")+
##   ylab("number of subjects")+
##   magnifytext(sze=16)
## 
## ## sd 250
## pow30<-pow40<-pow50<-pow60<-rep(NA,nsim)
## stddev <- 250
## 
## for(i in 1:nsim){
##   d<-rnorm(1,mean=d_mean,sd=d_sd)
##   pow30[i]<-power.t.test(delta=d,n=30,sd=stddev,type="one.sample",strict=TRUE)$power
##   pow40[i]<-power.t.test(delta=d,n=40,sd=stddev,type="one.sample",strict=TRUE)$power
##   pow50[i]<-power.t.test(delta=d,n=50,sd=stddev,type="one.sample",strict=TRUE)$power
##   pow60[i]<-power.t.test(delta=d,n=60,sd=stddev,type="one.sample",strict=TRUE)$power
## }
## 
## powerdistrn<-data.frame(nsubj=factor(rep(c(30,40,50,60),each=nsim)),power=c(pow30,pow40,pow50,pow60))
## 
## sd250<-ggplot(powerdistrn,
##        aes(x = power, y = nsubj, fill = nsubj, height = ..density..)) +
##   geom_density_ridges(scale = 4, stat = "density") +
##   scale_y_discrete(expand = c(0.01, 0)) +
##   #scale_x_continuous(expand = c(0.05, 0)) +
##   xlim(0.05,1)+
##   scale_fill_brewer(palette = "PuBuGn") +
##   theme_ridges() + theme(legend.position = "none")+
##   xlab("")+
##   ggtitle("Power estimates \n (sd=250)")+
##   ylab("")+
##   magnifytext(sze=16)
## 
## # sd 300
## pow30<-pow40<-pow50<-pow60<-rep(NA,nsim)
## stddev <- 300
## 
## for(i in 1:nsim){
##   d<-rnorm(1,mean=d_mean,sd=d_sd)
##   pow30[i]<-power.t.test(delta=d,n=30,sd=stddev,type="one.sample",strict=TRUE)$power
##   pow40[i]<-power.t.test(delta=d,n=40,sd=stddev,type="one.sample",strict=TRUE)$power
##   pow50[i]<-power.t.test(delta=d,n=50,sd=stddev,type="one.sample",strict=TRUE)$power
##   pow60[i]<-power.t.test(delta=d,n=60,sd=stddev,type="one.sample",strict=TRUE)$power
## }
## 
## powerdistrn<-data.frame(nsubj=factor(rep(c(30,40,50,60),each=nsim)),power=c(pow30,pow40,pow50,pow60))
## 
## sd300<-ggplot(powerdistrn,
##        aes(x = power, y = nsubj, fill = nsubj, height = ..density..)) +
##   geom_density_ridges(scale = 4, stat = "density") +
##   scale_y_discrete(expand = c(0.01, 0)) +
##   #scale_x_continuous(expand = c(0.05, 0)) +
##   xlim(0.05,1)+
##   scale_fill_brewer(palette = "PuBuGn") +
##   theme_ridges() + theme(legend.position = "none")+
##   xlab("")+
##   ggtitle("Power estimates \n (sd=300, \n effect ~ Normal(30,10))")+
##   ylab("number of subjects")+
##   magnifytext(sze=16)
## 
## multiplot(sd300,sd250,
##           sd200,sd150,cols=2)


## ----loadbf,echo=FALSE---------------------------------------------------
bayesfactors<-read.table("../data/bayesfactors.txt")


## ----agrmtattrn,eval=FALSE,echo=FALSE,warning=FALSE,message=FALSE--------
## ## Dillon et al 2013 Expt 1
## DillonE1<-read.table("../data/DillonE1.txt",header=T)
## ## Lago et al 2015 data (all expts):
## Lago<-read.csv("../data/Lago.csv",header=T)
## ##Wagers et al 2009 data (all expts):
## load("../data/Wagers.Rdata")
## load("../data/Tucker.RData")
## 
## ## Dillon E1:
## DillonE1$cond<-factor(DillonE1$cond)
## DillonE1Mism<-subset(DillonE1,fixationtype=="tt" & region==5 & cond%in%c(3,4) & value!="NA")
## DillonE1Mism$cond<-factor(DillonE1Mism$cond)
## DillonE1Mism$int<-ifelse(DillonE1Mism$cond==3,"low","high")
## DillonE1Mism$x<-ifelse(DillonE1Mism$cond==3,-1,1)
## dillonE1<-DillonE1Mism[,c(1,3,4,14,15)]
## dillonE1$expt<-factor("dillonE1")
## colnames(dillonE1)[3]<-"rt"
## 
## nsubj_dillonE1<-length(unique(dillonE1$subj))
## 
## ##Lago:
## dat<-Lago
## ## critical region: not used because published paper found
## ## significant effects in postcrit region only
## e1<-subset(dat,Experiment=="Experiment1" & Region=="06v1")
## e2<-subset(dat,Experiment=="Experiment2" & Region=="06aux")
## e3a<-subset(dat,Experiment=="Experiment3A" & Region=="06aux")
## e3b<-subset(dat,Experiment=="Experiment3B" & Region=="aux")
## 
## nsubj_lagoe1<-length(unique(e1$Subject))
## nsubj_lagoe2<-length(unique(e2$Subject))
## nsubj_lagoe3a<-length(unique(e3a$Subject))
## nsubj_lagoe3b<-length(unique(e3b$Subject))
## 
## 
## ## postcritical region:
## poste1<-subset(dat,Experiment=="Experiment1" & Region=="07prep")
## poste2<-subset(dat,Experiment=="Experiment2" & Region=="07adv")
## poste3a<-subset(dat,Experiment=="Experiment3A" & Region=="07a")
## poste3b<-subset(dat,Experiment=="Experiment3B" & Region=="a")
## 
## ##e1: a,b
## #-(a) Ungram , singular attractor (interference condition)
## #La *nota* que la chica escribieron en la clase alegró a su amiga
## #The note that the girl wrotepl during class cheered her friend up
## #-(b) Ungram , plural attractor (baseline condition)
## #Las *notas* que la chica escribieron en la clase alegraron a su amiga
## #The notes that the girl wrotepl during class cheered her friend up
## poste1<-subset(poste1,Condition%in%c("a","b"))
## poste1$Condition<-factor(poste1$Condition)
## poste1$x<-ifelse(poste1$Condition=="a",-1,1)
## poste1$int<-ifelse(poste1$Condition=="a","low","high")
## poste1<-poste1[,c(1,3,8,15,14)]
## poste1$expt<-factor("lagoE1")
## lagoE1<-poste1
## colnames(lagoE1)<-c("subj","item","rt","int","x","expt")
## 
## ## e2: c,d
## poste2<-subset(poste2,Condition%in%c("c","d"))
## poste2$Condition<-factor(poste2$Condition)
## poste2$x<-ifelse(poste2$Condition=="c",-1,1)
## poste2$int<-ifelse(poste2$Condition=="c","low","high")
## #head(poste2)
## poste2<-poste2[,c(1,3,8,15,14)]
## poste2$expt<-factor("lagoE2")
## lagoE2<-poste2
## colnames(lagoE2)<-c("subj","item","rt","int","x","expt")
## 
## ## e3a: e,f
## poste3a<-subset(poste3a,Condition%in%c("e","f"))
## poste3a$Condition<-factor(poste3a$Condition)
## #-(e) Ungram, singular attractor (interference condition)
## #La *nota* que la chica van a escribir en la clase alegrará a su amiga
## #The note that the girl are going to write during class will cheer her friend up
## #-(f) Ungram, plural attractor (baseline condition)
## #Las *notas* que la chica van a escribir en la clase alegrarán a su amiga
## #The notes that the girl are going to write during class will cheer her friend up
## #boxplot(RT~Condition,poste3a)
## poste3a$x<-ifelse(poste3a$Condition=="e",-1,1)
## poste3a$int<-ifelse(poste3a$Condition=="e","low","high")
## poste3a<-poste3a[,c(1,3,8,15,14)]
## poste3a$expt<-factor("lagoE3a")
## lagoE3a<-poste3a
## colnames(lagoE3a)<-c("subj","item","rt","int","x","expt")
## 
## ## e3b: e,f
## poste3b<-subset(poste3b,Condition%in%c("e","f"))
## poste3b$Condition<-factor(poste3b$Condition)
## #-(e) Ungram, singular attractor (baseline condition)
## #The player that the coach were always praising very enthusiastically decided to     leave the team
## #-(f) Ungram, plural attractor (interference condition)
## #The players that the coach were always praising very enthusiastically decided to     leave the team
## poste3b$x<-ifelse(poste3b$Condition=="e",-1,1)
## poste3b$int<-ifelse(poste3b$Condition=="e","low","high")
## poste3b<-poste3b[,c(1,3,8,15,14)]
## poste3b$expt<-factor("lagoE3b")
## lagoE3b<-poste3b
## colnames(lagoE3b)<-c("subj","item","rt","int","x","expt")
## 
## ## Wagers:
## E2postcrit<-subset(Experiment2,Region==7)
## nsubj_wagerse2<-length(unique(E2postcrit$Subj))
## #E2$intr.au<-ifelse(E2$rchead=="pl" & E2$gramm=="ungram",1/2,
## #                   ifelse(E2$rchead=="sg" & E2$gramm=="ungram",-1/2,
## #                          0))
## ## d (sing),h (plu)
## #unique(subset(E2postcrit,gramm=="ungram")$Condition)
## E2postcrit<-subset(E2postcrit,Condition%in%c("d","h"))
## E2postcrit$Condition<-factor(E2postcrit$Condition)
## E2postcrit$x<-ifelse(E2postcrit$Condition=="d",-1,1)
## E2postcrit$int<-ifelse(E2postcrit$Condition=="d","low","high")
## #colnames(E2postcrit)
## E2postcrit<-E2postcrit[,c(4,3,8,13,12)]
## E2postcrit$expt<-factor("wagersE2")
## wagersE2<-E2postcrit
## colnames(wagersE2)<-c("subj","item","rt","int","x","expt")
## 
## ## E3
## E3postcrit<-subset(Experiment3,Region==7)
## nsubj_wagerse3<-length(unique(E3postcrit$Subj))
## #E3crit$intr.au.pl<-ifelse(E3crit$gramm=="ungram" & E3crit$rcsubj=="sg" &
## #                            E3crit$rchead=="pl",1/2,
## #                         ifelse(E3crit$gramm=="ungram" & E3crit$rcsubj=="sg" &
## #                                   E3crit$rchead=="sg",-1/2,0))
## 
## #E3crit$intr.au.sg<-ifelse(E3crit$gramm=="ungram" & E3crit$rcsubj=="pl" &
## #                            E3crit$rchead=="sg",1/2,
## #                          ifelse(E3crit$gramm=="ungram" & E3crit$rcsubj=="pl" &
## #                                   E3crit$rchead=="pl",-1/2,0))
## 
## E3postcrit_pl<-subset(E3postcrit,gramm=="ungram" & rcsubj=="sg")
## E3postcrit_pl$Condition<-factor(E3postcrit_pl$Condition)
## E3postcrit_sg<-subset(E3postcrit,gramm=="ungram" & rcsubj=="pl")
## E3postcrit_sg$Condition<-factor(E3postcrit_sg$Condition)
## 
## #unique(E3postcrit_pl$Condition) ## b,f
## #unique(E3postcrit_sg$Condition) ## c,g
## 
## #head(subset(E3postcrit_sg,rchead=="sg"))
## #head(E3postcrit_pl)
## 
## ## plural:
## E3postcrit_pl$x<-ifelse(E3postcrit_pl$Condition=="b",-1,1)
## E3postcrit_pl$int<-ifelse(E3postcrit_pl$Condition=="b","low","high")
## E3postcrit_pl<-E3postcrit_pl[,c(4,3,8,15,14)]
## E3postcrit_pl$expt<-factor("wagersE3pl")
## colnames(E3postcrit_pl)<-c("subj","item","rt","int","x","expt")
## wagersE3pl<-E3postcrit_pl
## 
## ## singular:
## E3postcrit_sg$x<-ifelse(E3postcrit_sg$Condition=="c",-1,1)
## E3postcrit_sg$int<-ifelse(E3postcrit_sg$Condition=="c","low","high")
## E3postcrit_sg<-E3postcrit_sg[,c(4,3,8,15,14)]
## E3postcrit_sg$expt<-factor("wagersE3sg")
## colnames(E3postcrit_sg)<-c("subj","item","rt","int","x","expt")
## wagersE3sg<-E3postcrit_sg
## 
## ## E4
## E4postcrit<-subset(Experiment4,Region==8) ##
## nsubj_wagerse4<-length(unique(E4postcrit$Subj))
## #head(subset(Experiment4,Condition=="c"),n=10)
## #postcritical region
## #E4postcrit$intr.au<-ifelse(E4postcrit$gramm=="ungram" & E4postcrit$match=="match",-1/2,
## #                           ifelse(E4postcrit$gramm=="ungram" & E4postcrit$match=="mismatch",1/2,0))
## E4postcrit<-subset(E4postcrit,gramm=="ungram")
## E4postcrit$Condition<-factor(E4postcrit$Condition)
## E4postcrit$x<-ifelse(E4postcrit$Condition=="c",-1,1)
## E4postcrit$int<-ifelse(E4postcrit$Condition=="c","low","high")
## E4postcrit<-E4postcrit[,c(4,3,8,13,12)]
## E4postcrit$expt<-factor("wagersE4")
## colnames(E4postcrit)<-c("subj","item","rt","int","x","expt")
## wagersE4<-E4postcrit
## 
## # E5
## E5postcrit<-subset(Experiment5,Region==8) ##postcritical region
## nsubj_wagerse5<-length(unique(E5postcrit$Subj))
## E5postcrit<-subset(E5postcrit,gramm=="ungram")
## E5postcrit$Condition<-factor(E5postcrit$Condition)
## ## c,d
## E5postcrit$x<-ifelse(E5postcrit$Condition=="c",-1,1)
## E5postcrit$int<-ifelse(E5postcrit$Condition=="c","low","high")
## E5postcrit<-E5postcrit[,c(4,3,8,13,12)]
## colnames(E5postcrit)<-c("subj","item","rt","int","x")
## E5postcrit$expt<-factor("wagersE5")
## wagersE5<-E5postcrit
## 
## #head(wagersE5)
## 
## dat<-rbind(dillonE1,wagersE2,
##            lagoE1,lagoE2,
##            lagoE3a,lagoE3b,
##            wagersE2,
##            wagersE3pl,wagersE3sg,
##            wagersE4,wagersE5)
## dat$subj<-factor(paste(dat$expt,dat$subj,sep=""))
## dat$item<-factor(paste(dat$expt,dat$item,sep=""))
## #with(dat,tapply(subj,expt,function(x)length(unique(x))))


## ----agrmtattrn2,echo=FALSE,warning=FALSE,message=FALSE,eval=FALSE-------
## ## Dillon E1:
## stanDat<-createStanDat(d=dillonE1,
##                        rt=dillonE1$rt,
##                        form=as.formula("~ 1 + x"))
## #str(stanDat)
## DillonE1 <- stan(file = "StanModels/maxModelTargetMismatch.stan",
##              data = stanDat,
##              iter = 2000,
##              chains = 4)
## ## Int is Interference:
## pars<-c("Int","beta[2]","sigma_u[1]","sigma_u[2]","sigma_w[1]","sigma_w[2]","sigma_e")
## DillonE1_res<-stan_results(DillonE1,params=pars[1])
## 
## #with(dillonE1,tapply(rt,x,mean))
## 
## mDillonE1_lmer<-lmer(log(rt)~x+(1+x|subj)+(1+x|item),dillonE1)
## mDillonE1_lmer_res<-summary(mDillonE1_lmer)$coefficients[2,]
## 
## stanDat<-createStanDat(d=lagoE1,
##                        rt=lagoE1$rt,
##                        form=as.formula("~ 1 + x"))
## 
## LagoE1 <- stan(file = "StanModels/maxModelTargetMismatch.stan",
##                  data = stanDat,
##                  iter = 2000,
##                  chains = 4)
## LagoE1_res<-stan_results(LagoE1,params=pars[1])
## 
## mLagoE1_lmer<-lmer(log(rt)~x+(1+x||subj)+(1+x||item),lagoE1)
## mLagoE1_lmer_res<-summary(mLagoE1_lmer)$coefficients[2,]
## 
## 
## stanDat<-createStanDat(d=lagoE2,
##                        rt=lagoE2$rt,
##                        form=as.formula("~ 1 + x"))
## 
## LagoE2 <- stan(file = "StanModels/maxModelTargetMismatch.stan",
##                data = stanDat,
##                iter = 2000,
##                chains = 4)
## LagoE2_res<-stan_results(LagoE2,params=pars[1])
## 
## mLagoE2_lmer<-lmer(log(rt)~x+(1+x||subj)+(1+x||item),lagoE2)
## mLagoE2_lmer_res<-summary(mLagoE2_lmer)$coefficients[2,]
## 
## 
## stanDat<-createStanDat(d=lagoE3a,
##                           rt=lagoE3a$rt,
##                        form=as.formula("~ 1 + x"))
## 
## LagoE3a <- stan(file = "StanModels/maxModelTargetMismatch.stan",
##                data = stanDat,
##                iter = 2000,
##                chains = 4)
## 
## LagoE3a_res<-stan_results(LagoE3a,params=pars[1])
## 
## mLagoE3a_lmer<-lmer(log(rt)~x+(1|subj)+(1|item),lagoE3a)
## mLagoE3a_lmer_res<-summary(mLagoE3a_lmer)$coefficients[2,]
## 
## stanDat<-createStanDat(d=lagoE3b,
##                           rt=lagoE3b$rt,
##                        form=as.formula("~ 1 + x"))
## 
## LagoE3b <- stan(file = "StanModels/maxModelTargetMismatch.stan",
##                 data = stanDat,
##                 iter = 2000,
##                 chains = 4)
## 
## LagoE3b_res<-stan_results(LagoE3b,params=pars[1])
## 
## mLagoE3b_lmer<-lmer(log(rt)~x+(1+x||subj)+(1|item),lagoE3b)
## mLagoE3b_lmer_res<-summary(mLagoE3b_lmer)$coefficients[2,]
## 
## 
## stanDat<-createStanDat(d=wagersE2,
##                           rt=wagersE2$rt,
##                        form=as.formula("~ 1 + x"))
## 
## WagersE2 <- stan(file = "StanModels/maxModelTargetMismatch.stan",
##                 data = stanDat,
##                 iter = 2000,
##                 chains = 4)
## WagersE2_res<-stan_results(WagersE2,params=pars[1])
## 
## mWagersE2_lmer<-lmer(log(rt)~x+(1+x||subj)+(1+x||item),wagersE2)
## mWagersE2_lmer_res<-summary(mWagersE2_lmer)$coefficients[2,]
## 
## stanDat<-createStanDat(d=wagersE3pl,
##                           rt=wagersE3pl$rt,
##                        form=as.formula("~ 1 + x"))
## 
## WagersE3pl <- stan(file = "StanModels/maxModelTargetMismatch.stan",
##                  data = stanDat,
##                  iter = 2000,
##                  chains = 4)
## WagersE3pl_res<-stan_results(WagersE3pl,params=pars[1])
## 
## mWagersE3pl_lmer<-lmer(log(rt)~x+(1+x||subj)+(1+x||item),wagersE3pl)
## mWagersE3pl_lmer_res<-summary(mWagersE3pl_lmer)$coefficients[2,]
## 
## stanDat<-createStanDat(d=wagersE3sg,
##                           rt=wagersE3sg$rt,
##                           form=as.formula("~ 1 + x"))
## 
## WagersE3sg <- stan(file = "StanModels/maxModelTargetMismatch.stan",
##                    data = stanDat,
##                    iter = 2000,
##                    chains = 4)
## WagersE3sg_res<-stan_results(WagersE3sg,params=pars[1])
## 
## mWagersE3sg_lmer<-lmer(log(rt)~x+(1+x||subj)+(1+x||item),wagersE3sg)
## mWagersE3sg_lmer_res<-summary(mWagersE3sg_lmer)$coefficients[2,]
## 
## 
## stanDat<-createStanDat(d=wagersE4,
##                           rt=wagersE4$rt,
##                           form=as.formula("~ 1 + x"))
## 
## WagersE4 <- stan(file = "StanModels/maxModelTargetMismatch.stan",
##                    data = stanDat,
##                    iter = 2000,
##                    chains = 4)
## WagersE4_res<-stan_results(WagersE4,params=pars[1])
## 
## mWagersE4_lmer<-lmer(log(rt)~x+(1|subj),wagersE4)
## mWagersE4_lmer_res<-summary(mWagersE4_lmer)$coefficients[2,]
## 
## stanDat<-createStanDat(d=wagersE5,
##                        rt=wagersE5$rt,
##                        form=as.formula("~ 1 + x"))
## 
## WagersE5 <- stan(file = "StanModels/maxModelTargetMismatch.stan",
##                  data = stanDat,
##                  iter = 2000,
##                  chains = 4)
## WagersE5_res<-stan_results(WagersE5,params=pars[1])
## 
## mWagersE5_lmer<-lmer(log(rt)~x+(1+x|subj)+(1+x|item),wagersE5)
## mWagersE5_lmer_res<-summary(mWagersE5_lmer)$coefficients[2,]
## 
## ## posterior predicted data from Jaeger et al replication: See MethodsX paper:
## load("../model/au_predicted_meansD13rep.Rda")
## len_model<-length(au_predicted_means_rep)
## 
## agrmt_data<-data.frame(expt=factor(rep(c(1:10),
##                                          each=4000)),
##                        posterior=c(DillonE1_res[2][[1]][[1]],
##                        LagoE1_res[2][[1]][[1]],
##                        LagoE2_res[2][[1]][[1]],
##                        LagoE3a_res[2][[1]][[1]],
##                        LagoE3b_res[2][[1]][[1]],
##                        WagersE2_res[2][[1]][[1]],
##                        WagersE3pl_res[2][[1]][[1]],
##                        WagersE3sg_res[2][[1]][[1]],
##                        WagersE4_res[2][[1]][[1]],
##                        WagersE5_res[2][[1]][[1]]))
## 
## means<-sort(with(agrmt_data,tapply(posterior,expt,mean)))
## ordered_studies<-as.numeric(names(means))
## 
## model_pred<-data.frame(expt=rep("model",len_model),posterior=au_predicted_means_rep)
## head(model_pred)
## 
## data_model<-rbind(agrmt_data,model_pred)
## 
## lvls<-levels(data_model$expt)
## data_model$expt<-factor(data_model$expt,levels=lvls[c(11,ordered_studies
## )])
## 
## save(data_model,file="../data/data_model.Rda")


## ----tvals,echo=FALSE,eval=FALSE-----------------------------------------
## tvals<-c(mDillonE1_lmer_res[3],
## mLagoE1_lmer_res[3],
## mLagoE2_lmer_res[3],
## mLagoE3a_lmer_res[3],
## mLagoE3b_lmer_res[3],
## mWagersE2_lmer_res[3],
## mWagersE3pl_lmer_res[3],
## mWagersE3sg_lmer_res[3],
## mWagersE4_lmer_res[3],
## mWagersE5_lmer_res[3])
## tvals<-data.frame(t(tvals))
## colnames(tvals)<-1:10
## library(xtable)
## xtable(tvals)


## ----bayesfactors,echo=FALSE,eval=FALSE----------------------------------
## ## just compute BF for the two  strongest effects seen: Dillon et al E1, and Wagers et al E4.
## library(brms)
## 
## priorsUNINF <- c(set_prior('normal(0, 10)', class='Intercept'),
##             set_prior('normal(0, 1)', class='b'),
##             set_prior('normal(0, 1)', class='sd'),
##             set_prior('normal(0, 1)', class='sigma'),
##             set_prior('lkj(2)', class='cor'))
## 
## priors <- c(set_prior('normal(0, 10)', class='Intercept'),
##             set_prior('normal(0, 0.1)', class='b'),
##             set_prior('normal(0, 1)', class='sd'),
##             set_prior('normal(0, 1)', class='sigma'),
##             set_prior('lkj(2)', class='cor'))
## 
## priors0 <- c(set_prior('normal(0, 10)', class='Intercept'),
##             set_prior('normal(0, 1)', class='sd'),
##             set_prior('normal(0, 1)', class='sigma'),
##             set_prior('lkj(2)', class='cor'))
## 
## m1dillonE1UNINF<-brm(log(rt)~x + (1+x|subj)+(1+x|item),dillonE1,
## chains=3,
## iter=20000,
## warmup=2000,
## save_all_pars=TRUE,
## cores = 4,
## prior=priorsUNINF)
## 
## 
## m1dillonE1<-brm(log(rt)~x + (1+x|subj)+(1+x|item),dillonE1,
## chains=3,
## iter=20000,
## warmup=2000,
## save_all_pars=TRUE,
## cores=4,
## prior=priors)
## 
## m0dillonE1<-brm(log(rt)~ 1 + (1+x|subj)+(1+x|item),dillonE1,
## chains=3,
## iter=20000,
## warmup=2000,
## save_all_pars=TRUE,
## cores=4,
## prior=priors0)
## 
## ## first column for uninf, and second for inf
## bayesfactors<-matrix(rep(NA,20),ncol=2)
## 
## summary(m1dillonE1)
## 
## bayesfactors[1,1]<-bayes_factor(m1dillonE1UNINF,m0dillonE1)$bf
## bayesfactors[1,2]<-bayes_factor(m1dillonE1,m0dillonE1)$bf
## 
## ## wagers expt 3
## head(wagersE4)
## m1wagersE4UNINF<-brm(log(rt)~x + (1+x|subj)+(1+x|item),
##                      wagersE4,
## chains=3,
## iter=20000,
## warmup=2000,
## save_all_pars=TRUE,
## cores = 4,
## prior=priorsUNINF)
## 
## 
## m1wagersE4<-brm(log(rt)~x + (1+x|subj)+(1+x|item),wagersE4,
## chains=3,
## iter=20000,
## warmup=2000,
## save_all_pars=TRUE,
## cores=4,
## prior=priors)
## 
## summary(m1wagersE4)
## 
## m0wagersE4<-brm(log(rt)~ 1 + (1+x|subj)+(1+x|item),wagersE4,
## chains=3,
## iter=20000,
## warmup=2000,
## save_all_pars=TRUE,
## cores=4,
## prior=priors0)
## 
## bayesfactors[9,1]<-bayes_factor(m1wagersE4UNINF,m0wagersE4)$bf
## bayesfactors[9,2]<-bayes_factor(m1wagersE4,m0wagersE4)$bf
## write.table(bayesfactors,file="../data/bayesfactors.txt")


## ----jmvv2019rep,echo=FALSE,eval=FALSE-----------------------------------
## #save(crit,file="../data/JMVV2019replication.Rda")
## load("../data/JMVV2019replication.Rda")
## ## from JMVV_R1.Rnw:
## crit_tft_m2 <- createStanDat(d=subset(crit, TFT>0),
##                              form=as.formula("~1+Dep+Gram+DepxGram+Int_gram_refl+Int_gram_agr+Int_ungram_refl+Int_ungram_agr"),
##                              dv=subset(crit, TFT>0)$TFT)
## 
## # sample from posterior distribution.
## M2_tft <- stan(file = "StanModels/maxModel2.stan",
## data = crit_tft_m2,
## iter = 2000,
## chains = 4,
## control = list(adapt_delta=0.99))
## 
## DillonrepM2_res<-stan_results(M2_tft,params="Int_ungram_agr")
## 
## #str(DillonrepM2_res)
## #with(dillonE1,tapply(rt,x,mean))
## 
## dillonrep_posterior<-data.frame(expt=factor(rep(11,4000)),
##                                 posterior=DillonrepM2_res[2][[1]][[1]])
## 
## data_model_dillonrep<-rbind(data_model,dillonrep_posterior)
## save(data_model_dillonrep,file="../data/data_model_dillonrep.Rda")


## ----jmvv2019replmer,echo=FALSE,eval=FALSE-------------------------------
## mDillonrep_lmer<-lmer(log(TFT)~1+Dep+Gram+DepxGram+Int_gram_refl+Int_gram_agr+Int_ungram_refl+Int_ungram_agr+(1+Dep+Gram+DepxGram+Int_gram_refl+Int_gram_agr+Int_ungram_refl+Int_ungram_agr|subj)+(1+Dep+Gram+DepxGram+Int_gram_refl+Int_gram_agr+Int_ungram_refl+Int_ungram_agr|item),
##                       subset(crit, TFT>0))
## 
## mDillonrep_lmer_res<-summary(mDillonrep_lmer)$coefficients[8,]
## #  Estimate Std. Error    t value
## # -0.043589   0.024239  -1.798271


## ----ridgeplot,echo=FALSE------------------------------------------------
## precomputed above:
load("../data/data_model.Rda")
load("../data/data_model_dillonrep.Rda")

ggplot(data_model, 
       aes(x = posterior, y = expt, 
           fill = expt, height = ..density..)) +
  geom_density_ridges(scale = 4, stat = "density") +
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_fill_brewer(palette = "PuBuGn") +
  theme_ridges() + theme(legend.position = "none")+
  xlab("agreement attraction effect")+
  ylab("expt")+
  geom_vline(xintercept=0)+
  #geom_vline(xintercept=mean(process_A),linetype='dashed')+
  #geom_vline(xintercept=mean(process_B),linetype='dashed')+
  #ggtitle("An illustration of a race process")+
  magnifytext(sze=20)

