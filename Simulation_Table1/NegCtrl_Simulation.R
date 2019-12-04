###########################SIMULATION CODE###########################
### Xu Shi (shixu@umich.edu)
### Below is the code for Section 4 of the manuscript
### the simulation is done using seeds 1-4000
### Delta_1 is called GEST: biaesd if MR.wrong | MZ.wrong
### Delta_2 is called IPW: biaesd if MW.wrong | MZ.wrong
### Delta_3 is called OR: biaesd if MW.wrong | MR.wrong | MY.wrong
### MLE is called MIAO: biased if MW.wrong | MY.wrong
###########################SIMULATION CODE###########################
rm(list=ls())
wrong.i = as.numeric(Sys.getenv("arg1"))
myseed =  as.numeric(Sys.getenv("arg2"))
my.filepath = ""
source(paste0(my.filepath,"functions.R"))
n.rep=400#number of replications per seed
library(maxLik)
library(sandwich)
n=2000
all.wrong = rbind(c(F,F,F,F),c(T,F,F,F),c(F,T,F,F),c(F,F,T,F),c(F,F,T,T))
all.wrong = all.wrong[wrong.i,]
MW.wrong = all.wrong[1] ## E[W|AZX] wrong and deltaW in IPW misspecify
MR.wrong = all.wrong[2] ## R wrong
MZ.wrong = all.wrong[3] ## f(AZX) wrong
MY.wrong = all.wrong[4]  ## E[Y|Z=0,AX] wrong #this one never used
p=10 # 1 intercept, p-2 cov, 1 interaction
alpha = c(-0.1,rep(-0.1,(p-2)/2),rep(-0.1,(p-2)/2),2)/p
beta = c(-1,rep(-1/p,length.out=p-1))
dd = 0.1
tt = -0.2
ma=0#0.2
mb=0.2
mc=0.2
mr0=0#0.2
mr1=0.5 # mr1<1 mr0+mr1>=0
#mr1=0
cc=0.5
(me= ma+mb+mc+dd) ## +0.1 to make sure e-a>b+c>0
mf = 0#1-(me-ma)*(mr0+mr1)-cc ## -0.1 to make sure 1-(e-a)(r0+r1) > f which also gives 1-(b+c)(r0+r1) > f
(mt = me-ma-mb-mc)

rslt = NULL
for(i in 1:n.rep+(myseed-1)*n.rep){
  set.seed(i)
  tmp=try(do.one(n,p,alpha,beta,MW.wrong,MR.wrong,MY.wrong,MZ.wrong),silent=TRUE)
  if(!inherits(tmp, "try-error")){
    rslt=rbind(rslt,unlist(tmp$est))
  }
  if(i%%10==0){
    print(i)
  }
}
colnames(rslt)=names(tmp$est)

save(rslt,
     file=paste0(my.filepath,
                 "rslt_W_",MW.wrong,"_Z_",MZ.wrong,"_R_",MR.wrong,"_Y_",MY.wrong,"_n",n,"_nrep",n.rep,"_seed",myseed,".RData"))
