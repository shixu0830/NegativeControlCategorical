rm(list=ls())
mis.list = list(c(F,F,F,F),c(T,F,F,F),c(F,T,F,F),c(F,F,F,T),c(F,F,T,T))
n=2000;n.rep=400
my.filepath=""
mean.all = se.all =var.est.all= median.all = mse.ATE=meanbias.ATE=medianbias.ATE=percent.ATE=coverage.ATE=typeIerror=NULL
params.all = NULL
rslt.all = list()

for(i in 1:(length(mis.list))){
  mis.i = mis.list[[i]]
  MW.wrong =mis.i[1]  ## E[W|AZX] wrong and deltaW in IPW misspecify
  MR.wrong =mis.i[2]  ## R wrong
  MY.wrong =mis.i[3]  ## E[Y|Z=0,AX] wrong
  MZ.wrong =mis.i[4]  ## f(Z|A,X) wrong
  
  rslt.tt=rslt.mr=NULL
  for(myseed in 1:10){
    tmp=try(load(paste0(my.filepath,
                        "rslt_W_",MW.wrong,"_Z_",MZ.wrong,"_R_",MR.wrong,"_Y_",MY.wrong,"_n",n,"_nrep",n.rep,"_seed",myseed,".RData")),
            silent=TRUE)
    if(!inherits(tmp, "try-error")){
      rslt.tt = rbind(rslt.tt,rslt)
    }
  }
  rslt = rslt.tt
  ind.mr1 = grep("MR",colnames(rslt))[c(1,3)]
  ind.mr2 = grep("MR",colnames(rslt))[c(2,4)]
  if(i==3){
    rslt=rslt[,-ind.mr1]
  }else{
    rslt=rslt[,-ind.mr2]
  }
  
  ### stablize the result by removing outliers
  ind=NULL
  if(i==2){
    t.quantile=0.05
    ind = c(ind,unlist(apply(cbind(rslt[,grep("MIAO.var",colnames(rslt))[1]]),2,FUN=function(x){
      t=quantile(x,probs=c(t.quantile,1-t.quantile),na.rm=T)
      ind=which(x<t[1]|x>t[2])
      return(ind)
    })))
  }
  if(length(ind)>0){rslt=rslt[-ind,]}
  
  #### compute bias variance coverage 
  rslt.ATE = rslt[,grep(".ATE",colnames(rslt))]
  true.ATE = mean(rslt[,"ATE.true"],na.rm=T)
  mean.all = c(mean.all,apply(rslt.ATE,2,mean,na.rm=T))
  median.all = c(median.all,apply(rslt.ATE,2,median,na.rm=T))
  se.all = c( se.all, apply(rslt.ATE,2,var,na.rm=T) )
  var.est.all = c( var.est.all, (
    apply(rslt[,grep("var",colnames(rslt))],2,mean,na.rm=T)) )
  
  percent.ATE = c(percent.ATE,
                  apply(rslt.ATE,2,FUN=function(x){mean(x-true.ATE,na.rm=T)/true.ATE*100}))
  mse.ATE = c(mse.ATE,apply(rslt.ATE,2,FUN=function(x){mean((x-true.ATE)^2,na.rm=T)}))
  meanbias.ATE = c(meanbias.ATE,
                   apply(rslt.ATE,2,FUN=function(x){mean(x-true.ATE,na.rm=T)}))
  medianbias.ATE = c(medianbias.ATE,
                     apply(rslt.ATE,2,FUN=function(x){median(x-true.ATE,na.rm=T)}))
  t1=rslt.ATE+qnorm(0.975)*sqrt(rslt[,grep("var",colnames(rslt))])
  t2=rslt.ATE-qnorm(0.975)*sqrt(rslt[,grep("var",colnames(rslt))])
  t=mean(rslt[,"ATE.true"])>t2&mean(rslt[,"ATE.true"])<t1
  coverage.ATE=c(coverage.ATE,
                 apply(t,2,mean,na.rm=T))
  typeIerror = c(typeIerror,
                 apply(2*( 1-pnorm(  abs(rslt.ATE) / sqrt(rslt[,grep("var",colnames(rslt))])  ) ),2,FUN=function(x){mean(x<0.05)})
  )
}


### print table
round.text = "%.2f"
summary.rslt = cbind(meanbias=sprintf(round.text,meanbias.ATE*1000),
                     var=sprintf(round.text,se.all*1000),
                     est.var=sprintf(round.text,var.est.all*1000),
                     percent.bias=sprintf(round.text,percent.ATE),
                     mse=sprintf(round.text,mse.ATE*1000),
                     coverage=sprintf(round.text,coverage.ATE),
                     typeIerror=sprintf(round.text,typeIerror)
)
row.names(summary.rslt)=names(meanbias.ATE)
dd=5
summary.rslt=rbind(summary.rslt[1:dd,],rep("",ncol(summary.rslt)),
                   summary.rslt[(dd+1):(dd*2),],rep("",ncol(summary.rslt)),
                   summary.rslt[(dd*2+1):(dd*3),],rep("",ncol(summary.rslt)),
                   summary.rslt[(dd*3+1):(dd*4),],rep("",ncol(summary.rslt)),
                   summary.rslt[(dd*4+1):(dd*5),],rep("",ncol(summary.rslt))
)
if(dd==5){
  summary.rslt[dd+1+1,] = summary.rslt[1,]
  summary.rslt[2*(dd+1)+2,] = summary.rslt[2,]
  summary.rslt[3*(dd+1)+3,] = summary.rslt[3,]
  summary.rslt[2*(dd+1)+4,] = summary.rslt[3*(dd+1)+4,] = summary.rslt[4,]
}
summary.rslt = summary.rslt[-c(7,14,16,21,22),]
View(summary.rslt)
write.csv(summary.rslt,"simulation_results.csv")

