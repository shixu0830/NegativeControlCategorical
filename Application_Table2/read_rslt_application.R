###########################SIMULATION CODE###########################
### Xu Shi (shixu@umich.edu)
### Below is the code for generating Table 2
###########################SIMULATION CODE###########################
rm(list=ls())
mis.list = list(c(F,F,F,F),c(T,F,F,F),c(F,T,F,F),c(F,F,F,T),c(T,T,F,F)) #,c(F,F,T,T)
mypath=""
add = paste0(mypath,"rslt")
omit.age=F
rslt.all =NULL
for(i in 1:(length(mis.list))){
  mis.i = mis.list[[i]]
  MW.wrong =mis.i[1]  ## E[W|AZX] wrong and deltaW in IPW misspecify
  MR.wrong =mis.i[2]  ## R wrong
  MY.wrong =mis.i[3]  ## E[Y|Z=0,AX] wrong
  MZ.wrong =mis.i[4]
  
  load(file=paste0(add, MZ.wrong,"_W_",MW.wrong,"_R_",MR.wrong,"_Y_",MY.wrong,".RData"))
  rslt[,1:3]=rslt[,1:3]*1000
  rslt[,4:6]=rslt[,4:6]*1000^2
  
  true.ATE = rslt[1,1]
  true.naive = rslt[1,2]
  true.bias = rslt[1,3]
  
  ATE = rslt[2:nrow(rslt),1]
  naive=rslt[2:nrow(rslt),2]
  bias=rslt[2:nrow(rslt),3]
  sd.ATE=sqrt(rslt[2:nrow(rslt),1+3])
  sd.naive=sqrt(rslt[2:nrow(rslt),2+3])
  sd.bias=sqrt(rslt[2:nrow(rslt),3+3])
  
  bias.ATE = ATE-true.ATE
  bias.naive = naive-true.naive
  bias.bias = bias-true.bias
  
  percent_bias.ATE = (ATE-true.ATE)/true.ATE*100
  percent_bias.naive = (naive-true.naive)/true.naive*100
  percent_bias.bias = (bias-true.bias)/true.bias*100
  
  
  p.ATE = 2*( 1-pnorm(  abs(ATE) / sd.ATE ))
  p.naive = 2*( 1-pnorm(  abs(naive) / sd.naive ))
  p.bias = 2*( 1-pnorm(  abs(bias) / sd.bias ))
  
  L.ATE = ATE-2*sd.ATE
  L.naive = naive-2*sd.naive
  L.bias = bias-2*sd.bias
  
  U.ATE = ATE+2*sd.ATE
  U.naive = naive+2*sd.naive
  U.bias = bias+2*sd.bias
  
  print.rslt = cbind(
    ATE=paste0(sprintf("%.1f",ATE)," (",sprintf("%.1f",L.ATE),", ",sprintf("%.1f",U.ATE),")"), 
    percent_bias.ATE=sprintf("%.1f",percent_bias.ATE),
    p.ATE=sprintf("%.1f",p.ATE),
    
    naive=paste0(sprintf("%.1f",naive)," (",sprintf("%.1f",L.naive),", ",sprintf("%.1f",U.naive),")"), 
    percent_bias.naive=sprintf("%.1f",percent_bias.naive),
    p.naive=sprintf("%.1f",p.naive),
    
    bias=paste0(sprintf("%.1f",bias)," (",sprintf("%.1f",L.bias),", ",sprintf("%.1f",U.bias),")"), 
    percent_bias.bias=sprintf("%.1f",percent_bias.bias),
    p.bias=sprintf("%.1f",p.bias)
    
  )
  
  rslt.all = rbind(rslt.all,rep("",ncol(print.rslt)),print.rslt)
}
write.csv(rslt.all,file=paste0(add,"application_results.csv"))

