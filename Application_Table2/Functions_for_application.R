###########################APPLICATION CODE###########################
### Xu Shi (shixu@umich.edu)
### Below is the functions for Section 5 of the manuscript
### Delta_1 is called GEST: biaesd if MR.wrong | MZ.wrong
### Delta_2 is called IPW: biaesd if MW.wrong | MZ.wrong
### Delta_3 is called OR: biaesd if MW.wrong | MR.wrong | MY.wrong
### MLE is called MIAO: biased if MW.wrong | MY.wrong
###########################SIMULATION CODE###########################

U.ipw = function(par,A,Z,Y,W,Age,Sex,MZ.wrong,MW.wrong){
  n=length(A)
  ### A, Z
  A1=rep(1,n);A0=rep(0,n);A2=1-A
  Z1=rep(1,n);Z0=rep(0,n);
  data = data.frame(A0,A1,A2,A,Z0,Z1,Z,Age,Sex)
  
  X.A = model.matrix(A~Age*Sex)
  ### f(A,Z)
  if(MZ.wrong){
    if(omit.sex){
      AX.Z = model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Z",T,T))
      A0X.Z = model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Z",T,T))
      A1X.Z = model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Z",T,T))
      A2X.Z = model.matrix(data=data,object=get.formula(F,"A2",NULL,NULL,"Z",T,T))
      # par.z = lm(Z~A*Sex)$coefficients
    }else{
      AX.Z = model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Z",T,F))
      A0X.Z = model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Z",T,F))
      A1X.Z = model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Z",T,F))
      A2X.Z = model.matrix(data=data,object=get.formula(F,"A2",NULL,NULL,"Z",T,F))
      # par.z = lm(Z~A*(Age+Sex))$coefficients
    }
  }else{
    AX.Z = model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Z",F,F))
    A0X.Z = model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Z",F,F))
    A1X.Z = model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Z",F,F))
    A2X.Z = model.matrix(data=data,object=get.formula(F,"A2",NULL,NULL,"Z",F,F))
    # par.z = lm(Z~A*Age*Sex)$coefficients
  }
  
  par.z = par[1:ncol(AX.Z)] #par.z = lm(Z~A*Age*Sex)$coefficients
  par.a = par[(ncol(AX.Z)+1):(ncol(AX.Z)+ncol(X.A))] #par.a = lm(A~Age*Sex)$coefficients
  PZAX = (AX.Z%*%par.z); PZA1X = (A1X.Z%*%par.z); PZA0X = (A0X.Z%*%par.z); PZA2X = (A2X.Z%*%par.z);
  PAX = X.A%*%par.a; PZX = PZA1X*PAX+PZA0X*(1-PAX)
  PAZX = (PZA1X*Z+(1-PZA1X)*(1-Z)) * (PAX)  / (PZX*Z+(1-PZX)*(1-Z)) #equal to predict(lm(A~Z*Age*Sex))
  fAZX = PAZX*A+(1-PAZX)*(1-A)
  fZAX = PZAX*Z+(1-PZAX)*(1-Z)
  
  ###construct EE for par.w using par.a, par.z 
  Zc = Z-PZX #predict(lm(Z~Age*Sex)) ##Z.center E[Z|X]
  Ac = A-PAX #predict(lm(A~Age*Sex))##A.center E[A|X]
  AZc = A*Z-PZA1X*PAX#predict(lm(A*Z~Age*Sex))
  data = data.frame(A0,A1,A2,A,Z0,Z1,Z,Age,Sex,Ac,Zc,AZc)
  
  ## below are AZX.W delete Age, Sex, Age:Sex
  if(MW.wrong){
    if(omit.sex){
      # par.w_bl = lm(W~  Sex, data=data.frame(Age,Sex,W)[A==0&Z==0,])$coefficients
      X.W = model.matrix(W~Sex)
      AZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z","A*Z","W",T,T))
      AZ1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z1","A*Z1","W",T,T))
      AZ0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z0","A*Z0","W",T,T))
      A1ZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z","A1*Z","W",T,T))
      A0ZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z","A0*Z","W",T,T))
      A1Z1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z1","A1*Z1","W",T,T))
      A1Z0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z0","A1*Z0","W",T,T))
      A0Z1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z1","A0*Z1","W",T,T))
      A0Z0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z0","A0*Z0","W",T,T))
      
      A2Z1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A2","Z1","A2*Z1","W",T,T))
      A2Z0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A2","Z0","A2*Z0","W",T,T))
      AZX.W.center = model.matrix(data=data,object=get.formula(T,"Ac","Zc","AZc","W",T,T)) #each column should have mean zero
    }else{
      # par.w_bl = lm(W~  (Age+Sex), data=data.frame(Age,Sex,W)[A==0&Z==0,])$coefficients
      X.W = model.matrix(W~(Age+Sex))
      AZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z","A*Z","W",T,F))
      AZ1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z1","A*Z1","W",T,F))
      AZ0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z0","A*Z0","W",T,F))
      A1ZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z","A1*Z","W",T,F))
      A0ZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z","A0*Z","W",T,F))
      A1Z1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z1","A1*Z1","W",T,F))
      A1Z0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z0","A1*Z0","W",T,F))
      A0Z1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z1","A0*Z1","W",T,F))
      A0Z0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z0","A0*Z0","W",T,F))
      
      A2Z1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A2","Z1","A2*Z1","W",T,F))
      A2Z0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A2","Z0","A2*Z0","W",T,F))
      AZX.W.center = model.matrix(data=data,object=get.formula(T,"Ac","Zc","AZc","W",T,F)) #each column should have mean zero
    }
  }else{
    # par.w_bl = lm(W~  Age*Sex, data=data.frame(Age,Sex,W)[A==0&Z==0,])$coefficients
    X.W = model.matrix(W~Age*Sex)
    AZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z","A*Z","W",F,F))
    AZ1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z1","A*Z1","W",F,F))
    AZ0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z0","A*Z0","W",F,F))
    A1ZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z","A1*Z","W",F,F))
    A0ZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z","A0*Z","W",F,F))
    A1Z1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z1","A1*Z1","W",F,F))
    A1Z0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z0","A1*Z0","W",F,F))
    A0Z1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z1","A0*Z1","W",F,F))
    A0Z0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z0","A0*Z0","W",F,F))
    
    A2Z1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A2","Z1","A2*Z1","W",F,F))
    A2Z0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A2","Z0","A2*Z0","W",F,F))
    AZX.W.center = model.matrix(data=data,object=get.formula(T,"Ac","Zc","AZc","W",F,F)) #each column should have mean zero
  }
  
  par.w = par[(ncol(AX.Z)+ncol(X.A)+1):(ncol(AX.Z)+ncol(X.A)+ncol(AZX.W.DeleteCol))]
  D.AW = (A1ZX.W.DeleteCol%*%par.w)-(A0ZX.W.DeleteCol%*%par.w)
  D.ZW = (AZ1X.W.DeleteCol%*%par.w)-(AZ0X.W.DeleteCol%*%par.w);#mean(D.ZW);mean(eWAZ1X-eWAZ0X)
  D.ZWa2 = (A2Z1X.W.DeleteCol%*%par.w)-(A2Z0X.W.DeleteCol%*%par.w);#mean(D.ZWa2);mean(eWA2Z1X-eWA2Z0X)
  D.AWz1 = (A1Z1X.W.DeleteCol%*%par.w)-(A0Z1X.W.DeleteCol%*%par.w)
  D.AWz0 = (A1Z0X.W.DeleteCol%*%par.w)-(A0Z0X.W.DeleteCol%*%par.w)
  
  D.AY = (2*A-1)*Y/(fAZX)#mean(D.AY);mean(eYA1ZX-eYA0ZX)
  D.ZY = (2*Z-1)*Y/(fZAX) 
  
  ##need to reweight D.ZY to 1-A space
  fAX = PAX*A+(1-PAX)*(1-A); fA2X = PAX*A2+(1-PAX)*(1-A2)
  condE.D.AW = (D.AWz1*PZA2X+ D.AWz0*(1-PZA2X)) * fA2X/fAX
  #**********************************************
  delta = par[ncol(AX.Z)+ncol(X.A)+ncol(AZX.W.DeleteCol)+1]
  delta.naive = par[ncol(AX.Z)+ncol(X.A)+ncol(AZX.W.DeleteCol)+2]
  delta.bias = par[ncol(AX.Z)+ncol(X.A)+ncol(AZX.W.DeleteCol)+3]
  
  U1 = c(Z-(AX.Z%*%par.z)) *(AX.Z)
  U2 = c(A-(X.A%*%par.a)) *(X.A)
  U3 = c(W-(AZX.W.DeleteCol%*%par.w)) *(AZX.W.center)
  U4 = rep(delta,n) - (  D.AY - D.ZY /D.ZW *condE.D.AW  )
  U5 = rep(delta.naive,n) - (  D.AY )
  U6 = rep(delta.bias,n) - (  D.ZY /D.ZW *condE.D.AW  )
  U = cbind(U1,U2,U3,U4,U5,U6)
  return(U)
}
Delta.ipw = function(par,A,Z,Y,W,Age,Sex,MZ.wrong,MW.wrong){
  n=length(A)
  ### A, Z
  A1=rep(1,n);A0=rep(0,n);A2=1-A
  Z1=rep(1,n);Z0=rep(0,n);
  data = data.frame(A0,A1,A2,A,Z0,Z1,Z,Age,Sex)
  
  X.A = model.matrix(A~Age*Sex)
  ### f(A,Z)
  if(MZ.wrong){
    if(omit.sex){
      AX.Z = model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Z",T,T))
      A0X.Z = model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Z",T,T))
      A1X.Z = model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Z",T,T))
      A2X.Z = model.matrix(data=data,object=get.formula(F,"A2",NULL,NULL,"Z",T,T))
      par.z = lm(get.formula(F,"A",NULL,NULL,"Z",T,T))$coefficients
    }else{
      AX.Z = model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Z",T,F))
      A0X.Z = model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Z",T,F))
      A1X.Z = model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Z",T,F))
      A2X.Z = model.matrix(data=data,object=get.formula(F,"A2",NULL,NULL,"Z",T,F))
      par.z = lm(get.formula(F,"A",NULL,NULL,"Z",T,F))$coefficients
    }
  }else{
    AX.Z = model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Z",F,F))
    A0X.Z = model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Z",F,F))
    A1X.Z = model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Z",F,F))
    A2X.Z = model.matrix(data=data,object=get.formula(F,"A2",NULL,NULL,"Z",F,F))
    par.z = lm(get.formula(F,"A",NULL,NULL,"Z",F,F))$coefficients
  }
  PZAX = (AX.Z%*%par.z); PZA1X = (A1X.Z%*%par.z); PZA0X = (A0X.Z%*%par.z); PZA2X = (A2X.Z%*%par.z);
  fZAX = PZAX*Z+(1-PZAX)*(1-Z)
  
  par.a = lm(A~Age*Sex)$coefficients
  PAX = X.A%*%par.a; PZX = PZA1X*PAX+PZA0X*(1-PAX)
  PAZX = (PZA1X*Z+(1-PZA1X)*(1-Z)) * (PAX)  / (PZX*Z+(1-PZX)*(1-Z)) #equal to predict(lm(A~Z*Age*Sex))
  fAZX = PAZX*A+(1-PAZX)*(1-A)
  
  ###construct EE for par.w using par.a, par.z 
  Zc = Z-PZX #predict(lm(Z~Age*Sex)) ##Z.center E[Z|X]
  Ac = A-PAX #predict(lm(A~Age*Sex))##A.center E[A|X]
  AZc = A*Z-PZA1X*PAX#predict(lm(A*Z~Age*Sex))
  data = data.frame(A0,A1,A2,A,Z0,Z1,Z,Age,Sex,Ac,Zc,AZc)
  
  ## below are AZX.W delete Age, Sex, Age:Sex
  if(MW.wrong){
    if(omit.sex){
      # par.w_bl = lm(W~  Sex, data=data.frame(Age,Sex,W)[A==0&Z==0,])$coefficients
      X.W = model.matrix(W~Sex)
      AZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z","A*Z","W",T,T))
      AZ1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z1","A*Z1","W",T,T))
      AZ0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z0","A*Z0","W",T,T))
      A1ZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z","A1*Z","W",T,T))
      A0ZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z","A0*Z","W",T,T))
      A1Z1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z1","A1*Z1","W",T,T))
      A1Z0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z0","A1*Z0","W",T,T))
      A0Z1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z1","A0*Z1","W",T,T))
      A0Z0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z0","A0*Z0","W",T,T))
      
      A2Z1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A2","Z1","A2*Z1","W",T,T))
      A2Z0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A2","Z0","A2*Z0","W",T,T))
      AZX.W.center = model.matrix(data=data,object=get.formula(T,"Ac","Zc","AZc","W",T,T)) #each column should have mean zero
    }else{
      # par.w_bl = lm(W~  (Age+Sex), data=data.frame(Age,Sex,W)[A==0&Z==0,])$coefficients
      X.W = model.matrix(W~(Age+Sex))
      AZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z","A*Z","W",T,F))
      AZ1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z1","A*Z1","W",T,F))
      AZ0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z0","A*Z0","W",T,F))
      A1ZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z","A1*Z","W",T,F))
      A0ZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z","A0*Z","W",T,F))
      A1Z1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z1","A1*Z1","W",T,F))
      A1Z0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z0","A1*Z0","W",T,F))
      A0Z1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z1","A0*Z1","W",T,F))
      A0Z0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z0","A0*Z0","W",T,F))
      
      A2Z1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A2","Z1","A2*Z1","W",T,F))
      A2Z0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A2","Z0","A2*Z0","W",T,F))
      AZX.W.center = model.matrix(data=data,object=get.formula(T,"Ac","Zc","AZc","W",T,F)) #each column should have mean zero
    }
  }else{
    # par.w_bl = lm(W~  Age*Sex, data=data.frame(Age,Sex,W)[A==0&Z==0,])$coefficients
    X.W = model.matrix(W~Age*Sex)
    AZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z","A*Z","W",F,F))
    AZ1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z1","A*Z1","W",F,F))
    AZ0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z0","A*Z0","W",F,F))
    A1ZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z","A1*Z","W",F,F))
    A0ZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z","A0*Z","W",F,F))
    A1Z1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z1","A1*Z1","W",F,F))
    A1Z0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z0","A1*Z0","W",F,F))
    A0Z1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z1","A0*Z1","W",F,F))
    A0Z0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z0","A0*Z0","W",F,F))
    
    A2Z1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A2","Z1","A2*Z1","W",F,F))
    A2Z0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A2","Z0","A2*Z0","W",F,F))
    AZX.W.center = model.matrix(data=data,object=get.formula(T,"Ac","Zc","AZc","W",F,F)) #each column should have mean zero
  }
  obj = function(par.w){
    U = c(W-(AZX.W.DeleteCol%*%par.w)) *(AZX.W.center)
    return(sum(apply(U,2,sum)^2))
  }
  par.w=nlm(f=obj,p=rep(0.1,ncol(AZX.W.DeleteCol)),iterlim=1000,ndigit =20,gradtol=1e-15,stepmax=1000,steptol=1e-15)$estimate
  # rbind(c(par.w),lm(W~A*Z*Age*Sex)$coef[match(colnames(AZX.W.DeleteCol),names(lm(W~A*Z*Age*Sex)$coef))])
  D.AW = (A1ZX.W.DeleteCol%*%par.w)-(A0ZX.W.DeleteCol%*%par.w)
  D.ZW = (AZ1X.W.DeleteCol%*%par.w)-(AZ0X.W.DeleteCol%*%par.w);#mean(D.ZW);mean(eWAZ1X-eWAZ0X)
  D.ZWa2 = (A2Z1X.W.DeleteCol%*%par.w)-(A2Z0X.W.DeleteCol%*%par.w);#mean(D.ZWa2);mean(eWA2Z1X-eWA2Z0X)
  D.AWz1 = (A1Z1X.W.DeleteCol%*%par.w)-(A0Z1X.W.DeleteCol%*%par.w)
  D.AWz0 = (A1Z0X.W.DeleteCol%*%par.w)-(A0Z0X.W.DeleteCol%*%par.w)
  
  D.AY = (2*A-1)*Y/(fAZX)#mean(D.AY);mean(eYA1ZX-eYA0ZX)
  D.ZY = (2*Z-1)*Y/(fZAX) 
  
  ##need to reweight D.ZY to 1-A space
  fAX = PAX*A+(1-PAX)*(1-A); fA2X = PAX*A2+(1-PAX)*(1-A2)
  condE.D.AW = (D.AWz1*PZA2X+ D.AWz0*(1-PZA2X)) * fA2X/fAX
  # condE.D.AW =(D.AWz1*PZX+ D.AWz0*(1-PZX))/fAX - (D.AWz1*PZAX+ D.AWz0*(1-PZAX)) #is the same
  (delta=mean(  D.AY - D.ZY /D.ZW *condE.D.AW  ))
  (delta.naive=mean(  D.AY  ))
  (delta.bias=mean(  D.ZY /D.ZW *condE.D.AW  ))
  
  par = c(par.z,par.a,par.w,delta,delta.naive,delta.bias)
  return(par)
}

IPW = function(aa,zz,ww,yy,Age,Sex,MZ.wrong,MW.wrong){
  est=Delta.ipw(par=NA,A=aa,Z=zz,Y=yy,W=ww,Age=Age,Sex=Sex,MZ.wrong,MW.wrong)
  G.ipw = function(par,A,Z,Y,W,Age,Sex,MZ.wrong,MW.wrong){return(apply(U.ipw(par,A,Z,Y,W,Age,Sex,MZ.wrong,MW.wrong),2,sum))}
  bread=numDeriv::jacobian(func=G.ipw,x=est,A=aa,Z=zz,Y=yy,W=ww,Age=Age,Sex=Sex,MZ.wrong=MZ.wrong,MW.wrong=MW.wrong)
  meat.half=U.ipw(par=est,A=aa,Z=zz,Y=yy,W=ww,Age=Age,Sex=Sex,MZ.wrong=MZ.wrong,MW.wrong=MW.wrong)
  IF = meat.half%*%t(solve(-bread))
  delta = est[length(est)-2]
  delta.naive = est[length(est)-1]
  delta.bias = est[length(est)]
  (delta.var = sum(IF[,ncol(IF)-2]^2))
  (delta.naive.var = sum(IF[,ncol(IF)-1]^2))
  (delta.bias.var = sum(IF[,ncol(IF)]^2))
  
  return(
    list(
      ATE=as.numeric(delta),delta.naive=as.numeric(delta.naive),delta.bias=as.numeric(delta.bias),
      var=as.numeric(delta.var),delta.naive.var=as.numeric(delta.naive.var),delta.bias.var=as.numeric(delta.bias.var)
    )
  )
}

U.or = function(par,A,Z,Y,W,Age,Sex,MW.wrong,MR.wrong,MY.wrong){
  n=length(A)
  ### A, Z
  A1=rep(1,n);A0=rep(0,n);A2=1-A
  Z1=rep(1,n);Z0=rep(0,n);
  data = data.frame(A0,A1,A2,A,Z0,Z1,Z,Age,Sex)
  ####W
  if(MW.wrong){
    if(omit.sex){
      # par.w = lm(W~A*Z*Sex)$coefficients
      AZX.W = model.matrix(data=data,object=get.formula(F,"A","Z","A*Z","W",T,T))
      AZ1X.W = model.matrix(data=data,object=get.formula(F,"A","Z1","A*Z1","W",T,T))
      AZ0X.W = model.matrix(data=data,object=get.formula(F,"A","Z0","A*Z0","W",T,T))
      A1ZX.W = model.matrix(data=data,object=get.formula(F,"A1","Z","A1*Z","W",T,T))
      A0ZX.W = model.matrix(data=data,object=get.formula(F,"A0","Z","A0*Z","W",T,T))
      A1Z1X.W = model.matrix(data=data,object=get.formula(F,"A1","Z1","A1*Z1","W",T,T))
      A1Z0X.W = model.matrix(data=data,object=get.formula(F,"A1","Z0","A1*Z0","W",T,T))
      A0Z1X.W = model.matrix(data=data,object=get.formula(F,"A0","Z1","A0*Z1","W",T,T))
      A0Z0X.W = model.matrix(data=data,object=get.formula(F,"A0","Z0","A0*Z0","W",T,T))
    }else{
      # par.w = lm(W~A*Z*(Age+Sex))$coefficients
      AZX.W = model.matrix(data=data,object=get.formula(F,"A","Z","A*Z","W",T,F))
      AZ1X.W = model.matrix(data=data,object=get.formula(F,"A","Z1","A*Z1","W",T,F))
      AZ0X.W = model.matrix(data=data,object=get.formula(F,"A","Z0","A*Z0","W",T,F))
      A1ZX.W = model.matrix(data=data,object=get.formula(F,"A1","Z","A1*Z","W",T,F))
      A0ZX.W = model.matrix(data=data,object=get.formula(F,"A0","Z","A0*Z","W",T,F))
      A1Z1X.W = model.matrix(data=data,object=get.formula(F,"A1","Z1","A1*Z1","W",T,F))
      A1Z0X.W = model.matrix(data=data,object=get.formula(F,"A1","Z0","A1*Z0","W",T,F))
      A0Z1X.W = model.matrix(data=data,object=get.formula(F,"A0","Z1","A0*Z1","W",T,F))
      A0Z0X.W = model.matrix(data=data,object=get.formula(F,"A0","Z0","A0*Z0","W",T,F))
    }
    
  }else{
    # par.w = lm(W~A*Z*Age*Sex)$coefficients
    AZX.W = model.matrix(data=data,object=get.formula(F,"A","Z","A*Z","W",F,F))
    AZ1X.W = model.matrix(data=data,object=get.formula(F,"A","Z1","A*Z1","W",F,F))
    AZ0X.W = model.matrix(data=data,object=get.formula(F,"A","Z0","A*Z0","W",F,F))
    A1ZX.W = model.matrix(data=data,object=get.formula(F,"A1","Z","A1*Z","W",F,F))
    A0ZX.W = model.matrix(data=data,object=get.formula(F,"A0","Z","A0*Z","W",F,F))
    A1Z1X.W = model.matrix(data=data,object=get.formula(F,"A1","Z1","A1*Z1","W",F,F))
    A1Z0X.W = model.matrix(data=data,object=get.formula(F,"A1","Z0","A1*Z0","W",F,F))
    A0Z1X.W = model.matrix(data=data,object=get.formula(F,"A0","Z1","A0*Z1","W",F,F))
    A0Z0X.W = model.matrix(data=data,object=get.formula(F,"A0","Z0","A0*Z0","W",F,F))
  }
  par.w = par[1:ncol(AZ1X.W)]
  D.ZW = (AZ1X.W%*%par.w)-(AZ0X.W%*%par.w)
  D.ZWa1 = (A1Z1X.W%*%par.w)-(A1Z0X.W%*%par.w)
  D.ZWa0 = (A0Z1X.W%*%par.w)-(A0Z0X.W%*%par.w)
  D.AW = (A1ZX.W%*%par.w)-(A0ZX.W%*%par.w)
  ####Y
  if(MY.wrong){
    if(omit.sex){
      AZ0X.Y.DeletCol =  model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Y",T,T))
      A1Z0X.Y.DeletCol = model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Y",T,T))
      A0Z0X.Y.DeletCol = model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Y",T,T))
    }else{
      AZ0X.Y.DeletCol =  model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Y",T,F))
      A1Z0X.Y.DeletCol = model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Y",T,F))
      A0Z0X.Y.DeletCol = model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Y",T,F))
    }
  }else{
    AZ0X.Y.DeletCol =  model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Y",F,F))
    A1Z0X.Y.DeletCol = model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Y",F,F))
    A0Z0X.Y.DeletCol = model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Y",F,F))
  }
  par.y = par[(ncol(AZ1X.W)+1):(ncol(AZ1X.W)+ncol(AZ0X.Y.DeletCol))]
  ####R
  EYAZ0X = AZ0X.Y.DeletCol%*%par.y
  EWAZ0X = AZ0X.W%*%par.w
  if(MR.wrong){
    if(omit.sex){
      AX.R = model.matrix(Y~A*Sex)
      A0X.R = model.matrix(Y~A0*Sex)
      A1X.R = model.matrix(Y~A1*Sex)
      A2X.R = model.matrix(Y~A2*Sex)
    }else{
      AX.R = model.matrix(Y~A*(Age+Sex))
      A0X.R = model.matrix(Y~A0*(Age+Sex))
      A1X.R = model.matrix(Y~A1*(Age+Sex))
      A2X.R = model.matrix(Y~A2*(Age+Sex))
    }
  }else{
    AX.R = model.matrix(Y~A*Age*Sex)
    A0X.R = model.matrix(Y~A0*Age*Sex)
    A1X.R = model.matrix(Y~A1*Age*Sex)
    A2X.R = model.matrix(Y~A2*Age*Sex)
  }
  par.r = par[(ncol(AZ1X.W)+ncol(AZ0X.Y.DeletCol)+1):(ncol(AZ1X.W)+ncol(AZ0X.Y.DeletCol)+ncol(AX.R))]
  EYA1ZX = A1Z0X.Y.DeletCol%*%par.y + (A1X.R%*%par.r) * D.ZWa1 * Z
  EYA0ZX = A0Z0X.Y.DeletCol%*%par.y + (A0X.R%*%par.r) * D.ZWa0 * Z
  D.AY = EYA1ZX - EYA0ZX
  #********************************************************************
  delta = par[ncol(AZ1X.W)+ncol(AZ0X.Y.DeletCol)+ncol(AX.R)+1] 
  delta.naive = par[ncol(AZ1X.W)+ncol(AZ0X.Y.DeletCol)+ncol(AX.R)+2]
  delta.bias = par[ncol(AZ1X.W)+ncol(AZ0X.Y.DeletCol)+ncol(AX.R)+3]
  
  U1 = c(W-(AZX.W%*%par.w)) *(AZX.W)
  U2 = c(Y-(AZ0X.Y.DeletCol%*%par.y)) *(AZ0X.Y.DeletCol) *as.numeric(Z==0)
  U3 = c(  Y - EYAZ0X - (AX.R%*%par.r)*(W-EWAZ0X)  )*(AX.R)
  U4 = delta - (   D.AY - (A2X.R%*%par.r)*D.AW   )
  U5 = delta.naive - (   D.AY   )
  U6 = delta.bias - (   (A2X.R%*%par.r)*D.AW   )
  U = cbind(U1,U2,U3,U4,U5,U6)
  return(U)
}
Delta.or = function(par,A,Z,Y,W,Age,Sex,MW.wrong,MR.wrong,MY.wrong){
  n=length(A)
  ### A, Z
  A1=rep(1,n);A0=rep(0,n);A2=1-A
  Z1=rep(1,n);Z0=rep(0,n);
  data = data.frame(A0,A1,A2,A,Z0,Z1,Z,Age,Sex)
  ####W
  if(MW.wrong){
    if(omit.sex){
      par.w=lm(get.formula(F,"A","Z","A*Z","W",T,T))$coefficients
      AZX.W = model.matrix(data=data,object=get.formula(F,"A","Z","A*Z","W",T,T))
      AZ1X.W = model.matrix(data=data,object=get.formula(F,"A","Z1","A*Z1","W",T,T))
      AZ0X.W = model.matrix(data=data,object=get.formula(F,"A","Z0","A*Z0","W",T,T))
      A1ZX.W = model.matrix(data=data,object=get.formula(F,"A1","Z","A1*Z","W",T,T))
      A0ZX.W = model.matrix(data=data,object=get.formula(F,"A0","Z","A0*Z","W",T,T))
      A1Z1X.W = model.matrix(data=data,object=get.formula(F,"A1","Z1","A1*Z1","W",T,T))
      A1Z0X.W = model.matrix(data=data,object=get.formula(F,"A1","Z0","A1*Z0","W",T,T))
      A0Z1X.W = model.matrix(data=data,object=get.formula(F,"A0","Z1","A0*Z1","W",T,T))
      A0Z0X.W = model.matrix(data=data,object=get.formula(F,"A0","Z0","A0*Z0","W",T,T))
    }else{
      par.w=lm(get.formula(F,"A","Z","A*Z","W",T,F))$coefficients
      AZX.W = model.matrix(data=data,object=get.formula(F,"A","Z","A*Z","W",T,F))
      AZ1X.W = model.matrix(data=data,object=get.formula(F,"A","Z1","A*Z1","W",T,F))
      AZ0X.W = model.matrix(data=data,object=get.formula(F,"A","Z0","A*Z0","W",T,F))
      A1ZX.W = model.matrix(data=data,object=get.formula(F,"A1","Z","A1*Z","W",T,F))
      A0ZX.W = model.matrix(data=data,object=get.formula(F,"A0","Z","A0*Z","W",T,F))
      A1Z1X.W = model.matrix(data=data,object=get.formula(F,"A1","Z1","A1*Z1","W",T,F))
      A1Z0X.W = model.matrix(data=data,object=get.formula(F,"A1","Z0","A1*Z0","W",T,F))
      A0Z1X.W = model.matrix(data=data,object=get.formula(F,"A0","Z1","A0*Z1","W",T,F))
      A0Z0X.W = model.matrix(data=data,object=get.formula(F,"A0","Z0","A0*Z0","W",T,F))
    }
  }else{
    par.w=lm(get.formula(F,"A","Z","A*Z","W",F,F))$coefficients
    AZX.W = model.matrix(data=data,object=get.formula(F,"A","Z","A*Z","W",F,F))
    AZ1X.W = model.matrix(data=data,object=get.formula(F,"A","Z1","A*Z1","W",F,F))
    AZ0X.W = model.matrix(data=data,object=get.formula(F,"A","Z0","A*Z0","W",F,F))
    A1ZX.W = model.matrix(data=data,object=get.formula(F,"A1","Z","A1*Z","W",F,F))
    A0ZX.W = model.matrix(data=data,object=get.formula(F,"A0","Z","A0*Z","W",F,F))
    A1Z1X.W = model.matrix(data=data,object=get.formula(F,"A1","Z1","A1*Z1","W",F,F))
    A1Z0X.W = model.matrix(data=data,object=get.formula(F,"A1","Z0","A1*Z0","W",F,F))
    A0Z1X.W = model.matrix(data=data,object=get.formula(F,"A0","Z1","A0*Z1","W",F,F))
    A0Z0X.W = model.matrix(data=data,object=get.formula(F,"A0","Z0","A0*Z0","W",F,F))
  }
  D.ZW = (AZ1X.W%*%par.w)-(AZ0X.W%*%par.w)
  D.ZWa1 = (A1Z1X.W%*%par.w)-(A1Z0X.W%*%par.w)
  D.ZWa0 = (A0Z1X.W%*%par.w)-(A0Z0X.W%*%par.w)
  D.AW = (A1ZX.W%*%par.w)-(A0ZX.W%*%par.w)
  ####Y
  if(MY.wrong){
    if(omit.sex){
      par.y = lm(get.formula(F,"A",NULL,NULL,"Y",T,T), data=data.frame(A,Age,Sex,Y)[Z==0,])$coefficients
      AZ0X.Y.DeletCol =  model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Y",T,T))
      A1Z0X.Y.DeletCol = model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Y",T,T))
      A0Z0X.Y.DeletCol = model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Y",T,T))
    }else{
      par.y = lm(get.formula(F,"A",NULL,NULL,"Y",T,F), data=data.frame(A,Age,Sex,Y)[Z==0,])$coefficients
      AZ0X.Y.DeletCol =  model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Y",T,F))
      A1Z0X.Y.DeletCol = model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Y",T,F))
      A0Z0X.Y.DeletCol = model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Y",T,F))
    }
  }else{
    par.y = lm(get.formula(F,"A",NULL,NULL,"Y",F,F), data=data.frame(A,Age,Sex,Y)[Z==0,])$coefficients
    AZ0X.Y.DeletCol =  model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Y",F,F))
    A1Z0X.Y.DeletCol = model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Y",F,F))
    A0Z0X.Y.DeletCol = model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Y",F,F))
  }
  
  ####R
  EYAZ0X = AZ0X.Y.DeletCol%*%par.y
  EWAZ0X = AZ0X.W%*%par.w
  if(MR.wrong){
    if(omit.sex){
      AX.R = model.matrix(Y~A*Sex)
      A0X.R = model.matrix(Y~A0*Sex)
      A1X.R = model.matrix(Y~A1*Sex)
      A2X.R = model.matrix(Y~A2*Sex)
    }else{
      AX.R = model.matrix(Y~A*(Age+Sex))
      A0X.R = model.matrix(Y~A0*(Age+Sex))
      A1X.R = model.matrix(Y~A1*(Age+Sex))
      A2X.R = model.matrix(Y~A2*(Age+Sex))
    }
  }else{
    AX.R = model.matrix(Y~A*Age*Sex)
    A0X.R = model.matrix(Y~A0*Age*Sex)
    A1X.R = model.matrix(Y~A1*Age*Sex)
    A2X.R = model.matrix(Y~A2*Age*Sex)
  }
  obj = function(par.r){
    U = c(  (Y-EYAZ0X) - (AX.R%*%par.r)*(W-EWAZ0X)  )*(AX.R)
    return(sum(apply(U,2,sum)^2))
  }
  t= nlm(f=obj,p=rep(1,ncol(AX.R)),iterlim=1000,ndigit =20,gradtol=1e-15,stepmax=1000,steptol=1e-15)
  par.r = t$estimate
  # rbind(par.r,lm((EYAZ1X-EYAZ0X)/(EWAZ1X-EWAZ0X)~A*Age*Sex)$coef)
  EYA1ZX = A1Z0X.Y.DeletCol%*%par.y + (A1X.R%*%par.r) * D.ZWa1 * Z
  EYA0ZX = A0Z0X.Y.DeletCol%*%par.y + (A0X.R%*%par.r) * D.ZWa0 * Z
  D.AY = EYA1ZX - EYA0ZX
  delta = mean(   D.AY - (A2X.R%*%par.r)*D.AW   )
  delta.naive = mean(   D.AY   )
  delta.bias = mean(   (A2X.R%*%par.r)*D.AW   )
  
  par=c(par.w,par.y,par.r,delta,delta.naive,delta.bias)
  return(par)
}

OR = function(aa,zz,ww,yy,Age,Sex,MW.wrong,MR.wrong,MY.wrong){
  est=Delta.or(par=NA,A=aa,Z=zz,Y=yy,W=ww,Age=Age,Sex=Sex,MW.wrong=MW.wrong,MR.wrong=MR.wrong,MY.wrong=MY.wrong)
  G.or = function(par,A,Z,Y,W,Age,Sex,MW.wrong,MR.wrong,MY.wrong){return(apply(U.or(par,A,Z,Y,W,Age,Sex,MW.wrong,MR.wrong,MY.wrong),2,sum))}
  bread=numDeriv::jacobian(func=G.or,x=est,A=aa,Z=zz,Y=yy,W=ww,Age=Age,Sex=Sex,MW.wrong=MW.wrong,MR.wrong=MR.wrong,MY.wrong=MY.wrong)
  meat.half=U.or(par=est,A=aa,Z=zz,Y=yy,W=ww,Age=Age,Sex=Sex,MW.wrong=MW.wrong,MR.wrong=MR.wrong,MY.wrong=MY.wrong)
  IF = meat.half%*%t(solve(-bread))
  delta = est[length(est)-2]
  delta.var = sum(IF[,ncol(IF)-2]^2)
  delta.naive = est[length(est)-1]
  delta.naive.var = sum(IF[,ncol(IF)-1]^2)
  delta.bias = est[length(est)]
  delta.bias.var = sum(IF[,ncol(IF)]^2)
  
  return(
    list(
      ATE=as.numeric(delta),delta.naive=as.numeric(delta.naive),delta.bias=as.numeric(delta.bias),
      var=as.numeric(delta.var),delta.naive.var=as.numeric(delta.naive.var),delta.bias.var=as.numeric(delta.bias.var)
    )
  )
}

U.miao = function(par,A,Z,Y,W,Age,Sex,MW.wrong,MY.wrong){
  n=length(A)
  ### A, Z
  A1=rep(1,n);A0=rep(0,n);A2=1-A
  Z1=rep(1,n);Z0=rep(0,n);
  data = data.frame(A0,A1,A2,A,Z0,Z1,Z,Age,Sex)
  ### model data
  if(MW.wrong){
    if(omit.sex){
      AZX.W = model.matrix(data=data,object=get.formula(F,"A","Z","A*Z","W",T,T))
      A1ZX.W = model.matrix(data=data,object=get.formula(F,"A1","Z","A1*Z","W",T,T))
      A0ZX.W = model.matrix(data=data,object=get.formula(F,"A0","Z","A0*Z","W",T,T))
      
      A2Z1X.W = model.matrix(data=data,object=get.formula(F,"A2","Z1","A2*Z1","W",T,T))
      A2Z0X.W = model.matrix(data=data,object=get.formula(F,"A2","Z0","A2*Z0","W",T,T))
    }else{
      AZX.W = model.matrix(data=data,object=get.formula(F,"A","Z","A*Z","W",T,F))
      A1ZX.W = model.matrix(data=data,object=get.formula(F,"A1","Z","A1*Z","W",T,F))
      A0ZX.W = model.matrix(data=data,object=get.formula(F,"A0","Z","A0*Z","W",T,F))
      
      A2Z1X.W = model.matrix(data=data,object=get.formula(F,"A2","Z1","A2*Z1","W",T,F))
      A2Z0X.W = model.matrix(data=data,object=get.formula(F,"A2","Z0","A2*Z0","W",T,F))
    }
  }else{
    AZX.W = model.matrix(data=data,object=get.formula(F,"A","Z","A*Z","W",F,F))
    A1ZX.W = model.matrix(data=data,object=get.formula(F,"A1","Z","A1*Z","W",F,F))
    A0ZX.W = model.matrix(data=data,object=get.formula(F,"A0","Z","A0*Z","W",F,F))
    
    A2Z1X.W = model.matrix(data=data,object=get.formula(F,"A2","Z1","A2*Z1","W",F,F))
    A2Z0X.W = model.matrix(data=data,object=get.formula(F,"A2","Z0","A2*Z0","W",F,F))
  }
  
  if(MY.wrong){
    if(omit.sex){
      par.y = lm(Y~  A*Sex, data=data.frame(A,Age,Sex,Y)[Z==0,])$coefficients
      AZ0X.Y.DeletCol =  model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Y",T,T))
      A1Z0X.Y.DeletCol = model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Y",T,T))
      A0Z0X.Y.DeletCol = model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Y",T,T))
    }else{
      par.y = lm(Y~  A*(Age+Sex), data=data.frame(A,Age,Sex,Y)[Z==0,])$coefficients
      AZ0X.Y.DeletCol =  model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Y",T,F))
      A1Z0X.Y.DeletCol = model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Y",T,F))
      A0Z0X.Y.DeletCol = model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Y",T,F))
    }
  }else{
    par.y = lm(Y~  A*(Age+Sex), data=data.frame(A,Age,Sex,Y)[Z==0,])$coefficients
    AZ0X.Y.DeletCol =  model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Y",F,F))
    A1Z0X.Y.DeletCol = model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Y",F,F))
    A0Z0X.Y.DeletCol = model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Y",F,F))
  }
  
  if(MY.wrong){
    if(omit.sex){
      AZX.Y = model.matrix(data=data,object=get.formula(F,"A","Z","A*Z","Y",T,T))
      A0ZX.Y = model.matrix(data=data,object=get.formula(F,"A0","Z","A0*Z","Y",T,T))
      A1ZX.Y = model.matrix(data=data,object=get.formula(F,"A1","Z","A1*Z","Y",T,T))
      A2Z0X.Y = model.matrix(data=data,object=get.formula(F,"A2","Z0","A2*Z0","Y",T,T))
      A2Z1X.Y = model.matrix(data=data,object=get.formula(F,"A2","Z1","A3*Z1","Y",T,T))
    }else{
      AZX.Y = model.matrix(data=data,object=get.formula(F,"A","Z","A*Z","Y",T,F))
      A0ZX.Y = model.matrix(data=data,object=get.formula(F,"A0","Z","A0*Z","Y",T,F))
      A1ZX.Y = model.matrix(data=data,object=get.formula(F,"A1","Z","A1*Z","Y",T,F))
      A2Z0X.Y = model.matrix(data=data,object=get.formula(F,"A2","Z0","A2*Z0","Y",T,F))
      A2Z1X.Y = model.matrix(data=data,object=get.formula(F,"A2","Z1","A3*Z1","Y",T,F))
    }
  }else{
    AZX.Y = model.matrix(data=data,object=get.formula(F,"A","Z","A*Z","Y",F,F))
    A0ZX.Y = model.matrix(data=data,object=get.formula(F,"A0","Z","A0*Z","Y",F,F))
    A1ZX.Y = model.matrix(data=data,object=get.formula(F,"A1","Z","A1*Z","Y",F,F))
    A2Z0X.Y = model.matrix(data=data,object=get.formula(F,"A2","Z0","A2*Z0","Y",F,F))
    A2Z1X.Y = model.matrix(data=data,object=get.formula(F,"A2","Z1","A2*Z1","Y",F,F))
  }
  ### par
  par.w = par[1:ncol(AZX.W)]
  par.y = par[(ncol(AZX.W)+1):(ncol(AZX.W)+ncol(AZX.Y))]
  delta = par[(ncol(AZX.W)+ncol(AZX.Y)+1)]
  delta.naive = par[(ncol(AZX.W)+ncol(AZX.Y)+2)]
  delta.bias = par[(ncol(AZX.W)+ncol(AZX.Y)+3)]
  ### deltas
  D.AW = (A1ZX.W%*%par.w)-(A0ZX.W%*%par.w)
  D.ZWa2 = (A2Z1X.W%*%par.w)-(A2Z0X.W%*%par.w)
  D.AY = (A1ZX.Y%*%par.y)-(A0ZX.Y%*%par.y)
  D.ZY2 = (A2Z1X.Y%*%par.y)-(A2Z0X.Y%*%par.y)
  
  
  ### Estimating equation
  U1 = c(W-(AZX.W%*%par.w)) *(AZX.W)
  U2 = c(Y-(AZX.Y%*%par.y)) *(AZX.Y)
  U3 = delta - (   D.AY - D.ZY2/D.ZWa2*D.AW   )
  U4 = delta.naive - (   D.AY   )
  U5 = delta.bias - (   D.ZY2/D.ZWa2*D.AW   )
  
  U = cbind(U1,U2,U3,U4,U5)
  return(U)
}
Delta.miao = function(par,A,Z,Y,W,Age,Sex,MW.wrong,MY.wrong){
  n=length(A)
  ### A, Z
  A1=rep(1,n);A0=rep(0,n);A2=1-A
  Z1=rep(1,n);Z0=rep(0,n);
  data = data.frame(A0,A1,A2,A,Z0,Z1,Z,Age,Sex)
  ### model data
  if(MW.wrong){
    if(omit.sex){
      par.w = lm(get.formula(F,"A","Z","A*Z","W",T,T))$coefficients
      AZX.W = model.matrix(data=data,object=get.formula(F,"A","Z","A*Z","W",T,T))
      A1ZX.W = model.matrix(data=data,object=get.formula(F,"A1","Z","A1*Z","W",T,T))
      A0ZX.W = model.matrix(data=data,object=get.formula(F,"A0","Z","A0*Z","W",T,T))
      
      A2Z1X.W = model.matrix(data=data,object=get.formula(F,"A2","Z1","A2*Z1","W",T,T))
      A2Z0X.W = model.matrix(data=data,object=get.formula(F,"A2","Z0","A2*Z0","W",T,T))
    }else{
      par.w = lm(get.formula(F,"A","Z","A*Z","W",T,F))$coefficients
      AZX.W = model.matrix(data=data,object=get.formula(F,"A","Z","A*Z","W",T,F))
      A1ZX.W = model.matrix(data=data,object=get.formula(F,"A1","Z","A1*Z","W",T,F))
      A0ZX.W = model.matrix(data=data,object=get.formula(F,"A0","Z","A0*Z","W",T,F))
      
      A2Z1X.W = model.matrix(data=data,object=get.formula(F,"A2","Z1","A2*Z1","W",T,F))
      A2Z0X.W = model.matrix(data=data,object=get.formula(F,"A2","Z0","A2*Z0","W",T,F))
    }
  }else{
    par.w = lm(get.formula(F,"A","Z","A*Z","W",F,F))$coefficients
    AZX.W = model.matrix(data=data,object=get.formula(F,"A","Z","A*Z","W",F,F))
    A1ZX.W = model.matrix(data=data,object=get.formula(F,"A1","Z","A1*Z","W",F,F))
    A0ZX.W = model.matrix(data=data,object=get.formula(F,"A0","Z","A0*Z","W",F,F))
    
    A2Z1X.W = model.matrix(data=data,object=get.formula(F,"A2","Z1","A2*Z1","W",F,F))
    A2Z0X.W = model.matrix(data=data,object=get.formula(F,"A2","Z0","A2*Z0","W",F,F))
  }
  
  if(MY.wrong){
    if(omit.sex){
      AZX.Y = model.matrix(data=data,object=get.formula(F,"A","Z","A*Z","Y",T,T))
      A0ZX.Y = model.matrix(data=data,object=get.formula(F,"A0","Z","A0*Z","Y",T,T))
      A1ZX.Y = model.matrix(data=data,object=get.formula(F,"A1","Z","A1*Z","Y",T,T))
      A2Z0X.Y = model.matrix(data=data,object=get.formula(F,"A2","Z0","A2*Z0","Y",T,T))
      A2Z1X.Y = model.matrix(data=data,object=get.formula(F,"A2","Z1","A2*Z1","Y",T,T))
      par.y = lm(get.formula(F,"A","Z","A*Z","Y",T,T))$coefficients
    }else{
      AZX.Y = model.matrix(data=data,object=get.formula(F,"A","Z","A*Z","Y",T,F))
      A0ZX.Y = model.matrix(data=data,object=get.formula(F,"A0","Z","A0*Z","Y",T,F))
      A1ZX.Y = model.matrix(data=data,object=get.formula(F,"A1","Z","A1*Z","Y",T,F))
      A2Z0X.Y = model.matrix(data=data,object=get.formula(F,"A2","Z0","A2*Z0","Y",T,F))
      A2Z1X.Y = model.matrix(data=data,object=get.formula(F,"A2","Z1","A2*Z1","Y",T,F))
      par.y = lm(get.formula(F,"A","Z","A*Z","Y",T,F))$coefficients
    }
  }else{
    AZX.Y = model.matrix(data=data,object=get.formula(F,"A","Z","A*Z","Y",F,F))
    A0ZX.Y = model.matrix(data=data,object=get.formula(F,"A0","Z","A0*Z","Y",F,F))
    A1ZX.Y = model.matrix(data=data,object=get.formula(F,"A1","Z","A1*Z","Y",F,F))
    A2Z0X.Y = model.matrix(data=data,object=get.formula(F,"A2","Z0","A2*Z0","Y",F,F))
    A2Z1X.Y = model.matrix(data=data,object=get.formula(F,"A2","Z1","A2*Z1","Y",F,F))
    par.y = lm(get.formula(F,"A","Z","A*Z","Y",F,F))$coefficients
  }
  
  ### deltas
  D.AW = (A1ZX.W%*%par.w)-(A0ZX.W%*%par.w)
  D.ZWa2 = (A2Z1X.W%*%par.w)-(A2Z0X.W%*%par.w)
  D.AY = (A1ZX.Y%*%par.y)-(A0ZX.Y%*%par.y)
  D.ZY2 = (A2Z1X.Y%*%par.y)-(A2Z0X.Y%*%par.y)
  delta = mean(D.AY - D.ZY2/D.ZWa2*D.AW)
  delta.naive = mean(D.AY)
  delta.bias = mean(D.ZY2/D.ZWa2*D.AW)
  
  par = c(par.w,par.y,delta,delta.naive,delta.bias)
  return(par)
}

MIAO = function(aa,zz,ww,yy,Age,Sex,MW.wrong,MY.wrong){
  est=Delta.miao(par=NA,A=aa,Z=zz,Y=yy,W=ww,Age=Age,Sex=Sex,MW.wrong=MW.wrong,MY.wrong=MY.wrong)
  # t=nlm(f=U.miao,p=rep(1,35),A=aa,Z=zz,Y=yy,W=ww,Age=Age,Sex=Sex,MW.wrong=MW.wrong,MY.wrong=MY.wrong)
  G.miao = function(par,A,Z,Y,W,Age,Sex,MW.wrong,MY.wrong){return(apply(U.miao(par,A,Z,Y,W,Age,Sex,MW.wrong,MY.wrong),2,sum))}
  bread=numDeriv::jacobian(func=G.miao,x=est,A=aa,Z=zz,Y=yy,W=ww,Age=Age,Sex=Sex,MW.wrong=MW.wrong,MY.wrong=MY.wrong)
  meat.half=U.miao(par=est,A=aa,Z=zz,Y=yy,W=ww,Age=Age,Sex=Sex,MW.wrong=MW.wrong,MY.wrong=MY.wrong)
  IF = meat.half%*%t(solve(-bread)) #ginv
  delta = est[length(est)-2]
  delta.var = sum(IF[,ncol(IF)-2]^2)
  delta.naive = est[length(est)-1]
  delta.naive.var = sum(IF[,ncol(IF)-1]^2)
  delta.bias = est[length(est)]
  delta.bias.var = sum(IF[,ncol(IF)]^2)
  return(
    list(
      ATE=as.numeric(delta),delta.naive=as.numeric(delta.naive),delta.bias=as.numeric(delta.bias),
      var=as.numeric(delta.var),delta.naive.var=as.numeric(delta.naive.var),delta.bias.var=as.numeric(delta.bias.var)
    )
  )
}

U.gest = function(par,A,Z,Y,W,Age,Sex,MZ.wrong,MR.wrong){
  n=length(A)
  ### A, Z
  A1=rep(1,n);A0=rep(0,n);A2=1-A
  Z1=rep(1,n);Z0=rep(0,n);
  data = data.frame(A0,A1,A2,A,Z0,Z1,Z,Age,Sex)
  X.A = model.matrix(A~Age*Sex)
  par.a = par[1:ncol(X.A)]
  if(MZ.wrong){
    if(omit.sex){
      AX.Z = model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Z",T,T))
      A0X.Z = model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Z",T,T))
      A1X.Z = model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Z",T,T))
    }else{
      AX.Z = model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Z",T,F))
      A0X.Z = model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Z",T,F))
      A1X.Z = model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Z",T,F))
    }
  }else{
    AX.Z = model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Z",F,F))
    A0X.Z = model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Z",F,F))
    A1X.Z = model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Z",F,F))
  }
  par.z = par[(ncol(X.A)+1):(ncol(X.A)+ncol(AX.Z))]
  ### f(A,Z)
  A1=rep(1,n);A0=rep(0,n);A2=1-A
  PZAX = (AX.Z%*%par.z); PZA1X = (A1X.Z%*%par.z); PZA0X = (A0X.Z%*%par.z)
  PAX = X.A%*%par.a; PZX = PZA1X*PAX+PZA0X*(1-PAX)
  PAZX = (PZA1X*Z+(1-PZA1X)*(1-Z)) * (PAX)  / (PZX*Z+(1-PZX)*(1-Z))
  fAZX = PAZX*A+(1-PAZX)*(1-A)
  
  if(MR.wrong){
    if(omit.sex){
      AX.R = model.matrix(Y~A*Sex)
      A2X.R = model.matrix(Y~A2*Sex)
      A0X.R = model.matrix(Y~A0*Sex)
      A1X.R = model.matrix(Y~A1*Sex)
    }else{
      AX.R = model.matrix(Y~A*(Age+Sex))
      A2X.R = model.matrix(Y~A2*(Age+Sex))
      A0X.R = model.matrix(Y~A0*(Age+Sex))
      A1X.R = model.matrix(Y~A1*(Age+Sex))
    }
  }else{
    AX.R = model.matrix(Y~A*Age*Sex)
    A2X.R = model.matrix(Y~A2*Age*Sex)
    A0X.R = model.matrix(Y~A0*Age*Sex)
    A1X.R = model.matrix(Y~A1*Age*Sex)
  }
  
  par.r = par[(ncol(X.A)+ncol(AX.Z)+1):(ncol(X.A)+ncol(AX.Z)+ncol(AX.R))]
  
  D.AY = (2*A-1)*Y/(fAZX) #summary(D.AY-(eYA1ZX-eYA0ZX))
  D.AW = (2*A-1)*W/(fAZX) #summary(D.AW-(eWA1ZX-eWA0ZX))
  R2_cond_ZX = (A0X.R%*%par.r*PAZX)+(A1X.R%*%par.r*(1-PAZX))
  
  delta = par[ncol(X.A)+ncol(AX.Z)+ncol(AX.R)+1] #delta = mean(  D.AY - (R2_cond_ZX)*D.AW  )
  delta.naive = par[ncol(X.A)+ncol(AX.Z)+ncol(AX.R)+2]
  delta.bias = par[ncol(X.A)+ncol(AX.Z)+ncol(AX.R)+3]
  
  U1 = c(Z-(AX.Z%*%par.z)) *(AX.Z)
  U2 = c(A-(X.A%*%par.a)) *(X.A)
  U3 = c(   ( Z-(AX.Z%*%par.z) )*(Y-W*(AX.R%*%par.r))   )*AX.R
  U4 = rep(delta,n) - (  D.AY - (R2_cond_ZX)*D.AW  )
  U5 = rep(delta.naive,n) - (  D.AY  )
  U6 = rep(delta.bias,n) - ( (R2_cond_ZX)*D.AW  )
  U = cbind(U1,U2,U3,U4,U5,U6)
  return(U)
}
Delta.gest = function(par,A,Z,Y,W,Age,Sex,MZ.wrong,MR.wrong){
  n=length(A)
  ### A, Z
  A1=rep(1,n);A0=rep(0,n);A2=1-A
  Z1=rep(1,n);Z0=rep(0,n);
  data = data.frame(A0,A1,A2,A,Z0,Z1,Z,Age,Sex)
  X.A = model.matrix(A~Age*Sex)
  par.a = lm(A~Age*Sex)$coefficients
  if(MZ.wrong){
    if(omit.sex){
      AX.Z = model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Z",T,T))
      A0X.Z = model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Z",T,T))
      A1X.Z = model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Z",T,T))
      par.z = lm(get.formula(F,"A",NULL,NULL,"Z",T,T))$coefficients
    }else{
      AX.Z = model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Z",T,F))
      A0X.Z = model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Z",T,F))
      A1X.Z = model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Z",T,F))
      par.z = lm(get.formula(F,"A",NULL,NULL,"Z",T,F))$coefficients
    }
  }else{
    AX.Z = model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Z",F,F))
    A0X.Z = model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Z",F,F))
    A1X.Z = model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Z",F,F))
    par.z = lm(get.formula(F,"A",NULL,NULL,"Z",F,F))$coefficients
  }
  
  ### f(A,Z)
  A1=rep(1,n);A0=rep(0,n);A2=1-A
  PZAX = (AX.Z%*%par.z); PZA1X = (A1X.Z%*%par.z); PZA0X = (A0X.Z%*%par.z)
  PAX = X.A%*%par.a; PZX = PZA1X*PAX+PZA0X*(1-PAX)
  PAZX = (PZA1X*Z+(1-PZA1X)*(1-Z)) * (PAX)  / (PZX*Z+(1-PZX)*(1-Z))
  fAZX = PAZX*A+(1-PAZX)*(1-A)
  
  if(MR.wrong){
    if(omit.sex){
      AX.R = model.matrix(Y~A*Sex)
      A2X.R = model.matrix(Y~A2*Sex)
      A0X.R = model.matrix(Y~A0*Sex)
      A1X.R = model.matrix(Y~A1*Sex)
    }else{
      AX.R = model.matrix(Y~A*(Age+Sex))
      A2X.R = model.matrix(Y~A2*(Age+Sex))
      A0X.R = model.matrix(Y~A0*(Age+Sex))
      A1X.R = model.matrix(Y~A1*(Age+Sex))
    }
  }else{
    AX.R = model.matrix(Y~A*Age*Sex)
    A2X.R = model.matrix(Y~A2*Age*Sex)
    A0X.R = model.matrix(Y~A0*Age*Sex)
    A1X.R = model.matrix(Y~A1*Age*Sex)
  }
  
  obj = function(par.r){
    U = c(   ( Z-(AX.Z%*%par.z) )*(Y-W*(AX.R%*%par.r))   )*AX.R
    return(sum(apply(U,2,sum)^2))
  }
  par.r= nlm(f=obj,p=rep(0,ncol(AX.R)),iterlim=1000,ndigit =20,gradtol=1e-15,stepmax=1000,steptol=1e-15)$estimate
  #rbind(par.r,lm((EYAaZ1X-EYAaZ0X)/(EWAaZ1X-EWAaZ0X)~A*Age*Sex)$coef)
  D.AY = (2*A-1)*Y/(fAZX) #summary(D.AY-(eYA1ZX-eYA0ZX))
  D.AW = (2*A-1)*W/(fAZX) #summary(D.AW-(eWA1ZX-eWA0ZX))
  R2_cond_ZX = (A0X.R%*%par.r*PAZX)+(A1X.R%*%par.r*(1-PAZX))
  delta = mean(  D.AY - (R2_cond_ZX)*D.AW  )
  delta.naive = mean(  D.AY  )
  delta.bias = mean(  (R2_cond_ZX)*D.AW  )
  
  par=c(par.a,par.z,par.r,delta,delta.naive,delta.bias)
  return(par)
}

GEST = function(aa,zz,ww,yy,Age,Sex,MZ.wrong,MR.wrong){
  est=Delta.gest(par=NA,A=aa,Z=zz,Y=yy,W=ww,Age=Age,Sex=Sex,MZ.wrong=MZ.wrong,MR.wrong=MR.wrong)
  # t=nlm(f=U.gest,p=rep(1,23),A=aa,Z=zz,Y=yy,W=ww,Age=Age,Sex=Sex,MZ.wrong=MZ.wrong,MR.wrong=MR.wrong)
  G.gest = function(par,A,Z,Y,W,Age,Sex,MZ.wrong,MR.wrong){return(apply(U.gest(par,A,Z,Y,W,Age,Sex,MZ.wrong,MR.wrong),2,sum))}
  bread=numDeriv::jacobian(func=G.gest,x=est,A=aa,Z=zz,Y=yy,W=ww,Age=Age,Sex=Sex,MZ.wrong=MZ.wrong,MR.wrong=MR.wrong)
  meat.half=U.gest(par=est,A=aa,Z=zz,Y=yy,W=ww,Age=Age,Sex=Sex,MZ.wrong=MZ.wrong,MR.wrong=MR.wrong)
  IF = meat.half%*%t(solve(-bread)) #ginv
  delta = est[length(est)-2]
  delta.var = sum(IF[,ncol(IF)-2]^2)
  delta.naive = est[length(est)-1]
  delta.naive.var = sum(IF[,ncol(IF)-1]^2)
  delta.bias = est[length(est)]
  delta.bias.var = sum(IF[,ncol(IF)]^2)
  return(
    list(
      ATE=as.numeric(delta),delta.naive=as.numeric(delta.naive),delta.bias=as.numeric(delta.bias),
      var=as.numeric(delta.var),delta.naive.var=as.numeric(delta.naive.var),delta.bias.var=as.numeric(delta.bias.var)
    )
  )
}

U.mr = function(par,A,Z,Y,W,Age,Sex,MZ.wrong,MW.wrong,MR.wrong,MY.wrong){
  n=length(A)
  ### A, Z
  A1=rep(1,n);A0=rep(0,n);A2=1-A
  Z1=rep(1,n);Z0=rep(0,n);
  data = data.frame(A0,A1,A2,A,Z0,Z1,Z,Age,Sex)
  ### f(A,Z)
  if(MZ.wrong){
    if(omit.sex){
      AX.Z = model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Z",T,T))
      A0X.Z = model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Z",T,T))
      A1X.Z = model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Z",T,T))
      A2X.Z = model.matrix(data=data,object=get.formula(F,"A2",NULL,NULL,"Z",T,T))
      # par.z = lm(Z~A*Sex)$coefficients
    }else{
      AX.Z = model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Z",T,F))
      A0X.Z = model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Z",T,F))
      A1X.Z = model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Z",T,F))
      A2X.Z = model.matrix(data=data,object=get.formula(F,"A2",NULL,NULL,"Z",T,F))
      # par.z = lm(Z~A*(Age+Sex))$coefficients
    }
  }else{
    AX.Z = model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Z",F,F))
    A0X.Z = model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Z",F,F))
    A1X.Z = model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Z",F,F))
    A2X.Z = model.matrix(data=data,object=get.formula(F,"A2",NULL,NULL,"Z",F,F))
    # par.z = lm(Z~A*Age*Sex)$coefficients
  }
  par.z=par[1:ncol(AX.Z)]
  X.A = model.matrix(A~Age*Sex)
  # par.a = lm(A~Age*Sex)$coefficients
  par.a=par[(ncol(AX.Z)+1):(ncol(AX.Z)+ncol(X.A))]
  PZAX = (AX.Z%*%par.z); PZA1X = (A1X.Z%*%par.z); PZA0X = (A0X.Z%*%par.z); PZA2X = (A2X.Z%*%par.z);
  PAX = X.A%*%par.a; PZX = PZA1X*PAX+PZA0X*(1-PAX)
  PAZX = (PZA1X*Z+(1-PZA1X)*(1-Z)) * (PAX)  / (PZX*Z+(1-PZX)*(1-Z)) #equal to predict(lm(A~Z*Age*Sex))
  fAZX = PAZX*A+(1-PAZX)*(1-A)
  fZAX = PZAX*Z+(1-PZAX)*(1-Z)
  
  ###construct EE for par.w using par.a, par.z
  Zc = Z-PZX #predict(lm(Z~Age*Sex)) ##Z.center E[Z|X]
  Ac = A-PAX #predict(lm(A~Age*Sex))##A.center E[A|X]
  AZc = A*Z-PZA1X*PAX#predict(lm(A*Z~Age*Sex))
  data = data.frame(A0,A1,A2,A,Z0,Z1,Z,Age,Sex,Ac,Zc,AZc)
  
  ###E[W|A=0,Z=0,X]: mle
  if(MW.wrong){
    if(omit.sex){
      # par.w_bl = lm(W~  Sex, data=data.frame(Age,Sex,W)[A==0&Z==0,])$coefficients
      X.W = model.matrix(W~Sex)
      AZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z","A*Z","W",T,T))
      AZ1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z1","A*Z1","W",T,T))
      AZ0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z0","A*Z0","W",T,T))
      A1ZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z","A1*Z","W",T,T))
      A0ZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z","A0*Z","W",T,T))
      A1Z1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z1","A1*Z1","W",T,T))
      A1Z0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z0","A1*Z0","W",T,T))
      A0Z1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z1","A0*Z1","W",T,T))
      A0Z0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z0","A0*Z0","W",T,T))
      AZX.W.center = model.matrix(data=data,object=get.formula(T,"Ac","Zc","AZc","W",T,T)) #each column should have mean zero
    }else{
      # par.w_bl = lm(W~  (Age+Sex), data=data.frame(Age,Sex,W)[A==0&Z==0,])$coefficients
      X.W = model.matrix(W~(Age+Sex))
      AZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z","A*Z","W",T,F))
      AZ1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z1","A*Z1","W",T,F))
      AZ0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z0","A*Z0","W",T,F))
      A1ZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z","A1*Z","W",T,F))
      A0ZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z","A0*Z","W",T,F))
      A1Z1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z1","A1*Z1","W",T,F))
      A1Z0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z0","A1*Z0","W",T,F))
      A0Z1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z1","A0*Z1","W",T,F))
      A0Z0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z0","A0*Z0","W",T,F))
      AZX.W.center = model.matrix(data=data,object=get.formula(T,"Ac","Zc","AZc","W",T,F)) #each column should have mean zero
    }
  }else{
    # par.w_bl = lm(W~  Age*Sex, data=data.frame(Age,Sex,W)[A==0&Z==0,])$coefficients
    X.W = model.matrix(W~Age*Sex)
    AZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z","A*Z","W",F,F))
    AZ1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z1","A*Z1","W",F,F))
    AZ0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z0","A*Z0","W",F,F))
    A1ZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z","A1*Z","W",F,F))
    A0ZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z","A0*Z","W",F,F))
    A1Z1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z1","A1*Z1","W",F,F))
    A1Z0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z0","A1*Z0","W",F,F))
    A0Z1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z1","A0*Z1","W",F,F))
    A0Z0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z0","A0*Z0","W",F,F))
    AZX.W.center = model.matrix(data=data,object=get.formula(T,"Ac","Zc","AZc","W",F,F)) #each column should have mean zero
  }
  par.w_bl=par[(ncol(AX.Z)+ncol(X.A)+1):(ncol(AX.Z)+ncol(X.A)+ncol(X.W))]
  EWA0Z0X = X.W%*%par.w_bl
  par.w= par[(ncol(AX.Z)+ncol(X.A)+ncol(X.W)+1):(ncol(AX.Z)+ncol(X.A)+ncol(X.W)+ncol(AZ0X.W.DeleteCol))]
  
  EWAZ0X = EWA0Z0X+AZ0X.W.DeleteCol%*%par.w
  D.ZW = (AZ1X.W.DeleteCol%*%par.w)-(AZ0X.W.DeleteCol%*%par.w) #mean(D.ZW);mean(eWAZ1X-eWAZ0X)
  ###Y
  if(MY.wrong){
    if(omit.sex){
      AZ0X.Y.DeletCol =  model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Y",T,T))
      A1Z0X.Y.DeletCol = model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Y",T,T))
      A0Z0X.Y.DeletCol = model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Y",T,T))
    }else{
      AZ0X.Y.DeletCol =  model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Y",T,F))
      A1Z0X.Y.DeletCol = model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Y",T,F))
      A0Z0X.Y.DeletCol = model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Y",T,F))
    }
  }else{
    AZ0X.Y.DeletCol =  model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Y",F,F))
    A1Z0X.Y.DeletCol = model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Y",F,F))
    A0Z0X.Y.DeletCol = model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Y",F,F))
  }
  par.y=par[(ncol(AX.Z)+ncol(X.A)+ncol(X.W)+ncol(AZ0X.W.DeleteCol)+1):(ncol(AX.Z)+ncol(X.A)+ncol(X.W)+ncol(AZ0X.W.DeleteCol)+ncol(AZ0X.Y.DeletCol))]
  EYAZ0X = AZ0X.Y.DeletCol%*%par.y
  
  if(MR.wrong){
    if(omit.sex){
      AX.R =  model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Y",T,T))
      A0X.R =  model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Y",T,T))
      A1X.R =  model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Y",T,T))
      A2X.R =  model.matrix(data=data,object=get.formula(F,"A2",NULL,NULL,"Y",T,T))
    }else{
      AX.R =  model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Y",T,F))
      A0X.R =  model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Y",T,F))
      A1X.R =  model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Y",T,F))
      A2X.R =  model.matrix(data=data,object=get.formula(F,"A2",NULL,NULL,"Y",T,F))
    }
  }else{
    AX.R =  model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Y",F,F))
    A0X.R =  model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Y",F,F))
    A1X.R =  model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Y",F,F))
    A2X.R =  model.matrix(data=data,object=get.formula(F,"A2",NULL,NULL,"Y",F,F))
  }
  
  par.r=par[(ncol(AX.Z)+ncol(X.A)+ncol(X.W)+ncol(AZ0X.W.DeleteCol)+ncol(AZ0X.Y.DeletCol)+1):(ncol(AX.Z)+ncol(X.A)+ncol(X.W)+ncol(AZ0X.W.DeleteCol)+ncol(AZ0X.Y.DeletCol)+ncol(AX.R))]
  ####delta
  EYAZX = AZ0X.Y.DeletCol%*%par.y + (AX.R%*%par.r) * D.ZW * Z
  D.ZWa1 = (A1Z1X.W.DeleteCol%*%par.w)-(A1Z0X.W.DeleteCol%*%par.w)
  D.ZWa0 = (A0Z1X.W.DeleteCol%*%par.w)-(A0Z0X.W.DeleteCol%*%par.w)
  EYA1ZX = A1Z0X.Y.DeletCol%*%par.y + (A1X.R%*%par.r) * D.ZWa1 * Z
  EYA0ZX = A0Z0X.Y.DeletCol%*%par.y + (A0X.R%*%par.r) * D.ZWa0 * Z
  D.AY = EYA1ZX - EYA0ZX #D.AY = (A1ZX.Y%*%par.y)-(A0ZX.Y%*%par.y) is a wrong parameterization: only Miao can use this
  D_D.AY = (2*A-1)/fAZX*(Y-EYAZX)
  
  D.AWz1 = (A1Z1X.W.DeleteCol%*%par.w)-(A0Z1X.W.DeleteCol%*%par.w)
  D.AWz0 = (A1Z0X.W.DeleteCol%*%par.w)-(A0Z0X.W.DeleteCol%*%par.w)
  fAX = PAX*A+(1-PAX)*(1-A); fA2X = PAX*A2+(1-PAX)*(1-A2)
  condE.D.AW = (D.AWz1*PZA2X+ D.AWz0*(1-PZA2X)) * fA2X/fAX
  fAX = PAX*A+(1-PAX)*(1-A); fA2X = PAX*A2+(1-PAX)*(1-A2)
  odds.AX = fA2X/fAX
  D_R2 = (2*Z-1)/fZAX *(Y-EYAZ0X-(AX.R%*%par.r)*(W-EWAZ0X)) /D.ZW *condE.D.AW*odds.AX
  
  R2_cond_ZX = (A0X.R%*%par.r*PAZX)+(A1X.R%*%par.r*(1-PAZX))
  EWAZX = EWA0Z0X+AZX.W.DeleteCol%*%par.w
  D_D.AW = R2_cond_ZX* (2*A-1)/fAZX *(W-EWAZX)
  
  D.AW = (A1ZX.W.DeleteCol%*%par.w)-(A0ZX.W.DeleteCol%*%par.w)
  
  delta=par[ncol(AX.Z)+ncol(X.A)+ncol(X.W)+ncol(AZ0X.W.DeleteCol)+ncol(AZ0X.Y.DeletCol)+ncol(AX.R)+1]
  delta.naive=par[ncol(AX.Z)+ncol(X.A)+ncol(X.W)+ncol(AZ0X.W.DeleteCol)+ncol(AZ0X.Y.DeletCol)+ncol(AX.R)+2]
  delta.bias=par[ncol(AX.Z)+ncol(X.A)+ncol(X.W)+ncol(AZ0X.W.DeleteCol)+ncol(AZ0X.Y.DeletCol)+ncol(AX.R)+3]
  
  U1 = c(Z-(AX.Z%*%par.z)) *(AX.Z)
  U2 = c(A-(X.A%*%par.a)) *(X.A)
  U3 = c(W-(X.W%*%par.w_bl)) *(X.W) *as.numeric(A==0&Z==0)
  U4 = c(W-(AZX.W.DeleteCol%*%par.w) - EWA0Z0X) *(AZX.W.center)
  ### note that we only take E[W|Z=0,A=0,X,beta^{W0}_{mle}]
  U5 = c(Y-(AZ0X.Y.DeletCol%*%par.y)) *(AZ0X.Y.DeletCol) *as.numeric(Z==0)
  ### note that we only take E[Y|Z=0,A,X,beta_{mle}]
  # U6 = c(Y-EYAZ0X-(AX.R%*%par.r)*(W-EWAZ0X) )*(AX.R) ##this is not DR EE
  U6=c(Y-EYAZX)*AX.R
  U7 = delta - (   (  D.AY + D_D.AY  ) - (  (A2X.R%*%par.r)*D.AW + D_R2+D_D.AW  )   )
  U8 = delta.naive - (  (  D.AY + D_D.AY  )   )
  U9 = delta.bias - (  (  (A2X.R%*%par.r)*D.AW + D_R2+D_D.AW  )  )
  U=cbind(U1,U2,U3,U4,U5,U6,U7,U8,U9)
  
  return(U)
}
Delta.mr = function(par,A,Z,Y,W,Age,Sex,MZ.wrong,MW.wrong,MR.wrong,MY.wrong){
  n=length(A)
  ### A, Z
  A1=rep(1,n);A0=rep(0,n);A2=1-A
  Z1=rep(1,n);Z0=rep(0,n);
  data = data.frame(A0,A1,A2,A,Z0,Z1,Z,Age,Sex)
  ### f(A,Z)
  if(MZ.wrong){
    if(omit.sex){
      AX.Z = model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Z",T,T))
      A0X.Z = model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Z",T,T))
      A1X.Z = model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Z",T,T))
      A2X.Z = model.matrix(data=data,object=get.formula(F,"A2",NULL,NULL,"Z",T,T))
      par.z = lm(get.formula(F,"A",NULL,NULL,"Z",T,T))$coefficients
    }else{
      AX.Z = model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Z",T,F))
      A0X.Z = model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Z",T,F))
      A1X.Z = model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Z",T,F))
      A2X.Z = model.matrix(data=data,object=get.formula(F,"A2",NULL,NULL,"Z",T,F))
      par.z = lm(get.formula(F,"A",NULL,NULL,"Z",T,F))$coefficients
    }
  }else{
    AX.Z = model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Z",F,F))
    A0X.Z = model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Z",F,F))
    A1X.Z = model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Z",F,F))
    A2X.Z = model.matrix(data=data,object=get.formula(F,"A2",NULL,NULL,"Z",F,F))
    par.z = lm(get.formula(F,"A",NULL,NULL,"Z",F,F))$coefficients
  }
  X.A = model.matrix(A~Age*Sex)
  par.a = lm(A~Age*Sex)$coefficients
  PZAX = (AX.Z%*%par.z); PZA1X = (A1X.Z%*%par.z); PZA0X = (A0X.Z%*%par.z); PZA2X = (A2X.Z%*%par.z);
  PAX = X.A%*%par.a; PZX = PZA1X*PAX+PZA0X*(1-PAX)
  PAZX = (PZA1X*Z+(1-PZA1X)*(1-Z)) * (PAX)  / (PZX*Z+(1-PZX)*(1-Z)) #equal to predict(lm(A~Z*Age*Sex))
  fAZX = PAZX*A+(1-PAZX)*(1-A)
  fZAX = PZAX*Z+(1-PZAX)*(1-Z)
  
  ###construct EE for par.w using par.a, par.z 
  Zc = Z-PZX #predict(lm(Z~Age*Sex)) ##Z.center E[Z|X]
  Ac = A-PAX #predict(lm(A~Age*Sex))##A.center E[A|X]
  AZc = A*Z-PZA1X*PAX#predict(lm(A*Z~Age*Sex))
  data = data.frame(A0,A1,A2,A,Z0,Z1,Z,Age,Sex,Ac,Zc,AZc)
  
  ###E[W|A=0,Z=0,X]: mle
  if(MW.wrong){
    if(omit.sex){
      par.w_bl = lm(W~Sex, data=data.frame(Age,Sex,W)[A==0&Z==0,])$coefficients
      X.W = model.matrix(W~Sex)
      AZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z","A*Z","W",T,T))
      AZ1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z1","A*Z1","W",T,T))
      AZ0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z0","A*Z0","W",T,T))
      A1ZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z","A1*Z","W",T,T))
      A0ZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z","A0*Z","W",T,T))
      A1Z1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z1","A1*Z1","W",T,T))
      A1Z0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z0","A1*Z0","W",T,T))
      A0Z1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z1","A0*Z1","W",T,T))
      A0Z0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z0","A0*Z0","W",T,T))
      AZX.W.center = model.matrix(data=data,object=get.formula(T,"Ac","Zc","AZc","W",T,T)) #each column should have mean zero
    }else{
      par.w_bl = lm(W~(Age+Sex), data=data.frame(Age,Sex,W)[A==0&Z==0,])$coefficients
      X.W = model.matrix(W~(Age+Sex))
      AZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z","A*Z","W",T,F))
      AZ1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z1","A*Z1","W",T,F))
      AZ0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z0","A*Z0","W",T,F))
      A1ZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z","A1*Z","W",T,F))
      A0ZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z","A0*Z","W",T,F))
      A1Z1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z1","A1*Z1","W",T,F))
      A1Z0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z0","A1*Z0","W",T,F))
      A0Z1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z1","A0*Z1","W",T,F))
      A0Z0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z0","A0*Z0","W",T,F))
      AZX.W.center = model.matrix(data=data,object=get.formula(T,"Ac","Zc","AZc","W",T,F)) #each column should have mean zero
    }
  }else{
    par.w_bl = lm(W~Age*Sex, data=data.frame(Age,Sex,W)[A==0&Z==0,])$coefficients
    X.W = model.matrix(W~Age*Sex)
    AZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z","A*Z","W",F,F))
    AZ1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z1","A*Z1","W",F,F))
    AZ0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A","Z0","A*Z0","W",F,F))
    A1ZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z","A1*Z","W",F,F))
    A0ZX.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z","A0*Z","W",F,F))
    A1Z1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z1","A1*Z1","W",F,F))
    A1Z0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A1","Z0","A1*Z0","W",F,F))
    A0Z1X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z1","A0*Z1","W",F,F))
    A0Z0X.W.DeleteCol = model.matrix(data=data,object=get.formula(T,"A0","Z0","A0*Z0","W",F,F))
    AZX.W.center = model.matrix(data=data,object=get.formula(T,"Ac","Zc","AZc","W",F,F)) #each column should have mean zero
  }
  EWA0Z0X = X.W%*%par.w_bl
  obj = function(par.w){
    U = c(W-(AZX.W.DeleteCol%*%par.w) - EWA0Z0X) *(AZX.W.center)
    return(sum(apply(U,2,sum)^2))
  }
  t=nlm(f=obj,p=rep(0,ncol(AZX.W.DeleteCol)),iterlim=1000,ndigit =20,gradtol=1e-15,stepmax=1000,steptol=1e-15)
  
  par.w= t$estimate
  # rbind(colnames(AZX.W.DeleteCol),par.w);lm(W~A*Z*Age*Sex)$coefficients
  # system("say done")
  
  EWAZ0X = EWA0Z0X+AZ0X.W.DeleteCol%*%par.w
  D.ZW = (AZ1X.W.DeleteCol%*%par.w)-(AZ0X.W.DeleteCol%*%par.w) #mean(D.ZW);mean(eWAZ1X-eWAZ0X)
  ###Y
  if(MY.wrong){
    if(omit.sex){
      AZ0X.Y.DeletCol =  model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Y",T,T))
      A1Z0X.Y.DeletCol = model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Y",T,T))
      A0Z0X.Y.DeletCol = model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Y",T,T))
    }else{
      AZ0X.Y.DeletCol =  model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Y",T,F))
      A1Z0X.Y.DeletCol = model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Y",T,F))
      A0Z0X.Y.DeletCol = model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Y",T,F))
    }
  }else{
    AZ0X.Y.DeletCol =  model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Y",F,F))
    A1Z0X.Y.DeletCol = model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Y",F,F))
    A0Z0X.Y.DeletCol = model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Y",F,F))
  }
  obj = function(par.y){
    U = c(Y-(AZ0X.Y.DeletCol%*%par.y)) *(AZ0X.Y.DeletCol) *as.numeric(Z==0)
    return(sum(apply(U,2,sum)^2))
  }
  t=nlm(f=obj,p=rep(0,ncol(AZ0X.Y.DeletCol)),iterlim=1000,ndigit =20,gradtol=1e-15,stepmax=1000,steptol=1e-15)
  par.y = t$estimate
  # rbind(rbind(colnames(AZ0X.Y.DeletCol),par.y)[,par.y!=0],
  #       lm(Y~  A*Age*Sex, data=data.frame(A,Age,Sex,Y)[Z==0,])$coefficients)
  EYAZ0X = AZ0X.Y.DeletCol%*%par.y
  
  if(MR.wrong){
    if(omit.sex){
      AX.R =  model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Y",T,T))
      A0X.R =  model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Y",T,T))
      A1X.R =  model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Y",T,T))
      A2X.R =  model.matrix(data=data,object=get.formula(F,"A2",NULL,NULL,"Y",T,T))
    }else{
      AX.R =  model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Y",T,F))
      A0X.R =  model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Y",T,F))
      A1X.R =  model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Y",T,F))
      A2X.R =  model.matrix(data=data,object=get.formula(F,"A2",NULL,NULL,"Y",T,F))
    }
  }else{
    AX.R =  model.matrix(data=data,object=get.formula(F,"A",NULL,NULL,"Y",F,F))
    A0X.R =  model.matrix(data=data,object=get.formula(F,"A0",NULL,NULL,"Y",F,F))
    A1X.R =  model.matrix(data=data,object=get.formula(F,"A1",NULL,NULL,"Y",F,F))
    A2X.R =  model.matrix(data=data,object=get.formula(F,"A2",NULL,NULL,"Y",F,F))
  }
  
  obj = function(par.r){
    EYAZX = EYAZ0X + (AX.R%*%par.r) * D.ZW * Z
    U=c(Y-EYAZX)*AX.R
    # U=c(Y-EYAZ0X-(AX.R%*%par.r)*(W-EWAZ0X) )*(AX.R)##this is gest, not DR
    return(sum(apply(U,2,sum)^2))
  }
  t=nlm(f=obj,p=rep(0,ncol(AX.R)),iterlim=1000,ndigit =20,gradtol=1e-15,stepmax=1000,steptol=1e-15)
  par.r= t$estimate
  # rrrr=(predict(lm(Y~A*Z*Age*Sex),newdata=data.frame(Y,A,Z=1,Age,Sex))-EYAZ0X)/
    # (predict(lm(W~A*Z*Age*Sex),newdata=data.frame(W,A,Z=1,Age,Sex))-
       # predict(lm(W~A*Z*Age*Sex),newdata=data.frame(W,A,Z=0,Age,Sex)))
  # rbind(colnames(AX.R),par.r);lm(rrrr~A*Age*Sex)$coef
  # system("say done")
  
  ####delta
  EYAZX = AZ0X.Y.DeletCol%*%par.y + (AX.R%*%par.r) * D.ZW * Z
  D.ZWa1 = (A1Z1X.W.DeleteCol%*%par.w)-(A1Z0X.W.DeleteCol%*%par.w)
  D.ZWa0 = (A0Z1X.W.DeleteCol%*%par.w)-(A0Z0X.W.DeleteCol%*%par.w)
  EYA1ZX = A1Z0X.Y.DeletCol%*%par.y + (A1X.R%*%par.r) * D.ZWa1 * Z
  EYA0ZX = A0Z0X.Y.DeletCol%*%par.y + (A0X.R%*%par.r) * D.ZWa0 * Z
  D.AY = EYA1ZX - EYA0ZX #D.AY = (A1ZX.Y%*%par.y)-(A0ZX.Y%*%par.y) is a wrong parameterization: only Miao can use this
  D_D.AY = (2*A-1)/fAZX*(Y-EYAZX)
  
  D.AWz1 = (A1Z1X.W.DeleteCol%*%par.w)-(A0Z1X.W.DeleteCol%*%par.w)
  D.AWz0 = (A1Z0X.W.DeleteCol%*%par.w)-(A0Z0X.W.DeleteCol%*%par.w)
  fAX = PAX*A+(1-PAX)*(1-A); fA2X = PAX*A2+(1-PAX)*(1-A2)
  condE.D.AW = (D.AWz1*PZA2X+ D.AWz0*(1-PZA2X)) * fA2X/fAX
  fAX = PAX*A+(1-PAX)*(1-A); fA2X = PAX*A2+(1-PAX)*(1-A2)
  odds.AX = fA2X/fAX
  D_R2 = (2*Z-1)/fZAX *(Y-EYAZ0X-(AX.R%*%par.r)*(W-EWAZ0X)) /D.ZW *condE.D.AW*odds.AX
  
  R2_cond_ZX = (A0X.R%*%par.r*PAZX)+(A1X.R%*%par.r*(1-PAZX))
  EWAZX = EWA0Z0X+AZX.W.DeleteCol%*%par.w
  D_D.AW = R2_cond_ZX* (2*A-1)/fAZX *(W-EWAZX) 
  
  D.AW = (A1ZX.W.DeleteCol%*%par.w)-(A0ZX.W.DeleteCol%*%par.w)
  delta = mean(
    (  D.AY + D_D.AY  ) - (  (A2X.R%*%par.r)*D.AW + D_R2+D_D.AW  )
  )
  delta.naive = mean(D.AY + D_D.AY)
  delta.bias = mean((A2X.R%*%par.r)*D.AW + D_R2+D_D.AW)
  
  par=c(par.z,par.a,par.w_bl,par.w,par.y,par.r,delta,delta.naive,delta.bias)
  return(par)
}

MR = function(aa,zz,ww,yy,Age,Sex,MZ.wrong,MW.wrong,MR.wrong,MY.wrong){
  est=Delta.mr(par=NA,A=aa,Z=zz,Y=yy,W=ww,Age=Age,Sex=Sex,MZ.wrong=MZ.wrong,MW.wrong=MW.wrong,MR.wrong=MR.wrong,MY.wrong=MY.wrong)
  # rbind(est[(length(est)-2):length(est)],c(delta.true,mean((eYA1ZX-eYA0ZX)),mean((eYA2Z1X-eYA2Z0X)/(eWA2Z1X-eWA2Z0X)*(eWA1ZX-eWA0ZX))))
  G.mr = function(par,A,Z,Y,W,Age,Sex,MZ.wrong,MW.wrong,MR.wrong,MY.wrong){return(apply(U.mr(par,A,Z,Y,W,Age,Sex,MZ.wrong,MW.wrong,MR.wrong,MY.wrong),2,sum))}
  bread=numDeriv::jacobian(func=G.mr,x=est,A=aa,Z=zz,Y=yy,W=ww,Age=Age,Sex=Sex,MZ.wrong=MZ.wrong,MW.wrong=MW.wrong,MR.wrong=MR.wrong,MY.wrong=MY.wrong)
  meat.half=U.mr(par=est,A=aa,Z=zz,Y=yy,W=ww,Age=Age,Sex=Sex,MZ.wrong=MZ.wrong,MW.wrong=MW.wrong,MR.wrong=MR.wrong,MY.wrong=MY.wrong)
  IF = meat.half%*%t(solve(-bread)) #ginv
  delta = est[length(est)-2]
  delta.var = sum(IF[,ncol(IF)-2]^2)
  delta.naive = est[length(est)-1]
  delta.naive.var = sum(IF[,ncol(IF)-1]^2)
  delta.bias = est[length(est)]
  delta.bias.var = sum(IF[,ncol(IF)]^2)
  return(
    list(
      ATE=as.numeric(delta),delta.naive=as.numeric(delta.naive),delta.bias=as.numeric(delta.bias),
      var=as.numeric(delta.var),delta.naive.var=as.numeric(delta.naive.var),delta.bias.var=as.numeric(delta.bias.var)
    )
  )
}


get.formula = function(remove.baseline.cov=T,names.A = "A",names.Z = "Z",names.AZ = "A*Z",names.outcome = "Y",misspec=F,omit.sex=F){
  #first order, second order, and third order terms
  FO.Age = "I(Age)";  
  FO.Sex = "I(Sex)"
  SO.AgeSex = paste0("I(Age*Sex)")
  ### NULL means we are fitting baseline models such as #[W|A=0,Z=0,X] or E[Y|A,Z=0,X]
  if(is.null(names.A)){
    FO.A = NULL   
    SO.AAge = NULL
    SO.ASex = NULL
    TO.AAgeSex = NULL
  }else{
    FO.A = paste0("I(",names.A,")")
    SO.AAge = paste0("I(",names.A,"*Age",")")
    SO.ASex = paste0("I(",names.A,"*Sex",")")
    TO.AAgeSex = paste0("I(",names.A,"*Age*Sex",")")
  }
  if(is.null(names.Z)){
    FO.Z = NULL
    SO.ZAge = NULL
    SO.ZSex = NULL
    TO.ZAgeSex = NULL
  }else{
    FO.Z = paste0("I(",names.Z,")")
    SO.ZAge = paste0("I(",names.Z,"*Age",")")
    SO.ZSex = paste0("I(",names.Z,"*Sex",")")
    TO.ZAgeSex = paste0("I(",names.Z,"*Age*Sex",")")
  }
  if(is.null(names.AZ)){
    SO.AZ = NULL
    TO.AZAge = NULL
    TO.AZSex = NULL
    FO.AZAgeSex = NULL
  }else{
    SO.AZ = paste0("I(",names.AZ,")")
    TO.AZAge = paste0("I(",names.AZ,"*Age",")")#third order
    TO.AZSex = paste0("I(",names.AZ,"*Sex",")")
    FO.AZAgeSex = paste0("I(",names.AZ,"*Age*Sex",")")#forth order
  }
  
  
  if(remove.baseline.cov==T){
    all.names = c(
      FO.A,      FO.Z,    
      SO.AZ,     SO.AAge, SO.ASex,  SO.ZAge,    SO.ZSex,
      TO.AZAge,  TO.AZSex,TO.AAgeSex,TO.ZAgeSex,
      FO.AZAgeSex
    )
  }else{
    all.names = c(
      FO.A,      FO.Z,    FO.Age,   FO.Sex,
      SO.AZ,     SO.AAge, SO.ASex,  SO.ZAge,    SO.ZSex, SO.AgeSex,
      TO.AZAge,  TO.AZSex,TO.AAgeSex,TO.ZAgeSex,
      FO.AZAgeSex
    )
  }
  
  
  if(misspec==T){
    if(omit.sex==T){
      formula.string = paste(all.names[-grep("Sex",all.names)], collapse="+")
    }else{
      formula.string = paste(all.names[-intersect(grep("Age",all.names),grep("Sex",all.names))], collapse="+")
    }
  }else{
    formula.string = paste(all.names, collapse="+")
  }
  if(remove.baseline.cov==T){
    final.formula = paste0(names.outcome,"~-1+",formula.string)
  }else{
    final.formula = paste0(names.outcome,"~",formula.string)
  }
  
  return(as.formula(final.formula))
}
expit= function(x){1/(1+exp(-x))}
logit = function(x){log(x/(1-x))}
