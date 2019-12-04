###########################APPLICATION CODE###########################
### Xu Shi (shixu@umich.edu)
### Below is the code for Section 5 of the manuscript
### Data is available from Dr. Jennifer C. Nelson (Jen.Nelson@kp.org)
### with the permission of the Vaccine Safety Datalink.
### Below we mimic the real patient data using simulated data
### Delta_1 is called GEST: biaesd if MR.wrong | MZ.wrong
### Delta_2 is called IPW: biaesd if MW.wrong | MZ.wrong
### Delta_3 is called OR: biaesd if MW.wrong | MR.wrong | MY.wrong
### MLE is called MIAO: biased if MW.wrong | MY.wrong
###########################APPLICATION CODE###########################
rm(list=ls())
library(MASS)
n=2000
set.seed(27)
expit = function(x)  exp(x)/(1+exp(x))
Age=rbinom(n,size=1,p=0.5)
Sex=rbinom(n,size=1,p=0.5)
A=rbinom(n,size=1,p=expit(-1+0.5*Age+0.5*Sex+0.5*Age*Sex));A2=1-A
Z=rbinom(n,size=1,p=expit(-1+0.5*A+0.5*Age+0.5*Sex+0.5*A*Age+0.5*A*Sex+0.5*Age*Sex+0.5*A*Age*Sex))
W=rbinom(n,size=1,p=expit(-1+0.5*A+0.5*Z+0.5*Age+0.5*Sex+0.5*A*Z+0.5*A*Age+0.5*Z*Age+0.5*A*Sex+0.5*Z*Sex+0.5*Age*Sex+0.5*A*Z*Age+0.5*A*Z*Sex+0.5*A*Age*Sex+0.5*Z*Age*Sex+0.5*A*Z*Age*Sex))
Y=rbinom(n,size=1,p=expit(-1+0.5*A+0.5*Z+0.5*Age+0.5*Sex+0.5*A*Z+0.5*A*Age+0.5*Z*Age+0.5*A*Sex+0.5*Z*Sex+0.5*Age*Sex+0.5*A*Z*Age+0.5*A*Z*Sex+0.5*A*Age*Sex+0.5*Z*Age*Sex+0.5*A*Z*Age*Sex))
aa=A;zz=Z;ww=W;yy=Y
MZ.wrong=F
MW.wrong=F
MR.wrong=T
MY.wrong=F
omit.sex=T
source("Functions_for_application.R")

(GEST.est = round(unlist(GEST(aa=A,zz=Z,ww=W,yy=Y,Age=Age,Sex=Sex,MZ.wrong,MR.wrong)),5))
### IPW
(IPW.est = round(unlist(IPW(aa=A,zz=Z,ww=W,yy=Y,Age=Age,Sex=Sex,MZ.wrong,MW.wrong)),5))
### OR
(OR.est = round(unlist(OR(aa=A,zz=Z,ww=W,yy=Y,Age=Age,Sex=Sex,MW.wrong,MR.wrong,MY.wrong)),5))
### MIAO
(MIAO.est = round(unlist(MIAO(aa=A,zz=Z,ww=W,yy=Y,Age=Age,Sex=Sex,MW.wrong,MY.wrong)),5))
### GEST
(MR.est = round(unlist(MR(aa=A,zz=Z,ww=W,yy=Y,Age=Age,Sex=Sex,MZ.wrong,MW.wrong,MR.wrong,MY.wrong)),5))

system("say done")
rbind(GEST.est,IPW.est,OR.est,MIAO.est,MR.est)
