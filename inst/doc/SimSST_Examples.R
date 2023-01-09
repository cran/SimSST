## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(dplyr)
library(gamlss.dist)
library(MASS)
library(SimSST)
set.seed(1)

## -----------------------------------------------------------------------------
mySSTdata1 <- 
  simssfixed(
    pid = c("FNLN1","FNLN1"),
    block = c(1,2), 
    n = c(10,10), 
    m = c(4,4), 
    SSD.b = c(220,240),
    dist.go = c("ExG","ExG"),
    theta.go = as.matrix(rbind(c(440,90,90),c(440,90,90))),
    dist.stop = c("ExG","ExG"),
    theta.stop = as.matrix(rbind(c(120,80,70),c(120,80,70)))
  )
mySSTdata1 

## -----------------------------------------------------------------------------
mySSTdata2 <- 
  simsstrack(
    pid = c("FNLN1","FNLN1"),
    block = c(1,2), 
    n = c(10,10), 
    m = c(4,4), 
    SSD.b = c(220,240), 
    dist.go = c("ExG","ExG" ),
    theta.go = as.matrix(rbind(c(440,90,90),c(440,90,90))),
    dist.stop = c("ExG","ExG" ),
    theta.stop = as.matrix(rbind(c(120,80,70),c(120,80,70)))
  )
mySSTdata2

## -----------------------------------------------------------------------------
Datatemp2 <- mySSTdata2 
ss_presented <- recode(Datatemp2[,3], 'Stop' = "1", 'Go' = "0")
inhibited <- Datatemp2[,4]
ssd <- Datatemp2[,8]
rt <- Datatemp2[,5]
srrt <- Datatemp2[,7]
Data2 <- cbind.data.frame(ss_presented, inhibited, ssd, rt, srrt)
for(i in 1:20) if(Data2$inhibited[i]==0) Data2$rt[i] <- Data2$srrt[i]  
myBEESTSdata2 <- (Data2[,-5])[order(ss_presented),]
myBEESTSdata2

## -----------------------------------------------------------------------------
mySSTdata3 <- simssgen(
     pid = c("FNLN1","FNLN2","FNLN2"),
     block = c(1,1,2),
     n = c(50,100,150),
     m = c(10,20,30),
     SSD.b = c(200,220,240),
     dist.go = c("ExG","ExG","ExG"),
     theta.go = as.matrix(rbind(c(400,60,30),c(440,90,90),c(440,90,90))),
     dist.stop = c("ExG","ExG","ExG"),
     theta.stop = as.matrix(rbind(c(100,70,60),c(120,80,70),c(120,80,70))),
     rho = c(0.35,0.45,0.45),
     d = c(50,65,75))
mySSTdata3

