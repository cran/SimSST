#
# ./R/Updated-SST-Simulation-TrackingMethod.R
#
# Chel Hee Lee & Mohsen Soltanifar
# 2022-DEC-03
#

# @rdname simsstrack
simsstrack0 <- function(pid, n, m, SSD.b,
                                            dist.go, theta.go,
                                            dist.stop, theta.stop)
{

  SRRT0 <- -999
  index <- as.vector(matrix(0,nrow=1,ncol = n))
  id <- as.vector(matrix(pid,nrow=1,ncol = n))

  if(!(dist.go %in% c("ExG", "SW"))) stop("This is incorrect Go distribution!")
  if(!(dist.stop %in% c("ExG", "SW"))) stop("This is incorrect Stop distribution!")

  if(dist.go=="ExG" & dist.stop=="ExG") {
    GORT  <- round(gamlss.dist::rexGAUS(n, mu = theta.go[1], sigma = theta.go[2] , nu = theta.go[3]), digits = 1)
    SSRT  <- round(gamlss.dist::rexGAUS(n, mu = theta.stop[1], sigma = theta.stop[2], nu = theta.stop[3]), digits = 1)
  }

  if(dist.go == "ExG" & dist.stop == "SW") {
    GORT  <- round(gamlss.dist::rexGAUS(n, mu = theta.go[1], sigma = theta.go[2] , nu = theta.go[3]), digits = 1)
    SSRT  <- round(gamlss.dist::rIG(n, mu = theta.stop[1], sigma = theta.stop[2])+ theta.stop[3], digits = 1)
  }

  if(dist.go == "SW" & dist.stop == "ExG") {
    GORT  <- round(gamlss.dist::rIG(n, mu = theta.go[1], sigma = theta.go[2])+ theta.go[3], digits = 1)
    SSRT  <- round(gamlss.dist::rexGAUS(n, mu = theta.stop[1], sigma = theta.stop[2], nu = theta.stop[3]), digits = 1)
  }

  if(dist.go == "SW" & dist.stop == "SW") {
    GORT  <- round(gamlss.dist::rIG(n, mu = theta.go[1], sigma = theta.go[2])+ theta.go[3], digits = 1)
    SSRT  <- round(gamlss.dist::rIG(n, mu = theta.stop[1], sigma = theta.stop[2])+ theta.stop[3], digits = 1)
  }

  SSD <- as.vector(matrix(0 ,nrow=1,ncol = n))
  SSD[1] <- SSD.b
  SRRT <- as.vector(matrix(SRRT0,nrow=1,ncol=n))

  TrialType <- as.vector(matrix('Go',nrow=1,ncol = n))
  Inhibition <- as.vector(matrix('-999',nrow=1,ncol = n))

  ##
  ## First: Check the fist stop trial
  ##

  if (GORT[1] > SSRT[1]+SSD[1]) TrialType[1]<- 'Stop(Successful)'
  else TrialType[1]<- 'Stop(Failed)'

  if (GORT[1] > SSRT[1]+SSD[1]) Inhibition[1]<- '1'
  else Inhibition[1]<- '0'

  if (GORT[1] > SSRT[1]+SSD[1]) index[1]<-  +50
  else index[1]<- -50

  if (GORT[1] > SSRT[1]+SSD[1]) SSD[2]<-  SSD[1]+50
  else SSD[2]<- SSD[1]-50

  ##
  ## Second: Use information above to create a martingale of SSDs !
  ##

  for (i in 1:(m-1)) {
    if (index[i]== 50) SSD[i+1]<- SSD[i]+50
    if (index[i]== 50)  {
      if (GORT[i+1] > SSRT[i+1]+SSD[i+1]) TrialType[i+1]<- 'Stop(Successful)'
      else TrialType[i+1]<- 'Stop(Failed)'
    }
    if (index[i]== 50)  {
      if (GORT[i+1] > SSRT[i+1]+SSD[i+1]) index[i+1]<- +50
      else index[i+1]<- -50
    }

    if (index[i]== -50) SSD[i+1]<- SSD[i]-50
    if (index[i]== -50) {
      if (GORT[i+1] > SSRT[i+1]+SSD[i+1]) TrialType[i+1]<- 'Stop(Successful)'
      else TrialType[i+1]<- 'Stop(Failed)'
    }
    if (index[i]== -50) {
      if (GORT[i+1] > SSRT[i+1]+SSD[i+1]) index[i+1]<- +50
      else index[i+1]<- -50
    }

    if (index[i]== 50)  {
      if (GORT[i+1] > SSRT[i+1]+SSD[i+1]) Inhibition[i+1]<- '1'
      else Inhibition[i+1]<- '0'
    }

    if (index[i]== -50) {
      if (GORT[i+1] > SSRT[i+1]+SSD[i+1]) Inhibition[i+1]<- '1'
      else Inhibition[i+1]<- '0'
    }
  }

  ##
  ## Third: Update Final Values of SSD given values of Index !
  ##

  SSRTSSD <- as.vector(matrix(0,nrow=1,ncol =n))

  for (i in 1:n) SSRTSSD[i] <- SSRT[i] + SSD[i]

  ##
  ## Fourth: We clear unwanted data for Go trials !
  ##

  for (i in (m+1):n) {
    if (SSRT[i]!= 'NA') SSRT[i]= -999
    if (SSRTSSD[i]!= 'NA') SSRTSSD[i]= -999
    if (SSD[i] != 'NA') SSD[i]= -999
    if (index[i] != 'NA') index[i]= -999
  }

  SRRT <- ifelse(TrialType %in% c('Stop(Failed)'), GORT, SRRT0)
  GORT <- ifelse(TrialType %in% c('Stop(Successful)','Stop(Failed)'), -999, GORT)

  TrialTypee <- dplyr::recode(TrialType,
                              'Stop(Successful)' = "Stop",
                              'Stop(Failed)' = "Stop",
                              'Go' = "Go")
  MAT2 <- matrix(c(id , TrialTypee, Inhibition , GORT , SSRT ,SRRT , SSD ), nrow=length(GORT ))

  MAT22 <- matrix(NA,ncol = 7, nrow=length(GORT ))
  for (i in 1:m){
    MAT22[2*i-1,] <- MAT2[i,];
    MAT22[2*i,] <- MAT2[m+i,]
  }
  for(i in (2*m+1):n) MAT22[i,] <- MAT2[i,]

  colnames(MAT22) <- c( 'Participant.id', 'Trial','Inhibition', 'GORT', 'SSRT', 'SRRT','SSD');
  return(MAT22)
}


#' @rdname simsstrack
#' @title Simulating SSRT data using tracking method
#' @description Stop signal task data of go and stop trials is generated per participant. The tracking signal delay method with underlying exponentially modified Gaussian(ExG) or Shifted Wald(SW) distributions for each of go and stop process is applied. The output data can be converted to 'BEESTS' software input data enabling researchers to test and evaluate different distributional parameters of interest.
#' @param pid a character vector of size `b` of participant
#' @param block a numeric vector of size `b` blocks
#' @param SSD.b a numeric vector of size `b` of starting stop signal delay
#' @param n a numeric vector of size `b` of total number of trials
#' @param m a numeric vector of size `b` of total number of stops
#' @param dist.go a character vector of size `b` of distribution of go trials, either ExG or SW
#' @param dist.stop a character vector of size `b` of distribution of stop.trials, either ExG or SW
#' @param theta.go a numeric matrix of size `b` by columns `mu.go`, `sigma.go`, `tau.go`
#' @param theta.stop a numeric matrix of size `b` by columns `mu.stop`, `sigma.stop`, `tau.stop`
#' @returns a matrix with `sum(n)` rows and (8) columns
#'
#' @references
#' Gordon D. Logan. On the Ability to Inhibit Thought and Action: A User's Guide to the Stop Signal Paradigm. In D. Dagenbach, & T.H. Carr (Eds.), Inhibitory Process in Attention, Memory and Language. San Diego: Academic Press, 1994.
#'
#' Dora Matzke, Jonathon Love, Thomas V. Wiecki, Scott D. Brown, and et al. Release the BEESTS: Bayesian Estimation of Ex-Gaussian Stop Signal Reaction Times Distributions. Frontiers in Psychology, 4: Article 918, 2013.
#'
#' Mohsen Soltanifar. Stop Signal Reaction Times: New Estimations with Longitudinal, Bayesian and Time  Series based Methods, PhD Dissertation, Biostatistics Division, Dalla Lana School of Public Health, University of Toronto, Toronto, Canada, 2020.
#'
#' @examples
#' mySSTdata1 <- simsstrack(
#'     pid=c("John.Smith","Jane.McDonald","Jane.McDonald"), block=c(1,1,2),
#'     n=c(50,100,150), m=c(10,20,30),
#'     SSD.b=c(200,220,240), dist.go=c("ExG","ExG","ExG"),
#'     theta.go=as.matrix.data.frame(rbind(c(400,60,30),c(440,90,90),c(440,90,90))),
#'     dist.stop=c("ExG","ExG","ExG"),
#'     theta.stop=as.matrix.data.frame(rbind(c(100,70,60),c(120,80,70),c(120,80,70))))
#' mySSTdata1
#'
#' @export


# SstSimulatedTrackingMethod
simsstrack <- function(pid, block, n, m, SSD.b,
                                       dist.go, theta.go,
                                       dist.stop, theta.stop
                                       ) {

  b <- length(block)
  csn <- c(0, cumsum(n))
  M2 <- matrix(NA, nrow = sum(n), ncol = 8)

  for(i in 1:b){
    M2[c((csn[i]+1):csn[i+1]),1] <- block[i]
    M2[c((csn[i]+1):csn[i+1]),c(2:8)] <- simsstrack0(pid=pid[i], n=n[i], m=m[i], SSD.b=SSD.b[i], dist.go=dist.go[i], theta.go=theta.go[i,], dist.stop=dist.stop[i], theta.stop=theta.stop[i,])
  }

  M2[,c(1,2)] <- M2[,c(2,1)]
  M22 <- M2
  colnames(M22) <- c('Participant.id','Block','Trial','Inhibition','GORT','SSRT','SRRT','SSD')

  return(M22)
}

