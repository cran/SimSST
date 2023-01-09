[![R-CMD-check](https://github.com/imstatsbee/SimSST/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/imstatsbee/SimSST/actions/workflows/R-CMD-check.yaml)

README
================
Mohsen Soltanifar
2023-JAN-04

# SimSST

The goal of SimSST is to simulate stop signal task data based on fixed
ssd method and the tracking method.

## Installation

You can install the development version of SimSST with:

``` r
library(gamlss.dist)
library(dplyr)
library(MASS)
library(SimSST)
```

## Example: Simulation with fixed ssd method

This function takes in nine variables and produces a matrix of stop
signal task data based on fixed ssd method

``` r
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
```

## Example: Simulation with tracking method

This function takes in nine variables and produces a matrix of stop
signal task data based on tracking method

``` r
mySSTdata2 <- 
  simsstrack(
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
mySSTdata2 
```

## Example: Simulating correlated SST data using general tracking method

This function takes in eleven variables and produces a matrix of stop
signal task data based on the generalized tracking method.

```{r}
mySSTdata3 <- simssgen(
     pid = c("FNLN1", "FNLN2", "FNLN2"),
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
```
