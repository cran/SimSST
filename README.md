[![R-CMD-check](https://github.com/imstatsbee/SimSST/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/imstatsbee/SimSST/actions/workflows/R-CMD-check.yaml)

README
================
Mohsen Soltanifar
2022-DEC-05

<!-- README.md is generated from README.Rmd. Please edit that file -->

# SimSST

The goal of SimSST is to simulate stop signal task data based on fixed
ssd method and the tracking method.

## Installation

You can install the development version of SimSST with:

``` r
library(gamlss.dist)
library(dplyr)
# install.packages("SimSST")
```

## Example: Simulation with fixed ssd method

This function takes in nine variables and produces a matrix of stop
signal task data based on fixed ssd method

``` r
library(SimSST)

mySSTdata1 <- 
  simssfixed(
    pid=c("FNLN1","FNLN1"), 
    block = c(1,2),
    n=c(10,10), m=c(4,4), SSD.b=c(220,240),
    dist.go=c("ExG","ExG"),
    theta.go=as.matrix(rbind(c(440,90,90),c(440,90,90))),
    dist.stop=c("ExG","ExG"),
    theta.stop=as.matrix(rbind(c(120,80,70),c(120,80,70))))
mySSTdata1 
```

## Example: Simulation with tracking method

This function takes in nine variables and produces a matrix of stop
signal task data based on tracking method

``` r
library(SimSST)
mySSTdata1 <- 
  simsstrack(
    pid=c("FNLN1","FNLN1"), 
    block = c(1,2),
    n=c(10,10), m=c(4,4), SSD.b=c(220,240),
    dist.go=c("ExG","ExG"),
    theta.go=as.matrix(rbind(c(440,90,90),c(440,90,90))),
    dist.stop=c("ExG","ExG"),
    theta.stop=as.matrix(rbind(c(120,80,70),c(120,80,70))))
mySSTdata1 
```
