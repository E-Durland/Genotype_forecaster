###################################################################################
#This R script will walk you through a simple example with the genotype forecaster
###################################################################################

## Make sure you have ggplot2 and reshape2 before proceeding:
library(ggplot2)
library(reshape2)

## either load the core script or source to its path:
source("path/to/Gtype_forecaster.R")

## the script needs 3 vectors to run:
## 1) the days you want to estimate (sampling days)
## 2) The survival in your population on those days
## 3) Corresponding measures of minor allele frequency

## first, set up the days that were sampled:
s_days<-c(0,2,6,10,16,22)

## next, put in the observed survival rate for each sampling point:
surv<-c(1,0.64,0.3,0.25,0.19,0.04)

## The previous two vectors should be fixed for a single culture/experiment
## finally, input allele frequencies for a single locus:
freq<-c(0.50,	0.3,	0.4,	0.6,	06.,	0.6)

## the function will simulate one plausible scenario for changes in genotype frequency based on previous parameters
## with three inputs: 
## 1) the MAF vector ('freq' from above)
## 2) 'n_Gtyp' which is a numeric input (2-3) that accounts for the number of genotypes for this locus
## 3) Plotting (Y/N)

## One example run:
out<-Gtype_simulator(freq,3,10,"Y")

## you can find the raw data in the output:
head(out)

## From there you can try new frequencies or survival curves etc. 
## This function can easily be used in a loop to run on an entire dataset

#######################################################################################
## **NOTE** remember that this function essentially tells you what is possible given your inputs.
## It performs best (highest accuracy and precision) when mortality is high and changes big *and dynamic)
## It will also keep trying 
