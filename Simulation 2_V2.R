#################################################################################################################
### Title: ANCOM Sim 2 v2
### Current Analyst: Alexander Alsup
### Last Updated: (06/07/2024) 
### Notes: 
#################################################################################################################

  # Note: The function below clears the global environment
rm(list=ls())
# PACKAGES -----
library(ggplot2)
library(compositions)
library(tidyr)
library(ggfortify)
library(openxlsx)
library(readxl)
library(dplyr)

# GLOBAL ARGUMENTS ----
# Controls whether to save simulation results locally as an excel document
Outputs=TRUE
# Path for output document
Outputs_path = "Example Path"

# READ ME ------
## Simulation Parameters -----
# data.cells = Simulated Cell Proportions
# group = A vector
# ncell = The number of cell types used in simulation
# Nsims = Number of simulations to run per configuration
# lambda1 = The vector of Poisson coefficients used to generate count data for group 1
# lambda2 = The vector of Poisson coefficients used to generate count data for group 2
# W.star = The critical value of the ANCOM Test statistic
# sample.size = A vector of sample sizes to use in simulation

## Example Cell Count Ranges and Simulation Values ------
# Overall WBC Count: (3.08-7.83 10^9/L): 

# NEUTROPHILS - *ASSESSED FOR POWER*
  # Control: 1.5-4.08 x 10^9/L, Simulation Value: 2.8 (*10^9/L)
  # Mild Neutropenia: 1.0–1.5 (*10^9/L),  Simulation Value: 1.25 (*10^9/L)
  # Moderate Neutropenia: 0.5–1.0 (*10^9/L), Simulation Value: 0.75 (*10^9/L)
  # Severe Neutropenia: 0.2–0.5 (*10^9/L),  Simulation Value: 0.25 (*10^9/L)

# LYMPHOCYTES : 1.29-3.4 (*10^9/L), Simulation Value: 2.4 (*10^9/L)

# MONOCYTES (*10^9/L): 0.14-0.74, Simulation Value: 0.44 (*10^9/L)

# EOSINOPHILS (*10^9/L): 0.04-0.59, Simulation Value: 0.32 (*10^9/L) - *ASSESSED FOR FDR*

# BASOPHILS (*10^9/L): 0.01-0.07, Simulation Value: 0.04 (*10^9/L)

## SOURCES: 
# Normal White Blood Cell Ranges (Omuse, G. et. al, https://doi.org/10.1371/journal.pone.0198444 )
# Neutropenia Values (Newburger & Dale,  https://doi.org/10.1053/j.seminhematol.2013.06.010 ) 


# Simulation: One differentially abundant cell type -----
# FUNCTIONS ----
## Data Generation: ----
## Generates Cell Counts, Cell Proportions, and Patient Group Based on inputs
data.generation <- function(ncell=ncell,lambda1=lambda1,lambda2=lambda2,N=N){
  # Vector of "Cell Type" Names for naming
  cell.type <- paste("Cell_",seq(from=1,to=ncell),sep="")
  # Generate Count data
  pois1 <- sapply(1:ncell, function(x) rpois(n=N/2, lambda=lambda1[x]))+1
  pois2 <- sapply(1:ncell, function(x) rpois(n=N/2, lambda=lambda2[x]))+1
  data.cells <- rbind(pois1,pois2); colnames(data.cells) <- cell.type
  # Convert to Proportion data
  data.proportions <- data.cells/rowSums(data.cells)
  return(list(cell.type,data.cells,data.proportions))
}
## ANCOM Test: Simple T-test equivalent -----
ANCOM.test <- function(ncell=ncell,data.proportions=data.proportions,group=group,N=N){
  power_cell <- which(lambda1 != lambda2)
  fdr_cell <- 4
  model_test <- NULL
  model_pvals <- NULL
  model_pvals_corrected <- NULL
  for(i in 1:ncell){
    y1 <- log(data.proportions)
    y = apply(y1, 2, function(x) x - y1[,i])[,-i]  
    x <- as.factor(group) 
    lm <-summary(lm(formula = y ~ x, data = data.frame(X = x)))
    model_pvals <- c(unlist(lapply(coef(lm), function(x)x[2,4])))
    model_pvals_corrected <- p.adjust(model_pvals,"BH")
    W <- sum(model_pvals_corrected <= 0.05)
    model_test[i] <- ifelse(W>= W.star,1,0)   
  } 
  model_power <- sum(model_test[power_cell])
  model_fdr <- sum(model_test[(fdr_cell)])
  result <- data.frame("Sample_Size"=N,"Method"="ANCOM","Power"=model_power,"FDR"=model_fdr)
  return(result)
}
## Current Standard Test: T-test with Correction for multiple testing and log-transformed proportions -----
Standard.test <- function(data.proportions=data.proportions,group=group,ncell=ncell,N=N){
  power_cell <- which(lambda1 != lambda2)
  fdr_cell <- 4
  model_pvals <- NULL
  x <- as.factor(group)
  y1 <- log(data.proportions)
  for(i in 1:ncell){
    # Trad model
    lm <- summary(lm(formula = y1[,i] ~ x, data = data.frame(X = x)))
    model_pvals[i] <- coef(lm)[2,4]
  }   
  #model_pvals_corrected <- p.adjust(model_pvals,"BH")
  model_test <- ifelse(model_pvals<= 0.05,1,0) 
  model_power <- sum(model_test[power_cell])
  model_fdr <- sum(model_test[(fdr_cell)])
  result <- data.frame("Sample_Size"=N,"Method"="Standard","Power"=model_power,"FDR"=model_fdr)
  return(result)
}
## Poisson Count Test: A test using Cell Counts and a Poisson GLM ----
Poisson.test <- function(data.cells=data.cells,group=group,ncell=ncell,N=N){
  power_cell <- which(lambda1 != lambda2)
  fdr_cell <- 4
  model_pvals <- NULL
  x <- as.factor(group)
  y1 <- data.cells
  for(i in 1:ncell){
    lm <- summary(glm(formula = data.cells[,i] ~ x, data = data.frame(X = x), family = poisson(link = "log")))
    model_pvals[i] <- coef(lm)[2,4]
  }
  #model_pvals_corrected <- p.adjust(model_pvals,"BH")
  model_test <- ifelse(model_pvals<= 0.05,1,0) 
  model_power <- sum(model_test[power_cell])
  model_fdr <- sum(model_test[(fdr_cell)])
  result <- data.frame("Sample_Size"=N,"Method"="Poisson on Counts","Power"=model_power,"FDR"=model_fdr)
  return(result)
}
## Mann-Whitney U Test: A non-parametric test on Cell Proportions -----
MWU.test <- function(data.proportions=data.proportions,group=group,ncell=ncell,N=N){
  power_cell <- which(lambda1 != lambda2)
  fdr_cell <- 4
  model_pvals <- NULL
  x <- as.factor(group)
  y1 <- data.proportions
  for(i in 1:ncell){
    lm <- wilcox.test(data.proportions[,i] ~ x, exact = FALSE)
    model_pvals[i] <-lm$p.value
  }
  #model_pvals_corrected <- p.adjust(model_pvals,"BH")
  model_test <- ifelse(model_pvals<= 0.05,1,0) 
  model_power <- sum(model_test[power_cell])
  model_fdr <- sum(model_test[(fdr_cell)])
  result <- data.frame("Sample_Size"=N,"Method"="Mann-Whitney U","Power"=model_power,"FDR"=model_fdr)
  return(result)
}
## Simulation Function: Adjust for tests to be included -----
Simulation <- function(Nsims=Nsims,sample.size=sample.size,lambda1=lambda1,lambda2=lambda2,W.star=W.star,ncell=ncell){
  time_0 <- Sys.time()
  K_compare <- data.frame("Sample_Size"=NA,"Method"=NA,"Power"=NA,"FDR"=NA)
  # K Level - Cycling through Sample Size
  for(k in 1:length(sample.size)){
    N <- sample.size[k]
    group <- c(rep("A", N/2), rep("B", N/2))
    J_results <- data.frame("Sample_Size"=NA,"Method"=NA,"Power"=NA,"FDR"=NA)
    # J Level - Simulations
    for(j in 1:Nsims){
      data <- data.generation(ncell=ncell,lambda1=lambda1,lambda2=lambda2,N=N)
      cell.type <- data[[1]];data.cells <- data[[2]];data.proportions <- data[[3]]
      J_results <- bind_rows(J_results,
                             ANCOM.test(ncell=ncell,data.proportions=data.proportions,group=group,N=N),
                             Standard.test(ncell=ncell,data.proportions=data.proportions,group=group,N=N))
    }
    # Calculating Mean Power and FDR across J Simualtions
    J_results <- J_results %>%
      group_by(Method)%>%
      filter(!is.na(Power))%>%
      mutate(Power=mean(Power,na.rm=TRUE),
             FDR=mean(FDR,na.rm=TRUE))%>%
      ungroup()%>% distinct()
    # Appending results to K-Level Dataframe
    K_compare <- bind_rows(K_compare,J_results)
    # Posting Progress
    cat(paste0(round(k/length(sample.size)*100,1),"%","||"))
  }
  time_1 <- Sys.time()
  cat(paste0("Simulation Complete. Duration ",round(time_1-time_0,2)," minutes"))
  return(K_compare)
}

# SIMULATIONS ----
## Global Simulation Parameters -----
Nsims=800
ncell=5
sample.size=c(24,50,74,100,200,400,600,800,1000)
W.star=3

## Simulation: Mild Neutropenia -----
lambda1 <- c(2.8,2.4,0.44,0.32,0.04)
lambda2 <- c(1.25,2.4,0.44,0.32,0.04)
Neutropenia_Mild <- Simulation(Nsims=Nsims,sample.size=sample.size,lambda1=lambda1,lambda2=lambda2,W.star=W.star,ncell=ncell)%>%
  mutate(Neutropenia="Mild")

## Simulation: Moderate Neutropenia -----
lambda1 <- c(2.8,2.4,0.44,0.32,0.04)
lambda2 <- c(0.75,2.4,0.44,0.32,0.04)
Neutropenia_Moderate <- Simulation(Nsims=Nsims,sample.size=sample.size,lambda1=lambda1,lambda2=lambda2,W.star=W.star,ncell=ncell)%>%
  mutate(Neutropenia="Moderate")

## Simulation: Severe Neutropenia -----
lambda1 <- c(2.8,2.4,0.44,0.32,0.04)
lambda2 <- c(0.25,2.4,0.44,0.32,0.04)
Neutropenia_Severe <- Simulation(Nsims=Nsims,sample.size=sample.size,lambda1=lambda1,lambda2=lambda2,W.star=W.star,ncell=ncell)%>%
  mutate(Neutropenia="Severe")

## FINAL RESULTS ----
Full_Results <- bind_rows(Neutropenia_Mild,Neutropenia_Moderate,Neutropenia_Severe)%>%
  filter(!is.na(Method))
rm(Neutropenia_Mild,Neutropenia_Moderate,Neutropenia_Severe)


# SAVING RESULTS -----
if(Outputs==TRUE){
  openxlsx::write.xlsx(data.frame(Full_Results),paste0(Outputs_path,"Sim2_Results.xlsx"))
}
