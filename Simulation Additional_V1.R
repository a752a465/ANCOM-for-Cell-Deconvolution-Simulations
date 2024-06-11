#################################################################################################################
### Title: Additional Simulations V1
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
# fdr_cell = Specifies cell type(s) which are not differentially abundant

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
ANCOM.test <- function(ncell=ncell,data.proportions=data.proportions,group=group,N=N,fdr_cell=fdr_cell){
  power_cell <- which(lambda1 != lambda2)
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
## Multiple Testing ANCOM Test: Simple T-test equivalent -----
ANCOM.test.MT <- function(ncell=ncell,data.proportions=data.proportions,group=group,N=N,fdr_cell=fdr_cell){
  power_cell <- 1
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
Standard.test <- function(data.proportions=data.proportions,group=group,ncell=ncell,N=N,fdr_cell=fdr_cell){
  power_cell <- which(lambda1 != lambda2)
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
## Multiple Testing Current Standard Test: T-test with Correction for multiple testing and log-transformed proportions -----
Standard.test.MT <- function(data.proportions=data.proportions,group=group,ncell=ncell,N=N,fdr_cell=fdr_cell){
  power_cell <- 1
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
Poisson.test <- function(data.cells=data.cells,group=group,ncell=ncell,N=N,fdr_cell=fdr_cell){
  power_cell <- which(lambda1 != lambda2)
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
MWU.test <- function(data.proportions=data.proportions,group=group,ncell=ncell,N=N,fdr_cell=fdr_cell){
  power_cell <- which(lambda1 != lambda2)
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
Simulation <- function(Nsims=Nsims,sample.size=sample.size,lambda1=lambda1,lambda2=lambda2,W.star=W.star,ncell=ncell,fdr_cell=fdr_cell){
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
                             ANCOM.test(ncell=ncell,data.proportions=data.proportions,group=group,N=N,fdr_cell=fdr_cell),
                             Standard.test(ncell=ncell,data.proportions=data.proportions,group=group,N=N,fdr_cell=fdr_cell),
                             Poisson.test(data.cells=data.cells,group=group,ncell=ncell,N=N,fdr_cell=fdr_cell),
                             MWU.test(ncell=ncell,data.proportions=data.proportions,group=group,N=N,fdr_cell=fdr_cell))
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



## Simulation Function Family Wise Error Rate: Calculates the probability of at least 1 False Discovery -----
Simulation_FWER <- function(Nsims=Nsims,sample.size=sample.size,lambda1=lambda1,lambda2=lambda2,W.star=W.star,ncell=ncell,fdr_cell=fdr_cell){
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
                             ANCOM.test.MT(ncell=ncell,data.proportions=data.proportions,group=group,N=N,fdr_cell=fdr_cell),
                             Standard.test.MT(ncell=ncell,data.proportions=data.proportions,group=group,N=N,fdr_cell=fdr_cell))
    }
    # Calculating Mean Power and FDR across J Simualtions
    J_results <- J_results %>%
      group_by(Method)%>%
      filter(!is.na(Power))%>%
      mutate(FWER=case_when(FDR>=1~1,FDR==0~0,TRUE~NA))%>%
      mutate(Power=mean(Power,na.rm=TRUE),
             FDR=mean(FDR,na.rm=TRUE),
             FWER=mean(FWER,na.rm=TRUE))%>%
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


#################################################################################
# ADDITIONAL SIMULATION 1: False Discovery With a Large Number of Tests ----
#################################################################################
## Global Simulation Parameters -----
Nsims=200
ncell=150
sample.size=c(24,100,400)
W.star=120

### Simulation Parameters-----
set.seed(1200)
effect = 50
differential_celltypes <- c(1,2,3,4,5,51,52,53,54,55,101,102,103,104,105)
fdr_cell=which(!(seq(1,150,1)%in%differential_celltypes))
lambda1 <- c(rep(250,times=50),rep(500,times=50),rep(1000,times=50))
lambda2 <- lambda1
effects = rep(effect,15) * sample(c(1,1,-1),15,replace=TRUE)
lambda2[differential_celltypes] <- lambda2[differential_celltypes] + effects
### Simulation -----
Total_FWER <- Simulation_FWER(Nsims=Nsims,sample.size=sample.size,lambda1=lambda1,lambda2=lambda2,W.star=W.star,ncell=ncell,fdr_cell=fdr_cell)

## SAVING RESULTS -----
if(Outputs==TRUE){
  openxlsx::write.xlsx(data.frame(Total_FWER),paste0(Outputs_path,"Sim_Additional_FWER.xlsx"))
}

#################################################################################
# ADDITIONAL SIMULATION 2: ONE CELL TYPE IS MISSING IN ONE GROUP ----
#################################################################################
## Global Simulation Parameters -----
Nsims=400
ncell=12
sample.size=c(24,100,400)
W.star=9
fdr_cell=12
power_cell=1

## Generating Cell Types and Percentages -----
lambda2 = c(250,0  ,250,250,500,500,500,500,1000,1000,1000,1000)
sum_lambda2 = sum(lambda2)
sequence = seq(1,600,1)
totals = sequence+sum_lambda2
percentages = (sequence/totals)*100
options <- data.frame(sequence,totals,percentages)
### Sim: 0.1% -----
lambda1 = c(250,34,250,250,500,500,500,500,1000,1000,1000,1000) # Cell Type 2 comprises (34/6750+34)~0.1%
lambda2 = c(250,0  ,250,250,500,500,500,500,1000,1000,1000,1000) # Here, The Mean Count of Cell Type 2 is Zero
One_celltype_missing_01per <- Simulation(Nsims=Nsims,sample.size=sample.size,lambda1=lambda1,lambda2=lambda2,W.star=W.star,ncell=ncell,fdr_cell=fdr_cell)
One_celltype_missing_01per <- One_celltype_missing_01per %>% mutate(Missing_Cell_Percentage="0.1%")

### Sim: 2.5% -----
lambda1 = c(250,173,250,250,500,500,500,500,1000,1000,1000,1000) # Cell Type 2 comprises (173/6750+173)~2.5%
lambda2 = c(250,0  ,250,250,500,500,500,500,1000,1000,1000,1000) # Here, The Mean Count of Cell Type 2 is Zero
One_celltype_missing_25per <- Simulation(Nsims=Nsims,sample.size=sample.size,lambda1=lambda1,lambda2=lambda2,W.star=W.star,ncell=ncell,fdr_cell=fdr_cell)
One_celltype_missing_25per <- One_celltype_missing_25per %>% mutate(Missing_Cell_Percentage="2.5%")
### Sim: 5% -----
lambda1 = c(250,355,250,250,500,500,500,500,1000,1000,1000,1000) # Cell Type 2 comprises (355/6750+355)~5.0%
lambda2 = c(250,0  ,250,250,500,500,500,500,1000,1000,1000,1000) # Here, The Mean Count of Cell Type 2 is Zero
One_celltype_missing_5per <- Simulation(Nsims=Nsims,sample.size=sample.size,lambda1=lambda1,lambda2=lambda2,W.star=W.star,ncell=ncell,fdr_cell=fdr_cell)
One_celltype_missing_5per <- One_celltype_missing_5per %>% mutate(Missing_Cell_Percentage="5%")

### All Results ----
All_Results_onemissing <- bind_rows(One_celltype_missing_01per,One_celltype_missing_25per,One_celltype_missing_5per)%>%
  filter(!is.na(Sample_Size))

## SAVING RESULTS -----
if(Outputs==TRUE){
  openxlsx::write.xlsx(data.frame(All_Results_onemissing),paste0(Outputs_path,"Sim_Additional_OneMissing.xlsx"))
}

#################################################################################
# ADDITIONAL SIMULATION 3: POISSON AND MANN-WHITNEY U FOR SIMULATION 1 ----
#################################################################################

## Global Simulation Parameters -----
Nsims=400
ncell=15
sample.size=c(24,100,400)
W.star=9
fdr_cell=12

#######################################################
## Small Effect Size Simulations -----
#######################################################
### Low Abundance Cell Type -----
lambda1 = c(250,250,250,250,250,500,500,500,500,500,1000,1000,1000,1000,1000)
lambda2 = c(300,250,250,250,500,250,500,500,500,500,1000,1000,1000,1000,1000)
abund_low_effect_small <- Simulation(Nsims=Nsims,sample.size=sample.size,lambda1=lambda1,lambda2=lambda2,W.star=W.star,ncell=ncell,fdr_cell=fdr_cell)
abund_low_effect_small <- abund_low_effect_small %>% mutate(Abundance="Low Abund Cell")

### Moderate Abundance Cell Type ----
lambda1 = c(250,250,250,250,250,500,500,500,500,500,1000,1000,1000,1000,1000)
lambda2 = c(250,250,250,250,250,550,500,500,500,500,1000,1000,1000,1000,1000)
abund_mod_effect_small <- Simulation(Nsims=Nsims,sample.size=sample.size,lambda1=lambda1,lambda2=lambda2,W.star=W.star,ncell=ncell,fdr_cell=fdr_cell)
abund_mod_effect_small <- abund_mod_effect_small %>% mutate(Abundance="Moderate Abund Cell")

### High Abundance Cell Type ----
lambda1 = c(250,250,250,250,250,500,500,500,500,500,1000,1000,1000,1000,1000)
lambda2 = c(250,250,250,250,250,500,500,500,500,500,1050,1000,1000,1000,1000)
abund_high_effect_small <- Simulation(Nsims=Nsims,sample.size=sample.size,lambda1=lambda1,lambda2=lambda2,W.star=W.star,ncell=ncell,fdr_cell=fdr_cell)
abund_high_effect_small <- abund_high_effect_small %>% mutate(Abundance="High Abund Cell")

### All Results ----
Effect_small <- rbind(abund_low_effect_small,abund_mod_effect_small,abund_high_effect_small)%>%
  mutate(Effect="Small Effect")
rm(abund_low_effect_small,abund_mod_effect_small,abund_high_effect_small)

#######################################################
## Medium Effect Size Simulations -----
#######################################################
### Low Abundance Cell Type -----
lambda1 = c(250,250,250,250,250,500,500,500,500,500,1000,1000,1000,1000,1000)
lambda2 = c(350,250,250,250,250,500,500,500,500,500,1000,1000,1000,1000,1000)
abund_low_effect_medium <- Simulation(Nsims=Nsims,sample.size=sample.size,lambda1=lambda1,lambda2=lambda2,W.star=W.star,ncell=ncell,fdr_cell=fdr_cell)
abund_low_effect_medium <- abund_low_effect_medium %>% mutate(Abundance="Low Abund Cell")

### Moderate Abundance Cell Type ----
lambda1 = c(250,250,250,250,250,500,500,500,500,500,1000,1000,1000,1000,1000)
lambda2 = c(250,250,250,250,250,600,500,500,500,500,1000,1000,1000,1000,1000)
abund_mod_effect_medium <- Simulation(Nsims=Nsims,sample.size=sample.size,lambda1=lambda1,lambda2=lambda2,W.star=W.star,ncell=ncell,fdr_cell=fdr_cell)
abund_mod_effect_medium <- abund_mod_effect_medium %>% mutate(Abundance="Moderate Abund Cell")

### High Abundance Cell Type ----
lambda1 = c(250,250,250,250,250,500,500,500,500,500,1000,1000,1000,1000,1000)
lambda2 = c(250,250,250,250,250,500,500,500,500,500,1100,1000,1000,1000,1000)
abund_high_effect_medium <- Simulation(Nsims=Nsims,sample.size=sample.size,lambda1=lambda1,lambda2=lambda2,W.star=W.star,ncell=ncell,fdr_cell=fdr_cell)
abund_high_effect_medium <- abund_high_effect_medium %>% mutate(Abundance="High Abund Cell")

### All Results ----
Effect_medium <- rbind(abund_low_effect_medium,abund_mod_effect_medium,abund_high_effect_medium)%>%
  mutate(Effect="Medium Effect")
rm(abund_low_effect_medium,abund_mod_effect_medium,abund_high_effect_medium)

#######################################################
## Large Effect Size Simulations -----
#######################################################
### Low Abundance Cell Type -----
lambda1 = c(250,250,250,250,250,500,500,500,500,500,1000,1000,1000,1000,1000)
lambda2 = c(400,250,250,250,250,500,500,500,500,500,1000,1000,1000,1000,1000)
abund_low_effect_large <- Simulation(Nsims=Nsims,sample.size=sample.size,lambda1=lambda1,lambda2=lambda2,W.star=W.star,ncell=ncell,fdr_cell=fdr_cell)
abund_low_effect_large <- abund_low_effect_large %>% mutate(Abundance="Low Abund Cell")

### Moderate Abundance Cell Type ----
lambda1 = c(250,250,250,250,250,500,500,500,500,500,1000,1000,1000,1000,1000)
lambda2 = c(250,250,250,250,250,650,500,500,500,500,1000,1000,1000,1000,1000)
abund_mod_effect_large <- Simulation(Nsims=Nsims,sample.size=sample.size,lambda1=lambda1,lambda2=lambda2,W.star=W.star,ncell=ncell,fdr_cell=fdr_cell)
abund_mod_effect_large <- abund_mod_effect_large %>% mutate(Abundance="Moderate Abund Cell")
### High Abundance Cell Type ----
lambda1 = c(250,250,250,250,250,500,500,500,500,500,1000,1000,1000,1000,1000)
lambda2 = c(250,250,250,250,250,500,500,500,500,500,1150,1000,1000,1000,1000)
abund_high_effect_large <- Simulation(Nsims=Nsims,sample.size=sample.size,lambda1=lambda1,lambda2=lambda2,W.star=W.star,ncell=ncell,fdr_cell=fdr_cell)
abund_high_effect_large <- abund_high_effect_large %>% mutate(Abundance="High Abund Cell")
### All Results ----
Effect_large<- rbind(abund_low_effect_large,abund_mod_effect_large,abund_high_effect_large)%>%
  mutate(Effect="Large Effect")
rm(abund_low_effect_large,abund_mod_effect_large,abund_high_effect_large)

## FINAL RESULTS ----
Full_Results <- bind_rows(Effect_small,Effect_medium,Effect_large)%>%
  filter(!is.na(Method))
rm(Effect_small,Effect_medium,Effect_large)

# SAVING RESULTS -----
if(Outputs==TRUE){
  openxlsx::write.xlsx(data.frame(Full_Results),paste0(Outputs_path,"Sim_Additional_Sim1_Pois_MannWhitney.xlsx"))
}



