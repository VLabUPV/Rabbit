## ---------------------------
##
## Script name: runRabbit
##
## Purpose of script: User-friendly software to fit mixed linear models using Bayesian Inference
##
## Authors: Marina Martinez Alvaro and Cristina Casto Rebollo
##
## Date Created: 2023-02-12
##
## Email: mamaral9@upv.es; ciscasre@posgrado.upv.es
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

if ((!"pacman" %in% installed.packages())==TRUE){install.packages("pacman")}
pacman::p_load(brms,emmeans,tidybayes,tidyr,magrittr,HDInterval,crayon,ggplot2,gridExtra,gtable,grid,ggpubr,cowplot)

# library(brms)
# library(emmeans)
# library(dplyr)
# library(tidybayes)
# library(tidyr)
# library(magrittr)
# library(HDInterval)
# library(crayon)
# library(readxl)


#runRabbit <- function() {
  rm(list = ls())
## Read the data from file name ----------------------------------------------------------------
  file.name <- readline(sprintf("%s\n",green("Enter the name of the datafile with its extension (.csv or .xlsx) ")))

  while(!file.exists(file.name)==TRUE) {
    print(paste("DATA FILE",file.name,"NO FOUND",sep=" "))
    file.name <- readline(sprintf("%s\n",green("Enter the name of the datafile with its extension (.csv or .xlsx) ")))
  }
  
  Missing     <- readline(sprintf("%s\n",green("Has the data file missing values (Enter Yes=Y or No=N) ?  ")))
  while((Missing!= "y") && (Missing!= "n") && (Missing!= "Y") && (Missing!= "N") ==TRUE){ 
      print(paste("Weird response, try again!"))
      Missing     <- readline(sprintf("%s\n",green("Has the data file missing values (Enter Yes=Y or No=N) ?  ")))} 
        
  if(Missing == "Y"| Missing == "y"){
  MissingName <- readline(sprintf("%s\n",green("Please enter the missing value (if its a blank enter a space) ")))
  }
  
  if(sub(".*\\.", "", file.name) =="csv"){
    if (Missing == "Y"| Missing == "y"){data <- read.csv(file.name,header=T)} else {data <- read.csv(file.name,header=T, na.strings=MissingName)}}
  if(sub(".*\\.", "", file.name) =="xlsx"){
    if (Missing == "Y"| Missing == "y"){data <- data.frame(read_excel(file.name,col_names=T))} else {data <- data.frame(read_excel(file.name,col_names=T, na=MissingName))}}  

  ri=nrow(data)
  cat(green(paste0("The number of rows in the data file is ", ri)))
  
  
  
## Interactive parameter file --------------------------------------------------------------------
  nTrait <- as.numeric(readline(sprintf("%s\n",green("Enter the total number of traits  "))))
  hTrait<-NULL
  pTrait<-NULL
  cat(green(paste0(c("Help: the header of the datafile is ", colnames(data)),collapse=" ")))
  if (nTrait>1){
    TraitsAsString<-readline((sprintf("%s\n",green(paste("Do you want to enter the name of all Traits at once (Enter Yes=Y or No=N) ? ")))))
  }
  
  if (TraitsAsString == "Y"| TraitsAsString == "y"){
    cat(green("Enter the positions of the Traits in the datafile"))
    pTraitInString<- readline(sprintf("%s\n",green(paste("For example, enter c(1:3) for 1,2,3; c(1,4) for 1 and 4; or c(1:3,7) for 1,2,3 and 7 "))))
    pTrait<-eval(parse(text = pTraitInString))
    hTrait <- colnames(data)[pTrait]
    cat(green(paste0(c("Traits read are", hTrait),collapse=" ")))
  }else{
    for(n in 1:nTrait){
        hTrait[n]<- readline((sprintf("%s\n",green(paste("Enter the name of the Trait ", n," ")))))
          while (!hTrait[n]%in%colnames(data)==TRUE) {
            cat(green(paste0("Trait not found. Please check that spelling is correct")))
            hTrait[n]<- readline((sprintf("%s\n",green(paste("Enter the name of the Trait ", n," ")))))
          }
        pTrait[n] <- which(colnames(data)== hTrait[n])
     }
  }

  cat(green(paste0("Let's define the model. Remember, all traits will be analyzed with the same model!")))
  
  nTreatment  <- as.numeric(readline(sprintf("%s\n",green("Enter the number of treatments   "))))
  if(nTreatment != 0){
      hTreatment<-NULL
      pTreatment<-NULL
      nlevels_Treatment<-NULL
      
      cat(green(paste0(c("Help: the header of the datafile is ", colnames(data)),collapse=" ")))
      for(n in 1:nTreatment){
        hTreatment[n] <- readline(sprintf("%s\n",green(paste("Enter the name of Treatment ", n, " "))))
        while (!hTreatment[n]%in%colnames(data)==TRUE) {
          cat(green(paste0("Treatment not found. Please check that spelling is correct")))
          hTreatment[n] <- readline(sprintf("%s\n",green(paste("Enter the name of Treatment ", n, " "))))
        }
        pTreatment[n] <- which(colnames(data)== hTreatment[n])
        data[,pTreatment[n]]<-as.factor(data[,pTreatment[n]])
        nlevels_Treatment[n] <-nlevels(data[,pTreatment[n]])
        cat(green(paste0(c("The number of levels read in Treatment",hTreatment[n],"are: ", nlevels_Treatment[n]),collapse=" ")))
      }
      
      askCompare <- readline("Comparisons between treatment levels:   Enter DIFFERENCE=D or RATIO=R   ")
      while((askCompare!= "d") && (askCompare!= "D") && (askCompare!= "R") && (askCompare!= "r") ==TRUE){ 
        print(paste("Weird response, try again!"))
        askCompare <- readline("Comparisons between treatment levels:   Enter DIFFERENCE=D or RATIO=R   ")}
      
  }
  
  nNoise      <- as.numeric(readline(sprintf("%s\n",green("Enter the number of noise effects   "))))
  if(nNoise != 0){
    hNoise<-NULL
    pNoise<-NULL
    nlevels_Noise<-NULL
    cat(green(paste0(c("Help: the header of the datafile is ", colnames(data)),collapse=" ")))
    for(n in 1:nNoise){
      hNoise[n] <-readline(sprintf("%s\n",green(paste("Enter the name of the Noise effect", n, " "))))
      while (!hNoise[n]%in%colnames(data)==TRUE) {
        cat(green("Noise effect not found. Please check that spelling is correct"))
        hNoise[n] <-readline(sprintf("%s\n",green(paste("Enter the name of the Noise effect", n, " "))))
      }
        pNoise[n] <- which(colnames(data)== hNoise[n])
        data[,pNoise[n]]<-as.factor(data[,pNoise[n]])
        nlevels_Noise[n] <-nlevels(data[,pNoise[n]])
        cat(green(c("The number of levels read in Noise",hNoise[n],"are: ", nlevels_Noise[n]), collapse=" "))
      }
        
  }
  
  
  #askInterFix <- readline(sprintf("%s\n",green(paste("Do you want to consider any interactions of order 2 between fixed effects (Enter Yes=Y or No=N)  ? "))))
  #while((askInterFix!= "y") && (askInterFix!= "n") && (askInterFix!= "Y") && (askInterFix!= "N") ==TRUE){ 
  #  print(paste("Weird response, try again!"))
  #  askInterFix <- readline(sprintf("%s\n",green(paste("Do you want to consider any interactions of order 2 between fixed effects (Enter Yes=Y or No=N)  ? "))))} 
  askInterFix <- "n"
   if ((askInterFix=="Y") | (askInterFix == "y")){
    nInter <- as.numeric(readline(sprintf("%s\n",green("How many interactions do you want to consider? "))))
    hInter <-matrix(ncol=nInter, nrow=2)
    cInter <-matrix(ncol=nInter, nrow=2)
    nlevels_Interaction<-NULL
    for (n in 1:nInter){
      cat(green(paste0("Let's define Interaction ", n)))
          hInter[1,n] <-readline(sprintf("%s\n",green(paste("Enter the name of the first effect (noise or treatment) to be considered in the interaction", n, " "))))
              while (!hInter[1,n]%in%colnames(data)==TRUE) {
                cat(green(paste0("Effect not found. Please check that spelling is correct")))
                hInter[1,n] <-readline(sprintf("%s\n",green(paste("Enter the name of the first effect (noise or treatment) to be considered in the interaction", n, " "))))
              } 
              cInter[1,n] <-which(colnames(data)== hInter[1,n])
          
              hInter[2,n] <-readline(sprintf("%s\n",green(paste("Enter the name of the second effect (noise or treatment) to be considered in the interaction", n, " "))))
              
              while (!hInter[2,n]%in%colnames(data)==TRUE) {
                cat(green(paste0("Effect not found. Please check that spelling is correct")))
                hInter[2,n] <-readline(sprintf("%s\n",green(paste("Enter the name of the second effect (noise or treatment) to be considered in the interaction", n, " "))))
              }  
          
              cInter[2,n] <-which(colnames(data)== hInter[2,n])
          nlevels_Interaction[n]<-nlevels(as.factor(data[,cInter[1,n]]))*nlevels(as.factor(data[,cInter[2,n]]))
          cat(green(paste0("The number of levels for Interaction ", n, " is: ", nlevels_Interaction[n])))
  }}else{nInter <-0}
 
  
  nCov <- as.numeric(readline(sprintf("%s\n",green("Enter the number of covariates (0,1,2,...)   "))))
  if(nCov != 0){
    hCov<-NULL
    pCov<-NULL
    cat(green(paste0(c("Help: the header of the datafile is ", colnames(data)),collapse=" ")))
    for (n in 1:nCov){
    hCov[n] <-  readline(sprintf("%s\n",green(paste("Enter the name of Covariate ", n, " "))))
    while (!hCov[n]%in%colnames(data)==TRUE) {
      cat(green(paste0("Covariate not found. Please check that spelling is correct")))
      hCov[n] <-  readline(sprintf("%s\n",green(paste("Enter the name of Covariate ", n, " "))))
    }
      pCov[n] <- which(colnames(data)== hCov[n]) 
  }}

  nRand <- as.numeric(readline(sprintf("%s\n",green("Enter the number of random effects   "))))
  if(nRand != 0){
    hRand<-NULL
    pRand<-NULL
    cat(green(paste0(c("Help: the header of the datafile is ", colnames(data)),collapse=" ")))
    for (n in 1:nRand){
      hRand[n] <- readline(sprintf("%s\n",green(paste("Enter the name of the random effect ", n, " "))))
      while (!hRand[n]%in%colnames(data)==TRUE) {
        cat(green(paste0("Random effect not found. Please check that spelling is correct")))
        hRand[n] <- readline(sprintf("%s\n",green(paste("Enter the name of the random effect ", n, " "))))
      }
      pRand[n] <- which(colnames(data)== hRand[n])
  }}
  
  
  
## Descriptive analysis of traits, effects and contingency table -------------------------------------
  cat(green(paste0(c("See below the summary statitics of traits: ", hTrait),collapse=" ")))
  cat("\n")
  print(summary(data[,pTrait]))
  
  if(nCov != 0){
  cat(green(paste0(c("See below the summary statitics of covariates: ", hCov),collapse=" ")))
  cat("\n")
  print(summary(data[,pCov]))
  }
  
  if ((nTreatment+nNoise)>1) {
      print(paste0(c("Contingency table across effects")))
      cat("\n")
      fi=nTreatment+nNoise
      NTables=((fi*fi) - fi)/2
      
      if ((nNoise != 0) && (nTreatment !=0)) {FE<-data.frame(data[,c(pTreatment, pNoise)])
      } else if (nNoise == 0)  {
        FE<-data.frame(data[,c(pTreatment)])
      } else if (nTreatment == 0)  {
          FE<-data.frame(data[,c(pNoise)])} 
      
        for (m in 1:ncol(FE)){FE[,m] <- paste(names(FE)[m],FE[,m],sep="")}
        for (m in 1:ncol(FE)){
            for (j in (m+1):ncol(FE)) {
              if(j>ncol(FE)) break
              print(paste(colnames(FE)[m],colnames(FE)[j], sep=" vs "))
              print(table(FE[,m],FE[,j])) }
        }
  } else {
      print(paste0(c("Contingency tables cannot be created because there is none or only 1 effect")))}
    
  
## Define Priors --------------------------------------------------------------------------------------
#  set normal (0, 10*sd)
  
#  askPrior<-readline("Do you want to establish prior bounds? (Yes=Y, No=N) ")
#  if (askPrior =="Y" | askPrior =="y") {
#  ask_priorFixed<-readline("Do you want to ste up the priors for the fixed effects (1=Yes 0=No)" )
#    if (ask_priorFixed ==1) {
#    Which_priorFixed<-readline("Which fixed effect do you want to set up the prior for? " )
#    prior<-readline("Options: "normal(0,1)", "gaussian",")                
#    
#   lowerBound_Mean<-readline("Write the lower bound for the mean")
#    lowerBound_Mean<-readline("Write the upper bound for the mean")
#   lowerBound_Effects<-readline("Write one lower bound for all the effects (remember it should be negative)")
#    upperBound_Effects<-readline("Write one upper bound for all each effect (remember it should be positive)")
#    prior <- set_prior("normal(0,1)",class="b") #prior para todos los efectos fijos
#    ??set.prior
#  }
  
  
  
## Define MCMC  ----------------------------------------------------------------------------------------
  #askMCMC     <- readline(sprintf("%s\n",green("Do you want to establish the MCMC characteristics  (Enter Yes=Y or No=N) ?   ")))
  #while((askMCMC!= "y") && (askMCMC!= "Y") && (askMCMC!= "n") && (askMCMC!= "N") ==TRUE){ 
  #  print(paste("Weird response, try again!"))
  #  askMCMC     <- readline(sprintf("%s\n",green("Do you want to establish the MCMC characteristics  (Enter Yes=Y or No=N) ?   ")))}
  
  askMCMC <- "n" 
   if (askMCMC =="Y" | askMCMC =="y") {
      Seed   <- as.numeric(readline(sprintf("%s\n",green("Please enter a random seed   "))))
      chain  <- as.numeric(readline(sprintf("%s\n",green("How many chains do you want to launch   "))))
      iter   <- as.numeric(readline(sprintf("%s\n",green("Enter the chain length   "))))
      burnin <- as.numeric(readline(sprintf("%s\n",green("Enter the burn-in   "))))
      lag    <- as.numeric(readline(sprintf("%s\n",green("Enter the lag between samples   "))))
    } else { #Default parameters
      Seed = 1234  
      chain = 2
      iter = 30000
      burnin = 5000
      lag = 10
    }

  
## Define Inferences to be calculated from Posterior Chains -------------------------------------------------------
  probHPD     <- as.numeric(readline(sprintf("%s\n",green("Enter the probability for the HPD interval (for example 0.95)   "))))
  probK       <- as.numeric(readline(sprintf("%s\n",green("Enter the probability for the guaranteed value k (for example 0.80)   "))))
  
  if (nTreatment!=0){
  askProbRel  <- readline(sprintf("%s\n",green("Do you want to calculate probability of contrasts being greater than a relevant value for some traits (Enter Yes=Y or No=N)   ? ")))
  while((askProbRel!= "y") && (askProbRel!= "Y") && (askProbRel!= "n") && (askProbRel!= "N") ==TRUE){ 
  print(paste("Weird response, try again!"))
  askProbRel  <- readline(sprintf("%s\n",green("Do you want to calculate probability of contrasts being greater than a relevant value for some traits (Enter Yes=Y or No=N)   ? ")))}
  if (askProbRel =="Y" |askProbRel =="y") {
    rValue<-NULL
      for(n in 1:nTrait){
        if (nTrait>1){
            askValue <-readline(sprintf("%s\n",green(paste("Do you want to calculate it for Trait ",hTrait[n], " (Enter Yes=Y or No=N)  ? "))))  
            if (askValue =="Y" | askValue =="y") {
            rValue[n] <- as.numeric(readline(sprintf("%s\n",green(paste("Enter a relevant value for Trait ", hTrait[n], " (in case of ratio needs to be >1) ")))))
            }else{rValue[n]<-0}
        }else{rValue[n] <- as.numeric(readline(sprintf("%s\n",green(paste("Enter a relevant value for Trait ", hTrait[n], " (in case of ratio needs to be >1) ")))))}
        }}
      
  
    if (askProbRel =="Y" |askProbRel =="y") {       
      askProbSimil <-  readline(sprintf("%s\n",green("Do you want to calculate probability of similarity [-r, r] for some traits (Enter Yes=Y or No=N)   ? ")))
      }
  }   
  
    # askInterval <- readline("Do you want to calculate the probability of intervals [a,b] for some traits (Enter Yes=Y or No=N) ? ")
    # if (askInterval ==1) {
    #   rInterval<-matrix(ncol=nTrait, nrow=2)
    #     for(n in 1:nTrait){
    #     askInterval <-readline(paste("Do you want to calculate it for Trait ",n, " (Enter Yes=Y or No=N)  ? "))  
    #       if (askInterval ==1) {
    #         rInterval[1,n] <- as.numeric(readline(paste("Enter value a for Trait ", n, " ")))
    #         rInterval[2,n] <- as.numeric(readline(paste("Enter value b for Trait ", n, " ")))
    #       }else{
    #       rInterval[1,n]<-0
    #       rInterval[2,n]<-0}
    #     }  
    # }
  
  
  
  
## Write down the formula of the equation ------------------------------------------------------------------
 # remove(eq.T,eq.C,eq.I,eq.R,eq.N, eq.total)
  eq.T <- NULL #Part of the equation for Treatment
  if(nTreatment != 0){eq.T <- paste(hTreatment,collapse =" + ")}
  eq.N <- NULL #Part of the equation for Noise
  if(nNoise != 0){eq.N <- paste(hNoise,collapse =" + ")}
  eq.I <- NULL #Part of the equation for Interaction
  
  if(nInter != 0){
     n <- 1
    while(n <= nInter){
      eq.I[n] <- paste(hInter[1,n],":",hInter[2,n],sep ="")
      n <- n+1}
    eq.I <- paste(eq.I,collapse =" + ")
    }
  eq.C <- NULL 
  eq.C.name <- NULL
  if(nCov != 0){
    eq.C.name <- paste(paste("b*",hCov,sep=""),collapse =" + ")
    eq.C <- paste(hCov,collapse =" + ")
  }
  eq.R <- NULL #Part of the equation for Random
  eq.R.name<-NULL
  if(nRand != 0){
    n <- 1
    while(n <= nRand ){
      eq.R[n] <- paste("(1|",hRand[n],")",sep="")
      eq.R.name[n] <- paste("Random(",hRand[n],")",sep="")
      n <- n+1}
  eq.R <- paste(eq.R,collapse =" + ")
  }
  
  temp4<-c(eq.T,eq.N,eq.I,eq.C,eq.R)
  eq.total <- paste(temp4,collapse = " + ")
  
  temp4.1<-c(eq.T,eq.N,eq.I,eq.C.name,eq.R.name)
  eq.name <- paste(temp4.1,collapse = " + ")
  cat("\nModel equation for all Traits is :", red(paste("y = ",eq.name,sep="")))
  cat("\n")
  
## Run the model --------------------------------------------------------------------------------------------
  #remove(eq.multiple)
  # #if(nTrait > 1){
  # #  brm.equation <- list()
  # #  for (var in 1:nTrait) {
  #     brm.equation <- append(brm.equation,(formula(paste(names(data)[pTrait[var]],"~",eq.total,collapse = ""))))
  #   }
  #   eq.multiple <- mvbf(flist = brm.equation,rescor = F)
  #   }else{
  #   eq.multiple <- brmsformula(formula(paste(names(data)[pTrait[1]],"~",eq.total,collapse = ""))) #original [pRes[var]] y he quitado var
  #   }
  
  #Declare output files with results from all traits
  LSMeans.total  <-NULL
  Effects.total  <-NULL
  Contrasts.total<-NULL
  LSMeansNoise.total  <-NULL
  Covariate.total<-NULL
  SD_Random.total<-NULL
  SD_E.total     <-NULL
  MeanModel.total<-NULL
  SDTrait.total  <-NULL
  Inf_PM.total <-NULL
  Inf_PModel.total<-NULL
  Inf_PCov.total<-NULL
  Inf_PE.total<-NULL
  Inf_PC.total<-NULL  
  Inf_PRandom.total<-NULL  
  
for (u in 1:nTrait){
  cat("\n")
  cat(green(paste0(c("\nAnalysis for Trait ",hTrait[u]," in progress ")))) 
  cat("\n")
  eq.multiple <- brmsformula(formula(paste(names(data)[pTrait[u]],"~",eq.total,collapse = "")))
  model <- brm(eq.multiple,
             data    = data, 
             family  = gaussian(), 
             iter    = iter, 
             chains  = chain, 
             warmup  = burnin,
             thin    = lag,
             control = list(adapt_delta = 0.99),
             silent  = 2,
             refresh = 0,
             backend = "cmdstanr",
             threads = threading(5),
             seed    = 1234)
            

 
  
## Get Convergence ------------------------------------------------------------------------------------------
    conver<-data.frame(summary(model)$fixed)
    conver$convergence<-0
    conver$convergence<-ifelse(conver$Rhat<1.05 &conver$Rhat>0.95,"OK","FAIL")
    if ("FAIL"%in%conver$convergence==TRUE) {
    print(paste("Convergence not reached:",which(conver$convergence=="FAIL"),sep=" "), collapse = "")}
    
    
## Get Posterior chains of Means, effects, contrasts, mean,SD, VarE, VarR ---------------------------------------
  # To program LSMEANS with interactions: Estimable function to compute the means https://www.sfu.ca/sasdoc/sashtml/stat/chap30/sect39.htm
    
    LSMeans <- NULL
    Effects <- NULL
    Contrasts <-NULL
    LSMeansNoise <- NULL
    Covariate<-NULL
    MeanModel<-NULL
    SD_E<-NULL
    SD_Random<-NULL
    SDTrait<-NULL
    
    modelfit= as_draws_df(model)
    head(modelfit)
    colnames(modelfit)
    
    if (nTreatment!=0){
    for (i in 1:length(hTreatment)){
  
      #LSMeans of Treatments
          epred <- emmeans(model, hTreatment[i], epred = TRUE)
          iter_posterior <- gather_emmeans_draws(epred)
          iter_bylevels <- iter_posterior %>%
            pivot_wider(names_from = hTreatment[i] , values_from = ".value")
          iter_bylevels<-data.frame(iter_bylevels[-c(1:3)])
          name.frame <- paste(hTreatment[i],c(1:nlevels_Treatment[i]),sep="")
          colnames(iter_bylevels)<-c(paste(name.frame,hTrait[u],sep="."))
          LSMeans<-data.frame(c(LSMeans, iter_bylevels))
          dim(LSMeans)
          
      #Effects: Center LSMeans of Treatments to get what Old Rabbit calls "Effects"
          temp3<-iter_bylevels-rowMeans(iter_bylevels)
          Effects<-data.frame(c(Effects, temp3))
      
      #Contrasts of Treatments
          NContrasts=((nlevels_Treatment[i]*nlevels_Treatment[i]) - nlevels_Treatment[i])/2 
          k=0 #My counter
          temp1=matrix(0,nrow=NContrasts, ncol=3)
          colnames(temp1)=c("contrast","i_level1","j_level2")
              for (m in 1:nlevels_Treatment[i]){  
                for (j in (m+1):nlevels_Treatment[i]) { 
                  if(j>nlevels_Treatment[i]) break
                  k=k+1
                  if(askCompare=="D"|askCompare=="d"){
                    temp2<-as.matrix(iter_bylevels[,m]-iter_bylevels[,j])
                    colnames(temp2)<-paste(name.frame[m],"-",name.frame[j],".",hTrait[u],sep = "")}else{
                    temp2<-as.matrix(iter_bylevels[m]/iter_bylevels[j])
                    colnames(temp2)<-paste(name.frame[m],"/",name.frame[j],".",hTrait[u],sep="")}  
                  Contrasts<-cbind(Contrasts, temp2)
                  temp1[k,1]=k
                  temp1[k,2]=m
                  temp1[k,3]=j
                }}
            }
    
          #Save chains for all traits in the same file
          LSMeans.total <- data.frame(c(LSMeans.total,LSMeans))
          Effects.total <- data.frame(c(Effects.total,Effects))
          Contrasts.total <- cbind(Contrasts.total,Contrasts)
         
    } 
         
    # LSMeans of Noise (compute it just in case)
    if (nNoise!= 0) {
     for (i in 1:length(hNoise)){
          epred_Noise <- emmeans(model, hNoise[i], epred = TRUE)
          iter_posterior_Noise <- gather_emmeans_draws(epred_Noise)
          iter_bylevels_Noise <- iter_posterior_Noise %>%
            pivot_wider(names_from = hNoise[i] , values_from = ".value")
          iter_bylevels_Noise<-data.frame(iter_bylevels_Noise[-c(1:3)])
          name.frame <- paste(hNoise[i],c(1:nlevels_Noise[i]),sep="")
          colnames(iter_bylevels_Noise)<-c(paste(name.frame,hTrait[u],sep="."))
          LSMeansNoise<-data.frame(c(LSMeansNoise, iter_bylevels_Noise))
        }
      LSMeansNoise.total<-cbind(LSMeansNoise.total,LSMeansNoise)
    }
    
    
    # Covariates
        if (nCov!= 0) {
        for (i in 1:nCov){
          pos<-which(colnames(modelfit)==paste0("b_",hCov[i]))
          temp5<-modelfit[,pos]
          colnames(temp5)<-paste(colnames(modelfit[,pos]),hTrait[u], sep="_")
          Covariate<-data.frame(c(Covariate,temp5))
        }
            #Save chains for all traits in the same file
            Covariate.total <- data.frame(c(Covariate.total,Covariate))
        }   
        
    # Mean of the model
    if (nTreatment!=0){MeanModel=data.frame(rowMeans(LSMeans[1:nlevels_Treatment[1]])) #If there is a treatment, compute the mean from LSMeans
    } else if ((nTreatment==0)&&(nNoise!=0)){ #If there is no treatment but noise, compute the mean from LSMeansNoise
    MeanModel=data.frame(rowMeans(LSMeansNoise[1:nlevels_Noise[1]]))
    } else if ((nTreatment==0)&&(nNoise==0)&&(nCov!=0)) {
    #If there is no treatment or noise but there are covariates do it this way
      MeanModel<-NULL 
      for (i in 1:nrow(Covariate)){
          temp11=as.matrix(Covariate[i,])
          MeanModel[i]= mean(data[,pTrait[u]]- (temp11 %*% t(data[,pCov])) + modelfit$b_Intercept[i], na.rm=TRUE)}
        }
    } else if ((nTreatment==0)&&(nNoise==0)&&(nCov==0)&&nRand!=0) {
      MeanModel=data.frame(modelfit$b_Intercept) 
    }  
      MeanModel.total<-data.frame(c(MeanModel.total,MeanModel))
       
     
    # Residual Variance 
        pos2<-which(colnames(modelfit)=="sigma")
        SD_E<-modelfit[,pos2]
        SD_E.total<-data.frame(c(SD_E.total,SD_E))
    
            
    # Random effect Variance
    if (nRand!= 0) {
      SD_Random<-NULL
      for (i in 1:nRand){
        temp7=data.frame(apply(modelfit[,grep(paste0("r_", hRand[i]),colnames(modelfit))], 2, sd))
        colnames(temp7)<-paste("SD_Random",hRand[i], sep="_")
        SD_Random<-data.frame(c(SD_Random, temp7))
      }
      #Save chains for all traits in the same file
      SD_Random.total <- data.frame(c(SD_Random.total,SD_Random))
    }
      
    # SD of the trait
    if (nRand!= 0) {SDTrait=sum(SD_Random)+SD_E
      }else{SDTrait=SD_E}
    SDTrait.total<-data.frame(c(SDTrait.total,SDTrait))  
          
    
## Obtain Inferences from Posterior Chains  ------------------------------------------------------------------------------------ 
    ChainLength<-((iter-burnin)/lag)*chain
    
      #For Mean, SDTrait, SDE
       Inf_PModel<-matrix(ncol=4, nrow=3)
       colnames(Inf_PModel)<-c("Estimate","sd","HPD_1","HPD_2")
       rownames(Inf_PModel)<-c(paste("Mean",hTrait[u],sep="_"),paste("Sd",hTrait[u],sep="_"),paste("Residual_SD",hTrait[u],sep="_"))
       File=cbind(MeanModel, SDTrait, SD_E)

       for (i in 1:ncol(File)) {
         Inf_PModel[i,1]=median(File[,i])
         Inf_PModel[i,2]=sd(File[,i])
         Inf_PModel[i,3]=hdi(File[,i], credMass = probHPD)[1]
         Inf_PModel[i,4]=hdi(File[,i], credMass = probHPD)[2]
       }
       Inf_PModel.total <- data.frame(rbind(Inf_PModel.total,Inf_PModel))

       
      #For Covariates  
       if (nCov!= 0){
              Inf_PCov<-matrix(ncol=7, nrow=ncol(Covariate))
              colnames(Inf_PCov)<-c("Estimate","sd","HPD_1","HPD_2","P0","k","k_Guaranteed")  
              rownames(Inf_PCov)<-colnames(Covariate)
              for (i in 1:ncol(Covariate)) {
                
                # Median, mean, sd
                Inf_PCov[i,1]=median(Covariate[,i])
                Inf_PCov[i,2]=sd(Covariate[,i])
                
                # HPD, Prob 0
                Inf_PCov[i,3]=hdi(Covariate[,i], credMass = probHPD)[1] #HPD95% LOWER BOUND
                Inf_PCov[i,4]=hdi(Covariate[,i], credMass = probHPD)[2] #HPD95% UPPER BOUND
                if (Inf_PCov[i,1] > 0){P=length(which(Covariate[,i] > 0))} else {P=length(which(Covariate[,i] < 0))} #Prob 0
                Inf_PCov[i,5]=P/ChainLength
              
              # Guaranteed value 
                temp9=sort(Covariate[,i])
                Inf_PCov[i,6]=probK
                Inf_PCov[i,7]= temp9[round(ChainLength*(1-probK),0)]
                if (sign(Inf_PCov[i,1]) ==  sign(Inf_PCov[i,7])) {Inf_PCov[i,7]=Inf_PCov[i,7]} else {Inf_PCov[i,7]=NA}
                }   #End inferences from posterior chains of covariates
                
              #Save inferences for all traits in the same file
              Inf_PCov.total <- data.frame(rbind(Inf_PCov.total,Inf_PCov))
      } 

  if (nTreatment !=0){ 
        
        #For Effects
        Inf_PE<-matrix(ncol=7, nrow=ncol(Effects))
        colnames(Inf_PE)<-c("Estimate","sd","HPD_1","HPD_2","P0","k","k_Guaranteed")  
        rownames(Inf_PE)<-colnames(Effects)
        for (i in 1:ncol(Effects)) {
          
          # Median, mean, sd
          Inf_PE[i,1]=median(Effects[,i])
          Inf_PE[i,2]=sd(Effects[,i])
          
          # HPD, Prob 0
          Inf_PE[i,3]=hdi(Effects[,i], credMass = probHPD)[1] #HPD95% LOWER BOUND
          Inf_PE[i,4]=hdi(Effects[,i], credMass = probHPD)[2] #HPD95% UPPER BOUND
          if (Inf_PE[i,1] > 0){P=length(which(Effects[,i] > 0))} else {P=length(which(Effects[,i] < 0))} #Prob 0
          Inf_PE[i,5]=P/ChainLength
    
        # Guaranteed value 
          Inf_PE[i,6]=probK
          temp9=sort(Effects[,i])
          if (Inf_PE[i,1] > 0) {Inf_PE[i,7]= temp9[round(ChainLength*(1-probK),0)]} else {Inf_PE[i,7]= temp9[round(ChainLength*probK,0)]}
          if (sign(Inf_PE[i,1]) ==  sign(Inf_PE[i,7])) {Inf_PE[i,7]=Inf_PE[i,7]} else {Inf_PE[i,7]=NA}
      
          
        }   #End inferences from posterior chains of effects
    
        #Save inferences for all traits in the same file
        Inf_PE.total <- data.frame(rbind(Inf_PE.total,Inf_PE))
        
    
        #For LSMEANS
        Inf_PM<-matrix(ncol=7, nrow=ncol(LSMeans))
        colnames(Inf_PM)<-c("Estimate","sd","HPD_1","HPD_2","P0","k","k_Guaranteed")  
        rownames(Inf_PM)<-colnames(LSMeans)
        for (i in 1:ncol(LSMeans)) {
          
          # Median, mean, sd
          Inf_PM[i,1]=median(LSMeans[,i])
          Inf_PM[i,2]=sd(LSMeans[,i])
          
          # HPD, Prob 0
          Inf_PM[i,3]=hdi(LSMeans[,i], credMass = probHPD)[1] #HPD95% LOWER BOUND
          Inf_PM[i,4]=hdi(LSMeans[,i], credMass = probHPD)[2] #HPD95% UPPER BOUND
          if (Inf_PM[i,1] > 0){P=length(which(LSMeans[,i] > 0))} else {P=length(which(LSMeans[,i] < 0))} #Prob 0
          Inf_PM[i,5]=P/ChainLength
          
          # Guaranteed value 
          Inf_PM[i,6]=probK
          temp9=sort(LSMeans[,i])
          if (Inf_PM[i,1] > 0) {Inf_PM[i,7]= temp9[round(ChainLength*(1-probK),0)]} else {Inf_PM[i,7]= temp9[round(ChainLength*probK,0)]}
          if (sign(Inf_PM[i,1]) ==  sign(Inf_PM[i,7])) {Inf_PM[i,7]=Inf_PM[i,7]} else {Inf_PM[i,7]=NA}
          
          
        }   #End inferences from posterior chains of effects
        
        #Save inferences for all traits in the same file
        Inf_PM.total <- data.frame(rbind(Inf_PM.total,Inf_PM))
        
        
      #For Contrasts
        Inf_PC<-matrix(ncol=10, nrow=ncol(Contrasts))
        colnames(Inf_PC)<-c("Estimate","sd","HPD_1","HPD_2","P0","r","Pr","Psimil","k","k_Guaranteed")  
        rownames(Inf_PC)<-colnames(Contrasts)
        for (i in 1:ncol(Contrasts)) {
          hist(Contrasts[,i])
          # Median, mean, sd
          Inf_PC[i,1]=median(Contrasts[,i])
          Inf_PC[i,2]=sd(Contrasts[,i])
          
          # HPD, Prob 0
          Inf_PC[i,3]=hdi(Contrasts[,i], credMass = probHPD)[1] #HPD95% LOWER BOUND
          Inf_PC[i,4]=hdi(Contrasts[,i], credMass = probHPD)[2] #HPD95% UPPER BOUND
          if (askCompare=="D" | askCompare=="d"){
            if (Inf_PC[i,1] > 0) {P=length(which(Contrasts[,i] > 0))} else {P=length(which(Contrasts[,i] < 0))} #Prob 0
          }else{
            if (Inf_PC[i,1] > 1) {P=length(which(Contrasts[,i] > 1))} else {P=length(which(Contrasts[,i] < 1))} #Prob 1 
          }
          Inf_PC[i,5]=P/ChainLength
          
          # Probability of relevance 
          
          if(askProbRel=="Y"|askProbRel=="y"){ 
            Inf_PC[i,6]=rValue[u]
            if (askCompare=="D" | askCompare=="d"){
              if (Inf_PC[i,1] > 0) {Pr=length(which(Contrasts[,i] > rValue[u]))} else {Pr=length(which(Contrasts[,i] < -rValue[u]))}
            }else{
              if (Inf_PC[i,1] > 1) {Pr=length(which(Contrasts[,i] > rValue[u]))} else {Pr=length(which(Contrasts[,i] < 1/rValue[u]))}
            }
            Inf_PC[i,7]=Pr/ChainLength
          }
          
          # Probability of similarity
          if(askProbRel=="Y"|askProbRel=="y"){ 
          if(askProbSimil=="Y" | askProbSimil=="y"){
            if (askCompare=="D" | askCompare=="d"){
                 Inf_PC[i,8]=sum(Contrasts[,i] < rValue[u] & Contrasts[,i] > -rValue[u])/ChainLength} 
                 }else{
                 Inf_PC[i,8]=sum(Contrasts[,i] > 1/rValue[u] & Contrasts[,i] < rValue[u])/ChainLength
                 }
          }
          
         # Guaranteed value 
          temp6=sort(Contrasts[,i])
          Inf_PC[i,9]=probK
          if (askCompare=="D" | askCompare=="d"){
                if (Inf_PC[i,1] > 0) {Inf_PC[i,10]= temp6[round(ChainLength*(1-probK),0)]} else {Inf_PC[i,10]= temp6[round(ChainLength*probK,0)]}
              }else{
                if (Inf_PC[i,1] > 1) {Inf_PC[i,10]= temp6[round(ChainLength*(1-probK),0)]} else {Inf_PC[i,10]= temp6[round(ChainLength*probK,0)]}
              }
          if (askCompare=="D" | askCompare=="d"){ #Display only those K80 with the same sign as its median
                if (sign(Inf_PC[i,1]) ==  sign(Inf_PC[i,10])) {Inf_PC[i,10]=Inf_PC[i,10]} else {Inf_PC[i,10]=NA}
              }else{  
                if ((Inf_PC[i,1]>1)==(Inf_PC[i,10]>1)) {Inf_PC[i,10]=Inf_PC[i,10]} else {Inf_PC[i,10]=NA}
              }   
          
          #Interval [a,b]
          #if(askInterval==1){
          #        Inf_PC[i,10]=sum(temp6 > rInterval[1,u] & temp6 < rInterval[2,u])/ChainLength}
                
      } #End inferences from posterior chains of contrasts 
        
          #Save inferences for all traits in the same file
          Inf_PC.total <- data.frame(rbind(Inf_PC.total,Inf_PC))  
  }
       
    # For Random effect
      if (nRand!= 0){
        Inf_PRandom<-matrix(ncol=4, nrow=nRand)
        colnames(Inf_PRandom)<-c("SdRandomEstimate","Sd","HDP_1","HDP_2")
        rownames(Inf_PRandom)<-paste(hRand, hTrait[u],sep="_")
        for (i in 1:nRand){
          Inf_PRandom[i,1]=median(SD_Random[,i])
          Inf_PRandom[i,2]=sd(SD_Random[,i])
          Inf_PRandom[i,3]=hdi(SD_Random[,i], credMass = probHPD)[1]
          Inf_PRandom[i,4]=hdi(SD_Random[,i], credMass = probHPD)[2]
          }
        Inf_PRandom.total <- data.frame(rbind(Inf_PRandom.total,Inf_PRandom)) 
      }
      
cat("\n")     
cat(green(paste0(c("\nAnalysis for Trait ",hTrait[u]," done "))))  
cat("\n")

} #Repeat for every trait u
  
  
## Write down the outputs  ------------------------------------------------------------------------------------ 

    #Write Posterior Chains
    write.csv(MeanModel.total, file="PosteriorChain_Means.csv")
    write.csv(SD_E.total, file="PosteriorChain_ResidualVariance.csv")
    if(nRand != 0){
      write.csv(SD_Random.total, file="PosteriorChain_RandomEffVariances.csv")
    }
    if(nCov != 0){
      write.csv(Covariate.total, file="PosteriorChain_Covariate.csv")
    }
    if(nTreatment !=0){
    write.csv(LSMeans.total, file="PosteriorChain_Means.csv")
    write.csv(Effects.total, file="PosteriorChain_Effects.csv")
    write.csv(Contrasts.total, file="PosteriorChain_Contrasts.csv")
    }
    if(nNoise !=0){
      write.csv(LSMeansNoise.total, file="PosteriorChain_MeansNoise.csv")
    }
    
    #Write Inferences files
    write.csv(Inf_PModel.total , file=paste0("Results_Model.csv"))
    if(nRand != 0){
      write.csv(Inf_PRandom.total, file=paste0("Results_RandomEffects.csv"))
    }
    if(nCov != 0){
      write.csv(Inf_PCov.total, file=paste0("Results_Covariates.csv"))
    }
    if(nTreatment !=0){
    write.csv(Inf_PE.total, file=paste0("Results_Effects.csv"))
    write.csv(Inf_PM.total, file=paste0("Results_Means.csv"))
    write.csv(Inf_PC.total, file=paste0("Results_Contrasts.csv"))
    }
    
cat("\n")
cat(green("Progam finsihed!! :) "))
cat("\n")

    #Save Plots
    pdf.name <- readline(sprintf("%s\n",green("Enter the name for the pdf file containing the graphs: ")))
    pdf(paste(pdf.name,".pdf",sep=""), height=8, width=8)
    table.total <- NULL
    
    # for (m in 1:ncol(FE)){
    #   for (j in (m+1):ncol(FE)) {
    #     if(j>ncol(FE)) break
    #     table <- ggtexttable(table(FE[,m],FE[,j]),theme=ttheme("light"))
    #     table <- table %>%
    #       tab_add_title(text = paste(colnames(FE)[m],colnames(FE)[j], sep=" vs "), face = "bold",size = 12, padding = unit(2, "line"))
    #     grid.draw(table)
    #   }
    # }
    
    Contrasts.total <- data.frame(Contrasts.total)
    for (trait in 1:length(hTrait)) {
      stat <- t(matrix(round(summary(data[,pTrait[trait]])[1:6],2)))
      colnames(stat) <- c("Min","1st Qu.","Median","Mean","3rd Qu.","Max")
      table <- ggtexttable(stat, rows = NULL,theme=ttheme("light"))
      table <- table %>%
        tab_add_title(text = paste("Descriptive analyses for trait ",hTrait[trait],sep=""), face = "bold",size = 16, padding = unit(2, "line"))
      plot <- ggplot(data,aes(data[,hTrait[trait]]))
      hist <- plot +
        geom_histogram(stat = "bin",color="white",fill="skyblue")+
        theme_classic()+
        xlab(hTrait[trait])
      box<- plot +
        geom_boxplot(color="black",fill="skyblue")+
        theme_classic()+
        xlab(hTrait[trait])+
        coord_flip()
      total.plot <- ggarrange(hist,box,ncol=2,nrow=1)
      
      Contrasts.tmp1 <- Contrasts.total[,grep(hTrait[trait],names(Contrasts.total))]
      contrast.plot <- NULL
     
       for(treat in 1:dim(Contrasts.tmp1)[2]){
        if(askCompare=="D"|askCompare=="d"){
          contrast.tmp <- data.frame(Contrast=rep(names(Contrasts.tmp1)[treat],dim(Contrasts.tmp1)[1]),Difference = Contrasts.tmp1[,treat])
          contrast.plot <- rbind(contrast.plot,contrast.tmp)
          cPlot<- ggplot(contrast.plot,aes(y = Contrast, x = Difference, fill = after_stat(abs(x) < rValue[trait]))) +
            stat_halfeye() +
            theme_classic() +
            geom_vline(xintercept = c(-rValue[trait], rValue[trait]), linetype = "dashed") +
            scale_fill_manual(values = c("gray80", "skyblue"))+
            theme(legend.position="none")
        }
        if(askCompare=="R"|askCompare=="r"){
          contrast.tmp <- data.frame(Contrast=rep(names(Contrasts.tmp1)[treat],dim(Contrasts.tmp1)[1]),Ratio = Contrasts.tmp1[,treat])
          contrast.plot <- rbind(contrast.plot,contrast.tmp)
          colnames(contrast.plot) <- c("Contrast","Ratio")
          cPlot<- ggplot(contrast.plot,aes(y = Contrast, x = Ratio, fill = after_stat(x > rValue[trait] | x < 1/rValue[trait]))) +
            stat_halfeye() +
            theme_classic() +
            geom_vline(xintercept = c(1/rValue[trait], rValue[trait]), linetype = "dashed") +
            scale_fill_manual(values = c("gray80", "skyblue"))+
            theme(legend.position="none")
        }
      }
      x<- ggarrange(total.plot,table,ncol=1,nrow=2,heights=c(3,2), widths=c(2,1))
      grid.draw(x)
      grid.draw(cPlot)
    }
    dev.off()

#} #Close Rabbit



