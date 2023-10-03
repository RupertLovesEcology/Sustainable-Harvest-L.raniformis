

# for Geoff's request I will remove a set number of spawn tad adult perfectly from the population
# try 0 - 100 in 5s (lets me retain the 21 column recieving arrays)
## Remove everything
rm(list = ls())

# setup directories
source("C:/workspace/math0286/R/win-library/3.6/matrixOperators.r")

## libraries ####
# library(DescTools)
library(pracma)
library(deSolve)
library(reshape2)
library(data.table)
library(lhs)
library(snow)
library(doSNOW)
library(gbm)
library(foreach)
library(iterators)
library(parallel)
library(dismo)
library(ggplot2)
library(dplyr)

# beta distribution shape parameter estimator function
##Generates an alpha and a beta value to inform the beta distribution 
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta)) }



# call to remove adult frogs by number. Call using - n.mat[2:6,i+1] <- adultRemoval(adCollect.removal,n.mat[2:6,i+1]) 
adultRemoval <- function(collectionNo,demog) {
  if(collectionNo >= sum(demog)) { 
    result <- c(0,0,0,0,0)
  } else {
    adAge <- demog/sum(demog)
    remove <- floor(adAge*collectionNo)
    remain <- collectionNo - sum(remove)
    result <- demog - remove
    
    while (remain > 0) {
      adAge <- result
      adAge <- cumsum(adAge/sum(adAge))
      ageof <- runif(1,0,1)
      for (l in 1:5) {
        if (l == 1 && ageof > 0 && ageof <= adAge[1]) {
          result[l] <- result[l] - 1
          break
        } else if (ageof > adAge[l-1] && ageof <= adAge[l]) {
          result[l] <- result[l] - 1
        }
      }
      remain <- remain -1
    } 
  } 
  if (sum(result < 0) > 0) {stop("adultReduction function is fucked")}
  return(result)
}

# call to remove adult frogs by proportion. Call using - n.mat[2:6,i+1] <- adultReduction(adCollect.reduction,n.mat[2:6,i+1])
adultReduction <- function(collectionProp,demog) {
  collectionNo <- round(collectionProp*sum(demog))
  if(collectionNo >= sum(demog)) { 
    result <- c(0,0,0,0,0)
  } else {
    adAge <- demog/sum(demog)
    remove <- floor(adAge*collectionNo)
    remain <- collectionNo - sum(remove)
    result <- demog - remove
    
    while (remain > 0) {
      adAge <- result
      adAge <- cumsum(adAge/sum(adAge))
      ageof <- runif(1,0,1)
      for (l in 1:5) {
        if (l == 1 && ageof > 0 && ageof <= adAge[1]) {
          result[l] <- result[l] - 1
          break
        } else if (ageof > adAge[l-1] && ageof <= adAge[l]) {
          result[l] <- result[l] - 1
        }
      }
      remain <- remain -1
    } 
  } 
  if (sum(result < 0) > 0) {stop("adultReduction function is fucked")}
  return(result)  
}


## encapsulate the core metapopulation model as a function for the sustainable harvest model
sbf_sim <- function(input, dir.nm, rowNum) {
  
  ## assign all parameter values
  for (d in 1:ncol(input)) {assign(names(input)[d], input[,d])}
  # complementary log-log
  cloglog <- function(x) log(-log(1-x))
  
  ## Select the site and life stage to harvest #######
  siteList <- c("Hogwash")
  stageList <- c("tad")
  
  # assign iterations ####
  iteration <- 200
  
  ## Set up the code for the model
  # order of 'XXXXXX.K.init' is c(spawnK, tadpoleK, juvenileK, adultK)
  ## Set up the sites (K values)
  NapNap.K.init <- c(141,1417,570,570)
  Hogwash.K.init <- c(28,270,108,108) 
  Epsom.K.init <- c(40,412,165,165)
  
  ## set time limit for projection in 1-yr increments CAN PROBABLY CULL SOME OF THESE YEAR AND TIME LINES
  yr.now <- 2020
  yr.end <- 2020 + 86
  #************************
  ##Note I have constrained t to 86 years
  t <- (yr.end - yr.now)
  yrs <- seq(yr.now,yr.end,1)
  longev <- 5
  age.vec <- seq(0,longev,1)
  lage <- length(age.vec)
  sex.ratio <- 0.5
  stages <- lage
  
  ## set population storage matrices n.mat and the annual population change matrix pop.mat
  n.mat <- array(data = 0, dim = c(stages,(t + 1)))
  rownames(n.mat) <- c("0-1 Yrs","1-2 Yrs","2-3 Yrs","3-4 Yrs","4-5 Yrs","5+ Yrs")
  popmat <- matrix(0,nrow=stages,ncol=stages)
  colnames(popmat) <- age.vec[1:stages]
  rownames(popmat) <- c("Fecundity","Survival to 1 year","Survival to 2 year","Survival to 3 year","Survival to 4 year","Survival to 5 year")
  
  ## fertility data 
  clutch.size.NOTUSED <- c(1885,3893,2448,3090,3191,3644,4563) # L. raniformis
  clutch.sd <- sd(clutch.size.NOTUSED)
  prop.breeding <- c(0,rep(breedProp,5))
  fert.mn <- mean(clutch.size.lr)*prop.breeding
  
  #duration data (eggs and tadpoles)
  hatch.sd <- 0.06448
  
  ##the below uses tadpole survival TO METAMORPHOSIS figures from Bull (C.signifera) *note the beta distribution
  tadpole.mn.1 <- mean(c(.15,.26))
  tadpole.sd.1 <- ((.26-tadpole.mn.1)+(tadpole.mn.1-.15))/2/1.96
  tadpole.mn.2 <- mean(c(.07,.56))
  tadpole.sd.2 <- ((.56-tadpole.mn.1)+(tadpole.mn.1-.07))/2/1.96
  tadpole.sd <- sqrt(tadpole.sd.1^2 + tadpole.sd.2^2)
  
  tomet.dur.iter <- round(hatch.dur + tadpole.dur) 
  toad.dur.iter <- 365 - tomet.dur.iter
  
  #Adult annual survival from Turner et al 2022:Epsom=permanent,NapNap=semi-permanent,Hogwash=ephemeral,-(mean,sd)
  Hogwash.ad.s.yr.sd <- 0.08
  
  ##Calculate survivals
  #Probability of egg hatching and tadpole surviving to metamorphosis
  hatch <- rnorm(1, mean=hatch.pr, sd=hatch.sd)
  if (hatch > 1)  { hatch <- 1 }
  if(hatch < 0) { hatch <- 0 }
  tomet.s.iter <- rnorm(1,tadpole.mn,tadpole.sd)
  if(tomet.s.iter < 0) { tomet.s.iter <- 0 }
  if(tomet.s.iter > 1) { tomet.s.iter <- 1 }
  
  #placeholder
  site <- siteList[1]
  
  #calculate daily probability of adult survival
  toad.s.iter <- rnorm(1,mean = get(paste0(site,".ad.s.yr")),sd = get(paste0(site,".ad.s.yr.sd")))
  if(toad.s.iter < 0) { toad.s.iter <- 0 }
  if(toad.s.iter > 1) { toad.s.iter <- 1 }
  toad.daily.s.iter <- nthroot(toad.s.iter , 365)
  toad.s.season.iter <- toad.daily.s.iter ^ toad.dur.iter
  toad.s.iter <- tomet.s.iter * toad.s.season.iter                     
  
  #Create the survival vector, to adult survival then 4 adult survivals
  ad.s.vec.iter <- rep(NA,5)
  for (s in 1:5) {
    # ad.s.yr <- rnorm(1,mean = get(paste0(site,".ad.s.yr"))[1],sd = get(paste0(site,".ad.s.yr"))[2])
    ad.s.yr <- rnorm(1,mean = get(paste0(site,".ad.s.yr")),sd = get(paste0(site,".ad.s.yr.sd")))
    if(ad.s.yr < 0) { ad.s.yr <- 0 }
    if (ad.s.yr > 1) { ad.s.yr <- 1 }
    ad.s.vec.iter[s] <- ad.s.yr
  }
  surv.iter <- c(toad.s.iter, ad.s.vec.iter)
  init.vec <- rep(NA,6)
  dim(init.vec) <- c(6,1)
  
  ##Populate the Matix (popmat) and create a failure matrix(popmat.fail) for years with no breeding
  diag(popmat[2:(stages), ]) <- surv.iter[-stages]
  popmat[stages,stages] <- 0 # surv.mn[stages] 
  popmat[1,] <- fert.mn * sex.ratio
  popmat <- popmat.orig <- popmat ## save original matrix as popmat.orig
  
  # Create popmats (don't need all of these some are relicts from code development e.g.popmat.orig)
  popmat.fail <- popmat.current <- popmat
  popmat.fail[1,] <- 0 
  
  #create the egg density correcting function and variables
  ## Ignore the warning here, we truncate the selection range so it is not an issue
  surv.mult.egg.up <- 1.03
  surv.mult.egg.upmid <- 1.03
  surv.mult.egg.mid <- 0.94
  surv.mult.egg.lo <- 0.8
  surv.mult.egg.lo.lo <- 0.2
  surv.mult.egg.lo.lo.lo <- 0.05
  
  K.egg.up <- 1
  K.egg.upmid <- 0.98
  K.egg.mid <- 0.90
  K.egg.lo <- 0.7
  K.egg.lo.lo <- 0.3
  K.egg.lo.lo.lo <- 0.01
  
  K.egg.vec <- c(K.egg.up,K.egg.upmid, K.egg.mid,K.egg.lo, K.egg.lo.lo, K.egg.lo.lo.lo)
  surv.mult.egg.vec <- rev(c(surv.mult.egg.up, surv.mult.egg.upmid, surv.mult.egg.mid, surv.mult.egg.lo, surv.mult.egg.lo.lo, surv.mult.egg.lo.lo.lo))
  DD.dat <- data.frame(K.egg.vec, surv.mult.egg.vec)
  
  #the formula for the function
  SS<-getInitial(surv.mult.egg.vec~SSlogis(K.egg.vec,alpha,xmid,scale),data=DD.dat)
  fit.expd.egg <- nls(surv.mult.egg.vec ~ a/((exp((b-K.egg.vec)/c)) + 1),
                      data = DD.dat,
                      algorithm = "port",
                      start = c(a = as.numeric(SS["alpha"]), b = as.numeric(SS["xmid"]), c = as.numeric(SS["scale"])),
                      trace = TRUE,
                      nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
  
  eggfunc <- function(x) (coef(fit.expd.egg)[1])/((exp((coef(fit.expd.egg)[2] - x)/coef(fit.expd.egg)[3]))+ 1)
  intArea <- integrate(eggfunc,lower=0,upper=1.02)                    
  area1 <- (as.numeric(intArea[1])/1.02)                    
  s.mult.egg.iter <- as.numeric(coef(fit.expd.egg)[1])/((exp((coef(fit.expd.egg)[2] - 1.02)/coef(fit.expd.egg)[3]))+ 1)        
  
  ## invoke a density-feedback function on tadpole survival to year 1
  # density feedback survival multiplier for tadpoles hinges on the density of other tadpoles in the pond
  #form of the curve from "Effect of Stocking Density on the Survival and Growth of 
  # Hoplobatrachus occipitalis (Günther, 1858) (Amphibia: Dicroglossidae) 
  # of Tadpoles Reared in Ponds from Benin, Godome, 2018, International Journal of Aquaculture"  &
  # An Analysis of Density Effects and Predation in Bufo Americanus Tadpoles
  # from Brockelman 1969
  surv.mult.up <- 1.0
  surv.mult.upmid <- 0.58
  surv.mult.mid <- 0.19
  surv.mult.lo <- 0.10
  
  K.up <- 1
  K.upmid <- 0.83
  K.mid <- 0.45
  K.lo <- 0.3
  
  K.tad.vec <- c(K.up,K.upmid, K.mid,K.lo)
  surv.mult.tad.vec <- rev(c(surv.mult.up, surv.mult.upmid, surv.mult.mid, surv.mult.lo))
  plot(K.tad.vec, surv.mult.tad.vec, pch=19)
  
  # Bleasdale
  # y = (a + bx)^(-1/c)
  DD.dat <- data.frame(K.tad.vec, surv.mult.tad.vec)
  param.init <- c(-2.41e-01, 1.54, 1.17)
  fit.expd.tad <- nls(surv.mult.tad.vec ~ (a + (b*K.tad.vec))^(-1/c), 
                      data = DD.dat,
                      algorithm = "port",
                      start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                      trace = TRUE,      
                      nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
  plot(K.tad.vec, surv.mult.tad.vec, pch=19, xlab="K", ylab="reduction Tadpole survival to 1 yr")
  K.pred.tad.vec <- seq(K.lo,1,0.01)
  pred.surv.tad.mult <- (as.numeric(coef(fit.expd.tad)[1]) + (K.pred.tad.vec * as.numeric(coef(fit.expd.tad)[2])))^(-1/as.numeric(coef(fit.expd.tad)[3]))
  lines(K.pred.tad.vec, pred.surv.tad.mult, lty=2, lwd=1, col="red")
  
  ## invoke a density-feedback function on Juvenile survival from year 1 to year 2
  # density feedback survival multiplier for juveniles hinges on the density of themselves (but more strongly than adults)
  ### Make an approximation of the curve
  surv.mult.up <- 1.1
  surv.mult.upmid <- 1
  surv.mult.midmidup <- 0.5
  surv.mult.mid <- 0.3
  surv.mult.midlo <- 0.19
  surv.mult.lo <- 0.10
  
  K.up <- 1
  K.upmid <- 0.89
  K.mid <- 0.75
  K.midlo <- 0.6
  K.midlolo <- 0.3
  K.lo <- 0.1
  
  K.juv.vec <- c(K.up,K.upmid, K.mid,K.midlo,K.midlolo,K.lo)
  surv.mult.juv.vec <- rev(c(surv.mult.up, surv.mult.upmid,surv.mult.midmidup,surv.mult.mid,surv.mult.midlo,surv.mult.lo))
  plot(K.juv.vec, surv.mult.juv.vec, pch=19,main = "the curve I want to emulate")
  
  D.dat <- data.frame(K.juv.vec, surv.mult.juv.vec)
  
  # use the deSolve package to determine your starting parameters for the nls function
  # see more here     https://datascienceplus.com/first-steps-with-non-linear-regression-in-r/
  SS<-getInitial(surv.mult.juv.vec~SSlogis(K.juv.vec,alpha,xmid,scale),data=D.dat)
  
  #the formula for the function
  fit.expd.juv <- nls(surv.mult.juv.vec ~ a/((exp((b-K.juv.vec)/c)) + 1),
                      data = D.dat,
                      algorithm = "port",
                      start = c(a = as.numeric(SS["alpha"]), b = as.numeric(SS["xmid"]), c = as.numeric(SS["scale"])),
                      trace = TRUE,
                      nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
  
  fit.expd.adult <- nls(surv.mult.juv.vec ~ a/((exp((b-K.juv.vec)/c)) + 1),
                        data = D.dat,
                        algorithm = "port",
                        start = c(a = as.numeric(SS["alpha"]), b = as.numeric(SS["xmid"]), c = as.numeric(SS["scale"])),
                        trace = TRUE,
                        nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
  
  # Track the generation that the run went extinct 
  # row is %extraction col is iteration
  reductionTitles <- c("0%","5%","10%","15%","20%","25%","30%","35%","40%","45%","50%","55%","60%","65%","70%","75%","80%","85%","90%","95%","100%")
  reducList <- c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1)
  extinctionYearEgg <- matrix(data=85, nrow = length(reductionTitles),ncol= iteration)
  rownames(extinctionYearEgg) <- reductionTitles
  
  yrExt <- finalAdpop <- finalAllpop <- minAdpop <- minAllpop <- extinctionYearEgg   
  minAdpop[] <- minAllpop[] <- 999999
  yrExt[] <- finalAdpop[] <- finalAllpop[] <- NA
  
  #arrays for holding rRange 
  rRange <-  matrix(data=NA, nrow = 3, ncol = 21)
  r.mat <- matrix(data=NA,nrow =iteration,ncol = 84)
  rStoch <- c(rep(0,iteration))
  dryTracker <- 0
  wetTracker <- 0
  extCheck <- survCheck <- rep(0,21)
  
  ## And NOW for the actual model!  ##################
  
  # Start of the site loop
  for (S in 1:length(siteList)) {
    site <- siteList[S]
    
    # set the site values for K and the init ####
    current.K.init <- get(paste0(site,".K.init"))
    
    # convert the number of spawning masses to the number of female eggs (uses mean clutch size)
    current.K.init[1] <- current.K.init[1] * (mean(clutch.size.lr) * 0.5)
    
    ## Import the pregenerated datasets ####
    # make sure naming matches convention below
    wetdry <- fread(file = paste0("C:/Workspace/SustHarvRevSep22/wetdry",site,".csv"), header = F, sep = ",", dec = ".")
    wetdry <- wetdry[2:10001,]
    wetdry <- array(as.numeric(unlist(wetdry)), dim=c(iteration,99))
    
    startPops <- fread(file = paste0("C:/Workspace/SustHarvRevSep22/startPops",site,"2.csv"), header = T, sep = ",", dec = ".")
    startPops <- array(as.numeric(unlist(startPops)), dim=c(6,iteration))
    
    wetExtinct <- dryExtinct <- matrix(data = 0, nrow =6,ncol = length(reductionTitles))
    
    colnames(dryExtinct) <- colnames(wetExtinct) <- reductionTitles
    rownames(dryExtinct) <- c("D","DD","DDD","DDDD","DDDDD","DDDDDD")
    rownames(wetExtinct) <- c("W","WW","WWW","WWWW","WWWWW","6+ W")
    untracedExtinct <- rep(0,21)
    
    #start the stage of collection loop, 1-egg,2-tad,3-adult
    for (stage in 2:2) {  
      stg <- stageList[stage]
      
      # These variables determine the rate, type, and life stage of collection. retained separately to allow multiple approaches
      stageCollect <- c("eggCollect","tadCollect","adCollect")
      eggCollect <- tadCollect <- adCollect <- F
      eggCollect.removal <- eggCollect.reduction <- tadCollect.removal <- tadCollect.reduction <- adCollect.removal <- adCollect.reduction <- 0 
      assign(stageCollect[stage],T)
      
      #start the reduc loop
      #this will go from 0% reduction to 100% reduction in increments of 5% BUT this is the GSA so I remove no tadpoles
      for (r in 1:1) {
        reduc <- reducList[r]
        r.mat[] <- NA
        
        # start the iteration loop
        for (iter in 1:iteration) {
          #reset the population matrix except for the starting populations (n.mat[,1,])
          n.mat[] <- 0
          n.mat[1:6,1] <- startPops[1:6,iter]
          dryTracker <- 0
          wetTracker <- 0
          
          #retained as separate values in case I wish to re-examine multiple approaches on the same visits
          if (eggCollect == T) { eggCollect.reduction <- reduc }     
          if (tadCollect == T) { tadCollect.reduction <- reduc } 
          if (adCollect == T)  { adCollect.reduction <- reduc  } 
          
          ## The Innermost Loop: run the current projection set up for 85 years 
          for (i in 1:86) {
            # store minimum populations
            if (i >= 2) {
              if (sum(n.mat[2:6,i]) < minAdpop[r,iter]) { minAdpop[r,iter] <- sum(n.mat[2:6,i])  }
              if (sum(n.mat[1:6,i]) < minAllpop[r,iter]) { minAllpop[r,iter] <- sum(n.mat[1:6,i])   }
            }
            
            # if the population is extinct do some accounting
            if (sum(n.mat[,i]) == 0) { 
              #  minAdpop[r,iter] <- 0
              minAdpop[r,iter] <- 0
              minAllpop[r,iter] <- 0
              finalAdpop[r,iter] <- 0
              finalAllpop[r,iter] <- 0
              yrExt[r,iter] <- i
              extCheck[r] <- extCheck[r] + 1
              if (dryTracker > 0 && wetTracker == 0) {
                dryExtinct[dryTracker,r] <- dryExtinct[dryTracker,r] + 1
              } else if (wetTracker > 0 && dryTracker == 0) {
                wetExtinct[wetTracker,r] <- wetExtinct[wetTracker,r] + 1
              } else if (dryTracker == 0 && wetTracker == 0){
                untracedExtinct[r] <- untracedExtinct[r] + 1
              } else {
                stop("551 wetTracker and dryTracker, have I written this correctly? At the time of writing I cannot check my code for 2 x 0 or 2 x 1")
              }
              break 
            }  
            
            # if it is the 85th year then break here  
            if (i==85) { 
              finalAdpop[r,iter] <- sum(n.mat[2:6,i])
              finalAllpop[r,iter] <- sum(n.mat[1:6,i])
              survCheck[r] <-  survCheck[r] + 1
              break 
            }
            
            # resample the durations of egg to preCollection, Collection tometamorph and toadult
            tohatch.dur.iter <- round(hatch.dur)
            tomet.dur.iter <- round(hatch.dur + tadpole.dur)
            preCollect.dur.iter <- round(2/3 * tomet.dur.iter,digits = 0)
            postCollect.dur.iter <- tomet.dur.iter - preCollect.dur.iter
            toad.dur.iter <- (365 - tohatch.dur.iter) - tomet.dur.iter
            
            ##Calculate our survivals
            #likelyhood of egg hatching and tadpole surviving to metamorphosis
            hatch <- rnorm(1, mean=hatch.pr, sd=hatch.sd)
            if (hatch > 1)  { hatch <- 1 }
            if (hatch < 0)  { hatch <- 0 }
            tomet.s.iter <- rnorm(1,tadpole.mn,tadpole.sd)
            if(tomet.s.iter < 0) { tomet.s.iter <- 0 }
            if(tomet.s.iter > 1) { tomet.s.iter <- 1 }
            tomet.s.iter <- tomet.s.iter * hatch
            tomet.daily.s.iter <- nthroot(tomet.s.iter, tomet.dur.iter)
            preCollect.s.iter <- tomet.daily.s.iter ^ preCollect.dur.iter
            postCollect.s.iter <- tomet.daily.s.iter ^ postCollect.dur.iter
            
            toad.s.iter <- rnorm(1,mean = get(paste0(site,".ad.s.yr")),sd = get(paste0(site,".ad.s.yr.sd")))
            if (toad.s.iter < 0) { toad.s.iter <- 0 } 
            toad.daily.s.iter <- nthroot(toad.s.iter , 365)
            toad.s.iter <- toad.daily.s.iter ^ toad.dur.iter
            
            #Create the survival vector (popmat) for the year (density dependence not considered yet)
            # NOTE 0-1 surv is set to 1 and is applied after the matrix multiplication to allow for Collection of tadpoles
            ad.s.vec.iter <- rep(NA,5)
            for (s in 1:5) {
              # ad.s.yr <- rnorm(1,mean = get(paste0(site,".ad.s.yr"))[1],sd = get(paste0(site,".ad.s.yr"))[2])
              ad.s.yr <- rnorm(1,mean = get(paste0(site,".ad.s.yr")),sd = get(paste0(site,".ad.s.yr.sd")))
              if(ad.s.yr < 0) {  ad.s.yr <- 0 }
              if(ad.s.yr > 1) {  ad.s.yr <- 1 }
              ad.s.vec.iter[s] <- ad.s.yr
            }
            surv.iter <- c(1, ad.s.vec.iter)
            
            # fert.iter <- round((runif(stages, min=clutch.size.lr[1], max=clutch.size.lr[2])) * prop.breeding, 0)
            fert.iter <- round((rnorm(stages, mean=clutch.size.lr, sd=clutch.sd)) * prop.breeding, 0)
            popmat[1,] <- fert.iter * sex.ratio 
            for (ll in 1:stages) { if(surv.iter[ll] < 0) { surv.iter[ll] <- 0 }}
            diag(popmat[2:(stages), ]) <- surv.iter[-stages]
            
            popmat.current <- popmat
            
            # Implement density dependence effects for each of four densities and feed into the popmat.all[,,]
            #note density effect on egg laying is applied after matrix multiplication 
            #  density feedback for tadpoles to first year is strong and driven by the number of tadpoles in the cohort
            s.mult.iter.tad <- 1
            s.mult.iter.juv <- 1
            s.mult.iter.ad <- 1
            
            # calculate density dependence for tadpoles growing into year 1 adults
            # this is implemented after matrix multiplication to allow for collection midway
            K.rel.tad <- (n.mat[1,i]/current.K.init[2])
            if (!is.nan(K.rel.tad)) {
              if (K.rel.tad > 1.4)  { K.rel.tad <- 1.4 }
              if (K.rel.tad <= 1.4) {
                s.mult.iter.tad <- (as.numeric(coef(fit.expd.tad)[1]) + (K.rel.tad * as.numeric(coef(fit.expd.tad)[2])))^(-1/as.numeric(coef(fit.expd.tad)[3]))                       
              } 
            } 
            
            #  apply Density dependence reductions to each of the 0-1 suvival probabilities
            tad.Ddep.daily <- nthroot(s.mult.iter.tad, (365 - tohatch.dur.iter))
            preCollect.s.iterDD <- preCollect.s.iter * (tad.Ddep.daily^preCollect.dur.iter)
            postCollect.s.iterDD <- postCollect.s.iter * (tad.Ddep.daily^postCollect.dur.iter)
            toad.s.iterDD <- toad.s.iter * (tad.Ddep.daily^toad.dur.iter)
          
            # instil density dependence for juveniles (1 - 2 years) is driven by the  number of yr 1 present/competing per 
            K.rel.juv <- (n.mat[2,i]/current.K.init[3])
            if (!is.nan(K.rel.juv)) {
              if (K.rel.juv > 0.65)  { K.rel.juv <- 0.65 }
              if (K.rel.juv <= 0.65) {
                s.mult.iter.juv <- as.numeric(coef(fit.expd.juv)[1])/((exp((coef(fit.expd.juv)[2] - K.rel.juv)/coef(fit.expd.juv)[3]))+ 1)
                popmat.current[3,2] <- popmat.current[3,2] * s.mult.iter.juv
              }
            }
            
            # instill density dependence for adults  driven by the  number of yr 1s emerging from Berven 2009
            K.rel.adult <- (sum(n.mat[2:6,i])/current.K.init[4]) 
            if (!is.nan(K.rel.juv)) {
              if (K.rel.adult > 0.65)  { K.rel.adult <- 0.65 }
              if (K.rel.adult <= 0.65) {
                s.mult.iter.ad <- as.numeric(coef(fit.expd.adult)[1])/((exp((coef(fit.expd.adult)[2] - K.rel.adult)/coef(fit.expd.adult)[3]))+ 1) 
                for (adgens in 4:6) {
                  popmat.current[adgens,(adgens-1)] <- (popmat.current[adgens,(adgens - 1)] * s.mult.iter.ad)  
                }
              }
            }
            
            # set popmat.fail for use if this year is dry
            popmat.fail <- popmat.current
            popmat.fail[1,] <- 0 
            
            if (wetdry[iter,i] == 1) { 
              dryTracker <- 0 
              ifelse(wetTracker < 6, wetTracker <- wetTracker + 1, wetTracker <- 6)
            }
            if (wetdry[iter,i] == 0) {
              wetTracker <- 0
              ifelse(dryTracker < 6, dryTracker <- dryTracker + 1, dryTracker <- 6)
            }
            
            ## The matrix multiplication step
            if (wetdry[iter,i] == 1 && eggCollect == T) {
              # remove adults then produce spawn then iterate the survival (without breeding) then paste the spawn in
              matA <- n.mat[,i]
              spawnCollect <- 0
              if (eggCollect.removal > 0) { spawnCollect <- eggCollect.removal }
              if (eggCollect.reduction > 0) { spawnCollect <- (eggCollect.reduction*(sum(matA[2:6])))  }
              matA[2:6] <- adultRemoval(spawnCollect,matA[2:6]) 
              matA <- round(popmat.current %*% matA, digits = 0)
              n.mat[,i+1] <- round(popmat.fail %*% n.mat[,i], digits = 0)
              n.mat[1,i+1] <- matA[1]
            } else if (wetdry[iter,i] == 1 && eggCollect == F) {
              n.mat[,i+1] <- round(popmat.current %*% n.mat[,i], digits = 0)
            } else if  (wetdry[iter,i] == 0) {
              n.mat[,i+1] <- round(popmat.fail %*% n.mat[,i], digits = 0)
            } else {
              stop("hmm winFail is neither 1 or 0") 
            }
            
            n.mat[1,i+1] <- n.mat[1,i+1] * breedProp
            
            #  density feedback for eggs
            # I  use integration. i.e. early in the curve females will lay with 100% success. As it tends towards the pond limit 
            # successive females lay  with diminishing success 
            # above the egg limit for the pond, laying is possible but with a huge inhibition
            K.rel.egg <- (n.mat[ 1, i+1]/current.K.init[1]) 
            if (!is.nan(K.rel.egg) && (K.rel.egg > 0)) {
              if (K.rel.egg <= 1.02) {  
                intArea <-  integrate(eggfunc,lower=0,upper=K.rel.egg)
                area <- as.numeric(intArea[1])/K.rel.egg
                if (area >= 1) { area <- 1 }
                n.mat[ 1, i+1] <- (round(n.mat[ 1, i+1] * area))
              } else if (K.rel.egg > 1.02) {
                # if  > than the limit for the pond then calculate what happens to the first 102% of the pond limit (egg1) then apply the 102nd%ile inhibition on the remaining eggs(remain)
                egg1 <- (current.K.init[1] * area1)
                remain <- (n.mat[ 1, i+1] - current.K.init[1]) 
                n.mat[ 1, i+1] <- (round(egg1 + (s.mult.egg.iter * remain)))  
              } else { 
                stop("Crashed at line 1894: the egg conversion value K.rel.egg is misbehaving")
              }  
            }
            
            
            n.mat[n.mat < 0] <- 0 
            
            # apply partial survival, prior to tadpole collection (hatch and preCollect.s.iter)
            n.mat[2,i+1] <- round(n.mat[2,i+1] * hatch * preCollect.s.iterDD,digits = 0)
            
            # Apply tadpole collection here #### 
            if (tadCollect == T) {
              pretad <- n.mat[2,i+1]
              n.mat[2,i+1] <- n.mat[2,i+1] - tadCollect.removal
              n.mat[2,i+1] <- (n.mat[2,i+1] - (round(n.mat[2,i+1] * tadCollect.reduction)))
            }
            
            # apply post tadpole collection survival and adult life-stage survival
            n.mat[2,i+1] <- round(n.mat[2,i+1] * postCollect.s.iterDD * toad.s.iterDD, digits = 0)
            
            # apply adult collection here ####
            if (adCollect == T) { 
              if (adCollect.removal > 0) {  n.mat[2:6,i+1] <- adultRemoval(adCollect.removal,n.mat[2:6,i+1]) }
              if (adCollect.reduction > 0) {  n.mat[2:6,i+1] <- adultReduction(adCollect.reduction,n.mat[2:6,i+1]) }
            }
            
            ## save r for this iteration' stochastic matrix (for include relative change in mean instantaneous rate of population change (r))
            r.running <- log(sum(n.mat[2:6,i+1], na.rm=T) / sum(n.mat[2:6,i], na.rm=T))
            r.mat[iter,i] <- ifelse(r.running == -Inf, NA, r.running)
            
            # theoretically unnecessary but included for stability
            for (clean in 1:6) {
              if (is.nan(n.mat[clean,i+1])) {
                n.mat[clean,i+1] <- 0 
              }  } 
            
            #just a doublecheck
            if (i > 85) {
              stop("failed to break at the 85th year")
            }
            ##  Last line of the generation loop (86 years)
          }
          #  last line of the iterations loop
        }
        #last line of incrementing reduc loop 
      }  
      #last line of the stage loop
    }
    # last line of the site loop
  }
  
  dir.nm <- 'G:/My Drive/University Milestones/Sustainable harvest/Outputs/GSA Outputs/'
  # save
  aaa <- minAdpop[1,][minAdpop[1,] < 999999]
  input$meanminP <- sum(aaa)/length(aaa)
  input$minPsum <- sum(aaa)
  bbb <- aaa[aaa > 0]
  input$survminP <- sum(bbb)/length(bbb)
  input$PrExt <- sum(dryExtinct[,1] + wetExtinct[,1] + untracedExtinct[1], NA.rm = T)/100
  save.nm <- paste0('SustHarvGSA',sprintf("%09.0f", rowNum))
  assign(save.nm, input)
  save(list=save.nm,file=paste(dir.nm,save.nm,sep='/'))
  
  print("*******************")
  print(d) 
  print("*******************")
  minAdpop[] <- 999999
  dryExtinct[] <- wetExtinct[] <- untracedExtinct[] <- 0
  
} # end Latin hypercube loop





### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# Global Sensitivity Analysis
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
## parameter ranges
ranges <- list()
ranges$breedProp <- c(0.52, 1)
ranges$clutch.size.lr <- c(2000,5000)
ranges$hatch.dur <- c(1, 5)
ranges$hatch.pr <- c(0.5, 0.98)
ranges$tadpole.dur <- c(50, 90)
ranges$tadpole.mn <- c(0.03, 0.6)
ranges$Hogwash.ad.s.yr <- c(0.1, 0.4)

## create hypercube+
 nSamples <- 10000

lh <- data.frame(randomLHS(n=nSamples, k=length(ranges)))
names(lh) <- names(ranges)

## convert parameters to required scale
for (j in 1:ncol(lh)) {
  par <- names(lh)[j]
  lh[,par] <- qunif(lh[,j], min=ranges[[par]][1], max=ranges[[par]][2]) ## continuous
}



## number of iterations for each parameter set
lh$iter <- 1
dir.nm <- 'G:/My Drive/University Milestones/Sustainable harvest/Outputs/GSA Outputs/'

## uncomment to run in parallel
## Set up parallel processing (nproc is the number of processing cores to use)
# cores <- detectCores()
# nproc <- (cores - 2)
# cl.tmp = makeCluster(rep('localhost', nproc), type='SOCK')
# registerDoSNOW(cl.tmp)
# getDoParWorkers()
# res <- foreach(rowNum=1:nrow(lh),.verbose=T) %do% {sbf_sim(input=lh[rowNum,],dir.nm=dir.nm,rowNum=rowNum)}


lh <- as.matrix(lh)



# or run in series
for (rowNum in 1:nrow(lh)) { sbf_sim(input=lh[rowNum,],dir.nm=dir.nm,rowNum=rowNum) }








#########
## BRT Survival ####
#########
## retrieve results surv
res.nms <- list.files(dir.nm)
res.list <- lapply(res.nms, function(x) {load(paste(dir.nm,x,sep='/')) ; print(x) ; return(eval(as.name(x)))})
datsurv <- data.frame(rbindlist(res.list))
dat.nona <- data.frame(na.omit(datsurv[!is.infinite(rowSums(datsurv)),]))

colnames(dat.nona)

dat.nonaPrExt <- dat.nona[,c(-8,-9,-10)]
head(dat.nonaPrExt)

# the following might need to change
brt.fit <- gbm.step(dat.nonaPrExt, gbm.x = 1:7, gbm.y = 8, family="gaussian", n.trees = 500, max.trees=200000, 
                    tolerance = 0.0001, step.size = 50, learning.rate = 0.0001, bag.fraction=0.6, tree.complexity = 2)






dat.nonaPrExt
summary(brt.fit)
dim(vers.dat)[1]
D2 <- 100 * (brt.fit$cv.statistics$deviance.mean - brt.fit$self.statistics$mean.resid) / brt.fit$cv.statistics$deviance.mean
D2 # % deviance explained
gbm.plot(brt.fit)
gbm.plot.fits(brt.fit)

CV.cor <- 100 * brt.fit$cv.statistics$correlation.mean
CV.cor
CV.cor.se <- 100 *brt.fit$cv.statistics$correlation.se
CV.cor.se
print(c(CV.cor, CV.cor.se))

eq.sp.points <- 100
RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.points, ncol=7)
## output average predictions
for (p in 1:7) {
  RESP.val[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,1]
  RESP.pred[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,2]
}
RESP.val.dat <- as.data.frame(RESP.val)
colnames(RESP.val.dat) <- brt.fit$var.names
RESP.pred.dat <- as.data.frame(RESP.pred)
colnames(RESP.pred.dat) <- brt.fit$var.names
RESP.val.dat
RESP.pred.dat
