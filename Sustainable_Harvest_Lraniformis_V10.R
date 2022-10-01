### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
## Models to guide sustainable wild-harvest strategies for Litoria raniformis  
## model Version 10.0
## Rupert Mathwin
## October 2022
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
## 
## We model three contrasting populations of Litoria raniformis (NapNap Waterhole, Hogwash Bend and Bendigo Water Treatment Plant (herein "Epsom")
## Breeding is binary and predefined by a series of 0s or 1s - wetDry
##

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

# Custom functions  ####
# beta distribution shape parameter estimator function
##Generates an alpha and a beta value to inform the beta distribution 
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta)) }

# call to remove adult frogs by number. 
# Call using - n.mat[2:6,i+1] <- adultRemoval(adCollect.removal,n.mat[2:6,i+1]) 
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

# call to remove adult frogs by proportion. 
# Call using - n.mat[2:6,i+1] <- adultReduction(adCollect.reduction,n.mat[2:6,i+1])
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

# sites and stages
siteList <- c("NapNap","Hogwash","Epsom")
stageList <- c("egg","tad","adult")

# assign iterations ####
iteration <- 10000

## Set up the sites (K values)
NapNap.K.init <- c(141,1417,570,570)
Hogwash.K.init <- c(28,270,108,108) 
Epsom.K.init <- c(40,412,165,165)

# Adult annual survival from Turner et al 2022:
# Epsom=permanent,NapNap=semi-permanent,Hogwash=ephemeral,-(mean,sd)
NapNap.ad.s.yr <- c(0.25, 0.05)
Hogwash.ad.s.yr <- c(0.23,0.08)
Epsom.ad.s.yr <- c(0.28,0.06)

## set time limit for projection in 1-yr increments
yr.now <- 2020
#************************
#yr.end <- 2020 + round((40*gen.l), 0) # set projection end date

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
clutch.size.lr <- c(1885,3893,2448,3090,3191,3644,4563) # L. raniformis
prop.breeding <- c(0,rep(1,5))
fert.mn <- mean(clutch.size.lr)*prop.breeding

#duration data (eggs and tadpoles)
hatch.dur <- c(2,4)
tadpole.dur <- c(70,80) # duration (days) L. raniformis 23 deg

## survival data (not all used)
hatch.pr <- c(0.92575, 0.0644825) # hatch Prob from Christy PhD (L. aurea) - mn, SD
tadpole.mn.1 <- mean(c(.15,.26)) # wild tadpole survival from Bull (C.signifera) *note the beta distribution
tadpole.sd.1 <- ((.26-tadpole.mn.1)+(tadpole.mn.1-.15))/2/1.96
tadpole.mn.2 <- mean(c(.07,.56))
tadpole.sd.2 <- ((.56-tadpole.mn.1)+(tadpole.mn.1-.07))/2/1.96
tadpole.mn <- mean(c(tadpole.mn.1, tadpole.mn.2))
tadpole.sd <- sqrt(tadpole.sd.1^2 + tadpole.sd.2^2)
tp.s.alpha <- estBetaParams(tadpole.mn, tadpole.sd/10)$alpha
tp.s.beta <- estBetaParams(tadpole.mn, tadpole.sd/10)$beta

tomet.dur.iter <- round(sum(c(runif(1, hatch.dur[1], hatch.dur[2]), runif(1, tadpole.dur[1], tadpole.dur[2]))), 0)
toad.dur.iter <- 365 - tomet.dur.iter

tomet.dur.mn <- round(sum(c(mean(hatch.dur), mean(tadpole.dur))), 0)
toad.dur.mn <- 365 - tomet.dur.mn


#define adult survival (NapNap as a holder value)
site <- siteList[1]
ad.s.yr.mn <- get(paste0(site,".ad.s.yr"))[1]
ad.s.yr.sd <- get(paste0(site,".ad.s.yr"))[2]

ad.s.yr.alpha <- estBetaParams(ad.s.yr.mn, ad.s.yr.sd/10)$alpha
ad.s.yr.beta <- estBetaParams(ad.s.yr.mn, ad.s.yr.sd/10)$beta

#likelyhood of egg hatching and tadpole surviving to metamorphosis
hatch.val <- (rnorm(1,hatch.pr[1],hatch.pr[2]))
if (hatch.val > 1) { hatch.val <- 1 }
tomet.s.iter <- (rbeta(1, tp.s.alpha, tp.s.beta)) * hatch.val

#sample a daily probability of survival
toad.daily.s.iter <- nthroot(rbeta(1,ad.s.yr.alpha, ad.s.yr.beta) , 365)

toad.s.season.iter <- toad.daily.s.iter ^ toad.dur.mn
## change the above line to the below line just need the syntax to resample
## toad.s.season.iter <- prod(runif(toad.dur.mn, toad.daily.s.iter))
toad.s.iter <- tomet.s.iter * toad.s.season.iter                     

##Create the survival vector: to adult survival then 4 adult survivals
ad.s.vec.iter <- rep(NA,5)
for (s in 1:5) {
  ad.s.vec.iter[s] <- rbeta(1,ad.s.yr.alpha, ad.s.yr.beta)  }
surv.iter <- c(toad.s.iter, ad.s.vec.iter)

init.vec <- rep(NA,6)
dim(init.vec) <- c(6,1)

##Populate the Matix (popmat) and create a failure matrix(popmat.fail) for years with no breeding
## removed  diag(popmat[2:(stages), ]) <- surv.mn[-stages]
diag(popmat[2:(stages), ]) <- surv.iter[-stages]
popmat[stages,stages] <- 0 # surv.mn[stages] 
popmat[1,] <- fert.mn * sex.ratio
popmat.orig <- popmat ## save original matrix as popmat.orig
popmat <- popmat.orig

#Create init.vec (the starting poulation and age structure)
totalSurv <- sum(popmat[2:6,1:6])
popmat.current <- popmat
popmat.fail <- popmat
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

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
## And NOW for the actual model!  ##################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# Start of the site loop
 for (S in 1:length(siteList)) {
#for (S in 3:3) {
  site <- siteList[S]
  
  #define adult survival using Anna's values
  ad.s.yr.mn <- get(paste0(site,".ad.s.yr"))[1]
  ad.s.yr.sd <- get(paste0(site,".ad.s.yr"))[2]
  ad.s.yr.alpha <- estBetaParams(ad.s.yr.mn, ad.s.yr.sd/10)$alpha
  ad.s.yr.beta <- estBetaParams(ad.s.yr.mn, ad.s.yr.sd/10)$beta
  
  # set the site values for K and the init ####
  current.K.init <- get(paste0(site,".K.init"))
  
  # convert the number of spawning masses to the number of female eggs (uses mean clutch size)
  current.K.init[1] <- current.K.init[1] * (mean(clutch.size.lr) * 0.5)
  
  ## Import the pregenerated datasets ####
  wetdry <- fread(file = paste0("C:/Workspace/SustHarvRevSep22/wetdry",site,".csv"), header = F, sep = ",", dec = ".")
  wetdry <- wetdry[2:10001,]
  wetdry <- array(as.numeric(unlist(wetdry)), dim=c(iteration,99))
  
  startPops <- fread(file = paste0("C:/Workspace/SustHarvRevSep22/startPops",site,"2.csv"), header = T, sep = ",", dec = ".")
  startPops <- array(as.numeric(unlist(startPops)), dim=c(6,iteration))
  
  # create some storage martices
  wetExtinct <- dryExtinct <- matrix(data = 0, nrow =6,ncol = length(reductionTitles))
  colnames(dryExtinct) <- colnames(wetExtinct) <- reductionTitles
  rownames(dryExtinct) <- c("D","DD","DDD","DDDD","DDDDD","DDDDDD")
  rownames(wetExtinct) <- c("W","WW","WWW","WWWW","WWWWW","6+ W")
  untracedExtinct <- rep(0,21)
  
  #start the stage of collection loop, 1-egg,2-tad,3-adult
  for (stage in 1:length(stageList)) {  
    stg <- stageList[stage]
    cat("starting stage loop, stage is ", stage, "\n")
    
    # These variables determine the rate, type, and life stage of collection. retained sperately to allow mutiple approaches
    stageCollect <- c("eggCollect","tadCollect","adCollect")
    eggCollect <- tadCollect <- adCollect <- F
    eggCollect.removal <- eggCollect.reduction <- tadCollect.removal <- tadCollect.reduction <- adCollect.removal <- adCollect.reduction <- 0 
    assign(stageCollect[stage],T)
    
    #start the reduc loop
    #this will go from 0% reduction to 100% reduction in increments of 5%
    # havent used a while loop so I can reference r for columns in the survival array
    for (r in 1:length(reducList)) {
      reduc <- reducList[r]
      cat("stage is ", stage," r is (/21)", r, "\n")
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
              stop("line 433 wetTracker and dryTracker, have I written this correctly? Please double-check")
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
          tohatch.dur.iter <- round(runif(1,2,4),digits = 0)
          tomet.dur.iter <- round(sum(c(runif(1, hatch.dur[1], hatch.dur[2]), runif(1, tadpole.dur[1], tadpole.dur[2]))), 0)
          preCollect.dur.iter <- round(2/3 * tomet.dur.iter,digits = 0)
          postCollect.dur.iter <- tomet.dur.iter - preCollect.dur.iter
          toad.dur.iter <- (365 - tohatch.dur.iter) - tomet.dur.iter
          
          ##Calculate our survivals
          #likelyhood of egg hatching and tadpole surviving to metamorphosis
          hatch.val <- (rnorm(1,hatch.pr[1],hatch.pr[2]))
          if (hatch.val > 1) { hatch.val <- 1 }
          # hatch.s.iter <- runif(1,0.933,1)
          #  tomet.s.iter <- (rbeta(1, tp.s.alpha, tp.s.beta)) * (runif(1, min=0.933, max=1))
          tomet.s.iter <- (rbeta(1, tp.s.alpha, tp.s.beta)) * hatch.val
          tomet.daily.s.iter <- nthroot(tomet.s.iter, tomet.dur.iter)
          preCollect.s.iter <- tomet.daily.s.iter ^ preCollect.dur.iter
          postCollect.s.iter <- tomet.daily.s.iter ^ postCollect.dur.iter
          toad.daily.s.iter <- nthroot(rbeta(1,ad.s.yr.alpha, ad.s.yr.beta) , 365)
          toad.s.iter <- toad.daily.s.iter ^ toad.dur.iter
          
          #Create the survival vector (popmat) for the year (density dependence not considered yet)
          # NOTE 0-1 surv is set to 1 and is applied after the matrix multiplication to allow for Collection of tadpoles
          ad.s.vec.iter <- rep(NA,5)
          for (s in 1:5) {
            ad.s.vec.iter[s] <- rbeta(1,ad.s.yr.alpha, ad.s.yr.beta)  }
          surv.iter <- c(1, ad.s.vec.iter)
          
          fert.iter <- round((runif(stages, min=clutch.size.lr[1], max=clutch.size.lr[2])) * prop.breeding, 0)
          popmat[1,] <- fert.iter * sex.ratio 
          diag(popmat[2:(stages), ]) <- surv.iter[-stages]
          
          #note these are placeholder values for popmat.current, they are recalculated below not necessary but retained for safety
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
          preCollect.s.iter <- preCollect.s.iter * (tad.Ddep.daily^preCollect.dur.iter)
          postCollect.s.iter <- postCollect.s.iter * (tad.Ddep.daily^postCollect.dur.iter)
          toad.s.iter <- toad.s.iter * (tad.Ddep.daily^toad.dur.iter)
          
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
          
          # set popmat.fail for (use )only used if this year is dry, and does not support breeding)
          popmat.fail <- popmat.current
          popmat.fail[1,] <- 0 
          
          if (wetdry[iter,i] == 1) { 
            dryTracker <- 0 
            ifelse(wetTracker < 6, wetTracker <- wetTracker + 1, wetTracker <- 6)
          }
          if (wetdry[iter,i] == 0) {
            ifelse(dryTracker < 6, dryTracker <- dryTracker + 1, dryTracker <- 6)
            wetTracker <- 0
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
          
          
          # apply partial survival, prior to tadpole collection (hatch.val and preCollect.s.iter)
          n.mat[2,i+1] <- round(n.mat[2,i+1] *  hatch.val * preCollect.s.iter,digits = 0)
          
          # Apply tadpole collection here #### 
          if (tadCollect == T) {
            pretad <- n.mat[2,i+1]
            n.mat[2,i+1] <- n.mat[2,i+1] - tadCollect.removal
            n.mat[2,i+1] <- (n.mat[2,i+1] - (round(n.mat[2,i+1] * tadCollect.reduction)))
          }
          
          # apply post tadpole collection survival and adult life-stage survival
          n.mat[2,i+1] <- round(n.mat[2,i+1] * postCollect.s.iter * toad.s.iter, digits = 0)
          
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
          if (i > 85){
            stop("failed to break at the 85th year")
          }
          
          ##  Last line of the generation loop (86 years)
        }
        
        #  last line of the iterations loop
      }
      
      # Calulate the mean r for each iter then calculate the relatives for graphing
      for (i in 1:iteration) {
        rStoch[i] <- mean(r.mat[i,],na.rm = T)
      }
      r.mat[] <- 0
      rStoch[is.infinite(rStoch)]<-NA
      rRange[1,r] <- quantile(rStoch, probs=0.975, na.rm=T)
      rRange[2,r] <- mean(rStoch, na.rm=T)
      rRange[3,r] <- quantile(rStoch, probs=0.025, na.rm=T)
      rStoch[] <- 0
      
      #last line of incrementing reduc loop 
    }  
   
    #change below locations of course : ) 
    write.csv(dryExtinct,paste0("C:/Workspace/SustHarvRevSep22/",site,stg,"dryExta3NN2.csv"), row.names = T)
    write.csv(yrExt,paste0("C:/Workspace/SustHarvRevSep22/",site,stg,"yearExtincta3NN2.csv"), row.names = T)
    write.csv(minAllpop,paste0("C:/Workspace/SustHarvRevSep22/",site,stg,"minAllpopa3NN2.csv"), row.names = T)
    write.csv(minAdpop,paste0("C:/Workspace/SustHarvRevSep22/",site,stg,"minAdPopa3NN2.csv"), row.names = T)
    write.csv(finalAllpop,paste0("C:/Workspace/SustHarvRevSep22/",site,stg,"finalAllpopa3NN2.csv"), row.names = T)
    write.csv(finalAdpop,paste0("C:/Workspace/SustHarvRevSep22/",site,stg,"finalAdpopa3NN2.csv"), row.names = T)
    write.csv(survCheck,paste0("C:/Workspace/SustHarvRevSep22/",site,stg,"survChecka3NN2.csv"), row.names = T)
    write.csv(extCheck,paste0("C:/Workspace/SustHarvRevSep22/",site,stg,"extChecka3NN2.csv"), row.names = T)
    write.csv(untracedExtinct,paste0("C:/Workspace/SustHarvRevSep22/",site,stg,"untracedExtincta3NN2.csv"), row.names = T)
    write.csv(rRange,paste0("C:/Workspace/SustHarvRevSep22/",site,stg,"rstocha3NN2.csv"), row.names = T)
    
    # reset
    untracedExtinct[] <- extCheck[] <- survCheck[] <- dryExtinct[] <- 0
    minAdpop[] <- minAllpop[] <- 999999
    rRange[] <- yrExt[] <- finalAdpop[] <- finalAllpop[] <- NA  
    
    #  finished stage loop  
    #last line of the three stages of reduction (stage) 
  }
  
  # last line of the site loop
}



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Graphing output.csv from Sustainable_Harvest_Lraniformis_V10
# plot finishing is completed in another program
# R Mathwin October 2022
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###



## Remove everything
rm(list = ls())

# Libraries
library(ggplot2)
# to smooth the cdf plot data
library(tidyverse)
# for simple tick marks
library(Hmisc)

siteList <- c("NapNap","Hogwash","Epsom")
# siteList <- c("NapNap","Hogwash")
stgList <- c("egg","tad","adult")
reduction <- c("0%","5%","10%","15%","20%","25%","30%","35%","40%","45%","50%","55%","60%","65%","70%","75%","80%","85%","90%","95%","100%")
remNum <- c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1)

# LOAD ALL OF THE CSVS
for (st in 1:length(siteList)) {
  site <- siteList[st]
  assign(paste0(site,"startPops"),read.csv(paste0("C:/Workspace/SustHarvRevSep22/startPops",site,".csv")))   
  for (test in 1:1) {
    test <- "3NN2"
    for (stage in 1:3) {
      stg <- stgList[stage]
      assign(paste0(site,".dryExt.",stg),read.csv(paste0("C:/Workspace/SustHarvRevSep22/",site,stg,"dryExta",test,".csv"), row.names = 1))
      assign(paste0(site,".yrExt.",stg),read.csv(paste0("C:/Workspace/SustHarvRevSep22/",site,stg,"yearExtincta",test,".csv"), row.names = 1))
      assign(paste0(site,".minAllPop.",stg),read.csv(paste0("C:/Workspace/SustHarvRevSep22/",site,stg,"minAllpopa",test,".csv"), row.names = 1))
      assign(paste0(site,".minAdpop.",stg),read.csv(paste0("C:/Workspace/SustHarvRevSep22/",site,stg,"minAdpopa",test,".csv"), row.names = 1))
      assign(paste0(site,".finalAllPop.",stg),read.csv(paste0("C:/Workspace/SustHarvRevSep22/",site,stg,"finalAllPopa",test,".csv"), row.names = 1))
      assign(paste0(site,".finalAdPop.",stg),read.csv(paste0("C:/Workspace/SustHarvRevSep22/",site,stg,"finalAdPopa",test,".csv"), row.names = 1))
      assign(paste0(site,".survCheck.",stg),read.csv(paste0("C:/Workspace/SustHarvRevSep22/",site,stg,"survChecka",test,".csv"), row.names = 1))
      assign(paste0(site,".extCheck.",stg),read.csv(paste0("C:/Workspace/SustHarvRevSep22/",site,stg,"extChecka",test,".csv"), row.names = 1))
      assign(paste0(site,".untracedExtinct.",stg),read.csv(paste0("C:/Workspace/SustHarvRevSep22/",site,stg,"untracedExtincta",test,".csv"), row.names = 1))
      assign(paste0(site,".rStoch.",stg),read.csv(paste0("C:/Workspace/SustHarvRevSep22/",site,stg,"rstocha",test,".csv"), row.names = 1))
    }
  }  
}

# quick check that the extinction counts add up correctly
PrExtdf <- data.frame(reduction,remNum)
for (j in 1:length(siteList)) {
  site <- siteList[j]
  for (k in 1:length(stgList)) {
    stg <- stgList[k] 
    aa <- get(paste0(site,".extCheck.",stg))
    aa <- aa/10000
    PrExtdf <- cbind(PrExtdf, aa)
    colnames(PrExtdf)[ncol(PrExtdf)] <- paste0(site,stg,"PrExt")
  }
}
PrExtdf[PrExtdf==1] <- 0.995

# plot probability of extinction
# Epsom
plot(PrExtdf$remNum,PrExtdf$EpsomadultPrExt,type = "l",col = "red", xlab = "", ylab = "",
     xaxs = "i",yaxs = "i",axes = F,cex.axis = 1,ylim = c(0, 1),xlim = c(0, 1),lwd = 2.5,main = "Epsom PrExt")
axis(1,cex.axis = 1.5,lwd=2)
axis(2,las = 1,cex.axis = 1.5,lwd=2)
minor.tick(ny = 2, tick.ratio = 0.5)
lines(PrExtdf$remNum,PrExtdf$EpsomeggPrExt, type = "l", 
      col="chocolate4",lwd = 2.5)
lines(PrExtdf$remNum,PrExtdf$EpsomtadPrExt, type = "l", pch = 19, 
      col="blue",lwd = 2.5)

# NapNap
plot(PrExtdf$remNum,PrExtdf$NapNapadultPrExt,type = "l",col = "red", xlab = "", ylab = "",
     xaxs = "i",yaxs = "i",axes = F,cex.axis = 1,ylim = c(0, 1),xlim = c(0, 1),lwd = 2.5,main = "NapNap PrExt")
axis(1,cex.axis = 1.5,lwd=2)
axis(2,las = 1,cex.axis = 1.5,lwd=2)
minor.tick(ny = 2, tick.ratio = 0.5)
lines(PrExtdf$remNum,PrExtdf$NapNapeggPrExt, type = "l", 
      col="chocolate4",lwd = 2.5)
lines(PrExtdf$remNum,PrExtdf$NapNaptadPrExt, type = "l", pch = 19, 
      col="blue",lwd = 2.5)

# Hogwash
plot(PrExtdf$remNum,PrExtdf$HogwashadultPrExt,type = "l",col = "red", xlab = "", ylab = "",
     xaxs = "i",yaxs = "i",axes = F,cex.axis = 1,ylim = c(0, 1),xlim = c(0, 1),lwd = 2.5,main = "Hogwash")
axis(1,cex.axis = 1.5,lwd=2)
axis(2,las = 1,cex.axis = 1.5,lwd=2)
minor.tick(ny = 2, tick.ratio = 0.5)
lines(PrExtdf$remNum,PrExtdf$HogwasheggPrExt, type = "l", 
      col="chocolate4",lwd = 2.5)
lines(PrExtdf$remNum,PrExtdf$HogwashtadPrExt, type = "l", pch = 19, 
      col="blue",lwd = 2.5)

## add a datatable for the W/D/DD
# code is down below just bring it up 
# boxplot of yr extinct
holder <- rep(0,10000)
YrExtBoxdf <- data.frame(holder)
for (j in 1:length(siteList)) {
  site <- siteList[j]
  for (k in 1:length(stgList)) {
    stg <- stgList[k] 
    aa <- get(paste0(site,".yrExt.",stg))
    aa[aa==85] <- NA
    for (i in 1:21) {
      colnm <- remNum[i]
      bb <- as.numeric(aa[i,]) 
      YrExtBoxdf <- cbind(YrExtBoxdf, bb)
      colnames(YrExtBoxdf)[ncol(YrExtBoxdf)] <- colnm
    }
    YrExtBoxdf$holder <- NULL
    boxplot(YrExtBoxdf,data=YrExtBoxdf, na.rm = T, varwidth = T, outline = F,staplewex = T,
            main=paste0(site,stg,"YrExt"), xlab="", ylab="Year extinct")
    YrExtBoxdf <- data.frame(holder)
  }
}

# trying to remove the values < 0.025 % extinct
# boxplot of yr extinct
holder <- rep(0,10000)
holder2 <- rep(NA,10000)
YrExtBoxdf <- data.frame(holder)
for (j in 1:length(siteList)) {
  site <- siteList[j]
  for (k in 1:length(stgList)) {
    stg <- stgList[k] 
    aa <- get(paste0(site,".yrExt.",stg))
    aa[aa==85] <- NA
    for (i in 1:21) {
      colnm <- remNum[i]
      bb <- as.numeric(aa[i,])
      if (sum(is.na(bb)) < 9500) {
        YrExtBoxdf <- cbind(YrExtBoxdf, bb)
        colnames(YrExtBoxdf)[ncol(YrExtBoxdf)] <- colnm
        #might have to add a dummy to replace??
      } else {
        YrExtBoxdf <- cbind(YrExtBoxdf, holder2)
        colnames(YrExtBoxdf)[ncol(YrExtBoxdf)] <- colnm
      }
    }
    YrExtBoxdf$holder <- NULL
    boxplot(YrExtBoxdf,data=YrExtBoxdf, na.rm = T, varwidth = T, outline = F,staplewex = T, ylim = c(0,85),
            
            main=paste0(site,stg,"YrExt"), xlab="", ylab="Year extinct")
    YrExtBoxdf <- data.frame(holder)
  }
}

# this goes with the boxplots to be pasted together in powerpoint
# plot PrExt
# NapNap
plot(PrExtdf$remNum,PrExtdf$NapNapeggPrExt, type = "b", pch = 19
     ,col = "blue", xlab = "proportion removed", ylab = "probability of extinction (%)",main = "NapNap egg PrExt")
plot(PrExtdf$remNum,PrExtdf$NapNaptadPrExt, type = "b", pch = 19
     ,col = "chocolate4", xlab = "proportion removed", ylab = "probability of extinction (%)",main = "NapNap tad PrExt")
plot(PrExtdf$remNum,PrExtdf$NapNapadultPrExt, type = "b", pch = 19
     ,col = "red", xlab = "proportion removed", ylab = "probability of extinction (%)",main = "NapNap adult PrExt")

# Hogwash
plot(PrExtdf$remNum,PrExtdf$HogwasheggPrExt, type = "b", pch = 19
     ,col = "blue", xlab = "proportion removed", ylab = "probability of extinction (%)",main = "Hogwash PrExt egg")
plot(PrExtdf$remNum,PrExtdf$HogwashtadPrExt, type = "b", pch = 19
     ,col = "chocolate4", xlab = "proportion removed", ylab = "probability of extinction (%)",main = "Hogwash PrExt tad")
plot(PrExtdf$remNum,PrExtdf$HogwashadultPrExt, type = "b", pch = 19
     ,col = "red", xlab = "proportion removed", ylab = "probability of extinction (%)",main = "Hogwash PrExt adult")

# Epsom
plot(PrExtdf$remNum,PrExtdf$EpsomeggPrExt, type = "b", pch = 19
     ,col = "blue", xlab = "proportion removed", ylab = "probability of extinction (%)",main = "Epsom PrExt egg")
plot(PrExtdf$remNum,PrExtdf$EpsomtadPrExt, type = "b", pch = 19
     ,col = "chocolate4", xlab = "proportion removed", ylab = "probability of extinction (%)",main = "Epsom PrExt tad")
plot(PrExtdf$remNum,PrExtdf$EpsomadultPrExt, type = "b", pch = 19
     ,col = "red", xlab = "proportion removed", ylab = "probability of extinction (%)",main = "Epsom PrExt adult")

colnames(NapNap.dryExt.egg) <- reduction
colnames(NapNap.dryExt.tad) <- reduction
colnames(NapNap.dryExt.adult) <- reduction
colnames(Hogwash.dryExt.egg) <- reduction
colnames(Hogwash.dryExt.tad) <- reduction
colnames(Hogwash.dryExt.adult) <- reduction
colnames(Epsom.dryExt.egg) <- reduction
colnames(Epsom.dryExt.tad) <- reduction
colnames(Epsom.dryExt.adult) <- reduction

# create a holder and starting adult populations
holder <- Epsomstartadults <- NapNapstartadults <- Hogwashstartadults <- rep(0,10000)
for (i in 1:10000) {
  NapNapstartadults[i] <- sum(NapNapstartPops[2:6,i])
  Hogwashstartadults[i] <- sum(HogwashstartPops[2:6,i])
  Epsomstartadults[i] <- sum(EpsomstartPops[2:6,i])
}

# dry ext risk/site is the % of times it went extinct during that cycle
dryRisk <- cleanExt <- matrix(data=0, nrow = 3, ncol = 3)
stageRisk <- matrix(data=0, nrow = 3, ncol = 3)
colnames(cleanExt) <- siteList

colnames(stageRisk) <- rownames(cleanExt) <- stgList

colnames(dryRisk) <- siteList
rownames(stageRisk) <- rownames(dryRisk) <- c("W", "D","DD")
cleanExtsR <- stageRisk

# site, stage, removal, dryStage
allRisks <- array(data = 0, dim = 3*3*21*3)
dim(allRisks) <- c(3,3,21,3)

# just by site
for (i in 1:length(siteList)) {
  site <- siteList[i]
  for (j in 1:length(stgList)) {
    stg <- stgList[j]
    aa <- get(paste0(site,".extCheck.",stg))[1:21,1]
    bb <- get(paste0(site,".untracedExtinct.",stg))[1:21,1]
    cleanExt[j,i] <- (sum(aa) - sum(bb))
    cc <- get(paste0(site,".dryExt.",stg))[1,1:21]
    dd <- get(paste0(site,".dryExt.",stg))[2,1:21]
    ee <- (sum(cc) + sum(dd))
    dryRisk[2,i] <- ((sum(cc)/cleanExt[j,i])*100)
    dryRisk[3,i] <- ((sum(dd)/cleanExt[j,i])*100)
    dryRisk[1,i] <- 100 - (dryRisk[2,i] + dryRisk[3,i])
  }
}

# create stageRisk just for Hogwash and NapNap (because Epsom is only ever W)
for (i in 1:2) {
  site <- siteList[i]
  for (j in 1:length(stgList)) {
    stg <- stgList[j]
    aa <- sum(get(paste0(site,".untracedExtinct.",stg)))
    title <- paste0(site,".untracedExtinct.",stg)
    # quick check for untraced
    if (aa > 0) { cat("issue 1: ", title," is != 0 \n")  }
    bb <- get(paste0(site,".extCheck.",stg))[1:21,1]
    cc <- colSums(get(paste0(site,".dryExt.",stg))[,1:21])
    cleanExtsR[1,j] <- cleanExtsR[1,j] + sum(bb - cc)  # W
    cleanExtsR[2,j] <- cleanExtsR[2,j] + sum(get(paste0(site,".dryExt.",stg))[1,1:21]) # D
    cleanExtsR[3,j] <- cleanExtsR[3,j] + sum(get(paste0(site,".dryExt.",stg))[2,1:21]) # DD
  }
}

for (j in 1:length(stgList)) {
  stg <- stgList[j]
  for (i in 1:2) {
    site <- siteList[i]
    aa <- get(paste0(site,".extCheck.",stg))[1:21,1]
    bb <- get(paste0(site,".untracedExtinct.",stg))[1:21,1]
  }
  cleanExtsR[j,j] <- (sum(aa) - sum(bb))
  cc <- get(paste0(site,".dryExt.",stg))[1,1:21]
  dd <- get(paste0(site,".dryExt.",stg))[2,1:21]
  ee <- (sum(cc) + sum(dd))
  dryRisk[2,i] <- ((sum(cc)/cleanExt[j,i])*100)
  dryRisk[3,i] <- ((sum(dd)/cleanExt[j,i])*100)
  dryRisk[1,i] <- 100 - (dryRisk[2,i] + dryRisk[3,i])
}
cleanExt #count of extinctions that weren't immediate
dryRisk #times extinct during that wet/Dry/Drydry

#  this is the  one that I use OutsExt
OutsExt<- cleanExtsR
for (i in 1:3) {
  aa <- sum(OutsExt[,i])
  OutsExt[1,i] <- round(((OutsExt[1,i]/aa)*100),digits = 1)
  OutsExt[2,i] <- round(((OutsExt[2,i]/aa)*100),digits = 1) 
  OutsExt[3,i] <- round(((OutsExt[3,i]/aa)*100),digits = 1)
}
OutsExt
colSums(OutsExt)

# rStoch plot (plots seperated by site to enable fine tuning)

rstoch.df <- data.frame(remNum)
vals <- c("up","mn","low")

for (j in 1:length(siteList)) {
  site <- siteList[j] 
  for (i in 1:length(stgList)) {  
    stg <- stgList[i]
    aa <- get(paste0(site,".rStoch.",stg))[1:3,1:21] 
    for (k in 1:length(vals)) {
      rng <- vals[k]
      bb <- aa[k,]
      bb <- as.numeric(bb,decimal = 2)
      
      rstoch.df <- cbind(rstoch.df, bb)
      colnames(rstoch.df)[ncol(rstoch.df)] <- paste0(site,".rStoch.",stg,rng)
    }
  }
}

rstoch.df <- as.data.frame(rstoch.df)
rstoch.df[is.na(rstoch.df)] <- -100
rstoch.df <- rstoch.df[-c(1)]
rstoch.df[rstoch.df > 0.6] <- 0.599
rstoch.df[rstoch.df == -100] <- 100
rstoch.df[rstoch.df < -1.5] <- -1.499
rstoch.df[rstoch.df == 100] <- NA
rstoch.df <- cbind(rstoch.df, remNum)

TOP <-  0.6
BOT <- -1.5

#Hogwashegg
ifelse((min(abs(rstoch.df$Hogwash.rStoch.eggup)) < 0),vert <- (rstoch.df$remNum[which.min(abs(rstoch.df$Hogwash.rStoch.eggup))]) - 0.05,
       vert <- (rstoch.df$remNum[which.min(abs(rstoch.df$Hogwash.rStoch.eggup))]) + 0.05)
vert

HeggR <-  ggplot(data=rstoch.df, aes(x=remNum,y=Hogwash.rStoch.eggup), colour=Flow) + 
  # ggplot(data=NatExt) +
  theme_classic(base_size = 22) +
  theme(axis.line=element_line(color="black",size=1),axis.ticks.length=unit(.35, "lines"),
        axis.ticks = element_line(colour = "black", size = 1)) +
  scale_x_continuous(limits = c(0, 1),breaks=seq(0.1,1,0.1),  expand = c(0, 0)) +
  scale_y_continuous(limits = c(BOT, TOP),breaks=seq(BOT,TOP,0.5), expand = c(0, 0)) +
  ggtitle("Hogwash.rStoch.egg") +
  theme(plot.title = element_text(size=25, hjust=0.5)) +
  labs(y= "", x="") + 
  theme(aspect.ratio=.7) +
  theme(axis.text.x = element_text(size=30)) +
  theme(axis.text.y = element_text(size=30)) +
  geom_ribbon(aes(ymin = Hogwash.rStoch.egglow, ymax = Hogwash.rStoch.eggup), fill = "grey90") +
  geom_hline(yintercept=0, linetype="solid",color="black",size=1) +
  geom_vline(xintercept=vert, linetype="longdash",color="grey55",size=1) +
  geom_line(linetype="solid",color="chocolate4",size=1,aes(x=remNum,y=Hogwash.rStoch.eggmn), show.legend = TRUE) 

HeggR 


#Hogwashtad
ifelse((min(abs(rstoch.df$Hogwash.rStoch.tadup)) < 0),vert <- (rstoch.df$remNum[which.min(abs(rstoch.df$Hogwash.rStoch.tadup))]) - 0.05,
       vert <- (rstoch.df$remNum[which.min(abs(rstoch.df$Hogwash.rStoch.tadup))]) + 0.05)
vert

HtadR <-  ggplot(data=rstoch.df, aes(x=remNum,y=Hogwash.rStoch.tadup), colour=Flow) + 
  # ggplot(data=NatExt) +
  theme_classic(base_size = 22) +
  theme(axis.line=element_line(color="black",size=1),axis.ticks.length=unit(.35, "lines"),
        axis.ticks = element_line(colour = "black", size = 1)) +
  scale_x_continuous(limits = c(0, 1),breaks=seq(0.1,1,0.1),  expand = c(0, 0)) +
  scale_y_continuous(limits = c(BOT, TOP),breaks=seq(BOT,TOP,0.5), expand = c(0, 0)) +
  ggtitle("Hogwash.rStoch.tad") +
  theme(plot.title = element_text(size=25, hjust=0.5)) +
  labs(y= "", x="") + 
  theme(aspect.ratio=.7) +
  theme(axis.text.x = element_text(size=30)) +
  theme(axis.text.y = element_text(size=30)) +
  geom_ribbon(aes(ymin = Hogwash.rStoch.tadlow, ymax = Hogwash.rStoch.tadup), fill = "grey90") +
  geom_hline(yintercept=0, linetype="solid",color="black",size=1) +
  geom_vline(xintercept=vert, linetype="longdash",color="grey55",size=1) +
  geom_line(linetype="solid",color="blue",size=1,aes(x=remNum,y=Hogwash.rStoch.tadmn), show.legend = TRUE) 

HtadR 


#Hogwashadult
ifelse((min(abs(rstoch.df$Hogwash.rStoch.adultup)) < 0),vert <- (rstoch.df$remNum[which.min(abs(rstoch.df$Hogwash.rStoch.adultup))]) - 0.05,
       vert <- (rstoch.df$remNum[which.min(abs(rstoch.df$Hogwash.rStoch.adultup))]) + 0.05)
vert

HadultR <-  ggplot(data=rstoch.df, aes(x=remNum,y=Hogwash.rStoch.adultup), colour=Flow) + 
  # ggplot(data=NatExt) +
  theme_classic(base_size = 22) +
  theme(axis.line=element_line(color="black",size=1),axis.ticks.length=unit(.35, "lines"),
        axis.ticks = element_line(colour = "black", size = 1)) +
  scale_x_continuous(limits = c(0, 1),breaks=seq(0.1,1,0.1),  expand = c(0, 0)) +
  scale_y_continuous(limits = c(BOT, TOP),breaks=seq(BOT,TOP,0.5), expand = c(0, 0)) +
  ggtitle("Hogwash.rStoch.adult") +
  theme(plot.title = element_text(size=25, hjust=0.5)) +
  labs(y= "", x="") + 
  theme(aspect.ratio=.7) +
  theme(axis.text.x = element_text(size=30)) +
  theme(axis.text.y = element_text(size=30)) +
  geom_ribbon(aes(ymin = Hogwash.rStoch.adultlow, ymax = Hogwash.rStoch.adultup), fill = "grey90") +
  geom_hline(yintercept=0, linetype="solid",color="black",size=1) +
  #  geom_vline(xintercept=vert, linetype="longdash",color="blue",size=1) +
  geom_line(linetype="solid",color="red",size=1,aes(x=remNum,y=Hogwash.rStoch.adultmn), show.legend = TRUE) 

HadultR 


#NapNapegg
ifelse((min(abs(rstoch.df$NapNap.rStoch.eggup)) < 0),vert <- (rstoch.df$remNum[which.min(abs(rstoch.df$NapNap.rStoch.eggup))]) - 0.05,
       vert <- (rstoch.df$remNum[which.min(abs(rstoch.df$NapNap.rStoch.eggup))]) + 0.05)
vert

NeggR <-  ggplot(data=rstoch.df, aes(x=remNum,y=NapNap.rStoch.eggup), colour=Flow) + 
  # ggplot(data=NatExt) +
  theme_classic(base_size = 22) +
  theme(axis.line=element_line(color="black",size=1),axis.ticks.length=unit(.35, "lines"),
        axis.ticks = element_line(colour = "black", size = 1)) +
  scale_x_continuous(limits = c(0, 1),breaks=seq(0.1,1,0.1),  expand = c(0, 0)) +
  scale_y_continuous(limits = c(BOT, TOP),breaks=seq(BOT,TOP,0.5), expand = c(0, 0)) +
  ggtitle("NapNap.rStoch.egg") +
  theme(plot.title = element_text(size=25, hjust=0.5)) +
  labs(y= "", x="") + 
  theme(aspect.ratio=.7) +
  theme(axis.text.x = element_text(size=30)) +
  theme(axis.text.y = element_text(size=30)) +
  geom_ribbon(aes(ymin = NapNap.rStoch.egglow, ymax = NapNap.rStoch.eggup), fill = "grey90") +
  geom_hline(yintercept=0, linetype="solid",color="black",size=1) +
  geom_vline(xintercept=0.739, linetype="longdash",color="grey55",size=1) +
  geom_line(linetype="solid",color="chocolate4",size=1,aes(x=remNum,y=NapNap.rStoch.eggmn), show.legend = TRUE) 

NeggR 


#NapNaptad
ifelse((min(abs(rstoch.df$NapNap.rStoch.tadup)) < 0),vert <- (rstoch.df$remNum[which.min(abs(rstoch.df$NapNap.rStoch.tadup))]) - 0.05,
       vert <- (rstoch.df$remNum[which.min(abs(rstoch.df$NapNap.rStoch.tadup))]) + 0.05)
vert

NtadR <-  ggplot(data=rstoch.df, aes(x=remNum,y=NapNap.rStoch.tadup), colour=Flow) + 
  # ggplot(data=NatExt) +
  theme_classic(base_size = 22) +
  theme(axis.line=element_line(color="black",size=1),axis.ticks.length=unit(.35, "lines"),
        axis.ticks = element_line(colour = "black", size = 1)) +
  scale_x_continuous(limits = c(0, 1),breaks=seq(0.1,1,0.1),  expand = c(0, 0)) +
  scale_y_continuous(limits = c(BOT, TOP),breaks=seq(BOT,TOP,0.5), expand = c(0, 0)) +
  ggtitle("NapNap.rStoch.tad") +
  theme(plot.title = element_text(size=25, hjust=0.5)) +
  labs(y= "", x="") + 
  theme(aspect.ratio=.7) +
  theme(axis.text.x = element_text(size=30)) +
  theme(axis.text.y = element_text(size=30)) +
  geom_ribbon(aes(ymin = NapNap.rStoch.tadlow, ymax = NapNap.rStoch.tadup), fill = "grey90") +
  geom_hline(yintercept=0, linetype="solid",color="black",size=1) +
  geom_vline(xintercept=0.755, linetype="longdash",color="grey55",size=1) +
  geom_line(linetype="solid",color="blue",size=1,aes(x=remNum,y=NapNap.rStoch.tadmn), show.legend = TRUE) 

NtadR 


#NapNapadult
ifelse((min(abs(rstoch.df$NapNap.rStoch.adultup)) < 0),vert <- (rstoch.df$remNum[which.min(abs(rstoch.df$NapNap.rStoch.adultup))]) - 0.05,
       vert <- (rstoch.df$remNum[which.min(abs(rstoch.df$NapNap.rStoch.adultup))]) + 0.05)
vert

NadultR <-  ggplot(data=rstoch.df, aes(x=remNum,y=NapNap.rStoch.adultup), colour=Flow) + 
  # ggplot(data=NatExt) +
  theme_classic(base_size = 22) +
  theme(axis.line=element_line(color="black",size=1),axis.ticks.length=unit(.35, "lines"),
        axis.ticks = element_line(colour = "black", size = 1)) +
  scale_x_continuous(limits = c(0, 1),breaks=seq(0.1,1,0.1),  expand = c(0, 0)) +
  scale_y_continuous(limits = c(BOT, TOP),breaks=seq(BOT,TOP,0.5), expand = c(0, 0)) +
  ggtitle("NapNap.rStoch.adult") +
  theme(plot.title = element_text(size=25, hjust=0.5)) +
  labs(y= "", x="") + 
  theme(aspect.ratio=.7) +
  theme(axis.text.x = element_text(size=30)) +
  theme(axis.text.y = element_text(size=30)) +
  geom_ribbon(aes(ymin = NapNap.rStoch.adultlow, ymax = NapNap.rStoch.adultup), fill = "grey90") +
  geom_hline(yintercept=0, linetype="solid",color="black",size=1) +
  geom_vline(xintercept=0.648, linetype="longdash",color="grey55",size=1) +
  geom_line(linetype="solid",color="red",size=1,aes(x=remNum,y=NapNap.rStoch.adultmn), show.legend = TRUE) 

NadultR 

# and for Epsom
ifelse((min(abs(rstoch.df$Epsom.rStoch.eggup)) < 0),vert <- (rstoch.df$remNum[which.min(abs(rstoch.df$Epsom.rStoch.eggup))]) - 0.05,
       vert <- (rstoch.df$remNum[which.min(abs(rstoch.df$Epsom.rStoch.eggup))]) + 0.05)
vert

EeggR <-  ggplot(data=rstoch.df, aes(x=remNum,y=Epsom.rStoch.eggup), colour=Flow) + 
  # ggplot(data=NatExt) +
  theme_classic(base_size = 22) +
  theme(axis.line=element_line(color="black",size=1),axis.ticks.length=unit(.35, "lines"),
        axis.ticks = element_line(colour = "black", size = 1)) +
  scale_x_continuous(limits = c(0, 1),breaks=seq(0.1,1,0.1),  expand = c(0, 0)) +
  scale_y_continuous(limits = c(BOT, TOP),breaks=seq(BOT,TOP,0.5), expand = c(0, 0)) +
  ggtitle("Epsom.rStoch.egg") +
  theme(plot.title = element_text(size=25, hjust=0.5)) +
  labs(y= "", x="") + 
  theme(aspect.ratio=.7) +
  theme(axis.text.x = element_text(size=30)) +
  theme(axis.text.y = element_text(size=30)) +
  geom_ribbon(aes(ymin = Epsom.rStoch.egglow, ymax = Epsom.rStoch.eggup), fill = "grey90") +
  geom_hline(yintercept=0, linetype="solid",color="black",size=1) +
  geom_vline(xintercept=0.82, linetype="longdash",color="grey55",size=1) +
  geom_line(linetype="solid",color="chocolate4",size=1,aes(x=remNum,y=Epsom.rStoch.eggmn), show.legend = TRUE) 

EeggR 


#Epsomtad
ifelse((min(abs(rstoch.df$Epsom.rStoch.tadup)) < 0),vert <- (rstoch.df$remNum[which.min(abs(rstoch.df$Epsom.rStoch.tadup))]) - 0.05,
       vert <- (rstoch.df$remNum[which.min(abs(rstoch.df$Epsom.rStoch.tadup))]) + 0.05)
vert

EtadR <-  ggplot(data=rstoch.df, aes(x=remNum,y=Epsom.rStoch.tadup), colour=Flow) + 
  # ggplot(data=NatExt) +
  theme_classic(base_size = 22) +
  theme(axis.line=element_line(color="black",size=1),axis.ticks.length=unit(.35, "lines"),
        axis.ticks = element_line(colour = "black", size = 1)) +
  scale_x_continuous(limits = c(0, 1),breaks=seq(0.1,1,0.1),  expand = c(0, 0)) +
  scale_y_continuous(limits = c(BOT, TOP),breaks=seq(BOT,TOP,0.5), expand = c(0, 0)) +
  ggtitle("Epsom.rStoch.tad") +
  theme(plot.title = element_text(size=25, hjust=0.5)) +
  labs(y= "", x="") + 
  theme(aspect.ratio=.7) +
  theme(axis.text.x = element_text(size=30)) +
  theme(axis.text.y = element_text(size=30)) +
  geom_ribbon(aes(ymin = Epsom.rStoch.tadlow, ymax = Epsom.rStoch.tadup), fill = "grey90") +
  geom_hline(yintercept=0, linetype="solid",color="black",size=1) +
  geom_vline(xintercept=vert, linetype="longdash",color="grey55",size=1) +
  geom_line(linetype="solid",color="blue",size=1,aes(x=remNum,y=Epsom.rStoch.tadmn), show.legend = TRUE) 

EtadR 


#Epsomadult
ifelse((min(abs(rstoch.df$Epsom.rStoch.adultup)) < 0),vert <- (rstoch.df$remNum[which.min(abs(rstoch.df$Epsom.rStoch.adultup))]) - 0.05,
       vert <- (rstoch.df$remNum[which.min(abs(rstoch.df$Epsom.rStoch.adultup))]) + 0.05)
vert

EadultR <-  ggplot(data=rstoch.df, aes(x=remNum,y=Epsom.rStoch.adultup), colour=Flow) + 
  # ggplot(data=NatExt) +
  theme_classic(base_size = 22) +
  theme(axis.line=element_line(color="black",size=1),axis.ticks.length=unit(.35, "lines"),
        axis.ticks = element_line(colour = "black", size = 1)) +
  scale_x_continuous(limits = c(0, 1),breaks=seq(0.1,1,0.1),  expand = c(0, 0)) +
  scale_y_continuous(limits = c(BOT, TOP),breaks=seq(BOT,TOP,0.5), expand = c(0, 0)) +
  ggtitle("Epsom.rStoch.adult") +
  theme(plot.title = element_text(size=25, hjust=0.5)) +
  labs(y= "", x="") + 
  theme(aspect.ratio=.7) +
  theme(axis.text.x = element_text(size=30)) +
  theme(axis.text.y = element_text(size=30)) +
  geom_ribbon(aes(ymin = Epsom.rStoch.adultlow, ymax = Epsom.rStoch.adultup), fill = "grey90") +
  geom_hline(yintercept=0, linetype="solid",color="black",size=1) +
  geom_vline(xintercept=0.76, linetype="longdash",color="grey55",size=1) +
  geom_line(linetype="solid",color="red",size=1,aes(x=remNum,y=Epsom.rStoch.adultmn), show.legend = TRUE) 

EadultR 


### ### 
# smoothed ecdf plots (0,25,50 and 70% removal)
# have smoothed the outputs because the Hogwash ones are unreadable otherwise
# incs are the increments 1=no removal,6=25%,11=50%,16=75%
incCol <- c("blue","chocolate4","red")

#  manually rejigging each of the three sites 
# to reuse change j value form 1:1,2:2 and 3:3
# change 'adjust' value to change the "smoothness"
for (j in 3:3) {
  site <- siteList[j]
  for (k in 1:length(stgList)) {
    stg <- stgList[k] 
    group <- append(append(append(rep(1,10000),rep(6,10000)),rep(11,10000)),rep(16,10000))
    x <- append(append(append(as.numeric(get(paste0(site,".minAdpop.",stg))[1,]),as.numeric(get(paste0(site,".minAdpop.",stg))[6,])),
                       as.numeric(get(paste0(site,".minAdpop.",stg))[11,])),as.numeric(get(paste0(site,".minAdpop.",stg))[16,]))
    dat <- data.frame(x,group)
    # replace the holder values for extinction with minPop os 0
    dat$x[dat$x==999999] <- 0 # now dat is ready to go
    # Split the data by group and calculate the smoothed cumulative density for each group
    dens = split(dat, dat$group) %>% 
      map_df(function(d) {
        dens = density(d$x, adjust=2, from=min(dat$x) - 0.05*diff(range(dat$x)), 
                       to=max(dat$x) + 0.05*diff(range(dat$x)))
        data.frame(x=dens$x, y=dens$y, cd=cumsum(dens$y)/sum(dens$y), group=d$group[1])
      })
    dens$x[dens$x < 0] <- 0
    dens$cd[dens$cd  > 0.995] <- 0.995
    #subset to each removal proportion to help me graph
    x0 <- subset(dens, group == "1", select = c("x","cd"))
    x25 <- subset(dens, group == "6", select = c("x","cd"))
    x50 <- subset(dens, group == "11", select = c("x","cd"))    
    x75 <- subset(dens, group == "16", select = c("x","cd")) 
    plot(x0$x,x0$cd, type = "l",col=incCol[k] , xlab = "", ylab = "",
         xaxs = "i",xlim = c(0, 130),axes = F, yaxs = "i",ylim = c(0, 1),
         lwd = 2.5,main = paste0("ECDF of minimum adult population at ",site," removing ",stg,"s"))
    axis(1,at = seq(20,120, by = 20),cex.axis = 1.5, lwd=2)
    axis(2,las = 1,cex.axis = 1.5,lwd=2)
    #minor.tick(nx = 2, ny = 0, tick.ratio = 1)
    lines(x25$x,x25$cd, type = "l", col=incCol[k], lwd = 2.5)
    lines(x50$x,x50$cd, type = "l", pch = "", col=incCol[k] ,lwd = 2.5)
    lines(x75$x,x75$cd, type = "l", pch = "", col=incCol[k] ,lwd = 2.5)
    abline(h=0, col="black",lwd = 2.5)
  }
}

# Thanks for reading : ) R



