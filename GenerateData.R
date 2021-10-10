############ Guanbo Ph.D. Project 1: Overlapping Structural Penalization in Cox Models with Time-Dependent Covariates ##############

# The simulation below fits the data generating mechanism where only one trt is analyzed
library(survival)
library(PermAlgo)
library(dplyr)
library(spams)
library(ggplot2)



# generate trt process, m is the max follow-up time
trtonoff = function(m){
  trt.disc.time = round(runif(1,7,m),0) # at least on for a week
  trt.disc.duration = round(runif(1,7,m),0) # at least off for a week
  vec = c(rep(1, trt.disc.time), rep(0, trt.disc.duration))
  while (length(vec) <= m){
    trt.disc.time = round(runif(1,7,m),0) # at least on for a week
    trt.disc.duration = round(runif(1,7,m),0) # at least off for a week
    vec = append(vec, c(rep(1, trt.disc.time), rep(0, trt.disc.duration)))}
  return(vec[1:m])
}

# generate dose and adherence processes according to the trt process for each patient: two possible dose levels when trt is on; adherence is calculated as MPR, and categorized
generate.ind.dose.adherence = function(m, adherence.period, adherence.cutoff){
  ### DOSE
  dose.level = NULL
  trtonoff.0 = trtonoff(m)
  dose.level[which(trtonoff.0 == 0)]=0
  trtonoff.1 = c(1, head(trtonoff.0, -1))
  thelast1 = which(trtonoff.1 - trtonoff.0 == 1) - 1
  if (tail(trtonoff.0, 1) == 1) {thelast1 = c(thelast1, m)} 
  thefirst1 = c(1, which(trtonoff.1-trtonoff.0 == -1))
  trt.on.period = thelast1 - thefirst1+1
  
  k=1
  # for each trt on period
  while (k <= length(trt.on.period)){
  # generate duration of each dose
  duration = NULL; i = 1
  while ( sum(duration) <= trt.on.period[k] ){
    duration[i] = round(runif(1,1,trt.on.period[k]),0)
    i=i+1
  }
  dose.durations = c(head(duration,-1),trt.on.period[k]-sum(head(duration,-1)))
  
  start.low.or.high = round(runif(1,1,2),0)
  vec=NULL; i = 1
  while (i <= length(dose.durations)) {
    if (i %% 2 == 1){
      temp=rep(start.low.or.high,dose.durations[i])
    } else{
      temp = rep(3-start.low.or.high,dose.durations[i])
    }
    vec = c(vec,temp)
    i = i+1
  }
  dose.level[thefirst1[k]:thelast1[k]] = vec
  k = k+1
  }
  
  # generate dummy variables of dose
  dose.high = ifelse(dose.level == 2, 1, 0)
  dose.low = ifelse(dose.level ==1, 1, 0)
  
  ### ADHERENCE (more like mean possession ratio (MPR) in "Anticoagulants in Older Patients with Non-valvular Atrial Fibrillation after Intracranial Hemorrhage")
  # generate continuous adherence
  adherence.conti=NULL
  for (j in 1:m){
    if (j <= adherence.period) {adherence.conti[j] = mean(trtonoff.0[1:j])} else {adherence.conti[j] = mean(trtonoff.0[(j-adherence.period+1):j]) }
  }
  #generate binary adherence
  adherence.high = ifelse(adherence.conti < adherence.cutoff, 0, 1)
  return(data.frame(#trtonoff.0, dose.level, 
                    "D1" = dose.low, "D2" = dose.high, 
                    #adherence.conti, adherence.high
                    "I1" = dose.low * adherence.high, "I2" = dose.high * adherence.high))
}  

# generate n patients' dose and adherence processes ( (n*m) * p )
generate.DI = function(n, m, adherence.period, adherence.cutoff){
  X = list(); Y=NULL
  for (ID in 1:n) {
    X[[ID]] = generate.ind.dose.adherence (m, adherence.period, adherence.cutoff)
    Y=rbind(Y,cbind("ID" = rep(ID, m), as.data.frame(X[[ID]])))  }
  return(Y) 
  }

# generate.trt(m, adherence.period, adherence.cutoff)

# generate the other trt use for each individual
generate.ind.the.other.trt <- function(m){
  start <- round(runif(1,1,m/2),0) # individual start date, assume the trt must start before m/2
  duration <- 7 + 7*rpois(1,3) # assume the trt was taken by week, and must be taken for at least a week
  vec <- c(rep(0, start-1), rep(1, duration))
  while (length(vec) <= m){
    intermission <- 7 + 7*rpois(1,3) # at least off for a week
    duration <- 7 + 7*rpois(1,3) 
    vec <- append(vec, c(rep(0, intermission), rep(1, duration)))}
  return(vec[1:m])}

# generate n patients teh other trt use
generate.B = function(n, m){
  Y=NULL
  for (ID in 1:n) {
    X = generate.ind.the.other.trt (m)
    Y = c(Y,X) }
  return(Y) 
}

generate.cov = function(n, m, adherence.period, adherence.cutoff) {
  # time-independent covariates
  X1 = c(rep(rnorm(n), each=m))
  X2 = c(rep(rnorm(n), each=m))
  DDII = generate.DI(n, m, adherence.period, adherence.cutoff)
  B = generate.B(n, m)
  D1B = DDII$D1 * B
  D2B= DDII$D2 * B
  I1B = DDII$I1 * B
  I2B= DDII$I2 * B
  COV = data.frame(DDII, "B"=B, D1B, D2B, I1B, I2B, X1, X2) 
}

#set.seed(19842)
n=1000 # subjects
m=365 # days
adherence.cutoff = 0.8
adherence.period = 30

Xmat = as.matrix(generate.cov(n, m, adherence.period, adherence.cutoff))


# generate vectors of event and censoring times prior to calling the
# function for the algorithm
eventRandom <- round(rexp(n, 0.005)+14,0)
censorRandom <- round(runif(n, 50,700),0)

data <- permalgorithm(n, m, Xmat[,-1], 
                       XmatNames=c("D1", "D2", "I1","I2", "B","D1B","D2B","I1B","I2B", "X1", "X2"),
                      eventRandom = eventRandom, censorRandom=censorRandom, 
                      betas=c(
                        # D1, D2
                        log(2000), log(1900), 
                        # I1, I2
                        0, 0,
                        # B
                        0,
                        # D1B, D2B, I1B, I2B,
                        0, 0, 0, 0,
                        # X1, X2
                        log(400), 0),
                              groupByD=FALSE )
fit=coxph(Surv(Start, Stop, Event) ~ D1+D2+I1+I2+B+D1B+D2B+I1B+I2B+X1+X2,data)
fit









