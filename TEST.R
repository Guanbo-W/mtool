library(ggplot2)
library(directlabels)
library(spams)
library(Rcpp)
library(microbenchmark)
library(parallel)
library(foreach)
library(doParallel)

source("GenerateData.R")
source("ProxArgument.R")
#sourceCpp("GD.cpp")
source("CV.R")
######################################## IN R ################################################
# nrow = nrow(data)
# covariates=data[,c(6:dim(data)[2])]
# p=dim(data)[2]-6+1
# betas=rep(0,p)
# epsilon = 10^-3
# lambda = 0.1
# t = 3
# alpha = 0.8
# ####### likelihood
# log.lik=function(betas,covariates,data){
#   temp0=NULL
#   #event times TT (j=1,...,m)
#   TT=c(sort(unique(data[data$Event==1,]$Fup)))
#   for (t in TT) {
#     D = data$Id[which(data$Stop==t & data$Event==1)]
#     d = length(D)
#     X.D = as.matrix(covariates[which(data$Stop==t & data$Event==1),],ncol=p)    
#     censor = ifelse(data$Fup==data$Stop & data$Event==0, 1, 0)
#     X.R = as.matrix(covariates[which(data$Stop==t & censor==0),],ncol=p)
#     temp1 = colSums(X.D) %*% betas - d*log(sum(exp(X.R %*% betas)))
#     temp0 = c(temp0,temp1)
#   }
#   return(sum(temp0))
# }
# f.lik=function(betas, covariates, data){return(-t(log.lik(betas, covariates, data))/n )}
# ######## derivative
# derivative.log.lik=function(betas, covariates, data){
#   temp0=rep(NULL,p)
#   TT=c(sort(unique(data[data$Event==1,]$Fup)))
#   for (t in TT) {
#     D = data$Id[which(data$Stop==t & data$Event==1)]
#     d = length(D)
#     X.D = as.matrix(covariates[which(data$Stop==t & data$Event==1),],ncol=p)  
#     censor = ifelse(data$Fup==data$Stop & data$Event==0, 1, 0)
#     X.R = as.matrix(covariates[which(data$Stop==t & censor==0),],ncol=p)
#     temp1 = colSums(X.D) - d*(rowSums( t(X.R)  *  c(exp(X.R %*% betas))  ))/(sum(exp(X.R %*% betas)))
#     temp0 = rbind(temp0,temp1)
#   }
#   return(colSums(temp0))
# }
# f.der=function(betas, covariates, data){return(-t(derivative.log.lik(betas, covariates, data))/n )}
# ######### GD
# V.lambda = function(lambda, betas, covariates, data, t, alpha){
#   u = betas - t(f.der(betas, covariates, data)*t)
#   v = spams.proximalGraph(u, graph, return_val_loss = FALSE, lambda1=lambda*t, numThreads=num_threads, verbose=verbose, pos=pos, intercept=intercept, regul=regul)
#   # lambda = lambda * t ?
#   
#   # until convergence
#   while (sum(abs(betas-v))>epsilon){
#     # backtracking line search, calculate an appropriate step size t
#     while(
#       f.lik(v, covariates, data) >
#       f.lik(betas, covariates, data) + c(f.der(betas, covariates, data)) %*% (v-betas) + sum((v-betas)^2)/(2*t) 
#     )
#     {
#       t = alpha * t
#       u = betas - t(f.der(betas, covariates, data)*t)
#       v = spams.proximalGraph(u, graph, return_val_loss = FALSE, lambda1=lambda*t, numThreads=num_threads, verbose=verbose, pos=pos, intercept=intercept, regul=regul)
#     }
#     
#     # update, for the updated t and fixed lambda
#     betas = v
#     u = betas - t(f.der(betas, covariates, data)*t)
#     v = spams.proximalGraph(u, graph, return_val_loss = FALSE, lambda1=lambda*t, numThreads=num_threads, verbose=verbose, pos=pos, intercept=intercept, regul=regul)
#   }
#   return (v)}
# ######### RESULTS IN R
# f.lik(betas, covariates, data)
# f.der(betas, covariates, data)
# spams.proximalGraph(matrix(rep(1,p),p,1), graph, return_val_loss = FALSE, lambda1=lambda*t, numThreads=num_threads, verbose=verbose, pos=pos, intercept=intercept, regul=regul)
# betas=matrix(rep(0,p),p,1)
# ptm <- proc.time()
# c1=V.lambda(lambda, betas, covariates, data, t, alpha)
# proc.time() - ptm

######################################## IN Rcpp ################################################

p=dim(data)[2]-6+1
nrow = nrow(data)
covariates=matrix(unlist(data[,c(6:dim(data)[2])]),nrow = nrow, ncol = p)
betas = matrix(rep(0,p),nrow=p)
Id = data$Id
Event = data$Event
Fup = data$Fup
Start = data$Start
Stop = data$Stop
lam1 = 0.1
t = 3
epsilon = 10^-3
alpha = 0.8
grp = matrix(c(0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0,
                        1, 1, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0,
                        0, 1, 0, 1, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0),ncol = g,byrow = T)

grpV = matrix(c(1, 0, 0, 0, 0, 0, 0,
                            1, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 1, 0, 0, 0,
                            0, 0, 0, 1, 0, 0, 0,
                            0, 1, 0, 0, 0, 0, 0,
                            0, 0, 1, 0, 0, 0, 0,
                            0, 0, 1, 0, 0, 0, 0,
                            0, 0, 0, 0, 1, 0, 0,
                            0, 0, 0, 0, 1, 0, 0,
                            0, 0, 0, 0, 0, 1, 0,
                            0, 0, 0, 0, 0, 0, 1
),ncol = g,byrow = T)


lik(betas,covariates, Id, Event, Fup, Start, Stop)
der_lik(betas,covariates, Id, Event, Fup, Start, Stop)
proximalGraph(rep(1,p), grp, grpV, eta_g, regul, lam1*t)
# ptm <- proc.time()
# SurvGraphSelect(covariates, Id, Event, Fup, Start, Stop, grp, grpV, eta_g, regul, betas, t, alpha, epsilon, lam1)
# proc.time() - ptm


# cross-validation
# ptm <- proc.time()
# K = 5
# SurvGraphSelect_CV(covariates, Id, Event, Fup, Start, Stop, grp, grpV, eta_g, regul, betas, t, alpha, epsilon, lam1, K)
# proc.time() - ptm

# regularization path plot
lam1max = 0.2
lambda.ratio = 0.001

S_max = log(lam1max)
S_min = log(lam1max*lambda.ratio)
seq_S = seq(S_min, S_max, (S_max-S_min)/50)
seq_lambda = exp(seq_S)
length_lambda =length(seq_lambda)

coefs = NULL
errors = NULL
for(lam1 in seq_lambda){
  coef = SurvGraphSelect(covariates, Id, Event, Fup, Start, Stop, grp, grpV, eta_g, regul, betas, t, alpha, epsilon, lam1)
  coefs = cbind(coefs, coef)
}
coefs.col = matrix(t(coefs),ncol=1)
D1 = data.frame(coefs.col, var=as.character(rep(1:p, each= length_lambda)),lambda=rep(seq_lambda,p))

ggplot(D1, aes(x = lambda, y = coefs.col, color = var, group = var)) + 
  geom_line(position=position_dodge(width=0.01)) + labs(color = 'var')+ xlab(expression(lambda))+ylab("Coefficients")+
  geom_dl(aes(label = var), method = list("first.points"), cex = 0.02)



# CV error plot
K = 5
errors = error_uppers = error_lowers = NULL
# parallel computing
ptm <- proc.time()
cl = makeForkCluster(5)
registerDoParallel(cl)
error_on_each_fold = foreach(lam1 = seq_lambda, .combine=cbind) %dopar% {
  cv = SurvGraphSelect_CV(covariates, Id, Event, Fup, Start, Stop, grp, grpV, eta_g, regul, betas, t, alpha, epsilon, lam1, K)$error_on_each_fold
}
stopCluster(cl)
proc.time() - ptm
errors = apply(error_on_each_fold, 2, mean)
error_uppers = errors + apply(error_on_each_fold, 2, function(x) sd(x))
error_lowers = errors - apply(error_on_each_fold, 2, function(x) sd(x))

# ptm <- proc.time()
# for(lam1 in seq_lambda){
#   cv = SurvGraphSelect_CV(covariates, Id, Event, Fup, Start, Stop, grp, grpV, eta_g, regul, betas, t, alpha, epsilon, lam1, K)
# }
# proc.time() - ptm
# fast 2.5 times


D2 = data.frame(errors, error_uppers, error_lowers, lambda=seq_lambda)

lambda.min = seq_lambda[errors == min(errors)]
min_error.plus_one_se = error_uppers[errors == min(errors)]
lambda.1se = max(seq_lambda[errors<min_error.plus_one_se])

ggplot(D2, aes(x=seq_lambda, y=errors, colour="red")) + theme(legend.position = "none")+
  geom_point()+ xlab(expression(lambda)) + ylab("Mean (1SE) Neg Ave Log Lik")+
  geom_errorbar(aes(ymin=error_lowers, ymax=error_uppers), width=.01, colour="grey", position =position_dodge(0.1))+
  geom_vline(xintercept = lambda.min, linetype="dotted", 
             color = "blue", size=0.5) + geom_text(aes(lambda.min,error_lowers[seq_lambda==lambda.min],label="lambda.min"),size=3, colour="blue") +
  geom_vline(xintercept = lambda.1se, linetype="dotted", 
             color = "darkgreen", size=0.5) + geom_text(aes(lambda.1se,error_lowers[seq_lambda==lambda.min],label="lambda.1se"),size=3, colour="darkgreen")

COEF = coefs[,seq_lambda==lambda.1se]
COEF != 0
which(COEF != 0)
