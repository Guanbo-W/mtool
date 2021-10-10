
SurvGraphSelect_CV = function(covariates, Id, Event, Fup, Start, Stop, grp, grpV, eta_g, regul, betas, t, alpha, epsilon, lam1, K){
  error=NULL
  # sample splitting 
  n = length(unique(Id))
  cv.id = NULL
  for (j in 1: K){cv.id[Id %in% seq(n/K*(j-1)+1,n/K*j,1)] = j}
  for (i in 1:K){
    test.covariates = covariates[cv.id == i,]
    test.Id = Id[cv.id == i]
    test.Event = Event[cv.id == i]
    test.Fup = Fup[cv.id == i]
    test.Start = Start[cv.id == i]
    test.Stop = Stop[cv.id == i]
    
    train.covariates = covariates[cv.id != i,]
    train.Id = Id[cv.id != i]
    train.Event = Event[cv.id != i]
    train.Fup = Fup[cv.id != i]
    train.Start = Start[cv.id != i]
    train.Stop = Stop[cv.id != i]
    
    # negative averaged log partial likelihood as an error
    est.train = SurvGraphSelect(train.covariates, train.Id, train.Event, train.Fup, train.Start, train.Stop, grp, grpV, eta_g, regul, betas, t, alpha, epsilon, lam1)
    error[i] = lik(est.train ,test.covariates, test.Id, test.Event, test.Fup, test.Start, test.Stop)
  }
  results=list()
  results$error_on_each_fold = error
  results$error_average = mean(error)
  results$one_se = sd(error)
  results$upper = mean(error) + sd(error)
  results$lower = mean(error) - sd(error)
  return (results)}















