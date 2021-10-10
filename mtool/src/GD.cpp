//#include <RcppArmadillo.h>
//#include <math.h>
//#include <vector>
//#include <Rcpp.h>
//#include <algorithm>
//#include <numeric>
//#include <RcppArmadilloExtensions/sample.h>

#include "RcppArmadillo.h"
#include <iostream>
#include <string>
#include "spams.h"

//using namespace std;
//using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
arma::mat proximalGraph(
    arma::colvec U,
    arma::mat grp,
    arma::mat grpV,
    arma::colvec etaG,
    std::string regul,
    double lam1,
    double lam2 = 0.0,
    double lam3 = 0.0,
    int num_threads = -1,
    bool intercept = false,
    bool resetflow = false,
    bool verbose = false,
    bool pos = false,
    bool clever = true,
    bool eval = true,
    int size_group = 1,
    bool transpose = false) {
  
  // U dimensions
  int p = U.n_rows;
  
  // read in U and convert to spams::Matrix<double> alpha0
  double* ptrU = new double[p];
  for (int r = 0; r < p; r++) {
    ptrU[r] = U(r);
  }
  Matrix<double> alpha0(ptrU,p,1); 
  // alpha0.print("\n alpha0 \n");
  
  // grp dimensions
  int gr = grp.n_rows;
  int gc = grp.n_cols;
  int grc = int (gr * gc);
  
  // read in grp and convert to spams::Matrix<bool> grp_dense
  // then to spams::SpMatrix<bool> groups
  bool* ptrG = new bool[grc];
  for (int r = 0; r < gr; r++) {
    for (int c = 0; c < gc; c++) {
      ptrG[c * gr + r] = (grp(r, c) != 0.0);
    }
  }
  Matrix<bool> grp_dense(ptrG, gr, gc);
  SpMatrix<bool> groups;
  grp_dense.toSparse(groups);
  // groups.print("\n groups: \n");
  
  // grpV dimensions
  int gvr = grpV.n_rows;
  int gvc = grpV.n_cols;
  int gvrc = int (gvr * gvc);
  
  // read in grpV and convert to spams::Matrix<bool> grpV_dense
  // then to spams::SpMatrix<bool> groups_var
  bool* ptrGV = new bool[gvrc];
  for (int r = 0; r < gvr; r++) {
    for (int c = 0; c < gvc; c++) {
      ptrGV[c * gvr + r] = (grpV(r, c) != 0.0);
    }
  }
  Matrix<bool> grpV_dense(ptrGV, gvr, gvc);
  SpMatrix<bool> groups_var;
  grpV_dense.toSparse(groups_var);
  // groups_var.print("\n groups_var: \n");
  
  // read in etaG and convert to spams::Vector<int> eta_g
  int n_etaG = etaG.n_rows;
  
  double* ptrEG = new double[n_etaG];
  for (int i = 0; i < n_etaG; i++) {
    ptrEG[i] = etaG(i);
  }
  
  Vector<double> eta_g(ptrEG, n_etaG);
  // eta_g.print("\n eta_g \n");
  
  // read in regul and convert to char
  int l = regul.length();
  char* name_regul = new char[l];
  for (int i = 0; i < l+1; i++) {
    name_regul[i] = regul[i];
    // std::cout <<"\n char " << name_regul[i];
  }
  
  // Initialize alpha - proximal operator
  Matrix<double> alpha(p,1);
  alpha.setZeros();
  
  // alpha.print("\n alpha: \n");
  
  
  // call _proximalGraph
  _proximalGraph(&alpha0, &alpha,
                 &eta_g, &groups, &groups_var,
                 num_threads, lam1, lam2, lam3,
                 intercept, resetflow, name_regul,
                 verbose, pos, clever, eval,
                 size_group, transpose);
  
  // put updated alpha into arma::mat V to pass to R
  arma::mat V(p, 1);
  for (int r = 0; r < p; r++) {
    V.row(r) = alpha[r];
  }
  return V;
}















// [[Rcpp::export]]
double lik(arma::colvec betas, 
           arma::mat covariates, 
           arma::colvec Id, 
           arma::uvec Event, 
           arma::colvec Fup, 
           arma::colvec Start, 
           arma::uvec Stop){
  arma::colvec uniquId = unique(Id);
  int n = uniquId.n_rows;
  arma::uvec Event_Ind = (Event==1);
  arma::uvec Censor = (Fup==Stop) && (Event==0);
  arma::uvec Censor_Ind = (Censor==0);
  arma::colvec TT =  unique(Fup.rows(find(Event_Ind))); //event times TT (j=1,...,m)
  int n_TT = TT.n_rows;
  arma::colvec temp0(n_TT);
  for(int t = 0; t < n_TT; t++) {
    arma::colvec D = Id(find(Stop == TT(t) && Event_Ind));
    arma::mat X_D = covariates.rows(find(Stop==TT(t) && Event_Ind));
    arma::mat X_R = covariates.rows(find(Stop==TT(t) && Censor_Ind));
    arma::colvec EXP = exp(X_R * betas);
    temp0(t) = as_scalar(sum(X_D,0) * betas - (double) D.n_rows * log(sum(EXP)));
  }
  return -sum(temp0)/n;
  //return TT;
}


// [[Rcpp::export]]
arma::colvec der_lik(arma::colvec betas,
                     arma::mat covariates,
                     arma::colvec Id,
                     arma::uvec Event,
                     arma::colvec Fup,
                     arma::colvec Start,
                     arma::uvec Stop){
  arma::colvec uniquId = unique(Id);
  int n = uniquId.n_rows;
  arma::uvec Event_Ind = (Event==1);
  arma::uvec Censor = (Fup==Stop) && (Event==0);
  arma::uvec Censor_Ind = (Censor==0);
  arma::colvec TT = unique(Fup.rows(find(Event_Ind))); //event times TT (j=1,...,m)
  arma::colvec temp0(covariates.n_cols);
  temp0.zeros();
  for(int t = 0; t < TT.n_rows; t++) {
    arma::colvec D = Id(find(Stop == TT(t) && Event_Ind));
    arma::mat X_D = covariates.rows(find(Stop==TT(t) && Event_Ind));
    arma::mat X_R = covariates.rows(find(Stop==TT(t) && Censor_Ind));
    arma::colvec EXP = exp(X_R * betas);
    arma::colvec temp1(X_R.n_cols);
    temp1.zeros();
    for(int l = 0; l < X_R.n_rows; l++) {
      temp1 += trans(X_R.row(l) * EXP(l));
    }
    temp0 += trans(sum(X_D,0)) - (double) D.n_rows * temp1 / sum(EXP);
  }
  return -temp0/n;
}

// [[Rcpp::export]]
arma::colvec SurvGraphSelect(
    arma::mat covariates, 
    arma::colvec Id, 
    arma::uvec Event, 
    arma::colvec Fup, 
    arma::colvec Start, 
    arma::uvec Stop,
    arma::mat grp,
    arma::mat grpV,
    arma::colvec etaG,
    std::string regul,
    arma::colvec betas, 
    double t,
    double alpha,
    double epsilon,
    double lam1,
    double lam2 = 0.0,
    double lam3 = 0.0,
    int num_threads = -1,
    bool intercept = false,
    bool resetflow = false,
    bool verbose = false,
    bool pos = false,
    bool clever = true,
    bool eval = true,
    int size_group = 1,
    bool transpose = false) {
  arma::colvec U = betas - der_lik(betas, covariates, Id, Event, Fup, Start, Stop)*t;
  arma::mat V_mat = proximalGraph(U , grp , grpV , etaG , regul , lam1*t);
  arma::colvec V = vectorise(V_mat);
  while ( arma::as_scalar(sum(arma::abs(betas-V))) > epsilon ){
    while (lik(V, covariates, Id, Event, Fup, Start, Stop) > lik(betas, covariates, Id, Event, Fup, Start, Stop) + dot (der_lik(betas, covariates, Id, Event, Fup, Start, Stop) , V-betas) + sum(square(V-betas))/(2*t)){
      t = alpha * t;
      U = betas - der_lik(betas, covariates, Id, Event, Fup, Start, Stop)*t;
      V_mat = proximalGraph(U , grp , grpV , etaG , regul , lam1*t);
      arma::colvec V = vectorise(V_mat);
      } 
    betas = V;
    U = betas - der_lik(betas, covariates, Id, Event, Fup, Start, Stop)*t;
    V_mat = proximalGraph(U , grp , grpV , etaG , regul , lam1*t);
    arma::colvec V = vectorise(V_mat);
    } 
  return V;}



/*** R


*/
