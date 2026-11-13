#ifndef _medusa_MEDUSA_h
#define _medusa_MEDUSA_h

#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

arma::mat generate_feature_armadillo(const arma::vec& encoding, int K);
Rcpp::NumericMatrix generate_feature(Rcpp::NumericVector encoding, int K);

arma::mat deconvolve_nlreg(arma::mat BOLDobs, const arma::vec& kernel,
                           double nev_lr, double epsilon, double beta,
                           bool normalize, bool trim_kernel);

arma::vec sigmoid(const arma::vec& x, double beta);
arma::vec dsigmoid(const arma::vec& x, double beta);
double sigmoid(double x, double beta);
double dsigmoid(double x, double beta);

#endif
