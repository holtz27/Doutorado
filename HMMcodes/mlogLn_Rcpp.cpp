#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;
using namespace Rcpp;
#include<iostream>
#include<cmath>
#include <time.h>

// [[Rcpp::export]]
double mlogLn_Rcpp(arma::mat allprobs, arma::mat egamma, arma::rowvec foo, int n)
{
  double lscale=0;
  int i=0;
	double sumfoo;

	for(i=1;i<n;i++)
	{		
		foo=foo*egamma%allprobs.row(i);
		sumfoo=sum(foo);
		lscale=lscale+log(sumfoo);
		foo=foo/sumfoo;
	}

	return lscale;
}


