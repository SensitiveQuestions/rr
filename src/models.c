#include <string.h>
#include <stddef.h>
#include <stdio.h>      
#include <math.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <R.h>
#include "vector.h"
#include "subroutines.h"
#include "rand.h"


/**
   A Random Walk Metroplis Sampler for Binomial Logistic Regression 
   with Independent Normal Prior
   
   proposal distribution is the univariate normal whose mean is
   the current value and variance is given by the input. each
   parameter is updated one by one.
**/

void rrLogit(int *Y,        /* outcome variable: 0, 1 */
	     double **X,    /* (N x K) covariate matrix */
	     double p,     /* probability of answering truthfully */
	     double p1,    /* probability of answering Y = 1 */
	     double *beta,  /* K coefficients */
	     int n_samp,    /* # of obs */
	     int n_cov,     /* # of covariates, K */
	     double *beta0, /* K prior mean vector */
	     double **A0,   /* (K x K) prior precision */
	     double **Var, /* K proposal precision */
	     int n_gen,     /* # of MCMC draws */
	     int *counter   /* # of acceptance */
	     ) {
  
  int i, j, main_loop;
  double numer, denom, Xbeta, Xprop;
  double *prop = doubleArray(n_cov);

  for (main_loop = 0; main_loop < n_gen; main_loop++) {

    /** Sample from the proposal distribution **/
    rMVN(prop, beta, Var, n_cov);
    
    /** Calculating the ratio (log scale) **/
    /* prior */
    numer = dMVN(prop, beta0, A0, n_cov, 1);
    denom = dMVN(beta, beta0, A0, n_cov, 1);   
    
    /* likelihood */
    for (i = 0; i < n_samp; i++) {
      Xbeta = 0;
      Xprop = 0;
      for (j = 0; j < n_cov; j++) {
	Xprop += X[i][j]*prop[j];
	Xbeta += X[i][j]*beta[j];
      }
      if (Y[i] == 1) {
	denom += log(p / (1 + exp(-Xbeta)) + p1);
	numer += log(p / (1 + exp(-Xprop)) + p1);
      } else {
	denom += log(1 - p / (1 + exp(-Xbeta)) - p1);
	numer += log(1 - p / (1 + exp(-Xprop)) - p1);
      }
    }
      
    /** Rejection **/
    if (unif_rand() < fmin2(1.0, exp(numer-denom))) {
      counter[0]++;
      for (j = 0; j < n_cov; j++) {
	beta[j] = prop[j];
      }
    }
    
    
  }
  
  free(prop);
} /* end of rrLogit */


/**
   A Random Walk Metroplis Sampler for Binomial Logistic Mixed Effects 
   Regression with Independent Normal Prior and Normal random effects.
   
   proposal distribution for fixed effects is the normal whose mean is
   the current value and variance is given by the input. 

   proposal distribution for random effects is the multivariate normal
   whose mean is the current value and variance is given by the
   current value of covariance matrix multiplied by the input tuning
   parameter. 

**/

void rrLogitMixed(int *Y,          /* outcome variable: 0, 1 */
		  double **X,      /* (N x K) covariate matrix for
				      fixed effects */
		  double ***Z,     /* covariates for random effects 
				      organized by groups */
		  double p,        /* probability of answering
				      truthfully */
		  double p1,       /* probability of answering Y = 1 */
		  int *grp,        /* group indicator, 0, 1, ..., G-1 */
		  double *beta,    /* K coefficients for fixed effects */
		  double **gamma,  /* (G x L) matrix of random effects */
		  double **Psi,    /* LxL precision matrix for random effecs */
		  int n_samp,      /* # of obs */
		  int n_fixed,     /* # of fixed effects, K */
		  int n_random,    /* # of random effects, L */
		  int n_grp,       /* # of groups, G */
		  double *beta0,   /* K dimensional prior mean vector */
		  double **A0,     /* (K x K) prior precision */
		  int tau0,        /* prior df for Psi */
		  double **T0,     /* prior scale for Psi */
		  double **tune_fixed,  /* K proposal variance-covariance matrix */
		  double *tune_random, /* tuning constant for random effects of each group */
		  int n_gen,        /* # of MCMC draws */
		  int *acc_fixed,   /* # of acceptance for fixed effects */
		  int *acc_random   /* # of acceptance for random effects */
		  ) {
  
  int i, j, k, main_loop;
  int *vitemp = intArray(n_grp);
  double numer, denom;
  /* proposal values */
  double *beta1 = doubleArray(n_fixed);
  double *gamma1 = doubleArray(n_random);
  /* prior for gamma = 0 */
  double *gamma0 = doubleArray(n_random);
  /* data holders */
  double *Xbeta = doubleArray(n_samp);
  double *Xbeta1 = doubleArray(n_samp);
  double *Zgamma = doubleArray(n_samp);
  double *Zgamma1 = doubleArray(n_samp);
  /* matrix holders */
  double **mtemp = doubleMatrix(n_random, n_random);
  double **mtemp1 = doubleMatrix(n_random, n_random);

  for (j = 0; j < n_fixed; j++)
    beta1[j] = beta[j];

  for (j = 0; j < n_random; j++)
    gamma0[j] = 0;

  /** initializing Xbeta and Zgamma **/
  for (j = 0 ; j < n_grp; j++)
    vitemp[j] = 0;
  for (i = 0; i < n_samp; i++) {
    Xbeta[i] = 0; Zgamma[i] = 0;
    for (j = 0; j < n_fixed; j++) { 
      Xbeta[i] += X[i][j] * beta[j];
    }
    Xbeta1[i] = Xbeta[i];
    for (j = 0; j < n_random; j++) {
      /* Rprintf("%5g\n", Z[grp[i]][vitemp[grp[i]]][j]); */
      Zgamma[i] += Z[grp[i]][vitemp[grp[i]]][j]*gamma[grp[i]][j];
    }
    vitemp[grp[i]]++;
  }
  /* PintArray(grp, n_samp); */
  /* PdoubleArray(Zgamma, n_samp); */
  /* PdoubleMatrix(gamma, n_grp, n_random); */

  /** MCMC Sampler starts here **/
  for (main_loop = 0; main_loop < n_gen; main_loop++) {

    /** STEP 1: Update Fixed Effects **/
    rMVN(beta1, beta, tune_fixed, n_fixed);
    /* Rprintf("beta:");
    PdoubleArray(beta, n_fixed);
    Rprintf("beta1:");
    PdoubleArray(beta1, n_fixed); */
    numer = dMVN(beta1, beta0, A0, n_fixed, 1);
    denom = dMVN(beta, beta0, A0, n_fixed, 1);   
    /* Rprintf("numer: %5g\n", numer);
       Rprintf("denom: %5g\n", denom); */
    for (i = 0; i < n_samp; i++) {
      Xbeta1[i] = 0;
      for (j = 0; j < n_fixed; j++) {
	Xbeta1[i] += X[i][j] * beta1[j];
      }
      /* Rprintf("Zgamma: %5g\n", Zgamma[i]);*/
      if (Y[i] == 1) {
	denom += log(p / (1 + exp(-Xbeta[i]-Zgamma[i])) + p1);
	numer += log(p / (1 + exp(-Xbeta1[i]-Zgamma[i])) + p1);
      } else {
	denom += log(1 - p / (1 + exp(-Xbeta[i]-Zgamma[i])) - p1);
	numer += log(1 - p / (1 + exp(-Xbeta1[i]-Zgamma[i])) - p1);
      }
    }
    /* Rprintf("numer: %5g\n", numer);
       Rprintf("denom: %5g\n", denom); */
    if (unif_rand() < fmin2(1.0, exp(numer-denom))) {
      acc_fixed[0]++;
      for (j = 0; j < n_fixed; j++)
	beta[j] = beta1[j];
      for (i = 0; i < n_samp; i++) 
	Xbeta[i] = Xbeta1[i];
    } 
    /** STEP 2: Update Random Effects Given Fixed Effects **/
    for (j = 0; j < n_grp; j++) {
      for (i = 0; i < n_random; i++)
	for (k = 0; k < n_random; k++) {
	  if (i == k) {
	    mtemp[i][i] = tune_random[j];
	  } else {
	    mtemp[i][k] = 0;
	  }
	}
      rMVN(gamma1, gamma[j], mtemp, n_random);
      /* Calculating the ratio (log scale) */
      /* prior */
      numer = dMVN(gamma1, gamma0, Psi, n_random, 1);
      denom = dMVN(gamma[j], gamma0, Psi, n_random, 1); 
      /* likelihood for group j */
      for (k = 0 ; k < n_grp; k++)
	vitemp[k] = 0;
      for (i = 0; i < n_samp; i++) {
	if (grp[i] == j) {
	  Zgamma1[i] = 0;
	  for (k = 0; k < n_random; k++)
	    Zgamma1[i] += Z[j][vitemp[j]][k]*gamma1[k];
	  if (Y[i] == 1) {
	    denom += log(p / (1 + exp(-Xbeta[i]-Zgamma[i])) + p1);
	    numer += log(p / (1 + exp(-Xbeta[i]-Zgamma1[i])) + p1);
	  } else {
	    denom += log(1 - p / (1 + exp(-Xbeta[i]-Zgamma[i])) - p1);
	    numer += log(1 - p / (1 + exp(-Xbeta[i]-Zgamma1[i])) - p1);
	  }
	}
	vitemp[grp[i]]++;
      }
      if (unif_rand() < fmin2(1.0, exp(numer-denom))) {
	acc_random[j]++;
	for (k = 0; k < n_random; k++)
	  gamma[j][k] = gamma1[k];
	for (i = 0; i < n_samp; i++) {
	  if (grp[i] == j) {
	    Zgamma[i] = Zgamma1[i];
	  }      
	}
      }
    }
    
    /** STEP 3: Update Psi **/
    for (j = 0; j < n_random; j++)
      for (k = 0; k < n_random; k++)
	mtemp[j][k] = T0[j][k];
    for (i = 0; i < n_grp; i++)
      for (j = 0; j < n_random; j++)
	for (k = 0; k < n_random; k++)
	  mtemp[j][k] += gamma[i][j] * gamma[i][k];
    dinv(mtemp, n_random, mtemp1);
    rWish(Psi, mtemp1, (tau0+n_grp), n_random);
  }

  /* freeing memory */
  free(beta1);
  free(gamma1);
  free(gamma0);
  free(vitemp);
  free(Xbeta);
  free(Xbeta1);
  free(Zgamma);
  free(Zgamma1);
  FreeMatrix(mtemp, n_random);
  FreeMatrix(mtemp1, n_random);
} /* end of mixed effects logit */


