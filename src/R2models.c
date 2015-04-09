#include <string.h>
#include <stdio.h>      
#include <stdlib.h>
#include <math.h>
#include <Rmath.h>
#include <R.h>
#include "vector.h"
#include "subroutines.h"
#include "rand.h"
#include "models.h"


void R2rrLogit(int *Y,        /* outcome variable: 0, 1 */
	       double *dX,    /* (N x K) covariate matrix */
	       double *beta,  /* (K(J-1)) stacked coefficient vector */
	       double *p,     /* probability of answering truthfully */
	       double *p1,    /* probability of answering Y = 1 */
	       int *n_samp,   /* # of obs */
	       int *n_cov,    /* # of covariates, K */
	       double *beta0, /* (K(J-1)) prior mean vector */
	       double *dA0,   /* (K(J-1) x K(J-1)) prior precision */
	       double *Var,   /* K(J-1) proposal variances */
	       int *n_gen,     /* # of MCMC draws */
	       int *counter,  /* # of acceptance for each parameter */
	       int *verbose, /* want to print progress? */
	       double *store  /* Storage for beta */
	       ) {

  /* storage parameters and loop counters */
  int i, j, k, itemp, main_loop;  
  
  /* matrices */
  double **X = doubleMatrix(*n_samp, *n_cov+1);
  double **A0 = doubleMatrix(n_cov[0], n_cov[0]);
  double **V = doubleMatrix(n_cov[0], n_cov[0]);

  /* get random seed */
  GetRNGstate();

  /* packing the data */
  itemp = 0;
  for (j = 0; j < *n_cov; j++)
    for (i = 0; i < *n_samp; i++) 
      X[i][j] = dX[itemp++];

  /* packing the prior */
  itemp = 0; 
  for (k = 0; k < n_cov[0]; k++)
    for (j = 0; j < n_cov[0]; j++)
      A0[j][k] = dA0[itemp++];

  itemp = 0; 
  for (k = 0; k < n_cov[0]; k++)
    for (j = 0; j < n_cov[0]; j++)
      V[j][k] = Var[itemp++];

  counter[0] = 0;
  /* Gibbs Sampler! */
  itemp = 0;
  for(main_loop = 1; main_loop <= *n_gen; main_loop++) {

    rrLogit(Y, X, *p, *p1, beta, *n_samp, *n_cov, beta0, A0,
	    V, 1, counter);

    /* Storing the output */
    for (j = 0; j < *n_cov; j++)
      store[itemp++] = beta[j];

    if (*verbose) {
       Rprintf("acceptance rate: %5g\n", ((double) *counter / (double)
			    	       main_loop)); 
    }
    R_FlushConsole(); 
    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */

  PutRNGstate();

  /* freeing memory */
  FreeMatrix(X, *n_samp);
  FreeMatrix(A0, *n_cov);
  FreeMatrix(V, *n_cov);
}




void R2rrLogitMixed(int *Y,        /* outcome variable: 0, 1, ..., J-1 */
		    double *dX,    /* (N x K) covariate matrix for
				      fixed effects */
		    double *dZ,    /* covariates for random effects
				      organized by groups */
		    double *p,     /* probability of answering
				      truthfully */
		    double *p1,    /* probability of answering Y = 1 */
		    int *grp,      /* group indicator, 0, 1, ...,
					   G-1 */
		    double *beta,  /* (K(J-1)) stacked coefficient
					   vector for fixed effects */
		    double *dPsi,  /* LxL precision matrix for
				      random effecs for each equation */
		    int *n_samp,       /* # of obs */
		    int *n_fixed,      /* # of fixed effects, K */
		    int *n_random,     /* # of random effects, L */
		    int *n_grp,        /* # of groups, G */
		    int *max_samp_grp, /* max # of obs within each
					  group */
		    double *beta0,    /* (K(J-1)) prior mean vector */
		    double *dA0,      /* (K(J-1) x K(J-1)) prior
					 precision */
		    int *tau0,        /* prior df for Psi */
		    double *dT0,      /* prior scale for Psi */
		    double *tune_fixed,  /* K(J-1) proposal variances */
		    double *tune_random, /* tuning constant for random
					    effects of each random effect */
		    int *n_gen,        /* # of MCMC draws */
		    int *acc_fixed,    /* # of acceptance for fixed effects */
		    int *acc_random,   /* # of acceptance for random
					  effects */
		    int *verbose,  /* want to print progress? */
		    double *betaStore,
		    double *gammaStore,
		    double *PsiStore
		    ) {

   /* storage parameters and loop counters */
  int i, j, k, main_loop, itemp;  
  int *vitemp = intArray(*n_grp);
  int ibeta = 0, iPsi = 0, igamma = 0;

  /* matrices */
  double *gamma0 = doubleArray(*n_random);
  double **X = doubleMatrix(*n_samp, *n_fixed);
  double **gamma = doubleMatrix(*n_grp, *n_random);
  double **Psi = doubleMatrix(*n_random, *n_random);
  double **PsiInv = doubleMatrix(*n_random, *n_random);
  double **A0 = doubleMatrix(*n_fixed, *n_fixed);
  double **T0 = doubleMatrix(*n_random, *n_random);
  double **Tune = doubleMatrix(*n_fixed, *n_fixed);
  double ***Zgrp = doubleMatrix3D(*n_grp, *max_samp_grp, *n_random);

  /* get random seed */
  GetRNGstate();

  /* packing the data */
  itemp = 0; 
  for (j = 0; j < *n_fixed; j++)
    for (i = 0; i < *n_samp; i++) 
      X[i][j] = dX[itemp++];

  itemp = 0;
  for (j = 0; j < *n_grp; j++)
    vitemp[j] = 0;
  for (i = 0; i < *n_samp; i++) {
    for (j = 0; j < *n_random; j++) {
      Zgrp[grp[i]][vitemp[grp[i]]][j] = dZ[itemp++];
      /* Rprintf("%5g ", Zgrp[grp[i]][vitemp[grp[i]]][j]); */
    }
    /* Rprintf("\n"); */
    vitemp[grp[i]]++;
  }
  
  /* packing the prior */
  itemp = 0;
  for (k = 0; k < *n_random; k++)
    for (j = 0; j < *n_random; j++) {
      Psi[j][k] = dPsi[itemp];
      itemp++;
    }

  dinv(Psi, *n_random, PsiInv);

  itemp = 0;
  for (j = 0; j < *n_random; j++)
    gamma0[j] = 0;
  for (j = 0; j < *n_grp; j++)
    rMVN(gamma[j], gamma0, PsiInv, *n_random);

  itemp = 0; 
  for (k = 0; k < n_fixed[0]; k++)
    for (j = 0; j < n_fixed[0]; j++)
      A0[j][k] = dA0[itemp++];

  itemp = 0; 
  for (k = 0; k < n_fixed[0]; k++)
    for (j = 0; j < n_fixed[0]; j++)
      Tune[j][k] = tune_fixed[itemp++];

  itemp = 0; 
  for (k = 0; k < *n_random; k++)
    for (j = 0; j < *n_random; j++)
      T0[j][k] = dT0[itemp++];

  acc_fixed[0] = 0;
  for (j = 0; j < n_grp[0]; j++) 
    acc_random[j] = 0;

  /* Gibbs Sampler! */
  for(main_loop = 1; main_loop <= *n_gen; main_loop++) {

    rrLogitMixed(Y, X, Zgrp, *p, *p1, grp, beta, gamma, Psi, 
		 *n_samp, *n_fixed, *n_random, *n_grp,
		 beta0, A0, *tau0, T0, Tune, tune_random,
		 1, acc_fixed, acc_random);
    
    if (*verbose) {
      Rprintf("acceptance ratio for fixed effects:%5g\n", 
	      (double) acc_fixed[0] / (double) main_loop); 
      Rprintf("acceptance ratio for random effects:\n");
      for (j = 0; j < *n_grp; j++)
        Rprintf("%3g ", (double) acc_random[j] / (double) main_loop); 
      Rprintf("\n"); 
    }

    R_FlushConsole(); 
    /* Storing the output */
    for (j = 0; j < *n_fixed; j++)
      betaStore[ibeta++] = beta[j];
    for (j = 0; j < *n_grp; j++) 
      for (k = 0; k < *n_random; k++)
	gammaStore[igamma++] = gamma[j][k];
    for (j = 0; j < *n_random; j++) 
      for (k = j; k < *n_random; k++)
	PsiStore[iPsi++] = Psi[j][k]; 

    R_FlushConsole(); 
    R_CheckUserInterrupt();
  } /* end of Gibbs sampler */

  PutRNGstate();

  /* freeing memory */
  free(vitemp);
  free(gamma0);
  FreeMatrix(X, *n_samp);
  FreeMatrix(gamma, *n_grp);
  FreeMatrix(Psi, *n_random);
  FreeMatrix(PsiInv, *n_random);
  FreeMatrix(A0, *n_fixed);
  FreeMatrix(Tune, *n_fixed);
  FreeMatrix(T0, *n_random);
  Free3DMatrix(Zgrp, *n_grp, *max_samp_grp);
} 

