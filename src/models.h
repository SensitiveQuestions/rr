void rrLogit(int *Y, double **X, double p, double p1, double *beta, int n_samp, 
	     int n_cov, double *beta0, double **A0, double **InvVar, int n_gen, 
	     int *counter); 

void rrLogitMixed(int *Y, double **X, double ***Z, double p, double p1, int *grp, 
		  double *beta, double **gamma, double **Psi, int n_samp, int n_fixed,
		  int n_random, int n_grp, double *beta0, double **A0, int tau0,
		  double **T0, double **tune_fixed, double *tune_random, int n_gen,
		  int *acc_fixed, int *acc_random);

