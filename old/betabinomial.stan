//
// This Stan program defines a simple model for computing the posterior distribution
// for the Beta Binomial fully hierarchical model.
//
data {
	int<lower=0> N;         // Number of data (all features and cells)
	int<lower=0> J;         // # of features
	int<lower=0> n[N];      // total number of trials (long form)
	int<lower=0> y[N];      // number of successes (long form)
	int<lower=0> s[J];      // # of cells with obserations in each feature
	real a0_mu;             // Beta shape1 prior on proportions \mu
	real<lower=0> b0_mu;    // Beta shape2 prior on proportions \mu
	real a0_delta;          // Beta shape1 prior on population proportions \mu
	real<lower=0> b0_delta; // Beta shape2 prior on population proportions \mu
	real a0_rho;            // Beta shape1 prior on population dispersion \gamma
	real<lower=0> b0_rho;   // Beta shape2 prior on population dispersion \gamma
}

// Parameters to compute posterior, reparametrization of BB in terms of mean
// and overdispersion.
parameters {
	real<lower=0,upper=1> mu[J];       // Mean of BB
	real<lower=0,upper=1> gamma[J];    // Dispersion of BB
	real<lower=0,upper=1> delta;    // Population mean of BB
	real<lower=0,upper=1> rho; // Population dispersion of BB
	// real<lower=0,upper=1> theta[N];    // Proportions for conditioned Binomials
}

// Original parameters of the Beta distribution
transformed parameters {
	real<lower=0> alpha[J];
	real<lower=0> beta[J];
	real<lower=0> alpha_joint;
	real<lower=0> beta_joint;
	for (j in 1:J) {
	  alpha[j] = mu[j]/gamma[j] - mu[j];
	  beta[j] = (1 - mu[j])/gamma[j] + mu[j] - 1;
	}
	alpha_joint = delta/rho - delta;
	beta_joint = (1 - delta)/rho + delta - 1;
}

// Joint distribution for Beta Binomial model
model {
  int pos; // Counter to perform iteration only for cells with coverage
  pos = 1;

  delta ~ beta(a0_delta, b0_delta);         // Hyperprior on the population mean of BB
  rho ~ beta(a0_rho, b0_rho);               // Hyperprior on the population dispersion of BB
  mu ~ beta(a0_mu, b0_mu);                  // Hyperprior on proportions parameter
  gamma ~ beta(alpha_joint, beta_joint);    // Hyperprior on dispersion parameter
  // for (j in 1:J) {
  //   segment(theta, pos, s[j]) ~ beta(alpha[j], beta[j]); // prior for nuisance parameter \theta
  //   segment(y, pos, s[j]) ~ binomial(segment(n, pos, s[j]), segment(theta, pos, s[j])); // likelihood
  //   pos = pos + s[j];       // Increment counter
  // }
  // Beta-Binomial formulation, i.e. integrating out theta
  for (j in 1:J) {
    segment(y, pos, s[j]) ~ beta_binomial(segment(n, pos, s[j]), alpha[j], beta[j]); // likelihood
    pos = pos + s[j];     // Increment counter
  }
}
