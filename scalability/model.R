beta_binomial <- "
  data {
    int<lower=1> N;
    int<lower=1> s;
    int<lower=0> y[N];
  }
  parameters {
    real<lower=0, upper=1> theta;
  }
  model { 
    theta ~ beta(1,1);
    for (i in 1:N)
      y[i] ~ binomial(s,theta);
  }
"