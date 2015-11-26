data {
	int<lower=0> M;		// Number of samples
	int<lower=0> n[M];	// Success count in each sample
	int<lower=0> m[M];	// Failure count in each sample

	real<lower=0> PSS;	// Prior sample size
	real<lower=0> alpha;// Parameters of the gamma tau prior
	real<lower=0> beta;	// Parameters of the gamma tau prior
}

parameters {
	real x[M];			// logit-transformed success rate in each sample
	real mu;			// Underlying mean success rate (logit-transformed)
	real<lower=0> tau;	// Precision of success rate distribution
}

transformed parameters {
	real<lower=0,upper=1> r[M];	// Success rate in each sample
	for (i in 1:M)
		r[i] <- inv_logit(x[i]);
}

model {
	for (i in 1:M) {
		n[i] ~ binomial(n[i] + m[i], r[i]);
		x[i] ~ normal(mu, 1/tau);
	}

	tau ~ gamma(alpha, beta);
	mu ~ normal(0, 1/(tau*PSS));
}
