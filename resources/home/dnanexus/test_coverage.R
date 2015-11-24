library(plyr)
library(ggplot2)
library(abind)


betaBinomML = function(x, S)
{
    # Fit a Beta-binomial distribution:
    # x[i] ~ Binom(p[i]; S[i])
    # p[i] ~ Beta(a, b)
    # Returns the vector c(a, b) containing optimised values
    # for the beta parameters.
    
    # Sanity checks
    stopifnot(x <= S)
    stopifnot(S >= 0)

    if (all(S == 0))
        return(c(NA, NA))

    # Drop uninformative S=0 cases.  Not strictly necessary (doesn't affect
    # result), but reduces work.
    x = x[S != 0]
    S = S[S != 0]

    objective = function(logpars, .counts, .trials)
    {
        # Use a log transformation on the parameters to enforce positivity.
        a = exp(logpars[1])
        b = exp(logpars[2])
        logliks = lbeta(.counts+a,.trials-.counts+b)-lbeta(a,b)-lgamma(.counts+1)-lgamma(.trials-.counts+1)
        -sum(logliks)
    }

    # Use Nelder-Mead, as the derivatives are numerically tricky for large S.  For
    # reference, here is the derivation:
    # dbeta(x,a,b) := x^(a-1)*(1-x)^(b-1)/beta(a,b);
    # dbinom(x,r,S) := gamma(S+1)/(gamma(x+1)*gamma(S-x+1))*r^x*(1-r)^(S-x);
    # sum(log(integrate(dbinom(x[i],r[i],S)*dbeta(r[i],a,b), r[i], 0, 1)), i, 1, m);
    #   "Is "x[g18936]+a" positive, negative or zero?"positive;
    #   "Is "x[g18936]-S-b" positive, negative or zero?"negative;
    #   sum(log((beta(x[i]+a,S-x[i]+b)*gamma(S+1))/(beta(a,b)*gamma(x[i]+1)*gamma(S-x[i]+1))),i=1,1,m)
    # radcan(diff(sum(log((beta(x[i]+a,S-x[i]+b)*gamma(S+1))/(beta(a,b)*gamma(x[i]+1)*gamma(S-x[i]+1))),i,1,m),a,1));
    #   -sum(gamma(S+1)*psi[0](S+b+a)+(-psi[0](x[i]+a)-psi[0](b+a)+psi[0](a))*gamma(S+1),i=1,1,m)/gamma(S+1)
    # radcan(diff(sum(log((beta(x[i]+a,S-x[i]+b)*gamma(S+1))/(beta(a,b)*gamma(x[i]+1)*gamma(S-x[i]+1))),i,1,m),b,1));
    #   sum(gamma(S+1)*psi[0](S-x[i]+b)-gamma(S+1)*psi[0](S+b+a)+(psi[0](b+a)-psi[0](b))*gamma(S+1),i=1,1,m)/gamma(S+1)

    # test.a = seq(-5, 5, 0.5)
    # test.b = seq(-5, 5, 0.5)
    # test.obj = as.matrix(sapply(test.a, function(a) sapply(test.b, function(b) objective(c(a, b), x, S))))
    # dimnames(test.obj) = list(b = round(test.b, 2), a = round(test.a, 2))
    # print(test.obj-min(test.obj))

    # The likelihood appears well-behaved -- a uniform starting point seems to work well:
    start = c(0, 0)

    opt = optim(start, fn = objective, .counts = x, .trials = S, method = "Nelder-Mead")$par

    # In the single-class case (ie all x = 0, or all x = S), hard-threshold the optimal point 
    # to lie in the box [-4, 4] x [-4, 4] (approx 0.02 to 50 in linear units) for stability.
    if (all(x == 0) || all(x == S))
        opt = pmin(4, pmax(-4, opt))

    exp(opt)
}



replicatedBinomialCI = function(successes, failures, conf_level, model = c("bin", "betabin", "bayes", "combined"))
{
    # Estimate central tendency and confidence limits for the
    # success rate of a binomial process that has been measured
    # in multiple separate experiments.  Two models are 
    # currently implemented: Beta-binomial, and binomial.

    # The Beta-binomial supposes m experiments are performed, 
    # where in each experiment i, the measured variable x[i] is 
    # the number of successes in S[i] binomial trials, and the 
    # success rate in each experiment, r[i], is a random 
    # beta-distributed variable.  Symbolically,
    #  x[i] ~ Binom(r[i], S[i])
    #  r[i] ~ Beta(a, b)
    # The r[i]s allow experiments to have a slightly differing
    # success rate.

    # The binomial model is a straight wrapper around GLM, and
    # uses its profile functions to estimate CIs on the rate.

    # message(sprintf("replicatedBinomialCI(successes = %s, failures = %s, conf_level = %s, model = %s)", deparse(successes), deparse(failures), deparse(conf_level), deparse(model)))

    model = match.arg(model)

    result = list(est = NA, lcl = NA, ucl = NA, conf_level = conf_level)

    if (model == "betabin")
    {
        # Approach 1: ML on beta-binomial fit:
        fit = betaBinomML(successes, successes + failures)
        ci = qbeta(c(lower = (1-conf_level)/2, median = 0.5, upper = 1-(1-conf_level)/2), shape1 = fit[1], shape2 = fit[2])
        result$est = ci[["median"]]
        result$lcl = ci[["lower"]]
        result$ucl = ci[["upper"]]
    }
    else if (model == "bin")
    {
        # Approach 2: GLM
        data = cbind(successes = successes, failures = failures)
        if (!all(data[,"successes"] + data[,"failures"] == 0))
        {
            data = data[data[,"successes"] + data[,"failures"] != 0,,drop=FALSE]
            fit = glm(data ~ 1, family = binomial)
            mean_logit = coef(fit)[["(Intercept)"]]
            ci_logit = suppressMessages(confint(fit, parm = "(Intercept)", level = conf_level))
            result$est = exp(mean_logit) / (1 + exp(mean_logit))
            result$lcl = exp(ci_logit[[1]]) / (1 + exp(ci_logit[[1]]))
            result$ucl = exp(ci_logit[[2]]) / (1 + exp(ci_logit[[2]]))
        }
    }
    else if (model == "combined")
    {
        s = sum(successes)
        n = s + sum(failures)
        test = binom.test(s, n, conf.level = conf_level)
        result$est = test$estimate
        result$lcl = test$conf.int[[1]]
        result$ucl = test$conf.int[[2]]
    }
    else if (model == "bayes")
    {
        # Approach 3: Bayesian
        library(ars)

        # Prepare variables
        T = successes + failures
        n = successes
        M = length(successes)

        # TODO: add prior spec here
        S = 1
        alpha = 1
        beta = 1

        # Sample initial values
        tau = rgamma(1, shape = alpha, rate = beta)
        mu = rnorm(1, mean = 0, sd = 1/(tau*S))
        x = rnorm(M, mean = mu, sd = 1/tau)

        # print(paste(" Iter     Changed        mu       tau      ", paste(sprintf("x[%s]", 1:M), collapse = "      "), sep = ""))
        # print(sprintf("Start          NA    % 6.2f    % 6.2f    %s", mu, tau, paste(sprintf("% 6.2f", x), collapse = "    ")))

        # Conditional distribution functions (from Maxima, see code below)
        klncondPr_xj = function(xj, nj, Tj) (2*xj*(mu*tau+nj)-xj^2*tau-2*Tj*log(exp(xj)+1))/2
        dlncondPr_xj = function(xj, nj, Tj) (2*(mu*tau+nj)-2*xj*tau-(2*exp(xj)*Tj)/(exp(xj)+1))/2
        klncondPr_mu = function(mu)         (2*sum(x)*mu*tau-mu^2*tau*(S+M))/2
        dlncondPr_mu = function(mu)         (tau*(2*sum(x)-2*mu*M)-2*mu*tau*S)/2
        klncondPr_tau = function(tau)       (log(tau)*(M+2*alpha-1)-tau*(mu^2*S+mu^2*M-2*sum(x)*mu+sum(x)+2*beta))/2
        dlncondPr_tau = function(tau)       (-mu^2*S+M/tau-mu^2*M+(2*alpha)/tau-1/tau+2*sum(x)*mu-sum(x)-2*beta)/2

        mustore = c()

        # Run the Gibbs sampler
        for (i in 1:5000)
        {
            mu = ars(n = 1, f = klncondPr_mu, fprima = dlncondPr_mu)
            # print(sprintf("% 5d          mu    % 6.2f    % 6.2f    %s", i, mu, tau, paste(sprintf("% 6.2f", x), collapse = "    ")))

            tau = ars(n = 1, f = klncondPr_tau, fprima = dlncondPr_tau, x = c(1e-1, 1, 1e1), lb = TRUE, xlb = 1e-6)
            # print(sprintf("% 5d         tau    % 6.2f    % 6.2f    %s", i, mu, tau, paste(sprintf("% 6.2f", x), collapse = "    ")))

            for (j in 1:M)
            {
                x[j] = ars(n = 1, f = klncondPr_xj, fprima = dlncondPr_xj, nj = n[j], Tj = T[j])
                # print(sprintf("% 5d   %sx[%d]    % 6.2f    % 6.2f    %s", i, paste(rep(" ", 5-floor(log10(j))), collapse = ""), j, mu, tau, paste(sprintf("% 6.2f", x), collapse = "    ")))
            }

            if (i > 1000 && i %% 10 == 0)
            {
                # print(sprintf("% 5d   %sx[%d]    % 6.2f    % 6.2f    %s", i, paste(rep(" ", 5-floor(log10(j))), collapse = ""), j, mu, tau, paste(sprintf("% 6.2f", x), collapse = "    ")))
                mustore = c(mustore, mu)
            }
        }

        ilogit = function(x) exp(x)/(exp(x)+1)

        result$est = ilogit(median(mustore))
        result$lcl = ilogit(quantile(mustore, (1-conf_level)/2))
        result$ucl = ilogit(quantile(mustore, 1-(1-conf_level)/2))
    }

    result
}



# Test coverage probabilities for the CIs above.
# Model:
#   n[i] ~ binom(r[i]; T)
#   r[i] ~ beta(alpha, beta)
# The goal is to find a CI bracketing some MCT of the beta.

doCoverageTestIteration = function(params)
{
    T = as.integer(round(runif(1, params$Tmin, params$Tmax)))
    alpha = runif(1, params$alpha.min, params$alpha.max)
    beta = runif(1, 0, alpha)
    r = rbeta(params$M, alpha, beta)
    n = rbinom(params$M, T, r)

    cis = lapply(params$methods, function(method) replicatedBinomialCI(n, T - n, params$conf_level, model = method))

    target_statistics = list(
        mean = alpha / (alpha + beta),
        mode = (alpha - 1) / (alpha + beta - 2),
        median = qbeta(0.5, alpha, beta))

    coverage = laply(target_statistics, function(target_statistic) laply(cis, function(ci) c(
        under = ci[["lcl"]] > target_statistic, 
        over = ci[["ucl"]] < target_statistic)))

    dimnames(coverage) = list(statistic = names(target_statistics), method = params$method, call = c("under", "over"))

    coverage = abind(coverage, "good" = !(coverage[,,"under"] | coverage[,,"over"]))
    coverage
}


doCoverageTest = function(params, N = 1000, ...)
{
    result = laply(1:N, function(i) doCoverageTestIteration(params), ...)
    temp = dimnames(result)
    names(temp) = c("run", "statistic", "method", "call")
    dimnames(result) = temp
    result
}


summariseCoverageTest = function(test_results)
{
    table = adply(test_results, c(2,3), function(v) {
        failrate = mean(is.na(v[,"good"]))
        coverrate = mean(v[,"good"], na.rm = TRUE)
        userate = mean(v[,"good"] & !is.na(v[,"good"]))
        underrate = mean(v[,"under"], na.rm = TRUE)
        overrate = mean(v[,"over"], na.rm = TRUE)
        bias = underrate - overrate
        c(fail = failrate, cover = coverrate, use = userate, under = underrate, over = overrate, bias = bias)
    })
    table
}







set.seed(1234)
params.1 = list(Tmin = 5e4, Tmax = 5e4, M = 6, alpha.min = 1, alpha.max = 20, beta.min = 0, conf_level = 0.95, methods = c("bin", "betabin", "combined"))
test.1 = doCoverageTest(params.1, 5000, .progress = "text")
summariseCoverageTest(test.1)
#   statistic   method   fail      cover    use     under      over        bias
# 1      mean      bin 0.0032 0.02427769 0.0242 0.5048154 0.4709069  0.03390851
# 2      mode      bin 0.0032 0.01926164 0.0192 0.2548154 0.7259230 -0.47110754
# 3    median      bin 0.0032 0.02708668 0.0270 0.3792135 0.5936998 -0.21448636
# 4      mean  betabin 0.0000 0.98320000 0.9832 0.0156000 0.0012000  0.01440000
# 5      mode  betabin 0.0000 0.81480000 0.8148 0.0234000 0.1618000 -0.13840000
# 6    median  betabin 0.0000 0.99260000 0.9926 0.0044000 0.0030000  0.00140000
# 7      mean combined 0.0000 0.02420000 0.0242 0.5064000 0.4694000  0.03700000
# 8      mode combined 0.0000 0.01920000 0.0192 0.2550000 0.7258000 -0.47080000
# 9    median combined 0.0000 0.02980000 0.0298 0.3784000 0.5918000 -0.21340000

ggplot(summariseCoverageTest(test.1), aes(x = cover, y = bias, pch = method, colour = statistic)) + geom_point(size = 5)

# Summary:
#   * bin and combined have very poor coverage.
#   * betabin has acceptable coverage, though it tends to conservative


set.seed(1234)
params.2 = list(Tmin = 50, Tmax = 1000, M = 6, alpha.min = 1, alpha.max = 20, beta.min = 0, conf_level = 0.95, methods = c("bin", "betabin", "combined"))
test.2 = doCoverageTest(params.2, 5000, .progress = "text")
summariseCoverageTest(test.2)
#   statistic   method  fail     cover    use     under      over        bias
# 1      mean      bin 0.006 0.2879276 0.2862 0.3826962 0.3293763  0.05331992
# 2      mode      bin 0.006 0.2076459 0.2064 0.1782696 0.6140845 -0.43581489
# 3    median      bin 0.006 0.2796781 0.2780 0.2708249 0.4494970 -0.17867203
# 4      mean  betabin 0.000 0.9698000 0.9698 0.0244000 0.0058000  0.01860000
# 5      mode  betabin 0.000 0.7990000 0.7990 0.0262000 0.1748000 -0.14860000
# 6    median  betabin 0.000 0.9808000 0.9808 0.0110000 0.0082000  0.00280000
# 7      mean combined 0.000 0.2924000 0.2924 0.3822000 0.3254000  0.05680000
# 8      mode combined 0.000 0.2088000 0.2088 0.1766000 0.6146000 -0.43800000
# 9    median combined 0.000 0.2876000 0.2876 0.2672000 0.4452000 -0.17800000
ggplot(summariseCoverageTest(test.2), aes(x = cover, y = bias, pch = method, colour = statistic)) + geom_point(size = 5)

# Summary: broadly as above



