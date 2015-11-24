library(ggplot2)
library(scales)
library(plyr)


texquote = function(str) gsub("_", "\\\\_", sub("\\s+$", "", str))


formatC2 = function(x, ...)
{
    if (is.na(x))
        return("")
    else
        return(formatC(x, ...))
}


plotGenomeBreakdown = function(f_targ_of_wg, f_gold_of_targ,
    B = 0.3, D = 0.06, E = 0.03, 
    mar.left = 0.0, mar.right = 0.0, mar.top = 0.0, mar.bottom = 0.2,
    lwd.hatch = 1, lwd.box = 2, lwd.divide = 2, lwd.link = 1, 
    density.targ = NULL, angle.targ = 45, fill.targ = grey(0.8),
    density.giab = NULL, angle.giab = 135, fill.giab = grey(0.5))
{
    A = f_targ_of_wg
    C = f_gold_of_targ

    height = 1 - mar.top - mar.bottom
    width = 1 - mar.left - mar.right

    pars = par(no.readonly = TRUE)
    par(mar = c(0, 0, 0, 0))
    plot.new()

    # Top box hatched region (target)
    rect(mar.left, 1-(mar.top+B*height), mar.left+A*width, 1-mar.top, lwd = lwd.hatch, density = density.targ, angle = angle.targ, border = NA, col = fill.targ)
    # Top box divider
    segments(mar.left+A*width, 1-mar.top, mar.left+A*width, 1-(mar.top+B*height), lwd = lwd.divide)
    # Top box outer border
    rect(mar.left, 1-(mar.top+B*height), 1-mar.right, 1-mar.top, lwd = lwd.box)

    # Bottom box hatched region (target)
    rect(mar.left, mar.bottom+B*height, 1-mar.right, mar.bottom, lwd = lwd.hatch, density = density.targ, angle = angle.targ, border = NA, col = fill.targ)
    # Bottom box overshaded region (GiaB callable within target)
    rect(mar.left, mar.bottom+B*height, mar.left+C*width, mar.bottom, lwd = lwd.hatch, density = density.giab, angle = angle.giab, border = NA, col = fill.giab)
    # Bottom box divider
    segments(mar.left+C*width, mar.bottom+B*height, mar.left+C*width, mar.bottom, lwd = lwd.divide)
    # Bottom box outer border
    rect(mar.left, mar.bottom+B*height, 1-mar.right, mar.bottom, lwd = lwd.box)

    # Lines joining subsets
    segments(mar.left+c(0, A, A, 1)*width, mar.bottom+c(1-B-D, 1-B-D, 1-B-D-E, B+D+E)*height, mar.left+c(0, A, 1, 1)*width, mar.bottom+c(B+D, 1-B-D-E, B+D+E, B+D)*height, lwd = lwd.link)

    legend("bottom", legend = c("Genome", "Supplied region", "In gold standard"), fill = c(grey(1), fill.targ, fill.giab), horiz = TRUE, bty = "n")

    par(pars)
}


marginalizePerformance = function(perf_data, subset, vars, ...)
{
    ddply(perf_data[eval(subset, perf_data, parent.frame()),], vars, function(x) { 
        ntp = sum(x$ntp)
        nfn = sum(x$nfn)
        nfp = sum(x$nfp)
        ntn = sum(x$ntn)

        if (ntp + nfn == 0)
        {
            sens = NA
            sens.lci = NA
            sens.uci = NA
        }
        else
        {
            ci_test = binom.test(ntp, ntp + nfn, ...)
            sens = as.vector(ci_test$estimate)
            sens.lci = ci_test$conf.int[1]
            sens.uci = ci_test$conf.int[2]
        }

        c("sens" = sens, "sens.lci" = sens.lci, "sens.uci" = sens.uci, "n" = ntp + nfn + nfp + ntn, ntp = ntp, nfn = nfn, nfp = nfp, ntn = ntn)
    }, .drop = TRUE)
}


calcPerformanceStats = function(data, subset, vars, ..., model = "bin", conf_level = 0.95)
{
    marginalized_values = ldply(data, function(d) marginalizePerformance(d$class_subsets.performance_thresholded, subset, vars, ...))

    combined_stats = ddply(marginalized_values, vars, function(m) {
        sens_stats = replicatedBinomialCI(m$ntp, m$nfn, conf_level, model = model)
        fpr_stats = replicatedBinomialCI(m$nfp, m$n - m$nfp, conf_level, model = model)
        names(sens_stats) = paste("sens.", names(sens_stats), sep = "")
        names(fpr_stats) = paste("fpr.", names(fpr_stats), sep = "")
        unlist(c(sens_stats, fpr_stats))
    })

    combined_stats
}


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


replicatedBinomialCI = function(successes, failures, conf_level, model = c("bin", "betabin", "bayes"))
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
        for (i in 1:50000)
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

            if (i > 1000 && i %% 50 == 0)
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


replicatedBinomialCI(successes = c(0, 0, 4, 151), failures = c(0, 0, 0, 0), conf_level = 0.95, model = "bin")
replicatedBinomialCI(successes = c(0, 0, 4, 151), failures = c(0, 0, 0, 0), conf_level = 0.95, model = "betabin")
replicatedBinomialCI(successes = c(0, 0, 4, 151), failures = c(0, 0, 0, 0), conf_level = 0.95, model = "bayes")

# betaBinomML(c(0, 0, 4, 151), c(0, 0, 4, 151))

replicatedBinomialCI(successes = c(50, 70), failures = c(50, 30), conf_level = 0.95, model = "bin")
replicatedBinomialCI(successes = c(50, 7), failures = c(50, 3), conf_level = 0.95, model = "bin")
replicatedBinomialCI(successes = c(50, 70), failures = c(50, 30), conf_level = 0.95, model = "betabin")
replicatedBinomialCI(successes = c(50, 7), failures = c(50, 3), conf_level = 0.95, model = "betabin")
replicatedBinomialCI(successes = c(50, 70), failures = c(50, 30), conf_level = 0.95, model = "bayes")
replicatedBinomialCI(successes = c(50, 7), failures = c(50, 3), conf_level = 0.95, model = "bayes")


replicatedBinomialCI(successes = c(50, 50), failures = c(50, 50), conf_level = 0.95, model = "bin")
# 0.5 (0.431-0.569)
replicatedBinomialCI(successes = c(50, 50), failures = c(50, 50), conf_level = 0.95, model = "betabin")
# 0.4999 (0.4993-0.5005)
replicatedBinomialCI(successes = c(50, 50), failures = c(50, 50), conf_level = 0.95, model = "bayes")
# 0.5008 (0.239-0.7911)

replicatedBinomialCI(successes = c(50, 50, 50, 50, 50), failures = c(50, 50, 50, 50, 50), conf_level = 0.95, model = "bin")
# 0.5 (0.4563-0.5437)
replicatedBinomialCI(successes = c(50, 50, 50, 50, 50), failures = c(50, 50, 50, 50, 50), conf_level = 0.95, model = "betabin")
# 0.4999 (0.4993-0.5005)
replicatedBinomialCI(successes = c(50, 50, 50, 50, 50), failures = c(50, 50, 50, 50, 50), conf_level = 0.95, model = "bayes")
# 0.5012 (0.3732-0.636)



replicatedBinomialCI(successes = rep(5, 10), failures = rep(5, 10), conf_level = 0.95, model = "bayes")
# 0.5 (0.4111-0.6124)
replicatedBinomialCI(successes = rep(50, 10), failures = rep(50, 10), conf_level = 0.95, model = "bayes")
# 0.5 (0.4261-0.5708)



replicatedBinomialCI(successes = c(100, 100), failures = c(0, 0), conf_level = 0.95, model = "bin")
# 1 (NaN-NaN)
replicatedBinomialCI(successes = c(100, 100), failures = c(0, 0), conf_level = 0.95, model = "betabin")
# 1 (0.9969-1)
replicatedBinomialCI(successes = c(100, 100), failures = c(0, 0), conf_level = 0.95, model = "bayes")
# Fails

replicatedBinomialCI(successes = c(100, 100, 100, 100, 100), failures = c(0, 0, 0, 0, 0), conf_level = 0.95, model = "bin")
# 1 (NaN-NaN)
replicatedBinomialCI(successes = c(100, 100, 100, 100, 100), failures = c(0, 0, 0, 0, 0), conf_level = 0.95, model = "betabin")
# 1 (0.9969-1)
replicatedBinomialCI(successes = c(100, 100, 100, 100, 100), failures = c(0, 0, 0, 0, 0), conf_level = 0.95, model = "bayes")
# 0.8113 (0.5-0.9464)



# Maxima code to derive adaptive importance sampling steps 
# for Gibbs sampler:
#
# dbinom(k, n, p) := gamma(n+1)/(gamma(k+1)*gamma(n-k+1))*p^k*(1-p)^(n-k);
# dgamma(x, alpha, beta) := beta^alpha/gamma(alpha)*x^(alpha-1)*exp(-beta*x);
# dnorm(x, mu, tau) := sqrt(tau/2*%pi)*exp(-tau*(x-mu)^2/2);
# ilogit(x) := exp(x)/(exp(x)+1);
# countPr(i) := dbinom(n[i], T[i], ilogit(x[i]))*dnorm(x[i], mu, tau);
# jointPr : product(countPr(i), i, 1, M)*dnorm(mu, 0, S*tau)*dgamma(tau, alpha, beta);
# lnjointPr : radcan(log(jointPr));
# lncondPr_xj : facsum(ev(K1 + radcan(log(countPr(j))), logexpand = all), x[j]);
# lncondPr_mu : facsum(ev(K2 + lnjointPr, logexpand = all), mu);
# lncondPr_tau : facsum(ev(K3 + lnjointPr, logexpand = all), tau);
# declare (nounify(sum), linear);
# dlncondPr_xj : diff(lncondPr_xj, x[j], 1);
# # (2*(mu*tau+n[j])-2*x[j]*tau-(2*%e^x[j]*T[j])/(%e^x[j]+1))/2
# dlncondPr_mu : diff(lncondPr_mu, mu, 1);
# # (tau*(2*(sum(x[i],i=1,1,M))-2*mu*M)-2*mu*tau*S)/2
# dlncondPr_tau : diff(lncondPr_tau, tau, 1);
# # (-mu^2*S+M/tau-mu^2*M+(2*alpha)/tau-1/tau+2*(sum(x[i],i=1,1,M))*mu-(sum(x[i]^2,i=1,1,M))-2*beta)/2
# klncondPr_xj : facsum(integrate(dlncondPr_xj, x[j]), x[j]);
# # (2*x[j]*(mu*tau+n[j])-x[j]^2*tau-2*T[j]*log(%e^x[j]+1))/2
# klncondPr_mu : facsum(integrate(dlncondPr_mu, mu), mu);
# # (2*(sum(x[i],i=1,1,M))*mu*tau-mu^2*tau*(S+M))/2
# klncondPr_tau : facsum(integrate(dlncondPr_tau, tau), tau);
# # (log(tau)*(M+2*alpha-1)-tau*(mu^2*S+mu^2*M-2*(sum(x[i],i=1,1,M))*mu+(sum(x[i]^2,i=1,1,M))+2*beta))/2
#
# The klncondPr_var variables are K + log(Pr(var|others))
# The dlncondPr_var variables are diff(klncondPr_var, var)
