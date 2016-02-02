library(ggplot2)
library(scales)
library(plyr)
suppressPackageStartupMessages(library(lme4))

suppressPackageStartupMessages(library(robustbase)) # For Qn

suppressPackageStartupMessages(library(rstan))
suppressPackageStartupMessages(library(coda))       # For HPDinterval
STAN_LOGIT = stan_model("logit.stan")

set.seed(as.double(Sys.time()))


texquote = function(str) gsub("_", "\\\\_", sub("\\s+$", "", str))


formatC2 = function(x, ...)
{
    if (is.na(x))
        return("")
    else
        return(formatC(x, ...))
}


categoriseCallConfidence = function(data)
{
    # From https://kccg.garvan.org.au/confluence/display/DQMS/NATA-346+-+Editable+Copy+-+Appendix+1%3A+Analytical+Performance
    # Variants that pass VQSR are to be broken down into four sensitivity classes, based on variant type, size, and sequence depth around the variant:
    #   * Class H (High confidence):    99% <= s
    #   * Class C (Confident):          95% <= s < 99%
    #   * Class L (Limited confidence): 80% <= s < 95%
    #   * Class I (Indicative):         s < 80%
    # The criteria for each VQSR-passing variant were defined on the basis of the latest WGS performance report:
    #   * SNV:
    #      * Depth >= 20X:                  Class H
    #      * 20X > Depth >= 15X:            Class C
    #      * 15X > Depth >= 10X:            Class L
    #      * 10X > Depth:                   Class I
    #   * Ins, size [1,19] (short ins):
    #      * Depth >= 20X:                  Class C
    #      * 20X > Depth >= 15X:            Class L
    #      * 15X > Depth:                   Class I
    #   * Ins, size [20,inf) (long ins):    Class I
    #   * Del, size [1,19] (short del):
    #      * Depth >= 10X:                  Class C
    #      * 10X > Depth:                   Class I
    #   * Del, size [20,inf) (long del):    Class I
    #   * MNV:                              Class I

    llply(data, function (d) {
        perf = d$class_subsets.performance_thresholded
        perf$confidence_class = "Unclassified"
        perf$confidence_class[perf$muttype == "Subst" & perf$mutsize == "1" & perf$depth >= "20-24"] = "HighConfidence"
        perf$confidence_class[perf$muttype == "Subst" & perf$mutsize == "1" & perf$depth >= "15-19" & perf$depth < "20-24"] = "Confident"
        perf$confidence_class[perf$muttype == "Subst" & perf$mutsize == "1" & perf$depth >= "10-14" & perf$depth < "15-19"] = "LimitedConfidence"
        perf$confidence_class[perf$muttype == "Subst" & perf$mutsize == "1" & perf$depth < "10-14"] = "Indicative"
        perf$confidence_class[perf$muttype == "Subst" & perf$mutsize > "1"] = "Indicative"
        perf$confidence_class[perf$muttype == "Ins" & perf$mutsize < "20-29" & perf$depth >= "20-24"] = "Confident"
        perf$confidence_class[perf$muttype == "Ins" & perf$mutsize < "20-29" & perf$depth >= "15-19" & perf$depth < "20-24"] = "LimitedConfidence"
        perf$confidence_class[perf$muttype == "Ins" & perf$mutsize < "20-29" & perf$depth < "15-19"] = "Indicative"
        perf$confidence_class[perf$muttype == "Ins" & perf$mutsize >= "20-29"] = "Indicative"
        perf$confidence_class[perf$muttype == "Del" & perf$mutsize < "20-29" & perf$depth >= "10-14"] = "Confident"
        perf$confidence_class[perf$muttype == "Del" & perf$mutsize < "20-29" & perf$depth < "10-14"] = "Indicative"
        perf$confidence_class[perf$muttype == "Del" & perf$mutsize >= "20-29"] = "Indicative"
        perf$confidence_class[perf$muttype == "Other"] = "Indicative"
        perf$confidence_class = ordered(perf$confidence_class, levels = c("Unclassified", "Indicative", "LimitedConfidence", "Confident", "HighConfidence"))
        d$class_subsets.performance_thresholded = perf
        d
    })
}


relabelZyg = function(perfstats)
{
    if ("zyg" %in% colnames(perfstats))
        perfstats$zyg = ordered(c("RR" = "HomRef", "RA" = "Het", "AA" = "HomAlt", "AB" = "HetAlt")[as.character(perfstats$zyg)], levels = c("Hom", "Het", "HomAlt", "HetAlt"))
    perfstats
}


relabelMuttype = function(perfstats)
{
    if ("muttype" %in% colnames(perfstats))
        perfstats$muttype = ordered(c("Subst" = "SNVs", "Ins" = "Insertions", "Del" = "Deletions", "Other" = "Other", "None" = "None")[as.character(perfstats$muttype)], levels = c("SNVs", "Insertions", "Deletions", "Other", "None"))
    perfstats
}


dropSizeZero = function(perfstats)
{
    if ("mutsize" %in% colnames(perfstats))
        perfstats = perfstats[perfstats$mutsize != "0",]
    if ("mutsize_binned" %in% colnames(perfstats))
        perfstats = perfstats[perfstats$mutsize_binned != "0",]
    perfstats
}


dropDepthUnknown = function(perfstats)
{
    if ("depth" %in% colnames(perfstats))
        perfstats = perfstats[perfstats$depth != "Unknown",]
    perfstats
}


dropMuttypesOtherNone = function(perfstats)
{
    if ("muttype" %in% colnames(perfstats))
        perfstats = perfstats[!(perfstats$muttype %in% c("Other", "None")),]
    perfstats
}


generateSensitivityCountTable = function(data, subset, vars)
{
    tbl = ldply(data, function(x) marginalizePerformance(x$class_subsets.performance_thresholded, subset, vars))
    id_translator = unlist(llply(data, function(x) x$label))
    tbl$.id = id_translator[tbl$.id]

    giveTablePrettyHeaders(dropSizeZero(dropDepthUnknown(relabelZyg(tbl))))
}


generateSensitivitySummaryTable = function(data, subset, vars, ci_method, ci_level)
{
    stats = calcPerformanceStats(data, subset, vars, model = ci_method, conf_level = ci_level)
    stats$sens.ci = sprintf("%.4f-%.4f", stats$sens.lcl, stats$sens.ucl)

    giveTablePrettyHeaders(dropSizeZero(dropDepthUnknown(relabelZyg(stats))))
}


giveTablePrettyHeaders = function(tbl)
{
    # The table columns will be in this order.  Commented-out or missing
    # columns will not be kept in the output.
    colname_translator = c(
        ".id" = "Sample Label",
        "zyg" = "Zygosity",
        "depth" = "Depth",
        "mutsize" = "Size",
        "mutsize_binned" = "Size",
        "muttype" = "Variant class",
        "confidence_class" = "Confidence class",
        "ntp" = "True positives",
        "nfn" = "False negatives",
        # "n" = "Total variants",
        "sens" = "Sensitivity",
        "sens.est" = "Sensitivity",
        "sens.ci" = "CI"
    )

    tbl = tbl[,colnames(tbl) %in% names(colname_translator)]
    colname_translator = colname_translator[names(colname_translator) %in% colnames(tbl)]
    tbl = tbl[,names(colname_translator)]
    colnames(tbl) = colname_translator[colnames(tbl)]

    tbl
}


convertSensitivityCountTableToXTable = function(tbl, ...)
{
    # tbl is from generateSensitivityCountTable
    n = ncol(tbl)
    xtable(tbl, digits = c(rep(0, n), 4), display = c(rep("s", n - 2), "d", "d", "f"), ...)
}


convertSensitivitySummaryTableToXTable = function(tbl, total_var_counts = NULL, ci_heading = "CI", ...)
{
    # total_var_counts must be generated by a call like: total_var_counts = dropSizeZero(marginalizePerformance(data[[1]]$class_subsets.performance_thresholded, subset, vars))
    tbl[,"CI"] = paste("$", tbl[,"CI"], "$", sep = "")
    colnames(tbl)[colnames(tbl) == "CI"] = ci_heading
    nvars = ncol(tbl) - 2

    if (!is.null(total_var_counts))
    {
        # TODO: Should check that rows of tbl and total_var_counts match.  Not trivial
        # at this time, as columns are renamed, and factor codings are changed.
        old_colnames = colnames(tbl)
        n = ncol(tbl)
        tbl = cbind(tbl[,1:(n-2)], total_var_counts$ntp + total_var_counts$nfn, tbl[,(n-1):n])
        colnames(tbl) = c(old_colnames[1:(n-2)], "True variants", old_colnames[(n-1):n])
        n = ncol(tbl)
        xtbl = xtable(tbl, digits = c(rep(0, n-1), 4, 4), display = c(rep("s", n-2), "d", "f", "f"), ...)
    }
    else
    {
        n = ncol(tbl)
        xtbl = xtable(tbl, digits = c(rep(0, n-1), 4, 4), display = c(rep("s", n-1), "f", "f"), ...)
    }

    xtbl
}


generateSensitivitySummaryXTable = function(data, subset, vars, ci_method, ci_level, add_total_var_count = FALSE, ...)
{
    tbl = generateSensitivitySummaryTable(data, subset, vars, ci_method, ci_level)
    tbl[,"CI"] = paste("$", tbl[,"CI"], "$", sep = "")

    if (add_total_var_count)
    {
        counts = dropSizeZero(marginalizePerformance(data[[1]]$class_subsets.performance_thresholded, subset, vars))
        old_colnames = colnames(tbl)
        n = ncol(tbl)
        tbl = cbind(tbl[,1:(n-2)], counts$ntp + counts$nfn, tbl[,(n-1):n])
        colnames(tbl) = c(old_colnames[1:(n-2)], "True variants", old_colnames[(n-1):n])
        n = ncol(tbl)
        xtbl = xtable(tbl, digits = c(rep(0, n-1), 4, 4), display = c(rep("s", n-2), "d", "f", "f"), ...)
    }
    else
    {
        n = ncol(tbl)
        xtbl = xtable(tbl, digits = c(rep(0, n-1), 4, 4), display = c(rep("s", n-1), "f", "f"), ...)
    }

    xtbl
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


calcPerformanceStats = function(data, subset, vars, ..., model = "betabin", conf_level = 0.95)
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


replicatedBinomialCI = function(successes, failures, conf_level, model = c("betabin.ml", "binomial.ml", "logit.bayes", "none"))
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

    if (model == "betabin.ml")
    {
        # Approach 1: ML on beta-binomial fit:
        fit = betaBinomML(successes, successes + failures)
        ci = qbeta(c(lower = (1-conf_level)/2, median = 0.5, upper = 1-(1-conf_level)/2), shape1 = fit[1], shape2 = fit[2])
        result$est = ci[["median"]]
        result$lcl = ci[["lower"]]
        result$ucl = ci[["upper"]]
    }
    else if (model == "binomial.ml")
    {
        # The fit is flaky for unusual inputs, so try and sanitise 
        # them as much as possible, and handle errors gracefully.
        sel = successes + failures > 0
        successes = successes[sel]
        failures = failures[sel]
        if (length(successes) > 0)
        {
            result$est = sum(successes) / sum(successes + failures)
            data = data.frame(sample = factor(seq_along(successes)), successes = successes, failures = failures)
            if (diff(range(data$successes)) > 0)
            {
                fit = try(glmer(cbind(successes,failures) ~ 1 + (1 | sample), data = data, family = "binomial"))
                if (class(fit) != "try-error")
                {
                    ci = try(confint(fit, parm = "(Intercept)", level = conf_level, method = "profile"))
                    if (class(ci) != "try-error")
                    {
                        result$lcl = plogis(ci[[1]])
                        result$ucl = plogis(ci[[2]])
                    }
                }
            }
        }
    }
    else if (model == "logit.bayes")
    {
        # Prior:
        #   PSS: Prior 'sample size' / strength
        #   alpha, beta: parameters of gamma distribution prior for the precision tau.
        #   
        #   Derive alpha, beta from mean, sd of tau as:
        #     alpha = (mean/sd)^2     beta = mean/sd^2
        #   Reasonable starting values are PSS = 1, alpha = 5, beta = 10 (prior mass
        #   is concentrated at high and low rates), and PSS = 1, alpha = 5, beta = 5 
        #   (prior is approximately uniform, with slight peaks at the extremes).
        #   PSS is less obvious.  This could be tuned.

        stan_init_func = function(chain_id) { 
            old_seed = .Random.seed
            set.seed(old_seed + chain_id)
            init_r = successes / (successes + failures)
            init_r[is.na(init_r)] = 0.5
            init_r = init_r * 0.99 + 0.005
            init_x = qlogis(init_r)
            init_mu = mean(init_x)
            init_x = init_x + rnorm(length(init_x), 0, 0.1)
            init_tau = max(0.001, rnorm(1, 0.7, 0.2))
            set.seed(old_seed)
            list(x = init_x, mu = init_mu, tau = init_tau)
        }

        data = list(M = length(successes), PSS = 1, alpha = 5, beta = 10, n = successes, m = failures)

        # Original chatty code, always reported timing info.  See the silencing hack below.
        # stan.result = sampling(STAN_LOGIT, data = data, pars = c("mu", "tau", "r"), chains = 5, iter = 10000, thin = 10, init = stan_init_func, refresh = 0)
        
        # Hack to get around Stan being unavoidably chatty.  Bit of a hack, but we don't 
        # lose the stan exception warnings at least (are they in C code?)
        sink("/dev/null")
        stan.result = sampling(STAN_LOGIT, data = data, pars = c("mu", "tau", "r"), chains = 5, iter = 10000, thin = 10, init = stan_init_func, refresh = 0)
        sink()

        mu_samples = extract(stan.result, "mu")[[1]]
        mu_ci = HPDinterval(mcmc(as.vector(mu_samples)), prob = conf_level)
        result$est = plogis(median(mu_samples))
        result$lcl = plogis(mu_ci[[1]])
        result$ucl = plogis(mu_ci[[2]])
    }
    else if (model == "none")
        result$est = sum(successes) / (sum(successes) + sum(failures))

    result
}
