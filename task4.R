## =============================================================================
## ME317 Project 2026 -- Task 4: Dependence and copulas
## VineCopula implementation (no gsl / copula dependency)
##
## Assumes in the environment (from Tasks 1-3):
##   prices_clean : data.frame, date + 10 adjusted price columns
##   log_returns  : data.frame, date + 10 log-return columns
##
## Produces:
##   U                : pseudo-observations (uniform margins)
##   fit_table_pairs  : all 45 pairs x 5 families, log-lik and AIC
##   best_per_pair    : winning family per pair
##   cmp              : Gaussian-vs-t comparison per pair
##   tail_table       : model-implied vs empirical tail dependence
##   fit_table_vine   : 10-dimensional vine fits under family restrictions
##   var_compare      : portfolio VaR under each dependence model
##
## NOTE: not executed in this environment. Run section by section.
## =============================================================================

## ---- 0. Setup ----------------------------------------------------------------

set.seed(12345)

for (p in c("VineCopula", "MASS")) {
  if (!requireNamespace(p, quietly = TRUE))
    stop("Install ", p, ": install.packages('", p, "')")
}
library(VineCopula)
library(MASS)

stopifnot(exists("log_returns"), "date" %in% names(log_returns))

dir.create("results", showWarnings = FALSE)
dir.create("data",    showWarnings = FALSE)

IN_END <- as.Date("2020-12-31")
ALPHAS <- c(0.95, 0.99)

## VineCopula family codes. Restricting the headline analysis to the five
## families taught on the course keeps every number defensible.
##   1 = Gaussian, 2 = Student-t, 3 = Clayton, 4 = Gumbel, 5 = Frank
FAM_TAUGHT <- c(1, 2, 3, 4, 5)

## Rotations, used only as a robustness check in section 2c:
##   13 = survival Clayton (upper tail dependence, mirror of 3)
##   14 = survival Gumbel  (lower tail dependence, mirror of 4)
##   23/24/33/34 = 90 and 270 degree rotations, for negative dependence
FAM_EXTENDED <- c(1, 2, 3, 4, 5, 13, 14, 23, 24, 33, 34)

## Parameter counts by family code. Only family 2 has two parameters here.
TWO_PAR <- c(2)
npar_of <- function(f)
  ifelse(f == 0, 0L, ifelse(f %in% TWO_PAR, 2L, 1L))

## ---- 1. Data preparation -----------------------------------------------------

ret_dates <- as.Date(log_returns$date)
R_all     <- as.matrix(log_returns[, -1, drop = FALSE])
TICKERS   <- colnames(R_all)

in_idx <- ret_dates <= IN_END
R_in   <- R_all[in_idx, , drop = FALSE]
n      <- nrow(R_in)
d      <- ncol(R_in)

cat("In-sample returns for copula fitting: n = ", n, ", d = ", d, "\n", sep = "")
cat("Pairs: ", d * (d - 1) / 2, "\n\n", sep = "")

## Probability integral transform via ranks. Dividing by (n+1) rather than n
## keeps every value strictly inside (0,1), which is necessary because copula
## densities are undefined at the boundary. This is what copula::pobs() does;
## written out here so the script has no dependency on that package.
U <- apply(R_in, 2, function(x) rank(x, ties.method = "average") / (n + 1))
colnames(U) <- TICKERS

stopifnot(all(U > 0), all(U < 1))
cat("Pseudo-observations: range [", round(min(U), 4), ", ", round(max(U), 4),
    "]\n\n", sep = "")

## Note on estimation error: because the margins are estimated (ranks rather
## than known CDFs), this is the pseudo-maximum-likelihood estimator. Standard
## errors are therefore only asymptotically valid, and with n around 252 they
## should be treated as indicative.

## ---- 1b. Diagnostic: how much dependence is volatility clustering? ----------
## A copula fitted to raw returns cannot distinguish "these assets have tail
## dependence" from "these assets were all volatile in the same month". March
## 2020 raises volatility in every name at once, and the rank transform records
## that as copula dependence.
##
## Crude quantification: compare the average pairwise Kendall's tau on the full
## in-sample period against the same quantity with the crisis window excluded.
## A large gap means the fitted copula is substantially a crisis artefact.
## The proper fix is to filter each margin with a GARCH model and fit the copula
## to standardised residuals (McNeil, Frey & Embrechts, Ch. 9); that is beyond
## this brief, but the limitation should be stated in the report.

avg_tau <- function(X) {
  M <- cor(X, method = "kendall")
  mean(M[lower.tri(M)])
}
crisis_in <- ret_dates[in_idx] >= as.Date("2020-02-19") &
  ret_dates[in_idx] <= as.Date("2020-04-30")

cat("Average pairwise Kendall's tau\n")
cat("  full 2020            : ", round(avg_tau(R_in), 4), "\n", sep = "")
cat("  crisis window only   : ", round(avg_tau(R_in[crisis_in, ]), 4), "\n", sep = "")
cat("  2020 ex-crisis       : ", round(avg_tau(R_in[!crisis_in, ]), 4), "\n\n", sep = "")

## ---- 2. Copula fitting -------------------------------------------------------

## ---- 2a. Pairwise fits: the fair family comparison --------------------------
## In two dimensions every family has exactly one parameter (t has two), so
## differences in AIC reflect the SHAPE of the dependence -- symmetric,
## lower-tail-heavy, upper-tail-heavy -- rather than the number of correlations
## the model is allowed to fit. This is where the question "which family
## describes the dependence" can actually be answered.

pair_grid <- t(combn(TICKERS, 2))

## Log-likelihood computed directly as the sum of log copula densities at the
## pseudo-observations. Done explicitly rather than read off the fitted object
## so that the quantity being compared is unambiguous.
bicop_loglik <- function(u1, u2, fit) sum(log(BiCopPDF(u1, u2, obj = fit)))

fit_pair <- function(a, b) {
  u1 <- U[, a]; u2 <- U[, b]
  out <- do.call(rbind, lapply(FAM_TAUGHT, function(f) {
    ## max.df = 30 is BiCopEst's default ceiling on the t degrees of freedom.
    ## A fitted df at the ceiling means the data prefer a Gaussian-like copula
    ## and the t brings nothing -- worth checking for in the results.
    fit <- tryCatch(BiCopEst(u1, u2, family = f, method = "mle", max.df = 30),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    np <- npar_of(f)
    ll <- bicop_loglik(u1, u2, fit)
    data.frame(pair = paste(a, b, sep = "-"), stock_a = a, stock_b = b,
               family_code = f,
               family = BiCopName(f, short = FALSE),
               npar = np,
               loglik = ll,
               AIC = -2 * ll + 2 * np,
               BIC = -2 * ll + log(n) * np,
               par1 = fit$par,
               par2 = if (np > 1) fit$par2 else NA_real_,
               tau  = fit$tau,
               stringsAsFactors = FALSE)
  }))
  out
}

message("Fitting ", nrow(pair_grid), " pairs x ", length(FAM_TAUGHT), " families...")
fit_table_pairs <- do.call(rbind, lapply(seq_len(nrow(pair_grid)), function(i)
  fit_pair(pair_grid[i, 1], pair_grid[i, 2])))
rownames(fit_table_pairs) <- NULL

## Which family wins most often across the 45 pairs?
best_per_pair <- do.call(rbind, lapply(split(fit_table_pairs, fit_table_pairs$pair),
                                       function(g) g[which.min(g$AIC), ]))
rownames(best_per_pair) <- NULL

cat("\n--- Best-fitting family by AIC, counted over all 45 pairs ---\n")
print(sort(table(best_per_pair$family), decreasing = TRUE))

## Mean AIC by family: a family can win rarely but be a good all-round fit.
cat("\n--- Mean AIC by family across pairs (lower is better) ---\n")
print(round(sort(tapply(fit_table_pairs$AIC, fit_table_pairs$family, mean)), 1))

## Gaussian-vs-t comparison per pair. These are nested (t -> Gaussian as
## df -> infinity), so the AIC gap isolates the value of tail dependence alone,
## holding the correlation structure fixed. This is the cleanest single piece of
## evidence in the whole task.
gauss <- fit_table_pairs[fit_table_pairs$family_code == 1, ]
tcop  <- fit_table_pairs[fit_table_pairs$family_code == 2, ]
cmp   <- merge(gauss[, c("pair", "AIC")], tcop[, c("pair", "AIC", "par1", "par2")],
               by = "pair", suffixes = c("_gauss", "_t"))
names(cmp)[names(cmp) == "par1"] <- "rho_t"
names(cmp)[names(cmp) == "par2"] <- "df_t"
cmp$t_better_by <- cmp$AIC_gauss - cmp$AIC_t

cat("\nPairs where Student-t beats Gaussian on AIC: ",
    sum(cmp$t_better_by > 0), " of ", nrow(cmp), "\n", sep = "")
cat("Median fitted df: ", round(median(cmp$df_t, na.rm = TRUE), 2), "\n", sep = "")
cat("Pairs with fitted df at the max.df = 30 ceiling: ",
    sum(cmp$df_t > 29.9, na.rm = TRUE), "\n", sep = "")
cat("\n--- Ten pairs with the strongest preference for Student-t ---\n")
print(head(cmp[order(-cmp$t_better_by), ], 10), row.names = FALSE)

## ---- 2b. Ten-dimensional fits: regular vine ---------------------------------
## A regular vine decomposes the 10-dimensional copula density into
## d(d-1)/2 = 45 bivariate pair-copulas arranged in d-1 = 9 trees. Tree 1 holds
## 9 unconditional pair-copulas; every later tree holds pair-copulas of
## CONDITIONAL distributions.
##
## Why this replaces the exchangeable Archimedean fits: a 10-dimensional
## Clayton copula has ONE parameter shared by all 45 pairs, which for a
## portfolio spanning six sectors is an implausible restriction. A Clayton vine
## has one parameter PER EDGE, so the family restriction is tested without
## imposing that every pair share a dependence level.
##
## Fair-comparison design: the structure is selected ONCE using the
## unrestricted family set, then every restricted model is refitted on that
## SAME structure via RVineCopSelect. Selecting a new structure per family
## would confound the family comparison with the tree structure.

message("Selecting vine structure (Dissmann's maximum-spanning-tree algorithm)...")
RVM_free <- RVineStructureSelect(U, familyset = FAM_TAUGHT, type = 0,
                                 selectioncrit = "AIC", indeptest = FALSE,
                                 progress = FALSE)
M_fixed <- RVM_free$Matrix

vine_npar <- function(RVM) sum(npar_of(as.vector(RVM$family)))

vine_row <- function(label, RVM) {
  np <- vine_npar(RVM)
  ll <- RVineLogLik(U, RVM)$loglik
  data.frame(model = label, npar = np, loglik = ll,
             AIC = -2 * ll + 2 * np,
             BIC = -2 * ll + log(n) * np,
             stringsAsFactors = FALSE)
}

## Single-family vines, all on the fixed structure.
##   Gaussian vine: a vine with every pair-copula Gaussian is EXACTLY the
##     multivariate Gaussian copula, so this row is the direct replacement for
##     the unstructured Gaussian fit.
##   Student-t vine: each edge gets its own df, so this is a strict
##     GENERALISATION of the multivariate t copula (which forces a common df).
##     Worth stating precisely in the report rather than calling it "the t copula".
restricted <- list(
  "Gaussian vine"  = 1,
  "Student-t vine" = 2,
  "Clayton vine"   = 3,
  "Gumbel vine"    = 4,
  "Frank vine"     = 5
)

message("Refitting family restrictions on the fixed structure...")
vine_fits <- lapply(names(restricted), function(nm) {
  f <- restricted[[nm]]
  RVineCopSelect(U, familyset = f, Matrix = M_fixed,
                 selectioncrit = "AIC", indeptest = FALSE)
})
names(vine_fits) <- names(restricted)

vine_fits[["Mixed vine (family per edge)"]] <- RVM_free

fit_table_vine <- do.call(rbind, lapply(names(vine_fits), function(nm)
  vine_row(nm, vine_fits[[nm]])))

## Independence benchmark. The independence copula density is identically 1, so
## its log-likelihood is exactly 0 and its AIC is exactly 0 -- no fitting needed.
fit_table_vine <- rbind(fit_table_vine,
                        data.frame(model = "Independence", npar = 0, loglik = 0, AIC = 0, BIC = 0))

fit_table_vine <- fit_table_vine[order(fit_table_vine$AIC), ]
fit_table_vine$dAIC <- fit_table_vine$AIC - min(fit_table_vine$AIC)
rownames(fit_table_vine) <- NULL

cat("\n--- Ten-dimensional vine fits, ranked by AIC ---\n")
print(transform(fit_table_vine,
                loglik = round(loglik, 1), AIC = round(AIC, 1),
                BIC = round(BIC, 1), dAIC = round(dAIC, 1)),
      row.names = FALSE)

## Which families does the mixed vine actually select, and where?
## Tree 1 is the interesting row: those 9 edges are UNCONDITIONAL pair-copulas
## between actual stock returns, so their families and tail dependence are
## directly interpretable. Later trees describe conditional dependence and
## should not be read as statements about stock pairs.
## In VineCopula's lower-triangular layout, the entry in row i (below the
## diagonal) has d - i conditioning variables, so it sits in tree d - i + 1.
## Sanity check below: tree 1 must contain exactly d - 1 = 9 edges.
fam_mat <- RVM_free$family
sel <- which(fam_mat != 0, arr.ind = TRUE)

edge_info <- data.frame(
  tree   = d - sel[, "row"] + 1,
  family = vapply(fam_mat[sel], function(f) BiCopName(f, short = FALSE),
                  character(1)),
  stringsAsFactors = FALSE
)
stopifnot(sum(edge_info$tree == 1) == d - 1)

cat("\n--- Families selected by the mixed vine, by tree ---\n")
print(table(edge_info$tree, edge_info$family))

cat("\n--- Tree 1 only (unconditional pairs, directly interpretable) ---\n")
print(sort(table(edge_info$family[edge_info$tree == 1]), decreasing = TRUE))

## ---- 2c. Robustness: allowing rotated families ------------------------------
## Survival Clayton (13) has UPPER tail dependence and survival Gumbel (14) has
## LOWER tail dependence. If selection moves towards 14 and away from 4, that is
## evidence the dependence is concentrated in losses rather than gains -- which
## is the substantive question. Reported as a robustness check because the five
## taught families carry the headline result.

message("Robustness: re-selecting with rotated families allowed...")
RVM_ext <- RVineStructureSelect(U, familyset = FAM_EXTENDED, type = 0,
                                selectioncrit = "AIC", indeptest = FALSE,
                                progress = FALSE)
ext_row <- vine_row("Mixed vine (rotations allowed)", RVM_ext)
cat("\n--- Extended family set ---\n")
print(transform(ext_row, loglik = round(loglik, 1), AIC = round(AIC, 1),
                BIC = round(BIC, 1)), row.names = FALSE)
sel_e <- which(RVM_ext$family != 0, arr.ind = TRUE)
cat("Families selected (all trees):\n")
print(sort(table(vapply(RVM_ext$family[sel_e],
                        function(f) BiCopName(f, short = FALSE),
                        character(1))), decreasing = TRUE))

## ---- 3. Tail dependence ------------------------------------------------------

## Closed-form coefficients:
##   Gaussian : lambda_L = lambda_U = 0 for any rho < 1. Asymptotically
##              independent in both tails, however high the correlation.
##   Frank    : lambda_L = lambda_U = 0.
##   Clayton  : lambda_L = 2^(-1/theta), lambda_U = 0.
##   Gumbel   : lambda_L = 0, lambda_U = 2 - 2^(1/theta).
##   Student-t: lambda_L = lambda_U = 2 * pt(-sqrt((nu+1)(1-rho)/(1+rho)), nu+1),
##              symmetric and strictly positive for finite nu.
##
## BiCopPar2TailDep() returns these directly. The t formula is also written out
## explicitly and cross-checked below, since it is the one worth being able to
## reproduce from memory.

lambda_t <- function(rho, nu)
  2 * pt(-sqrt((nu + 1) * (1 - rho) / (1 + rho)), df = nu + 1)

.chk <- BiCopPar2TailDep(family = 2, par = 0.5, par2 = 5)$lower
stopifnot(abs(lambda_t(0.5, 5) - .chk) < 1e-8)
cat("\nClosed-form t tail-dependence formula cross-checked against VineCopula.\n")

## Empirical (non-parametric) estimator at threshold q:
##   lambda_L(q) = P(U <= q, V <= q) / q = C(q, q) / q
## This is a finite-q proxy for the asymptotic coefficient. With n around 252
## and q = 0.05 only ~13 observations sit below the threshold, so the estimator
## is extremely noisy -- reporting it across several q values shows whether it
## stabilises or is being driven by a couple of days.

lambda_emp    <- function(a, b, q) mean(U[, a] <= q & U[, b] <= q) / q
lambda_emp_up <- function(a, b, q) mean(U[, a] > 1 - q & U[, b] > 1 - q) / q

q_grid <- c(0.20, 0.15, 0.10, 0.05, 0.025)

## Focus pairs: within-sector and across-sector, matching Task 2.
focus <- list(c("PFE", "GILD"), c("MSFT", "AMZN"), c("AAL", "BA"),
              c("AAL", "MSFT"), c("MAR", "PFE"), c("BP", "MSFT"))

## combn() emits each pair in TICKERS order, so "AAL-MSFT" is stored as
## "MSFT-AAL" whenever MSFT precedes AAL in the data. Canonicalise the key so
## focus pairs can be written in either order.
pkey <- function(a, b) paste(c(a, b)[order(match(c(a, b), TICKERS))],
                             collapse = "-")
stopifnot(all(vapply(focus, function(pr) pkey(pr[1], pr[2]), character(1))
              %in% fit_table_pairs$pair))

tail_table <- do.call(rbind, lapply(focus, function(pr) {
  a <- pr[1]; b <- pr[2]; key <- pkey(a, b)
  g <- fit_table_pairs[fit_table_pairs$pair == key, ]
  stopifnot(nrow(g) == length(FAM_TAUGHT))
  best <- g[which.min(g$AIC), ]
  th <- BiCopPar2TailDep(family = best$family_code, par = best$par1,
                         par2 = if (is.na(best$par2)) 0 else best$par2)
  row <- data.frame(
    pair = key,
    kendall_tau = round(cor(R_in[, a], R_in[, b], method = "kendall"), 3),
    best_family = best$family,
    lambda_L_model = round(th$lower, 4),
    lambda_U_model = round(th$upper, 4),
    row.names = NULL
  )
  for (q in q_grid) row[[paste0("emp_L_", q)]] <- round(lambda_emp(a, b, q), 3)
  tf <- tcop[tcop$pair == key, ]
  row$lambda_t_implied <- round(lambda_t(tf$par1, tf$par2), 4)
  row
}))

cat("\n--- Tail dependence: model-implied vs empirical co-exceedance ---\n")
print(tail_table, row.names = FALSE)

cat("\nReminder: at q = 0.05 the empirical estimate is based on about ",
    round(0.05 * n), " observations below the threshold, so treat it as ",
    "indicative only.\n", sep = "")

## Pairwise t-copula tail dependence across all 45 pairs. Shows how much tail
## dependence the Gaussian copula sets to zero by assumption.
lam_all <- lambda_t(cmp$rho_t, cmp$df_t)
cat("\n--- Pairwise lower-tail dependence under the fitted t copulas ---\n")
print(round(summary(lam_all), 4))
cat("Gaussian copula implies exactly 0 for every one of these ",
    length(lam_all), " pairs.\n", sep = "")

## ---- 4. Why it matters: portfolio VaR under each dependence model -----------
## The interpretive claim -- that ignoring tail dependence understates portfolio
## risk -- is worth making quantitatively rather than in words. Procedure:
##   1. Fit a Student-t to each margin (consistent with Task 3's heavy-tail model).
##   2. Simulate from each fitted vine, holding the margins IDENTICAL.
##   3. Map to returns, aggregate with end-2020 weights, and read off VaR.
## Any difference between models is then attributable to dependence structure
## alone, since the marginal distributions are held fixed.

N_SIM <- 100000   # seed is set at the top, so this is reproducible

## Marginal t fits. Standardising before fitdistr keeps the optimiser well
## scaled; the df is scale-invariant and m, s rescale linearly.
marg <- lapply(seq_len(d), function(j) {
  x <- R_in[, j]; cs <- sd(x)
  f <- tryCatch(
    fitdistr(x / cs, densfun = "t",
             start = list(m = median(x / cs), s = IQR(x / cs) / 1.35, df = 5),
             lower = c(m = -Inf, s = 1e-6, df = 2.1)),
    error = function(e) NULL)
  if (is.null(f)) return(NULL)
  list(m = as.numeric(f$estimate["m"]) * cs,
       s = as.numeric(f$estimate["s"]) * cs,
       df = as.numeric(f$estimate["df"]))
})
names(marg) <- TICKERS
stopifnot(!any(vapply(marg, is.null, logical(1))))

cat("\n--- Fitted marginal t degrees of freedom ---\n")
print(round(vapply(marg, function(p) p$df, numeric(1)), 2))

## fitdistr's `s` is a SCALE, not a standard deviation: sd = s * sqrt(df/(df-2)).
## Only the scale enters the quantile function, so no conversion is needed here.
qmarg <- function(u, p) p$m + p$s * qt(u, df = p$df)

## Portfolio weights on the last in-sample day (share counts fixed since 2020,
## so weights have drifted from the initial 10% each).
P_in   <- as.matrix(prices_clean[, -1, drop = FALSE])[
  as.Date(prices_clean$date) <= IN_END, , drop = FALSE]
shares <- 1000 / P_in[1, ]
last_p <- P_in[nrow(P_in), ]
w      <- shares * last_p / sum(shares * last_p)

## RVineSim returns columns in the order of the vine's variable indexing, so
## the simulated matrix must be reordered to match TICKERS before applying the
## marginal quantile functions. Getting this wrong silently pairs the wrong
## margin with the wrong variable.
sim_from_vine <- function(RVM) {
  Us <- RVineSim(N_SIM, RVM)
  nm <- RVM$names
  ## If the fitted object carried no variable names, the vine's variable order
  ## is the column order of U, which is TICKERS. Fail loudly rather than
  ## silently pairing the wrong margin with the wrong variable.
  if (is.null(nm) || !all(TICKERS %in% nm)) {
    stopifnot(ncol(Us) == d)
    colnames(Us) <- TICKERS
  } else {
    colnames(Us) <- nm
  }
  Us[, TICKERS, drop = FALSE]
}

## One-day relative portfolio loss given a matrix of simulated log returns:
##   L = -sum_i w_i * (exp(r_i) - 1)
var_from_U <- function(Us, label) {
  Rs <- sapply(seq_len(d), function(j) qmarg(Us[, j], marg[[j]]))
  L  <- -as.numeric((exp(Rs) - 1) %*% w)
  v  <- quantile(L, ALPHAS)
  es <- sapply(seq_along(ALPHAS), function(k) mean(L[L >= v[k]]))
  data.frame(model = label, alpha = ALPHAS,
             VaR_pct = 100 * as.numeric(v),
             ES_pct  = 100 * as.numeric(es),
             row.names = NULL)
}

message("Simulating portfolio losses under each dependence model...")

sim_specs <- list(
  "Gaussian vine"               = vine_fits[["Gaussian vine"]],
  "Student-t vine"              = vine_fits[["Student-t vine"]],
  "Clayton vine"                = vine_fits[["Clayton vine"]],
  "Mixed vine (family per edge)" = vine_fits[["Mixed vine (family per edge)"]]
)

var_compare <- do.call(rbind, lapply(names(sim_specs), function(nm)
  var_from_U(sim_from_vine(sim_specs[[nm]]), nm)))

## Independence: the no-dependence benchmark, showing how much of the
## portfolio's risk comes from co-movement rather than from the margins.
## Simulated directly as independent uniforms -- no vine object required.
var_compare <- rbind(var_compare,
                     var_from_U(matrix(runif(N_SIM * d), N_SIM, d), "Independence"))

## Empirical benchmark from the actual realised in-sample losses.
L_emp <- -as.numeric((exp(R_in) - 1) %*% w)
var_compare <- rbind(var_compare,
                     data.frame(model = "Empirical (realised 2020)", alpha = ALPHAS,
                                VaR_pct = 100 * as.numeric(quantile(L_emp, ALPHAS)),
                                ES_pct = 100 * sapply(ALPHAS, function(a)
                                  mean(L_emp[L_emp >= quantile(L_emp, a)])),
                                row.names = NULL))

cat("\n--- One-day portfolio VaR / ES, identical margins, different copulas ---\n")
print(transform(var_compare[order(var_compare$alpha, var_compare$model), ],
                VaR_pct = round(VaR_pct, 3), ES_pct = round(ES_pct, 3)),
      row.names = FALSE)

## Understatement from assuming a Gaussian copula, at each level.
for (a in ALPHAS) {
  g  <- var_compare$VaR_pct[var_compare$model == "Gaussian vine"  & var_compare$alpha == a]
  t_ <- var_compare$VaR_pct[var_compare$model == "Student-t vine" & var_compare$alpha == a]
  if (length(g) && length(t_))
    cat(sprintf("At %.0f%%: Gaussian-vine VaR is %.1f%% below the t-vine VaR.\n",
                100 * a, 100 * (1 - g / t_)))
}

## ---- 5. Visual diagnostics ---------------------------------------------------

op <- par(no.readonly = TRUE)

## (a) Pseudo-observation scatter for the focus pairs. Point mass concentrating
##     in the lower-left corner is the visual signature of lower-tail dependence.
par(mfrow = c(2, 3), mar = c(4, 4, 2.6, 1))
for (pr in focus) {
  plot(U[, pr[1]], U[, pr[2]], pch = 19, cex = 0.5, col = rgb(0, 0, 0, 0.35),
       xlab = pr[1], ylab = pr[2], xlim = c(0, 1), ylim = c(0, 1),
       main = sprintf("%s-%s  (tau = %.2f)", pr[1], pr[2],
                      cor(R_in[, pr[1]], R_in[, pr[2]], method = "kendall")),
       cex.main = 0.95)
  abline(h = c(0.05, 0.95), v = c(0.05, 0.95), col = "red", lty = 3)
}
par(op)

## (b) Empirical lower-tail function against threshold, with the Gaussian
##     (zero) and fitted-t benchmarks. If the empirical curve stays well above
##     zero as q falls, the Gaussian copula is untenable.
par(mfrow = c(2, 3), mar = c(4, 4, 2.6, 1))
qq <- seq(0.02, 0.30, by = 0.01)
for (pr in focus) {
  key <- pkey(pr[1], pr[2])
  emp <- sapply(qq, function(q) lambda_emp(pr[1], pr[2], q))
  tf  <- tcop[tcop$pair == key, ]
  plot(qq, emp, type = "b", pch = 19, cex = 0.5, ylim = c(0, 1),
       xlab = "threshold q", ylab = expression(hat(lambda)[L](q)),
       main = paste(pr[1], pr[2], sep = "-"), cex.main = 0.95)
  abline(h = 0, col = "blue", lty = 2)                      # Gaussian
  if (nrow(tf)) abline(h = lambda_t(tf$par1, tf$par2), col = "red", lty = 2)
  legend("topright", c("empirical", "t copula", "Gaussian"),
         col = c("black", "red", "blue"), lty = c(1, 2, 2),
         pch = c(19, NA, NA), cex = 0.65, bty = "n")
}
par(op)

## (c) Optional: the selected first tree. Useful but costs a figure slot, so
##     only include if there is room within the page limit.
## plot(RVM_free, tree = 1, edge.labels = "family")

## ---- 6. Save -----------------------------------------------------------------

write.csv(fit_table_pairs, "results/task4_copula_fits_pairs.csv", row.names = FALSE)
write.csv(best_per_pair,   "results/task4_best_per_pair.csv",     row.names = FALSE)
write.csv(cmp,             "results/task4_gaussian_vs_t.csv",     row.names = FALSE)
write.csv(fit_table_vine,  "results/task4_vine_fits_10d.csv",     row.names = FALSE)
write.csv(tail_table,      "results/task4_tail_dependence.csv",   row.names = FALSE)
write.csv(var_compare,     "results/task4_var_by_copula.csv",     row.names = FALSE)

saveRDS(list(U = U, fit_table_pairs = fit_table_pairs,
             best_per_pair = best_per_pair, gaussian_vs_t = cmp,
             fit_table_vine = fit_table_vine, RVM_free = RVM_free,
             RVM_ext = RVM_ext, tail_table = tail_table,
             var_compare = var_compare, weights = w, marg = marg),
        "data/task4_results.rds")

message("Task 4 results written to results/ and data/task4_results.rds")