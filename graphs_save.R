dir.create("figures", showWarnings = FALSE)

## ---- Guardar figuras (ejecutar DESPUÉS del script principal) ----

# 1. Figura 1: Rebased index
png("figures/fig1_rebased_index.png", width = 1200, height = 600, res = 120)
P_idx <- sweep(P, 2, P[1, ], "/") * 100
par(mfrow = c(1, 1), mar = c(4, 4, 3, 6))
matplot(dates, P_idx, type = "l", lty = 1, lwd = 1.4,
        col = rainbow(length(TICKERS), v = 0.85),
        xlab = "", ylab = "Index (start = 100)", xaxt = "n",
        main = "Cumulative performance, rebased to 100")
axis.Date(1, at = seq(min(dates), max(dates), by = "3 months"), format = "%b %y")
abline(h = 100, col = "grey60", lty = 3)
legend("topleft", legend = TICKERS, col = rainbow(length(TICKERS), v = 0.85),
       lty = 1, lwd = 1.4, cex = 0.7, ncol = 2, bty = "n")
par(op)
dev.off()

# 2. Figura 2: Log-returns panel
png("figures/fig2_log_returns.png", width = 1200, height = 800, res = 120)
ylim_common <- range(R_log)
par(mfrow = c(5, 2), mar = c(2, 3.6, 1.8, 0.8), oma = c(0, 0, 2, 0))
for (tk in TICKERS) {
  plot(ret_dates, R_log[, tk], type = "l", col = "grey25", ylim = ylim_common,
       xlab = "", ylab = "log return", main = tk, cex.main = 0.95)
  abline(h = 0, col = "grey70")
}
mtext("Daily log returns (common scale)", outer = TRUE, cex = 1.0)
par(op)
dev.off()

# 3. Figura 3: QQ-plots
png("figures/fig3_qqplots.png", width = 1200, height = 800, res = 120)
par(mfrow = c(5, 2), mar = c(3.4, 3.4, 1.8, 0.8), mgp = c(2, 0.7, 0),
    oma = c(0, 0, 2, 0))
for (tk in TICKERS) {
  z <- as.numeric(scale(R_log[, tk]))
  qqnorm(z, main = tk, cex.main = 0.95, pch = 19, cex = 0.4,
         col = rgb(0, 0, 0, 0.45),
         xlab = "Theoretical quantiles", ylab = "Standardised return")
  abline(0, 1, col = "red", lwd = 1.3)
}
mtext("Normal QQ-plots of standardised daily log returns", outer = TRUE, cex = 1.0)
par(op)
dev.off()

# 4. Figura 4: Scatter plots (2x3)
png("figures/fig4_scatter.png", width = 1200, height = 800, res = 120)
within_pairs <- list(c("PFE", "GILD"), c("MSFT", "AMZN"), c("AAL", "BA"))
across_pairs <- list(c("AAL", "MSFT"), c("MAR", "PFE"), c("BP", "MSFT"))

scatter_pair <- function(a, b, tag) {
  x <- R_log[, a]; y <- R_log[, b]
  plot(x, y, pch = 19, cex = 0.45, col = rgb(0, 0, 0, 0.35),
       xlab = paste(a, "log return"), ylab = paste(b, "log return"),
       main = sprintf("%s vs %s (%s), rho = %.2f", a, b, tag, cor(x, y)),
       cex.main = 0.9)
  abline(h = 0, v = 0, col = "grey85")
  abline(lm(y ~ x), col = "red", lwd = 1.3)
}

par(mfrow = c(2, 3), mar = c(4, 4, 2.6, 1))
for (pr in within_pairs) scatter_pair(pr[1], pr[2], "within")
for (pr in across_pairs) scatter_pair(pr[1], pr[2], "across")
par(op)
dev.off()

message("Figures saved to figures/ directory")