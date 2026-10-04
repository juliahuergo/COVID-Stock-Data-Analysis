## =============================================================================
## ME317 Project 2026 -- Build a clean daily adjusted-price panel
##
## Output: prices_clean -- data.frame with `date` followed by one adjusted-price
##         column per ticker, in the order TPR, PFE, GILD, MSFT, AMZN, MAR,
##         AAL, BA, BP, AXP.
##
## NOTE: not executed here. Run each section and read the printed diagnostics
##       before using the output.
## =============================================================================

## ---- 0. Reproducibility ------------------------------------------------------

set.seed(12345)   # no randomness in this script, but fixed here so that any
# simulation sourced afterwards inherits the same seed

if (!requireNamespace("quantmod", quietly = TRUE)) {
  stop("Install quantmod first:  install.packages('quantmod')")
}

library(quantmod)   # getSymbols(), Ad()
library(xts)        # merge.xts, index(), coredata()
library(zoo)        # na.locf()

## ---- 1. Parameters -----------------------------------------------------------

## Ticker order here defines the column order of the final output.
TICKERS <- c("TPR", "PFE", "GILD", "MSFT", "AMZN",
             "MAR", "AAL", "BA", "BP", "AXP")

START <- as.Date("2020-01-01")
END   <- as.Date("2021-12-31")

## ---- 2. Download -------------------------------------------------------------
## Cached to disk on the first run. Yahoo rewrites its adjusted-close history
## retroactively whenever a new dividend or split occurs, so re-downloading
## later would silently change your numbers. Cite the download date in the
## report -- that, plus the cached file, is what makes the analysis reproducible.
##
## A slightly wider window is requested and then subset, because the exact
## inclusivity of getSymbols()'s `to` argument depends on what Yahoo returns for
## the boundary date.

dir.create("data", showWarnings = FALSE)
raw_path <- "data/raw_prices.rds"

if (!file.exists(raw_path)) {
  
  raw    <- list()
  failed <- character(0)
  
  for (tk in TICKERS) {
    message("Downloading ", tk, " ...")
    obj <- tryCatch(
      getSymbols(tk, src = "yahoo",
                 from = START - 10, to = END + 10,
                 auto.assign = FALSE),
      error = function(e) {
        warning("Download failed for ", tk, ": ", conditionMessage(e),
                call. = FALSE)
        NULL
      }
    )
    if (!is.null(obj)) raw[[tk]] <- obj else failed <- c(failed, tk)
    Sys.sleep(0.5)   # avoids intermittent Yahoo rate limiting
  }
  
  if (length(failed) > 0) {
    warning("Did not download: ", paste(failed, collapse = ", "), call. = FALSE)
  }
  
  attr(raw, "downloaded_on") <- Sys.Date()
  saveRDS(raw, raw_path)
  
} else {
  raw <- readRDS(raw_path)
  message("Loaded cached download from ", attr(raw, "downloaded_on"))
}

stopifnot(length(raw) > 0)

## ---- 3. Assemble the wide panel ---------------------------------------------
## Ad() extracts the adjusted close, which already accounts for splits and
## dividends. Cl() would show a split as a large fake one-day loss.

prices_raw <- do.call(merge, lapply(raw, Ad))   # outer join: union of all dates
colnames(prices_raw) <- names(raw)

## Enforce the requested column order (skipping any ticker that failed).
prices_raw <- prices_raw[, TICKERS[TICKERS %in% colnames(prices_raw)]]

## Restrict to the analysis window.
prices_raw <- prices_raw[paste0(format(START), "/", format(END))]

message("Raw panel: ", nrow(prices_raw), " dates x ", ncol(prices_raw), " tickers")

## ---- 4. Anomaly handling: zero or negative prices ---------------------------
## An adjusted close can never legitimately be <= 0. Such an observation is
## removed (set to NA) BEFORE filling, so it is repaired by the same
## carry-forward step that handles genuine gaps. Doing this after filling would
## preserve the corrupt value.

n_bad <- sum(prices_raw <= 0, na.rm = TRUE)
if (n_bad > 0) {
  warning(n_bad, " non-positive price(s) removed and carried forward",
          call. = FALSE)
  prices_raw[prices_raw <= 0] <- NA
} else {
  message("No zero or negative prices found.")
}

## ---- 5. Missing-value audit (before cleaning) -------------------------------

na_before <- colSums(is.na(prices_raw))
print(data.frame(ticker      = colnames(prices_raw),
                 n_missing   = as.integer(na_before),
                 pct_missing = round(100 * na_before / nrow(prices_raw), 3),
                 row.names   = NULL))

## Record which cells are synthetic before repairing them. A carried-forward
## price creates an artificial zero return, which biases volatility downward and
## pairwise correlations toward zero -- both matter for the later risk and
## dependence sections, so it is worth being able to quantify.
filled_flag <- is.na(prices_raw)

## ---- 6. Fill: LOCF, then backward fill for any leading gap ------------------
## Step 1: last observation carried forward.
## Step 2: any NA still remaining must sit at the start of a series (nothing
##         precedes it), so fill it from the next available price.

prices <- na.locf(prices_raw, na.rm = FALSE)                    # forward
n_still_na <- sum(is.na(prices))
if (n_still_na > 0) {
  message("Backward-filling ", n_still_na, " leading observation(s).")
  prices <- na.locf(prices, fromLast = TRUE, na.rm = FALSE)     # backward
}

stopifnot(!anyNA(prices))

## ---- 7. Final output: prices_clean ------------------------------------------

prices_clean <- data.frame(date = index(prices),
                           coredata(prices),
                           check.names = FALSE,
                           row.names   = NULL)

## ---- 8. Summary --------------------------------------------------------------

cat("\n=============== prices_clean SUMMARY ===============\n")
cat("Trading days     :", nrow(prices_clean), "\n")
cat("Date range       :", format(min(prices_clean$date)), "to",
    format(max(prices_clean$date)), "\n")
cat("Tickers          :", paste(colnames(prices_clean)[-1], collapse = ", "), "\n")
cat("Filled by LOCF   :", sum(filled_flag), "observations\n")
cat("----------------------------------------------------\n")
cat("Missing values per stock after cleaning:\n")
print(colSums(is.na(prices_clean[, -1, drop = FALSE])))
cat("====================================================\n\n")

str(prices_clean)
head(prices_clean, 3)
tail(prices_clean, 3)

## ---- 9. Corporate-action log (for the report's data section) ----------------
## Ad() already adjusts for splits and dividends, so no manual correction is
## needed -- but the brief asks you to DOCUMENT which actions occurred.
## Differencing raw-close returns against adjusted-close returns isolates them:
## a split appears as one large discrepancy, dividends as small recurring ones.
##
## Within a 2020-2021 window, check Pfizer in particular: the Upjohn/Viatris
## spin-off falls in November 2020, days from the vaccine efficacy announcement.
## Yahoo is less reliable on spin-offs than on splits, so confirm the dates
## against the company's investor-relations page rather than trusting this table.

corp_actions <- do.call(rbind, lapply(names(raw), function(tk) {
  x <- raw[[tk]][paste0(format(START), "/", format(END))]
  d <- na.omit(diff(log(Cl(x))) - diff(log(Ad(x))))
  hits <- d[abs(d) > 0.02]                       # threshold: tune and justify
  if (nrow(hits) == 0) return(NULL)
  data.frame(ticker = tk,
             date   = as.character(index(hits)),
             gap    = round(as.numeric(hits), 4),
             row.names = NULL)
}))
if (!is.null(corp_actions)) print(corp_actions)

## ---- 10. Save ----------------------------------------------------------------

saveRDS(list(prices       = prices,          # xts version, for modelling
             prices_clean = prices_clean,    # data.frame version
             filled_flag  = filled_flag,
             tickers      = colnames(prices),
             downloaded_on = attr(raw, "downloaded_on")),
        "data/prices_clean.rds")

write.csv(prices_clean, "data/prices_clean.csv", row.names = FALSE)

message("Saved data/prices_clean.rds and data/prices_clean.csv")