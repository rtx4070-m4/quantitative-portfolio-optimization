# ============================================================
# STOCK PORTFOLIO ANALYSIS - NSE & BSE
# Risk-Return Optimization using Maximum Sharpe Ratio
# Author: Tanuj
# Date: September 2026
# ============================================================

# Load required packages
library(quantmod)
library(PerformanceAnalytics)
library(quadprog)
library(dplyr)
library(tidyr)
library(ggplot2)
library(zoo)

# ------------------------------------------------------------
# Stock Universe (NSE + BSE, different industries)
# ------------------------------------------------------------
stock_list <- c(
  # Banking
  "HDFCBANK.NS", "ICICIBANK.NS", "SBIN.NS", "KOTAKBANK.NS", "AXISBANK.NS",
  # IT
  "TCS.NS", "INFY.NS", "WIPRO.NS", "HCLTECH.NS", "TECHM.NS",
  # Energy
  "RELIANCE.NS", "ONGC.NS", "BPCL.NS", "IOC.NS",
  # FMCG
  "HINDUNILVR.NS", "ITC.NS", "NESTLEIND.NS", "BRITANNIA.NS",
  # Pharma
  "SUNPHARMA.NS", "DRREDDY.NS", "CIPLA.NS", "DIVISLAB.NS",
  # Auto
  "MARUTI.NS", "M&M.NS", "BAJAJ-AUTO.NS", "HEROMOTOCO.NS",
  # Metals & Others
  "TATASTEEL.NS", "JSWSTEEL.NS", "HINDALCO.NS", "LT.NS", 
  "BAJFINANCE.NS", "ADANIENT.NS"
)

# ------------------------------------------------------------
# Function 1: Download Data
# ------------------------------------------------------------
download_stock_data <- function(tickers, years_back) {
  end_date   <- Sys.Date()
  start_date <- end_date - years_back * 365
  
  cat("Period:", as.character(start_date), "to", as.character(end_date), "\n")
  
  price_list <- list()
  success_count <- 0
  
  for (ticker in tickers) {
    tryCatch({
      data <- getSymbols(ticker, src = "yahoo", from = start_date, to = end_date,
                         auto.assign = FALSE, warnings = FALSE)
      adj <- Ad(data)
      colnames(adj) <- ticker
      
      if (nrow(adj) > 100) {
        price_list[[ticker]] <- adj
        success_count <- success_count + 1
        cat("✓", ticker, "- rows:", nrow(adj), "\n")
      }
    }, error = function(e) {
      cat("✗", ticker, "- failed\n")
    })
  }
  
  prices <- do.call(merge, price_list)
  prices <- na.omit(prices)
  
  cat("\nSuccessfully downloaded:", success_count, "stocks |", 
      ncol(prices), "stocks after cleaning |", nrow(prices), "days\n\n")
  
  return(prices)
}

# ------------------------------------------------------------
# Function 2: Calculate Returns
# Mathematical Formula: Daily Return = (P_t - P_{t-1}) / P_{t-1}
# ------------------------------------------------------------
calculate_returns <- function(prices) {
  returns <- CalculateReturns(prices, method = "discrete")
  returns <- na.omit(returns)
  return(returns)
}

# ------------------------------------------------------------
# Function 3: Analyze Stocks (Return, Risk, Sharpe)
# ------------------------------------------------------------
# mu     = mean(daily returns) × 252          → Annualized Return
# sigma  = sd(daily returns) × √252           → Annualized Risk
# Sharpe = (mu - Rf) / sigma
# ------------------------------------------------------------
analyze_stocks <- function(returns, risk_free_rate = 0.06) {
  expected_return <- colMeans(returns) * 252
  risk            <- apply(returns, 2, sd) * sqrt(252)
  sharpe          <- (expected_return - risk_free_rate) / risk
  risk_return     <- expected_return / risk
  daily_sharpe    <- (colMeans(returns) - risk_free_rate/252) / apply(returns, 2, sd)
  
  summary_df <- data.frame(
    Stock             = names(expected_return),
    Annual_Return_Pct = round(expected_return * 100, 2),
    Annual_Risk_Pct   = round(risk * 100, 2),
    Sharpe_Ratio      = round(sharpe, 3),
    Risk_Return_Ratio = round(risk_return, 3),
    Daily_Sharpe      = round(daily_sharpe, 4)
  )
  
  summary_df <- summary_df[order(-summary_df$Sharpe_Ratio), ]
  rownames(summary_df) <- NULL
  return(summary_df)
}

# ------------------------------------------------------------
# Function 4: Optimize Portfolio (Maximum Sharpe Ratio)
# Mathematical Model:
# Maximize  (w'μ - Rf) / √(w'Σw)
# Subject to: Σw = 1 , w ≥ 0
# Solved using Quadratic Programming
# ------------------------------------------------------------
optimize_portfolio <- function(returns, risk_free_rate = 0.06) {
  mu    <- colMeans(returns) * 252
  Sigma <- cov(returns) * 252
  n     <- length(mu)
  
  Dmat <- 2 * Sigma
  dvec <- rep(0, n)
  Amat <- cbind(rep(1, n), diag(n))
  bvec <- c(1, rep(0, n))
  
  sol <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 1)
  
  weights <- sol$solution
  names(weights) <- names(mu)
  weights[weights < 0.005] <- 0
  weights <- weights / sum(weights)
  
  port_ret  <- sum(weights * mu)
  port_risk <- as.numeric(sqrt(t(weights) %*% Sigma %*% weights))
  port_shp  <- (port_ret - risk_free_rate) / port_risk
  
  list(
    weights = round(weights, 4),
    return  = round(port_ret * 100, 2),
    risk    = round(port_risk * 100, 2),
    sharpe  = round(port_shp, 3)
  )
}

# ------------------------------------------------------------
# Main Analysis: Run for 10yr, 5yr, 3yr, 1yr
# ------------------------------------------------------------
time_frames <- c(10, 5, 3, 1)
all_results <- list()

for (yrs in time_frames) {
  cat("\n############################################################\n")
  cat("### ANALYZING", yrs, "YEAR PERIOD\n")
  cat("############################################################\n")
  
  prices  <- download_stock_data(stock_list, years_back = yrs)
  returns <- calculate_returns(prices)
  
  stock_summary <- analyze_stocks(returns)
  print(head(stock_summary, 10))
  
  top5 <- head(stock_summary$Stock, 5)
  cat("\n>>> Selected Top 5:", paste(top5, collapse = " | "), "\n")
  
  top5_returns <- returns[, top5, drop = FALSE]
  opt <- optimize_portfolio(top5_returns)
  
  cat("\nOptimized Weights:\n")
  print(opt$weights)
  cat("\nPortfolio Return:", opt$return, "% | Risk:", opt$risk, 
      "% | Sharpe:", opt$sharpe, "\n")
  
  all_results[[paste0(yrs, "yr")]] <- list(
    summary = stock_summary,
    top5    = top5,
    weights = opt$weights,
    metrics = c(Return = opt$return, Risk = opt$risk, Sharpe = opt$sharpe)
  )
}

# ------------------------------------------------------------
# Weight Comparison Table
# ------------------------------------------------------------
all_stocks <- unique(unlist(lapply(all_results, function(x) names(x$weights))))
weight_table <- data.frame(Stock = all_stocks)

for (tf in names(all_results)) {
  w <- all_results[[tf]]$weights
  weight_table[[tf]] <- w[match(all_stocks, names(w))]
}
weight_table[is.na(weight_table)] <- 0
print(weight_table)

# ------------------------------------------------------------
# Bar Chart of Weights
# ------------------------------------------------------------
weight_long <- weight_table %>%
  pivot_longer(cols = -Stock, names_to = "Period", values_to = "Weight")

ggplot(weight_long, aes(x = Stock, y = Weight, fill = Period)) +
  geom_col(position = "dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Optimized Portfolio Weights Across Time Frames",
       y = "Weight", x = "Stock") +
  scale_y_continuous(labels = scales::percent)

