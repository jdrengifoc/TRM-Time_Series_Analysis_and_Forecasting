rm(list = ls())
library(fUnitRoots)
# AR(p) -------------------------------------------------------------------
#' Standard error AR(p) series, with the parameterization and traditional
#' transformed series.
#' 
#' @param p Order of the AR. Must be a integer number.
#' @param n Sample size of the series.
#' @param phi Real coefficients of the AR. Predefined, are random numbers
#' between -1 and 1. Must be of length p.
#' @param x0 p initial values of the series. Predefined, are random numbers
#' between -1 and 1. Must be of length p.
#' @param alpha Tendency real parameter. Real number.
#' @return List that has the original AR(\code{p}) series with parameterization
#' \code{phi} and \code{alpha}, as well as an associated series. This series are:
#' x.tend (tendency stationary), x.int (differentiation stationary) and x.int.tend
#' (tendency and differentiation stationary)
#' @example 
#' # AR(1) with 1000 sample sizes.
#' p <- 1; n <- 1000
#' phi <- c(0.2, 0.4); alpha <- 0.03
#' x0 <- c(0.3, 0.5)
#' ar <- AR(p, n, phi, x0, alpha)
#' plot(ar$x)
#' plot(ar$x.int)
#' @example 
#' # AR(4) with 1000 sample sizes and default parameters.
#' ar <- AR(4, 1000)
#' plot(ar$x)
#' plot(ar$x.tend)
AR <- function(p, n, 
               phi = runif(p, min = -1, max = 1),
               x0 = runif(p, min = -1, max = 1),
               alpha = 0.01, delta = 0){
  # Initiallize series.
  x <- c(x0, rep(NA, n-p))
  
  # AR simulation for each observation.
  for (i in (p+1):n){
    I <- (i-p):(i-1)
    x[i] <- delta + sum(x[I]*phi) + rnorm(1)
  }
  # Organize results.
  results <- list(x = x,
                  params = list(phi = phi, alpha = alpha, delta = delta),
                  x.int = diffinv(x),
                  x.tend = alpha*1:n + x,
                  x.int.tend = alpha*1:(n+1) + diffinv(x)
  )
  return(results)
}

# ARMA(p,q) ---------------------------------------------------------------

#' Standard error ARMA(p, q) series, with the parameterization and traditional
#' transformed series.
#' 
#' @param p Order of the AR. Must be a integer number.
#' @param q Order of the MA. Must be a integer number.
#' @param n Sample size of the series.
#' @param phi AR coefficients. Predefined, are random numbers between -1/p and
#' 1/p, in the urge of guarantee stationarity. Must be of length p.
#' @param theta MA coefficients. Predefined, are random numbers between -1 and
#' 1. Must be of length p.
#' @param x0 p initial values of the series. Predefined, are random numbers
#' between -1 and 1. Must be of length p.
#' @param alpha Tendency real parameter.
#' @return List where \code{x} is the original ARMA(\code{p}, \code{q}) series,
#' \code{error} are the residuals, \code{params} are the parameters of the model
#' \code{phi}, \code{theta} \code{alpha} and \code{delta}. Furthermore, the list
#' has associated series. This series are: x.tend (tendency stationary), x.int 
#' (differentiation stationary) and x.int.tend (tendency and differentiation
#' stationary)
#' @example 
#' # AR(1) with 1000 sample sizes.
#' ar <- ARMA(1, 0, 1000, theta = 1)
#' plot(ar$x, type = 'l')
#' acf(ar$x)
#' pacf(ar$x)
#' @example 
#' # MA(4) with 100 sample sizes and default parameters.
#' ma <- ARMA(0, 4, 100)
#' plot(ma$x, type = 'l')
#' acf(ma$x)
#' pacf(ma$x)
#'  @example 
#' # ARMA(100, 100) with 1000 sample sizes and default parameters.
#' arma <- ARMA(100, 100, 1000)
#' plot(arma$x, type = 'l')
#' acf(arma$x)
#' pacf(arma$x)
ARMA <- function(p, q, n, 
                 phi = runif(p, min = -1, max = 1)/p,
                 theta = runif(q, min = -1, max = 1),
                 x0 = runif(p, min = -1, max = 1), alpha = 0.01, delta = 0)
{
  # Initiallize series.
  ifelse(p > 0, phi <- phi, phi <- 0)
  ifelse(q > 0, theta <- theta, theta <- 0)
  x <- c(x0, rep(0, n-p))
  # Asume initial errors equals to ZERO.
  fix_eps <- ifelse(q > p, q-p, 0)
  eps <- c(rep(0, fix_eps), rnorm(n))
  
  # AR simulation for each observation.
  for (i in (p+1):n){
    I_p <- (i-p):(i-1)
    I_q <- (i-q):(i-1) + fix_eps
    x[i] <- delta + sum(phi*x[I_p]) + eps[i+fix_eps] + sum(theta*eps[I_q]) 
  }
  # Organize results.
  results <- list(x = x, error = eps,
                  params = list(phi = phi, theta = theta,
                                alpha = alpha, delta = delta),
                  x.int = diffinv(x),
                  x.tend = alpha*1:n + x,
                  x.int.tend = alpha*1:(n+1) + diffinv(x)
  )
  return(results)
}


# Estimated autocorrelation functions -------------------------------------

#' Plot the estimated autocorrelation function and the estimated partial auto-
#' correlation function of a time series.
#' @param ts Time series.
#' @param lad_max The maximum number of lags to be plot, by default is 30.
#' @param sig_min Significance level of the estimations, by default is 0.05. All
#' the lags that are greater that it in absolute value are significant.
#' @return List with the acf and pacf information, we 
#' @import stats::acf, stats::pacf.
#' @example
#' Ease the idenfification of the order of AR.
#' ts <- ARMA(1,0,1000, phi = 0.4)
#' ecf <- plot.ecf(ts$x)
#' ecf$pacf$significant
#' @seealso  https://stats.stackexchange.com/questions/211628/how-is-the-confidence-interval-calculated-for-the-acf-function
plot.ecf <- function(ts, lag_max = 30, sig_min = 0.05){
  mine.acf <- acf(ts, lag.max = lag_max, plot = F)
  # Eliminate the lag 0, to ease the analysis.
  mine.acf$acf[1] <- NA
  mine.acf$lag[1] <- NA
  # Confidence intervals
  ci <- qnorm(1-sig_min/2)/sqrt(length(ts))
  
  # Plots.
  par(mfrow =  c(2, 1))
  plot(mine.acf, ci = 1-sig_min, main = "Función de autocorrelación estimada")
  # Significant lags for EACF.
  mine.acf$"significant" <- which(abs(mine.acf$acf)>= ci)-1
  mine.pacf <- pacf(ts, lag.max = lag_max, plot = F)
  plot(mine.pacf, ci = 1-sig_min, main = "Función de autocorrelación parcial estimada")
  # Significant lags for EPACF.
  mine.pacf$"significant" <- which(abs(mine.pacf$acf)>= ci)
  
  return(list("acf" = mine.acf, "pacf" = mine.pacf))
}
ar <- AR(2, 1000, delta = 1)
plot(ar$x, type = 'l')
mean(ar$x)


unit_roots <- function(ts){
  # Initialize.
  p_values <- NULL
  h0 <- c('No stationary', 'No stationary', 'Stationary')
  tests <- c('Augmented Dickey Fuller', 'Phillips-Perron', 'KPSS')
  
  adf <- tseries::adf.test(ts) # FALTA PARAM K
  p_values <- c(p_values, adf$p.value)
  
  pp <- tseries::pp.test(ts)
  p_values <- c(p_values, pp$p.value)
  
  kpss <- tseries::kpss.test(ts)
  p_values <- c(p_values, kpss$p.value)
  
  return(data.frame("Tests" = tests, "H0" = h0,
                    'p value' = p_values))
  
}
