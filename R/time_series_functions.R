rm(list = ls())
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


# Autocorrelation functions estimates -------------------------------------
## TO DO
# AC, KENDALL, SPEARMAN. Optimize.
# PAC, ROBUST. Optimize.
# ecf. Incorpore, default R option.

#' Plot the estimated autocorrelation function and the estimated partial auto-
#' correlation function of a time series.
#' @param ts Time series.
#' @param lag_max The maximum number of lags to be plot, by default is 30.
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

#' Estimate the autocorrelation (AC) of a time series.
#' @param x Time series.
#' @param k Order of the autocorrelation estimate. Supports non negative integer
#'  numbers. Must be lower than the length of the series.
#' @param type Estimate type of the AC, by default, is the unbiased AC of
#' Pearson \code{'pearson'}.
#' @return AC estimate.
#' @import stats
#' @example
#' k <-  1:30
#' x <- rnorm(1000)
#' sapply(k, ac, x = x)
#' @seealso https://en.wikipedia.org/wiki/Autocorrelation
ac <- function(x, k, type = 'pearson'){
  # For robustness.
  type <- tolower(type)
  
  if (type %in% c('pearson', 'p')) {
    # gamma_k / gamma_0
    rho_k <- ifelse(k==0, 1,
                    cov(x[-1:-k],x[(0:(k-1))-length(x)]) / var(x))
  } else {
    stop("The type of the estimated acf is not valid.")
  }
  return(rho_k)
}

#' Estimate the partial autocorrelation (PAC)of a time series, i.e. only the
#' direct effect of the k lags.
#' @param x Time series.
#' @param k Order of the autocorrelation estimate. Supports positive integer
#'  numbers. Must be lower than the length of the series.
#' @param alpha Significance level. Must be in the unit interval.
#' @param type PACF estimation type; by default, employs the linear regression
#' by OLS \code{lm}.
#' @return List with the PACF estimate (\code{pacf}) and the significance of
#' the estimate and the regression \code{pac_significance}.
#' model.
#' @import stats
#' @example
#' k <-  1:30
#' x <- rnorm(1000)
#' x.pacf <- unlist(sapply(k, pac, x = x)[1,])
#' x.pacf.significance <- unlist(sapply(k, pac, x = x)[2,])
#' @seealso https://en.wikipedia.org/wiki/Autocorrelation
pac <- function(x, k, alpha = 0.05, type = 'lm'){
  # For robustness.
  type <- tolower(type)
  
  if (type %in% c('lm', 'linear')) {
    n <- length(x) - k
    # Preallocate.
    X <- matrix(rep(NA, k*n), ncol = k)
    
    # For each lag.
    for (i in 1:k){
      X[,i] <- x[i+(1:n)-1]
    }
    
    # Linear regression.
    reg <- lm(x[-1:-k] ~ X)
    reg_sum <- summary(reg)
    
    # Estimate.
    pac <- reg$coefficients[2]
    # Significance.
    pac_pval <- reg_sum$coefficients[2, 4]
    f_stat <- reg_sum$fstatistic
    reg_pval <- pf(f_stat[[1]] , f_stat[[2]], f_stat[[3]], lower.tail = F)
    sig <- ifelse(pac_pval >= alpha | reg_pval >= alpha, F, T)
    
    
  } else {
    stop("The type of the estimated acf is not valid.")
  }
  return(list("pac" = pac, "pac_significance" = sig))
}

#' Estimate the partial autocorrelation (PAC)of a time series, i.e. only the
#' direct effect of the k lags.
#' @param x Time series.
#' @param k Order of the autocorrelation estimate. Supports atomic vector of 
#' positive integer numbers.
#' #' @param alpha Significance level. Must be in the unit interval.
#' @return List with a data frame that summarize the results of the estimates of
#' the ACF and the PACF, and the respectively ACF and PACF plots.
#' @import dplyr, ggplot2
#' @example
#' x <- rnorm(1000)
#' ecf(x)

ecf <- function(x, k =  1:30, alpha = 0.05){
  
  df <- data.frame(
    k, acf_pearson = sapply(k, ac, x = x),
    pacf = unlist(sapply(k, pac, x = x,  alpha = alpha)[1,]),
    ci = qnorm(1-alpha/2)/sqrt(length(x))
  ) %>% mutate(mycolor_acf = ifelse(ci <= abs(acf_pearson), T, F),
               mycolor_pacf = ifelse(ci <= abs(pacf) &
                                       unlist(sapply(k, pac, x = x,
                                                     alpha = alpha)[2,]), T, F))
  gg_acf <- ggplot(df, aes(k, acf_pearson)) +
    geom_hline(yintercept = c(0), color = 'black') +
    geom_segment(aes(x=k, xend=k, y=0, yend=acf_pearson, color=mycolor_acf),
                 size=0.75, alpha=0.9) +
    geom_point( aes(color = mycolor_acf), size=2) +
    geom_hline(yintercept = c(df$ci, -df$ci), linetype = "dashed") +
    labs(title = 'Estimated autocorrelation function', subtitle = 'Pearson',
         x = 'Lags', y = 'ACF',
         caption = TeX(sprintf(r'($\alpha = %g$.)', alpha))) +
    theme_light() +
    theme(
      legend.position = "none",
      panel.border = element_blank(),
      plot.margin = unit(c(1,1,1,1),"cm"),
      plot.title = element_text(size=16,face="bold"),
      axis.title = element_text(size=12, face='bold'),
      axis.text = element_text(size=11),
      plot.caption=element_text(size=12,colour="grey30",hjust=1),
      panel.grid.major = element_line(colour="grey70",size=0.35),
      panel.grid.minor =  element_line(colour = "grey90")
    ) +
    scale_x_continuous(breaks = seq(0, max(k), by = max(k)/10))
  
  gg_pacf <- ggplot(df, aes(k, pacf)) +
    geom_hline(yintercept = c(0), color = 'black') +
    geom_segment(aes(x=k, xend=k, y=0, yend=pacf, color=mycolor_pacf),
                 size=0.75, alpha=0.9) +
    geom_point( aes(color = mycolor_pacf), size=2) +
    geom_hline(yintercept = c(df$ci, -df$ci), linetype = "dashed") +
    labs(title = 'Estimated partial autocorrelation function',
         subtitle = 'Linear regression', x = 'Lags', y = 'PACF',
         caption = TeX(sprintf(r'($\alpha = %g$.)', alpha))) +
    theme_light() +
    theme(
      legend.position = "none",
      panel.border = element_blank(),
      plot.margin = unit(c(1,1,1,1),"cm"),
      plot.title = element_text(size=16,face="bold"),
      axis.title = element_text(size=12, face='bold'),
      axis.text = element_text(size=11),
      plot.caption=element_text(size=12,colour="grey30",hjust=1),
      panel.grid.major = element_line(colour="grey70",size=0.35),
      panel.grid.minor =  element_line(colour = "grey90")
    ) +
    scale_x_continuous(breaks = seq(0, max(k), by = max(k)/10))
  return(list('data' = df, 'plot_acf' = gg_acf, 'plot_pacf' = gg_pacf))
}
  
# Unit roots and stationarity tests -------------------------------------
##To develop
# 1. Other stationarity test (Leybourne-McCabe)
# 2. Modify all possible params of the functions.

#' Summary table of unit roots and stationatity tests.
#' @param ts Time series.
#' @return Dataframe of different units roots/stacionarity tests with their
#' associated null hypotesis (in terms of stationarity) and the p-value.
#' @import tseries
#' @example
#' ts <- ARMA(1,0,1000, phi = 0.4)
#' unit_roots(ts$x)
unit_roots <- function(ts){
  # Initialize.
  p_values <- NULL
  h0 <- c('No stationary', 'No stationary', 'Stationary')
  tests <- c('URT - Augmented Dickey Fuller', 'URT - Phillips-Perron',
             'ST - KPSS')
  
  # Augmented Dickey Fuller (1984)
  adf <- tseries::adf.test(ts) # FALTA PARAM K
  p_values <- c(p_values, adf$p.value)
  # Phillips-Perron (1988)
  pp <- tseries::pp.test(ts)
  p_values <- c(p_values, pp$p.value)
  # Kwiatkowski, Phhillips, Schimdt and Shain (1992)
  kpss <- tseries::kpss.test(ts)
  p_values <- c(p_values, kpss$p.value)
  
  return(data.frame("Tests" = tests, "H0" = h0,
                    'p value' = p_values))
  
}


