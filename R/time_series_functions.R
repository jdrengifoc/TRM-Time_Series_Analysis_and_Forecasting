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
plot.recf <- function(ts, lag_max = 30, alpha = 0.05){
  mine.acf <- acf(ts, lag.max = lag_max, plot = F)
  # Eliminate the lag 0, to ease the analysis.
  mine.acf$acf <- mine.acf$acf[-1]
  mine.acf$lag <- mine.acf$lag[-1]
  # Confidence intervals
  ci <- qnorm(1-alpha/2)/sqrt(length(ts))
 
  # Significant lags for EACF.
  mine.acf$"significant" <- which(abs(mine.acf$acf)>= ci)-1
  mine.pacf <- pacf(ts, lag.max = lag_max, plot = F)
 
  # Significant lags for EPACF.
  mine.pacf$"significant" <- which(abs(mine.pacf$acf)>= ci)
  df <- data.frame(
    k = 1:lag_max,
    acf = mine.acf$acf,
    pacf = mine.pacf$acf,
    ci = ci
  ) %>% mutate(mycolor_acf = ifelse(ci <= abs(acf), T, F),
               mycolor_pacf = ifelse(ci <= abs(pacf), T, F))
  # Plots.
  gg_acf <- ggplot(df, aes(k, acf)) +
    geom_hline(yintercept = c(0), color = 'black') +
    geom_segment(aes(x=k, xend=k, y=0, yend=acf, color=mycolor_acf),
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
    scale_x_continuous(breaks = seq(0, lag_max, by = lag_max/10))
  
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
    scale_x_continuous(breaks = seq(0, lag_max, by = lag_max/10))
  
  return(list('data' = df, 'plot_acf' = gg_acf, 'plot_pacf' = gg_pacf))
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
    # Confidence interval.
    # ci_lb <- confint(reg)[2,][1]
    # ci_ub <- confint(reg)[2,][2]
    
  } else {
    stop("The type of the estimated acf is not valid.")
  }
  return(list("pac" = pac, "pac_significance" = sig))
}

#' Estimate the partial autocorrelation (PAC)of a time series, i.e. only the
#' direct effect of the k lags. If the pacf is outside the confidence interval
#' and it is in red means that is not significant to the linear.
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
    k,
    acf_pearson = sapply(k, ac, x = x),
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
unit_roots <- function(ts, alpha = 0.05){
  # Initialize.
  stationary <- NULL
  p_values <- NULL
  h0 <- c('No stationary', 'No stationary', 'Stationary')
  tests <- c('URT - Augmented Dickey Fuller', 'URT - Phillips-Perron',
             'ST - KPSS')
  
  # Augmented Dickey Fuller (1984)
  adf <- tseries::adf.test(ts) # FALTA PARAM K
  p_values <- c(p_values, adf$p.value)
  stationary <- c(stationary, ifelse(last(p_values) <= alpha, 1, 0))
  # Phillips-Perron (1988)
  pp <- tseries::pp.test(ts)
  p_values <- c(p_values, pp$p.value)
  stationary <- c(stationary, ifelse(last(p_values) <= alpha, 1, 0))
  
  # Kwiatkowski, Phhillips, Schimdt and Shain (1992)
  kpss <- tseries::kpss.test(ts)
  p_values <- c(p_values, kpss$p.value)
  stationary <- c(stationary, ifelse(last(p_values) > alpha, 1, 0))
  
  return(data.frame("Tests" = tests, "H0" = h0,
                    'p value' = p_values, 'Stationary', stationary))
  
}



# Miscelaneous ------------------------------------------------------------

#' Plot xts data.
#' @param main Same as in base plot.
#' @param ylab Same as in base plot.
#' @param legend.names Same as in base plot.
plot_xts <- function(x, main = '', ylab = '', legend.names = NULL){
  plot.xts(x,
           col = 'blue4',
           main = main,
           ylab = ylab)
  addLegend(legend.loc = "topleft",
            legend.names = legend.names,
            col = NULL,
            bty = 'o',
            lty = 1,
            lwd = 2)
}

#' Box cox transformation.
#' @param X Time series.
#' @param lambda Transformation parameter.
bct <- function(X, lambda){
  if (lambda == 0){
    x_t <- log(X)
  } else{
    x_t <- (X ^ lambda - 1) / lambda
  }
  return(x_t)
}


#' Automatic Box-Cox transformation.
#' @param x Time serie.
#' @param alpha Significance of the lambda interval estimate.
#' @return List with the transformed timeseries, the optimal lambda with its
#' confidence interval and the plot of the transformed serie.
my_boxcox <- function(x, alpha = 0.05){
  # Delete nan.
  x <- na.omit(x)
  # Correct negative values with a displacement.
  if (any(x < 0)){
    x <- x - min(x) + 1
  }
  
  # Optimal lambda.
  b <- boxcox(lm(x ~ 1))
  lambda <- b$x[which.max(b$y)]
  # Confidence interval
  aux <- b$x[b$y > max(b$y) - 1/2 * qchisq(1-alpha,1)]
  ci <- aux[c(1, length(aux))]

  if (between(0,ci[1],ci[2])) {
    x <- bct(x, 0)
    msg <- 'Log transformation.'
  } else if (between(1, ci[1], ci[2])) {
    msg <- 'No transformation.'
  } else if (lambda > 1) {
    x <- bct(x, lambda)
    msg <- 'Warning. Increasing the variance. Has economic sense?'
  } else {
    x <- bct(x, lambda)
    msg <- 'Box-cox transformation.'
  }
  warning(msg)
  fig <- plot_xts(x, main = 'Box-Cox transformation',
                  legend.names = paste0('lambda = ', round(lambda,2)))
return(list('transformed' = x, 'lambda' = lambda, 'ci' = ci, 'plot' = fig))
}

#' Calculate the t_value of a stats::arima model.
#' @param model stats::arima model.
t_value <- function(model){
  se <- diag(model$var.coef)^0.5
  coef <- model$coef
  return(abs(coef) / se)
}

#' Calculate the p_value of a stats::arima model.
#' @param model stats::arima model.
p_value <- function(model){
  t_val <- t_value(model)
  p_val <- 2*(1 - pt(t_val, Inf))
  p_val[is.nan(p_val)] <- 1    # Manage NAN
  return(p_val)
}


#' Calculate the t_value of a simts::estimate model.
#' @param model simts::estimate model.
t_value_sarma <- function(model){
  se <- diag(model$mod$var.coef)^0.5
  coef <- model$mod$coef
  return(abs(coef) / se)
}

#' Calculate the p_value of a simts::estimate model.
#' @param model simts::estimate model.
p_value_sarma <- function(model){
  t_val <- t_value_sarma(model)
  p_val <- 2*(1 - pt(t_val, Inf))
  p_val[is.nan(p_val)] <- 1    # Manage NAN
  return(p_val)
}

#' Normalize data.
#' @param X Atomic vector, ts or xts.
normalize <- function(X){
  aux <- (X - min(X)) / diff(range(X))
  return(aux)
}


# Estimación y validación -------------------------------------------------

my_autoarima <- function(X_test, max_order = 3,  alpha = 0.05){
  bic <- Inf
  waitbar <- winProgressBar(title = "progress bar", min = 0,
                            max = max_order, width = 300)
  
  for (i in 0:max_order) {
    for (j in 0:max_order) {
      if (i + j > 0){
        aux <- arima(X_test, order = c(i,0,j), include.mean = F)
        
        if (all(p_value(aux) <= alpha)) {
          
          if (BIC(aux) < bic){
            bic <- BIC(aux)
            base <- aux
            order <- c(i, 0, j)
          }
        }
      }
      
      setWinProgressBar(waitbar, j,
                        title = paste( round(j/max_order*100, 0),
                                       "% done - ", i, "/", max_order))
    }
  }
  close(waitbar)
  return(base)
}


auto_sarma <- function(X_test, p = 1L, q = 0L, period = 12,
                       max_sorder = 3, alpha = 0.05){
  bic <- Inf

  n <- nrow(X_test)
  waitbar <- winProgressBar(title = "progress bar", min = 0,
                            max =   max_sorder, width = 300)
  
  for (i in 0:max_sorder) {
    for (j in 0:max_sorder) {
      if (i + j > 0){
        esp_sarma <- SARMA(ar = p, ma = q, sar = i, sma = j, s = period, sigma2 = sd(X_test))
        aux <- estimate(esp_sarma, X_test, demean = FALSE)
        
        if (all(p_value_sarma(aux) <= alpha)) {
          
          if (AIC(aux, k = log(n)) < bic){
            bic <- AIC(aux, k = log(n))
            sarma <- aux
            sorder <- c(i, j)
          }
        }
      }
      setWinProgressBar(waitbar, j,
                        title = paste( round(j/  max_sorder*100, 0),
                                       "% done - ", i, "/", max_sorder))
    }
  }
  close(waitbar)
  return(sarma)
}


auto_sarma_armax <- function(X_test, order = c(1L, 0L, 0L),
                             sorder = c(0L, 0L, 0L), s = NA,
                             EXO = NULL, include.mean = F, alpha = 0.05){
  if (!is.null(EXO)){
    bic <- Inf
    n_ex <- dim(EXO)[2]
    waitbar <- winProgressBar(title = "progress bar", min = 0,
                              max = n_ex, width = 300)
    for (k in 1:n_ex){
      combs <- combn(1:n_ex, k)
      for (i in 1:dim(combs)[2]){
        aux <- arima(X_test, order = order,
                     seasonal = list(order = sorder, period = s),
                     xreg = EXO[,combs[,i]], include.mean = include.mean)
        
        if (all(p_value(aux)<= alpha) & BIC(aux) < bic){
          bic <- BIC(aux)
          model <- aux
        }
        setWinProgressBar(waitbar, i,
                          title = paste(round(100*i/dim(combs)[2], 0),
                                        "% done - ", k, "/", n_ex))
      }
    }
    close(waitbar)
  }
  return(model)
}

auto_armax <- function(X_test, EXO = NULL, include.mean = F, alpha = 0.05){
  if (!is.null(EXO)){
    bic <- Inf
    n_ex <- dim(EXO)[2]
    waitbar <- winProgressBar(title = "progress bar", min = 0,
                              max = n_ex, width = 300)
    for (k in 1:n_ex){
      combs <- combn(1:n_ex, k)
      for (i in 1:dim(combs)[2]){
        aux <- auto.arima(X_test, xreg = EXO[,combs[,i]], allowmean = include.mean)
        
        if (all(p_value(aux)<= alpha) & BIC(aux) < bic){
          bic <- BIC(aux)
          armax <- aux
        }
        setWinProgressBar(waitbar, i,
                          title = paste(round(100*i/dim(combs)[2], 0),
                                        "% done - ", k, "/", n_ex))
      }
    }
    close(waitbar)
  }
  return(armax)
}

validation <- function(model, alpha = 0.05){
  val <- list('p_values' = p_value(model),
              'BIC' = BIC(model),
              'checkresiduals' = checkresiduals(model),
              'plot.recf' = plot.recf(model$residuals, alpha = alpha),
              'ecf' = ecf(model$residuals, alpha = alpha)
  )
  return(val)
}