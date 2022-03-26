pp <- ARMA(1,0, 1000, phi = 0.9, x0 = 25, delta = 10)
x <- pp$x
ecf(x, k = 1:500)

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

