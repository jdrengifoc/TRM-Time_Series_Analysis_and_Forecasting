rm(list = ls())
base_dir <- 'C:/Users/USUARIO_PC/OneDrive - Universidad EAFIT/2022-1/Econometrics/Windows/TRM-Time_Series_Analysis_and_Forecasting/'
source('requirements.R')
source('time_series_functions.R')

# Data -----------------------------------------------------------------------
DATA <- read_excel(paste0(base_dir,'Datasets/datasets.xlsx'), sheet = 'All',
                   col_types = c('date', rep("numeric", 7))) %>% 
  filter( Fecha >= as.Date('1999-01-01'))

data_xts <- xts(DATA[,2:8],
                order.by = as.Date(DATA$Fecha))

data_xts_diff <- data_xts[-1,-7]

for (i in 1:dim(data_xts_diff)[2]){
  data_xts_diff[,i] <- na.omit(diff(data_xts[,i]))
}

# Numeric statistics.
t(describe(data_xts))
print(xtable(t(describe(data_xts))), include.rownames = TRUE)
# TRM visualization.
plot_xts(X_test, ylab = 'COP / USD', legend.names = 'TRM Estacionaria')

# Statiocionarity ---------------------------------------------------------
# Testing.
x <- diff(data_xts$TRM)
unit_roots(x)
bc <- my_boxcox(x, 0.01)
unit_roots(bc$transformed)
# Selection.
X_test <- na.omit(diff(data_xts$TRM))
X_test <- X_test - mean(X_test)

# Identificación ----------------------------------------------------------
ecf(x, alpha = 0.01)
plot.recf(x, alpha = 0.01)

# ARMA --------------------------------------------------------------------

base <- my_autoarima(X_test, 9)
base1 <- auto.arima(X_test)


stargazer(base, base)

# ARIMAX ------------------------------------------------------------------
sarma_armax <- auto_sarma_armax(X_test, order = c(1L, 0L, 0L), sorder = c(2L, 0L, 1L),
                 s = 9, EXO = data_xts_diff[,2:6])

armax <- auto_armax(X_test, data_xts_diff[, 2:6])

# SARMA -------------------------------------------------------------------

sarma <- auto_sarma(X_test, p = 1L, q = 0L, period = 9, max_sorder = 9)

# Validation --------------------------------------------------------------
(model <- arima(X_test, order = c(1L, 0L, 0L),
               seasonal = list(order = c(0L, 0L, 0L), period = NA),
               xreg = data_xts_diff$OIL,
               include.mean = F))

validation(model)


