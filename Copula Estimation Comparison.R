library(copula)
library(qrmtools)
library(qrmdata)
library(rugarch)
library(xts)
#install.packages('qqtest')
library(qqtest)
library(ADGofTest)

### 1. Working with the data #################################################################
data('SP500_const')
#View(SP500_const)
stocks <- c('XOM','WMT','AAPL','GOOGL','MSFT') #intel and qualcomm are 1st two
S <- SP500_const['2012-01-01/2015-06-30', stocks]
#View(S)
plot.zoo(S)

acf_S <- acf(S[,1], plot = FALSE)
acf(S[,1], xlim = c(1.5, max(acf_S$lag)), ylim = c(-0.2,1)) 
# All of the individual stocks have significant ACF, which needs to be removed when working
# with copulas. So first we will difference the raw prices to stabilize the mean, and 
# take the logs to remove unconditional variance and make calculations easier afterwards.

X <- -returns(S)
#X <- diff(S) cannot just do this since as companies grow you would see bigger 
# fluctuations in the price, meaning that the time series will be heteroskedastic.

plot(X[,1])

plot.zoo(X)


acf_X <- acf(X[,1], plot = FALSE)
acf(X[,1], xlim = c(1.5, max(acf_X$lag)), ylim = c(-0.2,0.2)) 
# Now we see one or two spikes above our 95% threshold, which means that the differencing
# and the the log transformations did a good job, but now we need need to check how the 
# variance behaves. So the mean and the autocorrelation were dealt with, but now we need 
# to check for conditional variance (sigma sub t squared).

acf_X2 <- acf(X[,1]^2, plot = FALSE)
acf(X[,1]^2, xlim = c(1.5, max(acf_X2$lag)))
# So here we still see multiple spikes above the 95% threshold, meaning we need to run an 
# ARMA-GARCH model. (technically could just run a GARCH model)


n <- nrow(X)
d <- ncol(X)

### 2. Fit marginal ARMA(1,1) GARCH(1,1) Models #############################################
#So we already did the differencing with the returns function and doing ARIMA
# would do it again and then bring it back to regular prices after the model is fit,
# which we don't want.

# Define ARMA(1,1)-GARCH(1,1) specification with Normal innovations
uspec_norm <- ugarchspec(variance.model = list(model = 'sGARCH', garchOrder = c(1,1)),
                    mean.model     = list(armaOrder = c(1,1)),
                    distribution.model = 'norm') # abbrev. for normal t distr.
# baseline model

garch_fits_norm <- list()
Z_norm <- matrix(NA, nrow = n, ncol = d)
colnames(Z_norm) <- stocks


uspec_norm

for(j in 1:d) {
  garch_fits_norm[[j]] <- ugarchfit(spec = uspec_norm, data = X[,j], solver = 'hybrid')
  Z_norm[,j] <- residuals(garch_fits_norm[[j]], standardize = TRUE)
  # calculates the difference between the actual data and what the arma-garch model predicted
  # also divides those residuals by the standard deviation of a specific day
}

Z_norm
plot(X[,1])
plot.zoo(Z_norm)
# here we can see that there is no more volatility clustering after running the ARMA-GARCH 
# model on the data
acf_z_norm2 <- acf(Z_norm[,1]^2, plot = FALSE)
acf(Z_norm[,1]^2, xlim = c(1.5, max(acf_z_norm2$lag)), ylim = c(-0.2, 0.2))

# Success, it is white noise, so we can fit a copula!

U_norm <- pobs(Z_norm)
# So translate this white noise into uniform data.

### Marginal ARMA-GARCH model with Student-t innovations

uspec_std <- ugarchspec(variance.model = list(model = 'sGARCH', garchOrder = c(1,1)),
                        mean.model     = list(armaOrder = c(1,1)),
                        distribution.model = 'std')

garch_fits_std <- list()
Z_std <- matrix(NA, nrow = n, ncol = d)
colnames(Z_std) <- stocks
nu_margins <- numeric(d)

for(j in 1:d) {
  garch_fits_std[[j]] <- ugarchfit(spec = uspec_std, data = X[,j])
  Z_std[,j] <- residuals(garch_fits_std[[j]], standardize = TRUE)
  nu_margins[j] <- coef(garch_fits_std[[j]])['shape']
}

plot.zoo(Z_std)

# Now to rigorously check whether the conditional variance is not correlated, we run:
acf_z_std2 <- acf(Z_std[,1]^2, plot = FALSE)
acf(Z_std[,1]^2, xlim = c(1.5, max(acf_z_std2$lag)), ylim = c(-0.2,0.2))

# Success!! Now we can fit the copula!

U_std <- pobs(Z_std)

### 3. Fit copulas & Compare Model Fit ##################################################
fit_norm <- fitCopula(normalCopula(dim = d, dispstr = 'un'), data = U_std, method = 'mpl')

fit_t <- fitCopula(tCopula(dim = d, dispstr = 'un'), data = U_std, method = 'mpl')

fit_gumb <- fitCopula(gumbelCopula(dim = d), data = U_std, method = 'mpl')
# mpl is Maximum Pseudo-Likelihood

copula_comparison <- data.frame(
  Copula = c('Gaussian', 'Student-t','Gumbel'),
  LogLik = c(logLik(fit_norm), logLik(fit_t), logLik(fit_gumb)),
  AIC = c(AIC(fit_norm), AIC(fit_t), AIC(fit_gumb)))

print(copula_comparison)

best_copula <- fit_t@copula

### 4. Monte Carlo Simulation ##########################################################

n_sim <- 76 # number of days

N_rep <- 500 # number of simulations

sim_losses_norm <- matrix(NA, nrow = n_sim, ncol = N_rep)
sim_losses_std  <- matrix(NA, nrow = n_sim, ncol = N_rep)

set.seed(271)

for(b in 1:N_rep){
  # Generate joint uniform copula shocks
  U_sim <- rCopula(n_sim, copula = best_copula)
  
  # --- Model A: Normal Margins Simulation ---
  Z_sim_norm <- qnorm(U_sim) # this is where the difference happens
  X_sim_norm <- matrix(NA, nrow = n_sim, ncol = d)
  for(j in 1:d) {
    sim_obj_n <- ugarchsim(garch_fits_norm[[j]], n.sim = n_sim, 
                           custom.dist = list(name = 'sample', distfit = matrix(Z_sim_norm[,j], ncol = 1)))
    X_sim_norm[,j] <- fitted(sim_obj_n)
  }
  sim_losses_norm[,b] <- rowSums(X_sim_norm)
  
  # --- Model B: Student-t Margins Simulation ---
  Z_sim_std <- matrix(NA, nrow = n_sim, ncol = d)
  for(j in 1:d) {
    scale_factor <- sqrt((nu_margins[j] - 2) / nu_margins[j])
    Z_sim_std[,j] <- scale_factor * qt(U_sim[,j], df = nu_margins[j])
  }
  X_sim_std <- matrix(NA, nrow = n_sim, ncol = d)
  for(j in 1:d) {
    sim_obj_t <- ugarchsim(garch_fits_std[[j]], n.sim = n_sim, 
                           custom.dist = list(name = 'sample', distfit = matrix(Z_sim_std[,j], ncol = 1)))
    X_sim_std[,j] <- fitted(sim_obj_t)
  }
  sim_losses_std[,b] <- rowSums(X_sim_std)
}

# So here we created 500 independent Monte Carlo simulation paths for 76 days in the future

#View(sim_losses_norm)
plot(X_sim_norm[,1], type = 'l')

#View(sim_losses_std)

### 5. Forecast Metrics & Visualization ################################################
#View(sim_losses_std)
mean_loss <- rowMeans(sim_losses_std)
VaR_99_norm <- apply(sim_losses_norm, 1, function(path) quantile(path, 0.99))
VaR_99_std <- apply(sim_losses_std,  1, function(path) quantile(path, 0.99))

plot(1:n_sim, mean_loss, type = "l", col = "blue", lwd = 2,
     ylim = range(c(sim_losses_norm, sim_losses_std)),
     xlab = "Forecast Horizon (Days)", ylab = "Aggregated Portfolio Loss",
     main = "GARCH-Copula 99% VaR: Normal vs. Student-t Margins")
lines(1:n_sim, VaR_99_norm, col = "red", lwd = 2)
lines(1:n_sim, VaR_99_std,  col = "orange", lwd = 2)

## This result shows us that the t distribution margins will always have higher VaR, since
# the t distribution has fatter tails, meaning that the chance of an extreme event happening
# is higher than in the normal distribution. The normal distribution here acts as a baseline.

### 6. Fitting the Actual Data of the Forecasted Period ################################

# Extract the actual out-of-sample prices for the next 76 trading days following your training set
end_date_idx <- which(index(SP500_const) == index(last(S)))
S_actual <- SP500_const[end_date_idx:(end_date_idx + n_sim), stocks]
plot(S_actual)

# Compute actual realized portfolio losses using the same loss transformation (-returns)
X_actual <- -returns(S_actual)
nrow(X_actual)
actual_losses <- rowSums(X_actual)

# Count VaR exceptions (how many times realized losses breached your 99% Student-t VaR)
exceptions_std <- sum(actual_losses > VaR_99_std)
exceptions_norm <- sum(actual_losses > VaR_99_norm)

cat("--- Backtest Results (", n_sim, " Trading Days) ---\n")
cat("Student-t Margins 99% VaR Exceptions:", exceptions_std, "\n")
cat("Normal Margins 99% VaR Exceptions:   ", exceptions_norm, "\n")

### Enhanced Plot with Realized Losses ##############################################

plot(1:n_sim, mean_loss, type = "l", col = "blue", lwd = 2,
     ylim = range(c(sim_losses_std, actual_losses, VaR_99_std)),
     xlab = "Forecast Horizon (Days)", ylab = "Aggregated Portfolio Loss",
     main = "GARCH-Copula 99% VaR vs. Realized Out-of-Sample Losses")

# Add VaR threshold lines
lines(1:n_sim, VaR_99_std, col = "orange", lwd = 2)
lines(1:n_sim, VaR_99_norm, col = "red", lwd = 2)

# Plot the actual realized loss path
lines(1:n_sim, actual_losses, col = "black", lwd = 2, lty = 1)
length(actual_losses)
# Add a legend
legend("topleft", 
       legend = c("Simulated Mean Loss", "Student-t 99% VaR", "Normal 99% VaR", "Realized Actual Loss"),
       col = c("blue", "orange", "red", "black"), 
       lwd = 2, lty = c(1, 1, 2, 1), cex = 0.8)

### 7. Comprehensive Plot: Training History + Out-of-Sample Forecast & Actuals ###

# 1. Combine training losses with out-of-sample actual losses
train_losses <- rowSums(X)
full_losses  <- c(train_losses, actual_losses)
full_dates   <- c(index(X), index(X_actual))

# 2. Create padded forecast vectors to align with the full timeline length
n_train <- length(train_losses)
plot_VaR_std   <- c(rep(NA, n_train), VaR_99_std)
plot_VaR_norm  <- c(rep(NA, n_train), VaR_99_norm)
plot_mean_loss <- c(rep(NA, n_train), mean_loss)

# 3. Plot the complete timeline
plot(full_dates, full_losses, type = "l", col = "gray40", lwd = 1.5,
     ylim = range(c(full_losses, plot_VaR_std, plot_VaR_norm), na.rm = TRUE),
     xlab = "Date", ylab = "Aggregated Portfolio Loss",
     main = "GARCH-Copula 99% VaR: Training History vs. Out-of-Sample Forecast")

# Add the out-of-sample actual loss path
lines(index(X_actual), actual_losses, col = "black", lwd = 2)

# Add the forecasted risk bands and mean loss
lines(full_dates, plot_mean_loss, col = "blue", lwd = 2)
lines(full_dates, plot_VaR_std,  col = "orange", lwd = 2)
lines(full_dates, plot_VaR_norm, col = "red", lwd = 2)











