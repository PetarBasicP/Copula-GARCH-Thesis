### AR(1) process with mu = 1, phi = 0.9, and sigma.e = 1
# Need to start with a y1 value, so pick y1 = mu + epsilon1

phi = 0.9

# if we change this value to be negative, the process will fluctuate about the mean
# setting phi to be 0 gives us white noise centered around mu (in this case 1)

mu = 1 
sigma.e = 1 
n.obs = 250 
y = rep(0, n.obs)
set.seed(123) 
e = rnorm(n.obs, sd=sigma.e) 
y[1] = mu + e[1]
for (i in 2:n.obs) {
  y[i] = mu + phi*(y[i-1] - mu) + e[i] 
} 
head(y, 3)
plot(y, type = 'l')

### OR

y. = mu + arima.sim(model = list(ar = 0.9), n=250, n.start=1, start.innov=0)
plot(y.)

# So by changing the phi parameter we are able to generate cycles of data. This is
# useful for business cycles.