y <- c( 12.21,  9.57, 12.50, 10.55, 10.22, 10.67, 10.24, 11.43,  9.86, 12.25, 
        11.95,  9.41,  9.48,  9.48, 10.83, 11.84, 12.60, 10.04, 11.97, 10.87, 
        12.58, 12.01,  9.31, 10.38, 11.49, 12.32,  9.62, 12.42, 10.19, 11.78, 
        10.02, 12.40, 11.59, 10.05, 12.15, 10.76,  8.92,  9.75, 12.53,  9.41, 
        10.13, 10.24, 10.85, 12.05,  9.95, 12.00,  9.80, 11.41, 10.01, 11.88, 
        9.92, 11.50,  9.31, 12.22, 12.17, 11.40, 11.87, 10.47, 12.43, 10.37, 
        11.76, 10.58, 10.29,  9.02, 10.68, 10.29, 10.04, 12.20,  9.30, 12.13, 
        12.73, 11.68, 10.43,  9.78, 10.70, 12.33, 11.52, 10.42, 10.16,  9.64, 
        12.33, 11.97, 11.34, 12.78,  9.59,  9.19, 11.07, 12.24, 12.90,  9.61, 
        9.80,  9.45, 10.14,  9.66, 10.44, 12.09,  9.78, 11.94, 11.86, 10.15)

# summary statistics of the sample data
summary(y)
var(y); sd(y)
# plot histogram
hist(y, breaks=10, main="Histogram of amount of daily sales")

###############################################

B <- 5000 #Number of iterations to run burn in
n <- 100 # Number of samples
m <- 10 # thinning

# Initialization
alpha <- vector("numeric", B)
gamma <- vector("numeric", B)
psi <- vector("numeric", B)
p <- vector("numeric", B)

alpha[1] = 0
gamma[1] = 0
psi[1] = 1
p[1] = 0.4
# simulate labels
dis = rbinom(length(y), 1, p[1])

#Gibbs sampler
for (iter in 1:(B-1)){
  
  # assign labels
  y1 = y[dis == 1]
  y2 = y[dis == 0]
  
  # Update alpha
  meanalpha = (psi[iter]*(sum(y)-length(y1)*gamma[iter]))/(1/100+length(y)*psi[iter])
  varalpha = 1/(length(y)*psi[iter]+1/100)
  alpha[iter+1] <- rnorm(1, mean = meanalpha, sd=sqrt(varalpha))
  
  # update gamma
  meangamma = (psi[iter]*(sum(y1)-length(y1)*alpha[iter+1]))/(1/100+psi[iter]*length(y1))
  vargamma = 1/(length(y1)*psi[iter]+1/100)
  gamma[iter+1] <- rnorm(1, mean = meangamma, sd=sqrt(vargamma))
  
  # Update psi
  R0 <- y1 - alpha[iter+1] - gamma[iter+1]
  R1 <- y2 - alpha[iter+1]
  alphapsi = length(y)/2 + 0.01
  betapsi = 0.01 + sum(R0^2)/2 + sum(R1^2)/2
  psi[iter+1] <- rgamma(1, alphapsi, betapsi)
  
  # Update p
  alphap <- length(y1) + 0.4
  betap <- length(y2) + 0.6
  p[iter+1] = rbeta(1, alphap, betap)
  
  # Update d
  for(j in 1:length(y)) {
    R2 <- y[j] - alpha[iter+1] - gamma[iter+1]
    R3 <- y[j] - alpha[iter+1]
    pc = p[iter+1]*exp(-R2^2/2*psi[iter+1])/
      (p[iter+1]*exp(-R2^2/2*psi[iter+1])+
         (1-p[iter+1])*exp(-R3^2/2*psi[iter+1]))
    
    # simulate new d using updated pc
    dis[j] = rbinom(1, 1, pc)  
  }
}

# visualise convergence
par(mfrow=c(2,2))
plot(seq(1,B,1), alpha, type="l", xlab="Iterations", main="Convergence of alpha")
plot(seq(1,B,1), gamma, type="l", xlab="Iterations", main="Convergence of gamma")
plot(seq(1,B,1), psi, type="l", xlab="Iterations", main="Convergence of psi")
plot(seq(1,B,1), p, type="l", xlab="Iterations", main="Convergence of p")

# check distribution
hist(alpha[seq(B-(m*n)+1,B,m)], xlab="alpha", main="Histogram of alpha")
hist(gamma[seq(B-(m*n)+1,B,m)], xlab="gamma", main="Histogram of gamma")
hist(psi[seq(B-(m*n)+1,B,m)], xlab="psi", main="Histogram of psi")
hist(p[seq(B-(m*n)+1,B,m)], xlab="p", main="Histogram of p")

# check means and variances
mean(alpha); var(alpha)
mean(gamma); var(gamma)
mean(psi); var(psi)
mean(p); var(p)
