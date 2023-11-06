#install.packages("glue")
library(glue)

#loglikelihood function for gamma distribution
loglik.gamma = function(y, a, l) {
  n = length(y)
  r = n*a*log(l) - l*sum(y) + (a-1)*sum(log(y)) - n*log(gamma(a))
  return(r)
}

#score function for gamma distribution
score.gamma = function(y, a, l) {
  n = length(y)
  alpha = n*log(l) + sum(log(y)) - n*digamma(a)
  lambda = (n*a)/l - sum(y)
  r = c(alpha, lambda)
  
  names(r) = c("alpha", "lambda")
  return(r)
}

#fisher information matrix for gamma distribution
inf.gamma = function(y, a, l) {
  n = length(y)
  alpha = n*trigamma(a)
  cov = -n/l
  lambda = (n*a)/l^2
  r = matrix(c(alpha, cov, cov, lambda), 2, 2, byrow = T)
  colnames(r) = c("alpha", "lambda")
  rownames(r) = c("alpha", "lambda")
  return(r)
}


#let y be a sample from a gamma distribution with shape=alpha and rate=lambda
#let score.gamma be the score vector for alpha and lambda
#let inf.gamma be the fisher information matrix for alpha and lambda

#t0 is the vector of starting values of the algorithm
#t0 elements are chosen by estimates for alpha and lambda found by method of moments
#it can be shown that these estimates are
#alpha=(mean(y))^2/var(y)
#lambda=mean(y)/var(y)
#then t1 = t0 + solve(inf.gamma)%*%score.gamma

#let new t0=t1, start again until deltadev=abs(t1-t0)<eps

#first value of deltadev is initialized as deltadev>eps
#eps is fixed as small as you like

newton.gamma = function(y, eps = 0.00000001, plots=F) {
  #estimates with method of moments
  amm = (mean(y))^2/(var(y))
  lmm = mean(y)/var(y)
  t0 = c(amm, lmm)
  
  #matrix for saving results
  r = matrix(NA, 2, 1)
  rownames(r) = c("alpha", "lambda")
  colnames(r) = glue("step 0")
  r[, 1] = t0
  r
  
  
  #initialization
  
  deltadev = c(2, 2)
  i = 1
  
  #Algorithm
  while (abs(deltadev)[1] > eps & abs(deltadev)[2] > eps) {
    cat("
      ")
    
    t1 = t0 + solve(inf.gamma(y, t0[1], t0[2])) %*% score.gamma(y, t0[1], t0[2]) #t1 estimate
    colnames(t1) = glue("step {i}")
    print(t1)
    
    r = cbind(r, t1) #save t1 estimate in r
    
    deltadev = as.vector(abs(t1 - t0)) #calculate eps
    print(glue("eps for lambda = {deltadev[1]}"))
    print(glue("eps for alpha = {deltadev[2]}"))
    
    
    t0 = t1
    i = i + 1
    cat("

      ")
    
    if(i>100){break}
    
  }
  #
  
  #final estimates
  hat.alpha = r[1, length(r[1, ])]
  hat.lambda = r[2, length(r[2, ])]
  
  
  
  results = rbind(c(a, amm, hat.alpha), c(l, lmm, hat.lambda))
  colnames(results) = c("parameter", "moments", "algorithm")
  rownames(results) = c("alpha", "lambda")
  
  
  
  distance = cbind(abs(results[, 2] - results[, 1]), abs(results[, 3] - results[, 1]))
  colnames(distance) = c(" MM - trueparameters ", "NR - trueparameters")
  
  r = list(results = results, distance = distance)
  return(r)

  #plot of estimates within steps
  if(plots==T){
    par(mfrow = c(2, 1))
    plot(
      r[1, ],
      ylab = "alpha",
      xlab = "step",
      ylim = c(a - 0.2, a + 0.2),
      type = "b",
      lty = 3,
      main = "Red lines are true parameters value"
    )
    abline(h = a, col = "red")
    plot(
      r[2, ],
      ylab = "lambda",
      xlab = "step",
      ylim = c(l - 0.2, l + 0.2),
      type = "b",
      lty = 3
    )
    abline(h = l, col = "red")
    par(mfrow = c(1, 1))
  }  

}










####QUASI-NEWTON####
#In this approach we find estimates only for alpha

#We can find an explicit formula to estimate lambda with ML method
#hat.lambda(ML)= (n*a)/sum(y)

#We start the algorithm with estimates calculated by method of moment's
#and then compute Newton-Raphson algorithm only for alpha;
#with the alpha found by NR we calculate new lambda(ML) and the continue


quasi.newton.gamma = function(y, eps = 0.00000001, plots=F) {
  #estimates with method of moments
  amm = (mean(y))^2/(var(y))
  lmm = mean(y)/var(y)
  t0 = c(amm, lmm)
  
  #matrix for saving results
  r = matrix(NA, 2, 1)
  rownames(r) = c("alpha", "lambda")
  colnames(r) = glue("step 0")
  r[, 1] = t0
  r
  
  
  #initialization
  
  deltadev = c(2, 2)
  i = 1
  
  #Algorithm
  while (abs(deltadev)[1] > eps & abs(deltadev)[2] > eps) {
    cat("
      ")
    
    a1 = t0[1] + (score.gamma(y, t0[1], t0[2])[1]) / (inf.gamma(y, t0[1], t0[2])[1, 1]) #t1 estimate
    
    l1 = (n*a1)/sum(y)
    
    t1 = as.matrix(c(a1, l1), 2, 1)
    colnames(t1) = glue("step {i}")
    print(t1)
    
    r = cbind(r, t1) #save t1 estimate in r
    
    deltadev = as.vector(abs(t1 - t0)) #calculate eps
    print(glue("eps for lambda = {deltadev[1]}"))
    print(glue("eps for alpha = {deltadev[2]}"))
    
    
    t0 = t1
    i = i + 1
    cat("

      ")
    if (i>100) {break}
    
  }
  #
  
  
  #final estimates
  hat.alpha = r[1, length(r[1, ])]
  hat.lambda = r[2, length(r[2, ])]
  
  
  
  results = rbind(c(a, amm, hat.alpha), c(l, lmm, hat.lambda))
  colnames(results) = c("parameter", "moments", "Ql.algorithm")
  rownames(results) = c("alpha", "lambda")
  
  
  
  distance = cbind(abs(results[, 2] - results[, 1]), abs(results[, 3] -
                                                           results[, 1]))
  colnames(distance) = c(" MM - trueparameters ", "QL.NR - trueparameters")
  
  r = list(results = results, distance = distance)
  return(r)
  
  #plot of estimates within steps
  if(plots==T){
    par(mfrow = c(2, 1))
    plot(
      r[1, ],
      ylab = "alpha",
      xlab = "step",
      ylim = c(a - 0.2, a + 0.2),
      type = "b",
      lty = 3,
      main = "Red lines are true parameters value"
    )
    abline(h = a, col = "red")
    plot(
      r[2, ],
      ylab = "lambda",
      xlab = "step",
      ylim = c(l - 0.2, l + 0.2),
      type = "b",
      lty = 3
    )
    abline(h = l, col = "red")
    par(mfrow = c(1, 1))
  }
  
}



####simulation
eps = 0.000001
B = 100                           #number of samples
n = 100                          #sample dimension
a = runif(1, 1, 15)
a           #true alpha
l = runif(1, 1, 15)
l           #true lambda


#matrices to save results
parameter = matrix(c(a, l), B, 2, byrow = T)
colnames(parameter) = c("alpha", "lambda")
rownames(parameter) = c(1:B)

moments = matrix(NA, B, 2)
colnames(moments) = c("alphaMM", "lambdaMM")
rownames(moments) = c(1:B)

newton = matrix(NA, B, 2)
colnames(newton) = c("alphanewton", "lambdanewton")
rownames(newton) = c(1:B)

qnewton = matrix(NA, B, 2)
colnames(qnewton) = c("alphaquasi", "lambdaquasi")
rownames(qnewton) = c(1:B)


for (i in 1:B) {
  y = rgamma(n, shape = a, rate = l)     #sample
  
  nn = newton.gamma(y, eps)
  qn = quasi.newton.gamma(y, eps)
  
  moments[i, 1] = nn$results[1, 2]
  moments[i, 2] = nn$results[2, 2]
  
  newton[i, 1] = nn$results[1, 3]
  newton[i, 2] = nn$results[2, 3]
  qnewton[i, 1] = qn$results[1, 3]
  qnewton[i, 2] = qn$results[2, 3]
}


#all results
results = cbind(parameter, moments, newton, qnewton)
colnames(results) = c(
  "alpha",
  "lambda",
  "alphaMM",
  "lambdaMM",
  "alphanewton",
  "lambdanewton",
  "alphaquasi",
  "lambdaquasi"
)

head(results)


#distorsion of estimators
b2.moments = (colMeans(moments) - parameter[1, ])^2
b2.newton =  (colMeans(newton) - parameter[1, ])^2
b2.qnewton = (colMeans(qnewton) - parameter[1, ])^2
B2=rbind(b2.moments, b2.newton, b2.qnewton)
colnames(B2)=c("alpha", "lambda")
rownames(B2)=c("moments", "newton", "qnewton")
B2


#variance of estimators
var.moments = apply(moments, 2, var)
var.newton = apply(newton, 2, var)
var.qnewton = apply(qnewton, 2, var)
lirc=c(1/(n*trigamma(a)), 1/((n*a)/l^2))
V=rbind(var.moments, var.newton, var.qnewton, lirc)
colnames(V)=c("alpha", "lambda")
rownames(V)=c("moments", "newton", "qnewton", "lirc")
V


#MSE of estimators
MSE.moments = b2.moments + var.moments
MSE.newton = b2.newton + var.newton
MSE.qnewton = b2.qnewton + var.qnewton

MSE=rbind(MSE.moments, MSE.newton, MSE.qnewton)
colnames(MSE)=c("alpha", "lambda")
rownames(MSE)=c("moments", "newton", "qnewton")
MSE


s=seq(-4, 4, l=1000)

{par(mfrow=c(3,2))
hist(scale(moments[,1]), freq=F, main="amm"); lines(s, dnorm(s), col="blue")
hist(scale(moments[,2]), freq=F, main="lmm"); lines(s, dnorm(s), col="blue")
hist(scale(newton[,1]), freq=F, main="anewton"); lines(s, dnorm(s), col="blue")
hist(scale(newton[,2]), freq=F, main="lnewton"); lines(s, dnorm(s), col="blue")
hist(scale(qnewton[,1]), freq=F, main="lqnewton"); lines(s, dnorm(s), col="blue")
hist(scale(qnewton[,2]), freq=F, main="aqnewton"); lines(s, dnorm(s), col="blue")
par(mfrow=c(1,1))}


{par(mfrow=c(3,2))
  plot(ecdf(scale(moments[,1])), main="amm"); lines(s, pnorm(s), col="blue")
  plot(ecdf(scale(moments[,2])), main="lmm"); lines(s, pnorm(s), col="blue")
  plot(ecdf(scale(newton[,1])), main="anewton"); lines(s, pnorm(s), col="blue")
  plot(ecdf(scale(newton[,2])), main="lnewton"); lines(s, pnorm(s), col="blue")
  plot(ecdf(scale(qnewton[,1])), main="lqnewton"); lines(s, pnorm(s), col="blue")
  plot(ecdf((scale(qnewton[,2]))), main="aqnewton"); lines(s, pnorm(s), col="blue")
  par(mfrow=c(1,1))}
