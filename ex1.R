set.seed(123)
#data
weib.y = c(225, 171, 198, 189, 189, 135, 162, 135, 117, 162)
#####functions#####
t1 = function(data){ #sufficient statistic
  t = sum( log(data) )
  return(t)
}
nllik = function(th, data){ #log likelihood
  n = length(data)
  t = t1(data)
  l = -n*lgamma(th)+n*th*log(3)+(th-1)*t
  return(-l)
}
score = function(th, data){ #score function
  n = length(data)
  t = t1(data)
  s = -n*digamma(th)+n*log(3) + t
  return(s)
}
#### compute the mle from the score ####
# uniroot(nscore, c(1e-10,1000), data = weib.y)
mle_score = nleqslv::nleqslv(0.5, score, data = weib.y, jacobian = TRUE)
mle_score$x 
mle_score$fvec
mle_score$jac

#### compute the mle from the log-likelihood ####
#optimize(nllik, interval = c(0,1000), data = weib.y)
#optim(400, fn = nllik, method = 'SANN',  hessian = T, data = weib.y)
mle_llik = optim(400, fn = nllik, method = 'L-BFGS-B', 
                 lower = 0, upper = 1000, hessian = T, data = weib.y)
mle_llik$value
mle_llik$par
mle_llik$hessian

#### plot the relative log-likelihood for theta
#base plot
plot(function(x) -nllik(x, data = weib.y) + mle_llik$value, 0,1000, 
     ylab = 'relative log-likelihood', xlab = expression(theta))

# #teamggplot
tibble(theta = seq(0,1000,1)) %>%
  mutate(rel_ll = -nllik(theta, data = weib.y) + mle_llik$value) %>%
  ggplot(aes(x = theta, y = rel_ll))+
  geom_line() +
  labs(y = 'relative log-likelihood', x = expression(theta))+
  theme_minimal()
