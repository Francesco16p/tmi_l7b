library(tidyverse)
library(ggpubr)

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

thhat = mle_llik$par #mle
llik_thhat = -mle_llik$value #maximum loglikelihood
j_thhat = mle_llik$hessian #observed info in the maximum

#### plot the relative log-likelihood for theta ####
#base plot
plot(function(x) -nllik(x, data = weib.y) - llik_thhat, 0,1000, 
     ylab = 'relative log-likelihood', xlab = expression(theta))
# #teamggplot
tibble(theta = seq(0,1000,1)) %>%
  mutate(rel_ll = -nllik(theta, data = weib.y) - llik_thhat) %>%
  ggplot(aes(x = theta, y = rel_ll))+
  geom_line(alpha =.8)

#### wald confidence interval theta ####
conf = .9
th_wald_ci = c(
  thhat - qnorm(1 - (1 - conf)/2) *1/sqrt(j_thhat),
  thhat + qnorm(1 - (1 - conf)/2) *1/sqrt(j_thhat)
) %>% print()

#### deviance confidence interval theta ####
th_dev_ci = c(
  #left solution
  uniroot(function(x) -nllik(x, data = weib.y) - llik_thhat + qchisq(conf,1)/2, 
          c(480, thhat))$root,
  #right solution
  uniroot(function(x) -nllik(x, data = weib.y) - llik_thhat + qchisq(conf,1)/2, 
          c(thhat, 520))$root
)

#### results for theta parameterization ####
th_tab_ci = data.frame(
  ci = c('wald', 'deviance'),
  lower = c(th_wald_ci[1], th_dev_ci[1]),
  upper = c(th_wald_ci[2], th_dev_ci[2])
) %>% print()

p = tibble(theta = seq(0,1000,1)) %>%
  mutate(rel_ll = -nllik(theta, data = weib.y) - llik_thhat) %>%
  ggplot(aes(x = theta, y = rel_ll))+
  geom_line(alpha =.8) +
  geom_hline(yintercept = -qchisq(conf,1)/2, linetype = 'dashed', alpha =.5) +
  scale_y_continuous(limits = c(-10,0),
                     breaks = round(c(-10,-7.5,-5,-2.5, -qchisq(conf,1)/2, 0),1),
                     labels = c(-10,-7.5,-5,-2.5, expression(frac( {chi^2}[{paste(1,',',1-alpha)}],2)) , 0) ) +
  geom_vline(xintercept = thhat, linetype = 'dashed', alpha = .5) +
  scale_x_continuous(limits = c(460,530), 
                     breaks = round(c(460, 480, thhat ,520,540),0),
                     labels = c('460', '480', expression(hat(theta)),'520','540')) +
  geom_vline(aes(col='wald', xintercept = th_wald_ci[[1]]), size =1, linetype = 'dashed', alpha = 1) +
  geom_vline(aes(col='wald', xintercept = th_wald_ci[[2]]), size =1, linetype = 'dashed', alpha = 1) +
  geom_vline(aes(col='deviance', xintercept = th_dev_ci[[1]]), linetype = 'dashed', alpha = 1) +
  geom_vline(aes(col='deviance', xintercept = th_dev_ci[[2]]), linetype = 'dashed', alpha = 1) +
  labs(y = 'relative log-likelihood', x = '', col = 'Confidence Interval:',
       caption = substitute(paste(hat(theta), '=', nn), list(nn=round(thhat,1))))+
  theme_minimal() +
  theme(legend.position = 'top')
ggarrange(ggtexttable(th_tab_ci, rows = NULL),p, nrow = 2, heights = c(.5, 1.5))
#### w parameterization ####
w_of = function(theta){log(theta)}
th_of = function(w){exp(w)}
nllik_w = function(w, data){
  nllik(th = th_of(w), data = data)
}

#equivariance properties
what = w_of(thhat) #equivariance
llik_what = -nllik_w(what, data = weib.y)#equivariance
#compute numerically the observed information
j_what = optimHess(what, nllik_w, data = weib.y)
optimHess(th_of(what), nllik, data = weib.y) #equal to j_thhat

#### relative log-likelihood w ####
tibble(w = seq(5,7,0.001)) %>%
  mutate( rel_ll = -nllik(th=exp(w), data = weib.y) - llik_what) %>%
  ggplot(aes(x = w, y = rel_ll))+
  geom_line(alpha =.8)+
  geom_vline(xintercept = what) +
  geom_hline(yintercept = -qchisq(conf,1)/2, linetype = 'dashed', alpha =.5) +
  scale_y_continuous(limits = c(-5,0))+
  scale_x_continuous(limits = c(6.15,6.25))

  
#### wald confidence interval omega ####
conf = .9
w_wald_ci = c(
  what - qnorm(1 - (1 - conf)/2) *1/sqrt(j_what),
  what + qnorm(1 - (1 - conf)/2) *1/sqrt(j_what)
) %>% print()

#### deviance confidence interval omega ####
w_dev_ci = c(
  uniroot(function(x) -nllik_w(x, data = weib.y) - llik_what + qchisq(conf,1)/2, 
          c(6.175, what))$root,
  uniroot(function(x) -nllik_w(x, data = weib.y) - llik_what + qchisq(conf,1)/2, 
          c(what, 6.25))$root
  
) %>% print()
#### results for omega parameterization ####
w_tab_ci = data.frame(
  ci = c('wald', 'deviance'),
  lower = c(w_wald_ci[1], w_dev_ci[1]),
  upper = c(w_wald_ci[2], w_dev_ci[2])
) %>% print()

p = tibble(w = seq(5,7,0.001)) %>%
  mutate( rel_ll = -nllik(th=exp(w), data = weib.y) - llik_what) %>%
  ggplot(aes(x = w, y = rel_ll))+
  geom_line(alpha =.8)+
  geom_vline(xintercept = what, linetype = 'dashed', alpha =.5) +
  geom_hline(yintercept = -qchisq(conf,1)/2, linetype = 'dashed', alpha =.5) +
  scale_y_continuous(limits = c(-5,0),
                     breaks = round(c(-5, -3, -qchisq(conf,1)/2, 0),1),
                     labels = c(-5,-3,expression(frac( {chi^2}[{paste(1,',',1-alpha)}],2)) , 0) ) +
  scale_x_continuous(limits = c(6.15,6.25), 
                     breaks = round(c(6.15, what ,6.25),2),
                     labels = c('6.15', expression(hat(omega)),'6.25')) +
  geom_vline(aes(col='wald', xintercept = w_wald_ci[[1]]), size =1, linetype = 'dashed', alpha = 1) +
  geom_vline(aes(col='wald', xintercept = w_wald_ci[[2]]), size =1, linetype = 'dashed', alpha = 1) +
  geom_vline(aes(col='deviance', xintercept = w_dev_ci[[1]]), linetype = 'dashed', alpha = 1) +
  geom_vline(aes(col='deviance', xintercept = w_dev_ci[[2]]), linetype = 'dashed', alpha = 1) +
  labs(y = 'relative log-likelihood', x = '', col = 'Confidence Interval:', 
       caption = substitute(paste(hat(omega), '=', nn), list(nn=round(what,3))))+
  theme_minimal() +
  theme(legend.position = 'top') 
p
ggarrange(ggtexttable(w_tab_ci, rows = NULL),p, nrow = 2, heights = c(.5, 1.5))
