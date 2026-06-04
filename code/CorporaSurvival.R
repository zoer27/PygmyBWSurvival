#Estimating survival from Ovarian Corpora
#Corpora as a proxy for age--converted to age-since sexual maturity
#based on posterior of corpora accumulation rate estiamted in 
#Rand et al. (in review) Estimates of blue whale growth and reproduction parameters from historical data. Marine Mammal Science
#Code is based on brms catch-curve model (see EarplugSurvival.R) but was
#adapted in Stan so that it can explicitly incorporate the uncertainty from the 
#corpora accumulation rate
#
#Zoe Rand
#Last updated: June 4, 2026
#
library(tidyverse)
library(here)
library(posterior)
library(bayesplot)
library(brms)
library(rstan)

options(mc.cores = 4)


# Read in data ------------------------------------------------------------
Corp_py<-read_csv(here("data", "Branch2009_pygmy_corpora.csv"))

CorpFreq_py<-Corp_py %>% rowwise() %>% 
  mutate(Freq = sum(c_across(-corpora))) %>%
  rename(Corpora_count = corpora) %>%
  select(Corpora_count, Freq)


# Posterior for corpora rates ---------------------------------------------

rates<-read_csv(here("data", "rate_posterior.csv"))

#approximate distribution for prior
rsnorm<-rskew_normal(2000, 2.824, 0.482, 2)
dsnorm<-dskew_normal(rsnorm, 2.824, 0.482, 2)

supp_plot<-ggplot() + geom_histogram(aes(x = rates$P_rates, y = after_stat(density)), bins = 50, alpha = 0.5) +
  geom_density(aes(x = rsnorm, y = dsnorm), stat = "identity", color = "blue") +
  theme_classic() + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))+
  labs(x = "Corpora rate", y = "density")
supp_plot


# Fit Model -------------------------------------------------------------------

mu<-2.824
sd<-0.482
alpha<-2.00

#reparameterizing for stan
omega<-sd/(sqrt(1-((2/pi)*((alpha^2)/(1+alpha^2)))))
xi<-mu - (omega*(alpha/sqrt(1+alpha^2))*(sqrt(2/pi)))


#corpora_data
Corp_py_mod0<-CorpFreq_py%>% filter(Corpora_count > 0)

#corpora at highest frequency
Corp_py_mod0[which(Corp_py_mod0$Freq == max(Corp_py_mod0$Freq)),]


#starting the model at corpora count 2 so starting at 3
Corp_py_start<-3

Corp_mod<-Corp_py_mod0 %>% filter(Corpora_count >= Corp_py_start) %>% 
  add_row(Corpora_count = (max(Corp_py_mod0$Corpora_count) +1):(2*max(Corp_py_mod0$Corpora_count)), Freq = 0)

ggplot() + geom_point(aes(x = Corp_mod$Corpora_count, y = Corp_mod$Freq))

#using brms to get everything in the right format
mod<-brm(Freq ~ Corpora_count + (1|Corpora_count), data = Corp_mod, 
         family = poisson, 
         iter = 10, 
         warmup = 5)

#get data based on brms model
stan_dat<-standata(mod)

#add info for corpora rate prior
stan_dat$alpha <-alpha
stan_dat$omega<- omega
stan_dat$xi <-xi

#fit model with Stan
inits<-function(){
  list(Corpora_rate = 2.7, 
       b = array(0))
}

fit1<-stan(file = here("code", "CorporaSurvival.stan"), 
           data = stan_dat, 
           chains = 4, 
           warmup = 3000, 
           iter = 6000, 
           control = list(max_treedepth = 15, adapt_delta = 0.95), 
           init = inits, 
           seed = 134)


summary(fit1)


mcmc_trace(fit1, pars = c("b[1]", "b_Intercept", "sd_1[1]", "Corpora_rate"))
mcmc_dens(fit1, pars = c("b[1]", "b_Intercept", "sd_1[1]"))
mcmc_acf(fit1, pars = c("b[1]", "b_Intercept", "sd_1[1]"))
rhats<-rhat(fit1)
max(rhats, na.rm = TRUE)

#plotting
post<-as_draws(fit1)
post1<-merge_chains(post)

#making sure corpora_rate looks reasonable
cr_post<-subset_draws(post1, "Corpora_rate")
ggplot() + geom_density(aes(x = cr_post)) + 
  geom_density(aes(x = rsnorm, y = dsnorm), stat = "identity", color = "blue") + 
  labs(x = "Corpora_rate", y = "density")

#survival estimate
z_est<-subset_draws(post1, "b[1]")
head(z_est)
s_est<-exp(z_est)
summary(s_est)
ggplot() + geom_density(aes(x = s_est))


#save model fit for plotting
saveRDS(fit1, here("results", "corpora_survival_incorporated_posterior.RDS"))

