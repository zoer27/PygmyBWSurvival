#Estimating annual survival for pygmy blue whales 
#Using catch curve (Poisson GLMM)
#Model built in brms package
#Zoe Rand
#Last updated: June 4, 2026

library(tidyverse)
library(here)
library(brms)
library(bayesplot)
library(posterior)
options(mc.cores = parallel::detectCores())
options(brms.backend = "cmdstanr")


# Read in data ------------------------------------------------------------
#japanese data has catch year and information about aging accuracy
PyDat_J<-read_csv(here("data", "Earplug_pygmy_J.csv"))
summary(PyDat_J)
#only using "accurate ages" 
PyDat_J<-PyDat_J %>% filter(AccurateAge == "Yes")

#soviet data only has ages, lengths, and sex
#ages are based on 1 layer per year
#did not have whale stretching (and most whales are under 21m)
PyDat_S<-read_csv(here("data", "Earplug_pygmy_S.csv"))
summary(PyDat_S)
#combine into age frequency
AgeFreqS<- PyDat_S %>% group_by(Sex, Age) %>% summarise(n = n()) %>% mutate(Sex = ifelse(Sex == "Female", "F", "M"))
AgeFreqJ <- PyDat_J %>% group_by(Sex, Age) %>% summarise(n = n())

AgeFreq<-bind_rows(AgeFreqS, AgeFreqJ) %>% group_by(Sex, Age) %>% summarise(Freq = sum(n))
summary(AgeFreq)
ggplot(data = AgeFreq) + geom_col(aes(x = Age, y = Freq, fill = Sex), position = "dodge")

write_csv(AgeFreq, here("data", "Pygmy_AgeFreq_allearplugs.csv"))

# Model with all data sources and both sexes -------------------------------------------

#age frequency not seperated by sex
AgeFreq_mod1<-AgeFreq %>% ungroup() %>% group_by(Age) %>% summarise(Freq = sum(Freq))

#plot
ggplot(data = AgeFreq_mod1) + geom_col(aes(x = Age, y = Freq)) + theme_classic() + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))

#find age at max frequency
AgeFreq_mod1[which(AgeFreq_mod1$Freq == max(AgeFreq_mod1$Freq)),]
#age with maximum frequency is age 16, so starting model at age 17
Age_maxf<-17

#extending dataframe based on Millar 2015
AgeFreq_mod<-AgeFreq_mod1 %>% filter(Age >= Age_maxf) %>% add_row(Age = (max(AgeFreq_mod1$Age) +1):(2*max(AgeFreq_mod1$Age)), Freq = 0)
summary(AgeFreq_mod)

#plot
ggplot(AgeFreq_mod) + geom_point(aes(x = Age, y = Freq))

#run model
mod1<-brm(Freq ~ Age + (1|Age), data = AgeFreq_mod, family = poisson)
summary(mod1)

#save model results
saveRDS(mod1, here("results", "sex_aggregated_alldata.RDS"))

#diagnostics
mcmc_trace(mod1, pars = c("b_Intercept", "b_Age", "sd_Age__Intercept"))
mcmc_acf(mod1, pars = c("b_Intercept", "b_Age", "sd_Age__Intercept"))
mcmc_dens(mod1, pars = c("b_Intercept", "b_Age", "sd_Age__Intercept"))


#survival estimate
post<-as_draws(mod1)
post1<-merge_chains(post)
z_est<-subset_draws(post1, "b_Age")[[1]]
head(z_est$b_Age)
s_est<-exp(z_est$b_Age)
head(s_est)
s_med<-quantile(s_est, 0.5)

#plot
ggplot() + geom_density(aes(x = s_est), fill = "lightblue") + 
  geom_vline(aes(xintercept = s_med), color = "blue") +
  theme_classic() + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  labs(x = "Survival estimate") + theme(axis.text.y = element_blank(), 
                                        axis.ticks.y = element_blank())
summary(s_est)

#credible intervals
quantile(s_est, probs = c(0.025, 0.5, 0.975))

#plot posterior predictive with data
post_pred<-posterior_predict(mod1)
quant_post_pred<-apply(post_pred, 2, quantile, probs = c(0.025, 0.5, 0.975))
pp_tab<-as_tibble(t(quant_post_pred))
pp_tab$Age<-AgeFreq_mod$Age
p1<-ggplot()  + 
  geom_errorbar(data = pp_tab, aes(x = Age, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.5, color = "blue")+
  geom_point(data = pp_tab, aes(x = Age, y = `50%`), color = "blue", shape = 4) + 
  geom_point(data = AgeFreq_mod1, aes(x = Age, y = Freq)) + 
  theme_classic() + 
  coord_cartesian(xlim = c(0, 83)) + 
  labs(x = "Age", y = "Frequency")
p1

#plot posterior log scale
pred<-posterior_linpred(mod1, newdata = tibble(Age = seq(17,83)), re_formula = NA) #plot linear predictor without random effects
quant_pred<-apply(pred, 2, quantile, probs = c(0.025, 0.5, 0.975))
p_tab<-as_tibble(t(quant_pred))
p_tab$Age<-seq(17,83)

p2<-ggplot(p_tab) + 
  geom_point(data = AgeFreq_mod1, aes(x = Age, y = log(Freq))) +
  geom_ribbon(aes(x = Age, ymin = `2.5%`, ymax = `97.5%`), alpha = 0.5) +
  geom_line(aes(x = Age, y = `50%`), color = "blue") +
  coord_cartesian(xlim = c(0, 60)) +
  labs(y = "Log frequency") +
  theme_classic()

p2

#plotting together 
library(patchwork)
p2 + p1 + plot_annotation(tag_levels = "a",  tag_prefix = c("("), 
                          tag_suffix = c(")"))


# Model with just Soviet data (both sexes) ----------------------------------------

#model with only soviet data
AgeFreq_mod_S<- AgeFreqS %>% group_by(Age) %>% summarise(Freq = sum(n))
ggplot() + geom_col(data = AgeFreq_mod_S, aes(x = Age, y = Freq))

#find age at max frequency
which(AgeFreq_mod_S$Freq == max(AgeFreq_mod_S$Freq)) #this is 4 but doesn't really seem realisitc

which(AgeFreq_mod_S$Freq == sort(AgeFreq_mod_S$Freq, TRUE)[2]) #second highest is 2

which(AgeFreq_mod_S$Freq == sort(AgeFreq_mod_S$Freq, TRUE)[3]) #third highest is 14 so using age 15

#starting model at age 15

S_dat_mod<-AgeFreq_mod_S %>% filter(Age >= 15) %>% add_row(Age = (max(AgeFreq_mod_S$Age) +1):(2*max(AgeFreq_mod_S$Age)), Freq = 0)
summary(S_dat_mod)

ggplot(S_dat_mod) + geom_point(aes(x = Age, y = Freq))

#running model
mod2<-brm(Freq ~ Age + (1|Age), data = S_dat_mod, family = poisson)
summary(mod2)

saveRDS(mod2, here("results", "sex_aggregated_Sovietdat.RDS"))

#diagnostics
mcmc_trace(mod2, pars = c("b_Intercept", "b_Age", "sd_Age__Intercept"))
mcmc_acf(mod2, pars = c("b_Intercept", "b_Age", "sd_Age__Intercept"))
mcmc_dens(mod2, pars = c("b_Intercept", "b_Age", "sd_Age__Intercept"))

#survival estimate
post2<-as_draws(mod2)
post2<-merge_chains(post2)
z_est2<-subset_draws(post2, "b_Age")[[1]]
head(z_est2$b_Age)
s_est2<-exp(z_est2$b_Age)
head(s_est2)
ggplot() + geom_density(aes(x = s_est2), fill = "lightgreen") + 
  theme_classic() + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  labs(x = "Survival estimate") + theme(axis.text.y = element_blank(), 
                                        axis.ticks.y = element_blank())
summary(s_est2)
quantile(s_est2, probs = c(0.025, 0.5, 0.975))


# Model with just Japanese data (both sexes) --------------------------------------

AgeFreq_mod_J<- AgeFreqJ %>% group_by(Age) %>% summarise(Freq = sum(n))
ggplot() + geom_col(data = AgeFreq_mod_J, aes(x = Age, y = Freq))

#age at maximum frequency
which(AgeFreq_mod_J$Freq == max(AgeFreq_mod_J$Freq))
AgeFreq_mod_J[c(19,20,21),]
#starting model at age 22

J_mod_dat<-AgeFreq_mod_J %>% filter(Age >= 22) %>% add_row(Age = (max(AgeFreq_mod_J$Age) +1):(2*max(AgeFreq_mod_J$Age)), Freq = 0)
summary(J_mod_dat)

ggplot(J_mod_dat) + geom_point(aes(x = Age, y = Freq))

#running model
mod3<-brm(Freq ~ Age + (1|Age), data = J_mod_dat, family = poisson)
summary(mod3)

#save model results
saveRDS(mod3, here("results", "sex_aggregated_Japanesedat.RDS"))


#diagnostics
mcmc_trace(mod3, pars = c("b_Intercept", "b_Age", "sd_Age__Intercept"))
mcmc_acf(mod3, pars = c("b_Intercept", "b_Age", "sd_Age__Intercept"))
mcmc_dens(mod3, pars = c("b_Intercept", "b_Age", "sd_Age__Intercept"))

#survival estimate
post3<-as_draws(mod3)
post3<-merge_chains(post3)
z_est3<-subset_draws(post3, "b_Age")[[1]]
head(z_est3$b_Age)
s_est3<-exp(z_est3$b_Age)
head(s_est3)
ggplot() + geom_density(aes(x = s_est3), fill = "lavender") + 
  theme_classic() + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  labs(x = "Survival estimate") + theme(axis.text.y = element_blank(), 
                                        axis.ticks.y = element_blank())
summary(s_est3)

quantile(s_est3, probs = c(0.025, 0.5, 0.975))


# Compare data sources ----------------------------------------------------

ggplot() + geom_density(aes(x = s_est), fill = "lightblue", alpha = 0.5) + 
  geom_density(aes(x = s_est2), fill = "lightpink", alpha = 0.5) + 
  geom_density(aes(x = s_est3), fill = "lightgreen", alpha = 0.5) + 
  theme_classic() + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  labs(x = "Survival estimate") + theme(axis.text.y = element_blank(), 
                                        axis.ticks.y = element_blank())



# Sensitivity to age at model start ---------------------------------------

#function to run model given starting age
run_model<-function(start_age){
  AgeFreq_mod<-AgeFreq_mod1 %>% filter(Age >= start_age) %>% add_row(Age = (max(AgeFreq_mod1$Age) +1):(2*max(AgeFreq_mod1$Age)), Freq = 0)
  mod<-brm(Freq ~ Age + (1|Age), data = AgeFreq_mod, family = poisson)
  post<-as_draws(mod)
  post<-merge_chains(post)
  z_est<-subset_draws(post, "b_Age")[[1]]
  s_est<-exp(z_est$b_Age)
  return(s_est)
}

#plot data to get a sense of ages
ggplot(data = AgeFreq_mod1) + geom_point(aes(x = Age, y = Freq)) + theme_classic() + 
  scale_x_continuous(breaks = seq(1,80, by = 2), expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0.1, 0)) +
  geom_vline(aes(xintercept = 14))

#ages to test
age_test<-seq(14, 20, by = 1)

#survival estimates
surv_sens<-lapply(age_test, run_model)

#get median
med_surv<-sapply(surv_sens, function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
colnames(med_surv)<-as.character(age_test)

#save for plotting
saveRDS(med_surv, here("results", "survivalsensitivity.RDS"))

meds<-med_surv[2,]

#difference in survival 
(meds[1] - meds[length(meds)])/mean(c(meds[1], meds[length(meds)])) * 100


# Model for just females -----------------------------------------------------

#data for females
AgeFreqF<-AgeFreq %>% filter(Sex == "F")
AgeFreq_modF<-AgeFreqF %>% ungroup() %>% group_by(Age) %>% summarise(Freq = sum(Freq))
ggplot(data = AgeFreq_modF) + geom_col(aes(x = Age, y = Freq)) + theme_classic() + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))

#age at maximum frequency
AgeFreq_modF[which(AgeFreq_modF$Freq == max(AgeFreq_modF$Freq)),]

#both 12 and 14 have highest frequency so using 15
Age_maxf<-15

#extending dataframe based on Millar 2015
AgeFreq_mod_F<-AgeFreq_modF %>% filter(Age >= Age_maxf) %>% add_row(Age = (max(AgeFreq_modF$Age) +1):(2*max(AgeFreq_modF$Age)), Freq = 0)
summary(AgeFreq_mod_F)

ggplot(AgeFreq_mod_F) + geom_point(aes(x = Age, y = Freq))


#run model
mod4<-brm(Freq ~ Age + (1|Age), data = AgeFreq_mod_F, family = poisson)
summary(mod4)

#save results
saveRDS(mod4, here("results", "female_only_alldat.RDS"))

#diagnostics
mcmc_trace(mod4, pars = c("b_Intercept", "b_Age", "sd_Age__Intercept"))
mcmc_acf(mod4, pars = c("b_Intercept", "b_Age", "sd_Age__Intercept"))
mcmc_dens(mod4, pars = c("b_Intercept", "b_Age", "sd_Age__Intercept"))


#survival estimate
post<-as_draws(mod4)
post1<-merge_chains(post)
z_est<-subset_draws(post1, "b_Age")[[1]]
head(z_est$b_Age)
s_est<-exp(z_est$b_Age)
head(s_est)
s_med<-quantile(s_est, 0.5)
summary(s_est)


# Model for just males ----------------------------------------------------
#for males 
AgeFreqM<-AgeFreq %>% filter(Sex == "M")
AgeFreq_modM<-AgeFreqM %>% ungroup() %>% group_by(Age) %>% summarise(Freq = sum(Freq))
ggplot(data = AgeFreq_modM) + geom_col(aes(x = Age, y = Freq)) + theme_classic() + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))

#age at maximum frequency
AgeFreq_modM[which(AgeFreq_modM$Freq == max(AgeFreq_modM$Freq)),]

#16 so starting model at age 17
Age_maxf_M<-17

#extending dataframe based on Millar 2015
AgeFreq_mod_M<-AgeFreq_modM %>% filter(Age >= Age_maxf_M) %>% add_row(Age = (max(AgeFreq_modM$Age) +1):(2*max(AgeFreq_modM$Age)), Freq = 0)
summary(AgeFreq_mod_M)

ggplot(AgeFreq_mod_M) + geom_point(aes(x = Age, y = Freq))


#run model
mod5<-brm(Freq ~ Age + (1|Age), data = AgeFreq_mod_M, family = poisson)
summary(mod5)

#save results
saveRDS(mod5, here("results", "male_only_alldat.RDS"))

#diagnostics
mcmc_trace(mod5, pars = c("b_Intercept", "b_Age", "sd_Age__Intercept"))
mcmc_acf(mod5, pars = c("b_Intercept", "b_Age", "sd_Age__Intercept"))
mcmc_dens(mod5, pars = c("b_Intercept", "b_Age", "sd_Age__Intercept"))


#survival estimate
post<-as_draws(mod5)
post1<-merge_chains(post)
z_est<-subset_draws(post1, "b_Age")[[1]]
head(z_est$b_Age)
s_est<-exp(z_est$b_Age)
head(s_est)
s_med<-quantile(s_est, 0.5)
summary(s_est)

