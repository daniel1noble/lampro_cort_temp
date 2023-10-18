library(car)
library(ggplot2)
library(lme4)
library(lmerTest)
library(glmm)
library(MuMIn)
library (emmeans)
library(dplyr)
library(installr)
updateR()


Liz <- read.csv ("C:/Users/U1115575/Documents/lampro_cort_temp/Projects/Body Size/Pilot.csv")
View(Liz)

Liz$trt <- as.factor(Liz$trt)
Liz$hormone <- as.factor(Liz$hormone)
Liz$temp <- as.factor(Liz$temp)
Liz$clutch <- as.factor(Liz$clutch)

###Body conditions as mass/svl
Liz$condition_juv <- Liz$remeasure_mass / Liz$remeasure_SVL
Liz$condition_hatch <- Liz$bd_mass_orig / Liz$bd_svl_orig_mm


#Scaled mass calculations
##Ln transform data
Liz$lnbd_mass_orig <- log(Liz$bd_mass_orig)
Liz$lnbd_svl_orig_mm <- log(Liz$bd_svl_orig_mm)

#OLS regression of ln transformed data
fit <- lm(lnbd_mass_orig ~ lnbd_svl_orig_mm, data=Liz)
summary (fit)
plot(lnbd_mass_orig ~ lnbd_svl_orig_mm, data=Liz)
abline(fit)

##Pearsons correlation coefficient
PCor <- cor(Liz$lnbd_mass_orig, Liz$lnbd_svl_orig_mm)
mod1 <- lm(lnbd_mass_orig~lnbd_svl_orig_mm, data=Liz)

##OLS coefficient from summary
summary(mod1)

##Calculates mean tarsus for scaled mass equation
MeanSVL_orig <- mean(Liz$bd_svl_orig_mm)

SMCof <- 1.4234/PCor

##calculates scaled mass; plots scaled mass to residual mass
Liz$ScaledMass_orig <- Liz$bd_mass_orig * (MeanSVL_orig/Liz$bd_svl_orig_mm) ^ (SMCof)
plot (ScaledMass_orig ~ condition_hatch, data=Liz)



###treatments on body condition
conhatch_mod <- lm (condition_hatch ~ temp + hormone  + egg_mass, data = Liz)
Anova(conhatch_mod)
summary(conhatch_mod)

conjuv_mod <- lm(condition_juv ~ temp + hormone + days_posthatch, data = Liz)
anova(conjuv_mod)
summary(conjuv_mod)

emm_conjuv_mod <- emmeans (conjuv_mod, pairwise ~ hormone)
emm_conjuv_mod

Boxplot(Liz$condition_juv, Liz$hormone)

###range of days post-hatch
range(Liz$days_posthatch, na.rm=TRUE)
mean(Liz$days_posthatch, na.rm=TRUE)



##SVL corrected for days post hatching
###regression of SVL against days posthatch is highly significant 
##(F = 16.47, p<0.001, df = 1,103)
cor_SVL <- lm (remeasure_SVL ~ days_posthatch, data = Liz, na.action = na.exclude)
anova(cor_SVL)
summary(cor_SVL)

##plot SVL against days post hatch
plot (remeasure_SVL ~ days_posthatch, data = Liz, 
      xlab = "Days post hatch", ylab = "SVL (cm)")
abline (lm (remeasure_SVL ~ days_posthatch, data = Liz))

##save SVL residuals posthatch
Liz$SVL_residuals <- residuals(cor_SVL)


###tail length corrected for days post hatching
###regression of tail against days posthatch is highly significant 
##(F = 10.89, p=0.001, df = 1,103)
cor_tail <- lm (remeasure_tail ~ days_posthatch, data = Liz, na.action = na.exclude)
Anova(cor_tail)
summary(cor_tail)

##plot tail against days post hatch
plot (remeasure_tail ~ days_posthatch, data = Liz, 
      xlab = "Days post hatch", ylab = "tail length (cm)")
abline (lm (remeasure_tail ~ days_posthatch, data = Liz))


##save tail residuals posthatch
Liz$tail_residuals <- residuals (cor_tail)

###mass corrected for days post hatching
###regression of tail against days posthatch is highly significant 
##(F = 21.37, p<0.001, df = 1,103)
cor_mass <- lm (remeasure_mass ~ days_posthatch, data = Liz, na.action = na.exclude)
anova(cor_mass)
summary(cor_mass)

##plot mass against days post hatch
plot (remeasure_mass ~ days_posthatch, data = Liz, 
      xlab = "Days post hatch", ylab = "mass (g)")
abline (lm (remeasure_mass ~ days_posthatch, data = Liz))

##save mass residuals post hatch
Liz$mass_residuals <- residuals (cor_mass)

##Box plot of body measurement residuals (post hatch) by treatment group
Boxplot(Liz$mass_residual, Liz$hormone:Liz$temp)
Boxplot(Liz$SVL_residuals, Liz$temp:Liz$hormone)
Boxplot(Liz$tail_residuals, Liz$temp:Liz$hormone)

##box plot by of body measurement residuals (post hatch) by treatment group
Boxplot(Liz$mass_residual, Liz$hormone)
  order_mass <- with(Liz, reorder(hormone, Liz$mass_residual, median, na.rm=TRUE))
  Boxplot(Liz$mass_residual ~ order_mass, xlab="CORT treatment", ylab="Residual mass")
  

  
Boxplot(Liz$SVL_residual, Liz$temp)
Boxplot(Liz$tail_residual, Liz$hormone)
Boxplot(Liz$mass_residuals, Liz$temp:Liz$hormone)

###effects of treatments on mass
###temp has effect (F=4.44, p = 0.04, df = 1, 99)
###hormone has effect (F=3.82, p = 0.03, df=2, 99)
##no interaction effect (F=0.37, p = 0.69, df = 2, 99)
mass_mod <- lm (mass_residuals ~ temp + hormone + temp:hormone, data = Liz)
anova(mass_mod)
summary(mass_mod)

###Using days posthatch instead of residuals
mass_mod2 <- lm (remeasure_mass ~ temp + hormone + days_posthatch + temp:hormone, data = Liz)
anova(mass_mod2)
summary(mass_mod2)


##pairwise comparisons of hormone treatments on mass
######control - high (p=0.02, df=99); control-low (p=0.60, df=99); high - low (p=0.16, df=99)
emm_mass_mod_hormone <- emmeans (mass_mod2, pairwise ~ hormone)
emm_mass_mod_hormone

##Works...but the axes are backwards?
plot(emm_mass_mod_hormone, horizontal = FALSE, 
     comparisons = TRUE, 
     CIs = TRUE, 
     xlab = "Estimated marginal means", 
     ylab = "CORT treatment")


##pairwise comparisons of temp treatments on mass
######23-28 (p=0.052)
emm_mass_mod_temp <- emmeans (mass_mod, pairwise ~ temp)
emm_mass_mod_temp


###
###for mass juvenile data
data_summary <- aggregate(mass_residuals ~ temp, Liz,       # Create summary data
                          function(mass_residuals) c(mean = mean(mass_residuals),
                                                   se = sd(mass_residuals) / sqrt(length(mass_residuals))))
data_summary <- data.frame(group = data_summary[ , 1], data_summary$mass_residuals)
data_summary  


ggplot(data_summary, aes(y = mean, x = factor (group))) +
  geom_point(aes(colour = factor (group)), size = 8, stat = 'identity') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2,
                position=position_dodge(.9)) +
  
  labs(x="Incubation temperature (C)", y = "Residual mass")+
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.text.y = element_text(size = 30))+
  theme(axis.text.x = element_text(size = 30))+
  theme(axis.title = element_text(size=35))


###
###for mass juvenile data by CORT 


data_summary <- aggregate(mass_residuals ~ hormone, Liz,       # Create summary data
                          function(mass_residuals) c(mean = mean(mass_residuals),
                                                     se = sd(mass_residuals) / sqrt(length(mass_residuals))))
data_summary <- data.frame(group = data_summary[ , 1], data_summary$mass_residuals)
data_summary  





ggplot(data_summary, aes(y = mean, x = factor (group))) +
  geom_point(aes(size = 8, stat = 'identity')) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2,
                position=position_dodge(.9)) +
  
  labs(x="CORT treatment", y = "Residual mass")+
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.text.y = element_text(size = 30))+
  theme(axis.text.x = element_text(size = 30))+
  theme(axis.title = element_text(size=35))



###for juvenile SVL

data_summary <- aggregate(SVL_residuals ~ hormone, Liz,       # Create summary data
                          function(SVL_residuals) c(mean = mean(SVL_residuals),
                                                     se = sd(SVL_residuals) / sqrt(length(SVL_residuals))))
data_summary <- data.frame(group = data_summary[ , 1], data_summary$SVL_residuals)
data_summary  





ggplot(data_summary, aes(y = mean, x = factor (group))) +
  geom_point(aes(size = 8, stat = 'identity')) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2,
                position=position_dodge(.9)) +
  
  labs(x="CORT treatment", y = "Residual SVL")+
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.text.y = element_text(size = 30))+
  theme(axis.text.x = element_text(size = 30))+
  theme(axis.title = element_text(size=35))








##Works...but the axes are backwards?
plot(emm_mass_mod_temp, horizontal = FALSE, 
     comparisons = TRUE, 
     CIs = TRUE, 
     xlab = "Estimated marginal means", 
     ylab = "Temp treatment")

##take home message: high CORT lizards weigh less than conrol lizards


###effects of treatments on SVL
###temp has NO effect (F=2.74, p = 0.10, df = 1, 99)
###hormone has NO effect (F=2.00, p = 0.14, df=2, 99)
##no interaction effect (F=0.84, p = 0.44, df = 2, 99)
SVL_mod <- lm (SVL_residuals ~ temp + hormone + temp:hormone, data = Liz)
anova(SVL_mod)
summary(SVL_mod)


###Using days posthatch instead of residuals
svl_mod2 <- lm (remeasure_SVL ~ temp + hormone + days_posthatch + temp:hormone, data = Liz)
anova(svl_mod2)
summary(svl_mod2)

emm_svl_mod2 <- emmeans (svl_mod2, pairwise ~ temp)
emm_svl_mod2


###take home message for SVL: developmental treatment has no effect on SVL


###effects of treatments on tail
###temp has NO effect (F=2.29, p = 0.13, df = 1, 99)
###hormone has effect (F=3.91, p = 0.02, df=2, 99)
##no interaction effect (F=0.47, p = 0.62, df = 2, 99)
tail_mod <- lm (tail_residuals ~ temp + hormone + temp:hormone, data = Liz)
anova(tail_mod)
summary(tail_mod)

tail_mod2 <- lm (remeasure_tail ~ temp + hormone + days_posthatch, data = Liz)
anova(tail_mod2)
summary(tail_mod2)


##pairwise comparisons of hormone treatments on tail
######control - high (p=0.11, df=99); control-low (p=0.78, df=99); high - low (p=0.02, df=99)
##no need for pairwise comparisons of temp (only two)
##no need for pairwise comparisons of temp:hormone as the interaction term is NS
emm_tail_mod_2 <- emmeans (tail_mod, pairwise ~ temp)
emm_tail_mod_hormone

##Works...but the axes are backwards?
plot(emm_tail_mod_hormone, horizontal = FALSE, 
     comparisons = TRUE, 
     CIs = TRUE, 
     xlab = "Estimated marginal means", 
     ylab = "CORT treatment")


### take home message is that lizards that received low CORT had shorter tails than high CORT lizards
### no difference between high CORT and control which is midly concerning


######Now for hatch data

##box plot for body difference by hormone and temp
Boxplot(Liz$bd_tail_orig_mm, Liz$temp)
Boxplot(Liz$bd_svl_orig_mm, Liz$temp)
Boxplot(Liz$bd_mass_orig, Liz$temp)


###for mass hatchling data
data_summary <- aggregate(bd_mass_orig ~ hormone, Liz,       # Create summary data
                          function(bd_mass_orig) c(mean = mean(bd_mass_orig),
                                        se = sd(bd_mass_orig) / sqrt(length(bd_mass_orig))))
data_summary <- data.frame(group = data_summary[ , 1], data_summary$bd_mass_orig)
data_summary  


ggplot(data_summary, aes(y = mean, x = factor (group))) +
  geom_point(aes(size = 8, stat = 'identity')) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2,
                position=position_dodge(.9)) +

  labs(x="CORT treatment", y = "Mass (g)")+
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.text.y = element_text(size = 30))+
  theme(axis.text.x = element_text(size = 30))+
  theme(axis.title = element_text(size=35))

###effects of treatments on mass at hatching
###temp has effect (F=4.68, p = 0.03, df = 1, 120)
###hormone has effect (F=3.42, p = 0.04, df=2, 120)
##no interaction effect (F=0.09, p = 0.91, df = 2, 120)
hatch_mass_mod <- lm (bd_mass_orig ~ temp + hormone + egg_mass + temp:hormone, data = Liz)
anova(hatch_mass_mod)
summary(hatch_mass_mod)

summarise(Liz, n(), groups = "23")


##pairwise comparisons of hormone treatments on hatch mass
######control - high (p=0.04); control-low (p=0.15); high - low (p=0.76)
emm_hatch_mass_mod_temp <- emmeans (hatch_mass_mod, pairwise ~ temp)
emm_hatch_mass_mod_temp

##pairwise comparisons of hormone treatments on hatch mass
######control - high (p=0.04); control-low (p=0.15); high - low (p=0.76)
emm_hatch_mass_mod_hormone <- emmeans (hatch_mass_mod, pairwise ~ hormone)
emm_hatch_mass_mod_hormone

##Works...but the axes are backwards?
plot(emm_hatch_mass_mod_hormone, horizontal = FALSE, 
     comparisons = TRUE, 
     CIs = TRUE, 
     xlab = "Estimated marginal means of mass", 
     ylab = "CORT treatment")

##pairwise comparisons of temp treatments on hatch mass
######23 - 28 (p=0.04)
emm_hatch_mass_mod_temp <- emmeans (hatch_mass_mod, pairwise ~ temp)
emm_hatch_mass_mod_temp

##Works...but the axes are backwards?
plot(emm_hatch_mass_mod_temp, horizontal = FALSE, 
     comparisons = TRUE, 
     CIs = TRUE, 
     xlab = "Estimated marginal means of mass", 
     ylab = "Temp treatment")


###effects of treatments on svl at hatching
###temp has NO effect (F=0.14, p = 0.71, df = 1, 120)
###hormone has effect (F=5.69, p = 0.004, df=2, 120)
##no interaction effect (F=0.18, p = 0.18, df = 2, 120)
hatch_svl_mod <- lm (bd_svl_orig_mm ~ temp + hormone + egg_mass + temp:hormone, data = Liz)
anova(hatch_svl_mod)
summary(hatch_svl_mod)

##pairwise comparisons of hormone treatments on hatch svl
######control - high (p=0.006); control-low (p=0.90); high - low (p=0.02)
emm_hatch_svl_mod_hormone <- emmeans (hatch_svl_mod, pairwise ~ hormone)
emm_hatch_svl_mod_hormone


###for SVL hatchling data
data_summary <- aggregate(bd_svl_orig_mm ~ hormone, Liz,       # Create summary data
                          function(bd_svl_orig_mm) c(mean = mean(bd_svl_orig_mm),
                                                   se = sd(bd_svl_orig_mm) / sqrt(length(bd_svl_orig_mm))))
data_summary <- data.frame(group = data_summary[ , 1], data_summary$bd_svl_orig_mm)
data_summary  


ggplot(data_summary, aes(y = mean, x = factor (group))) +
  geom_point(aes(size = 8, stat = 'identity')) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2,
                position=position_dodge(.9)) +
  
  labs(x="CORT treatment", y = "SVL (mm)")+
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.text.y = element_text(size = 30))+
  theme(axis.text.x = element_text(size = 30))+
  theme(axis.title = element_text(size=35))








##Works...but the axes are backwards?
plot(emm_hatch_svl_mod_hormone, horizontal = FALSE, 
     comparisons = TRUE, 
     CIs = TRUE, 
     xlab = "Estimated marginal means", 
     ylab = "CORT treatment")

###effects of treatments on tail at hatching
###temp has NO effect (F=0.81, p = 0.37, df = 1, 120)
###hormone has NO effect (F=1.39, p = 0.25, df=2, 120)
##no interaction effect (F=0.82, p = 0.44, df = 2, 120)
hatch_tail_mod <- lm (bd_tail_orig_mm ~ temp + hormone + egg_mass + temp:hormone, data = Liz)
anova(hatch_tail_mod)
summary(hatch_tail_mod)


###Effect of treatments on days to hatch
###temp has effect (F=179.38, p <0.001, df = 1, 112)
###hormone has NO effect (F=0.15, p = 0.86, df=2, 112)
##no interaction effect (F=0.35, p = 0.71, df = 2, 112)
mod_days_to_hatch <- lm (days_to_hatch ~ temp + hormone + temp:hormone, data = Liz)
anova(mod_days_to_hatch)
summary(mod_days_to_hatch)





