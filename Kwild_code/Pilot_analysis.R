# WD, packages, data
pacman::p_load(dplyr, tidyverse, ggpubr, lme4, emmeans, car, lmerTest, MuMIn, glmm, installr, lubridate, performance)
setwd("~/Dropbox/lampro_cort_temp/")
liz_ondi <- read.csv ("Kwild_code/Pilot_cort_mortality_db.csv") %>%  # ondi data  database data
  mutate_at(c('trt', 'hormone', 'temp', 'clutch'), as.factor) %>% 
  rename(days_posthatch_juv1 = days_posthatch) %>% 
  mutate(remeasure_SVL = remeasure_SVL*10,
         remeasure_tail = remeasure_tail*10,
         measure_SVL_juv2 = measure_SVL_juv2*10,
         measure_tail_juv2 = measure_tail_juv2*10) # remeasurements for SVL & tail were in CM so transform

#determining days after hatch mortlity & days for remeasurements 2 (juv 2)
liz_ondi$remeasure_date <- dmy(liz_ondi$remeasure_date)
liz_ondi$hatch_date_ymd <- ymd(liz_ondi$hatch_date_ymd)
liz_ondi$mortality_date_dmy <- dmy(liz_ondi$mortality_date_dmy)
liz_ondi$measure_date_juv2 <- dmy(liz_ondi$measure_date_juv2)
liz_ondi$mortality_post_days <- difftime(liz_ondi$mortality_date_dmy,liz_ondi$hatch_date_ymd, units = "days")
liz_ondi$days_posthatch_juv2 <- difftime(liz_ondi$measure_date_juv2,liz_ondi$hatch_date_ymd, units = "days")
# naming final df and ordering for figures
Liz = mutate(liz_ondi, hormone = factor(hormone, 
                             levels = c("control",
                                        "low",
                                        "high"))) %>% 
               group_by(hormone)
################
#### 1.	What are the effects of developmental treatments (temp, cort, interaction) on time to hatch?
#### effect of temperature on incubation days - faster hatch warmer temps 
################
# mod ; interaction removed
hatch_dev_mod <- lm(days_to_hatch~ temp + hormone, data = Liz)
check_model(hatch_dev_mod)
hist(residuals(hatch_dev_mod))
# temperature effects days to hatch; no effect on hormone or temp x hormone interaction
summary(hatch_dev_mod)
anova(hatch_dev_mod)
# faster hatch warmer temps
hatch_dev_mod_emm <- emmeans(hatch_dev_mod, pairwise ~ temp)
plot(hatch_dev_mod_emm)


################
#### 2.	What are the effects of developmental treatments hatching success?
#### none: only 2 animals did not hatch this year, can't make comparison
################


################
#### 3.	What are the effects of developmental treatments on body size, mass, and condition (BCI)--HATCHLING
#### SVL:hormone effects body size - high hormone low svl
#### MASS: hormone & temp effects body size - high hormone lower mass; cooler temps higher mass
#### BCI: poor BCI warm temps
################
#### body size (SVL) - hatching
SVL_hatch_mod <- lm(bd_svl_orig_mm ~ temp + hormone, data = Liz)
check_model(SVL_hatch_mod)
hist(resid(SVL_hatch_mod))
# hormone effects body size; no effects on temp or temp x hormone interaction
summary(SVL_hatch_mod)
anova(SVL_hatch_mod)
# high hormone lower SVL
SVL_hatch_mod_emm <- emmeans(SVL_hatch_mod, pairwise ~ hormone)
plot(SVL_hatch_mod_emm)

#### body mass - hatching; interaction removed
Mass_hatch_mod <- lm(bd_mass_orig ~ temp + hormone, data = Liz)
check_model(Mass_hatch_mod)
hist(resid(Mass_hatch_mod))
# temperature and hormone effects mass; no interaction effects
summary(Mass_hatch_mod)
anova(Mass_hatch_mod)
# hormone pairwise comparison - control higher mass
Mass_hormone_hatch_mod_emm <- emmeans(Mass_hatch_mod, pairwise ~ hormone)
plot(Mass_hormone_hatch_mod_emm)
# temp pairwise comparison - cooler temps higher mass
Mass_temp_hatch_mod_emm <- emmeans(Mass_hatch_mod, pairwise ~ temp)
plot(Mass_temp_hatch_mod_emm)

#### BCI from residuals (OLS)
# hatchlings- residuals
condition_hatch_fit <- lm(bd_mass_orig ~ bd_svl_orig_mm, data = Liz, na.action=na.exclude) 
Liz$hatch_bc_residuals <- residuals(condition_hatch_fit, na.action=na.exclude) 
## hatchlings-anova on BCI residuals; interaction removed
condition_hatch_mod <- lm(hatch_bc_residuals ~ temp + hormone, data = Liz)
check_model(condition_hatch_mod)
hist(resid(condition_hatch_mod))
# temperature has effect on BCI; no effect on hormone, or temp x hormone interaction
summary(condition_hatch_mod)
anova(condition_hatch_mod)
# temp and BCI plot - poor BCI warm temps
condition_hatch_mod_emm <- emmeans(condition_hatch_mod, pairwise ~ temp)
plot(condition_hatch_mod_emm)


################
#### 4.	What are the effects of developmental treatments on juvenile body size, condition, and mass (at the 2 points measured after hatching)?
#### Juv1_SVL: effect of temperature- cooler temps larger svl
#### Juv2_SVL: no effects
#### Juv1_MASS: low mass with high hormones; low mass with high temperatures
#### Juv2_MASS: no effects on treatments; marginal effect hormones p 0.06
#### Juv1_BCI: no treatment effects
#### Juv2_BCI: no treatment effects
################
#### body size (SVL) - juvenile remeasure 1 (adding days since hatch because effect of temp on development time); interaction removed
SVL_Juv1_mod <- lm(remeasure_SVL ~ temp + hormone + scale(days_posthatch_juv1), data = Liz)
check_model(SVL_Juv1_mod)
hist(resid(SVL_Juv1_mod))
# effect on temperature and day of hatch on body size 
summary(SVL_Juv1_mod)
anova(SVL_Juv1_mod)
# plot - cooler temps larger svl
SVL_Juv1_mod_emm <- emmeans(SVL_Juv1_mod, pairwise ~ temp)
plot(SVL_Juv1_mod_emm)

#### body size (SVL) - juvenile remeasure 2; interaction removed
SVL_Juv2_mod <- lm(measure_SVL_juv2 ~ temp + hormone  + scale(days_posthatch_juv2), data = Liz)
check_model(SVL_Juv2_mod)
hist(resid(SVL_Juv2_mod))
# marginal hermonal effect on body size p = 0.06
summary(SVL_Juv2_mod)
anova(SVL_Juv2_mod)


#### body mass - juvenile remeasure 1; interaction removed
Mass_Juv1_mod <- lm(remeasure_mass ~ temp + hormone + scale(days_posthatch_juv1), data = Liz)
check_model(Mass_Juv1_mod)
hist(resid(Mass_Juv1_mod))
# effect on temp, hormone, and days since hatch; no effects on temp x hormone interaction
summary(Mass_Juv1_mod)
anova(Mass_Juv1_mod)
# hormone plot - low mass with high hormones
Mass_Juv1_hormone_mod_emm <- emmeans(Mass_Juv1_mod, pairwise ~ hormone)
plot(Mass_Juv1_hormone_mod_emm)
# temp plot - low mass with high temperatures
Mass_Juv1_temp_mod_emm <- emmeans(Mass_Juv1_mod, pairwise ~ temp)
plot(Mass_Juv1_temp_mod_emm)

#### body mass - juvenile remeasure 2
Mass_Juv2_mod <- lm(measure_mass_juv2 ~ temp + hormone  + scale(days_posthatch_juv2), data = Liz)
check_model(Mass_Juv2_mod)
hist(resid(Mass_Juv2_mod))
# no effects on treatments
summary(Mass_Juv2_mod)
anova(Mass_Juv2_mod)


### BCI Juv1 - residuals
condition_juv1_fit <- lm(remeasure_mass ~ remeasure_SVL + scale(days_posthatch_juv1), data = Liz, na.action=na.exclude) # fit the model for residuals
Liz$juv1_bc_residuals <- residuals(condition_juv1_fit, na.action=na.exclude) # Save the residual values
## Juvenile -anova on BCI residuals; interaction removed
condition_juv_mod <- lm(juv1_bc_residuals ~ temp + hormone  + scale(days_posthatch_juv1), data = Liz)
# no treatment effects on bci juv1
summary(condition_juv_mod)
anova(condition_juv_mod) 
check_model(condition_juv_mod)

### BCI Juv2 - residuals
condition_juv2_fit <- lm(measure_mass_juv2 ~ measure_SVL_juv2 + scale(days_posthatch_juv2), data = Liz, na.action=na.exclude) # fit the model for residuals
Liz$juv2_bc_residuals <- residuals(condition_juv2_fit, na.action=na.exclude) # Save the residual values
## Juvenile -anova on BCI residuals; interaction removed
condition_juv2_mod <- lm(juv2_bc_residuals ~ temp + hormone +scale(days_posthatch_juv2), data = Liz)
check_model(condition_juv2_mod)
hist(resid(condition_juv2_mod))
# no treatment effects on bci juv1
summary(condition_juv2_mod)
anova(condition_juv2_mod) 


################
#### 4.	developmental treatment & survival; Testing days to mortality across treatments
#### no differences between treatments
#### incubation temp has mortality post day
################
# survival fisher test- filter out missing animals
survival<- Liz %>% 
  filter(liz_status_db != "MISSING")
survival_table <- table(survival$liz_status_db, survival$hormone)
fisher.test(survival_table)
# figure
survival_final <- survival %>% 
  group_by(liz_status_db, hormone) %>%
  summarise(n = n()) %>% 
  mutate(freq = n/sum(n))

mortality.fig <- ggplot(survival_final, aes(fill=hormone, y=n, x=liz_status_db)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values=c("#999999", "#E69F00", "brown2"))+
  ylab("% of survival | mortality") +
  xlab("Dead or Alive")
mortality.fig

# days difference 
survival_dat <- Liz %>% 
  filter(liz_status_db == "DEAD") %>% 
  select(c("temp", "hormone", "mortality_post_days")) %>% 
  drop_na() 

survival_day_hormone <- kruskal.test(mortality_post_days ~ hormone, data = survival_dat)
survival_day_temp <- kruskal.test(mortality_post_days ~ temp, data = survival_dat)




################
#### 5.	Growthrate
#### Juv1_SVL: effect of temperature- cooler temps larger svl
#### Juv2_SVL: no effects
#### Juv1_MASS: low mass with high hormones; low mass with high temperatures
#### Juv2_MASS: no effects on treatments; marginal effect hormones p 0.06
#### Juv1_BCI: no treatment effects
#### Juv2_BCI: no treatment effects

################
# NOTE: Liz data - 
growth_data <- Liz %>% 
  select(bd_liz_id, temp, hormone, 
         hatch_date_ymd, bd_svl_orig_mm, bd_mass_orig,
         remeasure_date, remeasure_SVL, remeasure_mass,
         measure_date_juv2, measure_SVL_juv2, measure_mass_juv2) %>% 
  rename(ID = bd_liz_id, initial_date = hatch_date_ymd, inital_svl_mm = bd_svl_orig_mm, inital_mass_g=bd_mass_orig,
         juv1_date = remeasure_date, juv1_svl_mm = remeasure_SVL, juv1_mass_g = remeasure_mass,
         juv2_date = measure_date_juv2, juv2_svl_mm = measure_SVL_juv2, juv2_mass_g = measure_mass_juv2)
# calculating days difference
growth_data$inital_to_juv1 <- as.numeric(difftime(growth_data$juv1_date,growth_data$initial_date, units = "days"))
growth_data$juv1_to_juv2 <- as.numeric(difftime(growth_data$juv2_date,growth_data$juv1_date, units = "days"))  
growth_data$inital_to_juv2 <- as.numeric(difftime(growth_data$juv2_date,growth_data$initial_date, units = "days"))




# growth rate calculations: calculated growth rates by 
# dividing change in SVL (or mass) between initial, 
# juv 1, juv 2(final) measurements/by the total number of days elapsed

# Growth 1 -  initial to 1st remeasurement (Juv_1)
growth1_inital_juv1 <- growth_data %>% 
  group_by(ID, temp, hormone) %>% 
  summarise(day_SVL_growth1_inital_to_juv1 = ((juv1_svl_mm-inital_svl_mm)/inital_to_juv1),
            day_mass_growth1_inital_to_juv1 = ((juv1_mass_g-inital_mass_g)/inital_to_juv1))


# Growth 1 SVL test; interaction removed
growth1_SVL_mod <- lm(day_SVL_growth1_inital_to_juv1 ~hormone + temp, data = growth1_inital_juv1, na.action=na.exclude)
check_model(growth1_SVL_mod)
hist(resid(growth1_SVL_mod))
# treatment effect on temp but not growth 
summary(growth1_SVL_mod)
anova(growth1_SVL_mod) 
# temp and growth plot - faster growth rates in cooler temps
growth1_SVL_mod_emm <- emmeans(growth1_SVL_mod, pairwise ~ temp)
plot(growth1_SVL_mod_emm)

# Growth 1 Mass test; interaction removed
growth1_mass_mod <- lm(day_mass_growth1_inital_to_juv1 ~hormone + temp, data = growth1_inital_juv1, na.action=na.exclude) 
check_model(growth1_mass_mod)
hist(resid(growth1_mass_mod))
# treatment effect on temp but not growth 
summary(growth1_mass_mod)
anova(growth1_mass_mod) 
# temp and growth plot - faster growth rates in cooler temps
growth1_mass_mod_emm <- emmeans(growth1_mass_mod, pairwise ~ temp)
plot(growth1_mass_mod_emm)


# Growth 2 is 1st remeasurement (Juv_1) to 2nd remeasurement
growth2_juv1_juv2 <- growth_data %>% 
  group_by(ID, temp, hormone) %>% 
  summarise(day_SVL_growth2_juv1_to_juv2 = ((juv2_svl_mm - juv1_svl_mm)/juv1_to_juv2),
            day_mass_growth2_juv1_to_juv2 = ((juv2_mass_g-juv1_mass_g)/juv1_to_juv2))

# Growth 2 test SVL; interaction removed
growth2_SVL_mod <- lm(day_SVL_growth2_juv1_to_juv2 ~hormone + temp, data = growth2_juv1_juv2, na.action=na.exclude) 
check_model(growth2_SVL_mod)
hist(resid(growth2_SVL_mod))
# treatment effect on temp but not growth 
summary(growth2_SVL_mod)
anova(growth2_SVL_mod) 
# temp and growth plot - faster growth rates in warmer temps
growth2_SVL_mod_emm <- emmeans(growth2_SVL_mod, pairwise ~ temp)
plot(growth2_SVL_mod_emm)

# Growth 2 Mass test; interaction removed
growth2_mass_mod <- lm(day_mass_growth2_juv1_to_juv2 ~hormone + temp, data = growth2_juv1_juv2, na.action=na.exclude) 
check_model(growth2_mass_mod)
hist(resid(growth2_mass_mod))
# treatment effect on temp but not growth 
summary(growth2_mass_mod)
anova(growth2_mass_mod) 
# temp and growth plot - faster growth rates in warmer temps
growth2_mass_mod_emm <- emmeans(growth2_mass_mod, pairwise ~ temp)
plot(growth2_mass_mod_emm)


# Growth 3 - overall study growth rate - Initial to 2nd remeasurement 
growth3_study <- growth_data %>% 
  group_by(ID, temp, hormone) %>% 
  summarise(day_SVL_growth3_inital_to_juv2 = ((juv2_svl_mm-inital_svl_mm)/inital_to_juv2),
            day_mass_growth3_inital_to_juv2 = ((juv2_mass_g-inital_mass_g)/inital_to_juv2))
# Growth 3 test SVL; interaction removed
growth3_SVL_mod <- lm(day_SVL_growth3_inital_to_juv2 ~hormone + temp, data = growth3_study, na.action=na.exclude) 
check_model(growth3_SVL_mod)
hist(resid(growth3_SVL_mod))
# no treatment effect 
summary(growth3_SVL_mod)
anova(growth3_SVL_mod) 

# Growth 2 Mass test; interaction removed
growth3_mass_mod <- lm(day_mass_growth3_inital_to_juv2 ~hormone + temp, data = growth3_study, na.action=na.exclude) 
check_model(growth3_mass_mod)
hist(resid(growth3_mass_mod))
# no treatment effects 
summary(growth3_mass_mod)
anova(growth3_mass_mod) 


#####################
### Figures
#####################
hormones <- c("#999999", "#E69F00", "brown2")
temps <- c("Blue", "Red")
my_comparisons <- rev(list(c("control","low"),c("control","high-"),c("low","high")))
# Figure1 : effects of developmental treatments (temp, cort, interaction) on time to hatch: -	Effect on temperature: warmer temps faster development 
Fig1 <- ggboxplot(Liz, x = "temp", y = "days_to_hatch",
                   color = "temp", palette = temps,
                   short.panel.labs = FALSE,
                   font.label = list(size = 14, color = "black"),
                   ggtheme = theme_bw())+ 
  stat_compare_means(method = "anova", vjust = -.1, hjust = -3)+ 
  labs(x = "Temperature (°C)", y = "Incubation time (days)") +
  labs(color='Incubation Temperature')+
  scale_y_continuous(breaks=seq(25, 60, 5), limits = c(25,60))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
Fig1


# Figure 2: effects of developmental treatments on body size, condition, and mass at hatching;
# Hormone effect svl; mass hormone and temp effects body size; Temp effects BCI 
# MODS USED:  
# SVL_hatch_mod <- lm(bd_svl_orig_mm ~ temp + hormone, data = Liz) 
# Mass_hatch_mod <- lm(bd_mass_orig ~ temp + hormone, data = Liz) 
# condition_hatch_mod <- lm(hatch_bc_residuals ~ temp + hormone, data = Liz)
# A)
fig2.svl <- ggboxplot(Liz, x = "hormone", y = "bd_svl_orig_mm",
                      color = "hormone", palette = hormones,
                      short.panel.labs = FALSE,
                      font.label = list(size = 14, color = "black"),
                      ggtheme = theme_bw())+ 
  geom_jitter(aes(color = hormone), width = 0.2, alpha = 0.4) +  # Adding raw data points with color and transparency
  stat_compare_means(comparisons = list(c("control","low"), c("control","high"), c("low","high"))) +
  stat_compare_means(method = "anova")+
  labs(x = NULL, y = "Hatchling SVL") +
  labs(color='Hormone Treatment')+
  scale_y_continuous(breaks=seq(14, 22, 2), limits = c(14, 22))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
fig2.svl



# B)
fig2.bci <- ggboxplot(Liz, x = "temp", y = "hatch_bc_residuals",
                      color = "temp", palette = temps,
                      short.panel.labs = FALSE,
                      font.label = list(size = 14, color = "black"),
                      ggtheme = theme_bw()) + 
  geom_jitter(aes(color = temp), width = 0.2, alpha = 0.4) +  # Adding raw data points with color and transparency
  stat_compare_means(method = "anova", vjust = -.1, hjust = -3)+ 
  labs(x = "Temperature (°C)", y = "BCI") +
  labs(color='Incubation Treatment')+
  scale_y_continuous(breaks=seq(-0.06, 0.06, 0.02), limits = c(-0.06, 0.06))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
fig2.bci



# C)
Fig2.mass.h <- ggboxplot(Liz, x = "hormone", y = "bd_mass_orig",
                         color = "hormone", palette = hormones,
                         short.panel.labs = FALSE,
                         font.label = list(size = 14, color = "black"),
                         ggtheme = theme_bw()) +
  geom_jitter(aes(color = hormone), width = 0.2, alpha = 0.4) +  # Adding raw data points with color and transparency
  stat_compare_means(comparisons = list(c("control","low"), c("control","high"), c("low","high"))) +
  stat_compare_means(method = "anova") +
  labs(x = NULL, y = "Hatchling mass (g)") +
  labs(color='Hormone Treatment')+
  scale_y_continuous(breaks=seq(0.06, 0.20, .02), limits = c(0.06,0.20))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

Fig2.mass.h

# D)
Fig2.mass.t <- ggboxplot(Liz, x = "temp", y = "bd_mass_orig",
                         color = "temp", palette = temps,
                         short.panel.labs = FALSE,
                         font.label = list(size = 14, color = "black"),
                         ggtheme = theme_bw()) +
  geom_jitter(aes(color = temp), width = 0.2, alpha = 0.4) +  # Adding raw data points with color and transparency
  stat_compare_means(method = "anova") +
  labs(x = "Temperature (°C)", y = "Hatchling mass (g)") +
  labs(color='Temperature Treatment')+
  scale_y_continuous(breaks=seq(0.06, 0.18, .02), limits = c(0.06,0.18)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
Fig2.mass.t






ggplot(Liz, aes(x = hormone, y = bd_mass_orig, fill = hormone))+
  geom_violin(alpha=0.4, position = position_dodge(width = .9),size=1, color="black", trim = FALSE)+ 
  geom_point()+
  stat_compare_means(comparisons = list(c("control","low"), c("control","high"), c("low","high")))+
  stat_compare_means()+
  ggforce::geom_sina(alpha=0.4)

