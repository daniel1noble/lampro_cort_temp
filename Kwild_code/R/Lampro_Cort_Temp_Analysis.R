## Kwild Lampro Cort Temp Analysis
# WD, packages, data
pacman::p_load(dplyr, tidyverse, ggpubr, lme4, emmeans, car, lmerTest, MuMIn, glmm, installr, lubridate, performance, cowplot)


#########################
# Bring in data and merge
#########################
# 1) bring in data base data that has morphology and mortality data on across life stage: hatch, juv_1, juv_2
dat_hatch_juv1_2 <- read.csv ("Kwild_code/data/Cort_hatch_juv1_2.csv")

# 2) bring in data that has final morphology measurement and hormone data: "juv_3"
dat_juv3 <- read.csv(file = "Kwild_code/data/Cort_juv_3.csv") %>% 
  rename(sex = juv3_Sex)

# 3) final merge juv 3 data that contains hormonal results and final morph with dat_hatch_juv1_2 
merged_data <- merge(dat_hatch_juv1_2, dat_juv3, by.x = "Lizard_ID", 
                     by.y = "juv3_Lizard_ID", all = TRUE)


#########################
# "Issue data" from merge
#########################
# 1) find "new" values in Juv3: 
new_ids_juv3 <- dat_juv3$juv3_Lizard_ID[!dat_juv3$juv3_Lizard_ID %in% dat_hatch_juv1_2$Lizard_ID]
new_ids_juv3 # LD745_21???

# 2) Get IDs from dat_hatch_juv1_2 that are not in dat_juv3
missing_ids_juv3 <- dat_hatch_juv1_2$Lizard_ID[!dat_hatch_juv1_2$Lizard_ID %in% dat_juv3$juv3_Lizard_ID]
# Retrieve all rows from dat_hatch_juv1_2 where Lizard_ID is in missing_ids
missing_ids_dat <- dat_hatch_juv1_2[dat_hatch_juv1_2$Lizard_ID %in% missing_ids_juv3, ]
missing_ids_dat <- missing_ids_dat %>% 
  filter(liz_status_db == "ALIVE")
table(missing_ids_dat$liz_status_db) # 15 Dead # 21 Alive since juv 2 sample


#########################
# dates  and days post hatch days
#########################
# formatting dates
merged_data$hatch_date <- dmy(merged_data$hatch_date_dmy) # hatch
merged_data$juv1_measure_date <- dmy(merged_data$juv1_measure_date_dmy) # juv1
merged_data$juv2_measure_date <- dmy(merged_data$juv2_measure_date_dmy) # juv2
merged_data$juv3_measure_date <- dmy(merged_data$juv3_date) # juv3
merged_data$mortality_date <- dmy(merged_data$mortality_date_dmy) # mortality date - missing final

# days post hatch used to scale in models
# juv 1 days post hatch
merged_data$days_post_hatch_juv1 <- difftime(merged_data$juv1_measure_date,
                                           merged_data$hatch_date, units = "days")
# juv 2 days post hatch
merged_data$days_post_hatch_juv2 <- difftime(merged_data$juv2_measure_date,
                                           merged_data$hatch_date, units = "days")
# juv 3 days post hatch
merged_data$days_post_hatch_juv3 <- difftime(merged_data$juv3_measure_date,
                                           merged_data$hatch_date, units = "days")

# naming final df and ordering for figures
data_final <-  mutate(merged_data, hormone = factor(hormone, 
                                        levels = c("control",
                                                   "low",
                                                   "high"))) %>% 
  filter(!is.na(hormone)) %>% 
  group_by(hormone)


#########################
# data arrangement up for growth estimations that will be used later
#########################
# calculating days difference between measurments 
data_final$inital_to_juv1 <- as.numeric(difftime(data_final$juv1_measure_date,
                                                 data_final$hatch_date, 
                                                 units = "days"))
data_final$juv1_to_juv2 <- as.numeric(difftime(data_final$juv2_measure_date,
                                               data_final$juv1_measure_date, 
                                               units = "days"))  
data_final$juv2_to_juv3 <- as.numeric(difftime(data_final$juv3_measure_date,
                                               data_final$juv2_measure_date, 
                                               units = "days"))
data_final$inital_to_juv3 <- as.numeric(difftime(data_final$juv3_measure_date,
                                                 data_final$hatch_date, 
                                                 units = "days"))
########
# growth rate calculations: calculated growth rates by 
# dividing change in SVL (or mass) between initial, 
# juv 1, juv 2, juv 3(final) measurements/by the total number of days elapsed
data_final <- data_final %>% 
  mutate(day_SVL_growth1_inital_to_juv1 = ((juv1_measure_SVL_mm-hatch_svl_orig_mm)/inital_to_juv1), 
         day_mass_growth1_inital_to_juv1 = ((juv1_measure_mass_g-hatch_mass_orig)/inital_to_juv1),
         day_SVL_growth2_juv1_to_juv2 = ((juv2_measure_SVL_mm - juv1_measure_SVL_mm)/juv1_to_juv2),
         day_mass_growth2_juv1_to_juv2 = ((juv2_measure_mass_g-juv1_measure_mass_g)/juv1_to_juv2),
         day_SVL_growth3_juv2_to_juv3 = ((juv3_SVL_mm-juv2_measure_SVL_mm)/juv2_to_juv3),
         day_mass_growth3_juv2_to_juv3 = ((juv3_mass_g-juv2_measure_mass_g)/juv2_to_juv3), 
         day_SVL_growth4_inital_to_juv3 = ((juv3_SVL_mm-hatch_svl_orig_mm)/inital_to_juv3),
         day_mass_growth4_inital_to_juv3 = ((juv3_mass_g-hatch_mass_orig)/inital_to_juv3))



################
#### 1.	What are the effects of developmental treatments (temp, cort, interaction) on time to hatch?
#### effect of temperature on incubation days - faster hatch warmer temps 
################
# mod ; interaction removed
hatch_dev_mod <- lm(days_to_hatch~ temp + hormone, data = data_final)
check_model(hatch_dev_mod) # bimodal because there is an effect on temp to days of hatch

# temperature effects days to hatch; no effect on hormone or temp x hormone interaction
summary(hatch_dev_mod)
anova(hatch_dev_mod)
# faster hatch warmer temps
hatch_dev_mod_emm <- emmeans(hatch_dev_mod, pairwise ~ temp)
plot(hatch_dev_mod_emm)
saveRDS(hatch_dev_mod, "Kwild_code/models/hatch_dev_mod.RDS")


################
#### 2.	What are the effects of developmental treatments on body size, mass, and condition (BCI) of HATCHLINGS
#### SVL: hormone effects body size; no effects on temp or temp x hormone interaction
#### MASS: temperature and hormone effects mass; no interaction effects
#### BCI:  temperature has effect on BCI; no effect on hormone, or temp x hormone interaction

################
#### SVL; interaction removed
SVL_hatch_mod <- lm(hatch_svl_orig_mm ~ temp + hormone , data = data_final)
check_model(SVL_hatch_mod)
# hormone effects body size; no effects on temp or temp x hormone interaction
summary(SVL_hatch_mod)
anova(SVL_hatch_mod)
# high hormone lower SVL
SVL_hatch_mod_emm <- emmeans(SVL_hatch_mod, pairwise ~ hormone)
plot(SVL_hatch_mod_emm)
saveRDS(SVL_hatch_mod, "Kwild_code/models/SVL_hatch_mod.RDS")

#### MASS; interaction removed
Mass_hatch_mod <- lm(hatch_mass_orig ~ temp + hormone, data = data_final)
check_model(Mass_hatch_mod)
# temperature and hormone effects mass; no interaction effects
summary(Mass_hatch_mod)
anova(Mass_hatch_mod)
# hormone pairwise comparison - control higher mass
Mass_hormone_hatch_mod_emm <- emmeans(Mass_hatch_mod, pairwise ~ hormone)
plot(Mass_hormone_hatch_mod_emm)
# temp pairwise comparison - cooler temps higher mass
Mass_temp_hatch_mod_emm <- emmeans(Mass_hatch_mod, pairwise ~ temp)
plot(Mass_temp_hatch_mod_emm)
saveRDS(Mass_hatch_mod, "Kwild_code/models/MASS_hatch_mod.RDS")

#### BCI from residuals (OLS)
# hatchlings- residuals
condition_hatch_fit <- lm(hatch_mass_orig ~ hatch_svl_orig_mm, data = data_final, na.action=na.exclude) 
data_final$hatch_bc_residuals <- residuals(condition_hatch_fit, na.action=na.exclude) 
## hatchlings-anova on BCI residuals; interaction removed
condition_hatch_mod <- lm(hatch_bc_residuals ~ temp + hormone, data = data_final)
check_model(condition_hatch_mod)
# temperature has effect on BCI; no effect on hormone, or temp x hormone interaction
summary(condition_hatch_mod)
anova(condition_hatch_mod)
# temp and BCI plot - poor BCI warm temps
condition_hatch_mod_emm <- emmeans(condition_hatch_mod, pairwise ~ temp)
plot(condition_hatch_mod_emm)
saveRDS(Mass_hatch_mod, "Kwild_code/models/condition_hatch_mod.RDS")



################
#### 3.	What are the effects of developmental treatments on juvenile body size, condition, and mass (at the 3 points measured after hatching)?
#### Juv1_SVL: effect of temperature- cooler temps larger svl
#### Juv2_SVL: no effects - hormonal effect p = 0.07
#### Juv3_SVL: high hormone low body size; males smaller than females
#### Juv1_MASS: low mass with high hormones; low mass with high temperatures
#### Juv2_MASS: low mass with high temperatures
#### Juv3_MASS: low mass males, and low mass high temp
#### Juv1_BCI: no treatment effects
#### Juv2_BCI: no treatment effects
#### Juv3_BCI: Sex effect: male lower BCI than females - could be size difference

###############
# SVL ANALYSIS: Juv_1, Juv_2, Juv_3 (sex accounted for last measurment)
###############
#### JUV_1: SVL - adding days since hatch because effect of temp on development time); interaction removed
SVL_Juv1_mod <- lm(juv1_measure_SVL_mm ~ temp + hormone + scale(days_post_hatch_juv1), data = data_final)
check_model(SVL_Juv1_mod)
# effect on temperature and day of hatch on body size 
summary(SVL_Juv1_mod)
# plot - cooler temps larger svl
SVL_Juv1_mod_emm <- emmeans(SVL_Juv1_mod, pairwise ~ temp)
plot(SVL_Juv1_mod_emm)
saveRDS(SVL_Juv1_mod, "Kwild_code/models/SVL_Juv1_mod.RDS")

#### JUV_2: SVL; interaction removed
SVL_Juv2_mod <- lm(juv2_measure_SVL_mm ~ temp + hormone  + scale(days_post_hatch_juv2), data = data_final)
check_model(SVL_Juv2_mod)
# marginal hermonal effect on body size p = 0.06
summary(SVL_Juv2_mod)
anova(SVL_Juv2_mod)
saveRDS(SVL_Juv2_mod, "Kwild_code/models/SVL_Juv2_mod.RDS")

#### SVL - NOTE SEX ACCOUTNED; no interaction
SVL_Juv3_mod <- lm(juv3_SVL_mm ~ temp + hormone + sex + scale(days_post_hatch_juv3), data = data_final)
check_model(SVL_Juv3_mod)
# effect on hormones and sex
summary(SVL_Juv3_mod)
anova(SVL_Juv3_mod)
# plot - high hormone lower svl
SVL_Juv3_hormone_mod_emm <- emmeans(SVL_Juv3_mod, pairwise ~ hormone)
plot(SVL_Juv3_hormone_mod_emm)
# plot - males are smaller than girls
SVL_Juv3_sex_mod_emm <- emmeans(SVL_Juv3_mod, pairwise ~ sex)
plot(SVL_Juv3_sex_mod_emm)
saveRDS(SVL_Juv3_mod, "Kwild_code/models/SVL_Juv3_mod.RDS")


###############
# MASS ANALYSIS: Juv_1, Juv_2, Juv_3
###############
#### MASS: Juv_1 - interaction removed
Mass_Juv1_mod <- lm(juv1_measure_mass_g ~ temp + hormone + scale(days_post_hatch_juv1), data = data_final)
check_model(Mass_Juv1_mod)
# effect on temp, hormone, and days since hatch; no effects on temp x hormone interaction
summary(Mass_Juv1_mod)
anova(Mass_Juv1_mod)
# hormone plot - low mass with high hormones
Mass_Juv1_hormone_mod_emm <- emmeans(Mass_Juv1_mod, pairwise ~ hormone)
plot(Mass_Juv1_hormone_mod_emm)
# temp plot - low mass with high temperatures
Mass_Juv1_temp_mod_emm <- emmeans(Mass_Juv1_mod, pairwise ~ temp)
plot(Mass_Juv1_temp_mod_emm)
saveRDS(Mass_Juv1_mod, "Kwild_code/models/Mass_Juv1_mod.RDS")

#### MASS: Juv_2; interaction removed
Mass_Juv2_mod <- lm(juv2_measure_mass_g ~ temp + hormone  + scale(days_post_hatch_juv2), data = data_final)
check_model(Mass_Juv2_mod)
# temp effect
summary(Mass_Juv2_mod)
anova(Mass_Juv2_mod)
# temp plot - low mass with high temperatures
Mass_Juv2_temp_mod_emm <- emmeans(Mass_Juv2_mod, pairwise ~ temp)
plot(Mass_Juv2_temp_mod_emm)
saveRDS(Mass_Juv2_mod, "Kwild_code/models/Mass_Juv2_mod.RDS")

#### MASS: Juv_3; NOTE SEX IS ACCOUNTED FOR
Mass_Juv3_mod <- lm(juv3_mass_g ~ temp + hormone+ sex+ scale(days_post_hatch_juv2), data = data_final)
check_model(Mass_Juv3_mod)
# temp effect
summary(Mass_Juv3_mod)
anova(Mass_Juv3_mod)
# temp plot - low mass with high temperature
Mass_Juv3_temp_mod_emm <- emmeans(Mass_Juv3_mod, pairwise ~ temp)
plot(Mass_Juv3_temp_mod_emm)
# temp plot - males smaller than females
Sex_Juv3_temp_mod_emm <- emmeans(Mass_Juv3_mod, pairwise ~ sex)
plot(Sex_Juv3_temp_mod_emm)
saveRDS(Mass_Juv3_mod, "Kwild_code/models/Mass_Juv3_mod.RDS")


###############
# BCI ANALYSIS: Juv_1, Juv_2, Juv_3
###############
### BCI Juv1 - residuals
condition_juv1_fit <- lm(juv1_measure_mass_g ~ juv1_measure_SVL_mm + scale(days_post_hatch_juv1), data = data_final, na.action=na.exclude) # fit the model for residuals
data_final$juv1_bc_residuals <- residuals(condition_juv1_fit, na.action=na.exclude) # Save the residual values
## Juvenile -anova on BCI residuals; interaction removed
condition_juv_mod <- lm(juv1_bc_residuals ~ temp + hormone  + scale(days_post_hatch_juv1), data = data_final)
# no treatment effects on bci juv1
summary(condition_juv_mod)
anova(condition_juv_mod) 
check_model(condition_juv_mod)
# temp plot - 
condition_juv_mod_emm <- emmeans(condition_juv_mod, pairwise ~ temp)
plot(condition_juv_mod_emm)
saveRDS(condition_juv_mod, "Kwild_code/models/condition_juv_mod.RDS")

### BCI Juv2 - residuals
condition_juv2_fit <- lm(juv2_measure_mass_g ~ juv2_measure_SVL_mm + scale(days_post_hatch_juv2), data = data_final, na.action=na.exclude) # fit the model for residuals
data_final$juv2_bc_residuals <- residuals(condition_juv2_fit, na.action=na.exclude) # Save the residual values
## Juvenile -anova on BCI residuals; interaction removed
condition_juv2_mod <- lm(juv2_bc_residuals ~ temp + hormone +scale(days_post_hatch_juv2), data = data_final)
check_model(condition_juv2_mod)
# no treatment effects on bci juv1
summary(condition_juv2_mod)
anova(condition_juv2_mod) 
saveRDS(condition_juv2_mod, "Kwild_code/models/condition_juv2_mod.RDS")

### BCI Juv3 - residuals
condition_juv3_fit <- lm(juv3_mass_g ~ juv3_SVL_mm + scale(days_post_hatch_juv3), data = data_final, na.action=na.exclude) # fit the model for residuals
data_final$juv3_bc_residuals <- residuals(condition_juv3_fit, na.action=na.exclude) # Save the residual values
## Juvenile -anova on BCI residuals; NOTE SEX INCLUDED
condition_juv3_mod <- lm(juv3_bc_residuals ~ temp + hormone + sex + scale(days_post_hatch_juv3), data = data_final)
check_model(condition_juv3_mod)
# no treatment effects on bci juv1
summary(condition_juv3_mod)
anova(condition_juv3_mod) 
# sex plot - 
sex_condition_juv3_mod_emm <- emmeans(condition_juv3_mod, pairwise ~ sex)
plot(sex_condition_juv3_mod_emm)
saveRDS(condition_juv3_mod, "Kwild_code/models/condition_juv3_mod.RDS")


################
#### 4.	developmental treatment & survival; Testing days to mortality across treatments
#### no differences between treatments
################
# survival fisher test- filter out missing animals
survival<- data_final %>% 
  filter(liz_status_db != "MISSING")
survival_table <- table(survival$liz_status_db, survival$hormone)
fisher.test(survival_table)

# days difference 
survival_dat <- data_final %>% 
  filter(liz_status_db == "DEAD") %>% 
  select(c("temp", "hormone", "mortality_date_dmy")) %>% 
  drop_na() 

survival_day_hormone <- kruskal.test(mortality_date_dmy ~ hormone, data = survival_dat)
survival_day_temp <- kruskal.test(mortality_date_dmy ~ temp, data = survival_dat)


################
#### 5.	Growthrate
#### OVERALL SVL: effect of temperature- warm temps faster svl growth
#### OVERALL MASS: effect of temperature-warm temps slower mass growth
################
# SVL overall growth: hatch to juv_3; interaction removed
overall_SVL_growth_mod <- lm(day_SVL_growth4_inital_to_juv3 ~hormone + temp, 
                      data = data_final, na.action=na.exclude)
check_model(overall_SVL_growth_mod)
# treatment effect on temp but not hormone
summary(overall_SVL_growth_mod)
anova(overall_SVL_growth_mod) 
# temp and growth plot - faster growth rates in cooler temps
overall_SVL_growth_mod_emm <- emmeans(overall_SVL_growth_mod, pairwise ~ temp)
plot(overall_SVL_growth_mod_emm)
saveRDS(overall_SVL_growth_mod, "Kwild_code/models/overall_SVL_growth_mod.RDS")

# Mass test; interaction removed
growth4_mass_mod <- lm(day_SVL_growth4_inital_to_juv3 ~hormone + temp, 
                       data = data_final, na.action=na.exclude) 
check_model(growth4_mass_mod)
# treatment effect on temp but not hormone
summary(growth4_mass_mod)
anova(growth4_mass_mod) 
# temp and growth plot - faster growth rates in cooler temps
growth4_mass_mod_emm <- emmeans(growth4_mass_mod, pairwise ~ temp)
plot(growth4_mass_mod_emm)
saveRDS(growth4_mass_mod, "Kwild_code/models/growth4_mass_mod.RDS")


########################
# Testosterone, T4, Cort
########################
hormone_analysis_dat <- data_final %>% 
  filter(!is.na(juv3_CORT_Final_Hormone_ng_mL))


# cort differences
Cort_sex <- lmer(juv3_CORT_Final_Hormone_ng_mL ~Season + Sex + Season*Sex
                 + (1|Liz_ID), data = move.dat)
anova(Raw.Avg.mixed.lmer)  









#####################
### Figures
#####################
hormones <- c("#999999", "#E69F00", "brown2")
temps <- c("Blue", "Red")
my_comparisons <- rev(list(c("control","low"),c("control","high-"),c("low","high")))
# Figure1 : effects of developmental treatments (temp, cort, interaction) on time to hatch: -	Effect on temperature: warmer temps faster development 
Fig1 <- ggboxplot(data_final, x = "temp", y = "days_to_hatch",
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
# SVL_hatch_mod <- lm(hatch_svl_orig_mm ~ temp + hormone, data = Liz) 
# Mass_hatch_mod <- lm(hatch_mass_orig ~ temp + hormone, data = Liz) 
# condition_hatch_mod <- lm(hatch_bc_residuals ~ temp + hormone, data = Liz)
# A)
fig2_A <- ggboxplot(data_final, x = "hormone", y = "hatch_svl_orig_mm",
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
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none")




# B)
Fig2_B <- ggboxplot(data_final, x = "hormone", y = "hatch_mass_orig",
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
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = c(0.99, 0.99),  # Positions the legend inside the plot
        legend.justification = c("right", "top"), # Justifies the legend's position
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.box = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8))

Fig2 <- plot_grid(fig2_A, Fig2_B, ncol = 2, align = "h")

# Fig 3 Temperature
fig3_A <- ggboxplot(data_final, x = "temp", y = "hatch_bc_residuals",
                      color = "temp", palette = temps,
                      short.panel.labs = FALSE,
                      font.label = list(size = 14, color = "black"),
                      ggtheme = theme_bw()) + 
  geom_jitter(aes(color = temp), width = 0.2, alpha = 0.4) +  # Adding raw data points with color and transparency
  stat_compare_means(method = "anova", vjust = -.1, hjust = -2)+ 
  labs(x = "Temperature (°C)", y = "BCI") +
  labs(color='Incubation Treatment')+
  scale_y_continuous(breaks=seq(-0.06, 0.06, 0.02), limits = c(-0.06, 0.06))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



# D)
fig3_B <- ggboxplot(data_final, x = "temp", y = "hatch_mass_orig",
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




fig3 <- plot_grid(fig3_A, fig3_B, ncol = 2, align = "h")


ggplot(Liz, aes(x = hormone, y = bd_mass_orig, fill = hormone))+
  geom_violin(alpha=0.4, position = position_dodge(width = .9),size=1, color="black", trim = FALSE)+ 
  geom_point()+
  stat_compare_means(comparisons = list(c("control","low"), c("control","high"), c("low","high")))+
  stat_compare_means()+
  ggforce::geom_sina(alpha=0.4)



str(liz_ondi)
