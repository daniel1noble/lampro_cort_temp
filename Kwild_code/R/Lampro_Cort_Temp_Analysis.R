## Kwild Lampro Cort Temp Analysis
# WD, packages, data
pacman::p_load(dplyr, tidyverse, ggpubr, lme4, emmeans, car, lmerTest, MuMIn, glmm, installr, lubridate, performance)


#########################
# Bring in data and merge
#########################
# 1) bring in data base data that has morphology and mortality data on across life stage: hatch, juv_1, juv_2
dat_hatch_juv1_2 <- read.csv ("Kwild_code/data/Pilot_cort_hatch_juv1_2.csv")

# 2) bring in data that has final morphology measurement and hormone data: "juv_3"
dat_juv3 <- read.csv(file = "Kwild_code/data/Pilot_cort_juv_3.csv") %>% 
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
# dates and age 
#########################
# formatting dates
merged_data$hatch_date_ymd <- dmy(merged_data$hatch_date_dmy) # hatch
merged_data$remeasure_date <- dmy(merged_data$juv1_measure_date_dmy) # juv1
merged_data$measure_date_juv2 <- dmy(merged_data$juv2_measure_date_dmy) # juv2
merged_data$measure_date_juv3 <- dmy(merged_data$juv3_date) # juv3
merged_data$mortality_date_dmy <- dmy(merged_data$mortality_date_dmy) # mortality date - missing final data follow up
# Ages
merged_data$days_posthatch_juv1 <- difftime(merged_data$juv1_measure_date_dmy,
                                            merged_data$hatch_date_ymd, units = "days") # age to juv1
merged_data$days_posthatch_juv2 <- difftime(merged_data$juv2_measure_date_dmy,
                                            merged_data$hatch_date_ymd, units = "days") # age to juv2
merged_data$days_posthatch_juv3 <- difftime(merged_data$juv3_date,
                                            merged_data$hatch_date_ymd, units = "days") # age to juv3

# naming final df and ordering for figures
data_final <-  mutate(merged_data, hormone = factor(hormone, 
                                        levels = c("control",
                                                   "low",
                                                   "high"))) %>% 
  group_by(hormone)






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
#### Juv2_SVL: no effects
#### Juv3_SVL: hormonal effect and sex effect
#### Juv1_MASS: low mass with high hormones; low mass with high temperatures
#### Juv2_MASS: no effects on treatments; marginal effect hormones p 0.065
#### Juv3_MASS: XXXXXXXXXX
#### Juv1_BCI: no treatment effects
#### Juv2_BCI: no treatment effects
#### Juv3_BCI: XXXXXXXX
################
#### JUV_1: body size (SVL) - juvenile remeasure 1 (adding days since hatch because effect of temp on development time); interaction removed
SVL_Juv1_mod <- lm(juv1_measure_SVL_mm ~ temp + hormone + scale(days_posthatch_juv1), data = data_final)
check_model(SVL_Juv1_mod)
# effect on temperature and day of hatch on body size 
summary(SVL_Juv1_mod)
anova(SVL_Juv1_mod)
# plot - cooler temps larger svl
SVL_Juv1_mod_emm <- emmeans(SVL_Juv1_mod, pairwise ~ temp)
plot(SVL_Juv1_mod_emm)
saveRDS(SVL_Juv1_mod, "Kwild_code/models/SVL_Juv1_mod.RDS")

#### JUV_2: body size (SVL) - juvenile remeasure 2; interaction removed
SVL_Juv2_mod <- lm(juv2_measure_SVL_mm ~ temp + hormone  + scale(days_posthatch_juv2), data = data_final)
check_model(SVL_Juv2_mod)
# marginal hermonal effect on body size p = 0.06
summary(SVL_Juv2_mod)
anova(SVL_Juv2_mod)
saveRDS(SVL_Juv2_mod, "Kwild_code/models/SVL_Juv2_mod.RDS")

#### JUV_3: body size (SVL) - NOTE SEX ACCOUTNED FOR; no interaction
SVL_Juv3_mod <- lm(juv3_SVL_mm ~ temp + hormone + sex + scale(days_posthatch_juv3), data = data_final)
check_model(SVL_Juv3_mod)
# effect on hormones and sex
summary(SVL_Juv3_mod)
anova(SVL_Juv3_mod)
# plot - high hormone lower svl
SVL_Juv3_hormone_mod_emm <- emmeans(SVL_Juv3_mod, pairwise ~ hormone)
plot(SVL_Juv3_mod_emm)
# plot - males are smaller than girls
SVL_Juv3_sex_mod_emm <- emmeans(SVL_Juv3_mod, pairwise ~ sex)
plot(SVL_Juv3_sex_mod_emm)
saveRDS(SVL_Juv3_mod, "Kwild_code/models/SVL_Juv3_mod.RDS")





















#### Juv1_MASS - interaction removed
Mass_Juv1_mod <- lm(juv1_measure_mass_g ~ temp + hormone + scale(days_posthatch_juv1), data = data_final)
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

#### body mass - juvenile remeasure 2
Mass_Juv2_mod <- lm(measure_mass_juv2 ~ temp + hormone  + scale(days_posthatch_juv2), data = data_final)
check_model(Mass_Juv2_mod)
hist(resid(Mass_Juv2_mod))
# no effects on treatments
summary(Mass_Juv2_mod)
anova(Mass_Juv2_mod)


### BCI Juv1 - residuals
condition_juv1_fit <- lm(remeasure_mass ~ juv1_measure_SVL_mm + scale(days_posthatch_juv1), data = data_final, na.action=na.exclude) # fit the model for residuals
data_final$juv1_bc_residuals <- residuals(condition_juv1_fit, na.action=na.exclude) # Save the residual values
## Juvenile -anova on BCI residuals; interaction removed
condition_juv_mod <- lm(juv1_bc_residuals ~ temp + hormone  + scale(days_posthatch_juv1), data = data_final)
# no treatment effects on bci juv1
summary(condition_juv_mod)
anova(condition_juv_mod) 
check_model(condition_juv_mod)

### BCI Juv2 - residuals
condition_juv2_fit <- lm(measure_mass_juv2 ~ measure_SVL_juv2 + scale(days_posthatch_juv2), data = data_final, na.action=na.exclude) # fit the model for residuals
data_final$juv2_bc_residuals <- residuals(condition_juv2_fit, na.action=na.exclude) # Save the residual values
## Juvenile -anova on BCI residuals; interaction removed
condition_juv2_mod <- lm(juv2_bc_residuals ~ temp + hormone +scale(days_posthatch_juv2), data = data_final)
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
survival<- data_final %>% 
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
survival_dat <- data_final %>% 
  filter(liz_status_db == "DEAD") %>% 
  select(c("temp", "hormone", "mortality_post_days")) %>% 
  drop_na() 

survival_day_hormone <- kruskal.test(mortality_post_days ~ hormone, data = survival_dat)
survival_day_temp <- kruskal.test(mortality_post_days ~ temp, data = survival_dat)




