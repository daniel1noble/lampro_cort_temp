## Kwild Lampro Cort Temp Analysis
# WD, packages, data
pacman::p_load(dplyr, tidyverse, ggpubr, lme4, emmeans, here, plotly, performance, car)


library(lubridate)
library(Hmisc)
library(ppcor)
#########################
# Bring in data from samples and merge together
#########################
# 1) bring in data base data that has morphology and mortality data on across life stage: hatch, juv_1, juv_2
dat_hatch_juv1_2 <- read.csv ("Kwild_code/data/Cort_hatch_juv1_2_OC.csv")

# 2) bring in data that has final morphology measurement and hormone data: "juv_3"
dat_juv3 <- read.csv(file = "Kwild_code/data/Cort_juv_3.csv") %>% 
  rename(sex = juv3_Sex)

# 3) final merge juv 3 data that contains hormonal results and final morph with dat_hatch_juv1_2 
merged_data <- merge(dat_hatch_juv1_2, dat_juv3, by.x = "Lizard_ID", 
                     by.y = "juv3_Lizard_ID", all = TRUE)


#########################
# Merge checks
#########################
# 1) find "new" values in Juv3: # LD745_21???
new_ids_juv3 <- dat_juv3$juv3_Lizard_ID[!dat_juv3$juv3_Lizard_ID %in% dat_hatch_juv1_2$Lizard_ID]

# 2) Get IDs from dat_hatch_juv1_2 that are not in dat_juv3
missing_ids_juv3 <- dat_hatch_juv1_2$Lizard_ID[!dat_hatch_juv1_2$Lizard_ID %in% dat_juv3$juv3_Lizard_ID]
# Retrieve all rows from dat_hatch_juv1_2 where Lizard_ID is in missing_ids
missing_ids_dat <- dat_hatch_juv1_2[dat_hatch_juv1_2$Lizard_ID %in% missing_ids_juv3, ]
missing_ids_dat <- missing_ids_dat %>% 
  filter(liz_status_db == "ALIVE")
table(missing_ids_dat$liz_status_db) # 15 Dead # 21 Alive since juv 2 sample


############
# Rename columns
data_final <- merged_data %>% 
  rename(hatch_svl_mm = hatch_svl_orig_mm,
         hatch_mass_g = hatch_mass_orig, 
         juv1_SVL_mm = juv1_measure_SVL_mm,
         juv2_SVL_mm = juv2_measure_SVL_mm,
         juv1_mass_g = juv1_measure_mass_g,
         juv2_mass_g = juv2_measure_mass_g,)


#############
##Calculating age in days post hatching for juvenile 1 and 2 measurements and final adult measurement
####Using Julian dates

##Juvenile 1 measurements all collected on 03/06/2022 = julian date 154
##Juvenile 2 measurements all collected on 10/08/2022 = julian date 222
data_final$Juvenile1_Age <- 154 - data_final$Julian_hatch_date
data_final$Juvenile2_Age <- 222 - data_final$Julian_hatch_date

##Juvenile 3 measurements collected at variable dates; measurement date converted to Julian date and age calculated as difference
##Note that dates span 2022 - 2023
str(data_final$juv3_date)
data_final$juv3_date <- as.Date(data_final$juv3_date, "%d/%m/%y")
data_final$Juvenile3_Age <- (yday(data_final$juv3_date)) + (365- data_final$Julian_hatch_date)


###Mean, SD, and range of age of measurements post hatching
##Juvenile 1 = 105.7, 10.8, 85 - 123
##Juvenile 2 = 173.7, 10.8, 153 - 191
##Juvenile 3 (or adult) = 466.1, 12.4, 440 - 491
##note: SD for Juvenile 1 and Juvenile 2 are the same because these measurements were all collected on the same day
mean(data_final$Juvenile1_Age, na.rm = TRUE)
sd(data_final$Juvenile1_Age, na.rm = TRUE)
range(data_final$Juvenile1_Age, na.rm = TRUE)

mean(data_final$Juvenile2_Age, na.rm = TRUE)
sd(data_final$Juvenile2_Age, na.rm = TRUE)
range(data_final$Juvenile2_Age, na.rm = TRUE)

mean(data_final$Juvenile3_Age, na.rm = TRUE)
sd(data_final$Juvenile3_Age, na.rm = TRUE)
range(data_final$Juvenile3_Age, na.rm = TRUE)

##Ondi growth rate calculations
########
# growth rate calculations: calculated growth rates by between each time point (e.g., hatching to juvenile 1; juvenile 1 to juvenile 2) and overall
# dividing change body size (SVL and mass) by the total number of days elapsed
data_final <- data_final %>% 
  mutate(hatch_juv3_SVL_growth =  (juv3_SVL_mm-hatch_svl_mm)/Juvenile3_Age,
         hatch_juv3_MASS_growth = (juv3_mass_g-hatch_mass_g)/Juvenile3_Age,
         hatch_juv2_SVL_growth = (juv2_SVL_mm-hatch_svl_mm)/Juvenile2_Age,
         hatch_juv2_MASS_growth = (juv2_mass_g-hatch_mass_g)/Juvenile2_Age,
         hatch_juv1_SVL_growth = (juv1_SVL_mm-hatch_svl_mm)/Juvenile1_Age,
         hatch_juv1_MASS_growth = (juv1_mass_g-hatch_mass_g)/Juvenile1_Age, 
         juv1_juv2_MASS_growth = (juv2_mass_g - juv1_mass_g)/ (Juvenile2_Age - Juvenile1_Age), 
         juv1_juv2_SVL_growth = (juv2_SVL_mm - juv1_SVL_mm)/ (Juvenile2_Age - Juvenile1_Age), 
         juv2_juv3_MASS_growth = (juv3_mass_g - juv2_mass_g)/ (Juvenile3_Age - Juvenile2_Age),
         juv2_juv3_SVL_growth = (juv3_SVL_mm - juv2_SVL_mm)/ (Juvenile3_Age - Juvenile2_Age))

###checking variables
###Some variables need to be changed to factors so they are not treated as linear variables in analyses
data_final$temp <- as.factor(data_final$temp)
data_final$hormone <- as.factor(data_final$hormone)
data_final$clutch <- as.factor(data_final$clutch)
data_final$Plate_CORT_juv3 <- as.factor(data_final$Plate_CORT_juv3)
data_final$juv3_T4_plate <- as.factor(data_final$juv3_T4_plate)
data_final$juv3_Testosterone_plate <- as.factor(data_final$juv3_Testosterone_plate)
data_final$juv3_HandlingTime_sec <- as.numeric(data_final$juv3_HandlingTime_sec)
data_final$juv3_liver_time_sec <- as.numeric(data_final$juv3_liver_time_sec)
###############
####Body condition calculations
###residuals and scaled mass index calculated and produce equivalent statistical results
#Scaled mass calculations calculated separately for each age (hatching, juvenile1, juvenile 2, and juvenile 3)
##Calculates residual condition as mass divided by SVL
data_final <- data_final %>% 
  mutate(hatch_mass_svl =  hatch_mass_g / hatch_svl_mm, 
          juv1_mass_svl =  juv1_mass_g /  juv1_SVL_mm,
         juv2_mass_svl = juv2_mass_g / juv2_SVL_mm,
         juv3_mass_svl = juv3_mass_g / juv3_SVL_mm,
         log_hatch_mass = log(hatch_mass_g),
         log_hatch_SVL = log(hatch_svl_mm), 
         log_juv1_mass = log(juv1_mass_g),
         log_juv1_SVL = log(juv1_SVL_mm),
         log_juv2_mass = log(juv2_mass_g), 
         log_juv2_SVL = log(juv2_SVL_mm),
         log_juv3_mass = log(juv3_mass_g), 
         log_juv3_SVL = log(juv3_SVL_mm))
  
###Scaled mass for hatchlings; creates new variable 'ScaledMass_hatch' 
##Pearsons correlation coefficient hatching
PCor_hatch <- cor(data_final$log_hatch_mass, data_final$log_hatch_SVL, method = "pearson", use = "complete.obs")
#OLS regression of ln transformed data
mod_hatch <- lm(log_hatch_mass ~ log_hatch_SVL, data = data_final)
##OLS coefficient from summary (1.4263)
summary(mod_hatch)
##Calculates mean svl for scaled mass equation
MeanSVL_hatch <- mean(data_final$hatch_svl_mm, na.rm = TRUE)
SMCof_hatch <- 1.4263/PCor_hatch
##calculates scaled mass; plots scaled mass to residual mass
data_final$ScaledMass_hatch <- data_final$hatch_mass_g * (MeanSVL_hatch /data_final$hatch_svl_mm) ^ (SMCof_hatch)
plot (ScaledMass_hatch ~ hatch_mass_svl, data = data_final)
ols_scaled_mass_mod <- lm(ScaledMass_hatch ~ hatch_mass_svl, data = data_final)
summary(ols_scaled_mass_mod)
abline(ols_scaled_mass_mod)


##Scaled mass for juvenile 1; creates new variable 'ScaledMass_juv1' 
###There is one lizard with a small svl (22), but the measurements for svl for this lizard at Juvenile 2 are also small (svl = 25)
range(data_final$juv1_SVL_mm, na.rm = TRUE)
##Pearsons correlation coefficient hatching
PCor_juv1 <- cor(data_final$log_juv1_mass, data_final$log_juv1_SVL, method = "pearson", use = "complete.obs")
mod_juv1 <- lm(log_juv1_mass ~ log_juv1_SVL, data = data_final)
##OLS coefficient from summary (2.2480)
summary(mod_juv1)
##Calculates mean tarsus for scaled mass equation
MeanSVL_juv1 <- mean(data_final$juv1_SVL_mm, na.rm = TRUE)
SMCof_juv1 <- 2.2480/PCor_juv1
##calculates scaled mass; plots scaled mass to residual mass
data_final$ScaledMass_juv1 <- data_final$juv1_mass_g * (MeanSVL_juv1 /data_final$juv1_SVL_mm) ^ (SMCof_juv1)
plot (ScaledMass_juv1 ~ juv1_mass_svl, data = data_final)
ols_scaled_mass_mod <- lm(ScaledMass_juv1  ~ juv1_mass_svl, data = data_final)
summary(ols_scaled_mass_mod)
abline(ols_scaled_mass_mod)


###Scaled mass for juvenile 2; creates new variable 'ScaledMass_juv2' 
##Pearsons correlation coefficient hatching
PCor_juv2 <- cor(data_final$log_juv2_mass, data_final$log_juv2_SVL, method = "pearson", use = "complete.obs")
mod_juv2 <- lm(log_juv2_mass ~ log_juv2_SVL, data = data_final)
##OLS coefficient from summary (2.3591)
summary(mod_juv2)
##Calculates mean tarsus for scaled mass equation
MeanSVL_juv2 <- mean(data_final$juv2_SVL_mm, na.rm = TRUE)
SMCof_juv2 <- 2.3591/PCor_juv2
##calculates scaled mass; plots scaled mass to residual mass
data_final$ScaledMass_juv2 <- data_final$juv2_mass_g * (MeanSVL_juv2 /data_final$juv2_SVL_mm) ^ (SMCof_juv2)
plot (ScaledMass_juv2 ~ juv2_mass_svl, data = data_final)
ols_scaled_mass_mod <- lm(ScaledMass_juv2  ~ juv2_mass_svl, data = data_final)
summary(ols_scaled_mass_mod)
abline(ols_scaled_mass_mod)

###Scaled mass for juvenile 3; creates new variable 'ScaledMass_juv3' 
##Pearsons correlation coefficient hatching
PCor_juv3 <- cor(data_final$log_juv3_mass, data_final$log_juv3_SVL, method = "pearson", use = "complete.obs")
mod_juv3 <- lm(log_juv3_mass ~ log_juv3_SVL, data = data_final)
##OLS coefficient from summary (2.4744)
summary(mod_juv3)
##Calculates mean tarsus for scaled mass equation
MeanSVL_juv3 <- mean(data_final$juv3_SVL_mm, na.rm = TRUE)
SMCof_juv3 <- 2.4744/PCor_juv3
##calculates scaled mass; plots scaled mass to residual mass
data_final$ScaledMass_juv3 <- data_final$juv3_mass_g * (MeanSVL_juv3 /data_final$juv3_SVL_mm) ^ (SMCof_juv3)
plot (ScaledMass_juv3 ~ juv3_mass_svl, data = data_final)
ols_scaled_mass_mod <- lm(ScaledMass_juv3  ~ juv3_mass_svl, data = data_final)
summary(ols_scaled_mass_mod)
abline(ols_scaled_mass_mod)

################
#### 1.	What are the effects of developmental treatments (temp, cort, interaction) on time to hatch?
#### effect of temperature on incubation days - faster hatch warmer temps 
################
# mod ; interaction removed
hatch_dev_mod <- lm(days_to_hatch ~ temp + hormone, data = data_final)
check_model(hatch_dev_mod) # bimodal because there is an effect on temp to days of hatch

# temperature effects days to hatch; no effect on hormone or temp x hormone interaction
###Temperature: F1,115 = 184.64, p<0.001
###hormone: F2,115 = 0.13, p=0.88
Anova(hatch_dev_mod)
summary(hatch_dev_mod)

# faster hatch warmer temps
hatch_dev_mod_emm <- emmeans(hatch_dev_mod, pairwise ~ temp)
plot(hatch_dev_mod_emm)
saveRDS(hatch_dev_mod, "Kwild_code/models/hatch_dev_mod.RDS")


###Mean days to hatch: 28 - 30.93; 23 - 48.31; SD: 28 - 4.82; 23 - 8.38
tapply(data_final$days_to_hatch, data_final$temp, mean, na.rm = TRUE)
tapply(data_final$days_to_hatch, data_final$temp, sd, na.rm = TRUE)


################
#### 2.	What are the effects of developmental treatments on body size, mass, and condition (scaled mass) of HATCHLINGS
#### SVL: hormone effects body size; no effects on temp or temp x hormone interaction
#### MASS: temperature and hormone effects mass; no interaction effects
#### Condition:  temperature has effect on BCI; no effect on hormone, or temp x hormone interaction

################
#### SVL; hormone*temp interaction removed
SVL_hatch_mod <- lm(hatch_svl_mm ~ temp + hormone , data = data_final)
check_model(SVL_hatch_mod)
Residuals <- residuals(SVL_hatch_mod)
shapiro.test(Residuals)

hist(data_final$hatch_svl_mm)

# hormone effects body size; high CORT has larger SVL than low CORT and control; no effects on temp or temp x hormone interaction
##hormone: F2,123 = 5.68, p=0.004
##temp: F1,123 = 0.07, p = 0.79
Anova(SVL_hatch_mod)
summary(SVL_hatch_mod)

# high hormone lower SVL
##High - control: p = 0.006
##High - low: p = 0.02
##Low - control: p = 0.92
SVL_hatch_mod_emm <- emmeans(SVL_hatch_mod, pairwise ~ hormone)
plot(SVL_hatch_mod_emm)
saveRDS(SVL_hatch_mod, "Kwild_code/models/SVL_hatch_mod.RDS")


#### MASS; hormone*temp interaction removed
Mass_hatch_mod <- lm(hatch_mass_g ~ temp + hormone, data = data_final)
check_model(Mass_hatch_mod)

# temperature and hormone effects mass; no interaction effects
##temp: F1,123 = 4.08, p = 0.046
##hormone: F1,123 = 3.11, p = 0.048

Anova(Mass_hatch_mod)
summary(Mass_hatch_mod)

# hormone pairwise comparison - high CORT lower mass compared to control
##High - control: p = 0.047
##High - low: p =0.77
##Low - control = 0.19
Mass_hormone_hatch_mod_emm <- emmeans(Mass_hatch_mod, pairwise ~ hormone)
plot(Mass_hormone_hatch_mod_emm)

# temp pairwise comparison - cooler temps higher mass
Mass_temp_hatch_mod_emm <- emmeans(Mass_hatch_mod, pairwise ~ temp)
plot(Mass_temp_hatch_mod_emm)
saveRDS(Mass_hatch_mod, "Kwild_code/models/MASS_hatch_mod.RDS")

#### Condition
# OLS residuals
condition_hatch_fit <- lm(hatch_mass_g ~ hatch_svl_mm, data = data_final, na.action=na.exclude) 
data_final$hatch_bc_residuals <- residuals(condition_hatch_fit, na.action=na.exclude) 
## hatchlings-anova on BCI residuals; interaction removed
condition_hatch_mod <- lm(hatch_bc_residuals ~ temp + hormone, data = data_final)
check_model(condition_hatch_mod)
# temperature has effect on BCI; no effect on hormone, or temp x hormone interaction
Anova(condition_hatch_mod)
summary(condition_hatch_mod)
# temp and BCI plot - poor BCI warm temps
condition_hatch_mod_emm <- emmeans(condition_hatch_mod, pairwise ~ temp)
plot(condition_hatch_mod_emm)
saveRDS(condition_hatch_mod, "Kwild_code/models/condition_hatch_mod.RDS")

####Scaled mass at hatching; results are consistent with OLS model of condition
scaled_mass_hatch_mod <- lm(ScaledMass_hatch ~ temp + hormone, data = data_final)
check_model(scaled_mass_hatch_mod)

# temperature has effect on BCI; no effect on hormone, or temp x hormone interaction
##temp: F1,123 = 4.10, p = 0.045; cold temp in higher condition
##hormone: F1,123 = 2.63, p = 0.08
Anova(scaled_mass_hatch_mod)
summary(scaled_mass_hatch_mod)

scaled_mass_hatch_mod_emm <- emmeans(scaled_mass_hatch_mod, pairwise ~ temp)
plot(scaled_mass_hatch_mod_emm )


################
#### 3.	What are the effects of developmental treatments on juvenile body size, condition, and mass (at the 3 points measured after hatching)?
#### Juv1_SVL: effect of temperature- cooler temps larger svl
#### Juv1_MASS: low mass with high hormones; low mass with high temperatures
#### Juv1_BCI: no treatment effects
#### Juv2_SVL: no effects - hormonal effect p = 0.07
#### Juv2_MASS: low mass with high temperatures
#### Juv2_BCI: no treatment effects
#### Juv3_MASS: low mass males
#### Juv3_SVL: high hormone low body size; males smaller than females
#### Juv3_BCI: Sex effect: male lower BCI than females - could be size difference

###############
# SVL ANALYSIS: Juv_1, Juv_2, Juv_3 (sex accounted for last measurement)
###############
#### JUV_1: SVL - adding days since hatch because effect of temp on development time); interaction removed
SVL_Juv1_mod <- lm(juv1_SVL_mm ~ temp + hormone + scale(Juvenile1_Age), data = data_final)
check_model(SVL_Juv1_mod)
# effect on temperature and day of hatch on body size 
Anova(SVL_Juv1_mod)
summary(SVL_Juv1_mod)
# plot - cooler temps larger svl
SVL_Juv1_mod_emm <- emmeans(SVL_Juv1_mod, pairwise ~ temp)
plot(SVL_Juv1_mod_emm)
saveRDS(SVL_Juv1_mod, "Kwild_code/models/SVL_Juv1_mod.RDS")

#### JUV_2: SVL; interaction removed
SVL_Juv2_mod <- lm(juv2_SVL_mm ~ temp + hormone + scale(Juvenile2_Age), data = data_final)
check_model(SVL_Juv2_mod)
# marginal hermonal effect on body size p = 0.06
summary(SVL_Juv2_mod)
Anova(SVL_Juv2_mod)
SVL_Juv2_mod_emm <- emmeans(SVL_Juv2_mod, pairwise ~ temp)
plot(SVL_Juv2_mod_emm)
saveRDS(SVL_Juv2_mod, "Kwild_code/models/SVL_Juv2_mod.RDS")

#### SVL - NOTE SEX Accounted; no interaction
SVL_Juv3_mod <- lm(juv3_SVL_mm ~ temp + hormone + sex + scale(Juvenile3_Age), data = data_final)
check_model(SVL_Juv3_mod)
# effect on hormones and sex
Anova(SVL_Juv3_mod)
summary(SVL_Juv3_mod)

# plot - high hormone lower svl
SVL_Juv3_hormone_mod_emm <- emmeans(SVL_Juv3_mod, pairwise ~ hormone)
plot(SVL_Juv3_hormone_mod_emm)

# plot - males are smaller than females
SVL_Juv3_sex_mod_emm <- emmeans(SVL_Juv3_mod, pairwise ~ sex)
plot(SVL_Juv3_sex_mod_emm)
saveRDS(SVL_Juv3_mod, "Kwild_code/models/SVL_Juv3_mod.RDS")


###############
# MASS ANALYSIS: Juv_1, Juv_2, Juv_3
###############
#### MASS: Juv_1 - interaction removed
Mass_Juv1_mod <- lm(juv1_mass_g ~ temp + hormone + scale(Juvenile1_Age), data = data_final)
check_model(Mass_Juv1_mod)
# effect on temp, hormone, and days since hatch; no effects on temp x hormone interaction
summary(Mass_Juv1_mod)
Anova(Mass_Juv1_mod)
# hormone plot - low mass with high hormones
Mass_Juv1_hormone_mod_emm <- emmeans(Mass_Juv1_mod, pairwise ~ hormone)
plot(Mass_Juv1_hormone_mod_emm)
# temp plot - low mass with high temperatures
Mass_Juv1_temp_mod_emm <- emmeans(Mass_Juv1_mod, pairwise ~ temp)
plot(Mass_Juv1_temp_mod_emm)
saveRDS(Mass_Juv1_mod, "Kwild_code/models/Mass_Juv1_mod.RDS")

#### MASS: Juv_2; interaction removed
Mass_Juv2_mod <- lm(juv2_mass_g ~ temp + hormone + scale(Juvenile2_Age), data = data_final)
check_model(Mass_Juv2_mod)
# temp  effect
summary(Mass_Juv2_mod)
Anova(Mass_Juv2_mod)

# temp plot - low mass with high temperatures
Mass_Juv2_temp_mod_emm <- emmeans(Mass_Juv2_mod, pairwise ~ temp)
plot(Mass_Juv2_temp_mod_emm)
saveRDS(Mass_Juv2_mod, "Kwild_code/models/Mass_Juv2_mod.RDS")

#### MASS: Juv_3; NOTE SEX IS ACCOUNTED FOR
Mass_Juv3_mod <- lm(juv3_mass_g ~ hormone + temp + sex + scale(Juvenile3_Age), data = data_final)
check_model(Mass_Juv3_mod)
# sex effect
Anova(Mass_Juv3_mod)
summary(Mass_Juv3_mod)
# temp plot - males smaller than females
Sex_Juv3_temp_mod_emm <- emmeans(Mass_Juv3_mod, pairwise ~ hormone)
plot(Sex_Juv3_temp_mod_emm)
saveRDS(Mass_Juv3_mod, "Kwild_code/models/Mass_Juv3_mod.RDS")

###############
# Scaled mass ANALYSIS: Juv_1, Juv_2, Juv_3
###############
view(data_final)

## Juvenile 1-anova on scaled mass; interaction removed
condition_juv_mod <- lm(ScaledMass_juv1 ~ temp + hormone  + scale(Juvenile1_Age), data = data_final)
# no treatment effects on scaled mass juvenile 1
Anova(condition_juv_mod)
summary(condition_juv_mod)
check_model(condition_juv_mod)


### Juvenile 2 - scaled mass; temp*hormone interaction removed
condition_juv2_mod <- lm(ScaledMass_juv2 ~ temp + hormone +scale(Juvenile2_Age), data = data_final)
check_model(condition_juv2_mod)
# no treatment effects on bci juv1
Anova(condition_juv2_mod)
summary(condition_juv2_mod)
saveRDS(condition_juv2_mod, "Kwild_code/models/condition_juv2_mod.RDS")

### Juvenile 3 scaled mass, sex added; no treatment affects but males have lower condition
condition_juv3_mod <- lm(ScaledMass_juv3 ~ temp + hormone + sex + scale(Juvenile3_Age), data = data_final)
check_model(condition_juv3_mod)
# no treatment effects on bci juv1
Anova(condition_juv3_mod)
summary(condition_juv3_mod)
# sex plot - 
sex_condition_juv3_mod_emm <- emmeans(condition_juv3_mod, pairwise ~ sex)
plot(sex_condition_juv3_mod_emm)
saveRDS(condition_juv3_mod, "Kwild_code/models/condition_juv3_mod.RDS")


#########Kris code
# save final data with BCI residuals for figures and quarto doc
data_final <-  mutate(data_final, hormone = factor(hormone, 
                                                   levels = c("control","low","high"))) %>% 
  filter(!is.na(hormone)) %>%
  group_by(hormone)
write.csv(data_final, "Kwild_code/data/final_analysis_data_quarto.csv")
data_final <- read.csv(file = "Kwild_code/data/final_analysis_data_quarto.csv")

################
#### 4.	developmental treatment; Testing days to mortality across treatments
#### no differences between treatments. Mortality of duration of study
################
# survival fisher test- filter out missing animals
survival<- data_final %>% 
  filter(liz_status_db != "MISSING")
survival_table <- table(survival$liz_status_db, survival$hormone)
# analysis
survival_day_hormone <- kruskal.test(liz_status_db ~ hormone, data = survival)
survival_day_temp <- kruskal.test(liz_status_db ~ temp, data = survival)
# save analysis
saveRDS(survival_day_hormone, "Kwild_code/models/survival_day_hormone.RDS")
saveRDS(survival_day_temp, "Kwild_code/models/survival_day_temp.RDS")


################
#### 5.	Growthrate
#### CORT effects on mass for hatch and juvenile 1
#### Temperature effects on SVL for juvenile 1 and 2
#### Temperature effects on mass for hatch juvenile 1 and 2 
################
# SVL growth: hatch to juv_1; interaction removed
juvenile1_svl_growth_mod  <- lm(hatch_juv1_SVL_growth ~ hormone + temp + sex, 
                      data = data_final, na.action=na.exclude)
check_model(juvenile1_svl_growth_mod)
# treatment effect on temp but not hormone
Anova(juvenile1_svl_growth_mod)
summary(juvenile1_svl_growth_mod)
# temp and growth plot - faster growth rates in cooler temps
juvenile1_svl_growth_mod_emm <- emmeans(juvenile1_svl_growth_mod, pairwise ~ temp)
plot(juvenile1_svl_growth_mod_emm )
saveRDS(juvenile1_svl_growth_mod, "Kwild_code/models/overall_SVL_growth_mod.RDS")

# Mass growth - look at change in mass from hatching to juvenile 1; interaction removed
juvenile1_mass_growth_mod <- lm( hatch_juv1_MASS_growth ~ hormone + temp, 
                       data = data_final, na.action=na.exclude) 
check_model(juvenile1_mass_growth_mod)
# treatment effect on temp but not hormone
Anova(juvenile1_mass_growth_mod)
summary(juvenile1_mass_growth_mod)

#Temp no effect on change in mass between hatching and juvenile 1; CORT treatment affects growth; high CORT has slower growth than control
juvenile1_mass_growth_mod_emm <- emmeans(juvenile1_mass_growth_mod, pairwise ~ hormone)
plot(juvenile1_mass_growth_mod_emm)
saveRDS(juvenile1_mass_growth_mod, "Kwild_code/models/growth4_mass_mod.RDS")

##Growth rate for juvenile 1 to juvenile 2 
juvenile2_svl_growth_mod <- lm(juv1_juv2_SVL_growth ~hormone + temp, 
                          data = data_final, na.action=na.exclude)
Anova(juvenile2_svl_growth_mod)
summary(juvenile2_svl_growth_mod)

juvenile2_svl_growth_mod_emm <- emmeans(juvenile2_svl_growth_mod, pairwise ~ temp)
plot(juvenile2_svl_growth_mod_emm)

juvenile2_mass_growth_mod <- lm(juv1_juv2_MASS_growth ~ hormone + temp, 
                          data = data_final, na.action=na.exclude)
Anova(juvenile2_mass_growth_mod )
summary(juvenile2_mass_growth_mod )
juvenile2_mass_growth_mod_emm <- emmeans(juvenile2_mass_growth_mod, pairwise ~ temp)
plot(juvenile2_mass_growth_mod_emm )


###Growth rate from juvenile 2 to juvenile 3
juvenile3_svl_growth_mod <- lm(juv2_juv3_SVL_growth ~hormone + temp, 
                               data = data_final, na.action=na.exclude)
Anova(juvenile3_svl_growth_mod)
summary(juvenile3_svl_growth_mod)

juvenile3_svl_growth_mod_emm <- emmeans(juvenile3_svl_growth_mod, pairwise ~ temp)
plot(juvenile3_svl_growth_mod_emm )

juvenile3_mass_growth_mod <- lm(juv2_juv3_MASS_growth ~ hormone + temp, 
                                data = data_final, na.action=na.exclude)
Anova(juvenile3_mass_growth_mod)
summary(juvenile3_mass_growth_mod )
juvenile3_mass_growth_mod_emm <- emmeans(juvenile3_mass_growth_mod, pairwise ~ temp)
plot(juvenile3_mass_growth_mod_emm  )

###Overall growth hatch to juvenile 3; sex included
hatch_juv3_SVL_growth_mod <- lm(hatch_juv3_SVL_growth ~hormone + temp + sex, 
                               data = data_final, na.action=na.exclude)
Anova(hatch_juv3_SVL_growth_mod)
summary(hatch_juv3_SVL_growth_mod)

hatch_juv3_SVL_growth_mod_emm <- emmeans(hatch_juv3_SVL_growth_mod, pairwise ~ temp)
plot(hatch_juv3_SVL_growth_mod_emm  )

hatch_juv3_MASS_growth_mod <- lm(hatch_juv3_MASS_growth ~hormone + temp, 
                                data = data_final, na.action=na.exclude)
Anova(hatch_juv3_MASS_growth_mod)
summary(hatch_juv3_MASS_growth_mod)


########################
### 6: CORT, T4, testosterone
###OC changed >= to > to remove values of 4ul and below
########################
# data
cort_dat <- data_final %>%
  drop_na(juv3_CORT_Final_Hormone_ng_mL, temp, hormone, juv3_HandlingTime_sec) %>% 
  filter(juv3_Sample_volume_ul > 4 ) # remove values below 4 ul juv3_Sample_volume_ul

t4_dat <- data_final %>%
  drop_na(juv3_T4_corrected_ng_mL, temp, hormone, Lizard_ID) %>% 
  filter(juv3_Sample_volume_ul > 4 ) # remove values below 4 ul juv3_Sample_volume_ul

test_dat <- data_final %>%
  drop_na(juv3_Testosterone_Final_ng_ml, temp, hormone, Lizard_ID) %>% 
  filter(juv3_Sample_volume_ul > 4 ) # remove values below 4 ul juv3_Sample_volume_ul

plot(juv3_CORT_Final_Hormone_ng_mL ~ hormone, data = cort_dat)
plot(juv3_CORT_Final_Hormone_ng_mL ~ temp, data = cort_dat)

########
# CORT
######## 1) treatments- how does cort vary across temp and hormone treatment
###when cort plate is changed to nominal variable; hormone treatment effect is P=0.06

cort_development_mod <- lm(log(juv3_CORT_Final_Hormone_ng_mL) ~ hormone + temp + juv3_HandlingTime_sec + Plate_CORT_juv3 + sex, data = cort_dat)
check_model(cort_development_mod)
Anova(cort_development_mod)
summary(cort_development_mod)

hist(cort_dat$juv3_CORT_Final_Hormone_ng_mL)
Residuals <- residuals(cort_development_mod)
shapiro.test(Residuals)

###using gamma distribution; results qualitatively the same
cort_development_gam <- glm (juv3_CORT_Final_Hormone_ng_mL ~ hormone + temp + juv3_HandlingTime_sec + Plate_CORT_juv3 + sex, 
                             family = Gamma(link = "log"), data = cort_dat)

check_model(cort_development_gam)
Anova(cort_development_gam, type = "II", test = "F")
summary(cort_development_gam)
saveRDS(cort_development_mod, "Kwild_code/models/cort_development_mod.RDS")

###Save residuals from model that accounts for plate effects and handling time; handling time not significant so removed
cort_factors <- lm(log(juv3_CORT_Final_Hormone_ng_mL) ~ Plate_CORT_juv3, data = cort_dat)
check_model(cort_factors )
Anova(cort_factors )
summary(cort_factors )
residuals_CORT_factors <- residuals(cort_factors)

####Body size and condition at juvenile 3 related to CORT levels; 'handling time not included because it did not affect baseline CORT levels as above)
##CORT and scaled mass at juvenile 3; treatments not included because they have no effect on body condition at this time point
cort_body_condition_mod <- lm(ScaledMass_juv3 ~ (log(juv3_CORT_Final_Hormone_ng_mL)) + Plate_CORT_juv3 + sex + scale(Juvenile3_Age), data = cort_dat)
Anova(cort_body_condition_mod)
summary(cort_body_condition_mod)

###SVL with CORT; adult SVL positive assoication with CORT
cort_svl_juvenile3_mod <- lm(juv3_SVL_mm ~ (log(juv3_CORT_Final_Hormone_ng_mL)) +  Plate_CORT_juv3 + sex + scale(Juvenile3_Age), data = cort_dat)
Anova(cort_svl_juvenile3_mod)
summary(cort_svl_juvenile3_mod)

###Mass at juvenile 3 tested against CORT levels; CORT positive association with adult mass
cort_mass_juvenile3_mod <- lm(juv3_mass_g ~ (log(juv3_CORT_Final_Hormone_ng_mL)) +  Plate_CORT_juv3 + sex + scale(Juvenile3_Age), data = cort_dat)
Anova(cort_mass_juvenile3_mod)
summary(cort_mass_juvenile3_mod)


################
# T4
################ 
# 1) treatments -  does T4 vary across temp and hormone treatment: NS; plate and handling time are near significant
T4_temp_hormone_sex_mod <- lm(log(juv3_T4_corrected_ng_mL) ~hormone + temp + sex + juv3_T4_plate + juv3_HandlingTime_sec, data = t4_dat)
check_model(T4_temp_hormone_sex_mod)
Anova(T4_temp_hormone_sex_mod)
summary(T4_temp_hormone_sex_mod)
# sex differences: females have higher T4
T4_temp_hormone_sex_mod_emm <- emmeans(T4_temp_hormone_sex_mod, pairwise ~ sex)
plot(T4_temp_hormone_sex_mod_emm)
saveRDS(T4_temp_hormone_sex_mod, "Kwild_code/models/T4_temp_hormone_sex_mod.RDS")


# 2) SVL and associations with T4 and cort
T4_SVL_mod <- lm(juv3_SVL_mm  ~ ((log(juv3_T4_corrected_ng_mL))) + sex + scale(Juvenile3_Age), 
             data = t4_dat)
check_model(T4_SVL_mod )
summary(T4_SVL_mod )
Anova(T4_SVL_mod )

# 3) MASS and associations between T4
T4_mass_growth_mod <- lm(juv3_mass_g  ~ (log(juv3_T4_corrected_ng_mL)) + sex + scale(Juvenile3_Age), 
                        data = t4_dat)
check_model(T4_mass_growth_mod)
Anova(T4_mass_growth_mod)
summary(T4_mass_growth_mod)
saveRDS(T4_mass_growth_mod, "Kwild_code/models/T4_mass_growth_mod.RDS")


########################
###Testosterone

###Look at the effects of plate and handling time on testosterone levels; handling time affects T levels (p=0.02); plate is NS
Test_factor <- lm (juv3_Testosterone_Final_ng_ml ~ juv3_Testosterone_plate+ juv3_HandlingTime_sec, data = test_dat)
Anova(Test_factor)
summary(Test_factor)


###Does testosterone vary based on treatment; treats T levels as linear
Testosterone_mod <- lm((log(juv3_Testosterone_Final_ng_ml)) ~temp + hormone + juv3_HandlingTime_sec, data = test_dat)
check_model(Testosterone_mod)
Anova(Testosterone_mod)
summary(Testosterone_mod)


###Does mass or body size depend on testosterone? 
test_mass_mod <- lm(juv3_mass_g  ~ (log(juv3_Testosterone_Final_ng_ml)) + juv3_HandlingTime_sec + scale(Juvenile3_Age), 
                         data = test_dat)
Anova(test_mass_mod)
summary(test_mass_mod)


test_svl_mod <- lm(juv3_SVL_mm  ~ (log(juv3_Testosterone_Final_ng_ml)) + juv3_HandlingTime_sec  + scale(Juvenile3_Age), 
                    data = test_dat)
Anova(test_svl_mod)
summary(test_svl_mod)


test_condition_mod <- lm(ScaledMass_juv3  ~ (log(juv3_Testosterone_Final_ng_ml)) + juv3_HandlingTime_sec  + scale(Juvenile3_Age), 
                   data = test_dat)
Anova(test_condition_mod)
summary(test_condition_mod)


########################
### 7: Mito function
########################
# bring in data, rename, and save
###if looking at relationships with hormone levels - use cort_dat or equivalent to exclude samples with 
##low volumes
mito_dat <- data_final %>%
  dplyr::select(c("Lizard_ID", "clutch", "temp", "hormone", "Juvenile3_Age", "sex","hatch_juv3_MASS_growth",
                  "juv3_inject_time_sec", "juv3_liver_time_sec", "juv3_oroboros",
                  "juv3_T4_plate", "juv3_HandlingTime_sec", "juv3_T4_corrected_ng_mL", 
                  "juv3_CORT_Final_Hormone_ng_mL", "juv3_basal_corrected.pmol..sec.ng..",
                  "juv3_chamber","juv3_basal_corrected.pmol..sec.ng..", 
                  "juv3_adp_corrected.pmol..sec.ng..", "juv3_oligo_corrected.pmol..sec.ng..", 
                  "juv3_fccp_corrected.pmol..sec.ng..", "juv3_RCR.L.R.", "juv3_RCR.R.ETS.",
                  "juv3_oroboros_comments", "juv3_RCR.L.ETS.", "juv3_date", "Plate_CORT_juv3",
                  "juv3_Sample_volume_ul",
                  "juv3_mass_g", "juv3_SVL_mm", "juv2_juv3_MASS_growth", "ScaledMass_juv3",
                  "juv3_HandlingTime_sec", "juv2_juv3_SVL_growth", "juv3_chamber")) %>% 
  dplyr::rename(inject_time_sec = juv3_inject_time_sec, 
                chamber = juv3_chamber,
                juv3_basal_corrected = juv3_basal_corrected.pmol..sec.ng..,
                basal_corrected_pmol = juv3_basal_corrected.pmol..sec.ng..,
                oligo_corrected_pmol = juv3_oligo_corrected.pmol..sec.ng..,
                fccp_corrected_pmol = juv3_fccp_corrected.pmol..sec.ng..,
                adp_corrected_pmol =  juv3_adp_corrected.pmol..sec.ng..,
                T4_plate_ID = juv3_T4_plate,
                HandlingTime_sec = juv3_HandlingTime_sec,
                T4_corrected_ng_mL = juv3_T4_corrected_ng_mL, 
                CORT_Final_Hormone_ng_mL = juv3_CORT_Final_Hormone_ng_mL,
                RCR_L_R = juv3_RCR.L.R., 
                RCR_R_ETS = juv3_RCR.R.ETS.,
                RCR_L_ETS = juv3_RCR.L.ETS., 
                oroboros_comments = juv3_oroboros_comments) %>% 
  mutate(ID_and_Comments = paste(Lizard_ID, oroboros_comments, sep = ": "),
         RCR = adp_corrected_pmol/oligo_corrected_pmol) %>% 
  filter(Lizard_ID != c("LD736_21")) %>%
  filter(Lizard_ID != c("LD738_21")) %>%
  filter(Lizard_ID != c("LD829_21")) %>%
  filter(Lizard_ID != c("LD784_21"))

####LD736_21 and LD738_21 were first samples assayed; values are extremely high; filtered from analysis 

view(mito_dat)

###Checking factors that could affect oxygen consumption measurements - handling times, Oroboros, and chamber
###'chamber' is nested within Oroboros; there are three measurements of time (time to inject the drug, time to euthanize the animal
###'and time to collect liver, have opted to only use time to collect liver because it will encorporate all other time measurements and 
###'likely all time measures are correlated and so shouldn't all be going in the model anyway)
###all factors NS

#RCR - all NS
RCR_factors <- lm(RCR ~ juv3_liver_time_sec + juv3_oroboros/chamber, data = mito_dat)
Anova(RCR_factors)
summary(RCR_factors)
hist(mito_dat$RCR)

#RCR(L/ETS) - all NS
RCR_L_ETS_factors <- lm(RCR_L_ETS ~ juv3_liver_time_sec + juv3_oroboros/chamber, data = mito_dat)
Anova(RCR_L_ETS_factors)
summary(RCR_L_ETS_factors)
hist(mito_dat$RCR_L_ETS)

#RCR(R/ETS) - all NS
RCR_R_ETS_factors <- lm(RCR_R_ETS ~ juv3_liver_time_sec + juv3_oroboros/chamber, data = mito_dat)
Anova(RCR_R_ETS_factors)
summary(RCR_L_ETS_factors)
hist(mito_dat$RCR_R_ETS)

##RCR(L/R) - All factors NS
RCR_L_R_factors <- lm(RCR_L_R ~ juv3_liver_time_sec + juv3_oroboros/chamber, data = mito_dat)
Anova(RCR_L_R_factors)
summary(RCR_L_R_factors)
###There looks to be a strange value (IDLD829_21 that has really low basal- excluded for now but need to check raw oxygen files)
hist(mito_dat$RCR_L_R)
view(mito_dat)

###Basal; all factors NS
basal_factors <- lm(basal_corrected_pmol ~ juv3_liver_time_sec + juv3_oroboros, data = mito_dat)
Anova(basal_factors)
summary (basal_factors)
hist(mito_dat$basal_corrected_pmol)

##ADP - all factors NS
ADP_factors <- lm(adp_corrected_pmol ~ juv3_liver_time_sec + juv3_oroboros/chamber, data = mito_dat)
Anova(ADP_factors)
summary (ADP_factors)
#ADP - overdispersed?
hist(mito_dat$adp_corrected_pmol)

##Oligo - all factors NS but time to collect liver p = 0.06
Oligo_factors <- lm (oligo_corrected_pmol ~ juv3_liver_time_sec + juv3_oroboros/chamber, data = mito_dat)
Anova(Oligo_factors)
summary (Oligo_factors)
hist(mito_dat$oligo_corrected_pmol)

###FCCP - all factors = NS
FCCP_factors <- lm(fccp_corrected_pmol ~ juv3_liver_time_sec + juv3_oroboros/chamber, data = mito_dat) 
Anova(FCCP_factors)
summary (FCCP_factors)
hist(mito_dat$fccp_corrected_pmol)

#####Effects of treatments on mitochondrial function 
###Basal - all NS; interaction removed
Basal_treatment_mod <- lm(log(basal_corrected_pmol) ~ temp + hormone + sex, data = mito_dat)
Anova(Basal_treatment_mod)
summary(Basal_treatment_mod)

check_model(Basal_treatment_mod)
Mod_residuals <- residuals(Basal_treatment_mod, type = ("pearson"))
shapiro.test(Mod_residuals)

###checking gamma distribution; looks worse
basal_gamma_mod <- glm (basal_corrected_pmol ~ hormone + temp + sex, 
     family = Gamma(link = "log"), data = mito_dat)
Anova(basal_gamma_mod)
summary(basal_gamma_mod)
check_model(basal_gamma_mod)

###ADP - treatment effects; all NS, interaction removed
ADP_treatment_mod <- lm((log(adp_corrected_pmol)) ~ temp + hormone + sex, data = mito_dat)
Anova(ADP_treatment_mod)
summary(ADP_treatment_mod)

check_model(ADP_treatment_mod)
Mod_residuals <- residuals(ADP_treatment_mod, type = ("pearson"))
shapiro.test(Mod_residuals)

##Oligo mito by treatment; all factors NS
Oligo_treatment_mod <- lm((log(oligo_corrected_pmol)) ~ temp + hormone + sex, data = mito_dat)
Anova(Oligo_treatment_mod)
summary(Oligo_treatment_mod)

check_model(Oligo_treatment_mod)
Mod_residuals <- residuals(Oligo_treatment_mod, type = ("pearson"))
shapiro.test(Mod_residuals)

##FCCP mito by treatment; all factors NS 
FCCP_treatment_mod <- lm((log(fccp_corrected_pmol)) ~ temp + hormone + sex, data = mito_dat)
Anova(FCCP_treatment_mod)
summary(FCCP_treatment_mod)
check_model(FCCP_treatment_mod)

##RCR calculated as state 3 (ADP)/ state 4(olgio)
RCR_mod <- lm((log(RCR)) ~ temp + hormone + sex, data = mito_dat)
Anova(RCR_mod)
summary (RCR_mod)
check_model(RCR_mod)

# RCR (L/ETS): RCR_L_ETS effects by developmental environment and sex 
# interaction removed; 
RCR_L_ETS_treatment_mod <- lm((log(RCR_L_ETS)) ~ temp + hormone + sex, data = mito_dat)
Anova(RCR_L_ETS_treatment_mod)
summary(RCR_L_ETS_treatment_mod)
check_model(RCR_L_ETS_treatment_mod)

# RCR (R/ETS): RCR_R_ETS effects by developmental environment and sex 
# males have higher RCR(R/ETS); p = 0.03
RCR_R_ETS_sex_temp_hormone_mod <- lm((log(RCR_R_ETS)) ~ temp + hormone + sex, data = mito_dat)
Anova(RCR_R_ETS_sex_temp_hormone_mod)
summary(RCR_R_ETS_sex_temp_hormone_mod)
check_model(RCR_R_ETS_sex_temp_hormone_mod)

# RCR (L/R): treatment effects - all factors NS
RCR_L_R_treatment_mod <- lm((log(RCR_L_R)) ~ temp + hormone + sex, data = mito_dat)
Anova(RCR_L_R_treatment_mod)
summary(RCR_L_R_treatment_mod)
check_model(RRCR_L_R_treatment_mod)

###For models with hormone levels - need to filter n=7 data points from 
###samples with less than 7ul 
mito_hormone <- mito_dat %>%
  drop_na(CORT_Final_Hormone_ng_mL) %>% 
  filter(juv3_Sample_volume_ul > 4 ) # remove values below 4 ul juv3_Sample_volume_ul

###Associations between mitochondrial bioenergetics and hormones levels
view(mito_hormone)

##Basal - all factors NS
Basal_hormone_mod <- lm((log(basal_corrected_pmol)) ~ (log(CORT_Final_Hormone_ng_mL)) + (log(T4_corrected_ng_mL)) + sex + HandlingTime_sec, data = mito_hormone)
Anova(Basal_hormone_mod)
summary(Basal_hormone_mod)
check_model(Basal_hormone_mod)

###ADP - T4 positive association with ADP mito
ADP_hormone_mod <- lm((log(adp_corrected_pmol)) ~ (log(CORT_Final_Hormone_ng_mL)) + (log(T4_corrected_ng_mL)) + sex + HandlingTime_sec, data = mito_hormone)
Anova(ADP_hormone_mod)
summary(ADP_hormone_mod)
check_model(ADP_hormone_mod)

##oligo - all NS
Oligo_hormone_mod <- lm((log(oligo_corrected_pmol))  ~ (log(CORT_Final_Hormone_ng_mL)) + (log(T4_corrected_ng_mL)) + sex + HandlingTime_sec, data = mito_hormone)
Anova(Oligo_hormone_mod)
summary(Oligo_hormone_mod)
check_model(Oligo_hormone_mod)
Mod_residuals <- residuals(Oligo_treatment_mod, type = ("pearson"))
shapiro.test(Mod_residuals)

##FCCP - positive association with T4 (p = 0.04)
FCCP_hormone_mod <- lm((log(fccp_corrected_pmol)) ~ (log(CORT_Final_Hormone_ng_mL)) + (log(T4_corrected_ng_mL)) + sex + HandlingTime_sec, data = mito_hormone)
Anova(FCCP_hormone_mod)
summary(FCCP_hormone_mod)

RCR_hormone_mod <- lm((log(RCR)) ~  (log(CORT_Final_Hormone_ng_mL)) + (log(T4_corrected_ng_mL)) + sex + HandlingTime_sec, data = mito_hormone)
Anova(RCR_mod_cont)
summary (RCR_mod_cont)
Mod_residuals <- residuals(RCR_mod_cont, type = ("pearson"))
shapiro.test(Mod_residuals)

## RCR (L/ETS) - T4 negative association with RCR(L/ETS) - p =0.01
RCR_L_ETS_CORT_mod <- lm((log(RCR_L_ETS)) ~ (log(CORT_Final_Hormone_ng_mL)) + (log(T4_corrected_ng_mL)) + sex + HandlingTime_sec, data = mito_hormone)
Anova(RCR_L_ETS_CORT_mod)
summary(RCR_L_ETS_CORT_mod)
check_model(RCR_L_ETS_CORT_mod)
Mod_residuals <- residuals(RCR_L_ETS_CORT_mod, type = ("pearson"))
shapiro.test(Mod_residuals)

### RCR (L/R) - NS
RCR_L_R_hormone <- lm((log(RCR_L_R)) ~ (log(CORT_Final_Hormone_ng_mL)) + (log(T4_corrected_ng_mL)) + sex + HandlingTime_sec, data = mito_hormone)
Anova(RCR_L_R_hormone)
summary(RCR_L_R_hormone)
check_model(RCR_L_R_hormone)
Mod_residuals <- residuals(RCR_L_R_hormone, type = ("pearson"))
shapiro.test(Mod_residuals)

##RCR (R/ETS) - T4 negative assocation with RCR(R/ETS) p = 0.008
RCR_R_ETS_CORT_mod <- lm((log(RCR_R_ETS)) ~ (log(CORT_Final_Hormone_ng_mL)) + (log(T4_corrected_ng_mL)) + HandlingTime_sec + sex, data = mito_hormone)
Anova(RCR_R_ETS_CORT_mod)
summary(RCR_R_ETS_CORT_mod)
check_model(RCR_R_ETS_CORT_mod)

###Growth and body size and mito function 
###Basal - postive association with basal and growth (p = 0.02)
Growth_Basal <- lm(hatch_juv3_MASS_growth ~ (log(basal_corrected_pmol)) + sex, data = mito_dat)
Anova(Growth_Basal)
summary(Growth_Basal)

##ADP - positive association with ADP and growth (p = 0.04)
Growth_ADP <- lm(hatch_juv3_MASS_growth ~ (log(adp_corrected_pmol)) + sex, data = mito_dat)
Anova(Growth_ADP)
summary(Growth_ADP)

###Oligo - NS (p=0.07)
Growth_Oligo <- lm(hatch_juv3_MASS_growth ~ (log(oligo_corrected_pmol)) + sex, data = mito_dat)
Anova(Growth_Oligo)
summary(Growth_Oligo)

###FCCP NS p = 0.06
Growth_FCCP <- lm(hatch_juv3_MASS_growth ~ (log(fccp_corrected_pmol)) + sex, data = mito_dat)
Anova(Growth_FCCP)
summary(Growth_FCCP)

###NS for all RCR but sex sig for all with males growing slower
Growth_RCR <- lm(hatch_juv3_MASS_growth ~ (log(RCR)) + sex, data = mito_dat)
Anova(Growth_RCR)
summary(Growth_RCR)

Growth_RCR_L_R <- lm(hatch_juv3_MASS_growth ~ (log(RCR_L_R)) + sex, data = mito_dat)
Anova(Growth_RCR_L_R)
summary(Growth_RCR_L_R)

Growth_RCR_R_ETS <- lm(hatch_juv3_MASS_growth ~ (log(RCR_R_ETS)) + sex, data = mito_dat)
Anova(Growth_RCR_R_ETS)
summary(Growth_RCR_R_ETS)

Growth_RCR_L_ETS <- lm(hatch_juv3_MASS_growth ~ (log(RCR_L_ETS)) + sex, data = mito_dat)
Anova(Growth_RCR_L_ETS)
summary(Growth_RCR_L_ETS)

############# ############# ############# ############# ############# 
############# OLD CODE############# ############# ############# ############# 
############# ############# ############# ############# ############# ############# 
# Violin Plot of RCR_L_ETS by Hormone treatment
mito_dat <- mutate(mito_dat, hormone = factor(hormone, levels = c("control","low","high"))) %>%
  group_by(hormone)

# Function to calculate mean and standard error
mean_se <- function(x) {
  return(data.frame(y = mean(x), 
                    ymin = mean(x) - sd(x)/sqrt(length(x)), 
                    ymax = mean(x) + sd(x)/sqrt(length(x))))
}

# Creating the violin plot with data points overlaid, and means and standard errors
plot <- ggplot(data = mito_dat, aes(x = hormone, y = RCR_L_ETS)) +
  geom_violin(scale = "width", adjust = 1.5) +
  geom_jitter(aes(text = ID_and_Comments), width = 0.2, color = "blue") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "red") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "red") +
  labs(
    title = "Violin Plot of RCR(L/ETS) by Hormone Treatment",
    x = "Hormone Treatment",
    y = "RCR(L/ETS)") 

# Converting the ggplot object to a plotly object to enable interactive features
plot <- ggplotly(plot, tooltip = "text")

# Displaying the plot
plot



# Creating the violin plot with data points overlaid, and means and standard errors
hormones <- c("#999999", "#E69F00", "brown2")
temps <- c("Blue", "Red")
my_comparisons <- rev(list(c("control","low"),c("control","high-"),c("low","high")))

ggviolin(mito_dat, x = "hormone", y = "state2_state3",
         color = "hormone", palette = hormones,
         short.panel.labs = FALSE,
         font.label = list(size = 14, color = "black"),
         ggtheme = theme_bw())+ 
  geom_jitter(aes(color = hormone), width = 0.2, alpha = 0.6) +  # Adding raw data points with color and transparency
  
  # Adding mean and standard error
  stat_summary(
    mapping = aes(x = hormone, y = state2_state3),
    fun = "mean", geom = "point", size = 3, color = "black", alpha = 0.6) +
  stat_summary(
    mapping = aes(x = hormone, y = state2_state3),
    fun.data = "mean_se", geom = "errorbar", width = 0.1, color = "black", alpha = 0.6) +
  
  stat_compare_means(comparisons = list(c("control","low"), c("control","high"), c("low","high"))) +
  stat_compare_means(method = "anova")+
  labs(x = NULL, y = "state2_state3") +
  labs(color='Hormone Treatment')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none")

  





#####################
### Figures
#####################
hormones <- c("#999999", "#E69F00", "brown2")
temps <- c("Blue", "Red")
my_comparisons <- rev(list(c("control","low"),c("control","high-"),c("low","high")))

########## FIGURE 1: temperature and morphology: effects of developmental treatments (temp, cort, interaction) on time to hatch: -	Effect on temperature: warmer temps faster development 
fig1 <- ggviolin(data_final, x = "temp", y = "days_to_hatch",
                 color = "temp", palette = temps,
                 short.panel.labs = FALSE,
                 font.label = list(size = 14, color = "black"),
                 ggtheme = theme_bw()) + 
  geom_jitter(aes(color = temp), width = 0.2, alpha = 0.4) +  # Adding raw data points
  
  # Adding mean and standard error
  stat_summary(
    mapping = aes(x = temp, y = days_to_hatch),
    fun = "mean", geom = "point", size = 2, color = "black", alpha = 0.6) +
  stat_summary(
    mapping = aes(x = temp, y = days_to_hatch),
    fun.data = "mean_se", geom = "errorbar", width = 0.1, color = "black",  alpha = 0.6) +
  
  stat_compare_means(method = "anova", vjust = -.1, hjust = -3)+ 
  labs(x = "Temperature (°C)", y = "Incubation time (days)") +
  labs(color='Incubation Temperature')+
  scale_y_continuous(breaks=seq(25, 60, 5), limits = c(25,60))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


########### FIGURE 2: temperature and morphology effects of developmental 
# treatments on body size, condition, and mass at hatching;
fig2_A <- ggviolin(data_final, x = "hormone", y = "hatch_svl_mm",
                    color = "hormone", palette = hormones,
                    short.panel.labs = FALSE,
                    font.label = list(size = 14, color = "black"),
                    ggtheme = theme_bw())+ 
  geom_jitter(aes(color = hormone), width = 0.2, alpha = 0.6) +  # Adding raw data points with color and transparency
  
  # Adding mean and standard error
  stat_summary(
    mapping = aes(x = hormone, y = hatch_svl_mm),
    fun = "mean", geom = "point", size = 3, color = "black", alpha = 0.6) +
  stat_summary(
    mapping = aes(x = hormone, y = hatch_svl_mm),
    fun.data = "mean_se", geom = "errorbar", width = 0.1, color = "black", alpha = 0.6) +
  
  stat_compare_means(comparisons = list(c("control","low"), c("control","high"), c("low","high"))) +
  stat_compare_means(method = "anova")+
  labs(x = NULL, y = "Hatchling SVL") +
  labs(color='Hormone Treatment')+
  scale_y_continuous(breaks=seq(14, 22, 2), limits = c(14, 22))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none")
# B) hatchling mass * hormone
Fig2_B <- ggviolin(data_final, x = "hormone", y = "hatch_mass_g",
                    color = "hormone", palette = hormones,
                    short.panel.labs = FALSE,
                    font.label = list(size = 14, color = "black"),
                    ggtheme = theme_bw()) +
  geom_jitter(aes(color = hormone), width = 0.2, alpha = 0.4) +  # Adding raw data points with color and transparency
  
  # Adding mean and standard error
  stat_summary(
    mapping = aes(x = hormone, y = hatch_mass_g),
    fun = "mean", geom = "point", size = 3, color = "black", alpha = .6) +
  stat_summary(
    mapping = aes(x = hormone, y = hatch_mass_g),
    fun.data = "mean_se", geom = "errorbar", width = 0.1, color = "black", alpha = 0.6) +
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

fig2 <- plot_grid(fig2_A, Fig2_B, ncol = 2, align = "h")

########## FIGURE 3: temperature and morphology of hatchlings
fig3_A <- ggviolin(data_final, x = "temp", y = "hatch_mass_g",
                    color = "temp", palette = temps,
                    short.panel.labs = FALSE,
                    font.label = list(size = 14, color = "black"),
                    ggtheme = theme_bw()) +
  
  geom_jitter(aes(color = temp), width = 0.2, alpha = 0.4) + 
  stat_summary(
    mapping = aes(x = temp, y = hatch_mass_g),
    fun = "mean", geom = "point", size = 2, color = "black", alpha = .5) +
  stat_summary(
    mapping = aes(x = temp, y = hatch_mass_g),
    fun.data = "mean_se", geom = "errorbar", width = .1, color = "black", alpha = .5) +
  stat_compare_means(method = "anova") +
  labs(x = "Temperature (°C)", y = "Hatchling mass (g)") +
  labs(color='Temperature Treatment')+
  scale_y_continuous(breaks=seq(0.06, 0.18, .02), limits = c(0.06,0.18)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none")

# body condition and temp
fig3_B <- ggviolin(data_final, x = "temp", y = "",
                   color = "temp", palette = temps,
                   short.panel.labs = FALSE,
                   font.label = list(size = 14, color = "black"),
                   ggtheme = theme_bw()) + 
  geom_jitter(aes(color = temp), width = 0.2, alpha = 0.4) +  
  stat_summary(
    mapping = aes(x = temp, y = hatch_bc_residuals),
    fun = "mean", geom = "point", size = 2, color = "black", alpha = .5 ) +
  stat_summary(
    mapping = aes(x = temp, y = hatch_bc_residuals),
    fun.data = "mean_se", geom = "errorbar", width = 0.1, color = "black", alpha = .5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  stat_compare_means(method = "anova", vjust = -.1, hjust = 0)+ 
  labs(x = "Temperature (°C)", y = "BCI") +
  labs(color='Incubation Treatment')+
  scale_y_continuous(breaks=seq(-0.06, 0.06, 0.02), limits = c(-0.06, 0.06))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = c(0.99, 0.99),  # Positions the legend inside the plot
        legend.justification = c("right", "top"), # Justifies the legend's position
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.box = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8))

fig3 <- plot_grid(fig3_A, fig3_B, ncol = 2, align = "h")

#### FIGURE 4: Mass growth on final cort measurments 
fig4_A <- ggplot(cort_dat, aes(x = hatch_juv3_SVL_growth, y = juv3_CORT_Final_Hormone_ng_mL)) +
  geom_point(aes(color = as.factor(temp), shape = hormone), size = 2) + 
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = 2) +
  labs(x = "SVL growth rate (mm/d)", y = "CORT Final Hormone (ng/mL)", 
       color = "Temperature", shape = "Treatment") +
  scale_color_discrete(name = "Temperature") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = "none")
# B
fig4_B <- ggplot(cort_dat, aes(x = hatch_juv3_MASS_growth, 
                               y = juv3_CORT_Final_Hormone_ng_mL)) +
  geom_point(aes(color = as.factor(temp), shape = hormone), size = 2) + 
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = 2) +
  labs(x = "Mass growth rate (g/d)", y = NULL, 
       color = "Temperature", shape = "Treatment") +
  scale_color_discrete(name = "Temperature") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = c(0.99, 0.99),  # Positions the legend inside the plot
        legend.justification = c("right", "top"), # Justifies the legend's position
        legend.background = element_blank())

fig4 <- plot_grid(fig4_A, fig4_B, ncol = 2, align = "h")


########## FIGURE 5: cort values from cort treatments
fig5 <- ggviolin(data_final, x = "hormone", y = "juv3_CORT_Final_Hormone_ng_mL",
                   color = "hormone", palette = hormones,
                   short.panel.labs = FALSE,
                   font.label = list(size = 14, color = "black"),
                   ggtheme = theme_bw()) +
  geom_jitter(aes(color = hormone), width = 0.2, alpha = 0.4) +  #
  stat_summary(
    mapping = aes(x = hormone, y = juv3_CORT_Final_Hormone_ng_mL),
    fun = "mean", geom = "point", size = 3, color = "black", alpha = .6) +
  stat_summary(
    mapping = aes(x = hormone, y = juv3_CORT_Final_Hormone_ng_mL),
    fun.data = "mean_se", geom = "errorbar", width = 0.1, color = "black", alpha = 0.6) +
  labs(x = "Cort Environment Treatment", y = "Cort Final Hormone (ng/mL)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = c(0.99, 0.99),  # Positions the legend inside the plot
        legend.justification = c("right", "top"), # Justifies the legend's position
        legend.background = element_blank())


# Fig 6 T4 and cort effects on final growth contour plot
# another plot
s <- interp(x = newdata$juv3_T4_corrected_ng_mL, 
            y = newdata$juv3_CORT_Final_Hormone_ng_mL, 
            z = newdata$pred)
image.plot(s, xlab = "T4", ylab = "Cort", 
           las = 1, col = viridis(option = "magma", 50), 
           cex.axis = 1, axis.args = list(cex.axis = 1.5), cex.lab = 1.8)
contour(s, add = TRUE, col = "white")


fig1
fig2
fig3
fig4
fig5