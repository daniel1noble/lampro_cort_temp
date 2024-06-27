## Kwild Lampro Cort Temp Analysis
# WD, packages, data
pacman::p_load(dplyr, tidyverse, ggpubr, lme4, emmeans, here, plotly, performance, car,
               scales, Hmisc, ppcor, ordinal, cowplot)


library(lubridate)
library(Hmisc)
library(ppcor)
library(ordinal)
#########################
# Bring in data from samples and merge together
#########################
# 1) bring in data base data that has morphology and mortality data on across life stage: hatch, juv_1, juv_2
dat_hatch_juv1_2 <- read.csv ("Kwild_code/data/Cort_hatch_juv1_2_OC.csv")

# 2) bring in data that has final morphology measurement and hormone data: "juv_3"
dat_juv3 <- read.csv(file = "Kwild_code/data/Cort_juv_3_np.csv") %>% 
  rename(sex = juv3_Sex)

# 3) final merge juv 3 data that contains hormonal results and final morph with dat_hatch_juv1_2 
merged_data <- merge(dat_hatch_juv1_2, dat_juv3, by.x = "Lizard_ID", 
                     by.y = "juv3_Lizard_ID", all = TRUE)


#########################
# Merge checks
#########################
##Ondi: I don't know why this returns 21 alive since juv 2 sample when clearly there are many more alive than that
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

########
# growth rate calculations: calculated growth rate over lifespan 
# dividing change body size (SVL and mass) by the total number of days elapsed
data_final <- data_final %>% 
  mutate(hatch_juv3_SVL_growth =  (juv3_SVL_mm - hatch_svl_mm)/Juvenile3_Age,
         hatch_juv3_MASS_growth = (juv3_mass_g - hatch_mass_g)/Juvenile3_Age, 
         hatch_juv1_SVL_growth = (juv1_SVL_mm - hatch_svl_mm) / Juvenile1_Age, 
         hatch_juv1_mass_growth = (juv1_mass_g - hatch_mass_g) / Juvenile1_Age, 
         juv1_juv3_svl_growth = (juv3_SVL_mm - juv1_SVL_mm) / (Juvenile3_Age - Juvenile1_Age),
         juv1_juv3_mass_growth = (juv3_mass_g - juv1_mass_g)/ (Juvenile3_Age - Juvenile1_Age))
        
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
# temperature effects days to hatch; no effect on hormone or temp x hormone interaction
hatch_dev_mod <- lm(days_to_hatch ~ temp + hormone, data = data_final)
check_model(hatch_dev_mod) 
Anova(hatch_dev_mod)
summary(hatch_dev_mod)

# faster hatch warmer temps
hatch_dev_mod_emm <- emmeans(hatch_dev_mod, pairwise ~ temp)
plot(hatch_dev_mod_emm)
saveRDS(hatch_dev_mod, "Kwild_code/models/hatch_dev_mod.RDS")

###Mean days to hatch: 28 - 30.93; 23 - 48.31; SD: 28 - 4.82; 23 - 8.38
tapply(data_final$days_to_hatch, data_final$temp, mean, na.rm = TRUE)
tapply(data_final$days_to_hatch, data_final$temp, sd, na.rm = TRUE)

tapply(data_final$days_to_hatch, data_final$hormone, mean, na.rm = TRUE)
tapply(data_final$days_to_hatch, data_final$hormone, sd, na.rm = TRUE)

################
#### 2.	What are the effects of developmental treatments on body size, mass, and condition (scaled mass) of HATCHLINGS
#### SVL: hormone effects body size; no effects on temp or temp x hormone interaction
#### MASS: temperature and hormone effects mass; no interaction effects
#### Condition:  temperature has effect on BCI; no effect on hormone, or temp x hormone interaction

################
#### hatchling SVL
####hormone treatment affects svl; temp treatment NS; hormone*temp NS
SVL_hatch_mod <- lm(hatch_svl_mm ~ temp + hormone, data = data_final)
check_model(SVL_hatch_mod)
Anova(SVL_hatch_mod)
summary(SVL_hatch_mod)

# high hormone lower SVL
##High - control: p = 0.006
##High - low: p = 0.02
##Low - control: p = 0.91
SVL_hatch_mod_emm <- emmeans(SVL_hatch_mod, pairwise ~ hormone)
plot(SVL_hatch_mod_emm)
saveRDS(SVL_hatch_mod, "Kwild_code/models/SVL_hatch_mod.RDS")


#### MASS at hatching
####temp and hormone treatment affect mass; no interaction
Mass_hatch_mod <- lm(hatch_mass_g ~ temp + hormone, data = data_final)
check_model(Mass_hatch_mod)
Anova(Mass_hatch_mod)
summary(Mass_hatch_mod)

# hormone pairwise comparison - high CORT lower mass compared to control
##High - control: p = 0.02
##High - low: p =0.80
##Low - control = 0.06
Mass_hormone_hatch_mod_emm <- emmeans(Mass_hatch_mod, pairwise ~ hormone)
plot(Mass_hormone_hatch_mod_emm)

# temp pairwise comparison - cooler temps higher mass
Mass_temp_hatch_mod_emm <- emmeans(Mass_hatch_mod, pairwise ~ temp)
plot(Mass_temp_hatch_mod_emm)
saveRDS(Mass_hatch_mod, "Kwild_code/models/MASS_hatch_mod.RDS")

####Scaled mass at hatching; results are consistent with OLS model of condition
# temperature has effect on BCI; no effect on hormone (p=0.06), or temp x hormone interaction
scaled_mass_hatch_mod <- lm(ScaledMass_hatch ~ temp + hormone, data = data_final)
check_model(scaled_mass_hatch_mod)
Anova(scaled_mass_hatch_mod)
summary(scaled_mass_hatch_mod)

#28 temp in lower condition at hatch
scaled_mass_hatch_mod_emm <- emmeans(scaled_mass_hatch_mod, pairwise ~ temp)
plot(scaled_mass_hatch_mod_emm)


################
#### 3.	What are the effects of developmental treatments on juvenile body size and condition (at the 3 points measured after hatching)?
#### Juv1_SVL: effect of temperature- cooler temps larger svl
#### Juv1_MASS: low mass with high hormones; low mass with high temperatures
#### Juv1_scaled mass: no treatment effects
#### Juv3_MASS: low mass males
#### Juv3_SVL: high hormone low body size; males smaller than females
#### Juv3_scaled mass: Sex effect: male lower BCI than females - could be size difference

###############
# SVL ANALYSIS: Juv_1, Juv_3 (sex accounted for last measurement)
###############
#### JUV_1 - SVL: 28 temp has smaller SVL
SVL_Juv1_mod <- lm(juv1_SVL_mm ~ temp + hormone + scale(Juvenile1_Age), data = data_final)
check_model(SVL_Juv1_mod)
# effect on temperature and day of hatch on body size 
Anova(SVL_Juv1_mod)
summary(SVL_Juv1_mod)
# plot - cooler temps larger svl
SVL_Juv1_mod_emm <- emmeans(SVL_Juv1_mod, pairwise ~ temp)
plot(SVL_Juv1_mod_emm)
saveRDS(SVL_Juv1_mod, "Kwild_code/models/SVL_Juv1_mod.RDS")

#### JUV_3 (adult) SVL -hormone effect
SVL_Juv3_mod <- lm(juv3_SVL_mm ~ temp + hormone + sex + scale(Juvenile3_Age), data = data_final)
check_model(SVL_Juv3_mod)
# effect of hormones and sex
Anova(SVL_Juv3_mod)
summary(SVL_Juv3_mod)

# plot - high hormone lower svl than control
SVL_Juv3_hormone_mod_emm <- emmeans(SVL_Juv3_mod, pairwise ~ hormone)
plot(SVL_Juv3_hormone_mod_emm)

# plot - males are smaller than females
SVL_Juv3_sex_mod_emm <- emmeans(SVL_Juv3_mod, pairwise ~ sex)
plot(SVL_Juv3_sex_mod_emm)
saveRDS(SVL_Juv3_mod, "Kwild_code/models/SVL_Juv3_mod.RDS")


###############
# MASS ANALYSIS: Juv_1, Juv_3
###############
#### MASS: Juv_1 - effects of hormone and temp treatment 
Mass_Juv1_mod <- lm(juv1_mass_g ~ temp + hormone + scale(Juvenile1_Age), data = data_final)
check_model(Mass_Juv1_mod)
# effect on temp, hormone, and days since hatch; no effects on temp x hormone interaction
summary(Mass_Juv1_mod)
Anova(Mass_Juv1_mod)

# hormone plot - high CORT lower mass than control 
Mass_Juv1_hormone_mod_emm <- emmeans(Mass_Juv1_mod, pairwise ~ hormone)
plot(Mass_Juv1_hormone_mod_emm)

# temp plot - 28 temp has lower mass
Mass_Juv1_temp_mod_emm <- emmeans(Mass_Juv1_mod, pairwise ~ temp)
plot(Mass_Juv1_temp_mod_emm)
saveRDS(Mass_Juv1_mod, "Kwild_code/models/Mass_Juv1_mod.RDS")


#### MASS: Juv_3 sex effect
Mass_Juv3_mod <- lm(juv3_mass_g ~ hormone + temp + sex + scale(Juvenile3_Age), data = data_final)
check_model(Mass_Juv3_mod)
Anova(Mass_Juv3_mod)
summary(Mass_Juv3_mod)

# temp plot - males smaller than females
Sex_Juv3_temp_mod_emm <- emmeans(Mass_Juv3_mod, pairwise ~ hormone)
plot(Sex_Juv3_temp_mod_emm)
saveRDS(Mass_Juv3_mod, "Kwild_code/models/Mass_Juv3_mod.RDS")

###############
# Scaled mass ANALYSIS: Juv_1, Juv_3
###############
## Juvenile 1- no treatment effects on scaled mass juvenile 1
condition_juv_mod <- lm(ScaledMass_juv1 ~ temp + hormone+ scale(Juvenile1_Age), data = data_final)
Anova(condition_juv_mod)
summary(condition_juv_mod)
check_model(condition_juv_mod)

### Juvenile 3 scaled mass; no treatment affects; sex effect
condition_juv3_mod <- lm(ScaledMass_juv3 ~ temp + hormone + sex + scale(Juvenile3_Age), data = data_final)
check_model(condition_juv3_mod)
# no treatment effects on bci juv1
Anova(condition_juv3_mod)
summary(condition_juv3_mod)

# sex plot - males have lower condition at adulthood
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
###To get sample sizes
survival_table <- table(survival$liz_status_db, survival$temp)
survival_table <- table(survival$liz_status_db, survival$hormone)
# analysis
survival_day_hormone <- kruskal.test(liz_status_db ~ hormone, data = survival)
survival_day_temp <- kruskal.test(liz_status_db ~ temp, data = survival)
# save analysis
saveRDS(survival_day_hormone, "Kwild_code/models/survival_day_hormone.RDS")
saveRDS(survival_day_temp, "Kwild_code/models/survival_day_temp.RDS")


################
#### 5.	Growthrate
####hatching to juvenile 1 measurements
hatch_juv1_SVL_growth_mod <- lm(hatch_juv1_SVL_growth ~ hormone + temp, data = data_final)
Anova(hatch_juv1_SVL_growth_mod)
summary(hatch_juv1_SVL_growth_mod)

hatch_juv1_mass_growth_mod <- lm(hatch_juv1_mass_growth ~ hormone + temp, data = data_final)
Anova(hatch_juv1_mass_growth_mod)
summary(hatch_juv1_mass_growth_mod)
hatch_juv1_mass_growth_mod_emm <- emmeans(hatch_juv1_mass_growth_mod, pairwise ~ hormone)
plot(hatch_juv1_mass_growth_mod_emm)

tapply(data_final$hatch_juv1_SVL_growth, data_final$temp, mean, na.rm = TRUE)
tapply(data_final$hatch_juv1_SVL_growth, data_final$temp, sd, na.rm = TRUE)

tapply(data_final$hatch_juv1_mass_growth, data_final$temp, mean, na.rm = TRUE)
tapply(data_final$hatch_juv1_mass_growth, data_final$temp, sd, na.rm = TRUE)


#### Look at overall growth from hatch to adulthood
################
###Overall growth hatch to juvenile 3; sex included
hatch_juv3_SVL_growth_mod <- lm(hatch_juv3_SVL_growth ~hormone + temp + sex, data = data_final)
Anova(hatch_juv3_SVL_growth_mod)
summary(hatch_juv3_SVL_growth_mod)
check_model(hatch_juv3_SVL_growth_mod)

hatch_juv3_SVL_growth_mod_emm <- emmeans(hatch_juv3_SVL_growth_mod, pairwise ~ temp)
plot(hatch_juv3_SVL_growth_mod_emm)

hatch_juv3_MASS_growth_mod <- lm(hatch_juv3_MASS_growth ~ hormone + temp + sex, data = data_final)
Anova(hatch_juv3_MASS_growth_mod)
summary(hatch_juv3_MASS_growth_mod)


###growth between juvenile and adult
juv1_juv3_svl_growth_mod  <- lm(juv1_juv3_svl_growth ~ hormone + temp + sex, data = data_final)
Anova(juv1_juv3_svl_growth_mod)
summary(juv1_juv3_svl_growth_mod)

juv1_juv3_mass_growth_mod  <- lm(juv1_juv3_mass_growth ~ hormone + temp + sex, data = data_final)
Anova(juv1_juv3_mass_growth_mod)
summary(juv1_juv3_mass_growth_mod)

tapply(data_final$juv1_juv3_svl_growth, data_final$temp, mean, na.rm = TRUE)
tapply(data_final$juv1_juv3_svl_growth, data_final$temp, sd, na.rm = TRUE)

########################
### 6: CORT, T4, testosterone
###OC changed >= to > to remove values of 4ul and below

########
# CORT
# data
cort_dat <- data_final %>%
  drop_na(juv3_CORT_Final_Hormone_ng_mL, temp, hormone, juv3_HandlingTime_sec) %>% 
  filter(juv3_Sample_volume_ul > 4 ) # remove values below 4 ul juv3_Sample_volume_ul

plot(juv3_CORT_Final_Hormone_ng_mL ~ hormone, data = cort_dat)
plot(juv3_CORT_Final_Hormone_ng_mL ~ temp, data = cort_dat)

###Checking plate ID and handling time for effects on CORT; handling time = ns; plate ID is p=.007
cort_factors <- lm(log(juv3_CORT_Final_Hormone_ng_mL) ~ juv3_HandlingTime_sec + Plate_CORT_juv3, data = cort_dat)
Anova(cort_factors)
summary(cort_factors)


######## 1) treatments- how does cort vary across temp and hormone treatment; handling time not used because ns in model of CORT factors
cort_development_mod <- lm(log(juv3_CORT_Final_Hormone_ng_mL) ~ hormone + temp +  sex + scale(Juvenile3_Age) + Plate_CORT_juv3, data = cort_dat)
check_model(cort_development_mod)
Anova(cort_development_mod)
summary(cort_development_mod)
##low CORT treatment has higher baseline CORT than control treatment; but ns
cort_development_mod_emm <- emmeans(cort_development_mod, pairwise ~ hormone)
plot(cort_development_mod_emm)

###Visualize baseline CORT levels
Boxplot(cort_dat$juv3_CORT_Final_Hormone_ng_mL, cort_dat$hormone, na.action = na.exclude)

###summarise data for bar plot
sum.dat <- cort_dat  %>%
  
  group_by(hormone) %>%
  
  summarise(mean.value = mean(juv3_CORT_Final_Hormone_ng_mL, na.rm = TRUE),
            
            sd = sd(juv3_CORT_Final_Hormone_ng_mL, na.rm = TRUE),
            
            se=mean.value/sqrt(n()))



Bar.plot <- ggplot(sum.dat, aes(x=hormone, y=mean.value, fill=hormone)) +
  
  geom_bar(stat="identity", color="black",
           
           position=position_dodge()) +
  
  geom_errorbar(aes(ymin=mean.value -se, ymax=mean.value + se), width=.2,
                
                position=position_dodge(.9))+
  
  labs(title="Figure 3. Baselinecorticosterone levels in adults", y="CORT (pg/mg)", x="Corticosteorne treatment")+
  
  theme_classic()+
  
  
  scale_fill_manual(values = c("white", "grey46", "black"))

Bar.plot + theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
                 axis.title = element_text(size = 14), legend.text = element_text(size = 14), 
                 legend.title = element_text(size = 14))


### using a glm with gamma distribution; results equivalent
Gam_cort <- glm(juv3_CORT_Final_Hormone_ng_mL ~ hormone + temp + Plate_CORT_juv3 + sex + scale(Juvenile3_Age), family = Gamma (link = log), data = cort_dat)
Anova(Gam_cort)
summary(Gam_cort)
check_model(Gam_cort)
Gam_cort_emm <- emmeans(Gam_cort, pairwise ~ hormone)
plot(Gam_cort_emm)


##calculate residuals of plate ID and CORT levels
cort_residuals <- lm(juv3_CORT_Final_Hormone_ng_mL ~ Plate_CORT_juv3, data = cort_dat, na.action=na.exclude) # fit the model for residuals
cort_dat$cort_residuals <- residuals(cort_residuals, na.action=na.exclude) # Save the residual values


####Body size and condition at juvenile 3 related to CORT levels; 'handling time not included because it did not affect baseline CORT levels as above)
##CORT and scaled mass at juvenile 3; treatments not included because they have no effect on body condition at this time point
cort_body_condition_mod <- lm(ScaledMass_juv3 ~ cort_residuals + sex + scale(Juvenile3_Age) + sex, data = cort_dat)
Anova(cort_body_condition_mod)
summary(cort_body_condition_mod)
check_model(cort_body_condition_mod)

###SVL with CORT; adult SVL positive association with CORT
cort_svl_juvenile3_mod <- lm(juv3_SVL_mm ~ cort_residuals + sex + scale(Juvenile3_Age), data = cort_dat)
Anova(cort_svl_juvenile3_mod)
summary(cort_svl_juvenile3_mod)
check_model(cort_svl_juvenile3_mod)

###Mass at juvenile 3 tested against CORT levels; CORT positive association with adult mass
cort_mass_juvenile3_mod <- lm(juv3_mass_g ~ cort_residuals + sex + scale(Juvenile3_Age), data = cort_dat)
Anova(cort_mass_juvenile3_mod)
summary(cort_mass_juvenile3_mod)
check_model(cort_mass_juvenile3_mod)

################
# T4
################ 

t4_dat <- data_final %>%
  drop_na(juv3_T4_corrected_ng_mL, temp, hormone, Lizard_ID) %>% 
  filter(juv3_Sample_volume_ul > 4 ) # remove values below 4 ul juv3_Sample_volume_ul

###Checking plate ID and handling time for effects on T4; handling time = ns; plate ID is NS
T4_factors <- lm((log(juv3_T4_corrected_ng_mL)) ~ juv3_HandlingTime_sec + juv3_T4_plate, data = t4_dat)
Anova(T4_factors)
summary(T4_factors)

# 1) treatments -  does T4 vary across temp and hormone treatment: NS
T4_temp_hormone_sex_mod <- lm(log(juv3_T4_corrected_ng_mL) ~ hormone + temp + sex + scale(Juvenile3_Age), data = t4_dat)
check_model(T4_temp_hormone_sex_mod)
Anova(T4_temp_hormone_sex_mod)
summary(T4_temp_hormone_sex_mod)

# sex differences: females have higher T4
T4_temp_hormone_sex_mod_emm <- emmeans(T4_temp_hormone_sex_mod, pairwise ~ sex)
plot(T4_temp_hormone_sex_mod_emm)
saveRDS(T4_temp_hormone_sex_mod, "Kwild_code/models/T4_temp_hormone_sex_mod.RDS")


# 2) SVL and associations with T4 and cort
T4_SVL_mod <- lm(juv3_SVL_mm  ~ log(juv3_T4_corrected_ng_mL) + sex + scale(Juvenile3_Age), data = t4_dat)
check_model(T4_SVL_mod )
summary(T4_SVL_mod )
Anova(T4_SVL_mod )

# 3) MASS and associations between T4
T4_mass_growth_mod <- lm(juv3_mass_g  ~ log(juv3_T4_corrected_ng_mL) + sex + scale(Juvenile3_Age), 
                        data = t4_dat)
check_model(T4_mass_growth_mod)
Anova(T4_mass_growth_mod)
summary(T4_mass_growth_mod)
saveRDS(T4_mass_growth_mod, "Kwild_code/models/T4_mass_growth_mod.RDS")


########################
###Testosterone
test_dat <- data_final %>%
  drop_na(juv3_Testosterone_Final_ng_ml) %>% 
  filter(juv3_Sample_volume_ul > 4 ) # remove values below 4 ul juv3_Sample_volume_ul
test_dat$rank_test <- rank(test_dat$juv3_Testosterone_Final_ng_ml, ties.method = "first")


###Look at the effects of plate and handling time on testosterone levels; handling time affects T levels (p=0.02); plate is NS
Test_factor <- lm (log(juv3_Testosterone_Final_ng_ml) ~ juv3_Testosterone_plate+ juv3_HandlingTime_sec, data = test_dat)
Anova(Test_factor)
summary(Test_factor)

###Get residuals from model of T levels and handling time and save as new variable
test_residuals <- lm((log(juv3_Testosterone_Final_ng_ml)) ~ juv3_HandlingTime_sec, data = test_dat)
Anova(test_residuals)
summary(test_residuals)
test_dat$test_residuals <- residuals(test_residuals)


###Does testosterone vary based on treatment; use non-parametric test (t data determined from linear standard curve)
##No effect of treatments on adult T levels 
test_hormone <- kruskal.test(test_residuals ~ hormone, data = test_dat)
test_hormone_n <- table(test_dat$hormone)
test_temp <- kruskal.test(test_residuals ~ temp, data = test_dat)
test_temp_n <- table(test_dat$temp)

##correlation between testosterone and CORT levels
##T and CORT levels positively correlated
cor.test(test_dat$juv3_Testosterone_Final_ng_ml, test_dat$juv3_CORT_Final_Hormone_ng_mL, method = c("spearman"))

plot(test_dat$juv3_CORT_Final_Hormone_ng_mL, test_dat$juv3_Testosterone_Final_ng_ml)
cort_T_mod <- lm(juv3_Testosterone_Final_ng_ml ~ juv3_CORT_Final_Hormone_ng_mL, data = test_dat)
abline(cort_T_mod)

########################
### 7: Mito function
########################
# bring in data, rename, and save
###if looking at relationships with hormone levels - use cort_dat or equivalent to exclude samples with 
##low volumes
mito_dat <- data_final %>%
  dplyr::select(c("Lizard_ID", "clutch", "temp", "hormone", "Juvenile3_Age", "sex","hatch_juv3_MASS_growth", "hatch_juv3_SVL_growth",
                  "juv3_inject_time_sec", "juv3_liver_time_sec", "juv3_oroboros",
                  "juv3_T4_plate", "juv3_HandlingTime_sec", "juv3_T4_corrected_ng_mL", 
                  "juv3_CORT_Final_Hormone_ng_mL", "juv3_basal_corrected.pmol..sec.ng..",
                  "juv3_chamber","juv3_basal_corrected.pmol..sec.ng..", 
                  "juv3_adp_corrected.pmol..sec.ng..", "juv3_oligo_corrected.pmol..sec.ng..", 
                  "juv3_fccp_corrected.pmol..sec.ng..", "juv3_RCR.L.R.", "juv3_RCR.R.ETS.",
                  "juv3_oroboros_comments", "juv3_RCR.L.ETS.", "juv3_date", "Plate_CORT_juv3",
                  "juv3_Sample_volume_ul", "hatch_svl_mm", "juv3_Testosterone_Final_ng_ml",
                  "juv3_mass_g", "juv3_SVL_mm",  "ScaledMass_juv3",
                  "juv3_HandlingTime_sec",  "juv3_chamber")) %>% 
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
         OXPHOS = fccp_corrected_pmol - oligo_corrected_pmol,
         RCR = adp_corrected_pmol/oligo_corrected_pmol) %>% 

  filter(Lizard_ID != c("LD784_21"))
####LD784_21 had extremely high values and lab notes suggest should be filtered from analyses 

view(mito_dat)
###Checking factors that could affect oxygen consumption measurements - handling times, Oroboros, and chamber
###'chamber' is nested within Oroboros; there are three measurements of time (time to inject the drug, time to euthanize the animal
###'and time to collect liver, have opted to only use time to collect liver because it will incorporate all other time measurements and 
###'likely all time measures are correlated and so shouldn't all be going in the model anyway)
###all factors NS


#RCR - all NS
RCR_factors <- lm(RCR ~ juv3_liver_time_sec + juv3_oroboros/chamber, data = mito_dat)
Anova(RCR_factors)
summary(RCR_factors)
hist(mito_dat$RCR)

###Basal; all factors NS
basal_factors <- lm(basal_corrected_pmol ~ juv3_liver_time_sec + juv3_oroboros/chamber, data = mito_dat)
Anova(basal_factors)
summary (basal_factors)
hist(mito_dat$basal_corrected_pmol)

##ADP - all factors NS
ADP_factors <- lm(adp_corrected_pmol ~ juv3_liver_time_sec + juv3_oroboros/chamber, data = mito_dat)
Anova(ADP_factors)
summary (ADP_factors)

##Oligo - all factors NS but time to collect liver p = 0.06
Oligo_factors <- lm (oligo_corrected_pmol ~ juv3_liver_time_sec + juv3_oroboros/chamber, data = mito_dat)
Anova(Oligo_factors)
summary (Oligo_factors)
hist(mito_dat$oligo_corrected_pmol)


#####Effects of treatments on mitochondrial function 
mito_dat <- mito_dat %>% 
  mutate(mass_basal =  basal_corrected_pmol / juv3_mass_g,
         mass_ADP = adp_corrected_pmol / juv3_mass_g,
         mass_oligo = oligo_corrected_pmol / juv3_mass_g) %>%
  mutate(mass_RCR = mass_ADP / mass_oligo)

         
###Basal - interaction removed ns and removed; sex p = 0.005 males higher than females
Basal_treatment_mod <- lm(basal_corrected_pmol ~ temp + hormone + scale(Juvenile3_Age) + sex + juv3_mass_g + juv3_oroboros/chamber, data = mito_dat)
Anova(Basal_treatment_mod)
summary(Basal_treatment_mod)
check_model(Basal_treatment_mod)

Basal_treatment_mod_emm <- emmeans(Basal_treatment_mod, pairwise ~ sex)
plot(Basal_treatment_mod_emm)

Boxplot(mito_dat$mass_basal, mito_dat$sex, na.action = na.exclude)
saveRDS(Basal_treatment_mod, "Kwild_code/models/Basal_treatment_mod.RDS")

###ADP - treatment effects; interaction ns and removed; sex p = 0.045 males higher than females
ADP_treatment_mod <- lm(adp_corrected_pmol  ~ temp + hormone + sex + scale(Juvenile3_Age) + juv3_mass_g + juv3_oroboros/chamber, data = mito_dat)
Anova(ADP_treatment_mod)
summary(ADP_treatment_mod)
check_model(ADP_treatment_mod)

ADP_treatment_mod_emm <- emmeans(ADP_treatment_mod, pairwise ~ sex)
plot(ADP_treatment_mod_emm)

saveRDS(ADP_treatment_mod, "Kwild_code/models/ADP_treatment_mod.RDS")

##Oligo mito by treatment
Oligo_treatment_mod <- lm(oligo_corrected_pmol ~ temp + hormone + sex + scale(Juvenile3_Age) + juv3_mass_g + juv3_oroboros/chamber, data = mito_dat)
Anova(Oligo_treatment_mod)
summary(Oligo_treatment_mod)
check_model(Oligo_treatment_mod)

Oligo_treatment_mod_emm <- emmeans(Oligo_treatment_mod, pairwise ~ sex)
plot(Oligo_treatment_mod_emm)

saveRDS(Oligo_treatment_mod, "Kwild_code/models/Oligo_treatment_mod.RDS")

##RCR calculated as state 3 (ADP)/ state 4(olgio)
RCR_mod <- lm(RCR ~ temp + hormone + sex + scale(Juvenile3_Age) + juv3_mass_g + juv3_oroboros/chamber, data = mito_dat)
Anova(RCR_mod)
summary (RCR_mod)
check_model(RCR_mod)

saveRDS(RCR_mod, "Kwild_code/models/RCR_mod.RDS")

###For models with hormone levels - need to filter n=7 data points from 
###samples with less than 7ul 
mito_hormone <- mito_dat %>%
  drop_na(CORT_Final_Hormone_ng_mL, temp, hormone, HandlingTime_sec) %>% 
  filter(juv3_Sample_volume_ul > 4 ) # remove values below 4 ul juv3_Sample_volume_ul
view(mito_hormone)


##CORT and plate ID and handling time; NS
plot(mito_hormone$CORT_Final_Hormone_ng_mL, mito_hormone$HandlingTime_sec)
B_handling <- lm(CORT_Final_Hormone_ng_mL ~ HandlingTime_sec +  Plate_CORT_juv3, data = mito_hormone)
Anova(B_handling)
abline(B_handling)

CORT_residuals_mito <- lm(CORT_Final_Hormone_ng_mL ~ Plate_CORT_juv3, data = mito_hormone)
Anova(CORT_residuals_mito)

###T4 plate ID and handling time; NS
T4_residuals_mod <- lm(T4_corrected_ng_mL ~ HandlingTime_sec + T4_plate_ID, data = mito_hormone, na.action = na.exclude)
Anova(T4_residuals_mod)
summary(T4_residuals_mod)

###Associations between mitochondrial bioenergetics and hormones levels

mito_hormone <- mito_hormone %>% 
  mutate(mass_basal =  basal_corrected_pmol / juv3_mass_g,
         mass_ADP = adp_corrected_pmol / juv3_mass_g,
         mass_oligo = oligo_corrected_pmol / juv3_mass_g) %>%
  mutate(mass_RCR = mass_ADP / mass_oligo)

##Basal 
Basal_hormone_mod <- lm(mass_basal ~ log(CORT_Final_Hormone_ng_mL) + log(T4_corrected_ng_mL) + sex + juv3_oroboros/chamber, data = mito_hormone)
Anova(Basal_hormone_mod)
summary(Basal_hormone_mod)
check_model(Basal_hormone_mod)
Boxplot(mito_hormone$mass_basal, mito_hormone$sex, na.action = na.exclude)
saveRDS(Basal_hormone_mod, 'Kwild_code/models/Basal_hormone_mod.RDS')

###ADP -  T4 positive association with ADP mito 
ADP_hormone_mod <- lm(mass_ADP ~ log(CORT_Final_Hormone_ng_mL) + log(T4_corrected_ng_mL) + sex + juv3_oroboros/chamber, data = mito_hormone)
Anova(ADP_hormone_mod)
summary(ADP_hormone_mod)
check_model(ADP_hormone_mod)
saveRDS(ADP_hormone_mod, 'Kwild_code/models/ADP_hormone_mod.RDS')


##oligo - T4 near positive association
Oligo_hormone_mod <- lm(mass_oligo ~ log(CORT_Final_Hormone_ng_mL) + log(T4_corrected_ng_mL) + sex + juv3_oroboros/chamber, data = mito_hormone)
Anova(Oligo_hormone_mod)
summary(Oligo_hormone_mod)
check_model(Oligo_hormone_mod)
saveRDS(Oligo_hormone_mod, 'Kwild_code/models/Oligo_hormone_mod.RDS')


#RCR - 
RCR_hormone_mod <- lm(mass_RCR ~ log(CORT_Final_Hormone_ng_mL) + log(T4_corrected_ng_mL) + sex + juv3_oroboros/chamber, data = mito_hormone)
Anova(RCR_hormone_mod)
summary (RCR_hormone_mod)
check_model(RCR_hormone_mod)
saveRDS(RCR_hormone_mod, 'Kwild_code/models/RCR_hormone_mod.RDS')


###Growth and body size and mito function 
###Basal - postive association with basal and growth (p = 0.02)
Mass_basal <- lm(hatch_juv3_MASS_growth ~ basal_corrected_pmol + sex + log(CORT_Final_Hormone_ng_mL) + log(T4_corrected_ng_mL), data = mito_hormone)
Anova(Mass_basal)
summary(Mass_basal)
check_model(Mass_basal)
saveRDS(Mass_basal, "Kwild_code/models/Mass_basal.RDS")



SVL_basal <- lm(hatch_juv3_SVL_growth~ basal_corrected_pmol + sex + log(CORT_Final_Hormone_ng_mL) + log(T4_corrected_ng_mL), data = mito_hormone)
Anova(SVL_basal)
summary(SVL_basal)
check_model(SVL_basal)
saveRDS(SVL_basal, "Kwild_code/models/SVL_basal.RDS")


##ADP 
Mass_Growth_ADP <- lm(hatch_juv3_MASS_growth ~ adp_corrected_pmol + sex + log(CORT_Final_Hormone_ng_mL) + log(T4_corrected_ng_mL), data = mito_hormone)
Anova(Mass_Growth_ADP)
summary(Mass_Growth_ADP)
check_model(Mass_Growth_ADP)
saveRDS(Mass_Growth_ADP, "Kwild_code/models/Mass_Growth_ADP.RDS")

SVL_Growth_ADP <- lm(hatch_juv3_SVL_growth ~ adp_corrected_pmol + log(CORT_Final_Hormone_ng_mL) + log(T4_corrected_ng_mL) + sex, data = mito_hormone)
Anova(SVL_Growth_ADP)
summary(SVL_Growth_ADP)
check_model(SVL_Growth_ADP)
saveRDS(SVL_Growth_ADP, "Kwild_code/models/SVL_Growth_ADP.RDS")



###Oligo 
Mass_Growth_Oligo <- lm(hatch_juv3_MASS_growth ~ oligo_corrected_pmol + log(CORT_Final_Hormone_ng_mL) + log(T4_corrected_ng_mL) + sex, data = mito_hormone)
Anova(Mass_Growth_Oligo)
summary(Mass_Growth_Oligo)
saveRDS(Mass_Growth_Oligo, "Kwild_code/models/Mass_Growth_Oligo.RDS")


SVL_Growth_Oligo <- lm(hatch_juv3_SVL_growth ~ oligo_corrected_pmol + sex + log(CORT_Final_Hormone_ng_mL) + log(T4_corrected_ng_mL), data = mito_hormone)
Anova(SVL_Growth_Oligo)
summary(SVL_Growth_Oligo)
saveRDS(SVL_Growth_Oligo, "Kwild_code/models/SVL_Growth_Oligo.RDS")


###RCR
Mass_Growth_RCR <- lm (hatch_juv3_MASS_growth~ RCR + sex + log(CORT_Final_Hormone_ng_mL) + log(T4_corrected_ng_mL), data = mito_hormone)
Anova(Mass_Growth_RCR)
summary(Mass_Growth_RCR)
saveRDS(Mass_Growth_RCR, "Kwild_code/models/Mass_Growth_RCR.RDS")

SVL_Growth_RCR <- lm (hatch_juv3_SVL_growth ~ RCR + sex + log(CORT_Final_Hormone_ng_mL) + log(T4_corrected_ng_mL), data = mito_hormone)
Anova(SVL_Growth_RCR)
summary(SVL_Growth_RCR)
saveRDS(SVL_Growth_RCR, "Kwild_code/models/SVL_Growth_RCR.RDS")



cor_mito <- mito_dat %>%
  dplyr::select(c("oligo_corrected_pmol", "adp_corrected_pmol", "basal_corrected_pmol", "RCR")) %>%
  filter(if_all(c(oligo_corrected_pmol), complete.cases))

cor_mat <- rcorr(as.matrix(cor_mito))
cor_mat$r
cor_mat$P

view(cor_mito)

flatten_cor_mat <- function (cor_mat, pmat) {
  ut <- upper.tri(cor_mat)
  data.frame(
    row = rownames (cor_mat)[row(cor_mat)[ut]],
    column = rownames (cor_mat)[col(cor_mat)[ut]],
    cor = (cor_mat)[ut],
    P = pmat[ut]
  )
}

adjust_cor <- flatten_cor_mat(cor_mat$r, (p.adjust (cor_mat$P, method = 'bonferroni')))
print(adjust_cor, digits = 5)
                  
cor.test(cor_mito$oligo_corrected_pmol, cor_mito$adp_corrected_pmol, method = ("pearson"))
                  





############################################################
#### !!NEW Chris Freeze - analysis
### Chris asked about the scaling relationship between basal and OXPHOS #respiration and body mass
# Basal
Basal_treatment_mod <- lm(basal_corrected_pmol ~ temp + hormone + scale(Juvenile3_Age) + sex + juv3_mass_g + juv3_oroboros/chamber, 
                          data = mito_dat)
basal_mass_lm <- lm(basal_corrected_pmol ~ Juvenile3_Age, data = mito_dat)
# Get summary of the linear model
model_summary <- summary(basal_mass_lm)
intercept <- coef(basal_mass_lm)[1]
slope <- coef(basal_mass_lm)[2]
p_value <- model_summary$coefficients[2, 4]
r_squared <- model_summary$r.squared
# Create a scatter plot
plot(mito_dat$Juvenile3_Age, mito_dat$basal_corrected_pmol,
     main = "Relationship between Juvenile3 Age and Basal Corrected pmol",
     xlab = "Juvenile3 Age",
     ylab = "Basal Corrected pmol",
     pch = 19,
     col = "blue")
# Add the regression line
abline(basal_mass_lm, col = "red", lwd = 2)
# Annotate the plot with regression results
text(x = 440, 
     y = .7,
     labels = paste("Intercept:", 
                    round(intercept, 2), 
                    "\nSlope:", round(slope, 2), 
                    "\nP-value:", round(p_value, 4), 
                    "\nR-squared:", round(r_squared, 2)),
     pos = 4, cex = 0.8, col = "black")

######## ADP
ADP_treatment_mod <- lm(adp_corrected_pmol  ~ temp + hormone + sex + scale(Juvenile3_Age) + juv3_mass_g + juv3_oroboros/chamber, data = mito_dat)

OXPHOS_mass_lm <- lm(adp_corrected_pmol ~ Juvenile3_Age, data = mito_dat)
# Get summary of the linear model
model_summary <- summary(OXPHOS_mass_lm)
intercept <- coef(OXPHOS_mass_lm)[1]
slope <- coef(OXPHOS_mass_lm)[2]
p_value <- model_summary$coefficients[2, 4]
r_squared <- model_summary$r.squared
# Create a scatter plot
plot(mito_dat$Juvenile3_Age, mito_dat$adp_corrected_pmol,
     main = "Relationship between Juvenile3 Age and Basal Corrected pmol",
     xlab = "Juvenile3 Age",
     ylab = "Basal Corrected pmol",
     pch = 19,
     col = "blue")
# Add the regression line
abline(OXPHOS_mass_lm, col = "red", lwd = 2)
# Annotate the plot with regression results
text(x = 440, 
     y = 6,
     labels = paste("Intercept:", 
                    round(intercept, 2), 
                    "\nSlope:", round(slope, 2), 
                    "\nP-value:", round(p_value, 4), 
                    "\nR-squared:", round(r_squared, 2)),
     pos = 4, cex = 0.8, col = "black")










############# ############# ############# ############# ############# 
############# Kris Figures #############
###############################################################
####### !!!KRIS Figure 2 model - I know this is above but 
####### wanted to follow the code numbers you gave me
###############################################################
library(cowplot) # combines your plots into panels 
cort_dat$hormone <- factor(cort_dat$hormone, levels = c("control", "low", "high"))
scale_of_inference_dat <- data_final %>% dplyr::select(Lizard_ID, temp, hormone)
scale_of_inference <- scale_of_inference_dat %>%
  group_by(temp, hormone) %>%
  summarise(count = n()) %>%
  spread(hormone, count, fill = 0)


####################
###### FIGURE 2) MASS - Hormone 
######
### hatch model 
Mass_hatch_mod <- lm(hatch_mass_g ~ temp + hormone, data = data_final)

### Data from model for mean and wiskers 'cort_development_mod'  output
# hormone data - hatchling
Mass_hatch_mod_emm_hormone <- emmeans(Mass_hatch_mod, pairwise ~ hormone)
Mass_hatch_dat_emm_hormone <- as.data.frame(Mass_hatch_mod_emm_hormone$emmeans) %>% 
  mutate(Age = "hatchling",
         Treatment = "Hormone")
# sample size per group
Mass_hatch_dat_sum_hormone <- data_final %>%
  dplyr::select(Lizard_ID, hormone, hatch_mass_g) %>% 
  dplyr::filter(!is.na(hatch_mass_g)) %>%
  group_by(hormone) %>%
  summarise(n = n())
# df for plotting
mass_hatch_hormone <- left_join(Mass_hatch_dat_emm_hormone, 
                                Mass_hatch_dat_sum_hormone,
                                by = "hormone")


# temp data  - hatchling
Mass_hatch_mod_emm_temp <- emmeans(Mass_hatch_mod, pairwise ~ temp)
Mass_hatch_dat_emm_temp <- as.data.frame(Mass_hatch_mod_emm_temp$emmeans) %>% 
  mutate(Age = "hatchling",
         Treatment = "Temperature")
# sample size per group
Mass_hatch_dat_sum_temp <- data_final %>%
  dplyr::select(Lizard_ID, temp, hatch_mass_g) %>% 
  dplyr::filter(!is.na(hatch_mass_g)) %>%
  group_by(temp) %>%
  summarise(n = n())
# df for plotting
mass_hatch_temp<- left_join(Mass_hatch_dat_emm_temp, 
                            Mass_hatch_dat_sum_temp,
                            by = "temp") 

######
### JUV1 model 
Mass_Juv1_mod <- lm(juv1_mass_g ~ temp + hormone + scale(Juvenile1_Age), data = data_final)

# Hormone - Juvenile
Mass_Juv1_mod_emm_hormone <- emmeans(Mass_Juv1_mod, pairwise ~ hormone)
Mass_Juv1_emm_hormone <- as.data.frame(Mass_Juv1_mod_emm_hormone$emmeans) %>% 
  mutate(Age = "juvenile",
         Treatment = "Hormone")
# sample size per group
Mass_Juv1_dat_sum_hormone <- data_final %>%
  dplyr::select(Lizard_ID, hormone,temp,juv1_mass_g) %>% 
  dplyr::filter(!is.na(juv1_mass_g)) %>%
  group_by(hormone) %>%
  summarise(n = n())
# df for plotting
mass_juv_hormone <- left_join(Mass_Juv1_emm_hormone, 
                              Mass_Juv1_dat_sum_hormone,
                              by = "hormone")

# temp data  - juvenile
Mass_Juv1_mod_emm_temp <- emmeans(Mass_Juv1_mod, pairwise ~ temp)
Mass_Juv1_dat_emm_temp <- as.data.frame(Mass_Juv1_mod_emm_temp$emmeans) %>% 
  mutate(Age = "juvenile",
         Treatment = "Temperature")
# sample size per group
Mass_Juv1_dat_sum_temp <- data_final %>%
  dplyr::select(Lizard_ID, temp, juv1_mass_g) %>% 
  dplyr::filter(!is.na(juv1_mass_g)) %>%
  group_by(temp) %>%
  summarise(n = n())
# df for plotting
mass_juv_temp<- left_join(Mass_Juv1_dat_emm_temp, 
                          Mass_Juv1_dat_sum_temp,
                          by = "temp")


######
### MASS: Adult mod 
Mass_Juv3_mod <- lm(juv3_mass_g ~ hormone + temp + sex + 
                      scale(Juvenile3_Age), data = data_final)

# Data from model for mean and wiskers 'cort_development_mod'  output
Mass_Juv3_mod_emm_hormone <- emmeans(Mass_Juv3_mod, pairwise ~ hormone)
Mass_Juv3_emm_hormone <- as.data.frame(Mass_Juv3_mod_emm_hormone$emmeans) %>% 
  mutate(Age = "adult",
         Treatment = "Hormone")
# sample size per group
Mass_Juv3_dat_sum_hormone <- data_final %>%
  dplyr::select(Lizard_ID, hormone,temp,juv3_mass_g) %>% 
  dplyr::filter(!is.na(juv3_mass_g)) %>%
  group_by(hormone) %>%
  summarise(n = n())
# df for plotting
mass_adult_hormone <- left_join(Mass_Juv3_emm_hormone, 
                                Mass_Juv3_dat_sum_hormone,
                                by = "hormone")

# temp data  - Adult
Mass_Juv3_mod_emm_temp <- emmeans(Mass_Juv3_mod, pairwise ~ temp)
Mass_Juv3_dat_emm_temp <- as.data.frame(Mass_Juv3_mod_emm_temp$emmeans) %>% 
  mutate(Age = "adult",
         Treatment = "Temperature")
# sample size per group
Mass_Juv3_dat_sum_temp <- data_final %>%
  dplyr::select(Lizard_ID, temp, juv3_mass_g) %>% 
  dplyr::filter(!is.na(juv3_mass_g)) %>%
  group_by(temp) %>%
  summarise(n = n())
# df for plotting
mass_adult_temp<- left_join(Mass_Juv3_dat_emm_temp, 
                            Mass_Juv3_dat_sum_temp,
                            by = "temp")

#####
### final df for plotting
# hormone data
fig_2_plot_dat_hormone <- rbind(mass_hatch_hormone, 
                                mass_juv_hormone, 
                                mass_adult_hormone)%>% 
  dplyr::rename(Treatment_group = hormone)
# temperature data
fig_2_plot_dat_temp <- rbind(mass_hatch_temp, 
                             mass_juv_temp, 
                             mass_adult_temp) %>% 
  dplyr::rename(Treatment_group = temp)
# final plot data
figure_2_plot_dat <- rbind(fig_2_plot_dat_hormone, 
                           fig_2_plot_dat_temp) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

###########
#### !!Figure 2 plot
# reorder data
figure_2_plot_dat$Treatment_group <- factor(figure_2_plot_dat$Treatment_group, 
                                            levels = c("control", "low", "high",
                                                       '23', '28'))
figure_2_plot_dat$Age <- factor(figure_2_plot_dat$Age, 
                                levels = c("hatchling", 
                                           "juvenile", 
                                           "adult"))
###############
# 2A: hatchling plot - Temp
hatch_temp_data <- figure_2_plot_dat %>% filter(Age == 'hatchling',
                                                Treatment == 'Temperature')
figure_2A <- ggplot() +
  geom_errorbar(data = hatch_temp_data, 
                aes(x = Treatment_group, 
                    ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_point(data = hatch_temp_data, aes(x = Treatment_group, 
                                         y = emmean, 
                                         fill = Treatment_group), 
             shape = 22, color = "black", size = 5) +
  geom_text(data = hatch_temp_data, aes(x = Treatment_group, y = emmean + 
                                          SE + .00005, label = paste("n =", n)), 
            vjust = -0.5, size = 4) +
  labs(y = "Mass (g)", 
       x = NULL,
       fill = 'Treatment') +
  theme_classic() +
  scale_y_continuous(limits = c(0.112, 0.1285), breaks = seq(0.112, 0.128, by = 0.004)) +
  scale_fill_manual(values = c("dodgerblue", "tomato2")) +
  theme(axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_line(),  # Ensure ticks are kept
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14, face = "bold"),
        legend.position = c(.01, .45), # Position legend in the top left corner
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "transparent", color = NA), # Make legend background transparent
        axis.text.x = element_blank()) 
figure_2A

###############
# hatchling plot - hormone
hatch_hormone_data <- figure_2_plot_dat %>% filter(Age == 'hatchling',
                                                   Treatment == 'Hormone')
figure_2B <- ggplot() +
  geom_errorbar(data = hatch_hormone_data, 
                aes(x = Treatment_group, 
                    ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  
  geom_point(data = hatch_hormone_data, aes(x = Treatment_group, 
                                            y = emmean, 
                                            fill = Treatment_group), 
             shape = 22, color = "black", size = 5) +
  geom_text(data = hatch_hormone_data, aes(x = Treatment_group, y = emmean + 
                                             SE + .00005, label = paste("n =", n)), 
            vjust = -0.5, size = 4) +
  labs(y = "Mass (g)", 
       x = NULL,
       fill = 'Treatment') +
  theme_classic() +
  scale_y_continuous(limits = c(0.112, 0.1285), breaks = seq(0.112, 0.128, by = 0.004)) +
  scale_fill_manual(values = c("white", "lightgreen", "darkgreen")) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_line(),  # Ensure ticks are kept
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14, face = "bold"),
        legend.position = c(.01, .45), # Position legend in the top left corner
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "transparent", color = NA), # Make legend background transparent
        axis.text.x = element_blank())  # Adjust the anchor point of the legend

###############
# 2C: juvenile plot - Temp
juv_temp_data <- figure_2_plot_dat %>% filter(Age == 'juvenile',
                                              Treatment == 'Temperature')
figure_2C <- ggplot() +
  geom_errorbar(data = juv_temp_data, 
                aes(x = Treatment_group, 
                    ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_point(data = juv_temp_data, aes(x = Treatment_group, 
                                       y = emmean, 
                                       fill = Treatment_group), 
             shape = 22, color = "black", size = 5) +
  geom_text(data = juv_temp_data, aes(x = Treatment_group, y = emmean + 
                                        SE + .00005, label = paste("n =", n)), 
            vjust = -0.5, size = 4) +
  labs(y = "Mass (g)", 
       x = NULL,
       fill = 'Treatment') +
  theme_classic() +
  scale_y_continuous(limits = c(0.38, .55), breaks = seq(0.4, 0.55, by = 0.05)) +
  scale_fill_manual(values = c("dodgerblue", "tomato2")) +
  theme(axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_line(),  # Ensure ticks are kept
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14, face = "bold"),
        legend.position = 'none', # Position legend in the top left corner
        legend.justification = c(0, 1),
        axis.text.x = element_blank()) # Adjust the anchor point of the legend
figure_2C


###############
# juvenile plot - hormone
juv_hormone_data <- figure_2_plot_dat %>% filter(Age == 'juvenile',
                                                 Treatment == 'Hormone')
figure_2D <- ggplot() +
  geom_errorbar(data = juv_hormone_data, 
                aes(x = Treatment_group, 
                    ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_point(data = juv_hormone_data, aes(x = Treatment_group, 
                                          y = emmean, 
                                          fill = Treatment_group), 
             shape = 22, color = "black", size = 5) +
  labs(y = "Mass (g)", 
       x = NULL,
       fill = 'Treatment') +
  geom_text(data = juv_hormone_data, aes(x = Treatment_group, y = emmean + 
                                           SE + .00005, label = paste("n =", n)), 
            vjust = -0.5, size = 4) +
  theme_classic() +
  scale_y_continuous(limits = c(0.38, .55), breaks = seq(0.4, 0.55, by = 0.05)) +
  scale_fill_manual(values = c("white", "lightgreen", "darkgreen")) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_line(),  # Ensure ticks are kept
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14, face = "bold"),
        legend.position = 'none', # Position legend in the top left corner
        legend.justification = c(0, 1),
        axis.text.x = element_blank()) # Adjust the anchor point of the legend
figure_2D

###############
# 2E: adult plot - Temp
adult_temp_data <- figure_2_plot_dat %>% filter(Age == 'adult',
                                                Treatment == 'Temperature')
figure_2E <- ggplot() +
  geom_errorbar(data = adult_temp_data, 
                aes(x = Treatment_group, 
                    ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_point(data = adult_temp_data, aes(x = Treatment_group, 
                                         y = emmean, 
                                         fill = Treatment_group), 
             shape = 22, color = "black", size = 5) +
  geom_text(data = adult_temp_data, aes(x = Treatment_group, y = emmean + 
                                          SE + .00005, label = paste("n =", n)), 
            vjust = -0.5, size = 4) +
  labs(y = "Mass (g)", 
       x = NULL,
       fill = 'Treatment') +
  theme_classic() +
  scale_y_continuous(limits = c(1.4, 1.62), breaks = seq(1.4, 1.6, by = 0.05)) +
  scale_fill_manual(values = c("dodgerblue", "tomato2")) +
  theme(axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_line(),  # Ensure ticks are kept
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14, face = "bold"),
        legend.position = 'none', # Position legend in the top left corner
        legend.justification = c(0, 1),
        axis.text.x = element_text(size = 14, face = "bold")) # Adjust the anchor point of the legend
figure_2E


###############
# adult plot - hormone
adult_hormone_data <- figure_2_plot_dat %>% filter(Age == 'adult',
                                                   Treatment == 'Hormone')
figure_2F <- ggplot() +
  geom_errorbar(data = adult_hormone_data, 
                aes(x = Treatment_group, 
                    ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_point(data = adult_hormone_data, aes(x = Treatment_group, 
                                            y = emmean, 
                                            fill = Treatment_group), 
             shape = 22, color = "black", size = 5) +
  geom_text(data = adult_hormone_data, aes(x = Treatment_group, y = emmean + 
                                             SE + .00005, label = paste("n =", n)), 
            vjust = -0.5, size = 4) +
  labs(y = "Mass (g)", 
       x = NULL,
       fill = 'Treatment') +
  theme_classic() +
  scale_y_continuous(limits = c(1.4, 1.62), breaks = seq(1.4, 1.6, by = 0.05)) +
  scale_fill_manual(values = c("white", "lightgreen", "darkgreen")) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),  # Remove the x-axis title
        axis.ticks.x = element_line(),  # Ensure ticks are kept
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14, face = "bold"),
        legend.position = 'none',  # Position legend in the top left corner
        legend.justification = c(0, 1),
        axis.text.x = element_text(size = 14, face = "bold")) # Adjust the anchor point of the legend
figure_2F
##################
# Final figure
figure_2AB <- plot_grid(figure_2A, figure_2B, labels = c("A"), label_size = 16, label_fontface = "bold", ncol = 2)
figure_2CD <- plot_grid(figure_2C, figure_2D, labels = c("B"), label_size = 16, label_fontface = "bold", ncol = 2)
figure_2EF <- plot_grid(figure_2E, figure_2F, labels = c("C"), label_size = 16, label_fontface = "bold", ncol = 2)
# Combine all the rows into one figure
Figure_2 <- plot_grid(figure_2AB, figure_2CD, figure_2EF, ncol = 1, align = 'v')
# Display the combined plot
print(Figure_2)




##################################################
######### FIUGRE 2 SVL example
SVL_hatch_mod <- lm(hatch_svl_mm ~ temp + hormone, data = data_final)

### Data from model for mean and wiskers 'cort_development_mod'  output
# hormone data - hatchling
SVL_hatch_mod_emm_hormone <- emmeans(SVL_hatch_mod, pairwise ~ hormone)
SVL_hatch_dat_emm_hormone <- as.data.frame(SVL_hatch_mod_emm_hormone$emmeans) %>% 
  mutate(Age = "hatchling",
         Treatment = "Hormone")
# sample size per group
SVL_hatch_dat_sum_hormone <- data_final %>%
  dplyr::select(Lizard_ID, hormone, hatch_svl_mm) %>% 
  dplyr::filter(!is.na(hatch_svl_mm)) %>%
  group_by(hormone) %>%
  summarise(n = n())
# df for plotting
SVL_hatch_hormone <- left_join(SVL_hatch_dat_emm_hormone, 
                                SVL_hatch_dat_sum_hormone,
                                by = "hormone")


# temp data  - hatchling
SVL_hatch_mod_emm_temp <- emmeans(SVL_hatch_mod, pairwise ~ temp)
SVL_hatch_dat_emm_temp <- as.data.frame(SVL_hatch_mod_emm_temp$emmeans) %>% 
  mutate(Age = "hatchling",
         Treatment = "Temperature")
# sample size per group
SVL_hatch_dat_sum_temp <- data_final %>%
  dplyr::select(Lizard_ID, temp, hatch_svl_mm) %>% 
  dplyr::filter(!is.na(hatch_svl_mm)) %>%
  group_by(temp) %>%
  summarise(n = n())
# df for plotting
SVL_hatch_temp<- left_join(SVL_hatch_dat_emm_temp, 
                            SVL_hatch_dat_sum_temp,
                            by = "temp") 

######
### JUV1 model 
SVL_Juv1_mod <- lm(juv1_SVL_mm ~ temp + hormone + scale(Juvenile1_Age), data = data_final)

# Hormone - Juvenile
SVL_Juv1_mod_emm_hormone <- emmeans(SVL_Juv1_mod, pairwise ~ hormone)
SVL_Juv1_emm_hormone <- as.data.frame(SVL_Juv1_mod_emm_hormone$emmeans) %>% 
  mutate(Age = "juvenile",
         Treatment = "Hormone")
# sample size per group
SVL_Juv1_dat_sum_hormone <- data_final %>%
  dplyr::select(Lizard_ID, hormone,temp,juv1_SVL_mm) %>% 
  dplyr::filter(!is.na(juv1_SVL_mm)) %>%
  group_by(hormone) %>%
  summarise(n = n())
# df for plotting
SVL_juv_hormone <- left_join(SVL_Juv1_emm_hormone, 
                              SVL_Juv1_dat_sum_hormone,
                              by = "hormone")

# temp data  - juvenile
SVL_Juv1_mod_emm_temp <- emmeans(SVL_Juv1_mod, pairwise ~ temp)
SVL_Juv1_dat_emm_temp <- as.data.frame(SVL_Juv1_mod_emm_temp$emmeans) %>% 
  mutate(Age = "juvenile",
         Treatment = "Temperature")
# sample size per group
SVL_Juv1_dat_sum_temp <- data_final %>%
  dplyr::select(Lizard_ID, temp, juv1_SVL_mm) %>% 
  dplyr::filter(!is.na(juv1_SVL_mm)) %>%
  group_by(temp) %>%
  summarise(n = n())
# df for plotting
SVL_juv_temp<- left_join(SVL_Juv1_dat_emm_temp, 
                          SVL_Juv1_dat_sum_temp,
                          by = "temp")


######
### SVL: Adult mod 
SVL_Juv3_mod <- lm(juv3_SVL_mm ~ hormone + temp + sex + 
                      scale(Juvenile3_Age), data = data_final)

# Data from model for mean and wiskers 'cort_development_mod'  output
SVL_Juv3_mod_emm_hormone <- emmeans(SVL_Juv3_mod, pairwise ~ hormone)
SVL_Juv3_emm_hormone <- as.data.frame(SVL_Juv3_mod_emm_hormone$emmeans) %>% 
  mutate(Age = "adult",
         Treatment = "Hormone")
# sample size per group
SVL_Juv3_dat_sum_hormone <- data_final %>%
  dplyr::select(Lizard_ID, hormone,temp,juv3_SVL_mm) %>% 
  dplyr::filter(!is.na(juv3_SVL_mm)) %>%
  group_by(hormone) %>%
  summarise(n = n())
# df for plotting
SVL_adult_hormone <- left_join(SVL_Juv3_emm_hormone, 
                                SVL_Juv3_dat_sum_hormone,
                                by = "hormone")

# temp data  - Adult
SVL_Juv3_mod_emm_temp <- emmeans(SVL_Juv3_mod, pairwise ~ temp)
SVL_Juv3_dat_emm_temp <- as.data.frame(SVL_Juv3_mod_emm_temp$emmeans) %>% 
  mutate(Age = "adult",
         Treatment = "Temperature")
# sample size per group
SVL_Juv3_dat_sum_temp <- data_final %>%
  dplyr::select(Lizard_ID, temp, juv3_SVL_mm) %>% 
  dplyr::filter(!is.na(juv3_SVL_mm)) %>%
  group_by(temp) %>%
  summarise(n = n())
# df for plotting
SVL_adult_temp<- left_join(SVL_Juv3_dat_emm_temp, 
                            SVL_Juv3_dat_sum_temp,
                            by = "temp")

#####
### final df for plotting
# hormone data
fig_2_plot_dat_hormone <- rbind(SVL_hatch_hormone, 
                                SVL_juv_hormone, 
                                SVL_adult_hormone)%>% 
  dplyr::rename(Treatment_group = hormone)
# temperature data
fig_2_plot_dat_temp <- rbind(SVL_hatch_temp, 
                             SVL_juv_temp, 
                             SVL_adult_temp) %>% 
  dplyr::rename(Treatment_group = temp)
# final plot data
figure_2_plot_dat_SVL <- rbind(fig_2_plot_dat_hormone, 
                           fig_2_plot_dat_temp) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

###########
#### !!Figure 2 plot SVL
# reorder data
figure_2_plot_dat_SVL$Treatment_group <- factor(figure_2_plot_dat_SVL$Treatment_group, 
                                            levels = c("control", "low", "high",
                                                       '23', '28'))
figure_2_plot_dat_SVL$Age <- factor(figure_2_plot_dat_SVL$Age, 
                                levels = c("hatchling", 
                                           "juvenile", 
                                           "adult"))

###############
# 2A: hatchling plot - Temp
hatch_temp_data_SVL <- figure_2_plot_dat_SVL %>% filter(Age == 'hatchling',
                                                Treatment == 'Temperature')
figure_2A <- ggplot() +
  geom_errorbar(data = hatch_temp_data_SVL, 
                aes(x = Treatment_group, 
                    ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_point(data = hatch_temp_data_SVL, aes(x = Treatment_group, 
                                             y = emmean, 
                                             fill = Treatment_group), 
             shape = 22, color = "black", size = 5) +
  geom_text(data = hatch_temp_data_SVL, aes(x = Treatment_group, y = emmean + 
                                          SE + .00005, label = paste("n =", n)), 
            vjust = -0.5, size = 4) +
  labs(y = "SVL (mm)", 
       x = NULL,
       fill = 'Treatment') +
  theme_classic() +
  scale_y_continuous(limits = c(16.70, 17.8), breaks = seq(16.75, 17.75, by = .25)) +
  scale_fill_manual(values = c("dodgerblue", "tomato2")) +
  theme(axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_line(),  # Ensure ticks are kept
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14, face = "bold"),
        legend.position = c(.01, .45), # Position legend in the top left corner
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "transparent", color = NA), # Make legend background transparent
        axis.text.x = element_blank()) 
figure_2A

###############
# hatchling plot - hormone
hatch_hormone_data_SVL <- figure_2_plot_dat_SVL %>% filter(Age == 'hatchling',
                                                   Treatment == 'Hormone')
figure_2B <- ggplot() +
  geom_errorbar(data = hatch_hormone_data_SVL, 
                aes(x = Treatment_group, 
                    ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_point(data = hatch_hormone_data_SVL, aes(x = Treatment_group, 
                                                y = emmean, 
                                                fill = Treatment_group), 
             shape = 22, color = "black", size = 5) +
  geom_text(data = hatch_hormone_data_SVL, aes(x = Treatment_group, y = emmean + 
                                             SE + .00005, label = paste("n =", n)), 
            vjust = -0.5, size = 4) +
  labs(y = "SVL (mm)", 
       x = NULL,
       fill = 'Treatment') +
  theme_classic() +
  scale_y_continuous(limits = c(16.70, 17.8), breaks = seq(16.75, 17.75, by = .25)) +
  scale_fill_manual(values = c("white", "lightgreen", "darkgreen")) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_line(),  # Ensure ticks are kept
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14, face = "bold"),
        legend.position = c(.01, .45), # Position legend in the top left corner
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "transparent", color = NA), # Make legend background transparent
        axis.text.x = element_blank())  # Adjust the anchor point of the legend
figure_2B

###############
# 2C: juvenile plot - Temp
juv_temp_data_SVL <- figure_2_plot_dat_SVL %>% filter(Age == 'juvenile',
                                              Treatment == 'Temperature')
figure_2C <- ggplot() +
  geom_errorbar(data = juv_temp_data_SVL, 
                aes(x = Treatment_group, 
                    ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_point(data = juv_temp_data_SVL, aes(x = Treatment_group, 
                                           y = emmean, 
                                           fill = Treatment_group), 
             shape = 22, color = "black", size = 5) +
  geom_text(data = juv_temp_data_SVL, aes(x = Treatment_group, y = emmean + 
                                        SE + .00005, label = paste("n =", n)), 
            vjust = -0.5, size = 4) +
  labs(y = "SVL (mm)", 
       x = NULL,
       fill = 'Treatment') +
  theme_classic() +
  scale_y_continuous(limits = c(27.48, 31.4), 
                     breaks = seq(28.00, 31.00, by = 1),
                     labels = label_number(accuracy = 0.01))+
  scale_fill_manual(values = c("dodgerblue", "tomato2")) +
  theme(axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_line(),  # Ensure ticks are kept
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14, face = "bold"),
        legend.position = 'none', # Position legend in the top left corner
        legend.justification = c(0, 1),
        axis.text.x = element_blank()) # Adjust the anchor point of the legend
figure_2C


###############
# juvenile plot - hormone
juv_hormone_data_SVL <- figure_2_plot_dat_SVL %>% filter(Age == 'juvenile',
                                                 Treatment == 'Hormone')
figure_2D <- ggplot() +
  geom_errorbar(data = juv_hormone_data_SVL, 
                aes(x = Treatment_group, 
                    ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_point(data = juv_hormone_data_SVL, aes(x = Treatment_group, 
                                              y = emmean, 
                                              fill = Treatment_group), 
             shape = 22, color = "black", size = 5) +
  labs(y = "SVL (mm)", 
       x = NULL,
       fill = 'Treatment') +
  geom_text(data = juv_hormone_data_SVL, aes(x = Treatment_group, y = emmean + 
                                           SE + .00005, label = paste("n =", n)), 
            vjust = -0.5, size = 4) +
  theme_classic() +
  scale_y_continuous(limits = c(27.48, 31.4), breaks = seq(28.00, 31.00, by = 1)) +
  scale_fill_manual(values = c("white", "lightgreen", "darkgreen")) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_line(),  # Ensure ticks are kept
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14, face = "bold"),
        legend.position = 'none', # Position legend in the top left corner
        legend.justification = c(0, 1),
        axis.text.x = element_blank()) # Adjust the anchor point of the legend
figure_2D

###############
# 2E: adult plot - Temp
adult_temp_data_SVL <- figure_2_plot_dat_SVL %>% filter(Age == 'adult',
                                                Treatment == 'Temperature')
figure_2E <- ggplot() +
  geom_errorbar(data = adult_temp_data_SVL, 
                aes(x = Treatment_group, 
                    ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_point(data = adult_temp_data_SVL, aes(x = Treatment_group, 
                                             y = emmean, 
                                             fill = Treatment_group), 
             shape = 22, color = "black", size = 5) +
  geom_text(data = adult_temp_data_SVL, aes(x = Treatment_group, y = emmean + 
                                          SE + .00005, label = paste("n =", n)), 
            vjust = -0.5, size = 4) +
  labs(y = "SVL (mm)", 
       x = NULL,
       fill = 'Treatment') +
  theme_classic() +
  scale_y_continuous(limits = c(40.5, 42.75), breaks = seq(40.5, 42.5, by = 0.5),
                     labels = label_number(accuracy = 0.01)) +
  scale_fill_manual(values = c("dodgerblue", "tomato2")) +
  theme(axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_line(),  # Ensure ticks are kept
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14, face = "bold"),
        legend.position = 'none', # Position legend in the top left corner
        legend.justification = c(0, 1),
        axis.text.x = element_text(size = 14, face = "bold")) # Adjust the anchor point of the legend
figure_2E


###############
# adult plot - hormone
adult_hormone_data_SVL <- figure_2_plot_dat_SVL %>% filter(Age == 'adult',
                                                   Treatment == 'Hormone')
figure_2F <- ggplot() +
  geom_errorbar(data = adult_hormone_data_SVL, 
                aes(x = Treatment_group, 
                    ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_point(data = adult_hormone_data_SVL, aes(x = Treatment_group, 
                                                y = emmean, 
                                                fill = Treatment_group), 
             shape = 22, color = "black", size = 5) +
  geom_text(data = adult_hormone_data_SVL, aes(x = Treatment_group, y = emmean + 
                                             SE + .00005, label = paste("n =", n)), 
            vjust = -0.5, size = 4) +
  labs(y = "SVL (mm)", 
       x = NULL,
       fill = 'Treatment') +
  theme_classic() +
  scale_y_continuous(limits = c(40.5, 42.75), breaks = seq(40.5, 42.5, by = 0.5)) +
  scale_fill_manual(values = c("white", "lightgreen", "darkgreen")) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(),  # Remove the x-axis title
        axis.ticks.x = element_line(),  # Ensure ticks are kept
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14, face = "bold"),
        legend.position = 'none',  # Position legend in the top left corner
        legend.justification = c(0, 1),
        axis.text.x = element_text(size = 14, face = "bold")) # Adjust the anchor point of the legend
figure_2F
##################
# Final figure
figure_2AB_SVL <- plot_grid(figure_2A, figure_2B, labels = c("A"), label_size = 16, label_fontface = "bold", ncol = 2)
figure_2CD_SVL <- plot_grid(figure_2C, figure_2D, labels = c("B"), label_size = 16, label_fontface = "bold", ncol = 2)
figure_2EF_SVL <- plot_grid(figure_2E, figure_2F, labels = c("C"), label_size = 16, label_fontface = "bold", ncol = 2)
# Combine all the rows into one figure
Figure_2_SVL <- plot_grid(figure_2AB_SVL, figure_2CD_SVL, figure_2EF_SVL, ncol = 1, align = 'v')
# Display the combined plot
print(Figure_2_SVL)


###################
####### Figure 3
# Reorder the factor levels
cort_dat$hormone <- factor(cort_dat$hormone, levels = c("control", "low", "high"))

# Data from model for mean and wiskers 'cort_development_mod'  output
sum.dat.emm <- as.data.frame(cort_development_mod_emm$emmeans) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 2)))
# sample size per group
sum.dat <- cort_dat %>%
  group_by(hormone) %>%
  summarise(n = n())
# combine for plot
fig.3.emm.dat <- left_join(sum.dat.emm, sum.dat, by = "hormone")


# Plot
Figure_3 <- ggplot() +
  geom_errorbar(data = fig.3.emm.dat, aes(x = hormone, ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_text(data = fig.3.emm.dat, aes(x = hormone, y = emmean + SE + .01, label = paste("n =", n)), 
            vjust = -0.5, size = 4) +
  geom_point(data = fig.3.emm.dat, aes(x = hormone, y = emmean, fill = hormone), 
             shape = 22, color = "black", size = 5) +
  labs(y = "Log CORT", 
       x = "Corticosterone treatment") +
  theme_classic() +
  scale_fill_manual(values = c("white", "lightgreen", "darkgreen")) +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14),
        legend.position = c(0.02, 0.98), # Position legend in the top left corner
        legend.justification = c(0, 1)) # Adjust the anchor point of the legend
Figure_3




##############################################################
####### !!!KRIS Figure 4 model
mito_hormone$sex <- as.factor(mito_hormone$sex)

# Recode the levels of the 'sex' factor
Fig4_dat <- mito_hormone %>%
  mutate(sex = recode_factor(sex, "M" = "Male", "F" = "Female"))

######### Fig 4 A - BASAL
# Fit a linear model for each sex separately
Basal_male_mod <- lm(basal_corrected_pmol ~ juv3_mass_g, data = Fig4_dat %>% filter(sex == "Male"))
Basal_female_mod <- lm(basal_corrected_pmol ~ juv3_mass_g, data = Fig4_dat %>% filter(sex == "Female"))
slope_male_basal <- coef(Basal_male_mod)[2]
slope_female_basal <- coef(Basal_female_mod)[2]

# Plot
Fig_4A <- ggplot(Fig4_dat, aes(x = juv3_mass_g, y = basal_corrected_pmol)) +
  geom_point(aes(color = sex, shape = sex), size = 3.5) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  scale_shape_manual(values = c("Male" = 16, "Female" = 17)) +
  scale_color_manual(values = c("Male" = "steelblue", "Female" = 'darkred')) +
  annotate("text", x = 2.0, 
           y =.04, 
           label = paste0("Sex: P ", Fig_4A_pvalue, 
                          "\nSlope (Male): ", round(slope_male_basal, 2),
                          "\nSlope (Female): ", round(slope_female_basal, 2)),  
           size = 4, fontface = "bold") +
  theme_classic() +
  labs(x = "Mass (g)",
       y = "Basal Respirometry pmol O2/(sec*µg)",
       color = "Sex",
       shape = "Sex") +
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_line(), 
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14, face = "bold"),
        legend.position = c(0.05, 0.98),  
        legend.justification = c(0, 1))

Fig_4A

######
# Fig 4B ADP - OXPHOS 
# Fit a linear model for each sex separately
ADP_male_mod <- lm(adp_corrected_pmol ~ juv3_mass_g, data = Fig4_dat %>% filter(sex == "Male"))
ADP_female_mod <- lm(adp_corrected_pmol ~ juv3_mass_g, data = Fig4_dat %>% filter(sex == "Female"))
slope_male_adp <- coef(ADP_male_mod)[2]
slope_female_adp <- coef(ADP_female_mod)[2]

# Plot
Fig_4B <- ggplot(Fig4_dat, aes(x = juv3_mass_g, y = adp_corrected_pmol)) +
  geom_point(aes(color = sex, shape = sex), size = 3.5) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  scale_shape_manual(values = c("Male" = 16, "Female" = 17)) +
  scale_color_manual(values = c("Male" = "steelblue", "Female" = 'darkred')) +
  theme_classic() +
  annotate("text", x = 2.0, 
           y =.3, 
           label = paste0("Sex: P = ", Fig_4B_pvalue_adp, 
                          "\nSlope (Male): ", round(slope_male_adp, 2),
                          "\nSlope (Female): ", round(slope_female_adp, 2)),  
           size = 4, fontface = "bold") +
  theme_classic() +
  labs(x = "Mass (g)",
       y = "OXPHOS pmol O2/(sec*µg)",
       color = "Sex",
       shape = "Sex") +
  scale_y_continuous(limits = c(0, 6.5), 
                     breaks = seq(0, 6, by = 2),
                     labels = scales::label_number(accuracy = 0.1)) +
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_line(),  
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14, face = "bold"),
        legend.position = 'none',  
        legend.justification = c(0, 1))

Fig_4B

########################
###### Fig 4C Oligo
# Fit a linear model for each sex separately
Oligo_male_mod <- lm(oligo_corrected_pmol ~ juv3_mass_g, data = Fig4_dat %>% filter(sex == "Male"))
Oligo_female_mod <- lm(oligo_corrected_pmol ~ juv3_mass_g, data = Fig4_dat %>% filter(sex == "Female"))
slope_male_oligo <- coef(Oligo_male_mod)[2]
slope_female_oligo <- coef(Oligo_female_mod)[2]

# Plot
Fig_4C <- ggplot(Fig4_dat, aes(x = juv3_mass_g, y = oligo_corrected_pmol)) +
  geom_point(aes(color = sex, shape = sex), size = 3.5) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  scale_shape_manual(values = c("Male" = 16, "Female" = 17)) +
  scale_color_manual(values = c("Male" = 'steelblue', "Female" = 'darkred')) +
  theme_classic() +
  labs(x = "Mass (g)",
       y = "Leak pmol O2/(sec*µg)",
       color = "Sex",
       shape = "Sex") +
  annotate("text", x = 2.04, 
           y =.09, 
           label = paste0("Sex: P ", Fig_4C_pvalue_oligo, 
                          "\nSlope (Male): ", round(slope_male_oligo, 2),
                          "\nSlope (Female): ", round(slope_female_oligo, 2)),  
           size = 4, fontface = "bold") +
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        axis.title.x = element_text(size = 14), 
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14, face = "bold"),
        legend.position = 'none',  
        legend.justification = c(0, 1))

Fig_4C

########################
# Fig 4D - RCR
# Fit a linear model for each sex separately
RCR_male_mod <- lm(RCR ~ juv3_mass_g, data = Fig4_dat %>% filter(sex == "Male"))
RCR_female_mod <- lm(RCR ~ juv3_mass_g, data = Fig4_dat %>% filter(sex == "Female"))
slope_male_rcr <- coef(RCR_male_mod)[2]
slope_female_rcr <- coef(RCR_female_mod)[2]

# Plot
Fig_4D <- ggplot(Fig4_dat, aes(x = juv3_mass_g, y = RCR)) +
  geom_point(aes(color = sex, shape = sex), size = 3.5) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  scale_shape_manual(values = c("Male" = 16, "Female" = 17)) +
  scale_color_manual(values = c("Male" = "steelblue", "Female" = 'darkred')) +
  theme_classic() +
  labs(x = "Mass (g)",
       y = "RCR",
       color = "Sex",
       shape = "Sex") +
  annotate("text", x = 1.84, 
           y =4,  
           label = paste0("Sex: P = ", Fig_4D_pvalue_rcr, 
                          "\nSlope (Male): ", round(slope_male_rcr, 2),
                          "\nSlope (Female): ", round(slope_female_rcr, 2)),  
           size = 4, fontface = "bold") +
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        axis.title.x = element_text(size = 14), 
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14, face = "bold"),
        legend.position = 'none',  
        legend.justification = c(0, 1))

Fig_4D

################################
####### Figure 4 Final
Fig_4A <- Fig_4A + ggtitle("A") + theme(plot.title = element_text(hjust = -0.1, vjust = 1.5, face = 'bold'), 
                                        plot.margin = margin(15, 12, 20, 5.5, "pt"))
Fig_4B <- Fig_4B + ggtitle("B") + theme(plot.title = element_text(hjust = -0.1, vjust = 1.5, face = 'bold'), 
                                        plot.margin = margin(15, 12, 20, 5.5, "pt"))
Fig_4C <- Fig_4C + ggtitle("C") + theme(plot.title = element_text(hjust = -0.1, vjust = 1.5, face = 'bold'), 
                                        plot.margin = margin(15, 12, 20, 5.5, "pt"))
Fig_4D <- Fig_4D + ggtitle("D") + theme(plot.title = element_text(hjust = -0.1, vjust = 1.5, face = 'bold'), 
                                        plot.margin = margin(15, 12, 20, 5.5, "pt"))

# Combine the adjusted plots into a single figure
Figure_4 <- plot_grid(
  Fig_4A,
  Fig_4B,
  Fig_4C,
  Fig_4D,
  ncol = 2, nrow = 2, labels = NULL
)

Figure_4






######################
#### Figure 5 !!!Kris  
# Perform ANOVA
# Fit separate linear models for each sex to extract the slope coefficients
Mass_basal_male_mod <- lm(hatch_juv3_MASS_growth * 1000 ~ log(CORT_Final_Hormone_ng_mL), data = Fig5_dat %>% filter(sex == "Male"))
Mass_basal_female_mod <- lm(hatch_juv3_MASS_growth * 1000 ~ log(CORT_Final_Hormone_ng_mL), data = Fig5_dat %>% filter(sex == "Female"))
slope_male_mass_growth <- coef(Mass_basal_male_mod)[2]
slope_female_mass_growth <- coef(Mass_basal_female_mod)[2]

# Plot
Figure_5 <- ggplot(Fig5_dat, 
                   aes(x = log(CORT_Final_Hormone_ng_mL), 
                       y = hatch_juv3_MASS_growth * 1000)) +
  geom_point(aes(color = sex, shape = sex), size = 4) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  scale_shape_manual(values = c("Male" = 16, "Female" = 17)) +
  scale_color_manual(values = c("Male" = "steelblue", "Female" = 'darkred')) +
  theme_classic() +
  labs(x = "Log(CORT Final Hormone)",
       y = "Mass Growth Rate (mg/d)",
       color = "Sex",
       shape = "Sex") +
  annotate("text", x = 5.92,
           y = 2.3, 
           label = paste0("Sex: P =", Fig_5_pvalue_mass_growth, 
                          "\nSlope (Male): ", round(slope_male_mass_growth, 2),
                          "\nSlope (Female): ", round(slope_female_mass_growth, 2)),  
           size = 4, fontface = "bold") +
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        axis.title.x = element_text(size = 14), 
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14, face = "bold"),
        legend.position = c(.1, 0.98),   
        legend.justification = c(0, 1))
Figure_5

#########
# Final figures
Figure_2 #1100x900
Figure_2_SVL #1100x900
Figure_4 #1300x1000
Figure_5





















