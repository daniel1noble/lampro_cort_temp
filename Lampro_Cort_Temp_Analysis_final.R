################# Lampro Cort Temp Analysis
### 1. first dataset (Yolk_hormone.csv) is from the the egg hormone and analysis follows
#
### 2. second datasets (Cort_hatch_juv1_2_OC.csv & Cort_juv_3_np.csv) are from hatching, juv, adult 
# (morphometrics and bioenergetics) and analysis follows
### 3. all figures in MS follow analysis in order of MS
#

# packages
pacman::p_load(dplyr, tidyverse, ggpubr, lme4, emmeans, here, plotly, performance, car,
               scales, Hmisc, ppcor, ordinal, cowplot, lubridate, Hmisc, ppcor, ordinal,
               lmerTest)


#######################
# Egg hormone analysis
CORT <- read.csv ('data/Yolk_hormone.csv') %>%
  filter(set != c("1")) %>%
  filter(Sample.ID !=c("294"))
###egg 294 excluded because it has dissection note that says:
#"yolk sample dropped, found on bench one min later, a little dehydrated, exclude if hormone level is strange"
#and the value was strange

# data set up
 CORT$treatment <- as.factor(CORT$treatment)
     CORT$Plate <- as.factor(CORT$Plate)
CORT$CORT_value <- as.numeric(CORT$finalCORT.pg.mg.)
 CORT$treatment <- factor(CORT$treatment, 
                         levels = c("C_Topical", "CORT_5pg_Topical", "CORT_10pg_Topical"))

##linear model 
yolkB_mod <- lm(log(CORT_value) ~ treatment + Plate, data = CORT)
Anova(yolkB_mod)
summary(yolkB_mod)
emm_LogMod <- emmeans (yolkB_mod, pairwise ~ treatment)
colSums(model.matrix(yolkB_mod))
check_model(yolkB_mod)  
CORT$residuals <- residuals.lm(yolkB_mod, c(type = "pearson"))
shapiro.test(CORT$residuals)
#GAMma model
gam_yolk <- glm (CORT_value ~ treatment + Plate, family = Gamma (link = log), data = CORT)
Anova (gam_yolk)
summary (gam_yolk)
emm_gam <- emmeans (gam_yolk, pairwise ~ treatment)
check_model(gam_yolk)
saveRDS(gam_yolk, 'models/yolk_hormone_mod.RDS')

###Mean raw CORT levels
tapply(CORT$CORT_value, CORT$treatment, mean, na.rm = TRUE)
tapply(CORT$CORT_value, CORT$treatment, sd, na.rm = TRUE)
# quick plot
Boxplot(CORT$CORT_value, CORT$treatment, na.action = na.exclude)
Boxplot(CORT$LogCORT, CORT$treatment, na.action = na.exclude)


###########################################################################
###########################################################################
###########  Hatchling - Adult morphology, hormones, & bioenergetics
###########################################################################
###########################################################################

# 1) bring in data base data that has morphology and mortality data on across life stage: hatch, juv_1, juv_2
dat_hatch_juv1_2 <- read.csv ("data/Cort_hatch_juv.csv")

# 2) bring in data that has final morphology measurement and hormone data: "juv_3"
dat_juv3 <- read.csv(file = "data/Cort_adult.csv") %>% 
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
colnames(merged_data) <- gsub("juv3", "adult", colnames(merged_data))

# final dataframe
data_final <- merged_data %>% 
  rename(hatch_svl_mm = hatch_svl_orig_mm,
         hatch_mass_g = hatch_mass_orig, 
         juv1_SVL_mm = juv1_measure_SVL_mm,
         juv1_mass_g = juv1_measure_mass_g) 

str(data_final)


################################################################################
################# Adding columns for age, growth rate, julian day, body condition
################################################################################

#############
##Calculating age in days post hatching for juvenile 1 and 2 measurements and final adult measurement
####Using Julian dates
##Juvenile 1 measurements all collected on 03/06/2022 = julian date 154
data_final$Juvenile1_Age <- 154 - data_final$Julian_hatch_date
##Juvenile 2 measurements all collected on 10/08/2022 = julian date 222
data_final$Juvenile2_Age <- 222 - data_final$Julian_hatch_date
##adult measurements collected at variable dates; measurement date converted to Julian date and age calculated as difference
##Note that dates span 2022 - 2023
str(data_final$adult_date)
data_final$adult_date <- as.Date(data_final$adult_date, "%d/%m/%y")
data_final$Adult_Age <- (yday(data_final$adult_date)) + (365- data_final$Julian_hatch_date)

############
#### Summary statistics post hatching: mean, SD, and range of age of measurements post hatching
##Juvenile 1 = 105.7, 10.8, 85 - 123
##Juvenile 2 = 173.7, 10.8, 153 - 191
##Adult = 466.1, 12.4, 440 - 491
##note: SD for Juvenile 1 and Juvenile 2 are the same because these measurements were all collected on the same day
mean(data_final$Juvenile1_Age, na.rm = TRUE)
sd(data_final$Juvenile1_Age, na.rm = TRUE)
range(data_final$Juvenile1_Age, na.rm = TRUE)
mean(data_final$Juvenile2_Age, na.rm = TRUE)
sd(data_final$Juvenile2_Age, na.rm = TRUE)
range(data_final$Juvenile2_Age, na.rm = TRUE)
mean(data_final$Adult_Age, na.rm = TRUE)
sd(data_final$Adult_Age, na.rm = TRUE)
range(data_final$Adult_Age, na.rm = TRUE)

############
#### growth rate calculations: calculated growth rate over lifespan 
# dividing change body size (SVL and mass) by the total number of days elapsed
data_final <- data_final %>% 
  mutate(hatch_adult_SVL_growth =  (adult_SVL_mm - hatch_svl_mm)/Adult_Age,
         hatch_adult_MASS_growth = (adult_mass_g - hatch_mass_g)/Adult_Age, 
         hatch_juv1_SVL_growth = (juv1_SVL_mm - hatch_svl_mm) / Juvenile1_Age, 
         hatch_juv1_mass_growth = (juv1_mass_g - hatch_mass_g) / Juvenile1_Age, 
         juv1_adult_svl_growth = (adult_SVL_mm - juv1_SVL_mm) / (Adult_Age - Juvenile1_Age),
         juv1_adult_mass_growth = (adult_mass_g - juv1_mass_g)/ (Adult_Age - Juvenile1_Age))

###checking variables
###Some variables need to be changed to factors so they are not treated as linear variables in analyses
data_final$temp <- as.factor(data_final$temp)
data_final$hormone <- as.factor(data_final$hormone)
data_final$clutch <- as.factor(data_final$clutch)
data_final$Plate_CORT_adult <- as.factor(data_final$Plate_CORT_adult)
data_final$adult_T4_plate <- as.factor(data_final$adult_T4_plate)
data_final$adult_Testosterone_plate <- as.factor(data_final$adult_Testosterone_plate)
data_final$adult_HandlingTime_sec <- as.numeric(data_final$adult_HandlingTime_sec)
data_final$adult_liver_time_sec <- as.numeric(data_final$adult_liver_time_sec)


###############
####Body condition calculations
###residuals and scaled mass index calculated and produce equivalent statistical results
#Scaled mass calculations calculated separately for each age (hatching, juvenile, and juvenile 3 )
##Calculates residual condition as mass divided by SVL
data_final <- data_final %>% 
  mutate(hatch_mass_svl =  hatch_mass_g / hatch_svl_mm, 
         juv1_mass_svl =  juv1_mass_g /  juv1_SVL_mm,
         adult_mass_svl = adult_mass_g / adult_SVL_mm,
         log_hatch_mass = log(hatch_mass_g),
         log_hatch_SVL = log(hatch_svl_mm), 
         log_juv1_mass = log(juv1_mass_g),
         log_juv1_SVL = log(juv1_SVL_mm),
         log_adult_mass = log(adult_mass_g), 
         log_adult_SVL = log(adult_SVL_mm))

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

###Scaled mass for Adult ; creates new variable 'ScaledMass_adult' 
##Pearsons correlation coefficient hatching
PCor_adult <- cor(data_final$log_adult_mass, data_final$log_adult_SVL, method = "pearson", use = "complete.obs")
mod_adult <- lm(log_adult_mass ~ log_adult_SVL, data = data_final)
##OLS coefficient from summary (2.4744)
summary(mod_adult)
##Calculates mean tarsus for scaled mass equation
MeanSVL_adult <- mean(data_final$adult_SVL_mm, na.rm = TRUE)
SMCof_adult <- 2.4744/PCor_adult
##calculates scaled mass; plots scaled mass to residual mass
data_final$ScaledMass_adult <- data_final$adult_mass_g * (MeanSVL_adult /data_final$adult_SVL_mm) ^ (SMCof_adult)
plot (ScaledMass_adult ~ adult_mass_svl, data = data_final)
ols_scaled_mass_mod <- lm(ScaledMass_adult  ~ adult_mass_svl, data = data_final)
summary(ols_scaled_mass_mod)
abline(ols_scaled_mass_mod)




#################################################################################
#################################################################################
###### ###### ###### ######  ANALYSIS   ###### ###### ###### ###### ###### ###### 
#################################################################################
#################################################################################

interaction_table  <- data.frame("Variable" = NA, "Delta_AIC" = NA)

################
#### 1.	What are the effects of developmental treatments (temp, cort, interaction) on time to hatch?
#### effect of temperature on incubation days - faster hatch warmer temps 
################
# temperature effects days to hatch; no effect on hormone or temp x hormone interaction
#hatch_dev_mod <- lm(days_to_hatch ~ temp + hormone, data = data_final)
hatch_dev_mod_int <- lmer(days_to_hatch ~ temp + hormone + temp*hormone + (1 | clutch), 
                          data = data_final)
hatch_dev_mod <- lmer(days_to_hatch ~ temp + hormone + (1 | clutch), data = data_final)
# Update the models to use Maximum Likelihood
hatch_dev_mod_int_ml <- update(hatch_dev_mod_int, REML = FALSE)
hatch_dev_mod_ml <- update(hatch_dev_mod, REML = FALSE)
AIC(hatch_dev_mod_int_ml, hatch_dev_mod_ml) # go with hatch_dev_mod_ml
interaction_table[1,1]  <- "Incubation time (days)"
interaction_table[1,2] <- AIC(hatch_dev_mod_int_ml) - AIC(hatch_dev_mod_ml)
Anova(hatch_dev_mod)
summary(hatch_dev_mod)

# faster hatch warmer temps
hatch_dev_mod_emm <- emmeans(hatch_dev_mod, pairwise ~ temp)
plot(hatch_dev_mod_emm)
saveRDS(hatch_dev_mod, "models/hatch_dev_mod.RDS")

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
#SVL_hatch_mod <- lm(hatch_svl_mm ~ temp + hormone, data = data_final)
SVL_hatch_mod <- lmer(hatch_svl_mm ~ temp + hormone + (1 | clutch), data = data_final)
SVL_hatch_mod_ml <- update(SVL_hatch_mod, REML = FALSE)
SVL_hatch_mod_int <- lmer(hatch_svl_mm ~ temp + hormone +temp*hormone + (1 | clutch), 
                          data = data_final)
SVL_hatch_mod_int_ml <- update(SVL_hatch_mod_int, REML = FALSE)
AIC(SVL_hatch_mod_ml, SVL_hatch_mod_int_ml) # go SVL_hatch_mod_ml 
interaction_table[2,1]  <- "Hatchling Snout-vent Length (SVL) (mm)"
interaction_table[2,2] <- AIC(SVL_hatch_mod_int_ml) - AIC(SVL_hatch_mod_ml)
Anova(SVL_hatch_mod)
# high hormone lower SVL
SVL_hatch_mod_emm <- emmeans(SVL_hatch_mod, pairwise ~ hormone)
plot(SVL_hatch_mod_emm)
saveRDS(SVL_hatch_mod, "models/SVL_hatch_mod.RDS")


#### MASS at hatching
####temp and hormone treatment affect mass; no interaction
#Mass_hatch_mod <- lm(hatch_mass_g ~ temp + hormone, data = data_final)
# Initial mixed-effects model with random effect
Mass_hatch_mod <- lmer(hatch_mass_g ~ temp + hormone + (1 | clutch), data = data_final)
Mass_hatch_mod_ml <- update(Mass_hatch_mod, REML = FALSE)
# Mixed-effects model with interaction term
Mass_hatch_mod_int <- lmer(hatch_mass_g ~ temp + hormone + temp*hormone + (1 | clutch), 
                           data = data_final)
Mass_hatch_mod_int_ml <- update(Mass_hatch_mod_int, REML = FALSE)
# AIC comparison
AIC(Mass_hatch_mod_ml, Mass_hatch_mod_int_ml) # go with Mass_hatch_mod_ml
interaction_table[3,1]  <- "Hatchling Mass (g)"
interaction_table[3,2] <- AIC(Mass_hatch_mod_int_ml) - AIC(Mass_hatch_mod_ml)
# model checks
check_model(Mass_hatch_mod)
Anova(Mass_hatch_mod)
# hormone pairwise comparison - high CORT lower mass compared to control
Mass_hormone_hatch_mod_emm <- emmeans(Mass_hatch_mod, pairwise ~ hormone)
plot(Mass_hormone_hatch_mod_emm)
# temp pairwise comparison - cooler temps higher mass
Mass_temp_hatch_mod_emm <- emmeans(Mass_hatch_mod, pairwise ~ temp)
plot(Mass_temp_hatch_mod_emm)
saveRDS(Mass_hatch_mod, "models/MASS_hatch_mod.RDS")


####Scaled mass at hatching
#scaled_mass_hatch_mod <- lm(ScaledMass_hatch ~ temp + hormone, data = data_final)
# Initial mixed-effects model with random effect
scaled_mass_hatch_mod <- lmer(ScaledMass_hatch ~ temp + hormone + (1 | clutch), data = data_final)
scaled_mass_hatch_mod_ml <- update(scaled_mass_hatch_mod, REML = FALSE)

# Mixed-effects model with interaction term
scaled_mass_hatch_mod_int <- lmer(ScaledMass_hatch ~ temp + hormone + temp*hormone + (1 | clutch), data = data_final)
scaled_mass_hatch_mod_int_ml <- update(scaled_mass_hatch_mod_int, REML = FALSE)

AIC(scaled_mass_hatch_mod_ml, scaled_mass_hatch_mod_int_ml)

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
#### adult_MASS: low mass males
#### adult_SVL: high hormone low body size; males smaller than females
#### adult_scaled mass: Sex effect: male lower BCI than females - could be size difference

###############
# SVL ANALYSIS: Juv_1, adult (sex accounted for last measurement)
###############
#### JUV_1 - SVL: 28 temp has smaller SVL
#SVL_Juv1_mod <- lm(juv1_SVL_mm ~ temp + hormone + scale(Juvenile1_Age), data = data_final)
# Initial mixed-effects model with random effect
SVL_Juv1_mod <- lmer(juv1_SVL_mm ~ temp + hormone + scale(Juvenile1_Age) + (1 | clutch), data = data_final)
SVL_Juv1_mod_ml <- update(SVL_Juv1_mod, REML = FALSE)

#Mixed-effects model with interaction term
SVL_Juv1_mod_int <- lmer(juv1_SVL_mm ~ temp + hormone + scale(Juvenile1_Age) + temp*hormone + (1 | clutch), data = data_final)
SVL_Juv1_mod_int_ml <- update(SVL_Juv1_mod_int, REML = FALSE)

AIC(SVL_Juv1_mod_ml, SVL_Juv1_mod_int)
interaction_table[4,1]  <- "Juvenile SVL (mm)"
interaction_table[4,2] <- AIC(SVL_Juv1_mod_int_ml) - AIC(SVL_Juv1_mod_ml)
# effect on temperature and day of hatch on body size 
Anova(SVL_Juv1_mod)
summary(SVL_Juv1_mod)
# plot - cooler temps larger svl
SVL_Juv1_mod_emm <- emmeans(SVL_Juv1_mod, pairwise ~ temp)
plot(SVL_Juv1_mod_emm)
saveRDS(SVL_Juv1_mod, "models/SVL_Juv1_mod.RDS")

#### adult  SVL -hormone effect
# SVL_adult_mod <- lm(adult_SVL_mm ~ temp + hormone + sex + scale(Adult_Age), data = data_final)
# Initial mixed-effects model with random effect
SVL_adult_mod <- lmer(adult_SVL_mm ~ temp + hormone + sex + scale(Adult_Age) + (1 | clutch), data = data_final)
SVL_adult_mod_ml <- update(SVL_adult_mod, REML = FALSE)

# Mixed-effects model with interaction term
SVL_adult_mod_int <- lmer(adult_SVL_mm ~ temp + hormone + sex + scale(Adult_Age) + temp*hormone + (1 | clutch), data = data_final)
SVL_adult_mod_int_ml <- update(SVL_adult_mod_int, REML = FALSE)

# Compare models using AIC
AIC(SVL_adult_mod_ml, SVL_adult_mod_int_ml)
interaction_table[5,1]  <- "Adult SVL (mm)"
interaction_table[5,2] <- AIC(SVL_adult_mod_int_ml) - AIC(SVL_adult_mod_ml)

check_model(SVL_adult_mod)
# effect of hormones and sex
Anova(SVL_adult_mod)
summary(SVL_adult_mod)
# plot - high hormone lower svl than control
SVL_adult_hormone_mod_emm <- emmeans(SVL_adult_mod, pairwise ~ hormone)
plot(SVL_adult_hormone_mod_emm)
# plot - males are smaller than females
SVL_adult_sex_mod_emm <- emmeans(SVL_adult_mod, pairwise ~ sex)
plot(SVL_adult_sex_mod_emm)
saveRDS(SVL_adult_mod, "models/SVL_adult_mod.RDS")







###############
# MASS ANALYSIS: Juv_1, adult
###############
#### MASS: Juv_1 - effects of hormone and temp treatment 
#Mass_Juv1_mod <- lm(juv1_mass_g ~ temp + hormone + scale(Juvenile1_Age), data = data_final)
# Initial mixed-effects model with random effect
Mass_Juv1_mod <- lmer(juv1_mass_g ~ temp + hormone + scale(Juvenile1_Age) + (1 | clutch), data = data_final)
Mass_Juv1_mod_ml <- update(Mass_Juv1_mod, REML = FALSE)

# Mixed-effects model with interaction term
Mass_Juv1_mod_int <- lmer(juv1_mass_g ~ temp + hormone + scale(Juvenile1_Age) + temp*hormone + (1 | clutch), data = data_final)
Mass_Juv1_mod_int_ml <- update(Mass_Juv1_mod_int, REML = FALSE)

AIC(Mass_Juv1_mod_ml, Mass_Juv1_mod_int)
interaction_table[6,1]  <- "Juvenile Mass (g)"
interaction_table[6,2] <- AIC(Mass_Juv1_mod_int_ml) - AIC(Mass_Juv1_mod_ml)

check_model(Mass_Juv1_mod)
# summary
summary(Mass_Juv1_mod)
Anova(Mass_Juv1_mod)
# hormone plot - high CORT lower mass than control 
Mass_Juv1_hormone_mod_emm <- emmeans(Mass_Juv1_mod, pairwise ~ hormone)
plot(Mass_Juv1_hormone_mod_emm)
# temp plot - 28 temp has lower mass
Mass_Juv1_temp_mod_emm <- emmeans(Mass_Juv1_mod, pairwise ~ temp)
plot(Mass_Juv1_temp_mod_emm)
saveRDS(Mass_Juv1_mod, "models/Mass_Juv1_mod.RDS")

#### MASS: adult sex effect
# Mass_adult_mod <- lm(adult_mass_g ~ hormone + temp + sex + scale(Adult_Age), data = data_final)
# Initial mixed-effects model with random effect
Mass_adult_mod <- lmer(adult_mass_g ~ hormone + temp + sex + scale(Adult_Age) + (1 | clutch), data = data_final)
Mass_adult_mod_ml <- update(Mass_adult_mod, REML = FALSE)

# Mixed-effects model with interaction term
Mass_adult_mod_int <- lmer(adult_mass_g ~ hormone + temp + sex + scale(Adult_Age) + temp*hormone + (1 | clutch), data = data_final)
Mass_adult_mod_int_ml <- update(Mass_adult_mod_int, REML = FALSE)

# Compare models using AIC
AIC(Mass_adult_mod_ml, Mass_adult_mod_int_ml)
interaction_table[7,1]  <- "Adult Mass (g)"
interaction_table[7,2] <- AIC(Mass_adult_mod_int_ml) - AIC(Mass_adult_mod_ml)

check_model(Mass_adult_mod)
Anova(Mass_adult_mod)
summary(Mass_adult_mod)
# temp plot - males smaller than females
Sex_adult_temp_mod_emm <- emmeans(Mass_adult_mod, pairwise ~ hormone)
plot(Sex_adult_temp_mod_emm)
saveRDS(Mass_adult_mod, "models/Mass_adult_mod.RDS")




###############
# Scaled mass  ANALYSIS: Juv_1, adult
###############
#### Scaled mass : Juv_1 - effects of hormone and temp treatment 
####Scaled mass 
#scaled_ScaledMass_juv1_mod <- lm(ScaledMass_juv1 ~ temp + hormone, data = data_final)
# Initial mixed-effects model with random effect
scaledScaledMass_juv1_mod <- lmer(ScaledMass_juv1 ~ temp + hormone + (1 | clutch), data = data_final)
scaledScaledMass_juv1_mod_ml <- update(scaledScaledMass_juv1_mod, REML = FALSE)

# Mixed-effects model with interaction term
scaledScaledMass_juv1_mod_int <- lmer(ScaledMass_hatch ~ temp + hormone + temp*hormone + (1 | clutch), data = data_final)
scaledScaledMass_juv1_mod_int_ml <- update(scaledScaledMass_juv1_mod_int, REML = FALSE)

AIC(scaled_mass_hatch_mod_ml, scaled_mass_hatch_mod_int_ml)

summary(scaledScaledMass_juv1_mod)
Anova(scaledScaledMass_juv1_mod)


#### Scaled mass : Juv_1 - effects of hormone and temp treatment 
####Scaled mass 
#scaled_ScaledMass_juv1_mod <- lm(ScaledMass_adult ~ temp + hormone, data = data_final)
# Initial mixed-effects model with random effect
ScaledMass_adult_mod <- lmer(ScaledMass_adult ~ temp + hormone + (1 | clutch), data = data_final)
ScaledMass_adult_mod_ml <- update(ScaledMass_adult_mod, REML = FALSE)

# Mixed-effects model with interaction term
ScaledMass_adult_mod_int <- lmer(ScaledMass_hatch ~ temp + hormone + temp*hormone + (1 | clutch), data = data_final)
ScaledMass_adult_mod_int_ml <- update(ScaledMass_adult_mod_int, REML = FALSE)

AIC(scaled_mass_hatch_mod_ml, scaled_mass_hatch_mod_int_ml)
summary(ScaledMass_adult_mod)
Anova(ScaledMass_adult_mod)



# save final data with BCI residuals for figures and quarto doc
data_final <-  mutate(data_final, hormone = factor(hormone, 
                                                   levels = c("control","low","high"))) %>% 
  filter(!is.na(hormone)) %>%
  group_by(hormone)


################
#### 4.	developmental treatment; Testing days to mortality across treatments
#### no differences between treatments. Mortality of duration of study
################
# survival fisher test- filter out missing animals
survival<- data_final %>% 
  filter(liz_status_db != "MISSING")

# Create contingency tables
survival_table_temp <- table(survival$liz_status_db, survival$temp)
survival_table_hormone <- table(survival$liz_status_db, survival$hormone)

# Perform Fisher's Exact Test
survival_day_temp <- chisq.test(survival_table_temp)
survival_day_hormone <- chisq.test(survival_table_hormone)

# save analysis
saveRDS(survival_day_hormone, "models/survival_day_hormone.RDS")
saveRDS(survival_day_temp, "models/survival_day_temp.RDS")



################
#### 5.	Growthrate
####hatching to juvenile 1 measurements
# Initial mixed-effects model with random effect
hatch_juv1_SVL_growth_mod <- lmer(hatch_juv1_SVL_growth ~ hormone + temp + (1 | clutch), data = data_final)
hatch_juv1_SVL_growth_mod_ml <- update(hatch_juv1_SVL_growth_mod, REML = FALSE)
# Mixed-effects model with interaction term
hatch_juv1_SVL_growth_mod_int <- lmer(hatch_juv1_SVL_growth ~ hormone + temp + temp*hormone + (1 | clutch), data = data_final)
hatch_juv1_SVL_growth_mod_int_ml <- update(hatch_juv1_SVL_growth_mod_int, REML = FALSE)
# Compare models using AIC
AIC(hatch_juv1_SVL_growth_mod_ml, hatch_juv1_SVL_growth_mod_int_ml)
interaction_table[8,1]  <- "Hatchling to Juvenile Growth SVL (mm)"
interaction_table[8,2] <- AIC(hatch_juv1_SVL_growth_mod_int_ml) - AIC(hatch_juv1_SVL_growth_mod_ml)

Anova(hatch_juv1_SVL_growth_mod)
summary(hatch_juv1_SVL_growth_mod)


####### mass growth hatch to Juv1
# hatch_juv1_mass_growth_mod <- lm(hatch_juv1_mass_growth ~ hormone + temp, data = data_final)
# Initial mixed-effects model with random effect
hatch_juv1_mass_growth_mod <- lmer(hatch_juv1_mass_growth ~ hormone + temp + (1 | clutch), data = data_final)
hatch_juv1_mass_growth_mod_ml <- update(hatch_juv1_mass_growth_mod, REML = FALSE)
# Mixed-effects model with interaction term
hatch_juv1_mass_growth_mod_int <- lmer(hatch_juv1_mass_growth ~ hormone + temp + temp*hormone + (1 | clutch), data = data_final)
hatch_juv1_mass_growth_mod_int_ml <- update(hatch_juv1_mass_growth_mod_int, REML = FALSE)
# Compare models using AIC
AIC(hatch_juv1_mass_growth_mod_ml, hatch_juv1_mass_growth_mod_int_ml)
interaction_table[9,1]  <- "Hatchling to Juvenile Growth Mass (g)"
interaction_table[9,2] <- AIC(hatch_juv1_mass_growth_mod_int_ml) - AIC(hatch_juv1_mass_growth_mod_ml)

Anova(hatch_juv1_mass_growth_mod)

summary(hatch_juv1_mass_growth_mod)
hatch_juv1_mass_growth_mod_emm <- emmeans(hatch_juv1_mass_growth_mod, pairwise ~ hormone)
plot(hatch_juv1_mass_growth_mod_emm)
tapply(data_final$hatch_juv1_SVL_growth, data_final$temp, mean, na.rm = TRUE)
tapply(data_final$hatch_juv1_SVL_growth, data_final$temp, sd, na.rm = TRUE)
tapply(data_final$hatch_juv1_mass_growth, data_final$temp, mean, na.rm = TRUE)
tapply(data_final$hatch_juv1_mass_growth, data_final$temp, sd, na.rm = TRUE)





################
#### Growth from hatch to adulthood
# SVL Growth
# hatch_adult_SVL_growth_mod <- lm(hatch_adult_SVL_growth ~ hormone + temp + sex, data = data_final)
# Initial mixed-effects model with random effect
hatch_adult_SVL_growth_mod <- lmer(hatch_adult_SVL_growth ~ hormone + temp + sex + (1 | clutch), data = data_final)
hatch_adult_SVL_growth_mod_ml <- update(hatch_adult_SVL_growth_mod, REML = FALSE)
# Mixed-effects model with interaction term
hatch_adult_SVL_growth_mod_int <- lmer(hatch_adult_SVL_growth ~ hormone + temp + sex + temp*hormone + (1 | clutch), data = data_final)
hatch_adult_SVL_growth_mod_int_ml <- update(hatch_adult_SVL_growth_mod_int, REML = FALSE)
# Compare models using AIC
AIC(hatch_adult_SVL_growth_mod_ml, hatch_adult_SVL_growth_mod_int_ml)
interaction_table[10,1]  <- "Hatchling to Adult Growth SVL (mm)"
interaction_table[10,2] <- AIC(hatch_adult_SVL_growth_mod_int_ml) - AIC(hatch_adult_SVL_growth_mod_ml)

Anova(hatch_adult_SVL_growth_mod)
summary(hatch_adult_SVL_growth_mod)
check_model(hatch_adult_SVL_growth_mod)
hatch_adult_SVL_growth_mod_emm <- emmeans(hatch_adult_SVL_growth_mod, pairwise ~ temp)
plot(hatch_adult_SVL_growth_mod_emm)


######## mass growth : hatch to adult
# hatch_adult_MASS_growth_mod <- lm(hatch_adult_MASS_growth ~ hormone + temp + sex, data = data_final)
# Initial mixed-effects model with random effect
hatch_adult_MASS_growth_mod <- lmer(hatch_adult_MASS_growth ~ hormone + temp + sex + (1 | clutch), data = data_final)
hatch_adult_MASS_growth_mod_ml <- update(hatch_adult_MASS_growth_mod, REML = FALSE)
# Mixed-effects model with interaction term
hatch_adult_MASS_growth_mod_int <- lmer(hatch_adult_MASS_growth ~ hormone + temp + sex + temp*hormone + (1 | clutch), data = data_final)
hatch_adult_MASS_growth_mod_int_ml <- update(hatch_adult_MASS_growth_mod_int, REML = FALSE)
# Compare models using AIC
AIC(hatch_adult_MASS_growth_mod_ml, hatch_adult_MASS_growth_mod_int_ml)
interaction_table[11,1]  <- "Hatchling to Adult Growth Mass (g)"
interaction_table[11,2] <- AIC(hatch_adult_MASS_growth_mod_int_ml) - AIC(hatch_adult_MASS_growth_mod_ml)

Anova(hatch_adult_MASS_growth_mod)
summary(hatch_adult_MASS_growth_mod)





########################
### 6: CORT, T4, testosterone
### remove values of 4ul and below
########
# CORT
# data
cort_dat <- data_final %>%
  drop_na(adult_CORT_Final_Hormone_ng_mL, temp, hormone, adult_HandlingTime_sec) %>% 
  filter(adult_Sample_volume_ul > 4 ) # remove values below 4 ul juv3_Sample_volume_ul # remove values below 4 ul adult_Sample_volume_ul
# plot
plot(adult_CORT_Final_Hormone_ng_mL ~ hormone, data = cort_dat)
plot(adult_CORT_Final_Hormone_ng_mL ~ temp, data = cort_dat)

###Checking plate ID and handling time for effects on CORT; handling time = ns; plate ID is p=.007
cort_factors <- lmer(log(adult_CORT_Final_Hormone_ng_mL) ~ adult_HandlingTime_sec + Plate_CORT_adult +(1|clutch), data = cort_dat)
Anova(cort_factors)
summary(cort_factors)


######## 1) treatments- how does cort vary across temp and hormone treatment; handling time not used because ns in model of CORT factors
#cort_development_mod <- lm(log(adult_CORT_Final_Hormone_ng_mL) ~ hormone + temp +  sex + scale(Adult_Age) +Plate_CORT_adult, data = cort_dat)

# model used in MS where 
cort_development_mod <- lmer(log(adult_CORT_Final_Hormone_ng_mL) ~ hormone + temp +  sex + scale(Adult_Age) + Plate_CORT_adult + (1 | clutch), data = cort_dat)
check_model(cort_development_mod)
Anova(cort_development_mod)
summary(cort_development_mod)
##hormone
cort_development_mod_emm <- emmeans(cort_development_mod, pairwise ~ hormone)
plot(cort_development_mod_emm)

# sex
cort_development_sex_mod_emm <- emmeans(cort_development_mod, pairwise ~ sex)
plot(cort_development_sex_mod_emm)
#temperature
cort_temp_mod_emm <- emmeans(cort_development_mod, pairwise ~ temp)
plot(cort_temp_mod_emm)
saveRDS(cort_temp_mod_emm, 'models/cort_development_mod.rds')

###Visualize baseline CORT levels
Boxplot(cort_dat$adult_CORT_Final_Hormone_ng_mL, cort_dat$hormone, na.action = na.exclude)

##calculate residuals of plate ID and CORT levels
cort_residuals <- lm(adult_CORT_Final_Hormone_ng_mL ~ Plate_CORT_adult, data = cort_dat, na.action=na.exclude) # fit the model for residuals
cort_dat$cort_residuals <- residuals(cort_residuals, na.action=na.exclude) # Save the residual values


#### Body size posthoc model for discussion point where body mass is included
cort_development_mod_ondi <- lmer(log(adult_CORT_Final_Hormone_ng_mL) ~ hormone + temp + adult_mass_g +  sex + scale(Adult_Age) + Plate_CORT_adult + (1 | clutch), data = cort_dat)
Anova(cort_development_mod_ondi)




####Body size and condition at juvenile 3 related to CORT levels; 'handling time not included because it did not affect baseline CORT levels as above)
##CORT and scaled mass at juvenile 3; treatments not included because they have no effect on body condition at this time point
cort_body_condition_mod <- lm(ScaledMass_adult ~ cort_residuals + sex + scale(Adult_Age) + sex, data = cort_dat)
Anova(cort_body_condition_mod)
summary(cort_body_condition_mod)
check_model(cort_body_condition_mod)

###SVL with CORT; adult SVL positive association with CORT
cort_svl_juvenile3_mod <- lm(adult_SVL_mm ~ cort_residuals + sex + scale(Adult_Age), data = cort_dat)
Anova(cort_svl_juvenile3_mod)
summary(cort_svl_juvenile3_mod)
check_model(cort_svl_juvenile3_mod)

###Mass at juvenile 3 tested against CORT levels; CORT positive association with adult mass
cort_mass_juvenile3_mod <- lm(adult_mass_g ~ cort_residuals + sex + scale(Adult_Age), data = cort_dat)
Anova(cort_mass_juvenile3_mod)
summary(cort_mass_juvenile3_mod)
check_model(cort_mass_juvenile3_mod)



################
# T4
################ 
t4_dat <- data_final %>%
  drop_na(adult_T4_corrected_ng_mL, temp, hormone, Lizard_ID) %>% 
  filter(adult_Sample_volume_ul > 4 ) # remove values below 4 ul adult_Sample_volume_ul

###Checking plate ID and handling time for effects on T4; handling time = ns; plate ID is NS
T4_factors <- lmer((log(adult_T4_corrected_ng_mL)) ~ adult_HandlingTime_sec + adult_T4_plate + (1|clutch), data = t4_dat)
Anova(T4_factors)
summary(T4_factors)


# 1) treatments - 
T4_temp_hormone_sex_mod <- lmer(log(adult_T4_corrected_ng_mL) ~ hormone + temp + sex + scale(Adult_Age) + (1|clutch), data = t4_dat)
check_model(T4_temp_hormone_sex_mod)
Anova(T4_temp_hormone_sex_mod)
summary(T4_temp_hormone_sex_mod)
# sex differences: females have higher T4
T4_temp_hormone_sex_mod_emm <- emmeans(T4_temp_hormone_sex_mod, pairwise ~ sex)
plot(T4_temp_hormone_sex_mod_emm)
saveRDS(T4_temp_hormone_sex_mod, "models/T4_temp_hormone_sex_mod.RDS")

# 2) SVL and associations with T4 and cort
T4_SVL_mod <- lmer(adult_SVL_mm  ~ log(adult_T4_corrected_ng_mL) + sex + scale(Adult_Age) + (1|clutch), data = t4_dat)
check_model(T4_SVL_mod )
summary(T4_SVL_mod )
Anova(T4_SVL_mod )

# 3) MASS and associations between T4
T4_mass_growth_mod <- lm(adult_mass_g  ~ log(adult_T4_corrected_ng_mL) + sex + scale(Adult_Age), 
                         data = t4_dat)
check_model(T4_mass_growth_mod)
Anova(T4_mass_growth_mod)
summary(T4_mass_growth_mod)
saveRDS(T4_mass_growth_mod, "models/T4_mass_growth_mod.RDS")


########################
### Testosterone
test_dat <- data_final %>%
  drop_na(adult_Testosterone_Final_ng_ml) %>% 
  filter(adult_Sample_volume_ul > 4 ) # remove values below 4 ul adult_Sample_volume_ul
test_dat$rank_test <- rank(test_dat$adult_Testosterone_Final_ng_ml, ties.method = "first")

###Look at the effects of plate and handling time on testosterone levels; handling time affects T levels (p=0.04); plate is NS
Test_factor <- lmer(log(adult_Testosterone_Final_ng_ml) ~ adult_Testosterone_plate+ adult_HandlingTime_sec + (1|clutch), data = test_dat)
Anova(Test_factor)
summary(Test_factor)

###Get residuals from model of T levels and handling time and save as new variable
test_residuals <- lm((log(adult_Testosterone_Final_ng_ml)) ~ adult_HandlingTime_sec, data = test_dat)
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
cor.test(test_dat$adult_Testosterone_Final_ng_ml, test_dat$adult_CORT_Final_Hormone_ng_mL, method = c("spearman"))
plot(test_dat$adult_CORT_Final_Hormone_ng_mL, test_dat$adult_Testosterone_Final_ng_ml)
cort_T_mod <- lm(adult_Testosterone_Final_ng_ml ~ adult_CORT_Final_Hormone_ng_mL, data = test_dat)
abline(cort_T_mod)



########################
### 7: Mito function
########################
# bring in data, rename, and save
###if looking at relationships with hormone levels - use cort_dat or equivalent to exclude samples with 
##low volumes
mito_dat <- data_final %>%
  dplyr::select(c("Lizard_ID", "clutch", "temp", "hormone", "Adult_Age", "sex","hatch_adult_MASS_growth", "hatch_adult_SVL_growth",
                  "adult_inject_time_sec", "adult_liver_time_sec", "adult_oroboros",
                  "adult_T4_plate", "adult_HandlingTime_sec", "adult_T4_corrected_ng_mL", 
                  "adult_CORT_Final_Hormone_ng_mL", "adult_basal_corrected.pmol..sec.ng..",
                  "adult_chamber","adult_basal_corrected.pmol..sec.ng..", 
                  "adult_adp_corrected.pmol..sec.ng..", "adult_oligo_corrected.pmol..sec.ng..", 
                  "adult_fccp_corrected.pmol..sec.ng..", "adult_RCR.L.R.", "adult_RCR.R.ETS.",
                  "adult_oroboros_comments", "adult_RCR.L.ETS.", "adult_date", "Plate_CORT_adult",
                  "adult_Sample_volume_ul", "hatch_svl_mm", "adult_Testosterone_Final_ng_ml",
                  "adult_mass_g", "adult_SVL_mm",  "ScaledMass_adult",
                  "adult_HandlingTime_sec",  "adult_chamber")) %>% 
  dplyr::rename(inject_time_sec = adult_inject_time_sec, 
                chamber = adult_chamber,
                adult_basal_corrected = adult_basal_corrected.pmol..sec.ng..,
                basal_corrected_pmol = adult_basal_corrected.pmol..sec.ng..,
                oligo_corrected_pmol = adult_oligo_corrected.pmol..sec.ng..,
                fccp_corrected_pmol = adult_fccp_corrected.pmol..sec.ng..,
                adp_corrected_pmol =  adult_adp_corrected.pmol..sec.ng..,
                T4_plate_ID = adult_T4_plate,
                HandlingTime_sec = adult_HandlingTime_sec,
                T4_corrected_ng_mL = adult_T4_corrected_ng_mL, 
                CORT_Final_Hormone_ng_mL = adult_CORT_Final_Hormone_ng_mL,
                RCR_L_R = adult_RCR.L.R., 
                RCR_R_ETS = adult_RCR.R.ETS.,
                RCR_L_ETS = adult_RCR.L.ETS., 
                oroboros_comments = adult_oroboros_comments) %>% 
  mutate(ID_and_Comments = paste(Lizard_ID, oroboros_comments, sep = ": "),
         OXPHOS = fccp_corrected_pmol - oligo_corrected_pmol,
         RCR = adp_corrected_pmol/oligo_corrected_pmol) %>% 
  
  filter(Lizard_ID != c("LD784_21")) ####LD784_21 had extremely high values and lab notes suggest should be filtered from analyses 

mito_dat$cham_oro <- with(mito_dat, interaction(adult_oroboros, chamber))

##################### SANITY CHECKS BEFORE ANALYSIS #####################
###Checking factors that could affect oxygen consumption measurements - handling times, Oroboros, and chamber
###'chamber' is nested within Oroboros; there are three measurements of time (time to inject the drug, time to euthanize the animal
###'and time to collect liver, have opted to only use time to collect liver because it will incorporate all other time measurements and 
###'likely all time measures are correlated and so shouldn't all be going in the model anyway)
###all factors NS
#RCR - all NS
RCR_factors <- lm(RCR ~ adult_liver_time_sec + cham_oro, data = mito_dat)
Anova(RCR_factors)
summary(RCR_factors)
hist(mito_dat$RCR)
###Basal; all factors NS
basal_factors <- lm(basal_corrected_pmol ~ adult_liver_time_sec + adult_oroboros/chamber, data = mito_dat)
Anova(basal_factors)
summary (basal_factors)
hist(mito_dat$basal_corrected_pmol)
##ADP - all factors NS
ADP_factors <- lm(adp_corrected_pmol ~ adult_liver_time_sec + adult_oroboros/chamber, data = mito_dat)
Anova(ADP_factors)
summary (ADP_factors)
##Oligo - all factors NS but time to collect liver p = 0.06
Oligo_factors <- lm (oligo_corrected_pmol ~ adult_liver_time_sec + adult_oroboros/chamber, data = mito_dat)
Anova(Oligo_factors)
summary (Oligo_factors)
hist(mito_dat$oligo_corrected_pmol)
###########################################################

##################################################
#####Effects of treatments on mitochondrial function 
mito_dat <- mito_dat %>% 
  mutate(mass_basal =  basal_corrected_pmol / adult_mass_g,
         mass_ADP = adp_corrected_pmol / adult_mass_g,
         mass_oligo = oligo_corrected_pmol / adult_mass_g) %>%
  mutate(mass_RCR = mass_ADP / mass_oligo) %>% 
  rename(adult_age = Adult_Age)



###Basal - interaction removed ns and removed; sex p = 0.005 males higher than females
Basal_treatment_mod_int <- lmer(basal_corrected_pmol ~ temp*hormone + scale(adult_age) + sex + scale(adult_mass_g) + (1|clutch), data = mito_dat)
Basal_treatment_mod <- lmer(basal_corrected_pmol ~ temp + hormone + scale(adult_age) + sex + scale(adult_mass_g) + (1|clutch), data = mito_dat)
Anova(Basal_treatment_mod)
#interaction check
Basal_treatment_mod_int_ml <- update(Basal_treatment_mod_int, REML = FALSE)
Basal_treatment_mod_ml <- update(Basal_treatment_mod, REML = FALSE)
interaction_table[12,1]  <- "Basal respiration (pmol/sec/ng)"
interaction_table[12,2] <- AIC(Basal_treatment_mod_int_ml) - AIC(Basal_treatment_mod_ml)
# check model
check_model(Basal_treatment_mod)
Basal_treatment_mod_emm <- emmeans(Basal_treatment_mod, pairwise ~ sex)
plot(Basal_treatment_mod_emm)
Boxplot(mito_dat$mass_basal, mito_dat$sex, na.action = na.exclude)
saveRDS(Basal_treatment_mod, "models/Basal_treatment_mod.RDS")

###OXPHOS(ADP) - treatment effects;
ADP_treatment_mod_int <- lmer(adp_corrected_pmol  ~temp*hormone + sex + scale(adult_age) + scale(adult_mass_g)  + (1|clutch), data = mito_dat)
ADP_treatment_mod <- lmer(adp_corrected_pmol  ~temp + hormone + sex + scale(adult_age) + scale(adult_mass_g)  + (1|clutch), data = mito_dat)
Anova(ADP_treatment_mod)
summary(ADP_treatment_mod)
# interaction check
ADP_treatment_mod_int_ml <- update(ADP_treatment_mod_int, REML = FALSE)
ADP_treatment_mod_ml <- update(ADP_treatment_mod, REML = FALSE)
interaction_table[13,1]  <- "OXPHOS respiration (pmol/sec/ng)"
interaction_table[13,2] <- AIC(ADP_treatment_mod_int_ml) - AIC(ADP_treatment_mod_ml)
# check top mod
check_model(ADP_treatment_mod)
ADP_treatment_mod_emm <- emmeans(ADP_treatment_mod, pairwise ~ sex)
plot(ADP_treatment_mod_emm)
saveRDS(ADP_treatment_mod, "models/ADP_treatment_mod.RDS")

##LEAK (OLIGO) mito by treatment
Oligo_treatment_mod_int <- lmer(oligo_corrected_pmol ~ temp*hormone + sex + scale(adult_age) + scale(adult_mass_g)  + (1|clutch), data = mito_dat)
Oligo_treatment_mod <- lmer(oligo_corrected_pmol ~ temp + hormone + sex + scale(adult_age) + scale(adult_mass_g)  + (1|clutch), data = mito_dat)
Anova(Oligo_treatment_mod)
summary(Oligo_treatment_mod)
# check interaction
Oligo_treatment_mod_int_ml <- update(Oligo_treatment_mod_int, REML = FALSE)
Oligo_treatment_mod_ml <- update(Oligo_treatment_mod, REML = FALSE)
interaction_table[14,1]  <- "LEAK respiration (pmol/sec/ng)"
interaction_table[14,2] <- AIC(Oligo_treatment_mod_int_ml) - AIC(Oligo_treatment_mod_ml)
# check top model
check_model(Oligo_treatment_mod)
Oligo_treatment_mod_emm <- emmeans(Oligo_treatment_mod, pairwise ~ sex)
plot(Oligo_treatment_mod_emm)
saveRDS(Oligo_treatment_mod, "models/Oligo_treatment_mod.RDS")

##RCR calculated as state 3 (ADP)/ state 4(olgio)
RCR_mod_int <- lmer(RCR ~ temp*hormone + sex + scale(adult_age) + scale(adult_mass_g)  +(1|clutch), data = mito_dat)
RCR_mod <- lmer(RCR ~ temp + hormone + sex + scale(adult_age) + scale(adult_mass_g)  +(1|clutch), data = mito_dat)
Anova(RCR_mod)
RCR_mod_mod_emm <- emmeans(RCR_mod, pairwise ~ temp)
plot(RCR_mod_mod_emm)
summary (RCR_mod)
# interaction check
RCR_mod_int_ml <- update(RCR_mod_int, REML = FALSE)
RCR_mod_ml <- update(RCR_mod, REML = FALSE)
interaction_table[15,1]  <- "Respiratory Control Ratio (RCR)"
interaction_table[15,2] <- AIC(RCR_mod_int_ml) - AIC(RCR_mod_ml)
# check top model
check_model(RCR_mod)
saveRDS(RCR_mod, "models/RCR_treatment_mod.RDS")

# Write supplemental table comparing models with and without interactions. Note, interaction model AIC always subtracted from AIC of main effects model. Used maximum likelihood estimation for all models.
write.csv(interaction_table, "tables/interaction_table.csv")


###For models with hormone levels - need to filter n=7 data points from 
###samples with less than 7ul 
mito_hormone <- mito_dat %>%
  drop_na(CORT_Final_Hormone_ng_mL, temp, hormone, HandlingTime_sec) %>% 
  filter(adult_Sample_volume_ul > 4 ) # remove values below 4 ul adult_Sample_volume_ul


##CORT and plate ID and handling time; NS
plot(mito_hormone$CORT_Final_Hormone_ng_mL, mito_hormone$HandlingTime_sec)
B_handling <- lm(CORT_Final_Hormone_ng_mL ~ HandlingTime_sec +  Plate_CORT_adult, data = mito_hormone)
Anova(B_handling)
abline(B_handling)
CORT_residuals_mito <- lm(CORT_Final_Hormone_ng_mL ~ Plate_CORT_adult, data = mito_hormone)
Anova(CORT_residuals_mito)


###T4 plate ID and handling time; NS
T4_residuals_mod <- lm(T4_corrected_ng_mL ~ HandlingTime_sec + T4_plate_ID, data = mito_hormone, na.action = na.exclude)
Anova(T4_residuals_mod)
summary(T4_residuals_mod)



###############################################################
###Associations between mitochondrial bioenergetics and hormones levels
mito_hormone <- mito_hormone %>% 
  mutate(mass_basal =  basal_corrected_pmol / adult_mass_g,
         mass_ADP = adp_corrected_pmol / adult_mass_g,
         mass_oligo = oligo_corrected_pmol / adult_mass_g) %>%
  mutate(mass_RCR = mass_ADP / mass_oligo)

##Basal 
Basal_hormone_mod <- lmer(mass_basal ~ log(CORT_Final_Hormone_ng_mL) + log(T4_corrected_ng_mL) + sex + +(1|clutch), data = mito_hormone)
Anova(Basal_hormone_mod)
summary(Basal_hormone_mod)
check_model(Basal_hormone_mod)
Boxplot(mito_hormone$mass_basal, mito_hormone$sex, na.action = na.exclude)
saveRDS(Basal_hormone_mod, 'models/Basal_hormone_mod.RDS')


###ADP -  T4 positive association with ADP mito 
ADP_hormone_mod <- lmer(mass_ADP ~ log(CORT_Final_Hormone_ng_mL) + log(T4_corrected_ng_mL) + sex  +(1|clutch), data = mito_hormone)
Anova(ADP_hormone_mod)
summary(ADP_hormone_mod)
check_model(ADP_hormone_mod)
saveRDS(ADP_hormone_mod, 'models/ADP_hormone_mod.RDS')

##oligo - T4 near positive association
Oligo_hormone_mod <- lmer(mass_oligo ~ log(CORT_Final_Hormone_ng_mL) + log(T4_corrected_ng_mL) + sex + +(1|clutch), data = mito_hormone)
Anova(Oligo_hormone_mod)
summary(Oligo_hormone_mod)
check_model(Oligo_hormone_mod)
saveRDS(Oligo_hormone_mod, 'models/Oligo_hormone_mod.RDS')

#RCR - 
RCR_hormone_mod <- lmer(mass_RCR ~ log(CORT_Final_Hormone_ng_mL) + log(T4_corrected_ng_mL) + sex  +(1|clutch), data = mito_hormone)
Anova(RCR_hormone_mod)
summary (RCR_hormone_mod)
check_model(RCR_hormone_mod)
saveRDS(RCR_hormone_mod, 'models/RCR_hormone_mod.RDS')


##########################################
###Growth and body size and mito function  (cham_oro dropped for all)
###Basal - 
# mass
Mass_basal <- lmer(hatch_adult_MASS_growth ~ basal_corrected_pmol + sex + log(CORT_Final_Hormone_ng_mL) + log(T4_corrected_ng_mL)  +(1|clutch), data = mito_hormone)
Anova(Mass_basal)
summary(Mass_basal)
check_model(Mass_basal)
saveRDS(Mass_basal, "models/Mass_Growth_basal.RDS")
# SVL
SVL_basal <- lmer(hatch_adult_SVL_growth~ basal_corrected_pmol + sex + log(CORT_Final_Hormone_ng_mL) + log(T4_corrected_ng_mL)  +(1|clutch), data = mito_hormone)
Anova(SVL_basal)
check_model(SVL_basal)
saveRDS(SVL_basal, "models/SVL_Growth_basal.RDS")

##OXPHOS (ADP) 
# mass
Mass_Growth_ADP <- lmer(hatch_adult_MASS_growth ~ adp_corrected_pmol + sex + log(CORT_Final_Hormone_ng_mL) + log(T4_corrected_ng_mL)  +(1|clutch), data = mito_hormone)
Anova(Mass_Growth_ADP)
summary(Mass_Growth_ADP)
check_model(Mass_Growth_ADP)
saveRDS(Mass_Growth_ADP, "models/Mass_Growth_ADP.RDS")
#SVL
SVL_Growth_ADP <- lmer(hatch_adult_SVL_growth ~ adp_corrected_pmol + log(CORT_Final_Hormone_ng_mL) + log(T4_corrected_ng_mL) + sex +(1|clutch), data = mito_hormone)
Anova(SVL_Growth_ADP)
summary(SVL_Growth_ADP)
check_model(SVL_Growth_ADP)
saveRDS(SVL_Growth_ADP, "models/SVL_Growth_ADP.RDS")


###LEAK (OLIGO) 
#mass
Mass_Growth_Oligo <- lmer(hatch_adult_MASS_growth ~ oligo_corrected_pmol + log(CORT_Final_Hormone_ng_mL) + log(T4_corrected_ng_mL) + sex+(1|clutch), data = mito_hormone)
Anova(Mass_Growth_Oligo)
summary(Mass_Growth_Oligo)
saveRDS(Mass_Growth_Oligo, "models/Mass_Growth_Oligo.RDS")
#SVL
SVL_Growth_Oligo <- lmer(hatch_adult_SVL_growth ~ oligo_corrected_pmol + sex + log(CORT_Final_Hormone_ng_mL) + log(T4_corrected_ng_mL) +(1|clutch), data = mito_hormone)
Anova(SVL_Growth_Oligo)
summary(SVL_Growth_Oligo)
saveRDS(SVL_Growth_Oligo, "models/SVL_Growth_Oligo.RDS")


###RCR
#mass
Mass_Growth_RCR <- lmer(hatch_adult_MASS_growth~ RCR + sex + log(CORT_Final_Hormone_ng_mL) + log(T4_corrected_ng_mL) + +(1|clutch), data = mito_hormone)
Anova(Mass_Growth_RCR)
summary(Mass_Growth_RCR)
saveRDS(Mass_Growth_RCR, "models/Mass_Growth_RCR.RDS")
#SVL
SVL_Growth_RCR <- lmer(hatch_adult_SVL_growth ~ RCR + sex + log(CORT_Final_Hormone_ng_mL) + log(T4_corrected_ng_mL) +(1|clutch), data = mito_hormone)
Anova(SVL_Growth_RCR)
summary(SVL_Growth_RCR)
saveRDS(SVL_Growth_RCR, "models/SVL_Growth_RCR.RDS")



############# ############# ############# ############# ############# 
############# FINAL FIGURES  #############
# Figure 1 is timeline
###############################################################
####### Figure 2  Egg CORT (old figure 2)
###############################################################
gam_yolk <- readRDS('models/yolk_hormone_mod.RDS')
emm_gam <- emmeans (gam_yolk, pairwise ~ treatment)
# raw data
data_fig_2_raw <- CORT %>% 
  rename(hormone = treatment ) %>% 
  mutate(hormone = case_when(hormone == "C_Topical" ~ "control",
                             hormone == "CORT_5pg_Topical" ~ "low",
                             hormone == "CORT_10pg_Topical" ~ 'high',
                             TRUE ~ hormone))
# Data from model for mean and wiskers 'cort_development_mod'  output
data_fig_2_emm <- as.data.frame(emm_gam$emmeans) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>% 
  rename(hormone = treatment)
# sample size per group
fig.2.sum.dat <- CORT %>%
  group_by(treatment) %>% # hormone treatment
  summarise(n = n()) %>% 
  rename(hormone = treatment)
# combine for plot
data_fig_2_emm <- left_join(data_fig_2_emm, # emmean results
                           fig.2.sum.dat, # sample size per treatment
                           by = "hormone") %>% 
  mutate(hormone = case_when(hormone == "C_Topical" ~ "control",
                             hormone == "CORT_5pg_Topical" ~ "low",
                             hormone == "CORT_10pg_Topical" ~ 'high',
                             TRUE ~ hormone))

# Grouping for plot
data_fig_2_emm$hormone <- factor(data_fig_2_emm$hormone, levels = c("control", "low", "high"))
data_fig_2_raw$hormone <- factor(data_fig_2_raw$hormone, levels = c("control", "low", "high"))

# FIGURE 1
Figure_2<- ggplot() +
  geom_errorbar(data = data_fig_2_emm, 
                aes(x = hormone, ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_point(data = data_fig_2_emm,
             aes(x = hormone, y = emmean, fill = hormone), 
             shape = 23, color = "black", size = 3) +
  #geom_jitter(data = data_fig_2_raw, 
              #aes(x = hormone,
                  #y = log(CORT_value), color = hormone), 
              #width = 0.1, alpha = 0.6, size = 2) +
  geom_text(data = data_fig_2_emm, aes(x = hormone, 
                                      y = emmean + SE + .01, label = paste("n =", n),
                                      fontface = "bold"), 
            vjust = -0.5, size = 4) +
  labs(y = "Log CORT", 
       x = "Corticosterone treatment") +
  scale_y_continuous(limits = c(.5, 3.0), breaks = seq(0.5, 3.0, by = .5)) +
  #scale_y_continuous(limits = c(-., 3.3), breaks = seq(0.0, 3.0, by = 1.0)) +
  scale_fill_manual(values = c("control" = "white", 
                               "low" = "lightgreen", 
                               "high" = "darkgreen"), guide = "none") +  # Disable fill legend
  scale_color_manual(values = c("control" = "grey40", 
                                "low" = "lightgreen",
                                "high" = "darkgreen")) +
  theme_classic() +
  theme(axis.text.x = element_text(size= 20), 
        axis.text.y = element_text(size= 20),
        axis.title = element_text(size= 20), 
        legend.text = element_text(size= 20), 
        legend.title = element_text(size= 20),
        legend.position = c(.5, .5), # Position legend in the top left corner
        legend.justification = c(0, 1)) # Adjust the anchor point of the legend
Figure_2




###############################################################
####### Figure 3  growth (MASS & SVL) SVL will go in SUPP
###############################################################
library(cowplot) # combines your plots into panels 
cort_dat$hormone <- factor(cort_dat$hormone, levels = c("control", "low", "high"))

####################
###### FIGURE 2) MASS - Hormone 
######
### hatch model 
Mass_hatch_mod <- lmer(hatch_mass_g ~ temp + hormone + (1|clutch), data = data_final)

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
Mass_Juv1_mod <- lmer(juv1_mass_g ~ temp + hormone + scale(Juvenile1_Age) + (1|clutch), data = data_final)

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
Mass_adult_mod <- lmer(adult_mass_g ~ hormone + temp + sex + 
                       scale(Adult_Age) + (1|clutch), data = data_final)

# Data from model for mean and wiskers 'cort_development_mod'  output
Mass_adult_mod_emm_hormone <- emmeans(Mass_adult_mod, pairwise ~ hormone)
Mass_adult_emm_hormone <- as.data.frame(Mass_adult_mod_emm_hormone$emmeans) %>% 
  mutate(Age = "adult",
         Treatment = "Hormone")
# sample size per group
Mass_adult_dat_sum_hormone <- data_final %>%
  dplyr::select(Lizard_ID, hormone,temp,adult_mass_g) %>% 
  dplyr::filter(!is.na(adult_mass_g)) %>%
  group_by(hormone) %>%
  summarise(n = n())
# df for plotting
mass_adult_hormone <- left_join(Mass_adult_emm_hormone, 
                                Mass_adult_dat_sum_hormone,
                                by = "hormone")

# temp data  - Adult
Mass_adult_mod_emm_temp <- emmeans(Mass_adult_mod, pairwise ~ temp)
Mass_adult_dat_emm_temp <- as.data.frame(Mass_adult_mod_emm_temp$emmeans) %>% 
  mutate(Age = "adult",
         Treatment = "Temperature")
# sample size per group
Mass_adult_dat_sum_temp <- data_final %>%
  dplyr::select(Lizard_ID, temp, adult_mass_g) %>% 
  dplyr::filter(!is.na(adult_mass_g)) %>%
  group_by(temp) %>%
  summarise(n = n())
# df for plotting
mass_adult_temp<- left_join(Mass_adult_dat_emm_temp, 
                            Mass_adult_dat_sum_temp,
                            by = "temp")

#####
### final df for plotting
# hormone data
fig_3_plot_dat_hormone <- rbind(mass_hatch_hormone, 
                                mass_juv_hormone, 
                                mass_adult_hormone)%>% 
  dplyr::rename(Treatment_group = hormone)
# temperature data
fig_3_plot_dat_temp <- rbind(mass_hatch_temp, 
                             mass_juv_temp, 
                             mass_adult_temp) %>% 
  dplyr::rename(Treatment_group = temp)
# final plot data
figure_3_plot_dat <- rbind(fig_3_plot_dat_hormone, 
                           fig_3_plot_dat_temp) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

###########
#### !!Figure 3 plot
# reorder data
figure_3_plot_dat$Treatment_group <- factor(figure_3_plot_dat$Treatment_group, 
                                            levels = c("control", "low", "high",
                                                       '23', '28'))
figure_3_plot_dat$Age <- factor(figure_3_plot_dat$Age, 
                                levels = c("hatchling", 
                                           "juvenile", 
                                           "adult"))
###############
# 3A: hatchling plot - Temp
hatch_temp_data <- figure_3_plot_dat %>% 
  filter(Age == 'hatchling',
         Treatment == 'Temperature')

# raw data
hatch_temp_mass_raw <- data_final %>% 
  dplyr::select(Lizard_ID, hormone,temp,hatch_mass_g) %>% 
  dplyr::rename(Treatment_group = temp)
hatch_temp_mass_raw$Treatment_group <- factor(hatch_temp_mass_raw$Treatment_group, 
                                              levels = c('23', '28'))

#FIGURE 2A
figure_3A <- ggplot() +# Add color for the border of the shapes
  geom_errorbar(data = hatch_temp_data, 
                aes(x = Treatment_group, 
                    ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_point(data = hatch_temp_data, aes(x = Treatment_group, 
                                         y = emmean, 
                                         fill = Treatment_group), 
             shape = 23, color = "black", size =3) +
  geom_text(data = hatch_temp_data, aes(x = Treatment_group, y = emmean + 
                                        SE + .00005, label = paste("n =", n),
                                        fontface = "bold"), 
            vjust = -0.5, size = 4) +
  labs(y = "Mass (g)", 
       x = NULL,
       fill = 'Treatment') +
  scale_y_continuous(limits = c(.112, .128), breaks = seq(.112, .128, by = 0.004))+
  theme_classic() +
  scale_fill_manual(values = c("dodgerblue", "tomato2")) + 
  theme(axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_line(), 
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14, face = "bold", margin = margin(b = -1.5)), # Moves title closer to items
        legend.position = c(.01, .25), 
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA), 
        legend.key.size = unit(0.5, "lines"), 
        axis.text.x = element_blank())
figure_3A

###############
# hatchling plot - hormone
hatch_hormone_data <- figure_3_plot_dat %>% filter(Age == 'hatchling',
                                                   Treatment == 'Hormone')
# raw data
hatch_hormone_mass_raw <- data_final %>% 
  dplyr::select(Lizard_ID, hormone,temp,hatch_mass_g) %>% 
  dplyr::rename(Treatment_group = hormone)
hatch_hormone_mass_raw$Treatment_group <- factor(hatch_hormone_mass_raw$Treatment_group, 
                                                 levels = c("control", "low", "high"))

# FIGURE 2B
figure_3B <- ggplot() +
  geom_errorbar(data = hatch_hormone_data, 
                aes(x = Treatment_group, 
                    ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_point(data = hatch_hormone_data, aes(x = Treatment_group, 
                                            y = emmean, 
                                            fill = Treatment_group), 
             shape = 23, color = "black", size = 3) +
  geom_text(data = hatch_hormone_data, aes(x = Treatment_group, y = emmean + 
                                             SE + .00005, label = paste("n =", n)), 
            vjust = -0.5, size = 4, fontface = "bold") +
  labs(y = NULL, 
       x = NULL,
       fill = 'Treatment') +
  theme_classic() +
  scale_y_continuous(limits = c(.112, .128), breaks = seq(.112, .128, by = 0.004))+
  scale_fill_manual(values = c("white", "lightgreen", "darkgreen")) +
  theme(axis.title.y = element_blank(), 
        axis.ticks.y = element_line(),    # Keeps y-axis ticks
        axis.text.y = element_blank(),    # Removes y-axis text (numbers)
        axis.line.y = element_line(),     # Keeps y-axis line
        axis.title.x = element_blank(), 
        axis.ticks.x = element_line(),  
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14, 
                                    face = "bold", 
                                    margin = margin(b = -2.5)), # Moves title closer to items
        legend.position = c(.01, .25), 
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA), 
        legend.key.size = unit(0.5, "lines"), 
        axis.text.x = element_blank())  # Adjust the anchor point of the legend
figure_3B

###############
# 3C: juvenile plot - Temp
juv_temp_data <- figure_3_plot_dat %>% filter(Age == 'juvenile',
                                              Treatment == 'Temperature')
# raw data
juv_temp_data_raw <- data_final %>% 
  dplyr::select(Lizard_ID, hormone,temp,juv1_mass_g) %>% 
  dplyr::rename(Treatment_group = temp) %>% 
  dplyr::filter(!is.na(juv1_mass_g))
juv_temp_data_raw$Treatment_group <- factor(juv_temp_data_raw$Treatment_group, 
                                            levels = c('23', '28'))

# FIGURE 3C 
figure_3C <- ggplot()+  # Add color for the border of the shapes
  geom_errorbar(data = juv_temp_data, 
                aes(x = Treatment_group, 
                    ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_point(data = juv_temp_data, aes(x = Treatment_group, 
                                       y = emmean, 
                                       fill = Treatment_group), 
             shape = 23, color = "black", size = 3) +
  geom_text(data = juv_temp_data, aes(x = Treatment_group, y = emmean + 
                                        SE + .00005, label = paste("n =", n),
                                      fontface = "bold"), 
            vjust = -0.5, size = 4) +
  labs(y = "Mass (g)", 
       x = NULL,
       fill = 'Treatment') +
  theme_classic() +
  scale_y_continuous(limits = c(.38, .56), breaks = seq(0.4, 0.55, by = 0.05)) +
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
figure_3C


###############
# juvenile plot - hormone
juv_hormone_data <- figure_3_plot_dat %>% filter(Age == 'juvenile',
                                                 Treatment == 'Hormone')
# raw data
juv_hormone_mass_raw <- data_final %>% 
  dplyr::select(Lizard_ID, hormone,temp,juv1_mass_g) %>% 
  dplyr::rename(Treatment_group = hormone) %>% 
  dplyr::filter(!is.na(juv1_mass_g)) 
juv_hormone_mass_raw$Treatment_group <- factor(juv_hormone_mass_raw$Treatment_group, 
                                               levels = c("control", "low", "high"))
# FIGURE 4D
figure_3D <- ggplot() +  # Add color for the border of the shapes
  geom_errorbar(data = juv_hormone_data, 
                aes(x = Treatment_group, 
                    ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_point(data = juv_hormone_data, aes(x = Treatment_group, 
                                          y = emmean, 
                                          fill = Treatment_group), 
             shape = 23, color = "black", size = 3) +
  geom_text(data = juv_hormone_data, aes(x = Treatment_group, y = emmean + 
                                           SE + .00005, label = paste("n =", n)), 
            vjust = -0.5, size = 4) +
  labs(y = "Mass (g)", 
       x = NULL,
       fill = 'Treatment') +
  geom_text(data = juv_hormone_data, aes(x = Treatment_group, y = emmean + 
                                           SE + .00005, label = paste("n =", n),
                                         fontface = "bold"), 
            vjust = -0.5, size = 4) +
  theme_classic() +
  scale_y_continuous(limits = c(.38, .56), breaks = seq(0.4, 0.55, by = 0.05)) +
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
figure_3D

###############
# 3E: adult plot - Temp
adult_temp_data <- figure_3_plot_dat %>% filter(Age == 'adult',
                                                Treatment == 'Temperature')
# Raw data
adult_temp_mass_raw <- data_final %>% 
  dplyr::select(Lizard_ID, hormone,temp,adult_mass_g) %>% 
  dplyr::rename(Treatment_group = temp)%>% 
  dplyr::filter(!is.na(adult_mass_g)) 
adult_temp_mass_raw$Treatment_group <- factor(adult_temp_mass_raw$Treatment_group, 
                                              levels = c('23', '28'))

# FIGURE 3E
figure_3E <- ggplot() +  # Add color for the border of the shapes
  geom_errorbar(data = adult_temp_data, 
                aes(x = Treatment_group, 
                    ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_point(data = adult_temp_data, aes(x = Treatment_group, 
                                         y = emmean, 
                                         fill = Treatment_group), 
             shape = 23, color = "black", size = 3) +
  geom_text(data = adult_temp_data, aes(x = Treatment_group, y = emmean + 
                                          SE + .00005, label = paste("n =", n),
                                        fontface = "bold"), 
            vjust = -0.5, size = 4) +
  labs(y = "Mass (g)", 
       x = NULL,
       fill = 'Treatment') +
  theme_classic() +
  scale_y_continuous(limits = c(1.38, 1.61), breaks = seq(1.4, 1.6, by = 0.05)) +
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
figure_3E


###############
# adult plot - hormone
adult_hormone_data <- figure_3_plot_dat %>% filter(Age == 'adult',
                                                   Treatment == 'Hormone')
# Raw data
adult_hormone_mass_raw <- data_final %>% 
  dplyr::select(Lizard_ID, hormone,temp,adult_mass_g) %>% 
  dplyr::rename(Treatment_group = hormone) %>% 
  dplyr::filter(!is.na(adult_mass_g)) 
adult_hormone_mass_raw$Treatment_group <- factor(adult_hormone_mass_raw$Treatment_group, 
                                                 levels = c("control", "low", "high"))
# FIGURE 3F
figure_3F <- ggplot() +  # Add color for the border of the shapes
  geom_errorbar(data = adult_hormone_data, 
                aes(x = Treatment_group, 
                    ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_point(data = adult_hormone_data, aes(x = Treatment_group, 
                                            y = emmean, 
                                            fill = Treatment_group), 
             shape = 23, color = "black", size = 3) +
  geom_text(data = adult_hormone_data, aes(x = Treatment_group, y = emmean + 
                                             SE + .00005, label = paste("n =", n),
                                           fontface = "bold"), 
            vjust = -0.5, size = 4) +
  labs(y = "Mass (g)", 
       x = NULL,
       fill = 'Treatment') +
  theme_classic() +
  scale_y_continuous(limits = c(1.38, 1.61), breaks = seq(1.4, 1.6, by = 0.05)) +
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
figure_3F
##################
# Final figure
figure_3AB <- plot_grid(figure_3A, figure_3B, labels = c("A"), label_size = 16, label_fontface = "bold", ncol = 2)
figure_3CD <- plot_grid(figure_3C, figure_3D, labels = c("B"), label_size = 16, label_fontface = "bold", ncol = 2)
figure_3EF <- plot_grid(figure_3E, figure_3F, labels = c("C"), label_size = 16, label_fontface = "bold", ncol = 2)
# Combine all the rows into one figure
Figure_3 <- plot_grid(figure_3AB, figure_3CD, figure_3EF, ncol = 1, align = 'v')
# Display the combined plot
print(Figure_3)




############################################################################
########################## Figure 4 
############################################################################
########## FIGURE 4A
### MODEL INFO
cort_development_mod <- lmer(log(adult_CORT_Final_Hormone_ng_mL) ~ hormone + temp +  sex + scale(Adult_Age) + Plate_CORT_adult + (1 | clutch), data = cort_dat)
cort_development_mod_anova <- round(as.data.frame(Anova(cort_development_mod)), digits = 3)
# pvalues for plot
cort_development_mod_anova_CORT_Pvalue <- cort_development_mod_anova[1,3]
cort_development_mod_anova_Temp_Pvalue <- cort_development_mod_anova[2,3]
# emmeans pvalues for both
cort_development_mod_anova_CORT_emmeans <- as.data.frame(emmeans(cort_development_mod, pairwise ~ hormone)$contrasts) %>% mutate(across(where(is.numeric), ~ round(.x, 3)))
sum.dat.emm.hormone <- as.data.frame(emmeans(cort_development_mod, pairwise ~ hormone)$emmeans) %>% mutate(across(where(is.numeric), ~ round(.x, 3)))
cort_development_mod_anova_Temp_emmeans <- as.data.frame(emmeans(cort_development_mod, pairwise ~ temp)$contrasts) %>% mutate(across(where(is.numeric), ~ round(.x, 3)))
sum.dat.emm.temp <- as.data.frame(emmeans(cort_development_mod, pairwise ~ temp)$emmeans) %>% mutate(across(where(is.numeric), ~ round(.x, 3)))

### RAW DATA
cort_dat$hormone <- factor(cort_dat$hormone, levels = c("control", "low", "high"))
cort_dat$log_adult_CORT_Final_Hormone_ng_mL <- log(cort_dat$adult_CORT_Final_Hormone_ng_mL)
cort_dat_raw <- cort_dat %>% 
  dplyr::select(hormone, temp, log_adult_CORT_Final_Hormone_ng_mL) %>% 
  dplyr::filter(!is.na(log_adult_CORT_Final_Hormone_ng_mL))

########
# Figure 4A
# sample size per group
sum.dat.temp <- cort_dat %>%
  group_by(temp) %>%
  summarise(n = n())
# combine for plot
fig.4a.emm.dat <- left_join(sum.dat.temp, sum.dat.emm.temp, by = "temp")

# Plot
Figure_4_A <- ggplot() +
  #geom_jitter(data = cort_dat_raw, aes(x = hormone, 
  #y = log_adult_CORT_Final_Hormone_ng_mL, 
  #fill = hormone),
  #alpha = 0.3, 
  #width = 0.2, 
  #height = 0, 
  #shape = 21,  # Change shape to one that uses fill
  #color = "black") +
  geom_errorbar(data = fig.4a.emm.dat, aes(x = temp, ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_text(data = fig.4a.emm.dat, aes(x = temp, y = emmean + SE + .01, label = paste("n =", n),
                                       fontface = "bold"), 
            vjust = -0.5, size = 4) +
  geom_point(data = fig.4a.emm.dat, aes(x = temp, y = emmean, fill = temp), 
             shape = 23, color = "black", size = 3) +
  scale_y_continuous(limits = c(4.85, 5.6), breaks = seq(4.8, 5.6, by = .1)) +
  labs(y = "Log CORT", 
       x = "Temperature treatment",
       fill = "temperature") +  # Change legend title here
  annotate("text", 
           x = 1.85, 
           y = 5.55, 
           label = paste0("Temperature treatment: P = ", cort_development_mod_anova_Temp_Pvalue),  
           size = 4) +
  theme_classic() +
  scale_fill_manual(values = c("dodgerblue", "tomato2")) +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14),
        legend.position = c(0.02, 0.98), # Position legend in the top left corner
        legend.justification = c(0, 1), # Adjust the anchor point of the legend
        legend.background = element_rect(fill = "transparent", color = NA), # Make the legend background transparent
        legend.key = element_rect(fill = "transparent", color = NA))
Figure_4_A


########
# Figure 4B
# sample size per group
sum.dat.hormone <- cort_dat %>%
  group_by(hormone) %>%
  summarise(n = n())
# combine for plot
fig.4b.emm.dat <- left_join(sum.dat.hormone, sum.dat.emm.hormone, by = "hormone")

# Plot
Figure_4_B <- ggplot() +
  #geom_jitter(data = cort_dat_raw, aes(x = hormone, 
                                       #y = log_adult_CORT_Final_Hormone_ng_mL, 
                                       #fill = hormone),
              #alpha = 0.3, 
              #width = 0.2, 
              #height = 0, 
              #shape = 21,  # Change shape to one that uses fill
              #color = "black") +
  geom_errorbar(data = fig.4b.emm.dat, aes(x = hormone, ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_text(data = fig.4b.emm.dat, aes(x = hormone, y = emmean + SE + .01, label = paste("n =", n),
                                      fontface = "bold"), 
            vjust = -0.5, size = 4) +
  geom_point(data = fig.4b.emm.dat, aes(x = hormone, y = emmean, fill = hormone), 
             shape = 23, color = "black", size = 3) +
  scale_y_continuous(limits = c(4.85, 5.6), breaks = seq(4.8, 5.6, by = .1)) +
  labs(y = "Log CORT", 
       x = "Corticosterone treatment") +
  annotate("text", 
           x = 2.75, 
           y = 5.55, 
           label = paste0("CORT treatment: P = ", cort_development_mod_anova_CORT_Pvalue, 
                          "\ncontrol~low: P = ", cort_development_mod_anova_CORT_emmeans[1,6],
                          "\ncontrol~high: P = ", cort_development_mod_anova_CORT_emmeans[2,6],
                          "\nlow~high: P = ", cort_development_mod_anova_CORT_emmeans[3,6]),  
           size = 4)+
  theme_classic() +
  scale_fill_manual(values = c("white", "lightgreen", "darkgreen")) +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14),
        legend.position = c(0.02, 0.98), # Position legend in the top left corner
        legend.justification = c(0, 1), # Adjust the anchor point of the legend
        legend.background = element_rect(fill = "transparent", color = NA), # Make the legend background transparent
        legend.key = element_rect(fill = "transparent", color = NA))
Figure_4_B

######
# Final figure
Figure_4 <- plot_grid(Figure_4_A, Figure_4_B, 
                      labels = c("A", 'B'), 
                      label_size = 16, label_fontface = "bold", ncol = 2)



############################################################################
########################## Figure 5 !!KRIS
############################################################################
RCR_mod <- lmer(RCR ~ temp + hormone + sex + scale(adult_age) + scale(adult_mass_g)  +(1|clutch), 
                data = mito_dat)
RCR_mod_anova <- round(as.data.frame(Anova(RCR_mod)), digits = 3)
RCR_mod_anova_CORT_Pvalue <- RCR_mod_anova[1,3]
RCR_mod_anova_CORT_emmeans <- as.data.frame(emmeans(RCR_mod, pairwise ~ temp)$contrasts) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 3)))
rcr.sum.dat.emm <- as.data.frame(emmeans(RCR_mod, pairwise ~ temp)$emmeans) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 3)))


### RAW DATA & SAMPLE SIZE
rcr_raw <- mito_dat %>% 
  filter(!is.na(RCR)) %>% 
  dplyr::select(temp, RCR)
# sample size per group
rcr.sum.dat <- mito_dat %>%
  filter(!is.na(RCR)) %>% 
  group_by(temp) %>%
  summarise(n = n())
# combine for plot
fig.5.emm.dat <- left_join(rcr.sum.dat, rcr.sum.dat.emm, by = "temp")

#### Final figure
Figure_5 <- ggplot() +
  #geom_jitter(data = cort_dat_raw, aes(x = hormone, 
  #y = log_adult_CORT_Final_Hormone_ng_mL, 
  #fill = hormone),
  #alpha = 0.3, 
  #width = 0.2, 
  #height = 0, 
  #shape = 21,  # Change shape to one that uses fill
  #color = "black") +
  geom_errorbar(data = fig.5.emm.dat, aes(x = temp, ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_text(data = fig.5.emm.dat, aes(x = temp, y = emmean + SE + .01, label = paste("n =", n),
                                       fontface = "bold"), 
            vjust = -0.5, size = 4) +
  geom_point(data = fig.5.emm.dat, aes(x = temp, y = emmean, fill = temp), 
             shape = 23, color = "black", size = 3) +
  annotate("text", 
           x = .99, 
           y = 10.2, 
           label = "A",  
           size = 5, fontface = "bold") +
  annotate("text", 
           x = 1.99, 
           y = 8.1, 
           label = "B",  
           size = 5, fontface = "bold") +
  scale_y_continuous(limits = c(5.5, 10.3), breaks = seq(6, 10, by = 1)) +
  labs(y = "RCR", 
       x = "Temperature treatment",
       fill = "temperature") +  # Change legend title here
  annotate("text", 
           x = 2.2, 
           y = 10, 
           label = paste0("Temperature treatment: P = ", RCR_mod_anova_CORT_Pvalue),  
           size = 4) +
  theme_classic() +
  scale_fill_manual(values = c("dodgerblue", "tomato2")) +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14),
        legend.position = c(0.02, 0.98), # Position legend in the top left corner
        legend.justification = c(0, 1), # Adjust the anchor point of the legend
        legend.background = element_rect(fill = "transparent", color = NA), # Make the legend background transparent
        legend.key = element_rect(fill = "transparent", color = NA))


############################################################################
########################## Figure 6
############################################################################
mito_hormone$sex <- as.factor(mito_hormone$sex)
# Recode the levels of the 'sex' factor
Fig6_dat <- mito_hormone %>%
  mutate(sex = recode_factor(sex, "M" = "Male", "F" = "Female"))


######
########################
######  Fig 6 A - BASAL
# full model for sex pvalue
Basal_treatment_mod <- lmer(basal_corrected_pmol ~ temp + hormone + scale(adult_age) + sex + scale(adult_mass_g) + (1|clutch), data = mito_dat)
anova_result <- Anova(Basal_treatment_mod)
# Extract the p-value for the 'sex' variable
# Convert the ANOVA table to a data frame
anova_df <- as.data.frame(anova_result)
# Extract the p-value for 'sex'
# Save as a vector
Fig_6A_pvalue_sex <-  'P < 0.001'
Fig_6A_pvalue_CORT <- round(anova_df[5,3], digits = 3)


# Fit a linear model for each sex separately for slope
Basal_male_mod <- lm(basal_corrected_pmol ~ adult_mass_g, data = Fig6_dat %>% filter(sex == "Male"))
Basal_female_mod <- lm(basal_corrected_pmol ~ adult_mass_g, data = Fig6_dat %>% filter(sex == "Female"))
slope_male_basal <- coef(Basal_male_mod)[2]
slope_female_basal <- coef(Basal_female_mod)[2]

# Plot
Fig_6A <- ggplot(Fig6_dat, aes(x = adult_mass_g, y = basal_corrected_pmol)) +
  geom_point(aes(color = sex, shape = sex), size = 3.5) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  scale_shape_manual(values = c("Male" = 16, "Female" = 17)) +
  scale_color_manual(values = c("Male" = "steelblue", "Female" = 'darkred')) +
  annotate("text", x = 2.0, 
           y =.04, 
           label = paste0("Mass: P = ", Fig_6A_pvalue_CORT,
                          "\nSex: P = ", Fig_6A_pvalue_sex, 
                          "\nSlope (Male): ", round(slope_male_basal, 2),
                          "\nSlope (Female): ", round(slope_female_basal, 2)),  
           size = 4) +
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

Fig_6A



######
########################
######  Fig 6B ADP - OXPHOS 
# full model for sex pvalue
ADP_treatment_mod <- lmer(adp_corrected_pmol  ~temp + hormone + sex + scale(adult_age) + scale(adult_mass_g)  + (1|clutch), data = mito_dat)
# Perform ANOVA
anova_result_adp <- Anova(ADP_treatment_mod)
# Extract the p-value for the 'sex' variable
# Convert the ANOVA table to a data frame
anova_df_adp <- as.data.frame(anova_result_adp)
# Extract the p-value for 'sex'
Fig_6B_Sex_pvalue_adp <- round(anova_df_adp["sex", "Pr(>Chisq)"],3)
Fig_6B_Mass_pvalue_adp <- round(anova_df_adp["scale(adult_mass_g)", "Pr(>Chisq)"],3)
# Fit a linear model for each sex separately
ADP_male_mod <- lm(adp_corrected_pmol ~ adult_mass_g, data = Fig6_dat %>% filter(sex == "Male"))
ADP_female_mod <- lm(adp_corrected_pmol ~ adult_mass_g, data = Fig6_dat %>% filter(sex == "Female"))
slope_male_adp <- coef(ADP_male_mod)[2]
slope_female_adp <- coef(ADP_female_mod)[2]

# Plot
Fig_6B <- ggplot(Fig6_dat, aes(x = adult_mass_g, y = adp_corrected_pmol)) +
  geom_point(aes(color = sex, shape = sex), size = 3.5) +
  geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  scale_shape_manual(values = c("Male" = 16, "Female" = 17)) +
  scale_color_manual(values = c("Male" = "steelblue", "Female" = 'darkred')) +
  theme_classic() +
  annotate("text", x = 2.0, 
           y =.3, 
           label = paste0("Mass: P = ", Fig_6B_Mass_pvalue_adp, 
                          "\nSex: P = ", Fig_6B_Sex_pvalue_adp, 
                          "\nSlope (Male): ", round(slope_male_adp, 2),
                          "\nSlope (Female): ", round(slope_female_adp, 2)),  
           size = 4) +
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
Fig_6B




######
########################
###### Fig 4C Oligo
Oligo_treatment_mod <- lmer(oligo_corrected_pmol ~ temp + hormone + sex + scale(adult_age) + scale(adult_mass_g)  + (1|clutch), data = mito_dat)
# Perform ANOVA
anova_result_oligo <- Anova(Oligo_treatment_mod)
# Extract the p-value for the 'sex' variable
# Convert the ANOVA table to a data frame
anova_df_oligo <- as.data.frame(anova_result_oligo)
# Extract the p-value for 'sex' and mass
Fig_6C_Sex_pvalue_oligo <- "< 0.001"
Fig_6C_Mass_pvalue_oligo <- round(anova_df_oligo["scale(adult_mass_g)", "Pr(>Chisq)"],3)

# Fit a linear model for each sex separately
Oligo_male_mod <- lm(oligo_corrected_pmol ~ adult_mass_g, data = Fig6_dat %>% filter(sex == "Male"))
Oligo_female_mod <- lm(oligo_corrected_pmol ~ adult_mass_g, data = Fig6_dat %>% filter(sex == "Female"))
slope_male_oligo <- coef(Oligo_male_mod)[2]
slope_female_oligo <- coef(Oligo_female_mod)[2]

# Plot
Fig_6C <- ggplot(Fig6_dat, aes(x = adult_mass_g, y = oligo_corrected_pmol)) +
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
           label = paste0("Mass: P = ", Fig_6C_Mass_pvalue_oligo,
                          "\nSex: P ", Fig_6C_Sex_pvalue_oligo, 
                          "\nSlope (Male): ", round(slope_male_oligo, 2),
                          "\nSlope (Female): ", round(slope_female_oligo, 2)),  
           size = 4) +
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        axis.title.x = element_text(size = 14), 
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14, face = "bold"),
        legend.position = 'none',  
        legend.justification = c(0, 1))
Fig_6C

########################
# Fig 6D - RCR
RCR_mod <- lmer(RCR ~ temp + hormone + sex + scale(adult_age) + scale(adult_mass_g)  +(1|clutch), data = mito_dat)
# Perform ANOVA
anova_result_rcr <- Anova(RCR_mod)
# Extract the p-value for the 'sex' variable
# Convert the ANOVA table to a data frame
anova_df_rcr <- as.data.frame(anova_result_rcr)
# Extract the p-value for 'sex'
Fig_6D_Sex_pvalue_rcr <- round(anova_df_rcr["sex", "Pr(>Chisq)"],3)
Fig_6D_Mass_pvalue_rcr <- round(anova_df_rcr["scale(adult_mass_g)", "Pr(>Chisq)"],3)


# Plot
Fig_6D <- ggplot(Fig6_dat, aes(x = adult_mass_g, y = RCR)) +
  geom_point(aes(color = sex, shape = sex), size = 3.5) +
  #geom_smooth(method = "lm", aes(color = sex), se = FALSE) +
  scale_shape_manual(values = c("Male" = 16, "Female" = 17)) +
  scale_color_manual(values = c("Male" = "steelblue", "Female" = 'darkred')) +
  theme_classic() +
  labs(x = "Mass (g)",
       y = "RCR",
       color = "Sex",
       shape = "Sex") +
  annotate("text", x = 1.84, 
           y =4,  
           label = paste0("Mass: P = ", Fig_6D_Mass_pvalue_rcr, 
                          "\nSex: P = ", Fig_6D_Sex_pvalue_rcr), size = 4)+
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        axis.title.x = element_text(size = 14), 
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14, face = "bold"),
        legend.position = 'none',  
        legend.justification = c(0, 1))

Fig_6D

################################
####### Figure 6
Figure_6 <- plot_grid(
  Fig_6A + ggtitle("A") + theme(plot.title = element_text(hjust = -0.1, vjust = 1.5, face = 'bold')),
  Fig_6B + ggtitle("B") + theme(plot.title = element_text(hjust = -0.1, vjust = 1.5, face = 'bold')),
  Fig_6C + ggtitle("C") + theme(plot.title = element_text(hjust = -0.1, vjust = 1.5, face = 'bold')),
  Fig_6D + ggtitle("D") + theme(plot.title = element_text(hjust = -0.1, vjust = 1.5, face = 'bold')),
  ncol = 2, nrow = 2, labels = NULL)





################################
####### Figure 4 Final
Fig_6A <- Fig_6A + ggtitle("A") + theme(plot.title = element_text(hjust = -0.1, vjust = 1.5, face = 'bold'), 
                                        plot.margin = margin(15, 12, 20, 5.5, "pt"))
Fig_6B <- Fig_6B + ggtitle("B") + theme(plot.title = element_text(hjust = -0.1, vjust = 1.5, face = 'bold'), 
                                        plot.margin = margin(15, 12, 20, 5.5, "pt"))
Fig_6C <- Fig_6C + ggtitle("C") + theme(plot.title = element_text(hjust = -0.1, vjust = 1.5, face = 'bold'), 
                                        plot.margin = margin(15, 12, 20, 5.5, "pt"))
Fig_6D <- Fig_6D + ggtitle("D") + theme(plot.title = element_text(hjust = -0.1, vjust = 1.5, face = 'bold'), 
                                        plot.margin = margin(15, 12, 20, 5.5, "pt"))

# Combine the adjusted plots into a single figure
Figure_6 <- plot_grid(
  Fig_6A,
  Fig_6B,
  Fig_6C,
  Fig_6D,
  ncol = 2, nrow = 2, labels = NULL
)

Figure_6



############################################################################
########################## Figure 7
############################################################################
Fig7_dat <- mito_hormone %>%
  mutate(sex = recode_factor(sex, "M" = "Male", "F" = "Female"))
# full model
Mass_basal <- lmer(hatch_adult_MASS_growth ~ basal_corrected_pmol + sex + log(CORT_Final_Hormone_ng_mL) + log(T4_corrected_ng_mL) + (1|clutch), data = mito_hormone)
anova_result_mass_growth <- Anova(Mass_basal)
# Extract the p-value for the 'sex' variable
# Convert the ANOVA table to a data frame
anova_df_mass_growth <- as.data.frame(anova_result_mass_growth)
# Extract the p-value for 'sex'
sex_pvalue_mass_growth <- anova_df_mass_growth["sex", "Pr(>Chisq)"]
# Round the p-value
rounded_sex_pvalue_mass_growth <- round(sex_pvalue_mass_growth, 3)
# Save as a vector
Fig_7_pvalue_mass_growth <- '<0.001'
# sex 
Mass_basal_male_mod <- lm(hatch_adult_MASS_growth * 1000 ~ log(CORT_Final_Hormone_ng_mL), data = Fig7_dat %>% filter(sex == "Male"))
Mass_basal_female_mod <- lm(hatch_adult_MASS_growth * 1000 ~ log(CORT_Final_Hormone_ng_mL), data = Fig7_dat %>% filter(sex == "Female"))
slope_male_mass_growth <- coef(Mass_basal_male_mod)[2]
slope_female_mass_growth <- coef(Mass_basal_female_mod)[2]

# Plot
Figure_7 <- ggplot(Fig7_dat, 
                   aes(x = log(CORT_Final_Hormone_ng_mL), 
                       y = hatch_adult_MASS_growth * 1000)) +
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
           label = paste0("Sex: P =", Fig_7_pvalue_mass_growth, 
                          "\nSlope (Male): ", round(slope_male_mass_growth, 2),
                          "\nSlope (Female): ", round(slope_female_mass_growth, 2)),  
           size = 4) +
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        axis.title.x = element_text(size = 14), 
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14, face = "bold"),
        legend.position = c(.1, 0.98),   
        legend.justification = c(0, 1))
Figure_7




 
#########   #########   #########   #########   #########   #########   
#########   #########   #########   #########   #########   #########   
#########   #########   #########   #########   #########   #########   
#########   ######### FINAL FIGURES 
Figure_2 #1000x800 png 
Figure_3 #1100x1050
Figure_4 # 1000x850
Figure_5 #1000x800
Figure_6 #1300x1050
Figure_7 #1000x800
#########   #########   #########   #########   #########   #########   
#########   #########   #########   #########   #########   #########   
#########   #########   #########   #########   #########   ######### 










################################  ################################  ################################
################################  ################################  ################################
################################  ###   SUPPLEMENT FIGURES ###  ################################
##################################################
########### Fig 2 Raw
gam_yolk <- readRDS('models/yolk_hormone_mod.RDS')
emm_gam <- emmeans (gam_yolk, pairwise ~ treatment)
# raw data
data_fig_2_raw <- CORT %>% 
  rename(hormone = treatment ) %>% 
  mutate(hormone = case_when(hormone == "C_Topical" ~ "control",
                             hormone == "CORT_5pg_Topical" ~ "low",
                             hormone == "CORT_10pg_Topical" ~ 'high',
                             TRUE ~ hormone))
# Data from model for mean and wiskers 'cort_development_mod'  output
data_fig_2_emm <- as.data.frame(emm_gam$emmeans) %>% 
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>% 
  rename(hormone = treatment)
# sample size per group
fig.2.sum.dat <- CORT %>%
  group_by(treatment) %>% # hormone treatment
  summarise(n = n()) %>% 
  rename(hormone = treatment)
# combine for plot
data_fig_2_emm <- left_join(data_fig_2_emm, # emmean results
                            fig.2.sum.dat, # sample size per treatment
                            by = "hormone") %>% 
  mutate(hormone = case_when(hormone == "C_Topical" ~ "control",
                             hormone == "CORT_5pg_Topical" ~ "low",
                             hormone == "CORT_10pg_Topical" ~ 'high',
                             TRUE ~ hormone))

# Grouping for plot
data_fig_2_emm$hormone <- factor(data_fig_2_emm$hormone, levels = c("control", "low", "high"))
data_fig_2_raw$hormone <- factor(data_fig_2_raw$hormone, levels = c("control", "low", "high"))

# FIGURE 2
Figure_2_Raw<- ggplot() +
  geom_errorbar(data = data_fig_2_emm, 
                aes(x = hormone, ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_point(data = data_fig_2_emm,
             aes(x = hormone, y = emmean, fill = hormone), 
             shape = 23, color = "black", size = 3) +
  geom_jitter(data = data_fig_2_raw, 
              aes(x = hormone,
                  y = log(CORT_value), color = hormone), 
              width = 0.1, alpha = 0.6, size = 2) +
  geom_text(data = data_fig_2_emm, aes(x = hormone, 
                                       y = emmean + SE + .01, label = paste("n =", n),
                                       fontface = "bold"), 
            vjust = -0.5, size = 4) +
  labs(y = "Log CORT", 
       x = "Corticosterone treatment") +
  scale_y_continuous(limits = c(.5, 3.0), breaks = seq(0.5, 3.0, by = .5)) +
  #scale_y_continuous(limits = c(-., 3.3), breaks = seq(0.0, 3.0, by = 1.0)) +
  scale_fill_manual(values = c("control" = "white", 
                               "low" = "lightgreen", 
                               "high" = "darkgreen"), guide = "none") +  # Disable fill legend
  scale_color_manual(values = c("control" = "grey40", 
                                "low" = "lightgreen",
                                "high" = "darkgreen")) +
  theme_classic() +
  theme(axis.text.x = element_text(size= 20), 
        axis.text.y = element_text(size= 20),
        axis.title = element_text(size= 20), 
        legend.text = element_text(size= 20), 
        legend.title = element_text(size= 20),
        legend.position = c(.74, .35), # Position legend in the top left corner
        legend.justification = c(0, 1)) # Adjust the anchor point of the legend
Figure_2_Raw



##################### FIGURE 3 MASS RAW 
#FIGURE 3A
figure_3A <- ggplot() +
  geom_jitter(data = hatch_temp_mass_raw, aes(x = Treatment_group, 
                                              y = hatch_mass_g, 
                                              fill = Treatment_group),
              alpha = 0.3, 
              width = 0.2, 
              height = 0, 
              shape = 21,  # Change shape to one that uses fill
              color = "black") +  # Add color for the border of the shapes
  geom_errorbar(data = hatch_temp_data, 
                aes(x = Treatment_group, 
                    ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_point(data = hatch_temp_data, aes(x = Treatment_group, 
                                         y = emmean, 
                                         fill = Treatment_group), 
             shape = 23, color = "black", size =3) +
  geom_text(data = hatch_temp_data, aes(x = Treatment_group, y = emmean + 
                                          SE + .00005, label = paste("n =", n),
                                        fontface = "bold"), 
            vjust = -0.5, size = 4) +
  labs(y = "Mass (g)", 
       x = NULL,
       fill = 'Treatment') +
  scale_y_continuous(limits = c(0.07, 0.16), breaks = seq(0.08, .16, by = 0.02))+
  theme_classic() +
  scale_fill_manual(values = c("dodgerblue", "tomato2")) + 
  theme(axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_line(), 
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14, face = "bold", margin = margin(b = -1.5)), # Moves title closer to items
        legend.position = c(.01, .25), 
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA), 
        legend.key.size = unit(0.5, "lines"), 
        axis.text.x = element_blank())
figure_3A

###############
# hatchling plot - hormone
hatch_hormone_data <- figure_3_plot_dat %>% filter(Age == 'hatchling',
                                                   Treatment == 'Hormone')
# raw data
hatch_hormone_mass_raw <- data_final %>% 
  dplyr::select(Lizard_ID, hormone,temp,hatch_mass_g) %>% 
  dplyr::rename(Treatment_group = hormone)
hatch_hormone_mass_raw$Treatment_group <- factor(hatch_hormone_mass_raw$Treatment_group, 
                                                 levels = c("control", "low", "high"))

# FIGURE 2B
figure_3B <- ggplot() +
  geom_jitter(data = hatch_hormone_mass_raw, aes(x = Treatment_group, 
                                                 y = hatch_mass_g, 
                                                 fill = Treatment_group),
              alpha = 0.3, 
              width = 0.2, 
              height = 0, 
              shape = 21,  # Change shape to one that uses fill
              color = "black") +  # Add color for the border of the shapes
  geom_errorbar(data = hatch_hormone_data, 
                aes(x = Treatment_group, 
                    ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_point(data = hatch_hormone_data, aes(x = Treatment_group, 
                                            y = emmean, 
                                            fill = Treatment_group), 
             shape = 23, color = "black", size = 3) +
  geom_text(data = hatch_hormone_data, aes(x = Treatment_group, y = emmean + 
                                             SE + .00005, label = paste("n =", n)), 
            vjust = -0.5, size = 4, fontface = "bold") +
  labs(y = NULL, 
       x = NULL,
       fill = 'Treatment') +
  theme_classic() +
  scale_y_continuous(limits = c(0.07, 0.16), breaks = seq(0.08, .16, by = 0.02))+
  scale_fill_manual(values = c("white", "lightgreen", "darkgreen")) +
  theme(axis.title.y = element_blank(), 
        axis.ticks.y = element_line(),    # Keeps y-axis ticks
        axis.text.y = element_blank(),    # Removes y-axis text (numbers)
        axis.line.y = element_line(),     # Keeps y-axis line
        axis.title.x = element_blank(), 
        axis.ticks.x = element_line(),  
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14, 
                                    face = "bold", 
                                    margin = margin(b = -2.5)), # Moves title closer to items
        legend.position = c(.01, .25), 
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA), 
        legend.key.size = unit(0.5, "lines"), 
        axis.text.x = element_blank())  # Adjust the anchor point of the legend
figure_3B
###############
# 3C: juvenile plot - Temp
juv_temp_data <- figure_3_plot_dat %>% filter(Age == 'juvenile',
                                              Treatment == 'Temperature')
# raw data
juv_temp_data_raw <- data_final %>% 
  dplyr::select(Lizard_ID, hormone,temp,juv1_mass_g) %>% 
  dplyr::rename(Treatment_group = temp) %>% 
  dplyr::filter(!is.na(juv1_mass_g))
juv_temp_data_raw$Treatment_group <- factor(juv_temp_data_raw$Treatment_group, 
                                            levels = c('23', '28'))

# FIGURE 3C 
figure_3C <- ggplot()+
  geom_jitter(data = juv_temp_data_raw, aes(x = Treatment_group, 
                                            y = juv1_mass_g, 
                                            fill = Treatment_group),
              alpha = 0.3, 
              width = 0.2, 
              height = 0, 
              shape = 21,  # Change shape to one that uses fill
              color = "black") +  # Add color for the border of the shapes
  geom_errorbar(data = juv_temp_data, 
                aes(x = Treatment_group, 
                    ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_point(data = juv_temp_data, aes(x = Treatment_group, 
                                       y = emmean, 
                                       fill = Treatment_group), 
             shape = 23, color = "black", size = 3) +
  geom_text(data = juv_temp_data, aes(x = Treatment_group, y = emmean + 
                                        SE + .00005, label = paste("n =", n),
                                      fontface = "bold"), 
            vjust = -0.5, size = 4) +
  labs(y = "Mass (g)", 
       x = NULL,
       fill = 'Treatment') +
  theme_classic() +
  scale_y_continuous(limits = c(0.25, .74), breaks = seq(0.3, 0.7, by = 0.1)) +
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
figure_3C


###############
# juvenile plot - hormone
juv_hormone_data <- figure_3_plot_dat %>% filter(Age == 'juvenile',
                                                 Treatment == 'Hormone')
# raw data
juv_hormone_mass_raw <- data_final %>% 
  dplyr::select(Lizard_ID, hormone,temp,juv1_mass_g) %>% 
  dplyr::rename(Treatment_group = hormone) %>% 
  dplyr::filter(!is.na(juv1_mass_g)) 
juv_hormone_mass_raw$Treatment_group <- factor(juv_hormone_mass_raw$Treatment_group, 
                                               levels = c("control", "low", "high"))
# FIGURE 3D
figure_3D <- ggplot() +
  geom_jitter(data = juv_hormone_mass_raw, aes(x = Treatment_group, 
                                               y = juv1_mass_g, 
                                               fill = Treatment_group),
              alpha = 0.3, 
              width = 0.2, 
              height = 0, 
              shape = 21,  # Change shape to one that uses fill
              color = "black") +  # Add color for the border of the shapes
  geom_errorbar(data = juv_hormone_data, 
                aes(x = Treatment_group, 
                    ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_point(data = juv_hormone_data, aes(x = Treatment_group, 
                                          y = emmean, 
                                          fill = Treatment_group), 
             shape = 23, color = "black", size = 3) +
  geom_text(data = juv_hormone_data, aes(x = Treatment_group, y = emmean + 
                                           SE + .00005, label = paste("n =", n)), 
            vjust = -0.5, size = 4) +
  labs(y = "Mass (g)", 
       x = NULL,
       fill = 'Treatment') +
  geom_text(data = juv_hormone_data, aes(x = Treatment_group, y = emmean + 
                                           SE + .00005, label = paste("n =", n),
                                         fontface = "bold"), 
            vjust = -0.5, size = 4) +
  theme_classic() +
  scale_y_continuous(limits = c(0.25, .74), breaks = seq(0.3, 0.7, by = 0.1)) +
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
figure_3D

###############
# 3E: adult plot - Temp
adult_temp_data <- figure_3_plot_dat %>% filter(Age == 'adult',
                                                Treatment == 'Temperature')
# Raw data
adult_temp_mass_raw <- data_final %>% 
  dplyr::select(Lizard_ID, hormone,temp,adult_mass_g) %>% 
  dplyr::rename(Treatment_group = temp)%>% 
  dplyr::filter(!is.na(adult_mass_g)) 
adult_temp_mass_raw$Treatment_group <- factor(adult_temp_mass_raw$Treatment_group, 
                                              levels = c('23', '28'))

# FIGURE 3E
figure_3E <- ggplot() +
  geom_jitter(data = adult_temp_mass_raw, aes(x = Treatment_group, 
                                              y = adult_mass_g, 
                                              fill = Treatment_group),
              alpha = 0.3, 
              width = 0.2, 
              height = 0, 
              shape = 21,  # Change shape to one that uses fill
              color = "black") +  # Add color for the border of the shapes
  geom_errorbar(data = adult_temp_data, 
                aes(x = Treatment_group, 
                    ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_point(data = adult_temp_data, aes(x = Treatment_group, 
                                         y = emmean, 
                                         fill = Treatment_group), 
             shape = 23, color = "black", size = 3) +
  geom_text(data = adult_temp_data, aes(x = Treatment_group, y = emmean + 
                                          SE + .00005, label = paste("n =", n),
                                        fontface = "bold"), 
            vjust = -0.5, size = 4) +
  labs(y = "Mass (g)", 
       x = NULL,
       fill = 'Treatment') +
  theme_classic() +
  scale_y_continuous(limits = c(1.05, 2.25), breaks = seq(1.25, 2.25, by = 0.25)) +
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
figure_3E


###############
# adult plot - hormone
adult_hormone_data <- figure_3_plot_dat %>% filter(Age == 'adult',
                                                   Treatment == 'Hormone')
# Raw data
adult_hormone_mass_raw <- data_final %>% 
  dplyr::select(Lizard_ID, hormone,temp,adult_mass_g) %>% 
  dplyr::rename(Treatment_group = hormone) %>% 
  dplyr::filter(!is.na(adult_mass_g)) 
adult_hormone_mass_raw$Treatment_group <- factor(adult_hormone_mass_raw$Treatment_group, 
                                                 levels = c("control", "low", "high"))
# FIGURE 3F
figure_3F <- ggplot() +
  geom_jitter(data = adult_hormone_mass_raw, aes(x = Treatment_group, 
                                                 y = adult_mass_g, 
                                                 fill = Treatment_group),
              alpha = 0.3, 
              width = 0.2, 
              height = 0, 
              shape = 21,  # Change shape to one that uses fill
              color = "black") +  # Add color for the border of the shapes
  geom_errorbar(data = adult_hormone_data, 
                aes(x = Treatment_group, 
                    ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_point(data = adult_hormone_data, aes(x = Treatment_group, 
                                            y = emmean, 
                                            fill = Treatment_group), 
             shape = 23, color = "black", size = 3) +
  geom_text(data = adult_hormone_data, aes(x = Treatment_group, y = emmean + 
                                             SE + .00005, label = paste("n =", n),
                                           fontface = "bold"), 
            vjust = -0.5, size = 4) +
  labs(y = "Mass (g)", 
       x = NULL,
       fill = 'Treatment') +
  theme_classic() +
  scale_y_continuous(limits = c(1.05, 2.25), breaks = seq(1.25, 2.25, by = 0.25)) +
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
figure_3F
##################
# Final figure
figure_3AB <- plot_grid(figure_3A, figure_3B, labels = c("A"), label_size = 16, label_fontface = "bold", ncol = 2)
figure_3CD <- plot_grid(figure_3C, figure_3D, labels = c("B"), label_size = 16, label_fontface = "bold", ncol = 2)
figure_3EF <- plot_grid(figure_3E, figure_3F, labels = c("C"), label_size = 16, label_fontface = "bold", ncol = 2)
# Combine all the rows into one figure
Figure_3_Mass_Raw <- plot_grid(figure_3AB, figure_3CD, figure_3EF, ncol = 1, align = 'v')
# Display the combined plot
print(Figure_3)


######### FIUGRE 3 SVL example (SUPPLEMENT) RAW
SVL_hatch_mod <- lmer(hatch_svl_mm ~ temp + hormone + (1|clutch), data = data_final)
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
SVL_Juv1_mod <- lmer(juv1_SVL_mm ~ temp + hormone + scale(Juvenile1_Age) + (1|clutch), data = data_final)

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
SVL_adult_mod <- lm(adult_SVL_mm ~ hormone + temp + sex + 
                      scale(Adult_Age), data = data_final)

# Data from model for mean and wiskers 'cort_development_mod'  output
SVL_adult_mod_emm_hormone <- emmeans(SVL_adult_mod, pairwise ~ hormone)
SVL_adult_emm_hormone <- as.data.frame(SVL_adult_mod_emm_hormone$emmeans) %>% 
  mutate(Age = "adult",
         Treatment = "Hormone")
# sample size per group
SVL_adult_dat_sum_hormone <- data_final %>%
  dplyr::select(Lizard_ID, hormone,temp,adult_SVL_mm) %>% 
  dplyr::filter(!is.na(adult_SVL_mm)) %>%
  group_by(hormone) %>%
  summarise(n = n())
# df for plotting
SVL_adult_hormone <- left_join(SVL_adult_emm_hormone, 
                               SVL_adult_dat_sum_hormone,
                               by = "hormone")

# temp data  - Adult
SVL_adult_mod_emm_temp <- emmeans(SVL_adult_mod, pairwise ~ temp)
SVL_adult_dat_emm_temp <- as.data.frame(SVL_adult_mod_emm_temp$emmeans) %>% 
  mutate(Age = "adult",
         Treatment = "Temperature")
# sample size per group
SVL_adult_dat_sum_temp <- data_final %>%
  dplyr::select(Lizard_ID, temp, adult_SVL_mm) %>% 
  dplyr::filter(!is.na(adult_SVL_mm)) %>%
  group_by(temp) %>%
  summarise(n = n())
# df for plotting
SVL_adult_temp<- left_join(SVL_adult_dat_emm_temp, 
                           SVL_adult_dat_sum_temp,
                           by = "temp")

#####
### final df for plotting
# hormone data
fig_3_plot_dat_hormone <- rbind(SVL_hatch_hormone, 
                                SVL_juv_hormone, 
                                SVL_adult_hormone)%>% 
  dplyr::rename(Treatment_group = hormone)
# temperature data
fig_3_plot_dat_temp <- rbind(SVL_hatch_temp, 
                             SVL_juv_temp, 
                             SVL_adult_temp) %>% 
  dplyr::rename(Treatment_group = temp)
# final plot data
fig_3_plot_dat <- rbind(fig_3_plot_dat_hormone, 
                        fig_3_plot_dat_temp) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

###########
#### !!Figure 3 plot SVL
# reorder data
fig_3_plot_dat$Treatment_group <- factor(fig_3_plot_dat$Treatment_group, 
                                         levels = c("control", "low", "high",
                                                    '23', '28'))
fig_3_plot_dat$Age <- factor(fig_3_plot_dat$Age, 
                             levels = c("hatchling", 
                                        "juvenile", 
                                        "adult"))



###############
# 3A: hatchling plot - Temp
hatch_temp_data_SVL <- fig_3_plot_dat %>% filter(Age == 'hatchling',
                                                 Treatment == 'Temperature')

# raw data
hatch_temp_svl_raw <- data_final %>% 
  dplyr::select(Lizard_ID, hormone,temp,hatch_svl_mm) %>% 
  dplyr::rename(Treatment_group = temp)
hatch_temp_svl_raw$Treatment_group <- factor(hatch_temp_svl_raw$Treatment_group, 
                                             levels = c('23', '28'))

# Figure S3A
figure_3A <- ggplot() +
  geom_jitter(data = hatch_temp_svl_raw, aes(x = Treatment_group, 
                                             y = hatch_svl_mm, 
                                             fill = Treatment_group),
              alpha = 0.3, 
              width = 0.2, 
              height = 0, 
              shape = 21,  # Change shape to one that uses fill
              color = "black") +  # Add color for the border of the shapes
  geom_errorbar(data = hatch_temp_data_SVL, 
                aes(x = Treatment_group, 
                    ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_point(data = hatch_temp_data_SVL, aes(x = Treatment_group, 
                                             y = emmean, 
                                             fill = Treatment_group), 
             shape = 23, color = "black", size = 3) +
  geom_text(data = hatch_temp_data_SVL, aes(x = Treatment_group, y = emmean + 
                                              SE + .00005, label = paste("n =", n),
                                            fontface = "bold"), 
            vjust = -0.5, size = 4) +
  labs(y = "SVL (mm)", 
       x = NULL,
       fill = 'Treatment') +
  theme_classic() +
  scale_y_continuous(limits = c(13.9, 20.1), breaks = seq(14, 20, by = 2)) +
  scale_fill_manual(values = c("dodgerblue", "tomato2")) +
  theme(axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_line(), 
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14, face = "bold", margin = margin(b = -1.5)), # Moves title closer to items
        legend.position = c(.01, .25), 
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA), 
        legend.key.size = unit(0.5, "lines"), 
        axis.text.x = element_blank()) 
figure_3A

###############
# hatchling plot - hormone
hatch_hormone_data_SVL <- fig_3_plot_dat %>% filter(Age == 'hatchling',
                                                    Treatment == 'Hormone')
# raw data
hatch_hormone_svl_raw <- data_final %>% 
  dplyr::select(Lizard_ID, hormone,temp,hatch_svl_mm) %>% 
  dplyr::rename(Treatment_group = hormone)
hatch_hormone_svl_raw$Treatment_group <- factor(hatch_hormone_svl_raw$Treatment_group, 
                                                levels = c("control", "low", "high"))
# Figure S3B
figure_3B <- ggplot() +
  geom_jitter(data = hatch_hormone_svl_raw, aes(x = Treatment_group, 
                                                y = hatch_svl_mm, 
                                                fill = Treatment_group),
              alpha = 0.3, 
              width = 0.2, 
              height = 0, 
              shape = 21,  # Change shape to one that uses fill
              color = "black") +  # Add color for the border of the shapes
  geom_errorbar(data = hatch_hormone_data_SVL, 
                aes(x = Treatment_group, 
                    ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_point(data = hatch_hormone_data_SVL, aes(x = Treatment_group, 
                                                y = emmean, 
                                                fill = Treatment_group), 
             shape = 23, color = "black", size = 3) +
  geom_text(data = hatch_hormone_data_SVL, aes(x = Treatment_group, y = emmean + 
                                                 SE + .00005, label = paste("n =", n),
                                               fontface = "bold"), 
            vjust = -0.5, size = 4) +
  labs(y = "SVL (mm)", 
       x = NULL,
       fill = 'Treatment') +
  theme_classic() +
  scale_y_continuous(limits = c(13.9, 20.1), breaks = seq(14, 20, by = 2)) +
  scale_fill_manual(values = c("white", "lightgreen", "darkgreen")) +
  theme(axis.title.y = element_blank(), 
        axis.ticks.y = element_line(),    # Keeps y-axis ticks
        axis.text.y = element_blank(),    # Removes y-axis text (numbers)
        axis.line.y = element_line(),     # Keeps y-axis line
        axis.title.x = element_blank(), 
        axis.ticks.x = element_line(), 
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14, face = "bold", margin = margin(b = -2.5)), # Moves title closer to items
        legend.position = c(.01, .25), 
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA), 
        legend.key.size = unit(0.5, "lines"), 
        axis.text.x = element_blank())  # Adjust the anchor point of the legend
figure_3B

###############
# 3C: juvenile plot - Temp
juv_temp_data_SVL <- fig_3_plot_dat %>% filter(Age == 'juvenile',
                                               Treatment == 'Temperature')
# raw data
juv_temp_svl_raw <- data_final %>% 
  dplyr::select(Lizard_ID, hormone,temp,juv1_SVL_mm) %>% 
  dplyr::rename(Treatment_group = temp) %>% 
  dplyr::filter(!is.na(juv1_SVL_mm)) 
juv_temp_svl_raw$Treatment_group <- factor(juv_temp_svl_raw$Treatment_group, 
                                           levels = c('23', '28'))
# Figure s3C
figure_3C <- ggplot() +
  geom_jitter(data = juv_temp_svl_raw, aes(x = Treatment_group, 
                                           y = juv1_SVL_mm, 
                                           fill = Treatment_group),
              alpha = 0.3, 
              width = 0.2, 
              height = 0, 
              shape = 21,  # Change shape to one that uses fill
              color = "black") +  # Add color for the border of the shapes
  geom_errorbar(data = juv_temp_data_SVL, 
                aes(x = Treatment_group, 
                    ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_point(data = juv_temp_data_SVL, aes(x = Treatment_group, 
                                           y = emmean, 
                                           fill = Treatment_group), 
             shape = 23, color = "black", size = 3) +
  geom_text(data = juv_temp_data_SVL, aes(x = Treatment_group, y = emmean + 
                                            SE + .00005, label = paste("n =", n),
                                          fontface = "bold"), 
            vjust = -0.5, size = 4) +
  labs(y = "SVL (mm)", 
       x = NULL,
       fill = 'Treatment') +
  theme_classic() +
  scale_y_continuous(limits = c(22, 36.4), 
                     breaks = seq(24, 36, by = 2))+
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
juv_hormone_data_SVL <- fig_3_plot_dat %>% filter(Age == 'juvenile',
                                                  Treatment == 'Hormone')
# raw data
juv_hormone_svl_raw <- data_final %>% 
  dplyr::select(Lizard_ID, hormone,temp,juv1_SVL_mm) %>% 
  dplyr::rename(Treatment_group = hormone) %>% 
  dplyr::filter(!is.na(juv1_SVL_mm)) 
juv_hormone_svl_raw$Treatment_group <- factor(juv_hormone_svl_raw$Treatment_group, 
                                              levels = c("control", "low", "high"))

# Figure S3D
figure_3D <- ggplot() +
  geom_jitter(data = juv_hormone_svl_raw, aes(x = Treatment_group, 
                                              y = juv1_SVL_mm, 
                                              fill = Treatment_group),
              alpha = 0.3, 
              width = 0.2, 
              height = 0, 
              shape = 21,  # Change shape to one that uses fill
              color = "black") +  # Add color for the border of the shapes
  geom_errorbar(data = juv_hormone_data_SVL, 
                aes(x = Treatment_group, 
                    ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_point(data = juv_hormone_data_SVL, aes(x = Treatment_group, 
                                              y = emmean, 
                                              fill = Treatment_group), 
             shape = 23, color = "black", size = 3) +
  geom_text(data = juv_hormone_data_SVL, aes(x = Treatment_group, y = emmean + 
                                               SE + .00005, label = paste("n =", n)), 
            vjust = -0.5, size = 4) +
  labs(y = "SVL (mm)", 
       x = NULL,
       fill = 'Treatment') +
  geom_text(data = juv_hormone_data_SVL, aes(x = Treatment_group, y = emmean + 
                                               SE + .00005, label = paste("n =", n),
                                             fontface = "bold"), 
            vjust = -0.5, size = 4) +
  theme_classic() +
  scale_y_continuous(limits = c(22, 36.4), 
                     breaks = seq(24.00, 36.00, by = 2))+
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
figure_3D

###############
# 3E: adult plot - Temp
adult_temp_data_SVL <- fig_3_plot_dat %>% filter(Age == 'adult',
                                                 Treatment == 'Temperature')
# raw data
adult_temp_SVL_raw <- data_final %>% 
  dplyr::select(Lizard_ID, hormone,temp,adult_SVL_mm) %>% 
  dplyr::rename(Treatment_group = temp) %>% 
  dplyr::filter(!is.na(adult_SVL_mm)) 
adult_temp_SVL_raw$Treatment_group <- factor(adult_temp_SVL_raw$Treatment_group, 
                                             levels = c('23', '28'))
#FIGURE S3E
figure_3E <- ggplot() +
  geom_jitter(data = adult_temp_SVL_raw, aes(x = Treatment_group, 
                                             y = adult_SVL_mm, 
                                             fill = Treatment_group),
              alpha = 0.3, 
              width = 0.2, 
              height = 0, 
              shape = 21,  # Change shape to one that uses fill
              color = "black") +  # Add color for the border of the shapes
  geom_errorbar(data = adult_temp_data_SVL, 
                aes(x = Treatment_group, 
                    ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_point(data = adult_temp_data_SVL, aes(x = Treatment_group, 
                                             y = emmean, 
                                             fill = Treatment_group), 
             shape = 23, color = "black", size = 3) +
  geom_text(data = adult_temp_data_SVL, aes(x = Treatment_group, y = emmean + 
                                              SE + .00005, label = paste("n =", n),
                                            fontface = "bold"), 
            vjust = -0.5, size = 4)+
  labs(y = "SVL (mm)", 
       x = NULL,
       fill = 'Treatment') +
  theme_classic() +
  scale_y_continuous(limits = c(36.8, 48.2), breaks = seq(39, 48, by = 3)) +
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
figure_3E


###############
# adult plot - hormone
adult_hormone_data_SVL <- fig_3_plot_dat %>% filter(Age == 'adult',
                                                    Treatment == 'Hormone')
# raw data
adult_hormone_svl_raw <- data_final %>% 
  dplyr::select(Lizard_ID, hormone,temp,adult_SVL_mm) %>% 
  dplyr::rename(Treatment_group = hormone)%>% 
  dplyr::filter(!is.na(adult_SVL_mm))
adult_hormone_svl_raw$Treatment_group <- factor(adult_hormone_svl_raw$Treatment_group, 
                                                levels = c("control", "low", "high"))

# figure S3F
figure_3F <- ggplot() +
  geom_jitter(data = adult_hormone_svl_raw, aes(x = Treatment_group, 
                                                y = adult_SVL_mm, 
                                                fill = Treatment_group),
              alpha = 0.3, 
              width = 0.2, 
              height = 0, 
              shape = 21,  # Change shape to one that uses fill
              color = "black") +  # Add color for the border of the shapes
  geom_errorbar(data = adult_hormone_data_SVL, 
                aes(x = Treatment_group, 
                    ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_point(data = adult_hormone_data_SVL, aes(x = Treatment_group, 
                                                y = emmean, 
                                                fill = Treatment_group), 
             shape = 23, color = "black", size = 3) +
  geom_text(data = adult_hormone_data_SVL, aes(x = Treatment_group, y = emmean + 
                                                 SE + .00005, label = paste("n =", n),
                                               fontface = "bold"), 
            vjust = -0.5, size = 4) +
  labs(y = "SVL (mm)", 
       x = NULL,
       fill = 'Treatment') +
  theme_classic() +
  scale_y_continuous(limits = c(36.8, 48.2), breaks = seq(39, 48, by = 3),
                     labels = label_number(accuracy = 0.01)) +
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
figure_3F
##################
# Final figure
figure_3AB_SVL <- plot_grid(figure_3A, figure_3B, labels = c("A"), label_size = 16, label_fontface = "bold", ncol = 2)
figure_3CD_SVL <- plot_grid(figure_3C, figure_3D, labels = c("B"), label_size = 16, label_fontface = "bold", ncol = 2)
figure_3EF_SVL <- plot_grid(figure_3E, figure_3F, labels = c("C"), label_size = 16, label_fontface = "bold", ncol = 2)
# Combine all the rows into one figure
Figure_3_SVL_Raw <- plot_grid(figure_3AB_SVL, figure_3CD_SVL, figure_3EF_SVL, ncol = 1, align = 'v')


#################   #############      ##############  ###############
#################   #############      ##############  ###############
#################  Figure 4 RAW
Figure_4_A_Raw <- ggplot() +
  geom_jitter(data = cort_dat_raw, aes(x = temp, 
                                       y = log_adult_CORT_Final_Hormone_ng_mL, 
                                       fill = temp),
              alpha = 0.3, 
              width = 0.2, 
              height = 0, 
              shape = 21,  # Change shape to one that uses fill
              color = "black") +
  geom_errorbar(data = fig.4a.emm.dat, aes(x = temp, ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_text(data = fig.4a.emm.dat, aes(x = temp, y = emmean + SE + .01, label = paste("n =", n),
                                       fontface = "bold"), 
            vjust = -0.5, size = 4) +
  geom_point(data = fig.4a.emm.dat, aes(x = temp, y = emmean, fill = temp), 
             shape = 23, color = "black", size = 3) +
  scale_y_continuous(limits = c(3.5, 6.5), breaks = seq(3.5, 6.5, by = .5)) +
  labs(y = "Log CORT", 
       x = "Temperature treatment",
       fill = "temperature") +  # Change legend title here
  annotate("text", 
           x = 1.85, 
           y = 4,
           label = paste0("Temperature treatment: P = ", 
                          cort_development_mod_anova_Temp_Pvalue),  
           size = 4) +
  theme_classic() +
  scale_fill_manual(values = c("dodgerblue", "tomato2")) +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14),
        legend.position = c(0.02, 0.98), # Position legend in the top left corner
        legend.justification = c(0, 1), # Adjust the anchor point of the legend
        legend.background = element_rect(fill = "transparent", color = NA), # Make the legend background transparent
        legend.key = element_rect(fill = "transparent", color = NA))
Figure_4_A_Raw


########
# Figure 4B
# Plot
Figure_4_B_Raw <- ggplot() +
  geom_jitter(data = cort_dat_raw, aes(x = hormone, 
                                       y = log_adult_CORT_Final_Hormone_ng_mL, 
                                       fill = hormone),
              alpha = 0.3, 
              width = 0.2, 
              height = 0, 
              shape = 21,  # Change shape to one that uses fill
              color = "black") +
  geom_errorbar(data = fig.4b.emm.dat, aes(x = hormone, ymin = emmean - SE, ymax = emmean + SE), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_text(data = fig.4b.emm.dat, aes(x = hormone, y = emmean + SE + .01, label = paste("n =", n),
                                       fontface = "bold"), 
            vjust = -0.5, size = 4) +
  geom_point(data = fig.4b.emm.dat, aes(x = hormone, y = emmean, fill = hormone), 
             shape = 23, color = "black", size = 3) +
  scale_y_continuous(limits = c(3.5, 6.5), breaks = seq(3.5, 6.5, by = .5)) +
  labs(y = "Log CORT", 
       x = "Corticosterone treatment") +
  annotate("text", 
           x = 2.75, 
           y = 4.0, 
           label = paste0("CORT treatment: P = ", cort_development_mod_anova_CORT_Pvalue, 
                          "\ncontrol~low: P = ", cort_development_mod_anova_CORT_emmeans[1,6],
                          "\ncontrol~high: P = ", cort_development_mod_anova_CORT_emmeans[2,6],
                          "\nlow~high: P = ", cort_development_mod_anova_CORT_emmeans[3,6]),  
           size = 4)+
  theme_classic() +
  scale_fill_manual(values = c("white", "lightgreen", "darkgreen")) +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14), 
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14),
        legend.position = c(0.02, 0.98), # Position legend in the top left corner
        legend.justification = c(0, 1), # Adjust the anchor point of the legend
        legend.background = element_rect(fill = "transparent", color = NA), # Make the legend background transparent
        legend.key = element_rect(fill = "transparent", color = NA))
Figure_4_B_Raw
# Final figure
Figure_4_Raw <- plot_grid(Figure_4_A_Raw, Figure_4_B_Raw, 
                          labels = c("A", 'B'), 
                          label_size = 16, label_fontface = "bold", ncol = 2)

##########################
# Final Figures
Figure_2_Raw#1000x800 png 
Figure_3_Mass_Raw  #1100x1050
Figure_3_SVL_Raw #1100x1050
Figure_4_Raw #1000x850
