## Kwild Lampro Cort Temp Analysis
# WD, packages, data
pacman::p_load(dplyr, tidyverse, ggpubr, lme4, emmeans, here, plotly)


#########################
# Bring in data from samples and merge together
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


#########################
# dates  and days post hatch days
#########################
# formatting dates
merged_data$hatch_date <- dmy(merged_data$hatch_date_dmy) # hatch
merged_data$juv1_measure_date <- dmy(merged_data$juv1_measure_date_dmy) # juv1
merged_data$juv2_measure_date <- dmy(merged_data$juv2_measure_date_dmy) # juv2
merged_data$juv3_measure_date <- dmy(merged_data$juv3_date) # juv3
merged_data$mortality_date <- dmy(merged_data$mortality_date_dmy) # mortality date - missing final


#########################
# data arrangement up for growth estimations that will be used later
#########################
# calculating days difference between measurments 
merged_data$hatch_juv1_days <- as.numeric(difftime(merged_data$juv1_measure_date,
                                        merged_data$hatch_date, units = "days"))
merged_data$hatch_juv2_days <- as.numeric(difftime(merged_data$juv2_measure_date,
                                        merged_data$hatch_date, units = "days"))
# juv 3 days post hatch
merged_data$hatch_juv3_days <- as.numeric(difftime(merged_data$juv3_measure_date,
                                        merged_data$hatch_date, units = "days"))



############
# Rename columns
data_final <- merged_data %>% 
  rename(hatch_svl_mm = hatch_svl_orig_mm,
         hatch_mass_g = hatch_mass_orig, 
         juv1_SVL_mm = juv1_measure_SVL_mm,
         juv2_SVL_mm = juv2_measure_SVL_mm,
         juv1_mass_g = juv1_measure_mass_g,
         juv2_mass_g = juv2_measure_mass_g)



########
# growth rate calculations: calculated growth rates by 
# dividing change in SVL (or mass) between initial, 
# juv 1, juv 2, juv 3(final) measurements/by the total number of days elapsed
data_final <- data_final %>% 
  mutate(hatch_juv3_SVL_growth =  (juv3_SVL_mm-hatch_svl_mm)/hatch_juv3_days,
         hatch_juv3_MASS_growth = (juv3_mass_g-hatch_mass_g)/hatch_juv3_days)


################
#### 1.	What are the effects of developmental treatments (temp, cort, interaction) on time to hatch?
#### effect of temperature on incubation days - faster hatch warmer temps 
################
# mod ; interaction removed
hatch_dev_mod <- lm(days_to_hatch~ temp + hormone, data = data_final)
check_model(hatch_dev_mod) # bimodal because there is an effect on temp to days of hatch

# temperature effects days to hatch; no effect on hormone or temp x hormone interaction
summary(hatch_dev_mod)
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
SVL_hatch_mod <- lm(hatch_svl_mm ~ temp + hormone , data = data_final)
check_model(SVL_hatch_mod)
# hormone effects body size; no effects on temp or temp x hormone interaction
summary(SVL_hatch_mod)
# high hormone lower SVL
SVL_hatch_mod_emm <- emmeans(SVL_hatch_mod, pairwise ~ hormone)
plot(SVL_hatch_mod_emm)
saveRDS(SVL_hatch_mod, "Kwild_code/models/SVL_hatch_mod.RDS")

#### MASS; interaction removed
Mass_hatch_mod <- lm(hatch_mass_g ~ temp + hormone, data = data_final)
check_model(Mass_hatch_mod)
# temperature and hormone effects mass; no interaction effects
summary(Mass_hatch_mod)
# hormone pairwise comparison - control higher mass
Mass_hormone_hatch_mod_emm <- emmeans(Mass_hatch_mod, pairwise ~ hormone)
plot(Mass_hormone_hatch_mod_emm)
# temp pairwise comparison - cooler temps higher mass
Mass_temp_hatch_mod_emm <- emmeans(Mass_hatch_mod, pairwise ~ temp)
plot(Mass_temp_hatch_mod_emm)
saveRDS(Mass_hatch_mod, "Kwild_code/models/MASS_hatch_mod.RDS")

#### BCI from residuals (OLS)
# hatchlings- residuals
condition_hatch_fit <- lm(hatch_mass_g ~ hatch_svl_mm, data = data_final, na.action=na.exclude) 
data_final$hatch_bc_residuals <- residuals(condition_hatch_fit, na.action=na.exclude) 
## hatchlings-anova on BCI residuals; interaction removed
condition_hatch_mod <- lm(hatch_bc_residuals ~ temp + hormone, data = data_final)
check_model(condition_hatch_mod)
# temperature has effect on BCI; no effect on hormone, or temp x hormone interaction
summary(condition_hatch_mod)
# temp and BCI plot - poor BCI warm temps
condition_hatch_mod_emm <- emmeans(condition_hatch_mod, pairwise ~ temp)
plot(condition_hatch_mod_emm)
saveRDS(condition_hatch_mod, "Kwild_code/models/condition_hatch_mod.RDS")



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
# SVL ANALYSIS: Juv_1, Juv_2, Juv_3 (sex accounted for last measurment)
###############
#### JUV_1: SVL - adding days since hatch because effect of temp on development time); interaction removed
SVL_Juv1_mod <- lm(juv1_SVL_mm ~ temp + hormone + scale(hatch_juv1_days), data = data_final)
check_model(SVL_Juv1_mod)
# effect on temperature and day of hatch on body size 
summary(SVL_Juv1_mod)
# plot - cooler temps larger svl
SVL_Juv1_mod_emm <- emmeans(SVL_Juv1_mod, pairwise ~ temp)
plot(SVL_Juv1_mod_emm)
saveRDS(SVL_Juv1_mod, "Kwild_code/models/SVL_Juv1_mod.RDS")

#### JUV_2: SVL; interaction removed
SVL_Juv2_mod <- lm(juv2_SVL_mm ~ temp + hormone  + scale(hatch_juv2_days), data = data_final)
check_model(SVL_Juv2_mod)
# marginal hermonal effect on body size p = 0.06
summary(SVL_Juv2_mod)
anova(SVL_Juv2_mod)
saveRDS(SVL_Juv2_mod, "Kwild_code/models/SVL_Juv2_mod.RDS")

#### SVL - NOTE SEX ACCOUTNED; no interaction
SVL_Juv3_mod <- lm(juv3_SVL_mm ~ temp + hormone + sex + scale(hatch_juv3_days), data = data_final)
check_model(SVL_Juv3_mod)
# effect on hormones and sex
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
Mass_Juv1_mod <- lm(juv1_mass_g ~ temp + hormone + scale(hatch_juv1_days), data = data_final)
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
Mass_Juv2_mod <- lm(juv2_mass_g ~ temp + hormone  + scale(hatch_juv2_days), data = data_final)
check_model(Mass_Juv2_mod)
# temp  effect
summary(Mass_Juv2_mod)
# temp plot - low mass with high temperatures
Mass_Juv2_temp_mod_emm <- emmeans(Mass_Juv2_mod, pairwise ~ temp)
plot(Mass_Juv2_temp_mod_emm)
saveRDS(Mass_Juv2_mod, "Kwild_code/models/Mass_Juv2_mod.RDS")

#### MASS: Juv_3; NOTE SEX IS ACCOUNTED FOR
Mass_Juv3_mod <- lm(juv3_mass_g ~ temp + hormone+ sex + scale(hatch_juv3_days), data = data_final)
check_model(Mass_Juv3_mod)
# sex effect
summary(Mass_Juv3_mod)
# temp plot - males smaller than females
Sex_Juv3_temp_mod_emm <- emmeans(Mass_Juv3_mod, pairwise ~ sex)
plot(Sex_Juv3_temp_mod_emm)
saveRDS(Mass_Juv3_mod, "Kwild_code/models/Mass_Juv3_mod.RDS")


###############
# BCI ANALYSIS: Juv_1, Juv_2, Juv_3
###############
### BCI Juv1 - residuals
condition_juv1_fit <- lm(juv1_mass_g ~ juv1_SVL_mm + scale(hatch_juv1_days), data = data_final, na.action=na.exclude) # fit the model for residuals
data_final$juv1_bc_resid <- residuals(condition_juv1_fit, na.action=na.exclude) # Save the residual values
## Juvenile -anova on BCI residuals; interaction removed
condition_juv_mod <- lm(juv1_bc_resid ~ temp + hormone  + scale(hatch_juv1_days), data = data_final)
# no treatment effects on bci juv1
summary(condition_juv_mod)
check_model(condition_juv_mod)
# temp plot - 
condition_juv_mod_emm <- emmeans(condition_juv_mod, pairwise ~ temp)
plot(condition_juv_mod_emm)
saveRDS(condition_juv_mod, "Kwild_code/models/condition_juv_mod.RDS")

### BCI Juv2 - residuals
condition_juv2_fit <- lm(juv2_mass_g ~ juv2_SVL_mm + scale(hatch_juv2_days), data = data_final, na.action=na.exclude) # fit the model for residuals
data_final$juv2_bc_resid <- residuals(condition_juv2_fit, na.action=na.exclude) # Save the residual values
## Juvenile -anova on BCI residuals; interaction removed
condition_juv2_mod <- lm(juv2_bc_resid ~ temp + hormone +scale(hatch_juv2_days), data = data_final)
check_model(condition_juv2_mod)
# no treatment effects on bci juv1
summary(condition_juv2_mod)
saveRDS(condition_juv2_mod, "Kwild_code/models/condition_juv2_mod.RDS")

### BCI Juv3 - residuals
condition_juv3_fit <- lm(juv3_mass_g ~ juv3_SVL_mm + scale(hatch_juv2_days), data = data_final, na.action=na.exclude) # fit the model for residuals
data_final$juv3_bc_resid <- residuals(condition_juv3_fit, na.action=na.exclude) # Save the residual values
## Juvenile -anova on BCI residuals; NOTE SEX INCLUDED
condition_juv3_mod <- lm(juv3_bc_resid ~ temp + hormone + sex + scale(hatch_juv2_days), data = data_final)
check_model(condition_juv3_mod)
# no treatment effects on bci juv1
summary(condition_juv3_mod)
# sex plot - 
sex_condition_juv3_mod_emm <- emmeans(condition_juv3_mod, pairwise ~ sex)
plot(sex_condition_juv3_mod_emm)
saveRDS(condition_juv3_mod, "Kwild_code/models/condition_juv3_mod.RDS")


#########
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
#### OVERALL SVL: 
#### OVERALL MASS: effect of temperature-warm temps slower mass growth
################
# SVL overall growth: hatch to juv_3; interaction removed
overall_SVL_growth_mod <- lm(hatch_juv3_SVL_growth ~hormone + temp + sex+ scale(hatch_juv3_days), 
                      data = data_final, na.action=na.exclude)
check_model(overall_SVL_growth_mod)
# treatment effect on temp but not hormone
summary(overall_SVL_growth_mod)
# temp and growth plot - faster growth rates in cooler temps
overall_SVL_growth_mod_emm <- emmeans(overall_SVL_growth_mod, pairwise ~ temp)
plot(overall_SVL_growth_mod_emm)
saveRDS(overall_SVL_growth_mod, "Kwild_code/models/overall_SVL_growth_mod.RDS")

# Mass test; interaction removed
growth4_mass_mod <- lm( hatch_juv3_MASS_growth~hormone + temp + sex+ scale(hatch_juv3_days), 
                       data = data_final, na.action=na.exclude) 
check_model(growth4_mass_mod)
# treatment effect on temp but not hormone
summary(growth4_mass_mod)
# temp and growth plot - faster growth rates in cooler temps
growth4_mass_mod_emm <- emmeans(growth4_mass_mod, pairwise ~ temp)
plot(growth4_mass_mod_emm)
saveRDS(growth4_mass_mod, "Kwild_code/models/growth4_mass_mod.RDS")


########################
### 6: CORT, T4
# Cort: no differences across treatments but overall correlation with SVL and mass growth (high cort, high growth rate)
# T4: 
########################
# data
cort_dat <- data_final %>%
  drop_na(juv3_CORT_Final_Hormone_ng_mL, temp, hormone, juv3_HandlingTime_sec) %>% 
  filter(juv3_Sample_volume_ul >= 4 ) # remove values below 4 ul juv3_Sample_volume_ul

t4_dat <- data_final %>%
  drop_na(juv3_T4_corrected_ng_mL, temp, hormone, Lizard_ID) %>% 
  filter(juv3_Sample_volume_ul >= 4 ) # remove values below 4 ul juv3_Sample_volume_ul

test_dat <- data_final %>%
  drop_na(juv3_Testosterone_Final_ng_ml, temp, hormone, Lizard_ID) %>% 
  filter(juv3_Sample_volume_ul >= 4 ) # remove values below 4 ul juv3_Sample_volume_ul

########
# CORT
######## 1) treatments- how does cort vary across temp and hormone treatment
cort_development_mod <- lm(log(juv3_CORT_Final_Hormone_ng_mL) ~ hormone + temp + sex + juv3_HandlingTime_sec + scale(juv3_mass_g), data = cort_dat)
check_model(cort_development_mod)
summary(cort_development_mod)
saveRDS(cort_development_mod, "Kwild_code/models/cort_development_mod.RDS")

######## 2)  SVL growth - does cort correlate with growth rate and how does this vary across temp and hormone treatment
# high svl growth associated with high cort
cort_SVL_growth_mod <- lm(juv3_CORT_Final_Hormone_ng_mL ~ scale(hatch_juv3_SVL_growth) + temp + hormone + sex, data = cort_dat)
summary(cort_SVL_growth_mod)
saveRDS(cort_SVL_growth_mod, "Kwild_code/models/cort_SVL_growth_mod.RDS")
# plot 
ggplot(cort_dat, aes(x = hatch_juv3_SVL_growth, y = juv3_CORT_Final_Hormone_ng_mL)) +
  geom_point(aes(color = as.factor(temp), shape = hormone), size = 2) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(x = "SVL growth rate (mm/d)", y = "CORT Final Hormone (ng/mL)", 
       color = "Temperature", shape = "Treatment") +
  scale_color_discrete(name = "Temperature") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = c(0.99, 0.99),  # Positions the legend inside the plot
        legend.justification = c("right", "top"), # Justifies the legend's position
        legend.background = element_blank())

######## 3)  MASS growth - does cort correlate with mass growth rate and how does this vary across temp and hormone treatment
# high growth associated with high cort 
cort_mass_growth_mod <- lm(log(juv3_CORT_Final_Hormone_ng_mL) ~ scale(hatch_juv3_MASS_growth) + temp + hormone + sex, data = cort_dat)
summary(cort_mass_growth_mod)
saveRDS(cort_mass_growth_mod, "Kwild_code/models/cort_mass_growth_mod.RDS")
# plot
ggplot(cort_dat, aes(x = hatch_juv3_MASS_growth, y = juv3_CORT_Final_Hormone_ng_mL)) +
  geom_point(aes(color = as.factor(temp), shape = hormone), size = 2) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(x = "Mass growth rate (g/d)", y = "CORT Final Hormone (ng/mL)", 
       color = "Temperature", shape = "Treatment") +
  scale_color_discrete(name = "Temperature") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position = c(0.99, 0.99),  # Positions the legend inside the plot
        legend.justification = c("right", "top"), # Justifies the legend's position
        legend.background = element_blank())

######## 4) Body condition: NS
cort_BCI_growth <- lm(log(juv3_CORT_Final_Hormone_ng_mL) ~juv3_bc_resid, 
                      data = cort_dat)
summary(cort_BCI_growth)


################
# T4
################ 
# 1) treatments -  does T4 vary across temp and hormone treatment: NS
T4_temp_hormone_sex_mod <- lm(juv3_T4_corrected_ng_mL ~  temp + hormone + sex,
                          data = t4_dat)
check_model(T4_temp_hormone_sex_mod)
summary(T4_temp_hormone_sex_mod)
# sex differences: females have higher T4
T4_temp_hormone_sex_mod_emm <- emmeans(T4_temp_hormone_sex_mod, pairwise ~ sex)
plot(T4_temp_hormone_sex_mod_emm)
saveRDS(T4_temp_hormone_sex_mod, "Kwild_code/models/T4_temp_hormone_sex_mod.RDS")

# 2) SVL growth and the interaction between T4 and cort
T4_SVL_growth_mod <- lm(scale(hatch_juv3_SVL_growth) ~ scale(log(juv3_CORT_Final_Hormone_ng_mL))* 
                      scale(log(juv3_T4_corrected_ng_mL)) + temp, 
             data = t4_dat)
check_model(T4_SVL_growth_mod)
summary(T4_SVL_growth_mod)

# 3) MASS growth and the interaction between T4 and cort
T4_mass_growth_mod <- lm(scale(hatch_juv3_MASS_growth) ~ scale(log(juv3_CORT_Final_Hormone_ng_mL))* 
                          scale(log(juv3_T4_corrected_ng_mL)) + temp, 
                        data = t4_dat)
check_model(T4_mass_growth_mod)
summary(T4_mass_growth_mod)
saveRDS(T4_mass_growth_mod, "Kwild_code/models/T4_mass_growth_mod.RDS")
# PLOT: predicitions: growth and T4 and cort
newdata <- data.frame(`juv3_CORT_Final_Hormone_ng_mL` = seq(from = min(t4_dat$juv3_CORT_Final_Hormone_ng_mL), 
                                                            to = max(t4_dat$juv3_CORT_Final_Hormone_ng_mL),
                                                            length.out = 100),
                      `juv3_T4_corrected_ng_mL` = rep(seq(from = min(t4_dat$juv3_T4_corrected_ng_mL), 
                                                          to = max(t4_dat$juv3_T4_corrected_ng_mL)), 
                                                      each = 100), temp = 23)

newdata$pred <- predict(T4_mass_growth_mod, newdata = newdata)

# contour plot
s <- interp(x = newdata$juv3_T4_corrected_ng_mL, 
            y = newdata$juv3_CORT_Final_Hormone_ng_mL, 
            z = newdata$pred)

image.plot(s, xlab = "T4_ng_mL", ylab = "Cort_ng_mL", 
           las = 1, col = viridis(option = "magma", 50), 
           cex.axis = 1, axis.args = list(cex.axis = 1.5), cex.lab = 1.8,
           xlim = c(min(newdata$juv3_T4_corrected_ng_mL)+.012, 
                    max(newdata$juv3_T4_corrected_ng_mL)-.011),  # Explicitly defining x limits
           ylim = c(min(newdata$juv3_CORT_Final_Hormone_ng_mL)+8, 
                    max(newdata$juv3_CORT_Final_Hormone_ng_mL)))  # Explicitly defining y limits
contour(s, add = TRUE, col = "white")


######## 3) growth -  does T4 correlate with growth rate and how does this vary across temp and hormone treatment
# Fit a linear model to T4 the interaction between variables
T4_growth <- lm(juv3_T4_corrected_ng_mL ~ scale(hatch_juv3_SVL_growth) * temp * hormone, data = t4_dat)
summary(T4_growth)





########################
### 7: Mito function
########################
# bring in data, rename, and save
mito_dat <- data_final %>%
  dplyr::select(c("Lizard_ID", "temp", "hormone", "sex", "juv3_bc_resid",
                  "juv3_inject_time_sec", "juv3_liver_time_sec", "juv3_oroboros",
                  "juv3_T4_plate", "juv3_HandlingTime_sec", "juv3_T4_corrected_ng_mL", 
                  "juv3_CORT_Final_Hormone_ng_mL", "juv3_basal_corrected.pmol..sec.ng..",
                  "juv3_chamber","juv3_basal_corrected.pmol..sec.ng..", 
                  "juv3_adp_corrected.pmol..sec.ng..", "juv3_oligo_corrected.pmol..sec.ng..", 
                  "juv3_fccp_corrected.pmol..sec.ng..", "juv3_RCR.L.R.", "juv3_RCR.R.ETS.",
                  "juv3_oroboros_comments", "juv3_RCR.L.ETS.")) %>% 
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
         state2_state3 = basal_corrected_pmol/basal_corrected_pmol) %>% 
  filter(Lizard_ID != c("LD736_21, LD738_21"))
#write.csv(mito_dat, file = "Kwild_code/data/final_mito_dat_clean.csv")

# Scatter Plot of RCR(R/ETS) VS RCR(L/ETS
plot_ly(mito_dat, x = ~HandlingTime_sec, y = ~CORT_Final_Hormone_ng_mL, 
        text = ~ID_and_Comments, type = 'scatter', mode = 'markers') %>%
  layout(title = "Scatter Plot of HandlingTime_sec VS CORT_Final_Hormone_ng_mL ",
         xaxis = list(title = "HandlingTime_sec"),
         yaxis = list(title = "CORT_Final_Hormone_ng_mL")) %>%
  add_trace(marker = list(size = 10, opacity = 0.5), textposition = 'top left')



##########
# RCR (L/ETS): RCR_L_ETS effects by developmental environment and sex (sex:hormone) 
# with inject_time_sec as a covariate   
# interaction removed; no effect on covariate of injection time
RCR_L_ETS_sex_temp_hormone_mod <- glm(RCR_L_ETS ~ temp + hormone + sex  + inject_time_sec, data = mito_dat)
RCR_L_ETS_sex_temp_hormone_mod_summ <- summary(RCR_L_ETS_sex_temp_hormone_mod)
RCR_L_ETS_sex_temp_hormone_mod_summ <- as.data.frame(round(RCR_L_ETS_sex_temp_hormone_mod_summ$coefficients, 
                                                           digits = 2))
RCR_L_ETS_sex_temp_hormone_mod_anova <- as.data.frame(Anova(RCR_L_ETS_sex_temp_hormone_mod, type = 3)) %>% 
  mutate(across(where(is.numeric), round, digits = 2))
check_model(RCR_L_ETS_sex_temp_hormone_mod)
# No differences across treatments or sex or covariet
saveRDS(RCR_L_ETS_sex_temp_hormone_mod, "Kwild_code/models/RCR_L_ETS_sex_temp_hormone_mod.RDS")

#########
### RCR (L/ETS) VS CORT and T4
## RCR (L/ETS) VS CORT - NO EFFECTS
RCR_L_ETS_CORT_mod <- lm(RCR_L_ETS ~ CORT_Final_Hormone_ng_mL, data = mito_dat)
RCR_L_ETS_CORT_mod_summ <- summary(RCR_L_ETS_CORT_mod)
RCR_L_ETS_CORT_mod_summ <- as.data.frame(round(RCR_L_ETS_CORT_mod_summ$coefficients, 
                                               digits = 2))
RCR_L_ETS_CORT_mod_anova <- as.data.frame(Anova(RCR_L_ETS_CORT_mod, type = 3)) %>% 
  mutate(across(where(is.numeric), round, digits = 2))
check_model(RCR_L_ETS_CORT_mod)
saveRDS(RCR_L_ETS_CORT_mod, "Kwild_code/models/RCR_L_ETS_CORT_mod.RDS")

## RCR (R/ETS) VS T4 - EFFECT
RCR_L_ETS_T4_mod <- lm(RCR_L_ETS ~ T4_corrected_ng_mL, data = mito_dat)
RCR_L_ETS_T4_mod_summ <- summary(RCR_L_ETS_T4_mod)
RCR_L_ETS_T4_mod_summ <- as.data.frame(round(RCR_L_ETS_T4_mod_summ$coefficients, 
                                             digits = 2))
RCR_L_ETS_T4_mod_anova <- as.data.frame(Anova(RCR_L_ETS_T4_mod, type = 3)) %>% 
  mutate(across(where(is.numeric), round, digits = 2))
check_model(RCR_L_ETS_T4_mod)
plot(RCR_L_ETS_T4_mod_ano)
# No differences across treatments effects
saveRDS(RCR_L_ETS_T4_mod, "Kwild_code/models/RCR_L_ETS_T4_mod.RDS")


##########
# RCR (R/ETS): RCR_R_ETS effects by developmental environment and sex (sex:hormone) 
# with inject_time_sec as a covariate   
# interaction removed; no effect on covariate of injection time; Sex effects but no developmental effects
RCR_R_ETS_sex_temp_hormone_mod <- glm(RCR_R_ETS ~ temp + hormone + sex  + inject_time_sec , data = mito_dat)
RCR_R_ETS_sex_temp_hormone_mod_summ <- summary(RCR_R_ETS_sex_temp_hormone_mod)
RCR_R_ETS_sex_temp_hormone_mod_summ <- as.data.frame(round(RCR_R_ETS_sex_temp_hormone_mod_summ$coefficients, 
                                                           digits = 2))
RCR_R_ETS_sex_temp_hormone_mod_anova <- as.data.frame(Anova(RCR_R_ETS_sex_temp_hormone_mod, type = 3)) %>% 
  mutate(across(where(is.numeric), round, digits = 2))
check_model(RCR_R_ETS_sex_temp_hormone_mod)
# No differences across treatments or  covariate; but sex effects
RCR_R_ETS_sex_temp_hormone_mod_emm <- emmeans(RCR_R_ETS_sex_temp_hormone_mod, pairwise ~ sex)
plot(RCR_R_ETS_sex_temp_hormone_mod_emm)
saveRDS(RCR_R_ETS_sex_temp_hormone_mod, "Kwild_code/models/RCR_R_ETS_sex_temp_hormone_mod.RDS")

#########
### RCR (R/ETS) VS CORT and T4
## RCR (R/ETS) VS CORT - ONLY SEX DIFFERENCES
# with inject_time_sec as a covariate   
# interaction removed; no effect on covariate of injection time; Sex effects but no developmental effects
RCR_R_ETS_CORT_mod <- lm(RCR_R_ETS ~ CORT_Final_Hormone_ng_mL, data = mito_dat)
RCR_R_ETS_CORT_mod_summ <- summary(RCR_R_ETS_CORT_mod)
RCR_R_ETS_CORT_mod_summ <- as.data.frame(round(RCR_R_ETS_CORT_mod_summ$coefficients, 
                                                           digits = 2))
RCR_R_ETS_CORT_mod_anova <- as.data.frame(Anova(RCR_R_ETS_CORT_mod, type = 3)) %>% 
  mutate(across(where(is.numeric), round, digits = 2))
check_model(RCR_R_ETS_CORT_mod)
saveRDS(RCR_R_ETS_CORT_mod, "Kwild_code/models/RCR_R_ETS_CORT_mod.RDS")

## RCR (R/ETS) VS T4 - 
# with inject_time_sec as a covariate   
# interaction removed; no effect on covariate of injection time; Sex effects but no developmental effects
RCR_R_ETS_T4_mod <- lm(RCR_R_ETS ~ T4_corrected_ng_mL, data = mito_dat)
RCR_R_ETS_T4_mod_summ <- summary(RCR_R_ETS_T4_mod)
RCR_R_ETS_T4_mod_summ <- as.data.frame(round(RCR_R_ETS_T4_mod_summ$coefficients, 
                                               digits = 2))
RCR_R_ETS_T4_mod_anova <- as.data.frame(Anova(RCR_R_ETS_T4_mod, type = 3)) %>% 
  mutate(across(where(is.numeric), round, digits = 2))
check_model(RCR_R_ETS_T4_mod)
saveRDS(RCR_R_ETS_T4_mod, "Kwild_code/models/RCR_R_ETS_T4_mod.RDS")







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