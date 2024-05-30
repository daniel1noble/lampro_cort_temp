##########
##### Plots for Ondi talk
data_lf <- read.csv("C:/Users/U1115575/Documents/delicata_long_format.csv")

# morphology dataframe
morph_diffs <- data_final %>% 
  select(Liz_ID, hormone, temp, 
         hatch_mass_g, hatch_svl_mm, 
         juv1_mass_g, juv1_SVL_mm,
         juv2_mass_g, juv2_SVL_mm, 
         juv3_mass_g, juv3_SVL_mm,
         hatch_bc_residuals,
         juv1_bc_resid, 
         juv2_bc_resid,
         juv3_bc_resid) %>% 
  rename(Lizard_ID = Liz_ID,
         adult_mass_g = juv3_mass_g,
         adult_SVL_mm = juv3_SVL_mm,
         adult_bc_resid = juv3_bc_resid)

# Merge final sex df based on 'Lizard_ID'
mito_dat_sex <- read.csv(here("Kwild_code/data/final_mito_dat_clean.csv")) %>% 
  select(Lizard_ID, sex) 
morph_diffs_sex <- left_join(morph_diffs, 
                             mito_dat_sex, by = "Lizard_ID")


# reshape data for analysis and plots
morph_diffs_sex_long_format <- gather(morph_diffs_sex, 
                                      key = "test", value = "value", 
                                      hatch_mass_g, hatch_svl_mm, 
                                      juv1_mass_g, juv1_SVL_mm, 
                                      juv2_mass_g, juv2_SVL_mm, 
                                      adult_mass_g, adult_SVL_mm, 
                                      hatch_bc_residuals, 
                                      juv1_bc_resid, juv2_bc_resid, 
                                      adult_bc_resid) %>%
  mutate(Test_final = case_when(
    test %in% c("hatch_mass_g", "juv1_mass_g", "juv2_mass_g", "adult_mass_g") ~ "mass_g", # includes basking and inactive
    test %in% c("hatch_svl_mm", "juv1_SVL_mm", "juv2_SVL_mm", "adult_SVL_mm") ~ "svl_mm",
    test %in% c("hatch_bc_residuals", "juv1_bc_resid", "juv2_bc_resid", "adult_bc_resid") ~ "BC_resid"
  ), 
  hormone = factor(hormone, levels = c("control","low","high"))) %>% 
  filter(!is.na(hormone)) %>%
  group_by(hormone)

### filter by moph of interest
SVL_diffs <- morph_diffs_sex_long_format %>% filter(Test_final == "svl_mm")
mass_diffs <-morph_diffs_sex_long_format %>% filter(Test_final == "mass_g")
bc_diffs <- morph_diffs_sex_long_format %>% filter(Test_final == "BC_resid")

### Plot set up 
hormones <- c("#999999", "#E69F00", "brown2")
my_comparisons_CORT <- rev(list(c("control","low"),c("control","high-"),c("low","high")))



####################################################
####### Figures hormone and by CORT
#############
#### SVL PLOT: hormone
# plot set up
SVL_diffs$test <- factor(SVL_diffs$test, levels = c("hatch_svl_mm", 
                                                      "juv1_SVL_mm", 
                                                      "juv2_SVL_mm", 
                                                      "adult_SVL_mm"))

SVL_diffs_CORT_plot <- ggviolin(SVL_diffs, x = "hormone", y = "value",
                   color = "hormone", palette = hormones,
                   short.panel.labs = FALSE,
                   font.label = list(size = 12, color = "black"),
                   ggtheme = theme_bw()) + 
  geom_jitter(aes(color = hormone), width = 0.2, alpha = 0.6) +
  stat_summary(mapping = aes(x = hormone, y = value), fun = "mean", geom = "point", size = 3, color = "black", alpha = 0.6) +
  stat_summary(mapping = aes(x = hormone, y = value), fun.data = "mean_se", geom = "errorbar", width = 0.1, color = "black", alpha = 0.6) +
  stat_compare_means(comparisons = list(c("control", "low"), c("control", "high"), c("low", "high"))) +
  facet_wrap(~test, scales = "free_y", nrow =1) + # Free scales and ordered facets
  labs(x = NULL, y = "SVL (mm)", color = 'CORT Treatment') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none")

SVL_diffs_CORT_plot



#############
#### MASS PLOT:hormone
mass_diffs$test <- factor(mass_diffs$test, levels = c("hatch_mass_g", 
                                                      "juv1_mass_g", 
                                                      "juv2_mass_g", 
                                                      "adult_mass_g"))
mass_g_CORT_plot <- ggviolin(mass_diffs, x = "hormone", y = "value",
                   color = "hormone", palette = hormones,
                   short.panel.labs = FALSE,
                   font.label = list(size = 12, color = "black"),
                   ggtheme = theme_bw()) + 
  geom_jitter(aes(color = hormone), width = 0.2, alpha = 0.6) +
  stat_summary(mapping = aes(x = hormone, y = value), fun = "mean", geom = "point", size = 3, color = "black", alpha = 0.6) +
  stat_summary(mapping = aes(x = hormone, y = value), fun.data = "mean_se", geom = "errorbar", width = 0.1, color = "black", alpha = 0.6) +
  stat_compare_means(comparisons = list(c("control", "low"), c("control", "high"), c("low", "high"))) +
  facet_wrap(~test, scales = "free_y", nrow =1) + # Free scales and ordered facets
  labs(x = NULL, y = "mass (g)", color = 'CORT Treatment') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none")

mass_g_CORT_plot





#############
#### Bodycondition PLOT
bc_diffs$test <- factor(bc_diffs$test, levels = c("hatch_bc_residuals", 
                                                      "juv1_bc_resid", 
                                                      "juv2_bc_resid", 
                                                      "adult_bc_resid"))
bc_diffs_CORT_plot <- ggviolin(bc_diffs, x = "hormone", y = "value",
                   color = "hormone", palette = hormones,
                   short.panel.labs = FALSE,
                   font.label = list(size = 12, color = "black"),
                   ggtheme = theme_bw()) + 
  geom_jitter(aes(color = hormone), width = 0.2, alpha = 0.6) +
  stat_summary(mapping = aes(x = hormone, y = value), fun = "mean", geom = "point", size = 3, color = "black", alpha = 0.6) +
  stat_summary(mapping = aes(x = hormone, y = value), fun.data = "mean_se", geom = "errorbar", width = 0.1, color = "black", alpha = 0.6) +
  stat_compare_means(comparisons = list(c("control", "low"), 
                                        c("control", "high"), 
                                        c("low", "high"))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~test, scales = "free_y", nrow =1) + # Free scales and ordered facets
  labs(x = NULL, y = "Body Condition (BCI)", color = 'CORT Treatment') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none")

bc_diffs_CORT_plot



######################################################################
###################################
####### Plots by temperature
temp <- c("blue", "red")

SVL_diffs_temp_plot <- ggviolin(SVL_diffs, x = "temp", y = "value",
                                color = "temp", palette = temp,
                                short.panel.labs = FALSE,
                                font.label = list(size = 12, color = "black"),
                                ggtheme = theme_bw()) + 
  geom_jitter(aes(color = temp), width = 0.2, alpha = 0.6) +
  stat_summary(mapping = aes(x = temp, y = value), fun = "mean", geom = "point", size = 3, color = "black", alpha = 0.6) +
  stat_summary(mapping = aes(x = temp, y = value), fun.data = "mean_se", geom = "errorbar", width = 0.1, color = "black", alpha = 0.6) +
  stat_compare_means(comparisons = list(c("23", "28"))) +
  facet_wrap(~test, scales = "free_y", nrow =1) + # Free scales and ordered facets
  labs(x = NULL, y = "SVL (mm)", color = 'Incubation Treatment') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none")

SVL_diffs_temp_plot



#############
#### MASS PLOT:temp
mass_g_temp_plot <- ggviolin(mass_diffs, x = "temp", y = "value",
                             color = "temp", palette = temp,
                             short.panel.labs = FALSE,
                             font.label = list(size = 12, color = "black"),
                             ggtheme = theme_bw()) + 
  geom_jitter(aes(color = temp), width = 0.2, alpha = 0.6) +
  stat_summary(mapping = aes(x = temp, y = value), fun = "mean", geom = "point", size = 3, color = "black", alpha = 0.6) +
  stat_summary(mapping = aes(x = temp, y = value), fun.data = "mean_se", geom = "errorbar", width = 0.1, color = "black", alpha = 0.6) +
  stat_compare_means(comparisons = list(c("23", "28"))) +
  facet_wrap(~test, scales = "free_y", nrow =1) + # Free scales and ordered facets
  labs(x = NULL, y = "mass (g)", color = 'Incubation Treatment') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none")

mass_g_temp_plot

#############
#### Bodycondition PLOT
bc_diffs_temp_plot <- ggviolin(bc_diffs, x = "temp", y = "value",
                               color = "temp", palette = temp,
                               short.panel.labs = FALSE,
                               font.label = list(size = 12, color = "black"),
                               ggtheme = theme_bw()) + 
  geom_jitter(aes(color = temp), width = 0.2, alpha = 0.6) +
  stat_summary(mapping = aes(x = temp, y = value), fun = "mean", geom = "point", size = 3, color = "black", alpha = 0.6) +
  stat_summary(mapping = aes(x = temp, y = value), fun.data = "mean_se", geom = "errorbar", width = 0.1, color = "black", alpha = 0.6) +
  stat_compare_means(comparisons = list(c("23", "28"))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~test, scales = "free_y", nrow =1) + # Free scales and ordered facets
  labs(x = NULL, y = "Body Condition (BCI)", color = 'Incubation Treatment') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none")

bc_diffs_temp_plot


##############################################################################
##############################################################################
############# MITO DATA
#### 
mito_dat_sex <- read.csv(here("Kwild_code/data/final_mito_dat_clean.csv")) %>% 
  filter(RCR_L_R <20) # check lizard ID LD829_21 for RCR(L/R) value X8 higher than other values



#############
### basal_corrected_pmol : Hormone
basal_corrected_pmol_plot <- ggviolin(mito_dat_sex, x = "hormone",
                                      y = "basal_corrected_pmol",
                                color = "hormone", palette = hormones,
                                short.panel.labs = FALSE,
                                font.label = list(size = 12, color = "black"),
                                ggtheme = theme_bw()) + 
  geom_jitter(aes(color = hormone), width = 0.2, alpha = 0.6) +
  stat_summary(mapping = aes(x = hormone, y = basal_corrected_pmol), fun = "mean", geom = "point", size = 3, color = "black", alpha = 0.6) +
  stat_summary(mapping = aes(x = hormone, y = basal_corrected_pmol), fun.data = "mean_se", geom = "errorbar", width = 0.1, color = "black", alpha = 0.6) +
  stat_compare_means(comparisons = list(c("control", "low"), 
                                        c("control", "high"), 
                                        c("low", "high"))) +
  labs(x = NULL, y = "basal_corrected_pmol", color = 'CORT Treatment') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none")
basal_corrected_pmol_plot

#############
### oligo_corrected_pmol : Hormone
oligo_corrected_pmol_plot <- ggviolin(mito_dat_sex, x = "hormone",
                                      y = "oligo_corrected_pmol",
                                      color = "hormone", palette = hormones,
                                      short.panel.labs = FALSE,
                                      font.label = list(size = 12, color = "black"),
                                      ggtheme = theme_bw()) + 
  geom_jitter(aes(color = hormone), width = 0.2, alpha = 0.6) +
  stat_summary(mapping = aes(x = hormone, y = oligo_corrected_pmol), fun = "mean", geom = "point", size = 3, color = "black", alpha = 0.6) +
  stat_summary(mapping = aes(x = hormone, y = oligo_corrected_pmol), fun.data = "mean_se", geom = "errorbar", width = 0.1, color = "black", alpha = 0.6) +
  stat_compare_means(comparisons = list(c("control", "low"), 
                                        c("control", "high"), 
                                        c("low", "high"))) +
  labs(x = NULL, y = "oligo_corrected_pmol", color = 'CORT Treatment') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none")
oligo_corrected_pmol_plot

#############
### fccp_corrected_pmol : Hormone
fccp_corrected_pmol_plot <- ggviolin(mito_dat_sex, x = "hormone",
                                      y = "fccp_corrected_pmol",
                                      color = "hormone", palette = hormones,
                                      short.panel.labs = FALSE,
                                      font.label = list(size = 12, color = "black"),
                                      ggtheme = theme_bw()) + 
  geom_jitter(aes(color = hormone), width = 0.2, alpha = 0.6) +
  stat_summary(mapping = aes(x = hormone, y = fccp_corrected_pmol), fun = "mean", geom = "point", size = 3, color = "black", alpha = 0.6) +
  stat_summary(mapping = aes(x = hormone, y = fccp_corrected_pmol), fun.data = "mean_se", geom = "errorbar", width = 0.1, color = "black", alpha = 0.6) +
  stat_compare_means(comparisons = list(c("control", "low"), 
                                        c("control", "high"), 
                                        c("low", "high"))) +
  labs(x = NULL, y = "fccp_corrected_pmol", color = 'CORT Treatment') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none")
fccp_corrected_pmol_plot

#############
### RCR_L_R : Hormone
RCR_L_R_plot <- ggviolin(mito_dat_sex, x = "hormone",
                                      y = "RCR_L_R",
                                      color = "hormone", palette = hormones,
                                      short.panel.labs = FALSE,
                                      font.label = list(size = 12, color = "black"),
                                      ggtheme = theme_bw()) + 
  geom_jitter(aes(color = hormone), width = 0.2, alpha = 0.6) +
  stat_summary(mapping = aes(x = hormone, y = RCR_L_R), fun = "mean", geom = "point", size = 3, color = "black", alpha = 0.6) +
  stat_summary(mapping = aes(x = hormone, y = RCR_L_R), fun.data = "mean_se", geom = "errorbar", width = 0.1, color = "black", alpha = 0.6) +
  stat_compare_means(comparisons = list(c("control", "low"), 
                                        c("control", "high"), 
                                        c("low", "high"))) +
  labs(x = NULL, y = "RCR_L_R", color = 'CORT Treatment') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none")
RCR_L_R_plot





########################################################
########################################################
########################################################
################## FINAL PLOTS #########################
# Morphology - CORT treatment
SVL_diffs_CORT_plot
mass_g_CORT_plot
bc_diffs_CORT_plot
# Morphology - Temp treatment
SVL_diffs_temp_plot
mass_g_temp_plot
bc_diffs_temp_plot
# MITO data - CORT
basal_corrected_pmol_plot
oligo_corrected_pmol_plot
fccp_corrected_pmol_plot
RCR_L_R_plot