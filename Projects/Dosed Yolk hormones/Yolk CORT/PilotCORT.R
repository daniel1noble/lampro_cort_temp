## Install packages
intsall.packages("pacman")
pacman::p_load(tidyverse, readxl, GGally, car)

# Read data. It's found within a separate sheet. Probably better to avoid spaces in file names.
  CORT <- read_xlsx("./Projects/Dosed\ Yolk\ hormones/Yolk\ CORT/delicata_yolk_hormone_current.xlsx", sheet = "CORT\ treament")

# Have a look at the data
  View(CORT)
  str(CORT)

# Do some renaming and classification. Note we have an oredred factor
   CORT$treatment <- factor(CORT$treatment, 
                            levels = c("C_Topical", "CORT_5pg_Topical", "CORT_10pg_Topical"))
  CORT$CORT_value <- as.numeric(CORT$`final_CORT(pg/mg)`)
     CORT$LogCORT <- log10(CORT$CORT_value)

# Summarise data across treatments
   CORT %>% group_by(treatment) %>% summarise(mean_lgcort = mean(LogCORT, na.rm = TRUE),
                                                sd_lgcort = sd(LogCORT, na.rm = TRUE),
                                                        n = n())
   
# Look at the data plotted. Not too convincing to be honest....
   CORT %>% ggplot(aes(x = treatment, y = LogCORT, fill = treatment)) + geom_violin() + 
     geom_point() + labs(y = "log 10 Corticosterone Concentration", x = "Treatment") + 
     theme_bw() + theme(legend.position = 'none')
   
# Create a model. But, model assumes homogeneity of variance....
  yolkB_mod <- lm(CORT_value ~ treatment, data = CORT)
  summary(yolkB_mod)
  anova(yolkB_mod)

# Do a t.test to ditch homogeneity assumption, which is probably not right...
  new_dat <- CORT %>% filter(!treatment == "CORT_5pg_Topical")
  t.test(CORT_value ~ treatment, data = new_dat)

# Residuals don't look bad, but there is a clear outlier > 300. 
  Mod_residuals <- residuals.lm(yolkB_mod, type = ("pearson"))
  shapiro.test(Mod_residuals)
  hist(Mod_residuals)  
  leveneTest(CORT_value ~ treatment, data = new_dat)

# Lets find out what that residual is... Ok, so the really high one is actually form the 5pg treatment.
  CORT[which(Mod_residuals > 300),] %>% select(egg_id, CORT_value, treatment)

# Some extra plots.
Boxplot(CORT$CORT_value, CORT$treatment, na.action = na.exclude)
Boxplot(CORT$LogCORT, CORT$treatment, na.action = na.exclude)

# Try log transformed
  LogYolkB_mod <- lm(LogCORT ~ treatment, data = CORT)
  Anova(LogYolkB_mod)
  summary(LogYolkB_mod)
  
# t.test with log
  t.test(LogCORT ~ treatment, data = new_dat)
leveneTest(LogCORT ~ treatment, data = new_dat) # Although maybe not as bad with log data so can assume homogeneity of variance. 

## Could model with Gamma
  Gam_mod <- glm(CORT_value ~ treatment, family = Gamma, data = CORT)
  Anova(Gam_mod)
  summary(Gam_mod)
  emm_Gam_mod <- emmeans (Gam_mod, pairwise ~ treatment)