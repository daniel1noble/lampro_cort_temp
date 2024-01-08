CORT <- read.csv ("C:/Users/U1115575/Documents/delicata_yolk_hormone_pilot_second.csv")

View(CORT)

CORT$treatment <- as.factor(CORT$treatment)
CORT$CORT_value <- as.numeric(CORT$finalCORT.pg.mg.)
CORT$set <- as.factor(CORT$set)

CORT$treatment <- factor(CORT$treatment, 
                         levels = c("C_Topical", "CORT_5pg_Topical", "CORT_10pg_Topical"))


##linear model with clutch as a random effect
yolkB_mod <- lmer(LogCORT ~ treatment + (1|clutch), data = CORT)
Anova(yolkB_mod)
summary(yolkB_mod)
Mod_residuals <- residuals(yolkB_mod, type = ("pearson"))
shapiro.test(Mod_residuals)
        
yolkB_mod_emm  <- emmeans(yolkB_mod , pairwise ~ treatment)
plot(yolkB_mod_emm)




Boxplot(CORT$CORT_value, CORT$treatment, na.action = na.exclude)
Boxplot(CORT$LogCORT, CORT$treatment, na.action = na.exclude)



LogYolkB_mod <- lmer(LogCORT ~ treatment + set + (1|clutch), data = CORT)
Anova(LogYolkB_mod)
summary(LogYolkB_mod)
emm_LogMod <- emmeans (LogYolkB_mod, pairwise ~ treatment)



Gam_mod <- glmer(CORT_value ~ treatment + (1|clutch), family = Gamma, data = CORT)
Anova(Gam_mod)
summary(Gam_mod)
Gam_mod_emm <- emmeans(Gam_mod, pairwise ~ treatment)
plot(Gam_mod_emm)

Gam_mod_emm <- regrid(emmeans(Gam_mod, "treatment"), transform = "log")
confint(Gam_mod_emm, type = "CORT_value")
pairs(Gam_mod_emm, type = "CORT_value")

emm_Gam_mod <- emmeans (Gam_mod, pairwise ~ treatment)