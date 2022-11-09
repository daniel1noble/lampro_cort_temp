CORT <- read.csv ("C:/Users/U1115575/Documents/delicata_yolk_hormone_current.csv")

View(CORT)

CORT$treatment <- as.factor(CORT$treatment)
CORT$CORT_value <- as.numeric(CORT$final_CORT.pg.mg.)

CORT$LogCORT <- log10(CORT$CORT_value)

yolkB_mod <- lm(CORT_value ~ treatment, data = CORT)
anova(yolkB_mod)
summary(yolkB_mod)
Mod_residuals <- residuals.lm(yolkB_mod, type = ("pearson"))
shapiro.test(Mod_residuals)
        


Boxplot(CORT$CORT_value, CORT$treatment, na.action = na.exclude)
Boxplot(CORT$LogCORT, CORT$treatment, na.action = na.exclude)


LogYolkB_mod <- lm(LogCORT ~ treatment, data = CORT)
Anova(LogYolkB_mod)
summary(LogYolkB_mod)




Gam_mod <- glm(CORT_value ~ treatment, family = Gamma, data = CORT)
Anova(Gam_mod)
summary(Gam_mod)

emm_Gam_mod <- emmeans (Gam_mod, pairwise ~ treatment)