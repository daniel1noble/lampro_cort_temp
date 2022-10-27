CORT <- read.csv ("C:/Users/U1115575/Documents/delicata_yolk_hormone.csv")

View(CORT)

yolkB_mod <- lm(CORT_value ~ treatment, data = CORT)
anova(yolkB_mod)
summary(yolkB

        

CORT$CORT_value <- as.numeric(CORT$final_CORT.pg.mg.)
Boxplot(CORT$CORT_value, CORT$treatment, na.action = na.exclude)

CORT$LogCORT <- log10(CORT$CORT_value)

LogYolkB_mod <- lm(LogCORT ~ treatment, data = CORT)
nova(LogYolkB_mod)
summary(LogYolkB_mod)
