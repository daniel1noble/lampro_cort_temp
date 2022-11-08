T4Para <- read.csv ("C:/Users/U1115575/Documents/thyroid_para_test.csv")

View(T4Para)
T4Para$Type <- as.factor(T4Para$Type)

T4Para$LogT4 <- log10((T4Para$T4_ng.ml))
Para_Test <- lm(LogT4 ~ DilutionFactor + Type, data = T4Para)
anova(Para_Test)
summary(Para_Test)