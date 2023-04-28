CORTPara <- read.csv ("C:/Users/U1115575/Documents/CORT para.csv")

View(CORTPara)
CORTPara$Group <- as.factor(CORTPara$Group)


Para_TestB <- lm(Log.value ~ Factor + Group, data = CORTPara)
anova(Para_TestB)
summary(Para_TestB)