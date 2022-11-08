library(tidyverse)

T4_doses <- read.csv ("C:/Users/U1115575/Documents/T4_doses.csv")

View(T4_doses)
T4_doses$treatment <- as.factor(T4_doses$treatment)

##Box plots using non log transformed values
T4_doses$T4_value <- as.numeric(T4_doses$final_T4.ng.mg.)
Boxplot(T4_doses$T4_value, T4_doses$treatment, na.action = na.exclude)

#Box plot using log transformed T4 values
T4_doses$LogT4 <- log10(T4_dosePs$T4_value)
Boxplot(T4_doses$LogT4, T4_doses$treatment, na.action = na.exclude)


##linear regression with all data points and using log transformed hormone values
trt_mod1 <- lm(LogT4 ~ treatment, data = T4_doses)
anova(trt_mod1 )
summary(trt_mod1)


##Subset to remove values where CV >15.0 (n=12)
T4_data2 <- subset(T4_doses, CV<15)
View(T4_data2)

##linear regression with subset data points after high CV values removed and using log transformed hormone values
trt_mod2 <- lm(LogT4 ~ treatment, data = T4_data2)
anova(trt_mod2 )
summary(trt_mod2)

Boxplot(T4_data2$LogT4, T4_data2$treatment, na.action = na.exclude)


##Average values per treatment with sample sizes
aggregate(x = T4_data2$T4_value, 
          by = list(T4_data2$treatment), 
          FUN = mean)
