CORT <- read.csv ("C:/Users/U1115575/Documents/delicata_yolk_hormone_current09112022.csv")

View(CORT)

CORT$Treatment <- as.factor(CORT$Treatment)
CORT$CORT_value <- as.numeric(CORT$final_CORT.pg.mg.)

CORT$LogCORT <- log10(CORT$CORT_value)

Liz_order<- CORT %>% 
  filter(!is.na(Treatment))%>% mutate(Treatment = factor(Treatment,
                                        
                                        levels = c("Control",
                                                   
                                                   "Low CORT",
                                                   
                                                   "High CORT"))) %>%
  
  
  group_by(Treatment) 




###sumarise data for bar chart
sum.dat <- Liz_order  %>%
  
  group_by(Treatment) %>%
  
  summarise(mean.value = mean(CORT_value, na.rm = TRUE),
            
            sd = sd(CORT_value, na.rm = TRUE),
            
            se=mean.value/sqrt(n()))
            
            

            
##makes bar chart            
Bar.plot <- ggplot(sum.dat, aes(x=Treatment, y=mean.value, fill=Treatment)) +
              
              geom_bar(stat="identity", color="black",
                       
                       position=position_dodge()) +
              
              geom_errorbar(aes(ymin=mean.value -se, ymax=mean.value + se), width=.2,
                            
                            position=position_dodge(.9))+
              
              labs(title="Baller figure", y="CORT (pg/mg)", x="Treatment")+
              
              theme_classic()+
          
              
              scale_fill_manual(values = c("white", "grey46", "black"))

Bar.plot + theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20),
                 axis.title = element_text(size = 25), legend.text = element_text(size = 15), 
                 legend.title = element_text(size = 15))
      





yolkB_mod <- lmer(CORT_value ~ treatment + (1|clutch), data = CORT)
anova(yolkB_mod)
summary(yolkB_mod)
Mod_residuals <- residuals.lm(yolkB_mod, type = ("pearson"))
shapiro.test(Mod_residuals)
        


Boxplot(CORT$CORT_value, CORT$treatment, na.action = na.exclude)
Boxplot(CORT$LogCORT, CORT$treatment, na.action = na.exclude)


LogYolkB_mod <- lmer(LogCORT ~ Treatment + (1|clutch), data = CORT)
Anova(LogYolkB_mod)
summary(LogYolkB_mod)
emm_LogMod <- emmeans (LogYolkB_mod, pairwise ~ Treatment)



Gam_mod <- glmer(CORT_value ~ treatment + (1|enclosure), family = Gamma, data = CORT)
anova(Gam_mod)
summary(Gam_mod)

Gam_mod_emm <- regrid(emmeans(Gam_mod, "treatment"), transform = "log")
confint(Gam_mod_emm, type = "CORT_value")
pairs(Gam_mod_emm, type = "CORT_value")

emm_Gam_mod <- emmeans (Gam_mod, pairwise ~ treatment)