library(ggplot2)
library(lme4) ## estimated the multilevel model (random intercept)
library(dplyr)
library(broom)
library(sjstats) ## partial eta squared and cohens f effective size
library(lmerTest) 
library(emmeans) ## comperhensive anova output with p values
library(MuMIn) ## R for the model
library(bestNormalize)
library(lattice)
library(lmtest)
###___________________________________________________________


dm <- read.csv("data_04_22/filtr_Final_dat_04_22_NoChemTransform.csv")
##______________________M_O_D_E_L__________________________________________________________

## change the level of factor to have CONTROL as and INTERCEPT
dm$treat <- factor(dm$treat, levels = c(
  "control",
  "collar",
  "glue",
  "wax F",
  "wax C"))
levels(dm$treat)

## numerical variables as FACTOR
dm$year <- as.factor(dm$year)
dm$block <- as.factor(dm$block)


vymaz<-dm %>% group_by(year, treat) %>% 
  summarise(damageN = n())

## LINEAR MODELS
lm1 <-  lm(bxcx.x.t ~ treat + year + block + year*treat, data = dm)
lm2 <-  lm(bxcx.x.t ~ treat + year + block , data = dm)
lm3 <-  lm(bxcx.x.t ~ treat + year , data = dm) ## ***
anova( lm2,lm3) ## lm3


## Random effect - INTERCEPTS
lmm1.in1 <- lmer(bxcx.x.t ~ treat + (1|block),data = dm)
lmm1.in2 <-  lmer(bxcx.x.t ~ treat + treat*year + (1|block), data = dm) 
lmm1.in3 <-  lmer(bxcx.x.t ~ treat + year + (1|block), data = dm) ## ***
lmm1.in4 <- lmer(bxcx.x.t ~ treat + (1|block) + (1|year), data = dm)
lmm1.in5 <- lmer(bxcx.x.t ~ treat + (1|block) + (1|year), data = dm)
anova(lmm1.in1, lmm1.in2, lmm1.in3, lmm1.in4, lmm1.in5)
anova(lmm1.in3, lmm1.in4) ## lmm1.in3 ***
anova(lmm1.in3, lm3) ## lmm1.in3 ***
summary(lmm1.in3)
## EMMEANS - Tukey method for comparing a family
emmeans(lmm1.in3, pairwise~treat, type="response", adjustment = "none", confint=T)

 # dm <- filter(dm, damage<2000)
 # dm <- filter(dm, damage>20)
## Random effect - SLOPES
lmm2.sl1 <- lmer(bxcx.x.t ~ treat + (1|block) + (treat|year) + (treat|block), data = dm)
lmm2.sl2 <- lmer(bxcx.x.t ~ treat + (1|block) + (treat|year), data = dm) ## *** !!! Singular fit
lmm2.sl3 <- lmer(bxcx.x.t ~ treat + (treat|year/block), data = dm)
lmm2.sl4 <- lmer(bxcx.x.t ~ treat + (treat|block), data = dm)
anova(lmm2.sl1, lmm2.sl2, lmm2.sl3, lmm2.sl4)
anova(lmm2.sl1, lmm1.in3)
summary(lmm1.in1)
ranova(lmm1.in3)

## Residuals fitted vs pearson
plot(lmm1.in3)
qqmath(lmm1.in3, id = 0.05)


## EMMEANS - Tukey method for comparing a family
emmeans(lmm1.in3, pairwise~treat, type="response", adjustment = "none", confint=T)

## # Random effect chi square, p