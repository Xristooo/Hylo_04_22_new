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

# data <- read.csv("data_2021_Hylob/dataFinal.csv")
data <- read.csv("data_2021_Hylob/dataMORT_Final.csv")
names(data)
head(data)

data <- data %>%
  select(block, treat, year,notes, damage)

##__________________FILTER_______________________________________________________

# NO CHEMICAL, HYLASTES, OTHER DAMAGES
data1 <- data %>%
  filter(damage > 1 & treat!="chemical" & notes!="Hylastes" & notes!="other")
View(data1)

## CHECK the levels !=chemical, Hylastes, other
levels(factor(data1$notes))
levels(factor(data1$treat))

# # with CHEMICAL
# data2 <- data %>%
#   filter(damage > 1 & notes!="Hylastes" & notes!="other")
# View(data2)
# ## CHECK the levels != Hylastes, other
# levels(factor(data2$notes))
# 
# 
# # with CHEMICAL
# data3 <- data %>%
#   filter(notes!="Hylastes" & notes!="other")
# View(data3)
# ## CHECK the levels != Hylastes, other
# levels(factor(data2$notes))


##___________________TRANSFORMATION__________________________________________________________________
levels(factor(data1$treat))

bestNormalize(data1$damage) # 1. Box-Cox: 1.1486, 2. Yeo-Johnson: 1.1967
yeoT <- yeojohnson(data1$damage) 
bxcx <- boxcox(data1$damage)


## combine df, original + transform
## x original data - preserve to check when joining dfs,x.t transform
## x remove after combine frames and the original data match
# ordNormT<-data.frame(BN_objOrd$x.t, BN_objOrd$x) <- yeo performed better
bxcxT <-data.frame(bxcx$x.t, bxcx$x)
## Rename colums both are sam
colnames(bxcxT)[1] <-"bx_TRANS"
colnames(bxcxT)[2] <-"bx_Orig"
## made transformation for both types
datDamMod <- data.frame(data1, bxcxT)
# write.csv(datDamMod, "01_dam_NOchem.csv", row.names = F)

datDamMod$year <- as.factor(datDamMod$year)
datDamMod$block <- as.factor(datDamMod$block) 


View(datDamMod)



##______________________M_O_D_E_L__________________________________________________________

## !!!!!! this is important
## change the level of factor to have CONTROL as intercept
datDamMod$treat <- factor(datDamMod$treat, levels = c(
  "control",
  "collar",
  "glue",
  "wax F",
  "wax C"))


## MODEL - https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#should-i-treat-factor-xxx-as-fixed-or-random


## LINEAR MODELS

bxModel_00a <-  lm(bx_TRANS ~ treat + year + block + year*treat, data = datDamMod)
bxModel_00b <-  lm(bx_TRANS ~ treat + year + block , data = datDamMod)
bxModel_00c <-  lm(bx_TRANS ~ treat + year , data = datDamMod) ##4.179e-07 ***
anova(bxModel_00a, bxModel_00b,bxModel_00c)


## Random effect - INTERCEPTS
bxModel_01a <- lmer(bx_TRANS ~ treat + (1|block),data = datDamMod)
bxModel_01b.0 <-  lmer(bx_TRANS ~ treat + treat*year + (1|block), data = datDamMod) #anova(bxModel_01b.0,bxModel_01b.1) 01b.0<-0.09746.
bxModel_01b.1 <-  lmer(bx_TRANS ~ treat + year + (1|block), data = datDamMod)
bxModel_01c <-  lmer(bx_TRANS ~ treat + year + (1|block), data = datDamMod)
bxModel_01d <- lmer(bx_TRANS ~ treat + (1|block) + (1|year), data = datDamMod)
bxModel_01e <- lmer(bx_TRANS ~ treat + (1|block) + (1|year), data = datDamMod)
anova(bxModel_01b.0,bxModel_01b.1)

# datDamMod <- filter(datDamMod, damage<2000)

## Random effect - SLOPES
bxModel_03a <- lmer(bx_TRANS ~ treat + (1|block) + (treat|year) + (treat|block), data = datDamMod)
bxModel_03b <- lmer(bx_TRANS ~ treat + (1|block) + (treat|year), data = datDamMod)
bxModel_03c <- lmer(bx_TRANS ~ treat + (treat|year/block), data = datDamMod)
bxModel_03c1 <- lmer(bx_TRANS ~ treat + (treat|block), data = datDamMod)
## # Random effect chi square, p
ranova(bxModel)
ranef(bxModel_00a)

## Residuals fitted vs pearson
plot(bxModel_01c)
qqnorm(bxModel_03c)
## Model summary
summary(bxModel_01c)

## Model comparison
anova(bxModel_01b, bxModel_01c)
lrtest(bxModel_00c, bxModel_00b)

## EMMEANS - Tukey method for comparing a family
emmeans(bxModel_01c, pairwise~treat, type="response", adjustment = "none", confint=T)
confint(emBx, adjust = "none", level = 0.99) ## confidence intervals

##  standard fitted vs. residual plots
plot(bxModel_01b, type = c("p", "smooth"))
## scale-location plots
plot(bxModel_01b, sqrt(abs(resid(.))) ~ fitted(.),
     type = c("p", "smooth"))
## quantile-quantile plots (from lattice),
qqmath(bxModel_01b, id = 0.05)

plot(ranef(bxModel_03c))
plot(fixef(bxModel_03c))

## YEAR - ako fixed factor 
## there must be a reasonable number of random-effects levels (e.g. blocks) – more than 5 or 6 at a minimum. This is not surprising if you consider that random effects estimation is trying to estimate an among-block variance. For example, from Crawley (2002) p. 670:

## tabuľky
library(jtools) # summ function give an nicer output table
summ(bxModel_02) 
library(sjPlot)
tab_model(bxModel_02) ## library(sjtools)

##  model simulated
iqrvec <- sapply(simulate(bxModel_03c, 1000), IQR)
obsval <- IQR(datDamMod$bx_TRANS)
post.pred.p <- mean(obsval >= c(obsval, iqrvec))
summary(post.pred.p)



##______________Standardised residuals against BLOCKS______________________________
M1 <- lm(bx_TRANS~treat, data=datDamMod)
E1 <- rstandard(M1)
plot(E1 ~ datDamMod$block, ylab="Standardised residuals", xlab="blocks")
abline(0,0)
##______________Standardised residuals against YEAR______________________________
M2 <- lm(damage~treat, data=datDamMod)
E2 <- rstandard(M2)
plot(E2 ~ datDamMod$year, ylab="Standardised residuals", xlab="years")
abline(0,0)

## MODEL SELECTION
# bxModel npar    AIC    BIC  logLik deviance   Chisq Df Pr(>Chisq) 
# bxModel    8 1179.8 1212.6 -581.92   1163.8 22.0143  1  2.706e-06 ***
# mYeo1 <- lmer(bx_TRANS ~ treat + (1|block), data = datDamMod)
# mYeo2 <- lmer(bx_TRANS ~ treat +(1|year), data = datDamMod)
# mYeo3 <- lmer(bx_TRANS ~ treat +(1|year/block/treat) + (1|block/treat), data = datDamMod)
# mYeo4 <- lmer(bx_TRANS ~ treat +(1|year/treat) + (1|block), data = datDamMod)
# mYeo5 <- lmer(bx_TRANS ~ treat +(treat|year) + (1|block), data = datDamMod)
# 
# anova(bxModel,mYeo,mYeo1,mYeo2,mYeo3,mYeo4,mYeo5)

chem <- data %>%
  group_by(year, treat) %>%
  summarise(
    meanD = mean(damage),
    countD=n_distinct(factor(damage)),
    maxD=max(damage)
  )
ch_2018 <-  chem 




SUMARY_DATA <- data %>%
  group_by(year, treat)  %>%
  filter(site=="clearcut")%>%
  summarise(
    meanDam = mean(damage),
    medianDamage = median(damage),
    maxDamage = max(damage),
    totalDamage = sum(damage),
    count=n_distinct(unIDyearSite),
    SD = sd(damage),
    SE = SD/sqrt(count)
  )
View(SUMARY_DATA)

SUMARY_TOTAL <- data %>%
  group_by(year)  %>%
  filter(site=="clearcut")%>%
  summarise(
    meanDam = mean(damage),
    medianDamage = median(damage),
    maxDamage = max(damage),
    totalDamage = sum(damage),
    count=n_distinct(unIDyearSite),
    SD = sd(damage),
    SE = SD/sqrt(count)
  )
View(SUMARY_TOTAL)




names(data)
# 
# model1.1<- lmer(transDamage.x.t ~ treat +(1|year), data = data)
# model1.2<- lmer(transDamage.x.t ~ treat + (1|block), data = data)
# anova(model1, model1.1, model1.2)
# anova(model1)
# summary(model1)
# 
# 
# model2<- lmer(transDamage.x.t ~ treat +(1|year/treat) + (1|block/treat), data = data)
# model3<- lmer(transDamage.x.t ~ treat +(1|year/treat) + (1|block), data = data)
# model4<- lmer(transDamage.x.t ~ treat +(1+year|treat) + (1|block/treat), data = data)
# model5<- lmer(transDamage.x.t ~ treat + year + (1|block), data = data)
# model6<- lmer(transDamage.x.t ~ treat + (treat|year) + (1|block), data = data)
# 
# emmeansLMM1 <- emmeans(model1, pairwise~treat, type="response", adjustment = "none")
# emmeansLMM1
# summary(model1)
# anova(model1)
# coef(model1)
# 
# plot(ranef(model3))
# plot(fixef(model3))

# as_lmerModLmerTest(model1)
# qqplot(model3)

## no degrees of freedom
##model2<- lmer(damTrans ~ treat +(treat|year) + (1|block), data = data)

## https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#model-diagnostics
summary(model1, ddf="Kenward-Roger")



## Conditional F-tests are preferred for LMMs, if denominator degrees of freedom are known
## https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#model-diagnostics
anova(model1)
anova(model1, ddf="Kenward-Roger")

anova(model1, model2, model3, model4, model5)
eta_sq(model1, partial = T)

## Estimated Marginal Means - Tukey method
emmeans(model1, list(pairwise ~ treat), adjust ="tukey")
## Estimated Marginal Means - Bonferroni method
emmeans(model6, list(pairwise ~ treat), adjust ="tukey")

##______________________________________________________________
## GLMER
# treat as a random-effect intercept and Year as both
# a fixed- and random-effect slope
# data$year<-as.factor(data$year)
# data$block<-as.factor(data$block)
# data$notes<-as.factor(data$notes)
# 
# str(data)
# 
# 
# warnings()
# ## Block and year as a random effect
# m1<- glmer(damage ~ treat+(1|year), data = data, family=poisson, nAGQ =0)
# 
# 
# m2<- glmer(damage ~ treat+(1|year) , data = data, family=poisson)
# 
# m3<- glmer(damage ~ treat+(1+treat|year) + (1+treat|block), data = data, family=poisson)
# 
# anova(m1,m2,m3)
# 
# warnings()
# coef(model1)
# ranef(model1)
# fixef(model1)
# anova(model1)
# 
# ## Estimated Marginal Means - Tukey method
# emmeans(model1, list(pairwise ~ treat), adjust ="Tukey")
# ## Estimated Marginal Means - Bonferroni method
# emmeans(model1, list(pairwise ~ treat), adjust ="bonferroni")
# 
# summary(model1)
# 
# ## QQ residual plot
# qqnorm(resid(model1))
# qqline(resid(model1), col = "steelblue")
# 
# car::qqPlot(data$damage)
# 
# 
# ## Block and year as a random effect
# model11<- glmer(damage ~ treat+(1|year/treat), data = data, family=poisson)
# anova(model11)
# 
# 
# ## Estimated Marginal Means - Tukey method
# emmeans(model11, list(pairwise ~ treat), adjust ="tukey")
# ## Estimated Marginal Means - Bonferroni method
# emmeans(model11, list(pairwise ~ year), adjust ="bonferroni")
# 
# summary(model11)
# 
# ## QQ residual plot
# qqnorm(resid(model11))
# qqline(resid(model11), col = "steelblue")
# 
# 
# 
# ## Block as a random effect
# model2<- glmer(damage ~ treat+(1|block), data = data, family=poisson)
# anova(model2)
# eta_sq(model2, partial = T)
# 
# ## Estimated Marginal Means - Tukey method
# emmeans(model2, list(pairwise ~ treat), adjust ="tukey")
# ## Estimated Marginal Means - Bonferroni method
# emmeans(model2, list(pairwise ~ year), adjust ="bonferroni")
# 
# summary(model2)
# 
# ## QQ residual plot
# qqnorm(resid(model2))
# qqline(resid(model2), col = "steelblue")
# 
# 
# 
# ## Year as a random effect
# model3<- glmer(damage ~ treat+(1|year), data = data, family=poisson)
# anova(model3)
# eta_sq(model3, partial = T)
# 
# ## Estimated Marginal Means - Tukey method
# emmeans(model3, list(pairwise ~ treat), adjust ="tukey")
# ## Estimated Marginal Means - Bonferroni method
# emmeans(model3, list(pairwise ~ year), adjust ="bonferroni")
# 
# summary(model3)
# 
# ## QQ residual plot
# qqnorm(resid(model3))
# qqline(resid(model3), col = "steelblue")
# 
# 
# 
# 
# 
# 
# ## NO random effect
# model4<- glm(damage ~ treat, data = data, family=poisson)
# anova(model4)
# eta_sq(model4, partial = T)
# 
# ## Estimated Marginal Means - Tukey method
# emmeans(model4, list(pairwise ~ treat), adjust ="tukey")
# ## Estimated Marginal Means - Bonferroni method
# emmeans(model4, list(pairwise ~ year), adjust ="bonferroni")
# 
# summary(model4)
# 
# ## QQ residual plot
# qqnorm(resid(model4))
# qqline(resid(model4), col = "steelblue")
# 
# 
# 
# 
# ## NO random effect
# model5<- glmer(damage ~ treat + (treat | year), data = data, family = "poisson")
# anova(model5)
# eta_sq(model5, partial = T)
# 
# ## Estimated Marginal Means - Tukey method
# emmeans(model5, list(pairwise ~ treat), adjust ="tukey")
# ## Estimated Marginal Means - Bonferroni method
# emmeans(model5, list(pairwise ~ year), adjust ="bonferroni")
# 
# summary(model5)
# 
# ## QQ residual plot
# qqnorm(resid(model5))
# qqline(resid(model5), col = "steelblue")
# 
# 
# 
# 
# data1<- data
# data$uniqueID <- as.factor(data$uniqueID)
# 
# ## NO random effect
# model6<- glmer(damage ~ uniqueID + (1 | year), data = data1, family = "poisson")
# anova(model6)
# eta_sq(model6, partial = T)
# 
# ## Estimated Marginal Means - Tukey method
# emmeans(model6, list(pairwise ~ treat), adjust ="tukey")
# ## Estimated Marginal Means - Bonferroni method
# emmeans(model6, list(pairwise ~ year), adjust ="bonferroni")
# 
# summary(model6)
# 
# ## QQ residual plot
# qqnorm(resid(model6))
# 
# 
# 
# 
# 
##______________________________________________________________
## Simple linear model
# lm(damage ~ treat+0, data = data)
# 
# lm(damage ~ treat + factor(year), data = datDamMod)
# summary(lm(damage ~ treat, data = data))
# 
# summary(lm(damage ~ treat + factor(year), data = data))
# 
# summary(lm(damage ~ treat + factor(year)+factor(block), data = data))
# ## Simple linear model wrapped in anova
# anova(lm(damage ~ treat + factor(year), data = data))
##______________________________________________________________
## LMER
# treat as a random-effect intercept and Year as both
# a fixed- and random-effect slope
# 
# 
# broom::tidy(model1) %>%
#   filter(term == "treat")
# 
# 
# 
# ## Treat as a whole group indicator model1<-glmer(damage ~ treat + (1 | year)+}, data = data, family = "poisson")
# summary(model2)
# 
# 
# sumTab$year <- as.factor(sumTab$year) 
# sumTab$meanDamage <- as.integer(sumTab$meanDamage)
# 
# 
# model1Mena<-glmer(meanDamage ~ treat + (treat | year), data = sumTab, family = "poisson")
# summary(model1Mena)
# ?isSingular
# isSingular(x, tol = 1e-4)
# 
