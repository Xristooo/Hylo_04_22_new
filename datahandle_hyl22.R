## Read data, ns = new session
library(dplyr)
library(ggplot)
setwd("E:ROBOT/R_kody/2022_Hylobius_review/Hylob_04_22_new")
getwd()
## ns= new (data)set
ns <-read.csv("dataFinal.csv")
nsMort <-read.csv("dataMORT_Final.csv")
str(ns)

## Pivot table, to creat a summary Damage for each SEEDLING - kontingencna, damage za jednu sadenicu spolu, Pripojit k celkovym udajom ku kazdemu roku
## Sluzi to na vyradenie sadenic, ktore maju poskodenie menej ako 50 a sú mrtve
## tvorba unikatneho id pre kazdu site
ns$uniqueIDsite <- paste(ns$uniqueID, ns$site )
seedUNQ <- ns %>% group_by(uniqueIDsite) %>%
  summarise(sumDamage = sum(damage),
            sumMort = sum(mort_binar))
## Join the seedUNQ to ns dataset outer join in R using merge() 
ns_1<- merge(x=ns,y=seedUNQ,by="uniqueIDsite",all=TRUE)
View(ns_1)


## check how many seedlings 
DY <- ns_1 %>% group_by(year) %>% 
          summarise(damageN = n())

## Data for DAMAGE assessment - CLEARCUT
dmg <- filter(ns_1, site== "clearcut") 

dmg_Count2 <- dmg_02 %>% group_by(year, treat) %>% 
  summarise(damageN = n())

names(dmg)

## what is in column notes 
unique(dmg$notes)
## replace seedlings with damage > 60 &  notes="other"
dmg$notes[dmg$notes == "other"] <-  "Hylobius"

## NEW_COLUMN all the data where: damage = 0
dmg$ZeroDam <- ifelse(dmg$damage==0, 0,1)

# levels(factor(dmg$treat))
# dmg$damage <- if(dmg$treat=="wax F" && dmg$damage > 1900){
# dmg$damage==1800  }


###################### Final DAMAGE dataset#################################
dmg_01 <- filter(dmg, ZeroDam >0)
dmg<-dmg_01 ## rename for easier handeling in models
write.csv(dmg, "data_04_22/filtr_Final_dat_04_22.csv")


## Filter CHEMICAL and DAMAGE == 0
## Chemical nebudeme vyhodnocovat
dmg_02 <- dmg %>%
  filter(treat!="chemical")
dmg_02

dm <- dmg_02



##___________________TRANSFORMATION__________________________________________________________________
levels(factor(dm$treat))

bestNormalize(dm$damage) # 1. Box-Cox: 1.1486, 2. Yeo-Johnson: 1.1967
yeoT <- yeojohnson(dm$damage) 
bxcx <- boxcox(dm$damage)
arc <-arcsinh_x(dm$damage)


par(mfrow=c(2,2))
hist(bxcx$x.t, main = "Box-Cox transformation")
hist(arc$x.t, main = "arcsihn transformation")
hist(yeoT$x.t, main = "Yeo-Johnson transformation")
hist(dm_01$damage, main = "No transformation")

## combine df, original + transform
## x original data - preserve to check when joining dfs,x.t transform
## x remove after combine frames and the original data match
## Box-Cox
bxcxT <-data.frame(bxcx$x.t, bxcx$x)
## OrderNorm
arc <-data.frame(arc$x.t, arc$x)

## Include transformed 
dm_01 <- data.frame(dm, bxcxT)
dm_01 <- data.frame(dm_01, arc)



## Final data without CHEMICAL
write.csv(dm_01, "data_04_22/filtr_Final_dat_04_22_NoChemTransform.csv", row.names = F)








