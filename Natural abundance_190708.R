### 13C natural abundance in cyanos and plant leaves from C3 and C4 plants from 3 sites
# Eva Dettweiler-Robinson
# 2017-2-9
library(visreg)
library(reshape2)
library(data.table)
library(dplyr)
library(emmeans)
library(ggplot2)
library(lme4)

library(car)
# library(ggpubr)
library(MuMIn)
library(tidyr)
dodge <- position_dodge(width=0.25)
std <- function(x) sd(x)/sqrt(length(x))


# # CSI data --------------------------------------------------------------



setwd("/Users/eva_stricker/Documents/Post_Doc/Fungal_Loop/Natural_abundance_13C")


# bring in data
Data <- read.table("Natural_abundance_13C_190615.txt", header = T)
# Natural_abundance_13C_190615.txt"

str(Data) # check data for whether it read things in as numeric/factors/etc.

Data <- Data[which(Data$Species!="NaN" & !is.na(Data$d13C) & Data$remove != "yes"),]
# I don't know wtf ID1 is - it doesn't seem to perfectly match up to the combination of site-species-ID-distance

Data$Site= ordered( as.character(Data$Site), levels = c("MOAB",  "SEV",  "JOR"))
Data$Species= ordered( as.character(Data$Species), levels = c( "BOGR","PLJA", "BOER",   "ACHY","GUSA"  ))
# Data[Data$Site=="SEV","Site"] <- "Sevilleta"

# to make UTEP values match UNM values, based on what we know so far.
Data[which(Data$Where_processed == "UTEP"), "d13C_new"] <- Data[which(Data$Where_processed == "UTEP"), "d13C"]-2.5
Data[which(Data$Where_processed == "UNM"), "d13C_new"] <- Data[which(Data$Where_processed == "UNM"), "d13C"]
colnames(Data)

# make a useful variable that gives each plant a unique ID (rather than naming plants 1-12 from each site)
Data$unique <- paste(Data$Site, Data$Species, Data$ID, sep = ".")

write.csv(Data, "SI_data_200224.csv")

# CN ratio ----------------------------------------------------------------


Data_agN <- aggregate(Data$X.N, 
                      list(Data$Site,  Data$Species, Data$unique, Data$type, Data$distance), mean)

Data_agC <- aggregate(Data$X.C, 
                      list(Data$Site,  Data$Species, Data$unique, Data$type, Data$distance), mean)
Data_agC$x

CN_ag <- cbind (Data_agN, Data_agC[,6])
head(CN_ag)


names(CN_ag)[1] <- "Site"
names(CN_ag)[2] <- "Species"
names(CN_ag)[3] <- "unique"
names(CN_ag)[4] <- "type"
names(CN_ag)[5] <- "distance"
names(CN_ag)[6] <- "N"
names(CN_ag)[7] <- "C"

CN_ag$C

CN_ag_r <- CN_ag[which(CN_ag$type=="cyano" & CN_ag$distance != 10),]
plot(CN_ag_r$Site, CN_ag_r$N)
CN_ag_r$CN <- (CN_ag_r$C/12)/(CN_ag_r$N/14)
sort(CN_ag_r$CN)

CN_ag_r$Site.Sp <- paste(CN_ag_r$Site, CN_ag_r$Species, sep = ".")
hist(log(CN_ag_r$CN))
CN_ag_r$distance <- as.factor(CN_ag_r$distance)
moda <- lmer(CN ~ (Species+distance)^2 + (1|Site) +(1|unique), data = CN_ag_r, REML=FALSE)
hist(resid(moda))
Anova(moda, type = 3)
em<-emmeans(moda, ~Species)
pairs(em, adjust ="fdr")
visreg(moda, "distance", by = "Site.Sp")

CN_ag_r[CN_ag_r$CN>25,]

plot(factor(CN_ag_r$Site.Sp), CN_ag_r$CN)
### This doesn't make any sense - the site/species with more nostoc/scytonema have lower C:N??? 


# 13C ---------------------------------------------------------------------

Data[is.na(Data$d13C_new ),]

#### 13C
Data2 <- Data[!is.na(Data$d13C_new),]
head(Data2)
Data_ag <- aggregate(Data2$d13C_new, 
                     list(Data2$Site, Data2$Species, Data2$ID,  Data2$type, Data2$distance), mean)

names(Data_ag)[1] <- "Site"
names(Data_ag)[2] <- "Species"
names(Data_ag)[3] <- "ID"
names(Data_ag)[4] <- "type"
names(Data_ag)[5] <- "distance"
names(Data_ag)[6] <- "d13C"
# write.csv(Data_ag, "Data_ag.csv")

### How does %C relate to 13C
# Cy <- Data[Data$type == "cyano" & Data$distance == "0",]
# moda <- lmer(d13C ~ (X.C+Species+Site)^3 + (1|unique), data = Cy)
# Anova(moda, type =3)
# # average over our replicates to get a single value
# Data_ag <- aggregate(Data$d13C_new, 
#           list(Data$Site, Data$Species, Data$unique, Data$type, Data$distance), mean)
# # 316 observations
 head(Data_ag)
# 
# 

# 


# Q1 - distance = 0 vs plant ----------------------------------------------



# Q1 Do cynobacteria that are next to plants vary with plant 13C? 
# hypothesis is that if biocrusts take up plant carbon, they should be reflect variability in plant d13C.

# reshape the data to make analysis easier

Data_ag_0<- Data_ag[Data_ag$distance == "0",]
summary(Data_ag_0)

Data2 = melt(Data_ag_0, id.vars = c("Site", "Species", "ID", "distance", "d13C"), 
             measure.vars = c("type"))
head(Data2)

Data3 = dcast(Data2, Site +Species + ID + distance ~ value, fun.agg = mean, value.var = "d13C")
tail(Data3)

# subset to just the distance = 0 and type = leaf (which is coded with distance = 0 in the spreadsheet))
Dist0_full <- Data3[Data3$distance == "0",]
# to visualize them all together
plot(Dist0_full$cyano~ Dist0_full$leaf, col = as.numeric(Dist0_full$Species), pch = as.numeric(Dist0_full$Site), xlab = "plant d13C", ylab = "cyano d13C")
Dist0_full$sisp <- paste(Dist0_full$Site, Dist0_full$Species)

Dist0_full_nop<- Dist0_full
head(Dist0_full_nop)
Dist0_full_nop$sisp <- paste(Dist0_full_nop$Site, Dist0_full_nop$Species)



leaf_reg<- ggplot( Dist0_full2, aes(x=leaf, y=cyano)) +
  theme_bw() +   theme(  axis.title.x = element_text(vjust=-0.35),
                         axis.title.y = element_text(vjust=0.35) ,
                         strip.text.x = element_blank(),
                         axis.title = element_text(size = 12),
                         legend.position = "none",
                         legend.title=element_text(size=11),
                         axis.text = element_text(size = 12),
                         legend.text=element_text(size=10),
                         axis.ticks.x=element_blank(),
                         legend.key = element_rect(fill = "white"),
                         legend.background = element_rect(fill = "white"),
                         panel.grid.major = element_line(colour = "white"),
                         panel.grid.minor = element_blank()  ) +
  
  #   geom_line(stat = "identity", aes(group = Species), size = .5,  lty = 3) +
  scale_shape_manual(values=c( 17,19, 15))+
  scale_color_manual(values=c("blue", "grey","black",  "palegreen2", "seagreen"))+
  geom_point( aes(leaf, cyano, shape = Site, col = Species, alpha = .5) ,size = 1, data =Dist0_full2) +
  
  ylim (-29.5, -15)+
   xlim(-29.5, -15)+
    #geom_errorbar(data=Data4,  aes(y = mean13C_cyano,   ymin = mean13C_cyano-se13C_cyano, ymax = mean13C_cyano+se13C_cyano, width = 0), size = .6)+
  #geom_errorbarh(data=Data4,  aes(x = mean13C_leaf,   xmin = mean13C_leaf-se13C_leaf, xmax = mean13C_leaf+se13C_leaf, height = 0), size = .6)+
  # guides(fill = guide_legend(override.aes = list(linetype = 0)))+
  # geom_abline(slope = 1, intercept =0, lty = 3)+
  geom_smooth(method = "lm", se=F,aes( group = Species, col = Species)) +
  labs(x = expression("Observed leaf "*delta^{13}* "C (\u2030)") , y = expression("Observed cyanobacteria "*delta^{13}* "C (\u2030)") )+
  annotate("text", x = -29, y = -26, label = "e.")

# scale_x_discrete(expand=c(0.3,0), drop=FALSE)
leaf_reg


# scale_color_manual(values=c("blue","blue", "grey80", "grey80", "grey30","grey30",  "light green", "light green","dark green","dark green"))+
  

Dist0_full$sisp <- paste(Dist0_full$Site, Dist0_full$Species)

#update!!!
Dist0_full2<-Dist0_full [which(Dist0_full$sisp != "JOR BOER" & Dist0_full$sisp != "SEV PLJA"),]
head(Dist0_full2)
### note, model does NOT include JOR BOER or SEV PLJA because the plant values are not paired with the cyano values
mod1 <- lmer(cyano~leaf*Species +(1|Site), data = Dist0_full2)
Anova(mod1, type = 3)

emtrends(mod1, ~Species, var= "leaf")

Boer <- Dist0_full[Dist0_full$Species == "BOER"  ,]
mod1 <- lm(cyano ~ leaf*Site, data = Boer)
Anova(mod1,type=3)
library(emmeans)
emmeans(mod1, ~leaf*Site)



Bogr <- Dist0_full[Dist0_full$Species == "BOGR"  ,]
mod1 <- lm(cyano ~ leaf*Site, data = Bogr)
Anova(mod1,type=3)

# to update
PLJA <- Dist0_full[Dist0_full$Species == "PLJA"  & Dist0_full$Site != "SEV",]
mod1 <- lm(cyano ~ leaf, data = PLJA)
Anova(mod1,type=3)

GUSA <- Dist0_full[Dist0_full$Species == "GUSA"  ,]
mod1 <- lm(cyano ~ leaf*Site, data = GUSA)
Anova(mod1,type=3)


ACHY <- Dist0_full[Dist0_full$Species == "ACHY"  ,]
mod1 <- lm(cyano ~ leaf*Site, data = ACHY)
Anova(mod1,type=3)




detach("package:reshape2", unload=TRUE)
detach("package:data.table", unload=TRUE)
detach("package:dplyr", unload = TRUE)
library(data.table)
# library(reshape2)
library(dplyr)
Data_merge_ag <- Data_ag_0 %>% group_by(Site, Species, type, distance) %>%
  summarise(mean13C = mean(d13C), se13C=std(d13C))
head(Data_merge_ag)

Data4 = dcast(setDT(Data_merge_ag), Site + Species + distance ~ type, 
              fun.agg = mean, value.var = c('mean13C', 'se13C'))
Data4$sisp <- paste(Data4$Site, Data4$Species)


Data4$Species= ordered( as.character(Data4$Species), levels = c(  "BOGR","PLJA", "BOER","ACHY",  "GUSA" ))
Data4$Site= ordered( as.character(Data4$Site), levels = c("MOAB",  "SEV",  "JOR"))


### for plotting, we DO include JOR BOER and SEV PLJA for completeness sake

Data_5 <- Data4[which(Data4$sisp!="SEV PLJA" ),]
### To update with new data!!!
# 
# leaf_reg <- ggplot( Data_5, aes(x=mean13C_leaf, y=mean13C_cyano)) +
#   theme_bw() +   theme(  axis.title.x = element_text(vjust=-0.35),
#                          axis.title.y = element_text(vjust=0.35) ,
#                          strip.text.x = element_blank(),
#                          axis.title = element_text(size = 12),
#                          legend.position = "none",
#                          legend.title=element_text(size=11),
#                          axis.text = element_text(size = 12),
#                          legend.text=element_text(size=10),
#                          axis.ticks.x=element_blank(),
#                          legend.key = element_rect(fill = "white"),
#                          legend.background = element_rect(fill = "white"),
#                          panel.grid.major = element_line(colour = "white"),
#                          panel.grid.minor = element_blank()  ) +
#   
#   geom_line(stat = "identity", aes(group = Species), size = .5,  lty = 3) +
#   scale_shape_manual(values=c( 17, 19, 15))+
#   scale_color_manual(values=c("blue", "grey80", "grey30", "light green", "dark green"))+
#   
#   ylim (-26.7, -16.9)+ xlim(-28.7, -14.9)+
#   geom_point( aes(mean13C_leaf, mean13C_cyano, shape = Site, col = Species ) ,size = 3, data =Data_5) +
#   geom_errorbar(data=Data4,  aes(y = mean13C_cyano,   ymin = mean13C_cyano-se13C_cyano, ymax = mean13C_cyano+se13C_cyano, width = 0), size = .6)+
#   geom_errorbarh(data=Data4,  aes(x = mean13C_leaf,   xmin = mean13C_leaf-se13C_leaf, xmax = mean13C_leaf+se13C_leaf, height = 0), size = .6)+
#   
#   guides(fill = guide_legend(override.aes = list(linetype = 0)))+
# # geom_abline(slope = 1, intercept =0, lty = 3)+
#   labs(x = expression("Observed leaf "*delta^{13}* "C (\u2030)") , y = expression("Observed cyanobacteria "*delta^{13}* "C (\u2030)") )
# # scale_x_discrete(expand=c(0.3,0), drop=FALSE)
# leaf_reg
# 
cairo_pdf( "leaf_new_phyt.pdf", width = 4.4, height = 3 )
leaf_reg
dev.off()
  

# Q2 - by distance --------------------------------------------------------


# Q2 compare cyano 13C by distance from plant and species. 
# subset to only the cyano filament samples (exclude plants) and for now only look at 0 and 25

Data_ag$distance2 <- as.factor(Data_ag$distance)
Data_ag[Data_ag$Species=="PLJA",]
Cyano2<-Data_ag[which(Data_ag$type =="cyano" &  Data_ag$Species != "NaN" & Data_ag$distance != "10"),]
# check sample size 
Cyano2$Site.sp <- paste(Cyano2$Site, Cyano2$Species, sep = ".")
Cyano <- Cyano2[which(Cyano2$Site.sp!= "JOR.GUSA" & !is.na(Cyano2$d13C)),]
head(Cyano
     )

#### combine all the interspace values and just compare the species against those... I don't like it.
Cyano$microsite <- paste(Cyano$Species, Cyano$distance)
unique(Cyano$microsite)
Cyano[Cyano$distance=="25","microsite"]<- "I 25"
xtabs(~ID+Site,Cyano)
# modJenn<- lmer(d13C ~ microsite +(1|ID), data = Cyano)
# Anova(modJenn)
# em<-emmeans(modJenn, ~microsite)
# library(multcompView)
#   CLD(em, adjust="fdr")
# we need to account for the fact that these are measures of the same experimental unit- so we add a random effect that is the plant ID
# it will be useful to treat distance as a factor for plotting; because we have 2 points we still only use 1 df; same as if we treat it as continuous

Cyano$ID <- as.factor(Cyano$ID)

### 1/4/19 AFTER MEETING WITH JENN

# Q1. Does 13C cyano signature differ  by distance from diff plant species?

head(Cyano)
unique(Cyano$Site.sp)
Cyano[Cyano$Site.sp=="SEV.BOGR",]

Cyano$unique <- as.factor(paste(Cyano$Species, Cyano$ID))


mod1 <- lmer(d13C~(distance2*Species)  + (1|Site/unique), data = Cyano)
Anova(mod1, type = 3)

pairs(emmeans(mod1, ~distance2|Species), adjust = "fdr")

# Yes!
# going to move 10cm analysis to supplementals because it's incomplete... no sev plja
PLJAC <- Cyano[which(Cyano$Species == "PLJA" & Cyano$type =="cyano" & Cyano$distance!="10"&   !is.na(Cyano$d13C)),]
mod1 <- lmer(d13C~(distance2)*Site  +(1|ID), data = PLJAC)
Anova(mod1)
pairs(emmeans(mod1, ~distance2), adjust = "fdr") 

 PLJA10 <- Data_ag[which(Data_ag$Species == "PLJA" & Data_ag$Site!="SEV"& Data_ag$type =="cyano"  &  !is.na(Data_ag$d13C)),]
# summary(PLJA)
 mod10 <- lmer(d13C~(distance2)  +(1|ID), data = PLJA10)
 Anova(mod10, type = 3)
 lmerTest:: lmer(d13C~(distance2)  +(1|ID), data = PLJA10)
 anova(mod1)
# head(PLJA)
pairs(emmeans(mod1, ~distance2), adjust = "fdr")

def <- summary(emmeans(mod1, ~distance2*Site)) 
xyzplja <- as.data.frame(def[c("distance2","Site", "emmean", "SE", "df", "lower.CL", "upper.CL")])

# ggplot(PLJA, aes(x=distance2, y = d13C), group = Site, color = Site)+geom_point(aes(color = Site))

plja1<-ggplot( xyzplja, aes(x=distance2, y=emmean), group = Site) +
  theme_bw() +   theme(  axis.title.x = element_text(vjust=-0.35),
                         axis.title.y = element_text(vjust=0.35) ,
                         strip.text.x = element_blank(),
                         axis.title = element_text(size = 12),
                         legend.position = "none",
                         legend.title=element_text(size=11),
                         axis.text = element_text(size = 12),
                         legend.text=element_text(size=10),
                         axis.ticks.x=element_blank(),
                         legend.key = element_rect(fill = "white"),
                         legend.background = element_rect(fill = "white"),
                         panel.grid.major = element_line(colour = "white"),
                         panel.grid.minor = element_blank()  ) +
  geom_errorbar(position = dodge, data=xyzplja,  aes(y = emmean, group = Site,  ymin = lower.CL, ymax = upper.CL, width = 0.2), size = .8)+
  geom_line(position = dodge,stat = "identity", aes(group = Site), size = 1,  lty = 2) +
  scale_shape_manual(values=c( 17, 19), name="Site", breaks=c( "MOAB","JOR1"),    labels=c("MOAB", "JOR"))+
  ylim (-26.7, -15.9)+
  geom_point(position = dodge, aes(distance2, emmean, group = Site, shape = Site, size = distance2), col = "grey", data = xyzplja) +
#  guides(fill = guide_legend(override.aes = list(linetype = 0)))+
  labs(x = "Distance from plant (cm)",  y = expression("Cyanobacteria "*delta^{13}* "C (\u2030)") )+
  scale_x_discrete(expand=c(0.3,0), drop=FALSE)+
  annotate("text", x = 1.1, y = -17.2, label = "P. jamesii", fontface = 'italic')+
  annotate("text", x = .9, y = -26, label = "b.")



plja1


BOER <- Data_ag[which(Data_ag$Species == "BOER" & Data_ag$type =="cyano" &   Data_ag$distance != "10" & !is.na(Data_ag$d13C)),]

mod1 <- lmer(d13C~(distance2*Site)  +(1|ID), data = BOER)
Anova(mod1, type = 3)

def <- summary(emmeans(mod1, ~distance2|Site)) 
xyzboer <- as.data.frame(def[c("distance2","Site",  "emmean", "SE", "df", "lower.CL", "upper.CL")])
boer1<-ggplot( xyzboer, aes(x=distance2, y=emmean), group = Site) +
  theme_bw() +   theme(  axis.title.x = element_text(vjust=-0.35),
                         axis.title.y = element_text(vjust=0.35) ,
                         strip.text.x = element_blank(),
                         axis.title = element_text(size = 12),
                         legend.position = "none",
                         legend.title=element_text(size=11),
                         axis.text = element_text(size = 12),
                         legend.text=element_text(size=10),
                         axis.ticks.x=element_blank(),
                         legend.key = element_rect(fill = "white"),
                         legend.background = element_rect(fill = "white"),
                         panel.grid.major = element_line(colour = "white"),
                         panel.grid.minor = element_blank()  ) +
  geom_errorbar(position = dodge, data=xyzboer,  aes(y = emmean, group = Site,  ymin = lower.CL, ymax = upper.CL, width = 0.2), size = .8)+
  geom_line(position = dodge,stat = "identity", aes(group = Site), size = 1,  lty = 2) +
  scale_shape_manual(values=c(  19,15), name="Site", breaks=c( "SEV2","JOR1"),    labels=c("SEV", "JOR"))+
  ylim (-26.7, -15.9)+
  geom_point(position = dodge, aes(distance2, emmean, group = Site, shape = Site, size = distance2), col = "black", data = xyzboer) +
#  guides(fill = guide_legend(override.aes = list(linetype = 0)))+
  labs(x = "Distance from plant (cm)",  y = expression("Cyanobacteria "*delta^{13}* "C (\u2030)") )+
  scale_x_discrete(expand=c(0.3,0), drop=FALSE)+
  annotate("text", x = 1.1, y = -17.2, label = "B. eriopoda", fontface = 'italic')+
  annotate("text", x = .9, y = -26, label = "c.")



boer1




BOGR <- Cyano[which(Cyano$Species == "BOGR" & Cyano$type =="cyano"  &   Cyano$distance != "10"& !is.na(Cyano$d13C)),]

mod1 <- lmer(d13C~(distance2*Site)  +(1|ID), data = BOGR)
Anova(mod1, type = 3)


def <- summary(emmeans(mod1, ~distance2|Site)) 
xyzbogr <- as.data.frame(def[c("distance2","Site",  "emmean", "SE", "df", "lower.CL", "upper.CL")])


bogr1<-ggplot( xyzbogr, aes(x=distance2, y=emmean), group = Site) +
  theme_bw() +   theme(  axis.title.x = element_text(vjust=-0.35),
                         axis.title.y = element_text(vjust=0.35) ,
                         strip.text.x = element_blank(),
                         axis.title = element_text(size = 12),
                         legend.position = "none",
                         legend.title=element_text(size=11),
                         axis.text = element_text(size = 12),
                         legend.text=element_text(size=10),
                         axis.ticks.x=element_blank(),
                         legend.key = element_rect(fill = "white"),
                         legend.background = element_rect(fill = "white"),
                         panel.grid.major = element_line(colour = "white"),
                         panel.grid.minor = element_blank()  ) +
  geom_errorbar(position = dodge, data=xyzbogr,  aes(y = emmean, group = Site,  ymin = lower.CL, ymax = upper.CL, width = 0.2), size = .8)+
  geom_line(position = dodge,stat = "identity", aes(group = Site), size = 1,  lty = 2) +
   scale_shape_manual(values=c( 17, 19), name="Site", breaks=c( "MOAB2","SEV2"),    labels=c( "MOAB","SEV"))+
  ylim (-26.7, -15.9)+
  geom_point(position = dodge, aes(distance2, emmean, group = Site, shape = Site, size = distance2), col = "blue", data = xyzbogr) +
#  guides(fill = guide_legend(override.aes = list(linetype = 0)))+
  labs(x = "Distance from plant (cm)",  y = expression("Cyanobacteria "*delta^{13}* "C (\u2030)") )+
  scale_x_discrete(expand=c(0.3,0), drop=FALSE)+
  annotate("text", x = 1.1, y = -17.2, label = "B. gracilis", fontface = 'italic')+
  annotate("text", x = .9, y = -26, label = "a.")


bogr1



#sup
GUSA10 <- Data_ag[which(Data_ag$Species == "GUSA" & Data_ag$type =="cyano" & Data_ag$Site != "JOR" & !is.na(Data_ag$d13C)),]
mod1 <- lmer(d13C~(distance2*Site)  +(1|ID), data = GUSA)
Anova(mod1, type = 3)
pairs(emmeans(mod1, ~distance2|Site), adjust = "fdr") 


#main
# GUSA <- Data_ag[which(Data_ag$Species == "GUSA" & Data_ag$type =="cyano" & Data_ag$distance2!="10"& Data_ag$Site!="JOR"& !is.na(Data_ag$d13C)),]
# xtabs( ~distance2+Site,GUSA)
# GUSA[GUSA$Site=="JOR",]
# 
# mod1 <- lmer(d13C~(distance2*Site)  +(1|ID), data = GUSA)
# Anova(mod1, type = 3)
# pairs(emmeans(mod1, ~distance2), adjust = "fdr") 
# # 10cm and 25cm are the same

GUSA <- Data_ag[which(Data_ag$Species == "GUSA" & Data_ag$type =="cyano" & Data_ag$distance2!="10"&  !is.na(Data_ag$d13C)),]
mod1 <- lmer(d13C~(distance2*Site)  +(1|ID), data = GUSA)
Anova(mod1, type = 3)
pairs(emmeans(mod1, ~distance2), adjust = "fdr") 
# 10cm

def <- summary(emmeans(mod1, ~distance2|Site)) 
xyzgusa <- as.data.frame(def[c("distance2","Site",  "emmean", "SE", "df", "lower.CL", "upper.CL")])
# xyzgusa$species = "GUSA"
gusa1<-ggplot( xyzgusa, aes(x=distance2, y=emmean), group = Site) +
  theme_bw() +   theme(  axis.title.x = element_text(vjust=-0.35),
                         
                         axis.title.y = element_text(vjust=0.35) ,
                         strip.text.x = element_blank(),
                         axis.title = element_text(size = 12),
                         legend.position = "none",
                         legend.title=element_text(size=11),
                         axis.text = element_text(size = 12),
                         legend.text=element_text(size=10),
                         axis.ticks.x=element_blank(),
                         legend.key = element_rect(fill = "white"),
                         legend.background = element_rect(fill = "white"),
                         panel.grid.major = element_line(colour = "white"),
                         panel.grid.minor = element_blank(),
                         
  ) +
  geom_errorbar(position = dodge, data=xyzgusa,  aes(y = emmean, group = Site,  ymin = lower.CL, ymax = upper.CL, width = 0.2), size = .8)+
  geom_line(position = dodge,stat = "identity", aes(group = Site), size = 1,  lty = 2) +
   scale_shape_manual(values=c(17, 19, 15), name="Site", breaks=c("SEV2", "MOAB2", "JOR2"),
                     labels=c("MOAB","SEV", "JOR"))+
  ylim (-26.7, -15.9)+
  geom_point(position = dodge, 
             aes(distance2, emmean, group = Site, shape = Site, size = distance2), col = "seagreen", data = xyzgusa) +
 # guides(fill = guide_legend(override.aes = list(linetype = 0)))+
  labs(x = "Distance from plant (cm)",  y = expression("Cyanobacteria "*delta^{13}* "C (\u2030)") )+
  scale_x_discrete(expand=c(0.3,0), drop=FALSE)+
  annotate("text", x = 1.5, y = -17.2, label = "G. sarothrae", fontface = 'italic')+
annotate("text", x = .9, y = -26, label = "e.")

gusa1
# supp
ACHY <- Data_ag[which(Data_ag$Species == "ACHY" & Data_ag$type =="cyano" & !is.na(Data_ag$d13C)),]

mod1 <- lmer(d13C~(distance2*Site)  +(1|ID), data = ACHY)
Anova(mod1, type = 3)
pairs(emmeans(mod1, ~distance2|Site), adjust="fdr") 

# main
ACHY <- Data_ag[which(Data_ag$Species == "ACHY" & Data_ag$type =="cyano" & Data_ag$distance2!="10" & !is.na(Data_ag$d13C)),]

mod1 <- lmer(d13C~(distance2*Site)  +(1|ID), data = ACHY)
Anova(mod1, type = 3)
pairs(emmeans(mod1, ~distance2|Site), adjust="fdr") 
def <- summary(emmeans(mod1, ~distance2|Site)) 
xyzachy <- as.data.frame(def[c("distance2","Site",  "emmean", "SE", "df", "lower.CL", "upper.CL")])

achy1<-ggplot( xyzachy, aes(x=distance2, y=emmean), group = Site) +
  theme_bw() +   theme(  axis.title.x = element_text(vjust=-0.35),
                         axis.title.y = element_text(vjust=0.35) ,
                         strip.text.x = element_blank(),
                         axis.title = element_text(size = 12),
                         legend.position = "none",
                         legend.title=element_text(size=11),
                         axis.text = element_text(size = 12),
                         legend.text=element_text(size=10),
                         axis.ticks.x=element_blank(),
                         legend.key = element_rect(fill = "white"),
                         legend.background = element_rect(fill = "white"),
                         panel.grid.major = element_line(colour = "white"),
                         panel.grid.minor = element_blank(),
                         
  ) +
  geom_errorbar(position = dodge, data=xyzachy,  aes(y = emmean, group = Site,  ymin = lower.CL, ymax = upper.CL, width = 0.2), size = .8)+
  geom_line(position = dodge,stat = "identity", aes(group = Site), size = 1,  lty = 2) +
   scale_shape_manual(values=c(17, 19), name="Site", breaks=c("SEV2", "MOAB2"),
                      labels=c( "MOAB","SEV"))+
  ylim (-26.7, -15.9)+
  geom_point(position = dodge, aes(distance2, emmean, group = Site, shape = Site, size = distance2), col = "palegreen2", data = xyzachy) +
#  guides(fill = guide_legend(override.aes = list(linetype = 0)))+
  labs(x = "Distance from plant (cm)",  y = expression("Cyanobacteria "*delta^{13}* "C (\u2030)") )+
  scale_x_discrete(expand=c(0.3,0), drop=FALSE)+
  annotate("text", x = 1.5, y = -17.2, label = "A. hymenoides", fontface = 'italic')+
annotate("text", x = .9, y = -26, label = "d.")

achy1



gA <- ggplotGrob(boer1)
gB <- ggplotGrob(bogr1)
gC <- ggplotGrob(gusa1)
gD <- ggplotGrob(achy1)
gE <- ggplotGrob(plja1)
# 
# 
# 
gA$widths <- gB$widths
gC$widths <- gB$widths
gD$widths <- gB$widths
gE$widths <- gB$widths
# 
# gA$heights <- gB$heights
# gC$heights <- gB$heights
# 
# 
library(gridExtra)
library(grid)
# 
cairo_pdf(paste("Fig_sites", Sys.Date(), ".pdf", sep = ""), height = 5, width = 8.1)
grid.arrange( gB, gE, gA , gD, gC, ncol = 3) 
# # grid.arrange(gF, gG, gE, nrow=3)
dev.off()




# DOC ---------------------------------------------------------------------

DOC <- read.table(file.choose(), header = T)
# JORSEV_Spring_monsoon_2017_N_DOC_assays_180122.txt

master <- read.table(file.choose(), header = T)
# JORSEV_2017_Master_raceway_file_181226.txt

DOC_m <- left_join(DOC,master, by=c('Raceway_ID'='Raceway_ID'))

DOC_M2<- DOC_m[,c(1:2,4, 8:ncol(DOC_m))]

DOC_wide <- spread(DOC_M2, key = microsite, value =(initial))
head(DOC_wide)
plot(DOC_wide$crust_end, DOC_wide$plant_end, col = DOC_wide$site.sp)

DOC_M2$site.sp <- as.factor(paste(DOC_M2$Site, DOC_M2$Plant_species, sep = "."))

xtabs(~microsite+Site, DOC_M2)


DOC_M2_sub <- DOC_M2[DOC_M2$mesh_type=="none" , ]
xtabs(~microsite+Site, DOC_M2_sub)


mod1 <- lmer((initial)~(microsite*Plant_species)  +(1|Site)+(1|Raceway_ID), data = DOC_M2_sub)
Anova(mod1, type = 3)
emmeans(mod1, ~microsite)
(10.675258-9.498914)/9.498914
# 12% higher DOC in plant end


DOC_M2_sub$Site= ordered( as.character(DOC_M2_sub$Site), levels = c("MOAB",  "SEV",  "JOR"))
DOC_M2_sub$Site2 <- as.character(DOC_M2_sub$Site)

DOC_M2_sub[DOC_M2_sub$Site=="SEV", "Site2"]<- "Sevilleta"
DOC_M2_sub[DOC_M2_sub$Site=="JOR", "Site2"]<- "Jornada"

mod1 <- lm((initial)~(distance*Plant_species*Site2), data = DOC_M2_sub)
xyz <- summary(emmeans(mod1, ~distance*Plant_species*Site2))
xyz$Site2= ordered( as.character(xyz$Site2), levels = c("MOAB",  "Sevilleta",  "Jornada"))


DOC <- ggplot(xyz, aes (x = distance, y = emmean,width = .5, fill = Plant_species))+facet_grid(Plant_species~Site2)+
theme_bw() +  geom_bar(stat = "identity", position = position_dodge())+
  theme(  axis.title.x = element_text(vjust=-0.35),
                       axis.title.y = element_text(vjust=0.35) ,
                      
                       axis.title = element_text(size = 12),
                       legend.position = "none",
                       legend.title=element_text(size=11),
                       axis.text = element_text(size = 12),
                       legend.text=element_text(size=10),
                       axis.ticks.x=element_blank(),
                       legend.key = element_rect(fill = "white"),
                       legend.background = element_rect(fill = "white"),
                       panel.grid.major = element_line(colour = "white"),
                       panel.grid.minor = element_blank())+
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.2,
                position=position_dodge(.9))+
  labs(y = expression(paste("Dissolved organic C (", mu, "g g"^-1, "soil)")) , x = expression("Distance from plant (cm)") )+
  scale_fill_manual(values=c("grey30", "dark green"))+ylim(c(-2,28))

DOC

cairo_pdf("doc_microsite.pdf",width = 4, height = 3 )
DOC
dev.off()

DOC_M2_sub$distance <-as.character( DOC_M2_sub$microsite)

xtabs(~Plant_species + distance+Site,DOC_M2_sub)

DOC_M2_sub[DOC_M2_sub$distance=="crust_end","distance"] <- "25"
DOC_M2_sub[DOC_M2_sub$distance=="plant_end","distance"] <- "0"


DOC_M2_sub$distance <-as.factor( DOC_M2_sub$distance)
DOC_M2_sub$Site.sp <- paste(DOC_M2_sub$Site, DOC_M2_sub$Plant_species, sep = ".")
xtabs(~Site.sp + distance,DOC_M2_sub)


DOC_M2_ag <- DOC_M2_sub %>% group_by(Site, Plant_species, distance) %>%drop_na(initial)%>%
  summarise(meanDOC = mean(initial), sedoc=std(initial))

xtabs(~Site+Plant_species + distance,DOC_M2_ag)


Cyano_merge <- Cyano[, c(1,2, 6,7,8) ]

Cyano_merge_ag <- Cyano_merge %>% group_by(Site, Species, distance2) %>%drop_na(d13C)%>%
  summarise(mean13C = mean(d13C), se13C=std(d13C))
head(Cyano_merge_ag)

DOC_13C <- full_join(Cyano_merge_ag,DOC_M2_ag, by=c('Site'='Site',  'distance2'='distance', 'Species'='Plant_species'))
tail(DOC_13C)

mod1 <- lm(mean13C~meanDOC*Species, data = DOC_13C)
Anova(mod1, type = 3)
emmeans(mod1, ~Species, var = "meanDOC")

summary(mod1)
DOC_13C$sisp <- paste(DOC_13C$Site, DOC_13C$Species)


DOC_13C$Species= ordered( as.character(DOC_13C$Species),levels = c( "BOGR","PLJA", "BOER",   "ACHY","GUSA"  ))
DOC_13C$Site= ordered( as.character(DOC_13C$Site), levels = c("MOAB",  "SEV",  "JOR"))


DOC_reg <- ggplot( DOC_13C, aes(x=meanDOC, y=mean13C)) +
  theme_bw() +   theme(  axis.title.x = element_text(vjust=-0.35),
                         axis.title.y = element_text(vjust=0.35) ,
                         strip.text.x = element_blank(),
                         axis.title = element_text(size = 12),
                         legend.position = "none",
                         legend.title=element_text(size=11),
                         axis.text = element_text(size = 12),
                         legend.text=element_text(size=10),
                         axis.ticks.x=element_blank(),
                         legend.key = element_rect(fill = "white"),
                         legend.background = element_rect(fill = "white"),
                         panel.grid.major = element_line(colour = "white"),
                         panel.grid.minor = element_blank()  ) +
    
    geom_line(stat = "identity", aes(group = sisp), size = .5,  lty = 3) +
  scale_shape_manual(values=c( 17,19,  15))+
  scale_color_manual(values=c(  "light green", "dark green", "blue", "grey30", "dark green"))+
   ylim (-26.7, -16.9)+xlim(2,26)+
  geom_point( aes(meanDOC, mean13C, shape = Site, col = Species, size = distance2 ), data = DOC_13C) +
  geom_errorbar(data=DOC_13C,  aes(y = mean13C,   ymin = mean13C-se13C, ymax = mean13C+se13C, width = 0), size = .8)+
  geom_errorbarh(data=DOC_13C,  aes(x = meanDOC,   xmin = meanDOC-sedoc, xmax = meanDOC+sedoc, height = 0), size = .8)+

#   guides(fill = guide_legend(override.aes = list(linetype = 0)))+
  labs(x = expression(paste("Dissolved organic C (", mu, "g g"^-1, "soil)")) , y = expression("Observed cyanobacteria "*delta^{13}* "C (\u2030)") )
# scale_x_discrete(expand=c(0.3,0), drop=FALSE)
DOC_reg

cairo_pdf( "DOC_new_phyt.pdf", width = 4.4, height = 3 )
DOC_reg
dev.off()



mod1 <- lmer((initial)~(microsite*Plant_species)  +(1|Site)+(1|Raceway_ID), data = DOC_M2)
Anova(mod1, type =3)


# Chlorophyll -------------------------------------------------------------

Chlor<- read.table(file.choose(), header = T)
# JORSEVMOAB_chlor-2016_EDR_190502.txt
xtabs( ~Site+Species, Chlor)

head(Chlor)

Chlor_ag <- aggregate(Chlor$Chlor, 
                     list(Chlor$Site, Chlor$Site2,  Chlor$Species, Chlor$ID, Chlor$ID2, Chlor$Distance), mean)

names(Chlor_ag)[1] <- "Site"
names(Chlor_ag)[2] <- "Site2"
names(Chlor_ag)[3] <- "Species"
names(Chlor_ag)[4] <- "ID"
names(Chlor_ag)[5] <- "ID2"
names(Chlor_ag)[6] <- "Distance"
names(Chlor_ag)[7] <- "Chlor"


head(Chlor_ag)
Chlor_ag$Site.sp <- as.factor(paste(Chlor_ag$Site, Chlor_ag$Species, sep = "."))
Chlor_ag[Chlor_ag$Site.sp== "SEV.GUSA",]
Chlor_ag2<-Chlor_ag[ Chlor_ag$Distance != "10" ,]
Chlor_ag2$Distance <- as.factor(Chlor_ag2$Distance)

Cyano$Site2 <- as.character(Cyano$Site)
Cyano[Cyano$Site2=="JOR1","Site2"]<- "JOR"
Cyano[Cyano$Site2=="JOR2","Site2"]<- "JOR"
Cyano[Cyano$Site2=="SEV1","Site2"]<- "SEV"
Cyano[Cyano$Site2=="SEV2","Site2"]<- "SEV"
Cyano[Cyano$Site2=="MOAB1","Site2"]<- "MOAB"
Cyano[Cyano$Site2=="MOAB2","Site2"]<- "MOAB"
colnames(Cyano)
head(Cyano)



Cyano_merge <- Cyano[, c(2,3, 6,7,11) ]
Cyano_merge$Site.sp <- as.factor(paste(Cyano_merge$Site2, Cyano_merge$Species, sep = "."))

Cyano_merge[Cyano_merge$Site.sp == "SEV.GUSA",]
Chlor_ag2[Chlor_ag2$ID2 == "37",]
head(Chlor_ag2)

Chlor_ag2<-Chlor_ag2[,c(1,3,4,6,7,8)]
Chlor_ag2$ID2<-as.factor(Chlor_ag2$ID)
Cyano_merge$ID<-as.factor(Cyano_merge$ID)

Chlor_ag2[Chlor_ag2$Site.sp=="SEV.GUSA",]
head(Cyano_merge)
summary(Cyano_merge)
Alldata_sm <- left_join(Chlor_ag2,Cyano_merge, by=c('Site.sp'='Site.sp', 'Species' = 'Species', 'Site' = 'Site2', 
                                                    'ID2'=  'ID', 'Distance' = 'distance2'))

xtabs( ~Site+Species+Distance, Alldata_sm)


Alldata_sm[Alldata_sm$Site=="JOR",]

Alldata_sm[which(Alldata_sm$d13C<= -25),]
head(Alldata_sm)
Alldata_sm$unique<-paste(Alldata_sm$JOR, Alldata_sm$ID)
mod1 <- lmer(d13C ~ Chlor*Species  + (1|Site/unique), data=Alldata_sm)
Anova(mod1, type = 3)
# no interaction
summary(mod1)
emtrends(mod1, var= "Chlor")
emmeans(mod1, ~"Chlor", var = "Species")
Alldata_sm$Species= ordered( as.character(Alldata_sm$Species), levels = c("BOGR",  "PLJA",  "BOER", "ACHY", "GUSA"))

visreg(mod1,"Chlor", by = "Species",overlay = T, line=list(col=c("blue", "grey", "black", "light green", "dark green"), lwd=1),  points=list(col=c("blue", "grey", "black", "light green", "dark green")))


mod1 <- lmer(sqrt(Chlor)~(Distance*Species)  + (1|Site/ID2), data = Chlor_ag2)
Anova(mod1, type = 3)
# NO species x distance interaction; only effect of species.
xtabs(~Site.sp+Distance, Chlor_ag2)
cld.emmGrid(emmeans(mod1, ~Species))
mod1 <- lmer(sqrt(Chlor)~Distance*Site.sp+(1|ID2), data = Chlor_ag2)
.62^2
xyz <- summary(emmeans(mod1, ~Distance*Site.sp))
xyz$Site2 <- as.character(xyz$Site.sp)

xyz$Sp2 <- as.character(xyz$Site.sp)
xyz$mean <- xyz$emmean^2
xyz$up <- xyz$upper.CL^2
xyz$low <- xyz$lower.CL^2
xyz[xyz$Site.sp == "MOAB.ACHY", "Site2"]<- "Moab"
xyz[xyz$Site.sp == "MOAB.BOGR", "Site2"]<- "Moab"
xyz[xyz$Site.sp == "MOAB.PLJA", "Site2"]<- "Moab"
xyz[xyz$Site.sp == "MOAB.GUSA", "Site2"]<- "Moab"
xyz[xyz$Site.sp == "SEV.GUSA", "Site2"]<- "Sevilleta"
xyz[xyz$Site.sp == "SEV.ACHY", "Site2"]<- "Sevilleta"
xyz[xyz$Site.sp == "SEV.BOER", "Site2"]<- "Sevilleta"
xyz[xyz$Site.sp == "SEV.BOGR", "Site2"]<- "Sevilleta"
xyz[xyz$Site.sp == "SEV.PLJA", "Site2"]<- "Sevilleta"
xyz[xyz$Site.sp == "JOR.BOER", "Site2"]<- "Jornada"
xyz[xyz$Site.sp == "JOR.GUSA", "Site2"]<- "Jornada"

xyz[xyz$Site.sp == "MOAB.ACHY", "Sp2"]<- "ACHY"
xyz[xyz$Site.sp == "MOAB.BOGR", "Sp2"]<- "BOGR"
xyz[xyz$Site.sp == "MOAB.PLJA", "Sp2"]<- "PLJA"
xyz[xyz$Site.sp == "MOAB.GUSA", "Sp2"]<- "GUSA"
xyz[xyz$Site.sp == "SEV.GUSA", "Sp2"]<- "GUSA"
xyz[xyz$Site.sp == "SEV.ACHY", "Sp2"]<- "ACHY"
xyz[xyz$Site.sp == "SEV.BOER", "Sp2"]<- "BOER"
xyz[xyz$Site.sp == "SEV.BOGR", "Sp2"]<- "BOGR"
xyz[xyz$Site.sp == "SEV.PLJA", "Sp2"]<- "PLJA"
xyz[xyz$Site.sp == "JOR.BOER", "Sp2"]<- "BOER"
xyz[xyz$Site.sp == "JOR.GUSA", "Sp2"]<- "GUSA"

head(xyz)
xyz$Site2= ordered( as.character(xyz$Site2), levels = c("Moab",  "Sevilleta",  "Jornada"))
xyz$Sp2= ordered( as.character(xyz$Sp2), levels = c("BOGR",  "PLJA",  "BOER", "ACHY", "GUSA"))

Chl_bar <- ggplot(xyz, aes (x = Distance, y = emmean,width = .5, fill = Sp2))+facet_grid(Sp2~Site2)+
  theme_bw() +  geom_bar(stat = "identity", position = position_dodge())+
  theme(  axis.title.x = element_text(vjust=-0.35),
          axis.title.y = element_text(vjust=0.35) ,
          
          axis.title = element_text(size = 12),
          legend.position = "none",
          legend.title=element_text(size=11),
          axis.text = element_text(size = 12),
          legend.text=element_text(size=10),
          axis.ticks.x=element_blank(),
          legend.key = element_rect(fill = "white"),
          legend.background = element_rect(fill = "white"),
          panel.grid.major = element_line(colour = "white"),
          panel.grid.minor = element_blank())+
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.2,
                position=position_dodge(.9))+
  labs(y = expression(paste("Chlorophyll ", italic("a"), " (", mu, "g g"^-1, "soil)")) , x = expression("Distance from plant (cm)") )+
  scale_fill_manual(values=c("blue","grey80", "grey30","palegreen2","seagreen"   ))+ylim(c(-.3,4.5))

Chl_bar

cairo_pdf("chl_microsite.pdf",width = 4, height = 5 )
Chl_bar
dev.off()



detach("package:reshape2", unload=TRUE)
detach("package:data.table", unload=TRUE)
detach("package:dplyr", unload = TRUE)
library(data.table)
# library(reshape2)
library(dplyr)


 head(Alldata_sm)
melted <- Alldata_sm %>% group_by(Site.sp, Site, Species, Distance) %>%drop_na(d13C)%>%
   summarise(mean13C = mean(d13C), se13C=std(d13C), meanch = mean(Chlor), sechl= std(Chlor))


head(melted)
melted$Site.sp <- paste(melted$Site, melted$Species)

melted$color <- melted$Species
melted[which(melted$color == "ACHY"),"color"]<- "light green"
melted[which(melted$color == "BOER"),"color"]<- "black"
melted[which(melted$color == "BOGR"),"color"]<- "blue"
melted[which(melted$color == "GUSA"),"color"]<- "dark green"
melted[which(melted$color == "PLJA"),"color"]<- "grey"
sort(Alldata_sm$Site)
chlor_reg <- ggplot( Alldata_sm, aes(Chlor, d13C)) +
  theme_bw() +   theme(  axis.title.x = element_text(vjust=-0.35),
                         axis.title.y = element_text(vjust=0.35) ,
                         strip.text.x = element_blank(),
                         axis.title = element_text(size = 12),
                        legend.position ="none",
                         legend.title=element_text(size=11),
                         axis.text = element_text(size = 12),
                         legend.text=element_text(size = 10),
                         axis.ticks.x=element_blank(),
                         legend.key = element_rect(fill = "white"),
                         legend.background = element_rect(fill = "white"),
                         panel.grid.major = element_line(colour = "white"),
                         panel.grid.minor = element_blank()  ) +
    # geom_line(stat = "identity", aes(group = Site.sp), size = .5,  lty = 3) +
    scale_shape_manual(values=c(   15,17,19), name="Site", breaks=c( "MOAB", "SEV", "JOR"),    labels=c("MOAB" , "SEV","JOR"))+
scale_color_manual(values=c(  "palegreen2","black", "blue", "seagreen", "grey"), name="Species", breaks=c("BOER", "PLJA", "BOGR", "GUSA", "ACHY"),    labels=c("BOER", "PLJA", "BOGR", "GUSA", "ACHY"))+
 
ylim (-29.5, -15)+
    geom_point( aes(Chlor, d13C, shape = Site, col = Species, size = Distance , alpha = .5)) +
  # geom_errorbar( data=melted,  aes(y = mean13C, group = Site,  ymin = mean13C-se13C, ymax = mean13C+se13C, width = 0), size = .5)+
  # geom_errorbarh( data=melted,  aes(x = meanch, group = Site,  xmin = meanch-sechl, xmax = meanch+sechl, height = 0), size = .5)+
  # guides(fill = guide_legend(override.aes = list(linetype = 0)))+
  labs(x = expression(paste("Chlorophyll ", italic("a"), " (", mu, "g g"^-1, "soil)")) , y = expression("Observed cyanobacteria "*delta^{13}* "C (\u2030)") )+
  geom_smooth(method = "lm",  aes( color = Species), se=FALSE)+
# scale_x_continuous(limits = c(0, 16.4))+
  annotate("text", x = .9, y = -26, label = "d.")

chlor_reg





# check sample size 

xtabs(~Site.sp+Distance, data =Chlor_ag2)

cairo_pdf(paste("Fig_chl_", Sys.Date(), ".pdf", sep = ""), height = 3, width =4.3)
chlor_reg

dev.off()




# Root distribution -------------------------------------------------------
### stopped here 100719
Rootd <- read.table(file.choose(), header = T)
# JORSEVMOAB_root_dist_soil_weights_2017_EDR_180605.txt

head(Rootd)
summary(Rootd)
Rootd$Depth2 <- as.factor(Rootd$Depth)
Rootd$Depth2= ordered( as.character(Rootd$Depth2), levels = c("15",  "10",  "5"))
Rootd$Distance2<- as.factor(Rootd$Distance)

Rootd$sp_site <- (paste(Rootd$Species, Rootd$Site, sep = "."))
Rootd$sp_site= as.factor(ordered( as.character(Rootd$sp_site), levels = c("GUSA.M",  "GUSA.S",  "ACHY.M",
                                                               "ACHY.S", "PLJA.M", "BOGR.S", "BOER.S")))


Root05<- Rootd[which(Rootd$Depth2 == "5" & Rootd$Distance != "10" &  Rootd$ID!="F" & Rootd$ID != "H" & Rootd$ID != "A" & Rootd$ID != "I" & Rootd$ID !="D" & Rootd$ID!="E"),]
summary(Root05)
mod1 <- lmer(sqrt(Root_mass)~(Distance2*Species)  + (1|Site)+(1|ID), data = Root05)

Anova(mod1, type = 3)
#distance * sp interaction
pairs(emmeans (mod1, ~Distance2|Species), adjust = "fdr")
(-.163^2+.428^2)/.163^2
(-.177^2+.296^2)/.177^2
(-.145^2+.238^2)/.145^2

head(Root05)
Root05$Site <- as.character(Root05$Site2)
Root05[Root05$Site == "JOR1", "Site"]<- "JOR"
Root05[Root05$Site == "JOR2", "Site"]<- "JOR"
Root05[Root05$Site == "MOAB1", "Site"]<- "MOAB"
Root05[Root05$Site == "MOAB2", "Site"]<- "MOAB"
Root05[Root05$Site == "SEV1", "Site"]<- "SEV"
Root05[Root05$Site == "SEV2", "Site"]<- "SEV"

Root05$Site <- as.factor(Root05$Site)
head(Root05)
Root05$root_per_cm3<-Root05$Root_mass/14.25
plot(Root05$Root_mass, Root05$root_per_cm3)
root_ag <- Root05 %>% group_by(Site, Species, Distance) %>%
  summarise(meanmass = mean(root_per_cm3), semass=std(root_per_cm3))

root_ag$meanmass
head(root_ag)



head(Cyano)
Cyano_merge <- Cyano[, c(1,2,3, 6,7) ]
Cyano_merge$Site.sp <- paste(Cyano_merge$Site, Cyano_merge$Species, sep = ".")
head(Cyano_merge)
Cyano_merge_ag <- Cyano_merge %>% group_by(Site, Species, distance2) %>%drop_na(d13C)%>%
  summarise(mean13C = mean(d13C), se13C=std(d13C))
summary(Cyano_merge_ag)
summary(root_ag)
root_ag$Distance <- as.factor(root_ag$Distance)
root_13C <- full_join(Cyano_merge_ag,root_ag, by=c('Site'='Site',  'distance2'='Distance', 'Species'='Species'))
head(root_13C)
root_13C <- root_13C[which(root_13C$distance2 != "10"),]
root_13C$sisp <- paste(root_13C$Site, root_13C$Species, sep = ".")
### update with jor gusa?
# root_13Cno_sevplja <- root_13C[which(root_13C$sisp!= "SEV.PLJA"),]

mod1 <- lmer(mean13C~meanmass*Species+(1|Site), data = root_13C)
Anova(mod1, type = 3)
emmeans(mod1, ~meanmass)

emtrends(mod1, ~Species, var = "meanmass")
summary(mod1)

root_13C[root_13C$sisp=="MOAB.PLJA",]

root_reg <- ggplot( root_13C, aes(x=meanmass, y=mean13C)) +
  theme_bw() +   theme(  axis.title.x = element_text(vjust=-0.35),
                         axis.title.y = element_text(vjust=0.35) ,
                         strip.text.x = element_blank(),
                         axis.title = element_text(size = 12),
                         legend.position = "none",
                         legend.title=element_text(size=11),
                         axis.text = element_text(size = 12),
                         legend.text=element_text(size=10),
                         axis.ticks.x=element_blank(),
                         legend.key = element_rect(fill = "white"),
                         legend.background = element_rect(fill = "white"),
                         panel.grid.major = element_line(colour = "white"),
                         panel.grid.minor = element_blank()  ) +
  geom_smooth(method = "lm", aes(group=Species, color = Species), se=FALSE)+
 # geom_line(stat = "identity", aes(group = sisp), size = .5,  lty = 3) +
 scale_shape_manual(values=c( 15, 17, 19), name="Site", breaks=c( "MOAB", "SEV", "JOR"),    labels=c("MOAB" , "SEV","JOR"))+
  scale_color_manual(values=c(  "palegreen2","black", "blue", "seagreen", "grey"), name="Species", breaks=c("BOER", "PLJA", "BOGR", "GUSA", "ACHY"),    labels=c("BOER", "PLJA", "BOGR", "GUSA", "ACHY"))+
   ylim (-29.5, -15) + 
 # xlim  (-0.001 , 0.034)+
  geom_point( aes(meanmass, mean13C, shape = Site, col = Species, size = distance2, alpha = .5 ), data =root_13C) +
 # geom_errorbar(data=root_13C,  aes(y = mean13C,   ymin = mean13C-se13C, ymax = mean13C+se13C, width = 0), size = .5)+
#  geom_errorbarh(data=root_13C,  aes(x = meanmass,   xmin = meanmass-semass, xmax = meanmass+semass, height = 0), size = .5)+
 
 # guides(fill = guide_legend(override.aes = list(linetype = 0)))+
  labs(x = expression(paste("Roots (g cm"^-3, "soil)")) , y = expression("Observed cyanobacteria "*delta^{13}* "C (\u2030)") )+
  annotate("text", x = .001, y = -26, label = "f.")
# scale_x_discrete(expand=c(0.3,0), drop=FALSE)
root_reg




# cairo_pdf(paste("Root_reg_", Sys.Date(), ".pdf", sep = ""), height = 3, width =4.4)


# root_reg
# dev.off()


# 
 head(Root05)
 Root05$sp_site<- as.factor(paste(Root05$Site, Root05$Species, sep = "."))
 Root05$Distance<- as.factor(Root05$Distance)
 mod1 <- lmer(sqrt(root_per_cm3)~Distance*Species+(1|ID), data = Root05)
 Anova(mod1, type = 3)
# hist(resid(mod1))
 xyz <- summary(emmeans(mod1, ~Distance|Species))
# xyz$Site2 <- as.character(xyz$sp_site)

 (0.0407^2-.0612^2)/(.0612^2)
 
xyz$Sp2 <- as.character(xyz$sp_site)
xyz$mean <- xyz$emmean^2
xyz$up <- xyz$upper.CL^2
xyz$low <- xyz$lower.CL^2
xyz[xyz$sp_site == "MOAB.ACHY", "Site2"]<- "Moab"
xyz[xyz$sp_site == "MOAB.BOGR", "Site2"]<- "Moab"
xyz[xyz$sp_site == "MOAB.PLJA", "Site2"]<- "Moab"
xyz[xyz$sp_site == "MOAB.GUSA", "Site2"]<- "Moab"
xyz[xyz$sp_site == "SEV.GUSA", "Site2"]<- "Sevilleta"
xyz[xyz$sp_site == "SEV.ACHY", "Site2"]<- "Sevilleta"
xyz[xyz$sp_site == "SEV.BOER", "Site2"]<- "Sevilleta"
xyz[xyz$sp_site == "SEV.BOGR", "Site2"]<- "Sevilleta"
xyz[xyz$sp_site == "SEV.PLJA", "Site2"]<- "Sevilleta"
xyz[xyz$sp_site == "JOR.BOER", "Site2"]<- "Jornada"
xyz[xyz$sp_site == "JOR.GUSA", "Site2"]<- "Jornada"

xyz[xyz$sp_site == "MOAB.ACHY", "Sp2"]<- "ACHY"
xyz[xyz$sp_site == "MOAB.BOGR", "Sp2"]<- "BOGR"
xyz[xyz$sp_site == "MOAB.PLJA", "Sp2"]<- "PLJA"
xyz[xyz$sp_site == "MOAB.GUSA", "Sp2"]<- "GUSA"
xyz[xyz$sp_site == "SEV.GUSA", "Sp2"]<- "GUSA"
xyz[xyz$sp_site == "SEV.ACHY", "Sp2"]<- "ACHY"
xyz[xyz$sp_site == "SEV.BOER", "Sp2"]<- "BOER"
xyz[xyz$sp_site == "SEV.BOGR", "Sp2"]<- "BOGR"
xyz[xyz$sp_site == "SEV.PLJA", "Sp2"]<- "PLJA"
xyz[xyz$sp_site == "JOR.BOER", "Sp2"]<- "BOER"
xyz[xyz$sp_site == "JOR.GUSA", "Sp2"]<- "GUSA"

head(xyz)
xyz$Site2= ordered( as.character(xyz$Site2), levels = c("Moab",  "Sevilleta",  "Jornada"))
xyz$Sp2= ordered( as.character(xyz$Sp2), levels = c("BOGR",  "PLJA",  "BOER", "ACHY", "GUSA"))

Root_bar <- ggplot(xyz, aes (x = Distance, y = emmean,width = .5, fill = Sp2))+facet_grid(Sp2~Site2)+
  theme_bw() +  geom_bar(stat = "identity", position = position_dodge())+
  theme(  axis.title.x = element_text(vjust=-0.35),
          axis.title.y = element_text(vjust=0.35) ,
          
          axis.title = element_text(size = 12),
          legend.position = "none",
          legend.title=element_text(size=11),
          axis.text = element_text(size = 12),
          legend.text=element_text(size=10),
          axis.ticks.x=element_blank(),
          legend.key = element_rect(fill = "white"),
          legend.background = element_rect(fill = "white"),
          panel.grid.major = element_line(colour = "white"),
          panel.grid.minor = element_blank())+
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=.2,
                position=position_dodge(.9))+
  labs(y = expression(paste("Roots (g cm"^-3, "soil)"))  , x = expression("Distance from plant (cm)") )+
  scale_fill_manual(values=c("blue","light grey", "grey30","light green","dark green"   ))+ylim(ylim(c(-.02,.2))

Root_bar


cairo_pdf("root_microsite.pdf",width = 4, height = 5 )
Root_bar
dev.off()



 Root_all<- Rootd[ Rootd$Species != "PLJA" & Rootd$ID!="F" & Rootd$ID != "H" & Rootd$ID != "A" & Rootd$ID != "I" & Rootd$ID !="D" & Rootd$ID!="E",]
 xtabs(~Depth+Distance+Species,Root_all)
# mod1 <- lmer(sqrt(Root_mass)~(Distance2*Species*Depth)  + (1|Site2)+(1|ID), data = Root_all)
# Anova(mod1, type = 3)

# ggplot(Root_all, aes(x=factor(Distance), y= sqrt(Root_mass)))+geom_boxplot()+facet_grid(Depth~Species+Site)

Root05$Distance<-as.factor(Root05$Distance)
BOER <- Root05[Root05$Species=="BOER",]
BOER$Site

mod1 <- lm(sqrt(Root_mass)~(Distance*Site)  , data = BOER)
Anova(mod1, type = 3)
pairs(emmeans(mod1, ~Distance|Site),adjust="fdr")

PLJA <- Rootd[Rootd$Species=="PLJA",] # kaitlin to update!
mod1 <- lm(sqrt(Root_mass)~(Distance)  , data = PLJA)
Anova(mod1, type = 3)
pairs(emmeans(mod1, ~Distance|Site),adjust="fdr")


BOGR <- Rootd[Rootd$Species=="BOGR",] # Kristina to update!
mod1 <- lm(sqrt(Root_mass)~(Distance*Site)  , data = BOGR)
Anova(mod1, type = 3)
pairs(emmeans(mod1, ~Distance|Site),adjust="fdr")

ACHY <- Rootd[Rootd$Species=="ACHY",]
mod1 <- lm(sqrt(Root_mass)~(Distance*Site)  , data = ACHY)
Anova(mod1, type = 3)
pairs(emmeans(mod1, ~Distance|Site),adjust="fdr")


GUSA <- Rootd[Rootd$Species=="GUSA",]
mod1 <- lm(sqrt(Root_mass)~(Distance*Site)  , data = GUSA)
Anova(mod1, type = 3)




# Sequencing analysis -----------------------------------------------------

QPCR <- read.table(file.choose(), header = T)
# QPCR_results_nat_abund_190708.txt
QPCR$Distance <- as.factor(QPCR$Distance)
head(QPCR)
hist(log(QPCR$Qty_Mean))

plot(QPCR$Qty_Mean~factor(QPCR$Site3))
head(QPCR)
QPCR2<- QPCR[which(QPCR$Sample_ID!= "SPLJA0" & QPCR$Sample_ID!="SPLJA25"),]

#mod1<- lm(log(Qty_Mean)~Species*Distance, data = QPCR2)
#Anova(mod1, type =3)

mod1<- lm(log(Qty_Mean)~Species*Distance, data = QPCR)
Anova(mod1, type =3)
hist(resid(mod1))
AICc(mod1)
#88.81





# Composition_analysis ----------------------------------------------------

#comp<-read.table(file.choose(), header = T)
colnames(C_ag_t_qpcr)
head(C_ag_t_qpcr)
comp2 <- C_ag_t_qpcr[,c(1,2,9:12)]

comp <- spread(comp2, OTU, abund)
colnames(comp)
comp$Distance<- as.factor(comp$Distance)

head(Cyano_merge)
Cyano_merge$Site3 <- as.character(Cyano_merge$Site2)

Cyano_merge[Cyano_merge$Site2=="JOR", "Site3"]<- "J"
Cyano_merge[Cyano_merge$Site2=="MOAB", "Site3"]<- "C"
Cyano_merge[Cyano_merge$Site2=="SEV", "Site3"]<- "S"
Cyano_merge$Site3 <- as.factor(Cyano_merge$Site3)
colnames(Cyano_merge)
colnames(comp)
library(vegan)
comp$div <- ((diversity(comp[,c(5:22)])))

comp$J <- comp$div/log(specnumber(comp[,c(5:22)]))

comp$Srar<- rarefy(round(comp[,c(5:22)]), min(rowSums(comp[,c(5:22)])))

colnames(Cyano_merge)

Cyano_merge$Site3

comp$Site2 <- as.character(comp$Site)
comp[comp$Site=="C_BOGR_0", "Site2"] <- "C"
comp[comp$Site=="C_BOGR_25", "Site2"] <- "C"
comp[comp$Site=="S_BOGR_0", "Site2"] <- "S"
comp[comp$Site=="S_BOGR_25", "Site2"] <- "S"
comp[comp$Site=="C_PLJA_0", "Site2"] <- "C"
comp[comp$Site=="C_PLJA_25", "Site2"] <- "C"
comp[comp$Site=="S_PLJA_0", "Site2"] <- "S"
comp[comp$Site=="S_PLJA_25", "Site2"] <- "S"
comp[comp$Site=="C_GUSA_0", "Site2"] <- "C"
comp[comp$Site=="C_GUSA_25", "Site2"] <- "C"
comp[comp$Site=="S_GUSA_0", "Site2"] <- "S"
comp[comp$Site=="S_GUSA_25", "Site2"] <- "S"
comp[comp$Site=="C_ACHY_0", "Site2"] <- "C"
comp[comp$Site=="C_ACHY_25", "Site2"] <- "C"
comp[comp$Site=="S_ACHY_0", "Site2"] <- "S"
comp[comp$Site=="S_ACHY_25", "Site2"] <- "S"
comp[comp$Site=="S_BOER_0", "Site2"] <- "S"
comp[comp$Site=="S_BOER_25", "Site2"] <- "S"
comp[comp$Site=="J_GUSA_0", "Site2"] <- "J"
comp[comp$Site=="J_GUSA_25", "Site2"] <- "J"
comp[comp$Site=="J_BOER_0", "Site2"] <- "J"
comp[comp$Site=="J_BOER_25", "Site2"] <- "J"

head(Cyano_merge)
head(comp)

cyano_merg_ag <- aggregate(Cyano_merge$d13C, by=list(Site=Cyano_merge$Site3, Species = Cyano_merge$Species, Distance = Cyano_merge$distance2), FUN=mean)

head(comp)
head(cyano_merg_ag)
colnames(comp)
colnames(cyano_merg_ag)

comp_13 <- full_join(cyano_merg_ag, comp, by = c('Site' = 'Site2' , 'Distance' = 'Distance', 'Species' = 'Species' ))
colnames(comp_13)
head(comp_13)
comp_13$total <- rowSums(comp_13[,c(7:25)])
colnames(comp_13)
comp_13$site.sp <- paste(comp_13$Site, comp_13$Species)



library(tibble)
library(dplyr)
library(tibble)
colnames(comp)

# tots <- rowSums(comp_13[,c(5:22)])

hist(comp$div)
mod1<- lm((div)~Distance*Species, data = comp)
Anova(mod1, type = 3)
hist(resid(mod1))
emmeans(mod1,~ Species)
pairs(emmeans(mod1,~ Species), adjust = "fdr")

# hist(PCRcomp$J)
# mod1<- lm((J)~Species*Distance, data = PCRcomp)
# Anova(mod1, type = 3)
# hist(resid(mod1))

mod1<- lm((Srar)~Species*Distance, data = comp)
Anova(mod1, type = 3)
hist(resid(mod1))
emmeans(mod1,~ Species)
pairs(emmeans(mod1,~ Species), adjust = "fdr")







### update with PLJA
plot(comp_13$x ~ comp_13$div)
mod10 <- lm(x ~ Species *div, data = comp_13)
Anova(mod10, type = 3)

head(comp_13)
div_obs<-ggplot( comp_13, aes(x=div, y=x)) +
  theme_bw() +   theme(  axis.title.x = element_text(vjust=-0.35),
                         axis.title.y = element_text(vjust=0.35) ,
                         strip.text.x = element_blank(),
                         axis.title = element_text(size = 12),
                          legend.position = "none",
                         legend.title=element_text(size=11),
                         axis.text = element_text(size = 12),
                         legend.text=element_text(size=10),
                         axis.ticks.x=element_blank(),
                         legend.key = element_rect(fill = "white"),
                         legend.background = element_rect(fill = "white"),
                         panel.grid.major = element_line(colour = "white"),
                         panel.grid.minor = element_blank()  ) +
  # geom_errorbar(position = dodge, data=xyzbogr,  aes(y = emmean, group = Site,  ymin = lower.CL, ymax = upper.CL, width = 0.2), size = .8)+
  # geom_line(stat = "identity", aes(group = site.sp), size = .5,  lty = 3) +
    scale_shape_manual(values=c( 17, 15, 19), name="Site", breaks=c( "C", "S", "J"),    labels=c("M" , "S","J"))+
    scale_color_manual(values=c( "blue",  "grey","black",  "palegreen2","seagreen"), name="species", breaks=c("BOER", "PLJA", "BOGR", "GUSA", "ACHY"),    labels=c("BOER", "PLJA", "BOGR", "GUSA", "ACHY"))+
  annotate("text", x = .5, y = -26, label = "a.")+
  geom_smooth(method = "lm", aes(color = Species), se=FALSE)+
   ylim (-29.5, -15)+
  geom_point( aes(div, x, shape = Site, col = Species, size = Distance, alpha = .9 ), data = comp_13) +
 # guides(fill = guide_legend(override.aes = list(linetype = 0)))+
  labs(x = expression("Cyanobacteria diversity"),  y = expression("Observed cyanobacteria "*delta^{13}* "C (\u2030)") )
# scale_x_discrete(expand=c(0.3,0), drop=FALSE)
div_obs

#cairo_pdf(paste("Fig_div_pred_", Sys.Date(), ".pdf", sep = ""), height = 3, width =4.4)
#div_obs

#dev.off()



### update with PLJA
plot(comp_13$x ~ comp_13$div)
mod10 <- lm(x ~ Species *div, data = comp_13)
Anova(mod10, type = 3)
comp_13$Species
str(comp_13)
### update with PLJA
plot(comp_13$x ~ comp_13$Srar)
comp_13$Species
mod10 <- lm(x ~ Species *Srar, data = comp_13)
Anova(mod10, type = 3)
summary(mod10)
emtrends(mod10, ~Species, var = "Srar")


Srar_obs<-ggplot( comp_13, aes(x=Srar, y=x)) +
  theme_bw() +   theme(  axis.title.x = element_text(vjust=-0.35),
                         axis.title.y = element_text(vjust=0.35) ,
                         strip.text.x = element_blank(),
                         axis.title = element_text(size = 12),
                        legend.position = "none",
                         legend.title=element_text(size=11),
                         axis.text = element_text(size = 12),
                         legend.text=element_text(size=10),
                         axis.ticks.x=element_blank(),
                         legend.key = element_rect(fill = "white"),
                         legend.background = element_rect(fill = "white"),
                         panel.grid.major = element_line(colour = "white"),
                         panel.grid.minor = element_blank()  ) +
  # geom_errorbar(position = dodge, data=xyzbogr,  aes(y = emmean, group = Site,  ymin = lower.CL, ymax = upper.CL, width = 0.2), size = .8)+
 # geom_line(stat = "identity", aes(group = site.sp), size = .5,  lty = 3) +
  scale_shape_manual(values=c( 17, 15, 19), name="site", breaks=c( "C", "S", "J"),    labels=c("M" , "S","J"))+
  scale_color_manual(values=c( "blue",  "grey","black",  "palegreen2","seagreen"), name="species", breaks=c("BOER", "PLJA", "BOGR", "GUSA", "ACHY"),    labels=c("BOER", "PLJA", "BOGR", "GUSA", "ACHY"))+
  annotate("text", x = .9, y = -26, label = "b.")+
geom_smooth(method = "lm", aes(color = Species), se = FALSE)+
 ylim (-29.5, -15)+
  geom_point( aes(Srar, x, shape = Site, col = Species, size = Distance, alpha = .5 ), data = comp_13) +
 # guides(fill = guide_legend(override.aes = list(linetype = 0)))+
  labs(x = expression("Cyanobacteria richness"),  y = expression("Observed cyanobacteria "*delta^{13}* "C (\u2030)") )
# scale_x_discrete(expand=c(0.3,0), drop=FALSE)
Srar_obs


(comp[,1])
summary(comp)

samp2 <- comp[,c(5:22)]
rownames(samp2) <- comp[,1]
head(samp2)

library(labdsv)
comp.red <- vegtab(samp2, min=0.05*nrow(samp2))
head(comp.red)

comp.st <- wisconsin(comp.red)
head(comp.st)

comp.dist <- vegdist(comp.st, method = "bray")
comp$group <- paste(comp$Species, comp$Distance)
mod1 <- adonis(comp.dist~Species*Distance, data = comp, perm = 9999)
mod1
# can't do dispers on just 2 points
mod2 <- betadisper(comp.dist, comp$group)
anova(mod2)
(mod.HSD <- TukeyHSD(mod2))
plot(mod.HSD)
plot(mod2)
plot(mod2, ellipse = TRUE, hull = FALSE)
boxplot(mod2)
mds3 <- metaMDS(comp.dist,  k=3, perm=999)
stressplot(mds3)
mds2 <- metaMDS(comp.dist,  k=2, perm=999)

cairo_pdf(paste("NMDS", Sys.Date(), ".pdf"), width = 10, height = 8)
fig <- ordiplot(mds2 )
text(fig,"sites")
ordihull(mds2, comp$Species )
dev.off()

# nat_abund_cyano_community_180926.txt
# comp4$Sample_ID <- rownames(comp4)
head(QPCR)
head(comp4)
library (vegan)
library(dplyr)


# Predicted vs. observed d13C ---------------------------------------------
head(Cyano)
head(obpr)

Cyano$Site3<- paste(substr(as.character(Cyano$Site),1,1), Cyano$Species, Cyano$distance, sep = "_" )

unique(Cyano$Site3)
Cyano$Site4 <- Cyano$Site3
Cyano[Cyano$Site3 == "M_BOGR_0", "Site4"] <- "C_BOGR_0"
Cyano[Cyano$Site3 == "M_BOGR_25", "Site4"] <- "C_BOGR_25"
Cyano[Cyano$Site3 == "M_GUSA_0", "Site4"] <- "C_GUSA_0"
Cyano[Cyano$Site3 == "M_PLJA_0", "Site4"] <- "C_PLJA_0"
Cyano[Cyano$Site3 == "M_ACHY_0", "Site4"] <- "C_ACHY_0"
Cyano[Cyano$Site3 == "M_GUSA_25", "Site4"] <- "C_GUSA_25"
Cyano[Cyano$Site3 == "M_PLJA_25", "Site4"] <- "C_PLJA_25"
Cyano[Cyano$Site3 == "M_ACHY_25", "Site4"] <- "C_ACHY_25"
head(Cyano)
head(pred)

Cyano_ag <-  aggregate(Cyano$d13C, by=list(Site4= Cyano$Site4,Site=Cyano$Site, Species = Cyano$Species, distance = Cyano$distance), FUN=mean)


PredObs <- left_join(Cyano_ag, pred, by = c('Site4'='Site'))
head(PredObs)


# mod1 <- lmer(predicted~distance*observed+(1|site), data = PredObs)
mod1 <- lmer(x.x~x.y*Species+(1|Site), data = PredObs)
Anova(mod1, type = 3)


PredObs$site.sp <- as.factor(paste(PredObs$Site, PredObs$Species, sep = "."))
unique(PredObs$Species)
pred_obs<-ggplot( PredObs, aes(x=x.y, y=x.x)) +
  theme_bw() +   theme(  axis.title.x = element_text(vjust=-0.35),
                         axis.title.y = element_text(vjust=0.35) ,
                         strip.text.x = element_blank(),
                         axis.title = element_text(size = 12),
                         legend.position = "none",
                         legend.title=element_text(size=11),
                         axis.text = element_text(size = 12),
                         legend.text=element_text(size=10),
                         axis.ticks.x=element_blank(),
                         legend.key = element_rect(fill = "white"),
                         legend.background = element_rect(fill = "white"),
                         panel.grid.major = element_line(colour = "white"),
                         panel.grid.minor = element_blank()  ) +
  # geom_errorbar(position = dodge, data=xyzbogr,  aes(y = emmean, group = Site,  ymin = lower.CL, ymax = upper.CL, width = 0.2), size = .8)+
 #  geom_line(stat = "identity", aes(group = site.sp), size = .5,  lty = 3) +
  geom_smooth(method = "lm", aes(color = Species), se=FALSE)+
  scale_shape_manual(values=c(17,   19,15), name="Site", breaks=c( "MOAB", "SEV", "JOR"),    labels=c("M" , "S","J"))+
 scale_color_manual(values=c(  "blue","grey", "black", "palegreen2", "seagreen"), name="Species", breaks=c("BOER", "PLJA", "BOGR", "GUSA", "ACHY"),    labels=c("BOER", "PLJA", "BOGR", "GUSA", "ACHY"))+
  annotate("text", x = -17, y = -26, label = "c.")+
     ylim (-29.5, -15)+
  #+ xlim (-17.3, -15.3)+
  geom_point( aes(x.y, x.x, shape = Site, col = Species, size = distance, alpha = .5 ), data = PredObs) +
#  guides(fill = guide_legend(override.aes = list(linetype = 0)))+
  labs(x = expression("Predicted cyanobacteria "*delta^{13}* "C (\u2030)"),  y = expression("Observed cyanobacteria "*delta^{13}* "C (\u2030)") )
  # scale_x_discrete(expand=c(0.3,0), drop=FALSE)
pred_obs

 cairo_pdf(paste("Fig_obs_pred_", Sys.Date(), ".pdf", sep = ""), height = 3, width =4.4)
 pred_obs

 dev.off()
citation()

 ###### 

# Cyano_sequencing_2019 ---------------------------------------------------

Cyano_com <- read.csv(file.choose())
#otu-table-cyano_natabund_wtax_and_ontree_rem_doubletons_190709.csv
head(Cyano_com)
colnames(Cyano_com)

library(dplyr)
library(tidyr)

head(OTU2)

Cyano2 <- gather (OTU2, , "counts", 3:22)
head(Cyano2)
unique(Cyano2$Site)

Cyano2[Cyano2$Site == "M_BOGR_0", "Site"] <- "C_BOGR_0"
Cyano2[Cyano2$Site == "M_BOGR_25", "Site"] <- "C_BOGR_25"
Cyano2[Cyano2$Site == "M_GUSA_0", "Site"] <- "C_GUSA_0"
Cyano2[Cyano2$Site == "M_PLJA_0", "Site"] <- "C_PLJA_0"
Cyano2[Cyano2$Site == "M_ACHY_0", "Site"] <- "C_ACHY_0"
Cyano2[Cyano2$Site == "M_GUSA_25", "Site"] <- "C_GUSA_25"
Cyano2[Cyano2$Site == "M_PLJA_25", "Site"] <- "C_PLJA_25"
Cyano2[Cyano2$Site == "M_ACHY_25", "Site"] <- "C_ACHY_25"

unique(Cyano2$Site)

#recode(Cyano2$Site, M_BOGR_0 = "C_BOGR_0", M_BOGR_25 = "C_BOGR_25")

#Cyano2 %>%
#  mutate(Site=replace(Site, Site=="M_BOGR_0", "C_BOGR_0")) %>%
#  mutate(Site=replace(Site, Site=="M_BOGR_25", "C_BOGR_25")) %>%
#  mutate(Site=replace(Site, Site=="M_GUSA_0", "C_GUSA_0")) %>%
#  mutate(Site=replace(Site, Site=="M_GUSA_25", "C_GUSA_25")) %>%
#  mutate(Site=replace(Site, Site=="M_PLJA_0", "C_PLJA_0)")) %>%
#  mutate(Site=replace(Site, Site=="M_PLJA_25", "C_PLJA_25)")) %>%
#  mutate(Site=replace(Site, Site=="M_ACHY_0", "C_ACHY_0)")) %>%
#  mutate(Site=replace(Site, Site=="M_ACHY_25", "C_ACHY_25)")) 
  

Cyano2$Site= as.factor(ordered( as.character(Cyano2$Site),
                          levels = c("C_BOGR_0","C_BOGR_25",
                                     "S_BOGR_0","S_BOGR_25",
                                     "C_PLJA_0","C_PLJA_25",
                                     "S_PLJA_0","S_PLJA_25",
                                     "S_BOER_0","S_BOER_25",
                                     "J_BOER_0","J_BOER_25",
                                     "C_ACHY_0","C_ACHY_25",
                                     "S_ACHY_0","S_ACHY_25",
                                     "C_GUSA_0","C_GUSA_25",
                                     "S_GUSA_0","S_GUSA_25",
                                     "J_GUSA_0","J_GUSA_25"
                          )))

unique(Cyano2$Site)

# aggregate seq numbers by phyllum
cyano_ag_p <- aggregate(Cyano2$counts, by=list(Site=Cyano2$Site, Phylla = Cyano2$p), FUN=sum)
head(cyano_ag_p)

# total sequence numbers per sample
cyano_ag_tot <- aggregate(Cyano2$counts, by=list(Site=Cyano2$Site), FUN=sum)
head(cyano_ag_tot)


Cyano_ag_pt <- full_join(cyano_ag_p, cyano_ag_tot, by = c('Site' = 'Site' ))
head(Cyano_ag_pt)
unique(Cyano_ag_pt$Site)
plot(as.factor(Cyano_ag_pt$Site), Cyano_ag_pt$x.y)

Cyano_ag_pt$rel_abun <- Cyano_ag_pt$x.x/Cyano_ag_pt$x.y


head(QPCR)
fillers <- c("Unassigned",
             "Verrucomicrobia",

             "Thermi",
             "Tenericutes",
             
             "Proteobacteria",
             "Planctomycetes",
             "Parvarchaeota",
             "Nitrospirae",
             
             "Gemmatimonadetes",
             "Fusobacteria",
             "Firmicutes",
             "Fibrobacteres",
             "Euryarchaeota",
             "Elusimicrobia",
             "Cyanobacteria",
             "Crenarchaeota",
             "Chloroflexi",
             "Chlorobi",
             "Chlamydiae",
             "Bacteroidetes",
             "Armatimonadetes",
             "Actinobacteria",
             "Acidobacteria")

Cyano_ag_pt$PHYLUM <- factor(Cyano_ag_pt$Phylla, fillers)
sort(unique(Cyano_ag_pt$Site))



Cyano_ag_pt$Site= ordered( as.character(Cyano_ag_pt$Site),
                          levels = c("C_BOGR_0","C_BOGR_25",
                                     "S_BOGR_0","S_BOGR_25",
                                     "C_PLJA_0","C_PLJA_25",
                                     "S_PLJA_0","S_PLJA_25",
                                     "S_BOER_0","S_BOER_25",
                                     "J_BOER_0","J_BOER_25",
                                     "C_ACHY_0","C_ACHY_25",
                                     "S_ACHY_0","S_ACHY_25",
                                     "C_GUSA_0","C_GUSA_25",
                                     "S_GUSA_0","S_GUSA_25"
                                  
                          ))




colour <- c("azure3",
            "lightblue2",
            "lightcyan2",
            "gray70",
            "gray48",
            "gray61",
            "gray25",
            "lightsteelblue4",
            
            "steelblue",
            "darkcyan",
            "cyan3",
            "lightblue4",
            "lightsteelblue",
            "gray50",
            "olivedrab3",
            "lightsteelblue3",
            "darkslategray4",
            "lightcyan4",
            "cornflowerblue",
            "skyblue4",
            "aliceblue",
            "gray79",
            "lightblue")

head(Cyano_ag_pt)
QPCR$Site3 <- paste(QPCR$Site2, QPCR$Species, QPCR$Distance, sep = "_")

QPCR[QPCR$Site3 == "M_BOGR_0", "Site3"] <- "C_BOGR_0"
QPCR[QPCR$Site3 == "M_BOGR_25", "Site3"] <- "C_BOGR_25"
QPCR[QPCR$Site3 == "M_GUSA_0", "Site3"] <- "C_GUSA_0"
QPCR[QPCR$Site3 == "M_PLJA_0", "Site3"] <- "C_PLJA_0"
QPCR[QPCR$Site3 == "M_ACHY_0", "Site3"] <- "C_ACHY_0"
QPCR[QPCR$Site3 == "M_GUSA_25", "Site3"] <- "C_GUSA_25"
QPCR[QPCR$Site3 == "M_PLJA_25", "Site3"] <- "C_PLJA_25"
QPCR[QPCR$Site3 == "M_ACHY_25", "Site3"] <- "C_ACHY_25"


P_ag_t_qpcr <- full_join(Cyano_ag_pt, QPCR, by = c('Site' = 'Site3' ))
head(P_ag_t_qpcr)
P_ag_t_qpcr$abund <- P_ag_t_qpcr$rel_abun*P_ag_t_qpcr$Qty_Mean

P_ag_t_qpcr$Site.y<- as.character(P_ag_t_qpcr$Site.y)
P_ag_t_qpcr[which(P_ag_t_qpcr$Site.y=="MOAB"),"Site.y"]<- "COL"

P_ag_t_qpcr$xlab <- paste(P_ag_t_qpcr$Site.y, P_ag_t_qpcr$Distance, sep = " ")

 P_ag_t_qpcr$OTU <- factor(P_ag_t_qpcr$OTU, fillers)
P_ag_t_qpcr$xlab= as.factor(ordered( as.character(P_ag_t_qpcr$xlab),
                                     levels = c("COL 0","COL 25",
                                                "SEV 0","SEV 25",
                                                "JOR 0","JOR 25"
                                     )))

P_ag_t_qpcr$Species= as.factor(ordered( as.character(P_ag_t_qpcr$Species),
                                        levels = c("BOGR","PLJA",
                                                   "BOER","ACHY",
                                                   "GUSA"
                                        )))

P_ag_t_qpcr$lab <- P_ag_t_qpcr$Species
P_ag_t_qpcr$lab<- recode_factor(P_ag_t_qpcr$lab, 'BOGR' = "B. gracilis", 'BOER' = "B. eriopoda",
                                'PLJA' = "P. jamesii", 'ACHY' = "A. hymenoides", 'GUSA' = "G. sarothrae")
P_ag_t_qpcr$lab = as.factor(ordered( as.character(P_ag_t_qpcr$lab),
                                     levels = c("B. gracilis","P. jamesii",
                                                "B. eriopoda","A. hymenoides",
                                                "G. sarothrae")))
# 
# 
# P_ag_t_qpcr$Site= ordered( as.character(P_ag_t_qpcr$Site),
#                            levels = c("C_BOGR_0","C_BOGR_25",
#                                       "S_BOGR_0","S_BOGR_25",
#                                       "C_PLJA_0","C_PLJA_25",
#                                       "S_PLJA_0","S_PLJA_25",
#                                       "S_BOER_0","S_BOER_25",
#                                       "J_BOER_0","J_BOER_25",
#                                       "C_ACHY_0","C_ACHY_25",
#                                       "S_ACHY_0","S_ACHY_25",
#                                       "C_GUSA_0","C_GUSA_25",
#                                       "S_GUSA_0","S_GUSA_25"
#                                       
#                            ))


#Hacemos el gráfico de barras apiladas
phyl<-ggplot(P_ag_t_qpcr,aes(x = xlab, y = abund/1000000, fill = PHYLUM, order = -as.numeric(PHYLUM))) + 
  geom_bar(stat="identity") + facet_grid(.~lab)+
  ylab(label = expression("N 16SrRNA gene copies x 10"^{6})) +
  scale_y_continuous() +
  scale_fill_manual(values=colour, fillers) +
 # scale_x_discrete(limits=c("CONTROL 13", "CONTROL 21", "CONTROL 4", "CONTROL 6", "CONTROL 30", "CONTROL 26", "CONTROL 16", "CONTROL 10", "CONTROL 23", "CONTROL 1","", "DELAY 19", "DELAY 29", "DELAY 15", "DELAY 25", "DELAY 28", "DELAY 8", "DELAY 32", "DELAY 2", "DELAY 7", "DELAY 12","", "DROUGHT 20", "DROUGHT 5", "DROUGHT 27", "DROUGHT 31", "DROUGHT 9", "DROUGHT 24", "DROUGHT 18", "DROUGHT 11", "DROUGHT 3", "DROUGHT 14"), name="Treatments") +
  theme_bw() +
  theme(strip.text = element_text( face = "italic", size = 14))+
  theme(panel.border = element_rect(colour = "white", linetype="blank")) +
  theme(legend.key = element_blank()) +
  theme(panel.grid.minor = element_line(colour = NA)) +
  theme(panel.grid.major = element_line(colour = NA)) +
  theme(axis.title.x = element_text(colour="black", size=30, family="Times New Roman"),
        axis.text.x  = element_text(angle=0, vjust=0.5, size=30, family="Times New Roman")) +
  theme(axis.title.y = element_text(colour="black", size=30, family="Times New Roman"),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=30, family="Times New Roman")) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(colour="black", face="italic", size = 30, family="Times New Roman")) +
  theme(axis.line = element_line(colour = "black", size = 0.5, linetype = "solid")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x = "Site, Distance from plant")

cairo_pdf("Phyllum.pdf", height = 8, width = 20)
phyl

dev.off()







Only_cyano <- Cyano2[Cyano2$p == "Cyanobacteria" | Cyano2$p == "Chloroflexi",]
head(Only_cyano)


cyano_ag_c <- aggregate(Only_cyano$counts, by=list(Site=Only_cyano$Site, OTU = Only_cyano$OTU), FUN=sum)
head(cyano_ag_c)
cyano_ag_c[cyano_ag_c$OTU  == "Unassigned",]
# For predicted - use relative to CYANO only!
 cyano_ag_ctot <- aggregate(Only_cyano$counts, by=list(Site=Only_cyano$Site), FUN=sum)
head(cyano_ag_ctot)




Cyanoonly_ag_ct <- left_join(cyano_ag_c, cyano_ag_ctot, by = c('Site' = 'Site' ))
head(Cyanoonly_ag_ct)

Cyanoonly_ag_ct$rel_ab <- Cyanoonly_ag_ct$x.x/Cyanoonly_ag_ct$x.y





# for QPCR, use relative to aLL sequences
Cyanoonly_ag_t <- left_join(cyano_ag_c, cyano_ag_tot, by = c('Site' = 'Site' ))
head(Cyanoonly_ag_t)
plot(as.factor(Cyanoonly_ag_t$Site), Cyanoonly_ag_t$x.x)
Cyanoonly_ag_t[Cyanoonly_ag_t$OTU  == "Unassigned",]

Cyanoonly_ag_t$rel_abun <- Cyanoonly_ag_t$x.x/Cyanoonly_ag_t$x.y

sort(unique(Cyanoonly_ag_t$OTU))
head(Cyanoonly_ag_t)


Cyanoonly_ag_t$Site= as.factor(ordered( as.character(Cyanoonly_ag_t$Site),
                                levels = c("C_BOGR_0","C_BOGR_25",
                                           "S_BOGR_0","S_BOGR_25",
                                           "C_PLJA_0","C_PLJA_25",
                                           "S_PLJA_0","S_PLJA_25",
                                           "S_BOER_0","S_BOER_25",
                                           "J_BOER_0","J_BOER_25",
                                           "C_ACHY_0","C_ACHY_25",
                                           "S_ACHY_0","S_ACHY_25",
                                           "C_GUSA_0","C_GUSA_25",
                                           "S_GUSA_0","S_GUSA_25"
                                       )))




QPCR$Site3 <- paste(QPCR$Site2, QPCR$Species, QPCR$Distance, sep = "_")

QPCR[QPCR$Site3 == "M_BOGR_0", "Site3"] <- "C_BOGR_0"
QPCR[QPCR$Site3 == "M_BOGR_25", "Site3"] <- "C_BOGR_25"
QPCR[QPCR$Site3 == "M_GUSA_0", "Site3"] <- "C_GUSA_0"
QPCR[QPCR$Site3 == "M_PLJA_0", "Site3"] <- "C_PLJA_0"
QPCR[QPCR$Site3 == "M_ACHY_0", "Site3"] <- "C_ACHY_0"
QPCR[QPCR$Site3 == "M_GUSA_25", "Site3"] <- "C_GUSA_25"
QPCR[QPCR$Site3 == "M_PLJA_25", "Site3"] <- "C_PLJA_25"
QPCR[QPCR$Site3 == "M_ACHY_25", "Site3"] <- "C_ACHY_25"


C_ag_t_qpcr <- full_join(Cyanoonly_ag_t, QPCR, by = c('Site' = 'Site3' ))
head(C_ag_t_qpcr)

C_ag_t_qpcr$abund <- C_ag_t_qpcr$rel_abun*C_ag_t_qpcr$Qty_Mean
head(sort(C_ag_t_qpcr$abund, decreasing = T),20)

unique(C_ag_t_qpcr$OTU)

C_ag_t_qpcr[is.na(C_ag_t_qpcr$OTU),]

fillers <- c("Microcoleus_vaginatus",
             "Microcoleus_steenstrupii_complex", # new
             "Microcoleus_steenstrupii_Clade_I",
           #  "CLADE 2 - M. steenstrupii",
            # "CLADE 3 - M. steenstrupii",
             "Microcoleus_steenstrupii_Clade_IV",
            # "CLADE 5 - M. steenstrupii",
            # "CLADE 6 - M. steenstrupii",
             "Microcoleus_steenstrupii_Clade_VII",
             "Microcoleus_paludosus",
           "Microcoleus_chthonoplastes", #new
             "Leptolyngbya",
           "Lyngbya",#new
           #  "Phormidium",
           #  "Trichocoleus",
             "Chroococcidiopsis",
           "Crinalium", #new
             "Scytonema",
             "Nostoc",
           "Fischerella", #new
           "Calothrix", #new
           "Cephalothrix",#new
             "Tolypothrix",
           "Vampirovibrio", #new
             "Unassigned")



  colour <- c("skyblue3",
              "green", #new
              "palegreen2",
             # "olivedrab3",
            #  "darkgreen",
              "springgreen",
            #  "springgreen4",
            #  "yellowgreen",
              "seagreen2",
              "darkolivegreen2",
              "springgreen4",#new
              "mediumpurple",
             "purple",#new
            #  "darkgoldenrod4",
            #  "indianred2",
              "indianred4",
            "red",#new
              "darkgoldenrod2",
              "darkgoldenrod3",
            "khaki", #new
            "yellow3",#new
            "yellow",#new
         
            "lavender",#new
            "gray36",
              "gray53"
  )

head(C_ag_t_qpcr)
C_ag_t_qpcr$Site.y<- as.character(C_ag_t_qpcr$Site.y)
C_ag_t_qpcr[which(C_ag_t_qpcr$Site.y=="MOAB"),"Site.y"]<- "COL"

C_ag_t_qpcr$xlab <- paste(C_ag_t_qpcr$Site.y, C_ag_t_qpcr$Distance, sep = " ")

C_ag_t_qpcr$OTU <- factor(C_ag_t_qpcr$OTU, fillers)
C_ag_t_qpcr$xlab= as.factor(ordered( as.character(C_ag_t_qpcr$xlab),
                                     levels = c("COL 0","COL 25",
                                                "SEV 0","SEV 25",
                                                "JOR 0","JOR 25"
                                     )))

C_ag_t_qpcr$Species= as.factor(ordered( as.character(C_ag_t_qpcr$Species),
                                     levels = c("BOGR","PLJA",
                                                "BOER","ACHY",
                                                "GUSA"
                                     )))

C_ag_t_qpcr$lab <- C_ag_t_qpcr$Species
C_ag_t_qpcr$lab<- recode_factor(C_ag_t_qpcr$lab, 'BOGR' = "B. gracilis", 'BOER' = "B. eriopoda",
              'PLJA' = "P. jamesii", 'ACHY' = "A. hymenoides", 'GUSA' = "G. sarothrae")
C_ag_t_qpcr$lab = as.factor(ordered( as.character(C_ag_t_qpcr$lab),
                                        levels = c("B. gracilis","P. jamesii",
                                                   "B. eriopoda","A. hymenoides",
                                                   "G. sarothrae")))
# C_ag_t_qpcr$Site= as.factor(ordered( as.character(C_ag_t_qpcr$Site),
#                                         levels = c("C_BOGR_0","C_BOGR_25",
#                                                    "S_BOGR_0","S_BOGR_25",
#                                                    "C_PLJA_0","C_PLJA_25",
#                                                    "S_PLJA_0","S_PLJA_25",
#                                                    "S_BOER_0","S_BOER_25",
#                                                    "J_BOER_0","J_BOER_25",
#                                                    "C_ACHY_0","C_ACHY_25",
#                                                    "S_ACHY_0","S_ACHY_25",
#                                                    "C_GUSA_0","C_GUSA_25",
#                                                    "S_GUSA_0","S_GUSA_25"
#                                         )))
# 

head(C_ag_t_qpcr)

plot_otu<-ggplot(C_ag_t_qpcr,aes(x = xlab, y = abund/1000000, fill = OTU, order = -as.numeric(OTU))) + 
  geom_bar(stat="identity") + facet_grid(.~lab)+
  ylab(label = expression("N 16SrRNA gene copies x 10"^{6})) +
  scale_fill_manual(values=colour, fillers) +
  # scale_x_discrete(limits=c("CONTROL 13", "CONTROL 21", "CONTROL 4", "CONTROL 6", "CONTROL 30", "CONTROL 26", "CONTROL 16", "CONTROL 10", "CONTROL 23", "CONTROL 1","", "DELAY 19", "DELAY 29", "DELAY 15", "DELAY 25", "DELAY 28", "DELAY 8", "DELAY 32", "DELAY 2", "DELAY 7", "DELAY 12","", "DROUGHT 20", "DROUGHT 5", "DROUGHT 27", "DROUGHT 31", "DROUGHT 9", "DROUGHT 24", "DROUGHT 18", "DROUGHT 11", "DROUGHT 3", "DROUGHT 14"), name="Treatments") +
  theme_bw() +
  theme(strip.text = element_text( face = "italic", size = 14))+
  theme(panel.border = element_rect(colour = "white", linetype="blank")) +
  theme(legend.key = element_blank()) +
  theme(panel.grid.minor = element_line(colour = NA)) +
  theme(panel.grid.major = element_line(colour = NA)) +
  theme(axis.title.x = element_text(colour="black", size=30, family="Times New Roman"),
        axis.text.x  = element_text(angle=0, vjust=0.5, size=30, family="Times New Roman")) +
  theme(axis.title.y = element_text(colour="black", size=30, family="Times New Roman"),
        axis.text.y  = element_text(angle=45, vjust=0.2, size=30, family="Times New Roman")) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(colour="black", face="italic", size = 30, family="Times New Roman")) +
  theme(axis.line = element_line(colour = "black", size = 0.2, linetype = "solid")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x = "Site, Distance from plant")
plot_otu

?facet_grid

cairo_pdf(paste("OTU", Sys.Date(),".pdf"), height = 10, width = 20)
plot_otu



dev.off()


write.csv(Cyanoonly_ag_t, paste("Rel_abund", Sys.Date(), ".csv", sep = ""))

#### To get predicted

pure <- read.csv(file.choose())
# pure_culture_191107.csv
head(pure)

pure_ag<- aggregate(pure$d13C, by=list(Site=pure$Site, OTU=pure$OTU), FUN=mean)

pure_ag_J<-pure_ag[which(pure_ag$Site == "JOR"),]

head(Cyanoonly_ag_ct)
sort(unique(pure$OTU))
sort(unique(Cyanoonly_ag_ct$OTU))
obpr<- left_join(Cyanoonly_ag_ct, pure_ag_J, by=c('OTU'='OTU'))

head(obpr)

obpr$pred <- obpr$rel_ab * obpr$x

pred <- aggregate(obpr$pred, by=list(Site=obpr$Site.x), FUN=sum)




### To do: compare  % C for pure vs. field

# except that vanessa didn't rarefy... 
Torarefy<- read.table(file.choose(), header = T)
library(GUniFrac)

nat_abun_rarefied <- Rarefy(Torarefy, depth = min(rowSums(Torarefy)))

head(nat_abun_rarefied$otu.tab.rff)

a



#### Figure 3


gA <- ggplotGrob(div_obs)
gB <- ggplotGrob(Srar_obs)
gC <- ggplotGrob(pred_obs)
gD <- ggplotGrob(chlor_reg)
gE <- ggplotGrob(leaf_reg)
gF <- ggplotGrob(root_reg)
# 
# 
# 
gA$widths <- gB$widths
gC$widths <- gB$widths
gD$widths <- gB$widths
gE$widths <- gB$widths
gF$widths <- gB$widths
# 
# gA$heights <- gB$heights
# gC$heights <- gB$heights
# 
# 
library(gridExtra)
library(grid)
library(ggsn)
# 
cairo_pdf(paste("Fig_regs", Sys.Date(), ".pdf", sep = ""), height = 8, width = 7)
plot_grid(gA, gB,gC,gD,gE,gF, ncol = 2)
# grid.arrange(gA,gB,gC, gD,gE,gF, ncol = 2) 
# # grid.arrange(gF, gG, gE, nrow=3)
dev.off()

library(cowplot)

