library(lme4)
library(tidyverse)
library(car)
# library(broom.lm)
library(emmeans)

dat <- read.csv("Combined_Final_Data_05_14.csv")

dim(dat)
str(dat)


#first pass trying to look at some smaller subsets. 
# only24 <- filter(dat, Temperature == 24)
# 
# str(only24)
# 
# #we didnt want to have to worry about modeling out the parent of origin effect with this(at first)
# sameparents <- only24[only24$Maternal_Allele == only24$Paternal_Allele,]
# 
# dim(sameparents)
# 
# #releving this so that wt is the base. 
# sameparents$Maternal_Allele <- as.factor(sameparents$Maternal_Allele)
# sameparents$Maternal_Allele <- relevel(sameparents$Maternal_Allele, ref = "wt")
# 
# #I also don't want the Hybrids. I think it is messing with my model
# #need to think more about these random effects. 
# #I made an assumption about how this was set up for the biological reps
# #I don't really want the three way interaction. but the ^2 thing didn't work.  
# #
# mod1 <- lmer(Wing_Area ~ Maternal_Allele*F1_Background*Sex
#              + (1|Block) 
#              +  (1|Maternal_Allele:F1_Background:Replicate), 
#              data = sameparents[sameparents$F1_Background != "HYBRID",])
# 
# #I don't have a great intuition for this and if it is right. 
# summary(mod1)
# 
# 
# #Nothing looks too wild to me here?
# #Anova(mod1)
# 
# #Now I want to get the contrasts. 
# #fist to look at this. 
# #are these error bars too large?
# plot(emmeans(mod1, ~Maternal_Allele | F1_Background))
# 
# formod <- pairs(emmeans(mod1, ~Maternal_Allele | F1_Background))


#trying it again but better I want to think more clearly about what I am doing with this. Trying to fit the entire model and then talk to Ian. 

#Things in the data: temp, maternal allele/background, paternal allele/background, F1 background (3 levels), block (3), replicate (not all blocks have all replicates), sex, indiv, wing area and some semiquant measures. 
with(dat, table(Block, Replicate))
#Mostly all within a single block. 
with(dat, table(Block, Replicate, F1_Background, Maternal_Allele))

with(dat, table(Maternal_Allele, Paternal_Allele, F1_Background))

#I really dont think I want to mess with compound hets here. Ask ID? Not sure how to even fit that effect and what it would mean for what I care about
sameparents <- dat[dat$Maternal_Allele == dat$Paternal_Allele,]
 
dim(sameparents)
with(sameparents, table(Block, Replicate, F1_Background, Maternal_Allele))
#Not a ton of hybrid data actually. I don't think this is worth trying to estimate. 
with(sameparents, table(F1_Background, Maternal_Allele))


#releving this so that wt is the base. 
sameparents$Maternal_Allele <- as.factor(sameparents$Maternal_Allele)
sameparents$Maternal_Allele <- relevel(sameparents$Maternal_Allele, ref = "wt")

#Is there a strong sex* background effect to worry about? I left it out for now. 
#I don't want this temp interaction at all. 
mod2 <- lmer(Wing_Area ~ Maternal_Allele*F1_Background + Sex + Temperature
                          + (1|Block)
                          +  (1|Maternal_Allele:F1_Background:Replicate),
                          data = sameparents[sameparents$F1_Background != "HYBRID",])

summary(mod2)

#What if I drop those replicate (vial) level effects? 
#The answer here is almost the same. 
mod3 <- lmer(Wing_Area ~ Maternal_Allele*F1_Background  + Sex + Temperature
             + (1|Block),
             data = sameparents[sameparents$F1_Background != "HYBRID",])

summary(mod3)

#Does this make sense? I think so. 
plot(emmeans(mod2, ~ Maternal_Allele | F1_Background))

formod <- pairs(emmeans(mod2, ~Maternal_Allele | F1_Background))

head(formod)

dat.out <- data.frame(formod)
dat.wt <- dat.out[grep("^wt", dat.out$contrast),]

#pulling the allele to make things easier later. 
dat.wt$mutant.allele <- substring(dat.wt$contrast, 6)
dat.wt$mutant.allele

#saving for Arteen 
#One thing to note is that these are all positive, when they should be neg? 
write.csv(dat.wt2, "wingSize_modelContrasts_forRNAseqMod_update.csv", row.names = FALSE)

#I don't think this is a GREAT idea but what if we average over sex?




