library(lme4)
library(tidyverse)
library(car)
# library(broom.lm)
library(emmeans)
library(glmmTMB)

dat <- read.csv("Combined_Final_Data_05_14.csv")

dim(dat)
str(dat)

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

#log transforming size
sameparents$log_wingSize <- log2(sameparents$Wing_Area)

str(sameparents)

sameparents[,2:9] <- lapply(sameparents[,2:9], as.factor)
str(sameparents)


####This is what we should actually use based on a convo with ID about this data. 
real.mod.semiQ <- lm(Semiqt_Num ~ Maternal_Allele*F1_Background, 
                     data = only24.m[only24.m$F1_Background != "HYBRID",])

summary(real.mod.semiQ)

#Extract these values and send these to Arteen. 
plot(emmeans(real.mod.semiQ, ~Maternal_Allele | F1_Background))


forAT <- data.frame(emmeans(real.mod.semiQ, ~Maternal_Allele | F1_Background))

write.csv(forAT, "final_alleleEffects_semiQscale.csv", row.names = FALSE)

forAT2 <- data.frame(emmeans(real.mod.semiQ, ~Maternal_Allele))
write.csv(forAT2, "final_alleleEffects_semiQscale_noBackground.csv", row.names = FALSE)

