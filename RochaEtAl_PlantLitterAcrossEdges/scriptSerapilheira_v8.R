library(mgcv) # For adjusting generalized additive models
library(gamm4) # For adjusting generalized additive mixed models
library(lme4) # For adjusting generalized linear mixed models
library(bbmle) # For model selection
library(statmod)
library(gplots) # For some figures
library(ggplot2) # For some other figures
library(PMCMRplus) # For comparisons among fractions
library(tidyverse) # For data reshaping


# Set the working directory and load the data

setwd("C:/Users/55759/Desktop/Artigo_litter_mestrado/analises")
setwd("e:/Pavel/Profissional/Pesquisa/MyPapers_Working/Janaine_2019_EILitterUna/v6_PlantEcol/")

dir.main <- getwd()

dados <- read.table("dadosSerapilheira.txt", header=T, sep="\t")

str(dados)

# Not currently used code for exploratory figures

# par(mfrow=c(2,2), mar=c(3,3,2,2), oma=c(2,2,2,2))
# plot(Mass_leaves ~ Distance, data=dados, pch=21, bg=Transect, xlab="", 
# 	ylab="", main="Folhas")
# abline(v=0, col="gray", lty=3)
# plot(Mass_fwd ~ Distance, data=dados, pch=21, bg=Transect, xlab="", 
# 	ylab="", main="Galhos")
# abline(v=0, col="gray", lty=3)
# plot(Mass_misc ~ Distance, data=dados, pch=21, bg=Transect, xlab="", 
# 	ylab="", main="Miscelânia")
# abline(v=0, col="gray", lty=3)
# plot(Mass_total ~ Distance, data=dados, pch=21, bg=Transect, xlab="", 
# 	ylab="", main="Total")
# abline(v=0, col="gray", lty=3)
# mtext(side=1, text="Distância", outer=T)
# mtext(side=2, text="Biomassa (g)", outer=T)

# Calculate total amount of litter per class for each plot

sum.leaves <- aggregate(Mass_leaves ~ Transect + Distance, FUN=sum, data=dados)
sum.fwd <- aggregate(Mass_fwd ~ Transect + Distance, FUN=sum, data=dados)
sum.misc <- aggregate(Mass_misc ~ Transect + Distance, FUN=sum, data=dados)
sum.total <- aggregate(Mass_total ~ Transect + Distance, FUN=sum, data=dados)

dados.sum <- merge(sum.leaves, sum.fwd)
dados.sum <- merge(dados.sum, sum.misc)
dados.sum <- merge(dados.sum, sum.total)

str(dados.sum)

# Divide the values by 100 to convert from grams per meter to tons per hectare.

dados.sum$Mass_leaves <- dados.sum$Mass_leaves/100
dados.sum$Mass_fwd <- dados.sum$Mass_fwd/100
dados.sum$Mass_misc <- dados.sum$Mass_misc/100
dados.sum$Mass_total <- dados.sum$Mass_total/100

str(dados.sum)

# Some more exploratory figures - not currently used
# par(mfrow=c(2,2), mar=c(3,3,2,2), oma=c(2,2,2,2))
# plot(Mass_leaves ~ Distance, data=dados.sum, pch=21, 
# 	bg=Transect, xlab="", ylab="", main="Folhas")
# abline(v=0, col="gray", lty=3)
# plot(Mass_fwd ~ Distance, data=dados.sum, pch=21, 
# 	bg=Transect, xlab="", ylab="", main="Galhos")
# abline(v=0, col="gray", lty=3)
# plot(Mass_misc ~ Distance, data=dados.sum, pch=21, 
# 	bg=Transect, xlab="", ylab="", main="Miscelânia")
# abline(v=0, col="gray", lty=3)
# plot(Mass_total ~ Distance, data=dados.sum, pch=21, 
# 	bg=Transect, xlab="", ylab="", main="Total")
# abline(v=0, col="gray", lty=3)
# mtext(side=1, text="Distância", outer=T)
# mtext(side=2, text="Biomassa total (t/ha)", outer=T)

# First step of data analysis - adjusting models for the entire transects
# In this analyses, the transects are treated as a whole.
# We adjust the following models:
# Null: no effects of side (burnt vs unburnt) or distance
# Cat_edgeAsForest: difference between burnt and unburnt; edge classifed
##  as forest
# Cat_edgeAsBurnt: same as above but edge classified as burnt
# GAM: generalized additive model, testing for a gradual change along the
##  transects
# GamCat_edgeAsForest: gradual change as well as abrupt differences between
##  the forest and burnt area, edge classifed as forest
# GamCat_edgeAsBurnt: same as above but edge classified as burnt

# Descriptive analysis and comparisons between litter portions

dados.sum.comp <- dados.sum
dados.sum.comp$Collector <- with(dados.sum.comp, paste("T",Transect,
	"_Dist",formatC(Distance, flag="+", zero.print=TRUE),
	sep=""))
dados.sum.comp.long <- pivot_longer(data=dados.sum.comp,
	cols=Mass_leaves:Mass_total, names_to="Fraction",
	values_to="Biomass")
dados.sum.comp.long$Location <- as.factor(with(dados.sum.comp.long,
	ifelse(Distance < 0, "Burnt", ifelse(Distance == 0, "Edge", "Forest"))))


# Calculate descriptive statistics

aggregate(Biomass ~ Fraction + Location, data=dados.sum.comp.long,
	FUN=function(x) round(mean(x),3))

aggregate(Biomass ~ Fraction + Location, data=dados.sum.comp.long,
	FUN=function(x) round(sd(x),3))

aggregate(Biomass ~ Fraction + Location, data=dados.sum.comp.long,
	FUN=length)

# Comparisons among fractions

dados.sum.comp.long2 <- subset(dados.sum.comp.long, Fraction != "Mass_total")

# Fire

friedman.test(Biomass~Fraction|Collector, 
	data=subset(dados.sum.comp.long2, Distance<0))

frdAllPairsSiegelTest(y=subset(dados.sum.comp.long2, Distance<0)$Biomass,
	groups=subset(dados.sum.comp.long2, Distance<0)$Fraction,
	blocks=subset(dados.sum.comp.long2, Distance<0)$Collector,
	p.adjust = "bonferroni")

# Forest

friedman.test(Biomass~Fraction|Collector, 
	data=subset(dados.sum.comp.long2, Distance>0))

frdAllPairsSiegelTest(y=subset(dados.sum.comp.long2, Distance>0)$Biomass,
	groups=subset(dados.sum.comp.long2, Distance>0)$Fraction,
	blocks=subset(dados.sum.comp.long2, Distance>0)$Collector,
	p.adjust = "bonferroni")


# Edge

friedman.test(Biomass~Fraction|Collector, 
	data=subset(dados.sum.comp.long2, Distance==0))

frdAllPairsSiegelTest(y=subset(dados.sum.comp.long2, Distance==0)$Biomass,
	groups=subset(dados.sum.comp.long2, Distance==0)$Fraction,
	blocks=subset(dados.sum.comp.long2, Distance==0)$Collector,
	p.adjust = "bonferroni")



# Model selection


# Create columns showing whether each plot is forst or burnt area
dados.sum$Side_edgeAsForest <- as.factor(ifelse(dados.sum$Distance < 0, 
	"Burnt", "Forest"))
dados.sum$Side_edgeAsBurnt <- as.factor(ifelse(dados.sum$Distance < 1, 
	"Burnt", "Forest"))

# Models for leaf litter

mod.leaves.null <- glmer(Mass_leaves ~ 1 + (1 | Transect), data=dados.sum, 
	family=Gamma)
mod.leaves.cat_edgeAsForest <- glmer(Mass_leaves ~ Side_edgeAsForest + 
	(1|Transect), data=dados.sum, family=Gamma)
mod.leaves.cat_edgeAsBurnt <- glmer(Mass_leaves ~ Side_edgeAsBurnt + 
	(1|Transect), data=dados.sum, family=Gamma)
mod.leaves.gam <- gamm4(Mass_leaves ~ s(Distance, fx=F, k=4), 
	random=~(1|Transect), data=dados.sum, family=Gamma)
mod.leaves.gamcat_edgeAsForest <- gamm4(Mass_leaves ~ s(Distance, fx=F, k=4) + 
	Side_edgeAsForest, random=~(1|Transect), data=dados.sum, family=Gamma)
mod.leaves.gamcat_edgeAsBurnt <- gamm4(Mass_leaves ~ s(Distance, fx=F, k=4) + 
	Side_edgeAsBurnt, random=~(1|Transect), data=dados.sum, family=Gamma)

AICctab(mod.leaves.null, mod.leaves.cat_edgeAsForest, 
	mod.leaves.cat_edgeAsBurnt, mod.leaves.gam$mer, 
	mod.leaves.gamcat_edgeAsForest$mer, mod.leaves.gamcat_edgeAsBurnt$mer)

# The best model is the categorical model with edge classified as burnt.

# Models for fine woody debris (FWD)

# First adjust a null model to selected the optimal parameter p of the
##  Tweedie distribution

mod.fwd.foo <- gam(Mass_fwd ~ 1, data=dados.sum, family=tw)

summary(mod.fwd.foo)

fam <- summary(mod.fwd.foo)$family$family
fam <- strsplit(fam, "=")[[1]][2]
fam <- gsub("\\)","",fam)
fam <- as.numeric(fam)

mod.fwd.null <- gamm4(Mass_fwd ~ 1, random=~(1|Transect), data=dados.sum, 
	family=Tweedie(p=fam))
mod.fwd.cat_edgeAsForest <- gamm4(Mass_fwd ~ Side_edgeAsForest, 
	random=~(1|Transect), data=dados.sum, family=Tweedie(p=fam))
mod.fwd.cat_edgeAsBurnt <- gamm4(Mass_fwd ~ Side_edgeAsBurnt, 
	random=~(1|Transect), data=dados.sum, family=Tweedie(p=fam))
mod.fwd.gam <- gamm4(Mass_fwd ~ s(Distance, fx=F, k=4), 
	random=~(1|Transect), data=dados.sum, family=Tweedie(p=fam))
mod.fwd.gamcat_edgeAsForest <- gamm4(Mass_fwd ~ s(Distance, fx=F, k=4) + 
	Side_edgeAsForest, random=~(1|Transect), data=dados.sum, 
	family=Tweedie(p=fam))
mod.fwd.gamcat_edgeAsBurnt <- gamm4(Mass_fwd ~ s(Distance, fx=F, k=4) + 
	Side_edgeAsBurnt, random=~(1|Transect), data=dados.sum, 
	family=Tweedie(p=fam))

AICctab(mod.fwd.null$mer, mod.fwd.cat_edgeAsForest$mer, 
	mod.fwd.cat_edgeAsBurnt$mer, mod.fwd.gam$mer, 
	mod.fwd.gamcat_edgeAsForest$mer, mod.fwd.gamcat_edgeAsBurnt$mer)


# The best nodel is the null model

### Misc ###

# Determine the parameter of the Tweedie distribution

mod.misc.foo <- gam(Mass_misc ~ 1, data=dados.sum, family=tw)

fam <- summary(mod.misc.foo)$family$family
fam <- strsplit(fam, "=")[[1]][2]
fam <- gsub("\\)","",fam)
fam <- as.numeric(fam)

mod.misc.null <- gamm4(Mass_misc ~ 1, random=~(1|Transect), data=dados.sum, 
	family=Tweedie(p=fam))
mod.misc.cat_edgeAsForest <- gamm4(Mass_misc ~ Side_edgeAsForest, 
	random=~(1|Transect), data=dados.sum, family=Tweedie(p=fam))
mod.misc.cat_edgeAsBurnt <- gamm4(Mass_misc ~ Side_edgeAsBurnt, 
	random=~(1|Transect), data=dados.sum, family=Tweedie(p=fam))
mod.misc.gam <- gamm4(Mass_misc ~ s(Distance, fx=F, k=4), 
	random=~(1|Transect), data=dados.sum, family=Tweedie(p=fam))
mod.misc.gamcat_edgeAsForest <- gamm4(Mass_misc ~ s(Distance, fx=F, k=4) + 
	Side_edgeAsForest, random=~(1|Transect), data=dados.sum, 
	family=Tweedie(p=fam))
mod.misc.gamcat_edgeAsBurnt <- gamm4(Mass_misc ~ s(Distance, fx=F, k=4) + 
	Side_edgeAsBurnt, random=~(1|Transect), data=dados.sum, 
	family=Tweedie(p=fam))

AICctab(mod.misc.null$mer, mod.misc.cat_edgeAsForest$mer, 
	mod.misc.cat_edgeAsBurnt$mer, mod.misc.gam$mer, 
	mod.misc.gamcat_edgeAsForest$mer, mod.misc.gamcat_edgeAsBurnt$mer,
	sort=FALSE)

# The best model is the categorical model, with edge classified as burnt

# Models for total leaf litter


mod.total.null <- glmer(Mass_total ~ 1 + (1 | Transect), data=dados.sum, 
	family=Gamma)
mod.total.cat_edgeAsForest <- glmer(Mass_total ~ Side_edgeAsForest + 
	(1|Transect), data=dados.sum, family=Gamma)
mod.total.cat_edgeAsBurnt <- glmer(Mass_total ~ Side_edgeAsBurnt + 
	(1|Transect), data=dados.sum, family=Gamma)
mod.total.gam <- gamm4(Mass_total ~ s(Distance, fx=F, k=4), 
	random=~(1|Transect), data=dados.sum, family=Gamma)
mod.total.gamcat_edgeAsForest <- gamm4(Mass_total ~ s(Distance, fx=F, k=4) + 
	Side_edgeAsForest, random=~(1|Transect), data=dados.sum, family=Gamma)
mod.total.gamcat_edgeAsBurnt <- gamm4(Mass_total ~ s(Distance, fx=F, k=4) + 
	Side_edgeAsBurnt, random=~(1|Transect), data=dados.sum, family=Gamma)

AICctab(mod.total.null, mod.total.cat_edgeAsForest, 
	mod.total.cat_edgeAsBurnt, mod.total.gam$mer, 
	mod.total.gamcat_edgeAsForest$mer, mod.total.gamcat_edgeAsBurnt$mer,
	sort=FALSE)

# The best model is the categorial model with edge classifed as burnt.

# Thus, there is less leaf litter, miscelaneous litter and total litter 
##  in the burnt area, with edge being more similar to the burnt area than
##  to the unburnt forest. There is no pattern in FWD.

# Calculate mean values per side, for the figures

mean.leaves.burnt <- mean(subset(dados.sum, Side_edgeAsBurnt=="Burnt")$
	Mass_leaves)
mean.leaves.forest <- mean(subset(dados.sum, Side_edgeAsBurnt=="Forest")$
	Mass_leaves)

mean.fwd <- mean(dados.sum$Mass_fwd)

mean.misc.burnt <- mean(subset(dados.sum, Side_edgeAsBurnt=="Burnt")$
	Mass_misc)
mean.misc.forest <- mean(subset(dados.sum, Side_edgeAsBurnt=="Forest")$
	Mass_misc)

mean.total.burnt <- mean(subset(dados.sum, Side_edgeAsBurnt=="Burnt")$
	Mass_total)
mean.total.forest <- mean(subset(dados.sum, Side_edgeAsBurnt=="Forest")$
	Mass_total)


# Figure with model selection results


png(filename="Figure2_modSel_v6.png", height=20, width=20, res=300, unit="cm")
par(mfrow=c(2,2), mar=c(3,3,2,2), oma=c(2,2,2,2))
plot(Mass_leaves ~ Distance, data=dados.sum, pch=16, 
	col=alpha("black",0.5), xlab="", ylab="", main="a. Leaf litter")
abline(v=0, col="gray", lty=3, lwd=2)
segments(x0=-150, x1=5, y0=mean.leaves.burnt, y1=mean.leaves.burnt, 
	col="black")
segments(x0=5, x1=150, y0=mean.leaves.forest, y1=mean.leaves.forest, 
	col="black")
plot(Mass_fwd ~ Distance, data=dados.sum, pch=16, col=alpha("black",0.5), 
	xlab="", ylab="", main="b. Fine woody debris (FWD)")
abline(v=0, col="gray", lty=3, lwd=2)
segments(x0=-150, x1=150, y0=mean.fwd, y1=mean.fwd, col="black")
plot(Mass_misc ~ Distance, data=dados.sum, pch=16, col=alpha("black",0.5), 
	xlab="", ylab="", main="c. Miscelaneous")
abline(v=0, col="gray", lty=3, lwd=2)
segments(x0=-150, x1=5, y0=mean.misc.burnt, y1=mean.misc.burnt, col="black")
segments(x0=5, x1=150, y0=mean.misc.forest, y1=mean.misc.forest, col="black")
plot(Mass_total ~ Distance, data=dados.sum, pch=16, col=alpha("black",0.5), 
	xlab="", ylab="", main="d. Total")
abline(v=0, col="gray", lty=3, lwd=2)
segments(x0=-150, x1=5, y0=mean.total.burnt, y1=mean.total.burnt, col="black")
segments(x0=5, x1=150, y0=mean.total.forest, y1=mean.total.forest, col="black")
mtext(side=1, text="Distance (m)", outer=T)
mtext(side=2, text="Annual litterfall (t/ha)", outer=T)
dev.off()

# Second step of the analysis - monthly analysis.
# The objective here is to see whether these patterns - the best models
##  selected - is something consisted throughout the year.
# For this, we performed this same model selection for each month separately
##  and synthesized their results.

colnames(dados)
unique(dados$Month)

dados$Month <- ordered(dados$Month, levels=c("Sep", "Oct", "Nov", "Dec", 
	"Jan","Feb","Mar", "Apr", "May", "Jun", "Jul", "Aug"))
dados$Side_edgeAsForest <- as.factor(ifelse(dados$Distance < 0, 
	"Burnt", "Forest"))
dados$Side_edgeAsBurnt <- as.factor(ifelse(dados$Distance < 1, 
	"Burnt", "Forest"))

# Convert into t/ha

dadost <- dados


dadost$Mass_leaves <- dadost$Mass_leaves/100
dadost$Mass_fwd <- dadost$Mass_fwd/100
dadost$Mass_misc <- dadost$Mass_misc/100
dadost$Mass_total <- dadost$Mass_total/100

str(dadost)

months <- levels(dadost$Month)

results.modSel <- list()

for(i in 1:12) {
	print(paste(i, "of 12"))
	results.modSel[[i]] <- list()
	names(results.modSel)[i] <- months[i]
	dadost.i <- subset(dadost, Month==months[[i]])
	mod.leaves.foo <- gam(Mass_leaves ~ 1, data=dadost.i, family=tw)
	fam <- summary(mod.leaves.foo)$family$family
	fam <- strsplit(fam, "=")[[1]][2]
	fam <- gsub("\\)","",fam)
	fam <- as.numeric(fam)
	mod.leaves.null <- gamm4(Mass_leaves ~ 1, random=~(1|Transect), 
		data=dadost.i, family=Tweedie(p=fam))
	mod.leaves.cat_edgeAsForest <- gamm4(Mass_leaves ~ Side_edgeAsForest, 
		random=~(1|Transect), data=dadost.i, family=Tweedie(fam))
	mod.leaves.cat_edgeAsBurnt <- gamm4(Mass_leaves ~ Side_edgeAsBurnt, 
		random=~(1|Transect), data=dadost.i, family=Tweedie(fam))
	mod.leaves.gam <- gamm4(Mass_leaves ~ s(Distance, fx=F, k=4), 
		random=~(1|Transect), data=dadost.i, family=Tweedie(fam))
	mod.leaves.gamcat_edgeAsForest <- gamm4(Mass_leaves ~ 
		s(Distance, fx=F, k=4) + Side_edgeAsForest, random=~(1|Transect), 
		data=dadost.i, family=Tweedie(fam))
	mod.leaves.gamcat_edgeAsBurnt <- gamm4(Mass_leaves ~ 
		s(Distance, fx=F, k=4) + Side_edgeAsBurnt, random=~(1|Transect), 
		data=dadost.i, family=Tweedie(fam))
	results.modSel[[i]][[1]] <- AICctab(mod.leaves.null$mer, 
		mod.leaves.cat_edgeAsForest$mer, mod.leaves.cat_edgeAsBurnt$mer, 
		mod.leaves.gam$mer, mod.leaves.gamcat_edgeAsForest$mer, 
		mod.leaves.gamcat_edgeAsBurnt$mer)
	
	mod.fwd.foo <- gam(Mass_fwd ~ 1, data=dadost.i, family=tw)
	fam <- summary(mod.fwd.foo)$family$family
	fam <- strsplit(fam, "=")[[1]][2]
	fam <- gsub("\\)","",fam)
	fam <- as.numeric(fam)
	mod.fwd.null <- gamm4(Mass_fwd ~ 1, random=~(1|Transect), data=dadost.i, 
		family=Tweedie(p=fam))
	mod.fwd.cat_edgeAsForest <- gamm4(Mass_fwd ~ Side_edgeAsForest, 
		random=~(1|Transect), data=dadost.i, family=Tweedie(fam))
	mod.fwd.cat_edgeAsBurnt <- gamm4(Mass_fwd ~ Side_edgeAsBurnt, 
		random=~(1|Transect), data=dadost.i, family=Tweedie(fam))
	mod.fwd.gam <- gamm4(Mass_fwd ~ s(Distance, fx=F, k=4), 
		random=~(1|Transect), data=dadost.i, family=Tweedie(fam))
	mod.fwd.gamcat_edgeAsForest <- gamm4(Mass_fwd ~ s(Distance, fx=F, k=4) + 
		Side_edgeAsForest, random=~(1|Transect), data=dadost.i, 
		family=Tweedie(fam))
	mod.fwd.gamcat_edgeAsBurnt <- gamm4(Mass_fwd ~ s(Distance, fx=F, k=4) + 
		Side_edgeAsBurnt, random=~(1|Transect), data=dadost.i, 
		family=Tweedie(fam))
	results.modSel[[i]][[2]] <- AICctab(mod.fwd.null$mer, 
		mod.fwd.cat_edgeAsForest$mer, mod.fwd.cat_edgeAsBurnt$mer, 
		mod.fwd.gam$mer, mod.fwd.gamcat_edgeAsForest$mer, 
		mod.fwd.gamcat_edgeAsBurnt$mer)

	mod.misc.foo <- gam(Mass_misc ~ 1, data=dadost.i, family=tw)
	#str(summary(mod.misc.foo))
	fam <- summary(mod.misc.foo)$family$family
	fam <- strsplit(fam, "=")[[1]][2]
	fam <- gsub("\\)","",fam)
	fam <- as.numeric(fam)
	mod.misc.null <- gamm4(Mass_misc ~ 1, random=~(1|Transect), data=dadost.i, 
		family=Tweedie(p=fam))
	mod.misc.cat_edgeAsForest <- gamm4(Mass_misc ~ Side_edgeAsForest, 
		random=~(1|Transect), data=dadost.i, family=Tweedie(fam))
	mod.misc.cat_edgeAsBurnt <- gamm4(Mass_misc ~ Side_edgeAsBurnt, 
		random=~(1|Transect), data=dadost.i, family=Tweedie(fam))
	mod.misc.gam <- gamm4(Mass_misc ~ s(Distance, fx=F, k=4), 
		random=~(1|Transect), data=dadost.i, family=Tweedie(fam))
	mod.misc.gamcat_edgeAsForest <- gamm4(Mass_misc ~ s(Distance, fx=F, k=4) + 
		Side_edgeAsForest, random=~(1|Transect), data=dadost.i, 
		family=Tweedie(fam))
	mod.misc.gamcat_edgeAsBurnt <- gamm4(Mass_misc ~ s(Distance, fx=F, k=4) + 
		Side_edgeAsBurnt, random=~(1|Transect), data=dadost.i, 
		family=Tweedie(fam))
	results.modSel[[i]][[3]] <- AICctab(mod.misc.null$mer, 
		mod.misc.cat_edgeAsForest$mer, mod.misc.cat_edgeAsBurnt$mer, 
		mod.misc.gam$mer, mod.misc.gamcat_edgeAsForest$mer, 
		mod.misc.gamcat_edgeAsBurnt$mer)

	mod.total.foo <- gam(Mass_total ~ 1, data=dadost.i, family=tw)
	#str(summary(mod.total.foo))
	fam <- summary(mod.total.foo)$family$family
	fam <- strsplit(fam, "=")[[1]][2]
	fam <- gsub("\\)","",fam)
	fam <- as.numeric(fam)
	mod.total.null <- gamm4(Mass_total ~ 1, random=~(1|Transect), 
		data=dadost.i, family=Tweedie(p=fam))
	mod.total.cat_edgeAsForest <- gamm4(Mass_total ~ Side_edgeAsForest, 
		random=~(1|Transect), data=dadost.i, family=Tweedie(fam))
	mod.total.cat_edgeAsBurnt <- gamm4(Mass_total ~ Side_edgeAsBurnt, 
		random=~(1|Transect), data=dadost.i, family=Tweedie(fam))
	mod.total.gam <- gamm4(Mass_total ~ s(Distance, fx=F, k=4), 
		random=~(1|Transect), data=dadost.i, family=Tweedie(fam))
	mod.total.gamcat_edgeAsForest <- gamm4(Mass_total ~ 
		s(Distance, fx=F, k=4) + Side_edgeAsForest, random=~(1|Transect), 
		data=dadost.i, family=Tweedie(fam))
	mod.total.gamcat_edgeAsBurnt <- gamm4(Mass_total ~ s(Distance, fx=F, k=4) + 
		Side_edgeAsBurnt, random=~(1|Transect), data=dadost.i, 
		family=Tweedie(fam))
	results.modSel[[i]][[4]] <- AICctab(mod.total.null$mer, 
		mod.total.cat_edgeAsForest$mer, mod.total.cat_edgeAsBurnt$mer, 
		mod.total.gam$mer, mod.total.gamcat_edgeAsForest$mer, 
		mod.total.gamcat_edgeAsBurnt$mer)

	names(results.modSel[[i]]) <- c("Leaves", "FWD", "Misc", "Total")

}
results.modSel

# In which was there more of each fraction per month?

aggregate(Mass_leaves ~ Side_edgeAsBurnt + Month, data=dadost, FUN=mean)
aggregate(Mass_fwd ~ Side_edgeAsBurnt + Month, data=dadost, FUN=mean)
aggregate(Mass_misc ~ Side_edgeAsBurnt + Month, data=dadost, FUN=mean)
aggregate(Mass_misc ~ Side_edgeAsForest + Month, data=dadost, FUN=mean)
aggregate(Mass_total ~ Side_edgeAsBurnt + Month, data=dadost, FUN=mean)
aggregate(Mass_total ~ Side_edgeAsForest + Month, data=dadost, FUN=mean)


# Select the best model for each month and fraction

bestModels <- matrix(nrow=12, ncol=4)
rownames(bestModels) <- months
colnames(bestModels) <- c("Leaves", "FWD", "Misc", "Total")

for(i in 1:length(months)) {
	for(j in 1:4) {
		results.i.j <- as.data.frame(results.modSel[[i]][[j]])
		rownames(results.i.j) <- gsub("mod.leaves.", "",rownames(results.i.j))
		rownames(results.i.j) <- gsub("mod.fwd.", "",rownames(results.i.j))
		rownames(results.i.j) <- gsub("mod.misc.", "",rownames(results.i.j))
		rownames(results.i.j) <- gsub("mod.total.", "",rownames(results.i.j))
		rownames(results.i.j) <- gsub("\\$mer", "",rownames(results.i.j))
		results.foo <- subset(results.i.j, dAICc<=2)
		df.min <- min(results.foo$df)
		best.model <- rownames(subset(results.foo, df==df.min))[1]
		bestModels[i,j] <- best.model
	}
}

bestModels

View(bestModels)



### Separate analyses for the burnt and unburnt sides, to check whether
###  there are edge effects that were not detected in the model selection
###  for being at only one side of the edge


### Separate the data

dados.sum.forest <- subset(dados.sum, Side_edgeAsForest == "Forest")
dados.sum.burnt <- subset(dados.sum, Side_edgeAsBurnt == "Burnt")

### Generalized additive models for the forest

mod.leaves.gam.forest <- gamm4(Mass_leaves ~ s(Distance, fx=F, k=4), 
	random=~(1|Transect), 
	data=dados.sum.forest, family=Gamma)

mod.fwd.foo <- gam(Mass_fwd ~ 1, data=dados.sum.forest, family=tw)
fam <- summary(mod.fwd.foo)$family$family
fam <- strsplit(fam, "=")[[1]][2]
fam <- gsub("\\)","",fam)
fam <- as.numeric(fam)
mod.fwd.gam.forest <- gamm4(Mass_fwd ~ s(Distance, fx=F, k=4), 
	random=~(1|Transect), 
	data=dados.sum.forest, family=Tweedie(fam))

mod.misc.foo <- gam(Mass_misc ~ 1, data=dados.sum.forest, family=tw)
fam <- summary(mod.misc.foo)$family$family
fam <- strsplit(fam, "=")[[1]][2]
fam <- gsub("\\)","",fam)
fam <- as.numeric(fam)
mod.misc.gam.forest <- gamm4(Mass_misc ~ s(Distance, fx=F, k=4), 
	random=~(1|Transect), 
	data=dados.sum.forest, family=Tweedie(p=fam))

mod.total.gam.forest <- gamm4(Mass_total ~ s(Distance, fx=F, k=4), 
	random=~(1|Transect), 
	data=dados.sum.forest, family=Gamma)


summary(mod.leaves.gam.forest$gam) # Significant
summary(mod.fwd.gam.forest$gam) # Non-significant
summary(mod.misc.gam.forest$gam) # Non-significant
summary(mod.total.gam.forest$gam) # Significant


p.gam.leaves.forest <- summary(mod.leaves.gam.forest$gam)$s.table[4]
p.gam.fwd.forest <- summary(mod.fwd.gam.forest$gam)$s.table[4]
p.gam.misc.forest <- summary(mod.misc.gam.forest$gam)$s.table[4]
p.gam.total.forest <- summary(mod.total.gam.forest$gam)$s.table[4]

newdistances <- 0:150

pred.leaves.gam.forest <- predict(mod.leaves.gam.forest$gam, 
	newdata=list(Distance=newdistances), 
	type="response")
pred.fwd.gam.forest <- predict(mod.fwd.gam.forest$gam, 
	newdata=list(Distance=newdistances), 
	type="response")
pred.misc.gam.forest <- predict(mod.misc.gam.forest$gam, 
	newdata=list(Distance=newdistances), 
	type="response")
pred.total.gam.forest <- predict(mod.total.gam.forest$gam, 
	newdata=list(Distance=newdistances), 
	type="response")

png(filename="figure3_GAM_forest.png", height=20, 
	width=20, res=300, unit="cm")
par(mfrow=c(2,2), mar=c(3,3,2,2), oma=c(2,2,2,2))
plot(Mass_leaves ~ Distance, data=dados.sum.forest, pch=16, 
	col=alpha("black",0.5), xlab="", ylab="", main="a. Leaf litter")
lines(pred.leaves.gam.forest ~ newdistances, lwd=1.5,lty=1)
legend("topright", legend=paste("p = ", formatC(p.gam.leaves.forest, digits=3,
	format="f")), bty="n")
plot(Mass_fwd ~ Distance, data=dados.sum.forest, pch=16, 
	col=alpha("black",0.5), xlab="", ylab="", main="b. Fine woody debris (FWD)")
lines(pred.fwd.gam.forest ~ newdistances, lwd=1.5, lty=2)
legend("topright", legend=paste("p = ", formatC(p.gam.fwd.forest, digits=3,
	format="f")), bty="n")
plot(Mass_misc ~ Distance, data=dados.sum.forest, pch=16, 
	col=alpha("black",0.5), xlab="", ylab="", main="c. Miscelaneous")
lines(pred.misc.gam.forest ~ newdistances, lwd=1.5, lty=2)
legend("topright", legend=paste("p = ", formatC(p.gam.misc.forest, digits=3,
	format="f")), bty="n")
plot(Mass_total ~ Distance, data=dados.sum.forest, pch=16, 
	col=alpha("black",0.5), xlab="", ylab="", main="d. Total")
legend("topright", legend=paste("p = ", formatC(p.gam.total.forest, digits=3,
	format="f")), bty="n")
lines(pred.total.gam.forest ~ newdistances, lwd=1.5, lty=1)
mtext(side=1, text="Distance from edge (m)", outer=T)
mtext(side=2, text="Annual litterfall (t/ha)", outer=T)
dev.off()


### Generalized addtive models for the burnt area

dados.sum.burnt$Distance <- abs(dados.sum.burnt$Distance)

mod.leaves.gam.burnt <- gamm4(Mass_leaves ~ s(Distance, fx=F, k=4), 
	random=~(1|Transect), 
	data=dados.sum.burnt, family=Gamma)

mod.fwd.foo <- gam(Mass_fwd ~ 1, data=dados.sum.burnt, family=tw)
fam <- summary(mod.fwd.foo)$family$family
fam <- strsplit(fam, "=")[[1]][2]
fam <- gsub("\\)","",fam)
fam <- as.numeric(fam)
mod.fwd.gam.burnt <- gamm4(Mass_fwd ~ s(Distance, fx=F, k=4), 
	random=~(1|Transect), 
	data=dados.sum.burnt, family=Tweedie(p=1.72))

mod.misc.foo <- gam(Mass_misc ~ 1, data=dados.sum.burnt, family=tw)
fam <- summary(mod.misc.foo)$family$family
fam <- strsplit(fam, "=")[[1]][2]
fam <- gsub("\\)","",fam)
fam <- as.numeric(fam)
mod.misc.gam.burnt <- gamm4(Mass_misc ~ s(Distance, fx=F, k=4), 
	random=~(1|Transect), 
	data=dados.sum.burnt, family=Tweedie(p=1.659))

mod.total.gam.burnt <- gamm4(Mass_total ~ s(Distance, fx=F, k=4), 
	random=~(1|Transect), 
	data=dados.sum.burnt, family=Gamma)

summary(mod.leaves.gam.burnt$gam) # Non-significant
summary(mod.fwd.gam.burnt$gam)  # Non-significant
summary(mod.misc.gam.burnt$gam) # Non-significant
summary(mod.total.gam.burnt$gam) # Non-significant


pred.leaves.gam.burnt <- predict(mod.leaves.gam.burnt$gam, 
	newdata=list(Distance=newdistances), 
	type="response")
pred.fwd.gam.burnt <- predict(mod.fwd.gam.burnt$gam, 
	newdata=list(Distance=newdistances), 
	type="response")
pred.misc.gam.burnt <- predict(mod.misc.gam.burnt$gam, 
	newdata=list(Distance=newdistances), 
	type="response")
pred.total.gam.burnt <- predict(mod.total.gam.burnt$gam, 
	newdata=list(Distance=newdistances), 
	type="response")



# png(filename="EIlitter_results_burnt_v6.png", 
# 	height=20, width=20, res=300, unit="cm")
# par(mfrow=c(2,2), mar=c(3,3,2,2), oma=c(2,2,2,2))
# plot(Mass_leaves ~ Distance, data=dados.sum.burnt, pch=16, 
# 	col=alpha("black",0.5), xlab="", ylab="", main="a. Leaf litter")
# lines(pred.leaves.gam.burnt ~ newdistances, lwd=1.5,lty=2)
# plot(Mass_fwd ~ Distance, data=dados.sum.burnt, pch=16, 
# 	col=alpha("black",0.5), xlab="", ylab="", main="b. Fine woody debris (fwd)")
# lines(pred.fwd.gam.burnt ~ newdistances, lwd=1.5, lty=2)
# plot(Mass_misc ~ Distance, data=dados.sum.burnt, pch=16, 
# 	col=alpha("black",0.5), xlab="", ylab="", main="c. Miscelaneous")
# lines(pred.misc.gam.burnt ~ newdistances, lwd=1.5, lty=2)
# plot(Mass_total ~ Distance, data=dados.sum.burnt, pch=16, 
# 	col=alpha("black",0.5), xlab="", ylab="", main="d. Total")
# lines(pred.total.gam.burnt ~ newdistances, lwd=1.5, lty=2)
# mtext(side=1, text="Distance (m)", outer=T)
# mtext(side=2, text="Annual litter production (t/ha)", outer=T)
# dev.off()


# Now, separate analyses for the burnt and unburnt sides for each month


results.gam.forest <- list()
results.gam.burnt <- list()

newdistances.forest <- 0:150
newdistances.burnt <- -150:0

setwd("GAMMperMonth")

for(i in 1:12) {
	print(paste(i, "of 12"))
	# Separate the data
	results.gam.forest[[i]] <- list()
	names(results.gam.forest)[i] <- months[i]
	results.gam.burnt[[i]] <- list()
	names(results.gam.burnt)[i] <- months[i]
	dados.i <- subset(dadost, Month==months[[i]])
	dados.i.forest <- subset(dados.i, Side_edgeAsForest == "Forest")
	dados.i.burnt <- subset(dados.i, Side_edgeAsBurnt == "Burnt")



	png(filename=paste("Results_GAM_perSide_",formatC(i, width=2, flag=0),"_",
		months[i], ifelse(months[i] %in% c("Sep", "Oct", "Nov", "Dec"),
			2017, 2018), ".png",sep=""), 
		height=30, width=25, unit="cm", res=300)
	par(mfcol=c(4,2), mar=c(4,4,3,3), oma=c(3,3,2,2))

	# Models for burnt area

	mod.leaves.burnt.foo <- gam(Mass_leaves ~ 1, data=dados.i.burnt, 
		family=tw)
	fam <- summary(mod.leaves.burnt.foo)$family$family
	fam <- strsplit(fam, "=")[[1]][2]
	fam <- gsub("\\)","",fam)
	fam <- as.numeric(fam)
	results.gam.burnt[[i]][[1]] <- gamm4(Mass_leaves ~ s(Distance, fx=F, k=4), 
		random=~(1|Transect), data=dados.i.burnt, family=Tweedie(fam))

	predicted.leaves.burnt <- predict(results.gam.burnt[[i]][[1]]$gam,
		type="response", newdata=list(Distance=newdistances.burnt))
	p.leaves.burnt <- summary(results.gam.burnt[[i]][[1]]$gam)$s.table[,4]

	plot(Mass_leaves ~ Distance, data=dados.i.burnt, pch=16, 
		bg=alpha("black", 0.5), xlab="", ylab="", 
		main="Burnt area - leaves",
		ylim=c(0,max(dados.i$Mass_leaves)))
	lines(predicted.leaves.burnt ~ newdistances.burnt, type="l",
		lty=ifelse(p.leaves.burnt < 0.05, 1, 2), lwd=1.5)
	text(x=0, y=max(dados.i$Mass_leaves),
		labels=paste("p=", formatC(p.leaves.burnt, 
		digits=2, format="g"), sep=""), pos=2)	

	mod.fwd.burnt.foo <- gam(Mass_fwd ~ 1, data=dados.i.burnt, 
		family=tw)
	fam <- summary(mod.fwd.burnt.foo)$family$family
	fam <- strsplit(fam, "=")[[1]][2]
	fam <- gsub("\\)","",fam)
	fam <- as.numeric(fam)
	results.gam.burnt[[i]][[2]] <- gamm4(Mass_fwd ~ s(Distance, fx=F, k=4), 
		random=~(1|Transect), data=dados.i.burnt, family=Tweedie(fam))

	predicted.fwd.burnt <- predict(results.gam.burnt[[i]][[2]]$gam,
		type="response", newdata=list(Distance=newdistances.burnt))
	p.fwd.burnt <- summary(results.gam.burnt[[i]][[2]]$gam)$s.table[,4]

	plot(Mass_fwd ~ Distance, data=dados.i.burnt, pch=16, 
		bg=alpha("black", 0.5), xlab="", ylab="", 
		main="Burnt area - fwd",
		ylim=c(0,max(dados.i$Mass_fwd)))
	lines(predicted.fwd.burnt ~ newdistances.burnt, type="l",
		lty=ifelse(p.fwd.burnt < 0.05, 1, 2), lwd=1.5)
	text(x=0, y=max(dados.i$Mass_fwd),
		labels=paste("p=", formatC(p.fwd.burnt, 
		digits=2, format="g"), sep=""), pos=2)	

	mod.misc.burnt.foo <- gam(Mass_misc ~ 1, data=dados.i.burnt,
	 family=tw)
	#str(summary(mod.misc.foo))
	fam <- summary(mod.misc.burnt.foo)$family$family
	fam <- strsplit(fam, "=")[[1]][2]
	fam <- gsub("\\)","",fam)
	fam <- as.numeric(fam)
	results.gam.burnt[[i]][[3]] <- gamm4(Mass_misc ~ s(Distance, fx=F, k=4), 
		random=~(1|Transect), data=dados.i.burnt, family=Tweedie(fam))

	predicted.misc.burnt <- predict(results.gam.burnt[[i]][[3]]$gam,
		type="response", newdata=list(Distance=newdistances.burnt))
	p.misc.burnt <- summary(results.gam.burnt[[i]][[3]]$gam)$s.table[,4]

	plot(Mass_misc ~ Distance, data=dados.i.burnt, pch=16, 
		bg=alpha("black", 0.5), xlab="", ylab="", 
		main="Burnt area - misc",
		ylim=c(0,max(dados.i$Mass_misc)))
	lines(predicted.misc.burnt ~ newdistances.burnt, type="l",
		lty=ifelse(p.misc.burnt < 0.05, 1, 2), lwd=1.5)
	text(x=0, y=max(dados.i$Mass_misc),
		labels=paste("p=", formatC(p.misc.burnt, 
		digits=2, format="g"), sep=""), pos=2)		
	
	mod.total.burnt.foo <- gam(Mass_total ~ 1, data=dados.i.burnt, 
		family=tw)
	#str(summary(mod.total.foo))
	fam <- summary(mod.total.burnt.foo)$family$family
	fam <- strsplit(fam, "=")[[1]][2]
	fam <- gsub("\\)","",fam)
	fam <- as.numeric(fam)
	results.gam.burnt[[i]][[4]] <- gamm4(Mass_total ~ s(Distance, fx=F, k=4), 
		random=~(1|Transect), data=dados.i.burnt, family=Tweedie(fam))

	predicted.total.burnt <- predict(results.gam.burnt[[i]][[4]]$gam,
		type="response", newdata=list(Distance=newdistances.burnt))
	p.total.burnt <- summary(results.gam.burnt[[i]][[4]]$gam)$s.table[,4]

	plot(Mass_total ~ Distance, data=dados.i.burnt, pch=16, 
		bg=alpha("black", 0.5), xlab="", ylab="", 
		main="Burnt area - total",
		ylim=c(0,max(dados.i$Mass_total)))
	lines(predicted.total.burnt ~ newdistances.burnt, type="l",
		lty=ifelse(p.total.burnt < 0.05, 1, 2), lwd=1.5)
	text(x=0, y=max(dados.i$Mass_total),
		labels=paste("p=", formatC(p.total.burnt, 
		digits=2, format="g"), sep=""), pos=2)		

	mtext(side=3, text=paste(formatC(i, width=2, flag=0), "-",
		months[i], ifelse(months[i] %in% c("Sep", "Oct", "Nov", "Dec"),
			2017, 2018), sep=" "), outer=T, cex=1.5)
	mtext(text="Distance from edge (m)", outer=T, side=1)
	mtext(text="Mass (t/ha)", outer=T, side=2)
	
	names(results.gam.burnt[[i]]) <- c("Leaves", "FWD", "Misc", "Total")


	# Models for forest

	mod.leaves.forest.foo <- gam(Mass_leaves ~ 1, data=dados.i.forest, 
		family=tw)
	fam <- summary(mod.leaves.forest.foo)$family$family
	fam <- strsplit(fam, "=")[[1]][2]
	fam <- gsub("\\)","",fam)
	fam <- as.numeric(fam)
	results.gam.forest[[i]][[1]] <- gamm4(Mass_leaves ~ s(Distance, fx=F, k=4), 
		random=~(1|Transect), data=dados.i.forest, family=Tweedie(fam))

	predicted.leaves.forest <- predict(results.gam.forest[[i]][[1]]$gam,
		type="response", newdata=list(Distance=newdistances.forest))
	p.leaves.forest <- summary(results.gam.forest[[i]][[1]]$gam)$s.table[,4]

	plot(Mass_leaves ~ Distance, data=dados.i.forest, pch=16, 
		bg=alpha("black", 0.5), xlab="", ylab="", 
		main="Unburnt forest - leaves",
		ylim=c(0,max(dados.i$Mass_leaves)))
	lines(predicted.leaves.forest ~ newdistances.forest, type="l",
		lty=ifelse(p.leaves.forest < 0.05, 1, 2), lwd=1.5)
	text(x=150, y=max(dados.i$Mass_leaves),
		labels=paste("p=", formatC(p.leaves.forest, 
		digits=2, format="g"), sep=""), pos=2)

	mod.fwd.forest.foo <- gam(Mass_fwd ~ 1, data=dados.i.forest, 
		family=tw)
	fam <- summary(mod.fwd.forest.foo)$family$family
	fam <- strsplit(fam, "=")[[1]][2]
	fam <- gsub("\\)","",fam)
	fam <- as.numeric(fam)
	results.gam.forest[[i]][[2]] <- gamm4(Mass_fwd ~ s(Distance, fx=F, k=4), 
		random=~(1|Transect), data=dados.i.forest, family=Tweedie(fam))

	predicted.fwd.forest <- predict(results.gam.forest[[i]][[2]]$gam,
		type="response", newdata=list(Distance=newdistances.forest))
	p.fwd.forest <- summary(results.gam.forest[[i]][[2]]$gam)$s.table[,4]

	plot(Mass_fwd ~ Distance, data=dados.i.forest, pch=16, 
		bg=alpha("black", 0.5), xlab="", ylab="", 
		main="Unburnt forest - fwd",
		ylim=c(0,max(dados.i$Mass_fwd)))
	lines(predicted.fwd.forest ~ newdistances.forest, type="l",
		lty=ifelse(p.fwd.forest < 0.05, 1, 2), lwd=1.5)
	text(x=150, y=max(dados.i$Mass_fwd),
		labels=paste("p=", formatC(p.fwd.forest, 
		digits=2, format="g"), sep=""), pos=2)

	mod.misc.forest.foo <- gam(Mass_misc ~ 1, data=dados.i.forest,
	 family=tw)
	#str(summary(mod.misc.foo))
	fam <- summary(mod.misc.forest.foo)$family$family
	fam <- strsplit(fam, "=")[[1]][2]
	fam <- gsub("\\)","",fam)
	fam <- as.numeric(fam)
	results.gam.forest[[i]][[3]] <- gamm4(Mass_misc ~ s(Distance, fx=F, k=4), 
		random=~(1|Transect), data=dados.i.forest, family=Tweedie(fam))
	
	predicted.misc.forest <- predict(results.gam.forest[[i]][[3]]$gam,
		type="response", newdata=list(Distance=newdistances.forest))
	p.misc.forest <- summary(results.gam.forest[[i]][[3]]$gam)$s.table[,4]

	plot(Mass_misc ~ Distance, data=dados.i.forest, pch=16, 
		bg=alpha("black", 0.5), xlab="", ylab="", 
		main="Unburnt forest - misc",
		ylim=c(0,max(dados.i$Mass_misc)))
	lines(predicted.misc.forest ~ newdistances.forest, type="l",
		lty=ifelse(p.misc.forest < 0.05, 1, 2), lwd=1.5)
	text(x=150, y=max(dados.i$Mass_misc),
		labels=paste("p=", formatC(p.misc.forest, 
		digits=2, format="g"), sep=""), pos=2)

	mod.total.forest.foo <- gam(Mass_total ~ 1, data=dados.i.forest, 
		family=tw)
	#str(summary(mod.total.foo))
	fam <- summary(mod.total.forest.foo)$family$family
	fam <- strsplit(fam, "=")[[1]][2]
	fam <- gsub("\\)","",fam)
	fam <- as.numeric(fam)
	results.gam.forest[[i]][[4]] <- gamm4(Mass_total ~ s(Distance, fx=F, k=4), 
		random=~(1|Transect), data=dados.i.forest, family=Tweedie(fam))

	predicted.total.forest <- predict(results.gam.forest[[i]][[4]]$gam,
		type="response", newdata=list(Distance=newdistances.forest))
	p.total.forest <- summary(results.gam.forest[[i]][[4]]$gam)$s.table[,4]

	plot(Mass_total ~ Distance, data=dados.i.forest, pch=16, 
		bg=alpha("black", 0.5), xlab="", ylab="", 
		main="Unburnt forest - total",
		ylim=c(0,max(dados.i$Mass_total)))
	lines(predicted.total.forest ~ newdistances.forest, type="l",
		lty=ifelse(p.total.forest < 0.05, 1, 2), lwd=1.5)
	text(x=150, y=max(dados.i$Mass_total),
		labels=paste("p=", formatC(p.total.forest, 
		digits=2, format="g"), sep=""), pos=2)

	names(results.gam.forest[[i]]) <- c("Leaves", "FWD", "Misc", "Total")


	dev.off()

}

> citation("car")                                                                                                                                                                                                                             
                                                                                                                                                                                                                                              
To cite the car package in publications use:                                                                                                                                                                                                  
                                                                                                                                                                                                                                              
  John Fox and Sanford Weisberg (2019). An {R} Companion to Applied                                                                                                                                                                           
  Regression, Third Edition. Thousand Oaks CA: Sage. URL:                                                                                                                                                                                     
  https://socialsciences.mcmaster.ca/jfox/Books/Companion/                                                                                                                                                                                    
                                                                                                                                                                                                                                              
A BibTeX entry for LaTeX users is                                                                                                                                                                                                             
                                                                                                                                                                                                                                              
  @Book{,                                                                                                                                                                                                                                     
    title = {An {R} Companion to Applied Regression},                                                                                                                                                                                         
    edition = {Third},                                                                                                                                                                                                                        
    author = {John Fox and Sanford Weisberg},                                                                                                                                                                                                 
    year = {2019},                                                                                                                                                                                                                            
    publisher = {Sage},                                                                                                                                                                                                                       
    address = {Thousand Oaks {CA}},                                                                                                                                                                                                           
    url = {https://socialsciences.mcmaster.ca/jfox/Books/Companion/},                                                                                                                                                                         
  }                                                                                                                                                                                                                                           
                                                                                                                                                                                                                                              
> citation("ape")                                                                                                                                                                                                                             
                                                                                                                                                                                                                                              
To cite ape in a publication please use:                                                                                                                                                                                                      
                                                                                                                                                                                                                                              
  Paradis E, Schliep K (2019). "ape 5.0: an environment for modern                                                                                                                                                                            
  phylogenetics and evolutionary analyses in R." _Bioinformatics_,                                                                                                                                                                            
  *35*, 526-528. doi:10.1093/bioinformatics/bty633                                                                                                                                                                                            
  <https://doi.org/10.1093/bioinformatics/bty633>.                                                                                                                                                                                            
                                                                                                                                                                                                                                              
A BibTeX entry for LaTeX users is                                                                                                                                                                                                             
                                                                                                                                                                                                                                              
  @Article{,                                                                                                                                                                                                                                  
    title = {ape 5.0: an environment for modern phylogenetics and evolutionary analyses in {R}},                                                                                                                                              
    author = {Emmanuel Paradis and Klaus Schliep},                                                                                                                                                                                            
    journal = {Bioinformatics},                                                                                                                                                                                                               
    year = {2019},                                                                                                                                                                                                                            
    volume = {35},                                                                                                                                                                                                                            
    pages = {526-528},                                                                                                                                                                                                                        
    doi = {10.1093/bioinformatics/bty633},                                                                                                                                                                                                    
  }                                                                                                                                                                                                                                           
                                                                                                                                                                                                                                              
ape is evolving quickly, so you may want to cite its version number                                                                                                                                                                           
(found with 'library(help = ape)' or 'packageVersion("ape")').                                                                                                                                                                                
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
> citation("spdep")                                                                                                                                                                                                                           
Error in citation("spdep") : there is no package called 'spdep'                                                                                                                                                                               
> citation("psych")                                                                                                                                                                                                                           
Error in citation("psych") : there is no package called 'psych'                                                                                                                                                                               
> citation("plyr")                                                                                                                                                                                                                            
                                                                                                                                                                                                                                              
To cite plyr in publications use:                                                                                                                                                                                                             
                                                                                                                                                                                                                                              
  Hadley Wickham (2011). The Split-Apply-Combine Strategy for Data                                                                                                                                                                            
  Analysis. Journal of Statistical Software, 40(1), 1-29. URL                                                                                                                                                                                 
  https://www.jstatsoft.org/v40/i01/.                                                                                                                                                                                                         
                                                                                                                                                                                                                                              
A BibTeX entry for LaTeX users is                                                                                                                                                                                                             
                                                                                                                                                                                                                                              
  @Article{,                                                                                                                                                                                                                                  
    title = {The Split-Apply-Combine Strategy for Data Analysis},                                                                                                                                                                             
    author = {Hadley Wickham},                                                                                                                                                                                                                
    journal = {Journal of Statistical Software},                                                                                                                                                                                              
    year = {2011},                                                                                                                                                                                                                            
    volume = {40},                                                                                                                                                                                                                            
    number = {1},                                                                                                                                                                                                                             
    pages = {1--29},                                                                                                                                                                                                                          
    url = {https://www.jstatsoft.org/v40/i01/},                                                                                                                                                                                               
  }                                                                                                                                                                                                                                           
                                                                                                                                                                                                                                              
> citation("modEVA")                                                                                                                                                                                                                          
Error in citation("modEVA") : there is no package called 'modEVA'                                                                                                                                                                             
> citation("mgcv")                                                                                                                                                                                                                            
                                                                                                                                                                                                                                              
2011 for generalized additive model method; 2016 for beyond exponential                                                                                                                                                                       
family; 2004 for strictly additive GCV based model method and basics of                                                                                                                                                                       
gamm; 2017 for overview; 2003 for thin plate regression splines.                                                                                                                                                                              
                                                                                                                                                                                                                                              
  Wood, S.N. (2011) Fast stable restricted maximum likelihood and                                                                                                                                                                             
  marginal likelihood estimation of semiparametric generalized linear                                                                                                                                                                         
  models. Journal of the Royal Statistical Society (B) 73(1):3-36                                                                                                                                                                             
                                                                                                                                                                                                                                              
  Wood S.N., N. Pya and B. Saefken (2016) Smoothing parameter and model                                                                                                                                                                       
  selection for general smooth models (with discussion). Journal of the                                                                                                                                                                       
  American Statistical Association 111:1548-1575.                                                                                                                                                                                             
                                                                                                                                                                                                                                              
  Wood, S.N. (2004) Stable and efficient multiple smoothing parameter                                                                                                                                                                         
  estimation for generalized additive models. Journal of the American                                                                                                                                                                         
  Statistical Association. 99:673-686.                                                                                                                                                                                                        
                                                                                                                                                                                                                                              
  Wood, S.N. (2017) Generalized Additive Models: An Introduction with R                                                                                                                                                                       
  (2nd edition). Chapman and Hall/CRC.                                                                                                                                                                                                        
                                                                                                                                                                                                                                              
  Wood, S.N. (2003) Thin-plate regression splines. Journal of the Royal                                                                                                                                                                       
  Statistical Society (B) 65(1):95-114.                                                                                                                                                                                                       
                                                                                                                                                                                                                                              
To see these entries in BibTeX format, use 'print(<citation>,                                                                                                                                                                                 
bibtex=TRUE)', 'toBibtex(.)', or set                                                                                                                                                                                                          
'options(citation.bibtex.max=999)'.                                                                                                                                                                                                           
                                                                                                                                                                                                                                              
> citation("MASS")                                                                                                                                                                                                                            
                                                                                                                                                                                                                                              
To cite the MASS package in publications use:                                                                                                                                                                                                 
                                                                                                                                                                                                                                              
  Venables, W. N. & Ripley, B. D. (2002) Modern Applied Statistics with                                                                                                                                                                       
  S. Fourth Edition. Springer, New York. ISBN 0-387-95457-0                                                                                                                                                                                   
                                                                                                                                                                                                                                              
A BibTeX entry for LaTeX users is                                                                                                                                                                                                             
                                                                                                                                                                                                                                              
  @Book{,                                                                                                                                                                                                                                     
    title = {Modern Applied Statistics with S},                                                                                                                                                                                               
    author = {W. N. Venables and B. D. Ripley},                                                                                                                                                                                               
    publisher = {Springer},                                                                                                                                                                                                                   
    edition = {Fourth},                                                                                                                                                                                                                       
    address = {New York},                                                                                                                                                                                                                     
    year = {2002},                                                                                                                                                                                                                            
    note = {ISBN 0-387-95457-0},                                                                                                                                                                                                              
    url = {https://www.stats.ox.ac.uk/pub/MASS4/},                                                                                                                                                                                            
  }                                                                                                                                                                                                                                           
                                                                                                                                                                                                                                              
> citation("visreg")                                                                                                                                                                                                                          
Error in citation("visreg") : there is no package called 'visreg'                                                                                                                                                                             
> citation("ggplot2")                                                                                                                                                                                                                         
                                                                                                                                                                                                                                              
To cite ggplot2 in publications, please use                                                                                                                                                                                                   
                                                                                                                                                                                                                                              
  H. Wickham. ggplot2: Elegant Graphics for Data Analysis.                                                                                                                                                                                    
  Springer-Verlag New York, 2016.                                                                                                                                                                                                             
                                                                                                                                                                                                                                              
A BibTeX entry for LaTeX users is                                                                                                                                                                                                             
                                                                                                                                                                                                                                              
  @Book{,                                                                                                                                                                                                                                     
    author = {Hadley Wickham},                                                                                                                                                                                                                
    title = {ggplot2: Elegant Graphics for Data Analysis},                                                                                                                                                                                    
    publisher = {Springer-Verlag New York},                                                                                                                                                                                                   
    year = {2016},                                                                                                                                                                                                                            
    isbn = {978-3-319-24277-4},                                                                                                                                                                                                               
    url = {https://ggplot2.tidyverse.org},                                                                                                                                                                                                    
  }                                                                                                                                                                                                                                           
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                              
> citation("bbmle")                                                                                                                                                                                                                           
                                                                                                                                                                                                                                              
To cite package 'bbmle' in publications use:                                                                                                                                                                                                  
                                                                                                                                                                                                                                              
  Bolker B, R Development Core Team (2022). _bbmle: Tools for General                                                                                                                                                                         
  Maximum Likelihood Estimation_. R package version 1.0.25,                                                                                                                                                                                   
  <https://CRAN.R-project.org/package=bbmle>.                                                                                                                                                                                                 
                                                                                                                                                                                                                                              
A BibTeX entry for LaTeX users is                                                                                                                                                                                                             
                                                                                                                                                                                                                                              
  @Manual{,                                                                                                                                                                                                                                   
    title = {bbmle: Tools for General Maximum Likelihood Estimation},                                                                                                                                                                         
    author = {Ben Bolker and {R Development Core Team}},                                                                                                                                                                                      
    year = {2022},                                                                                                                                                                                                                            
    note = {R package version 1.0.25},                                                                                                                                                                                                        
    url = {https://CRAN.R-project.org/package=bbmle},                                                                                                                                                                                         
  }                                                                                                                                                                                                                                           
                                                                                                                                                                                                                                              
      
### Check at which months there were significant relations for each fraction
###  at each side.



pvalues.gam.forest <- matrix(nrow=12, ncol=4)
rownames(pvalues.gam.forest) <- months
colnames(pvalues.gam.forest) <- c("Leaves", "FWD", "Misc", "Total")

pvalues.gam.burnt <- matrix(nrow=12, ncol=4)
rownames(pvalues.gam.burnt) <- months
colnames(pvalues.gam.burnt) <- c("Leaves", "FWD", "Misc", "Total")

for(i in 1:length(months)) {
	for(j in 1:4) {
		results.i.j.forest <- results.gam.forest[[i]][[j]]
		p.i.j.forest <- summary(results.i.j.forest$gam)$s.table[,4]
		pvalues.gam.forest[i,j] <- p.i.j.forest
		results.i.j.burnt <- results.gam.burnt[[i]][[j]]
		p.i.j.burnt <- summary(results.i.j.burnt$gam)$s.table[,4]
		pvalues.gam.burnt[i,j] <- p.i.j.burnt
	}
}


formatC(pvalues.gam.forest,3,format="g")
formatC(pvalues.gam.burnt,3, format="g")

# How many significant?

apply(pvalues.gam.forest, 2, function(x) sum(x<0.05))
apply(pvalues.gam.burnt, 2, function(x) sum(x<0.05))


# Figure for montly patterns for leaves in forest


newdistances.forest <- 0:150
newdistances.burnt <- -150:0

setwd(dir.main)

png(filename="Figure4_GAMM_signif_leaves.png", height=25, width=25,
	res=300, unit="cm")
par(mfrow=c(4,3), mar=c(3,3,2,2), oma=c(4,4,2,2))

il <- 1 #index for letters

months.long <- c("January", "February", "March", "April", "May", "June",
	"July", "August", "September", "October", "November", "December")

for(i in 1:12) {
	print(paste(i, "of 12"))
	# Separate the data
	dados.i <- subset(dadost, Month==months[[i]])
	dados.i.forest <- subset(dados.i, Side_edgeAsForest == "Forest")

	# Models for forest

	mod.leaves.forest.foo <- gam(Mass_leaves ~ 1, data=dados.i.forest, 
		family=tw)
	fam <- summary(mod.leaves.forest.foo)$family$family
	fam <- strsplit(fam, "=")[[1]][2]
	fam <- gsub("\\)","",fam)
	fam <- as.numeric(fam)
	gamm.i <- gamm4(Mass_leaves ~ s(Distance, fx=F, k=4), 
		random=~(1|Transect), data=dados.i.forest, family=Tweedie(fam))

	predicted.leaves.forest <- predict(gamm.i$gam,
		type="response", newdata=list(Distance=newdistances.forest))
	p.leaves.forest <- summary(gamm.i$gam)$s.table[,4]

	plot(Mass_leaves ~ Distance, data=dados.i.forest, pch=16, 
		type="n",
		xlab="", ylab="", 
		main=paste(letters[il],". ", months.long[i], sep="")
		)
	#abline(v=0, col="gray")
	points(Mass_leaves ~ Distance, data=dados.i.forest, pch=16, 
		bg=alpha("black", 0.5))
	lines(predicted.leaves.forest ~ newdistances.forest, type="l",
		lty=ifelse(p.leaves.forest < 0.05, 1, 2), lwd=1.5)
	text(x=150, y=max(dados.i.forest$Mass_leaves)-
		max(dados.i.forest$Mass_leaves)/10,
		labels=paste("p=", formatC(p.leaves.forest, 
		digits=2, format="g"), sep=""), pos=2)
	il <- il+1

	mtext(side=1, text="Distance from edge (m)", outer=T)
	mtext(side=2, text="Monthly leaf litter production (t/ha)", outer=T)



}

dev.off()


### Another version - standardizing the y scale


png(filename="Figure4_GAMM_signif_leaves_v2.png", height=25, width=25,
	res=300, unit="cm")
par(mfrow=c(4,3), mar=c(3,3,2,2), oma=c(4,4,2,2))

il <- 1 #index for letters

months.long <- c("January", "February", "March", "April", "May", "June",
	"July", "August", "September", "October", "November", "December")

for(i in 1:12) {
	print(paste(i, "of 12"))
	# Separate the data
	dados.i <- subset(dadost, Month==months[[i]])
	dados.i.forest <- subset(dados.i, Side_edgeAsForest == "Forest")

	# Models for forest

	mod.leaves.forest.foo <- gam(Mass_leaves ~ 1, data=dados.i.forest, 
		family=tw)
	fam <- summary(mod.leaves.forest.foo)$family$family
	fam <- strsplit(fam, "=")[[1]][2]
	fam <- gsub("\\)","",fam)
	fam <- as.numeric(fam)
	gamm.i <- gamm4(Mass_leaves ~ s(Distance, fx=F, k=4), 
		random=~(1|Transect), data=dados.i.forest, family=Tweedie(fam))

	predicted.leaves.forest <- predict(gamm.i$gam,
		type="response", newdata=list(Distance=newdistances.forest))
	p.leaves.forest <- summary(gamm.i$gam)$s.table[,4]

	plot(Mass_leaves ~ Distance, data=dados.i.forest, pch=16, 
		type="n", 
		ylim=range(subset(dadost, Side_edgeAsForest == "Forest")$Mass_leaves),
		xlab="", ylab="", 
		main=paste(letters[il],". ", months.long[i], sep="")
		)
	#abline(v=0, col="gray")
	points(Mass_leaves ~ Distance, data=dados.i.forest, pch=16, 
		bg=alpha("black", 0.5))
	lines(predicted.leaves.forest ~ newdistances.forest, type="l",
		lty=ifelse(p.leaves.forest < 0.05, 1, 2), lwd=1.5)
	text(x=150, y=max(dados.i.forest$Mass_leaves)-
		max(dados.i.forest$Mass_leaves)/10,
		labels=paste("p=", formatC(p.leaves.forest, 
		digits=2, format="g"), sep=""), pos=2)
	il <- il+1

	mtext(side=1, text="Distance from edge (m)", outer=T)
	mtext(side=2, text="Monthly leaf litter production (t/ha)", outer=T)



}

dev.off()


