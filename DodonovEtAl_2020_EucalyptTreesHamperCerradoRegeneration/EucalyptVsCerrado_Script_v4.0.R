#############################################################
### This is the script for the analyses used in           ###
### "Remaining eucalypt trees may hamper woody plant      ###
###   regeneration in a Neotropical savanna"              ###
### by Pavel Dodonov, Andreza L. Braga,                   ###
###  Maria José Dias Sales, Rafael de Oliveira Xavier     ###
###           and Dalva M. Silva-Matos                    ###
###                                                       ###
### Script by Pavel Dodonov                               ###
### No rights reserved except for those reserved by the   ###
###  publisher.                                           ###
### Contact: pdodonov@gmail.com                           ###
#############################################################

library(mgcv)
library(car)
library(vegan)

setwd("/home/pavel/Profissional/Pesquisa/MyPapers-Working/EucalyptVsCerrado")

# First step - check if summary statistics are affected by eucalypt

data.univar <- read.table("dadosALB_univar.txt", header=T, sep="\t")
data.multivar <- read.table("dadosALB_multivar.txt", header=T, sep="\t", row.names=1)

# Response variables: Nind, S, Schvin, Micalb, Micleg, Piprot, Strobo, BAtot, Shannon, Hmean.

# Explanatory variables: eucalypts in the plot, all eucalypts, invasive grasses.


# Distribution of the explanatory variables
par(mfrow=c(2,2))
plot(data.univar$Neuc_tot)
plot(data.univar$Aeuc_tot)
plot(data.univar$Cgram)

# There's one point with a very large number of eucalypts...

# Solution: remove this plot!
plot.remove <- which(data.univar$Neuc_tot == max(data.univar$Neuc_tot))
data.univar2 <- data.univar[-plot.remove,]

par(mfrow=c(2,2))
plot(data.univar2$Neuc_tot)
plot(data.univar2$Aeuc_tot)
plot(data.univar2$Cgram)

# The distribution seems OK now.

names(data.univar)

### Explanatory variables: Aeuc_tot, logNeuc_tot
### Vegetation response variables: Nind, BAtot, Hmean, S, Shannon; Cgram;
### Population response variables: Schvin, Micalb, Micleg, Piprot, Strobo

### Exploratory plots - Aeuc_tot

data.univar2$Cgram_sqrt <- sqrt(asin(data.univar2$Cgram/100))

data.univar2$Cmelinis_sqrt <- sqrt(asin(data.univar2$Cmelinis/100))

### Calculate mean basal area per plot

data.univar2$BAmean <- data.univar2$BAtot / data.univar2$Nind


var.names <- c("Nind", "BAtot", "BAmean", "Hmean", "S", "Shannon",  "Cgram_sqrt", "Schvin", "Micalb", "Micleg", "Piprot", "Strobo")
var.names.long <- c("Total abundance", "Total basal area", "Mean basal area", "Mean height", "Richness", "Diversity", "Invasive grasses", "Schefflera vinosa", "Miconia albicans", "Miconia ligustroides", "Piptocarpha rotundifolia", "Stryphnodendron obovatum")

setwd("/home/pavel/Profissional/Pesquisa/MyPapers-Working/EucalyptVsCerrado/resubmitted4")

png(filename="eucalypts_expl_Neuc_tot_2.png", height=20, width=20, res=300, unit="cm")

par(mfrow=c(4,3), mar=c(3,3,2,2), oma=c(3,3,2,2))
for(i in 1:12) {
	y.temp <- data.univar2[,var.names[i]]
	plot(y.temp ~ data.univar2[,"Neuc_tot"], xlab="", ylab="", main=var.names.long[i])
}
mtext(side=1, text="Number of eucalypts", outer=T)

dev.off()


png(filename="eucalypts_expl_Aeuc_tot_2.png", height=20, width=20, res=300, unit="cm")

par(mfrow=c(4,3), mar=c(3,3,2,2), oma=c(3,3,2,2))
for(i in 1:12) {
	y.temp <- data.univar2[,var.names[i]]
	plot(y.temp ~ data.univar2[,"Aeuc_tot"], xlab="", ylab="", main=var.names.long[i])
}
mtext(side=1, text="Basal area of eucalypts", outer=T)

dev.off()


### Adjusting GLMs - for number of eucalypts and their basal area as explanatory variables


result.Neuc_tot <- matrix(nrow=12, ncol=3)
row.names(result.Neuc_tot) <- var.names
colnames(result.Neuc_tot) <- c("EDF", "Deviance_explained", "p-value")

Neuc_tot.new <- seq(min(data.univar2$Neuc_tot), max(data.univar2$Neuc_tot), length.out=100)

predict.Neuc_tot <- matrix(nrow=100, ncol=12)
colnames(predict.Neuc_tot) <- var.names


result.Aeuc_tot <- matrix(nrow=12, ncol=3)
row.names(result.Aeuc_tot) <- var.names
colnames(result.Aeuc_tot) <- c("EDF", "Deviance_explained", "p-value")

Aeuc_tot.new <- seq(min(data.univar$Aeuc_tot), max(data.univar$Aeuc_tot), length.out=100)

predict.Aeuc_tot <- matrix(nrow=100, ncol=12)
colnames(predict.Aeuc_tot) <- var.names


for(i in 1:12) {
	y.temp <- data.univar2[,var.names[i]]
	gam.temp <- gam(y.temp ~ s(Aeuc_tot, fx=F, k=3), data=data.univar2, family=ifelse(is.integer(y.temp),quasipoisson,ifelse(any(y.temp==0),gaussian,Gamma)))
	edf.temp <- summary(gam.temp)$s.table[1,1]
	dev.temp <- summary(gam.temp)$dev.expl
	sig.temp <- summary(gam.temp)$s.table[1,4]
	result.Aeuc_tot[i,] <- c(edf.temp, dev.temp, sig.temp)
	predict.Aeuc_tot[,i] <- predict(gam.temp, newdata=list(Aeuc_tot=Aeuc_tot.new), type="response")
}




for(i in 1:12) {
	y.temp <- data.univar2[,var.names[i]]
	gam.temp <- gam(y.temp ~ s(Neuc_tot, fx=F, k=3), data=data.univar2, family=ifelse(is.integer(y.temp),quasipoisson,ifelse(any(y.temp==0),gaussian,Gamma)))
	edf.temp <- summary(gam.temp)$s.table[1,1]
	dev.temp <- summary(gam.temp)$dev.expl
	sig.temp <- summary(gam.temp)$s.table[1,4]
	result.Neuc_tot[i,] <- c(edf.temp, dev.temp, sig.temp)
	predict.Neuc_tot[,i] <- predict(gam.temp, newdata=list(Neuc_tot=Neuc_tot.new), type="response")
}


result.Aeuc_tot
result.Neuc_tot


### Number of eucalypts has consistently greater explanation

y.names.long <- c("Number of individuals", "Area (cm^2)", "Area (cm^2)", "Height (cm)", "Number of species", "Shannon index", "Cover (%)", "Number of individuals", "Number of individuals", "Number of individuals", "Number of individuals", "Number of individuals")

png(filename="eucalypts_figure_result_final.png", height=20, width=20, res=300, unit="cm")

par(mfrow=c(4,3), mar=c(3,4,2,2), oma=c(3,3,2,2))
for(i in 1:12) {
	y.temp <- data.univar2[,var.names[i]]
	if(i < 7) plot(y.temp ~ data.univar2[,"Neuc_tot"], xlab="", ylab="", main=paste(letters[i],") ",var.names.long[i], sep=""), xaxt="s")
	if(i == 7) {
		plot(y.temp ~ data.univar2[,"Neuc_tot"], xlab="", ylab="", main=paste(letters[i],") ",var.names.long[i], sep=""), xaxt="s", yaxt="n")
		axis(side=2, at=sqrt(asin(c(0, 10, 50, 100)/100)), labels=c(0, 10, 50, 100))
	}
	if(i == 8) plot(y.temp ~ data.univar2[,"Neuc_tot"], xlab="", ylab="", main=expression(paste(bold("h)"), " ", bolditalic("Schefflera vinosa"))), xaxt="s")
	if(i == 9) plot(y.temp ~ data.univar2[,"Neuc_tot"], xlab="", ylab="", main=expression(paste(bold("i)"), " ", bolditalic("Miconia albicans"))), xaxt="s")
	if(i == 10) plot(y.temp ~ data.univar2[,"Neuc_tot"], xlab="", ylab="", main=expression(paste(bold("j)"), " ", bolditalic("Miconia ligustroides"))), xaxt="s")
	if(i == 11) plot(y.temp ~ data.univar2[,"Neuc_tot"], xlab="", ylab="", main=expression(paste(bold("k)"), " ", bolditalic("Piptocarpha rotundifolia"))), xaxt="s")
	if(i == 12) plot(y.temp ~ data.univar2[,"Neuc_tot"], xlab="", ylab="", main=expression(paste(bold("l)"), " ", bolditalic("Stryphnodendron obovatum"))), xaxt="s")
	mtext(text=ifelse(i!=2 & i!=3,y.names.long[i],expression(paste("Area (",cm^2,")"))), side=2, line=2, cex=0.8)
	line.col <- ifelse(result.Neuc_tot[i,"p-value"] < 0.05, "black", "gray60")
	line.type <- ifelse(result.Neuc_tot[i,"p-value"] < 0.10, 1, 2)
	lines(predict.Neuc_tot[,i] ~ Neuc_tot.new, col=line.col, lty=line.type, lwd=2)
	text(x=max(data.univar2[,"Neuc_tot"])-2, y=max(y.temp)-diff(range(y.temp))/10, labels=paste(round(result.Neuc_tot[i,"Deviance_explained"]*100),"% dev. expl."), cex=1)
	#axis(side=1, at=log(c(0,1,3,7,14,25)+1), lab=c(0,1,3,7,14,25))
}
mtext(side=1, text="Number of eucalypts", outer=T)

dev.off()


### CCA
# Remove species that occur in less than five plots
Nplots <- apply(data.multivar, 2, function(x) sum(x>0) )
data.cca <- data.multivar[,Nplots>3]
data.cca2 <- data.cca[data.univar$Neuc_tot<26,]
ncol(data.cca2)

data.cca2.presabs <- data.cca2
data.cca2.presabs[data.cca2>1] <- 1

### CCA for presence/absence data

cca.abund <- cca(data.cca2 ~ Neuc_tot + Aeuc_tot, data=data.univar2)
cca.abund.anova <- anova(cca.abund, permutations=how(nperm=9999))
cca.abund.anova

cca.presabs <- cca(data.cca2.presabs ~ Neuc_tot + Aeuc_tot, data=data.univar2)
cca.presabs.anova <- anova(cca.presabs, permutations=how(nperm=9999))
cca.presabs.anova

plot(cca.abund, display=c("species", "bp"))
plot(cca.presabs, display=c("species", "bp"))


scores.abund <- scores(cca.abund, display=c("species","bp"))
scores.presabs <- scores(cca.presabs, display=c("species","bp"))

png(filename="resubmitted3/EucalyptsVsCerrado_SM4.png", height=15, width=30, unit="cm", res=300)
par(mfrow=c(1,2), mar=c(4,4,2,2), oma=c(1,1,1,1))
plot(scores.abund$species[,2] ~ scores.abund$species[,1], xlab="CCA 1", ylab="CCA 2", main="Abundance", type="n", xlim=c(-0.9, 1.1))
text(labels=attr(scores.abund$species, "dimnames")[[1]], x=scores.abund$species[,1], y=scores.abund$species[,2], adj=0.7, cex=0.8)
segments(x0=0, y0=0, x1=scores.abund$biplot[,1], y1=scores.abund$biplot[,2], col="darkred")
text(labels=c("Eucalypt\nnumber", "Eucalypt\nbasal area"), y=scores.abund$biplot[,2], x=scores.abund$biplot[,1], cex=0.8, col="darkred", adj=0.8)

plot(scores.presabs$species[,2] ~ scores.presabs$species[,1], xlab="CCA 1", ylab="CCA 2", main="Composition", type="n", xlim=c(-0.9, 1.0), ylim=c(-0.7,0.4))
text(labels=attr(scores.presabs$species, "dimnames")[[1]], x=scores.presabs$species[,1], y=scores.presabs$species[,2], adj=0.5, cex=0.8)
segments(x0=0, y0=0, x1=scores.presabs$biplot[,1], y1=scores.presabs$biplot[,2], col="darkred")
text(labels=c("Eucalypt\nnumber", "Eucalypt\nbasal area"), y=scores.presabs$biplot[,2], x=scores.presabs$biplot[,1], cex=0.8, col="darkred", adj=0.8)
dev.off()


