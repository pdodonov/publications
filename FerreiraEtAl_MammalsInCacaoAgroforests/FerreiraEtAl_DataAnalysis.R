### Data analysis and visualization for the manuscript
### "Mammals in cacao agroforests: implications of management intensification
### in two contrasting landscapes", by Aluane Silva Ferreiraa, 
### Carlos A. Peresb, Pavel Dodonov, Eduardo Mariano-Neto, Deborah Faria
### and Camila Righetto Cassanoa

# 0 - Adjust settings to better show model results:

options("width"=120)

# 1 - Load the required packages
library(bbmle) # for model selection
library(car) # for VIF calculation
library(lme4) # for model selection
library(MASS) # for model selection
library(gstat) # for spatial analyses
library(sp) # for spatial analyses
library(nlme) # for mixed-effects models
library(MuMIn) # for multi-model inference
library(yarrr) # for some data visualization
library(scales) # for transparency
library(vegan) # for multivariate analyses

# 2 - Load the data and perform some data exploration

setwd("e:/Pavel/Profissional/Pesquisa/MyPapers_Working/Aluane_2022_Cabrucas/")

# Don't forget to change the above line to your own working directory.

# Univariate data

dados_cac <- read.csv2("FerreiraEtAl_data.csv", stringsAsFactors=T)
str(dados_cac)


# Data exploration and collinearity

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "purple", ...)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))^2
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}


pairs(dados_cac[,-c(1:32)], upper.panel=panel.smooth, lower.panel=panel.cor, 
  diag.panel=panel.hist)

pairs(dados_cac[,-c(1:32,34,37)], upper.panel=panel.smooth, 
  lower.panel=panel.cor, diag.panel=panel.hist)

# There appear to be some outliers for two variables.
# Transforming these variables into square-root and log:

par(mfrow=c(3,3), mar=c(4,4,2,2))
plot(richness_40 ~ density_nat, ylab = "Richness", xlab="Native trees number", 
  col=dados_cac$l_, dados_cac)
plot(1:1,type="n",xlab="",ylab="", xaxt="n", yaxt="n", bty="n")
plot(1:1,type="n",xlab="",ylab="", xaxt="n", yaxt="n", bty="n")
plot(richness_40 ~ BA_frutexo, ylab = "Richness", xlab="Basal area (frut_exo)", 
  col=dados_cac$l_, dados_cac)
plot(richness_40 ~ sqrt(BA_frutexo), ylab = "Richness", 
  xlab="sqrt Basal area (frut_exo)", col=dados_cac$l_, dados_cac)
plot(richness_40 ~ log(BA_frutexo+1), ylab = "Richness", 
  xlab="log Basal area (frut_exo)", col=dados_cac$l_, dados_cac)
plot(richness_40 ~ dograte_30, ylab = "Richness", xlab="Dog rate", 
  col=dados_cac$l_, dados_cac)
plot(richness_40 ~ sqrt(dograte_30), ylab = "Richness", xlab="sqrt Dog rate", 
  col=dados_cac$l_, dados_cac)
plot(richness_40 ~ log10(dograte_30+1), ylab = "Richness", xlab="log Dog rate", 
  col=dados_cac$l_, dados_cac)


# log transformation works better for dagrate, sqrt for BA_frutexo

dados_cac$BA_frutexo_sqrt <- sqrt(dados_cac$BA_frutexo)
dados_cac$dograte_30_log <- log10(dados_cac$dograte_30+1)

# Testing collinearity, using GVIF/2Df due to the presence
#  of categorical variables

mod.teste <- glm(richness_40 ~ l_ + density_nat + BA_frutexo_sqrt + 
  dograte_30_log, data=dados_cac)
vif(mod.teste) 

# All VIF values are smaller than two - no collinearity

#########################################################################
### Analyses for total richness
#########################################################################

# Full model

mod.riq.tot <- glm(richness_40 ~ l_ + density_nat + BA_frutexo_sqrt + 
  dograte_30_log + l_:density_nat + l_:BA_frutexo_sqrt + l_:dograte_30_log, 
  data=dados_cac)

# The variables are as follows:
# l_ - landscape
# density_nat - density of native trees
# BA_frutexo_sqrt - basal area of fructiferous non-native (exotic) trees
#  (square-root transformed)
# dograte_30_log - records of dogs, log-transformed.

# Checking residual distribution

resid.mod.glm <- resid(mod.riq.tot)

hist(resid.mod.glm)
shapiro.test(resid.mod.glm)

# The residuals appear to be normally distributed.

pred.mod.glm <- predict(mod.riq.tot, type="response")
par(mfrow=c(3,2), mar=c(2,2,2,2))
plot(mod.riq.tot)
plot(resid.mod.glm, main="Residuals in order")
hist(resid.mod.glm)
dev.off()

# No serious deviations from model assumptions.

par(mfrow=c(2,2), mar = c(2,2,2,2))
plot(resid.mod.glm ~dados_cac$l_, main="Resid ~ X") #ver linearidade
plot(resid.mod.glm ~dados_cac$density_nat, main="Resid ~ X")
plot(resid.mod.glm ~sqrt(dados_cac$BA_frutexo), main="Resid ~ X")
plot(resid.mod.glm ~log10(dados_cac$dograte_30+1), main="Resid ~ X")
# No serious deviations from model assumptions.
dev.off()

# Check for autocorrelation 
resid.mod.glm = data.frame(resid.mod.glm, dados_cac$x, dados_cac$y)
names(resid.mod.glm) = c("resid", "x", "y")
coordinates(resid.mod.glm) = c("x","y")
bubble(resid.mod.glm, "resid", col=c("black", "red"))
vario.mod.glm = variogram(resid~1, resid.mod.glm)
plot(vario.mod.glm)
# There is evidence for autocorrelation.

# GLS (Generalized Least Squares) to correct for autocorrelation, 
#  with different autocorrelation structures.

mod.riq.tot.gls <- gls(richness_40 ~ l_ + density_nat + BA_frutexo_sqrt + 
  dograte_30_log + l_:density_nat + l_:BA_frutexo_sqrt + l_:dograte_30_log, 
  data=dados_cac, cor=corLin(form = ~x+y|l_), method="ML")
## false convergence

mod.riq.tot.gls <- gls(richness_40 ~ l_ + density_nat + BA_frutexo_sqrt + 
  dograte_30_log + l_:density_nat + l_:BA_frutexo_sqrt + l_:dograte_30_log, 
  data=dados_cac, cor=corSpher(form = ~x+y|l_), method="ML")
plot(nlme:::Variogram(mod.riq.tot.gls, form = ~x+y, resType = "normalized")) 

mod.riq.tot.gls <- gls(richness_40 ~ l_ + density_nat + BA_frutexo_sqrt + 
  dograte_30_log + l_:density_nat + l_:BA_frutexo_sqrt + l_:dograte_30_log, 
  data=dados_cac, cor=corExp(form = ~x+y|l_), method="ML")
plot(nlme:::Variogram(mod.riq.tot.gls, form = ~x+y, resType = "normalized")) 

mod.riq.tot.gls <- gls(richness_40 ~ l_ + density_nat + BA_frutexo_sqrt + 
  dograte_30_log + l_:density_nat + l_:BA_frutexo_sqrt + l_:dograte_30_log, 
  data=dados_cac, cor=corGaus(form = ~x+y|l_), method="ML")
plot(nlme:::Variogram(mod.riq.tot.gls, form = ~x+y, resType = "normalized")) ## Feioso

# To maintain consistency with the other analyses below, the corSpher
#  structure was selected.

mod.riq.tot.gls <- gls(richness_40 ~ l_ + density_nat + BA_frutexo_sqrt + 
  dograte_30_log + l_:density_nat + l_:BA_frutexo_sqrt + l_:dograte_30_log, 
  data=dados_cac, cor=corSpher(form = ~x+y|l_), method="ML")

# Which variables are important?

mod.riq.tot.gls.dredge <- dredge(mod.riq.tot.gls)

mod.riq.tot.gls.dredge

### Best model: dog_rate_30_log. 

bestModel.riq.tot <- gls(richness_40 ~ dograte_30_log, data=dados_cac, 
  cor=corSpher(form = ~x+y|l_), method="ML")

summary(bestModel.riq.tot)

plot(richness_40 ~ dograte_30_log, data=dados_cac)
abline(bestModel.riq.tot)


# The number of species decreases as dog records increase.


# As the analyses below follow the same steps, they will not be explained
#  in detail.

#########################################################################
# Analyses for richness of sensitive species
#########################################################################

mod.riq.sens <- glm(richness_sens40~ l_ + density_nat + BA_frutexo_sqrt + 
  dograte_30_log + l_:density_nat + l_:BA_frutexo_sqrt + 
  l_:dograte_30_log,data=dados_cac)


summary(mod.riq.sens)
resid.mod.glm.s <- resid(mod.riq.sens, type="deviance")
pred.mod.glm.s <- predict(mod.riq.sens, type="response")
par(mfrow=c(3,2), mar = c(2,2,2,2))
plot(mod.riq.sens) 
plot(resid.mod.glm.s, main="Deviance residuals ~ Predicted")
hist(resid.mod.glm.s)
dev.off()
shapiro.test(resid.mod.glm.s)
par(mfrow=c(2,2), mar = c(2,2,2,2))
plot(resid.mod.glm.s ~dados_cac$l_, main="Resid ~ X") #ver linearidade
plot(resid.mod.glm.s ~dados_cac$density_nat, main="Resid ~ X")
plot(resid.mod.glm.s ~log(dados_cac$BA_frutexo_sqrt), main="Resid ~ X")
plot(resid.mod.glm.s ~dados_cac$dograte_30_log, main="Resid ~ X")
dev.off()
### No serious deviations from the model assumptions

# Checking for autocorrelation
resid.mod.glm.s = data.frame(resid.mod.glm.s, dados_cac$x, dados_cac$y)
str(resid.mod.glm.s)
names(resid.mod.glm.s) = c("resid", "x", "y")
coordinates(resid.mod.glm.s) = c("x","y")
bubble(resid.mod.glm.s, "resid", col=c("black", "red"))
vario.mod.gls.s.var = variogram(resid~1, resid.mod.glm.s)
plot(vario.mod.gls.s.var)

# There is evidence of autocorrelation.

mod.riq.sens.gls.autocor <- gls(richness_sens40~ l_ + density_nat + 
  BA_frutexo_sqrt + dograte_30_log +  l_:density_nat + l_:BA_frutexo_sqrt + 
  l_:dograte_30_log, 
  data=dados_cac, cor=corLin(form = ~x+y|l_), method="ML") 
# False convergence

mod.riq.sens.gls.autocor <- gls(richness_sens40~ l_ + density_nat + 
  BA_frutexo_sqrt + dograte_30_log +  l_:density_nat + l_:BA_frutexo_sqrt + 
  l_:dograte_30_log, 
  data=dados_cac, cor=corExp(form = ~x+y|l_), method="ML") 
plot(nlme:::Variogram(mod.riq.sens.gls.autocor, form = ~x+y, 
  resType = "normalized"))

mod.riq.sens.gls.autocor <- gls(richness_sens40~ l_ + density_nat + 
  BA_frutexo_sqrt + dograte_30_log +  l_:density_nat + l_:BA_frutexo_sqrt + 
  l_:dograte_30_log, 
  data=dados_cac, cor=corSpher(form = ~x+y|l_), method="ML")
plot(nlme:::Variogram(mod.riq.sens.gls.autocor, form = ~x+y, 
  resType = "normalized"))

mod.riq.sens.gls.autocor <- gls(richness_sens40~ l_ + density_nat + 
  BA_frutexo_sqrt + dograte_30_log +  l_:density_nat + l_:BA_frutexo_sqrt + 
  l_:dograte_30_log, 
  data=dados_cac, cor=corGaus(form = ~x+y|l_), method="ML")
plot(nlme:::Variogram(mod.riq.sens.gls.autocor, form = ~x+y, 
  resType = "normalized"))

### For consistency, corSpher was chosen.

mod.riq.sens.gls.autocor <- gls(richness_sens40~ l_ + density_nat + 
  BA_frutexo_sqrt + dograte_30_log +  l_:density_nat + l_:BA_frutexo_sqrt + 
  l_:dograte_30_log, data=dados_cac, cor=corSpher(form = ~x+y|l_), method="ML")

### Which variables are important?

mod.riq.sens.gls.autocor.dredge <- dredge(mod.riq.sens.gls.autocor)
mod.riq.sens.gls.autocor.dredge

### MThe best model includes only the difference between the two landscapes.

bestModel.riq.sens <- gls(richness_sens40~ l_, 
  data=dados_cac, cor=corSpher(form = ~x+y|l_), method="ML")

plot(richness_sens40 ~ l_, data=dados_cac, range=0)
points(richness_sens40 ~ jitter(as.numeric(l_)), data=dados_cac, pch=21,
	bg="white")



#########################################################################
# Analysis for richness of insensitive species
#########################################################################

mod.riq.ns <- glm(richness_ns40~ l_ + density_nat + BA_frutexo_sqrt + 
	dograte_30_log +  l_:density_nat + l_:BA_frutexo_sqrt + l_:dograte_30_log, 
	data=dados_cac)

summary(mod.riq.ns)
resid.mod.gls.ns <- resid(mod.riq.ns, type="deviance")
pred.mod.glm.ns <- predict(mod.riq.ns, type="response") #response p ele ajustar ao preditor linear
par(mfrow=c(3,2), mar = c(2,2,2,2))
plot(mod.riq.ns) 
plot(resid.mod.gls.ns, main="Deviance residuals ~ Predicted")
hist(resid.mod.gls.ns)
shapiro.test(resid.mod.gls.ns)
par(mfrow=c(2,2), mar = c(2,2,2,2))
plot(resid.mod.gls.ns ~dados_cac$l_, main="Resid ~ X") #ver linearidade
plot(resid.mod.gls.ns ~dados_cac$density_nat, main="Resid ~ X")
plot(resid.mod.gls.ns ~log(dados_cac$BA_frutexo_sqrt), main="Resid ~ X")
plot(resid.mod.gls.ns ~dados_cac$dograte_30_log, main="Resid ~ X")
dev.off()
# No serious deviations from model assumptions.

resid.mod.gls.ns = data.frame(resid.mod.gls.ns, dados_cac$x, dados_cac$y)
str(resid.mod.gls.ns)
names(resid.mod.gls.ns) = c("resid", "x", "y")
coordinates(resid.mod.gls.ns) = c("x","y")
bubble(resid.mod.gls.ns, "resid", col=c("black", "red"))
vario.mod.gls.s.var = variogram(resid~1, resid.mod.gls.ns)
plot(vario.mod.gls.s.var, pch=16, cex=1.2)
# There is evidence of autocorrelation.

mod.riq.ns.gls.autocor <- gls(richness_ns40~ l_ + density_nat + 
	BA_frutexo_sqrt + dograte_30_log +  l_:density_nat + l_:BA_frutexo_sqrt + 
	l_:dograte_30_log, data=dados_cac, cor=corLin(form = ~x+y|l_))
plot(nlme:::Variogram(mod.riq.ns.gls.autocor, form = ~x+y, 
	resType = "normalized")) 

mod.riq.ns.gls.autocor <- gls(richness_ns40~ l_ + density_nat + 
	BA_frutexo_sqrt + dograte_30_log +  l_:density_nat + l_:BA_frutexo_sqrt + 
	l_:dograte_30_log, data=dados_cac, cor=corExp(form = ~x+y|l_))
plot(nlme:::Variogram(mod.riq.ns.gls.autocor, form = ~x+y, 
	resType = "normalized")) 

mod.riq.ns.gls.autocor <- gls(richness_ns40~ l_ + density_nat + 
	BA_frutexo_sqrt + dograte_30_log +  l_:density_nat + l_:BA_frutexo_sqrt + 
	l_:dograte_30_log, data=dados_cac, cor=corSpher(form = ~x+y|l_))
plot(nlme:::Variogram(mod.riq.ns.gls.autocor, form = ~x+y, 
	resType = "normalized")) 

mod.riq.ns.gls.autocor <- gls(richness_ns40~ l_ + density_nat + 
	BA_frutexo_sqrt + dograte_30_log +  l_:density_nat + l_:BA_frutexo_sqrt + 
	l_:dograte_30_log, data=dados_cac, cor=corGaus(form = ~x+y|l_))
plot(nlme:::Variogram(mod.riq.ns.gls.autocor, form = ~x+y, 
	resType = "normalized"))

mod.riq.ns.gls.autocor <- gls(richness_ns40~ l_ + density_nat + 
	BA_frutexo_sqrt + dograte_30_log +  l_:density_nat + l_:BA_frutexo_sqrt + 
	l_:dograte_30_log, data=dados_cac, cor=corExp(form = ~x+y|l_))
plot(nlme:::Variogram(mod.riq.ns.gls.autocor, form = ~x+y, 
	resType = "normalized")) 

# For consistency, corSpher was chosen.

mod.riq.ns.gls.autocor <- gls(richness_ns40~ l_ + density_nat + 
	BA_frutexo_sqrt + dograte_30_log +  l_:density_nat + l_:BA_frutexo_sqrt + 
	l_:dograte_30_log, data=dados_cac, cor=corSpher(form = ~x+y|l_), 
	method="ML") 

# Which variables are important?

mod.riq.ns.gls.autocor.dredge <- dredge(mod.riq.ns.gls.autocor)
mod.riq.ns.gls.autocor.dredge

# The null model has a dAICc < 2 and was chosen.

bestModel.riq.ns <- gls(richness_ns40~ 1, 
	data=dados_cac, cor=corSpher(form = ~x+y|l_), 
	method="ML") 



#########################################################################
# Analyses for abundance
#########################################################################

# Correct abundance by the sampling effort

par(mfrow=c(2,2), mar=c(4,4,2,2))
plot(abundance ~ effort, data=dados_cac)
plot(ab_sens ~ effort, data=dados_cac)
plot(ab_notsens ~ effort, data=dados_cac)

# Standardize per sampling effort (divide)
dados_cac$abundance_st <- dados_cac$abundance / dados_cac$effort
dados_cac$ab_sens_st <- dados_cac$abundance / dados_cac$effort
dados_cac$ab_notsens_st <- dados_cac$ab_notsens / dados_cac$effort

par(mfrow=c(2,2), mar=c(4,4,2,2))
plot(abundance_st ~ effort, data=dados_cac)
plot(ab_sens_st ~ effort, data=dados_cac)
plot(ab_notsens_st ~ effort, data=dados_cac)


#########################################################################
# Total abundance
#########################################################################

mod.ab.tot <- glm(abundance_st ~ l_ + density_nat + BA_frutexo_sqrt + 
	dograte_30_log +  l_:density_nat + l_:BA_frutexo_sqrt + l_:dograte_30_log, 
	data=dados_cac)

resid.abundance <- resid(mod.ab.tot, type="deviance")
pred.glm.abundance <- predict(mod.ab.tot, type="response") 
par(mfrow=c(5,2), mar = c(4,4,2,2))
plot(mod.ab.tot) 
plot(resid.abundance, main="Deviance residuals ~ Predicted")
hist(resid.abundance, breaks=20)
plot(resid.abundance ~dados_cac$l_, main="Resid ~ X") #ver linearidade
plot(resid.abundance ~dados_cac$density_nat, main="Resid ~ X")
plot(resid.abundance ~log(dados_cac$BA_frutexo), main="Resid ~ X")
plot(resid.abundance ~dados_cac$dograte_30, main="Resid ~ X")
dev.off()
shapiro.test(resid.abundance)

# No serious deviations from model assumptions.
# There seems to be some heteroscedasticity but not too serious.

resid.abundance = data.frame(resid.abundance, dados_cac$x, dados_cac$y)
str(resid.abundance)
names(resid.abundance) = c("resid", "x", "y")
coordinates(resid.abundance) = c("x","y")
bubble(resid.abundance, "resid", col=c("black", "red"))
vario.abundance = variogram(resid~1, resid.abundance)
plot(vario.abundance, pch=16, cex=1.2)

# There is evidence of autocorrelation.

mod.ab.tot.autocor <- gls(abundance_st ~ l_ + density_nat + 
	BA_frutexo_sqrt + dograte_30_log +  l_:density_nat + l_:BA_frutexo_sqrt + 
	l_:dograte_30_log, data=dados_cac, cor=corLin(form = ~x+y|l_), method="ML")
plot(nlme:::Variogram(mod.ab.tot.autocor, form = ~x+y, 
	resType = "normalized")) 

mod.ab.tot.autocor <- gls(abundance_st ~ l_ + density_nat + 
	BA_frutexo_sqrt + dograte_30_log +  l_:density_nat + l_:BA_frutexo_sqrt + 
	l_:dograte_30_log, data=dados_cac, cor=corSpher(form = ~x+y|l_), method="ML")
plot(nlme:::Variogram(mod.ab.tot.autocor, form = ~x+y, 
	resType = "normalized")) 

mod.ab.tot.autocor <- gls(abundance_st ~ l_ + density_nat + 
	BA_frutexo_sqrt + dograte_30_log +  l_:density_nat + l_:BA_frutexo_sqrt + 
	l_:dograte_30_log, data=dados_cac, cor=corGaus(form = ~x+y|l_), method="ML")
plot(nlme:::Variogram(mod.ab.tot.autocor, form = ~x+y, 
	resType = "normalized")) 

mod.ab.tot.autocor <- gls(abundance_st ~ l_ + density_nat + 
	BA_frutexo_sqrt + dograte_30_log +  l_:density_nat + l_:BA_frutexo_sqrt + 
	l_:dograte_30_log, data=dados_cac, cor=corExp(form = ~x+y|l_), method="ML")
plot(nlme:::Variogram(mod.ab.tot.autocor, form = ~x+y, 
	resType = "normalized")) 

# The different models are similar, so corSpher was used for consistency.

mod.ab.tot.autocor <- gls(abundance_st ~ l_ + density_nat + 
	BA_frutexo_sqrt + dograte_30_log +  l_:density_nat + l_:BA_frutexo_sqrt + 
	l_:dograte_30_log, data=dados_cac, cor=corSpher(form = ~x+y|l_), method="ML")

# Which variables are important?

mod.ab.tot.autocor.dredge <- dredge(mod.ab.tot.autocor)
mod.ab.tot.autocor.dredge

# The best model has only differences between the landscapes.

bestModel.ab.tot <- gls(abundance_st ~ l_, 
	data=dados_cac, cor=corSpher(form = ~x+y|l_), method="ML")

#########################################################################
# Abundance of sensitive species
#########################################################################

mod.ab.sens <- glm(ab_sens_st ~ l_ + density_nat + BA_frutexo_sqrt + 
	dograte_30_log +  l_:density_nat + l_:BA_frutexo_sqrt + l_:dograte_30_log, 
	data=dados_cac)


resid.ab_sens <- resid(mod.ab.sens, type="deviance")
pred.ab_sens <- predict(mod.ab.sens, type="response") #response p ele ajustar ao preditor linear
par(mfrow=c(5,2), mar = c(4,4,2,2))
plot(mod.ab.sens) 
plot(resid.ab_sens, main="Deviance residuals ~ Predicted")
hist(resid.ab_sens, breaks=20)
plot(resid.ab_sens ~dados_cac$l_, main="Resid ~ X") #ver linearidade
plot(resid.ab_sens ~dados_cac$density_nat, main="Resid ~ X")
plot(resid.ab_sens ~log(dados_cac$BA_frutexo), main="Resid ~ X")
plot(resid.ab_sens ~dados_cac$dograte_30, main="Resid ~ X")
dev.off()
shapiro.test(resid.ab_sens)

# No serious deviations from model assumptions.

resid.ab_sens = data.frame(resid.ab_sens, dados_cac$x, dados_cac$y)
str(resid.ab_sens)
names(resid.ab_sens) = c("resid", "x", "y")
coordinates(resid.ab_sens) = c("x","y")
bubble(resid.ab_sens, "resid", col=c("black", "red"))
vario.ab_sens = variogram(resid~1, resid.ab_sens)
plot(vario.ab_sens, pch=16, cex=1.2)

# There is evidence of autocorrelation.

mod.ab.sens.autocor <- gls(ab_sens_st ~ l_ + density_nat + BA_frutexo_sqrt + 
	dograte_30_log +  l_:density_nat + l_:BA_frutexo_sqrt + l_:dograte_30_log, 
	data=dados_cac, cor=corLin(form = ~x+y|l_), method="ML")
plot(nlme:::Variogram(mod.ab.sens.autocor, form = ~x+y, 
	resType = "normalized")) 

mod.ab.sens.autocor <- gls(ab_sens_st ~ l_ + density_nat + BA_frutexo_sqrt + 
	dograte_30_log +  l_:density_nat + l_:BA_frutexo_sqrt + l_:dograte_30_log, 
	data=dados_cac, cor=corSpher(form = ~x+y|l_), method="ML") 
plot(nlme:::Variogram(mod.ab.sens.autocor, form = ~x+y, 
	resType = "normalized"))

mod.ab.sens.autocor <- gls(ab_sens_st ~ l_ + density_nat + BA_frutexo_sqrt + 
	dograte_30_log +  l_:density_nat + l_:BA_frutexo_sqrt + l_:dograte_30_log, 
	data=dados_cac, cor=corGaus(form = ~x+y|l_), method="ML")
plot(nlme:::Variogram(mod.ab.sens.autocor, form = ~x+y, 
	resType = "normalized"))

mod.ab.sens.autocor <- gls(ab_sens_st ~ l_ + density_nat + BA_frutexo_sqrt + 
	dograte_30_log +  l_:density_nat + l_:BA_frutexo_sqrt + l_:dograte_30_log, 
	data=dados_cac, cor=corExp(form = ~x+y|l_), method="ML") 
plot(nlme:::Variogram(mod.ab.sens.autocor, form = ~x+y, 
	resType = "normalized"))

# For consistency, corSpher was selected.

mod.ab.sens.autocor <- gls(ab_sens_st ~ l_ + density_nat + BA_frutexo_sqrt + 
	dograte_30_log +  l_:density_nat + l_:BA_frutexo_sqrt + l_:dograte_30_log, 
	data=dados_cac, cor=corSpher(form = ~x+y|l_), method="ML")

# Which variables are important?

mod.ab.sens.autocor.dredge <- dredge(mod.ab.sens.autocor)
mod.ab.sens.autocor.dredge

# The best model has only differences between the two landscapes.

bestModel.ab.sens <- gls(ab_sens_st ~ l_, 
	data=dados_cac, cor=corSpher(form = ~x+y|l_), method="ML")

#########################################################################
# Abundance of insensitive species
#########################################################################

mod.ab.ns <- glm(ab_notsens_st ~ l_ + density_nat + BA_frutexo_sqrt + 
	dograte_30_log + l_:density_nat + l_:BA_frutexo_sqrt + l_:dograte_30_log, 
	data=dados_cac)

resid.ab.ns <- resid(mod.ab.ns, type="deviance")
pred.ab.ns <- predict(mod.ab.ns, type="response") #response p ele ajustar ao preditor linear
par(mfrow=c(5,2), mar = c(4,4,2,2))
plot(mod.ab.ns) 
plot(resid.ab.ns, main="Deviance residuals ~ Predicted")
hist(resid.ab.ns, breaks=20)
plot(resid.ab.ns ~dados_cac$l_, main="Resid ~ X") #ver linearidade
plot(resid.ab.ns ~dados_cac$density_nat, main="Resid ~ X")
plot(resid.ab.ns ~log(dados_cac$BA_frutexo), main="Resid ~ X")
plot(resid.ab.ns ~dados_cac$dograte_30, main="Resid ~ X")
dev.off()
shapiro.test(resid.ab.ns)

# No serious deviations from model assumptions.

resid.ab.ns = data.frame(resid.ab.ns, dados_cac$x, dados_cac$y)
str(resid.ab.ns)
names(resid.ab.ns) = c("resid", "x", "y")
coordinates(resid.ab.ns) = c("x","y")
bubble(resid.ab.ns, "resid", col=c("black", "red"))
vario.ab.ns_sqrt.gls.s.var = variogram(resid~1, resid.ab.ns)
plot(vario.ab.ns_sqrt.gls.s.var, pch=16, cex=1.2)

# There is evidence of autocorrelation.

mod.ab.ns.autocor <- gls(ab_notsens_st ~ l_ + density_nat + BA_frutexo_sqrt + 
	dograte_30_log + l_:density_nat + l_:BA_frutexo_sqrt + l_:dograte_30_log, 
	data=dados_cac, cor=corLin(form = ~x+y|l_), method="ML")
# False convergence

mod.ab.ns.autocor <- gls(ab_notsens_st ~ l_ + density_nat + BA_frutexo_sqrt + 
	dograte_30_log + l_:density_nat + l_:BA_frutexo_sqrt + l_:dograte_30_log, 
	data=dados_cac, cor=corSpher(form = ~x+y|l_), method="ML")
plot(nlme:::Variogram(mod.ab.ns.autocor, form = ~x+y, 
	resType = "normalized"))

mod.ab.ns.autocor <- gls(ab_notsens_st ~ l_ + density_nat + BA_frutexo_sqrt + 
	dograte_30_log + l_:density_nat + l_:BA_frutexo_sqrt + l_:dograte_30_log, 
	data=dados_cac, cor=corGaus(form = ~x+y|l_), method="ML")
plot(nlme:::Variogram(mod.ab.ns.autocor, form = ~x+y, 
	resType = "normalized"))

mod.ab.ns.autocor <- gls(ab_notsens_st ~ l_ + density_nat + BA_frutexo_sqrt + 
	dograte_30_log + l_:density_nat + l_:BA_frutexo_sqrt + l_:dograte_30_log, 
	data=dados_cac, cor=corExp(form = ~x+y|l_), method="ML")
plot(nlme:::Variogram(mod.ab.ns.autocor, form = ~x+y, 
	resType = "normalized"))

# corSpher seems to work well enough.

mod.ab.ns.autocor <- gls(ab_notsens_st ~ l_ + density_nat + BA_frutexo_sqrt + 
	dograte_30_log + l_:density_nat + l_:BA_frutexo_sqrt + l_:dograte_30_log, 
	data=dados_cac, cor=corSpher(form = ~x+y|l_), method="ML")

# Which variables are important?

mod.ab.ns.autocor.dredge <- dredge(mod.ab.ns.autocor)
mod.ab.ns.autocor.dredge

# Best model contains only landscape.

bestModel.ab.ns <- gls(ab_notsens_st ~ l_, 
	data=dados_cac, cor=corSpher(form = ~x+y|l_), method="ML")

# Summarizing the best models:

summary(bestModel.riq.tot)
summary(bestModel.riq.sens)
summary(bestModel.riq.ns)
summary(bestModel.ab.tot)
summary(bestModel.ab.sens)
summary(bestModel.ab.ns)


# Total richness - effect of dog records
# Senstive richness - differences between the two landscapes
# insensitive richness - null model
# Total abundace - differences between the two landscapes
# Sensitive abundance - differences between the two landscapes
# insensitive abundance - differences between the two landscapes


#########################################################################
# Data visualization
#########################################################################

# Figure 3 - Comparing the six response variables between the two landscapes

pointcolor <- "gray30"

# Unused script, for jitterplot with boxplot:
#
# par(mfrow=c(3,2), mar=c(3,3,2,2), oma=c(3,3,2,2))
# plot(richness_40 ~ l_, data=dados_cac, range=0, xaxt="n", 
# 	ylab="Number of species")
# points(richness_40 ~ jitter(as.numeric(l_)), data=dados_cac, pch=21, 
# 	bg=pointcolor)
# title("a. Total richness", adj=0)
# 
# plot(richness_sens40 ~ l_, data=dados_cac, range=0, xaxt="n", 
# 	ylab="Number of species")
# points(richness_sens40 ~ jitter(as.numeric(l_)), data=dados_cac, pch=21, 
# 	bg=pointcolor)
# title("b. Richness of sensitive species", adj=0)
# text(x=1.5, y=max(dados_cac$richness_sens40), labels="*", cex=2)
# 
# plot(richness_ns40 ~ l_, data=dados_cac, range=0, xaxt="n", 
# 	ylab="Number of species")
# points(richness_ns40 ~ jitter(as.numeric(l_)), data=dados_cac, pch=21, 
# 	bg=pointcolor)
# title("c. Richness of insensitive species", adj=0)
# 
# plot(abundance_st ~ l_, data=dados_cac, range=0, xaxt="n", 
# 	ylab="Standardized abundance")
# points(abundance_st ~ jitter(as.numeric(l_)), data=dados_cac, pch=21, 
# 	bg=pointcolor)
# title("e. Total standardized abundance", adj=0)
# text(x=1.5, y=max(dados_cac$abundance_st), labels="*", cex=2)
# 
# plot(ab_sens_st ~ l_, data=dados_cac, range=0, xaxt="n", 
# 	ylab="Standardized abundance")
# points(ab_sens_st ~ jitter(as.numeric(l_)), data=dados_cac, pch=21, 
# 	bg=pointcolor)
# title("f. Standardized abundance of sensitive species", adj=0)
# text(x=1.5, y=max(dados_cac$ab_sens_st), labels="*", cex=2)
# 
# plot(ab_notsens_st ~ l_, data=dados_cac, range=0, xaxt="n", 
# 	ylab="Standardized abundance")
# points(ab_notsens_st ~ jitter(as.numeric(l_)), data=dados_cac, pch=21, 
# 	bg=pointcolor)
# title("g. Standardized abundance of insensitive species", adj=0)
# text(x=1.5, y=max(dados_cac$ab_notsens_st), labels="*", cex=2)


png(filename="figure3_landscapes.png", height=25, width=20, 
	unit="cm", res=300)
par(mfcol=c(3,2), mar=c(4,4,2,2), oma=c(3,3,2,2))
pirateplot(richness_40 ~ l_, data=dados_cac, 
	pal="gray", jitter.val=0.05, inf.f.o=0, inf.b.o=0, avg.line.lwd=1, 
	bar.f.o=0,bean.b.o=1, bean.f.o=0.8, 
	bean.f.col=c("#5ab4ac", "#d8b365"), 
	bean.b.col="black", 
	point.col="black", gl.lwd=0, point.o=0.6,gl.lty=5, gl.col="white", 
	ylab="Rarefied number of species", xaxt="n",
	xlab="", point.cex=1, yaxt="n", cex.lab=1.4, bty="l")
axis(side=1,at=c(1,2), labels=c("Forest-\ndominated", "Cacao-\ndominated"),
  tick=FALSE, cex.axis=1.2, line=1.1)
axis(side=2, tck=0.02)
axis(side=1, tck=0.02, labels=c("",""), at=c(1,2))
title(expression(italic("a. Total richness")), 
  adj=0, cex.main=1.4)

pirateplot(richness_sens40 ~ l_, data=dados_cac, 
	pal="gray", jitter.val=0.05, inf.f.o=0, inf.b.o=0, avg.line.lwd=1, 
	bar.f.o=0,bean.b.o=1, bean.f.o=0.8, 
	bean.f.col=c("#5ab4ac", "#d8b365"), 
	bean.b.col="black", 
	point.col="black", gl.lwd=0, point.o=0.6,gl.lty=5, gl.col="white", 
	ylab="Rarefied number of species", xaxt="n",
	xlab="", point.cex=1, yaxt="n", cex.lab=1.4, bty="l")
axis(side=1,at=c(1,2), labels=c("Forest-\ndominated", "Cacao-\ndominated"),
  tick=FALSE, cex.axis=1.2, line=1.1)
axis(side=2, tck=0.02)
axis(side=1, tck=0.02, labels=c("",""), at=c(1,2))
title(expression(italic("b. Richness of sensitive species")), 
  adj=0, cex.main=1.4)
text(x=1.5, y=max(dados_cac$richness_sens40), labels="*", cex=2)

pirateplot(richness_ns40 ~ l_, data=dados_cac, 
	pal="gray", jitter.val=0.05, inf.f.o=0, inf.b.o=0, avg.line.lwd=1, 
	bar.f.o=0,bean.b.o=1, bean.f.o=0.8, 
	bean.f.col=c("#5ab4ac", "#d8b365"), 
	bean.b.col="black", 
	point.col="black", gl.lwd=0, point.o=0.6,gl.lty=5, gl.col="white", 
	ylab="Rarefied number of species", xaxt="n",
	xlab="", point.cex=1, yaxt="n", cex.lab=1.4, bty="l")
axis(side=1,at=c(1,2), labels=c("Forest-\ndominated", "Cacao-\ndominated"),
  tick=FALSE, cex.axis=1.2, line=1.1)
axis(side=2, tck=0.02)
axis(side=1, tck=0.02, labels=c("",""), at=c(1,2))
title(expression(italic("c. Richness of insensitive species")), 
  adj=0, cex.main=1.4)

pirateplot(abundance_st ~ l_, data=dados_cac, 
	pal="gray", jitter.val=0.05, inf.f.o=0, inf.b.o=0, avg.line.lwd=1, 
	bar.f.o=0,bean.b.o=1, bean.f.o=0.8, 
	bean.f.col=c("#5ab4ac", "#d8b365"), 
	bean.b.col="black", 
	point.col="black", gl.lwd=0, point.o=0.6,gl.lty=5, gl.col="white", 
	ylab="Standardized abundance", xaxt="n",
	xlab="", point.cex=1, yaxt="n", cex.lab=1.4, bty="l")
axis(side=1,at=c(1,2), labels=c("Forest-\ndominated", "Cacao-\ndominated"),
  tick=FALSE, cex.axis=1.2, line=1.1)
axis(side=2, tck=0.02)
axis(side=1, tck=0.02, labels=c("",""), at=c(1,2))
title(expression(italic("d. Total abundance")),
  adj=0, cex.main=1.4)
text(x=1.5, y=max(dados_cac$abundance_st), labels="*", cex=2)

pirateplot(ab_sens_st ~ l_, data=dados_cac, 
	pal="gray", jitter.val=0.05, inf.f.o=0, inf.b.o=0, avg.line.lwd=1, 
	bar.f.o=0,bean.b.o=1, bean.f.o=0.8, 
	bean.f.col=c("#5ab4ac", "#d8b365"), 
	bean.b.col="black", 
	point.col="black", gl.lwd=0, point.o=0.6,gl.lty=5, gl.col="white", 
	ylab="Standardized abundance", xaxt="n",
	xlab="", point.cex=1, yaxt="n", cex.lab=1.4, bty="l")
axis(side=1,at=c(1,2), labels=c("Forest-\ndominated", "Cacao-\ndominated"),
  tick=FALSE, cex.axis=1.2, line=1.1)
axis(side=2, tck=0.02)
axis(side=1, tck=0.02, labels=c("",""), at=c(1,2))
title(expression(italic("e. Abundance of sensitive species")), 
  adj=0, cex.main=1.4)
text(x=1.5, y=max(dados_cac$ab_sens_st), labels="*", cex=2)

pirateplot(ab_notsens_st ~ l_, data=dados_cac, 
	pal="gray", jitter.val=0.05, inf.f.o=0, inf.b.o=0, avg.line.lwd=1, 
	bar.f.o=0,bean.b.o=1, bean.f.o=0.8, 
	bean.f.col=c("#5ab4ac", "#d8b365"), 
	bean.b.col="black", 
	point.col="black", gl.lwd=0, point.o=0.6,gl.lty=5, gl.col="white", 
	ylab="Standardized abundance", xaxt="n",
	xlab="", point.cex=1, yaxt="n", cex.lab=1.4, bty="l")
axis(side=1,at=c(1,2), labels=c("Forest-\ndominated", "Cacao-\ndominated"),
  tick=FALSE, cex.axis=1.2, line=1.1)
axis(side=2, tck=0.02)
axis(side=1, tck=0.02, labels=c("",""), at=c(1,2))
title(expression(italic("f. Abundance of insensitive species")), 
  adj=0, cex.main=1.4)
text(x=1.5, y=max(dados_cac$ab_notsens_st), labels="*", cex=2)

mtext(side=1, text="Landscape", outer=T, line=-1)

dev.off()


# Figure 4 - Relation with dog records

dogsNew <- seq(min(dados_cac$dograte_30_log), max(dados_cac$dograte_30_log),
	length=100)

dogsPred <- predict(bestModel.riq.tot, type="response", se.fit=TRUE,
	newdata=list(dograte_30_log=dogsNew))

png(filename="figure4_richnessVSdograte.png", height=17, width=17,
	unit="cm", res=300)
plot(richness_40 ~ dograte_30_log, data=dados_cac, 
	xlab="Dog record rate", 
  ylab="Rarefied species richness", xaxt="n",
	pch=21,  type="n", bty="l", tck=0.02)
polygon(x=c(dogsNew,dogsNew[100:1]),
	y=c(dogsPred$fit-1.96*dogsPred$se.fit,
		dogsPred$fit[100:1] + 1.96*dogsPred$se.fit[100:1]), 
	border=NA, col="gray90")
lines(dogsPred$fit ~ dogsNew)
points(richness_40 ~ dograte_30_log, data=dados_cac, 
	pch=21, bg=ifelse(l_=="HFC ","darkgreen","orange"))
axis(side=1, at=log10(c(0,1,2,5,10,20)+1), labels=c(0,1,2,5,10,20), tck=0.02)
legend(pch=21, pt.bg=c("darkgreen", "orange"), 
	legend=c("Forest-dominated landscape","Cacao-dominated landscape"), 
  x=log10(5), y=9,
	bty="o", cex=0.9)
dev.off()



#########################################################################
# Multivariate analyses
#########################################################################

# Preparing the data

dados.expl <- dados_cac[,c("l_", "density_nat", "BA_frutexo_sqrt", 
	"dograte_30_log")]
effort <- dados_cac$effort
dados.resp <- dados_cac[,c("Ckuhlii", "Cpaca", "Cthous", "Dnovemcinctus", 
	"Daurita", "Dleporina", "Ebarbara", "Esexcinctus", "Lchrysomelas", 
	"Lwiedii", "Llongicaudis", "M.americana", "Mgouazoubira", "Nnasua", 
	"Ptajacu", "Pcancrivorus", "Pconcolor", "Pyagouarondi", 
	"Sbrasiliensis", "Ttetradactyla")]
dados.resp <- as.matrix(dados.resp)

### Remove species found in one or two sites

freq <- function(x) {
  result <- sum(x>0)
  return(result)
}


species.frequencies <- apply(dados.resp, 2, freq)
species.frequencies

species.keep <- species.frequencies >= 3
species.keep

dados.resp2 <- dados.resp[,species.keep]

str(dados.resp2)


# Divide by sampling effort

dados.resp3 <- dados.resp2

for(i in 1:ncol(dados.resp2)) {
  dados.resp3[,i] <- dados.resp2[,i] / effort
}

# Remove rows with no species in them

test.rows <- apply(dados.resp3, 1, sum)
rows.keep <- test.rows > 0
dados.resp4 <- dados.resp3[rows.keep, ]
dados.expl4 <- dados.expl[rows.keep, ]

# Canonical Correspondence Analysis

cca.abund <- cca(dados.resp4 ~ l_ + density_nat + BA_frutexo_sqrt + 
	dograte_30_log, data=dados.expl4)

cca.abund

# First two CCA axes: 0.21, 0.067. Constrained variance: 17% explained.

plot(cca.abund, display=c("sp", "bp"))

# Stepwise removal of explanatory variables

set.seed(12724) # This number was fully arbitrary.

cca.anova1 <- anova(cca.abund, by="margin", permutations=how(nperm=9999))
cca.anova1

# Remove frutexo
cca.abund2 <- cca(dados.resp4 ~ density_nat + l_ + dograte_30_log, data=dados.expl4)
cca.anova2 <- anova(cca.abund2, by="margin", permutations=how(nperm=9999))
cca.anova2

# Remove l_
cca.abund3 <- cca(dados.resp4 ~ density_nat + dograte_30_log, data=dados.expl4)
cca.anova3 <- anova(cca.abund3, by="margin", permutations=how(nperm=9999))
cca.anova3

# Dograte marginally significant (p=0.0921)


plot(cca.abund3, display=c("sp", "bp"))


# Data visualization

spnames.4 <- c("Calkuh", "Cunpac", "Certho", "Dasnov", "Didaur", "Eirbar",
	"Leochr", "Leowie", "Mazgou", "Nasnas", "Procan", "Tamtet")

scores.abund3 <- scores(cca.abund3, display=c("sp", "wa", "bp"), scaling=3)

png(filename="figure5_CCA.png", height=20, width=20, unit="cm", res=300)
plot(scores.abund3$sites[,2] ~ scores.abund3$sites[,1],	type="p", pch=21, 
	bg = alpha(ifelse(dados.expl4$l_=="HFC ", "#5ab4ac", "#d8b365"),0.9),
	xlab="CCA 1 (11%)", ylab= "CCA 2 (03%)", bty="l")
segments(x0=c(0,0), y0=c(0,0), x1=scores.abund3$biplot[,1]*2, 
	y1=scores.abund3$biplot[,2]*2, col=c("black", "black"), lwd=2)
text(x=scores.abund3$biplot[,1]*2, y=scores.abund3$biplot[,2]*2,
	labels=c(expression(italic("Native trees\n(p=0.0004)")), 
    expression(italic("Dog record rate\n(p=0.0921)"))), 
	col="black", pos=c(1,4))
text(x=scores.abund3$species[,1], y=scores.abund3$species[,2],
	labels=spnames.4, cex=0.8)
legend(pch=21, pt.bg=c("#5ab4ac", "#d8b365"), 
  legend=c("Forest-dominated landscape",
	"Cacao-dominated landscape"), x="bottomright")
dev.off()

#########################################################################
# The end.
#########################################################################