library(mgcv)
library(gamm4)
library(lme4)
library(bbmle)
library(statmod)
library(scales)

setwd("/home/pavel/Profissional/Pesquisa/MyPapers_working/Janaine_2019_EILitterUna")

setwd("e:\\Pavel\\Profissional\\Pesquisa\\MyPapers_Working\\Janaine_2019_EILitterUna\\")


dados <- read.table("dadosSerapilheira.txt", header=T, sep="\t")

str(dados)

### Figuras exploratórias iniciais
par(mfrow=c(2,2), mar=c(3,3,2,2), oma=c(2,2,2,2))
plot(Peso_folhas ~ Distancia, data=dados, pch=21, bg=Transecto, xlab="", ylab="", main="Folhas")
abline(v=0, col="gray", lty=3)
plot(Peso_galhos ~ Distancia, data=dados, pch=21, bg=Transecto, xlab="", ylab="", main="Galhos")
abline(v=0, col="gray", lty=3)
plot(Peso_misc ~ Distancia, data=dados, pch=21, bg=Transecto, xlab="", ylab="", main="Miscelânia")
abline(v=0, col="gray", lty=3)
plot(Peso_total ~ Distancia, data=dados, pch=21, bg=Transecto, xlab="", ylab="", main="Total")
abline(v=0, col="gray", lty=3)
mtext(side=1, text="Distância", outer=T)
mtext(side=2, text="Biomassa (g)", outer=T)

### Calcular média de cada variável por parcela

mean.folhas <- aggregate(Peso_folhas ~ Transecto + Distancia, FUN=mean, data=dados)
mean.galhos <- aggregate(Peso_galhos ~ Transecto + Distancia, FUN=mean, data=dados)
mean.misc <- aggregate(Peso_misc ~ Transecto + Distancia, FUN=mean, data=dados)
mean.total <- aggregate(Peso_total ~ Transecto + Distancia, FUN=mean, data=dados)

dados.mean <- merge(mean.folhas, mean.galhos)
dados.mean <- merge(dados.mean, mean.misc)
dados.mean <- merge(dados.mean, mean.total)

str(dados.mean)

### Figuras com base na média
par(mfrow=c(2,2), mar=c(3,3,2,2), oma=c(2,2,2,2))
plot(Peso_folhas ~ Distancia, data=dados.mean, pch=21, bg=Transecto, xlab="", ylab="", main="Folhas")
abline(v=0, col="gray", lty=3)
plot(Peso_galhos ~ Distancia, data=dados.mean, pch=21, bg=Transecto, xlab="", ylab="", main="Galhos")
abline(v=0, col="gray", lty=3)
plot(Peso_misc ~ Distancia, data=dados.mean, pch=21, bg=Transecto, xlab="", ylab="", main="Miscelânia")
abline(v=0, col="gray", lty=3)
plot(Peso_total ~ Distancia, data=dados.mean, pch=21, bg=Transecto, xlab="", ylab="", main="Total")
abline(v=0, col="gray", lty=3)
mtext(side=1, text="Distância", outer=T)
mtext(side=2, text="Biomassa média (g)", outer=T)

### Análise dos dados
### Modelos:
### Null: sem efeito do local ou distância
### Cat_edgeAsForest: diferença entre área queimada e não queimada, borda como floresta
### Cat_edgeAsBurnt: diferença entre área queimada e não queimada, borda como área queimada
### Gradient: modelo aditivo, com uma variação gradual ao longo do transecto
### CatGrad_edgeAsForest: diferença entre área queimada e não queimada + modelo aditivo, borda como floresta
### CatGrad_edgeAsBurnt: diferença entre área queimada e não queimada + modelo aditivo, borda como área queimada

dados.mean$Ambiente_edgeAsForest <- as.factor(ifelse(dados.mean$Distancia < 0, "Burnt", "Forest"))
dados.mean$Ambiente_edgeAsBurnt <- as.factor(ifelse(dados.mean$Distancia < 1, "Burnt", "Forest"))

### Folhas ###

mod.folhas.null <- glmer(Peso_folhas ~ 1 + (1 | Transecto), data=dados.mean, family=Gamma)
mod.folhas.cat_edgeAsForest <- glmer(Peso_folhas ~ Ambiente_edgeAsForest + (1|Transecto), data=dados.mean, family=Gamma)
mod.folhas.cat_edgeAsBurnt <- glmer(Peso_folhas ~ Ambiente_edgeAsBurnt + (1|Transecto), data=dados.mean, family=Gamma)
mod.folhas.gam <- gamm4(Peso_folhas ~ s(Distancia, fx=F, k=4), random=~(1|Transecto), data=dados.mean, family=Gamma)
mod.folhas.gamcat_edgeAsForest <- gamm4(Peso_folhas ~ s(Distancia, fx=F, k=4) + Ambiente_edgeAsForest, random=~(1|Transecto), data=dados.mean, family=Gamma)
mod.folhas.gamcat_edgeAsBurnt <- gamm4(Peso_folhas ~ s(Distancia, fx=F, k=4) + Ambiente_edgeAsBurnt, random=~(1|Transecto), data=dados.mean, family=Gamma)

AICctab(mod.folhas.null, mod.folhas.cat_edgeAsForest, mod.folhas.cat_edgeAsBurnt, mod.folhas.gam$mer, mod.folhas.gamcat_edgeAsForest$mer, mod.folhas.gamcat_edgeAsBurnt$mer)

### Categorical model selected; edge as burnt.

### Galhos ###

### Select the optimal p for the Tweedie family: a null model, without random factors.

mod.galhos.foo <- gam(Peso_galhos ~ 1, data=dados.mean, family=tw)

summary(mod.galhos.foo)

# p = 1.72


mod.galhos.null <- gamm4(Peso_galhos ~ 1, random=~(1|Transecto), data=dados.mean, family=Tweedie(p=1.72))
mod.galhos.cat_edgeAsForest <- gamm4(Peso_galhos ~ Ambiente_edgeAsForest, random=~(1|Transecto), data=dados.mean, family=Tweedie(p=1.72))
mod.galhos.cat_edgeAsBurnt <- gamm4(Peso_galhos ~ Ambiente_edgeAsBurnt, random=~(1|Transecto), data=dados.mean, family=Tweedie(p=1.72))
mod.galhos.gam <- gamm4(Peso_galhos ~ s(Distancia, fx=F, k=4), random=~(1|Transecto), data=dados.mean, family=Tweedie(p=1.72))
mod.galhos.gamcat_edgeAsForest <- gamm4(Peso_galhos ~ s(Distancia, fx=F, k=4) + Ambiente_edgeAsForest, random=~(1|Transecto), data=dados.mean, family=Tweedie(p=1.72))
mod.galhos.gamcat_edgeAsBurnt <- gamm4(Peso_galhos ~ s(Distancia, fx=F, k=4) + Ambiente_edgeAsBurnt, random=~(1|Transecto), data=dados.mean, family=Tweedie(p=1.72))

AICctab(mod.galhos.null$mer, mod.galhos.cat_edgeAsForest$mer, mod.galhos.cat_edgeAsBurnt$mer, mod.galhos.gam$mer, mod.galhos.gamcat_edgeAsForest$mer, mod.galhos.gamcat_edgeAsBurnt$mer)


# Null model selected

### Misc ###

mod.misc.foo <- gam(Peso_misc ~ 1, data=dados.mean, family=tw)

summary(mod.misc.foo)

# p = 1.659

mod.misc.null <- gamm4(Peso_misc ~ 1, random=~(1|Transecto), data=dados.mean, family=Tweedie(p=1.659))
mod.misc.cat_edgeAsForest <- gamm4(Peso_misc ~ Ambiente_edgeAsForest, random=~(1|Transecto), data=dados.mean, family=Tweedie(p=1.659))
mod.misc.cat_edgeAsBurnt <- gamm4(Peso_misc ~ Ambiente_edgeAsBurnt, random=~(1|Transecto), data=dados.mean, family=Tweedie(p=1.659))
mod.misc.gam <- gamm4(Peso_misc ~ s(Distancia, fx=F, k=4), random=~(1|Transecto), data=dados.mean, family=Tweedie(p=1.659))
mod.misc.gamcat_edgeAsForest <- gamm4(Peso_misc ~ s(Distancia, fx=F, k=4) + Ambiente_edgeAsForest, random=~(1|Transecto), data=dados.mean, family=Tweedie(p=1.659))
mod.misc.gamcat_edgeAsBurnt <- gamm4(Peso_misc ~ s(Distancia, fx=F, k=4) + Ambiente_edgeAsBurnt, random=~(1|Transecto), data=dados.mean, family=Tweedie(p=1.659))

AICctab(mod.misc.null$mer, mod.misc.cat_edgeAsForest$mer, mod.misc.cat_edgeAsBurnt$mer, mod.misc.gam$mer, mod.misc.gamcat_edgeAsForest$mer, mod.misc.gamcat_edgeAsBurnt$mer)

# Categorical - edge as burnt selected.

### Total ###


mod.total.null <- glmer(Peso_total ~ 1 + (1 | Transecto), data=dados.mean, family=Gamma)
mod.total.cat_edgeAsForest <- glmer(Peso_total ~ Ambiente_edgeAsForest + (1|Transecto), data=dados.mean, family=Gamma)
mod.total.cat_edgeAsBurnt <- glmer(Peso_total ~ Ambiente_edgeAsBurnt + (1|Transecto), data=dados.mean, family=Gamma)
mod.total.gam <- gamm4(Peso_total ~ s(Distancia, fx=F, k=4), random=~(1|Transecto), data=dados.mean, family=Gamma)
mod.total.gamcat_edgeAsForest <- gamm4(Peso_total ~ s(Distancia, fx=F, k=4) + Ambiente_edgeAsForest, random=~(1|Transecto), data=dados.mean, family=Gamma)
mod.total.gamcat_edgeAsBurnt <- gamm4(Peso_total ~ s(Distancia, fx=F, k=4) + Ambiente_edgeAsBurnt, random=~(1|Transecto), data=dados.mean, family=Gamma)

AICctab(mod.total.null, mod.total.cat_edgeAsForest, mod.total.cat_edgeAsBurnt, mod.total.gam$mer, mod.total.gamcat_edgeAsForest$mer, mod.total.gamcat_edgeAsBurnt$mer)

# Categorical - edge as burnt.

### Thus, there is less leaf litter, miscelaneous litter and total litter in the burnt area. Edge is similar to the burnt area.

### Calcular médias para área queimada e floresta, com borda como área queimada

mean.folhas.burnt <- mean(subset(dados.mean, Ambiente_edgeAsBurnt=="Burnt")$Peso_folhas)
mean.folhas.forest <- mean(subset(dados.mean, Ambiente_edgeAsBurnt=="Forest")$Peso_folhas)

mean.galhos <- mean(dados.mean$Peso_galhos)

mean.misc.burnt <- mean(subset(dados.mean, Ambiente_edgeAsBurnt=="Burnt")$Peso_misc)
mean.misc.forest <- mean(subset(dados.mean, Ambiente_edgeAsBurnt=="Forest")$Peso_misc)

mean.total.burnt <- mean(subset(dados.mean, Ambiente_edgeAsBurnt=="Burnt")$Peso_total)
mean.total.forest <- mean(subset(dados.mean, Ambiente_edgeAsBurnt=="Forest")$Peso_total)


### Figura com os resultados do modelo

png(filename="EIlitter_results_grayscale.png", height=20, width=20, res=300, unit="cm")
par(mfrow=c(2,2), mar=c(3,3,2,2), oma=c(2,2,2,2))
plot(Peso_folhas ~ Distancia, data=dados.mean, pch=16, col=alpha("black",0.5), xlab="", ylab="", main="a. Leaf litter")
abline(v=0, col="gray", lty=3, lwd=2)
segments(x0=-150, x1=5, y0=mean.folhas.burnt, y1=mean.folhas.burnt, col="black")
segments(x0=5, x1=150, y0=mean.folhas.forest, y1=mean.folhas.forest, col="black")
plot(Peso_galhos ~ Distancia, data=dados.mean, pch=16, col=alpha("black",0.5), xlab="", ylab="", main="b. Fine woody debris (FWD)")
abline(v=0, col="gray", lty=3, lwd=2)
segments(x0=-150, x1=150, y0=mean.galhos, y1=mean.galhos, col="black")
plot(Peso_misc ~ Distancia, data=dados.mean, pch=16, col=alpha("black",0.5), xlab="", ylab="", main="c. Miscelaneous")
abline(v=0, col="gray", lty=3, lwd=2)
segments(x0=-150, x1=5, y0=mean.misc.burnt, y1=mean.misc.burnt, col="black")
segments(x0=5, x1=150, y0=mean.misc.forest, y1=mean.misc.forest, col="black")
plot(Peso_total ~ Distancia, data=dados.mean, pch=16, col=alpha("black",0.5), xlab="", ylab="", main="d. Total")
abline(v=0, col="gray", lty=3, lwd=2)
segments(x0=-150, x1=5, y0=mean.total.burnt, y1=mean.total.burnt, col="black")
segments(x0=5, x1=150, y0=mean.total.forest, y1=mean.total.forest, col="black")
mtext(side=1, text="Distance (m)", outer=T)
mtext(side=2, text="Mean biomass (g)", outer=T)
dev.off()

### Análise dos dados
