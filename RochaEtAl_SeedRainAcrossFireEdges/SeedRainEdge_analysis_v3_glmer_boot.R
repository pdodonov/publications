library(mgcv)
library(nlme)
library(bbmle)
library(lme4)
library(MASS)
library(gamm4)

setwd("/home/pavel/Profissional/Pesquisa/MyPapers/2018_Janaine_EISeedRain")

dados <- read.table("SeedRainEdge_data_all.txt", header=T, sep="\t")

### Add vegetation

dados$Vegetation_edgeAsForest <- as.factor(ifelse(dados$Distance<0,"Fire","Forest"))
dados$Vegetation_edgeAsFire <- as.factor(ifelse(dados$Distance<1,"Fire","Forest"))

# Replace zero weight with 0.0001

dados$Weight_all[dados$Weight_all==0 & !is.na(dados$Weight_all)] <- 0.0001
dados$Weight_nzoo[dados$Weight_nzoo==0 & !is.na(dados$Weight_nzoo)] <- 0.0001



### t-test and boostrap to compare between burnt and uburnt areas
dados.noEdge <- subset(dados, Distance != 0) # Remove the edge
dados.Fire <- subset(dados, Distance < 0)
dados.Forest <- subset(dados, Distance > 0)


str(dados.noEdge)
dados.noEdge$Vegetation <- ifelse(dados.noEdge$Distance < 0, "Fire", "Forest")

x <- dados.noEdge$Abundance_tot
y <- dados.noEdge$Vegetation
group <- dados.noEdge$Transect

t.test.restRand <- function(x, y, group, Nperm=1E4) { #function to perform a t-test with restricted randomizations
	groups <- unique(group)
	y.levels <- unique(y)
	mean.1.real <- mean(x[y==y.levels[1]], na.rm=TRUE)
	mean.2.real <- mean(x[y==y.levels[2]], na.rm=TRUE)
	diff.real <- mean.1.real - mean.2.real
	diff.perm <- numeric(Nperm)
	diff.perm[1] <- diff.real
	for(i in 2:Nperm) {
		for(j in 1:length(groups)) {
			y.foo <- subset(y, group == groups[j])
			y.foo <- sample(y.foo)
			if(j == 1) {
				y.rand <- y.foo
				} else {
					y.rand <- c(y.rand, y.foo)
			}
		}
		mean.1.perm <- mean(x[y.rand==y.levels[1]], na.rm=TRUE)
		mean.2.perm <- mean(x[y.rand==y.levels[2]], na.rm=TRUE)
		diff.perm[i] <- mean.1.perm - mean.2.perm
	}
	signif <- sum(abs(diff.perm) >= abs(diff.real))/Nperm
	result <- c(mean.1.real, mean.2.real, diff.real, signif)
	names(result) <- c(y.levels, "diff.real", "signif")
	return(result)
}

boostrap.restricted <- function(x, group, Nboot=1E4) {
	groups <- unique(group)
	mean.boot <- numeric(Nboot)
	mean.real <- mean(x, na.rm=TRUE)
	mean.boot[1] <- mean.real
	for(i in 2:Nboot) {
		for(j in 1:length(groups)) {
			x.foo <- subset(x, group == groups[j])
			x.foo <- sample(x.foo, replace=TRUE)
			if(j == 1) {
				x.boot <- x.foo
				} else {
					x.boot <- c(x.boot, x.foo)
			}
		}
		mean.boot[i] <- mean(x.boot, na.rm=TRUE)
	}
	result <- quantile(mean.boot, c(0.025, 0.975))
	return(result)
}



### Comparing burnt area and unburnt forest

t.test.Abundance_tot <- t.test.restRand(x=dados.noEdge$Abundance_tot, y=dados.noEdge$Vegetation, group=dados.noEdge$Transect)
t.test.Richness_mean <- t.test.restRand(x=dados.noEdge$Richness_mean, y=dados.noEdge$Vegetation, group=dados.noEdge$Transect)
t.test.Richness_mean_zoo <- t.test.restRand(x=dados.noEdge$Richness_mean_zoo, y=dados.noEdge$Vegetation, group=dados.noEdge$Transect)
t.test.Richness_mean_nzoo <- t.test.restRand(x=dados.noEdge$Richness_mean_nzoo, y=dados.noEdge$Vegetation, group=dados.noEdge$Transect)
t.test.Richness_tot <- t.test.restRand(x=dados.noEdge$Richness_tot, y=dados.noEdge$Vegetation, group=dados.noEdge$Transect)
t.test.Richness_tot_zoo <- t.test.restRand(x=dados.noEdge$Richness_tot_zoo, y=dados.noEdge$Vegetation, group=dados.noEdge$Transect)
t.test.Richness_tot_nzoo <- t.test.restRand(x=dados.noEdge$Richness_tot_nzoo, y=dados.noEdge$Vegetation, group=dados.noEdge$Transect)
t.test.Abundance_tot_zoo <- t.test.restRand(x=dados.noEdge$Abundance_tot_zoo, y=dados.noEdge$Vegetation, group=dados.noEdge$Transect)
t.test.Abundance_tot_nzoo <- t.test.restRand(x=dados.noEdge$Abundance_tot_nzoo, y=dados.noEdge$Vegetation, group=dados.noEdge$Transect)
t.test.Width_all <- t.test.restRand(x=dados.noEdge$Width_all, y=dados.noEdge$Vegetation, group=dados.noEdge$Transect)
t.test.Weight_all <- t.test.restRand(x=dados.noEdge$Weight_all, y=dados.noEdge$Vegetation, group=dados.noEdge$Transect)
t.test.Width_zoo <- t.test.restRand(x=dados.noEdge$Width_zoo, y=dados.noEdge$Vegetation, group=dados.noEdge$Transect)
t.test.Weight_nzoo <- t.test.restRand(x=dados.noEdge$Weight_nzoo, y=dados.noEdge$Vegetation, group=dados.noEdge$Transect)


# t.test.Abundance_tot
# t.test.Richness_mean
# t.test.Richness_mean_zoo
# t.test.Richness_mean_nzoo
# t.test.Richness_tot
# t.test.Richness_tot_zoo
# t.test.Richness_tot_nzoo
# t.test.Abundance_tot_zoo
# t.test.Abundance_tot_nzoo
# t.test.Width_all
# t.test.Weight_all
# t.test.Width_zoo
# t.test.Weight_nzoo


### Combine into a matrix
t.test.all <- matrix(c(t.test.Abundance_tot, t.test.Richness_mean, t.test.Richness_mean_zoo, t.test.Richness_mean_nzoo, t.test.Richness_tot, t.test.Richness_tot_zoo, t.test.Richness_tot_nzoo, t.test.Abundance_tot_zoo, t.test.Abundance_tot_nzoo, t.test.Width_all, t.test.Weight_all, t.test.Width_zoo, t.test.Weight_nzoo), ncol=4, byrow=T)
colnames(t.test.all) <- names(t.test.Abundance_tot)
row.names(t.test.all) <- c("Abundance_tot", "Richness_mean", "Richness_mean_zoo", "Richness_mean_nzoo", "Richness_tot", "Richness_tot_zoo", "Richness_tot_nzoo", "Abundance_tot_zoo", "Abundance_tot_nzoo", "Width_all", "Weight_all", "Width_zoo", "Weight_nzoo")



bootstrap.Abundance_tot.Fire <- boostrap.restricted (x=dados.Fire$Abundance_tot, group=dados.Fire$Transect)
bootstrap.Richness_mean.Fire <- boostrap.restricted (x=dados.Fire$Richness_mean, group=dados.Fire$Transect)
bootstrap.Richness_mean_zoo.Fire <- boostrap.restricted (x=dados.Fire$Richness_mean_zoo, group=dados.Fire$Transect)
bootstrap.Richness_mean_nzoo.Fire <- boostrap.restricted (x=dados.Fire$Richness_mean_nzoo, group=dados.Fire$Transect)
bootstrap.Richness_tot.Fire <- boostrap.restricted (x=dados.Fire$Richness_tot, group=dados.Fire$Transect)
bootstrap.Richness_tot_zoo.Fire <- boostrap.restricted (x=dados.Fire$Richness_tot_zoo, group=dados.Fire$Transect)
bootstrap.Richness_tot_nzoo.Fire <- boostrap.restricted (x=dados.Fire$Richness_tot_nzoo, group=dados.Fire$Transect)
bootstrap.Abundance_tot_zoo.Fire <- boostrap.restricted (x=dados.Fire$Abundance_tot_zoo, group=dados.Fire$Transect)
bootstrap.Abundance_tot_nzoo.Fire <- boostrap.restricted (x=dados.Fire$Abundance_tot_nzoo, group=dados.Fire$Transect)
bootstrap.Width_all.Fire <- boostrap.restricted (x=dados.Fire$Width_all, group=dados.Fire$Transect)
bootstrap.Weight_all.Fire <- boostrap.restricted (x=dados.Fire$Weight_all, group=dados.Fire$Transect)
bootstrap.Width_zoo.Fire <- boostrap.restricted (x=dados.Fire$Width_zoo, group=dados.Fire$Transect)
bootstrap.Weight_nzoo.Fire <- boostrap.restricted (x=dados.Fire$Weight_nzoo, group=dados.Fire$Transect)

bootstrap.Abundance_tot.Forest <- boostrap.restricted (x=dados.Forest$Abundance_tot, group=dados.Forest$Transect)
bootstrap.Richness_mean.Forest <- boostrap.restricted (x=dados.Forest$Richness_mean, group=dados.Forest$Transect)
bootstrap.Richness_mean_zoo.Forest <- boostrap.restricted (x=dados.Forest$Richness_mean_zoo, group=dados.Forest$Transect)
bootstrap.Richness_mean_nzoo.Forest <- boostrap.restricted (x=dados.Forest$Richness_mean_nzoo, group=dados.Forest$Transect)
bootstrap.Richness_tot.Forest <- boostrap.restricted (x=dados.Forest$Richness_tot, group=dados.Forest$Transect)
bootstrap.Richness_tot_zoo.Forest <- boostrap.restricted (x=dados.Forest$Richness_tot_zoo, group=dados.Forest$Transect)
bootstrap.Richness_tot_nzoo.Forest <- boostrap.restricted (x=dados.Forest$Richness_tot_nzoo, group=dados.Forest$Transect)
bootstrap.Abundance_tot_zoo.Forest <- boostrap.restricted (x=dados.Forest$Abundance_tot_zoo, group=dados.Forest$Transect)
bootstrap.Abundance_tot_nzoo.Forest <- boostrap.restricted (x=dados.Forest$Abundance_tot_nzoo, group=dados.Forest$Transect)
bootstrap.Width_all.Forest <- boostrap.restricted (x=dados.Forest$Width_all, group=dados.Forest$Transect)
bootstrap.Weight_all.Forest <- boostrap.restricted (x=dados.Forest$Weight_all, group=dados.Forest$Transect)
bootstrap.Width_zoo.Forest <- boostrap.restricted (x=dados.Forest$Width_zoo, group=dados.Forest$Transect)
bootstrap.Weight_nzoo.Forest <- boostrap.restricted (x=dados.Forest$Weight_nzoo, group=dados.Forest$Transect)



#bootstrap.Abundance_tot.Fire
#bootstrap.Richness_mean.Fire
#bootstrap.Richness_mean_zoo.Fire
#bootstrap.Richness_mean_nzoo.Fire
#bootstrap.Richness_tot.Fire
#bootstrap.Richness_tot_zoo.Fire
#bootstrap.Richness_tot_nzoo.Fire
#bootstrap.Abundance_tot_zoo.Fire
#bootstrap.Abundance_tot_nzoo.Fire
#bootstrap.Width_all.Fire
#bootstrap.Weight_all.Fire
#bootstrap.Width_zoo.Fire
#bootstrap.Weight_nzoo.Fire
#bootstrap.Abundance_tot.Forest
#bootstrap.Richness_mean.Forest
#bootstrap.Richness_mean_zoo.Forest
#bootstrap.Richness_mean_nzoo.Forest
#bootstrap.Richness_tot.Forest
#bootstrap.Richness_tot_zoo.Forest
#bootstrap.Richness_tot_nzoo.Forest
#bootstrap.Abundance_tot_zoo.Forest
#bootstrap.Abundance_tot_nzoo.Forest
#bootstrap.Width_all.Forest
#bootstrap.Weight_all.Forest
#bootstrap.Width_zoo.Forest
#bootstrap.Weight_nzoo.Forest


bootstrap.Fire.all <- matrix(c(bootstrap.Abundance_tot.Fire, bootstrap.Richness_mean.Fire, bootstrap.Richness_mean_zoo.Fire, bootstrap.Richness_mean_nzoo.Fire, bootstrap.Richness_tot.Fire, bootstrap.Richness_tot_zoo.Fire, bootstrap.Richness_tot_nzoo.Fire, bootstrap.Abundance_tot_zoo.Fire, bootstrap.Abundance_tot_nzoo.Fire, bootstrap.Width_all.Fire, bootstrap.Weight_all.Fire, bootstrap.Width_zoo.Fire, bootstrap.Weight_nzoo.Fire), ncol=2, byrow=T)
colnames(bootstrap.Fire.all) <- names(bootstrap.Abundance_tot.Fire)
row.names(bootstrap.Fire.all) <- c("Abundance_tot", "Richness_mean", "Richness_mean_zoo", "Richness_mean_nzoo", "Richness_tot", "Richness_tot_zoo", "Richness_tot_nzoo", "Abundance_tot_zoo", "Abundance_tot_nzoo", "Width_all", "Weight_all", "Width_zoo", "Weight_nzoo")


bootstrap.Forest.all <- matrix(c(bootstrap.Abundance_tot.Forest, bootstrap.Richness_mean.Forest, bootstrap.Richness_mean_zoo.Forest, bootstrap.Richness_mean_nzoo.Forest, bootstrap.Richness_tot.Forest, bootstrap.Richness_tot_zoo.Forest, bootstrap.Richness_tot_nzoo.Forest, bootstrap.Abundance_tot_zoo.Forest, bootstrap.Abundance_tot_nzoo.Forest, bootstrap.Width_all.Forest, bootstrap.Weight_all.Forest, bootstrap.Width_zoo.Forest, bootstrap.Weight_nzoo.Forest), ncol=2, byrow=T)
colnames(bootstrap.Forest.all) <- names(bootstrap.Abundance_tot.Forest)
row.names(bootstrap.Forest.all) <- c("Abundance_tot", "Richness_mean", "Richness_mean_zoo", "Richness_mean_nzoo", "Richness_tot", "Richness_tot_zoo", "Richness_tot_nzoo", "Abundance_tot_zoo", "Abundance_tot_nzoo", "Width_all", "Weight_all", "Width_zoo", "Weight_nzoo")





### Adjust the models...

Nvar <- ncol(dados)-2-2
Nmodels <- 6

response.nb <- matrix(ncol=Nmodels, nrow=Nvar)

colnames(response.nb) <- c("Null", "Cat_edgeAsForest", "Cat_edgeAsFire", "Gam", "GamCat_edgeAsForest", "GamCat_edgeAsFire")
row.names(response.nb) <- colnames(dados)[-c(1,2,16,17)]


### Transform mean values into absolute numbers

# 01 "Abundance_tot"          
# 02 "Richness_mean"          
# 03 "Richness_mean_zoo"      
# 04 "Richness_mean_nzoo"     
# 05 "Richness_tot"           
# 06 "Richness_tot_zoo"       
# 07 "Richness_tot_nzoo"      
# 08 "Abundance_tot_zoo"      
# 09 "Abundance_tot_nzoo"     
# 10 "Width_all"              
# 11 "Weight_all"             
# 12 "Width_zoo"              
# 13 "Weight_nzoo"            

for(i in 1:Nvar) {
	print(i)
	y <- dados[,i+2]
	dados.use <- subset(dados, !is.na(y))
	y.use <- y[!is.na(y)]
	if(i %in% c(1,5,6,7,8,9)) {
		mod.i.null <- glmer.nb(y.use ~ 1 + (1 | Transect), data=dados.use, method="ML", control=glmerControl(tolPwrss=0.00001, optCtrl = list(maxfun=1E4)))
		mod.i.cat_edgeAsForest <- glmer.nb(y.use ~ Vegetation_edgeAsForest + (1|Transect), data=dados.use, method="ML", control=glmerControl(tolPwrss=0.00001, optCtrl = list(maxfun=1E4)))
		mod.i.cat_edgeAsFire <- glmer.nb(y.use ~ Vegetation_edgeAsFire + (1|Transect), data=dados.use, method="ML", control=glmerControl(tolPwrss=0.00001, optCtrl = list(maxfun=1E4)))

		theta.i.null <- getME(mod.i.null, "glmer.nb.theta")
		theta.i.cat_edgeAsForest <- getME(mod.i.cat_edgeAsForest, "glmer.nb.theta")
		theta.i.cat_edgeAsFire <- getME(mod.i.cat_edgeAsFire, "glmer.nb.theta")

		mod.i.gam <- gamm4(y.use ~ s(Distance, fx=F, k=4), random=~(1|Transect), data=dados.use, method="ML", family=negbin(theta=theta.i.null), control=glmerControl(tolPwrss=0.00001, optCtrl = list(maxfun=1E4)))
		mod.i.gamcat_edgeAsForest <- gamm4(y.use ~ s(Distance, fx=F, k=4) + Vegetation_edgeAsForest, random=~(1|Transect), data=dados.use, method="ML", family=negbin(theta=theta.i.cat_edgeAsForest), control=glmerControl(tolPwrss=0.00001, optCtrl = list(maxfun=1E4)))
		mod.i.gamcat_edgeAsFire <- gamm4(y.use ~ s(Distance, fx=F, k=4) + Vegetation_edgeAsFire, random=~(1|Transect), data=dados.use, method="ML", family=negbin(theta=theta.i.cat_edgeAsFire), control=glmerControl(tolPwrss=0.00001, optCtrl = list(maxfun=1E4)))
	}
	if(i %in% c(2,3,4)) {
		mod.i.null <- lmer(y.use ~ 1 + (1 | Transect), data=dados.use, REML=F)
		mod.i.cat_edgeAsForest <- lmer(y.use ~ Vegetation_edgeAsForest + (1|Transect), data=dados.use, REML=F)
		mod.i.cat_edgeAsFire <- lmer(y.use ~ Vegetation_edgeAsFire + (1|Transect), data=dados.use, REML=F)
		mod.i.gam <- gamm4(y.use ~ s(Distance, fx=F, k=4), random=~(1|Transect), data=dados.use, REML=F, family=gaussian)
		mod.i.gamcat_edgeAsForest <- gamm4(y.use ~ s(Distance, fx=F, k=4) + Vegetation_edgeAsForest, random=~(1|Transect), data=dados.use, REML=F, family=gaussian)
		mod.i.gamcat_edgeAsFire <- gamm4(y.use ~ s(Distance, fx=F, k=4) + Vegetation_edgeAsFire, random=~(1|Transect), data=dados.use, REML=F, family=gaussian)

	}
	if(i %in% c(10,11,12,13)) {
		mod.i.null <- glmer(y.use ~ 1 + (1 | Transect), data=dados.use, method="ML", control=glmerControl(tolPwrss=1E-6, optCtrl = list(maxfun=1E4)), family=Gamma)
		mod.i.cat_edgeAsForest <- glmer(y.use ~ Vegetation_edgeAsForest + (1|Transect), data=dados.use, method="ML", control=glmerControl(tolPwrss=1E-6, optCtrl = list(maxfun=1E4)), family=Gamma)
		mod.i.cat_edgeAsFire <- glmer(y.use ~ Vegetation_edgeAsFire + (1|Transect), data=dados.use, method="ML", control=glmerControl(tolPwrss=1E-6, optCtrl = list(maxfun=1E4)), family=Gamma)
		mod.i.gam <- gamm4(y.use ~ s(Distance, fx=F, k=4), random=~(1|Transect), data=dados.use, method="ML", family=Gamma, control=glmerControl(tolPwrss=1E-6, optCtrl = list(maxfun=1E4)))
		mod.i.gamcat_edgeAsForest <- gamm4(y.use ~ s(Distance, fx=F, k=4) + Vegetation_edgeAsForest, random=~(1|Transect), data=dados.use, method="ML", family=Gamma, control=glmerControl(tolPwrss=1E-6, optCtrl = list(maxfun=1E4)))
		mod.i.gamcat_edgeAsFire <- gamm4(y.use ~ s(Distance, fx=F, k=4) + Vegetation_edgeAsFire, random=~(1|Transect), data=dados.use, method="ML", family=Gamma, control=glmerControl(tolPwrss=1E-6, optCtrl = list(maxfun=1E4)))

	}
	
	AICc.i<- AICctab(mod.i.null, mod.i.cat_edgeAsForest, mod.i.cat_edgeAsFire, mod.i.gam$mer, mod.i.gamcat_edgeAsForest$mer, mod.i.gamcat_edgeAsFire$mer, nobs=nrow(dados.use), sort=F)
	response.nb[i,] <- AICc.i$dAICc
}


# There was an error in the last iteration... Readjust the model, lowering the tolerance in tolPwrss
i <- 13
y <- dados[,i+2]
dados.use <- subset(dados, !is.na(y))
y.use <- y[!is.na(y)]
mod.i.null <- glmer(y.use ~ 1 + (1 | Transect), data=dados.use, method="ML", control=glmerControl(tolPwrss=1E-6, optCtrl = list(maxfun=1E4)), family=Gamma)
mod.i.cat_edgeAsForest <- glmer(y.use ~ Vegetation_edgeAsForest + (1|Transect), data=dados.use, method="ML", control=glmerControl(tolPwrss=1E-6, optCtrl = list(maxfun=1E4)), family=Gamma)
mod.i.cat_edgeAsFire <- glmer(y.use ~ Vegetation_edgeAsFire + (1|Transect), data=dados.use, method="ML", control=glmerControl(tolPwrss=1E-6, optCtrl = list(maxfun=1E4)), family=Gamma)
mod.i.gam <- gamm4(y.use ~ s(Distance, fx=F, k=4), random=~(1|Transect), data=dados.use, method="ML", family=Gamma, control=glmerControl(tolPwrss=1E-4, optCtrl = list(maxfun=1E4)))
mod.i.gamcat_edgeAsForest <- gamm4(y.use ~ s(Distance, fx=F, k=4) + Vegetation_edgeAsForest, random=~(1|Transect), data=dados.use, method="ML", family=Gamma, control=glmerControl(tolPwrss=1E-2, optCtrl = list(maxfun=1E4)))
mod.i.gamcat_edgeAsFire <- gamm4(y.use ~ s(Distance, fx=F, k=4) + Vegetation_edgeAsFire, random=~(1|Transect), data=dados.use, method="ML", family=Gamma, control=glmerControl(tolPwrss=1E-2, optCtrl = list(maxfun=1E4))) # Convergence failed


AICc.i<- AICctab(mod.i.null, mod.i.cat_edgeAsForest, mod.i.cat_edgeAsFire, mod.i.gam$mer, mod.i.gamcat_edgeAsForest$mer, nobs=nrow(dados.use), sort=F)
response.nb[i,1:5] <- AICc.i$dAICc


write.table(round(response.nb,2), file="resultados_AICc_nb.txt", quote=F, sep="\t")


# Selected models:
# 01 Abundance_mean - Cat_edgeAsFire
# 02 Abundance_mean_zoo - Cat_edgeAsFire
# 03 Abundance_mean_nzoo - Cat_edgeAsFire
# 04 Richness_mean  - Cat_edgeAsForest
# 05 Richness_mean_zoo - Cat_edgeAsForest
# 06 Richness_mean_nzoo - Cat_edgeAsForest
# 07 Richness_tot - Cat_edgeAsForest
# 08 Richness_tot_zoo - Cat_edgeAsForest
# 09 Richness_tot_nzoo - Cat_edgeAsForest
# 10 Width_all - Cat_edgeAsForest
# 11 Weight_all - Cat_edgeAsForest
# 12 Width_zoo - Cat_edgeAsFire
# 13 Weight_nzoo - Cat_edgeAsForest


# Making plots

# Figura 2 - Models
dist.test <- -150:150
Ndists <- length(dist.test)

mean.null <- rep(15, Ndists)
mean.cat <- ifelse(dist.test<=-10,16,20)
mean.cat2 <- ifelse(dist.test<=10,16,20)
mean.gam.orig <- 20/(1+exp(0.05*(-dist.test)))
gam.foo <- gam(mean.gam.orig ~ s(dist.test, fx=F, k=4))
mean.gam <- predict(gam.foo, type="response")
mean.gam.cat1 <- mean.gam + mean.cat*2
mean.gam.cat2 <- mean.gam + mean.cat2*2


var.null <- rep(3, Ndists)



fig1.width=26
fig1.height=13

#null.null

png(filename="SeedRainEdge_fig2_mods_a_null.png", height=fig1.height, width=fig1.width, unit="cm", res=300)
plot(mean.null~dist.test, xlab="Distance (m)", ylab="", xaxt="n", yaxt="n", type="l", ylim=range(c(mean.null+2*var.null, mean.null-2*var.null)), main="", cex.lab=1.8, lwd=4)
lines(I(mean.null+var.null) ~ dist.test, lty=2, lwd=4)
lines(I(mean.null-var.null) ~ dist.test, lty=2, lwd=4)
#axis(side=1, at=c(-150, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 150))
mtext(side=2,text="Seed rain", line=2, cex=1.8)
abline(v=0, col="gray", lwd=4)
dev.off()

#cat
png(filename="SeedRainEdge_fig2_mods_b_cat1.png", height=fig1.height, width=fig1.width, unit="cm", res=300)
plot(mean.cat~dist.test, xlab="Distance (m)", ylab="", xaxt="n", yaxt="n", type="n", ylim=range(c(mean.cat+2*var.null, mean.cat-2*var.null)), main="", cex.lab=1.8)
lines(I(mean.cat)[dist.test<=-10] ~ dist.test[dist.test<=-10], lty=1, lwd=4)
lines(I(mean.cat)[dist.test>-10] ~ dist.test[dist.test>-10], lty=1, lwd=4)
lines(I(mean.cat+var.null)[dist.test<=-10] ~ dist.test[dist.test<=-10], lty=2, lwd=4)
lines(I(mean.cat+var.null)[dist.test>-10] ~ dist.test[dist.test>-10], lty=2, lwd=4)
lines(I(mean.cat-var.null)[dist.test<=-10] ~ dist.test[dist.test<=-10], lty=2, lwd=4)
lines(I(mean.cat-var.null)[dist.test>-10] ~ dist.test[dist.test>-10], lty=2, lwd=4)
#axis(side=1, at=c(-150, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 150))
abline(v=0, col="gray", lwd=4)
mtext(side=2,text="Seed rain", line=2, cex=1.8)
dev.off()


png(filename="SeedRainEdge_fig2_mods_c_cat2.png", height=fig1.height, width=fig1.width, unit="cm", res=300)
plot(mean.cat2~dist.test, xlab="Distance (m)", ylab="", xaxt="n", yaxt="n", type="n", ylim=range(c(mean.cat2+2*var.null, mean.cat2-2*var.null)), main="", cex.lab=1.8)
lines(I(mean.cat2)[dist.test<=+10] ~ dist.test[dist.test<=+10], lty=1, lwd=4)
lines(I(mean.cat2)[dist.test>+10] ~ dist.test[dist.test>+10], lty=1, lwd=4)
lines(I(mean.cat2+var.null)[dist.test<=+10] ~ dist.test[dist.test<=+10], lty=2, lwd=4)
lines(I(mean.cat2+var.null)[dist.test>+10] ~ dist.test[dist.test>+10], lty=2, lwd=4)
lines(I(mean.cat2-var.null)[dist.test<=+10] ~ dist.test[dist.test<=+10], lty=2, lwd=4)
lines(I(mean.cat2-var.null)[dist.test>+10] ~ dist.test[dist.test>+10], lty=2, lwd=4)
#axis(side=1, at=c(-150, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 150))
abline(v=0, col="gray", lwd=4)
mtext(side=2,text="Seed rain", line=2, cex=1.8)
dev.off()


#gam.null
png(filename="SeedRainEdge_fig2_mods_d_gam.png", height=fig1.height, width=fig1.width, unit="cm", res=300)
plot(mean.gam~dist.test, xlab="Distance (m)", ylab="", xaxt="n", yaxt="n", type="l", ylim=range(c(mean.gam+2*var.null, mean.gam-2*var.null)), main="", cex.lab=1.8, lwd=4)
lines(I(mean.gam+var.null) ~ dist.test, lty=2, lwd=4)
lines(I(mean.gam+var.null) ~ dist.test, lty=2, lwd=4)
lines(I(mean.gam-var.null) ~ dist.test, lty=2, lwd=4)
lines(I(mean.gam-var.null) ~ dist.test, lty=2, lwd=4)
#axis(side=1, at=c(-150, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 150))
abline(v=0, col="gray", lwd=4)
mtext(side=2,text="Seed rain", line=2, cex=1.8)
dev.off()




#gam.cat
png(filename="SeedRainEdge_fig2_mods_e_gamcat.png", height=fig1.height, width=fig1.width, unit="cm", res=300)
plot(mean.gam.cat1~dist.test, xlab="Distance (m)", ylab="", xaxt="n", yaxt="n", type="n", ylim=range(c(mean.gam.cat1+2*var.null, mean.gam.cat1-2*var.null)), main="", cex.lab=1.8)
lines(I(mean.gam.cat1)[dist.test<=-10] ~ dist.test[dist.test<=-10], lty=1, lwd=4)
lines(I(mean.gam.cat1)[dist.test>-10] ~ dist.test[dist.test>-10], lty=1, lwd=4)
lines(I(mean.gam.cat1+var.null)[dist.test<=-10] ~ dist.test[dist.test<=-10], lty=2, lwd=4)
lines(I(mean.gam.cat1+var.null)[dist.test>-10] ~ dist.test[dist.test>-10], lty=2, lwd=4)
lines(I(mean.gam.cat1-var.null)[dist.test<=-10] ~ dist.test[dist.test<=-10], lty=2, lwd=4)
lines(I(mean.gam.cat1-var.null)[dist.test>-10] ~ dist.test[dist.test>-10], lty=2, lwd=4)
#axis(side=1, at=c(-150, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 150))
abline(v=0, col="gray", lwd=4)
mtext(side=2,text="Seed rain", line=2, cex=1.8)
dev.off()


#gam.cat2
png(filename="SeedRainEdge_fig2_mods_f_gamcat2.png", height=fig1.height, width=fig1.width, unit="cm", res=300)
plot(mean.gam.cat2~dist.test, xlab="Distance (m)", ylab="", xaxt="n", yaxt="n", type="n", ylim=range(c(mean.gam.cat2+2*var.null, mean.gam.cat2-2*var.null)), main="", cex.lab=1.8)
lines(I(mean.gam.cat2)[dist.test<=+10] ~ dist.test[dist.test<=+10], lty=1, lwd=4)
lines(I(mean.gam.cat2)[dist.test>+10] ~ dist.test[dist.test>+10], lty=1, lwd=4)
lines(I(mean.gam.cat2+var.null)[dist.test<=+10] ~ dist.test[dist.test<=+10], lty=2, lwd=4)
lines(I(mean.gam.cat2+var.null)[dist.test>+10] ~ dist.test[dist.test>+10], lty=2, lwd=4)
lines(I(mean.gam.cat2-var.null)[dist.test<=+10] ~ dist.test[dist.test<=+10], lty=2, lwd=4)
lines(I(mean.gam.cat2-var.null)[dist.test>+10] ~ dist.test[dist.test>+10], lty=2, lwd=4)
#axis(side=1, at=c(-150, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 150))
abline(v=0, col="gray", lwd=2)
mtext(side=2,text="Seed rain", line=2, cex=1.8)
dev.off()





### Making the figures...

xcoord1.edgeAsFire <- c(-150,5)
xcoord2.edgeAsFire <- c(5,150)
newdata.edgeAsFire <- matrix(ncol=2, nrow=14)
newdata.edgeAsFire[,1] <- sort(c(1:7, 1:7))
newdata.edgeAsFire[,2] <- rep(c("Fire","Forest"),7)
colnames(newdata.edgeAsFire) <- c("Transect", "Vegetation_edgeAsFire")
newdata.edgeAsFire <- apply(newdata.edgeAsFire, c(1,2), as.factor)
newdata.edgeAsFire <- as.data.frame(newdata.edgeAsFire)

xcoord1.edgeAsForest <- c(-150,-5)
xcoord2.edgeAsForest <- c(-5,150)
newdata.edgeAsForest <- matrix(ncol=2, nrow=14)
newdata.edgeAsForest[,1] <- sort(c(1:7, 1:7))
newdata.edgeAsForest[,2] <- rep(c("Fire","Forest"),7)
colnames(newdata.edgeAsForest) <- c("Transect", "Vegetation_edgeAsForest")
newdata.edgeAsForest <- apply(newdata.edgeAsForest, c(1,2), as.factor)
newdata.edgeAsForest <- as.data.frame(newdata.edgeAsForest)


### Figure 3 - Abundance
# Selected models:
# 01 Abundance_mean - Cat_edgeAsFire
# 02 Abundance_mean_zoo - Cat_edgeAsFire
# 03 Abundance_mean_nzoo - Cat_edgeAsFire


mod.final.Abundance_tot <- glmer.nb(Abundance_tot ~ Vegetation_edgeAsFire + (1|Transect), data=dados, method="ML", control=glmerControl(tolPwrss=0.00001, optCtrl=list(maxfun=1E4)))	
mod.final.Abundance_tot_zoo <- glmer.nb(Abundance_tot_zoo ~ Vegetation_edgeAsFire + (1|Transect), data=dados, method="ML", control=glmerControl(tolPwrss=0.00001, optCtrl=list(maxfun=1E4)))
mod.final.Abundance_tot_nzoo <- glmer.nb(Abundance_tot_nzoo ~ Vegetation_edgeAsFire + (1|Transect), data=dados, method="ML", control=glmerControl(tolPwrss=0.00001, optCtrl=list(maxfun=1E4)))


pred.Abundance_tot <- predict(mod.final.Abundance_tot, newdata=newdata.edgeAsFire, type="response")
mean1.Abundance_tot <- mean(pred.Abundance_tot[c(T,F)])
mean2.Abundance_tot <- mean(pred.Abundance_tot[c(F,T)])

pred.Abundance_tot_zoo <- predict(mod.final.Abundance_tot_zoo, newdata=newdata.edgeAsFire, type="response")
mean1.Abundance_tot_zoo <- mean(pred.Abundance_tot_zoo[c(T,F)])
mean2.Abundance_tot_zoo <- mean(pred.Abundance_tot_zoo[c(F,T)])

pred.Abundance_tot_nzoo <- predict(mod.final.Abundance_tot_nzoo, newdata=newdata.edgeAsFire, type="response")
mean1.Abundance_tot_nzoo <- mean(pred.Abundance_tot_nzoo[c(T,F)])
mean2.Abundance_tot_nzoo <- mean(pred.Abundance_tot_nzoo[c(F,T)])

png(filename="SeedRainEdge_fig4_abundance.png", height=20, width=10, unit="cm", res=300)
par(mfrow=c(3,1), mar=c(2,2,2,2), oma=c(4,4,2,2))
plot(Abundance_tot ~ Distance, data=dados, main="a. Abundance - all seeds", ylim=range(dados$Abundance_tot))
abline(v=0, col="gray40")
segments(x0=xcoord1.edgeAsFire[1], x1=xcoord1.edgeAsFire[2], y0=mean1.Abundance_tot, y1=mean1.Abundance_tot)
segments(x0=xcoord2.edgeAsFire[1], x1=xcoord2.edgeAsFire[2], y0=mean2.Abundance_tot, y1=mean2.Abundance_tot)
plot(Abundance_tot_zoo ~ Distance, data=dados, main="b. Abundance - zoochoric", ylim=range(dados$Abundance_tot))
abline(v=0, col="gray40")
segments(x0=xcoord1.edgeAsFire[1], x1=xcoord1.edgeAsFire[2], y0=mean1.Abundance_tot_zoo, y1=mean1.Abundance_tot_zoo)
segments(x0=xcoord2.edgeAsFire[1], x1=xcoord2.edgeAsFire[2], y0=mean2.Abundance_tot_zoo, y1=mean2.Abundance_tot_zoo)
plot(Abundance_tot_nzoo ~ Distance, data=dados, main="c. Abundance - non-zoochoric", ylim=range(dados$Abundance_tot))
abline(v=0, col="gray40")
segments(x0=xcoord1.edgeAsFire[1], x1=xcoord1.edgeAsFire[2], y0=mean1.Abundance_tot_nzoo, y1=mean1.Abundance_tot_nzoo)
segments(x0=xcoord2.edgeAsFire[1], x1=xcoord2.edgeAsFire[2], y0=mean2.Abundance_tot_nzoo, y1=mean2.Abundance_tot_nzoo)
mtext(side=1, text="Distance along transect (m)", outer=T, line=1)
mtext(side=2, text="Number of seeds", outer=T, line=1)
dev.off()

# 04 Richness_mean  - Cat_edgeAsFire
# 05 Richness_mean_zoo - Cat_edgeAsFire
# 06 Richness_mean_nzoo - Cat_edgeAsFire
# 07 Richness_tot - Cat_edgeAsForest
# 08 Richness_tot_zoo - Cat_edgeAsForest
# 09 Richness_tot_nzoo - Cat_edgeAsForest


mod.final.Richness_mean <- lmer(Richness_mean ~ Vegetation_edgeAsFire + (1|Transect), data=dados.use, REML=F)
mod.final.Richness_mean_zoo <- lmer(Richness_mean_zoo ~ Vegetation_edgeAsFire + (1|Transect), data=dados.use, REML=F)
mod.final.Richness_mean_nzoo <- lmer(Richness_mean_nzoo ~ Vegetation_edgeAsFire + (1|Transect), data=dados.use, REML=F)
		
mod.final.Richness_tot <- glmer.nb(Richness_tot ~ Vegetation_edgeAsForest + (1|Transect), data=dados, method="ML", control=glmerControl(tolPwrss=0.00001, optCtrl=list(maxfun=1E4)))	
mod.final.Richness_tot_zoo <- glmer.nb(Richness_tot_zoo ~ Vegetation_edgeAsForest + (1|Transect), data=dados, method="ML", control=glmerControl(tolPwrss=0.00001, optCtrl=list(maxfun=1E4)))
mod.final.Richness_tot_nzoo <- glmer.nb(Richness_tot_nzoo ~ Vegetation_edgeAsForest + (1|Transect), data=dados, method="ML", control=glmerControl(tolPwrss=0.00001, optCtrl=list(maxfun=1E4)))


pred.Richness_mean <- predict(mod.final.Richness_mean, newdata=newdata.edgeAsFire, type="response")
mean1.Richness_mean <- mean(pred.Richness_mean[c(T,F)])
mean2.Richness_mean <- mean(pred.Richness_mean[c(F,T)])
pred.Richness_mean_zoo <- predict(mod.final.Richness_mean_zoo, newdata=newdata.edgeAsFire, type="response")
mean1.Richness_mean_zoo <- mean(pred.Richness_mean_zoo[c(T,F)])
mean2.Richness_mean_zoo <- mean(pred.Richness_mean_zoo[c(F,T)])
pred.Richness_mean_nzoo <- predict(mod.final.Richness_mean_nzoo, newdata=newdata.edgeAsFire, type="response")
mean1.Richness_mean_nzoo <- mean(pred.Richness_mean_nzoo[c(T,F)])
mean2.Richness_mean_nzoo <- mean(pred.Richness_mean_nzoo[c(F,T)])

pred.Richness_tot <- predict(mod.final.Richness_tot, newdata=newdata.edgeAsForest, type="response")
mean1.Richness_tot <- mean(pred.Richness_tot[c(T,F)])
mean2.Richness_tot <- mean(pred.Richness_tot[c(F,T)])
pred.Richness_tot_zoo <- predict(mod.final.Richness_tot_zoo, newdata=newdata.edgeAsForest, type="response")
mean1.Richness_tot_zoo <- mean(pred.Richness_tot_zoo[c(T,F)])
mean2.Richness_tot_zoo <- mean(pred.Richness_tot_zoo[c(F,T)])
pred.Richness_tot_nzoo <- predict(mod.final.Richness_tot_nzoo, newdata=newdata.edgeAsForest, type="response")
mean1.Richness_tot_nzoo <- mean(pred.Richness_tot_nzoo[c(T,F)])
mean2.Richness_tot_nzoo <- mean(pred.Richness_tot_nzoo[c(F,T)])

png(filename="SeedRainEdge_fig5_richness.png", height=20, width=15, unit="cm",res=300)

par(mfcol=c(3,2), mar=c(3,3,2,2), oma=c(3,3,2,2))

plot(Richness_mean ~ Distance, data=dados, main="a. Mean richness - all seeds", ylim=range(dados$Richness_mean))
abline(v=0, col="gray40")
segments(x0=xcoord1.edgeAsFire[1], x1=xcoord1.edgeAsFire[2], y0=mean1.Richness_mean, y1=mean1.Richness_mean)
segments(x0=xcoord2.edgeAsFire[1], x1=xcoord2.edgeAsFire[2], y0=mean2.Richness_mean, y1=mean2.Richness_mean)

plot(Richness_mean_zoo ~ Distance, data=dados, main="b. Mean richness - zoochoric", ylim=range(dados$Richness_mean))
abline(v=0, col="gray40")
segments(x0=xcoord1.edgeAsFire[1], x1=xcoord1.edgeAsFire[2], y0=mean1.Richness_mean_zoo, y1=mean1.Richness_mean_zoo)
segments(x0=xcoord2.edgeAsFire[1], x1=xcoord2.edgeAsFire[2], y0=mean2.Richness_mean_zoo, y1=mean2.Richness_mean_zoo)


plot(Richness_mean_nzoo ~ Distance, data=dados, main="c. Mean richness - non-zoochoric", ylim=range(dados$Richness_mean))
abline(v=0, col="gray40")
segments(x0=xcoord1.edgeAsFire[1], x1=xcoord1.edgeAsFire[2], y0=mean1.Richness_mean_nzoo, y1=mean1.Richness_mean_nzoo)
segments(x0=xcoord2.edgeAsFire[1], x1=xcoord2.edgeAsFire[2], y0=mean2.Richness_mean_nzoo, y1=mean2.Richness_mean_nzoo)


plot(Richness_tot ~ Distance, data=dados, main="d. Total richness - all seeds", ylim=range(dados$Richness_tot))
abline(v=0, col="gray40")
segments(x0=xcoord1.edgeAsForest[1], x1=xcoord1.edgeAsForest[2], y0=mean1.Richness_tot, y1=mean1.Richness_tot)
segments(x0=xcoord2.edgeAsForest[1], x1=xcoord2.edgeAsForest[2], y0=mean2.Richness_tot, y1=mean2.Richness_tot)

plot(Richness_tot_zoo ~ Distance, data=dados, main="e. Total richness - zoochoric", ylim=range(dados$Richness_tot))
abline(v=0, col="gray40")
segments(x0=xcoord1.edgeAsForest[1], x1=xcoord1.edgeAsForest[2], y0=mean1.Richness_tot_zoo, y1=mean1.Richness_tot_zoo)
segments(x0=xcoord2.edgeAsForest[1], x1=xcoord2.edgeAsForest[2], y0=mean2.Richness_tot_zoo, y1=mean2.Richness_tot_zoo)


plot(Richness_tot_nzoo ~ Distance, data=dados, main="f. Total richness - non-zoochoric", ylim=range(dados$Richness_tot))
abline(v=0, col="gray40")
segments(x0=xcoord1.edgeAsForest[1], x1=xcoord1.edgeAsForest[2], y0=mean1.Richness_tot_nzoo, y1=mean1.Richness_tot_nzoo)
segments(x0=xcoord2.edgeAsForest[1], x1=xcoord2.edgeAsForest[2], y0=mean2.Richness_tot_nzoo, y1=mean2.Richness_tot_nzoo)

mtext(side=1, text="Distance along transect (m)", outer=T, line=1)
mtext(side=2, text="Number of species", outer=T, line=1)

dev.off()

# 10 Width_all - Cat_edgeAsForest
# 11 Weight_all - Cat_edgeAsForest
# 12 Width_zoo - Cat_edgeAsFire
# 13 Weight_nzoo - Cat_edgeAsForest

mod.final.Width_all <- glmer(Width_all ~ Vegetation_edgeAsForest + (1|Transect), data=dados.use, method="ML", control=glmerControl(tolPwrss=1E-6, optCtrl = list(maxfun=1E4)), family=Gamma)
mod.final.Weight_all <- glmer(Weight_all ~ Vegetation_edgeAsFire + (1|Transect), data=dados.use, method="ML", control=glmerControl(tolPwrss=1E-6, optCtrl = list(maxfun=1E4)), family=Gamma)
mod.final.Width_zoo <- glmer(Width_zoo ~ Vegetation_edgeAsForest + (1|Transect), data=dados.use, method="ML", control=glmerControl(tolPwrss=1E-6, optCtrl = list(maxfun=1E4)), family=Gamma)
mod.final.Weight_nzoo <- glmer(Weight_nzoo ~ Vegetation_edgeAsForest + (1|Transect), data=dados.use, method="ML", control=glmerControl(tolPwrss=1E-6, optCtrl = list(maxfun=1E4)), family=Gamma)

pred.Width_all <- predict(mod.final.Width_all, newdata=newdata.edgeAsForest, type="response")
mean1.Width_all <- mean(pred.Width_all[c(T,F)])
mean2.Width_all <- mean(pred.Width_all[c(F,T)])

pred.Weight_all <- predict(mod.final.Weight_all, newdata=newdata.edgeAsFire, type="response")
mean1.Weight_all <- mean(pred.Weight_all[c(T,F)])
mean2.Weight_all <- mean(pred.Weight_all[c(F,T)])

pred.Width_zoo <- predict(mod.final.Width_zoo, newdata=newdata.edgeAsForest, type="response")
mean1.Width_zoo <- mean(pred.Width_zoo[c(T,F)])
mean2.Width_zoo <- mean(pred.Width_zoo[c(F,T)])

pred.Weight_nzoo <- predict(mod.final.Weight_nzoo, newdata=newdata.edgeAsForest, type="response")
mean1.Weight_nzoo <- mean(pred.Weight_nzoo[c(T,F)])
mean2.Weight_nzoo <- mean(pred.Weight_nzoo[c(F,T)])

png(filename="SeedRainEdge_fig6_size.png", height=20, width=20, unit="cm", res=300)

par(mfrow=c(2,2), mar=c(3,3,2,2), oma=c(3,3,2,2))

plot(Width_all ~ Distance, data=dados, main="a. Mean width - all seeds", ylim=range(c(dados$Width_all, dados$Width_zoo), na.rm=T), ylab="")
abline(v=0, col="gray40")
segments(x0=xcoord1.edgeAsForest[1], x1=xcoord1.edgeAsForest[2], y0=mean1.Width_all, y1=mean1.Width_all)
segments(x0=xcoord2.edgeAsForest[1], x1=xcoord2.edgeAsForest[2], y0=mean2.Width_all, y1=mean2.Width_all)

plot(Width_zoo ~ Distance, data=dados, main="b. Mean width - zoochoric", ylim=range(c(dados$Width_all, dados$Width_zoo), na.rm=T))
abline(v=0, col="gray40")
segments(x0=xcoord1.edgeAsFire[1], x1=xcoord1.edgeAsFire[2], y0=mean1.Width_zoo, y1=mean1.Width_zoo)
segments(x0=xcoord2.edgeAsFire[1], x1=xcoord2.edgeAsFire[2], y0=mean2.Width_zoo, y1=mean2.Width_zoo)


plot(Weight_all ~ Distance, data=dados, main="c. Mean weight - all seeds", ylim=range(c(dados$Weight_all, dados$Weight_nzoo), na.rm=T), ylab="")
abline(v=0, col="gray40")
segments(x0=xcoord1.edgeAsForest[1], x1=xcoord1.edgeAsForest[2], y0=mean1.Weight_all, y1=mean1.Weight_all)
segments(x0=xcoord2.edgeAsForest[1], x1=xcoord2.edgeAsForest[2], y0=mean2.Weight_all, y1=mean2.Weight_all)

plot(Weight_nzoo ~ Distance, data=dados, main="d. Mean weight - non-zoochoric", ylim=range(c(dados$Weight_all, dados$Weight_nzoo), na.rm=T))
abline(v=0, col="gray40")
segments(x0=xcoord1.edgeAsForest[1], x1=xcoord1.edgeAsForest[2], y0=mean1.Weight_nzoo, y1=mean1.Weight_nzoo)
segments(x0=xcoord2.edgeAsForest[1], x1=xcoord2.edgeAsForest[2], y0=mean2.Weight_nzoo, y1=mean2.Weight_nzoo)

mtext(side=1, text="Distance along transect (m)", outer=T, line=1)
mtext(side=2, text=c("Weight (g)", "Width (mm)"), at=c(0.24,0.75), outer=T, line=0.5)

dev.off()


### Calculate mean values per side (forest VS fire)...

dados.fire <- subset(dados, Distance < 0)
dados.forest <- subset(dados, Distance > 0)

means.fire <- apply(dados.fire[,-c(1,2,16,17)], 2, mean, na.rm=T)
means.forest <- apply(dados.forest[,-c(1,2,16,17)], 2, mean, na.rm=T)

means.ratio <- means.forest/means.fire

### Make histograms for seed size...

dados.size <- read.table("dados_originais_R.txt", header=T)

dados.size$Tamanho <- apply(dados.size[,c("Altura","Largura")],1,min)

foo.peso <- hist(dados.size$Peso, breaks=42, plot=F)
foo.tamanho <- hist(dados.size$Tamanho, breaks=42, plot=F)
foo.peso.nzoo <- hist(subset(dados.size, Sindrome==0)$Peso, breaks=foo.peso$breaks, plot=F)
foo.tamanho.zoo <- hist(subset(dados.size, Sindrome==1)$Tamanho, breaks=foo.tamanho$breaks, plot=F)


foo.peso$counts <- log10(foo.peso$counts+1) 
foo.tamanho$counts <- log10(foo.tamanho$counts+1)
foo.peso.nzoo$counts <- log10(foo.peso.nzoo$counts+1)
foo.tamanho.zoo$counts <- log10(foo.tamanho.zoo$counts+1)


png(filename="SeedRainEdge_fig3_sizeHistogams.png", height=20, width=20, unit="cm", res=300)
par(mfrow=c(2,2), mar=c(4,3,2,2), oma=c(1,3,1,1))

plot(foo.tamanho, main="a. Seed width", xlab="Width (mm)", ylab="", yaxt="n", xlim=range(dados.size$Tamanho))
axis(side=2, at=log10(c(0,1,10,100,1000,10000)+1), labels=c(0,1,10,100,1000,10000), las=1)
plot(foo.tamanho.zoo, main="b. Seed width - zoochoric", xlab="Width (mm)", ylab="", yaxt="n", xlim=range(dados.size$Tamanho))
axis(side=2, at=log10(c(0,1,10,100,1000,10000)+1), labels=c(0,1,10,100,1000,10000), las=1)
plot(foo.peso, main="c. Seed weight", xlab="Weight (g)", ylab="", yaxt="n", xlim=range(dados.size$Peso))
axis(side=2, at=log10(c(0,1,10,100,1000,10000)+1), labels=c(0,1,10,100,1000,10000), las=1)
plot(foo.peso.nzoo, main="d. Seed weight - non-zoochoric", xlab="Weight (g)", ylab="", yaxt="n", xlim=range(dados.size$Peso))
axis(side=2, at=log10(c(0,1,10,100,1000,10000)+1), labels=c(0,1,10,100,1000,10000), las=1)

mtext(side=2, text="Number of seeds", outer=T, line=1)
dev.off()



