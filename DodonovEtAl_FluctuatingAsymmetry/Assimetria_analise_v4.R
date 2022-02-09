panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, use="complete.obs")^2
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste(prefix, txt, sep = "")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}


setwd("/home/pavel/Profissional/Pesquisa/MyPapers_working/Pavel_2011_Asymmetry")

data.inds=read.table("assimetria_dados_ind.txt", sep="", dec=".", header=T, stringsAsFactors=TRUE)
str(data.inds)
pairs(data.inds[data.inds$Medidor=="Andreza",-c(1,2)], lower.panel=panel.cor)

#Hmedia e Hmax altamente correlacionadas - vamos usar apenas Hmax. (correlações abaixo de 0.7)

###Análise stepwise, separadamante pra cada medidor e pra AFtot e AFrel

#FAtot

result.FAtot.Andreza=step(lm(AFtot_abs~Dens+ABtotal+ABmedia+Hmax, 
	data=subset(data.inds, Medidor=="Andreza")), dir="both", trace=5)


result.FAtot.Luiz=step(lm(AFtot_abs~Dens+ABtotal+ABmedia+Hmax, 
	data=subset(data.inds, Medidor=="Luiz")), dir="both", trace=5)


result.FAtot.Pavel=step(lm(AFtot_abs~Dens+ABtotal+ABmedia+Hmax, 
	data=subset(data.inds, Medidor=="Pavel")), dir="both", trace=5)


#FArel


result.FArel.Andreza=step(lm(AFrel_abs~Dens+ABtotal+ABmedia+Hmax, 
	data=subset(data.inds, Medidor=="Andreza")), dir="both", trace=5)


result.FArel.Luiz=step(lm(AFrel_abs~Dens+ABtotal+ABmedia+Hmax, 
	data=subset(data.inds, Medidor=="Luiz")), dir="both", trace=5)


result.FArel.Pavel=step(lm(AFrel_abs~Dens+ABtotal+ABmedia+Hmax, 
	data=subset(data.inds, Medidor=="Pavel")), dir="both", trace=5)


summary(result.FAtot.Andreza)
summary(result.FAtot.Luiz)
summary(result.FAtot.Pavel)
summary(result.FArel.Andreza)
summary(result.FArel.Luiz)
summary(result.FArel.Pavel)


###Por folha

library(nlme)
library(mgcv)
library(quantreg)

setwd("/home/pavel/Profissional/Pesquisa/MyPapers_working/Pavel_2011_Asymmetry")

dados.folha=read.table("assimetria_dados_folha.txt", header=T, sep="", dec=".", stringsAsFactors=TRUE)
plot(FAtot_abs~Atot, data=dados.folha, col=Medidor)
plot(FArel~Atot, data=dados.folha, col=Medidor)

result.folha.Andreza=lme(FArel_abs~Atot, random=~1|as.factor(Planta), data=subset(dados.folha, Medidor=="Andreza"))
result.folha.Luiz=lme(FArel_abs~Atot, random=~1|as.factor(Planta), data=subset(dados.folha, Medidor=="Luiz"))
result.folha.Pavel=lme(FArel_abs~Atot, random=~1|as.factor(Planta), data=subset(dados.folha, Medidor=="Pavel"))

summary(lm(FArel_abs~Atot, data=subset(dados.folha, Medidor=="Andreza")))

summary(result.folha.Andreza)
summary(result.folha.Luiz)
summary(result.folha.Pavel)
png(filename="assimetria_fig2.png", res=300, height=15, width=15, unit="cm")
plot(FArel_abs~Atot, data=subset(dados.folha, Medidor=="Andreza"), 
	xlab=expression(paste("Total leaf area (cm"^"2",")",sep="")), ylab="Absolute FArel")
abline(coef=c(0.05910578, -0.00023007), col="grey", lwd=2)
dev.off()





result.folha.Andreza=lme(FArel_abs~Atot, random=~1|as.factor(Planta), data=subset(dados.folha, Medidor=="Andreza"))
result.folha.Luiz=lme(FArel_abs~Atot, random=~1|as.factor(Planta), data=subset(dados.folha, Medidor=="Luiz"))
result.folha.Pavel=lme(FArel_abs~Atot, random=~1|as.factor(Planta), data=subset(dados.folha, Medidor=="Pavel"))

summary(lm(FArel_abs~Atot, data=subset(dados.folha, Medidor=="Pavel")))

summary(result.folha.Andreza)
summary(result.folha.Luiz)
summary(result.folha.Pavel)




###Quantilic regression
library(quantreg)

# Function for quantile regression with restricted randomization, to assess significance while accounting for non-independnce of the leaves on the same plant

#y = subset(dados.folha, Medidor=="Andreza")$FArel_abs
#x = subset(dados.folha, Medidor=="Andreza")$Atot
#block = subset(dados.folha, Medidor=="Andreza")$Planta
#quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9)

rq.restrand <- function(y, x, quantiles, block, Nperm=5000) {
	quantreg.orig <- rq(y ~ x, tau=quantiles)
	coefs.orig <- quantreg.orig$coefficients[2,]
	levels.block <- unique(block)
	Nlevels <- length(levels.block)
	coefs.perm <- matrix(nrow=Nperm, ncol=length(quantiles))
	coefs.perm[1,] <- coefs.orig
	for(i in 2:Nperm) {
		y.perm <- numeric(length(y))
		for(j in 1:Nlevels) {
			y.now <- y[block==levels.block[j]]
			y.now <- sample(y.now)
			if(j==1) y.perm <- y.now else y.perm <- c(y.perm, y.now)
		}
		quantreg.now <- rq(y.perm ~ x, tau=quantiles)
		coefs.now <- quantreg.now$coefficients[2,]
		coefs.perm[i,] <- coefs.now
		print(i)
	}
	calculate.signif <- function(x) sum(abs(x) >= abs(x[1]))/length(x)
	significance <- apply (coefs.perm, 2, calculate.signif)
	result <- data.frame(quantiles = quantiles, coef = coefs.orig, signif = significance)
	return(result)
}

#rm(x)
#rm(y)
#rm(quantiles)
#rm(block)

quantreg.restrand.Andreza <- rq.restrand(y=subset(dados.folha, Medidor=="Andreza")$FArel_abs, x=subset(dados.folha, Medidor=="Andreza")$Atot, block=subset(dados.folha, Medidor=="Andreza")$Planta, quantiles=c(0.1, 0.25, 0.5, 0.75, 0.9), Nperm=5000)
quantreg.restrand.Luiz <- rq.restrand(y=subset(dados.folha, Medidor=="Luiz")$FArel_abs, x=subset(dados.folha, Medidor=="Luiz")$Atot, block=subset(dados.folha, Medidor=="Luiz")$Planta, quantiles=c(0.1, 0.25, 0.5, 0.75, 0.9), Nperm=5000)
quantreg.restrand.Pavel <- rq.restrand(y=subset(dados.folha, Medidor=="Pavel")$FArel_abs, x=subset(dados.folha, Medidor=="Pavel")$Atot, block=subset(dados.folha, Medidor=="Pavel")$Planta, quantiles=c(0.1, 0.25, 0.5, 0.75, 0.9), Nperm=5000)



quantreg.Andreza=rq(FArel_abs~Atot, tau=c(0.1, 0.25, 0.5, 0.75, 0.9), data=subset(dados.folha, Medidor=="Andreza"))
quantreg.Luiz=rq(FArel_abs~Atot, tau=c(0.1, 0.25, 0.5, 0.75, 0.9), data=subset(dados.folha, Medidor=="Luiz"))
quantreg.Pavel=rq(FArel_abs~Atot, tau=c(0.1, 0.25, 0.5, 0.75, 0.9), data=subset(dados.folha, Medidor=="Pavel"))


summary(quantreg.Andreza)
quantreg.restrand.Andreza

summary(quantreg.Luiz)
quantreg.restrand.Luiz

summary(quantreg.Pavel)
quantreg.restrand.Pavel




png(filename="assimetria_fig2.png", res=300, height=15, width=15, unit="cm")
plot(FArel_abs~Atot, data=subset(dados.folha, Medidor=="Andreza"), 
	xlab=expression(paste("Total leaf area (cm"^"2",")",sep="")), ylab="Absolute FArel")
abline(coef=c(0.05910578, -0.00023007), col="grey20", lwd=2)
abline(rq(FArel_abs~Atot, tau=0.90, data=subset(dados.folha, Medidor=="Andreza")), lty=2, lwd=2, col="grey20")
dev.off()


###Agora uma abordagem conservadora, tentando reduzir o erro de medida

dados.folha.corrected=aggregate(dados.folha[,c("FAtot","Atot")], by=list(Planta=dados.folha$Planta, Folha=dados.folha$Folha), mean)
dados.folha.corrected$FArel=dados.folha.corrected$FAtot/dados.folha.corrected$Atot
dados.folha.corrected$FArel_abs=abs(dados.folha.corrected$FArel)

lme.folha.corrected=lme(FArel_abs~Atot, random=~1|as.factor(Planta), data=dados.folha.corrected)
summary(lme.folha.corrected)
quantreg.corrected=rq(FArel_abs~Atot, tau=c(0.1, 0.25, 0.5, 0.75, 0.9), dados.folha.corrected)
summary(quantreg.corrected)
quantreg.restrand.mean <- rq.restrand(y=dados.folha.corrected$FArel_abs, x=dados.folha.corrected$Atot, block=dados.folha.corrected$Planta, quantiles=c(0.1, 0.25, 0.5, 0.75, 0.9), Nperm=5000)
quantreg.restrand.mean



png(filename="assimetria_fig2.png", res=300, height=15, width=15, unit="cm")
plot(FArel_abs~Atot, data=dados.folha.corrected, 
	xlab=expression(paste("Total leaf area (cm"^"2",")",sep="")), ylab="Absolute FArel")
abline(coef=c(0.05722946, -0.00019177), col="grey20", lwd=2)
abline(rq(FArel_abs~Atot, tau=0.90, data=dados.folha.corrected), lty=2, lwd=2, col="grey20")
dev.off()

dados.folha.corrected$FAtot_abs=abs(dados.folha.corrected$FAtot)

dados.ind.mean=aggregate(dados.folha.corrected[,c("FAtot_abs","FArel_abs")], by=list(Planta=dados.folha.corrected$Planta), mean)

dados.ind.mean=subset(dados.ind.mean,Planta!=62)


dados.ind=data.inds[1:69,]

dados.ind.corrected=cbind(dados.ind.mean, dados.ind)


result.FAtot.mean=step(lm(FAtot_abs~Dens+ABtotal+ABmedia+Hmax, 
	data=dados.ind.corrected), dir="both", trace=5)

summary(result.FAtot.mean)


result.FArel.mean=step(lm(FArel_abs~Dens+ABtotal+ABmedia+Hmax, 
	data=dados.ind.corrected), dir="both", trace=5)

summary(result.FArel.mean)


###Calculando o intervalo de confiança para AFtot (não a absoluta!)

bootstrapT <- function (x, n=10000, alpha=0.05) { # requires a numeric vector as input
	#Critical values are calculated using bootstrap-t (Manly 1998, p. 56)
	Teta.boot=numeric(0)
	SEM.boot=numeric(0)
	Tb=numeric(0) 
	Teta.est=mean(x)
	SEM.est=sd(x)/sqrt(length(x))
	for (i in 1:n) {
		x.boot=sample(x,size=length(x),replace=TRUE,)
		Tb[i]=(mean(x.boot)-Teta.est)/(sd(x.boot)/sqrt(length(x)))
		}
	hist(Tb)
	#Calculating the alpha critical T values
	Tsup=quantile(Tb,alpha/2,na.rm=TRUE)
	Tinf=quantile(Tb,1-alpha/2,na.rm=TRUE)
	lim.inf=Teta.est-Tinf*SEM.est
	lim.sup=Teta.est-Tsup*SEM.est
	c(lim.inf,lim.sup)
	}



dados.folha=read.table("assimetria_dados_folha.txt", header=T, sep="", dec=".")

FAtot.media=aggregate(dados.folha$FAtot, by=list(Planta=dados.folha$Planta, Folha=dados.folha$Folha), mean)
names(FAtot.media)[3]="FAtot"

CI.Andreza=bootstrapT(subset(dados.folha, Medidor=="Andreza")$FAtot)
CI.Luiz=bootstrapT(subset(dados.folha, Medidor=="Luiz")$FAtot)
CI.Pavel=bootstrapT(subset(dados.folha, Medidor=="Pavel")$FAtot)
CI.media=bootstrapT(FAtot.media$FAtot)

rbind(CI.Andreza, CI.Luiz, CI.Pavel, CI.media)

write.table(FAtot.media, file="FAtot.media", sep="\t")

png(filename="assimetria_fig2.png", res=300, height=20, width=20, unit="cm")
hist(subset(dados.folha, Medidor=="Andreza")$FAtot, freq=F, breaks=50,
	xlab=expression(paste("Right area - left area (cm"^"2",")",sep="")), axes=F, main="")
axis(side=1, at=c(-12.5,-10,-7.5,-5,-2.5,0,2.5,5,7.5,10,12.5), line=-0.8)
axis(side=2, at=c(0,0.05,0.10,0.15,0.20,0.25), line=-0.8)

curve(dnorm(x,mean=0.23843, sd=2.8356)*0.61112, add=T, lwd=1.5) 
curve(dnorm(x,mean=0.099606, sd=1.159)*0.3888, add=T, lwd=1.5) 

dev.off()

