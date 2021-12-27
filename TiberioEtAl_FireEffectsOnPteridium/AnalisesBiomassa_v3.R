setwd("/home/pavel/Profissional/Pesquisa/MyPapers-Revising/PteridiumBrasilia")

library(yarrr)

dados=read.table("dados_biomassa.txt", sep="", header=T)
Coletas=unique(dados$Coleta)
Areas=unique(dados$Area)

dados$Coleta.jitter=dados$Coleta+runif(160,-0.2,0.2)

Nperm=5000



plot(Densidade~Coleta.jitter, col=Area, data=dados)
plot(Altura~Coleta.jitter, col=Area, data=dados)
plot(Biomassa~Coleta.jitter, col=Area, data=dados)
plot(Serrapilheira~Coleta.jitter, col=Area, data=dados)


# Diferença entre os anos na área queimada:

dados.area1=subset(dados,Area==1)

dens.coleta.area1=aov(Densidade~Coleta, data=dados.area1)
dens.coleta.area1.orig=summary(dens.coleta.area1)[[1]]$"F value"[1]
dens.coleta.area1.perm=numeric(Nperm)
dens.coleta.area1.perm[1]=dens.coleta.area1.orig

alt.coleta.area1=aov(Altura~Coleta, data=dados.area1)
alt.coleta.area1.orig=summary(alt.coleta.area1)[[1]]$"F value"[1]
alt.coleta.area1.perm=numeric(Nperm)
alt.coleta.area1.perm[1]=alt.coleta.area1.orig

biom.coleta.area1=aov(Biomassa~Coleta, data=dados.area1)
biom.coleta.area1.orig=summary(biom.coleta.area1)[[1]]$"F value"[1]
biom.coleta.area1.perm=numeric(Nperm)
biom.coleta.area1.perm[1]=biom.coleta.area1.orig

serap.coleta.area1=aov(Serrapilheira~Coleta, data=dados.area1)
serap.coleta.area1.orig=summary(serap.coleta.area1)[[1]]$"F value"[1]
serap.coleta.area1.perm=numeric(Nperm)
serap.coleta.area1.perm[1]=serap.coleta.area1.orig


for(i in 2:Nperm) {
	dados.temp=dados.area1
	dados.temp$Coleta=sample(dados.temp$Coleta)
	
	dens.foo=aov(Densidade~Coleta, data=dados.temp)
	alt.foo=aov(Altura~Coleta, data=dados.temp)
	biom.foo=aov(Biomassa~Coleta, data=dados.temp)
	serap.foo=aov(Serrapilheira~Coleta, data=dados.temp)

	
	dens.coleta.area1.perm[i]=summary(dens.foo)[[1]]$"F value"[1]
	alt.coleta.area1.perm[i]=summary(alt.foo)[[1]]$"F value"[1]
	biom.coleta.area1.perm[i]=summary(biom.foo)[[1]]$"F value"[1]
	serap.coleta.area1.perm[i]=summary(serap.foo)[[1]]$"F value"[1]
	if (i %% 100 == 0) print (i)

	}



signif.dens.coleta.area1=sum(dens.coleta.area1.perm >= dens.coleta.area1.orig) / Nperm		
signif.alt.coleta.area1=sum(alt.coleta.area1.perm >= alt.coleta.area1.orig) / Nperm		
signif.biom.coleta.area1=sum(biom.coleta.area1.perm >= biom.coleta.area1.orig) / Nperm		
signif.serap.coleta.area1=sum(serap.coleta.area1.perm >= serap.coleta.area1.orig) / Nperm		

signif.dens.coleta.area1
signif.alt.coleta.area1
signif.biom.coleta.area1
signif.serap.coleta.area1




### Diferença entre os anos na área controle:


dados.area2=subset(dados,Area==2)

dens.coleta.area2=aov(Densidade~Coleta, data=dados.area2)
dens.coleta.area2.orig=summary(dens.coleta.area2)[[1]]$"F value"[1]
dens.coleta.area2.perm=numeric(Nperm)
dens.coleta.area2.perm[1]=dens.coleta.area2.orig

alt.coleta.area2=aov(Altura~Coleta, data=dados.area2)
alt.coleta.area2.orig=summary(alt.coleta.area2)[[1]]$"F value"[1]
alt.coleta.area2.perm=numeric(Nperm)
alt.coleta.area2.perm[1]=alt.coleta.area2.orig

biom.coleta.area2=aov(Biomassa~Coleta, data=dados.area2)
biom.coleta.area2.orig=summary(biom.coleta.area2)[[1]]$"F value"[1]
biom.coleta.area2.perm=numeric(Nperm)
biom.coleta.area2.perm[1]=biom.coleta.area2.orig

serap.coleta.area2=aov(Serrapilheira~Coleta, data=dados.area2)
serap.coleta.area2.orig=summary(serap.coleta.area2)[[1]]$"F value"[1]
serap.coleta.area2.perm=numeric(Nperm)
serap.coleta.area2.perm[1]=serap.coleta.area2.orig


for(i in 2:Nperm) {
	dados.temp=dados.area2
	dados.temp$Coleta=sample(dados.temp$Coleta)
	
	dens.foo=aov(Densidade~Coleta, data=dados.temp)
	alt.foo=aov(Altura~Coleta, data=dados.temp)
	biom.foo=aov(Biomassa~Coleta, data=dados.temp)
	serap.foo=aov(Serrapilheira~Coleta, data=dados.temp)

	
	dens.coleta.area2.perm[i]=summary(dens.foo)[[1]]$"F value"[1]
	alt.coleta.area2.perm[i]=summary(alt.foo)[[1]]$"F value"[1]
	biom.coleta.area2.perm[i]=summary(biom.foo)[[1]]$"F value"[1]
	serap.coleta.area2.perm[i]=summary(serap.foo)[[1]]$"F value"[1]
	if (i %% 100 == 0) print (i)

	}



signif.dens.coleta.area2=sum(dens.coleta.area2.perm >= dens.coleta.area2.orig) / Nperm		
signif.alt.coleta.area2=sum(alt.coleta.area2.perm >= alt.coleta.area2.orig) / Nperm		
signif.biom.coleta.area2=sum(biom.coleta.area2.perm >= biom.coleta.area2.orig) / Nperm		
signif.serap.coleta.area2=sum(serap.coleta.area2.perm >= serap.coleta.area2.orig) / Nperm		

signif.dens.coleta.area2
signif.alt.coleta.area2
signif.biom.coleta.area2
signif.serap.coleta.area2



### Comparações par-a-par entre os anos - área queimada

Ncoletas=4

coletas=Coletas



dens.coleta.area1.pairwise=matrix(ncol=4, nrow=4)
alt.coleta.area1.pairwise=matrix(ncol=4, nrow=4)
biom.coleta.area1.pairwise=matrix(ncol=4, nrow=4)
serap.coleta.area1.pairwise=matrix(ncol=4, nrow=4)

for(i in 1:(Ncoletas-1)) {
	coleta.a=coletas[i]
	for(j in (i+1):Ncoletas) {
		coleta.b=coletas[j]
		dados.a=subset(dados.area1, Coleta==coleta.a)
		dados.b=subset(dados.area1, Coleta==coleta.b)
		dens.a=mean(dados.a$Densidade)
		dens.b=mean(dados.b$Densidade)
		dens.ab=abs(dens.a-dens.b)

		length.a=nrow(dados.a)
		length.b=nrow(dados.b)

		alt.a=mean(dados.a$Altura)
		alt.b=mean(dados.b$Altura)
		alt.ab=abs(alt.a-alt.b)

		biom.a=mean(dados.a$Biomassa)
		biom.b=mean(dados.b$Biomassa)
		biom.ab=abs(biom.a-biom.b)

		serap.a=mean(dados.a$Serrapilheira)
		serap.b=mean(dados.b$Serrapilheira)
		serap.ab=abs(serap.a-serap.b)


		join.dens=c(dados.a$Densidade, dados.b$Densidade)
		join.alt=c(dados.a$Altura, dados.b$Altura)
		join.biom=c(dados.a$Biomassa, dados.b$Biomassa)
		join.serap=c(dados.a$Serrapilheira, dados.b$Serrapilheira)

		dens.ab.perm=numeric(Nperm)
		alt.ab.perm=numeric(Nperm)
		biom.ab.perm=numeric(Nperm)
		serap.ab.perm=numeric(Nperm)

		dens.ab.perm[1]=dens.ab
		alt.ab.perm[1]=alt.ab
		biom.ab.perm[1]=biom.ab
		serap.ab.perm[1]=serap.ab

		for (k in 2:Nperm) {
			sorteio=sample(1:(length.a+length.b))
			sorteio.a=sorteio[1:length.a]
			sorteio.b=sorteio[(length.a+1):length(sorteio)]

			dens.a.perm=mean(join.dens[sorteio.a])
			dens.b.perm=mean(join.dens[sorteio.b])
			dens.ab.perm[k]=abs(dens.a.perm-dens.b.perm)

			alt.a.perm=mean(join.alt[sorteio.a])
			alt.b.perm=mean(join.alt[sorteio.b])
			alt.ab.perm[k]=abs(alt.a.perm-alt.b.perm)

			biom.a.perm=mean(join.biom[sorteio.a])
			biom.b.perm=mean(join.biom[sorteio.b])
			biom.ab.perm[k]=abs(biom.a.perm-biom.b.perm)

			serap.a.perm=mean(join.serap[sorteio.a])
			serap.b.perm=mean(join.serap[sorteio.b])
			serap.ab.perm[k]=abs(serap.a.perm-serap.b.perm)

			if(k %% 500 == 0) print(c(i,j,k))
			}
		dens.coleta.area1.pairwise[i,j]=sum(dens.ab.perm >= dens.ab) / Nperm
		alt.coleta.area1.pairwise[i,j]=sum(alt.ab.perm >= alt.ab) / Nperm
		biom.coleta.area1.pairwise[i,j]=sum(biom.ab.perm >= biom.ab) / Nperm
		serap.coleta.area1.pairwise[i,j]=sum(serap.ab.perm >= serap.ab) / Nperm
		}
	}

dens.coleta.area1.pairwise
alt.coleta.area1.pairwise
biom.coleta.area1.pairwise
serap.coleta.area1.pairwise








### Comparações par-a-par entre os anos - área controle


Ncoletas=4

coletas=Coletas



dens.coleta.area2.pairwise=matrix(ncol=4, nrow=4)
alt.coleta.area2.pairwise=matrix(ncol=4, nrow=4)
biom.coleta.area2.pairwise=matrix(ncol=4, nrow=4)
serap.coleta.area2.pairwise=matrix(ncol=4, nrow=4)

for(i in 1:(Ncoletas-1)) {
	coleta.a=coletas[i]
	for(j in (i+1):Ncoletas) {
		coleta.b=coletas[j]
		dados.a=subset(dados.area2, Coleta==coleta.a)
		dados.b=subset(dados.area2, Coleta==coleta.b)
		dens.a=mean(dados.a$Densidade)
		dens.b=mean(dados.b$Densidade)
		dens.ab=abs(dens.a-dens.b)

		length.a=nrow(dados.a)
		length.b=nrow(dados.b)

		alt.a=mean(dados.a$Altura)
		alt.b=mean(dados.b$Altura)
		alt.ab=abs(alt.a-alt.b)

		biom.a=mean(dados.a$Biomassa)
		biom.b=mean(dados.b$Biomassa)
		biom.ab=abs(biom.a-biom.b)

		serap.a=mean(dados.a$Serrapilheira)
		serap.b=mean(dados.b$Serrapilheira)
		serap.ab=abs(serap.a-serap.b)


		join.dens=c(dados.a$Densidade, dados.b$Densidade)
		join.alt=c(dados.a$Altura, dados.b$Altura)
		join.biom=c(dados.a$Biomassa, dados.b$Biomassa)
		join.serap=c(dados.a$Serrapilheira, dados.b$Serrapilheira)

		dens.ab.perm=numeric(Nperm)
		alt.ab.perm=numeric(Nperm)
		biom.ab.perm=numeric(Nperm)
		serap.ab.perm=numeric(Nperm)

		dens.ab.perm[1]=dens.ab
		alt.ab.perm[1]=alt.ab
		biom.ab.perm[1]=biom.ab
		serap.ab.perm[1]=serap.ab

		for (k in 2:Nperm) {
			sorteio=sample(1:(length.a+length.b))
			sorteio.a=sorteio[1:length.a]
			sorteio.b=sorteio[(length.a+1):length(sorteio)]

			dens.a.perm=mean(join.dens[sorteio.a])
			dens.b.perm=mean(join.dens[sorteio.b])
			dens.ab.perm[k]=abs(dens.a.perm-dens.b.perm)

			alt.a.perm=mean(join.alt[sorteio.a])
			alt.b.perm=mean(join.alt[sorteio.b])
			alt.ab.perm[k]=abs(alt.a.perm-alt.b.perm)

			biom.a.perm=mean(join.biom[sorteio.a])
			biom.b.perm=mean(join.biom[sorteio.b])
			biom.ab.perm[k]=abs(biom.a.perm-biom.b.perm)

			serap.a.perm=mean(join.serap[sorteio.a])
			serap.b.perm=mean(join.serap[sorteio.b])
			serap.ab.perm[k]=abs(serap.a.perm-serap.b.perm)

			if(k %% 500 == 0) print(c(i,j,k))
			}
		dens.coleta.area2.pairwise[i,j]=sum(dens.ab.perm >= dens.ab) / Nperm
		alt.coleta.area2.pairwise[i,j]=sum(alt.ab.perm >= alt.ab) / Nperm
		biom.coleta.area2.pairwise[i,j]=sum(biom.ab.perm >= biom.ab) / Nperm
		serap.coleta.area2.pairwise[i,j]=sum(serap.ab.perm >= serap.ab) / Nperm
		}
	}

dens.coleta.area2.pairwise
alt.coleta.area2.pairwise
biom.coleta.area2.pairwise
serap.coleta.area2.pairwise


### Making nice beanplots

setwd("/home/pavel/Profissional/Pesquisa/MyPapers-Revising/PteridiumBrasilia/resubmitted2")

png(filename="PteridiumBSB_fig01_corrigida4.png", height=21, width=15, res=300, unit="cm")

par(mfrow=c(4,2), mar=c(2,2,1,1), oma=c(3,6,3,5))

pirateplot(Densidade ~ Coleta, data=subset(dados, Area==1), pal="gray", jitter.val=0.05, inf.f.o=0, inf.b.o=0, avg.line.lwd=1, bar.f.o=0,bean.b.o=1, bean.f.o=0.8, bean.f.col="gray80", bean.b.col="black", point.col="black", gl.lwd=0, point.o=0.6,gl.lty=5, gl.col="white", ylab="", xlab="", ylim=range(dados$Densidade)+c(0,max(dados$Densidade/20)), point.cex=1, xaxt="n", yaxt="n")
text(labels=c("a","b","c","b"), x=c(1,2,3,4), y=max(dados$Densidade)+max(dados$Densidade)/20, cex=1.3)
axis(side=2, at=c(0,10,20,30), las=2)
abline(v=1.5, col="red", lwd=2, lty=2)

pirateplot(Densidade ~ Coleta, data=subset(dados, Area==2), pal="gray", jitter.val=0.05, inf.f.o=0, inf.b.o=0, avg.line.lwd=1, bar.f.o=0,bean.b.o=1, bean.f.o=0.8, bean.f.col="white", bean.b.col="black", point.col="black", gl.lwd=0, point.o=0.6,gl.lty=5, gl.col="white", ylab="", xlab="", ylim=range(dados$Densidade)+c(0,max(dados$Densidade/20)), point.cex=1, xaxt="n", yaxt="n")
text(labels=c("a","b","b","b"), x=c(1,2,3,4), y=max(dados$Densidade)+max(dados$Densidade)/20, cex=1.3)

pirateplot(Altura ~ Coleta, data=subset(dados, Area==1), pal="gray", jitter.val=0.05, inf.f.o=0, inf.b.o=0, avg.line.lwd=1, bar.f.o=0,bean.b.o=1, bean.f.o=0.8, bean.f.col="gray80", bean.b.col="black", point.col="black", gl.lwd=0, point.o=0.6,gl.lty=5, gl.col="white", ylab="", xlab="", ylim=range(dados$Altura)+c(0,max(dados$Altura/20)), point.cex=1, xaxt="n", yaxt="n")
text(labels=c("a","a","b","b"), x=c(1,2,3,4), y=max(dados$Altura)+max(dados$Altura)/20, cex=1.3)
axis(side=2, at=c(0,1,2,3), las=2)
abline(v=1.5, col="red", lwd=2, lty=2)

pirateplot(Altura ~ Coleta, data=subset(dados, Area==2), pal="gray", jitter.val=0.05, inf.f.o=0, inf.b.o=0, avg.line.lwd=1, bar.f.o=0,bean.b.o=1, bean.f.o=0.8, bean.f.col="gray80", bean.b.col="black", point.col="black", gl.lwd=0, point.o=0.6,gl.lty=5, gl.col="white", ylab="", xlab="", ylim=range(dados$Altura)+c(0,max(dados$Altura/20)), point.cex=1, xaxt="n", yaxt="n")
text(labels=c("a","b","b","b"), x=c(1,2,3,4), y=max(dados$Altura)+max(dados$Altura)/20, cex=1.3)


pirateplot(Biomassa ~ Coleta, data=subset(dados, Area==1), pal="gray", jitter.val=0.05, inf.f.o=0, inf.b.o=0, avg.line.lwd=1, bar.f.o=0,bean.b.o=1, bean.f.o=0.8, bean.f.col="gray80", bean.b.col="black", point.col="black", gl.lwd=0, point.o=0.6,gl.lty=5, gl.col="white", ylab="", xlab="", ylim=range(dados$Biomassa)+c(0,max(dados$Biomassa/20)), point.cex=1, xaxt="n", yaxt="n")
text(labels=c("a","b","b","b"), x=c(1,2,3,4), y=max(dados$Biomassa)+max(dados$Biomassa)/20, cex=1.3)
axis(side=2, at=c(0,500,1000,1500), las=2)
abline(v=1.5, col="red", lwd=2, lty=2)

pirateplot(Biomassa ~ Coleta, data=subset(dados, Area==2), pal="gray", jitter.val=0.05, inf.f.o=0, inf.b.o=0, avg.line.lwd=1, bar.f.o=0,bean.b.o=1, bean.f.o=0.8, bean.f.col="gray80", bean.b.col="black", point.col="black", gl.lwd=0, point.o=0.6,gl.lty=5, gl.col="white", ylab="", xlab="", ylim=range(dados$Biomassa)+c(0,max(dados$Biomassa/20)), point.cex=1, xaxt="n", yaxt="n")
text(labels=c("a","a","a","a"), x=c(1,2,3,4), y=max(dados$Biomassa)+max(dados$Biomassa)/20, cex=1.3)


pirateplot(Serrapilheira ~ Coleta, data=subset(dados, Area==1), pal="gray", jitter.val=0.05, inf.f.o=0, inf.b.o=0, avg.line.lwd=1, bar.f.o=0,bean.b.o=1, bean.f.o=0.8, bean.f.col="gray80", bean.b.col="black", point.col="black", gl.lwd=0, point.o=0.6,gl.lty=5, gl.col="white", ylab="", xlab="", ylim=c(0,max(dados$Serrapilheira))+c(0,max(dados$Serrapilheira/20)), point.cex=1, yaxt="n")
axis(side=2, at=c(0,500,1000,1500,2000),las=2)
text(labels=c("a","a","a,b","a,b"), x=c(1,2,3,4), y=max(dados$Serrapilheira)+max(dados$Serrapilheira)/20, cex=1.3)
abline(v=1.5, col="red", lwd=2, lty=2)

pirateplot(Serrapilheira ~ Coleta, data=subset(dados, Area==2), pal="gray", jitter.val=0.05, inf.f.o=0, inf.b.o=0, avg.line.lwd=1, bar.f.o=0,bean.b.o=1, bean.f.o=0.8, bean.f.col="gray80", bean.b.col="black", point.col="black", gl.lwd=0, point.o=0.6,gl.lty=5, gl.col="white", ylab="", xlab="", ylim=c(0,max(dados$Serrapilheira))+c(0,max(dados$Serrapilheira/20)), point.cex=1, yaxt="n")
text(labels=c("a","b","b","b"), x=c(1,2,3,4), y=max(dados$Serrapilheira)+max(dados$Serrapilheira)/20, cex=1.3)


mtext(side=2, text=c("Litter \n biomass (g)", "Standing \n biomass (g)", "Height (m)", "Number of \n fronds"), at=c(0, 0.25, 0.5, 0.75)+0.125, outer=T, line=2.5)

mtext(side=3, text=c("Burnt area", "Unburnt area"), at=c(0.25, 0.75), outer=T)

mtext(side=1, text="Sampling period", outer=T)

dev.off()

