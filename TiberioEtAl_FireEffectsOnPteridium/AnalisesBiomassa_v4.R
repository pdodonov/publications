setwd("/home/pavel/Profissional/Pesquisa/MyPapers_working/Fernanda_2012_PteridiumBrasilia")

library(yarrr)

dados=read.table("dados_biomassa.txt", sep="", header=T)
Coletas=unique(dados$Coleta)
Areas=unique(dados$Area)

dados$Coleta.jitter=dados$Coleta+runif(160,-0.2,0.2)

dados$Coleta <- as.factor(dados$Coleta)

dados$Area <- as.factor(dados$Area)

Nperm=5000



plot(Densidade~Coleta.jitter, col=Area, data=dados)
plot(Altura~Coleta.jitter, col=Area, data=dados)
plot(Biomassa~Coleta.jitter, col=Area, data=dados)
plot(Serrapilheira~Coleta.jitter, col=Area, data=dados)

# Diferença entre as áreas antes da queimada:


dados.coleta1=subset(dados,Coleta==1)

dens.area.coleta1=aov(Densidade~Area, data=dados.coleta1)
dens.area.coleta1.orig=summary(dens.area.coleta1)[[1]]$"F value"[1]
dens.area.coleta1.perm=numeric(Nperm)
dens.area.coleta1.perm[1]=dens.area.coleta1.orig

alt.area.coleta1=aov(Altura~Area, data=dados.coleta1)
alt.area.coleta1.orig=summary(alt.area.coleta1)[[1]]$"F value"[1]
alt.area.coleta1.perm=numeric(Nperm)
alt.area.coleta1.perm[1]=alt.area.coleta1.orig

biom.area.coleta1=aov(Biomassa~Area, data=dados.coleta1)
biom.area.coleta1.orig=summary(biom.area.coleta1)[[1]]$"F value"[1]
biom.area.coleta1.perm=numeric(Nperm)
biom.area.coleta1.perm[1]=biom.area.coleta1.orig

serap.area.coleta1=aov(Serrapilheira~Area, data=dados.coleta1)
serap.area.coleta1.orig=summary(serap.area.coleta1)[[1]]$"F value"[1]
serap.area.coleta1.perm=numeric(Nperm)
serap.area.coleta1.perm[1]=serap.area.coleta1.orig


for(i in 2:Nperm) {
	dados.temp=dados.coleta1
	dados.temp$Area=sample(dados.temp$Area)
	
	dens.foo=aov(Densidade~Area, data=dados.temp)
	alt.foo=aov(Altura~Area, data=dados.temp)
	biom.foo=aov(Biomassa~Area, data=dados.temp)
	serap.foo=aov(Serrapilheira~Area, data=dados.temp)

	
	dens.area.coleta1.perm[i]=summary(dens.foo)[[1]]$"F value"[1]
	alt.area.coleta1.perm[i]=summary(alt.foo)[[1]]$"F value"[1]
	biom.area.coleta1.perm[i]=summary(biom.foo)[[1]]$"F value"[1]
	serap.area.coleta1.perm[i]=summary(serap.foo)[[1]]$"F value"[1]
	if (i %% 100 == 0) print (i)

	}



signif.dens.area.coleta1=sum(dens.area.coleta1.perm >= dens.area.coleta1.orig) / Nperm		
signif.alt.area.coleta1=sum(alt.area.coleta1.perm >= alt.area.coleta1.orig) / Nperm		
signif.biom.area.coleta1=sum(biom.area.coleta1.perm >= biom.area.coleta1.orig) / Nperm		
signif.serap.area.coleta1=sum(serap.area.coleta1.perm >= serap.area.coleta1.orig) / Nperm		

signif.dens.area.coleta1
signif.alt.area.coleta1
signif.biom.area.coleta1
signif.serap.area.coleta1

### > signif.dens.area.coleta1
### [1] 2e-04
### > signif.alt.area.coleta1
### [1] 2e-04
### > signif.biom.area.coleta1
### [1] 0.3072
### > signif.serap.area.coleta1
### [1] 4e-04


# Testando se há interação entre ano e área:


dados.resid=dados  #Calculando resíduos

dados.resid$Densidade=resid(aov(Densidade~Area+Coleta, data=dados))
dados.resid$Altura=resid(aov(Altura~Area+Coleta, data=dados))
dados.resid$Biomassa=resid(aov(Biomassa~Area+Coleta, data=dados))
dados.resid$Serrapilheira=resid(aov(Serrapilheira~Area+Coleta, data=dados))



dens.interacao=aov(Densidade~Area*Coleta, data=dados.resid)
dens.interacao.orig=summary(dens.interacao)[[1]]$"F value"[3]
dens.interacao.perm=numeric(Nperm)
dens.interacao.perm[1]=dens.interacao.orig

alt.interacao=aov(Altura~Area*Coleta, data=dados.resid)
alt.interacao.orig=summary(alt.interacao)[[1]]$"F value"[3]
alt.interacao.perm=numeric(Nperm)
alt.interacao.perm[1]=alt.interacao.orig

biom.interacao=aov(Biomassa~Area*Coleta, data=dados.resid)
biom.interacao.orig=summary(biom.interacao)[[1]]$"F value"[3]
biom.interacao.perm=numeric(Nperm)
biom.interacao.perm[1]=biom.interacao.orig

serap.interacao=aov(Serrapilheira~Area*Coleta, data=dados.resid)
serap.interacao.orig=summary(serap.interacao)[[1]]$"F value"[3]
serap.interacao.perm=numeric(Nperm)
serap.interacao.perm[1]=serap.interacao.orig

for(i in 2:Nperm) {
	dados.temp=dados.resid
	ordem.temp = sample(1:nrow(dados.resid))
	dados.temp[,3:7]=dados.temp[ordem.temp,3:7]
	#dados.temp$Area=sample(dados.temp$Area)
	#dados.temp$Coleta=sample(dados.temp$Coleta)
	
	dens.foo=aov(Densidade~Area*Coleta, data=dados.temp)
	alt.foo=aov(Altura~Area*Coleta, data=dados.temp)
	biom.foo=aov(Biomassa~Area*Coleta, data=dados.temp)
	serap.foo=aov(Serrapilheira~Area*Coleta, data=dados.temp)

	
	dens.interacao.perm[i]=summary(dens.foo)[[1]]$"F value"[3]
	alt.interacao.perm[i]=summary(alt.foo)[[1]]$"F value"[3]
	biom.interacao.perm[i]=summary(biom.foo)[[1]]$"F value"[3]
	serap.interacao.perm[i]=summary(serap.foo)[[1]]$"F value"[3]
	if (i %% 100 == 0) print (i)

	}


signif.dens.interacao=sum(dens.interacao.perm >= dens.interacao.orig) / Nperm		
signif.alt.interacao=sum(alt.interacao.perm >= alt.interacao.orig) / Nperm		
signif.biom.interacao=sum(biom.interacao.perm >= biom.interacao.orig) / Nperm		
signif.serap.interacao=sum(serap.interacao.perm >= serap.interacao.orig) / Nperm		

signif.dens.interacao
signif.alt.interacao
signif.biom.interacao
signif.serap.interacao


### > signif.dens.interacao
### [1] 0.0476
### > signif.alt.interacao
### [1] 2e-04
### > signif.biom.interacao
### [1] 4e-04
### > signif.serap.interacao
### [1] 2e-04


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


### > signif.dens.coleta.area1
### [1] 2e-04
### > signif.alt.coleta.area1
### [1] 2e-04
### > signif.biom.coleta.area1
### [1] 2e-04
### > signif.serap.coleta.area1
### [1] 0.006



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

### > signif.dens.coleta.area2
### [1] 4e-04
### > signif.alt.coleta.area2
### [1] 0.0016
### > signif.biom.coleta.area2
### [1] 0.1118
### > signif.serap.coleta.area2
### [1] 2e-04

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


### > dens.coleta.area1.pairwise
###      [,1]  [,2]  [,3]   [,4]
### [1,]   NA 4e-04 2e-04 0.0034
### [2,]   NA    NA 3e-03 0.7642
### [3,]   NA    NA    NA 0.0008
### [4,]   NA    NA    NA     NA
### > alt.coleta.area1.pairwise
###      [,1]  [,2]   [,3]   [,4]
### [1,]   NA 0.017 0.0002 0.0002
### [2,]   NA    NA 0.0018 0.0008
### [3,]   NA    NA     NA 0.7608
### [4,]   NA    NA     NA     NA
### > biom.coleta.area1.pairwise
###      [,1]  [,2]   [,3]   [,4]
### [1,]   NA 2e-04 0.0002 0.0002
### [2,]   NA    NA 0.7232 0.6152
### [3,]   NA    NA     NA 0.4064
### [4,]   NA    NA     NA     NA
### > serap.coleta.area1.pairwise
###      [,1]  [,2]   [,3]   [,4]
### [1,]   NA 6e-04 0.0972 0.0152
### [2,]   NA    NA 0.1274 0.1952
### [3,]   NA    NA     NA 0.6302
### [4,]   NA    NA     NA     NA






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


### > dens.coleta.area2.pairwise
###      [,1]   [,2]   [,3]   [,4]
### [1,]   NA 0.0052 0.0004 0.0004
### [2,]   NA     NA 0.4964 0.5706
### [3,]   NA     NA     NA 0.8884
### [4,]   NA     NA     NA     NA
### > alt.coleta.area2.pairwise
###      [,1]   [,2]   [,3]   [,4]
### [1,]   NA 0.0014 0.0014 0.0014
### [2,]   NA     NA 0.7080 0.5920
### [3,]   NA     NA     NA 0.8632
### [4,]   NA     NA     NA     NA
### > biom.coleta.area2.pairwise
###      [,1]   [,2]   [,3]   [,4]
### [1,]   NA 0.9306 0.2824 0.1632
### [2,]   NA     NA 0.3576 0.1750
### [3,]   NA     NA     NA 0.0138
### [4,]   NA     NA     NA     NA
### > serap.coleta.area2.pairwise
###      [,1]  [,2]   [,3]   [,4]
### [1,]   NA 2e-04 0.0002 0.0002
### [2,]   NA    NA 0.0954 0.0266
### [3,]   NA    NA     NA 0.6780
### [4,]   NA    NA     NA     NA



### Making nice beanplots

setwd("/home/pavel/Profissional/Pesquisa/MyPapers_working/Fernanda_2012_PteridiumBrasilia/resubmitted3")

png(filename="PteridiumBSB_fig02.png", height=21, width=15, res=300, unit="cm")

par(mfrow=c(4,2), mar=c(2,2,1,1), oma=c(3,6,3,5))

pirateplot(Densidade ~ Coleta, data=subset(dados, Area==1), pal="gray", jitter.val=0.05, inf.f.o=0, inf.b.o=0, avg.line.lwd=1, bar.f.o=0,bean.b.o=1, bean.f.o=0.8, bean.f.col="gray80", bean.b.col="black", point.col="black", gl.lwd=0, point.o=0.6,gl.lty=5, gl.col="white", ylab="", xlab="", ylim=range(dados$Densidade)+c(0,max(dados$Densidade/20)), point.cex=1, xaxt="n", yaxt="n")
text(labels=c("a","b","c","b"), x=c(1,2,3,4), y=max(dados$Densidade)+max(dados$Densidade)/20, cex=1.3)
axis(side=2, at=c(0,10,20,30), las=2)
abline(v=1 + 2/9, col="red", lwd=2, lty=2)

pirateplot(Densidade ~ Coleta, data=subset(dados, Area==2), pal="gray", jitter.val=0.05, inf.f.o=0, inf.b.o=0, avg.line.lwd=1, bar.f.o=0,bean.b.o=1, bean.f.o=0.8, bean.f.col="white", bean.b.col="black", point.col="black", gl.lwd=0, point.o=0.6,gl.lty=5, gl.col="white", ylab="", xlab="", ylim=range(dados$Densidade)+c(0,max(dados$Densidade/20)), point.cex=1, xaxt="n", yaxt="n")
text(labels=c("a","b","b","b"), x=c(1,2,3,4), y=max(dados$Densidade)+max(dados$Densidade)/20, cex=1.3)

pirateplot(Altura ~ Coleta, data=subset(dados, Area==1), pal="gray", jitter.val=0.05, inf.f.o=0, inf.b.o=0, avg.line.lwd=1, bar.f.o=0,bean.b.o=1, bean.f.o=0.8, bean.f.col="gray80", bean.b.col="black", point.col="black", gl.lwd=0, point.o=0.6,gl.lty=5, gl.col="white", ylab="", xlab="", ylim=range(dados$Altura)+c(0,max(dados$Altura/20)), point.cex=1, xaxt="n", yaxt="n")
text(labels=c("a","a","b","b"), x=c(1,2,3,4), y=max(dados$Altura)+max(dados$Altura)/20, cex=1.3)
axis(side=2, at=c(0,1,2,3), las=2)
abline(v=1 + 2/9, col="red", lwd=2, lty=2)

pirateplot(Altura ~ Coleta, data=subset(dados, Area==2), pal="gray", jitter.val=0.05, inf.f.o=0, inf.b.o=0, avg.line.lwd=1, bar.f.o=0,bean.b.o=1, bean.f.o=0.8, bean.f.col="gray80", bean.b.col="black", point.col="black", gl.lwd=0, point.o=0.6,gl.lty=5, gl.col="white", ylab="", xlab="", ylim=range(dados$Altura)+c(0,max(dados$Altura/20)), point.cex=1, xaxt="n", yaxt="n")
text(labels=c("a","b","b","b"), x=c(1,2,3,4), y=max(dados$Altura)+max(dados$Altura)/20, cex=1.3)


pirateplot(Biomassa ~ Coleta, data=subset(dados, Area==1), pal="gray", jitter.val=0.05, inf.f.o=0, inf.b.o=0, avg.line.lwd=1, bar.f.o=0,bean.b.o=1, bean.f.o=0.8, bean.f.col="gray80", bean.b.col="black", point.col="black", gl.lwd=0, point.o=0.6,gl.lty=5, gl.col="white", ylab="", xlab="", ylim=range(dados$Biomassa)+c(0,max(dados$Biomassa/20)), point.cex=1, xaxt="n", yaxt="n")
text(labels=c("a","b","b","b"), x=c(1,2,3,4), y=max(dados$Biomassa)+max(dados$Biomassa)/20, cex=1.3)
axis(side=2, at=c(0,500,1000,1500), las=2)
abline(v=1 + 2/9, col="red", lwd=2, lty=2)

pirateplot(Biomassa ~ Coleta, data=subset(dados, Area==2), pal="gray", jitter.val=0.05, inf.f.o=0, inf.b.o=0, avg.line.lwd=1, bar.f.o=0,bean.b.o=1, bean.f.o=0.8, bean.f.col="gray80", bean.b.col="black", point.col="black", gl.lwd=0, point.o=0.6,gl.lty=5, gl.col="white", ylab="", xlab="", ylim=range(dados$Biomassa)+c(0,max(dados$Biomassa/20)), point.cex=1, xaxt="n", yaxt="n")
text(labels=c("a","a","a","a"), x=c(1,2,3,4), y=max(dados$Biomassa)+max(dados$Biomassa)/20, cex=1.3)


pirateplot(Serrapilheira ~ Coleta, data=subset(dados, Area==1), pal="gray", jitter.val=0.05, inf.f.o=0, inf.b.o=0, avg.line.lwd=1, bar.f.o=0,bean.b.o=1, bean.f.o=0.8, bean.f.col="gray80", bean.b.col="black", point.col="black", gl.lwd=0, point.o=0.6,gl.lty=5, gl.col="white", ylab="", xlab="", ylim=c(0,max(dados$Serrapilheira))+c(0,max(dados$Serrapilheira/20)), point.cex=1, yaxt="n")
axis(side=2, at=c(0,500,1000,1500,2000),las=2)
text(labels=c("a","b","a,b","a,b"), x=c(1,2,3,4), y=max(dados$Serrapilheira)+max(dados$Serrapilheira)/20, cex=1.3)
abline(v=1 + 2/9, col="red", lwd=2, lty=2)

pirateplot(Serrapilheira ~ Coleta, data=subset(dados, Area==2), pal="gray", jitter.val=0.05, inf.f.o=0, inf.b.o=0, avg.line.lwd=1, bar.f.o=0,bean.b.o=1, bean.f.o=0.8, bean.f.col="gray80", bean.b.col="black", point.col="black", gl.lwd=0, point.o=0.6,gl.lty=5, gl.col="white", ylab="", xlab="", ylim=c(0,max(dados$Serrapilheira))+c(0,max(dados$Serrapilheira/20)), point.cex=1, yaxt="n")
text(labels=c("a","b","b","b"), x=c(1,2,3,4), y=max(dados$Serrapilheira)+max(dados$Serrapilheira)/20, cex=1.3)


mtext(side=2, text=c("Litter \n biomass (g/m²)", "Standing \n biomass (g/m²)", "Height (m)", "Number of \n fronds"), at=c(0, 0.25, 0.5, 0.75)+0.125, outer=T, line=2.5)

mtext(side=3, text=c("Burnt site", "Unburnt site"), at=c(0.25, 0.75), outer=T)

mtext(side=1, text="Sampling period", outer=T)

dev.off()

