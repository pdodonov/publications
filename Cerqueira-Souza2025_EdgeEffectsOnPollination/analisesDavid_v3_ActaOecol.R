voteCount <- function(x, Nrand=1E4, toPlot=FALSE) {
	observed <- mean(x)
	random <- numeric(Nrand)
	random[1] <- observed
	for(i in 2:Nrand) {
		random[i] <- mean(sample(c(-1,0,1), size=length(x), replace=T))
	}
	if(toPlot) {
		hist(random)
		abline(v=observed, col="red", lwd=2)
	}
	signif <- sum(abs(random) >= abs(observed))/Nrand
	result <- c(observed, signif)
	names(result) <- c("MeanEI", "Significance")
	return(result)
}

sensitivity <- function(x, study, show.iteration=FALSE) {
	studies <- sort(unique(study))
	Nstudies <- length(studies)
	result <- matrix(nrow=Nstudies+1, ncol=2)
	colnames(result) <- c("MeanEI", "Significance")
	result[1,] <- voteCount(x)
	rownames(result) <- 1:nrow(result)
	rownames(result)[1] <- "None removed"
	for (i in 1:Nstudies) {
		study.i <- studies[i]
		x.i <- x[study != study.i]
		result[i+1,] <- voteCount(x.i)
		row.names(result)[1+i] <- paste(study.i, "removed")
		if(show.iteration) print(i)
	}
	return(result)
}


library(tidyverse)

setwd("e:/Pavel/Profissional/Orientacoes/MES_ori/2020_David/Artigo_David_EIPollination/v2024-01/submission3_ActaOecologica/newAnalyses/")

dados <- read.table("dados_david_paraanalise.csv", stringsAsFactors=TRUE,
	sep="\t", header=TRUE)
dados$Estudo <- as.character(dados$Estudo)
str(dados)


foo <- strsplit(as.character(dados$Estudo), " ")
ano <- numeric(length(foo))
for(i in 1:length(foo)) {
	bar <- foo[[i]][length(foo[[i]])]
	ano[i] <- as.numeric(bar)
}
ano

dados$Ano <- ano


### Studies per year per region

dados.YR <- dados[,c("Estudo", "Ano", "Regiao")]
dados.YR <- unique(dados.YR)
nrow(dados.YR)

dados.YR.N <- aggregate(Estudo ~ Ano + Regiao, data=dados.YR, FUN=length)

dados.YR.N <- as.data.frame(pivot_wider(dados.YR.N, id_cols = Ano, 
	names_from=Regiao, 
	values_from=Estudo))
dados.YR.N <- replace_na(dados.YR.N, list(Temperada=0, Tropical=0))

dados.YR.N <- dados.YR.N[order(dados.YR.N$Ano),]


dados.YR.N$CumulTemperada <- cumsum(dados.YR.N$Temperada)
dados.YR.N$CumulTropical <- cumsum(dados.YR.N$Tropical)

dados.YR.N$Total <- dados.YR.N$Temperada + dados.YR.N$Tropical
str(dados.YR.N)

dados.YR.N$CumulTotal <- cumsum(dados.YR.N$Total)


color.temperate <- "#998ec3"
color.tropical <- "#f1a340"

png(filename="Fig2_PapersPerYear.png", height=20, width=20, unit="cm", res=300)
plot(CumulTemperada ~ Ano, 
	data=dados.YR.N, type="s", lwd=2, 
	col=color.temperate, xlab="Year", ylab="Cumulative number of papers",
	ylim=c(0,max(dados.YR.N$CumulTemperada)))
lines(CumulTropical ~ Ano, data=dados.YR.N, type="s", lwd=2,
	col=color.tropical)
legend(x=2000, y=20, fill=c(color.tropical, color.temperate),
	legend=c("Tropical region", "Temperate region"))
dev.off()


#####################################################################

### Which subsets to analyze:
# abundance
# abundance - temperate region
# abundance - tropical region
# abundance – agricultural edges
# abundance - linear opening edges
# abundance - Diptera
# abundance - Hymenoptera
# abundance - Hymenoptera in the temperate region
# richness
# richness - temperate region
# richness - tropical region
# richness - agricultural edges
# richness - linear opening edges
# richness - Hymenoptera
# visitation
# visitation - temperate region
# visitation - tropical region
# visitation - agricultural edges
# visitation - linear opening edges
# visitation - Hymenoptera
# production
# production - temperate region
# production - tropical region
# production - agricultural edges
# production - linear opening edges
# production - plants visited by Hymenoptera
# production - plants visited by Hymenoptera and Diptera
# production - plants visited by Hymenoptera, Diptera and Lepidoptera.


dados.abundancia <- subset(dados, Abundancia==1)
length(unique(dados.abundancia$Estudo))
nrow(dados.abundancia) 
table(dados.abundancia$Resposta)
sensitivity(x = dados.abundancia$Resposta, study=dados.abundancia$Estudo)

dados.abund.temper <- subset(dados, Abundancia==1 & Regiao=="Temperada")
length(unique(dados.abund.temper$Estudo))
nrow(dados.abund.temper) 
table(dados.abund.temper$Resposta)
sensitivity(x = dados.abund.temper$Resposta, study=dados.abund.temper$Estudo)

dados.abund.trop <- subset(dados, Abundancia==1 & Regiao=="Tropical")
length(unique(dados.abund.trop$Estudo))
nrow(dados.abund.trop) 
table(dados.abund.trop$Resposta)
sensitivity(x = dados.abund.trop$Resposta, study=dados.abund.trop$Estudo)

dados.abund.agric <- subset(dados, Abundancia == 1 & Agropecuaria_ocupaca == 1)
length(unique(dados.abund.agric$Estudo))
nrow(dados.abund.agric) 
table(dados.abund.agric$Resposta)
sensitivity(x = dados.abund.agric$Resposta, study=dados.abund.agric$Estudo)

dados.abund.linear <- subset(dados, Abundancia == 1 & Aberturas_lineares == 1)
length(unique(dados.abund.linear$Estudo))
nrow(dados.abund.linear) 
table(dados.abund.linear$Resposta)
sensitivity(x = dados.abund.linear$Resposta, study=dados.abund.linear$Estudo)

dados.abund.hymenoptera <- subset(dados, Abundancia == 1 & 
	Inseto_ordem == "Hymenoptera")
length(unique(dados.abund.hymenoptera$Estudo))
nrow(dados.abund.hymenoptera) 
table(dados.abund.hymenoptera$Resposta)
sensitivity(x = dados.abund.hymenoptera$Resposta, 
	study=dados.abund.hymenoptera$Estudo)

dados.abund.hymenoptera.temper <- subset(dados, Abundancia == 1 & 
	Inseto_ordem == "Hymenoptera" & Regiao == "Temperada")
length(unique(dados.abund.hymenoptera.temper$Estudo))
nrow(dados.abund.hymenoptera.temper) 
table(dados.abund.hymenoptera.temper$Resposta)
sensitivity(x = dados.abund.hymenoptera.temper$Resposta, 
	study=dados.abund.hymenoptera.temper$Estudo)

dados.riqueza <- subset(dados, Riqueza==1)
length(unique(dados.riqueza$Estudo))
nrow(dados.riqueza) 
table(dados.riqueza$Resposta)
sensitivity(x = dados.riqueza$Resposta, study=dados.riqueza$Estudo)

dados.riqueza.temper <- subset(dados, Riqueza==1 & Regiao == "Temperada")
length(unique(dados.riqueza.temper$Estudo))
nrow(dados.riqueza.temper) 
table(dados.riqueza.temper$Resposta)
sensitivity(x = dados.riqueza.temper$Resposta, study=dados.riqueza.temper$Estudo)

dados.riqueza.trop <- subset(dados, Riqueza==1 & Regiao == "Tropical")
length(unique(dados.riqueza.trop$Estudo))
nrow(dados.riqueza.trop) 
table(dados.riqueza.trop$Resposta)
sensitivity(x = dados.riqueza.trop$Resposta, study=dados.riqueza.trop$Estudo)

dados.riqueza.agric <- subset(dados, Riqueza==1 & 
	Agropecuaria_ocupaca == 1)
length(unique(dados.riqueza.agric$Estudo))
nrow(dados.riqueza.agric) 
table(dados.riqueza.agric$Resposta)
sensitivity(x = dados.riqueza.agric$Resposta, study=dados.riqueza.agric$Estudo)

dados.riqueza.linear <- subset(dados, Riqueza==1 & 
	Aberturas_lineares == 1)
length(unique(dados.riqueza.linear$Estudo))
nrow(dados.riqueza.linear) 
table(dados.riqueza.linear$Resposta)
sensitivity(x = dados.riqueza.linear$Resposta, study=dados.riqueza.linear$Estudo)


dados.riqueza.hymenoptera <- subset(dados, Riqueza==1 & 
	Inseto_ordem == "Hymenoptera")
length(unique(dados.riqueza.hymenoptera$Estudo))
nrow(dados.riqueza.hymenoptera) 
table(dados.riqueza.hymenoptera$Resposta)
sensitivity(x = dados.riqueza.hymenoptera$Resposta, study=dados.riqueza.hymenoptera$Estudo)

dados.visit <- subset(dados, Visitacao == 1)
length(unique(dados.visit$Estudo))
nrow(dados.visit) 
table(dados.visit$Resposta)
sensitivity(x = dados.visit$Resposta, study=dados.visit$Estudo)

dados.visit.temper <- subset(dados, Visitacao == 1 & Regiao=="Temperada")
length(unique(dados.visit.temper$Estudo))
nrow(dados.visit.temper) 
table(dados.visit.temper$Resposta)
sensitivity(x = dados.visit.temper$Resposta, study=dados.visit.temper$Estudo)

dados.visit.trop <- subset(dados, Visitacao == 1 & Regiao=="Tropical")
length(unique(dados.visit.trop$Estudo))
nrow(dados.visit.trop) 
table(dados.visit.trop$Resposta)
sensitivity(x = dados.visit.trop$Resposta, study=dados.visit.trop$Estudo)

dados.visit.agric <- subset(dados, Visitacao == 1 & 
	Agropecuaria_ocupaca==1)
length(unique(dados.visit.agric$Estudo))
nrow(dados.visit.agric) 
table(dados.visit.agric$Resposta)
sensitivity(x = dados.visit.agric$Resposta, study=dados.visit.agric$Estudo)

dados.visit.linear <- subset(dados, Visitacao == 1 & 
	Aberturas_lineares==1)
length(unique(dados.visit.linear$Estudo))
nrow(dados.visit.linear) 
table(dados.visit.linear$Resposta)
sensitivity(x = dados.visit.linear$Resposta, study=dados.visit.linear$Estudo)

dados.visit.hymenoptera <- subset(dados, Visitacao == 1 & 
	Inseto_ordem=="Hymenoptera")
length(unique(dados.visit.hymenoptera$Estudo))
nrow(dados.visit.hymenoptera) 
table(dados.visit.hymenoptera$Resposta)
sensitivity(x = dados.visit.hymenoptera$Resposta, study=dados.visit.hymenoptera$Estudo)

dados.visit.temper.hymenoptera <- subset(dados, Visitacao == 1 & 
	Inseto_ordem=="Hymenoptera" & Regiao=="Temperada")
length(unique(dados.visit.temper.hymenoptera$Estudo))
nrow(dados.visit.temper.hymenoptera) 
table(dados.visit.temper.hymenoptera$Resposta)
sensitivity(x = dados.visit.temper.hymenoptera$Resposta, 
	study=dados.visit.temper.hymenoptera$Estudo)

dados.visit.temper.hymenoptera <- subset(dados, Visitacao == 1 & 
	Inseto_ordem=="Hymenoptera" & Regiao=="Temperada")
length(unique(dados.visit.temper.hymenoptera$Estudo))
nrow(dados.visit.temper.hymenoptera) 
table(dados.visit.temper.hymenoptera$Resposta)
sensitivity(x = dados.visit.temper.hymenoptera$Resposta, 
	study=dados.visit.temper.hymenoptera$Estudo)

dados.visit.temper.agric.hymenoptera <- subset(dados, Visitacao == 1 & 
	Inseto_ordem=="Hymenoptera" & Regiao=="Temperada" &
	Agropecuaria_ocupaca==1)
length(unique(dados.visit.temper.agric.hymenoptera$Estudo))
nrow(dados.visit.temper.agric.hymenoptera) 
table(dados.visit.temper.agric.hymenoptera$Resposta)
sensitivity(x = dados.visit.temper.agric.hymenoptera$Resposta, 
	study=dados.visit.temper.agric.hymenoptera$Estudo)

dados.prod <- subset(dados, (Producao_frutos==1|Producao_sementes==1))
length(unique(dados.prod$Estudo))
nrow(dados.prod) 
table(dados.prod$Resposta)
sensitivity(x = dados.prod$Resposta, 
	study=dados.prod$Estudo)

dados.prod.temper <- subset(dados, (Producao_frutos==1|Producao_sementes==1) &
	Regiao=="Temperada")
length(unique(dados.prod.temper$Estudo))
nrow(dados.prod.temper) 
table(dados.prod.temper$Resposta)
sensitivity(x = dados.prod.temper$Resposta, 
	study=dados.prod.temper$Estudo)

dados.prod.trop <- subset(dados, (Producao_frutos==1|Producao_sementes==1) &
	Regiao=="Tropical")
length(unique(dados.prod.trop$Estudo))
nrow(dados.prod.trop) 
table(dados.prod.trop$Resposta)
sensitivity(x = dados.prod.trop$Resposta, 
	study=dados.prod.trop$Estudo)

dados.prod.agric <- subset(dados, (Producao_frutos==1|Producao_sementes==1) &
	Agropecuaria_ocupaca==1)
length(unique(dados.prod.agric$Estudo))
nrow(dados.prod.agric) 
table(dados.prod.agric$Resposta)
sensitivity(x = dados.prod.agric$Resposta, 
	study=dados.prod.agric$Estudo)

dados.prod.linear <- subset(dados, (Producao_frutos==1|Producao_sementes==1) &
	Aberturas_lineares==1)
length(unique(dados.prod.linear$Estudo))
nrow(dados.prod.linear) 
table(dados.prod.linear$Resposta)
sensitivity(x = dados.prod.linear$Resposta, 
	study=dados.prod.linear$Estudo)

dados.prod.hymenoptera <- subset(dados, 
	(Producao_frutos==1|Producao_sementes==1) &
	Inseto_ordem == "Hymenoptera")
length(unique(dados.prod.hymenoptera$Estudo))
nrow(dados.prod.hymenoptera) 
table(dados.prod.hymenoptera$Resposta)
sensitivity(x = dados.prod.hymenoptera$Resposta, 
	study=dados.prod.hymenoptera$Estudo)

dados.prod.hymenoptera_diptera <- subset(dados, 
	(Producao_frutos==1|Producao_sementes==1) &
	Inseto_ordem == "Hymenoptera e Diptera")
length(unique(dados.prod.hymenoptera_diptera$Estudo))
nrow(dados.prod.hymenoptera_diptera) 
table(dados.prod.hymenoptera_diptera$Resposta)
sensitivity(x = dados.prod.hymenoptera_diptera$Resposta, 
	study=dados.prod.hymenoptera_diptera$Estudo)

dados.prod.hymenoptera_diptera_lepidoptera <- subset(dados, 
	(Producao_frutos==1|Producao_sementes==1) &
	Inseto_ordem == "Hymenoptera, Diptera, Lepidoptera")
length(unique(dados.prod.hymenoptera_diptera_lepidoptera$Estudo))
nrow(dados.prod.hymenoptera_diptera_lepidoptera) 
table(dados.prod.hymenoptera_diptera_lepidoptera$Resposta)
sensitivity(x = dados.prod.hymenoptera_diptera_lepidoptera$Resposta, 
	study=dados.prod.hymenoptera_diptera_lepidoptera$Estudo)

dados.prod.hymenoptera_diptera_lepidoptera <- subset(dados, 
	(Producao_frutos==1|Producao_sementes==1) &
	Inseto_ordem == "Hymenoptera, Diptera, Lepidoptera")
length(unique(dados.prod.hymenoptera_diptera_lepidoptera$Estudo))
nrow(dados.prod.hymenoptera_diptera_lepidoptera) 
table(dados.prod.hymenoptera_diptera_lepidoptera$Resposta)
sensitivity(x = dados.prod.hymenoptera_diptera_lepidoptera$Resposta, 
	study=dados.prod.hymenoptera_diptera_lepidoptera$Estudo)

dados.crossing <- subset(dados, Taxa_cruzamento==1)
length(unique(dados.crossing$Estudo))
nrow(dados.crossing) 
table(dados.crossing$Resposta)
sensitivity(x = dados.crossing$Resposta, 
	study=dados.crossing$Estudo)

dados.crossing.temper <- subset(dados, Taxa_cruzamento==1 &
	Regiao == "Temperada")
length(unique(dados.crossing.temper$Estudo))
nrow(dados.crossing.temper) 
table(dados.crossing.temper$Resposta)
sensitivity(x = dados.crossing.temper$Resposta, 
	study=dados.crossing.temper$Estudo)



observacoesPorFamilia <- table(dados$Planta_familia, dados$Estudo)
estudosPorFamilia <- observacoesPorFamilia > 0
NestudosPorFamilia <- apply(estudosPorFamilia, 1, FUN=sum)
NestudosPorFamilia

observacoesPorOrdem <- table(dados$Inseto_ordem, dados$Estudo)
estudosPorOrdem <- observacoesPorOrdem > 0
NestudosPorOrdem <- apply(estudosPorOrdem, 1, FUN=sum)
NestudosPorOrdem


