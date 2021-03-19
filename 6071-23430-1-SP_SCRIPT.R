Supplementary Material 

Effects of land conversion on ant functional diversity: a global review[W1]

Roberta de Jesus Santos, Pavel Dodonov, Jacques Hubert Charles Delabie


Appendix 2. Statistical analysis script



### Queremos saber como o tipo de resposta depende de:
# Tipo de uso da terra convertido;
# Tipo de vegetação;
# Grupo funcional;
# Zona climática;

# De forma mais fácil e didático analisamos uma coisa por vez.

##### Tipo de resposta em função de tipo de vegetação
setwd() 

dados <- read.table("6071-23430-1-SP_DATABASE.txt", header=T, sep="\t")

response.vegetation <- aggregate(Response ~ NativeHabitat, data=dados, FUN=mean)

##### Mas serão significativamente diferente de zero?
Nperm <- 10000

response.vegetation.MC <- matrix(ncol=Nperm, nrow=length(unique(dados$NativeHabitat)))
row.names(response.vegetation.MC) <- response.vegetation[,1]
response.vegetation.MC[,1] <- response.vegetation[,2]

for(i in 2:Nperm) {
	response.rand <- sample(c(-1,0,1), replace=T, size=nrow(dados))
	response.vegetation.rand <- aggregate(response.rand ~ NativeHabitat, data=dados, FUN=mean)
	response.vegetation.MC[,i] <- response.vegetation.rand[,2]
	print(i)
}

### Plotando os gráficos
par(mfrow=c(3,2), mar=c(3,3,2,2))
for(i in 1:6) {
	hist(response.vegetation.MC[i,], main=rownames(response.vegetation.MC)[i], col="pink")
	abline(v=response.vegetation.MC[i,1], col="purple", lwd=2)
}

### Calculando significância...
response.vegetation$Signif <- NA
for(i in 1:6) {
	foo <- sum(abs(response.vegetation.MC[i,]) >= abs(response.vegetation.MC[i,1]))/Nperm
	response.vegetation$Signif[i] <- foo
}

response.vegetation

#########################################################################
###  Tipo de resposta em função do tipo de hábitat alterado

response.disturbed <- aggregate(Response ~ DisturbedHabitat, data=dados, FUN=mean)
response.disturbed.MC <- matrix(ncol=Nperm, nrow=length(unique(dados$DisturbedHabitat)))
row.names(response.disturbed.MC) <- response.disturbed[,1]
response.disturbed.MC[,1] <- response.disturbed[,2]

for(i in 2:Nperm) {
	response.rand <- sample(c(-1,0,1), replace=T, size=nrow(dados))
	response.disturbed.rand <- aggregate(response.rand ~ DisturbedHabitat, data=dados, FUN=mean)
	response.disturbed.MC[,i] <- response.disturbed.rand[,2]
	print(i)
}
nrow(response.disturbed.MC)

par(mfrow=c(3,2), mar=c(3,3,2,2))
for(i in 1:4) {
	hist(response.disturbed.MC[i,], main=rownames(response.disturbed.MC)[i], col="pink")
	abline(v=response.disturbed.MC[i,1], col="purple", lwd=2)
}

response.disturbed$Signif <- NA
for(i in 1:4) {
	foo <- sum(abs(response.disturbed.MC[i,]) >= abs(response.disturbed.MC[i,1]))/Nperm
	response.disturbed$Signif[i] <- foo
}

response.disturbed

###############################################################
###  Tipo de resposta em função da zona climática
response.climaticz <- aggregate(Response ~ ClimaticZone, data=dados, FUN=mean)

### É diferente de zero???
response.climaticz.MC <- matrix(ncol=Nperm, nrow=length(unique(dados$ClimaticZone)))
row.names(response.climaticz.MC) <- response.climaticz[,1]
response.climaticz.MC[,1] <- response.climaticz[,2]

for(i in 2:Nperm) {
	response.rand <- sample(c(-1,0,1), replace=T, size=nrow(dados))
	response.climaticz.rand <- aggregate(response.rand ~ ClimaticZone, data=dados, FUN=mean)
	response.climaticz.MC[,i] <- response.climaticz.rand[,2]
	print(i)
}
nrow(response.climaticz.MC)

par(mfrow=c(3,2), mar=c(3,3,2,2))
for(i in 1:2) {
	hist(response.climaticz.MC[i,], main=rownames(response.climaticz.MC)[i], col="pink")
	abline(v=response.climaticz.MC[i,1], col="purple", lwd=2)
}

response.climaticz$Signif <- NA
for(i in 1:2) {

	foo <- sum(abs(response.climaticz.MC[i,]) >= abs(response.climaticz.MC[i,1]))/Nperm
	response.climaticz$Signif[i] <- foo
}

response.climaticz


###############################################################
### Para os grupos funcionais

response.fg <- aggregate(Response ~ FunctionalGroup, data=dados, FUN=mean)
response.fg.MC <- matrix(ncol=Nperm, nrow=length(unique(dados$FunctionalGroup)))
row.names(response.fg.MC) <- response.fg[,1]
response.fg.MC[,1] <- response.fg[,2]

for(i in 2:Nperm) {
	response.rand <- sample(c(-1,0,1), replace=T, size=nrow(dados))
	response.fg.rand <- aggregate(response.rand ~ FunctionalGroup, data=dados, FUN=mean)
	response.fg.MC[,i] <- response.fg.rand[,2]
	print(i)
}
nrow(response.fg.MC)

par(mfrow=c(4,4), mar=c(3,3,2,2))
for(i in 1:13) {
	hist(response.fg.MC[i,], main=rownames(response.fg.MC)[i], col="pink")
	abline(v=response.fg.MC[i,1], col="purple", lwd=2)
}

response.fg$Signif <- NA
for(i in 1:13) {
	foo <- sum(abs(response.fg.MC[i,]) >= abs(response.fg.MC[i,1]))/Nperm
	response.fg$Signif[i] <- foo
}

response.fg


########################################################################
###Como os GFs respondem a conversão de habitat nativos nas diferentes zonas climáticas?######

##############################################################
###Criando SUBCONJUNTOS para ZC Tropical e Temperada###

data_trop <- subset(dados, ClimaticZone == "Tropical")
str (data_trop)

data_temp <- subset(dados, ClimaticZone == "Temperate")
str (data_temp)

###############################################################
###Resposta dos grupos funcionais em função da ZC tropical

response.fgTrop <- aggregate(Response ~ FunctionalGroup, data=data_trop, FUN=mean)
response.fgTrop.MC <- matrix(ncol=Nperm, nrow=length(unique(data_trop$FunctionalGroup)))
row.names(response.fgTrop.MC) <- response.fgTrop[,1]
response.fgTrop.MC[,1] <- response.fgTrop[,2]

for(i in 2:Nperm) {
	response.rand <- sample(c(-1,0,1), replace=T, size=nrow(data_trop))
	response.fgTrop.rand <- aggregate(response.rand ~ FunctionalGroup, data=data_trop, FUN=mean)
	response.fgTrop.MC[,i] <- response.fgTrop.rand[,2]
	print(i)
}
nrow(response.fgTrop.MC)

par(mfrow=c(4,4), mar=c(3,3,2,2))
for(i in 1:13) {
	hist(response.fgTrop.MC[i,], main=rownames(response.fgTrop.MC)[i], col="pink")
	abline(v=response.fgTrop.MC[i,1], col="purple", lwd=2)
}

response.fgTrop$Signif <- NA
for(i in 1:13) {
	foo <- sum(abs(response.fgTrop.MC[i,]) >= abs(response.fgTrop.MC[i,1]))/Nperm
	response.fgTrop$Signif[i] <- foo
}

response.fgTrop


############################################################################
### Resposta dos grupos funcionais em função da ZC temperada

response.fgtem <- aggregate(Response ~ FunctionalGroup, data=data_temp, FUN=mean)
response.fgtem.MC <- matrix(ncol=Nperm, nrow=length(unique(data_temp$FunctionalGroup)))
row.names(response.fgtem.MC) <- response.fgtem[,1]
response.fgtem.MC[,1] <- response.fgtem[,2]

for(i in 2:Nperm) {
	response.rand <- sample(c(-1,0,1), replace=T, size=nrow(data_temp))
	response.fgtem.rand <- aggregate(response.rand ~ FunctionalGroup, data=data_temp, FUN=mean)
	response.fgtem.MC[,i] <- response.fgtem.rand[,2]
	print(i)
}
nrow(response.fgtem.MC)

par(mfrow=c(4,4), mar=c(3,3,2,2))
for(i in 1:9) {
	hist(response.fgtem.MC[i,], main=rownames(response.fgtem.MC)[i], col="pink")
	abline(v=response.fgtem.MC[i,1], col="purple", lwd=2)
}

response.fgtem$Signif <- NA
for(i in 1:9) {
	foo <- sum(abs(response.fgtem.MC[i,]) >= abs(response.fgtem.MC[i,1]))/Nperm
	response.fgtem$Signif[i] <- foo

}
response.fgtem
[W1]Conferir título 
