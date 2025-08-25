setwd("e:/Pavel/Profissional/Orientacoes/MES_ori/2020_David/Artigo_David_EIPollination/v2024-01/submission3_ActaOecologica/newAnalyses/")


dados <- read.table("dados_barrasEmpilhadas.txt", sep="\t", header=TRUE)
str(dados)

dados.abund <- subset(dados, Variable=="Abundance")
matriz.abund <- dados.abund[,c("Negative","Neutral","Positive", "Mean", 
	"P.value",
	"MeanJacknife", "PvalueJacknife")]
row.names(matriz.abund) <- dados.abund$Subgroup
matriz.abund <- as.matrix(matriz.abund)
matriz.abund <- matriz.abund[nrow(matriz.abund):1, ]

cores <- c("#cbc9e2", "#9e9ac8", "#6a51a3")[3:1]

rownames(matriz.abund)

rownames(matriz.abund) <- c(
	"Hymenoptera visitors within the\ntemperate region",
	"Hymenoptera visitors",
	"Edges with\nlinear disturbances",
	"Large opening edges",
	"Edges within the\ntropical region",
	"Edges within the\ntemperate region",
	"All observations")

png(filename="../fig5_abundance.png", height=20, width=25, unit="cm", res=300)

par(mar=c(6,13,3,3))
foo <- barplot(t(matriz.abund[,1:3]), legend.text=NULL, horiz=T, las=1, 
	xlab="Number of observations", xlim=c(0,42), xaxt="n",
	main="Pollinator abundance", col=cores)
nvalues <- length(foo)
means<-as.numeric(dados.abund$Mean)
meanJacknife <- as.numeric(dados.abund$MeanJacknife)
axis(side=1, at=c(0,5,10,15,20,25, 30))
legend(legend=c("Negative edge effects", "Neutral edge effects", 
	"Positive edge effects"), fill=cores, 
	x=sum(matriz.abund[1,]), y=1.2)
text(rep(sum(matriz.abund["All observations",]),nvalues), 
	y=foo[nvalues:1], 
	labels=paste("mean = ", 
		formatC(matriz.abund[nvalues:1,"Mean"], digits=2, format="f"),
		paste(" (", formatC(matriz.abund[nvalues:1,"MeanJacknife"],digits=2, format="f"),
		 ")", sep=""),
		"\np = ", formatC(matriz.abund[nvalues:1,"P.value"], digits=4, format="f"), 
		sep="",
		paste(" (", formatC(matriz.abund[nvalues:1,"PvalueJacknife"],digits=4, format="f"),
		 ")", sep="")),
	font=ifelse(matriz.abund[nvalues:1,"PvalueJacknife"] < 0.05, 2, 1),
	pos=4)

dev.off()

#------------------------------------------------------------------------------

dados.rich <- subset(dados, Variable=="Richness")
matriz.rich <- dados.rich[,c("Negative","Neutral","Positive", "Mean", 
	"P.value", "MeanJacknife", "PvalueJacknife")]
row.names(matriz.rich) <- dados.rich$Subgroup
matriz.rich <- as.matrix(matriz.rich)
matriz.rich <- matriz.rich[nrow(matriz.rich):1, ]

rownames.rich.old <- rownames(matriz.rich)

t(t(rownames.rich.old))

rownames(matriz.rich) <- c(
	"Hymenoptera visitors",
	"Edges with\nlinear disturbances",
	"Large opening edges",
	"Edges within the\ntropical region",
	"Edges within the\ntemperate region",
	"All observations"
	)


png(filename="../fig4_richness.png", height=20, width=25, unit="cm", res=300)

par(mar=c(6,13,3,3))
foo <- barplot(t(matriz.rich[,1:3]), legend.text=NULL, horiz=T, las=1, 
	xlab="Number of observations", xlim=c(0,35), xaxt="n",
	main="Pollinator richness", col=cores)
nvalues <- length(foo)
means=as.numeric(dados.rich$Mean)
axis(side=1, at=c(0,5,10,15,20))
legend(legend=c("Negative edge effects", "Neutral edge effects",
	"Positive edge effects"), fill=cores, 
	x=9, y=foo[2]+0.5)
text(rep(sum(matriz.rich["All observations",]),nvalues), 
	y=foo[nvalues:1], 
	labels=paste("mean = ", 
		formatC(matriz.rich[nvalues:1,"Mean"], digits=2, format="f"),
		paste(" (", formatC(matriz.rich[nvalues:1,"MeanJacknife"],digits=2, format="f"),
		 ")", sep=""),
		"\np = ", formatC(matriz.rich[nvalues:1,"P.value"], digits=4, format="f"), 
		sep="",
		paste(" (", formatC(matriz.rich[nvalues:1,"PvalueJacknife"],digits=4, format="f"),
		 ")", sep="")),
	font=ifelse(matriz.rich[nvalues:1,"PvalueJacknife"] < 0.05, 2, 1),
	pos=4)
dev.off()



#------------------------------------------------------------------------------

dados.visit <- subset(dados, Variable=="Visitation")
matriz.visit <- dados.visit[,c("Negative","Neutral","Positive", "Mean", 
	"P.value", "MeanJacknife", "PvalueJacknife")]
row.names(matriz.visit) <- dados.visit$Subgroup
matriz.visit <- as.matrix(matriz.visit)
matriz.visit <- matriz.visit[nrow(matriz.visit):1, ]

rownames.visit.old <- rownames(matriz.visit)
t(t(rownames.visit.old))


rownames(matriz.visit) <- c(
	"Hymenoptera visitors at\nagricultural edges within the\ntemperate region",
	"Hymenoptera visitors within\nthe temperate region",
	"Hymenoptera visitors",
	"Edges with\nlinear disturbances",
	"Large opening edges",
	"Edges within the\ntropical region",
	"Edges within the\ntemperate region",
	"All observations")

png(filename="../fig6_visitation.png", height=20, width=27, unit="cm", res=300)

par(mar=c(6,15,3,3))
foo <- barplot(t(matriz.visit[,1:3]), legend.text=NULL, horiz=T, las=1, 
	xlab="Number of observations", xlim=c(0,50), xaxt="n",
	main="Visitation rate", col=cores)
nvalues <- length(foo)
means=as.numeric(dados.visit$Mean)
axis(side=1, at=c(0,5,10,15,20,25, 30, 35, 40))
legend(legend=c("Negative edge effects", "Neutral edge effects",
	"Positive edge effects"), fill=cores, 
	x=20, y=1.5)
text(rep(sum(matriz.visit["All observations",]),nvalues), 
	y=foo[nvalues:1], 
	labels=paste("mean = ", 
		formatC(matriz.visit[nvalues:1,"Mean"], digits=2, format="f"),
		paste(" (", formatC(matriz.visit[nvalues:1,"MeanJacknife"],digits=2, format="f"),
		 ")", sep=""),
		"\np = ", formatC(matriz.visit[nvalues:1,"P.value"], digits=4, format="f"), 
		sep="",
		paste(" (", formatC(matriz.visit[nvalues:1,"PvalueJacknife"],digits=4, format="f"),
		 ")", sep="")),
	font=ifelse(matriz.visit[nvalues:1,"PvalueJacknife"] < 0.05, 2, 1),
	pos=4)

dev.off()

#------------------------------------------------------------------------------


dados.product <- subset(dados, Variable=="Production")
matriz.product <- dados.product[,c("Negative","Neutral","Positive", "Mean", 
	"P.value", "MeanJacknife", "PvalueJacknife")]
row.names(matriz.product) <- dados.product$Subgroup
matriz.product <- as.matrix(matriz.product)
matriz.product <- matriz.product[nrow(matriz.product):1, ]

rownames.product.old <- rownames(matriz.product)

t(t(rownames.product.old))

rownames(matriz.product) <- c(
	"Plants visited by\nHymenoptera, Diptera and Lepidoptera",
	"Plants visited by\nHymenoptera and Diptera",
	"Plants visited by\nHymenoptera",
	"Edges with\nlinear disturbances",
	"Large opening edges",
	"Edges within the\ntropical region",
	"Edges within the\ntemperate region",
	"All observations"
	)


png(filename="../fig7_production.png", height=20, width=25, unit="cm", res=300)

par(mar=c(6,16,3,3))
foo <- barplot(t(matriz.product[,1:3]), legend.text=NULL, horiz=T, las=1, 
	xlab="Number of observations", xlim=c(0,60), xaxt="n",
	main="Fruit and seed production", col=cores)
nvalues <- length(foo)
means=as.numeric(dados.product$Mean)
axis(side=1, at=c(0,5,10,15,20,25,30,35))
legend(legend=c("Negative edge effects", "Neutral edge effects",
	"Positive edge effects"), fill=cores, 
	x=10, y=3.5)
text(rep(sum(matriz.product["All observations",]),nvalues), 
	y=foo[nvalues:1], 
	labels=paste("mean = ", 
		formatC(matriz.product[nvalues:1,"Mean"], digits=2, format="f"),
		paste(" (", formatC(matriz.product[nvalues:1,"MeanJacknife"],digits=2, format="f"),
		 ")", sep=""),
		"\np = ", formatC(matriz.product[nvalues:1,"P.value"], digits=4, format="f"), 
		sep="",
		paste(" (", formatC(matriz.product[nvalues:1,"PvalueJacknife"],digits=4, format="f"),
		 ")", sep="")),
	font=ifelse(matriz.product[nvalues:1,"PvalueJacknife"] < 0.05, 2, 1),
	pos=4)

dev.off()



#----------------------------------------------------------------------
# Será que fica bom tudo em uma figura só?

png(filename="figure8_combined.png", height=35, width=47, unit="cm", res=300)

par(mfrow=c(2,2), mar=c(3,16,3,3), oma=c(3,3,2,2))
foo <- barplot(t(matriz.rich[,1:3]), legend.text=NULL, horiz=T, las=1, 
	xlab="", main="a. Richness", 
	xlim=c(0,45), xaxt="n")
nvalues <- length(foo)
means=as.numeric(dados.rich$Mean)
axis(side=1)
legend(legend=c("Negative", "Neutral", "Positive"), fill=grey.colors(n=3), 
	x=20, y=5)
text(x=37, 
	y=foo[nvalues:1], 
	labels=paste("mean = ", 
		formatC(matriz.rich[nvalues:1,4], digits=3, format="f"),
		"\np = ", formatC(matriz.rich[nvalues:1,5], digits=4, format="f"), sep=""),
	pos=4)

foo <- barplot(t(matriz.abund[,1:3]), legend.text=NULL, horiz=T, las=1, 
	xlab="", main="b. Abundance", 
	xlim=c(0,45), xaxt="n")
nvalues <- length(foo)
means=as.numeric(dados.abund$Mean)
axis(side=1)
text(x=37, 
	y=foo[nvalues:1], 
	labels=paste("mean = ", 
		formatC(matriz.abund[nvalues:1,4], digits=3, format="f"),
		"\np = ", formatC(matriz.abund[nvalues:1,5], digits=4, format="f"), sep=""),
	pos=4)

foo <- barplot(t(matriz.visit[,1:3]), legend.text=NULL, horiz=T, las=1, 
	xlab="", main="c. Visitation", 
	xlim=c(0,45), xaxt="n")
nvalues <- length(foo)
means=as.numeric(dados.visit$Mean)
axis(side=1)
text(x=37, 
	y=foo[nvalues:1], 
	labels=paste("mean = ", 
		formatC(matriz.visit[nvalues:1,4], digits=3, format="f"),
		"\np = ", formatC(matriz.visit[nvalues:1,5], digits=4, format="f"), sep=""),
	pos=4)

foo <- barplot(t(matriz.product[,1:3]), legend.text=NULL, horiz=T, las=1, 
	xlab="", main="d. Seed production", 
	xlim=c(0,45), xaxt="n")
nvalues <- length(foo)
means=as.numeric(dados.product$Mean)
axis(side=1)
text(x=37, 
	y=foo[nvalues:1], 
	labels=paste("mean = ", 
		formatC(matriz.product[nvalues:1,4], digits=3, format="f"),
		"\np = ", formatC(matriz.product[nvalues:1,5], digits=4, format="f"), sep=""),
	pos=4)

mtext(side=1, text="Number of observations", outer=TRUE, cex=1.5)

dev.off()


