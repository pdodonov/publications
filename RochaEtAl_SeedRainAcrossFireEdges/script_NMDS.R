library(vegan)
setwd("/home/pavel/Profissional/Pesquisa/MyPapers-Working/SeedRainEdge")
dados <- read.table("dados_morfotipo.txt", header=T)

dados.expl <- dados[,1:2]
dados.resp <- as.matrix(dados[,-c(1,2)])
dados.expl$Ambiente <- ifelse(dados.expl$Distancia < 0, "Fire", "Forest")
dados.expl$Ambiente[dados.expl$Distancia==0] <- "Edge"
dados.expl$Ambiente <- as.factor(dados.expl$Ambiente)

# Remove species with less than five occurrence

freq <- apply(dados.resp, 2, function(x) sum(x>0))

dados.resp2 <- dados.resp[,freq>5]

# Remove sites with no species
freq2 <- apply(dados.resp2, 1, sum)

dados.resp3 <- dados.resp2[freq2>0,]

dados.expl2 <- dados.expl[freq2>0,]

# NMDS

NMDS.bray <- metaMDS(dados.resp3, try=100, trymax=100)
NMDS.bray # Stress = 0.20
scores.bray <- scores(NMDS.bray, display="sites")
plot(scores.bray[,2] ~ scores.bray[,1], pch=21, bg=as.numeric(dados.expl2$Ambiente))

NMDS.jaccard <- metaMDS(dados.resp3, distance="jaccard", try=100, trymax=100)
scores.jaccard <- scores(NMDS.jaccard, display="sites")
plot(scores.jaccard[,2] ~ scores.jaccard[,1], pch=21, bg=as.numeric(dados.expl2$Ambiente))

png(filename="SeedRainEdge_NMDS_bray.png", height=20, width=20, unit="cm", res=300)
plot(scores.bray[,2] ~ scores.bray[,1], pch=21, bg=as.numeric(dados.expl2$Ambiente), main="Bray-Curtis")
dev.off()

### Can't resist a PERMANOVA, now can I?

permanova.bray <- adonis(dados.resp3 ~ dados.expl2$Ambiente, permutations=4999)
permanova.bray




