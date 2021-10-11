# The package wmtsa, essential for these analyses, is no longer available at CRAN. In ortder to install it, use the following script:

#library(devtools)
#install_url('https://cran.r-project.org/src/contrib/Archive/wmtsa/wmtsa_2.0-3.tar.gz')

### Load required packages
library(vegan) # used in the function to calculate structural diversity
library(cluster) # used in the function to calculate structural diversity
library(devtools)
library(wmtsa) # required for the wavelet analysis
library(scales) # to make plots with transparency
library(ecodist) # for cluster analysis required to calculate the structural diversity measure
          

# Study areas forest-tundra ecotone and tundra, near Churchill, MB, Canada. The ecotone has two transects separated by a lake. The tundra has three transects separated by lakes and with a road crossing one of them.
# Objective 1a: To assess the spatial scales of structural diversity
# Objective 1b: To compare the spatial scales of structural diveristy between the environments
# Objective 1c: To compare the spatial scales of structural diveristy between structural diversity components. 
# Objective 2: To assess whether structural diversity may be explained by environmental factors, namely soil pH, microtopography, and distance from edges.

### Load custom functions

setwd("/home/pavel/Profissional/Pesquisa/MyPapers/2014_Pavel_StrDivTundra")

source("scripts/StrDivTundra_functions.R")

### Load the data

# Data on the quantity of each structural element in each quadrat

data.all <- read.table("data/StrDivTundra_data.txt", header=T, sep="\t", stringsAsFactors=F)

str(data.all)
unique(data.all$Transect)

# A description of the characteristics of the different structural elements

data.elements <- read.table("data/StrDivTundra_elements.txt", header=T, sep="\t", stringsAsFactors=T)

str(data.elements)

# A description of which elements are to be considered for the different types ("layers") of structural diversity

data.layers <- read.table("data/StrDivTundra_layers.txt", header=T, sep="\t")
str(data.layers, stringsAsFactors=T)

### Remove "tree root", as it had only two occurrences, and remove the overhanging elements because they will not be used in this study.
data.all <- subset(data.all, Element_detail!="Tree_root")
data.all <- subset(data.all, Overhanging!="y")
data.all <- data.all[,-which(colnames(data.all)=="Overhanging")]

### Remove the state of decay of some variables, as it will not be used for the analyses

data.all <- data.all[,-7]

### Combine the data with the elements' characteristics


### Combine data when the same element is observed more than once in a quadrat
elements.paste <- as.character(apply(data.all[,-4], 1, paste, collapse="_"))
# Which combinations occur more than once?
elements.table <- table(elements.paste)
elements.repeated <- names(elements.table)[elements.table>1]
lines.repeated <- which(elements.paste %in% elements.repeated)
data.all.repeated <- data.all[lines.repeated,]
data.all.repeated <- data.all.repeated[order(data.all.repeated$Transect, data.all.repeated$Distance, data.all.repeated$Element_detail),]
for(i in 1:length(elements.repeated)){
	foo <- which(elements.paste %in% elements.repeated[i])
	data.all[foo[1],]$Quantity <- max(data.all[foo,]$Quantity)
	data.all[foo[-1],]$Quantity <- NA
}
data.all <- subset(data.all, !is.na(Quantity))

### Checking if it worked...
elements.paste <- as.character(apply(data.all[,-4], 1, paste, collapse="_"))
elements.table <- table(elements.paste)
elements.repeated <- names(elements.table)[elements.table>1]
elements.repeated # no elements are repeated any longer. :-)


# Testing if all quadrats are present (no quadrat numbers should be absent)...
Nquad <- numeric()
Distmax <- numeric()
transects <- unique(data.all$Transect)
for(i in 1:5) {
	foo <- subset(data.all, Transect==transects[i])
	Nquad[i] <- length(unique(foo$Distance))
	Distmax[i] <- max(foo$Distance)
}
names(Distmax) <- transects
Nquad
Distmax

# All numbers of quadrats are equal to maximum distances (i.e. to the transect lengths), thus all quadrats are present.


### Reorder the elements

elements <- sort(unique(data.all$Element_detail))

t(t(elements))

#### Reorder
# [6,] Ground_soil         
# [4,] Ground_gravel       
# [5,] Ground_rock         
#[16,] Litter_conifer      
#[17,] Litter_graminoid    
#[18,] Litter_herbaceous   
#[19,] Litter_lichen       
#[20,] Litter_twig         
# [1,] Deadwood_log        
# [2,] Deadwood_shrub      
# [3,] Deadwood_snag       
#[11,] Lichen_crustose     
#[12,] Lichen_folious      
#[13,] Lichen_fruticose1   
#[14,] Lichen_fruticose2   
#[15,] Litter_broadleaf    
# [9,] Herbaceous_horsetail
# [8,] Herbaceous_graminoid
# [7,] Herbaceous_forb     
#[10,] Herbaceous_other    
#[22,] Moss_sphagnum-like  
#[21,] Moss_other          
#[23,] Shrub_prostrate     
#[24,] Shrub_standing      
#[25,] Tree_layering       
#[26,] Tree_skirt          
#[27,] Tree_standing   

elements <- elements[c(6,4,5,16,17,18,19,20,1,2,3,11,12,13,14,15,9,8,7,10,22,21,23,24,25,26,27)]

t(t(elements))

### Create some some additional objects

Nelements <- length(elements)

transects <- unique(data.all$Transect)
Ntransects <- length(transects)

### Calculate structural diversity...

layer.names <- colnames(data.layers)[-c(1,2)]

### First organize as a list... To afterwards transform into a dataframe! 

lista.StrDiv <- list()

Nindices <- 2
Nlayers <- 5
Ntransects <- length(transects)

elements.used <- list()


for(i in 1:Nlayers) { #first layer: layer
	layer.i <- layer.names[i]
	elements.i <- data.layers$Element[data.layers[,layer.i]==1]
	data.i <- subset(data.all, as.character(Element_detail) %in% as.character(elements.i))
	lista.StrDiv[[i]] <- list()
	names(lista.StrDiv)[i] <- layer.i
	elements.used[[i]] <- list()
	names(elements.used)[[i]] <- layer.i
	for(j in 1:Ntransects) { #second layer: transect
		data.i.j <- subset(data.i, Transect==transects[j])
		# StrS: Number of structural elements
		StrS.i.j <- as.data.frame(table(data.i.j$Distance))
		colnames(StrS.i.j) <- c("Distance", "StrS")
		# StrD: Our new metric of structural diversity
		data.StrD <- merge(data.i.j, data.elements, sort=F, all.x=T, all.y=F)
		data.StrD.calc <- data.StrD[,c("Layer", "Element", "Standing", "Live", "Woody", "Diameter", "Height_tot")]
		# Keep only columns for which there are no NA values, as these values are likely to introduce artefacts into the analyses
		all.values.present <- apply(data.StrD.calc, 2, function(x) !any(is.na(x)))
		data.StrD.calc <- data.StrD.calc[,all.values.present]
		StrD.weights <- rep(1, ncol(data.StrD.calc))
		StrD.i.j <- numeric()
		distances.i.j <- sort(unique(data.StrD$Distance))
		l <- 1
		for(k in distances.i.j) {
			dist.i.j.k <- k
			data.StrD.calc.k <- data.StrD.calc[data.StrD$Distance==dist.i.j.k,]
			StrD.i.j[l] <- ifelse(nrow(data.StrD.calc.k)==1,0,StrD.calc(data.StrD.calc.k, w=StrD.weights))
			l <- l+1
		}
		StrD.i.j <- data.frame(Distance=distances.i.j,StrD=StrD.i.j)
		test <- all(StrS.i.j[,1] == StrD.i.j[,1])
		if(test == F) print(paste("Something wrong at i=",i, ", and j=",j, sep=""))
		StrSD.i.j <- cbind(StrS.i.j,StrD.i.j)
		StrSD.i.j <- StrSD.i.j[,c(1,2,4)]
		lista.StrDiv[[i]][[j]] <- StrSD.i.j
		elements.used[[i]][[j]] <- colnames(data.StrD.calc)
		names(lista.StrDiv[[i]])[j] <- as.character(transects[j])
		names(elements.used[[i]])[[j]] <- as.character(transects[j])
		print(c(i,j))
	}
}

### How are these indices correlated? I.e. does it matter to use a more complex index, or do they all give the same information?

StrDiv.corr <- matrix(nrow=Nlayers, ncol=Ntransects)
colnames(StrDiv.corr) <- transects
rownames(StrDiv.corr) <- layer.names

for(i in 1:Nlayers) {
	for(j in 1:Ntransects) {
		StrDiv.i.j <- lista.StrDiv[[i]][[j]]
		determination <- cor(StrDiv.i.j[,c(2,3)])[2]^2
		StrDiv.corr[i,j] <- determination
	}
}

StrDiv.corr

StrDiv.corr <- round(StrDiv.corr, 2)

range(StrDiv.corr)


### R2 values between 0.54 and 0.95 between StrS and StrD --> These two metrics apparently have some complementarity.


### Convert list to data.frame...
Nquadrats <- nrow(unique(data.all[,c("Transect","Distance")]))
StrDiv <- matrix(0,nrow=Nquadrats*Nlayers,ncol=3+Nindices)
StrDiv <- as.data.frame(StrDiv)
StrDiv[,1:2] <- unique(data.all[,c("Transect","Distance")])
StrDiv[,3] <- rep(layer.names, each=Nquadrats)
colnames(StrDiv) <- c("Transect", "Distance", "Layer", "StrS", "StrD")

for(i in 1:Nlayers) {
	layer.i <- layer.names[i]
	for(j in 1:Ntransects) {
		transect.j <- transects[j]
		data.i.j <- lista.StrDiv[[i]][[j]]
		StrDiv[StrDiv$Transect == transect.j & StrDiv$Layer == layer.i & StrDiv$Distance %in% data.i.j$Distance,]$StrS <- data.i.j$StrS
		StrDiv[StrDiv$Transect == transect.j & StrDiv$Layer == layer.i & StrDiv$Distance %in% data.i.j$Distance,]$StrD <- data.i.j$StrD
	}
}

write.table(StrDiv, file="data/StrDivTundra_data_StrDiv.txt", sep="\t", row.names=F, col.names=T, quote=F)


# Figure showing correlations

StrDiv.noZero <- subset(StrDiv, StrS>0)

layer.names.plot <- layer.names
layer.names.plot[5] <- "Live plants"

transects.plot <- c("Ecotone 1", "Ecotone 2", "Tundra 1",  "Tundra 2",  "Tundra 3")

cex.outer <- 1.1

png(filename="figures_paper/StrDivTundra_figure2_corr.png", height=25, width=25, unit="cm", res=300)
par(mfrow=c(Nlayers, Ntransects), mar=c(3,2,2,2), oma=c(5,6,4,4))
for(i in 1:Nlayers) {
	for(j in 1:Ntransects) {
		StrDiv.i.j <- subset(StrDiv.noZero, Layer==layer.names[i] & Transect == transects[j])
		plot(StrD ~ StrS, data=StrDiv.i.j, bg=alpha("black",0.3), pch=21, xlab="", ylab="")
		if(j == 1) mtext(text=layer.names.plot[i], side=2, line=4.3, col="gray40", cex=cex.outer)
		if(i == 1) mtext(text=transects.plot[j], side=3, line=1, col="gray40", cex=cex.outer)
		text(x=min(StrDiv.i.j$StrS), y=max(StrDiv.i.j$StrD), labels=formatC(paste("R²=",round(cor(StrDiv.i.j$StrS, StrDiv.i.j$StrD)^2, 2), sep=""),digits=2), adj=c(0,1), cex=1.2)
	}
}
mtext(side=1, outer=T, text="Number of structural elements (StrS)", line=0, cex=1)
mtext(side=2, outer=T, text="Dissimilarity-based structural diversity (StrD)", line=0.5, cex=1)
mtext(side=2, outer=T, text="Layer", line=4.5, col="gray40", cex=cex.outer)
mtext(side=3, outer=T, text="Transect", line=2, col="gray40", cex=cex.outer)
dev.off()




# Figures showing how structural diversity varies along the transects

props <- Distmax/sum(Distmax)
props2 <- numeric(length=length(props)*2)
props2[c(F,T)] <- props
props2[c(T,F)] <- 0.04

png(filename="figures_paper/StrDivTundra_figure3_StrS.png", height=15, width=20, unit="cm", res=300)

layout(mat=matrix(1:50, nrow=Nlayers, byrow=T), widths=props2)
par(mar=c(1,0,2,0), oma=c(5,6,4,4))

for(i in 1:Nlayers) {
	data.i <- subset(StrDiv, Layer == layer.names[i])
	for(j in 1:Ntransects) {
		data.i.j <- subset (data.i, Transect == transects[j])
		plot(0,type="n",xaxt="n", yaxt="n", xlab="", ylab="",bty="n")
		plot(StrS ~ Distance, data=data.i.j, type="l", lwd=1, ylim=range(data.i$StrS), yaxt=ifelse(j==1,"s","n"))
		if(j == 1) mtext(text=layer.names.plot[i], side=2, line=4.3, col="gray40", cex=cex.outer)
		if(i == 1) mtext(text=transects.plot[j], side=3, line=1, col="gray40", cex=cex.outer)

	}
}
	
mtext(side=2, outer=T, text="Number of structural elements (StrS)", line=1, cex=1)
mtext(side=1, outer=T, text="Distance along transect (m)", line=1.5, cex=1)
mtext(side=2, outer=T, text="Layer", line=4.5, col="gray40", cex=cex.outer)
mtext(side=3, outer=T, text="Transect", line=2, col="gray40", cex=cex.outer)

dev.off()



png(filename="figures_paper/StrDivTundra_figure4_StrD.png", height=15, width=20, unit="cm", res=300)

layout(mat=matrix(1:50, nrow=Nlayers, byrow=T), widths=props2)
par(mar=c(1,0,2,0), oma=c(5,6,4,4))

for(i in 1:Nlayers) {
	data.i <- subset(StrDiv, Layer == layer.names[i])
	for(j in 1:Ntransects) {
		data.i.j <- subset (data.i, Transect == transects[j])
		plot(0,type="n",xaxt="n", yaxt="n", xlab="", ylab="",bty="n")
		plot(StrD ~ Distance, data=data.i.j, type="l", lwd=1, ylim=range(data.i$StrD), yaxt=ifelse(j==1,"s","n"))
		if(j == 1) mtext(text=layer.names.plot[i], side=2, line=4.3, col="gray40", cex=cex.outer)
		if(i == 1) mtext(text=transects.plot[j], side=3, line=1, col="gray40", cex=cex.outer)

	}
}
	
mtext(side=2, outer=T, text="Dissimilarity-based structural diversity (StrD)", line=1, cex=1)
mtext(side=1, outer=T, text="Distance along transect (m)", line=1.5, cex=1)
mtext(side=2, outer=T, text="Layer", line=4.5, col="gray40", cex=cex.outer)
mtext(side=3, outer=T, text="Transect", line=2, col="gray40", cex=cex.outer)

dev.off()






### Question 2a: What are the scales of spatial patterns for the different measures of structural diversity?

### Question 2a.i: How do these patterns deviate from complete serial randomness (CSR)?

# Step 1: Simulate CSR null models for the structural diversity values:
StrDiv.CSR <- list()
#progress <- txtProgressBar(min=1, max=Nlayers * Ntransects * 2, initial=1, char="#", style=3)
k <- 1
for(i in 1:Nlayers) {
	for(j in 1:Ntransects) {
		foo <- subset(StrDiv, Layer==layer.names[i] & Transect==transects[j])
		StrS.CSR <- sim.CSR(foo$StrS, print.loop=F, Nsim=1E4)
		StrDiv.CSR[[k]] <- StrS.CSR
		names(StrDiv.CSR)[k] <- paste(layer.names[i], transects[j], "StrS", sep="_")
		#setTxtProgressBar(progress, k)
		k <- k+1
		print(c(i,j,1))
		StrD.CSR <- sim.CSR(foo$StrD, print.loop=F, Nsim=1E4)
		StrDiv.CSR[[k]] <- StrD.CSR
		names(StrDiv.CSR)[k] <- paste(layer.names[i], transects[j], "StrD", sep="_")
		k <- k+1
		print(c(i,j,2))
	}
}

# Step 2: calculate wavelet scale variance on these data!
### Unless your computer is really really fast, this will take a while.
date.start <- date()														
StrDiv.wav <- list()																																
#progress <- txtProgressBar(min=1, max=length(StrDiv.AR1), initial=1, char="^", style=3)	
proc.time.start <- proc.time()
for(i in 1:length(StrDiv.CSR)) {
	foo <- StrDiv.CSR[[i]]
	bar <- wavCWTvarCIs(foo, pos.var=F, CIquant=c(0.95), print.loop=F, effect.size=F, test=F, scale.per.section=F, scale.max=floor(nrow(foo)/2))
	StrDiv.wav[[i]] <- bar
	print(i)
#	setTxtProgressBar(progress, i)
}
proc.time.spent <- proc.time() - proc.time.start
date.end <- date()
names(StrDiv.wav) <- names(StrDiv.CSR)

save(StrDiv.wav, file="data/StrDivTundra_image_StrDiv_wav.RData")
save.image(file="data/StrDivTundra_image_all.RData")


load("data/StrDivTundra_image_all.RData") # In case you quit R after calculating the wavelets


# Transform the list into a data.frame

Nlines <- numeric()
for(i in 1:length(StrDiv.wav)) {
	Nlines[i] <- length(StrDiv.wav[[i]][[1]])
}
lines.tot <- sum(Nlines)

StrDiv.wav.df <- data.frame(foo=rep(names(StrDiv.wav)[1],length(StrDiv.wav[[1]]$Scales)),  Scales=StrDiv.wav[[1]]$Scales, ScaleVar_Obs = StrDiv.wav[[1]]$ScaleVar_Obs, ScaleVar_CI = StrDiv.wav[[1]]$ScaleVar_CI[1,]) 

for(i in 2:length(StrDiv.wav)) {
	foo <- data.frame(foo=rep(names(StrDiv.wav)[i],length(StrDiv.wav[[i]]$Scales)),  Scales=StrDiv.wav[[i]]$Scales, ScaleVar_Obs = StrDiv.wav[[i]]$ScaleVar_Obs, ScaleVar_CI = StrDiv.wav[[i]]$ScaleVar_CI[1,])
	StrDiv.wav.df <- rbind(StrDiv.wav.df, foo)
}

bar <- strsplit(as.character(StrDiv.wav.df$foo),"_")
bar <- matrix(unlist(bar), ncol=3, byrow=T)

StrDiv.wav.df$Layer <- bar[,1]
StrDiv.wav.df$Transect <- bar[,2]
StrDiv.wav.df$Index <- bar[,3]

StrDiv.wav.df <- StrDiv.wav.df[,c("Layer", "Transect", "Index", "Scales", "ScaleVar_Obs", "ScaleVar_CI")]

# Remove NA
StrDiv.wav.df <- subset(StrDiv.wav.df, !is.na(ScaleVar_Obs))

# Make plots - scale variance, for StrS
scale.range <- range(StrDiv.wav.df$Scales)

png(filename="figures_supplementary/StrDivTundra_StrDivWav_ScaleVar_StrS.png", height=20, width=20, unit="cm", res=300)
par(mfrow=c(Nlayers, Ntransects), mar=c(2,2,1,1), oma=c(5,8,3,3))
for(i in 1:Nlayers) {
	for(j in 1:Ntransects) {
		foo <- subset(StrDiv.wav.df, Layer==layer.names[i] & Transect==transects[j] & Index=="StrS")
		plot(ScaleVar_Obs ~ Scales, type="l", lwd=1, xlim=scale.range, data=foo, main="")
		lines(ScaleVar_CI ~ Scales, lty=2, lwd=1, data=foo, col="gray40")
	}
}
mtext(side=1, text="Scale (m)", outer=T, line=1.5)
mtext(side=2, text="Scale variance", outer=T, line=3.5)
mtext(side=2, outer=T, text=c("Ground", "Herbaceous", "Woody", "Deadwood", "Live \n plants"), at=c(0.90, 0.70, 0.50, 0.30, 0.10), line=1.5)
mtext(side=3, outer=T, text=c("Ecotone 1", "Ecotone 2", "Tundra 1", "Tundra 2", "Tundra 3"), at=c(0.10, 0.30, 0.50, 0.70, 0.90), line=0)
dev.off()


png(filename="figures_supplementary/StrDivTundra_StrDivWav_ScaleVar_StrD.png", height=20, width=20, unit="cm", res=300)
par(mfrow=c(Nlayers, Ntransects), mar=c(2,2,1,1), oma=c(5,8,3,3))
for(i in 1:Nlayers) {
	for(j in 1:Ntransects) {
		foo <- subset(StrDiv.wav.df, Layer==layer.names[i] & Transect==transects[j] & Index=="StrD")
		plot(ScaleVar_Obs ~ Scales, type="l", lwd=1, xlim=scale.range, data=foo, main="")
		lines(ScaleVar_CI ~ Scales, lty=1, lwd=1, data=foo, col="gray40")
	}
}
mtext(side=1, text="Scale (m)", outer=T, line=1.5)
mtext(side=2, text="Scale variance", outer=T, line=3.5)
mtext(side=2, outer=T, text=c("Ground", "Herbaceous", "Woody", "Deadwood", "Live \n plants"), at=c(0.90, 0.70, 0.50, 0.30, 0.10), line=1.5)
mtext(side=3, outer=T, text=c("Ecotone 1", "Ecotone 2", "Tundra 1", "Tundra 2", "Tundra 3"), at=c(0.10, 0.30, 0.50, 0.70, 0.90), line=0)
dev.off()





### Check which scales are significant...

scale.max <- aggregate(Scales~Layer+Transect+Index, data=StrDiv.wav.df, FUN=max)
colnames(scale.max)[4] <- "Scale_max"
StrDiv.wav.df.signif <- subset(StrDiv.wav.df, ScaleVar_Obs > ScaleVar_CI)
StrDiv.wav.df.signif

scales.signif.StrDiv <- aggregate(Scales ~ Layer+Transect+Index, data=StrDiv.wav.df.signif, FUN=paste, collapse="_")
colnames(scales.signif.StrDiv)[4] <- "Scales_signif"

scales.signifRange.StrDiv <- aggregate(Scales ~ Layer+Transect+Index, data=StrDiv.wav.df.signif, FUN=range)
scales.signifRange.StrDiv$Scales_signif_min <- scales.signifRange.StrDiv$Scales[,1]
scales.signifRange.StrDiv$Scales_signif_max <- scales.signifRange.StrDiv$Scales[,2]
scales.signifRange.StrDiv <- scales.signifRange.StrDiv[,-4]

# Which scale correspondos to maximum scale variance?

which.scale.max <- function(x, scales) {
	value <- max(x)
	index <- which(x == value)[1]
	return(scales[index])
}

# Which scales correspond to peaks?

which.scale.peak <- function(x, scales, as.char=F) {
	if(length(x)==0) return(NA)
	if(length(x)==1) return(scales)
	if(length(x)==2) {
		test <- logical(length(x))
		test[1] <- x[1] > x[2]
		test[length(x)] <- x[length(x)] > x[length(x)-1]
		if(as.char) return(paste(scales[test], collapse="_")) else return(scales[test])
	}
	if(length(x)>2) {
		test <- logical(length(x))
		test[1] <- x[1] > x[2]
		for(i in 2:(length(x)-1)) {
			test[i] <- x[i] > x[i-1] & x[i] > x[i+1]
		}
		test[length(x)] <- x[length(x)] > x[length(x)-1]
		if(as.char) return(paste(scales[test], collapse="_")) else return(scales[test])
	}
}

foo.trans <- character()
foo.layer <- character()
foo.index <- character()
foo.scale <- character()
l <- 1
for(i in 1:Nlayers) {
	for(j in 1:Ntransects) {
		for(k in 1:2) {
			data.foo <- subset(StrDiv.wav.df.signif, Layer==layer.names[i] & Transect==transects[j] & Index==c("StrD", "StrS")[k])
			foo.layer[l] <- layer.names[i]
			foo.trans[l] <- as.character(transects)[j]
			foo.index[l] <- c("StrD", "StrS")[k]
			foo.scale[l] <- which.scale.max(data.foo$ScaleVar_Obs, data.foo$Scales)
			l <- l+1
		}
	}
}

scales.max <- data.frame(Layer=foo.layer, Transect=foo.trans, Index=foo.index, Scale.max=foo.scale)

fooo.trans <- character()
fooo.layer <- character()
fooo.index <- character()
fooo.scale <- character()
l <- 1
for(i in 1:Nlayers) {
	for(j in 1:Ntransects) {
		for(k in 1:2) {
			data.fooo <- subset(StrDiv.wav.df.signif, Layer==layer.names[i] & Transect==transects[j] & Index==c("StrD", "StrS")[k])
			fooo.layer[l] <- layer.names[i]
			fooo.trans[l] <- as.character(transects)[j]
			fooo.index[l] <- c("StrD", "StrS")[k]
			fooo.scale[l] <- which.scale.peak(data.fooo$ScaleVar_Obs, data.fooo$Scales, as.char=T)
			l <- l+1
		}
	}
}

scales.peak <- data.frame(Layer=fooo.layer, Transect=fooo.trans, Index=fooo.index, Scale.peak=fooo.scale)


bar <- merge(scale.max, scales.signifRange.StrDiv, all.x=T)
bar <- merge(bar, scales.signif.StrDiv, all.x=T)
bar <- merge(bar, scales.max, all.x=T)
bar <- merge(bar, scales.peak, all.x=T)

scales.signif.StrDiv.df <- bar

write.csv(x=scales.signif.StrDiv.df, file="StrDivTundra_significantScales_StrDiv.csv")


###################################################################################################
### Second part - Relation of structural diversity with the environment
library(car) # to calculate VIF
library(MASS) # to use glmmPQL
library(MuMIn) # to calculate r2 values with 
library(scales) # for transparency in plots
library(nlme) # for mixed effects models

rm(list=ls())


# Functions taken from the help file for pairs:

panel.hist <- function(x, ...) {
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(usr[1:2], 0, 1.5) )
         h <- hist(x, plot = FALSE)
         breaks <- h$breaks; nB <- length(breaks)
         y <- h$counts; y <- y/max(y)
         rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(0, 1, 0, 1))
         r <- abs(cor(x, y)^2)
         txt <- format(c(r, 0.123456789), digits = digits)[1]
         txt <- paste0(prefix, txt)
         if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
         text(0.5, 0.5, txt, cex = cex.cor * r)
}

setwd("/home/pavel/Profissional/Pesquisa/MyPapers/2014_Pavel_StrDivTundra")

data.StrDiv <- read.table("data/StrDivTundra_data_StrDiv.txt", header=T, stringsAsFactors=F)
data.expl <- read.table("data/StrDivTundra_data_explanatory.txt", header=T, stringsAsFactor=F)
data.full <- merge(data.StrDiv, data.expl, all.x=T, all.y=T)

str(data.StrDiv)
str(data.expl)
str(data.full)

### Remove rows with no pH data

data.full <- subset(data.full, !is.na(pH))
str(data.full)

### Remove disturbances (i.e. roads)

data.full <- subset(data.full, Disturbance=="No")
str(data.full)

### Check distribution of each explanatory variable
par(mfrow=c(2,2))
hist(data.full$Altitude)
hist(data.full$DistFromLake)
hist(data.full$pH)
hist(data.full$Microtopography)

### The are two plots with high Microtopography values. On closer examination, they appear to correspond to embankments near one edge of the Tundra1 transect and it seems better to remove them as outliers.

data.full <- subset(data.full, Microtopography < 30)

### Re-check distribution of each explanatory variable
par(mfrow=c(2,2))
hist(data.full$Altitude)
hist(data.full$DistFromLake)
hist(data.full$pH)
hist(data.full$Microtopography)

### Although there are still plots with higher microtopography than the others, there are no clear ecological reasons to remove them.

### Reclassify the levels of Transect

### Check how many quadrats have structural elements of each type, in each transect

StrS.exists <- aggregate(StrS ~ Transect + Layer, data=data.full, FUN=function(x) sum(x>0)/length(x))
StrD.exists <- aggregate(StrD ~ Transect + Layer, data=data.full, FUN=function(x) sum(x>0)/length(x))

StrS.exists
StrD.exists

# I think there's a sufficient number of non-zero values to permit the analyses...

# Check which transects are present

unique(data.full$Transect)

# Transform the Layer column into factor

data.full$Layer <- as.factor(data.full$Layer)


### Perform separate analyses for tundra and ecotone

data.full.tundra <- subset(data.full, Transect %in% c("Tundra1", "Tundra2"))
data.full.ecotone <- subset(data.full, Transect %in% c("Ecotone1", "Ecotone2"))

data.full.tundra$Transect <- as.factor(data.full.tundra$Transect)
data.full.ecotone$Transect <- as.factor(data.full.ecotone$Transect)


### Exploratory figures
layer.names <- as.character(levels(data.full$Layer))
Nlayers <- length(layer.names)

colors.trans <- alpha(c("Red", "Orange", "Blue", "Purple"), 0.5)

png(filename="figures_exploratory/StrDivTundra_StrSvsEnvironment_Tundra.png", height=30, width=20, unit="cm", res=300)

par(mfrow=c(5,4), mar=c(4,3,3,2), oma=c(2,2,2,3))
for(i in 1:Nlayers) {
	data.temp <- subset(data.full.tundra, Layer==layer.names[i])
	plot(StrS ~ Altitude, data=data.temp, pch=21, bg=colors.trans[as.numeric(Transect)+2], main=paste(layer.names[i], " StrS vs \n Altitude", sep=""))
	plot(StrS ~ DistFromLake, data=data.temp, pch=21, bg=colors.trans[as.numeric(Transect)+2], main=paste(layer.names[i], " StrS vs \n Distance", sep=""))
	plot(StrS ~ pH, data=data.temp, pch=21, bg=colors.trans[as.numeric(Transect)+2], main=paste(layer.names[i], " StrS vs \n pH", sep=""))
	plot(StrS ~ Microtopography, data=data.temp, pch=21, bg=colors.trans[as.numeric(Transect)+2], main=paste(layer.names[i], " StrS vs \n Microtopography", sep=""))
	mtext(text="StrS", side=2, outer=T)
	mtext(text="StrS - Tundra", side=3, outer=T)
}
dev.off()

png(filename="figures_exploratory/StrDivTundra_StrDvsEnvironment_Tundra.png", height=30, width=20, unit="cm", res=300)

par(mfrow=c(5,4), mar=c(4,3,3,2), oma=c(2,2,2,3))
for(i in 1:Nlayers) {
	data.temp <- subset(data.full.tundra, Layer==layer.names[i])
	plot(StrD ~ Altitude, data=data.temp, pch=21, bg=colors.trans[as.numeric(Transect)+2], main=paste(layer.names[i], " StrD vs \n Altitude", sep=""))
	plot(StrD ~ DistFromLake, data=data.temp, pch=21, bg=colors.trans[as.numeric(Transect)+2], main=paste(layer.names[i], " StrD vs \n Distance", sep=""))
	plot(StrD ~ pH, data=data.temp, pch=21, bg=colors.trans[as.numeric(Transect)+2], main=paste(layer.names[i], " StrD vs \n pH", sep=""))
	plot(StrD ~ Microtopography, data=data.temp, pch=21, bg=colors.trans[as.numeric(Transect)+2], main=paste(layer.names[i], " StrD vs \n Microtopography", sep=""))
	mtext(text="StrD", side=2, outer=T)
	mtext(text="StrD - Tundra ", side=3, outer=T)
}
dev.off()


png(filename="figures_exploratory/StrDivTundra_StrSvsEnvironment_Ecotone.png", height=30, width=20, unit="cm", res=300)

par(mfrow=c(5,4), mar=c(4,3,3,2), oma=c(2,2,2,3))
for(i in 1:Nlayers) {
	data.temp <- subset(data.full.ecotone, Layer==layer.names[i])
	plot(StrS ~ Altitude, data=data.temp, pch=21, bg=colors.trans[as.numeric(Transect)], main=paste(layer.names[i], " StrS vs \n Altitude", sep=""))
	plot(StrS ~ DistFromLake, data=data.temp, pch=21, bg=colors.trans[as.numeric(Transect)], main=paste(layer.names[i], " StrS vs \n Distance", sep=""))
	plot(StrS ~ pH, data=data.temp, pch=21, bg=colors.trans[as.numeric(Transect)], main=paste(layer.names[i], " StrS vs \n pH", sep=""))
	plot(StrS ~ Microtopography, data=data.temp, pch=21, bg=colors.trans[as.numeric(Transect)], main=paste(layer.names[i], " StrS vs \n Microtopography", sep=""))
	mtext(text="StrS", side=2, outer=T)
	mtext(text="StrS - Ecotone", side=3, outer=T)
}
dev.off()



png(filename="figures_exploratory/StrDivTundra_StrDvsEnvironment_Ecotone.png", height=30, width=20, unit="cm", res=300)

par(mfrow=c(5,4), mar=c(4,3,3,2), oma=c(2,2,2,3))
for(i in 1:Nlayers) {
	data.temp <- subset(data.full.ecotone, Layer==layer.names[i])
	plot(StrD ~ Altitude, data=data.temp, pch=21, bg=colors.trans[as.numeric(Transect)], main=paste(layer.names[i], " StrD vs \n Altitude", sep=""))
	plot(StrD ~ DistFromLake, data=data.temp, pch=21, bg=colors.trans[as.numeric(Transect)], main=paste(layer.names[i], " StrD vs \n Distance", sep=""))
	plot(StrD ~ pH, data=data.temp, pch=21, bg=colors.trans[as.numeric(Transect)], main=paste(layer.names[i], " StrD vs \n pH", sep=""))
	plot(StrD ~ Microtopography, data=data.temp, pch=21, bg=colors.trans[as.numeric(Transect)], main=paste(layer.names[i], " StrD vs \n Microtopography", sep=""))
	mtext(text="StrD", side=2, outer=T)
	mtext(text="StrD - Ecotone", side=3, outer=T)
}
dev.off()



### Create objects for layer and vegetation type

data.full.tundra.ground <- subset(data.full.tundra, Layer=="Ground")
data.full.tundra.deadwood <- subset(data.full.tundra, Layer=="Deadwood")
data.full.tundra.herbaceous <- subset(data.full.tundra, Layer=="Herbaceous")
data.full.tundra.woody <- subset(data.full.tundra, Layer=="Woody")
data.full.tundra.liveplants <- subset(data.full.tundra, Layer=="LivePlants")
data.full.ecotone.ground <- subset(data.full.ecotone, Layer=="Ground")
data.full.ecotone.deadwood <- subset(data.full.ecotone, Layer=="Deadwood")
data.full.ecotone.herbaceous <- subset(data.full.ecotone, Layer=="Herbaceous")
data.full.ecotone.woody <- subset(data.full.ecotone, Layer=="Woody")
data.full.ecotone.liveplants <- subset(data.full.ecotone, Layer=="LivePlants")

### Check for colinearity in the tundra and the ecotone

mod.temp.tundra <- lm(StrS ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.tundra.ground)
vif(mod.temp.tundra) # all VIF values < 3. No variable has to be removed.
mod.temp.ecotone <- lm(StrS ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.ecotone.ground)
vif(mod.temp.ecotone) # all VIF values < 2. No variable has to be removed.


###### Adjusting the models... using glmmPQL from MASS package, with stepwise selection based on p-values, and r.squaredGLMM from MuMIn package to calculate r2 of the final model. Critical p-value of 0.01. For the pseudo-R2: use marginal (m) values and the trigamma method when available.

### Tundra - StrS - ground
mod.StrS.tundra.ground.glmm.0 <- glmmPQL(StrS ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.tundra.ground, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrS.tundra.ground.glmm.0) # remove DistFromLake
mod.StrS.tundra.ground.glmm.1 <- glmmPQL(StrS ~ Altitude + pH + Microtopography, data=data.full.tundra.ground, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrS.tundra.ground.glmm.1) # remove pH
mod.StrS.tundra.ground.glmm.2 <- glmmPQL(StrS ~ Altitude + Microtopography, data=data.full.tundra.ground, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrS.tundra.ground.glmm.2) # remove Altitude
mod.StrS.tundra.ground.glmm.3 <- glmmPQL(StrS ~ Microtopography, data=data.full.tundra.ground, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrS.tundra.ground.glmm.3) # no significant variables
r.squaredGLMM(mod.StrS.tundra.ground.glmm.3)
mod.StrS.tundra.ground.glmm.final <- glmmPQL(StrS ~ 1, data=data.full.tundra.ground, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)

### Tundra - StrD - ground
mod.StrD.tundra.ground.glmm.0 <- glmmPQL(StrD ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.tundra.ground, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.tundra.ground.glmm.0) # remove DistFromLake
mod.StrD.tundra.ground.glmm.1 <- glmmPQL(StrD ~ Altitude + pH + Microtopography, data=data.full.tundra.ground, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.tundra.ground.glmm.1) # remove Altitude
mod.StrD.tundra.ground.glmm.2 <- glmmPQL(StrD ~ pH + Microtopography, data=data.full.tundra.ground, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.tundra.ground.glmm.2) # remove pH
mod.StrD.tundra.ground.glmm.3 <- glmmPQL(StrD ~ Microtopography, data=data.full.tundra.ground, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.tundra.ground.glmm.3) # Nothing significant
r.squaredGLMM(mod.StrD.tundra.ground.glmm.3)


### Ecotone - StrS - ground
mod.StrS.ecotone.ground.glmm.0 <- glmmPQL(StrS ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.ecotone.ground, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrS.ecotone.ground.glmm.0) # remove pH
mod.StrS.ecotone.ground.glmm.1 <- glmmPQL(StrS ~ Altitude + DistFromLake + Microtopography, data=data.full.ecotone.ground, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrS.ecotone.ground.glmm.1) # remove Altitude
mod.StrS.ecotone.ground.glmm.2 <- glmmPQL(StrS ~ DistFromLake + Microtopography, data=data.full.ecotone.ground, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrS.ecotone.ground.glmm.2) # DistFromLake and Microtopography significant
mod.StrS.ecotone.ground.glmm.final <- mod.StrS.ecotone.ground.glmm.2
r.squaredGLMM(mod.StrS.ecotone.ground.glmm.final) # R2 = 0.09. So, there is a significant relationship but it doesn't really mean much :-)

### Ecotone - StrD - ground
mod.StrD.ecotone.ground.glmm.0 <- glmmPQL(StrD ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.ecotone.ground, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.ecotone.ground.glmm.0) # remove Altitude
mod.StrD.ecotone.ground.glmm.1 <- glmmPQL(StrD ~ DistFromLake + pH + Microtopography, data=data.full.ecotone.ground, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.ecotone.ground.glmm.1) # remove pH
mod.StrD.ecotone.ground.glmm.2 <- glmmPQL(StrD ~ DistFromLake +  Microtopography, data=data.full.ecotone.ground, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.ecotone.ground.glmm.2) # remove DistFromLake
mod.StrD.ecotone.ground.glmm.3 <- glmmPQL(StrD ~  Microtopography, data=data.full.ecotone.ground, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.ecotone.ground.glmm.3) # Microtopography significant
r.squaredGLMM(mod.StrD.ecotone.ground.glmm.3)

### Tundra - StrS - Herbaceous
mod.StrS.tundra.herbaceous.glmm.0 <- glmmPQL(StrS ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.tundra.herbaceous, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrS.tundra.herbaceous.glmm.0) # remove Altitude
mod.StrS.tundra.herbaceous.glmm.1 <- glmmPQL(StrS ~ DistFromLake + pH + Microtopography, data=data.full.tundra.herbaceous, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrS.tundra.herbaceous.glmm.1) # remove DistFromLake
mod.StrS.tundra.herbaceous.glmm.2 <- glmmPQL(StrS ~ pH + Microtopography, data=data.full.tundra.herbaceous, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrS.tundra.herbaceous.glmm.2) # all significant
mod.StrS.tundra.herbaceous.glmm.final <- mod.StrS.tundra.herbaceous.glmm.2
r.squaredGLMM(mod.StrS.tundra.herbaceous.glmm.final) #

### Tundra - StrD - Herbaceous
mod.StrD.tundra.herbaceous.glmm.0 <- glmmPQL(StrD ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.tundra.herbaceous, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.tundra.herbaceous.glmm.0) # remove Atitude
mod.StrD.tundra.herbaceous.glmm.1 <- glmmPQL(StrD ~ DistFromLake + pH + Microtopography, data=data.full.tundra.herbaceous, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.tundra.herbaceous.glmm.1) # remove DistFromLake
mod.StrD.tundra.herbaceous.glmm.2 <- glmmPQL(StrD ~ pH + Microtopography, data=data.full.tundra.herbaceous, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.tundra.herbaceous.glmm.2) # remove pH
mod.StrD.tundra.herbaceous.glmm.3 <- glmmPQL(StrD ~ Microtopography, data=data.full.tundra.herbaceous, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.tundra.herbaceous.glmm.3) # remove pH
mod.StrD.tundra.herbaceous.glmm.final <- mod.StrD.tundra.herbaceous.glmm.3
r.squaredGLMM(mod.StrD.tundra.herbaceous.glmm.final) # 

### Ecotone - StrS - Herbaceous
mod.StrS.ecotone.herbaceous.glmm.0 <- glmmPQL(StrS ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.ecotone.herbaceous, family=poisson, correlation=corAR1(form=~Distance|
	Transect), random=~1|Transect)
summary(mod.StrS.ecotone.herbaceous.glmm.0) # remove Altitude
mod.StrS.ecotone.herbaceous.glmm.1 <- glmmPQL(StrS ~ DistFromLake + pH + Microtopography, data=data.full.ecotone.herbaceous, family=poisson, correlation=corAR1(form=~Distance|
	Transect), random=~1|Transect)
summary(mod.StrS.ecotone.herbaceous.glmm.1) # remove pH
mod.StrS.ecotone.herbaceous.glmm.2 <- glmmPQL(StrS ~ DistFromLake + Microtopography, data=data.full.ecotone.herbaceous, family=poisson, correlation=corAR1(form=~Distance|
	Transect), random=~1|Transect)
summary(mod.StrS.ecotone.herbaceous.glmm.2) # remove Microtopography
mod.StrS.ecotone.herbaceous.glmm.3 <- glmmPQL(StrS ~ DistFromLake, data=data.full.ecotone.herbaceous, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrS.ecotone.herbaceous.glmm.3) # not significant
mod.StrS.ecotone.herbaceous.glmm.final <-  mod.StrS.ecotone.herbaceous.glmm.3
r.squaredGLMM(mod.StrS.ecotone.herbaceous.glmm.final) # R2 of 0.02 - not meaningful.

### Ecotone - StrD - Herbaceous
mod.StrD.ecotone.herbaceous.glmm.0 <- glmmPQL(StrD ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.ecotone.herbaceous, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.ecotone.herbaceous.glmm.0) # remove Altitude
mod.StrD.ecotone.herbaceous.glmm.1 <- glmmPQL(StrD ~ DistFromLake + pH + Microtopography, data=data.full.ecotone.herbaceous, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.ecotone.herbaceous.glmm.1) # remove pH
mod.StrD.ecotone.herbaceous.glmm.2 <- glmmPQL(StrD ~ DistFromLake + Microtopography, data=data.full.ecotone.herbaceous, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.ecotone.herbaceous.glmm.2) # remove Microtopography
mod.StrD.ecotone.herbaceous.glmm.3 <- glmmPQL(StrD ~  DistFromLake, data=data.full.ecotone.herbaceous, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.ecotone.herbaceous.glmm.3) # nothing significant
r.squaredGLMM(mod.StrD.ecotone.herbaceous.glmm.3)


### Tundra - StrS - Woody
mod.StrS.tundra.woody.glmm.0 <- glmmPQL(StrS ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.tundra.woody, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrS.tundra.woody.glmm.0) # remove Altitude
mod.StrS.tundra.woody.glmm.1 <- glmmPQL(StrS ~ DistFromLake + pH + Microtopography, data=data.full.tundra.woody, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrS.tundra.woody.glmm.1) # remove Microtopography
mod.StrS.tundra.woody.glmm.2 <- glmmPQL(StrS ~ DistFromLake + pH, data=data.full.tundra.woody, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrS.tundra.woody.glmm.2) # remove DistFromLake
mod.StrS.tundra.woody.glmm.3 <- glmmPQL(StrS ~ pH, data=data.full.tundra.woody, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrS.tundra.woody.glmm.3) # nothing significant
r.squaredGLMM(mod.StrS.tundra.woody.glmm.3) # r2 of 0.02 - not meaningful

### Tundra - StrD - Woody
mod.StrD.tundra.woody.glmm.0 <- glmmPQL(StrD ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.tundra.woody, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.tundra.woody.glmm.0) # remove Microtopography
mod.StrD.tundra.woody.glmm.1 <- glmmPQL(StrD ~ Altitude + DistFromLake + pH, data=data.full.tundra.woody, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.tundra.woody.glmm.1) # remove Altitude
mod.StrD.tundra.woody.glmm.2 <- glmmPQL(StrD ~ DistFromLake + pH, data=data.full.tundra.woody, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.tundra.woody.glmm.2) # remove DistFromLake
mod.StrD.tundra.woody.glmm.3 <- glmmPQL(StrD ~ pH, data=data.full.tundra.woody, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.tundra.woody.glmm.3) # nothing significant
r.squaredGLMM(mod.StrD.tundra.woody.glmm.3)

### Ecotone - StrS - Woody
mod.StrS.ecotone.woody.glmm.0 <- glmmPQL(StrS ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.ecotone.woody, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrS.ecotone.woody.glmm.0) # remove Altitude
mod.StrS.ecotone.woody.glmm.1 <- glmmPQL(StrS ~ DistFromLake + pH + Microtopography, data=data.full.ecotone.woody, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrS.ecotone.woody.glmm.1) # remove DistFromLake
mod.StrS.ecotone.woody.glmm.2 <- glmmPQL(StrS ~ pH + Microtopography, data=data.full.ecotone.woody, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrS.ecotone.woody.glmm.2) # remove microtopography
mod.StrS.ecotone.woody.glmm.3 <- glmmPQL(StrS ~ pH, data=data.full.ecotone.woody, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrS.ecotone.woody.glmm.3) #
mod.StrS.ecotone.woody.glmm.final <- mod.StrS.ecotone.woody.glmm.3
r.squaredGLMM(mod.StrS.ecotone.woody.glmm.final) # r squared of 0.07.

### Ecotone - StrD - Woody
mod.StrD.ecotone.woody.glmm.0 <- glmmPQL(StrD ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.ecotone.woody, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.ecotone.woody.glmm.0) # Remove altitude
mod.StrD.ecotone.woody.glmm.1 <- glmmPQL(StrD ~ DistFromLake + pH + Microtopography, data=data.full.ecotone.woody, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.ecotone.woody.glmm.1) # Remove DistFromLake
mod.StrD.ecotone.woody.glmm.2 <- glmmPQL(StrD ~ pH + Microtopography, data=data.full.ecotone.woody, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.ecotone.woody.glmm.2) # Remove Microtopography
mod.StrD.ecotone.woody.glmm.3 <- glmmPQL(StrD ~ pH, data=data.full.ecotone.woody, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.ecotone.woody.glmm.3) # pH significant
mod.StrD.ecotone.woody.glmm.final <- mod.StrD.ecotone.woody.glmm.3
r.squaredGLMM(mod.StrD.ecotone.woody.glmm.final) # R2 of 0.09

### Tundra - StrS - Deadwood
mod.StrS.tundra.deadwood.glmm.0 <- glmmPQL(StrS ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.tundra.deadwood, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrS.tundra.deadwood.glmm.0) # Remove Altitude
mod.StrS.tundra.deadwood.glmm.1 <- glmmPQL(StrS ~ DistFromLake + pH + Microtopography, data=data.full.tundra.deadwood, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrS.tundra.deadwood.glmm.1) # Remove Microtopography
mod.StrS.tundra.deadwood.glmm.2 <- glmmPQL(StrS ~ DistFromLake + pH, data=data.full.tundra.deadwood, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrS.tundra.deadwood.glmm.2) # Remove DistFromLake
mod.StrS.tundra.deadwood.glmm.3 <- glmmPQL(StrS ~ pH, data=data.full.tundra.deadwood, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrS.tundra.deadwood.glmm.3) # pH not significant
r.squaredGLMM(mod.StrS.tundra.deadwood.glmm.3)




### Tundra - StrD - Deadwood
mod.StrD.tundra.deadwood.glmm.0 <- glmmPQL(StrD ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.tundra.deadwood, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.tundra.deadwood.glmm.0) # Remove Altitude
mod.StrD.tundra.deadwood.glmm.1 <- glmmPQL(StrD ~ DistFromLake + pH + Microtopography, data=data.full.tundra.deadwood, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.tundra.deadwood.glmm.1) # Remove Microtopography
mod.StrD.tundra.deadwood.glmm.2 <- glmmPQL(StrD ~ DistFromLake + pH, data=data.full.tundra.deadwood, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.tundra.deadwood.glmm.2) # remove pH
mod.StrD.tundra.deadwood.glmm.3 <- glmmPQL(StrD ~ DistFromLake, data=data.full.tundra.deadwood, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.tundra.deadwood.glmm.3) # non significant (p>0.01)
r.squaredGLMM(mod.StrD.tundra.deadwood.glmm.3)


### Ecotone - StrS - Deadwood
mod.StrS.ecotone.deadwood.glmm.0 <- glmmPQL(StrS ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.ecotone.deadwood, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrS.ecotone.deadwood.glmm.0) # Remove DistFromLake
mod.StrS.ecotone.deadwood.glmm.1 <- glmmPQL(StrS ~ Altitude + pH + Microtopography, data=data.full.ecotone.deadwood, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrS.ecotone.deadwood.glmm.1) # Remove Altitude
mod.StrS.ecotone.deadwood.glmm.2 <- glmmPQL(StrS ~ pH + Microtopography, data=data.full.ecotone.deadwood, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrS.ecotone.deadwood.glmm.2) # pH and microtopography significant
mod.StrS.ecotone.deadwood.glmm.final <- mod.StrS.ecotone.deadwood.glmm.2
r.squaredGLMM(mod.StrS.ecotone.deadwood.glmm.final) # R2 of 0.04


### Ecotone - StrD - Deadwood
mod.StrD.ecotone.deadwood.glmm.0 <- glmmPQL(StrD ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.ecotone.deadwood, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.ecotone.deadwood.glmm.0) # Remove DistFromLake
mod.StrD.ecotone.deadwood.glmm.1 <- glmmPQL(StrD ~ Altitude + pH + Microtopography, data=data.full.ecotone.deadwood, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.ecotone.deadwood.glmm.1) # Remove Altitude
mod.StrD.ecotone.deadwood.glmm.2 <- glmmPQL(StrD ~ pH + Microtopography, data=data.full.ecotone.deadwood, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.ecotone.deadwood.glmm.2) # Remove microtopography
mod.StrD.ecotone.deadwood.glmm.3 <- glmmPQL(StrD ~ pH, data=data.full.ecotone.deadwood, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.ecotone.deadwood.glmm.3) # All significant
mod.StrD.ecotone.deadwood.glmm.final <- mod.StrD.ecotone.deadwood.glmm.3
r.squaredGLMM(mod.StrD.ecotone.deadwood.glmm.final) # R2 of 0.09

### Tundra - StrS - LivePlants
mod.StrS.tundra.liveplants.glmm.0 <- glmmPQL(StrS ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.tundra.liveplants, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrS.tundra.liveplants.glmm.0) # Remove Altitude
mod.StrS.tundra.liveplants.glmm.1 <- glmmPQL(StrS ~ DistFromLake + pH + Microtopography, data=data.full.tundra.liveplants, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrS.tundra.liveplants.glmm.1) # Remove DistFromLake
mod.StrS.tundra.liveplants.glmm.2 <- glmmPQL(StrS ~ pH + Microtopography, data=data.full.tundra.liveplants, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrS.tundra.liveplants.glmm.2) # pH and Micropography significant
mod.StrS.tundra.liveplants.glmm.final <- mod.StrS.tundra.liveplants.glmm.2
r.squaredGLMM(mod.StrS.tundra.liveplants.glmm.final) # R2 of 0.04

### Tundra - StrD - LivePlants
mod.StrD.tundra.liveplants.glmm.0 <- glmmPQL(StrD ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.tundra.liveplants, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.tundra.liveplants.glmm.0) # Remove DistFromaLake
mod.StrD.tundra.liveplants.glmm.1 <- glmmPQL(StrD ~ Altitude + pH + Microtopography, data=data.full.tundra.liveplants, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.tundra.liveplants.glmm.1) # Remove Altitude
mod.StrD.tundra.liveplants.glmm.2 <- glmmPQL(StrD ~ pH + Microtopography, data=data.full.tundra.liveplants, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.tundra.liveplants.glmm.2) # Remove microtopography
mod.StrD.tundra.liveplants.glmm.3 <- glmmPQL(StrD ~ pH, data=data.full.tundra.liveplants, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.tundra.liveplants.glmm.3) # pH and Microtopography significant
mod.StrD.tundra.liveplants.glmm.final <- mod.StrD.tundra.liveplants.glmm.2
r.squaredGLMM(mod.StrD.tundra.liveplants.glmm.final) # R2 of 0.10


### Ecotone - StrS - LivePlants
mod.StrS.ecotone.liveplants.glmm.0 <- glmmPQL(StrS ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.ecotone.liveplants, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrS.ecotone.liveplants.glmm.0) # Remove Altitude
mod.StrS.ecotone.liveplants.glmm.1 <- glmmPQL(StrS ~ DistFromLake + pH + Microtopography, data=data.full.ecotone.liveplants, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrS.ecotone.liveplants.glmm.1) # remove pH
mod.StrS.ecotone.liveplants.glmm.2 <- glmmPQL(StrS ~ DistFromLake + Microtopography, data=data.full.ecotone.liveplants, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrS.ecotone.liveplants.glmm.2) # remove microtopography
mod.StrS.ecotone.liveplants.glmm.3 <- glmmPQL(StrS ~ DistFromLake, data=data.full.ecotone.liveplants, family=poisson, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrS.ecotone.liveplants.glmm.3) # remove microtopography
r.squaredGLMM(mod.StrS.ecotone.liveplants.glmm.3) # R2 of 0.05

### Ecotone - StrD - LivePlants
mod.StrD.ecotone.liveplants.glmm.0 <- glmmPQL(StrD ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.ecotone.liveplants, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.ecotone.liveplants.glmm.0) # Remove DistFromLake
mod.StrD.ecotone.liveplants.glmm.1 <- glmmPQL(StrD ~ Altitude + pH + Microtopography, data=data.full.ecotone.liveplants, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.ecotone.liveplants.glmm.1) # Remove pH
mod.StrD.ecotone.liveplants.glmm.2 <- glmmPQL(StrD ~ Altitude + Microtopography, data=data.full.ecotone.liveplants, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.ecotone.liveplants.glmm.2) # Remove altitude
mod.StrD.ecotone.liveplants.glmm.3 <- glmmPQL(StrD ~ Microtopography, data=data.full.ecotone.liveplants, family=gaussian, correlation=corAR1(form=~Distance|Transect), random=~1|Transect)
summary(mod.StrD.ecotone.liveplants.glmm.3) # Altitude and Microtopography significant
mod.StrD.ecotone.liveplants.glmm.final <- mod.StrD.ecotone.liveplants.glmm.3
r.squaredGLMM(mod.StrD.ecotone.liveplants.glmm.final) # R2 of 0.06


### Present results as table: response variables as rows; explanatory variables as columns; coefficients of the final model in the cells; R2 in the final column.





### Figures for supplementary material


### Figures showing the distribution and height of different elements along the transects

for(i in 1:Ntransects) {
	data.all.i <- subset(data.all, Transect==transects[[i]])
	distmax.i <- max(data.all.i$Distance)
	foo <- rep(1:Nelements, length.out=distmax.i)
	png(filename=paste("figures_supplementary/StrDivTundra_Elements_Quantity_", transects[i], ".png", sep=""), height=20, width=20, unit="cm", res=300)
	par(mar=c(5,11,3,3))
	plot(foo ~ c(1:distmax.i), type="n", ylab="", yaxt="n", xlab="Distance (m)", main=paste("Quantity - ",transects[i], sep=""))
	axis(side=2, at=Nelements:1, labels=elements, las=1, tick=F)
	abline(h=1:Nelements, col="gray70")
	for(j in 1:Nelements) {
		data.all.j <- subset(data.all.i, Element_detail==elements[j], select=c("Element_detail","Distance", "Quantity"))
		y.temp <- rep(Nelements+1-j, nrow(data.all.j))
		points(y.temp ~ data.all.j$Distance, pch="I", cex=data.all.j$Quantity/2)
	}
	dev.off()
}


for(i in 1:Ntransects) {
	data.all.i <- subset(data.all, Transect==transects[[i]])
	distmax.i <- max(data.all.i$Distance)
	foo <- rep(1:Nelements, length.out=distmax.i)
	png(filename=paste("figures_supplementary/StrDivTundra_Elements_Height_", transects[i], ".png", sep=""), height=20, width=20, unit="cm", res=300)
	par(mar=c(5,11,3,3))
	plot(foo ~ c(1:distmax.i), type="n", ylab="", yaxt="n", xlab="Distance (m)", main=paste("Height - ", transects[i], sep=""))
	axis(side=2, at=Nelements:1, labels=elements, las=1, tick=F)
	abline(h=1:Nelements, col="gray70")
	for(j in 1:Nelements) {
		data.all.j <- subset(data.all.i, Element_detail==elements[j], select=c("Element_detail","Distance", "Height_tot"))
		y.temp <- rep(Nelements+1-j, nrow(data.all.j))
		points(y.temp ~ data.all.j$Distance, pch="I", cex=data.all.j$Height_tot/2)
	}
	dev.off()
}





