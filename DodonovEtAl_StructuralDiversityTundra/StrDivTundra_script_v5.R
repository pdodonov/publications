# The package wmtsa, essential for these analyses, is no longer available at CRAN. In ortder to install it, use the following script:

#library(devtools)
#install_url('https://cran.r-project.org/src/contrib/Archive/wmtsa/wmtsa_2.0-3.tar.gz')

# Study areas forest-tundra ecotone and tundra, near Churchill, MB, Canada. The ecotone has two transects separated by a lake. The tundra has three transects separated by lakes and with a road crossing one of them.
# Objective 1: To assess whether structural diversity may be explained by environmental factors, namely soil pH, microtopography, and distance from edges.
# Objective 2a: To assess the spatial scales of structural diversity
# Objective 2b: To compare the spatial scales of structural diveristy between the environments
# Objective 2c: To compare the spatial scales of structural diveristy between structural diversity components. 


##### First part - scales of spatial pattern


### Load required packages
library(vegan) # used in the function to calculate structural diversity
library(cluster) # used in the function to calculate structural diversity
library(devtools)
library(wmtsa) # required for the wavelet analysis
library(scales) # to make plots with transparency
library(ecodist) # for cluster analysis required to calculate the structural diversity measure
          




### Load custom functions

setwd("/home/pavel/Profissional/Pesquisa/MyPapers_working/Pavel_2014_StrDivTundra")

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
#library(MASS) # to use glmmPQL - no longer used
library(MuMIn) # used mainly for the dredge function
library(scales) # for transparency in plots
#library(nlme) # for mixed effects models - no longer used
#library(bbmle) # for AIC - no longer used
library(geepack) # for generalized estimation equations
library(rmarkdown) # to export some of the results as PDF
library(rsq) # to calculate pseudo-R2 values

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

setwd("/home/pavel/Profissional/Pesquisa/MyPapers_working/Pavel_2014_StrDivTundra")

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

# The code below is marked as comments in order to not remake the figures every time it is run.

#png(filename="figures_exploratory/StrDivTundra_StrSvsEnvironment_Tundra.png", height=30, width=20, unit="cm", res=300)
#
#par(mfrow=c(5,4), mar=c(4,3,3,2), oma=c(2,2,2,3))
#for(i in 1:Nlayers) {
#   data.temp <- subset(data.full.tundra, Layer==layer.names[i])
#   plot(StrS ~ Altitude, data=data.temp, pch=21, bg=colors.trans[as.numeric(Transect)+2], main=paste(layer.names[i], " StrS vs \n Altitude", sep=""))
#   plot(StrS ~ DistFromLake, data=data.temp, pch=21, bg=colors.trans[as.numeric(Transect)+2], main=paste(layer.names[i], " StrS vs \n Distance", sep=""))
#   plot(StrS ~ pH, data=data.temp, pch=21, bg=colors.trans[as.numeric(Transect)+2], main=paste(layer.names[i], " StrS vs \n pH", sep=""))
#   plot(StrS ~ Microtopography, data=data.temp, pch=21, bg=colors.trans[as.numeric(Transect)+2], main=paste(layer.names[i], " StrS vs \n Microtopography", sep=""))
#   mtext(text="StrS", side=2, outer=T)
#   mtext(text="StrS - Tundra", side=3, outer=T)
#}
#dev.off()
#
#png(filename="figures_exploratory/StrDivTundra_StrDvsEnvironment_Tundra.png", height=30, width=20, unit="cm", res=300)
#
#par(mfrow=c(5,4), mar=c(4,3,3,2), oma=c(2,2,2,3))
#for(i in 1:Nlayers) {
#   data.temp <- subset(data.full.tundra, Layer==layer.names[i])
#   plot(StrD ~ Altitude, data=data.temp, pch=21, bg=colors.trans[as.numeric(Transect)+2], main=paste(layer.names[i], " StrD vs \n Altitude", sep=""))
#   plot(StrD ~ DistFromLake, data=data.temp, pch=21, bg=colors.trans[as.numeric(Transect)+2], main=paste(layer.names[i], " StrD vs \n Distance", sep=""))
#   plot(StrD ~ pH, data=data.temp, pch=21, bg=colors.trans[as.numeric(Transect)+2], main=paste(layer.names[i], " StrD vs \n pH", sep=""))
#   plot(StrD ~ Microtopography, data=data.temp, pch=21, bg=colors.trans[as.numeric(Transect)+2], main=paste(layer.names[i], " StrD vs \n Microtopography", sep=""))
#   mtext(text="StrD", side=2, outer=T)
#   mtext(text="StrD - Tundra ", side=3, outer=T)
#}
#dev.off()
#
#
#png(filename="figures_exploratory/StrDivTundra_StrSvsEnvironment_Ecotone.png", height=30, width=20, unit="cm", res=300)
#
#par(mfrow=c(5,4), mar=c(4,3,3,2), oma=c(2,2,2,3))
#for(i in 1:Nlayers) {
#   data.temp <- subset(data.full.ecotone, Layer==layer.names[i])
#   plot(StrS ~ Altitude, data=data.temp, pch=21, bg=colors.trans[as.numeric(Transect)], main=paste(layer.names[i], " StrS vs \n Altitude", sep=""))
#   plot(StrS ~ DistFromLake, data=data.temp, pch=21, bg=colors.trans[as.numeric(Transect)], main=paste(layer.names[i], " StrS vs \n Distance", sep=""))
#   plot(StrS ~ pH, data=data.temp, pch=21, bg=colors.trans[as.numeric(Transect)], main=paste(layer.names[i], " StrS vs \n pH", sep=""))
#   plot(StrS ~ Microtopography, data=data.temp, pch=21, bg=colors.trans[as.numeric(Transect)], main=paste(layer.names[i], " StrS vs \n Microtopography", sep=""))
#   mtext(text="StrS", side=2, outer=T)
#   mtext(text="StrS - Ecotone", side=3, outer=T)
#}
#dev.off()
#
#
#
#png(filename="figures_exploratory/StrDivTundra_StrDvsEnvironment_Ecotone.png", height=30, width=20, unit="cm", res=300)
#
#par(mfrow=c(5,4), mar=c(4,3,3,2), oma=c(2,2,2,3))
#for(i in 1:Nlayers) {
#   data.temp <- subset(data.full.ecotone, Layer==layer.names[i])
#   plot(StrD ~ Altitude, data=data.temp, pch=21, bg=colors.trans[as.numeric(Transect)], main=paste(layer.names[i], " StrD vs \n Altitude", sep=""))
#   plot(StrD ~ DistFromLake, data=data.temp, pch=21, bg=colors.trans[as.numeric(Transect)], main=paste(layer.names[i], " StrD vs \n Distance", sep=""))
#   plot(StrD ~ pH, data=data.temp, pch=21, bg=colors.trans[as.numeric(Transect)], main=paste(layer.names[i], " StrD vs \n pH", sep=""))
#   plot(StrD ~ Microtopography, data=data.temp, pch=21, bg=colors.trans[as.numeric(Transect)], main=paste(layer.names[i], " StrD vs \n Microtopography", sep=""))
#   mtext(text="StrD", side=2, outer=T)
#   mtext(text="StrD - Ecotone", side=3, outer=T)
#}
#dev.off()



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


###### Models: adjusted via Generalized Estimation Equations with ar1 autocorrelation stucture, per substransect. After adjusting the full model (without interactions), compare between all possible models using dredge. Afterwards, select the simples model among the ones with dQIC < 2 and, in addition, calculate relative variable importance using the importance funtion, which sums the wQIC values of all the models containing the explanatory variable being considered. QIC is used because the models with autocorrelation do not permit calculating AIC.

### Tundra - StrS - ground
mod.StrS.tundra.ground.gee.full <- geeglm(StrS ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.tundra.ground, family=poisson(link="log"), id=Transect, corstr="ar1", na.action=na.fail)
dredge.StrS.tundra.ground.gee <- MuMIn::dredge(mod.StrS.tundra.ground.gee.full)
dredge.StrS.tundra.ground.gee
# The simplest model is the null model, with a dQIC of 0..79


### Tundra - StrD - ground
mod.StrD.tundra.ground.gee.full <- geeglm(StrD ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.tundra.ground, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)
dredge.StrD.tundra.ground.gee <- MuMIn::dredge(mod.StrD.tundra.ground.gee.full)
dredge.StrD.tundra.ground.gee
# Two simplest models: one with distance from lake, one with microtopography
mod.StrD.tundra.ground.gee.DistFromLake <- geeglm(StrD ~ DistFromLake, data=data.full.tundra.ground, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)

### Ecotone - StrS - ground
mod.StrS.ecotone.ground.gee.full <- geeglm(StrS ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.ecotone.ground, family=poisson, id=Transect, corstr="ar1", na.action=na.fail)
dredge.StrS.ecotone.ground.gee <- MuMIn::dredge(mod.StrS.ecotone.ground.gee.full)
dredge.StrS.ecotone.ground.gee 
# Simplest model: distance from lake and microtopography

### Ecotone - StrD - ground
mod.StrD.ecotone.ground.gee.full <- geeglm(StrD ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.ecotone.ground, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)
dredge.StrD.ecotone.ground.gee <- MuMIn::dredge(mod.StrD.ecotone.ground.gee.full)
dredge.StrD.ecotone.ground.gee 
# Simplest model: microtopography




### Tundra - StrS - herbaceous
mod.StrS.tundra.herbaceous.gee.full <- geeglm(StrS ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.tundra.herbaceous, family=poisson, id=Transect, corstr="ar1", na.action=na.fail)
dredge.StrS.tundra.herbaceous.gee <- MuMIn::dredge(mod.StrS.tundra.herbaceous.gee.full)
dredge.StrS.tundra.herbaceous.gee
# Simplest model: distance from lake


### Tundra - StrD - herbaceous
mod.StrD.tundra.herbaceous.gee.full <- geeglm(StrD ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.tundra.herbaceous, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)
dredge.StrD.tundra.herbaceous.gee <- MuMIn::dredge(mod.StrD.tundra.herbaceous.gee.full)
dredge.StrD.tundra.herbaceous.gee
# Simplest model: distance from lake, microtopography, pH



### Ecotone - StrS - herbaceous
mod.StrS.ecotone.herbaceous.gee.full <- geeglm(StrS ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.ecotone.herbaceous, family=poisson, id=Transect, corstr="ar1", na.action=na.fail)
dredge.StrS.ecotone.herbaceous.gee <- MuMIn::dredge(mod.StrS.ecotone.herbaceous.gee.full)
dredge.StrS.ecotone.herbaceous.gee 
# Simplest model: null model

### Ecotone - StrD - herbaceous
mod.StrD.ecotone.herbaceous.gee.full <- geeglm(StrD ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.ecotone.herbaceous, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)
dredge.StrD.ecotone.herbaceous.gee <- MuMIn::dredge(mod.StrD.ecotone.herbaceous.gee.full)
dredge.StrD.ecotone.herbaceous.gee 
# Simplest model: null


### Tundra - StrS - woody
mod.StrS.tundra.woody.gee.full <- geeglm(StrS ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.tundra.woody, family=poisson, id=Transect, corstr="ar1", na.action=na.fail)
dredge.StrS.tundra.woody.gee <- MuMIn::dredge(mod.StrS.tundra.woody.gee.full)
dredge.StrS.tundra.woody.gee
# The simplest model contains pH


### Tundra - StrD - woody
mod.StrD.tundra.woody.gee.full <- geeglm(StrD ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.tundra.woody, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)
dredge.StrD.tundra.woody.gee <- MuMIn::dredge(mod.StrD.tundra.woody.gee.full)
dredge.StrD.tundra.woody.gee
# The simplest model is the null model


### Ecotone - StrS - woody
mod.StrS.ecotone.woody.gee.full <- geeglm(StrS ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.ecotone.woody, family=poisson, id=Transect, corstr="ar1", na.action=na.fail)
dredge.StrS.ecotone.woody.gee <- MuMIn::dredge(mod.StrS.ecotone.woody.gee.full)
dredge.StrS.ecotone.woody.gee 
# Simplest model: microtopography and pH

### Ecotone - StrD - woody
mod.StrD.ecotone.woody.gee.full <- geeglm(StrD ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.ecotone.woody, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)
dredge.StrD.ecotone.woody.gee <- MuMIn::dredge(mod.StrD.ecotone.woody.gee.full)
dredge.StrD.ecotone.woody.gee 
# Simplest model: null


### Tundra - StrS - deadwood
mod.StrS.tundra.deadwood.gee.full <- geeglm(StrS ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.tundra.deadwood, family=poisson, id=Transect, corstr="ar1", na.action=na.fail)
dredge.StrS.tundra.deadwood.gee <- MuMIn::dredge(mod.StrS.tundra.deadwood.gee.full)
dredge.StrS.tundra.deadwood.gee
# The simplest model contains distance from lake and pH


### Tundra - StrD - deadwood
mod.StrD.tundra.deadwood.gee.full <- geeglm(StrD ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.tundra.deadwood, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)
dredge.StrD.tundra.deadwood.gee <- MuMIn::dredge(mod.StrD.tundra.deadwood.gee.full)
dredge.StrD.tundra.deadwood.gee
# Three simplest models: one with distance from lake, one with altitude and one with pH

### Ecotone - StrS - deadwood
mod.StrS.ecotone.deadwood.gee.full <- geeglm(StrS ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.ecotone.deadwood, family=poisson, id=Transect, corstr="ar1", na.action=na.fail)
dredge.StrS.ecotone.deadwood.gee <- MuMIn::dredge(mod.StrS.ecotone.deadwood.gee.full)
dredge.StrS.ecotone.deadwood.gee 
# Simplest model: microtopography and pH

### Ecotone - StrD - deadwood
mod.StrD.ecotone.deadwood.gee.full <- geeglm(StrD ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.ecotone.deadwood, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)
dredge.StrD.ecotone.deadwood.gee <- MuMIn::dredge(mod.StrD.ecotone.deadwood.gee.full)
dredge.StrD.ecotone.deadwood.gee 
# Simplest model: null

### Tundra - StrS - liveplants
mod.StrS.tundra.liveplants.gee.full <- geeglm(StrS ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.tundra.liveplants, family=poisson, id=Transect, corstr="ar1", na.action=na.fail)
dredge.StrS.tundra.liveplants.gee <- MuMIn::dredge(mod.StrS.tundra.liveplants.gee.full)
dredge.StrS.tundra.liveplants.gee

### Tundra - StrD - liveplants
mod.StrD.tundra.liveplants.gee.full <- geeglm(StrD ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.tundra.liveplants, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)
dredge.StrD.tundra.liveplants.gee <- MuMIn::dredge(mod.StrD.tundra.liveplants.gee.full)
dredge.StrD.tundra.liveplants.gee
# Simplest model: distance from lake

### Ecotone - StrS - liveplants
mod.StrS.ecotone.liveplants.gee.full <- geeglm(StrS ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.ecotone.liveplants, family=poisson, id=Transect, corstr="ar1", na.action=na.fail)
dredge.StrS.ecotone.liveplants.gee <- MuMIn::dredge(mod.StrS.ecotone.liveplants.gee.full)
dredge.StrS.ecotone.liveplants.gee 
# Simplest model: altitude and pH

### Ecotone - StrD - liveplants
mod.StrD.ecotone.liveplants.gee.full <- geeglm(StrD ~ Altitude + DistFromLake + pH + Microtopography, data=data.full.ecotone.liveplants, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)
dredge.StrD.ecotone.liveplants.gee <- MuMIn::dredge(mod.StrD.ecotone.liveplants.gee.full)
dredge.StrD.ecotone.liveplants.gee 
# Simplest model: pH

render(input="scripts/StrDivTundra_script_v5_dredge_v2.Rmd", output_file="../resubmitted/resultsDredge.pdf")
getwd()




### Now, to try and calculate pseudo-R2 values... For the pseudo-R2: use marginal (m) values and the trigamma method when available.
mod.StrS.tundra.ground.gee.final <- geeglm(StrS ~ 1, data=data.full.tundra.ground, family=poisson(link="log"), id=Transect, corstr="ar1", na.action=na.fail)
mod.StrS.tundra.herbaceous.gee.final <- geeglm(StrS ~ DistFromLake, data=data.full.tundra.herbaceous, family=poisson, id=Transect, corstr="ar1", na.action=na.fail)
mod.StrS.tundra.woody.gee.final <- geeglm(StrS ~ pH, data=data.full.tundra.woody, family=poisson, id=Transect, corstr="ar1", na.action=na.fail)
mod.StrS.tundra.deadwood.gee.final <- geeglm(StrS ~ DistFromLake + pH, data=data.full.tundra.deadwood, family=poisson, id=Transect, corstr="ar1", na.action=na.fail)
mod.StrS.tundra.liveplants.gee.final <- geeglm(StrS ~ DistFromLake, data=data.full.tundra.liveplants, family=poisson, id=Transect, corstr="ar1", na.action=na.fail)


rsq(mod.StrS.tundra.ground.gee.final)
rsq(mod.StrS.tundra.herbaceous.gee.final)
rsq(mod.StrS.tundra.woody.gee.final)
rsq(mod.StrS.tundra.deadwood.gee.final)
rsq(mod.StrS.tundra.liveplants.gee.final)

mod.StrD.tundra.ground.gee.final1 <- geeglm(StrD ~ DistFromLake, data=data.full.tundra.ground, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)

mod.StrD.tundra.ground.gee.null <- geeglm(StrD ~ 1, data=data.full.tundra.ground, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)

mod.StrD.tundra.ground.gee.final2 <- geeglm(StrD ~ Microtopography, data=data.full.tundra.ground, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)
mod.StrD.tundra.herbaceous.gee.final <- geeglm(StrD ~ DistFromLake + Microtopography + pH, data=data.full.tundra.herbaceous, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)
mod.StrD.tundra.woody.gee.final <- geeglm(StrD ~ 1, data=data.full.tundra.woody, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)
mod.StrD.tundra.deadwood.gee.final1 <- geeglm(StrD ~ DistFromLake, data=data.full.tundra.deadwood, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)
mod.StrD.tundra.deadwood.gee.final2 <- geeglm(StrD ~ Altitude, data=data.full.tundra.deadwood, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)
mod.StrD.tundra.deadwood.gee.final3 <- geeglm(StrD ~ pH, data=data.full.tundra.deadwood, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)
mod.StrD.tundra.liveplants.gee.final <- geeglm(StrD ~ DistFromLake, data=data.full.tundra.liveplants, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)

rsq(mod.StrD.tundra.ground.gee.final1)
rsq(mod.StrD.tundra.ground.gee.final2)
rsq(mod.StrD.tundra.herbaceous.gee.final)
rsq(mod.StrD.tundra.woody.gee.final1)
rsq(mod.StrD.tundra.deadwood.gee.final1)
rsq(mod.StrD.tundra.deadwood.gee.final2)
rsq(mod.StrD.tundra.deadwood.gee.final3)
rsq(mod.StrD.tundra.liveplants.gee.final)





mod.StrS.ecotone.ground.gee.final <- geeglm(StrS ~  DistFromLake + Microtopography, data=data.full.ecotone.ground, family=poisson, id=Transect, corstr="ar1", na.action=na.fail)
mod.StrS.ecotone.herbaceous.gee.final <- geeglm(StrS ~ 1, data=data.full.ecotone.herbaceous, family=poisson, id=Transect, corstr="ar1", na.action=na.fail)
mod.StrS.ecotone.woody.gee.final <- geeglm(StrS ~ pH + Microtopography, data=data.full.ecotone.woody, family=poisson, id=Transect, corstr="ar1", na.action=na.fail)
mod.StrS.ecotone.deadwood.gee.final <- geeglm(StrS ~ pH + Microtopography, data=data.full.ecotone.deadwood, family=poisson, id=Transect, corstr="ar1", na.action=na.fail)
mod.StrS.ecotone.liveplants.gee.final <- geeglm(StrS ~ Altitude + pH, data=data.full.ecotone.liveplants, family=poisson, id=Transect, corstr="ar1", na.action=na.fail)

rsq(mod.StrS.ecotone.ground.gee.final)
rsq(mod.StrS.ecotone.herbaceous.gee.final)
rsq(mod.StrS.ecotone.woody.gee.final)
rsq(mod.StrS.ecotone.deadwood.gee.final)
rsq(mod.StrS.ecotone.liveplants.gee.final)



mod.StrD.ecotone.ground.gee.final <- geeglm(StrD ~ Microtopography, data=data.full.ecotone.ground, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)
mod.StrD.ecotone.herbaceous.gee.final <- geeglm(StrD ~ 1, data=data.full.ecotone.herbaceous, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)
mod.StrD.ecotone.woody.gee.final <- geeglm(StrD ~ 1, data=data.full.ecotone.woody, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)
mod.StrD.ecotone.deadwood.gee.final <- geeglm(StrD ~ 1, data=data.full.ecotone.deadwood, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)
mod.StrD.ecotone.liveplants.gee.final <- geeglm(StrD ~ pH, data=data.full.ecotone.liveplants, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)

rsq(mod.StrD.ecotone.ground.gee.final)
rsq(mod.StrD.ecotone.herbaceous.gee.final)
rsq(mod.StrD.ecotone.woody.gee.final)
rsq(mod.StrD.ecotone.deadwood.gee.final)
rsq(mod.StrD.ecotone.liveplants.gee.final)



### Figure - relation with explanatory variables

# Tundra - StrS: 5 figures
# Herbaceous ~ Distance; Woody ~ pH; Deadwood ~ Distance + pH; Live plants ~ Distance
# Tundra - StrD: 9 figures
# Ground ~ Distance; Ground ~ Microtopography; Herbaceous ~ Distance + Microtopography + pH; Deadwood ~ Distance; Deadwood ~ Elevation; Deadwood ~ pH; Live plants ~ Distance



# Just make figure for StrS and make it figure .... so

library(scales)


newdata.tundra.DistFromLake <- seq(min(data.full.tundra$DistFromLake), max(data.full.tundra$DistFromLake), length.out=100)
newdata.tundra.pH <- seq(min(data.full.tundra$pH), max(data.full.tundra$pH), length.out=100)

predict.tundra.herbaceous.StrS.DistFromLake <- predict(mod.StrS.tundra.herbaceous.gee.final, newdata=list(DistFromLake=newdata.tundra.DistFromLake), type="response")
predict.tundra.woody.StrS.pH <- predict(mod.StrS.tundra.woody.gee.final, newdata=list(pH=newdata.tundra.pH), type="response")
predict.tundra.deadwood.StrS.DistFromLake <- predict(mod.StrS.tundra.deadwood.gee.final, newdata=list(DistFromLake = newdata.tundra.DistFromLake, pH = rep(mean(data.full.tundra$pH), 100)), type="response")
predict.tundra.deadwood.StrS.pH <- predict(mod.StrS.tundra.deadwood.gee.final, newdata=list(pH = newdata.tundra.pH, DistFromLake = rep(mean(data.full.tundra$DistFromLake), 100)), type="response")
predict.tundra.liveplants.StrS.DistFromLake <- predict(mod.StrS.tundra.liveplants.gee.final, newdata=list(DistFromLake = newdata.tundra.DistFromLake), type="response")


setwd("/home/pavel/Profissional/Pesquisa/MyPapers_working/Pavel_2014_StrDivTundra/resubmitted")

png(filename="StrDivTundra_figure4_expl_tundra.png", height=25, width=25, unit="cm", res=300)

par(mfrow = c(5,4), mar=c(1,1,2,1), oma=c(6,5,2,2))
plot(StrS ~ Altitude, data=data.full.tundra.ground, xlab="", ylab="", yaxt="s", las=1, xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="A)", outer=F, adj=0, cex=1.1)
plot(StrS ~ DistFromLake, data=data.full.tundra.ground, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="B)", outer=F, adj=0, cex=1.1)
plot(StrS ~ Microtopography, data=data.full.tundra.ground, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="C)", outer=F, adj=0, cex=1.1)
plot(StrS ~ pH, data=data.full.tundra.ground, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="D)", outer=F, adj=0, cex=1.1)

plot(StrS ~ Altitude, data=data.full.tundra.herbaceous, xlab="", ylab="", yaxt="s", las=1, xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="E)", outer=F, adj=0, cex=1.1)
plot(StrS ~ DistFromLake, data=data.full.tundra.herbaceous, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="F)", outer=F, adj=0, cex=1.1)
lines(predict.tundra.herbaceous.StrS.DistFromLake ~ newdata.tundra.DistFromLake, lwd=2, col="blue")
plot(StrS ~ Microtopography, data=data.full.tundra.herbaceous, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="G)", outer=F, adj=0, cex=1.1)
plot(StrS ~ pH, data=data.full.tundra.herbaceous, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="H)", outer=F, adj=0, cex=1.1)


plot(StrS ~ Altitude, data=data.full.tundra.woody, xlab="", ylab="", yaxt="s", las=1, xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="I)", outer=F, adj=0, cex=1.1)
plot(StrS ~ DistFromLake, data=data.full.tundra.woody, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="J)", outer=F, adj=0, cex=1.1)
plot(StrS ~ Microtopography, data=data.full.tundra.woody, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="K)", outer=F, adj=0, cex=1.1)
plot(StrS ~ pH, data=data.full.tundra.woody, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="L)", outer=F, adj=0, cex=1.1)
lines(predict.tundra.woody.StrS.pH ~ newdata.tundra.pH, lwd=2, col="blue")

plot(StrS ~ Altitude, data=data.full.tundra.deadwood, xlab="", ylab="", yaxt="n", las=1, xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="M)", outer=F, adj=0, cex=1.1)
axis(side=2, las=1, at=c(0,1,2,3,4))
plot(StrS ~ DistFromLake, data=data.full.tundra.deadwood, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="N)", outer=F, adj=0, cex=1.1)
lines(predict.tundra.deadwood.StrS.DistFromLake ~ newdata.tundra.DistFromLake, lwd=2, col="blue")
plot(StrS ~ Microtopography, data=data.full.tundra.deadwood, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="O)", outer=F, adj=0, cex=1.1)
plot(StrS ~ pH, data=data.full.tundra.deadwood, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="P)", outer=F, adj=0, cex=1.1)
lines(predict.tundra.deadwood.StrS.pH ~ newdata.tundra.pH, lwd=2, col="blue")

plot(StrS ~ Altitude, data=data.full.tundra.liveplants, xlab="Elevation (m)", ylab="", las=1, yaxt="s", xaxt="s", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="Q)", outer=F, adj=0, cex=1.1)
plot(StrS ~ DistFromLake, data=data.full.tundra.liveplants, xlab="Distance from/n lakes (m)", ylab="", yaxt="n", xaxt="s", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="R)", outer=F, adj=0, cex=1.1)
lines(predict.tundra.liveplants.StrS.DistFromLake ~ newdata.tundra.DistFromLake, lwd=2, col="blue")
plot(StrS ~ Microtopography, data=data.full.tundra.liveplants, xlab="Microtopography (cm²)", ylab="", yaxt="n", xaxt="s", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="S)", outer=F, adj=0, cex=1.1)
plot(StrS ~ pH, data=data.full.tundra.liveplants, xlab="pH", ylab="", yaxt="n", xaxt="s", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="T)", outer=F, adj=0, cex=1.1)

mtext(side=1, at=0.125 + 0.25*c(0,1,2,3), text=c("Elevation (m)", "Distance from\nlakes (m)", "Variance in\nmicrotopography (cm²)", "pH"), outer=T, line=3, padj=0)
mtext(side=2, at=c(0.1, 0.3, 0.5, 0.7, 0.9), text=c("Live plants", "Deadwood", "Woody", "Herbaceous", "Ground"), outer=T, line=2)
mtext(side=2, text="Number of structural elements", outer=T, line=3.5)

dev.off()





newdata.ecotone.DistFromLake <- seq(min(data.full.ecotone$DistFromLake), max(data.full.ecotone$DistFromLake), length.out=100)
newdata.ecotone.pH <- seq(min(data.full.ecotone$pH), max(data.full.ecotone$pH), length.out=100)
newdata.ecotone.Altitude <- seq(min(data.full.ecotone$Altitude), max(data.full.ecotone$Altitude), length.out=100)
newdata.ecotone.Microtopography <- seq(min(data.full.ecotone$Microtopography), max(data.full.ecotone$Microtopography), length.out=100)



predict.ecotone.ground.StrS.Microtopography <- predict(mod.StrS.ecotone.ground.gee.final, newdata=list(Microtopography = newdata.ecotone.Microtopography, DistFromLake = rep(mean(data.full.ecotone$DistFromLake), 100)), type="response")
predict.ecotone.ground.StrS.DistFromLake <- predict(mod.StrS.ecotone.ground.gee.final, newdata=list(DistFromLake = newdata.ecotone.DistFromLake, Microtopography = rep(mean(data.full.ecotone$Microtopography), 100)), type="response")

predict.ecotone.woody.StrS.Microtopography <- predict(mod.StrS.ecotone.woody.gee.final, newdata=list(Microtopography = newdata.ecotone.Microtopography, pH = rep(mean(data.full.ecotone$pH), 100)), type="response")
predict.ecotone.woody.StrS.pH <- predict(mod.StrS.ecotone.woody.gee.final, newdata=list(pH = newdata.ecotone.pH, Microtopography = rep(mean(data.full.ecotone$Microtopography), 100)), type="response")

predict.ecotone.deadwood.StrS.Microtopography <- predict(mod.StrS.ecotone.deadwood.gee.final, newdata=list(Microtopography = newdata.ecotone.Microtopography, pH = rep(mean(data.full.ecotone$pH), 100)), type="response")
predict.ecotone.deadwood.StrS.pH <- predict(mod.StrS.ecotone.deadwood.gee.final, newdata=list(pH = newdata.ecotone.pH, Microtopography = rep(mean(data.full.ecotone$Microtopography), 100)), type="response")

predict.ecotone.liveplants.StrS.Altitude <- predict(mod.StrS.ecotone.liveplants.gee.final, newdata=list(Altitude = newdata.ecotone.Altitude, pH = rep(mean(data.full.ecotone$pH), 100)), type="response")
predict.ecotone.liveplants.StrS.pH <- predict(mod.StrS.ecotone.liveplants.gee.final, newdata=list(pH = newdata.ecotone.pH, Altitude = rep(mean(data.full.ecotone$Altitude), 100)), type="response")


setwd("/home/pavel/Profissional/Pesquisa/MyPapers_working/Pavel_2014_StrDivTundra/resubmitted")

png(filename="StrDivTundra_figure5_expl_ecotone.png", height=25, width=25, unit="cm", res=300)

par(mfrow = c(5,4), mar=c(1,1,2,1), oma=c(6,5,2,2))
plot(StrS ~ Altitude, data=data.full.ecotone.ground, xlab="", ylab="", yaxt="s", las=1, xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="A)", outer=F, adj=0, cex=1.1)
plot(StrS ~ DistFromLake, data=data.full.ecotone.ground, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="B)", outer=F, adj=0, cex=1.1)
lines(predict.ecotone.ground.StrS.DistFromLake ~ newdata.ecotone.DistFromLake, lwd=2, col="blue")
plot(StrS ~ Microtopography, data=data.full.ecotone.ground, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="C)", outer=F, adj=0, cex=1.1)
lines(predict.ecotone.ground.StrS.Microtopography ~ newdata.ecotone.Microtopography, lwd=2, col="blue")
plot(StrS ~ pH, data=data.full.ecotone.ground, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="D)", outer=F, adj=0, cex=1.1)

plot(StrS ~ Altitude, data=data.full.ecotone.herbaceous, xlab="", ylab="", yaxt="s", las=1, xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="E)", outer=F, adj=0, cex=1.1)
plot(StrS ~ DistFromLake, data=data.full.ecotone.herbaceous, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="F)", outer=F, adj=0, cex=1.1)
plot(StrS ~ Microtopography, data=data.full.ecotone.herbaceous, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="G)", outer=F, adj=0, cex=1.1)
plot(StrS ~ pH, data=data.full.ecotone.herbaceous, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="H)", outer=F, adj=0, cex=1.1)


plot(StrS ~ Altitude, data=data.full.ecotone.woody, xlab="", ylab="", yaxt="s", las=1, xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="I)", outer=F, adj=0, cex=1.1)
plot(StrS ~ DistFromLake, data=data.full.ecotone.woody, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="J)", outer=F, adj=0, cex=1.1)
plot(StrS ~ Microtopography, data=data.full.ecotone.woody, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="K)", outer=F, adj=0, cex=1.1)
lines(predict.ecotone.woody.StrS.Microtopography ~ newdata.ecotone.Microtopography, lwd=2, col="blue")
plot(StrS ~ pH, data=data.full.ecotone.woody, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="L)", outer=F, adj=0, cex=1.1)
lines(predict.ecotone.woody.StrS.pH ~ newdata.ecotone.pH, lwd=2, col="blue")

plot(StrS ~ Altitude, data=data.full.ecotone.deadwood, xlab="", ylab="", yaxt="n", las=1, xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="M)", outer=F, adj=0, cex=1.1)
axis(side=2, las=1, at=c(0,1,2,3,4))
plot(StrS ~ DistFromLake, data=data.full.ecotone.deadwood, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="N)", outer=F, adj=0, cex=1.1)
plot(StrS ~ Microtopography, data=data.full.ecotone.deadwood, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="O)", outer=F, adj=0, cex=1.1)
lines(predict.ecotone.deadwood.StrS.Microtopography ~ newdata.ecotone.Microtopography, lwd=2, col="blue")
plot(StrS ~ pH, data=data.full.ecotone.deadwood, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="P)", outer=F, adj=0, cex=1.1)
lines(predict.ecotone.deadwood.StrS.pH ~ newdata.ecotone.pH, lwd=2, col="blue")

plot(StrS ~ Altitude, data=data.full.ecotone.liveplants, xlab="Elevation (m)", ylab="", las=1, yaxt="s", xaxt="s", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="Q)", outer=F, adj=0, cex=1.1)
lines(predict.ecotone.liveplants.StrS.Altitude ~ newdata.ecotone.Altitude, lwd=2, col="blue")
plot(StrS ~ DistFromLake, data=data.full.ecotone.liveplants, xlab="Distance from/n lakes (m)", ylab="", yaxt="n", xaxt="s", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="R)", outer=F, adj=0, cex=1.1)
plot(StrS ~ Microtopography, data=data.full.ecotone.liveplants, xlab="Microtopography (cm²)", ylab="", yaxt="n", xaxt="s", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="S)", outer=F, adj=0, cex=1.1)
plot(StrS ~ pH, data=data.full.ecotone.liveplants, xlab="pH", ylab="", yaxt="n", xaxt="s", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="T)", outer=F, adj=0, cex=1.1)
lines(predict.ecotone.liveplants.StrS.pH ~ newdata.ecotone.pH, lwd=2, col="blue")

mtext(side=1, at=0.125 + 0.25*c(0,1,2,3), text=c("Elevation (m)", "Distance from\nlakes (m)", "Variance in\nmicrotopography (cm²)", "pH"), outer=T, line=3, padj=0)
mtext(side=2, at=c(0.1, 0.3, 0.5, 0.7, 0.9), text=c("Live plants", "Deadwood", "Woody", "Herbaceous", "Ground"), outer=T, line=2)
mtext(side=2, text="Number of structural elements", outer=T, line=3.5)

dev.off()





### Figures for supplementary material (creating them here before placing them in the markdown file):

### StrD - tundra


mod.StrD.tundra.ground.gee.final1 <- geeglm(StrD ~ DistFromLake, data=data.full.tundra.ground, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)
mod.StrD.tundra.ground.gee.final2 <- geeglm(StrD ~ Microtopography, data=data.full.tundra.ground, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)
mod.StrD.tundra.herbaceous.gee.final <- geeglm(StrD ~ DistFromLake + Microtopography + pH, data=data.full.tundra.herbaceous, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)
mod.StrD.tundra.woody.gee.final <- geeglm(StrD ~ 1, data=data.full.tundra.woody, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)
mod.StrD.tundra.deadwood.gee.final1 <- geeglm(StrD ~ DistFromLake, data=data.full.tundra.deadwood, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)
mod.StrD.tundra.deadwood.gee.final2 <- geeglm(StrD ~ Altitude, data=data.full.tundra.deadwood, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)
mod.StrD.tundra.deadwood.gee.final3 <- geeglm(StrD ~ pH, data=data.full.tundra.deadwood, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)
mod.StrD.tundra.liveplants.gee.final <- geeglm(StrD ~ DistFromLake, data=data.full.tundra.liveplants, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)

newdata.tundra.DistFromLake <- seq(min(data.full.tundra$DistFromLake), max(data.full.tundra$DistFromLake), length.out=100)
newdata.tundra.pH <- seq(min(data.full.tundra$pH), max(data.full.tundra$pH), length.out=100)
newdata.tundra.Altitude <- seq(min(data.full.tundra$Altitude), max(data.full.tundra$Altitude), length.out=100)
newdata.tundra.Microtopography <- seq(min(data.full.tundra$Microtopography), max(data.full.tundra$Microtopography), length.out=100)


predict.tundra.ground.StrD.DistFromLake <- predict(mod.StrD.tundra.ground.gee.final1, newdata=list(DistFromLake=newdata.tundra.DistFromLake), type="response")
predict.tundra.ground.StrD.Microtopography <- predict(mod.StrD.tundra.ground.gee.final2, newdata=list(Microtopography=newdata.tundra.Microtopography), type="response")
predict.tundra.herbaceous.StrD.DistFromLake <- predict(mod.StrD.tundra.herbaceous.gee.final, newdata=list(DistFromLake=newdata.tundra.DistFromLake, Microtopography = rep(mean(data.full.tundra$Microtopography), 100), pH = rep(mean(data.full.tundra$pH), 100)), type="response")
predict.tundra.herbaceous.StrD.Microtopography <- predict(mod.StrD.tundra.herbaceous.gee.final, newdata=list(Microtopography=newdata.tundra.Microtopography, DistFromLake = rep(mean(data.full.tundra$DistFromLake), 100), pH = rep(mean(data.full.tundra$pH), 100)), type="response")
predict.tundra.herbaceous.StrD.pH <- predict(mod.StrD.tundra.herbaceous.gee.final, newdata=list(pH=newdata.tundra.pH, Microtopography = rep(mean(data.full.tundra$Microtopography), 100), DistFromLake = rep(mean(data.full.tundra$DistFromLake), 100)), type="response")
predict.tundra.deadwood.StrD.DistFromLake <- predict(mod.StrD.tundra.deadwood.gee.final1, newdata=list(DistFromLake=newdata.tundra.DistFromLake), type="response")
predict.tundra.deadwood.StrD.Altitude <- predict(mod.StrD.tundra.deadwood.gee.final2, newdata=list(Altitude=newdata.tundra.Altitude), type="response")
predict.tundra.deadwood.StrD.pH <- predict(mod.StrD.tundra.deadwood.gee.final3, newdata=list(pH=newdata.tundra.pH), type="response")
predict.tundra.liveplants.StrD.DistFromLake <- predict(mod.StrD.tundra.liveplants.gee.final, newdata=list(DistFromLake=newdata.tundra.DistFromLake), type="response")

par(mfrow = c(5,4), mar=c(1,1,2,1), oma=c(6,5,2,2))
plot(StrD ~ Altitude, data=data.full.tundra.ground, xlab="", ylab="", yaxt="s", las=1, xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="A)", outer=F, adj=0, cex=1.1)
plot(StrD ~ DistFromLake, data=data.full.tundra.ground, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="B)", outer=F, adj=0, cex=1.1)
lines(predict.tundra.ground.StrD.DistFromLake ~ newdata.tundra.DistFromLake, lwd=2, col="blue")
plot(StrD ~ Microtopography, data=data.full.tundra.ground, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="C)", outer=F, adj=0, cex=1.1)
lines(predict.tundra.ground.StrD.Microtopography ~ newdata.tundra.Microtopography, lwd=2, col="blue")
plot(StrD ~ pH, data=data.full.tundra.ground, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="D)", outer=F, adj=0, cex=1.1)

plot(StrD ~ Altitude, data=data.full.tundra.herbaceous, xlab="", ylab="", yaxt="s", las=1, xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="E)", outer=F, adj=0, cex=1.1)
plot(StrD ~ DistFromLake, data=data.full.tundra.herbaceous, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="F)", outer=F, adj=0, cex=1.1)
lines(predict.tundra.herbaceous.StrD.DistFromLake ~ newdata.tundra.DistFromLake, lwd=2, col="blue")
plot(StrD ~ Microtopography, data=data.full.tundra.herbaceous, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="G)", outer=F, adj=0, cex=1.1)
lines(predict.tundra.herbaceous.StrD.Microtopography ~ newdata.tundra.Microtopography, lwd=2, col="blue")
plot(StrD ~ pH, data=data.full.tundra.herbaceous, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="H)", outer=F, adj=0, cex=1.1)
lines(predict.tundra.herbaceous.StrD.pH ~ newdata.tundra.pH, lwd=2, col="blue")


plot(StrD ~ Altitude, data=data.full.tundra.woody, xlab="", ylab="", yaxt="s", las=1, xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="I)", outer=F, adj=0, cex=1.1)
plot(StrD ~ DistFromLake, data=data.full.tundra.woody, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="J)", outer=F, adj=0, cex=1.1)
plot(StrD ~ Microtopography, data=data.full.tundra.woody, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="K)", outer=F, adj=0, cex=1.1)
plot(StrD ~ pH, data=data.full.tundra.woody, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="L)", outer=F, adj=0, cex=1.1)

plot(StrD ~ Altitude, data=data.full.tundra.deadwood, xlab="", ylab="", yaxt="n", las=1, xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="M)", outer=F, adj=0, cex=1.1)
axis(side=2, las=1, at=c(0,1,2,3,4))
lines(predict.tundra.deadwood.StrD.Altitude ~ newdata.tundra.Altitude, lwd=2, col="blue")
plot(StrD ~ DistFromLake, data=data.full.tundra.deadwood, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="N)", outer=F, adj=0, cex=1.1)
lines(predict.tundra.deadwood.StrD.DistFromLake ~ newdata.tundra.DistFromLake, lwd=2, col="blue")
plot(StrD ~ Microtopography, data=data.full.tundra.deadwood, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="O)", outer=F, adj=0, cex=1.1)
plot(StrD ~ pH, data=data.full.tundra.deadwood, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="P)", outer=F, adj=0, cex=1.1)
lines(predict.tundra.deadwood.StrD.pH ~ newdata.tundra.pH, lwd=2, col="blue")

plot(StrD ~ Altitude, data=data.full.tundra.liveplants, xlab="Elevation (m)", ylab="", las=1, yaxt="s", xaxt="s", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="Q)", outer=F, adj=0, cex=1.1)
plot(StrD ~ DistFromLake, data=data.full.tundra.liveplants, xlab="Distance from/n lakes (m)", ylab="", yaxt="n", xaxt="s", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="R)", outer=F, adj=0, cex=1.1)
lines(predict.tundra.liveplants.StrD.DistFromLake ~ newdata.tundra.DistFromLake, lwd=2, col="blue")
plot(StrD ~ Microtopography, data=data.full.tundra.liveplants, xlab="Microtopography (cm²)", ylab="", yaxt="n", xaxt="s", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="S)", outer=F, adj=0, cex=1.1)
plot(StrD ~ pH, data=data.full.tundra.liveplants, xlab="pH", ylab="", yaxt="n", xaxt="s", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="T)", outer=F, adj=0, cex=1.1)

mtext(side=1, at=0.125 + 0.25*c(0,1,2,3), text=c("Elevation (m)", "Distance from\nlakes (m)", "Variance in\nmicrotopography (cm²)", "pH"), outer=T, line=3, padj=0)
mtext(side=2, at=c(0.1, 0.3, 0.5, 0.7, 0.9), text=c("Live plants", "Deadwood", "Woody", "Herbaceous", "Ground"), outer=T, line=2)
mtext(side=2, text="Number of structural elements", outer=T, line=3.5)



### StrD - Ecotone
mod.StrD.ecotone.ground.gee.final <- geeglm(StrD ~ Microtopography, data=data.full.ecotone.ground, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)
mod.StrD.ecotone.liveplants.gee.final <- geeglm(StrD ~ pH, data=data.full.ecotone.liveplants, family=gaussian, id=Transect, corstr="ar1", na.action=na.fail)


newdata.ecotone.DistFromLake <- seq(min(data.full.ecotone$DistFromLake), max(data.full.ecotone$DistFromLake), length.out=100)
newdata.ecotone.pH <- seq(min(data.full.ecotone$pH), max(data.full.ecotone$pH), length.out=100)
newdata.ecotone.Altitude <- seq(min(data.full.ecotone$Altitude), max(data.full.ecotone$Altitude), length.out=100)
newdata.ecotone.Microtopography <- seq(min(data.full.ecotone$Microtopography), max(data.full.ecotone$Microtopography), length.out=100)


predict.ecotone.ground.StrD.Microtopography <- predict(mod.StrD.ecotone.ground.gee.final, newdata=list(Microtopography=newdata.ecotone.Microtopography), type="response")
predict.ecotone.liveplants.StrD.pH <- predict(mod.StrD.ecotone.liveplants.gee.final, newdata=list(pH=newdata.ecotone.pH), type="response")




par(mfrow = c(5,4), mar=c(1,1,2,1), oma=c(6,5,2,2))
plot(StrD ~ Altitude, data=data.full.ecotone.ground, xlab="", ylab="", yaxt="s", las=1, xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="A)", outer=F, adj=0, cex=1.1)
plot(StrD ~ DistFromLake, data=data.full.ecotone.ground, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="B)", outer=F, adj=0, cex=1.1)
plot(StrD ~ Microtopography, data=data.full.ecotone.ground, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="C)", outer=F, adj=0, cex=1.1)
lines(predict.ecotone.ground.StrD.Microtopography ~ newdata.ecotone.Microtopography, lwd=2, col="blue")
plot(StrD ~ pH, data=data.full.ecotone.ground, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="D)", outer=F, adj=0, cex=1.1)

plot(StrD ~ Altitude, data=data.full.ecotone.herbaceous, xlab="", ylab="", yaxt="s", las=1, xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="E)", outer=F, adj=0, cex=1.1)
plot(StrD ~ DistFromLake, data=data.full.ecotone.herbaceous, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="F)", outer=F, adj=0, cex=1.1)
plot(StrD ~ Microtopography, data=data.full.ecotone.herbaceous, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="G)", outer=F, adj=0, cex=1.1)
plot(StrD ~ pH, data=data.full.ecotone.herbaceous, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="H)", outer=F, adj=0, cex=1.1)


plot(StrD ~ Altitude, data=data.full.ecotone.woody, xlab="", ylab="", yaxt="s", las=1, xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="I)", outer=F, adj=0, cex=1.1)
plot(StrD ~ DistFromLake, data=data.full.ecotone.woody, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="J)", outer=F, adj=0, cex=1.1)
plot(StrD ~ Microtopography, data=data.full.ecotone.woody, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="K)", outer=F, adj=0, cex=1.1)
plot(StrD ~ pH, data=data.full.ecotone.woody, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="L)", outer=F, adj=0, cex=1.1)

plot(StrD ~ Altitude, data=data.full.ecotone.deadwood, xlab="", ylab="", yaxt="n", las=1, xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="M)", outer=F, adj=0, cex=1.1)
axis(side=2, las=1, at=c(0,1,2,3,4))
plot(StrD ~ DistFromLake, data=data.full.ecotone.deadwood, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="N)", outer=F, adj=0, cex=1.1)
plot(StrD ~ Microtopography, data=data.full.ecotone.deadwood, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="O)", outer=F, adj=0, cex=1.1)
plot(StrD ~ pH, data=data.full.ecotone.deadwood, xlab="", ylab="", yaxt="n", xaxt="n", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="P)", outer=F, adj=0, cex=1.1)

plot(StrD ~ Altitude, data=data.full.ecotone.liveplants, xlab="Elevation (m)", ylab="", las=1, yaxt="s", xaxt="s", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="Q)", outer=F, adj=0, cex=1.1)
plot(StrD ~ DistFromLake, data=data.full.ecotone.liveplants, xlab="Distance from/n lakes (m)", ylab="", yaxt="n", xaxt="s", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="R)", outer=F, adj=0, cex=1.1)
plot(StrD ~ Microtopography, data=data.full.ecotone.liveplants, xlab="Microtopography (cm²)", ylab="", yaxt="n", xaxt="s", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="S)", outer=F, adj=0, cex=1.1)
plot(StrD ~ pH, data=data.full.ecotone.liveplants, xlab="pH", ylab="", yaxt="n", xaxt="s", pch=21, bg=alpha("black", 0.1)); mtext(side=3, text="T)", outer=F, adj=0, cex=1.1)
lines(predict.ecotone.liveplants.StrD.pH ~ newdata.ecotone.pH, lwd=2, col="blue")


mtext(side=1, at=0.125 + 0.25*c(0,1,2,3), text=c("Elevation (m)", "Distance from\nlakes (m)", "Variance in\nmicrotopography (cm²)", "pH"), outer=T, line=3, padj=0)
mtext(side=2, at=c(0.1, 0.3, 0.5, 0.7, 0.9), text=c("Live plants", "Deadwood", "Woody", "Herbaceous", "Ground"), outer=T, line=2)
mtext(side=2, text="Number of structural elements", outer=T, line=3.5)









### Figures not used in the end


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










