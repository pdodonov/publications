##################################################################
### This is the code use for data analysis of the manuscript   ###
### Spatial pattern of invasive and native graminoids          ###   
###          in the Brazilian cerrado                          ###
### by Pavel Dodonov, Rafael Xavier, Karen Harper              ###
###             and Dalva Silva-Matos.                         ###
### Contact: pdodonov@gmail.com                                ###
### No rights reserved except those reserved by the publisher. ###
### Updated versions of the functions may be available at      ###
###            github.com/pdodonov                             ###
##################################################################

### Below, parts I, II and III refer to pre-analysis data organization and some plots;
### Parts 1-5 refer specifically to the manuscript's questions.

### Note: I chose to write a script for each transect and sometimes each response variables instead of automating this procedure using loops. I made this decision because having multiple loops within loops increases the possibilities of making a mistake somewhere. On the other hand, the text editor I use - Sublime Text - makes it quite easy to quickly replace parts of object names in a block of code.

### Another note: I designed this script to be run by parts, i.e. R can be closed and the computer my be turned off after each analysis without necessarily saving the workspace. Workspace is saved in two places only for safety.

### The script has three parts: part I for the libraries, part II for the objects, part III for the exploratory figures, parts 1-4 for data analysis, and part 5 for organizing the results of some analyses. Parts I and II must always be run before running parts 1-4. Part 2 does not depend on part 1, part 3 does not depend on parts 1 or 2 and so on.

### Part I - set up the packages and user-defined functions
library(wmtsa) # Used for wavelet analysis
library(mgcv) # Used for the additive models
library(bbmle) # Used to calculate AIC values.

count <- function(x) {
  res <- sum(x>0)
  return(res)
}

freq <- function(x) {
  res <- sum(x>0) / length(x)
  return(res)
}

test.firebreak <- function(x) any(x==c("firebreak","Firebreak","firebreak_regenerating","railroad"))


calc.freq <- function(x) return(sum(x>0)/length(x))

mean.present <- function(x) {
  foo <- mean(x[x>0])
  if(is.na(foo)) foo <- 0
  return(foo)
}

rep.patch <- function(patches) {
  patch.lengths <- patches$End-patches$Start+1
  data.patch <- rep(as.character(patches[1,1]), patch.lengths[1] )
  for (i in 2:nrow(patches)) {
    data.patch <- c(data.patch, rep(as.character(patches[i,1]), patch.lengths[i] ) )
    }
  data.patch <- as.factor(data.patch)
  data.patch.char <- as.character(data.patch)
  return(data.patch.char)
}

replace.cover <- function(x) {
  foo <- x
  foo[x==1] <- 0.0625
  foo[x==2] <- 0.1875
  foo[x==3] <- 0.375
  foo[x==4] <- 0.625
  foo[x>=5] <- 0.875
  return(foo)
  }

gam.modSel <- function(resp, expl) { # x data.frame with one response variable; y has the explanatory variables.
  foo <- cbind(resp,expl)
  DistElev <- gam(resp ~ s(DistToFirebreak, fx=F, k=5) + s(Elevation, fx=F, k=5) + Vegetation, family=binomial, data=foo)
  Dist <- gam(resp ~ s(DistToFirebreak, fx=F, k=5) + Vegetation, family=binomial, data=foo)
  Elev <- gam(resp ~ s(Elevation, fx=F, k=5) + Vegetation, family=binomial, data=foo)
  null <- gam(resp ~ Vegetation, family=binomial, data=foo)
  AIC.resp <- AICtab(DistElev, Dist, Elev, null, weights=T, sort=F)
  wAIC.Dist <- sum(AIC.resp$weight[c(1,2)])
  wAIC.Elev <- sum(AIC.resp$weight[c(1,3)])
  return(list(AICctab = AIC.resp, wAIC.Dist = wAIC.Dist, wAIC.Elev = wAIC.Elev))
}

signif.perm.above <- function(x) {
  signif <- sum(x >= x[1]) / length(x)
  return(signif)
}

signif.perm.below <- function(x) {
  signif <- sum(x <= x[1]) / length(x)
  return(signif)
}

MC1h <- function(x, Nperm=5000, keep.obs=T, print.loop=T, which.loop=100) {
  Nquad <- length(x)
  data.sim <- matrix(ncol=Nperm, nrow=Nquad)
  if(keep.obs) data.sim[,1] <- x
  #Create transition matrix used in the MC1 simulations
  foo <- c(x[2:Nquad], NA)
  bar <- c(NA, x[1:(Nquad-1)])
  x.trans <- data.frame(x, foo, bar)
  foo <- c(x.trans$x, x.trans$x)
  bar <- c(x.trans$foo, x.trans$bar)  
  transmatrix <- table(foo,bar)
  sums <- apply(transmatrix, 1, sum)
  transmatrix <- transmatrix / sums
  for (i in ifelse(keep.obs, 2, 1):Nperm) {
    foo <- numeric(Nquad)
    foo.section <- x.trans$x
    foo.trans <- transmatrix
    foo.classes <- colnames(foo.trans)
    first <- sample(length(foo),1)
    foo[first] <- sample(foo.section,1)
    #Going backward
    if (first != 1 ) {      
      for(k in first:2) {
        ref <- as.character(foo[k])
        trans.row <- foo.trans[ref,]
        trans.row <- cumsum(trans.row)
        random <- runif(1)
        test <- random > trans.row
        foo[k-1] <- foo.classes[sum(test)+1]
      }
    }
    if (first != Nquad) {
      for(k in first:(length(foo)-1)) {
        ref <- as.character(foo[k])
        trans.row <- foo.trans[ref,]
        trans.row <- cumsum(trans.row)
        random <- runif(1)
        test <- random > trans.row
        foo[k+1] <- foo.classes[sum(test)+1]
      }
    }
    foobar <- as.numeric(foo)
    if (print.loop) {
      if(i%%which.loop == 0) print(i)
    }
    data.sim[,i] <- foobar
  }
  return(data.sim)
}

MC1s <- function(x, sections, Nperm=5000, keep.obs=T, print.loop=T, which.loop=100) {
  Nquad <- length(x)
  data.sim <- matrix(ncol=Nperm, nrow=Nquad)
  if(keep.obs) data.sim[,1] <- x
  section.lengths <- sections$End-sections$Start+1
  data.section <- rep(as.character(sections[1,1]), section.lengths[1])
  for (i in 2:nrow(sections)) {
    data.section <- c(data.section, rep(as.character(sections[i,1]), section.lengths[i]))
    }
  data.section <- as.factor(data.section)
  data.section.char <- as.character(data.section)
  foo <- c(x[2:Nquad], NA)
  bar <- c(NA, x[1:(Nquad-1)])
  x.trans <- data.frame(x, foo, bar)
  Nsections <- nrow(sections)  
  sectiontypes <- unique(data.section.char)
  Nsectiontypes <- length(sectiontypes)
  x.trans$loc <- character(Nquad)
  data.section.char.next <- c(data.section.char[2:Nquad],NA)
  data.section.char.prev <- c(NA, data.section.char[1:Nquad-1])
  #State whether is quadrat is located in the middle of a transect section or between sections; in this latter case it may not be included in the transition matrix.
  x.trans$loc <- ifelse(data.section.char == data.section.char.next & data.section.char == data.section.char.prev, "middle", ifelse(data.section.char != data.section.char.prev, "first", ifelse(data.section.char != data.section.char.next, "last", "error")))
  x.trans$loc[1] <- "first"
  x.trans$loc[Nquad] <- "last"
  #Make a transition matrix per section type
  trans.matrices <- list()
  for(i in 1:Nsectiontypes) {
    sectiontype <- sectiontypes[i]
    foobar.rows <- (1:Nquad)[data.section.char == sectiontype]
    foobar <- x.trans[foobar.rows,]
    foobar [foobar$loc == "last", 2] <- NA
    foobar[foobar$loc == "first", 3] <- NA
    foo <- c(foobar$x, foobar$x)
    bar <- c(foobar$foo, foobar$bar)  
    trans.matrices[[i]] <- table(foo,bar)
    sums <- apply(trans.matrices[[i]], 1, sum)
    trans.matrices[[i]] <- trans.matrices[[i]] / sums
  }    
  names(trans.matrices) = sectiontypes
  #Simulate the data, separately for each section
  for (i in ifelse(keep.obs, 2, 1):Nperm) {
    foobar = numeric()
    for (j in 1:Nsections) {
      foo.section <- as.character(sections[j,1])
      foo.trans <- trans.matrices[[foo.section]]
      foo.classes <- colnames(foo.trans)
      foo <- numeric((sections[j,3]-sections[j,2]+1) )
      first <- sample(length(foo),1)
      foo[first] <- sample(x[sections[j,2]:sections[j,3]],1)
      #Going backward
        if (first != 1 ) {      
          for(k in first:2) {
            ref <- as.character(foo[k])
            trans.row <- foo.trans[ref,]
            trans.row <- cumsum(trans.row)
            random <- runif(1)
            test <- random > trans.row
            foo[k-1] <- foo.classes[sum(test)+1]
          }
        }
        if (first != length(foo)) {
          for(k in first:(length(foo)-1)) {
            ref <- as.character(foo[k])
            trans.row <- foo.trans[ref,]
            trans.row <- cumsum(trans.row)
            random <- runif(1)
            test <- random > trans.row
            foo[k+1] <- foo.classes[sum(test)+1]
          }
        }
        if (j == 1) foobar = foo else foobar = c(foobar,foo)
      }
    data.sim[,i]=as.numeric(foobar)
    if (print.loop) {
      if(i%%which.loop == 0) print(i)
    }
  }
  return(data.sim)
}

wavCWTvarCIs <- function(x, significance=T, sections=NULL, make.plot=F, zero.padding=T, remove.COI=T, test=F, scale.max=75, scale.var=T, pos.var=T, CI=T, CIquant=c(0.025, 0.975), effect.size=T,  keep.all=F, print.loop=T, which.loop=1, save.cwt=F, file.title="my_beautiful_wavelet", wav.template="gaussian2", scale.per.section=T) {
  #This function does not save the original wavelet output, retaining only scale and position variance. Such functionality may be added in the future. For now, the keep.all,save.cwt and file.title arguments are not fully functional.
  #This may not work for wavelets other than "gaussian2" (Mexican hat) and "haar" (Haar).
  data.main <- x
  Nperm <- ifelse(test,5,ncol(data.main))
  dimnames(data.main) <- NULL
  data.foo <- data.main[,1]
  Nquad <- length(data.foo)
  if (scale.per.section) {
    section.lengths <- sections$End-sections$Start+1
    data.section <- rep(as.character(sections[1,1]), section.lengths[1])
    for (k in 2:nrow(sections)) {
      data.section <- c(data.section, rep(as.character(sections[k,1]), section.lengths[k]))
    }
    data.section <- as.factor(data.section)
    section.types <- unique(sections[,1])
    Nsection.types <- length(section.types)
  }
  full.length <- Nquad
  half.length <- round(Nquad/2)
  if(zero.padding) data.foo <- c(rep(0,half.length),data.foo,rep(0,half.length)) #adds zeroes to the extremities of the transect
  wav.orig <- wavCWT(data.foo,scale.range=c(1,scale.max),wavelet=wav.template)
  scales <- attr(wav.orig,"scale")
  Nscales <- length(scales)
  #Zero-padding adds zeroes to the extremities of the transect, to a length equal to half of the transect
  if(zero.padding) {  
    #Save the attributes of the original wavelet transform
    wav.scale <- attr(wav.orig,"scale")
    wav.time <- 1:full.length
    wav.wavelet <- attr(wav.orig,"wavelet")
    wav.series <- attr(wav.orig,"series")
    wav.sampling.interval <- attr(wav.orig,"sampling.interval")
    wav.series.name <- attr(wav.orig,"series.name")
    wav.n.sample <- length(wav.time)
    wav.n.scale <- attr(wav.orig,"n.scale")
    wav.filter.arg <- attr(wav.orig,"filter.arg")
    #Remove parts of the wavelet transform corresponding to the artifficially added zeroes
    wav.orig <- wav.orig[(half.length+1):(half.length+full.length),]
    #Return the attributes to the corrected wavelet transform
    class(wav.orig) <- "wavCWT"
    attr(wav.orig,"scale") <- wav.scale
    attr(wav.orig,"time") <- wav.time
    attr(wav.orig,"wavelet") <- wav.wavelet
    attr(wav.orig,"series") <- wav.series
    attr(wav.orig,"sampling.interval") <- wav.sampling.interval
    attr(wav.orig,"series.names") <- wav.series.name
    attr(wav.orig,"n.sample") <- wav.n.sample
    attr(wav.orig,"n.scale") <- wav.n.scale
    attr(wav.orig,"filter.arg") <- wav.filter.arg
  }
  #The Cone Of Influence - COI - is the part of the wavelet transform that is affected by the transect's extremities, i.e. by quadrats for which no data has been collected. Inference made with these numbers is unreliable. The area affected is greater at larger scales, and therefore the COI has the shape of a, well, cone. Here I used an approach similar to that used by the PASSaGE software, which removes two quadrats for each scale for the MH wavelet and 1 quadrat for the Haar wavelet, from each side of the transect.  
  if(remove.COI) {
    for(k in 1:length(scales)) {
      index1 <- 1:(scales[k]*ifelse(wav.template=="gaussian2",2,1))
      index2 <- nrow(wav.orig)-index1+1
      indices <- c(index1,index2)
      wav.orig[indices,k] <- NA
    }
  }
  if(scale.var) wav.orig.scale <- apply(as.matrix(wav.orig)^2,2,mean, na.rm=T)
  if(pos.var) wav.orig.pos <- apply(as.matrix(wav.orig)^2,1,mean, na.rm=T)
  if(scale.per.section) {
    scale.var.section <- list()
    for (k in 1:Nsection.types) {
      quadrats.use <- (1:Nquad)[data.section==section.types[k]]
      scale.var.section[[k]] <- apply(as.matrix(wav.orig[quadrats.use,])^2,2,mean,na.rm=T)
    }
    names(scale.var.section) <- section.types
  }
  if(make.plot) {
    png(paste(file.title,".png",sep=""),res=300,width=20,height=20,unit="cm")
    plot(wav.orig)
    dev.off()
  }
  if(scale.var) {
    data.scale <- matrix(ncol=ncol(data.main),nrow=Nscales)
    data.scale[,1] <- wav.orig.scale
  }
  if(pos.var){
    data.pos <- matrix(ncol=ncol(data.main),nrow=length(wav.orig.pos))
    data.pos[,1] <- wav.orig.pos
  }
  if(scale.per.section){
    data.scale.section <- list()
    for (k in 1:Nsection.types) {
      data.scale.section[[k]] <- matrix(ncol=ncol(data.main), nrow=Nscales)
      data.scale.section[[k]][,1] <- scale.var.section[[k]]
    }
    names(data.scale.section) <- section.types
  }
  #The data.scale.section object stores the scale variance for each scale (scales are rows).
  if(save.cwt) write.table(wav.orig, paste(file.title,"_orig.txt",sep=""), sep=" ", dec=".", row.names=F, col.names=F)
  #Significance calculation based on the randomized data. Set significance=F if no randomized values for significance calculation are supplied.
  if(significance) {
    for(j in 2:Nperm) {
      #The analysis of the randomized data starts here
      foo <- data.main[,j]
      full.length <- length(foo)
      half.length <- round(length(foo)/2)
      if(zero.padding) foo <- c(rep(0,half.length),foo,rep(0,half.length))
      foo.wav <- wavCWT(foo, scale.range=c(1,scale.max), wavelet=wav.template)
      if(zero.padding) {
        wav.scale <- attr(foo.wav,"scale")
        wav.time <- 1:full.length
        wav.wavelet <- attr(foo.wav,"wavelet")
        wav.series <- attr(foo.wav,"series")
        wav.sampling.interval <- attr(foo.wav,"sampling.interval")
        wav.series.name <- attr(foo.wav,"series.name")
        wav.n.sample <- length(wav.time)
        wav.n.scale <- attr(foo.wav,"n.scale")
        wav.filter.arg <- attr(foo.wav,"filter.arg")
        #Remove the uninformative wavelet transform values
        foo.wav <- foo.wav[(half.length+1):(half.length+full.length),]
        #Return the attributes of the wavelet transform
        class(foo.wav) <- "wavCWT"
        attr(foo.wav,"scale") <- wav.scale
        attr(foo.wav,"time") <- wav.time
        attr(foo.wav,"wavelet") <- wav.wavelet
        attr(foo.wav,"series") <- wav.series
        attr(foo.wav,"sampling.interval") <- wav.sampling.interval
        attr(foo.wav,"series.names") <- wav.series.name
        attr(foo.wav,"n.sample") <- wav.n.sample
        attr(foo.wav,"n.scale") <- wav.n.scale
        attr(foo.wav,"filter.arg") <- wav.filter.arg
      }
      if(remove.COI) {
        for(k in 1:length(scales)) {
          index1 <- 1:(scales[k]*ifelse(wav.template=="gaussian2",2,1))
          index2 <- nrow(wav.orig)-index1+1
          indices <- c(index1,index2)
          foo.wav[indices,k] <- NA
        }
      }
      if(all(c(keep.all, save.cwt))) write.table(foo.wav,paste(file.title,"_",formatC(j,digits=3,format="d",flag=0),".txt",sep=""),sep=" ", dec=".", row.names=F, col.names=F)
      if(scale.var) foo.scale <- apply(as.matrix(foo.wav)^2,2,mean,na.rm=T)
      if(pos.var) foo.pos <- apply(as.matrix(foo.wav)^2,1,mean,na.rm=T)
      if(scale.per.section) {
        foo.var.section <- list()
        for (k in 1:Nsection.types) {
          quadrats.use <- (1:Nquad)[data.section==section.types[k]]
          foo.var.section[[k]] <- apply(as.matrix(foo.wav[quadrats.use,])^2,2,mean,na.rm=T)
        }
      }
      if(scale.var) data.scale[,j] <- foo.scale
      if(pos.var) data.pos[,j] <- foo.pos
      if(scale.per.section){
        for (k in 1:Nsection.types) {
          data.scale.section[[k]][,j] <- foo.var.section[[k]]
        }
      }
      if (print.loop) {
        if(j%%which.loop == 0) print(j)
      }
    }
    #This is the end of the loop. Each iteration calculated the wavelet transform on a randomized dataset and calculates scale variance, position variance, and scale variance per section (if specified by the user).
    #Now calculation of confidence intervals and effect sizes.    
    if (CI) { #Confidence intervals      
      if(pos.var) pos.var.CI <- t(apply(data.pos, 1, quantile, CIquant, na.rm=T))
      if(scale.var) scale.var.CI <- t(apply(data.scale, 1, quantile, CIquant, na.rm=T))
      if(scale.per.section) {
        scale.section.CI <- list()
        for (k in 1:Nsection.types) {
          scale.section.CI[[k]] <- t(apply(data.scale.section[[k]], 1, quantile, CIquant, na.rm=T))
        }
        names(scale.section.CI) <- section.types
      }
    }
    if(effect.size) {
      if(pos.var) {
        pos.var.avg <- apply(data.pos, 1, mean, na.rm=T)
        pos.var.sd <- apply(data.pos, 1, sd, na.rm=T)
        pos.var.effect <- (wav.orig.pos - pos.var.avg) / pos.var.sd
      }
      if(scale.var) {
        scale.var.avg <- apply(data.scale, 1, mean, na.rm=T)
        scale.var.sd <- apply(data.scale, 1, sd, na.rm=T)
        scale.var.effect <- (wav.orig.scale - scale.var.avg) / scale.var.sd
      }
      if(scale.per.section) {
        scale.section.avg <- list()
        scale.section.sd <- list()
        scale.section.effect <- list()
        for (k in 1:Nsection.types) {
          scale.section.avg[[k]] <- apply(data.scale.section[[k]], 1, mean, na.rm=T)
          scale.section.sd[[k]] <- apply(data.scale.section[[k]], 1, sd, na.rm=T)
          scale.section.effect[[k]] <- (scale.var.section[[k]] - scale.section.avg[[k]]) / scale.section.sd[[k]]
        }
        names(scale.section.avg) <- section.types
        names(scale.section.sd) <- section.types
        names(scale.section.effect) <- section.types
      }
    }
    #Now, export the final object
    result <- list()
    k <- 1
    if (pos.var) {
      result[[k]] <- wav.orig.pos
      names(result)[k] <- "PosVar_Obs"
      k <- k+1
      if(CI) {
        result[[k]] <- pos.var.CI
        names(result)[k] <- "PosVar_CI" 
        k <- k+1
      }
      if(effect.size) { 
        result[[k]] <- data.frame(pos.var.avg, pos.var.sd, pos.var.effect)
        names(result[[k]]) <- c("Mean", "SD", "Effect_size")
        names(result)[k]  <- "PosVar_Effect"
        k <- k+1
      }
    }
    if(any(c(scale.var, scale.per.section))) { 
      result[[k]] <- scales
      names(result)[[k]] <- "Scales"
      k <- k+1
    }
    if (scale.var) {
      result[[k]] <- wav.orig.scale
      names(result)[k] <- "ScaleVar_Obs"
      k <- k+1
      if(CI) { 
        result[[k]] <- scale.var.CI
        names(result)[k] <- "ScaleVar_CI" 
        k <- k+1
      }
      if(effect.size) {
        result[[k]] <- list(scale.var.avg, scale.var.sd, scale.var.effect)
        names(result[[k]]) <- c("Mean", "SD", "Effect_size")
        names(result)[k]  <- "ScaleVar_Effect"
        k <- k+1
      }
    }
    if (scale.per.section) {
      result[[k]] <- scale.var.section
      names(result)[k] <- "ScaleSection_Obs"
      k <- k+1
      if(CI) {
        result[[k]] <-scale.section.CI 
        names(result)[k] <- "ScaleSection_CI" 
        k <- k+1
      }
      if(effect.size) {
        result[[k]] <- list(scale.section.avg, scale.section.sd, scale.section.effect)
        names(result[[k]]) <- c("Mean", "SD", "Effect_size")
        names(result)[k]  <- "ScaleSection_Effect"
        k <- k+1
      }
    }
  } else {
    result <- list()
    k <- 1
    if(pos.var) {
      result[[k]] <- wav.orig.pos
      names(result)[k] <- "PosVar_Obs"
      k <- k+1
    }
    if (scale.var) {
    result[[k]] <- wav.orig.scale
      names(result)[k] <- "ScaleVar_Obs" 
      k <- k+1
    }
    if (scale.per.section) {
      result[[k]] <- scale.var.section
      names(result)[k] <- "ScaleSection_Obs"
      k <- k+1
    }
  }
  return(result)
}




wavCWTvarCIs.bivar <- function(x, y, significance=T, sections=NULL, zero.padding=T, remove.COI=T, test=F, scale.max=75, scale.var=T, pos.var=T, CI=T, CIquant=c(0.025, 0.975), effect.size=T, keep.all=F, print.loop=T, which.loop=1, save.cwt=F, file.title="my_beautiful_wavelet", wav.template="gaussian2", scale.per.section=T) {
  #This function does not save the original wavelet output, retaining only scale and position variance. Such functionality may be added in the future. For now, the keep.all, save.cwt and file.title arguments are not fully functional.
  #This may not work for wavelets other than "gaussian2" ("Mexican hat") and "haar" (Haar).
  if (scale.per.section) {
    section.lengths <- sections$End-sections$Start+1
    data.section <- rep(as.character(sections[1,1]), section.lengths[1])
    for (k in 2:nrow(sections)) {
      data.section <- c(data.section, rep(as.character(sections[k,1]), section.lengths[k]))
    }
    data.section <- as.factor(data.section)
    section.types <- unique(sections[,1])
    Nsection.types <- length(section.types)
  }
  data.main.x <- x
  data.main.y <- y
  Nperm <- ifelse(test,5,ncol(data.main.x))
  dimnames(data.main.x) <- NULL
  dimnames(data.main.y) <- NULL
  data.foo.x <- data.main.x[,1]
  data.foo.y <- data.main.y[,1]
  Nquad <- length(data.foo.x)
  full.length <- Nquad
  half.length <- round(Nquad/2)
  # Wavelet for x
  if(zero.padding) {
    data.foo.x <- c(rep(0,half.length),data.foo.x,rep(0,half.length)) #adds zeroes to the extremities of the transect
  }
  wav.orig.x <- wavCWT(data.foo.x,scale.range=c(1,scale.max),wavelet=wav.template)
  scales <- attr(wav.orig.x,"scale")
  Nscales <- length(scales)
  #Zero-padding adds zeroes to the extremities of the transect, to a length equal to half of the transect
  if(zero.padding) {
    #Save the attributes of the original wavelet transform
    wav.scale <- attr(wav.orig.x,"scale")
    wav.time <- 1:full.length
    wav.wavelet <- attr(wav.orig.x,"wavelet")
    wav.series <- attr(wav.orig.x,"series")
    wav.sampling.interval <- attr(wav.orig.x,"sampling.interval")
    wav.series.name <- attr(wav.orig.x,"series.name")
    wav.n.sample <- length(wav.time)
    wav.n.scale <- attr(wav.orig.x,"n.scale")
    wav.filter.arg <- attr(wav.orig.x,"filter.arg")
    #Remove parts of the wavelet transform corresponding to the artifficially added zeroes
    wav.orig.x <- wav.orig.x[(half.length+1):(half.length+full.length),]
    #Return the attributes to the corrected wavelet transform
    class(wav.orig.x) <- "wavCWT"
    attr(wav.orig.x,"scale") <- wav.scale
    attr(wav.orig.x,"time") <- wav.time
    attr(wav.orig.x,"wavelet") <- wav.wavelet
    attr(wav.orig.x,"series") <- wav.series
    attr(wav.orig.x,"sampling.interval") <- wav.sampling.interval
    attr(wav.orig.x,"series.names") <- wav.series.name
    attr(wav.orig.x,"n.sample") <- wav.n.sample
    attr(wav.orig.x,"n.scale") <- wav.n.scale
    attr(wav.orig.x,"filter.arg") <- wav.filter.arg
  }
  #The Cone Of Influence - COI - is the part of the wavelet transform that is affected by the transect's extremities, i.e. by quadrats for which no data has been collected. Inference made with these numbers is unreliable. The area affected is greater at larger scales, and therefore the COI has the shape of a, well, cone. Here I used an approach similar to that used by the PASSaGE software, which removes two quadrats for each scale for the MH wavelet and 1 quadrat for the Haar wavelet, from each side of the transect.
  if(remove.COI) {
    for(k in 1:length(scales)) {
      index1 <- 1:(scales[k]*ifelse(wav.template=="gaussian2",2,1))
      index2 <- nrow(wav.orig.x)-index1+1
      indices <- c(index1,index2)
      wav.orig.x[indices,k] <- NA
    }
  }
  # Wavelet for y 
  if(zero.padding) {
    data.foo.y <- c(rep(0,half.length),data.foo.y,rep(0,half.length)) #adds zeroes to the extremities of the transect
  }
  wav.orig.y <- wavCWT(data.foo.y,scale.range=c(1,scale.max),wavelet=wav.template)
  scales <- attr(wav.orig.y,"scale") # The scales should be the same for both wavelets.
  if (any(attr(wav.orig.x,"scale") != attr(wav.orig.y,"scale"))) stop("Scales differ - something's wrong")
  Nscales <- length(scales)
  #Zero-padding adds zeroes to the extremities of the transect, to a length equal to half of the transect
  if(zero.padding) {
    #Save the attributes of the original wavelet transform
    wav.scale <- attr(wav.orig.y,"scale")
    wav.time <- 1:full.length
    wav.wavelet <- attr(wav.orig.y,"wavelet")
    wav.series <- attr(wav.orig.y,"series")
    wav.sampling.interval <- attr(wav.orig.y,"sampling.interval")
    wav.series.name <- attr(wav.orig.y,"series.name")
    wav.n.sample <- length(wav.time)
    wav.n.scale <- attr(wav.orig.y,"n.scale")
    wav.filter.arg <- attr(wav.orig.y,"filter.arg")
    #Remove parts of the wavelet transform corresponding to the artifficially added zeroes
    wav.orig.y <- wav.orig.y[(half.length+1):(half.length+full.length),]
    #Return the attributes to the corrected wavelet transform
    class(wav.orig.y) <- "wavCWT"
    attr(wav.orig.y,"scale") <- wav.scale
    attr(wav.orig.y,"time") <- wav.time
    attr(wav.orig.y,"wavelet") <- wav.wavelet
    attr(wav.orig.y,"series") <- wav.series
    attr(wav.orig.y,"sampling.interval") <- wav.sampling.interval
    attr(wav.orig.y,"series.names") <- wav.series.name
    attr(wav.orig.y,"n.sample") <- wav.n.sample
    attr(wav.orig.y,"n.scale") <- wav.n.scale
    attr(wav.orig.y,"filter.arg") <- wav.filter.arg
  }
  #The Cone Of Influence - COI - is the part of the wavelet transform that is affected by the transect's extremities, i.e. by quadrats for which no data has been collected. Inference made with these numbers is unreliable. The area affected is greater at larger scales, and therefore the COI has the shape of a, well, cone. Here I used an approach similar to that used by the PASSaGE software, which removes two quadrats for each scale for the MH wavelet and 1 quadrat for the Haar wavelet, from each side of the transect.
  if(remove.COI) {
    for(k in 1:length(scales)) {
      index1 <- 1:(scales[k]*ifelse(wav.template=="gaussian2",2,1))
      index2 <- nrow(wav.orig.y)-index1+1
      indices <- c(index1,index2)
      wav.orig.y[indices,k] <- NA
    }
  }

  wav.covar <- as.matrix(wav.orig.x) * as.matrix(wav.orig.y)
  if(scale.var) wav.bivar.scale <- apply(wav.covar,2,mean, na.rm=T)
  if(pos.var) wav.bivar.pos <- apply(wav.covar,1,mean, na.rm=T)
  if(scale.per.section) {
    scale.var.section <- list()
    for (k in 1:Nsection.types) {
      quadrats.use <- (1:Nquad)[data.section==section.types[k]]
      scale.var.section[[k]] <- apply(wav.covar[quadrats.use,],2,mean,na.rm=T)
    }
    names(scale.var.section) <- section.types
  }
  if(scale.var) {
    data.scale <- matrix(ncol=ncol(data.main.x),nrow=Nscales)
    data.scale[,1] <- wav.bivar.scale
  }
  if(pos.var){
    data.pos <- matrix(ncol=ncol(data.main.x),nrow=length(wav.bivar.pos))
    data.pos[,1] <- wav.bivar.pos
  }
  if(scale.per.section){
    data.scale.section <- list()
    for (k in 1:Nsection.types) {
      data.scale.section[[k]] <- matrix(ncol=ncol(data.main.x), nrow=Nscales)
      data.scale.section[[k]][,1] <- scale.var.section[[k]]
    }
    names(data.scale.section) <- section.types
  }
  #The data.scale.section object stores the scale variance for each scale (scales are rows).
  #if(save.cwt) write.table(wav.orig.x, paste(file.title,"_x_orig.txt",sep=""), sep=" ", dec=".", row.names=F, col.names=F)
  #if(save.cwt) write.table(wav.orig.y, paste(file.title,"_y_orig.txt",sep=""), sep=" ", dec=".", row.names=F, col.names=F)
  #Significance calculation based on the randomized data. Set significance=F if no randomized values for significance calculation are supplied.
  if(significance) {
    for(j in 2:Nperm) {
      #The analysis of the randomized data starts here
      foo.x <- data.main.x[,j]
      full.length <- length(foo.x)
      half.length <- round(length(foo.x)/2)
      if(zero.padding) foo.x <- c(rep(0,half.length),foo.x,rep(0,half.length))
      foo.x.wav <- wavCWT(foo.x, scale.range=c(1,scale.max), wavelet=wav.template)
      if(zero.padding) {
        wav.scale <- attr(foo.x.wav,"scale")
        wav.time <- 1:full.length
        wav.wavelet <- attr(foo.x.wav,"wavelet")
        wav.series <- attr(foo.x.wav,"series")
        wav.sampling.interval <- attr(foo.x.wav,"sampling.interval")
        wav.series.name <- attr(foo.x.wav,"series.name")
        wav.n.sample <- length(wav.time)
        wav.n.scale <- attr(foo.x.wav,"n.scale")
        wav.filter.arg <- attr(foo.x.wav,"filter.arg")
        #Remove the uninformative wavelet transform values
        foo.x.wav <- foo.x.wav[(half.length+1):(half.length+full.length),]
        #Return the attributes of the wavelet transform
        class(foo.x.wav) <- "wavCWT"
        attr(foo.x.wav,"scale") <- wav.scale
        attr(foo.x.wav,"time") <- wav.time
        attr(foo.x.wav,"wavelet") <- wav.wavelet
        attr(foo.x.wav,"series") <- wav.series
        attr(foo.x.wav,"sampling.interval") <- wav.sampling.interval
        attr(foo.x.wav,"series.names") <- wav.series.name
        attr(foo.x.wav,"n.sample") <- wav.n.sample
        attr(foo.x.wav,"n.scale") <- wav.n.scale
        attr(foo.x.wav,"filter.arg") <- wav.filter.arg
      }
      if(remove.COI) {
        for(k in 1:length(scales)) {
          index1 <- 1:(scales[k]*ifelse(wav.template=="gaussian2",2,1))
          index2 <- nrow(wav.orig.x)-index1+1
          indices <- c(index1,index2)
          foo.x.wav[indices,k] <- NA
        }
      }
      #if(all(c(keep.all, save.cwt))) write.table(foo.wav,paste(file.title,"_x_",formatC(j,digits=3,format="d",flag=0),".txt",sep=""),sep=" ", dec=".", row.names=F, col.names=F)

      foo.y <- data.main.y[,j]
      full.length <- length(foo.y)
      half.length <- round(length(foo.y)/2)
      if(zero.padding) foo.y <- c(rep(0,half.length),foo.y,rep(0,half.length))
      foo.y.wav <- wavCWT(foo.y, scale.range=c(1,scale.max), wavelet=wav.template)
      if(zero.padding) {
        wav.scale <- attr(foo.y.wav,"scale")
        wav.time <- 1:full.length
        wav.wavelet <- attr(foo.y.wav,"wavelet")
        wav.series <- attr(foo.y.wav,"series")
        wav.sampling.interval <- attr(foo.y.wav,"sampling.interval")
        wav.series.name <- attr(foo.y.wav,"series.name")
        wav.n.sample <- length(wav.time)
        wav.n.scale <- attr(foo.y.wav,"n.scale")
        wav.filter.arg <- attr(foo.y.wav,"filter.arg")
        #Remove the uninformative wavelet transform values
        foo.y.wav <- foo.y.wav[(half.length+1):(half.length+full.length),]
        #Return the attributes of the wavelet transform
        class(foo.y.wav) <- "wavCWT"
        attr(foo.y.wav,"scale") <- wav.scale
        attr(foo.y.wav,"time") <- wav.time
        attr(foo.y.wav,"wavelet") <- wav.wavelet
        attr(foo.y.wav,"series") <- wav.series
        attr(foo.y.wav,"sampling.interval") <- wav.sampling.interval
        attr(foo.y.wav,"series.names") <- wav.series.name
        attr(foo.y.wav,"n.sample") <- wav.n.sample
        attr(foo.y.wav,"n.scale") <- wav.n.scale
        attr(foo.y.wav,"filter.arg") <- wav.filter.arg
      }
      if(remove.COI) {
        for(k in 1:length(scales)) {
          index1 <- 1:(scales[k]*ifelse(wav.template=="gaussian2",2,1))
          index2 <- nrow(wav.orig.y)-index1+1
          indices <- c(index1,index2)
          foo.y.wav[indices,k] <- NA
        }
      }
      #if(all(c(keep.all, save.cwt))) write.table(foo.wav,paste(file.title,"_y_",formatC(j,digits=3,format="d",flag=0),".txt",sep=""),sep=" ", dec=".", row.names=F, col.names=F)

      foo.covar <- as.matrix(foo.x.wav) * as.matrix(foo.y.wav)

      if(scale.var) foo.scale <- apply(foo.covar,2,mean,na.rm=T)
      if(pos.var) foo.pos <- apply(foo.covar,1,mean,na.rm=T)
      if(scale.per.section) {
        foo.var.section <- list()
        for (k in 1:Nsection.types) {
          quadrats.use <- (1:Nquad)[data.section==section.types[k]]
          foo.var.section[[k]] <- apply(as.matrix(foo.covar[quadrats.use,]),2,mean,na.rm=T)
        }
      }
      if(scale.var) data.scale[,j] <- foo.scale
      if(pos.var) data.pos[,j] <- foo.pos
      if(scale.per.section){
        for (k in 1:Nsection.types) {
          data.scale.section[[k]][,j] <- foo.var.section[[k]]
        }
      }
      if (print.loop) {
        if(j%%which.loop == 0) print(j)
      }
    }
    #This is the end of the loop. Each iteration calculated the wavelet transform on a randomized dataset and calculated scale covariance, position covariance, and scale covariance per section (if specified by the user).
    #Now comes the calculation of confidence intervals and effect sizes.
    if (CI) { #Confidence intervals
      if(pos.var) pos.var.CI <- t(apply(data.pos, 1, quantile, CIquant, na.rm=T))
      if(scale.var) scale.var.CI <- t(apply(data.scale, 1, quantile, CIquant, na.rm=T))
      if(scale.per.section) {
        scale.section.CI <- list()
        for (k in 1:Nsection.types) {
          scale.section.CI[[k]] <- t(apply(data.scale.section[[k]], 1, quantile, CIquant, na.rm=T))
        }
        names(scale.section.CI) <- section.types
      }
    }
    if(effect.size) {
      if(pos.var) {
        pos.var.avg <- apply(data.pos, 1, mean, na.rm=T)
        pos.var.sd <- apply(data.pos, 1, sd, na.rm=T)
        pos.var.effect <- (wav.orig.pos - pos.var.avg) / pos.var.sd
      }
      if(scale.var) {
        scale.var.avg <- apply(data.scale, 1, mean, na.rm=T)
        scale.var.sd <- apply(data.scale, 1, sd, na.rm=T)
        scale.var.effect <- (wav.orig.scale - scale.var.avg) / scale.var.sd
      }
      if(scale.per.section) {
        scale.section.avg <- list()
        scale.section.sd <- list()
        scale.section.effect <- list()
        for (k in 1:Nsection.types) {
          scale.section.avg[[k]] <- apply(data.scale.section[[k]], 1, mean, na.rm=T)
          scale.section.sd[[k]] <- apply(data.scale.section[[k]], 1, sd, na.rm=T)
          scale.section.effect[[k]] <- (scale.var.section[[k]] - scale.section.avg[[k]]) / scale.section.sd[[k]]
        }
        names(scale.section.avg) <- section.types
        names(scale.section.sd) <- section.types
        names(scale.section.effect) <- section.types
      }
    }
    #Now, export the final object
    result <- list()
    k <- 1
    if (pos.var) {
      result[[k]] <- wav.bivar.pos
      names(result)[k] <- "PosVar_Obs"
      k <- k+1
      if(CI) {
        result[[k]] <- pos.var.CI
        names(result)[k] <- "PosVar_CI" 
        k <- k+1
      }
      if(effect.size) { 
        result[[k]] <- data.frame(pos.var.avg, pos.var.sd, pos.var.effect)
        names(result[[k]]) <- c("Mean", "SD", "Effect_size")
        names(result)[k] <- "PosVar_Effect"
        k <- k+1
      }
    }
    if(any(c(scale.var, scale.per.section))) { 
      result[[k]] <- scales
      names(result)[[k]] <- "Scales"
      k <- k+1
    }
    if (scale.var) {
      result[[k]] <- wav.bivar.scale
      names(result)[k] <- "ScaleVar_Obs"
      k <- k+1
      if(CI) { 
        result[[k]] <- scale.var.CI
        names(result)[k] <- "ScaleVar_CI" 
        k <- k+1
      }
      if(effect.size) {
        result[[k]] <- list(scale.var.avg, scale.var.sd, scale.var.effect)
        names(result[[k]]) <- c("Mean", "SD", "Effect_size")
        names(result)[k] <- "ScaleVar_Effect"
        k <- k+1
      }
    }
    if (scale.per.section) {
      result[[k]] <- scale.var.section
      names(result)[k] <- "ScaleSection_Obs"
      k <- k+1
      if(CI) {
        result[[k]] <-scale.section.CI 
        names(result)[k] <- "ScaleSection_CI" 
        k <- k+1
      }
      if(effect.size) {
        result[[k]] <- list(scale.section.avg, scale.section.sd, scale.section.effect)
        names(result[[k]]) <- c("Mean", "SD", "Effect_size")
        names(result)[k] <- "ScaleSection_Effect"
        k <- k+1
      }
    }
  } else {
    result <- list()
    k <- 1
    if(pos.var) {
      result[[k]] <- wav.bivar.pos
      names(result)[k] <- "PosVar_Obs"
      k <- k+1
    }
    if (scale.var) {
    result[[k]] <- wav.bivar.scale
      names(result)[k] <- "ScaleVar_Obs" 
      k <- k+1
    }
    if (scale.per.section) {
      result[[k]] <- scale.var.section
      names(result)[k] <- "ScaleSection_Obs"
      k <- k+1
    }
  }
  return(result)
}


### Part II - Insert and organize the data

setwd("/home/pavel/Profissional/Pesquisa/MyPapers-Revising/EIgrasses/forR")
data.itirapina <- read.table("itirapina_data.txt", header=T, sep="\t")
data.ufscar1 <- read.table("ufscar1_data.txt", header=T, sep="\t")
data.ufscar2 <- read.table("ufscar2_data.txt", header=T, sep="\t")
patches.itirapina <- read.table("itirapina_patches.txt", header=T, sep="\t")
patches.ufscar1 <- read.table("ufscar1_patches.txt", header=T, sep="\t")
patches.ufscar2 <- read.table("ufscar2_patches.txt", header=T, sep="\t")
elevation.itirapina <- read.table("itirapina_elevation.txt", header=T)
elevation.ufscar1 <- read.table("ufscar1_elevation.txt", header=T)
elevation.ufscar2 <- read.table("ufscar2_elevation.txt", header=T)

data.itirapina <- replace.cover(data.itirapina)
data.ufscar1 <- replace.cover(data.ufscar1)
data.ufscar2 <- replace.cover(data.ufscar2)


# Prepare some objects

patches.itirapina.rep <- rep.patch(patches=patches.itirapina)
patches.ufscar1.rep <- rep.patch(patches=patches.ufscar1)
patches.ufscar2.rep <- rep.patch(patches=patches.ufscar2)

patch.lengths.itirapina <- patches.itirapina$End-patches.itirapina$Start+1
patches2.itirapina <- rep(as.character(patches.itirapina[1,1]), patch.lengths.itirapina[1] )
for (i in 2:nrow(patches.itirapina)) {
  patches2.itirapina <- c(patches2.itirapina, rep(as.character(patches.itirapina[i,1]), patch.lengths.itirapina[i] ) )
  }
dist.itirapina <- 1:nrow(data.itirapina)
dist.firebreaks.itirapina <- numeric()
tests <- logical()
for(i in 1:nrow(patches.itirapina)) tests[i]=test.firebreak(patches.itirapina[i,1])
patches.firebreak.itirapina <- patches.itirapina[tests,]
dist.firebreak.itirapina <- numeric()
for (i in 1:nrow(patches.firebreak.itirapina)) {
  foo <- patches.firebreak.itirapina[i,2]:patches.firebreak.itirapina[i,3]
  if (i==1) dist.firebreak.itirapina <- foo else dist.firebreak.itirapina <- c(dist.firebreak.itirapina, foo)
}
dist.toFirebreak.itirapina <- numeric()
for(i in 1:nrow(data.itirapina)) {
  foo <- dist.itirapina[i]
  if (any(foo==dist.firebreak.itirapina)) dist.toFirebreak.itirapina[i]=NA else dist.toFirebreak.itirapina[i]=min(abs(foo-dist.firebreak.itirapina))
}
data.itirapina$Dist <- dist.itirapina
data.itirapina$DistToFirebreak <- dist.toFirebreak.itirapina
data.itirapina$Vegetation <- as.factor(patches2.itirapina)
data.itirapina$Elevation <- elevation.itirapina[,1]


patch.lengths.ufscar1 <- patches.ufscar1$End-patches.ufscar1$Start+1
patches2.ufscar1 <- rep(as.character(patches.ufscar1[1,1]), patch.lengths.ufscar1[1] )
for (i in 2:nrow(patches.ufscar1)) {
  patches2.ufscar1 <- c(patches2.ufscar1, rep(as.character(patches.ufscar1[i,1]), patch.lengths.ufscar1[i] ) )
  }
dist.ufscar1 <- 1:nrow(data.ufscar1)
dist.firebreaks.ufscar1 <- numeric()
tests <- logical()
for(i in 1:nrow(patches.ufscar1)) tests[i]=test.firebreak(patches.ufscar1[i,1])
patches.firebreak.ufscar1 <- patches.ufscar1[tests,]
dist.firebreak.ufscar1 <- numeric()
for (i in 1:nrow(patches.firebreak.ufscar1)) {
  foo <- patches.firebreak.ufscar1[i,2]:patches.firebreak.ufscar1[i,3]
  if (i==1) dist.firebreak.ufscar1 <- foo else dist.firebreak.ufscar1 <- c(dist.firebreak.ufscar1, foo)
}
dist.toFirebreak.ufscar1 <- numeric()
for(i in 1:nrow(data.ufscar1)) {
  foo <- dist.ufscar1[i]
  if (any(foo==dist.firebreak.ufscar1)) dist.toFirebreak.ufscar1[i]=NA else dist.toFirebreak.ufscar1[i]=min(abs(foo-dist.firebreak.ufscar1))
}
data.ufscar1$Dist <- dist.ufscar1
data.ufscar1$DistToFirebreak <- dist.toFirebreak.ufscar1
data.ufscar1$Vegetation <- as.factor(patches2.ufscar1)
data.ufscar1$Elevation <- elevation.ufscar1[,1]

patch.lengths.ufscar2 <- patches.ufscar2$End-patches.ufscar2$Start+1
patches2.ufscar2 <- rep(as.character(patches.ufscar2[1,1]), patch.lengths.ufscar2[1] )
for (i in 2:nrow(patches.ufscar2)) {
  patches2.ufscar2 <- c(patches2.ufscar2, rep(as.character(patches.ufscar2[i,1]), patch.lengths.ufscar2[i] ) )
  }
dist.ufscar2 <- 1:nrow(data.ufscar2)
dist.firebreaks.ufscar2 <- numeric()
tests <- logical()
for(i in 1:nrow(patches.ufscar2)) tests[i]=test.firebreak(patches.ufscar2[i,1])
patches.firebreak.ufscar2 <- patches.ufscar2[tests,]
dist.firebreak.ufscar2 <- numeric()
for (i in 1:nrow(patches.firebreak.ufscar2)) {
  foo <- patches.firebreak.ufscar2[i,2]:patches.firebreak.ufscar2[i,3]
  if (i==1) dist.firebreak.ufscar2 <- foo else dist.firebreak.ufscar2 <- c(dist.firebreak.ufscar2, foo)
}
dist.toFirebreak.ufscar2 <- numeric()
for(i in 1:nrow(data.ufscar2)) {
  foo <- dist.ufscar2[i]
  if (any(foo==dist.firebreak.ufscar2)) dist.toFirebreak.ufscar2[i]=NA else dist.toFirebreak.ufscar2[i]=min(abs(foo-dist.firebreak.ufscar2))
}
data.ufscar2$Dist <- dist.ufscar2
data.ufscar2$DistToFirebreak <- dist.toFirebreak.ufscar2
data.ufscar2$Vegetation <- as.factor(patches2.ufscar2)
data.ufscar2$Elevation <- elevation.ufscar2[,1]

# Remove firebreaks
data2.itirapina <- subset(data.itirapina, !is.na(DistToFirebreak))
data2.itirapina$Vegetation <- as.factor(as.character(data2.itirapina$Vegetation))
data2.ufscar1 <- subset(data.ufscar1, !is.na(DistToFirebreak))
data2.ufscar1$Vegetation <- as.factor(as.character(data2.ufscar1$Vegetation))
data2.ufscar2 <- subset(data.ufscar2, !is.na(DistToFirebreak))
data2.ufscar2$Vegetation <- as.factor(as.character(data2.ufscar2$Vegetation))

# Differentiate between response and explanatory variables

data2.itirapina.r <- data2.itirapina[,1:4]
data2.itirapina.e <- data2.itirapina[,6:8]
data2.ufscar1.r <- data2.ufscar1[,1:4]
data2.ufscar1.e <- data2.ufscar1[,6:8]
data2.ufscar2.r <- data2.ufscar2[,1:4]
data2.ufscar2.e <- data2.ufscar2[,6:8]


### Part III - Make the general figures: altimetric profile and grass cover along the transects.

patches.itirapina.plot <- patches.itirapina
patches.ufscar1.plot <- patches.ufscar1
patches.ufscar2.plot <- patches.ufscar2


patches.itirapina.plot$Color <- c("white","gray80","white","gray80","white","gray40","gray30","gray40","gray80","white","gray80")
patches.ufscar1.plot$Color <- c("gray80","white","gray80","white","gray60","white","gray60","white","gray60","white")
patches.ufscar2.plot$Color <- c("white","gray60","white","gray60","white","gray40","gray30","gray60","white")

alt.itirapina <- data.itirapina$Elevation
alt.ufscar1 <- data.ufscar1$Elevation
alt.ufscar2 <- data.ufscar2$Elevation


range.itirapina <- range(alt.itirapina)
range.ufscar1 <- range(alt.ufscar1)
range.ufscar2 <- range(alt.ufscar2)

min.itirapina <- min(alt.itirapina)
min.ufscar1 <- min(alt.ufscar1)
min.ufscar2 <- min(alt.ufscar2)

alt.diff <- max(c(diff(range(alt.ufscar1)), diff(range(alt.ufscar2)),diff(range(alt.itirapina))))

range2.itirapina <- c(min.itirapina, min.itirapina+alt.diff)
range2.ufscar1 <- c(min.ufscar1, min.ufscar1+alt.diff)
range2.ufscar2 <- c(min.ufscar2, min.ufscar2+alt.diff)

Nquad.itirapina <- nrow(data.itirapina)
Nquad.ufscar1 <- nrow(data.ufscar1)
Nquad.ufscar2 <- nrow(data.ufscar2)

props <- c(Nquad.itirapina, Nquad.ufscar1, Nquad.ufscar2) / sum(c(Nquad.itirapina, Nquad.ufscar1, Nquad.ufscar2))
props2 <- c(0.04,props[1],0.04,props[2],0.04,props[3])

setwd("/home/pavel/Profissional/Pesquisa/MyPapers-Working/EIgrasses/figures")
png(filename="EIgrasses_fig1e_2.png", height=8, width=37, unit="cm", res=300)                                   
layout(mat=matrix(c(1,2,3,4,5,6),nrow=1), widths=props2)  ### Define the size of each plot as proportional to transect length
par(mar=c(3,0,3,0), oma=c(2,3,1,1))

plot(0,type="n",xaxt="n", yaxt="n", xlab="", ylab="",bty="n")

plot(alt.itirapina, type="n", xlab="", ylab="", cex.lab=1.4, xaxt="n", ylim=range2.itirapina, xaxs="r", bty="o", main="e) Itirapina - I1", cex.main=1.6)
for(i in 1:nrow(patches.itirapina.plot)){
  polygon(x=c(patches.itirapina.plot[i,2],patches.itirapina.plot[i,2],patches.itirapina.plot[i,3]+1,patches.itirapina.plot[i,3]+1),y=c(range2.itirapina,range2.itirapina[c(2,1)]),col=patches.itirapina.plot[i,4], border=NA)
}
axis(side=1, at=c(0,200, 400, 600))
points(alt.itirapina, type="l", lwd=2)

plot(0,type="n",xaxt="n", yaxt="n", xlab="", ylab="",bty="n")

plot(alt.ufscar1, type="n", xlab="", ylab="", xaxt="n", ylim=range2.ufscar1, xaxs="r", bty="o", main="f) São Carlos - S1", cex.main=1.6)
for(i in 1:nrow(patches.ufscar1.plot)){
  polygon(x=c(patches.ufscar1.plot[i,2],patches.ufscar1.plot[i,2],patches.ufscar1.plot[i,3]+1,patches.ufscar1.plot[i,3]+1),y=c(range2.ufscar1,range2.ufscar1[c(2,1)]),col=patches.ufscar1.plot[i,4], border=NA)
}
points(alt.ufscar1, type="l", lwd=2)
axis(side=1, at=c(0,400,800,1200))

plot(0,type="n",xaxt="n", yaxt="n", xlab="", ylab="",bty="n")

plot(alt.ufscar2, type="n", xlab="", ylab="", xaxt="n", ylim=range2.ufscar2, xaxs="r", bty="o", main="g) São Carlos - S2", cex.main=1.6)
for(i in 1:nrow(patches.ufscar2.plot)){
  polygon(x=c(patches.ufscar2.plot[i,2],patches.ufscar2.plot[i,2],patches.ufscar2.plot[i,3]+1,patches.ufscar2.plot[i,3]+1),y=c(range2.ufscar2,range2.ufscar2[c(2,1)]),col=patches.ufscar2.plot[i,4], border=NA)
}
points(alt.ufscar2, type="l", lwd=2)
axis(side=1, at=c(0,100,200))

mtext(side=1, text="Distance along transect (m)", outer=T, line=1, cex=1.3)
mtext(side=2, text="Elevation (m a. s. l.)", outer=T, line=1, cex=1.3)

dev.off()





setwd("/home/pavel/Profissional/Pesquisa/MyPapers-Revising/EIgrasses/Submitted3")

png(filename="EIgrasses_fig2_distribution_subm3.png", height=20, width=22, unit="cm",res=300)

layout(mat=matrix(1:24, nrow=4, ncol=6), widths=props2) # Define the width of each plot as proportional to transect length
par(mar=c(1,0,1,0), oma=c(3,6,3.5,1))

plot(0,type="n",xaxt="n", yaxt="n", xlab="", ylab="",bty="n")
plot(0,type="n",xaxt="n", yaxt="n", xlab="", ylab="",bty="n")
plot(0,type="n",xaxt="n", yaxt="n", xlab="", ylab="",bty="n")
plot(0,type="n",xaxt="n", yaxt="n", xlab="", ylab="",bty="n")


plot(data.itirapina$Urochloa, type="n", xlab="", ylab="", xaxt="n", ylim=c(0,1), xaxs="r", yaxt="n")
for(i in 1:nrow(patches.itirapina.plot)){
  polygon(x=c(patches.itirapina.plot[i,2],patches.itirapina.plot[i,2],patches.itirapina.plot[i,3]+1,patches.itirapina.plot[i,3]+1),y=c(0,1,1,0),col=patches.itirapina.plot[i,4], border=NA)
}
points(data.itirapina$Urochloa, type="h")
mtext(side=3,text="a)",adj=0)
axis(side=2, at=c(0,0.25, 0.50, 0.75, 1), labels=c(0,25,50,75,100), las=1)

plot(data.itirapina$Melinis, type="n", xlab="", ylab="", xaxt="n", ylim=c(0,1), xaxs="r", yaxt="n")
for(i in 1:nrow(patches.itirapina.plot)){
  polygon(x=c(patches.itirapina.plot[i,2],patches.itirapina.plot[i,2],patches.itirapina.plot[i,3]+1,patches.itirapina.plot[i,3]+1),y=c(0,1,1,0),col=patches.itirapina.plot[i,4], border=NA)
}
points(data.itirapina$Melinis, type="h")
mtext(side=3,text="b)",adj=0)
axis(side=2, at=c(0,0.25, 0.50, 0.75, 1), labels=c(0,25,50,75,100), las=1)


plot(data.itirapina$Grasses, type="n", xlab="", ylab="", xaxt="n", ylim=c(0,1), xaxs="r", yaxt="n")
for(i in 1:nrow(patches.itirapina.plot)){
  polygon(x=c(patches.itirapina.plot[i,2],patches.itirapina.plot[i,2],patches.itirapina.plot[i,3]+1,patches.itirapina.plot[i,3]+1),y=c(0,1,1,0),col=patches.itirapina.plot[i,4], border=NA)
}
points(data.itirapina$Grasses, type="h")
mtext(side=3,text="c)",adj=0)
axis(side=2, at=c(0,0.25, 0.50, 0.75, 1), labels=c(0,25,50,75,100), las=1)

plot(data.itirapina$Sedges, type="n", xlab="", ylab="", xaxt="n", ylim=c(0,1), xaxs="r", yaxt="n")
for(i in 1:nrow(patches.itirapina.plot)){
  polygon(x=c(patches.itirapina.plot[i,2],patches.itirapina.plot[i,2],patches.itirapina.plot[i,3]+1,patches.itirapina.plot[i,3]+1),y=c(0,1,1,0),col=patches.itirapina.plot[i,4], border=NA)
}
points(data.itirapina$Sedges, type="h")
axis(side=1)
mtext(side=3,text="d)",adj=0)
axis(side=2, at=c(0,0.25, 0.50, 0.75, 1), labels=c(0,25,50,75,100), las=1)


plot(0,type="n",xaxt="n", yaxt="n", xlab="", ylab="",bty="n")
plot(0,type="n",xaxt="n", yaxt="n", xlab="", ylab="",bty="n")
plot(0,type="n",xaxt="n", yaxt="n", xlab="", ylab="",bty="n")
plot(0,type="n",xaxt="n", yaxt="n", xlab="", ylab="",bty="n")


par(yaxt="n")
plot(data.ufscar1$Urochloa, type="n", xlab="", ylab="", xaxt="n", ylim=c(0,1), xaxs="r")
for(i in 1:nrow(patches.ufscar1.plot)){
  polygon(x=c(patches.ufscar1.plot[i,2],patches.ufscar1.plot[i,2],patches.ufscar1.plot[i,3]+1,patches.ufscar1.plot[i,3]+1),y=c(0,1,1,0),col=patches.ufscar1.plot[i,4], border=NA)
}
points(data.ufscar1$Urochloa, type="h")
mtext(side=3,text="e)",adj=0)

plot(data.ufscar1$Melinis, type="n", xlab="", ylab="", xaxt="n", ylim=c(0,1), xaxs="r")
for(i in 1:nrow(patches.ufscar1.plot)){
  polygon(x=c(patches.ufscar1.plot[i,2],patches.ufscar1.plot[i,2],patches.ufscar1.plot[i,3]+1,patches.ufscar1.plot[i,3]+1),y=c(0,1,1,0),col=patches.ufscar1.plot[i,4], border=NA)
}
points(data.ufscar1$Melinis, type="h")
mtext(side=3,text="f)",adj=0)

plot(data.ufscar1$Grasses, type="n", xlab="", ylab="", xaxt="n", ylim=c(0,1), xaxs="r")
for(i in 1:nrow(patches.ufscar1.plot)){
  polygon(x=c(patches.ufscar1.plot[i,2],patches.ufscar1.plot[i,2],patches.ufscar1.plot[i,3]+1,patches.ufscar1.plot[i,3]+1),y=c(0,1,1,0),col=patches.ufscar1.plot[i,4], border=NA)
}
points(data.ufscar1$Grasses, type="h")
mtext(side=3,text="g)",adj=0)

plot(data.ufscar1$Sedges, type="n", xlab="", ylab="", xaxt="n", ylim=c(0,1), xaxs="r")
for(i in 1:nrow(patches.ufscar1.plot)){
  polygon(x=c(patches.ufscar1.plot[i,2],patches.ufscar1.plot[i,2],patches.ufscar1.plot[i,3]+1,patches.ufscar1.plot[i,3]+1),y=c(0,1,1,0),col=patches.ufscar1.plot[i,4], border=NA)
}
points(data.ufscar1$Sedges, type="h")
axis(side=1)
mtext(side=3,text="h)",adj=0)


plot(0,type="n",xaxt="n", yaxt="n", xlab="", ylab="",bty="n")
plot(0,type="n",xaxt="n", yaxt="n", xlab="", ylab="",bty="n")
plot(0,type="n",xaxt="n", yaxt="n", xlab="", ylab="",bty="n")
plot(0,type="n",xaxt="n", yaxt="n", xlab="", ylab="",bty="n")


plot(data.ufscar2$Urochloa, type="n", xlab="", ylab="", xaxt="n", ylim=c(0,1), xaxs="r")
for(i in 1:nrow(patches.ufscar2.plot)){
  polygon(x=c(patches.ufscar2.plot[i,2],patches.ufscar2.plot[i,2],patches.ufscar2.plot[i,3]+1,patches.ufscar2.plot[i,3]+1),y=c(0,1,1,0),col=patches.ufscar2.plot[i,4], border=NA)
}
points(data.ufscar2$Urochloa, type="h")
mtext(side=3,text="i)",adj=0)

plot(data.ufscar2$Melinis, type="n", xlab="", ylab="", xaxt="n", ylim=c(0,1), xaxs="r")
for(i in 1:nrow(patches.ufscar2.plot)){
  polygon(x=c(patches.ufscar2.plot[i,2],patches.ufscar2.plot[i,2],patches.ufscar2.plot[i,3]+1,patches.ufscar2.plot[i,3]+1),y=c(0,1,1,0),col=patches.ufscar2.plot[i,4], border=NA)
}
points(data.ufscar2$Melinis, type="h")
mtext(side=3,text="j)",adj=0)

plot(data.ufscar2$Grasses, type="n", xlab="", ylab="", xaxt="n", ylim=c(0,1), xaxs="r")
for(i in 1:nrow(patches.ufscar2.plot)){
  polygon(x=c(patches.ufscar2.plot[i,2],patches.ufscar2.plot[i,2],patches.ufscar2.plot[i,3]+1,patches.ufscar2.plot[i,3]+1),y=c(0,1,1,0),col=patches.ufscar2.plot[i,4], border=NA)
}
points(data.ufscar2$Grasses, type="h")
mtext(side=3,text="k)",adj=0)

plot(data.ufscar2$Sedges, type="n", xlab="", ylab="", xaxt="n", ylim=c(0,1), xaxs="r")
for(i in 1:nrow(patches.ufscar2.plot)){
  polygon(x=c(patches.ufscar2.plot[i,2],patches.ufscar2.plot[i,2],patches.ufscar2.plot[i,3]+1,patches.ufscar2.plot[i,3]+1),y=c(0,1,1,0),col=patches.ufscar2.plot[i,4], border=NA)
}
points(data.ufscar2$Sedges, type="h")
axis(side=1)
mtext(side=3,text="l)",adj=0)

mtext(side=1, outer=T, text="Distance along the transect (m)", line=1.5)
mtext(side=2, outer=T, text="Cover (%)", line=4)
mtext(side=2, outer=T, text=expression("sedges", "grasses", italic("minutiflora"), italic("decumbens")),
  at=c(0.125, 0.375, 0.625, 0.875), line=1, adj=0.5)
mtext(side=2, outer=T, text=expression("Native", "Native", italic("Melinis"), italic("Urochloa")),
  at=c(0.125, 0.375, 0.625, 0.875), line=2.5, adj=0.5)

mtext(side=3, text=c("Itirapina - I1","São Carlos - S1", "São Carlos - S2"), at=c(0.15, 0.62, 0.94), outer=T, line=1.5)

dev.off()



### Part 1 - What is the variation in cover between vegetation types?

# Simulate the null models

Nperm <- 5000

rand.MC1h.itirapina <- lapply(data.itirapina[,1:4], MC1h, which.loop=1000, Nperm=Nperm)
rand.MC1h.ufscar1 <- lapply(data.ufscar1[,1:4], MC1h, which.loop=1000, Nperm=Nperm)
rand.MC1h.ufscar2 <- lapply(data.ufscar2[,1:4], MC1h, which.loop=1000, Nperm=Nperm)

patch.types.itirapina <- sort(unique(patches.itirapina.rep))
patch.types.ufscar1 <- sort(unique(patches.ufscar1.rep))
patch.types.ufscar2 <- sort(unique(patches.ufscar2.rep))

### Calculate 95% confidence intervals for each null model for each vegegation type (or patch.type)

freqs.itirapina.rand <- list()
covers.itirapina.rand <- list()
for(i in 1:4) {
	freqs.itirapina.rand[[i]] <- matrix(nrow=length(patch.types.itirapina), ncol=Nperm)
	covers.itirapina.rand[[i]] <- matrix(nrow=length(patch.types.itirapina), ncol=Nperm)	
	row.names(freqs.itirapina.rand[[i]]) <- patch.types.itirapina
	row.names(covers.itirapina.rand[[i]]) <- patch.types.itirapina
	for(j in 1:Nperm) {
		freqs.itirapina.rand[[i]][,j] <- aggregate(rand.MC1h.itirapina[[i]][,j], by=list(section=patches.itirapina.rep), FUN=calc.freq)$x
		covers.itirapina.rand[[i]][,j] <- aggregate(rand.MC1h.itirapina[[i]][,j], by=list(section=patches.itirapina.rep), FUN=mean.present)$x
		if(j%%500==0) print(c(i,j))
	}
}
freqs.itirapina <- list()
covers.itirapina <- list()
for(i in 1:4) {
	freqs.itirapina[[i]] <- matrix(nrow=length(unique(patches.itirapina.rep)), ncol=3)
	colnames(freqs.itirapina[[i]]) <- c("Observed", "CI_0.025", "CI_0.975")
	row.names(freqs.itirapina[[i]]) <- patch.types.itirapina
	freqs.itirapina[[i]][,1] <- aggregate(data.itirapina[,i],by=list(section=patches.itirapina.rep),FUN=calc.freq)$x
	freqs.itirapina[[i]][,2:3] <- t(apply(freqs.itirapina.rand[[i]], 1, quantile, probs=c(0.025, 0.975)))
	
	covers.itirapina[[i]] <- matrix(nrow=length(unique(patches.itirapina.rep)), ncol=3)
	colnames(covers.itirapina[[i]]) <- c("Observed", "CI_0.025", "CI_0.975")
	row.names(covers.itirapina[[i]]) <- patch.types.itirapina
	covers.itirapina[[i]][,1] <- aggregate(data.itirapina[,i],by=list(section=patches.itirapina.rep),FUN=mean.present)$x
	covers.itirapina[[i]][,2:3] <- t(apply(covers.itirapina.rand[[i]], 1, quantile, probs=c(0.025, 0.975)))
}
names(freqs.itirapina) <- colnames(data.itirapina)[1:4]
names(covers.itirapina) <- colnames(data.itirapina)[1:4]
freqs.itirapina.m <- rbind(freqs.itirapina[[1]],freqs.itirapina[[2]],freqs.itirapina[[3]],freqs.itirapina[[4]])
covers.itirapina.m <- rbind(covers.itirapina[[1]],covers.itirapina[[2]],covers.itirapina[[3]],covers.itirapina[[4]])
freqs.itirapina.df <- data.frame(freqs.itirapina.m)
covers.itirapina.df <- data.frame(covers.itirapina.m)
freqs.itirapina.df$Patch_type <- row.names(freqs.itirapina.m)
covers.itirapina.df$Patch_type <- row.names(covers.itirapina.m)
freqs.itirapina.df$Variable <- rep(names(freqs.itirapina), each=length(unique(patches.itirapina.rep)))
covers.itirapina.df$Variable <- rep(names(covers.itirapina), each=length(unique(patches.itirapina.rep)))
freqs.itirapina.df <- freqs.itirapina.df[,c(5,4,1,2,3)]
covers.itirapina.df <- covers.itirapina.df[,c(5,4,1,2,3)]

freqs.ufscar1.rand <- list()
covers.ufscar1.rand <- list()
for(i in 1:4) {
	freqs.ufscar1.rand[[i]] <- matrix(nrow=length(patch.types.ufscar1), ncol=Nperm)
	covers.ufscar1.rand[[i]] <- matrix(nrow=length(patch.types.ufscar1), ncol=Nperm)	
	row.names(freqs.ufscar1.rand[[i]]) <- patch.types.ufscar1
	row.names(covers.ufscar1.rand[[i]]) <- patch.types.ufscar1
	for(j in 1:Nperm) {
		freqs.ufscar1.rand[[i]][,j] <- aggregate(rand.MC1h.ufscar1[[i]][,j], by=list(section=patches.ufscar1.rep), FUN=calc.freq)$x
		covers.ufscar1.rand[[i]][,j] <- aggregate(rand.MC1h.ufscar1[[i]][,j], by=list(section=patches.ufscar1.rep), FUN=mean.present)$x
		if(j%%500==0) print(c(i,j))
	}
}
freqs.ufscar1 <- list()
covers.ufscar1 <- list()
for(i in 1:4) {
	freqs.ufscar1[[i]] <- matrix(nrow=length(unique(patches.ufscar1.rep)), ncol=3)
	colnames(freqs.ufscar1[[i]]) <- c("Observed", "CI_0.025", "CI_0.975")
	row.names(freqs.ufscar1[[i]]) <- patch.types.ufscar1
	freqs.ufscar1[[i]][,1] <- aggregate(data.ufscar1[,i],by=list(section=patches.ufscar1.rep),FUN=calc.freq)$x
	freqs.ufscar1[[i]][,2:3] <- t(apply(freqs.ufscar1.rand[[i]], 1, quantile, probs=c(0.025, 0.975)))
	
	covers.ufscar1[[i]] <- matrix(nrow=length(unique(patches.ufscar1.rep)), ncol=3)
	colnames(covers.ufscar1[[i]]) <- c("Observed", "CI_0.025", "CI_0.975")
	row.names(covers.ufscar1[[i]]) <- patch.types.ufscar1
	covers.ufscar1[[i]][,1] <- aggregate(data.ufscar1[,i],by=list(section=patches.ufscar1.rep),FUN=mean.present)$x
	covers.ufscar1[[i]][,2:3] <- t(apply(covers.ufscar1.rand[[i]], 1, quantile, probs=c(0.025, 0.975)))
}
names(freqs.ufscar1) <- colnames(data.ufscar1)[1:4]
names(covers.ufscar1) <- colnames(data.ufscar1)[1:4]
freqs.ufscar1.m <- rbind(freqs.ufscar1[[1]],freqs.ufscar1[[2]],freqs.ufscar1[[3]],freqs.ufscar1[[4]])
covers.ufscar1.m <- rbind(covers.ufscar1[[1]],covers.ufscar1[[2]],covers.ufscar1[[3]],covers.ufscar1[[4]])
freqs.ufscar1.df <- data.frame(freqs.ufscar1.m)
covers.ufscar1.df <- data.frame(covers.ufscar1.m)
freqs.ufscar1.df$Patch_type <- row.names(freqs.ufscar1.m)
covers.ufscar1.df$Patch_type <- row.names(covers.ufscar1.m)
freqs.ufscar1.df$Variable <- rep(names(freqs.ufscar1), each=length(unique(patches.ufscar1.rep)))
covers.ufscar1.df$Variable <- rep(names(covers.ufscar1), each=length(unique(patches.ufscar1.rep)))
freqs.ufscar1.df <- freqs.ufscar1.df[,c(5,4,1,2,3)]
covers.ufscar1.df <- covers.ufscar1.df[,c(5,4,1,2,3)]

freqs.ufscar2.rand <- list()
covers.ufscar2.rand <- list()
for(i in 1:4) {
	freqs.ufscar2.rand[[i]] <- matrix(nrow=length(patch.types.ufscar2), ncol=Nperm)
	covers.ufscar2.rand[[i]] <- matrix(nrow=length(patch.types.ufscar2), ncol=Nperm)	
	row.names(freqs.ufscar2.rand[[i]]) <- patch.types.ufscar2
	row.names(covers.ufscar2.rand[[i]]) <- patch.types.ufscar2
	for(j in 1:Nperm) {
		freqs.ufscar2.rand[[i]][,j] <- aggregate(rand.MC1h.ufscar2[[i]][,j], by=list(section=patches.ufscar2.rep), FUN=calc.freq)$x
		covers.ufscar2.rand[[i]][,j] <- aggregate(rand.MC1h.ufscar2[[i]][,j], by=list(section=patches.ufscar2.rep), FUN=mean.present)$x
		if(j%%500==0) print(c(i,j))
	}
}
freqs.ufscar2 <- list()
covers.ufscar2 <- list()
for(i in 1:4) {
	freqs.ufscar2[[i]] <- matrix(nrow=length(unique(patches.ufscar2.rep)), ncol=3)
	colnames(freqs.ufscar2[[i]]) <- c("Observed", "CI_0.025", "CI_0.975")
	row.names(freqs.ufscar2[[i]]) <- patch.types.ufscar2
	freqs.ufscar2[[i]][,1] <- aggregate(data.ufscar2[,i],by=list(section=patches.ufscar2.rep),FUN=calc.freq)$x
	freqs.ufscar2[[i]][,2:3] <- t(apply(freqs.ufscar2.rand[[i]], 1, quantile, probs=c(0.025, 0.975)))
	
	covers.ufscar2[[i]] <- matrix(nrow=length(unique(patches.ufscar2.rep)), ncol=3)
	colnames(covers.ufscar2[[i]]) <- c("Observed", "CI_0.025", "CI_0.975")
	row.names(covers.ufscar2[[i]]) <- patch.types.ufscar2
	covers.ufscar2[[i]][,1] <- aggregate(data.ufscar2[,i],by=list(section=patches.ufscar2.rep),FUN=mean.present)$x
	covers.ufscar2[[i]][,2:3] <- t(apply(covers.ufscar2.rand[[i]], 1, quantile, probs=c(0.025, 0.975)))
}
names(freqs.ufscar2) <- colnames(data.ufscar2)[1:4]
names(covers.ufscar2) <- colnames(data.ufscar2)[1:4]
freqs.ufscar2.m <- rbind(freqs.ufscar2[[1]],freqs.ufscar2[[2]],freqs.ufscar2[[3]],freqs.ufscar2[[4]])
covers.ufscar2.m <- rbind(covers.ufscar2[[1]],covers.ufscar2[[2]],covers.ufscar2[[3]],covers.ufscar2[[4]])
freqs.ufscar2.df <- data.frame(freqs.ufscar2.m)
covers.ufscar2.df <- data.frame(covers.ufscar2.m)
freqs.ufscar2.df$Patch_type <- row.names(freqs.ufscar2.m)
covers.ufscar2.df$Patch_type <- row.names(covers.ufscar2.m)
freqs.ufscar2.df$Variable <- rep(names(freqs.ufscar2), each=length(unique(patches.ufscar2.rep)))
covers.ufscar2.df$Variable <- rep(names(covers.ufscar2), each=length(unique(patches.ufscar2.rep)))
freqs.ufscar2.df <- freqs.ufscar2.df[,c(5,4,1,2,3)]
covers.ufscar2.df <- covers.ufscar2.df[,c(5,4,1,2,3)]

# Save the results
setwd("/home/pavel/Profissional/Pesquisa/MyPapers-Working/EIgrasses/results")

write.table(freqs.itirapina.df, file="part1_freqs_itirapina.txt", row.names=F, quote=F)
write.table(covers.itirapina.df, file="part1_covers_itirapina.txt", row.names=F, quote=F)
write.table(freqs.ufscar1.df, file="part1_freqs_ufscar1.txt", row.names=F, quote=F)
write.table(covers.ufscar1.df, file="part1_covers_ufscar1.txt", row.names=F, quote=F)
write.table(freqs.ufscar2.df, file="part1_freqs_ufscar2.txt", row.names=F, quote=F)
write.table(covers.ufscar2.df, file="part1_covers_ufscar2.txt", row.names=F, quote=F)


###Part 2 - Is cover by elevation and by distance to firebreak edges? 

# Calculate distance to firebreaks and organize the data

# Adjust the gam models and calculate the relative importance of each variable (wAICv)

# GAM results:

test <- gam(Urochloa ~ s(DistToFirebreak, fx=F, k=5) + Vegetation, data=data2.itirapina, family=binomial)
test2 <- gam(Urochloa ~ Vegetation, data=data2.itirapina, family=binomial)

summary(test)$dev.expl
summary(test)$edf



# Simulate the null models
Nperm=5000
data.itirapina.rand <- lapply(data.itirapina[,1:4], MC1s, sections=patches.itirapina, which.loop=1000, Nperm=Nperm)
data.ufscar1.rand <- lapply(data.ufscar1[,1:4], MC1s, sections=patches.ufscar1, which.loop=1000, Nperm=Nperm)
data.ufscar2.rand <- lapply(data.ufscar2[,1:4], MC1s, sections=patches.ufscar2, which.loop=1000, Nperm=Nperm)

# Remove firebreaks from the null models

for(i in 1:4) {
	data.itirapina.rand[[i]] <- data.itirapina.rand[[i]][!is.na(data.itirapina$DistToFirebreak),]
	data.ufscar1.rand[[i]] <- data.ufscar1.rand[[i]][!is.na(data.ufscar1$DistToFirebreak),]
	data.ufscar2.rand[[i]] <- data.ufscar2.rand[[i]][!is.na(data.ufscar2$DistToFirebreak),]
}

elev.itirapina <- data2.itirapina.e$Elevation
dist.itirapina <- data2.itirapina.e$DistToFirebreak
veg.itirapina <- data2.itirapina$Vegetation
elev.ufscar1 <- data2.ufscar1.e$Elevation
dist.ufscar1 <- data2.ufscar1.e$DistToFirebreak
veg.ufscar1 <- data2.ufscar1$Vegetation
elev.ufscar2 <- data2.ufscar2.e$Elevation
dist.ufscar2 <- data2.ufscar2.e$DistToFirebreak
veg.ufscar2 <- data2.ufscar2$Vegetation
gamDev.dist.itirapina <- matrix(NA,nrow=4, ncol=Nperm)
gamDev.dist.ufscar1 <- matrix(NA,nrow=4, ncol=Nperm)
gamDev.dist.ufscar2 <- matrix(NA,nrow=4, ncol=Nperm)
gamDev.elev.itirapina <- matrix(NA,nrow=4, ncol=Nperm)
gamDev.elev.ufscar1 <- matrix(NA,nrow=4, ncol=Nperm)
gamDev.elev.ufscar2 <- matrix(NA,nrow=4, ncol=Nperm)

row.names(gamDev.dist.itirapina) <- row.names(gamDev.elev.itirapina) <- row.names(gamDev.dist.ufscar1) <- row.names(gamDev.elev.ufscar1) <- row.names(gamDev.dist.ufscar2) <- row.names(gamDev.elev.ufscar2) <- c("Urochloa", "Melinis", "Grasses", "Sedges")
for(i in 1:Nperm) {
	foo <- data.itirapina.rand$Urochloa[,i]
	bar <- gam(foo ~ s(dist.itirapina, fx=F, k=5) + veg.itirapina, family=binomial)
	foobar <- summary(bar)$dev.expl
	gamDev.dist.itirapina[1,i] <- foobar
	foo <- data.itirapina.rand$Melinis[,i]
	bar <- gam(foo ~ s(dist.itirapina, fx=F, k=5) + veg.itirapina, family=binomial)
	foobar <- summary(bar)$dev.expl
	gamDev.dist.itirapina[2,i] <- foobar
	foo <- data.itirapina.rand$Grasses[,i]
	bar <- gam(foo ~ s(dist.itirapina, fx=F, k=5) + veg.itirapina, family=binomial)
	foobar <- summary(bar)$dev.expl
	gamDev.dist.itirapina[3,i] <- foobar
	foo <- data.itirapina.rand$Sedges[,i]
	bar <- gam(foo ~ s(dist.itirapina, fx=F, k=5) + veg.itirapina, family=binomial)
	foobar <- summary(bar)$dev.expl
	gamDev.dist.itirapina[4,i] <- foobar
	foo <- data.ufscar1.rand$Urochloa[,i]
	bar <- gam(foo ~ s(dist.ufscar1, fx=F, k=5) + veg.ufscar1, family=binomial)
	foobar <- summary(bar)$dev.expl
	gamDev.dist.ufscar1[1,i] <- foobar
	foo <- data.ufscar1.rand$Melinis[,i]
	bar <- gam(foo ~ s(dist.ufscar1, fx=F, k=5) + veg.ufscar1, family=binomial)
	foobar <- summary(bar)$dev.expl
	gamDev.dist.ufscar1[2,i] <- foobar
	foo <- data.ufscar1.rand$Grasses[,i]
	bar <- gam(foo ~ s(dist.ufscar1, fx=F, k=5) + veg.ufscar1, family=binomial)
	foobar <- summary(bar)$dev.expl
	gamDev.dist.ufscar1[3,i] <- foobar
	foo <- data.ufscar1.rand$Sedges[,i]
	bar <- gam(foo ~ s(dist.ufscar1, fx=F, k=5) + veg.ufscar1, family=binomial)
	foobar <- summary(bar)$dev.expl
	gamDev.dist.ufscar1[4,i] <- foobar
	foo <- data.ufscar2.rand$Urochloa[,i]
	bar <- gam(foo ~ s(dist.ufscar2, fx=F, k=5) + veg.ufscar2, family=binomial)
	foobar <- summary(bar)$dev.expl
	gamDev.dist.ufscar2[1,i] <- foobar
	foo <- data.ufscar2.rand$Melinis[,i]
	bar <- gam(foo ~ s(dist.ufscar2, fx=F, k=5) + veg.ufscar2, family=binomial)
	foobar <- summary(bar)$dev.expl
	gamDev.dist.ufscar2[2,i] <- foobar
	foo <- data.ufscar2.rand$Grasses[,i]
	bar <- gam(foo ~ s(dist.ufscar2, fx=F, k=5) + veg.ufscar2, family=binomial)
	foobar <- summary(bar)$dev.expl
	gamDev.dist.ufscar2[3,i] <- foobar
	foo <- data.ufscar2.rand$Sedges[,i]
	bar <- gam(foo ~ s(dist.ufscar2, fx=F, k=5) + veg.ufscar2, family=binomial)
	foobar <- summary(bar)$dev.expl
	gamDev.dist.ufscar2[4,i] <- foobar

	foo <- data.itirapina.rand$Urochloa[,i]
	bar <- gam(foo ~ s(elev.itirapina, fx=F, k=5) + veg.itirapina, family=binomial)
	foobar <- summary(bar)$dev.expl
	gamDev.elev.itirapina[1,i] <- foobar
	foo <- data.itirapina.rand$Melinis[,i]
	bar <- gam(foo ~ s(elev.itirapina, fx=F, k=5) + veg.itirapina, family=binomial)
	foobar <- summary(bar)$dev.expl
	gamDev.elev.itirapina[2,i] <- foobar
	foo <- data.itirapina.rand$Grasses[,i]
	bar <- gam(foo ~ s(elev.itirapina, fx=F, k=5) + veg.itirapina, family=binomial)
	foobar <- summary(bar)$dev.expl
	gamDev.elev.itirapina[3,i] <- foobar
	foo <- data.itirapina.rand$Sedges[,i]
	bar <- gam(foo ~ s(elev.itirapina, fx=F, k=5) + veg.itirapina, family=binomial)
	foobar <- summary(bar)$dev.expl
	gamDev.elev.itirapina[4,i] <- foobar
	foo <- data.ufscar1.rand$Urochloa[,i]
	bar <- gam(foo ~ s(elev.ufscar1, fx=F, k=5) + veg.ufscar1, family=binomial)
	foobar <- summary(bar)$dev.expl
	gamDev.elev.ufscar1[1,i] <- foobar
	foo <- data.ufscar1.rand$Melinis[,i]
	bar <- gam(foo ~ s(elev.ufscar1, fx=F, k=5) + veg.ufscar1, family=binomial)
	foobar <- summary(bar)$dev.expl
	gamDev.elev.ufscar1[2,i] <- foobar
	foo <- data.ufscar1.rand$Grasses[,i]
	bar <- gam(foo ~ s(elev.ufscar1, fx=F, k=5) + veg.ufscar1, family=binomial)
	foobar <- summary(bar)$dev.expl
	gamDev.elev.ufscar1[3,i] <- foobar
	foo <- data.ufscar1.rand$Sedges[,i]
	bar <- gam(foo ~ s(elev.ufscar1, fx=F, k=5) + veg.ufscar1, family=binomial)
	foobar <- summary(bar)$dev.expl
	gamDev.elev.ufscar1[4,i] <- foobar
	foo <- data.ufscar2.rand$Urochloa[,i]
	bar <- gam(foo ~ s(elev.ufscar2, fx=F, k=5) + veg.ufscar2, family=binomial)
	foobar <- summary(bar)$dev.expl
	gamDev.elev.ufscar2[1,i] <- foobar
	foo <- data.ufscar2.rand$Melinis[,i]
	bar <- gam(foo ~ s(elev.ufscar2, fx=F, k=5) + veg.ufscar2, family=binomial)
	foobar <- summary(bar)$dev.expl
	gamDev.elev.ufscar2[2,i] <- foobar
	foo <- data.ufscar2.rand$Grasses[,i]
	bar <- gam(foo ~ s(elev.ufscar2, fx=F, k=5) + veg.ufscar2, family=binomial)
	foobar <- summary(bar)$dev.expl
	gamDev.elev.ufscar2[3,i] <- foobar
	foo <- data.ufscar2.rand$Sedges[,i]
	bar <- gam(foo ~ s(elev.ufscar2, fx=F, k=5) + veg.ufscar2, family=binomial)
	foobar <- summary(bar)$dev.expl
	gamDev.elev.ufscar2[4,i] <- foobar
	print(i)
}


# Make histograms
setwd("/home/pavel/Profissional/Pesquisa/MyPapers-Working/EIgrasses/figs_gam_hist")
png(filename="gam_hist_itirapina.png", width=20, height=30, unit="cm", res=300)
par(mfrow=c(4,2), mar=c(3,3,2,2))
hist(gamDev.dist.itirapina[1,], main="Urochloa_dist")
abline(v=gamDev.dist.itirapina[1,1], lwd=2, col="red")
hist(gamDev.elev.itirapina[1,], main="Urochloa_elev")
abline(v=gamDev.elev.itirapina[1,1], lwd=2, col="red")
hist(gamDev.dist.itirapina[2,], main="Melinis_dist")
abline(v=gamDev.dist.itirapina[2,1], lwd=2, col="red")
hist(gamDev.elev.itirapina[2,], main="Melinis_elev")
abline(v=gamDev.elev.itirapina[2,1], lwd=2, col="red")
hist(gamDev.dist.itirapina[3,], main="Grasses_dist")
abline(v=gamDev.dist.itirapina[3,1], lwd=2, col="red")
hist(gamDev.elev.itirapina[3,], main="Grasses_elev")
abline(v=gamDev.elev.itirapina[3,1], lwd=2, col="red")
hist(gamDev.dist.itirapina[4,], main="Sedges_dist")
abline(v=gamDev.dist.itirapina[4,1], lwd=2, col="red")
hist(gamDev.elev.itirapina[4,], main="Sedges_elev")
abline(v=gamDev.elev.itirapina[4,1], lwd=2, col="red")
dev.off()

png(filename="gam_hist_ufscar1.png", width=20, height=30, unit="cm", res=300)
par(mfrow=c(4,2), mar=c(3,3,2,2))
hist(gamDev.dist.ufscar1[1,], main="Urochloa_dist")
abline(v=gamDev.dist.ufscar1[1,1], lwd=2, col="red")
hist(gamDev.elev.ufscar1[1,], main="Urochloa_elev")
abline(v=gamDev.elev.ufscar1[1,1], lwd=2, col="red")
hist(gamDev.dist.ufscar1[2,], main="Melinis_dist")
abline(v=gamDev.dist.ufscar1[2,1], lwd=2, col="red")
hist(gamDev.elev.ufscar1[2,], main="Melinis_elev")
abline(v=gamDev.elev.ufscar1[2,1], lwd=2, col="red")
hist(gamDev.dist.ufscar1[3,], main="Grasses_dist")
abline(v=gamDev.dist.ufscar1[3,1], lwd=2, col="red")
hist(gamDev.elev.ufscar1[3,], main="Grasses_elev")
abline(v=gamDev.elev.ufscar1[3,1], lwd=2, col="red")
hist(gamDev.dist.ufscar1[4,], main="Sedges_dist")
abline(v=gamDev.dist.ufscar1[4,1], lwd=2, col="red")
hist(gamDev.elev.ufscar1[4,], main="Sedges_elev")
abline(v=gamDev.elev.ufscar1[4,1], lwd=2, col="red")
dev.off()

png(filename="gam_hist_ufscar2.png", width=20, height=30, unit="cm", res=300)
par(mfrow=c(4,2), mar=c(3,3,2,2))
hist(gamDev.dist.ufscar2[1,], main="Urochloa_dist")
abline(v=gamDev.dist.ufscar2[1,1], lwd=2, col="red")
hist(gamDev.elev.ufscar2[1,], main="Urochloa_elev")
abline(v=gamDev.elev.ufscar2[1,1], lwd=2, col="red")
hist(gamDev.dist.ufscar2[2,], main="Melinis_dist")
abline(v=gamDev.dist.ufscar2[2,1], lwd=2, col="red")
hist(gamDev.elev.ufscar2[2,], main="Melinis_elev")
abline(v=gamDev.elev.ufscar2[2,1], lwd=2, col="red")
hist(gamDev.dist.ufscar2[3,], main="Grasses_dist")
abline(v=gamDev.dist.ufscar2[3,1], lwd=2, col="red")
hist(gamDev.elev.ufscar2[3,], main="Grasses_elev")
abline(v=gamDev.elev.ufscar2[3,1], lwd=2, col="red")
hist(gamDev.dist.ufscar2[4,], main="Sedges_dist")
abline(v=gamDev.dist.ufscar2[4,1], lwd=2, col="red")
hist(gamDev.elev.ufscar2[4,], main="Sedges_elev")
abline(v=gamDev.elev.ufscar2[4,1], lwd=2, col="red")
dev.off()

# Calculate significance
signif.dist.itirapina <- apply(gamDev.dist.itirapina,1,signif.perm.above)
signif.elev.itirapina <- apply(gamDev.elev.itirapina,1,signif.perm.above)
signif.dist.ufscar1 <- apply(gamDev.dist.ufscar1,1,signif.perm.above)
signif.elev.ufscar1 <- apply(gamDev.elev.ufscar1,1,signif.perm.above)
signif.dist.ufscar2 <- apply(gamDev.dist.ufscar2,1,signif.perm.above)
signif.elev.ufscar2 <- apply(gamDev.elev.ufscar2,1,signif.perm.above)
signif.all <- rbind(signif.dist.itirapina, signif.elev.itirapina, signif.dist.ufscar1, signif.elev.ufscar1, signif.dist.ufscar2, signif.elev.ufscar2)
setwd("/home/pavel/Profissional/Pesquisa/MyPapers-Working/EIgrasses/results")
write.table(signif.all, file="part2_gamSignif.txt", row.names=T, quote=F, sep="\t")


# Adjust the models and make the figures
# Itirapina: Grasses - elev
# Ufscar1: Urochloa - elev, Sedges - dist
# Ufscar2: Grasses - dist, elev (marginally significant)
gam.itirapina.Grasses.elev <- gam(Grasses ~ s(Elevation, fx=F, k=5) + Vegetation, data=data2.itirapina, family=binomial)
new.itirapina.elev.ecotone <- list(Elevation=seq(min(subset(data2.itirapina,Vegetation=="ecotone",Elevation)), max(subset(data2.itirapina,Vegetation=="ecotone",Elevation)), length.out=50), Vegetation=rep("ecotone",50))
new.itirapina.elev.forest <- list(Elevation=seq(min(subset(data2.itirapina,Vegetation=="forest",Elevation)), max(subset(data2.itirapina,Vegetation=="forest",Elevation)), length.out=50), Vegetation=rep("forest",50))
new.itirapina.elev.grassland <- list(Elevation=seq(min(subset(data2.itirapina,Vegetation=="grassland",Elevation)), max(subset(data2.itirapina,Vegetation=="grassland",Elevation)), length.out=50), Vegetation=rep("grassland",50))
new.itirapina.elev.grassland_invaded <- list(Elevation=seq(min(subset(data2.itirapina,Vegetation=="grassland_invaded",Elevation)), max(subset(data2.itirapina,Vegetation=="grassland_invaded",Elevation)), length.out=50), Vegetation=rep("grassland_invaded",50))
pred.itirapina.Grasses.elev.ecotone <- predict(gam.itirapina.Grasses.elev, newdata=new.itirapina.elev.ecotone, type="response")
pred.itirapina.Grasses.elev.forest <- predict(gam.itirapina.Grasses.elev, newdata=new.itirapina.elev.forest, type="response")
pred.itirapina.Grasses.elev.grassland <- predict(gam.itirapina.Grasses.elev, newdata=new.itirapina.elev.grassland, type="response")
pred.itirapina.Grasses.elev.grassland_invaded <- predict(gam.itirapina.Grasses.elev, newdata=new.itirapina.elev.grassland_invaded, type="response")

gam.ufscar1.Urochloa.elev <- gam(Urochloa ~ s(Elevation, fx=F, k=5) + Vegetation, data=data2.ufscar1, family=binomial)
new.ufscar1.elev.cerrado <- list(Elevation=seq(min(subset(data2.ufscar1,Vegetation=="cerrado",Elevation)), max(subset(data2.ufscar1,Vegetation=="cerrado",Elevation)), length.out=50), Vegetation=rep("cerrado",50))
new.ufscar1.elev.cerrado_regenerating <- list(Elevation=seq(min(subset(data2.ufscar1,Vegetation=="cerrado_regenerating",Elevation)), max(subset(data2.ufscar1,Vegetation=="cerrado_regenerating",Elevation)), length.out=50), Vegetation=rep("cerrado_regenerating",50))
new.ufscar1.elev.grassland_invaded <- list(Elevation=seq(min(subset(data2.ufscar1,Vegetation=="grassland_invaded",Elevation)), max(subset(data2.ufscar1,Vegetation=="grassland_invaded",Elevation)), length.out=50), Vegetation=rep("grassland_invaded",50))
pred.ufscar1.Urochloa.elev.cerrado <- predict(gam.ufscar1.Urochloa.elev, newdata=new.ufscar1.elev.cerrado, type="response")
pred.ufscar1.Urochloa.elev.cerrado_regenerating <- predict(gam.ufscar1.Urochloa.elev, newdata=new.ufscar1.elev.cerrado_regenerating, type="response")
pred.ufscar1.Urochloa.elev.grassland_invaded <- predict(gam.ufscar1.Urochloa.elev, newdata=new.ufscar1.elev.grassland_invaded, type="response")


gam.ufscar1.Sedges.dist <- gam(Sedges ~ s(DistToFirebreak, fx=F, k=5) + Vegetation, data=data2.ufscar1, family=binomial)
new.ufscar1.dist.cerrado <- list(DistToFirebreak=seq(min(subset(data2.ufscar1,Vegetation=="cerrado",DistToFirebreak)), max(subset(data2.ufscar1,Vegetation=="cerrado",DistToFirebreak)), length.out=50), Vegetation=rep("cerrado",50))
new.ufscar1.dist.cerrado_regenerating <- list(DistToFirebreak=seq(min(subset(data2.ufscar1,Vegetation=="cerrado_regenerating",DistToFirebreak)), max(subset(data2.ufscar1,Vegetation=="cerrado_regenerating",DistToFirebreak)), length.out=50), Vegetation=rep("cerrado_regenerating",50))
new.ufscar1.dist.grassland_invaded <- list(DistToFirebreak=seq(min(subset(data2.ufscar1,Vegetation=="grassland_invaded",DistToFirebreak)), max(subset(data2.ufscar1,Vegetation=="grassland_invaded",DistToFirebreak)), length.out=50), Vegetation=rep("grassland_invaded",50))
pred.ufscar1.Sedges.dist.cerrado <- predict(gam.ufscar1.Sedges.dist, newdata=new.ufscar1.dist.cerrado, type="response")
pred.ufscar1.Sedges.dist.cerrado_regenerating <- predict(gam.ufscar1.Sedges.dist, newdata=new.ufscar1.dist.cerrado_regenerating, type="response")
pred.ufscar1.Sedges.dist.grassland_invaded <- predict(gam.ufscar1.Sedges.dist, newdata=new.ufscar1.dist.grassland_invaded, type="response")



gam.ufscar2.Grasses.elev <- gam(Grasses ~ s(Elevation, fx=F, k=5) + Vegetation, data=data2.ufscar2, family=binomial)
new.ufscar2.elev.cerrado <- list(Elevation=seq(min(subset(data2.ufscar2,Vegetation=="cerrado",Elevation)), max(subset(data2.ufscar2,Vegetation=="cerrado",Elevation)), length.out=50), Vegetation=rep("cerrado",50))
new.ufscar2.elev.cerrado_regenerating <- list(Elevation=seq(min(subset(data2.ufscar2,Vegetation=="cerrado_regenerating",Elevation)), max(subset(data2.ufscar2,Vegetation=="cerrado_regenerating",Elevation)), length.out=50), Vegetation=rep("cerrado_regenerating",50))
new.ufscar2.elev.ecotone <- list(Elevation=seq(min(subset(data2.ufscar2,Vegetation=="ecotone",Elevation)), max(subset(data2.ufscar2,Vegetation=="ecotone",Elevation)), length.out=50), Vegetation=rep("ecotone",50))
new.ufscar2.elev.forest <- list(Elevation=seq(min(subset(data2.ufscar2,Vegetation=="forest",Elevation)), max(subset(data2.ufscar2,Vegetation=="forest",Elevation)), length.out=50), Vegetation=rep("forest",50))
pred.ufscar2.Grasses.elev.cerrado <- predict(gam.ufscar2.Grasses.elev, newdata=new.ufscar2.elev.cerrado, type="response")
pred.ufscar2.Grasses.elev.cerrado_regenerating <- predict(gam.ufscar2.Grasses.elev, newdata=new.ufscar2.elev.cerrado_regenerating, type="response")
pred.ufscar2.Grasses.elev.ecotone <- predict(gam.ufscar2.Grasses.elev, newdata=new.ufscar2.elev.ecotone, type="response")
pred.ufscar2.Grasses.elev.forest <- predict(gam.ufscar2.Grasses.elev, newdata=new.ufscar2.elev.forest, type="response")

gam.ufscar2.Grasses.dist <- gam(Grasses ~ s(DistToFirebreak, fx=F, k=5) + Vegetation, data=data2.ufscar2, family=binomial)
new.ufscar2.dist.cerrado <- list(DistToFirebreak=seq(min(subset(data2.ufscar2,Vegetation=="cerrado",DistToFirebreak)), max(subset(data2.ufscar2,Vegetation=="cerrado",DistToFirebreak)), length.out=50), Vegetation=rep("cerrado",50))
new.ufscar2.dist.cerrado_regenerating <- list(DistToFirebreak=seq(min(subset(data2.ufscar2,Vegetation=="cerrado_regenerating",DistToFirebreak)), max(subset(data2.ufscar2,Vegetation=="cerrado_regenerating",DistToFirebreak)), length.out=50), Vegetation=rep("cerrado_regenerating",50))
new.ufscar2.dist.ecotone <- list(DistToFirebreak=seq(min(subset(data2.ufscar2,Vegetation=="ecotone",DistToFirebreak)), max(subset(data2.ufscar2,Vegetation=="ecotone",DistToFirebreak)), length.out=50), Vegetation=rep("ecotone",50))
new.ufscar2.dist.forest <- list(DistToFirebreak=seq(min(subset(data2.ufscar2,Vegetation=="forest",DistToFirebreak)), max(subset(data2.ufscar2,Vegetation=="forest",DistToFirebreak)), length.out=50), Vegetation=rep("forest",50))
pred.ufscar2.Grasses.dist.cerrado <- predict(gam.ufscar2.Grasses.dist, newdata=new.ufscar2.dist.cerrado, type="response")
pred.ufscar2.Grasses.dist.cerrado_regenerating <- predict(gam.ufscar2.Grasses.dist, newdata=new.ufscar2.dist.cerrado_regenerating, type="response")
pred.ufscar2.Grasses.dist.ecotone <- predict(gam.ufscar2.Grasses.dist, newdata=new.ufscar2.dist.ecotone, type="response")
pred.ufscar2.Grasses.dist.forest <- predict(gam.ufscar2.Grasses.dist, newdata=new.ufscar2.dist.forest, type="response")

setwd("/home/pavel/Profissional/Pesquisa/MyPapers-Revising/EIgrasses/Submitted4")

png(filename="EIgrasses_fig4_GAM_final_subm3.png", res=300, unit="cm", height=20, width=20)

par(mfrow=c(3,2), mar=c(2,2,3,2), oma=c(4,4,1,1))

plot(Grasses ~ Elevation, data=data2.itirapina, col="gray50", ylim=c(0,1), main=expression("a. Native grasses - I1"), yaxt="n", cex.main=1.5)
lines(pred.itirapina.Grasses.elev.ecotone ~ new.itirapina.elev.ecotone$Elevation, lwd=2)
lines(pred.itirapina.Grasses.elev.forest ~ new.itirapina.elev.forest$Elevation, lwd=2)
lines(pred.itirapina.Grasses.elev.grassland ~ new.itirapina.elev.grassland$Elevation, lwd=2)
lines(pred.itirapina.Grasses.elev.grassland_invaded ~ new.itirapina.elev.grassland_invaded$Elevation, lwd=2)
axis(side=2, at=c(0,0.25, 0.50, 0.75, 1), labels=c(0,25,50,75,100), las=1)

plot(0, xaxt="n", yaxt="n", bty="n", type="n", xlab="", ylab="")

plot(Urochloa ~ Elevation, data=data2.ufscar1, col="gray50", ylim=c(0,1), main=expression(paste("b. ", italic("Urochloa decumbens"), " - S1")), yaxt="n", cex.main=1.5)
lines(pred.ufscar1.Urochloa.elev.cerrado ~ new.ufscar1.elev.cerrado$Elevation, lwd=2)
lines(pred.ufscar1.Urochloa.elev.cerrado_regenerating ~ new.ufscar1.elev.cerrado_regenerating$Elevation, lwd=2)
lines(pred.ufscar1.Urochloa.elev.grassland_invaded ~ new.ufscar1.elev.grassland_invaded$Elevation, lwd=2)
axis(side=2, at=c(0,0.25, 0.50, 0.75, 1), labels=c(0,25,50,75,100), las=1)

plot(Sedges ~ DistToFirebreak, data=data2.ufscar1, col="gray50", ylim=c(0,1), main=expression("d. Native sedges - S1"), yaxt="n", cex.main=1.5)
lines(pred.ufscar1.Sedges.dist.cerrado ~ new.ufscar1.dist.cerrado$DistToFirebreak, lwd=2)
lines(pred.ufscar1.Sedges.dist.cerrado_regenerating ~ new.ufscar1.dist.cerrado_regenerating$DistToFirebreak, lwd=2)
lines(pred.ufscar1.Sedges.dist.grassland_invaded ~ new.ufscar1.dist.grassland_invaded$DistToFirebreak, lwd=2)
axis(side=2, at=c(0,0.25, 0.50, 0.75, 1), labels=c(0,25,50,75,100), las=1)

plot(Grasses ~ Elevation, data=data2.ufscar2, col="gray50", ylim=c(0,1), main=expression("c. Native grasses - S2"), yaxt="n", cex.main=1.5)
lines(pred.ufscar2.Grasses.elev.cerrado ~ new.ufscar2.elev.cerrado$Elevation, lwd=2)
lines(pred.ufscar2.Grasses.elev.cerrado_regenerating ~ new.ufscar2.elev.cerrado_regenerating$Elevation, lwd=2)
lines(pred.ufscar2.Grasses.elev.ecotone ~ new.ufscar2.elev.ecotone$Elevation, lwd=2)
lines(pred.ufscar2.Grasses.elev.forest ~ new.ufscar2.elev.forest$Elevation, lwd=2)
axis(side=2, at=c(0,0.25, 0.50, 0.75, 1), labels=c(0,25,50,75,100), las=1)

plot(Grasses ~ DistToFirebreak, data=data2.ufscar2, col="gray50", ylim=c(0,1), main=expression("e. Native grasses - S2"), yaxt="n", cex.main=1.5)
lines(pred.ufscar2.Grasses.dist.cerrado ~ new.ufscar2.dist.cerrado$DistToFirebreak, lwd=2)
lines(pred.ufscar2.Grasses.dist.cerrado_regenerating ~ new.ufscar2.dist.cerrado_regenerating$DistToFirebreak, lwd=2)
lines(pred.ufscar2.Grasses.dist.ecotone ~ new.ufscar2.dist.ecotone$DistToFirebreak, lwd=2)
lines(pred.ufscar2.Grasses.dist.forest ~ new.ufscar2.dist.forest$DistToFirebreak, lwd=2)
axis(side=2, at=c(0,0.25, 0.50, 0.75, 1), labels=c(0,25,50,75,100), las=1)

mtext(text="Cover (%)", side=2, outer=T, line=1.5)
mtext(text=c("Elevation (m a. s. l.)", "Distance to firebreak (m)"), at=c(0.25, 0.75), side=1, outer=T, line=1.5)

dev.off()

# Get the other useful information...
gam.itirapina.Urochloa.elev <- gam(Urochloa ~ s(Elevation, fx=F, k=5) + Vegetation, data=data2.itirapina, family=binomial)
gam.itirapina.Urochloa.dist <- gam(Urochloa ~ s(DistToFirebreak, fx=F, k=5) + Vegetation, data=data2.itirapina, family=binomial)
gam.itirapina.Melinis.elev <- gam(Melinis ~ s(Elevation, fx=F, k=5) + Vegetation, data=data2.itirapina, family=binomial)
gam.itirapina.Melinis.dist <- gam(Melinis ~ s(DistToFirebreak, fx=F, k=5) + Vegetation, data=data2.itirapina, family=binomial)
gam.itirapina.Grasses.elev <- gam(Grasses ~ s(Elevation, fx=F, k=5) + Vegetation, data=data2.itirapina, family=binomial)
gam.itirapina.Grasses.dist <- gam(Grasses ~ s(DistToFirebreak, fx=F, k=5) + Vegetation, data=data2.itirapina, family=binomial)
gam.itirapina.Sedges.elev <- gam(Sedges ~ s(Elevation, fx=F, k=5) + Vegetation, data=data2.itirapina, family=binomial)
gam.itirapina.Sedges.dist <- gam(Sedges ~ s(DistToFirebreak, fx=F, k=5) + Vegetation, data=data2.itirapina, family=binomial)

gam.ufscar1.Urochloa.elev <- gam(Urochloa ~ s(Elevation, fx=F, k=5) + Vegetation, data=data2.ufscar1, family=binomial)
gam.ufscar1.Urochloa.dist <- gam(Urochloa ~ s(DistToFirebreak, fx=F, k=5) + Vegetation, data=data2.ufscar1, family=binomial)
gam.ufscar1.Melinis.elev <- gam(Melinis ~ s(Elevation, fx=F, k=5) + Vegetation, data=data2.ufscar1, family=binomial)
gam.ufscar1.Melinis.dist <- gam(Melinis ~ s(DistToFirebreak, fx=F, k=5) + Vegetation, data=data2.ufscar1, family=binomial)
gam.ufscar1.Grasses.elev <- gam(Grasses ~ s(Elevation, fx=F, k=5) + Vegetation, data=data2.ufscar1, family=binomial)
gam.ufscar1.Grasses.dist <- gam(Grasses ~ s(DistToFirebreak, fx=F, k=5) + Vegetation, data=data2.ufscar1, family=binomial)
gam.ufscar1.Sedges.elev <- gam(Sedges ~ s(Elevation, fx=F, k=5) + Vegetation, data=data2.ufscar1, family=binomial)
gam.ufscar1.Sedges.dist <- gam(Sedges ~ s(DistToFirebreak, fx=F, k=5) + Vegetation, data=data2.ufscar1, family=binomial)

gam.ufscar2.Urochloa.elev <- gam(Urochloa ~ s(Elevation, fx=F, k=5) + Vegetation, data=data2.ufscar2, family=binomial)
gam.ufscar2.Urochloa.dist <- gam(Urochloa ~ s(DistToFirebreak, fx=F, k=5) + Vegetation, data=data2.ufscar2, family=binomial)
gam.ufscar2.Melinis.elev <- gam(Melinis ~ s(Elevation, fx=F, k=5) + Vegetation, data=data2.ufscar2, family=binomial)
gam.ufscar2.Melinis.dist <- gam(Melinis ~ s(DistToFirebreak, fx=F, k=5) + Vegetation, data=data2.ufscar2, family=binomial)
gam.ufscar2.Grasses.elev <- gam(Grasses ~ s(Elevation, fx=F, k=5) + Vegetation, data=data2.ufscar2, family=binomial)
gam.ufscar2.Grasses.dist <- gam(Grasses ~ s(DistToFirebreak, fx=F, k=5) + Vegetation, data=data2.ufscar2, family=binomial)
gam.ufscar2.Sedges.elev <- gam(Sedges ~ s(Elevation, fx=F, k=5) + Vegetation, data=data2.ufscar2, family=binomial)
gam.ufscar2.Sedges.dist <- gam(Sedges ~ s(DistToFirebreak, fx=F, k=5) + Vegetation, data=data2.ufscar2, family=binomial)

summary(gam.itirapina.Urochloa.elev)
summary(gam.itirapina.Urochloa.dist)
summary(gam.itirapina.Melinis.elev)
summary(gam.itirapina.Melinis.dist)
summary(gam.itirapina.Grasses.elev)
summary(gam.itirapina.Grasses.dist)
summary(gam.itirapina.Sedges.elev)
summary(gam.itirapina.Sedges.dist)
summary(gam.ufscar1.Urochloa.elev)
summary(gam.ufscar1.Urochloa.dist)
summary(gam.ufscar1.Melinis.elev)
summary(gam.ufscar1.Melinis.dist)
summary(gam.ufscar1.Grasses.elev)
summary(gam.ufscar1.Grasses.dist)
summary(gam.ufscar1.Sedges.elev)
summary(gam.ufscar1.Sedges.dist)
summary(gam.ufscar2.Urochloa.elev)
summary(gam.ufscar2.Urochloa.dist)
summary(gam.ufscar2.Melinis.elev)
summary(gam.ufscar2.Melinis.dist)
summary(gam.ufscar2.Grasses.elev)
summary(gam.ufscar2.Grasses.dist)
summary(gam.ufscar2.Sedges.elev)
summary(gam.ufscar2.Sedges.dist)

### Part 3 - What are the main scales of spatial pattern?

# Simulate the null models...

rand.itirapina <- lapply(data.itirapina[,1:4], MC1s, sections=patches.itirapina)
rand.ufscar1 <- lapply(data.ufscar1[,1:4], MC1s, sections=patches.ufscar1)
rand.ufscar2 <- lapply(data.ufscar2[,1:4], MC1s, sections=patches.ufscar2)


# Apply the wavelet :-)

wav.itirapina <- lapply(rand.itirapina, wavCWTvarCIs, significance=T, sections=patches.itirapina, make.plot=F, zero.padding=T, remove.COI=T, test=F, scale.max=75, scale.var=T, pos.var=F, CI=T, CIquant=c(0.95), effect.size=F,  keep.all=F, print.loop=T, which.loop=100, save.cwt=F, file.title="my_beautiful_wavelet", wav.template="gaussian2", scale.per.section=T)

wav.ufscar1 <- lapply(rand.ufscar1, wavCWTvarCIs, significance=T, sections=patches.ufscar1, make.plot=F, zero.padding=T, remove.COI=T, test=F, scale.max=75, scale.var=T, pos.var=F, CI=T, CIquant=c(0.95), effect.size=F,  keep.all=F, print.loop=T, which.loop=100, save.cwt=F, file.title="my_beautiful_wavelet", wav.template="gaussian2", scale.per.section=T)

wav.ufscar2 <- lapply(rand.ufscar2, wavCWTvarCIs, significance=T, sections=patches.ufscar2, make.plot=F, zero.padding=T, remove.COI=T, test=F, scale.max=75, scale.var=T, pos.var=F, CI=T, CIquant=c(0.95), effect.size=F,  keep.all=F, print.loop=T, which.loop=100, save.cwt=F, file.title="my_beautiful_wavelet", wav.template="gaussian2", scale.per.section=T)


save.image(file="wav_univar.RData")

# Organize as something prettier

setwd("/home/pavel/Profissional/Pesquisa/MyPapers-Working/EIgrasses/results")

names.itirapina <- numeric(15)
names.itirapina[1:3] <- c("Scales", "Total_Obs", "Total_CI")
names.itirapina[c(4,6,8,10,12,14)] <- paste(names(wav.itirapina$Urochloa$ScaleSection_Obs),"Obs",sep="_")
names.itirapina[c(5,7,9,11,13,15)] <- paste(names(wav.itirapina$Urochloa$ScaleSection_CI),"CI",sep="_")

for(i in 1:length(wav.itirapina$Urochloa$ScaleSection_CI)) {
  wav.itirapina$Urochloa$ScaleSection_CI[[i]] <- t(wav.itirapina$Urochloa$ScaleSection_CI[[i]])
}
wav.itirapina.Urochloa <- matrix(NA,ncol=15, nrow=50)
colnames(wav.itirapina.Urochloa) <- names.itirapina
wav.itirapina.Urochloa <- as.data.frame(wav.itirapina.Urochloa)
wav.itirapina.Urochloa$Scales <- wav.itirapina$Urochloa$Scales
wav.itirapina.Urochloa$Total_Obs <- wav.itirapina$Urochloa$ScaleVar_Obs
wav.itirapina.Urochloa$Total_CI <- wav.itirapina$Urochloa$ScaleVar_CI[1,]
wav.itirapina.Urochloa[c(4,6,8,10,12,14)] <- as.data.frame(wav.itirapina$Urochloa$ScaleSection_Obs)
wav.itirapina.Urochloa[c(5,7,9,11,13,15)] <- as.data.frame(wav.itirapina$Urochloa$ScaleSection_CI)
write.table(wav.itirapina.Urochloa, file="part3_wavScales_itirapina_Urochloa.txt", quote=F, row.names=F, sep="\t")

for(i in 1:length(wav.itirapina$Melinis$ScaleSection_CI)) {
  wav.itirapina$Melinis$ScaleSection_CI[[i]] <- t(wav.itirapina$Melinis$ScaleSection_CI[[i]])
}
wav.itirapina.Melinis <- matrix(NA,ncol=15, nrow=50)
colnames(wav.itirapina.Melinis) <- names.itirapina
wav.itirapina.Melinis <- as.data.frame(wav.itirapina.Melinis)
wav.itirapina.Melinis$Scales <- wav.itirapina$Melinis$Scales
wav.itirapina.Melinis$Total_Obs <- wav.itirapina$Melinis$ScaleVar_Obs
wav.itirapina.Melinis$Total_CI <- wav.itirapina$Melinis$ScaleVar_CI[1,]
wav.itirapina.Melinis[c(4,6,8,10,12,14)] <- as.data.frame(wav.itirapina$Melinis$ScaleSection_Obs)
wav.itirapina.Melinis[c(5,7,9,11,13,15)] <- as.data.frame(wav.itirapina$Melinis$ScaleSection_CI)
write.table(wav.itirapina.Melinis, file="part3_wavScales_itirapina_Melinis.txt", quote=F, row.names=F, sep="\t")

for(i in 1:length(wav.itirapina$Grasses$ScaleSection_CI)) {
  wav.itirapina$Grasses$ScaleSection_CI[[i]] <- t(wav.itirapina$Grasses$ScaleSection_CI[[i]])
}
wav.itirapina.Grasses <- matrix(NA,ncol=15, nrow=50)
colnames(wav.itirapina.Grasses) <- names.itirapina
wav.itirapina.Grasses <- as.data.frame(wav.itirapina.Grasses)
wav.itirapina.Grasses$Scales <- wav.itirapina$Grasses$Scales
wav.itirapina.Grasses$Total_Obs <- wav.itirapina$Grasses$ScaleVar_Obs
wav.itirapina.Grasses$Total_CI <- wav.itirapina$Grasses$ScaleVar_CI[1,]
wav.itirapina.Grasses[c(4,6,8,10,12,14)] <- as.data.frame(wav.itirapina$Grasses$ScaleSection_Obs)
wav.itirapina.Grasses[c(5,7,9,11,13,15)] <- as.data.frame(wav.itirapina$Grasses$ScaleSection_CI)
write.table(wav.itirapina.Grasses, file="part3_wavScales_itirapina_Grasses.txt", quote=F, row.names=F, sep="\t")

for(i in 1:length(wav.itirapina$Sedges$ScaleSection_CI)) {
  wav.itirapina$Sedges$ScaleSection_CI[[i]] <- t(wav.itirapina$Sedges$ScaleSection_CI[[i]])
}
wav.itirapina.Sedges <- matrix(NA,ncol=15, nrow=50)
colnames(wav.itirapina.Sedges) <- names.itirapina
wav.itirapina.Sedges <- as.data.frame(wav.itirapina.Sedges)
wav.itirapina.Sedges$Scales <- wav.itirapina$Sedges$Scales
wav.itirapina.Sedges$Total_Obs <- wav.itirapina$Sedges$ScaleVar_Obs
wav.itirapina.Sedges$Total_CI <- wav.itirapina$Sedges$ScaleVar_CI[1,]
wav.itirapina.Sedges[c(4,6,8,10,12,14)] <- as.data.frame(wav.itirapina$Sedges$ScaleSection_Obs)
wav.itirapina.Sedges[c(5,7,9,11,13,15)] <- as.data.frame(wav.itirapina$Sedges$ScaleSection_CI)
write.table(wav.itirapina.Sedges, file="part3_wavScales_itirapina_Sedges.txt", quote=F, row.names=F, sep="\t")



names.ufscar1 <- numeric(11)
names.ufscar1[1:3] <- c("Scales", "Total_Obs", "Total_CI")
names.ufscar1[c(4,6,8,10)] <- paste(names(wav.ufscar1$Urochloa$ScaleSection_Obs),"Obs",sep="_")
names.ufscar1[c(5,7,9,11)] <- paste(names(wav.ufscar1$Urochloa$ScaleSection_CI),"CI",sep="_")

for(i in 1:length(wav.ufscar1$Urochloa$ScaleSection_CI)) {
  wav.ufscar1$Urochloa$ScaleSection_CI[[i]] <- t(wav.ufscar1$Urochloa$ScaleSection_CI[[i]])
}
wav.ufscar1.Urochloa <- matrix(NA,ncol=11, nrow=50)
colnames(wav.ufscar1.Urochloa) <- names.ufscar1
wav.ufscar1.Urochloa <- as.data.frame(wav.ufscar1.Urochloa)
wav.ufscar1.Urochloa$Scales <- wav.ufscar1$Urochloa$Scales
wav.ufscar1.Urochloa$Total_Obs <- wav.ufscar1$Urochloa$ScaleVar_Obs
wav.ufscar1.Urochloa$Total_CI <- wav.ufscar1$Urochloa$ScaleVar_CI[1,]
wav.ufscar1.Urochloa[c(4,6,8,10)] <- as.data.frame(wav.ufscar1$Urochloa$ScaleSection_Obs)
wav.ufscar1.Urochloa[c(5,7,9,11)] <- as.data.frame(wav.ufscar1$Urochloa$ScaleSection_CI)
write.table(wav.ufscar1.Urochloa, file="part3_wavScales_ufscar1_Urochloa.txt", quote=F, row.names=F, sep="\t")

for(i in 1:length(wav.ufscar1$Melinis$ScaleSection_CI)) {
  wav.ufscar1$Melinis$ScaleSection_CI[[i]] <- t(wav.ufscar1$Melinis$ScaleSection_CI[[i]])
}
wav.ufscar1.Melinis <- matrix(NA,ncol=11, nrow=50)
colnames(wav.ufscar1.Melinis) <- names.ufscar1
wav.ufscar1.Melinis <- as.data.frame(wav.ufscar1.Melinis)
wav.ufscar1.Melinis$Scales <- wav.ufscar1$Melinis$Scales
wav.ufscar1.Melinis$Total_Obs <- wav.ufscar1$Melinis$ScaleVar_Obs
wav.ufscar1.Melinis$Total_CI <- wav.ufscar1$Melinis$ScaleVar_CI[1,]
wav.ufscar1.Melinis[c(4,6,8,10)] <- as.data.frame(wav.ufscar1$Melinis$ScaleSection_Obs)
wav.ufscar1.Melinis[c(5,7,9,11)] <- as.data.frame(wav.ufscar1$Melinis$ScaleSection_CI)
write.table(wav.ufscar1.Melinis, file="part3_wavScales_ufscar1_Melinis.txt", quote=F, row.names=F, sep="\t")

for(i in 1:length(wav.ufscar1$Grasses$ScaleSection_CI)) {
  wav.ufscar1$Grasses$ScaleSection_CI[[i]] <- t(wav.ufscar1$Grasses$ScaleSection_CI[[i]])
}
wav.ufscar1.Grasses <- matrix(NA,ncol=11, nrow=50)
colnames(wav.ufscar1.Grasses) <- names.ufscar1
wav.ufscar1.Grasses <- as.data.frame(wav.ufscar1.Grasses)
wav.ufscar1.Grasses$Scales <- wav.ufscar1$Grasses$Scales
wav.ufscar1.Grasses$Total_Obs <- wav.ufscar1$Grasses$ScaleVar_Obs
wav.ufscar1.Grasses$Total_CI <- wav.ufscar1$Grasses$ScaleVar_CI[1,]
wav.ufscar1.Grasses[c(4,6,8,10)] <- as.data.frame(wav.ufscar1$Grasses$ScaleSection_Obs)
wav.ufscar1.Grasses[c(5,7,9,11)] <- as.data.frame(wav.ufscar1$Grasses$ScaleSection_CI)
write.table(wav.ufscar1.Grasses, file="part3_wavScales_ufscar1_Grasses.txt", quote=F, row.names=F, sep="\t")

for(i in 1:length(wav.ufscar1$Sedges$ScaleSection_CI)) {
  wav.ufscar1$Sedges$ScaleSection_CI[[i]] <- t(wav.ufscar1$Sedges$ScaleSection_CI[[i]])
}
wav.ufscar1.Sedges <- matrix(NA,ncol=11, nrow=50)
colnames(wav.ufscar1.Sedges) <- names.ufscar1
wav.ufscar1.Sedges <- as.data.frame(wav.ufscar1.Sedges)
wav.ufscar1.Sedges$Scales <- wav.ufscar1$Sedges$Scales
wav.ufscar1.Sedges$Total_Obs <- wav.ufscar1$Sedges$ScaleVar_Obs
wav.ufscar1.Sedges$Total_CI <- wav.ufscar1$Sedges$ScaleVar_CI[1,]
wav.ufscar1.Sedges[c(4,6,8,10)] <- as.data.frame(wav.ufscar1$Sedges$ScaleSection_Obs)
wav.ufscar1.Sedges[c(5,7,9,11)] <- as.data.frame(wav.ufscar1$Sedges$ScaleSection_CI)
write.table(wav.ufscar1.Sedges, file="part3_wavScales_ufscar1_Sedges.txt", quote=F, row.names=F, sep="\t")


names.ufscar2 <- numeric(13)
names.ufscar2[1:3] <- c("Scales", "Total_Obs", "Total_CI")
names.ufscar2[c(4,6,8,10,12)] <- paste(names(wav.ufscar2$Urochloa$ScaleSection_Obs),"Obs",sep="_")
names.ufscar2[c(5,7,9,11,13)] <- paste(names(wav.ufscar2$Urochloa$ScaleSection_CI),"CI",sep="_")

for(i in 1:length(wav.ufscar2$Urochloa$ScaleSection_CI)) {
  wav.ufscar2$Urochloa$ScaleSection_CI[[i]] <- t(wav.ufscar2$Urochloa$ScaleSection_CI[[i]])
}
wav.ufscar2.Urochloa <- matrix(NA,ncol=13, nrow=50)
colnames(wav.ufscar2.Urochloa) <- names.ufscar2
wav.ufscar2.Urochloa <- as.data.frame(wav.ufscar2.Urochloa)
wav.ufscar2.Urochloa$Scales <- wav.ufscar2$Urochloa$Scales
wav.ufscar2.Urochloa$Total_Obs <- wav.ufscar2$Urochloa$ScaleVar_Obs
wav.ufscar2.Urochloa$Total_CI <- wav.ufscar2$Urochloa$ScaleVar_CI[1,]
wav.ufscar2.Urochloa[c(4,6,8,10,12)] <- as.data.frame(wav.ufscar2$Urochloa$ScaleSection_Obs)
wav.ufscar2.Urochloa[c(5,7,9,11,13)] <- as.data.frame(wav.ufscar2$Urochloa$ScaleSection_CI)
write.table(wav.ufscar2.Urochloa, file="part3_wavScales_ufscar2_Urochloa.txt", quote=F, row.names=F, sep="\t")

for(i in 1:length(wav.ufscar2$Melinis$ScaleSection_CI)) {
  wav.ufscar2$Melinis$ScaleSection_CI[[i]] <- t(wav.ufscar2$Melinis$ScaleSection_CI[[i]])
}
wav.ufscar2.Melinis <- matrix(NA,ncol=13, nrow=50)
colnames(wav.ufscar2.Melinis) <- names.ufscar2
wav.ufscar2.Melinis <- as.data.frame(wav.ufscar2.Melinis)
wav.ufscar2.Melinis$Scales <- wav.ufscar2$Melinis$Scales
wav.ufscar2.Melinis$Total_Obs <- wav.ufscar2$Melinis$ScaleVar_Obs
wav.ufscar2.Melinis$Total_CI <- wav.ufscar2$Melinis$ScaleVar_CI[1,]
wav.ufscar2.Melinis[c(4,6,8,10,12)] <- as.data.frame(wav.ufscar2$Melinis$ScaleSection_Obs)
wav.ufscar2.Melinis[c(5,7,9,11,13)] <- as.data.frame(wav.ufscar2$Melinis$ScaleSection_CI)
write.table(wav.ufscar2.Melinis, file="part3_wavScales_ufscar2_Melinis.txt", quote=F, row.names=F, sep="\t")

for(i in 1:length(wav.ufscar2$Grasses$ScaleSection_CI)) {
  wav.ufscar2$Grasses$ScaleSection_CI[[i]] <- t(wav.ufscar2$Grasses$ScaleSection_CI[[i]])
}
wav.ufscar2.Grasses <- matrix(NA,ncol=13, nrow=50)
colnames(wav.ufscar2.Grasses) <- names.ufscar2
wav.ufscar2.Grasses <- as.data.frame(wav.ufscar2.Grasses)
wav.ufscar2.Grasses$Scales <- wav.ufscar2$Grasses$Scales
wav.ufscar2.Grasses$Total_Obs <- wav.ufscar2$Grasses$ScaleVar_Obs
wav.ufscar2.Grasses$Total_CI <- wav.ufscar2$Grasses$ScaleVar_CI[1,]
wav.ufscar2.Grasses[c(4,6,8,10,12)] <- as.data.frame(wav.ufscar2$Grasses$ScaleSection_Obs)
wav.ufscar2.Grasses[c(5,7,9,11,13)] <- as.data.frame(wav.ufscar2$Grasses$ScaleSection_CI)
write.table(wav.ufscar2.Grasses, file="part3_wavScales_ufscar2_Grasses.txt", quote=F, row.names=F, sep="\t")

for(i in 1:length(wav.ufscar2$Sedges$ScaleSection_CI)) {
  wav.ufscar2$Sedges$ScaleSection_CI[[i]] <- t(wav.ufscar2$Sedges$ScaleSection_CI[[i]])
}
wav.ufscar2.Sedges <- matrix(NA,ncol=13, nrow=50)
colnames(wav.ufscar2.Sedges) <- names.ufscar2
wav.ufscar2.Sedges <- as.data.frame(wav.ufscar2.Sedges)
wav.ufscar2.Sedges$Scales <- wav.ufscar2$Sedges$Scales
wav.ufscar2.Sedges$Total_Obs <- wav.ufscar2$Sedges$ScaleVar_Obs
wav.ufscar2.Sedges$Total_CI <- wav.ufscar2$Sedges$ScaleVar_CI[1,]
wav.ufscar2.Sedges[c(4,6,8,10,12)] <- as.data.frame(wav.ufscar2$Sedges$ScaleSection_Obs)
wav.ufscar2.Sedges[c(5,7,9,11,13)] <- as.data.frame(wav.ufscar2$Sedges$ScaleSection_CI)
write.table(wav.ufscar2.Sedges, file="part3_wavScales_ufscar2_Sedges.txt", quote=F, row.names=F, sep="\t")


# Supplementary figures
setwd("/home/pavel/Profissional/Pesquisa/MyPapers-Working/EIgrasses/figures")

scales.itirapina <- wav.itirapina$Urochloa$Scales
png(filename="SM3_scales_itirapina.png", height=35, width=30, unit="cm", res=300)
par(mfrow=c(5,4), mar=c(3,3,2,2), oma=c(3,3,1,1))
for(k in c(2,6,10,12,14)) {
  name.temp <- paste(names(wav.itirapina.Urochloa)[k],"Urochloa",sep="-")
  name.temp <- gsub("_Obs", "", name.temp)
  foo <- wav.itirapina.Urochloa[,k]
  bar <- wav.itirapina.Urochloa[,k+1]
  plot(foo ~ scales.itirapina, type="l", main=name.temp)
  lines(bar ~ scales.itirapina, col="gray30")

  name.temp <- paste(names(wav.itirapina.Melinis)[k],"Melinis",sep="-")
  name.temp <- gsub("_Obs", "", name.temp)
  foo <- wav.itirapina.Melinis[,k]
  bar <- wav.itirapina.Melinis[,k+1]
  plot(foo ~ scales.itirapina, type="l", main=name.temp)
  lines(bar ~ scales.itirapina, col="gray30")

  name.temp <- paste(names(wav.itirapina.Grasses)[k],"Grasses",sep="-")
  name.temp <- gsub("_Obs", "", name.temp)
  foo <- wav.itirapina.Grasses[,k]
  bar <- wav.itirapina.Grasses[,k+1]
  plot(foo ~ scales.itirapina, type="l", main=name.temp)
  lines(bar ~ scales.itirapina, col="gray30")

  name.temp <- paste(names(wav.itirapina.Sedges)[k],"Sedges",sep="-")
  name.temp <- gsub("_Obs", "", name.temp)
  foo <- wav.itirapina.Sedges[,k]
  bar <- wav.itirapina.Sedges[,k+1]
  plot(foo ~ scales.itirapina, type="l", main=name.temp)
  lines(bar ~ scales.itirapina, col="gray30")
}
mtext(side=1, text="Scale (m)", outer=T)
mtext(side=2, text="Scale variance", outer=T)
dev.off()


scales.ufscar1 <- wav.ufscar1$Urochloa$Scales
png(filename="SM3_scales_ufscar1.png", height=35, width=30, unit="cm", res=300)
par(mfrow=c(4,4), mar=c(3,3,2,2), oma=c(3,3,1,1))
for(k in c(2,4,8,10)) {
  name.temp <- paste(names(wav.ufscar1.Urochloa)[k],"Urochloa",sep="-")
  name.temp <- gsub("_Obs", "", name.temp)
  foo <- wav.ufscar1.Urochloa[,k]
  bar <- wav.ufscar1.Urochloa[,k+1]
  plot(foo ~ scales.ufscar1, type="l", main=name.temp)
  lines(bar ~ scales.ufscar1, col="gray30")

  name.temp <- paste(names(wav.ufscar1.Melinis)[k],"Melinis",sep="-")
  name.temp <- gsub("_Obs", "", name.temp)
  foo <- wav.ufscar1.Melinis[,k]
  bar <- wav.ufscar1.Melinis[,k+1]
  plot(foo ~ scales.ufscar1, type="l", main=name.temp)
  lines(bar ~ scales.ufscar1, col="gray30")

  name.temp <- paste(names(wav.ufscar1.Grasses)[k],"Grasses",sep="-")
  name.temp <- gsub("_Obs", "", name.temp)
  foo <- wav.ufscar1.Grasses[,k]
  bar <- wav.ufscar1.Grasses[,k+1]
  plot(foo ~ scales.ufscar1, type="l", main=name.temp)
  lines(bar ~ scales.ufscar1, col="gray30")

  name.temp <- paste(names(wav.ufscar1.Sedges)[k],"Sedges",sep="-")
  name.temp <- gsub("_Obs", "", name.temp)
  foo <- wav.ufscar1.Sedges[,k]
  bar <- wav.ufscar1.Sedges[,k+1]
  plot(foo ~ scales.ufscar1, type="l", main=name.temp)
  lines(bar ~ scales.ufscar1, col="gray30")
}
mtext(side=1, text="Scale (m)", outer=T)
mtext(side=2, text="Scale variance", outer=T)
dev.off()



scales.ufscar2 <- wav.ufscar2$Urochloa$Scales
png(filename="SM3_scales_ufscar2.png", height=35, width=30, unit="cm", res=300)
par(mfrow=c(5,4), mar=c(3,3,2,2), oma=c(3,3,1,1))
for(k in c(2,6,8,10,12)) {
  name.temp <- paste(names(wav.ufscar2.Urochloa)[k],"Urochloa",sep="-")
  name.temp <- gsub("_Obs", "", name.temp)
  foo <- wav.ufscar2.Urochloa[,k]
  bar <- wav.ufscar2.Urochloa[,k+1]
  plot(foo ~ scales.ufscar2, type="l", main=name.temp)
  lines(bar ~ scales.ufscar2, col="gray30")

  name.temp <- paste(names(wav.ufscar2.Melinis)[k],"Melinis",sep="-")
  name.temp <- gsub("_Obs", "", name.temp)
  foo <- wav.ufscar2.Melinis[,k]
  bar <- wav.ufscar2.Melinis[,k+1]
  plot(foo ~ scales.ufscar2, type="l", main=name.temp)
  lines(bar ~ scales.ufscar2, col="gray30")

  name.temp <- paste(names(wav.ufscar2.Grasses)[k],"Grasses",sep="-")
  name.temp <- gsub("_Obs", "", name.temp)
  foo <- wav.ufscar2.Grasses[,k]
  bar <- wav.ufscar2.Grasses[,k+1]
  plot(foo ~ scales.ufscar2, type="l", main=name.temp)
  lines(bar ~ scales.ufscar2, col="gray30")

  name.temp <- paste(names(wav.ufscar2.Sedges)[k],"Sedges",sep="-")
  name.temp <- gsub("_Obs", "", name.temp)
  foo <- wav.ufscar2.Sedges[,k]
  bar <- wav.ufscar2.Sedges[,k+1]
  plot(foo ~ scales.ufscar2, type="l", main=name.temp)
  lines(bar ~ scales.ufscar2, col="gray30")
}
mtext(side=1, text="Scale (m)", outer=T)
mtext(side=2, text="Scale variance", outer=T)
dev.off()


### Part 4 - At what scales are the invasive and native graminoids negatively or positively related?

# Simulate the null models...
Nperm=5000

rand.itirapina <- lapply(data.itirapina[,1:4], MC1s, sections=patches.itirapina, Nperm=Nperm)
rand.ufscar1 <- lapply(data.ufscar1[,1:4], MC1s, sections=patches.ufscar1, Nperm=Nperm)
rand.ufscar2 <- lapply(data.ufscar2[,1:4], MC1s, sections=patches.ufscar2, Nperm=Nperm)

# Apply the wavelet :-)

wav2.itirapina.UrochloaVsGrasses <- wavCWTvarCIs.bivar(x=rand.itirapina$Urochloa, y=rand.itirapina$Grasses, significance=T, sections=patches.itirapina, zero.padding=T, remove.COI=T, test=F, scale.max=75, scale.var=T, pos.var=F, CI=T, CIquant=c(0.05), effect.size=F,  keep.all=F, print.loop=T, which.loop=100, save.cwt=F, wav.template="gaussian2", scale.per.section=T)
wav2.itirapina.UrochloaVsSedges <- wavCWTvarCIs.bivar(x=rand.itirapina$Urochloa, y=rand.itirapina$Sedges, significance=T, sections=patches.itirapina, zero.padding=T, remove.COI=T, test=F, scale.max=75, scale.var=T, pos.var=F, CI=T, CIquant=c(0.05), effect.size=F,  keep.all=F, print.loop=T, which.loop=100, save.cwt=F, wav.template="gaussian2", scale.per.section=T)
wav2.itirapina.MelinisVsGrasses <- wavCWTvarCIs.bivar(x=rand.itirapina$Melinis, y=rand.itirapina$Grasses, significance=T, sections=patches.itirapina, zero.padding=T, remove.COI=T, test=F, scale.max=75, scale.var=T, pos.var=F, CI=T, CIquant=c(0.05), effect.size=F,  keep.all=F, print.loop=T, which.loop=100, save.cwt=F, wav.template="gaussian2", scale.per.section=T)
wav2.itirapina.MelinisVsSedges <- wavCWTvarCIs.bivar(x=rand.itirapina$Melinis, y=rand.itirapina$Sedges, significance=T, sections=patches.itirapina, zero.padding=T, remove.COI=T, test=F, scale.max=75, scale.var=T, pos.var=F, CI=T, CIquant=c(0.05), effect.size=F,  keep.all=F, print.loop=T, which.loop=100, save.cwt=F, wav.template="gaussian2", scale.per.section=T)
wav2.itirapina.UrochloaVsMelinis <- wavCWTvarCIs.bivar(x=rand.itirapina$Urochloa, y=rand.itirapina$Melinis, significance=T, sections=patches.itirapina, zero.padding=T, remove.COI=T, test=F, scale.max=75, scale.var=T, pos.var=F, CI=T, CIquant=c(0.025, 0.975), effect.size=F,  keep.all=F, print.loop=T, which.loop=100, save.cwt=F, wav.template="gaussian2", scale.per.section=T)

wav2.ufscar1.UrochloaVsGrasses <- wavCWTvarCIs.bivar(x=rand.ufscar1$Urochloa, y=rand.ufscar1$Grasses, significance=T, sections=patches.ufscar1, zero.padding=T, remove.COI=T, test=F, scale.max=75, scale.var=T, pos.var=F, CI=T, CIquant=c(0.05), effect.size=F,  keep.all=F, print.loop=T, which.loop=100, save.cwt=F, wav.template="gaussian2", scale.per.section=T)
wav2.ufscar1.UrochloaVsSedges <- wavCWTvarCIs.bivar(x=rand.ufscar1$Urochloa, y=rand.ufscar1$Sedges, significance=T, sections=patches.ufscar1, zero.padding=T, remove.COI=T, test=F, scale.max=75, scale.var=T, pos.var=F, CI=T, CIquant=c(0.05), effect.size=F,  keep.all=F, print.loop=T, which.loop=100, save.cwt=F, wav.template="gaussian2", scale.per.section=T)
wav2.ufscar1.MelinisVsGrasses <- wavCWTvarCIs.bivar(x=rand.ufscar1$Melinis, y=rand.ufscar1$Grasses, significance=T, sections=patches.ufscar1, zero.padding=T, remove.COI=T, test=F, scale.max=75, scale.var=T, pos.var=F, CI=T, CIquant=c(0.05), effect.size=F,  keep.all=F, print.loop=T, which.loop=100, save.cwt=F, wav.template="gaussian2", scale.per.section=T)
wav2.ufscar1.MelinisVsSedges <- wavCWTvarCIs.bivar(x=rand.ufscar1$Melinis, y=rand.ufscar1$Sedges, significance=T, sections=patches.ufscar1, zero.padding=T, remove.COI=T, test=F, scale.max=75, scale.var=T, pos.var=F, CI=T, CIquant=c(0.05), effect.size=F,  keep.all=F, print.loop=T, which.loop=100, save.cwt=F, wav.template="gaussian2", scale.per.section=T)
wav2.ufscar1.UrochloaVsMelinis <- wavCWTvarCIs.bivar(x=rand.ufscar1$Urochloa, y=rand.ufscar1$Melinis, significance=T, sections=patches.ufscar1, zero.padding=T, remove.COI=T, test=F, scale.max=75, scale.var=T, pos.var=F, CI=T, CIquant=c(0.025, 0.975), effect.size=F,  keep.all=F, print.loop=T, which.loop=100, save.cwt=F, wav.template="gaussian2", scale.per.section=T)

wav2.ufscar2.UrochloaVsGrasses <- wavCWTvarCIs.bivar(x=rand.ufscar2$Urochloa, y=rand.ufscar2$Grasses, significance=T, sections=patches.ufscar2, zero.padding=T, remove.COI=T, test=F, scale.max=75, scale.var=T, pos.var=F, CI=T, CIquant=c(0.05), effect.size=F,  keep.all=F, print.loop=T, which.loop=100, save.cwt=F, wav.template="gaussian2", scale.per.section=T)
wav2.ufscar2.UrochloaVsSedges <- wavCWTvarCIs.bivar(x=rand.ufscar2$Urochloa, y=rand.ufscar2$Sedges, significance=T, sections=patches.ufscar2, zero.padding=T, remove.COI=T, test=F, scale.max=75, scale.var=T, pos.var=F, CI=T, CIquant=c(0.05), effect.size=F,  keep.all=F, print.loop=T, which.loop=100, save.cwt=F, wav.template="gaussian2", scale.per.section=T)
wav2.ufscar2.MelinisVsGrasses <- wavCWTvarCIs.bivar(x=rand.ufscar2$Melinis, y=rand.ufscar2$Grasses, significance=T, sections=patches.ufscar2, zero.padding=T, remove.COI=T, test=F, scale.max=75, scale.var=T, pos.var=F, CI=T, CIquant=c(0.05), effect.size=F,  keep.all=F, print.loop=T, which.loop=100, save.cwt=F, wav.template="gaussian2", scale.per.section=T)
wav2.ufscar2.MelinisVsSedges <- wavCWTvarCIs.bivar(x=rand.ufscar2$Melinis, y=rand.ufscar2$Sedges, significance=T, sections=patches.ufscar2, zero.padding=T, remove.COI=T, test=F, scale.max=75, scale.var=T, pos.var=F, CI=T, CIquant=c(0.05), effect.size=F,  keep.all=F, print.loop=T, which.loop=100, save.cwt=F, wav.template="gaussian2", scale.per.section=T)
wav2.ufscar2.UrochloaVsMelinis <- wavCWTvarCIs.bivar(x=rand.ufscar2$Urochloa, y=rand.ufscar2$Melinis, significance=T, sections=patches.ufscar2, zero.padding=T, remove.COI=T, test=F, scale.max=75, scale.var=T, pos.var=F, CI=T, CIquant=c(0.025, 0.975), effect.size=F,  keep.all=F, print.loop=T, which.loop=100, save.cwt=F, wav.template="gaussian2", scale.per.section=T)


save.image(file="wav_bivar.RData")

# Organize as something prettier and save the output as .txt files
# Once again, I could have used a loop, but this way is simpler.

setwd("/home/pavel/Profissional/Pesquisa/MyPapers-Working/EIgrasses/results")

names.itirapina.a <- numeric(15)
names.itirapina.a[1:3] <- c("Scales", "Total_Obs", "Total_CI_inf")
names.itirapina.a[c(4,6,8,10,12,14)] <- paste(names(wav2.itirapina.UrochloaVsGrasses$ScaleSection_Obs),"Obs",sep="_")
names.itirapina.a[c(5,7,9,11,13,15)] <- paste(names(wav2.itirapina.UrochloaVsGrasses$ScaleSection_CI),"CI_inf",sep="_")
names.itirapina.b <- numeric(22)
names.itirapina.b[1:4] <- c("Scales", "Total_Obs", "Total_CI_inf", "Total_CI_sup")
names.itirapina.b[c(5,8,11,14,17,20)] <- paste(names(wav2.itirapina.UrochloaVsGrasses$ScaleSection_Obs),"Obs",sep="_")
names.itirapina.b[c(6,9,12,15,18,21)] <- paste(names(wav2.itirapina.UrochloaVsGrasses$ScaleSection_CI),"CI_inf",sep="_")
names.itirapina.b[c(7,10,13,16,19,22)] <- paste(names(wav2.itirapina.UrochloaVsGrasses$ScaleSection_CI),"CI_sup",sep="_")

for(i in 1:length(wav2.itirapina.UrochloaVsGrasses$ScaleSection_CI)) {
  wav2.itirapina.UrochloaVsGrasses$ScaleSection_CI[[i]] <- t(wav2.itirapina.UrochloaVsGrasses$ScaleSection_CI[[i]])
}
for(i in 1:length(wav2.itirapina.UrochloaVsSedges$ScaleSection_CI)) {
  wav2.itirapina.UrochloaVsSedges$ScaleSection_CI[[i]] <- t(wav2.itirapina.UrochloaVsSedges$ScaleSection_CI[[i]])
}
for(i in 1:length(wav2.itirapina.MelinisVsGrasses$ScaleSection_CI)) {
  wav2.itirapina.MelinisVsGrasses$ScaleSection_CI[[i]] <- t(wav2.itirapina.MelinisVsGrasses$ScaleSection_CI[[i]])
}
for(i in 1:length(wav2.itirapina.MelinisVsSedges$ScaleSection_CI)) {
  wav2.itirapina.MelinisVsSedges$ScaleSection_CI[[i]] <- t(wav2.itirapina.MelinisVsSedges$ScaleSection_CI[[i]])
}

wav2m.itirapina.UrochloaVsGrasses <- matrix(NA,ncol=15, nrow=50)
colnames(wav2m.itirapina.UrochloaVsGrasses) <- names.itirapina.a
wav2m.itirapina.UrochloaVsGrasses <- as.data.frame(wav2m.itirapina.UrochloaVsGrasses)
wav2m.itirapina.UrochloaVsGrasses$Scales <- wav2.itirapina.UrochloaVsGrasses$Scales
wav2m.itirapina.UrochloaVsGrasses$Total_Obs <- wav2.itirapina.UrochloaVsGrasses$ScaleVar_Obs
wav2m.itirapina.UrochloaVsGrasses$Total_CI_inf <- wav2.itirapina.UrochloaVsGrasses$ScaleVar_CI[1,]
wav2m.itirapina.UrochloaVsGrasses[c(4,6,8,10,12,14)] <- as.data.frame(wav2.itirapina.UrochloaVsGrasses$ScaleSection_Obs)
wav2m.itirapina.UrochloaVsGrasses[c(5,7,9,11,13,15)] <- as.data.frame(wav2.itirapina.UrochloaVsGrasses$ScaleSection_CI)
write.table(wav2m.itirapina.UrochloaVsGrasses, file="part4_wav2Scales_itirapina_UrochloaVsGrasses.txt", quote=F, row.names=F, sep="\t")

wav2m.itirapina.UrochloaVsSedges <- matrix(NA,ncol=15, nrow=50)
colnames(wav2m.itirapina.UrochloaVsSedges) <- names.itirapina.a
wav2m.itirapina.UrochloaVsSedges <- as.data.frame(wav2m.itirapina.UrochloaVsSedges)
wav2m.itirapina.UrochloaVsSedges$Scales <- wav2.itirapina.UrochloaVsSedges$Scales
wav2m.itirapina.UrochloaVsSedges$Total_Obs <- wav2.itirapina.UrochloaVsSedges$ScaleVar_Obs
wav2m.itirapina.UrochloaVsSedges$Total_CI_inf <- wav2.itirapina.UrochloaVsSedges$ScaleVar_CI[1,]
wav2m.itirapina.UrochloaVsSedges[c(4,6,8,10,12,14)] <- as.data.frame(wav2.itirapina.UrochloaVsSedges$ScaleSection_Obs)
wav2m.itirapina.UrochloaVsSedges[c(5,7,9,11,13,15)] <- as.data.frame(wav2.itirapina.UrochloaVsSedges$ScaleSection_CI)
write.table(wav2m.itirapina.UrochloaVsSedges, file="part4_wav2Scales_itirapina_UrochloaVsSedges.txt", quote=F, row.names=F, sep="\t")

wav2m.itirapina.MelinisVsGrasses <- matrix(NA,ncol=15, nrow=50)
colnames(wav2m.itirapina.MelinisVsGrasses) <- names.itirapina.a
wav2m.itirapina.MelinisVsGrasses <- as.data.frame(wav2m.itirapina.MelinisVsGrasses)
wav2m.itirapina.MelinisVsGrasses$Scales <- wav2.itirapina.MelinisVsGrasses$Scales
wav2m.itirapina.MelinisVsGrasses$Total_Obs <- wav2.itirapina.MelinisVsGrasses$ScaleVar_Obs
wav2m.itirapina.MelinisVsGrasses$Total_CI_inf <- wav2.itirapina.MelinisVsGrasses$ScaleVar_CI[1,]
wav2m.itirapina.MelinisVsGrasses[c(4,6,8,10,12,14)] <- as.data.frame(wav2.itirapina.MelinisVsGrasses$ScaleSection_Obs)
wav2m.itirapina.MelinisVsGrasses[c(5,7,9,11,13,15)] <- as.data.frame(wav2.itirapina.MelinisVsGrasses$ScaleSection_CI)
write.table(wav2m.itirapina.MelinisVsGrasses, file="part4_wav2Scales_itirapina_MelinisVsGrasses.txt", quote=F, row.names=F, sep="\t")

wav2m.itirapina.MelinisVsSedges <- matrix(NA,ncol=15, nrow=50)
colnames(wav2m.itirapina.MelinisVsSedges) <- names.itirapina.a
wav2m.itirapina.MelinisVsSedges <- as.data.frame(wav2m.itirapina.MelinisVsSedges)
wav2m.itirapina.MelinisVsSedges$Scales <- wav2.itirapina.MelinisVsSedges$Scales
wav2m.itirapina.MelinisVsSedges$Total_Obs <- wav2.itirapina.MelinisVsSedges$ScaleVar_Obs
wav2m.itirapina.MelinisVsSedges$Total_CI_inf <- wav2.itirapina.MelinisVsSedges$ScaleVar_CI[1,]
wav2m.itirapina.MelinisVsSedges[c(4,6,8,10,12,14)] <- as.data.frame(wav2.itirapina.MelinisVsSedges$ScaleSection_Obs)
wav2m.itirapina.MelinisVsSedges[c(5,7,9,11,13,15)] <- as.data.frame(wav2.itirapina.MelinisVsSedges$ScaleSection_CI)
write.table(wav2m.itirapina.MelinisVsSedges, file="part4_wav2Scales_itirapina_MelinisVsSedges.txt", quote=F, row.names=F, sep="\t")

wav2m.itirapina.UrochloaVsMelinis <- matrix(NA,ncol=22, nrow=50)
colnames(wav2m.itirapina.UrochloaVsMelinis) <- names.itirapina.b
wav2m.itirapina.UrochloaVsMelinis <- as.data.frame(wav2m.itirapina.UrochloaVsMelinis)
wav2m.itirapina.UrochloaVsMelinis <- as.data.frame(wav2m.itirapina.UrochloaVsMelinis)
wav2m.itirapina.UrochloaVsMelinis$Scales <- wav2.itirapina.UrochloaVsMelinis$Scales
wav2m.itirapina.UrochloaVsMelinis$Total_Obs <- wav2.itirapina.UrochloaVsMelinis$ScaleVar_Obs
wav2m.itirapina.UrochloaVsMelinis$Total_CI_inf <- wav2.itirapina.UrochloaVsMelinis$ScaleVar_CI[,1]
wav2m.itirapina.UrochloaVsMelinis$Total_CI_sup <- wav2.itirapina.UrochloaVsMelinis$ScaleVar_CI[,2]
wav2m.itirapina.UrochloaVsMelinis[c(5,8,11,14,17,20)] <- as.data.frame(wav2.itirapina.UrochloaVsMelinis$ScaleSection_Obs)
wav2m.itirapina.UrochloaVsMelinis[c(6,7,9,10,12,13,15,16,18,19,21,22)] <- as.data.frame(wav2.itirapina.UrochloaVsMelinis$ScaleSection_CI)
write.table(wav2m.itirapina.UrochloaVsMelinis, file="part4_wav2Scales_itirapina_UrochloaVsMelinis.txt", quote=F, row.names=F, sep="\t")

names.ufscar2.a <- numeric(13)
names.ufscar2.a[1:3] <- c("Scales", "Total_Obs", "Total_CI_inf")
names.ufscar2.a[c(4,6,8,10,12)] <- paste(names(wav2.ufscar2.UrochloaVsGrasses$ScaleSection_Obs),"Obs",sep="_")
names.ufscar2.a[c(5,7,9,11,13)] <- paste(names(wav2.ufscar2.UrochloaVsGrasses$ScaleSection_CI),"CI_inf",sep="_")
names.ufscar2.b <- numeric(19)
names.ufscar2.b[1:4] <- c("Scales", "Total_Obs", "Total_CI_inf", "Total_CI_sup")
names.ufscar2.b[c(5,8,11,14,17)] <- paste(names(wav2.ufscar2.UrochloaVsGrasses$ScaleSection_Obs),"Obs",sep="_")
names.ufscar2.b[c(6,9,12,15,18)] <- paste(names(wav2.ufscar2.UrochloaVsGrasses$ScaleSection_CI),"CI_inf",sep="_")
names.ufscar2.b[c(7,10,13,16,19)] <- paste(names(wav2.ufscar2.UrochloaVsGrasses$ScaleSection_CI),"CI_sup",sep="_")

for(i in 1:length(wav2.ufscar2.UrochloaVsGrasses$ScaleSection_CI)) {
  wav2.ufscar2.UrochloaVsGrasses$ScaleSection_CI[[i]] <- t(wav2.ufscar2.UrochloaVsGrasses$ScaleSection_CI[[i]])
}
for(i in 1:length(wav2.ufscar2.UrochloaVsSedges$ScaleSection_CI)) {
  wav2.ufscar2.UrochloaVsSedges$ScaleSection_CI[[i]] <- t(wav2.ufscar2.UrochloaVsSedges$ScaleSection_CI[[i]])
}
for(i in 1:length(wav2.ufscar2.MelinisVsGrasses$ScaleSection_CI)) {
  wav2.ufscar2.MelinisVsGrasses$ScaleSection_CI[[i]] <- t(wav2.ufscar2.MelinisVsGrasses$ScaleSection_CI[[i]])
}
for(i in 1:length(wav2.ufscar2.MelinisVsSedges$ScaleSection_CI)) {
  wav2.ufscar2.MelinisVsSedges$ScaleSection_CI[[i]] <- t(wav2.ufscar2.MelinisVsSedges$ScaleSection_CI[[i]])
}

wav2m.ufscar2.UrochloaVsGrasses <- matrix(NA,ncol=13, nrow=50)
colnames(wav2m.ufscar2.UrochloaVsGrasses) <- names.ufscar2.a
wav2m.ufscar2.UrochloaVsGrasses <- as.data.frame(wav2m.ufscar2.UrochloaVsGrasses)
wav2m.ufscar2.UrochloaVsGrasses$Scales <- wav2.ufscar2.UrochloaVsGrasses$Scales
wav2m.ufscar2.UrochloaVsGrasses$Total_Obs <- wav2.ufscar2.UrochloaVsGrasses$ScaleVar_Obs
wav2m.ufscar2.UrochloaVsGrasses$Total_CI_inf <- wav2.ufscar2.UrochloaVsGrasses$ScaleVar_CI[1,]
wav2m.ufscar2.UrochloaVsGrasses[c(4,6,8,10,12)] <- as.data.frame(wav2.ufscar2.UrochloaVsGrasses$ScaleSection_Obs)
wav2m.ufscar2.UrochloaVsGrasses[c(5,7,9,11,13)] <- as.data.frame(wav2.ufscar2.UrochloaVsGrasses$ScaleSection_CI)
write.table(wav2m.ufscar2.UrochloaVsGrasses, file="part4_wav2Scales_ufscar2_UrochloaVsGrasses.txt", quote=F, row.names=F, sep="\t")

wav2m.ufscar2.UrochloaVsSedges <- matrix(NA,ncol=13, nrow=50)
colnames(wav2m.ufscar2.UrochloaVsSedges) <- names.ufscar2.a
wav2m.ufscar2.UrochloaVsSedges <- as.data.frame(wav2m.ufscar2.UrochloaVsSedges)
wav2m.ufscar2.UrochloaVsSedges$Scales <- wav2.ufscar2.UrochloaVsSedges$Scales
wav2m.ufscar2.UrochloaVsSedges$Total_Obs <- wav2.ufscar2.UrochloaVsSedges$ScaleVar_Obs
wav2m.ufscar2.UrochloaVsSedges$Total_CI_inf <- wav2.ufscar2.UrochloaVsSedges$ScaleVar_CI[1,]
wav2m.ufscar2.UrochloaVsSedges[c(4,6,8,10,12)] <- as.data.frame(wav2.ufscar2.UrochloaVsSedges$ScaleSection_Obs)
wav2m.ufscar2.UrochloaVsSedges[c(5,7,9,11,13)] <- as.data.frame(wav2.ufscar2.UrochloaVsSedges$ScaleSection_CI)
write.table(wav2m.ufscar2.UrochloaVsSedges, file="part4_wav2Scales_ufscar2_UrochloaVsSedges.txt", quote=F, row.names=F, sep="\t")

wav2m.ufscar2.MelinisVsGrasses <- matrix(NA,ncol=13, nrow=50)
colnames(wav2m.ufscar2.MelinisVsGrasses) <- names.ufscar2.a
wav2m.ufscar2.MelinisVsGrasses <- as.data.frame(wav2m.ufscar2.MelinisVsGrasses)
wav2m.ufscar2.MelinisVsGrasses$Scales <- wav2.ufscar2.MelinisVsGrasses$Scales
wav2m.ufscar2.MelinisVsGrasses$Total_Obs <- wav2.ufscar2.MelinisVsGrasses$ScaleVar_Obs
wav2m.ufscar2.MelinisVsGrasses$Total_CI_inf <- wav2.ufscar2.MelinisVsGrasses$ScaleVar_CI[1,]
wav2m.ufscar2.MelinisVsGrasses[c(4,6,8,10,12)] <- as.data.frame(wav2.ufscar2.MelinisVsGrasses$ScaleSection_Obs)
wav2m.ufscar2.MelinisVsGrasses[c(5,7,9,11,13)] <- as.data.frame(wav2.ufscar2.MelinisVsGrasses$ScaleSection_CI)
write.table(wav2m.ufscar2.MelinisVsGrasses, file="part4_wav2Scales_ufscar2_MelinisVsGrasses.txt", quote=F, row.names=F, sep="\t")

wav2m.ufscar2.MelinisVsSedges <- matrix(NA,ncol=13, nrow=50)
colnames(wav2m.ufscar2.MelinisVsSedges) <- names.ufscar2.a
wav2m.ufscar2.MelinisVsSedges <- as.data.frame(wav2m.ufscar2.MelinisVsSedges)
wav2m.ufscar2.MelinisVsSedges$Scales <- wav2.ufscar2.MelinisVsSedges$Scales
wav2m.ufscar2.MelinisVsSedges$Total_Obs <- wav2.ufscar2.MelinisVsSedges$ScaleVar_Obs
wav2m.ufscar2.MelinisVsSedges$Total_CI_inf <- wav2.ufscar2.MelinisVsSedges$ScaleVar_CI[1,]
wav2m.ufscar2.MelinisVsSedges[c(4,6,8,10,12)] <- as.data.frame(wav2.ufscar2.MelinisVsSedges$ScaleSection_Obs)
wav2m.ufscar2.MelinisVsSedges[c(5,7,9,11,13)] <- as.data.frame(wav2.ufscar2.MelinisVsSedges$ScaleSection_CI)
write.table(wav2m.ufscar2.MelinisVsSedges, file="part4_wav2Scales_ufscar2_MelinisVsSedges.txt", quote=F, row.names=F, sep="\t")

wav2m.ufscar2.UrochloaVsMelinis <- matrix(NA,ncol=19, nrow=50)
colnames(wav2m.ufscar2.UrochloaVsMelinis) <- names.ufscar2.b
wav2m.ufscar2.UrochloaVsMelinis <- as.data.frame(wav2m.ufscar2.UrochloaVsMelinis)
wav2m.ufscar2.UrochloaVsMelinis <- as.data.frame(wav2m.ufscar2.UrochloaVsMelinis)
wav2m.ufscar2.UrochloaVsMelinis$Scales <- wav2.ufscar2.UrochloaVsMelinis$Scales
wav2m.ufscar2.UrochloaVsMelinis$Total_Obs <- wav2.ufscar2.UrochloaVsMelinis$ScaleVar_Obs
wav2m.ufscar2.UrochloaVsMelinis$Total_CI_inf <- wav2.ufscar2.UrochloaVsMelinis$ScaleVar_CI[,1]
wav2m.ufscar2.UrochloaVsMelinis$Total_CI_sup <- wav2.ufscar2.UrochloaVsMelinis$ScaleVar_CI[,2]
wav2m.ufscar2.UrochloaVsMelinis[c(5,8,11,14,17)] <- as.data.frame(wav2.ufscar2.UrochloaVsMelinis$ScaleSection_Obs)
wav2m.ufscar2.UrochloaVsMelinis[c(6,7,9,10,12,13,15,16,18,19)] <- as.data.frame(wav2.ufscar2.UrochloaVsMelinis$ScaleSection_CI)
write.table(wav2m.ufscar2.UrochloaVsMelinis, file="part4_wav2Scales_ufscar2_UrochloaVsMelinis.txt", quote=F, row.names=F, sep="\t")


names.ufscar1.a <- numeric(11)
names.ufscar1.a[1:3] <- c("Scales", "Total_Obs", "Total_CI_inf")
names.ufscar1.a[c(4,6,8,10)] <- paste(names(wav2.ufscar1.UrochloaVsGrasses$ScaleSection_Obs),"Obs",sep="_")
names.ufscar1.a[c(5,7,9,11)] <- paste(names(wav2.ufscar1.UrochloaVsGrasses$ScaleSection_CI),"CI_inf",sep="_")
names.ufscar1.b <- numeric(16)
names.ufscar1.b[1:4] <- c("Scales", "Total_Obs", "Total_CI_inf", "Total_CI_sup")
names.ufscar1.b[c(5,8,11,14)] <- paste(names(wav2.ufscar1.UrochloaVsGrasses$ScaleSection_Obs),"Obs",sep="_")
names.ufscar1.b[c(6,9,12,15)] <- paste(names(wav2.ufscar1.UrochloaVsGrasses$ScaleSection_CI),"CI_inf",sep="_")
names.ufscar1.b[c(7,10,13,16)] <- paste(names(wav2.ufscar1.UrochloaVsGrasses$ScaleSection_CI),"CI_sup",sep="_")

for(i in 1:length(wav2.ufscar1.UrochloaVsGrasses$ScaleSection_CI)) {
  wav2.ufscar1.UrochloaVsGrasses$ScaleSection_CI[[i]] <- t(wav2.ufscar1.UrochloaVsGrasses$ScaleSection_CI[[i]])
}
for(i in 1:length(wav2.ufscar1.UrochloaVsSedges$ScaleSection_CI)) {
  wav2.ufscar1.UrochloaVsSedges$ScaleSection_CI[[i]] <- t(wav2.ufscar1.UrochloaVsSedges$ScaleSection_CI[[i]])
}
for(i in 1:length(wav2.ufscar1.MelinisVsGrasses$ScaleSection_CI)) {
  wav2.ufscar1.MelinisVsGrasses$ScaleSection_CI[[i]] <- t(wav2.ufscar1.MelinisVsGrasses$ScaleSection_CI[[i]])
}
for(i in 1:length(wav2.ufscar1.MelinisVsSedges$ScaleSection_CI)) {
  wav2.ufscar1.MelinisVsSedges$ScaleSection_CI[[i]] <- t(wav2.ufscar1.MelinisVsSedges$ScaleSection_CI[[i]])
}

wav2m.ufscar1.UrochloaVsGrasses <- matrix(NA,ncol=11, nrow=50)
colnames(wav2m.ufscar1.UrochloaVsGrasses) <- names.ufscar1.a
wav2m.ufscar1.UrochloaVsGrasses <- as.data.frame(wav2m.ufscar1.UrochloaVsGrasses)
wav2m.ufscar1.UrochloaVsGrasses$Scales <- wav2.ufscar1.UrochloaVsGrasses$Scales
wav2m.ufscar1.UrochloaVsGrasses$Total_Obs <- wav2.ufscar1.UrochloaVsGrasses$ScaleVar_Obs
wav2m.ufscar1.UrochloaVsGrasses$Total_CI_inf <- wav2.ufscar1.UrochloaVsGrasses$ScaleVar_CI[1,]
wav2m.ufscar1.UrochloaVsGrasses[c(4,6,8,10)] <- as.data.frame(wav2.ufscar1.UrochloaVsGrasses$ScaleSection_Obs)
wav2m.ufscar1.UrochloaVsGrasses[c(5,7,9,11)] <- as.data.frame(wav2.ufscar1.UrochloaVsGrasses$ScaleSection_CI)
write.table(wav2m.ufscar1.UrochloaVsGrasses, file="part4_wav2Scales_ufscar1_UrochloaVsGrasses.txt", quote=F, row.names=F, sep="\t")

wav2m.ufscar1.UrochloaVsSedges <- matrix(NA,ncol=11, nrow=50)
colnames(wav2m.ufscar1.UrochloaVsSedges) <- names.ufscar1.a
wav2m.ufscar1.UrochloaVsSedges <- as.data.frame(wav2m.ufscar1.UrochloaVsSedges)
wav2m.ufscar1.UrochloaVsSedges$Scales <- wav2.ufscar1.UrochloaVsSedges$Scales
wav2m.ufscar1.UrochloaVsSedges$Total_Obs <- wav2.ufscar1.UrochloaVsSedges$ScaleVar_Obs
wav2m.ufscar1.UrochloaVsSedges$Total_CI_inf <- wav2.ufscar1.UrochloaVsSedges$ScaleVar_CI[1,]
wav2m.ufscar1.UrochloaVsSedges[c(4,6,8,10)] <- as.data.frame(wav2.ufscar1.UrochloaVsSedges$ScaleSection_Obs)
wav2m.ufscar1.UrochloaVsSedges[c(5,7,9,11)] <- as.data.frame(wav2.ufscar1.UrochloaVsSedges$ScaleSection_CI)
write.table(wav2m.ufscar1.UrochloaVsSedges, file="part4_wav2Scales_ufscar1_UrochloaVsSedges.txt", quote=F, row.names=F, sep="\t")

wav2m.ufscar1.MelinisVsGrasses <- matrix(NA,ncol=11, nrow=50)
colnames(wav2m.ufscar1.MelinisVsGrasses) <- names.ufscar1.a
wav2m.ufscar1.MelinisVsGrasses <- as.data.frame(wav2m.ufscar1.MelinisVsGrasses)
wav2m.ufscar1.MelinisVsGrasses$Scales <- wav2.ufscar1.MelinisVsGrasses$Scales
wav2m.ufscar1.MelinisVsGrasses$Total_Obs <- wav2.ufscar1.MelinisVsGrasses$ScaleVar_Obs
wav2m.ufscar1.MelinisVsGrasses$Total_CI_inf <- wav2.ufscar1.MelinisVsGrasses$ScaleVar_CI[1,]
wav2m.ufscar1.MelinisVsGrasses[c(4,6,8,10)] <- as.data.frame(wav2.ufscar1.MelinisVsGrasses$ScaleSection_Obs)
wav2m.ufscar1.MelinisVsGrasses[c(5,7,9,11)] <- as.data.frame(wav2.ufscar1.MelinisVsGrasses$ScaleSection_CI)
write.table(wav2m.ufscar1.MelinisVsGrasses, file="part4_wav2Scales_ufscar1_MelinisVsGrasses.txt", quote=F, row.names=F, sep="\t")

wav2m.ufscar1.MelinisVsSedges <- matrix(NA,ncol=11, nrow=50)
colnames(wav2m.ufscar1.MelinisVsSedges) <- names.ufscar1.a
wav2m.ufscar1.MelinisVsSedges <- as.data.frame(wav2m.ufscar1.MelinisVsSedges)
wav2m.ufscar1.MelinisVsSedges$Scales <- wav2.ufscar1.MelinisVsSedges$Scales
wav2m.ufscar1.MelinisVsSedges$Total_Obs <- wav2.ufscar1.MelinisVsSedges$ScaleVar_Obs
wav2m.ufscar1.MelinisVsSedges$Total_CI_inf <- wav2.ufscar1.MelinisVsSedges$ScaleVar_CI[1,]
wav2m.ufscar1.MelinisVsSedges[c(4,6,8,10)] <- as.data.frame(wav2.ufscar1.MelinisVsSedges$ScaleSection_Obs)
wav2m.ufscar1.MelinisVsSedges[c(5,7,9,11)] <- as.data.frame(wav2.ufscar1.MelinisVsSedges$ScaleSection_CI)
write.table(wav2m.ufscar1.MelinisVsSedges, file="part4_wav2Scales_ufscar1_MelinisVsSedges.txt", quote=F, row.names=F, sep="\t")

wav2m.ufscar1.UrochloaVsMelinis <- matrix(NA,ncol=16, nrow=50)
colnames(wav2m.ufscar1.UrochloaVsMelinis) <- names.ufscar1.b
wav2m.ufscar1.UrochloaVsMelinis <- as.data.frame(wav2m.ufscar1.UrochloaVsMelinis)
wav2m.ufscar1.UrochloaVsMelinis <- as.data.frame(wav2m.ufscar1.UrochloaVsMelinis)
wav2m.ufscar1.UrochloaVsMelinis$Scales <- wav2.ufscar1.UrochloaVsMelinis$Scales
wav2m.ufscar1.UrochloaVsMelinis$Total_Obs <- wav2.ufscar1.UrochloaVsMelinis$ScaleVar_Obs
wav2m.ufscar1.UrochloaVsMelinis$Total_CI_inf <- wav2.ufscar1.UrochloaVsMelinis$ScaleVar_CI[,1]
wav2m.ufscar1.UrochloaVsMelinis$Total_CI_sup <- wav2.ufscar1.UrochloaVsMelinis$ScaleVar_CI[,2]
wav2m.ufscar1.UrochloaVsMelinis[c(5,8,11,14)] <- as.data.frame(wav2.ufscar1.UrochloaVsMelinis$ScaleSection_Obs)
wav2m.ufscar1.UrochloaVsMelinis[c(6,7,9,10,12,13,15,16)] <- as.data.frame(wav2.ufscar1.UrochloaVsMelinis$ScaleSection_CI)
write.table(wav2m.ufscar1.UrochloaVsMelinis, file="part4_wav2Scales_ufscar1_UrochloaVsMelinis.txt", quote=F, row.names=F, sep="\t")


# Supplementary figures
setwd("/home/pavel/Profissional/Pesquisa/MyPapers-Working/EIgrasses/figures")

scales2.itirapina <- wav2.itirapina.UrochloaVsGrasses$Scales
png(filename="SM4_scalesBivar_itirapina.png", height=35, width=30, unit="cm", res=300)
par(mfrow=c(5,5), mar=c(3,3,2,2), oma=c(3,3,1,1))
for(i in 1:5) {
	k <- c(2,6,10,12,14)[i]
  name.temp <- paste(names(wav2m.itirapina.UrochloaVsGrasses)[k],"UrochloaVsGrasses",sep="-")
  name.temp <- gsub("_Obs", "", name.temp)
  foo <- wav2m.itirapina.UrochloaVsGrasses[,k]
  bar <- wav2m.itirapina.UrochloaVsGrasses[,k+1]
  plot(foo ~ scales2.itirapina, type="l", main=name.temp, cex.main=1.0)
  lines(bar ~ scales2.itirapina, col="gray30", lty=2, lwd=2)

  name.temp <- paste(names(wav2m.itirapina.UrochloaVsSedges)[k],"UrochloaVsSedges",sep="-")
  name.temp <- gsub("_Obs", "", name.temp)
  foo <- wav2m.itirapina.UrochloaVsSedges[,k]
  bar <- wav2m.itirapina.UrochloaVsSedges[,k+1]
  plot(foo ~ scales2.itirapina, type="l", main=name.temp, cex.main=1.0)
  lines(bar ~ scales2.itirapina, col="gray30", lty=2, lwd=2)

  name.temp <- paste(names(wav2m.itirapina.MelinisVsGrasses)[k],"MelinisVsGrasses",sep="-")
  name.temp <- gsub("_Obs", "", name.temp)
  foo <- wav2m.itirapina.MelinisVsGrasses[,k]
  bar <- wav2m.itirapina.MelinisVsGrasses[,k+1]
  plot(foo ~ scales2.itirapina, type="l", main=name.temp, cex.main=1.0)
  lines(bar ~ scales2.itirapina, col="gray30", lty=2, lwd=2)

  name.temp <- paste(names(wav2m.itirapina.MelinisVsSedges)[k],"MelinisVsSedges",sep="-")
  name.temp <- gsub("_Obs", "", name.temp)
  foo <- wav2m.itirapina.MelinisVsSedges[,k]
  bar <- wav2m.itirapina.MelinisVsSedges[,k+1]
  plot(foo ~ scales2.itirapina, type="l", main=name.temp, cex.main=1.0)
  lines(bar ~ scales2.itirapina, col="gray30", lty=2, lwd=2)

	k <- c(2,8,14,17,20)[i]

  name.temp <- paste(names(wav2m.itirapina.UrochloaVsMelinis)[k],"UrochloaVsMelinis",sep="-")
  name.temp <- gsub("_Obs", "", name.temp)
  foo <- wav2m.itirapina.UrochloaVsMelinis[,k]
  bar <- wav2m.itirapina.UrochloaVsMelinis[,k+1]
  bar2 <- wav2m.itirapina.UrochloaVsMelinis[,k+2]
  plot(foo ~ scales2.itirapina, type="l", main=name.temp, cex.main=1.0)
  lines(bar ~ scales2.itirapina, col="gray30", lty=2, lwd=2)
  lines(bar2 ~ scales2.itirapina, col="gray30", lwd=2)
}
mtext(side=1, text="Scale (m)", outer=T)
mtext(side=2, text="Scale variance", outer=T)
dev.off()


scales2.ufscar2 <- wav2.ufscar2.UrochloaVsGrasses$Scales
png(filename="SM4_scalesBivar_ufscar2.png", height=35, width=30, unit="cm", res=300)
par(mfrow=c(5,5), mar=c(3,3,2,2), oma=c(3,3,1,1))
for(i in 1:5) {
	k <- c(2,6,8,10,12)[i]
  name.temp <- paste(names(wav2m.ufscar2.UrochloaVsGrasses)[k],"UrochloaVsGrasses",sep="-")
  name.temp <- gsub("_Obs", "", name.temp)
  foo <- wav2m.ufscar2.UrochloaVsGrasses[,k]
  bar <- wav2m.ufscar2.UrochloaVsGrasses[,k+1]
  plot(foo ~ scales2.ufscar2, type="l", main=name.temp, cex.main=1.0)
  lines(bar ~ scales2.ufscar2, col="gray30", lty=2, lwd=2)

  name.temp <- paste(names(wav2m.ufscar2.UrochloaVsSedges)[k],"UrochloaVsSedges",sep="-")
  name.temp <- gsub("_Obs", "", name.temp)
  foo <- wav2m.ufscar2.UrochloaVsSedges[,k]
  bar <- wav2m.ufscar2.UrochloaVsSedges[,k+1]
  plot(foo ~ scales2.ufscar2, type="l", main=name.temp, cex.main=1.0)
  lines(bar ~ scales2.ufscar2, col="gray30", lty=2, lwd=2)

  name.temp <- paste(names(wav2m.ufscar2.MelinisVsGrasses)[k],"MelinisVsGrasses",sep="-")
  name.temp <- gsub("_Obs", "", name.temp)
  foo <- wav2m.ufscar2.MelinisVsGrasses[,k]
  bar <- wav2m.ufscar2.MelinisVsGrasses[,k+1]
  plot(foo ~ scales2.ufscar2, type="l", main=name.temp, cex.main=1.0)
  lines(bar ~ scales2.ufscar2, col="gray30", lty=2, lwd=2)

  name.temp <- paste(names(wav2m.ufscar2.MelinisVsSedges)[k],"MelinisVsSedges",sep="-")
  name.temp <- gsub("_Obs", "", name.temp)
  foo <- wav2m.ufscar2.MelinisVsSedges[,k]
  bar <- wav2m.ufscar2.MelinisVsSedges[,k+1]
  plot(foo ~ scales2.ufscar2, type="l", main=name.temp, cex.main=1.0)
  lines(bar ~ scales2.ufscar2, col="gray30", lty=2, lwd=2)

	k <- c(2,8,11,14,17)[i]

  name.temp <- paste(names(wav2m.ufscar2.UrochloaVsMelinis)[k],"UrochloaVsMelinis",sep="-")
  name.temp <- gsub("_Obs", "", name.temp)
  foo <- wav2m.ufscar2.UrochloaVsMelinis[,k]
  bar <- wav2m.ufscar2.UrochloaVsMelinis[,k+1]
  bar2 <- wav2m.ufscar2.UrochloaVsMelinis[,k+2]
  plot(foo ~ scales2.ufscar2, type="l", main=name.temp, cex.main=1.0)
  lines(bar ~ scales2.ufscar2, col="gray30", lty=2, lwd=2)
  lines(bar2 ~ scales2.ufscar2, col="gray30", lwd=2)
}
mtext(side=1, text="Scale (m)", outer=T)
mtext(side=2, text="Scale variance", outer=T)
dev.off()



scales2.ufscar1 <- wav2.ufscar1.UrochloaVsGrasses$Scales
png(filename="SM4_scalesBivar_ufscar1.png", height=35, width=30, unit="cm", res=300)
par(mfrow=c(4,5), mar=c(3,3,2,2), oma=c(3,3,1,1))
for(i in 1:4) {
	k <- c(2,4,8,10)[i]
  name.temp <- paste(names(wav2m.ufscar1.UrochloaVsGrasses)[k],"UrochloaVsGrasses",sep="-")
  name.temp <- gsub("_Obs", "", name.temp)
  foo <- wav2m.ufscar1.UrochloaVsGrasses[,k]
  bar <- wav2m.ufscar1.UrochloaVsGrasses[,k+1]
  plot(foo ~ scales2.ufscar1, type="l", main=name.temp, cex.main=1.0)
  lines(bar ~ scales2.ufscar1, col="gray30", lty=2, lwd=2)

  name.temp <- paste(names(wav2m.ufscar1.UrochloaVsSedges)[k],"UrochloaVsSedges",sep="-")
  name.temp <- gsub("_Obs", "", name.temp)
  foo <- wav2m.ufscar1.UrochloaVsSedges[,k]
  bar <- wav2m.ufscar1.UrochloaVsSedges[,k+1]
  plot(foo ~ scales2.ufscar1, type="l", main=name.temp, cex.main=1.0)
  lines(bar ~ scales2.ufscar1, col="gray30", lty=2, lwd=2)

  name.temp <- paste(names(wav2m.ufscar1.MelinisVsGrasses)[k],"MelinisVsGrasses",sep="-")
  name.temp <- gsub("_Obs", "", name.temp)
  foo <- wav2m.ufscar1.MelinisVsGrasses[,k]
  bar <- wav2m.ufscar1.MelinisVsGrasses[,k+1]
  plot(foo ~ scales2.ufscar1, type="l", main=name.temp, cex.main=1.0)
  lines(bar ~ scales2.ufscar1, col="gray30", lty=2, lwd=2)

  name.temp <- paste(names(wav2m.ufscar1.MelinisVsSedges)[k],"MelinisVsSedges",sep="-")
  name.temp <- gsub("_Obs", "", name.temp)
  foo <- wav2m.ufscar1.MelinisVsSedges[,k]
  bar <- wav2m.ufscar1.MelinisVsSedges[,k+1]
  plot(foo ~ scales2.ufscar1, type="l", main=name.temp, cex.main=1.0)
  lines(bar ~ scales2.ufscar1, col="gray30", lty=2, lwd=2)

	k <- c(2,5,11,14)[i]

  name.temp <- paste(names(wav2m.ufscar1.UrochloaVsMelinis)[k],"UrochloaVsMelinis",sep="-")
  name.temp <- gsub("_Obs", "", name.temp)
  foo <- wav2m.ufscar1.UrochloaVsMelinis[,k]
  bar <- wav2m.ufscar1.UrochloaVsMelinis[,k+1]
  bar2 <- wav2m.ufscar1.UrochloaVsMelinis[,k+2]
  plot(foo ~ scales2.ufscar1, type="l", main=name.temp, cex.main=1.0)
  lines(bar ~ scales2.ufscar1, col="gray30", lty=2, lwd=2)
  lines(bar2 ~ scales2.ufscar1, col="gray30", lwd=2)
}
mtext(side=1, text="Scale (m)", outer=T)
mtext(side=2, text="Scale variance", outer=T)
dev.off()


### Part 5 - Check which univariate and bivariate scales are significant

setwd("/home/pavel/Profissional/Pesquisa/MyPapers-Working/EIgrasses/results")

files.univar <- list.files(pattern="part3")
files.bivar <- list.files(pattern="part4")

wav.univar <- list()
for(i in 1:length(files.univar)) {
  wav.univar[[i]] <- read.table(files.univar[i], sep="\t", header=T)
}
foo <- gsub(pattern=c("part3_wavScales_"), replacement=c(""), x=files.univar)
foo <- gsub(pattern=c(".txt"), replacement=c(""), x=foo)
names(wav.univar) <- foo

wav.bivar <- list()
for(i in 1:length(files.bivar)) {
  wav.bivar[[i]] <- read.table(files.bivar[i], sep="\t", header=T)
}
foo <- gsub(pattern=c("part4_wav2Scales_"), replacement=c(""), x=files.bivar)
foo <- gsub(pattern=c(".txt"), replacement=c(""), x=foo)
names(wav.bivar) <- foo

wav.bivar1 <- wav.bivar[-c(4,9,14)]
wav.bivar2 <- wav.bivar[c(4,9,14)]

scales <- wav.univar[[1]]$Scales
signif.univar <- list()
for(i in 1:length(wav.univar)) {
  signif.univar[[i]] <- list()
  foo <- wav.univar[[i]]
  Nuse <- (ncol(foo)-1)/2
  for(j in 1:Nuse) {
    k <- j*2
    test <- foo[,k] > foo[,k+1]
    bar <- scales[test]
    bar <- bar[!is.na(bar)]
    bar <- paste(bar, collapse=" ")
    signif.univar[[i]][[j]] <- bar
  }
  names(signif.univar[[i]]) <- names(foo)[(1:Nuse)*2]
}
names(signif.univar) <- names(wav.univar)

signif.bivar1 <- list()
for(i in 1:length(wav.bivar1)) {
  signif.bivar1[[i]] <- list()
  foo <- wav.bivar1[[i]]
  Nuse <- (ncol(foo)-1)/2
  for(j in 1:Nuse) {
    k <- j*2
    test <- foo[,k] < foo[,k+1]
    bar <- scales[test]
    bar <- bar[!is.na(bar)]
    bar <- paste(bar, collapse=" ")
    signif.bivar1[[i]][[j]] <- bar
  }
  names(signif.bivar1[[i]]) <- names(foo)[(1:Nuse)*2]
}
names(signif.bivar1) <- names(wav.bivar1)

signif.bivar2 <- list()
for(i in 1:length(wav.bivar2)) {
  signif.bivar2[[i]] <- list()
  foo <- wav.bivar2[[i]]
  Nuse <- (ncol(foo)-1)/3
  for(j in 1:Nuse) {
    k <- 2 + (j-1)*3
    test1 <- foo[,k] < foo[,k+1]
    test2 <- foo[,k] > foo[,k+2]
    bar1 <- scales[test1]
    bar1 <- bar1[!is.na(bar1)]
    bar1 <- paste(bar1, collapse=" ")
    bar2 <- scales[test2]
    bar2 <- bar2[!is.na(bar2)]
    bar2 <- paste(bar2, collapse=" ")
    signif.bivar2[[i]][[j]] <- paste(bar1, bar2, sep=" ; ")
  }
  names(signif.bivar2[[i]]) <- names(foo)[2 + ((1:Nuse)-1)*3]
}
names(signif.bivar2) <- names(wav.bivar2)

























