replace.diameters <- function(x) {
  foo <- x
  foo[x==1] <- 0.03
  foo[x==2] <- 1.55
  foo[x==3] <- 3.75
  foo[x==4] <- 10
  foo[x==5] <- 22.5
  return(foo)
}

shannon <- function(x) {
  foo <- sum(x)
  Pi <- x/foo
  result <- -sum(Pi * log(Pi))
  return(result)
}

StrD.calc <- function(x, w) { # x is a matrix
  dists <- daisy(x, metric="gower", weights=w)
  clust <- agnes(dists, diss=1, method="average", keep.diss=F, keep.data=F)
  result <- treeheight(as.hclust(clust))
  return(result)
}


sum.cover <- function(x) {
	Nx <- length(x)
	if (Nx == 1) result <- x
	if (Nx == 2) {
		if (x[1] == x[2]) result <- x[1]+1 else result <- x[1]
	}
	if (Nx == 3) {
		foo <- sort(x, decreasing=T)
		if (any (foo[-1] == foo[1])) result <- foo[1]+1 else result <- foo[1]
	}
	if(Nx == 4) {
		foo <- sort(x, decreasing=T)
		if(all(foo[-1] == foo[1])) {
			result <- foo[1]+2
			} else {
				if(any(foo[-1] == foo[1])) result <- foo[1]+1 else result <- foo[1]
		}
	}
	if(result>5) return(5) else return(result)
}



sim.CSR <- function(x, Nsim=5000, print.loop=F) {
  data.sim <- matrix(nrow=length(x), ncol=Nsim)
  data.sim[,1] <- x
  for(i in 2:Nsim) {
    data.sim[,i] <- sample(x)
    if(print.loop) if(i %% 500 == 0) print(i)
  }
  return(data.sim)
}


sim.AR1 <- function(x, Nsim=5000, print.loop=T) {
  xprev <- c(NA, x[1:(length(x)-1)])
  xnext <- c(x[2:length(x)],NA)
  mod1 <- lm(x ~ xprev)
  mod2 <- lm(x ~ xnext)
  slope <- mean(c(coef(mod1)[2],coef(mod2)[2]))
  intercept <- mean(c(coef(mod1)[1],coef(mod2)[1]))
  error <- mean(c(sd(resid(mod1)),sd(resid(mod2))))
  data.sim <- matrix(nrow=length(x), ncol=Nsim)
  data.sim[,1] <- x
  for(i in 2:Nsim) {
    sim.now <- numeric(length(x))
    k.start <- sample(1:length(x),1)
    sim.now[k.start] <- x[k.start]
    if(k.start > 1) {
      for(k in (k.start-1):1) {
        x.foo <- sim.now[k+1]
        sim.now[k] <- rnorm(1, mean=intercept + slope*x.foo, sd=error)
        if(sim.now[k] < 0) sim.now[k] <- 0
        if(sim.now[k] > max(x)) sim.now[k] <- max(x)
      }
    }
    if(k.start < length(x)) {
      for(k in (k.start+1) : length(x)) {
        x.foo <- sim.now[k-1]
        sim.now[k] <- rnorm(1, mean=intercept + slope*x, sd=error)
        if(sim.now[k] < 0) sim.now[k] <- 0
        if(sim.now[k] > max(x)) sim.now[k] <- max(x)
      }
    }
    data.sim[,i] <- sim.now
    if(print.loop) if(i %% 500 == 0) print(i)
  }
  return(data.sim)
}


sim.MC1 <- function(x, Nperm=5000, keep.obs=T, print.loop=T, which.loop=100) {
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



wavCWTvarCIs <- function(x, significance=T, sections=NULL, make.plot=F, zero.padding=T, remove.COI=T, test=F, scale.max=75, scale.var=T, pos.var=T, CI=T, CIquant=c(0.025, 0.975), effect.size=T,  keep.all=F, print.loop=T, which.loop=1, save.cwt=F, file.title="my_beautiful_wavelet", wav.template="gaussian2", scale.per.section=T) {
  #This function does not save the original wavelet output, retaining only scale and position variance. Such functionality may be added in the future. For now, the keep.all, save.cwt and file.title arguments are not fully functional.
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


SVD <- function(mat1,mat2) { ###Function by T. Rouyer
  R <- mat1%*%t(mat2)
  e1 <- svd(R)
  coef1<-e1$d
  square <- coef1*coef1
  sq <- square/sum(square)
  ### LEADINGS MODES
  vect1 <- e1$u
  vect2 <- e1$v
  ### PROJECTION: LEADINGS PATTERNS
  A <- (t(vect1)%*%mat1)
  B <- (t(vect2)%*%mat2)
  return(list(A,B,vect1,vect2,sq))
}

distance5 <- function(A,B) {###Function by T. Rouyer
  L <- length(A)
  d <- abs((A[1:(L-1)]-B[1:(L-1)])-(A[2:L]-B[2:L]))
  return(d)
}

stand <- function(v) {###Function by T. Rouyer
  d <- v+abs(min(v))
  d <- d/max(d)
  return(d)
}

M_cor4 <- function(PW2,seuil) {###Function by T. Rouyer
  T <- matrix(NA,length(PW2),length(PW2))
  F <- matrix(NA,length(PW2),length(PW2))
  SQ <- matrix(NA,length(PW2),length(PW2))
  N <- matrix(NA,length(PW2),length(PW2))
  for (j in 1:(length(PW2)-1)) {
    for (i in (j+1):length(PW2)) {
      mat1 <- PW2[[i]]
      mat2 <- PW2[[j]]

      x <- mat1
      y <- mat2
      ### MCA TEMPORELLE
      S <- SVD(x,y)
      A <- S[[1]];B <- S[[2]];vect1 <- S[[3]];vect2 <- S[[4]]
      sq <- S[[5]]
      nbaxes <- which(cumsum(sq)>=seuil)[1]
      ### LEADINGS PATTERNS STANDARDISES
      lpat <- 0;v <- 0
      for (n in 1:nbaxes) {
        Anorm <- stand(A[n,])
        Bnorm <- stand(B[n,])
        Vect1norm <- stand(vect1[,n])
        Vect2norm <- stand(vect2[,n])
        lpat <- lpat+sq[n]*sum(atan(distance5(Anorm,Bnorm)))
        v <- v+sq[n]*sum(atan(distance5(Vect1norm,Vect2norm)))
      }
      lpat <- lpat/sum(sq[1:nbaxes]);v <- v/sum(sq[1:nbaxes])
      ###
      
      T[i,j] <- lpat
      F[i,j] <- v
      SQ[i,j] <- sum(sq[1:nbaxes])
      N[i,j] <- nbaxes
      print(c(j,i))     
    }

 
  } 
  return(list(T,F,SQ,N))
} 


wavCWT2 <- function(x, zero.padding=T, remove.COI=T, scale.max=75, wav.template="gaussian2") {
	data.foo=x
	full.length=length(data.foo)
	half.length=round(length(data.foo)/2)
	if(zero.padding) data.foo=c(rep(0,half.length),data.foo,rep(0,half.length/2))
	wav.orig=wavCWT(data.foo,scale.range=c(1,scale.max),wavelet="gaussian2")
	scales=attr(wav.orig,"scale")
	
	if(zero.padding) {  
		wav.scale=attr(wav.orig,"scale")
		wav.time=1:full.length
		wav.wavelet=attr(wav.orig,"wavelet")
		wav.series=attr(wav.orig,"series")
		wav.sampling.interval=attr(wav.orig,"sampling.interval")
		wav.series.name=attr(wav.orig,"series.name")
		wav.n.sample=length(wav.time)
		wav.n.scale=attr(wav.orig,"n.scale")
		wav.filter.arg=attr(wav.orig,"filter.arg")
		wav.orig=wav.orig[(half.length+1):(half.length+full.length),]
		 
		class(wav.orig)="wavCWT"
		attr(wav.orig,"scale")=wav.scale
		attr(wav.orig,"time")=wav.time
		attr(wav.orig,"wavelet")=wav.wavelet
		attr(wav.orig,"series")=wav.series
		attr(wav.orig,"sampling.interval")=wav.sampling.interval
		attr(wav.orig,"series.names")=wav.series.name
		attr(wav.orig,"n.sample")=wav.n.sample
		attr(wav.orig,"n.scale")=wav.n.scale
		attr(wav.orig,"filter.arg")=wav.filter.arg
	}
	if(remove.COI) {
		for(k in 1:length(scales)) {
			index1=1:(scales[k]*2)
			index2=nrow(wav.orig)-index1+1
			indices=c(index1,index2)
			wav.orig[indices,k]=NA
		}
	}
	return(wav.orig)
}

