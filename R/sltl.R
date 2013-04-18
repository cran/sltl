

sltl <-function(x,c.window,c.degree=1,Fm=6,Fs=0,searchQ=FALSE,ic=c("BIC","AIC"),
			type=c("daily","monthly"), filling=FALSE, npass=1)
{
	original <- x
	nx <- length(x)
	xt <- 1:nx

    if((type=="daily")&(filling==TRUE)) x <- filling.d(x)
    if((type=="monthly")&(filling==TRUE)) x <- filling.m(x)

    if (is.timeSeries(x)&&frequency(x)==12) type <- "monthly" else type <- match.arg(type)
    if (type=="daily") s <- 365.25 else s <- 12
    if (! (s%in%c(12,365.25))) 
	stop("Error: s must be 12 or 365.25")
    if (min(Fm,Fs) < 0) 
	stop("Error: both Fm and Fs must be >= 0")
    if (max(Fm,Fs) > 6) 
	stop("Error: Fm and Fs are restricted to be <= 6")
    if (max(Fm,Fs) > s/2) 
	stop("error: Fm or Fs setting too large")

	ic <- match.arg(ic)
	if (length(x) < s) stop("error: need at least one year (366 days or 12 months)")

if(is.timeSeries(x)){
	t <- time(x)
	x <- as.numeric(x)
	}
dt <- function(x,c.window,c.degree){
	loess(x~c(1:length(x)),degree=c.degree,span=c.window/length(x))$fitted
	}
dss <- function(x,Fm,Fs,s,ic){
 		 if (Fm == 0) {
 		   zds <- x-mean(x)
 		   } else if (Fm > 0) {
 		   X <- matrix(rep(1,(2*Fm+1)*nx), ncol=2*Fm+1)
  		  jj <- 2
  		  for (j in 1:Fm){
  		    xx <- j*(2*pi/s)*xt
   		   X[,jj] <- sin(xx)
   		   X[,jj+1] <- cos(xx)
  		    jj <- jj+2
  		    }
  		  zds <- resid(lm.fit(x=X, y=x))
		}
zds
}
trd1 <- dt(x,c.window=c.window[1],c.degree=c.degree[1])
if(searchQ){
      maxFm <- Fm
      maxFs <- Fs
      m <- matrix(numeric(4*((maxFm+1)*(maxFs+1))), ncol=4)
      colnames(m) <- c("Fm", "Fs", "p", ic)
      i <- 0
      for (iFm in 0:maxFm){
        for (iFs in 0:maxFs) {
          i <- i+1
 		 ICBest <- AdmissPenalty <- Inf
		  pHat <- zds <- NA
		zds <- dss(x,Fm=iFm,Fs=iFs,s=s,ic=ic)
		AdmissPenalty <- 0
		ans <- SelectModel(zds, lag.max = 20, Criterion = ic, Best=2)
		pHat <- ans[1,1]
   		ICBest <- ans[1,2]
 		 if (ic == "BIC") parPenalty <- log(nx)*2*(iFm+iFs) else
  	     parPenalty <- 4*(iFm+iFs) 
 		 BestIC <- ICBest + parPenalty + AdmissPenalty
 		 out <- c(iFm, iFs, pHat, BestIC)
 		 if (ic == "BIC") names(out) <- c("Fm", "Fs", "p", "BIC") else
 		       names(out) <- c("Fm", "Fs", "p", "AIC")
  		m[i,] <- out
 	   }
	}
      rownames(m) <- rep(" ", nrow(m))
      ind <- which.min(m[,4])
      rownames(m)[ind] <- "*" 
      FmOpt <- m[ind,1]
      FsOpt <- m[ind,2]
	Fms <- c(FmOpt,FsOpt)
	} else if(!searchQ) Fms <- c(Fm,Fs)
if(class(original)=="timeSeries" & frequency(original)==12) s <- 12 else s <- 365.25
irr1 <- dss(x-trd1,Fm=Fms[1],Fs=Fms[2],s=s,ic=ic)
ssn1 <- x - trd1 - irr1

E <- 1
maxit <- 0
trd2 <- trd1
ssn2 <- ssn1
irr2 <- irr1
while(E >= 10^-6 & maxit<5)
	{
	maxit <- maxit + 1
	trd1 <- trd2
	ssn1 <- ssn2
	irr1 <- irr2
	trd2 <- dt(x-ssn1,c.window=c.window[1],c.degree=c.degree[1])
	irr2 <- dss(x-trd2,Fm=Fms[1],Fs=Fms[2],s=s,ic=ic)
	ssn2 <- x - trd2 - irr2
	E <- (sum((irr1-irr2)^2))/(sum(irr2^2))
	} 
if(maxit==5) stop("Warning: MAXIT is > 5")

n.lt <- length(c.window)
if(n.lt==1) 
{
comps <- 
if(filling) timeSeries(cbind(data=original,filled_data=x,seasonal=ssn2,trend=trd2,remainder=irr2),charvec=t)
else timeSeries(cbind(data=x,seasonal=ssn2,trend=trd2,remainder=irr2),charvec=t)
}
if(n.lt>1)
{
lt.comp <- matrix(0,nrow=length(x),ncol=n.lt)
colnames(lt.comp) <- paste("long_term", 1:n.lt, sep="")
for(j in 1:npass)
{
	for(k in 1:n.lt)
	{
	temp <- x - ssn2 - if(n.lt==2) lt.comp[,-k] else apply(lt.comp[,-k],1,sum)
	model <- loess(temp~c(1:length(x)), degree=c.degree[k], span=c.window[k]/length(x))
	lt.comp[,k] <- model$fitted
	}
	irr <- model$resid
}
if(filling) comps <- timeSeries(cbind(data=original,filled_data=x,seasonal=ssn2,lt.comp,remainder=irr),charvec=t)
if(!filling) comps <- timeSeries(cbind(data=x,seasonal=ssn2,lt.comp,remainder=irr),charvec=t)
}

res <- list(call=match.call(),time.series=comps,Fms =
as.numeric(c(Fm=Fms[1],Fs=Fms[2])),window = c.window, degree=c.degree)
class(res) <- "sltl"
res
}

plot.sltl <- function(x, labels = colnames(X),
		     set.pars = list(mar = c(0, 6, 0, 6), oma = c(6, 0, 4, 0),
				     tck = -0.01, mfrow = c(nplot, 1)),
		     main = NULL, ...)
{
    X <- x$time.series
    nplot <- ncol(X)
    range.bars = TRUE
    col.range = "light gray"
    if(range.bars)
	mx <- min(apply(rx <- apply(X,2, range,na.rm=TRUE), 2, diff,na.rm=TRUE))
    dev.hold(); on.exit(dev.flush())
    if(length(set.pars)) {
	oldpar <- do.call("par", as.list(names(set.pars)))
	on.exit(par(oldpar), add = TRUE)
	do.call("par", set.pars)
    }
    for(i in 1L:nplot) {
	plot(X[, i], type = if(i < nplot) "l" else "h",
	     xlab = "", ylab = "", axes = FALSE, xaxt = "n") #, ...)
	if(range.bars) {
	    dx <- 1/64 * diff(ux <- par("usr")[1L:2])
	    y <- mean(rx[,i])
	    rect(ux[2L] - dx, y + mx/2, ux[2L] - 0.4*dx, y - mx/2,
		 col = col.range, xpd = TRUE)
	}
	if(i == 1 && !is.null(main))
	    title(main, line = 2, outer = par("oma")[3L] > 0)
	if(i == nplot) abline(h=0)
	box()
	right <- i %% 2 == 0
	axis(2, labels = !right)
	axis(4, labels = right)
	axis(1, labels = i == nplot)
	mtext(labels[i], side = 2, 3)
    }
    mtext("time", side = 1, line = 3)
    invisible()
}

summary.sltl <- function(object, ...){
	cat(" Call:\n ")
	dput(object$call, control=NULL)
	cat("\n Time Series Components: \n")
	print(summary(object$time.series[,-1]))
	cat("\n Other Components:\n")
	cat("window : ")
	print(object$window)
	cat("degree : ")
	print(object$degree)
	cat("Fm : ")
	print(object$Fms[1])
	cat("Fs : ")
	print(object$Fms[2])
	invisible(object)
	}



filling.m <- function(x)          
{
c1 <- median(x,na.rm=TRUE)
Sm.hat <- numeric(12)
Sm.hat[1:12] <- NA
for(i in 1:12){
	a <- x[month(x)==i]
	Sm.hat[i] <- median(a-c1,na.rm=TRUE)
	}
Sm <- Sm.hat - mean(Sm.hat,na.rm=TRUE)
irr <- numeric(length(x))
for(i in 1:12){
	ind <- which(month(x)==i)
	irr[ind] <- x[ind]-c1-Sm[i]
	}
qt <- quantile(irr,na.rm=TRUE)[c(2,4)]
step <- 1.5 * (qt[2] - qt[1])
far.out1 <- which(irr < as.vector(qt[1]-2*step))
far.out2 <- which(irr > as.vector(qt[2]+2*step))
if(length(far.out1)>0) irr[far.out1] <- as.vector(qt[1]-2*step)
if(length(far.out2)>0) irr[far.out2] <- as.vector(qt[2]+2*step)

desnl <- c1 + irr

year <- year(x)
c2 <- numeric(max(year)-min(year)+1)
for(i in min(year):max(year)){
	j <- i - min(year) + 1
	index <- which((year)==i)
	if(sum(x[index],na.rm=TRUE)==0)
	c2[j] <- mean(desnl[year==i-1|year==i+1],na.rm=TRUE) else
	c2[j] <- mean(desnl[year==i],na.rm=TRUE)
	}

Sm2 <- numeric(12)
Sm2[1:12] <- NA
z.c2 <- x - c2[year-min(year)+1]
for(m in 1:12){
	mth <- z.c2[month(z.c2)==m]
	Sm2[m] <- median(mth,na.rm=TRUE)
	}
Sm2 <- Sm2 - mean(Sm2,na.rm=TRUE)

irr <- z.c2 - Sm2[month(x)]

qt <- quantile(irr,na.rm=TRUE)[c(2,4)]
step <- 1.5 * (qt[2] - qt[1])
far.out1 <- which(irr < as.vector(qt[1]-2*step))
far.out2 <- which(irr > as.vector(qt[2]+2*step))
if(length(far.out1)>0) irr[far.out1] <- as.vector(qt[1]-2*step)
if(length(far.out2)>0) irr[far.out2] <- as.vector(qt[2]+2*step)

c3 <- c2 + mean(irr,na.rm=TRUE)

ind <- which(is.na(x))
out <- x
out[ind] <- c3[year[ind]-min(year)+1] + Sm2[month(out)[ind]]
out
}


filling.d <- function(x){
c1 <- median(x,na.rm=TRUE)
time.dummy <- timeDate("1972-01-01",format="%Y-%m-%d")
Sm.hat <- numeric(366)
Sm.hat[1:366] <- NA
for(i in 1:366){
	ind <- which(month(x)==month(time.dummy) & day(x)==day(time.dummy))
	if(!all(is.na(x[ind]))) Sm.hat[i] <- median(x[ind]-c1,na.rm=TRUE)
	time.dummy <- time.dummy + days(1)
	}
Sm <- Sm.hat - mean(Sm.hat,na.rm=TRUE)
tt <- time(x)
year(tt) <- year(tt)-1
month(tt) <- 12
day(tt) <- 31
irr <- x - c1 - Sm[time(x)-tt]

qt <- quantile(irr,na.rm=TRUE)[c(2,4)]
step <- 1.5 * (qt[2] - qt[1])
far.out1 <- which(irr < as.vector(qt[1]-2*step))
far.out2 <- which(irr > as.vector(qt[2]+2*step))
if(length(far.out1)>0) irr[far.out1] <- as.vector(qt[1]-2*step)
if(length(far.out2)>0) irr[far.out2] <- as.vector(qt[2]+2*step)

desnl <- c1 + irr

M <- max(year(x))
m <- min(year(x))
c2 <- numeric(M-m+1)
for(i in m:M){
	j <- i - m + 1
	ind <- which(year(x)==i)
	if(all(is.na(x[ind])))	c2[j] <- mean(desnl[year(x)==i-1|year(x)==i+1],na.rm=TRUE) 
	else c2[j] <- mean(desnl[year(x)==i],na.rm=TRUE)
}

Sm2 <- numeric(366)
Sm2[1:366] <- NA
Sm_2 <- x - c2[year(x)-m+1]
time.dummy <- timeDate("1972-01-01",format="%Y-%m-%d")
for(i in 1:366){
	ind <- which(month(x)==month(time.dummy) & day(x)==day(time.dummy))
	if(!all(is.na(x[ind]))) Sm2[i] <- median(Sm_2[ind],na.rm=TRUE)
	time.dummy <- time.dummy + days(1)
	}
Sm2 <- Sm2 - mean(Sm2,na.rm=TRUE)

irr <- Sm_2 - Sm2[time(x)-tt]

qt <- quantile(irr,na.rm=TRUE)[c(2,4)]
step <- 1.5 * (qt[2] - qt[1])
far.out1 <- which(irr < as.vector(qt[1]-2*step))
far.out2 <- which(irr > as.vector(qt[2]+2*step))
if(length(far.out1)>0) irr[far.out1] <- as.vector(qt[1]-2*step)
if(length(far.out2)>0) irr[far.out2] <- as.vector(qt[2]+2*step)

c3 <- c2 + mean(irr,na.rm=TRUE)

ind <- which(is.na(x))
out <- x
out[ind] <- c3[year(x)[ind]-m+1] + Sm2[time(x)[ind]-tt[ind]]
out 
}
