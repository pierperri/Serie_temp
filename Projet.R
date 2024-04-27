require(zoo) #format de serie temporelle pratique et facile d'utilisation (mais plus volumineux)
require(tseries) #diverses fonctions sur les series temporelles
require(stats)

path <- "C:\\Users\\pierr\\OneDrive\\Documents\\ENSAE\\2A - S2\\Serie temp\\Projet"
setwd(path)

datafile <- "valeurs_mensuelles.csv"
data <- read.csv(datafile,sep=";")
T <- length(data)
xm <- as.numeric(zoo(data[-(1:3), 2]))

plot(xm)

#On enlève la tendance qui parait quadratique

trend <- 1:length(xm)
trend2 <- trend^2
lt <- lm(xm ~ trend2) #
summary(lt) #
r <- lt$residuals #
par(mfrow=c(1,2))
plot(r)

#On teste UR 

pp.test(r)

#On peut rejeter l'hypothèse d'une UR

x <- r

acf(x, 50)
pacf(x)

#On peut dire que le p max (pour AR) semble à 6 et q max (pour MA) à 31

pmax = 6; qmax = 31

#On utilise la fonction qui teste tous les p et q en dessous de pmax et qmax

Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}

signif <- function(estim){ #fonction de test des significations individuelles des coefficients
  coef <- estim$coef
  se <- sqrt(diag(estim$var.coef))
  t <- coef/se
  pval <- (1-pnorm(abs(t)))*2
  return(rbind(coef,se,pval))
}

modelchoice <- function(p,q,data=x, k=24){
  estim <- try(arima(data, c(p,0,q),optim.control=list(maxit=20000)))
  if (class(estim)=="try-error") return(c("p"=p,"q"=q,"arsignif"=NA,"masignif"=NA,"resnocorr"=NA, "ok"=NA))
  arsignif <- if (p==0) NA else signif(estim)[3,p]<=0.05
  masignif <- if (q==0) NA else signif(estim)[3,p+q]<=0.05
  resnocorr <- sum(Qtests(estim$residuals,24,length(estim$coef)-1)[,2]<=0.05,na.rm=T)==0
  checks <- c(arsignif,masignif,resnocorr)
  ok <- as.numeric(sum(checks,na.rm=T)==(3-sum(is.na(checks))))
  return(c("p"=p,"q"=q,"arsignif"=arsignif,"masignif"=masignif,"resnocorr"=resnocorr,"ok"=ok))
}

## fonction pour estimer et verifier tous les arima(p,q) avec p<=pmax et q<=max
armamodelchoice <- function(pmax,qmax){
  pqs <- expand.grid(0:pmax,0:qmax)
  t(apply(matrix(1:dim(pqs)[1]),1,function(row) {
    p <- pqs[row,1]; q <- pqs[row,2]
    cat(paste0("Computing ARMA(",p,",",q,") \n"))
    modelchoice(p,q)
  }))
}

armamodels <- armamodelchoice(pmax,qmax)

selec <- armamodels[armamodels[,"ok"]==1&!is.na(armamodels[,"ok"]),] #modeles bien ajustes et valides
selec

pqs <- apply(selec,1,function(row) list("p"=as.numeric(row[1]),"q"=as.numeric(row[2])))
names(pqs) <- paste0("arma(",selec[,1],",",selec[,2],")")
models <- lapply(pqs, function(pq) arima(r,c(pq[["p"]],0,pq[["q"]])))
vapply(models, FUN.VALUE=numeric(2), function(m) c("AIC"=AIC(m),"BIC"=BIC(m)))
