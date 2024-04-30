require(zoo) #format de serie temporelle pratique et facile d'utilisation (mais plus volumineux)
require(tseries) #diverses fonctions sur les series temporelles
require(stats)

path <- ""
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

#Ok donc on a une certaine liste de modèles qui marchent, je la copie là : 
#p <- c(1, 2, 5, 3, 4, 5, 3, 6, 2, 5, 6, 0, 2, 5, 5, 6, 0, 1, 2, 4, 5, 6, 5, 4, 5, 6, 0, 1, 3)
#q <- c(12, 13, 13, 14, 15, 19, 21, 22, 23, 23, 24, 25, 25, 25, 27, 27, 28, 28, 28, 28, 28, 28, 29, 30, 30, 30, 31, 31, 31)
#arsignif <- c(1, NA, 1, 1, NA, NA, 1, 1, 1, 1, NA, NA, 1, 1, 1, 1, NA, 1, 1, 1, 1, NA, 1, 1, 1, 1, NA, 1, 1)
#masignif <- c(1, NA, 1, 1, NA, NA, 1, 1, 1, 1, NA, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, NA, 1, 1, 1, 1, 1, 1, 1)
#resnocorr <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#ok <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#selec2 <- data.frame(p, q, arsignif, masignif, resnocorr, ok)

pqs <- apply(selec,1,function(row) list("p"=as.numeric(row[1]),"q"=as.numeric(row[2])))
names(pqs) <- paste0("arma(",selec[,1],",",selec[,2],")")
models <- lapply(pqs, function(pq) arima(r,c(pq[["p"]],0,pq[["q"]])))
vapply(models, FUN.VALUE=numeric(2), function(m) c("AIC"=AIC(m),"BIC"=BIC(m)))

#De même, liste de AIC et BIC
#modeles <- c("arma(1,12)", "arma(2,13)", "arma(5,13)", "arma(3,14)", "arma(4,15)", "arma(5,19)", "arma(3,21)", "arma(6,22)", 
#"arma(2,23)", "arma(5,23)", "arma(6,24)", "arma(0,25)", "arma(2,25)", "arma(5,25)", "arma(5,27)", "arma(6,27)",
#"arma(0,28)", "arma(1,28)", "arma(2,28)", "arma(4,28)", "arma(5,28)", "arma(6,28)", "arma(5,29)", "arma(4,30)",
#"arma(5,30)", "arma(6,30)", "arma(0,31)", "arma(1,31)", "arma(3,31)")
#AIC <- c(3670.189, 3674.189, 3676.161, 3661.907, 3665.845, 3674.263, 3680.553, 3681.299, 
#3684.639, 3675.041, 3681.351, 3686.902, 3688.032, 3690.622, 3679.842, 3681.371,
#3688.092, 3688.123, 3679.831, 3685.073, 3683.305, 3678.110, 3680.398, 3677.711,
#3681.543, 3671.954, 3679.818, 3679.391, 3683.972)
#BIC <- c(3730.431, 3742.464, 3756.484, 3738.214, 3750.185, 3778.683, 3784.973, 3801.784,
#3793.075, 3795.526, 3809.868, 3795.338, 3804.500, 3819.139, 3816.392, 3821.937,
#3808.577, 3812.624, 3808.349, 3821.623, 3823.870, 3822.692, 3824.979, 3822.292,
#3830.141, 3824.567, 3812.351, 3815.941, 3828.554)
#tableau_aic_bic <- data.frame(modeles, AIC, BIC)

#On minimise l'AIC avec arma(3,14) et le BIC avec arma(1,12)

selected_indices <- grep("arma\\(1,12\\)|arma\\(3,14\\)", names(models))
selected_models <- models[selected_indices]

library(ggplot2)

# Prédire les deux prochaines valeurs des résidus
pred_1_12 <- predict(selected_models[["arma(1,12)"]], n.ahead = 2)
pred_3_14 <- predict(selected_models[["arma(3,14)"]], n.ahead = 2)

# Prédire la tendance pour les prochaines valeurs
future_trend <- predict(lt, newdata = data.frame(trend2 = (length(xm) + 1:2)^2))

# Prédictions pour la série d'origine (xm) en ajoutant les résidus prédits à la tendance prédite
pred_xm_1_12 <- future_trend + pred_1_12$pred
pred_xm_3_14 <- future_trend + pred_3_14$pred

# Calculer les intervalles de confiance pour xm
ci_lower_1_12 <- pred_xm_1_12 - 1.96 * pred_1_12$se
ci_upper_1_12 <- pred_xm_1_12 + 1.96 * pred_1_12$se
ci_lower_3_14 <- pred_xm_3_14 - 1.96 * pred_3_14$se
ci_upper_3_14 <- pred_xm_3_14 + 1.96 * pred_3_14$se

# Créer un data frame pour le tracé
time <- c(1:length(xm), length(xm) + 1:length(pred_xm_1_12))
data <- data.frame(
  time = time,
  values_1_12 = c(xm, pred_xm_1_12),
  lower_ci_1_12 = c(rep(NA, length(xm)), ci_lower_1_12),
  upper_ci_1_12 = c(rep(NA, length(xm)), ci_upper_1_12),
  values_3_14 = c(xm, pred_xm_3_14),
  lower_ci_3_14 = c(rep(NA, length(xm)), ci_lower_3_14),
  upper_ci_3_14 = c(rep(NA, length(xm)), ci_upper_3_14)
)

# Tracer la série d'origine avec intervalles de confiance
p_1_12 <- ggplot(data, aes(x = time, y = values_1_12)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci_1_12, ymax = upper_ci_1_12), alpha = 0.2) +
  labs(title = "Prédictions avec Région de Confiance à 95% (ARMA(1,12))",
       x = "Temps",
       y = "Série Originale")

print(p_1_12)

# Tracer la série d'origine avec le modèle ARMA(2,14)
p_3_14 <- ggplot(data, aes(x = time, y = values_3_14)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci_3_14, ymax = upper_ci_3_14), alpha = 0.2) +
  labs(title = "Prédictions avec Région de Confiance à 95% (ARMA(3,14))",
       x = "Temps",
       y = "Série Originale")

print(p_3_14)





