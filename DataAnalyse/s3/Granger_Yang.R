
##################################################################################
set.seed(42) 
borehole <- function(x)
{
  ##########################################################################
  #
  # BOREHOLE FUNCTION
  #
  ##########################################################################
  #
  # OUTPUT AND INPUT:
  #
  # y  = water flow rate
  # x = c(rw, riw, r, Tu, Hu, Tum, Hum, Tlm, Hlm, Tl, Hl, L, Kw)
  #
  ##########################################################################
  
  xx <- matrix(x, ncol=13)
  
  rw  <- xx[,1]
  riw <- xx[,2]
  r   <- xx[,3]
  Tu  <- xx[,4]
  Hu  <- xx[,5]
  Tum <- xx[,6]
  Hum <- xx[,7]
  Tlm <- xx[,8]
  Hlm <- xx[,9]
  Tl  <- xx[,10]
  Hl  <- xx[,11]
  L   <- xx[,12]
  Kw  <- xx[,13]
  
  frac1 <- 2 * pi * Tu * (Hu-Hl)
  frac11 <- (Tum - Tlm) * (Hum - Hlm)
  frac2 <- frac1 / frac11
  
  frac2a <- 2*L*Tu / (log(r/rw)*rw^2*Kw)
  frac2b <- Tu / Tl
  frac2 <- log(r/rw) * (1+frac2a+frac2b)
  
  y <- frac1 / frac2
  return(y)
}

##################################################################################

EchantBorehole <- function(N){
  
  # Description de la fonction :
  # Cette fonction genere un échantillon (pour le modele borhole) de taille N.
  #
  # Ces 13 variables sont supposées statistiquement independantes
  #
  # Entrées de la fonction :
  # - taille = taille de l'échantillon
  #
  # Description des réponses : 
  # - X = matrice de taille N x 13
  
  X = matrix(NA, N, 13)
  
  X[,1] <- rnorm(N, 0.1, 0.015)
  X[,2] <- rnorm(N, 0.05, 0.01)
  X[,3] <- rlnorm(N, 7.71, 1.0056)
  X[,4] <- runif(N, 63100, 116000)
  X[,5] <- runif(N, 1000, 1100)
  X[,6] <- runif(N, 6310, 11600)
  X[,7] <- runif(N, 900, 1000)
  X[,8] <- runif(N, 631, 1160)
  X[,9] <- runif(N, 800, 900)
  X[,10] <- runif(N, 63.1, 116)
  X[,11] <- runif(N, 700, 800)
  X[,12] <- runif(N, 1120, 1680)
  X[,13] <- runif(N, 3000, 12000)
  
  colnames(X) <- c("rw","riw","r","Tu","Hu","Tum","Hum","Tlm","Hlm","Tl","Hl","L","Kw") # noms de variables aleatoires
  
  return(X)
  
}#end function


#1.

#Valeurs minimales
val_min <- c(0.1 - 3.0*0.015, 0.05 - 3.0*0.01, exp(7.71 - 3.0*1.0056), 
          63100, 1000, 6310, 900, 631, 800, 63.1, 700, 1120, 3000)

#Valeurs maximales
val_max <- c(0.1 + 3.0*0.015, 0.05 + 3.0*0.01, exp(7.71 + 3.0*1.0056), 
          116000, 1100, 11600, 1000, 1160, 900, 116, 800, 1680, 12000)

#Valeurs moyennes
val_moy <- c(0.1, 0.05, exp(7.71 + (1.0056^2) / 2), 
          mean(c(63100, 116000)), mean(c(1000, 1100)), 
          mean(c(6310, 11600)), mean(c(900, 1000)), 
          mean(c(631, 1160)), mean(c(800, 900)), 
          mean(c(63.1, 116)), mean(c(700, 800)), 
          mean(c(1120, 1680)), mean(c(3000, 12000)))


#################################

#2.
#a.
echantillon <- EchantBorehole(1000)
resultats <- apply(echantillon, 1, function(x) borehole(as.numeric(x)))

moyEchant <- mean(resultats)
varEchant <- var(resultats)

cat("Moyenne : ",moyEchant)
cat("Variance : ",varEchant)
hist(resultats, main="Histogramme du modèle Borehole", xlab="Valeur de sortie",col="lightgreen")

#b.
borehole_moyenne <- borehole(val_moy)
cat("Moyenne de la fonction sur l’échantillon :", moyEchant, "\n")
cat("Valeur de la fonction évaluée en la moyenne des entrées :", borehole_moyenne, "\n")

diff <- abs(moyEchant - borehole_moyenne)
cat("Différence :", diff)

#c.
#Calcul du quantile 95% de la sortie
q95 <- quantile(resultats, 0.95)
cat("Quantile d'ordre 95% de la sortie : ", q95)


#Estimation de Monte-Carlo
monte_carlo_inter <- function(N, iter = 10000) {
  resultats_mc_global <- numeric(iter)
  
  for (i in 1:iter) {
    #nouveau échantillon
    echantillon_mc <- EchantBorehole(N)
    resultats_mc <- apply(echantillon_mc, 1, function(x) borehole(as.numeric(x)))
    
    resultats_mc_global[i] <- quantile(resultats_mc, 0.95)
  }
  
  #Calcul les quantiles à 2.5% et 97.5%
  inter_conf <- quantile(resultats_mc_global, c(0.025, 0.975))
  return(inter_conf)
}

inter_conf_95 <- monte_carlo_inter(1000, iter = 100)
cat("Intervalle de confiance à 95% : [", inter_conf_95[1], "; ", inter_conf_95[2], "]")

#d.
#Calcul de la proba que le débit soit supérieur à 250 m^3/an
debit_sup_250 <- function(N) {
  echantillon <- EchantBorehole(N)
  res <- apply(echantillon, 1, function(x) borehole(as.numeric(x)))
  return(mean(res > 250))
}


#Initialisation
taille <- 2 * 10^6
incr <- 2 * 10^6   
err_target <- 0.10

probas <- c()
errs <- c()
tailles <- c(taille)
iter <- 0  
iter_max <- 10

while (iter < iter_max) {
  prob <- debit_sup_250(taille)
  probas <- c(probas, prob)
  
  #Si l'erreur relative est inférieure à err_target
  if (length(probas) > 1) {
    erreur_relative <- abs(probas[length(probas)] - probas[length(probas) - 1]) / probas[length(probas) - 1]
    errs <- c(errs, erreur_relative)
    
    if (erreur_relative < err_target) {
      cat("Erreur relative : ", erreur_relative, "\n")
      cat("Taille nécessaire : ", taille, "\n")
      break
    }
  }
  
  taille <- taille + incr
  tailles <- c(tailles,taille)
  iter <- iter + 1
}


cat(errs,"\n")
cat(tailles)

plot(tailles, probas, type = "b", 
     xlab = "Taille de l'échantillon", ylab = "Probabilité (débit > 250)",
     main = "Convergence de la probabilité", 
     col = "green4", pch = 16, lwd = 2)

if (length(errs) > 0) { #si on a plus qu'une erreur
  plot(seq(taille_init + incr, taille, by = incr), errs, type = "b", 
       xlab = "Taille de l'échantillon", ylab = "Erreur relative",
       main = "Évolution de l'erreur relative",
       col = "darkgreen", pch = 16, lwd = 2, log = "y") 
  
  abline(h = 0.1, col = "red", lty = 2, lwd = 2)  
}

#################################

#2.
#a. methode de Morris

install.packages("sensitivity")
library(sensitivity)


morris_design <- list(
  type = "oat",      
  levels = 6,        
  grid.jump = 3      
)

binf <- sapply(1:13, function(i) min(val_min[i], val_max[i]))
bsup <- sapply(1:13, function(i) max(val_min[i], val_max[i]))
names(binf) <- names(bsup) <- colnames(EchantBorehole(1))

morris_result <- morris(
  model = borehole,
  factors = colnames(EchantBorehole(1)),
  r = 10,   
  design = morris_design,
  binf = binf,
  bsup = bsup
)
plot(morris_result, main = "Morris sensibilite",cex.lab = 1.5,cex = 1.5)


mu_ini <- apply(morris_result$ee, 2, function(x) mean(abs(x)))
non_influential <- names(mu_ini)[mu_ini < 0.1*max(mu_ini)]
cat("Coût total d'évaluation:", 4*(13+1), "\n",
    "Variables non influentes (à fixer):\n",
    paste(non_influential, collapse = ", "), "\n")


#b.i Indices bases sur la correlation/regression

input_mc <- EchantBorehole(1000)
output_mc <- borehole(input_mc)
input_filt <- input_mc[, !colnames(input_mc) %in% non_influential]
output_filt <- output_mc

# scatterplots
pairs(cbind(output = output_filt, input_filt),
      main = "scatters entre la sortie et les entrées",
      col = "blue",
      pch = 19, cex = 0.8)


#b.ii

lm_model <- lm(output_filt ~ ., data = data.frame(input_filt))
src <- coef(lm_model)[-1] * apply(input_filt, 2, sd)/sd(output_filt)
src_sq <- src^2 / sum(src^2)
cat("nom des variables:", colnames(input_filt), "\n",
    "SRC²:\n",
    sort(round(src_sq, 3), decreasing = TRUE), "\n")



#c. Indices de Sobol 
sobol_ind <- sobol2007(
  model = borehole,
  X1 = EchantBorehole(2000),
  X2 = EchantBorehole(2000),
  nboot = 60
)
par(lwd = 2)  

plot(sobol_ind, main = "Sobol Indices")
