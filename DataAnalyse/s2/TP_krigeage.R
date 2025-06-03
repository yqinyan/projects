################################################################################
########################## Atelier de B. Iooss #################################
######################## MIM 1 Cote d'Azur 2024 ################################
################################################################################

# Etude de précipitations en Suisse (Swiss rainfall)
# Ref: Christensen, O.F., Diggle, P.J. and Ribeiro Jr, P.J. (2001) Analysing positive-valued spatial data: the transformed Gaussian model. In Monestiez, P., Allard, D. and Froidevaux (eds), GeoENV III - Geostatistics for environmental applications. Quantitative Geology and Geostatistics, Kluwer Series, 11, 287–298

# Ce TP est inspiré de : Manuel d’analyse spatiale (rapport INSEE), chapitre 5 de Jean-Michel Floch

graphics.off()
rm(list=ls())

# 3 bases dans le package geoR :
# — sic.100 : échantillon de 100 observations qui pourront servir à effectuer les interpolations ;
# — sic.367 : observations non incluses dans l’échantillon qui permettront de comparer les estimations et les observations ;
# — sic.all : ensemble.
install.packages("geoR")
library(geoR)

help(SIC)
summary(sic.all)

################################################################################
# 1 Visualisation des donnees

x11() ; points(sic.100, borders=sic.borders,col="green")

x11()
points(sic.100, borders=sic.borders,col="green")
points(sic.367, borders=sic.borders,col="red",add=TRUE)

x11() ; plot.geodata(sic.100,bor=sic.borders) # quartiles in the first plot are divided following blue, green, yellow and red

################################################################################
# 2 Analyse variographique

help(variog)

# variogram cloud
vario.c <- variog(sic.100, op="cloud")
#binned variogram and stores the cloud
vario.bc <- variog(sic.100, bin.cloud=TRUE)
# smoothed variogram
vario.s <- variog(sic.100, op="sm", band=10)
# binned variogram
vario.b <- variog(sic.100)

x11() ; par(mfrow=c(2,2))
plot(vario.c, main="variogram cloud")
plot(vario.bc, bin.cloud=TRUE, main="clouds for binned variogram")
plot(vario.s, main="smoothed variogram", ylim=c(0,25000), xlim=c(0,300))
plot(vario.b, main="binned variogram", ylim=c(0,25000), xlim=c(0,300)) 

# variogrammes directionnels
help(variog4)
vario4 <- variog4(sic.100)
x11() ; plot(vario4,same=FALSE)

################################################################################
# 3 Ajustement d'un modèle de variogramme

vario.ex<- variog(sic.100,option="bin")
x11() ; plot(vario.ex) # nous donne une idee des parametres de variance et de portee

#exemple de fit
help(variofit)
vario.sphe <- variofit(vario.ex,cov.model= "spher", ini.cov.pars=c(14522,69))
print(vario.sphe)
vario.exp <- variofit(vario.ex,cov.model= "exp", ini.cov.pars=c(15000,100))
print(vario.exp)
x11() ; par(mfrow=c(1,2))
plot(vario.ex,main="Sphérique")
lines.variomodel(cov.model="sphe",cov.pars=c(vario.sphe$cov.pars[1],vario.sphe$cov.pars[2]), nugget=0, max.dist=350)
plot(vario.ex,main="Exponentiel")
lines.variomodel(cov.model="exp",cov.pars=c(vario.exp$cov.pars[1],vario.exp$cov.pars[2]), nugget=0, max.dist=350)

# autre methode de fit avec likfit (max de vraisemblance)

# ajustement a l'oeil
x11() ; par(mfrow=c(2,2), mar=c(3,3,1,1), mgp =c (2,1,0))
plot(vario.ex,main="Sphérique")
lines.variomodel(cov.model="sphe",cov.pars=c(15000,100), nugget=0, max.dist=350)
plot(vario.ex,main="Exponentiel")
lines.variomodel(cov.model="exp",cov.pars=c(15000,100), nugget=0, max.dist=350)
plot(vario.ex,main="Exponentiel")
lines.variomodel(cov.model="exp",cov.pars=c(15000,30), nugget=0, max.dist=350)
plot(vario.ex,main="Exponentiel avec pépite")
lines.variomodel(cov.model="exp",cov.pars=c(10000,100), nugget=5000, max.dist=350)

################################################################################
# 4 krigeage

vario.ex <- variog(sic.100, bin.cloud=TRUE)
x11() ; plot(vario.ex,main="")
lines.variomodel(cov.model="spher",cov.pars=c(15000,50), nug=0,max.dist=300)

pred.grid <- expand.grid(seq(0,350, l=51), seq(0,220, l=51))
rgb.palette <- colorRampPalette(c("blue", "lightblue", "orange", "red"),space = "rgb")

help(krige.conv)
kc <- krige.conv(sic.100, loc = pred.grid, krige=krige.control(cov.model="spherical",cov.pars=c(14600,69)))
x11() ; image(kc, loc = pred.grid,col =rgb.palette(20) ,xlab="Coord X", ylab="Coord Y",borders=sic.borders,main="Estimation")
x11() ; image(kc, kc$krige.var, loc = pred.grid,col=rgb.palette(20), xlab="Coord X",ylab="Coord Y",borders=sic.borders, main="Variance de krigeage")
# names(kc)

# validation du krigeage

kc1<- krige.conv(sic.100, loc = sic.100$coords, krige=krige.control(cov.model="spherical",cov.pars=c(15000,70)))
kc2<- krige.conv(sic.100, loc = sic.367$coords, krige=krige.control(cov.model="spherical",cov.pars=c(15000,70)))
x11() ; par(mfrow=c(1,2))
plot(sic.100$data,kc1$predict,xlab="Observe",ylab="Predite", main="Echantillon")
abline(a=0,b=1,col="red")
plot(sic.367$data,kc2$predict, xlab="Observe",ylab="Predite", main="Tests")
abline(a=0,b=1,col="red")

# erreur de test (mean square error)

MSE <- mean((sic.367$data-kc2$predict)^2)
print(MSE)
Q2 <- 1 - MSE/var(sic.367$data)
print(Q2)

################################################################################
# 5 Planification d'experiences
install.packages("DiceDesign")
library(DiceDesign)

###########
# Methode a)

R <- 1e4 # nb de repetitions
critere_opt <- 0
for (i in 1:R){
  sel <- sample(1:467, 45) # selection de 45 numeros d'observations parmi 467 (taille de sic.all)
  coords <- sic.all$coords[sel,]
  critere <- mindist(coords)
  if (critere > critere_opt){
    critere_opt <- critere
    sel_opt <- sel
  }
}

# Visualisation des resultats et calcul de l'erreur de test

# sic.sel contient le resultat de la selection des 45 stations
sic.sel <- sic.all
class(sic.sel) <- "geodata"
sic.sel$coords <- sic.all$coords[sel_opt,]
sic.sel$data <- sic.all$data[sel_opt]
sic.sel$covariate <- sic.all$covariate[[1]][sel_opt]

# sic.reste contient toutes les autres stations
sic.reste <- sic.all
class(sic.reste) <- "geodata"
sic.reste$coords <- sic.all$coords[-sel_opt,]
sic.reste$data <- sic.all$data[-sel_opt]
sic.reste$covariate <- sic.all$covariate[[1]][-sel_opt]

x11() ; points(sic.sel, borders=sic.borders,col="red")
x11() ; points(sic.all, borders=sic.borders,col="black")

kc<- krige.conv(sic.sel, loc = pred.grid, krige=krige.control(cov.model="spher",cov.pars=c(15000,50)))
x11() ; image(kc, loc = pred.grid,col=rgb.palette(20) ,xlab="Coord X", ylab="Coord Y",borders=sic.borders,main="Estimation")
x11() ; image(kc, kc$krige.var, loc = pred.grid,col=rgb.palette(20), xlab="Coord X",ylab="Coord Y",borders=sic.borders, main="Variance de krigeage")

kc1<- krige.conv(sic.sel, loc = sic.sel$coords, krige=krige.control(cov.model="spherical",cov.pars=c(15000,50)))
kc2<- krige.conv(sic.sel, loc = sic.reste$coords, krige=krige.control(cov.model="spherical",cov.pars=c(15000,50)))
x11() ; par(mfrow=c(1,2))
plot(sic.sel$data,kc1$predict,xlab="Observe",ylab="Predite", main="Echantillon")
abline(a=0,b=1,col="red")
plot(sic.reste$data,kc2$predict, xlab="Observe",ylab="Predite", main="Tests")
abline(a=0,b=1,col="red")

MSE <- mean((sic.reste$data-kc2$predict)^2)
print(MSE)
Q2 <- 1 - MSE/var(sic.reste$data)
print(Q2)

#############
# Methode b)

r <- 28 # rayon du cercle de suppression des points
x <- wspDesign(sic.all$coords, dmin = r, init = "center")$design
print(dim(x)[1])

# A faire : trouver le r optimal de maniere automatique (avec une boucle 'while' par exemple)

# on recupere les indices dans sic.all de la selection des 45 stations
sel <- NULL
for (i in 1:dim(x)[1]) sel <- c(sel, which(x[i,1] == sic.all$coords[,1] & x[i,2] == sic.all$coords[,2]))
sel_opt <- as.vector(sel)

# Visualisation des resultats et calcul de l'erreur de test

# sic.sel contient le resultat de la selection des 45 stations
sic.sel <- sic.all
class(sic.sel) <- "geodata"
sic.sel$coords <- sic.all$coords[sel_opt,]
sic.sel$data <- sic.all$data[sel_opt]
sic.sel$covariate <- sic.all$covariate[[1]][sel_opt]

# sic.reste contient toutes les autres stations
sic.reste <- sic.all
class(sic.reste) <- "geodata"
sic.reste$coords <- sic.all$coords[-sel_opt,]
sic.reste$data <- sic.all$data[-sel_opt]
sic.reste$covariate <- sic.all$covariate[[1]][-sel_opt]

x11() ; points(sic.sel, borders=sic.borders,col="red")
x11() ; points(sic.all, borders=sic.borders,col="black")

kc<- krige.conv(sic.sel, loc = pred.grid, krige=krige.control(cov.model="spher",cov.pars=c(15000,50)))
x11() ; image(kc, loc = pred.grid,col=rgb.palette(20) ,xlab="Coord X", ylab="Coord Y",borders=sic.borders,main="Estimation")
x11() ; image(kc, kc$krige.var, loc = pred.grid,col=rgb.palette(20), xlab="Coord X",ylab="Coord Y",borders=sic.borders, main="Variance de krigeage")

kc1<- krige.conv(sic.sel, loc = sic.sel$coords, krige=krige.control(cov.model="spherical",cov.pars=c(15000,50)))
kc2<- krige.conv(sic.sel, loc = sic.reste$coords, krige=krige.control(cov.model="spherical",cov.pars=c(15000,50)))
x11() ; par(mfrow=c(1,2))
plot(sic.sel$data,kc1$predict,xlab="Observe",ylab="Predite", main="Echantillon")
abline(a=0,b=1,col="red")
plot(sic.reste$data,kc2$predict, xlab="Observe",ylab="Predite", main="Tests")
abline(a=0,b=1,col="red")

MSE <- mean((sic.reste$data-kc2$predict)^2)
print(MSE)
Q2 <- 1 - MSE/var(sic.reste$data)
print(Q2)


###############
# Methode c)
install.packages("FNN")
library(FNN)

# 1er indice
sel <- sample(1:467, 1)
sic.sel <- sic.all
class(sic.sel) <- "geodata"
sic.sel$coords <- matrix(sic.all$coords[sel,], ncol=2)
sic.reste <- sic.all
class(sic.reste) <- "geodata"
sic.reste$coords <- sic.all$coords[-sel,]

#indices suivants
for (i in 1:44){
  j <- which.max(knnx.dist(sic.sel$coords, sic.reste$coords, k=1)) # indice dans sic.reste
  sel <- c(sel, which(sic.reste$coords[j,1] == sic.all$coords[,1] & sic.reste$coords[j,2] == sic.all$coords[,2]))
  sic.sel$coords <- sic.all$coords[sel,]
  sic.reste$coords <- sic.all$coords[-sel,]
}

sel_opt <- as.vector(sel)

# Visualisation des resultats et calcul de l'erreur de test

# sic.sel contient le resultat de la selection des 45 stations
sic.sel <- sic.all
class(sic.sel) <- "geodata"
sic.sel$coords <- sic.all$coords[sel_opt,]
sic.sel$data <- sic.all$data[sel_opt]
sic.sel$covariate <- sic.all$covariate[[1]][sel_opt]

# sic.reste contient toutes les autres stations
sic.reste <- sic.all
class(sic.reste) <- "geodata"
sic.reste$coords <- sic.all$coords[-sel_opt,]
sic.reste$data <- sic.all$data[-sel_opt]
sic.reste$covariate <- sic.all$covariate[[1]][-sel_opt]

x11() ; points(sic.sel, borders=sic.borders,col="red")
x11() ; points(sic.all, borders=sic.borders,col="black")

kc<- krige.conv(sic.sel, loc = pred.grid, krige=krige.control(cov.model="spher",cov.pars=c(15000,50)))
x11() ; image(kc, loc = pred.grid,col=rgb.palette(20) ,xlab="Coord X", ylab="Coord Y",borders=sic.borders,main="Estimation")
x11() ; image(kc, kc$krige.var, loc = pred.grid,col=rgb.palette(20), xlab="Coord X",ylab="Coord Y",borders=sic.borders, main="Variance de krigeage")

kc1<- krige.conv(sic.sel, loc = sic.sel$coords, krige=krige.control(cov.model="spherical",cov.pars=c(15000,70)))
kc2<- krige.conv(sic.sel, loc = sic.reste$coords, krige=krige.control(cov.model="spherical",cov.pars=c(15000,70)))
x11() ; par(mfrow=c(1,2))
plot(sic.sel$data,kc1$predict,xlab="Observe",ylab="Predite", main="Echantillon")
abline(a=0,b=1,col="red")
plot(sic.reste$data,kc2$predict, xlab="Observe",ylab="Predite", main="Tests")
abline(a=0,b=1,col="red")

MSE <- mean((sic.reste$data-kc2$predict)^2)
print(MSE)
Q2 <- 1 - MSE/var(sic.reste$data)
print(Q2)
