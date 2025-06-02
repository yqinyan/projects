getwd()
setwd('E:\\2023_2024\\M1\\S1\\Aanalyse de données\\Base de données\\Projet')

# Installation des librairies
install.packages("rpart")
library(rpart)
install.packages('rpart.plot')
library(rpart.plot)
install.packages("C50")
library(C50)
install.packages("tree")
library(tree)
install.packages("ROCR")
library(ROCR)
install.packages("ggplot2")
library(ggplot2)
install.packages("regclass")
library(regclass)
install.packages("cluster")
library(cluster)
install.packages("fpc")
library(fpc)
install.packages("dbscan")
library(dbscan)
install.packages("tsne")
library(tsne)
install.packages("questionr")
library(questionr)


# Chargement des donnees
fraudulent <- read.csv('Data_Projet_1.csv',header = T,sep = ',', dec = '.', stringsAsFactors = TRUE)
str(fraudulent)


# Visualisation des données et pré-traitement des données

names(fraudulent) #Affichage des noms des variables 
#Affichage des effectives
table(fraudulent$fraudulent)
table(fraudulent$gender)
table(fraudulent$claim_area)
table(fraudulent$police_report)
table(fraudulent$claim_type)

# Affichage graphique sectoriel
pie(table(fraudulent$fraudulent), main = "Répartition des classes de fraudulent")
pie(table(fraudulent$gender), main = "Répartition des classes de gender")
pie(table(fraudulent$claim_area), main = "Répartition des classes de claim_area")
pie(table(fraudulent$police_report), main = "Répartition des classes de police_report")
pie(table(fraudulent$claim_type), main = "Répartition des classes de claim_type")

#Histogrammes d'effectifs
qplot(fraudulent, data = fraudulent, main = 'fraudulent', fill=fraudulent,xlab="Class de fraudulent", ylab="Nombre d'instances") + theme(text = element_text(size=16))
qplot(gender, data = fraudulent, main = 'gender', fill=fraudulent)
qplot(claim_area, data = fraudulent, main = 'claim_area', fill=fraudulent)
qplot(police_report, data = fraudulent, main = 'police_report', fill=fraudulent)
qplot(claim_type, data = fraudulent, main = 'claim_type', fill=fraudulent)

# Nuages de points
qplot(claim_type, police_report, data = fraudulent, color=fraudulent)+geom_jitter(height = 0.4)
qplot(claim_area, police_report, data = fraudulent, color=fraudulent)+geom_jitter(height = 0.4, width = 0.4)
qplot(claim_area, claim_type, data = fraudulent, color=fraudulent)+geom_jitter(height = 0.4, width = 0.4)


# Boites a moustaches
boxplot(claim_amount~fraudulent,data = fraudulent, col=c("tomato","darkturquoise"), main = "claim_amount selon fraudulent", xlab = 'fraudulent',ylab = 'claim_amount')
tapply(fraudulent$claim_amount,fraudulent$fraudulent, summary) 

boxplot(days_to_incident~fraudulent,data = fraudulent, col=c("tomato","darkturquoise"), main = "days_to_incident selon fraudulent", xlab = 'fraudulent',ylab = 'days_to_incident')
tapply(fraudulent$claim_amount,fraudulent$fraudulent, summary) 

# Tables de contingence de variables dicrètes
table(fraudulent$claim_area,fraudulent$fraudulent)
table(fraudulent$claim_type,fraudulent$fraudulent)
table(fraudulent$police_report,fraudulent$fraudulent)

# Tables de contingences en proportions
options(digits = 2)
prop.table(table(fraudulent$claim_area,fraudulent$fraudulent))
prop.table(table(fraudulent$claim_type,fraudulent$fraudulent))
prop.table(table(fraudulent$police_report,fraudulent$fraudulent))

# Tables de contingences en pourcentage
prop.table(table(fraudulent$claim_area,fraudulent$fraudulent,dnn=c("claim_area", "fraudulent")))*100
prop.table(table(fraudulent$claim_type,fraudulent$fraudulent,dnn=c("claim_type", "fraudulent")))*100                                                       
prop.table(table(fraudulent$police_report,fraudulent$fraudulent,dnn=c("police_rep", "fraudulent")))*100


#Calcul de la corrélation de Pearson
chisq.test(fraudulent$claim_area, fraudulent$claim_type) #p-value = 0.3618
chisq.test(fraudulent$claim_area, fraudulent$police_report) #p-value = 0.7502
chisq.test(fraudulent$police_report, fraudulent$claim_type) #p-value < 2e-16 donc ne sont pas corrélées  

#Tests de dépendance entre valeurs de variables discrètes
rprop(table(fraudulent$claim_area, fraudulent$fraudulent))
rprop(table(fraudulent$claim_type, fraudulent$fraudulent))
rprop(table(fraudulent$police_report, fraudulent$fraudulent))

chisq.residuals(rprop(table(fraudulent$claim_area, fraudulent$fraudulent)))
chisq.residuals(rprop(table(fraudulent$claim_type, fraudulent$fraudulent)))
chisq.residuals(rprop(table(fraudulent$police_report, fraudulent$fraudulent)))
# les valeurs de police_report indiquae une dépendance significative

# Graphiques en mosaïque des proportions de cooccurrences de valeurs de variables discrètes
mosaic(fraudulent$claim_area~fraudulent$fraudulent, data=fraudulent, color=TRUE)
mosaic(fraudulent$claim_type~fraudulent$fraudulent, data=fraudulent, color=TRUE) 
mosaic(fraudulent$police_report~fraudulent$fraudulent, data=fraudulent, color=TRUE) 


#l'ensemble d'apprentissage
#claim_id et customer_id ne constituant pas une information utile pour l'objectif poursuivi
apprenti <- fraudulent[1:750,]
apprenti <- subset(apprenti, select = -claim_id)
apprenti <- subset(apprenti, select = -customer_id)
table(apprenti$fraudulent)
#l'emsemble de test
test <- fraudulent[751:1100,]
#les caractéristiques principales des data frames
summary(apprenti)
summary(test)

# Apprentissage d'un arbre de décision 
tree1 <- rpart(fraudulent~., apprenti)
tree2 <- C5.0(fraudulent~., apprenti)
tree3 <- tree(fraudulent~., apprenti)
# Représentation graphique des arbres de décision 
plot(tree1)
text(tree1, pretty = 0)
prp(tree1, cex = 0.75)

plot(tree2, type="simple")

plot(tree3)
text(tree3, pretty = 3)
#Test des arbres 
test_tree1 <- predict(tree1, test, type="class")
test_tree2 <- predict(tree2, test, type="class")
test_tree3 <- predict(tree3, test, type="class")

table(test$fraudulent)
table(test_tree1)
table(test_tree2)
table(test_tree3)

#Calcul des taux de succès
test$Prediction1 <- test_tree1  #creations une colonnes de prediction par tree1
test$Prediction2 <- test_tree2  #creations une colonnes de prediction par tree2
test$Prediction3 <- test_tree3  #creations une colonnes de prediction par tree3

nbr_succes1 <- nrow(test[test$fraudulent==test$Prediction1,])
nbr_succes1 #267
nbr_succes2 <- nrow(test[test$fraudulent==test$Prediction2,])
nbr_succes2 #270
nbr_succes3 <- nrow(test[test$fraudulent==test$Prediction3,])
nbr_succes3 #275

taux_succes1 <- nbr_succes1/nrow(test)
taux_succes1  # 76.28%
taux_succes2 <- nbr_succes2/nrow(test)
taux_succes2  # 77.14%
taux_succes3 <- nbr_succes3/nrow(test)
taux_succes3  # 78.57%

# d'apres les calculs des taux de succes des trois predictions on trouve que le classifieur tree est plus performant

#Précision des prédictions
#Affichez la liste des probabilités de prédiction
proba_tree1 <- predict(tree1, test, type = "prob") 
print(proba_tree1)
proba_tree2 <- predict(tree2, test, type = "prob") 
print(proba_tree2)
proba_tree3 <- predict(tree3, test, type = "vector") 
print(proba_tree3)
#Affichage de la liste des probabilités de prédiction fraudulent = Yes pour les 3 
print(proba_tree1[,2])
print(proba_tree2[,2])
print(proba_tree3[,2])
#Affichage de la liste des probabilités de prédiction fraudulent = No
print(proba_tree1[,1])
print(proba_tree12[,1])
print(proba_tree3[,1])

#Construction de data frame pour afficher des statistiques concernant ces probabilités de prédictions

#par tree1
df_result1 <-  data.frame(test$fraudulent, test_tree1, proba_tree1[,2], proba_tree1[,1])
colnames(df_result1) = list('class','prediction','fraud_yes','fraud_no')
View(df_result1)

#par tree2
df_result2 <-  data.frame(test$fraudulent, test_tree2, proba_tree2[,2], proba_tree2[,1])
colnames(df_result2) = list('class','prediction','fraud_yes','fraud_no')
View(df_result2)

#par tree3
df_result3 <-  data.frame(test$fraudulent, test_tree3, proba_tree3[,2], proba_tree3[,1])
colnames(df_result3) = list('class','prediction','fraud_yes','fraud_no')
View(df_result3)

#afficher les quartiles et la moyenne des probabilités des prédiction  pour la classe fraudulent = Yes
summary(df_result1[df_result1$prediction=="Yes", "fraud_yes"]) #la moyenne des probabilités est 71.69%
summary(df_result2[df_result2$prediction=="Yes", "fraud_yes"]) #la moyenne des probabilités est 68.03%
summary(df_result3[df_result3$prediction=="Yes", "fraud_yes"]) #la moyenne des probabilités est 100%

# Test de paramètres de l’apprentissage d’arbres de décision rpart()
coef_Gini_10 <- rpart(fraudulent~., apprenti, parms = list(split = "gini"), control = rpart.control(minbucket = 10))
coef_Gini_5 <- rpart(fraudulent~., apprenti, parms = list(split = "gini"), control = rpart.control(minbucket = 5))
info_gain_10 <- rpart(fraudulent~., apprenti, parms = list(split = "information"), control = rpart.control(minbucket = 10))
info_gain_5 <- rpart(fraudulent~., apprenti, parms = list(split = "information"), control = rpart.control(minbucket = 5))

plot(coef_Gini_10)
text(coef_Gini_10, pretty = 0)

plot(coef_Gini_5)  #différent que les 3 autres
text(coef_Gini_5, pretty = 0)

plot(info_gain_10)  
text(info_gain_10, pretty = 0)

plot(info_gain_5)
text(info_gain_5, pretty = 0)
# on obtien 2 arbres de décision différents seulement parmi les 4.
# claim_amount n’a pas eu d’impact sur l’apprentissage.

#Appliquation de l’arbre à l'ensemble de test
test_tree11 <- predict(info_gain_10, test, type="class") 
print(test_tree11)
table(test_tree11)

test_tree22 <- predict(info_gain_5, test, type = "class")
print(test_tree22)
table(test_tree22)

test_tree33 <- predict(coef_Gini_10, test, type = "class")
print(test_tree33)
table(test_tree33)

test_tree44 <- predict(coef_Gini_5, test, type = "class")
print(test_tree44)
table(test_tree44)

test$Tree1 <- test_tree11
test$Tree2 <- test_tree22
test$Tree3 <- test_tree33
test$Tree4 <- test_tree44
View(test [,c("fraudulent","Tree1","Tree2","Tree3","Tree4")])

# taux de succes
taux_succes11 <- length(test[test$fraudulent==test$Tree1,'claim_id']) / nrow(test)
taux_succes22 <- length(test[test$fraudulent==test$Tree2,'claim_id']) / nrow(test)
taux_succes33 <- length(test[test$fraudulent==test$Tree3,'claim_id']) / nrow(test)
taux_succes44 <- length(test[test$fraudulent==test$Tree4,'claim_id']) / nrow(test)

taux_succes11 #77.14%
taux_succes22 #77.14%
taux_succes33 #77.14%
taux_succes44 #76.28%

# Test de paramètres de l’apprentissage d’arbres de décision C5.0()

c50_minCa10_Fal <- C5.0(fraudulent~., apprenti, control = C5.0Control(minCases = 10, noGlobalPruning = FALSE))
c50_minCa10_Tur <- C5.0(fraudulent~., apprenti, control = C5.0Control(minCases = 10, noGlobalPruning = TRUE))
c50_minCa5_Fal <- C5.0(fraudulent~., apprenti, control = C5.0Control(minCases = 5, noGlobalPruning = FALSE))
c50_minCa5_Tur <- C5.0(fraudulent~., apprenti, control = C5.0Control(minCases = 5, noGlobalPruning = TRUE))

plot(c50_minCa10_Fal, type = "simple")

plot(c50_minCa10_Tur, type = "simple")

plot(c50_minCa5_Fal, type = "simple")

plot(c50_minCa5_Tur, type = "simple")

#le paramètre police_report qui n'a eu d’impact sur l’apprentissage


test_c50_tree1 <- predict(c50_minCa10_Fal, test, type="class") 
print(test_c50_tree1)
table(test_c50_tree1)

test_c50_tree2 <- predict(c50_minCa10_Tur, test, type = "class")
print(test_c50_tree2)
table(test_c50_tree2)

test_c50_tree3 <- predict(c50_minCa5_Fal, test, type = "class")
print(test_c50_tree3)
table(test_c50_tree3)

test_c50_tree4 <- predict(c50_minCa5_Tur, test, type = "class")
print(test_c50_tree4)
table(test_c50_tree4)

test$c50Tree1 <- test_c50_tree1
test$c50Tree2 <- test_c50_tree2
test$c50Tree3 <- test_c50_tree3
test$c50Tree4 <- test_c50_tree4
View(test [,c("fraudulent","c50Tree1","c50Tree2","c50Tree3","c50Tree4")])

taux_c50_succes1 <- length(test[test$fraudulent==test$c50Tree1,"claim_id"]) / nrow(test)
taux_c50_succes2 <- length(test[test$fraudulent==test$c50Tree2,"claim_id"]) / nrow(test)
taux_c50_succes3 <- length(test[test$fraudulent==test$c50Tree3,"claim_id"]) / nrow(test)
taux_c50_succes4 <- length(test[test$fraudulent==test$c50Tree4,"claim_id"]) / nrow(test)

taux_c50_succes1 #79.14%
taux_c50_succes2 #77.42%
taux_c50_succes3 #77.14%
taux_c50_succes4 #76.57%

# Test de paramètres de l’apprentissage d’arbres de décision tree()

tree_dev10 <- tree(fraudulent~., apprenti, split = "deviance", control = tree.control(nrow(apprenti),mincut = 10))
tree_dev5 <- tree(fraudulent~., apprenti, split = "deviance", control = tree.control(nrow(apprenti),mincut = 5))
tree_gini10 <- tree(fraudulent~., apprenti, split = "gini", control = tree.control(nrow(apprenti),mincut = 10))
tree_gini5 <- tree(fraudulent~., apprenti, split = "gini", control = tree.control(nrow(apprenti),mincut = 5))

plot(tree_dev10)
text(tree_dev10,pretty = 3)

plot(tree_dev5)
text(tree_dev5,pretty = 3)

plot(tree_gini10)
text(tree_gini10,pretty = 3)

plot(tree_gini5)
text(tree_gini5,pretty = 3)

#pas de paramètre qui n'a eu d’impact sur l’apprentissage.

test_tr_tree1 <- predict(tree_dev10, test, type="class") 
print(test_tr_tree1)
table(test_tr_tree1)

test_tr_tree2 <- predict(tree_dev5, test, type = "class")
print(test_tr_tree2)
table(test_tr_tree2)

test_tr_tree3 <- predict(tree_gini10, test, type = "class")
print(test_tr_tree3)
table(test_tr_tree3)

test_tr_tree4 <- predict(tree_gini5, test, type = "class")
print(test_tr_tree4)
table(test_tr_tree4)

test$tree1 <- test_tr_tree1
test$tree2 <- test_tr_tree2
test$tree3 <- test_tr_tree3
test$tree4 <- test_tr_tree4
View(test [,c("fraudulent","tree1","tree2","tree3","tree4")])

taux_tr_succes1 <- length(test[test$fraudulent==test$tree1,"claim_id"]) / nrow(test)
taux_tr_succes2 <- length(test[test$fraudulent==test$tree2,"claim_id"]) / nrow(test)
taux_tr_succes3 <- length(test[test$fraudulent==test$tree3,"claim_id"]) / nrow(test)
taux_tr_succes4 <- length(test[test$fraudulent==test$tree4,"claim_id"]) / nrow(test)
options(digits = 5)
taux_tr_succes1 #78.57%
taux_tr_succes2 #78.57%
taux_tr_succes3 #72.00%
taux_tr_succes4 #70.57%

#calcul des matrices de confusion et mesures d’évaluation 

mc_tree1 <- table(test$fraudulent,test_tree1)
print(mc_tree1)
VP1 <- mc_tree1[2,2]
FP1 <- mc_tree1[1,2]
VN1 <- mc_tree1[1,1]
FN1 <- mc_tree1[2,1]
print(VP1)
print(FP1)
print(VN1)
print(FN1)

mc_tree2 <- table(test$fraudulent,test_tree2)
print(mc_tree2)
VP2 <- mc_tree2[2,2]
FP2 <- mc_tree2[1,2]
VN2 <- mc_tree2[1,1]
FN2 <- mc_tree2[2,1]
print(VP2)
print(FP2)
print(VN2)
print(FN2)

mc_tree3 <- table(test$fraudulent,test_tree3)
print(mc_tree3)
VP3 <- mc_tree3[2,2]
FP3 <- mc_tree3[1,2]
VN3 <- mc_tree3[1,1]
FN3 <- mc_tree3[2,1]
print(VP3)
print(FP3)
print(VN3)
print(FN3)

# taux de  succès

Rappel1 <- mc_tree1[2,2]/(mc_tree1[2,2]+mc_tree1[2,1])
Specificit1 <-mc_tree1[1,1]/(mc_tree1[1,1]+mc_tree1[1,2])
Precision1 <- mc_tree1[2,2]/(mc_tree1[1,2]+mc_tree1[2,2])
TauxN1 <- mc_tree1[1,1]/(mc_tree1[1,1]+mc_tree1[2,1])

Rappel2 <- mc_tree2[2,2]/(mc_tree2[2,2]+mc_tree2[2,1])
Specificit2 <-mc_tree2[1,1]/(mc_tree2[1,1]+mc_tree2[1,2])
Precision2 <- mc_tree2[2,2]/(mc_tree2[1,2]+mc_tree2[2,2])
TauxN2 <- mc_tree2[1,1]/(mc_tree2[1,1]+mc_tree2[2,1])

Rappel3 <- mc_tree3[2,2]/(mc_tree3[2,2]+mc_tree3[2,1])
Specificit3 <-mc_tree3[1,1]/(mc_tree3[1,1]+mc_tree3[1,2])
Precision3 <- mc_tree3[2,2]/(mc_tree3[1,2]+mc_tree3[2,2])
TauxN3 <- mc_tree3[1,1]/(mc_tree3[1,1]+mc_tree3[2,1])

print(Rappel1)  #30.48%
print(Specificit1) #90.29%
print(Precision1) #49.01%
print(TauxN1) #80.93%

print(Rappel2) #24.39%
print(Specificit2) #93.28%
print(Precision2) #52.63%
print(TauxN2) #80.12%

print(Rappel3) #8.53%
print(Specificit3) #100%
print(Precision3) #100%
print(TauxN3)  #78.13%
# le classifieur le plus performant parmi les trois arbres de décision 
#selon le critère d'évaluation du Rappel est tree1

# le classifieur le plus performant parmi les trois arbres de décision 
#selon le critère d'évaluation de la Spécificité est tree3

# le classifieur le plus performant parmi les trois arbres de décision 
#selon le critère d'évaluation de la Précision est tree3

# le classifieur le plus performant parmi les trois arbres de décision 
#selon le critère d'évaluation du Taux de Vrais Négatifs est tree1

#Calcul des courbes ROC
roc_pred1 <- prediction(proba_tree1[,2], test$fraudulent)
roc_pref1 <- performance(roc_pred1, "tpr", "fpr")
plot(roc_pref1, col = 'green')
legend("bottomright",legend='tree1',col = "green", lty = 1)

roc_pred2 <- prediction(proba_tree2[,2], test$fraudulent)
roc_pref2 <- performance(roc_pred2, "tpr", "fpr")
plot(roc_pref2, col = 'orange',add = TRUE)
legend("right",legend='tree2',col = "orange", lty = 1)

roc_pred3 <- prediction(proba_tree3[,2], test$fraudulent)
roc_pref3 <- performance(roc_pred3, "tpr", "fpr")
plot(roc_pref3, col ='black',add = TRUE)
legend("topleft",legend='tree3',col = "black", lty = 1)
#d'après la comparaision des courbes ROC, tree3 est le plus performant

#Calcul des indices AUC
auc_tree1 <- performance(roc_pred1,"auc")
str(auc_tree1)
attr(auc_tree1,"y.values")  #61.82%

auc_tree2 <- performance(roc_pred2,"auc")
str(auc_tree2)
attr(auc_tree2,"y.values")  #58.22%

auc_tree3 <- performance(roc_pred3,"auc")
str(auc_tree3)
attr(auc_tree3,"y.values")  #63.02%
#d'après le calcule des aucs des trees, le plus performant est tree3



# matrice de distance
dmatrix <- daisy(fraudulent)
summary(dmatrix)


#Clustering par partitionnement
#le clustering par K-means
for (k in 4:7){
  km <- kmeans(dmatrix, k)
  print(table(km$cluster, fraudulent$fraudulent))
  print(qplot(km$cluster, data=fraudulent, fill=fraudulent,main = paste0('k =',k)))
} 

# k=6 est le plus performant
#le clustering par K-means en fixant le nombre de clusters à 6
km6 <- kmeans(dmatrix, 6)
fraudulent_km6 <- data.frame(fraudulent, km6$cluster)
View(fraudulent_km6)
table(km6$cluster, fraudulent$fraudulent)

qplot(km6$cluster, data = fraudulent, fill = fraudulent,main='km6$cluster selon fraudulent')

qplot(police_report, km6$cluster, data = fraudulent, color = fraudulent)+geom_jitter(height = 0.3, width = 0.3)
qplot(claim_area,km6$cluster, data = fraudulent, color = fraudulent)+geom_jitter(height = 0.3, width = 0.3)


#Clustering hiérarchique par agglomération
agn <- agnes(dmatrix)
plot(agn)
rect.hclust(agn, k = 4, border = 'red')
plot(agn)
for (k in 3:8) {
  print(rect.hclust(agn, k, border = 'red'))
  print(agnk <- cutree(agn, k))
  print(table(agnk,fraudulent$fraudulent))
  print(qplot(agnk, data = fraudulent, fill = fraudulent, main = paste0('k=',k)))
} #k=8


agn8 <- cutree(agn ,k =8)
View(agn8)
table(agn8, fraudulent$fraudulent) #table de contingence
qplot(agn8, data = fraudulent,fill = fraudulent)


#Clustering hiérarchique par division

dia <- diana(dmatrix)
plot(dia)
rect.hclust(dia, k =3, border = 'red')

plot(dia)
for (k in 3:8) {
  print(rect.hclust(dia, k, border = 'red'))
  print(diak <- cutree(dia,k))
  print(table(diak, fraudulent$fraudulent))
  print(qplot(diak,data=fraudulent,fill=fraudulent, main=paste0('k=',k)))
} #k=3

dia3 <- cutree(dia, k =3)
View(dia3)
table(dia3, fraudulent$fraudulent)
qplot(dia3, data = fraudulent,fill=fraudulent, main='dia3/fraudulent')
# 2 clusters la proportion d'instances de même classe dans chaque cluster est maximale


# Clustering basé sur la densité
# x pour une distance de voisinage (eps) de 0.125 et une taille minimale de voisinage (minPts) de 3

# avec la fonction dbscan
dbs <- dbscan(dmatrix, eps = 0.125, minPts = 3)
table(dbs$cluster,fraudulent$fraudulent)
qplot(as.factor(dbs$cluster),data = fraudulent,fill=fraudulent)
#Évaluation interne des clusters par t-SNE
tsne_out <- tsne(dmatrix, k=2)
tsne_out <- data.frame(tsns_out)
qplot(tsne_out[,1], tsne_out[,2], col=as.factor(dbs$cluster))

for (k in 3:8) {
  dbsm <- dbscan(dmatrix, eps = 0.125, minPts = k)
  print(table(dbsm$cluster,fraudulent$fraudulent))
  print(qplot(as.factor(dbsm$cluster),data = fraudulent,fill=fraudulent, main = paste0('k=',k)))
}  #k=7

# avec la fonction optics
for (k in 3:8) {
  opt <- optics(dmatrix, eps = 0.125, minPts = k)
  print(table(opt$cluster,fraudulent$fraudulent))
  print(qplot(as.factor(opt$cluster),data = fraudulent,fill=fraudulent, main = paste0('k=',k)))
}  # ne marche pas

# avec la fonction hdbscan
for (k in 3:8) {
  hdb <- hdbscan(dmatrix, minPts = k)
  print(table(hdb$cluster,fraudulent$fraudulent))
  print(qplot(as.factor(hdb$cluster),data = fraudulent,fill=fraudulent, main = paste0('k=',k)))
}  #k=4

# la fonction hdbscan est la plus performante
options(digits = 6)


# Condclusion: avec les analyes faites, on trouve le classifieur tree qui est plus performant 

# Application du classifieur aux données à prédire.
pred <- read.csv('Data_Projet_1_New.csv',header = T,sep = ',', dec = '.', stringsAsFactors = TRUE)
pred_tree3 <- predict(tree3, pred, type="class")
pred$fraudulent_predict <- pred_tree3
View(pred)
table(pred$fraudulent_predict)
proba <- predict(tree3, pred, type = "vector") 
print(proba)

df_result <-  data.frame(pred_tree3, proba[,2], proba[,1])
colnames(df_result) = list('prediction','fraud_yes','fraud_no')
View(df_result)
table(df_result$prediction)
summary(df_result[df_result$prediction=="Yes", "fraud_yes"])

min(df_result$fraud_yes)
min(df_result$fraud_no)

max(df_result$fraud_yes)
max(df_result$fraud_no)

mean(df_result$fraud_yes)
mean(df_result$fraud_no)

#Enregistrement des prédictions dans le fichier resultat.csv
write.table(pred,"resultats.csv", sep = "\t", dec = ".", row.names = F)
