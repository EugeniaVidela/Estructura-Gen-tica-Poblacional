# Se presenta el código utilizado para la implementación de los algoritmos de agrupamiento: UPGMA, k-means y Método Bayesiano Sstructure. Estos tres algoritmos se presentan como funciones independientes y que, internamente, calcula los índices de validación del número de grupo: CH, Dunn, Silueta y Conectividad. El código está paralelizado utilizando la función mclapply de la librería "parallel".

# Directorio del Escenario de Simulacion
path<-"/home/evidela/E4/Simulaciones/" 

# Librerías
############

library(parallel)
library(vegan)
library(stats)
library(fpc)
library(caret)
library(clValid)
library(pastecs)
library(LEA)


# Función UPGMA
###############

UPGMA<-function(m, n=15, k, Simulaciones){
  Simulaciones<-list.files(path)
  Base<-read.table(paste(path, Simulaciones[m], sep="/"))[,-c(1:2)]   # Le saco la columna de asignacion y la de fenotipo
  VA<-read.table(paste(path, Simulaciones[m], sep="/"))[,1]
  
  Base<-t(apply(Base, 1, Codif))
  Base[is.na(Base)]<-0
  D<-vegdist(as.matrix(Base), method="jaccard", binary=TRUE)
  
  # Metodo: Vector de Asignacion
  V<-matrix(c(rep(0,nrow(Base)*(n-1))),nrow=nrow(Base),ncol=(n-1), dimnames = list(1:nrow(Base), c(paste0("k=", 2:n))))
  for (l in 2:n) V[,l-1]=as.vector(cutree(hclust(D, method="average"), k=l))
  
  # Indices
  Ind<-data.frame((matrix(NA, nrow=34, ncol=(n-1))))
  for (j in 1:(n-1)) Ind[,j] = as.data.frame(as.matrix(cluster.stats(D,V[,j])))
  names(Ind)<-c(paste0("k=", 2:n))
  rownames(Ind)<-rownames(as.data.frame(as.matrix(cluster.stats(D,V[,1]))))
  C=c()
  for(j in 1:(n-1)) C[j]=connectivity(D,V[,j])
  IndexC<-rbind(Ind, Conectividad=C)
  Index<-IndexC[c(21,26,29,35),]
  rownames(Index)<- c("Silhouette", "Dunn", "Ch", "Conectividad")
  
  # Matriz de Confusion
  MCo<-confusionMatrix(as.factor(V[,(k-1)]), as.factor(VA))
  a<-apply(MCo$table,2,which.max)
  a[duplicated(apply(MCo$table,2,which.max))]<-seq(1,k)[!(1:k%in%a)]
  MC<-MCo$table[a,]
  
  # Proporcion de Mala Clasificación
  diag(MCo$table[a,])<-c(rep(0, ncol(MCo$table[a,])))
  PMA<-sum(MCo$table[a,])/length(VA) # Proporcion de mala asignacion
  return(list(V, Index, MC, PMA))
  
}
UP<-mclapply(1:100, UPGMA, k=2, mc.cores = 8)
save(UP, file="/home/evidela/Metodos/E4/UP.RData")


# Función k-means
#################

KMEANS<-function(m, n=15, k, Simulaciones){
  Simulaciones<-list.files(path)
  Base<-read.table(paste(path, Simulaciones[m], sep="/"))[,-c(1:2)]   # Le saco la columna de asignacion y la de fenotipo
  VA<-read.table(paste(path, Simulaciones[m], sep="/"))[,1]
  
  Base<-t(apply(Base, 1, Codif))
  Base[is.na(Base)]<-0
  D<-vegdist(as.matrix(Base), method="jaccard", binary=TRUE)
  
  # Metodo: Vector de Asignacion
  V<-matrix(c(rep(0,nrow(Base)*(n-1))),nrow=nrow(Base),ncol=(n-1), dimnames = list(1:nrow(Base), c(paste0("k=", 2:n))))
  for (l in 2:n) V[,l-1]=as.vector((kmeans(Base,l)$cluster))
  
  # Indices
  Ind<-data.frame((matrix(NA, nrow=34, ncol=(n-1))))
  for (j in 1:(n-1)) Ind[,j] = as.data.frame(as.matrix(cluster.stats(D,V[,j])))
  names(Ind)<-c(paste0("k=", 2:n))
  rownames(Ind)<-rownames(as.data.frame(as.matrix(cluster.stats(D,V[,1]))))
  C=c()
  for(j in 1:(n-1)) C[j]=connectivity(D,V[,j])
  Index<-IndexC[c(21,26,29,35),]
  rownames(Index)<- c("Silhouette", "Dunn", "Ch", "Conectividad")
  
  
  # Matriz de Confusion
  MCo<-confusionMatrix(as.factor(V[,(k-1)]), as.factor(VA))
  a<-apply(MCo$table,2,which.max)
  a[duplicated(apply(MCo$table,2,which.max))]<-seq(1,k)[!(1:k%in%a)]
  MC<-MCo$table[a,]
  
  # Proporcion de Mala Clasificación
  diag(MCo$table[a,])<-c(rep(0, ncol(MCo$table[a,])))
  PMA<-sum(MCo$table[a,])/length(VA) # Proporcion de mala asignacion
  return(list(V, Index, MC, PMA))
  
}

KM<-mclapply(1:100, KMEANS, k=2, mc.cores = 8)
save(KM, file="/home/evidela/Metodos/E4/KM.RData")


# Función MBS
#############

MBS<-function(m, n=15, k, Simulaciones){
  Simulaciones<-list.files(path)
  Base<-read.table(paste(path, Simulaciones[m], sep="/"))[,-c(1:2)]   # Le saco la columna de asignacion y la de fenotipo
  VA<-read.table(paste(path, Simulaciones[m], sep="/"))[,1]
  Base<-t(apply(Base, 1, Codif))
  Base[is.na(Base)]<-0
  D<-vegdist(as.matrix(Base), method="jaccard", binary=TRUE)
  
  colnames(Base) <- NULL
  rownames(Base) <- NULL
  
  write.geno(Base, paste0(paste(path2,"base",m, sep=""),".geno"))
  
  best.k <- snmf(input.file = paste0(paste(path2,"base",m, sep=""),".geno"), K = 2:n, project = "force", entropy = T)
  
  Qlist<-list(0)
  for (r in 1:(n-1)) Qlist[[r]]<-Q(best.k, K=(r+1), run=1)
  
  # Metodo: Vector de Asignacion
  V<-matrix(c(rep(0,nrow(Base)*(n-1))),nrow=nrow(Base),ncol=(n-1), dimnames = list(1:nrow(Base), c(paste0("k=", 2:n))))
  for (l in 2:n) V[,l-1]=as.vector(apply(Qlist[[l-1]], 1, function(x)which.max(x)))
  
  # Indices
  Ind<-data.frame((matrix(NA, nrow=34, ncol=(n-1))))
  for (j in 1:(n-1)) Ind[,j] = as.data.frame(as.matrix(cluster.stats(D,V[,j])))
  names(Ind)<-c(paste0("k=", 2:n))
  rownames(Ind)<-rownames(as.data.frame(as.matrix(cluster.stats(D,V[,1]))))
  C=c()
  for(j in 1:(n-1)) C[j]=connectivity(D,V[,j])
  IndexC<-rbind(Ind, Conectividad=C)
  Index<-IndexC[c(21,26,29,35),]
  rownames(Index)<- c("Silhouette", "Dunn", "Ch", "Conectividad")
  
  
  # Matriz de Confusion
  MCo<-confusionMatrix(as.factor(V[,(k-1)]), as.factor(VA))
  a<-apply(MCo$table,2,which.max)
  a[duplicated(apply(MCo$table,2,which.max))]<-seq(1,k)[!(1:k%in%a)]
  MC<-MCo$table[a,]
  
  # Proporcion de Mala Asignacion
  diag(MCo$table[a,])<-c(rep(0, ncol(MCo$table[a,])))
  PMA<-sum(MCo$table[a,])/length(VA) # Proporcion de mala asignacion
  return(list(Qlist, V, Index, MC, PMA))
  
}

St<-mclapply(1:50, MBS, k=2, mc.cores = 8)
save(St, file="/home/evidela/Metodos/E4/MBS.RData")
