# Se presenta el código utilizado para simular 100 réplicas del escenario de simulación 2. Cada réplica cuenta con n=1000 individuos, p=80000 marcadores SNPs y k=2 subpoblaciones con baja divergencia genética. El código está paralelizado utilizando la función mclapply de la librería "parallel".

# Librerías
###########
library("xbreed")
library("parallel")
library("StAMPP")
library("spDataLarge")

# Argumentos definidos a priori
###############################

args <- commandArgs(TRUE)

sim_start <-1
sim_end <- 100
if (length(args) == 2) {
  sim_start <- as.numeric(args[[1]])
  sim_end <- as.numeric(args[[2]])
}
cat(paste0("Simular desde ", sim_start, " hasta ", sim_end, "\n"))
n_cores <- 8
if (is.na(n_cores)) {
  n_cores <- Sys.getenv("SLURM_CPUS_PER_TASK")
  if (n_cores == "") {
    n_cores <- max(1, detectCores() - 1)
  }
}
cat(paste0("Corriendo en NCORES: ", n_cores, "\n"))

# Genoma
#########

genome <- data.frame(matrix(NA, nrow = 10, ncol = 6))
names(genome) <- c("chr", "len", "nmrk", "mpos", "nqtl", "qpos")
genome$chr <- c(1:10)
genome$len <- rep(8000,10)
genome$nmrk <- rep(7900,10)
genome$mpos <- rep("even", 10)
genome$nqtl <- c(171,130,153,112,129,66,60,90,63,106)
genome$qpos <- rep("even", 10)


# Marco de datos de Selección
#############################

Selection <- data.frame(matrix(NA, nrow = 2, ncol = 3))
names(Selection) <- c("Number", "type", "value")
Selection$Number[1:2] <- c(140, 140)

input <- list(genome = genome, Selection = Selection)
rm(list = c("genome", "Selection"))


# Funcion que simulará una población según argumentos indicados
###############################################################

Pob <- function(s1, v1,  hpsize, ng, h2, d2, phen_var, mutr, laf, input) {
  genome <- input$genome
  Selection <- input$Selection
  # Poblacion Historica
  historical <- make_hp(
    hpsize = hpsize, ng = ng, h2 = h2, d2 = d2, phen_var = phen_var, 
    genome = genome, mutr = mutr, laf = laf
  )
  # Poblacion Simulada
  Breed_A_Male_fndrs <- data.frame(number = Selection$Number[1], select = s1, value = v1)
  Breed_A_Female_fndrs <- data.frame(number = Selection$Number[2], select = s1, value = v1)
  Selection$type[1:2] <- c(s1, s1)
  Selection$value[1:2] <- c(v1, v1)
  Breed_A <- sample_hp(
    hp_out = historical, Male_founders = Breed_A_Male_fndrs,
    Female_founders = Breed_A_Female_fndrs,
    ng = 5, Selection = Selection,
    litter_size = 3, Display = TRUE
  )
  P<-Breed_A$output[[6]]$data$phen
  A <- (Breed_A$output[[6]]$sequ)
  A <- A[, -c(1, 2)]
  B <- matrix(rep(0, nrow(A) * ncol(A) / 2), nrow = nrow(A), ncol(A) / 2)
  for (j in 1:nrow(A)) {
    for (i in 1:(ncol(A) / 2)) {
      B[j, i] <- as.numeric(paste0(A[j, c(i, i + 1)], collapse = ""))
    }
  }
  
  BSimu<- cbind(P, B)
  BSimu
}

# Especificación de los argumentos
##################################

Sim <- function(i, input) {
  cat(i, "\n")
  Sim1<-Pob(s1='phen', v1="h", hpsize=1000, ng=200, h2=0.2, d2=0.1, phen_var=10000, mutr=0.001, laf=0.5, input=input)
  Sim2<-Pob(s1='phen', v1="l", hpsize=1000, ng=200, h2=0.2, d2=0.1, phen_var=10000, mutr=0.001, laf=0.5, input=input)
  
  
  VA<-factor(c(rep(1, nrow(Sim1)), rep(2, nrow(Sim2))))
  BaseP <- rbind(Sim1, Sim2)
  
  
  BCodif<-t(apply(BaseP[,-c(1:2)], 1, Codif))
  colnames(BCodif)<-c(paste0("snp", c(1:ncol(BCodif))))
  Base<-data.frame(Sample=c(1:nrow(BCodif)), Pop=VA, Ploydi=c(rep(2,nrow(BCodif))), Format=(rep("BiA",nrow(BCodif))), BCodif)
  BaseFreq<-stamppConvert(Base,"r")
  fst<-round(stamppFst(BaseFreq, 100, 95, 2)$Fsts, digits = 4)
  fstMin<-min(fst, na.rm=TRUE )
  fstMax<-max(fst, na.rm=TRUE)
  fstProm<-round(mean(fst, na.rm = TRUE), digits=4)
  fstDE<-round(sd(fst, na.rm = TRUE), digits=4)
  
  Base <- cbind(VA, BaseP)
  write.table(Base, paste0(paste(i, "Simulacion_K2","FstMax",fstMax,"FstMin",fstMin, "FstProm",fstProm, "FstDe",fstDE, sep = "_"), ".txt"),
              sep = "\t", eol = "\n", dec = ".", row.names = F, col.names = F
  )
  write.table(fst, paste0(paste(i, "fst", sep = "_"), ".txt"))
}

# Paralelización
################

library("parallel")
Tiempo <- system.time({
  Simul<-mclapply(sim_start:sim_end, function(i,input2) {
    Sim(i, input = input2)
  },input, mc.cores = n_cores)
})
Tiempo
