source("code_ML-DynSBM.R")
set.seed(1234)
library(latex2exp)
library(ggplot2)
library(randnet)


k=6
T=5
pi = rep(1/k,k)
n_values <- c(50,75,100,125,150,175,200,300,500,1000,1500,2000)
nb_essai <- 10

create_matrix_commu <- function(k) {
  a <- 0.4 
  mat <- matrix(a, nrow = k, ncol = k)
  diag(mat) <- runif(k, min = 0.6, max = 1) 
  return(mat)
}

create_matrix_anti_commu <- function(k) {
  a <- runif(1,min=0.6,max=1) 
  mat <- matrix(a, nrow = k, ncol = k)
  diag(mat) <-0.4 
  return(mat)
}

create_matrix <- function(k) {
  prem <- create_matrix_commu(k/2)
  deu <- create_matrix_anti_commu(k/2)
  reste <- matrix(0.2, nrow = k/2, ncol = k/2) 
  mat <- rbind(cbind(prem, reste), cbind(reste, deu))
  return(mat)
}

temps_KT <- c()
temps_PML <- c()
temps_BHMC <- c()
temps_NCV <- c()
k_KT <- c()
k_PML <- c()
k_BHMC <- c()
k_NCV <- c()
P <- lapply(1:T, function(x) create_matrix(k)) 
for(n in n_values){
  print(paste("n=", n))
  t_KT <- 0
  t_PML <- 0
  t_BHMC <- 0
  t_NCV <- 0
  for(i in 1:nb_essai){
    print(paste("essai ", i))
    G <- MLsbm(n,pi,P,T)
    temps_exe_KT <- system.time({
      k_chap_KT <- k_chap_ML(G)
      print(paste("k_chapeau KT = ",k_chap_KT))
    })
    k_KT <- c(k_KT, k_chap_KT)
    t_KT <- t_KT+ temps_exe_KT["elapsed"]
    temps_exe_PML <- system.time({
      k_chap_PML <- maj(lapply(1:T, function(t) LRBIC(as_adjacency_matrix(G[[t]]), 15, lambda = NULL, model = "SBM")$SBM.K))
      print(paste("k_chapeau PML = ",k_chap_PML))
    })
    k_PML <- c(k_PML, k_chap_PML)
    t_PML <- t_PML+ temps_exe_PML["elapsed"]
    temps_exe_BHMC <- system.time({
      k_chap_BHMC <-maj(lapply(1:T, function(t) BHMC.estimate(as_adjacency_matrix(G[[t]]), K.max = 15)$K))
      print(paste("k_chapeau BHMC = ",k_chap_BHMC))
    })
    k_BHMC <- c(k_BHMC, k_chap_BHMC)
    t_BHMC <-t_BHMC+ temps_exe_BHMC["elapsed"]
    temps_exe_NCV <- system.time({
      k_chap_NCV_texte <-maj(lapply(1:T, function(t) NCV.select(as_adjacency_matrix(G[[t]]), 15, cv = 3)$l2.model))
      k_chap_NCV <- substr(k_chap_NCV_texte, 5, nchar(k_chap_NCV_texte))
      print(paste("k_chapeau NCV = ",k_chap_NCV))
    })
    k_NCV <- c(k_NCV, k_chap_NCV)
    t_NCV <- t_NCV + temps_exe_NCV["elapsed"]
  }
  temps_KT <- c(temps_KT, t_KT/nb_essai)
  temps_PML <- c(temps_PML, t_PML/nb_essai)
  temps_BHMC <- c(temps_BHMC, t_BHMC/nb_essai)
  temps_NCV <- c(temps_NCV, t_NCV/nb_essai)
}


data <- data.frame(
  n = n_values,
  time_KT = temps_KT,
  time_PML = temps_PML,
  time_BHMC = temps_BHMC,
  time_NCV = temps_NCV
)

ggplot(data, aes(x = n)) +
  geom_line(aes(y = time_KT, color = "KT", linetype = "KT"), linewidth = 1) +
  geom_line(aes(y = time_PML, color = "PML", linetype = "PML"), linewidth = 1) +
  geom_line(aes(y = time_BHMC, color = "BHMC", linetype = "BHMC"), linewidth = 1) +
  geom_line(aes(y = time_NCV, color = "NCV", linetype = "NCV"), linewidth = 1) +
  labs(
    x = "Number of nodes (n) per layer", 
    y = "Execution time (seconds)"
  ) +
  scale_color_manual(
    name = "Method", 
    values = c("KT" = "red", "PML" = "blue", "BHMC" = "green", "NCV" = "yellow")
  ) +
  scale_linetype_manual(
    name = "Method",
    values = c("KT" = "solid", "PML" = "dashed", "BHMC" = "dotted", "NCV" = "dotdash")
  ) +
  theme_minimal()

print(data)




