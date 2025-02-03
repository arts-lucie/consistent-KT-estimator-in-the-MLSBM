source("code_ML-DynSBM.R")
set.seed(1234)
library(latex2exp)
library(ggplot2)
library(randnet)


k=6
T=5
pi = rep(1/k,k)
n_values <- c(50,75,100,125,150,175,200,300,500,700)
nb_essai <- 100


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

res_KT <- c()
res_PML <- c()
res_BHMC <- c()
res_NCV <- c()
k_KT <- c()
P <- lapply(1:T, function(x) create_matrix(k))
for(n in n_values){
  print(paste("n=", n))
  bon_KT <- 0
  bon_PML <-0
  bon_BHMC <-0
  bon_NCV <- 0
  for (i in 1:nb_essai){
    print(paste("Essai", i))
    G <- MLsbm(n,pi,P,T)
    k_chap_KT <- k_chap_ML(G)
    print(paste("k_chapeau = ",k_chap_KT))
    k_KT <- c(k_KT, k_chap_KT)
    if (k_chap_KT == k) { bon_KT <- bon_KT + 1 }
    for(t in 1:T){
      k_chap_PML <- LRBIC(as_adjacency_matrix(G[[t]]), 15, lambda = NULL, model = "SBM")$SBM.K
      print(paste("k_chapeau PML = ",k_chap_PML))
      if (k_chap_PML == k) { bon_PML <- bon_PML + 1 }
      k_chap_BHMC <-BHMC.estimate(as_adjacency_matrix(G[[t]]), K.max = 15)$K
      print(paste("k_chapeau BHMC = ",k_chap_BHMC))
      if (k_chap_BHMC == k) { bon_BHMC <- bon_BHMC + 1 }
      k_chap_NCV_texte <-NCV.select(as_adjacency_matrix(G[[t]]), 15, cv = 3)$l2.model
      k_chap_NCV <- substr(k_chap_NCV_texte, 5, nchar(k_chap_NCV_texte))
      print(paste("k_chapeau NCV = ",k_chap_NCV))
      if (k_chap_NCV == k) { bon_NCV <- bon_NCV + 1 }
    }
  }
  res_KT <- c(res_KT, bon_KT / nb_essai)
  res_PML <- c(res_PML, bon_PML / (nb_essai*T))
  res_BHMC <- c(res_BHMC, bon_BHMC / (nb_essai*T))
  res_NCV <- c(res_NCV, bon_NCV / (nb_essai*T))
}



data <- data.frame(
  n = n_values,
  bon_KT = res_KT,
  bon_PML = res_PML,
  bon_BHMC = res_BHMC,
  bon_NCV = res_NCV
)


ggplot(data, aes(x = n)) +
  geom_line(aes(y = bon_KT, color = "KT", linetype = "KT"), linewidth = 1) +
  geom_line(aes(y = bon_PML, color = "PML", linetype = "PML"), linewidth = 1) +
  geom_line(aes(y = bon_BHMC, color = "BHMC", linetype = "BHMC"), linewidth = 1) +
  geom_line(aes(y = bon_NCV, color = "NCV", linetype = "NCV"), linewidth = 1) +
  labs(
    x = "Number of nodes (n) per layer", 
    y = "Accuracy"
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






