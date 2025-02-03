source("code_ML-DynSBM.R")
set.seed(1234)
library(latex2exp)
library(ggplot2)
library(randnet)


T=4
k=3
pi = rep(1/k,k)
n <- 300
nb_essai <- 100


rho <- seq(0.05,0.45,by=0.05)

S <- lapply(1:T, function(t) {
  a <-runif(1, min = -0.1, max = 0.1)
  M <- matrix(1,k,k)
  diag(M) <- rep(2,k)
  M <- M+a
  return(M)})

res_KT <- c()
res_PML <- c()
res_BHMC <- c()
res_NCV <- c()


for (r in rho){
  print(paste("rho =", r))
  bon_KT <- 0
  bon_PML <-0
  bon_BHMC <-0
  bon_NCV <- 0
  for (i in 1:nb_essai){
    print(paste("Essai", i))
    rS <- lapply(1:T, function(t) r*S[[t]])
    G <- MLsbm(n,pi,rS,T)
    k_chap_KT <- k_chap_ML(G)
    print(paste("k_chapeau = ",k_chap_KT))
    if (k_chap_KT == k) { bon_KT <- bon_KT + 1 }
    for(t in 1:T){
      k_chap_PML <- LRBIC(as_adjacency_matrix(G[[t]]), 15, lambda = NULL, model = "SBM")$SBM.K
      if (k_chap_PML == k) { bon_PML <- bon_PML + 1 }
      k_chap_BHMC <-BHMC.estimate(as_adjacency_matrix(G[[t]]), K.max = 15)$K
      if (k_chap_BHMC == k) { bon_BHMC <- bon_BHMC + 1 }
      k_chap_NCV_texte <-NCV.select(as_adjacency_matrix(G[[t]]), 15, cv = 3)$l2.model
      k_chap_NCV <- substr(k_chap_NCV_texte, 5, nchar(k_chap_NCV_texte))
      if (k_chap_NCV == k) { bon_NCV <- bon_NCV + 1 }
    }
  }
  res_KT <- c(res_KT, bon_KT / nb_essai)
  res_PML <- c(res_PML, bon_PML / (nb_essai*T))
  res_BHMC <- c(res_BHMC, bon_BHMC / (nb_essai*T))
  res_NCV <- c(res_NCV, bon_NCV / (nb_essai*T))
}




data <- data.frame(
  r = rho,
  bon_KT = res_KT,
  bon_PML = res_PML,
  bon_BHMC = res_BHMC,
  bon_NCV = res_NCV
)


print(data)










