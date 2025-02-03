source("code_ML-DynSBM.R")
set.seed(1234)
library(latex2exp)
library(ggplot2)
library(randnet)


k=6
T=c(1,4,9,16)
pi = rep(1/k,k)
n_values <- c(60,120,180,240,300,360,420,480,540,600,660,720,780,840,900,960,1020)
nb_essai <- 100

create_matrix_commu <- function(k) {
  a <- runif(1, min = 0, max = 0.1)
  mat <- matrix(a, nrow = k, ncol = k)
  diag(mat) <- runif(k, min = 0.7, max = 1) 
  return(mat)
}

res_KT <- c()
for (t in T){
  print(paste("T =", t))
  P_commu <- lapply(1:t, function(x) create_matrix_commu(k))
  for(n in n_values){
    print(paste("n=", n))
    bon_KT <- 0
    for (i in 1:nb_essai){
      print(paste("Essai", i))
      G <- MLsbm(n/(t^(1/2)),pi,P_commu,t)
      k_chap_KT <- k_chap_ML(G)
      print(paste("k_chapeau = ",k_chap_KT))
      if (k_chap_KT == k) { bon_KT <- bon_KT + 1 }
    }
    res_KT <- c(res_KT, bon_KT / nb_essai)
  }
}


data <- data.frame(
  n = n_values,
  bon_KT_1 = res_KT[1:17],
  bon_KT_4 = res_KT[18:34],
  bon_KT_9 = res_KT[35:51],
  bon_KT_16 = res_KT[52:68]
)

ggplot(data, aes(x = n)) +
  geom_line(aes(y = bon_KT_1, color = "1", linetype = "1"), linewidth = 1) +
  geom_line(aes(y = bon_KT_4, color = "4", linetype = "4"), linewidth = 1) +
  geom_line(aes(y = bon_KT_9, color = "9", linetype = "9"), linewidth = 1) +
  geom_line(aes(y = bon_KT_16, color = "16", linetype = "16"), linewidth = 1) +
  labs(x = TeX(r'(total number of possible edges ($n^2T$))'), y = "Accuracy") +
  scale_color_manual(name = "T", 
                     values = c("1" = "red", "4" = "blue", "9" = "green", "16" = "yellow"),
                     breaks = c("1", "4", "9", "16")) +
  scale_linetype_manual(name = "T", 
                        values = c("1" = "solid", "4" = "dashed", "9" = "dotted", "16" = "dotdash"),
                        breaks = c("1", "4", "9", "16")) +
  theme_minimal()


