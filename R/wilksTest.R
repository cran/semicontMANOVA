#############################################################
#	Wilks statistic
#	Authors: E. Sabbioni, C. Agostinelli and A. Farcomeni
#	Maintainer e-mail: elena.sabbioni@polito.it
#	Date: 25 September 2023
#	Version: 0.1
#	Copyright (C) 2023 E. Sabbioni, C. Agostinelli and A. Farcomeni
#############################################################

# -------- TEST ------------------------------------------------------
# Valori in iput: 
# logLik0: log-likelihood sotto H0
# logLik:  log-likelihood sotto nessuna asusnzione
# p: numeri di geni
# K: numero di gruppi
#
# Valori di output: 
# statistic: valore della statistica test
# p.value: p-value
# parameter: gradi di liberta' della chi-quadrato
  
wilksTest <- function(logLik0, logLik, p, K) {
  D <- -2 * (logLik0 - logLik)
  # pvalue, calcolato assumendo che D abbia distribuzione chi-quadrato
  pval <- 1 - pchisq(D, 2*p*(K-1))
  result <- list(statistic=D, parameter=2*p*(K - 1), p.value=pval)
  return(result)
}
