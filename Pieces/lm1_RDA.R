#***************************************************************************************#
#                           Construction site2                                            ####
#***************************************************************************************#

###DO RDA converge ALPH and BETA##########
library(elasticnet)
library(mvtnorm)

lm2_RDA <- function(X,Y,lambda,nonzero) {
  ###########################  
  #1. Preperation of the data
  ###########################  
  
  Y.mat <- as.matrix(Y)
  Yc <- scale(Y.mat, scale = FALSE)
  
  X.mat <- as.matrix(X)
  Xcr <- scale(X.mat)
  
  
  p <- ncol(Xcr)
  q <- ncol(Yc)
  
  CRT = matrix( c(1.000), nrow = 1, ncol = 1, byrow = TRUE)
  
  #   For initalization, let:
  #
  #         ETA^(0) = y_1 + y_2 + .... + y_q = Y*BETA_hat^(0) ( = Y_i)
  #           were BETA^(0) = [1,1....,1]
  # an y-weight column vector of q elements
  BETA = matrix( rep(1, q), nrow = q, ncol = 1, byrow = TRUE)
  # an x-weight column vector of p elements
  ALPH = matrix( rep(1, p), nrow = p, ncol = 1, byrow = TRUE)
  
  CRT = 1
  
  ##Peanlization algorithm
  calculateVectorEnet = function(x,y, lambda, nonzero)
  {
    colnames(x) = paste("x", 1:ncol(x),sep=".")
    epsilon.enet = enet(x, y ,lambda=lambda, max.steps=nonzero)[[4]]
    
    #reorder??
    tmp = rep(0,ncol(x))
    
    #Empty vector of regression coefficients
    tmp[colnames(x) %in% colnames(epsilon.enet)] = epsilon.enet[nonzero,]
    
    names(tmp) <- colnames(Xcr)
    
    tmp
  }
  
  #LOOP
  while(CRT > 0.00000000000001) {
    # and ^(i) denotes the ith iteration
    
    #BETA
    
    ETA = Yc %*% BETA
    XI = Xcr %*% ALPH
    #ETA
    #XI
    
    # For the value of ahat^(0) (ALPHA) and hence XIhat^(0), regress ETAhat^(0) on X jointly to get ALPH,
    
    # lm without penalization
    #ALPH_0 = as.matrix(lm(ETA ~ 0+Xcr)$coefficients)
    
    ALPH_0 = as.matrix(calculateVectorEnet(Xcr,ETA,lambda,nonzero))
    
    #ALPH_0
    #ALPH
    
    #         and then compute the value of XI^(1):           
    #           XIhat^(1) = SUM_running to p where t=1 ( ahat_t^(0) *x_t )
    XI = Xcr %*% ALPH_0
    #XI
    #XI
    
    #         and normalize XIhat^(1) such that
    #           t(XIhat^(1))*XIhat^(1) = t(ahat^(0)) t(X)*X*ahat^(0) = 1
    #           that is it variance is 1
    XI = scale(XI)
    
    #is t(ahat^(0)) t(X)*X*ahat^(0) = 1 /normalized?
    t(ALPH_0) %*%  t(Xcr) %*% Xcr %*% ALPH_0
    
    #normalize ALPH
    #ALPH4 = as.numeric(ALPH) / (sqrt(t(ALPH) %*%  t(Xcr) %*% Xcr %*% ALPH )) //<- this is not the correct form
    ALPH4_0 = as.numeric(ALPH_0) %*% solve(sqrt(t(ALPH_0) %*%  cov(Xcr,Xcr) %*% ALPH_0 ))
    
    t(ALPH4_0) %*%  t(Xcr) %*% Xcr %*% ALPH4_0
    
    ALPH_0 = as.matrix(ALPH4_0)
    
    #         For the value BETAhat^(1) and hence ETAhat^(1), regress y1,y2 ... yq separately on XIhat^(1),
    #
    #           y1 = BETA_1 * XIhat^(1) + Epsilon_1
    #           .
    
    #           yq = BETA_q * XIhat^(1) + Epsilon_q
    BETA_0 = solve(t(XI)%*%XI) %*% t(XI) %*% Yc
    BETA_0 = t(as.matrix(BETA_0))
    
    #normalize BETA
    BETA_0 = as.numeric(BETA_0) %*% solve(sqrt(t(BETA_0) %*%  cov(Yc,Yc) %*% BETA_0 ))
    BETA_0 = as.matrix(BETA_0)
    
    #BETA
    #         and then compute the vaule of ETAhat^(1),
    #           
    #           ETAhat^(1) = SUM_running to q where k=1 ( BETAhat_k^(1) * y_k)
    ETA = Yc %*% BETA_0
    
    CRT = sum((ALPH - ALPH_0)^2, (BETA - BETA_0)^2);
    #print(CRT)
    
    ALPH = ALPH_0
    BETA = BETA_0
    
  }
  
  #covariance matrices
  
  n = nrow(Yc)
  
#  XXMA = (1/(n-1)) * t(Xcr) %*% Xcr
#  XYMA = (1/(n-1)) * t(Xcr) %*% Yc
#  YYMA = (1/(n-1)) * t(Yc) %*% Yc
  
  #XLDG = X-LOADINGS / eigenvectors / site scores
#  XLDG = XXMA %*% ALPH
  
#  t(XLDG) %*% XLDG
  
  #YWGT = Y-WEIGHTS / ginverse
#  YWGT = BETA %*% solve(sqrt(t(BETA) %*% YYMA %*% BETA))
  
  #YLDG = Y-LOADINGS
#  YLDG = YYMA %*% YWGT
  
  #PATH = PATH COEFFICIENT / correlation between XI and ETA
#  PATH = t(ALPH) %*% XYMA %*% YWGT;
  
  #REDN = REDUNDANCY INDEX: redundancy of y variables
#  REDN = (sum(BETA^2))/4
  
  #************************
  
#  t(BETA) %*% t(XYMA) %*% solve(XXMA) %*% XYMA %*% BETA
  
#  PATH
  
  
  #***************************
  
result <- list(ALPH,       #ALPH = X-WEIGHTS
               BETA    #BETA = Y-CROSSLOADINGS
            
)

names(result) <- c("ALPH", 
                   "BETA"
                   
)

result
  
}

#***************************************************************************************#
#                           RUN ALGORITHM                                            ####
#***************************************************************************************#


book_data <- read.csv("~/R/DATA/table10_5.csv", header = TRUE)

book_data.m <- as.matrix(book_data)

Y <- book_data.m[,2:6]
X <- book_data.m[,7:9]

X
Y

lm2results <- lm2_RDA(X,Y,0.1,3)

lm2results

lm2results$YLDG_Y_LOADINGS
lm2results$XLDG_X_LOADINGS_site_scores

f <- lm(results$ALPH ~lm2results$ALPH)
f

library(vegan)

veganRDA <- rda(Y,X)

veganRDA$CCA$biplot
summary(veganRDA)

plot(veganRDA$CCA$biplot)

veganRDA$CA$v

veganRDA$CCA$v
sum(abs(veganRDA$CCA$wa[,1]))/nrow(veganRDA$CCA$wa)

#***************************************************************************************#
#######Play around data##########
#***************************************************************************************#
n = 400
hoogCor = 50  
lageCor = 10

covMatrix = matrix(rep(hoogCor ,8100), ncol=90)

diag(covMatrix ) = 100

hogeCor = rmvnorm(n, rep(0,90), covMatrix )

covMatrix = matrix(rep(lageCor ,400), ncol=20)
diag(covMatrix ) = 100
x1=x2=c()

for(i in 1:20)
  x1 = cbind(x1, rmvnorm(n, rep(0,20), covMatrix ))

for(i in 1:50)
  x2 = cbind(x2, rmvnorm(n, rep(0,20), covMatrix ))


x1[,1:10] = hogeCor[,1:10]
x2[,1:10] = hogeCor[,11:20]

dim(x1);dim(x2)

lm2results <- lm2_RDA(x1,x2,0.1,100)

lm2results$ALPH
lm2results$BETA
plot(lm2results$BETA)
plot(lm2results$ALPH)

#***************************************************************************************#
#######DELETE ME##########
#***************************************************************************************#
#***************************************
"""
    x = Xcr
    y = ETA
    lambda = 0.5
    nonzero = 2
    
    colnames(x) = paste("x", 1:ncol(x),sep=".")
    colnames(x)
    epsilon.enet = enet(x, y ,lambda=lambda, max.steps=nonzero)[[4]]
    epsilon.enet
    
    #reorder??
    tmp = rep(0,ncol(x))
    
    #Empty vector of regression coefficients
    tmp[colnames(x) %in% colnames(epsilon.enet)] = epsilon.enet[nonzero,]
    
    names(tmp) <- colnames(Xcr)
    
    tmp
    """
#***************************************

M = matrix( c(1, 2, 3, 4,
                 5, 6, 7, 8,
                 9, 10, 11, 12,
                 13, 14, 15,16), nrow = 4, ncol = 4, byrow = TRUE)
M

scale(M)

M = matrix( c(1, 2
            ),nrow = 2, ncol = 1, byrow = TRUE)
M
M = scale(M)

mean(M)

sd(M)


result <- list(ALPH,       #ALPH = X-WEIGHTS
               BETA,    #BETA = Y-CROSSLOADINGS
               XLDG,    #XLDG = X-LOADINGS, site scores = Biplot scores for constraining variables"
               YWGT,    #YWGT = Y-WEIGHTS
               YLDG,    #"YLDG = Y-LOADINGS"
               PATH,    #"PATH = PATH COEFFICIENT / correlation between XI and ETA"
               REDN    #REDN = REDUNDANCY INDEX: redundancy of y variables
)

names(result) <- c("ALPH", 
                   "BETA", 
                   "XLDG_X_LOADINGS_site_scores", 
                   "YWGT_Y_WEIGHTS",
                   "YLDG_Y_LOADINGS",
                   "PATH_PATH_COEFFICIENT",
                   "REDN_REDUNDANCY_INDEX"
)

result
