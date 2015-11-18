#######################################
#Simple iterative RDA after Fornell (1998)
#
#Fornell, Claes, Donald W. Barclay, and Byong-Duk Rhee. 
#"A model and simple iterative algorithm for redundancy analysis." 
#Multivariate behavioral research 23.3 (1988): 349-360.


########################
####Hardcoded matrices####
########################

#INPUT DATA X'X; 
XXMA = matrix( c(1.000, 0.800, 0.140, 0.060,
                 0.800, 1.000, 0.060, 0.140,
                 0.140, 0.060, 1.000, 0.800,
                 0.060, 0.140, 0.800, 1.000), nrow = 4, ncol = 4, byrow = TRUE)

#INPUT DATA X'Y; 
XYMA = matrix( c(-0.003, 0.265, 0.404, 0.723,
                 0.062, 0.203, 0.709, 0.461,
                 0.422, 0.714, -0.142, -0.012,
                 0.710, 0.440, 0.089,-0.037), nrow = 4, ncol = 4, byrow = TRUE)

#INPUT DATA Y'Y;
YYMA = matrix( c(-0.003, 0.265, 0.404, 0.723,
                 0.062, 0.203, 0.709, 0.461,
                 0.422, 0.714, -0.142, -0.012,
                 0.710, 0.440, 0.089,-0.037), nrow = 4, ncol = 4, byrow = TRUE)

# INITIAL VALUES OF;ALPH, BETA, & CRT;
ALPH = matrix( c(1.000, 1.000, 1.000, 1.000), nrow = 4, ncol = 1, byrow = TRUE)
BETA = matrix( c(1.000, 1.000, 1.000, 1.000), nrow = 4, ncol = 1, byrow = TRUE)
CRT = matrix( c(1.000), nrow = 1, ncol = 1, byrow = TRUE)



#START Partial least square regression
while(CRT > 0.0000000001){
  
  
  ETA = XYMA %*% BETA %*% (solve(sqrt(t(BETA) %*% XYMA %*% BETA)))
  
  AP = solve(XXMA) %*% ETA;
  AQ = AP %*% solve(sqrt(t(AP) %*% XXMA %*% AP));
  
  BT = t(XYMA) %*% AQ;
  
  #Sum of squares of all elements
  CRT = sum((ALPH - AQ)^2, (BETA - BT)^2);
  
  ALPH = AQ;
  
  BETA = BT;
  
}

BETA
CRT


########################################################################################
##########################################AS A FUNCTION########################
#########################################################################################

simple_RDA <- function(Y,X) {
  
  ###########################  
  #1. Preperation of the data
  ###########################  
  
  Y.mat <- as.matrix(Y)
  Yc <- scale(Y.mat, scale = FALSE)
  
  X.mat <- as.matrix(X)
  Xcr <- scale(X.mat)
  
  
  n <- nrow(Yc)
  
  
  #covariance matrices
  XXMA = (1/(n-1)) * t(Xcr) %*% Xcr
  XYMA = (1/(n-1)) * t(Yc) %*% Xcr
  YYMA = (1/(n-1)) * t(Yc) %*% Yc
  
  
  # INITIAL VALUES OF;ALPH, BETA, & CRT;
  
  initials = rep(1,6)
  
  # an x-weight column vector of p elements
  ALPH = matrix( rep(1, ncol(X)), nrow = ncol(X), ncol = 1, byrow = TRUE)
  
  # an y-weight column vector of q elements
  BETA = matrix( rep(1, ncol(Y)), nrow = ncol(Y), ncol = 1, byrow = TRUE)
  
  
  CRT = matrix( c(1.000), nrow = 1, ncol = 1, byrow = TRUE)
  
  while(CRT > 0.0000000001){
    
    
    ETA = XYMA %*% BETA %*% (solve(sqrt(t(BETA) %*% XYMA %*% BETA)))
    
    AP = solve(XXMA) %*% ETA;
    AQ = AP %*% solve(sqrt(t(AP) %*% XXMA %*% AP));
    
    BT = t(XYMA) %*% AQ;
    
    #Sum of squares of all elements
    CRT = sum((ALPH - AQ)^2, (BETA - BT)^2);
    
    ALPH = AQ;
    
    BETA = BT;
    
  }
  
  Yhat = Xcr %*% BETA
  
  S <- cov(Yhat)
  
  # Eigenvalue decomposition of covariance matrix of Yhat
  eigenS <- eigen(S)
  
  #Number of canonical axes
  ka <- length(which(eigenS$values > 0.0000000001))
  
  #Eigenvalues of canonical axes
  ev <- eigenS$values[1:ka]
  
  print(BETA)
  print(CRT)
  
  print(ev)
  
  
}

book_data <- read.csv("~/R/DATA/table10_5.csv", header = TRUE)

book_data.m <- as.matrix(book_data)

Y <- book_data.m[,2:6]
X <- book_data.m[,7:9]

simple_RDA(Y,X)

library(vegan)
vegan_rda <- rda(Y,X)
vegan_rda

########################################################################################
##########################################CONSTRUCTION SITE########################
#########################################################################################

book_data <- read.csv("~/R/DATA/table10_5.csv", header = TRUE)

book_data.m <- as.matrix(book_data)

Y <- book_data.m[,2:6]
X <- book_data.m[,7:9]

Y.mat <- as.matrix(Y)
Yc <- scale(Y.mat, scale = FALSE)

X.mat <- as.matrix(X)
Xcr <- scale(X.mat)


n <- nrow(Yc)
n
p <- ncol(Xcr)
p
q <- ncol(Yc)
q

#covariance matrices
XXMA = (1/(n-1)) * t(Xcr) %*% Xcr
XYMA = (1/(n-1)) * t(Xcr) %*% Yc
YYMA = (1/(n-1)) * t(Yc) %*% Yc

XXMA
XYMA
YYMA

# an Y-weight column vector of p elements
BETA = matrix( rep(1, ncol(Y)), nrow = ncol(Y), ncol = 1, byrow = TRUE)

BETA
XYMA

ETA = XYMA %*% BETA %*% (solve(sqrt(t(BETA) %*% XYMA %*% BETA)))
