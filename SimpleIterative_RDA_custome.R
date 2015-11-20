#*******************************************************************************#
#Simple iterative RDA after Fornell (1998)                                      #####
#                                                                               #
#Fornell, Claes, Donald W. Barclay, and Byong-Duk Rhee.                         #
#"A model and simple iterative algorithm for redundancy analysis."              #
#Multivariate behavioral research 23.3 (1988): 349-360.                         #
#*******************************************************************************#


#***************************************************************************************#
#                           SIMPLE RDA AS A FUNCTION                                    ####
#***************************************************************************************#

simple_RDA <- function(X,Y) {
  
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
  XYMA = (1/(n-1)) * t(Xcr) %*% Yc
  YYMA = (1/(n-1)) * t(Yc) %*% Yc
  
  p <- ncol(XXMA)
  q <- ncol(YYMA)
  
  # an x-weight column vector of p elements
  ALPH = matrix( rep(1, p), nrow = p, ncol = 1, byrow = TRUE)
  # an y-weight column vector of q elements
  BETA = matrix( rep(1, q), nrow = q, ncol = 1, byrow = TRUE)
  
  CRT = matrix( c(1.000), nrow = 1, ncol = 1, byrow = TRUE)
  
  ETA = XYMA %*% BETA %*% (solve(sqrt(t(BETA) %*% YYMA %*% BETA)))

  
  #START Partial least square regression
  while(CRT > 0.0000000001){
    
    ETA = XYMA %*% BETA %*% (solve(sqrt(t(BETA) %*% YYMA %*% BETA)))
    
    AP = solve(XXMA) %*% ETA;
    AQ = AP %*% solve(sqrt(t(AP) %*% XXMA %*% AP));
    
    BT = t(XYMA) %*% AQ;
    
    #Sum of squares of all elements
    CRT = sum((ALPH - AQ)^2, (BETA - BT)^2);
    
    ALPH = AQ;
    
    BETA = BT;
    
  }
  
  #XLDG = X-LOADINGS
  XLDG = XXMA %*% ALPH
  
  #YWGT = Y-WEIGHTS
  YWGT = BETA %*% solve(sqrt(t(BETA) %*% YYMA %*% BETA))
  
  #YLDG = Y-LOADINGS
  YLDG = YYMA %*% YWGT
  
  
  #PATH = PATH COEFFICIENT / correlation between XI and ETA
  PATH = t(ALPH) %*% XYMA %*% YWGT;
  
  #REDN = REDUNDANCY INDEX: redundancy of y variables
  REDN = (sum(BETA^2))/4
  
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
                     "YWGT_Y-WEIGHTS",
                     "YLDG_Y-LOADINGS",
                     "PATH_PATH_COEFFICIENT",
                     "REDN_REDUNDANCY_INDEX"
                      )
  result

}

book_data <- read.csv("~/R/DATA/table10_5.csv", header = TRUE)

book_data.m <- as.matrix(book_data)

Y <- book_data.m[,2:6]
X <- book_data.m[,7:9]

results <- simple_RDA(X,Y)
results$ALPH
results$BETA
results$XLDG_X_LOADINGS_site_scores

library(vegan)
vegan_rda <- rda(Y,X)
summary(vegan_rda)
vegan_rda$CCA$biplot

#***************************************************************************************#
#                           Simulate X and Y from XX, XY and YY                         ####
#***************************************************************************************#
"""

x = rnorm(1000,mean=4,2)

mean(x)
sd

library(mvtnorm)


Sigma = diag(2)
Sigma[1,2] = 0.3
Sigma[2,1] = 0.3

X = rmvnorm(100, c(0,1), sigma=Sigma)
apply(X,2,mean)
var(X)



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
YYMA = matrix( c(1.000, 0.400, 0.200, 0.000,
                 0.400, 1.000, 0.000, 0.200,
                 0.200, 0.000, 1.000, 0.400,
                 0.000, 0.200, 0.400, 1.000), nrow = 4, ncol = 4, byrow = TRUE)


SIGMA = diag(8)
SIGMA[SIGMA==1] = 0
SIGMA[1:4,1:4] = XXMA
SIGMA[5:8,5:8] = YYMA
SIGMA[5:8,1:4] = XYMA
SIGMA[1:4,5:8] = t(XYMA)

SIGMA
xy_sim = rmvnorm(1000000, rep(0,8), sigma=SIGMA )

var(X)
var(X,Y)

xy_sim[,1:4]
X = xy_sim[,1:4]
Y = xy_sim[,5:8]

simple_RDA(X,Y)
"""
