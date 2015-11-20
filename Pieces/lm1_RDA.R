
#***************************************************************************************#
#                           Construction site                                            ####
#***************************************************************************************#

book_data <- read.csv("~/R/DATA/table10_5.csv", header = TRUE)

book_data.m <- as.matrix(book_data)

Y <- book_data.m[,2:6]
X <- book_data.m[,7:9]

simple_RDA(X,Y)
lm_RDA(X,Y)

###DO RDA##########


lm_RDA <- function(X,Y) {
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
  
  #LOOP
  while(CRT > 0.0000000001){
    #           and ^(i) denotes the ith iteration
    
    ETA = Yc %*% BETA
    ETA
    
    #         For the value of ahat^(0) (ALPHA) and hence XIhat^(0), regress ETAhat^(0) on X jointly to get ALPH,
    
    ALPH_0 = lm(ETA ~ 0+Xcr)
    
    ALPH_0 = ALPH_0$coefficients
    ALPH_0 = as.matrix(ALPH_0)
    ALPH_0
    
    #         and then compute the value of XI^(1):           
    #           XIhat^(1) = SUM_running to p where t=1 ( ahat_t^(0) *x_t )
    XI = Xcr %*% ALPH_0
    
    #         and normalize XIhat^(1) such that
    #           t(XIhat^(1))*XIhat^(1) = t(ahat^(0)) t(X)*X*ahat^(0) = 1
    #           that is it variance is 1
    XI = XI / sqrt(sum(XI^2))
    t(XI) %*% XI
    #SCALE ALPHA?
    
      
    #         For the value BETAhat^(1) and hence ETAhat^(1), regress y1,y2 ... yq separately on XIhat^(1),
    #
    #           y1 = BETA_1 * XIhat^(1) + Epsilon_1
    #           .
    
    #           yq = BETA_q * XIhat^(1) + Epsilon_q
    BETA_0 = solve(t(XI)%*%XI) %*% t(XI) %*% Yc
    BETA_0 = t(as.matrix(BETA_0))
    BETA_0
    #         and then compute the vaule of ETAhat^(1),
    #           
    #           ETAhat^(1) = SUM_running to q where k=1 ( BETAhat_k^(1) * y_k)
    ETA = Yc %*% BETA_0
    
    ALPH_0  
    BETA
    BETA_0
    
    CRT = sum((ALPH - ALPH_0)^2, (BETA - BETA_0)^2);
    
    ALPH = ALPH_0
    
  
    
    BETA = BETA_0
  }
  
  
  ALPH4 = as.numeric(ALPH) /   sqrt(t(ALPH) %*%  t(Xcr) %*% Xcr %*% ALPH )

  t(ALPH4) %*%  t(Xcr) %*% Xcr %*% ALPH4
   
  ALPH = ALPH4
  
  
  
  #covariance matrices
  
  n = nrow(Yc)
  
  XXMA = (1/(n-1)) * t(Xcr) %*% Xcr
  XYMA = (1/(n-1)) * t(Xcr) %*% Yc
  YYMA = (1/(n-1)) * t(Yc) %*% Yc

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
