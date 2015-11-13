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

