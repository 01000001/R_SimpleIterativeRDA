#######################################
#Simple iterative RDA after Fornell (1998)
#
#Fornell, Claes, Donald W. Barclay, and Byong-Duk Rhee. 
#"A model and simple iterative algorithm for redundancy analysis." 
#Multivariate behavioral research 23.3 (1988): 349-360.

###Denotation from Fornell (1998)####
# XI:                 an X variate
# ETA:                an Y variate
# x_t, t = 1,...p:    x variables
# y_i, i = 1,...q:    y variables
# ALPH:               an x-weight column vector of p elements
# BETA:               a y-weight column vecotr of q elements
# Lambda:             an x-loading row vector of p elements
# v:                  an y-loading row vector of q elements
# mu:                 a residual
# Epsilon:            an error term column vector of q elements
# b:                  correlation between XI and ETA

####Iterative algoright described from Fornell (1998)####
# The iterative algorithm proceeds as follows:
#
#   For initalization, let:
#
#         ETA^(0) = y_1 + y_2 + .... + y_q = Y*BETA_hat^(0) ( = Y_i)
#           were BETA^(0) = [1,1....,1]
#           and ^(i) denotes the ith iteration
#         
#         For the value of ahat^(0) (ALPH) and hence XIhat^(0), regress ETAhat^(0) on X jointly to get ALPH,
#
#           ETAhat^(0) = a_1*x_1 + a_2*x_2 + .... a_p*x_p + mu (residual)
#
#
#         and then compute the value of XI^(1):           
#           XIhat^(1) = SUM_running to p where t=1 ( ahat_t^(0) *x_t )
#
#         and normalize XIhat^(1) such that
#           t(XIhat^(1))*XIhat^(1) = t(ahat^(0)) t(X)*X*ahat^(0) = 1
#
#         For the value BETAhat^(1) and hence ETAhat^(1), regress y1,y2 ... yq separately on XIhat^(1),
#
#           y1 = BETA_1 * XIhat^(1) + Epsilon_1
#           .
#           .
#           . 
#           .
#           .
#           yq = BETA_q * XIhat^(1) + Epsilon_q
#
#         and then compute the vaule of ETAhat^(1),
#           
#           ETAhat^(1) = SUM_running to q where k=1 ( BETAhat_k^(1) * y_k)
#
#         and compare XIhat^(1) such that,
#           
#           Compare XIhat and XI, and ETAhat and ETA. If XIhat and ETAhat are not equal to XI and ETA
#           (within some chosen convergence criteria), repeat the process for the next values of XIhat
#           and ETAhat with the current value of ETAhat instead of the initilized value.
#
#           If the differences between XIhat(i) and XIhat(i-1) and ETAhat(i) and ETAhat(i-1) are within
#           the chosen convergence criteria, regress XIhat(i)  and ETAhat(i) for the value of b.


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
YYMA = matrix( c(1.000, 0.400, 0.200, 0.000,
                 0.400, 1.000, 0.000, 0.200,
                 0.200, 0.000, 1.000, 0.400,
                 0.000, 0.200, 0.400, 1.000), nrow = 4, ncol = 4, byrow = TRUE)

p <- ncol(XXMA)
q <- ncol(YYMA)

# an x-weight column vector of p elements
ALPH = matrix( rep(1, p), nrow = p, ncol = 1, byrow = TRUE)
               
# an y-weight column vector of q elements
BETA = matrix( rep(1, q), nrow = q, ncol = 1, byrow = TRUE)

CRT = matrix( c(1.000), nrow = 1, ncol = 1, byrow = TRUE)


ETA = XYMA %*% BETA %*% (solve(sqrt(t(BETA) %*% YYMA %*% BETA)))

ETA

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

#ALPH = X-WEIGHTS
ALPH

#BETA = Y-CROSSLOADINGS
BETA

#XLDG = X-LOADINGS
XLDG = XXMA %*% ALPH
XLDG

#YWGT = Y-WEIGHTS
YWGT = BETA %*% solve(sqrt(t(BETA) %*% YYMA %*% BETA))
YWGT

#YLDG = Y-LOADINGS
YLDG = YYMA %*% YWGT
YLDG

#PATH = PATH COEFFICIENT / correlation between XI and ETA
PATH = t(ALPH) %*% XYMA %*% YWGT;
PATH

#REDN = REDUNDANCY INDEX: redundancy of y variables
REDN = (sum(BETA^2))/4
REDN
