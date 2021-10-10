#Define groups among covariates
# g1={D1, D2, D1B, D2B}
# g2={B, D1B, D2B，I1B, I2B}
# g3={D1B, D2B}
# g4={I1, I2, I1B, I2B}
# g5={I1B, I2B}



# g1={1, 2, 6, 7}
# g2={5, 6, 7, 8, 9}
# g3={6, 7}
# g4={3, 4, 8, 9}
# g5={8, 9}



g=5
eta_g = as.vector(rep(1,g),mode='double') # weights of the penalties
groups = as(matrix(as.vector(c(0, 0, 0, 0, 0, 
                               0, 0, 0, 0, 0,
                               1, 1, 0, 0, 0,
                               0, 0, 0, 0, 0,
                               0, 1, 0, 1, 0),mode='logical'),ncol = g,byrow = T),'CsparseMatrix')

groups_var = as(matrix(as.vector(c(1, 0, 0, 0, 0,
                                   1, 0, 0, 0, 0,
                                   0, 0, 0, 1, 0,
                                   0, 0, 0, 1, 0, 
                                   0, 1, 0, 0, 0, 
                                   0, 0, 1, 0, 0, 
                                   0, 0, 1, 0, 0, 
                                   0, 0, 0, 0, 1, 
                                   0, 0, 0, 0, 1
),
mode='logical'),ncol = g,byrow = T),'CsparseMatrix')
graph = list('eta_g'= eta_g,'groups' = groups,'groups_var' = groups_var)
regul='graph'
#lambda=1
num_threads = -1 # all cores (-1 by default)
verbose = FALSE   # verbosity, false by default
pos = FALSE       # can be used with all the other regularizations
intercept = FALSE # can be used with all the other regularizations 



# # simplest case
# # #Define groups in covariates
# # g1={D1, D2, D1B, D2B}
# # g2={I1, I2, B, I1B, I2B}
# # g3={X1}
# # g4={X2}
# 
# 
# # g1={1, 2, 6, 7}
# # g2={3, 4, 5, 8, 9}
# # g3={10}
# # g4={11}
# g=4
# eta_g = as.vector(rep(1,g),mode='double') # weights of the penalties
# groups = as(matrix(as.vector(c(0, 0, 0, 0, 
#                                0, 0, 0, 0, 
#                                0, 0, 0, 0,
#                                0, 0, 0, 0
# ),mode='logical'),ncol = g,byrow = T),'CsparseMatrix')
# groups_var = as(matrix(as.vector(c(1, 0, 0, 0,
#                                    1, 0, 0, 0,  
#                                    0, 1, 0, 0, 
#                                    0, 1, 0, 0,  
#                                    0, 1, 0, 0, 
#                                    1, 0, 0, 0, 
#                                    1, 0, 0, 0, 
#                                    0, 1, 0, 0, 
#                                    0, 1, 0, 0, 
#                                    0, 0, 0, 1,  
#                                    0, 0, 0, 1
# ),
# mode='logical'),ncol = g,byrow = T),'CsparseMatrix')
# graph = list('eta_g'= eta_g,'groups' = groups,'groups_var' = groups_var)







