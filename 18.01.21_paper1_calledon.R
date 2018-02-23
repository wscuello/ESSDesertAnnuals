# SUMMARY: THIS IS CALLED ON BY HIREC-OCT2616-EXTESS.R
# THIS HAS BEEN PETE TRIMMER APPROVED

###################################################################################
# Purpose of this file is to simplify "everything" from HIREC-NODATA-714-2016 #####
# HIREC-REPORT-714-2016.Rnw                                                   #####
###################################################################################

#####################################
##expressions (dr/dg) and (fitness)##
#####################################

fit_func = function(Y, N, g, s, a){
  return(g*Y/(1+g*a*N) + (1-g)*s)
}

fit_phenoplast_func = function(Y,N,x,s,a,E){ # E is the environmental cues
  g = inverse_logit((x[1]+E*x[2]))
  return(g*Y/(1+g*a*N) + (1-g)*s)
}

ddg_func = function(y,n,g,s,a){ # This is the expression in the E(dr(g,gtilde)/dgtilde)
  z = y/(1+a*g*n)
  return((z-s)/(z*g+(1-g)*s))
}

ddg_phenoplast_1_func = function(y,n,x,s,a,E){
  g = inverse_logit(x[1]+E*x[2])
  z = y/(1+a*g*n)
  chain1 = (z-s)/(z*g+(1-g)*s)
  chain2 = exp(-(x[1]+E*x[2]))/((1+exp(-(x[1]+E*x[2])))^2)
  return(chain1*chain2)
}

ddg_phenoplast_2_func = function(y,n,x,s,a,E){
  g = inverse_logit(x[1]+E*x[2])
  z = y/(1+a*g*n)
  chain1 = (z-s)/(z*g+(1-g)*s)
  chain2 = E*exp(-(x[1]+E*x[2]))/((1+exp(-(x[1]+E*x[2])))^2)
  return(chain1*chain2)
}

#######################################
##end express. (dr/dg) and (fitness)###
#######################################

############             ##############
############             ##############

#######################################
## MAIN FUNCTIONS                   ###
#######################################

# Notes: Never deleted the N0 from the function -- 
# it's okay with just being 1... very useless function argument, but I 
# didn't want to go through all the functions changing it throughout.

ESS_func = function(Ys,N0,s,a){ # For autocorrelation, keep track of the environmental factor it uses
  h=function(g){mean_drdg2(Ys,N0,g,s,a)}
  ESS.val = 1
  if (h(0) < 0 && h(1) < 0){
    ESS.val = 0
  }else if( h(0)>0 && h(1)>0){ # totally forgot to put this in before
    ESS.val = 1 
  }else if (h(0) > 0 && h(1) < 0){
    out = uniroot(h, c(0,1)); ESS.val = out$root
  }else if (h(0) < 0 && h(1) > 0){
    print('Minimum Acquired')
  }
  return(ESS.val)
}


mean_drdg2 = function(Ys,N0,g,s,a){ # Outputs the mean (dr/dg) expression; this is for uniroot to use.
  Tf = length(Ys) 
  N = N0
  run_sum = 0 # this isn't actually the mean, it's the running sum
  count = 0
  cut_off = 2/4
  for (i in 1:(Tf-1)){
    if (i <= cut_off*Tf){
      Y = Ys[i+1]
      fitness = fit_func(Y,N,g,s,a) # create the fitness
      N=N*fitness
    }else{
      Y = Ys[i+1]
      expression = ddg_func(Y,N,g,s,a)
      run_sum = run_sum + expression
      fitness = fit_func(Y,N,g,s,a) # create the fitness
      N=N*fitness
      count = count+1
    }
  }
  return(run_sum/count) # sum divided by count gets sample mean
} 

############################################
# Now to some Paper 1 stuff
############################################


#############################################
# Paper 3 stuff
#############################################

ESS_func_two_fracs = function(Ys, N0, s, a, Es){ # x is the tuple (a,b) -- the phenotypic traits
  h_phenoplast = function(x){c(F1 = mean_drdg_phenoplast1(Ys,N0,x,s,a,Es), F2 = mean_drdg_phenoplast2(Ys, N0, x, s, a, Es))} # some_function is a place-holder for now
  ess_tuple = multiroot(f = h_phenoplast, start = c(.25,.25))$root
  return(ess_tuple)
}

mean_drdg_phenoplast1 = function(Ys, N0, x, s, a, Es){
  Tf = length(Ys)
  N = N0
  run_sum = 0
  count = 0
  cut_off = 2/4 # what fraction of Tf to start averaging
  for (i in 1:(Tf-1)){
    if (i <= cut_off*Tf){ # don't average, but iterate forward with the model
      Y = Ys[i+1]
      E = Es[i+1] # grab the next environmental cue 
      fitness = fit_phenoplast_func(Y,N,x,s,a,E) # create the fitness
      N=N*fitness
    }else{ # start averaging the first phenotypic plasticity term
      Y = Ys[i+1]
      E = Es[i+1]
      expression = ddg_phenoplast_1_func(Y,N,x,s,a,E)
      run_sum = run_sum + expression
      fitness = fit_phenoplast_func(Y,N,x,s,a,E) # create the fitness
      N=N*fitness
      count = count+1
    }
  }
  return(run_sum/count) # sum divided by count gets sample mean
}

mean_drdg_phenoplast2 = function(Ys, N0, x, s, a, Es){
  Tf = length(Ys)
  N = N0
  run_sum = 0
  count = 0
  cut_off = 1/10
  for (i in 1:(Tf-1)){
    if (i <= cut_off*Tf){
      Y = Ys[i+1]
      E = Es[i+1]
      fitness = fit_phenoplast_func(Y,N,x,s,a,E) # create the fitness
      N=N*fitness
    }else{
      Y = Ys[i+1]
      E = Es[i+1]
      expression = ddg_phenoplast_2_func(Y,N,x,s,a,E)
      run_sum = run_sum + expression
      fitness = fit_phenoplast_func(Y,N,x,s,a,E) # create the fitness
      N=N*fitness
      count = count+1
    }
  }
  return(run_sum/count) # sum divided by count gets sample mean
}

#############################################
# back to other stuff needed for both papers
#############################################
auto_cor_func = function(Ps, rho){
  # Note the Ps are of form exp(base*sigma + mu)
  # rho is the autocorrelation constant
  # Pdata is the output
  
  Tf = length(Ps)
  Pdata = array(NA,Tf)
  mu = mean(log(Ps))
  lPdata = log(Ps)-mu # subtracting off the mean
  Pdata[1] = exp(lPdata[1]+mu)
  V = lPdata[1]
  for (i in 1:(Tf-1)){
    Z = lPdata[i+1]
    V = rho*V+sqrt(1-rho^2)*Z
    P = exp(V+mu)
    Pdata[i+1]=P
  }
  return(Pdata)
}


Ns_func = function(Ys, N0, g, s, a){# input: Ydata, N0 -- starting pop, g -- germination rate, s -- survivor rate, a -- intraspecific competitive factor
  Tf = length(Ys)
  Ns = Tf
  Ns[1] = N0
  for (i in 1:(Tf-1)){
    Y = Ys[i+1] # initiate with the *next* environmental factor 
    fitness = fit_func(Y,Ns[i],g,s,a) # create the fitness
    # Ns[i+1]=log(Ns[i])+log(fitness)
    # Ns[i+1]=exp(Ns[i+1])
    # print(Ns[i])
    Ns[i+1] = Ns[i]*fitness
  }
  Ns_and_Ys = list(pops = Ns, envs = Ys)
  return(Ns_and_Ys)
}

####################################################################
# These are the vectorized counterparts of Ns_func and auto_cor_func
####################################################################

multi_Ns_func = function(Ymat, N0.vec, g, s, a){# Input: Matrix of environmental factors, vector of initial population(s), vector of corresponding ESS value(s), survival rate, and intraspecific competition
  t = length(Ymat[,1]) # Ys should be a matrix of environmental stuff
  kk = length(N0.vec)
  Nmat = matrix(NA, nrow = t, ncol = kk)
  Nmat[1,] = N0.vec
  for (i in 1:(t-1)){
    Y = Ymat[i+1,] # initiate with the *next* environmental factor 
    fitness = fit_func(Y,Nmat[i,],g,s,a) # create the fitness
    Nmat[i+1,]=Nmat[i,]*fitness
  }
  return(Nmat)
}


multi_auto_cor_func = function(m1, P.mat, rho){# first row of P.mat is already part of an autocorrelated sequence
  mu = m1 # this should be the THEORETICAL mean of P.mat... or should it?
  Tf = length(P.mat[,1])
  lPdata.mat = log(P.mat) - mu
  Pnew.mat = matrix(NA, nrow = Tf, ncol = dim(P.mat)[2])
  Pnew.mat[1,] = P.mat[1,] #exp(lPdata.mat[1,] + mu)
  V = lPdata.mat[1,]
  for (i in 1:(Tf-1)){
    Z = lPdata.mat[i+1,]
    V = rho*V+sqrt(1-rho^2)*Z
    P = exp(V+mu)
    Pnew.mat[i+1,]=P
  }
  return(Pnew.mat)
}

########################################
## For creating a new lognormal
########################################

Pnew_func = function(base, m1, sd1, sd2){
  m2 = (sd1^2/2)*(1-sd2^2)
  P2 = exp((m1+m2)+sd1*sd2*base)
  return(P2)
}

pos_old_stdev_new_mean_lognormal_func = function(b, sd1, m1){
  inside_argument = (1/b^2)*(exp(sd1^2)-1) + 1
  sd2 = sqrt(log(inside_argument))
  m2 = log(b)+m1 + (sd1^2-sd2^2)/2
  list_of_m2_sd2 = list(new_mean = m2, new_sd = sd2)
  return(list_of_m2_sd2)
}

neg_old_stdev_new_mean_lognormal_func = function(b, sd1, m1){
  inside_argument = (1/b^2)*(exp(sd1^2)-1) + 1
  sd2 = -sqrt(log(inside_argument))
  m2 = log(b)+m1 + (sd1^2-sd2^2)/2
  list_of_m2_sd2 = list(new_mean = m2, new_sd = sd2)
  return(list_of_m2_sd2)
}

new_stdev_old_mean_lognormal_func = function(sd_mult_factor, sd1, m1){
  inside_argument = sd_mult_factor^2*(exp(sd1^2)-1) + 1
  sd2 = sqrt(log(inside_argument))
  m2 = m1 + (sd1^2-sd2^2)/2
  list_of_m2_sd2 = list(new_mean = m2, new_sd = sd2)
  return(list_of_m2_sd2)
}

new_stdev_old_mean_lognormal_func_neg = function(sd_mult_factor, sd1, m1){
  inside_argument = sd_mult_factor^2*(exp(sd1^2)-1) + 1
  sd2 = -sqrt(log(inside_argument))
  m2 = m1 + (sd1^2-sd2^2)/2
  list_of_m2_sd2 = list(new_mean = m2, new_sd = sd2)
  return(list_of_m2_sd2)
}

new_stdev_old_mean_lognormal_func_pos = function(sd.diff, sd1, m1){
  inside_argument = sd.diff*(exp(sd1^2)-1) + 1
  sd2 = sqrt(log(inside_argument))
  m2 = m1 + (sd1^2-sd2^2)/2
  list_of_m2_sd2 = list(new_mean = m2, new_sd = sd2)
  return(list_of_m2_sd2)
}

fitness_mat_producer = function(N){ # takes in matrix of whose rows are updates of previous row N(t+1) = N(t)*f(N(t), Y(t+1))
  nrows = dim(N)[1]
  ncols = dim(N)[2]
  fitness_mat = array(NA, c(nrows, ncols))
  for (m in 2:nrows){
    fitness_mat[m,] = (N[m,])/(N[(m-1),]) 
  }
  return(fitness_mat[2:nrows,])
}

sd.real = function(x){
  m = mean(x)
  integrand = (x-m)^2
  return(sqrt(mean(integrand)))
}

dead_func = function(M, thresh){
  len_M = length(M[,1])
  for (t in 1:(len_M-1)){
    indices = which(M[t,] <= thresh)
    M[t,indices] = thresh
    M[t:(len_M), indices] = thresh
  }
  return(M)
}

zeros_and_ones_func = function(log_ydata, y_hyperplane){ # any log_data that shows you didn't replace yourself set it to 0, otherwise set it to 1
  y = array(NA, length(log_ydata))
  y[which(log_ydata <= y_hyperplane)] = 0
  y[which(log_ydata > y_hyperplane)] = 1
  return(y)
}


probability_logistic_coeffs_func = function(b0, b1, x){
  return(exp(b0+b1*x)/(1+exp(b0+b1*x)))
}

one_or_zero_given_b0_and_b1_func = function(b0,b1,x){
  pi0 = probability_logistic_coeffs_func(b0,b1,x) 
  return(rbinom(1,1,pi0))
  # return(round(pi0))
}

hockey_stick_func = function(xvals, intercept, m, thresh, power){
  yvals = array(NA, length(xvals))
  for (j in 1:length(xvals)){
    if (xvals[j] < thresh){
      yvals[j] = intercept
    }else{
      yvals[j] = m*(xvals[j]-thresh)^power + intercept
    }
  }
  return(yvals)
}

matrix_hockey_stick_func = function(XMAT, intercept, m, thresh, power){
  YMAT = array(NA, dim(XMAT))
  for (j in 1:dim(YMAT)[2]){
    YMAT[,j] = hockey_stick_func(XMAT[,j], intercept, m, thresh, power)
  }
  return(YMAT)
}

rigid_s_func = function(xvals, intercept, m, thresh1, thresh2, power){
  yvals = array(NA, length(xvals))
  for (j in 1:length(xvals)){
    if (xvals[j] >= thresh2){
      yvals[j] = (hockey_stick_func(thresh2, intercept, m, thresh1, power))
    }else{
      yvals[j] = (hockey_stick_func(xvals[j], intercept, m, thresh1, power))
    }
  }
  return(yvals)
}

broken_hockey_stick_func = function(xvals, intercept, m, thresh, alpha){ # xvals, intercept, m, thresh, alpha
  yvals = array(NA, length(xvals))
  for (j in 1:length(xvals)){
    if (xvals[j] < thresh){
      yvals[j] = intercept
    }else{
      yvals[j] = m*(xvals[j]) + alpha
    }
  }
  return(yvals)
}

power_func = function(xvals, intercept, slope, power){
  yvals = array(NA, length(xvals))
  for (j in 1:length(xvals)){
    yvals[j] = intercept + slope*xvals[j]^power
  }
  return(yvals)
}

logit = function(x){return(log(x/(1-x)))} # creates a logit function (only works from values 0 < x < 1)
inverse_logit = function(x){return(1/(1+exp(-x)))} # inverses logit values
standardize = function(x){return((x-mean(x))/sd.real(x))}
add_or_subtract_a_little_from_zeros_and_ones = function(x){
  if (x == 0){
    return(x + .0001)
  }else if (x==1){
    return(x - .0001)
  }else{
    return(x)
  }
}

