# Dependencies: 18.01.21_paper1_calledon.R

rm(list=ls())

require(MASS)
require(fitdistrplus)
require(RColorBrewer)
require(leaps)
require(fmsb)
source('18.01.21_paper1_calledon.R')

# NAME: WILLIAM CUELLO
# PURPOSE: FOR PAPER IN COLLABORATION WITH SEBASTIAN SCHREIBER, ANDY SIH, JENNIFER GREMER, PETE TRIMMER.

observed_ess_values_gremer_2014 = c(.648424684, .542108916, .463880691, .422518228, .419136362, .71601903, .672763181, .3281005, .088460635, .176970447) # Observed ESS values taken from Jenny
jenny_baseline_ess_values = c(.7431519, .5518486, .5274235, .4028515, .535545, .487864, .6176489, .5204265, .2183197, .492692) # Jenny's predicted values

s_mobe = c(.593, .635, .102,.430) # survivorships for mobe, plin, plpa, etc.  
s_plin = c(.813, .836, .168,.713)
s_plpa = c(.788, .814, .160,.676)
s_pere = c(.836, .858, .173,.731)
s_stmi = c(.715, .743, .145,.616)
s_erte = c(.777, .799, .170,.722)
s_erci = c(.731, .766, .132,.560)
s_erla = c(.795, .823, .153,.646)
s_evmu = c(.907, .916, .214,.904)
s_scba = c(.778, .807, .150,.635)
s_mat = rbind(s_mobe, s_plin, s_plpa, s_pere, s_stmi, s_erte, s_erci, s_erla, s_evmu, s_scba) # placing them in a matrix
colnames(s_mat) = c('Fall', 'Winter', 'Summer_new', 'Summer_old')

Tf = 1*10^4 #6*10^4 # (Tf) precipitation values to be generated for ESS 
base = rnorm(Tf, mean = 0, sd=1) # Generate (Tf) N(0,1) distributed values
delta = 1e-2 # Just changed the extinction threshold to be a fraction of (to be calculated) the mean of an individual species' run
len_of_sdvec_bsvec = 5 # YOU MUST CHOOSE AN ODD (# > 1) length of vector for percentage differences in SD and MEAN
a_comp_vec = c(.0163, .0380, .0069, .0167, .0156, .0532, .0351, .0114, .0119, .0112) # Competitive Params found from Venable & Gremer
num_species = length(a_comp_vec)
s_vec = as.numeric(s_mat[,2]*s_mat[,4]) # DONT TOUCH -- Note that when only using s_mat[,4] we get eerily close values to Jenny's estimates when sampling from her yield values...
s_vec_new = as.numeric(s_mat[,3])
amt = Tf/2 # From the (Tf)-long ESS run, we grab the last (amt) of population densities and yields
t.mat = 100 # We push each of the (amt) of population densities forward with NEWLY generated Precipitation Values

new_string_of_names = c() # creating an array to start appending to in the next snippet
new_string_of_names[1] = "mobe"
new_string_of_names[2] = "plin" 
new_string_of_names[3] = "plpa"
new_string_of_names[4] = "pere"
new_string_of_names[5] = "stmi"
new_string_of_names[6] = "erte"
new_string_of_names[7] = "erci"
new_string_of_names[8] = "erla"
new_string_of_names[9] = "evmu"
new_string_of_names[10] = "scba"




# Loading in Venable Data
beginning_year = 1984
end_year = 2013
table_precipitation = read.csv('LTERDATA1982_2016.csv')
table_precipitation_jenny_2005_2016 = read.csv('daily_precip_clean_all.csv')
array_of_years = as.character(2006:end_year) # because Venable had 1984-2005
array_of_months =  c('-01-', '-02-', '-03-')
array_of_first_days = c('01', '02', '03', '04', '05', '06', '07', '08', '09')
array_of_all_days = c(array_of_first_days, as.character(10:27))
leap_years = as.character(c(1988, 1992, 1996, 2000, 2004, 2008, 2012, 2016, 2020))
running_precip_sum = 0
running_year_count = 1
vector_of_precipitations = array(NA, length(array_of_years))

for (prep_year in array_of_years){
  for (prep_month in array_of_months){
    for (prep_day in array_of_all_days){
      current_date = paste(prep_year, prep_month, sep='')
      current_date = paste(current_date, prep_day, sep='')
      table_row = which(table_precipitation$Date == current_date)
      running_precip_sum = running_precip_sum + is.numeric(table_precipitation[which(table_precipitation$Date == current_date),]$Precip..total.mm.)*table_precipitation[which(table_precipitation$Date == current_date),]$Precip..total.mm.
    }
  }
  vector_of_precipitations[running_year_count] = running_precip_sum
  running_year_count = running_year_count + 1
  running_precip_sum = 0
}

head(table_precipitation)

vector_of_precipitations = vector_of_precipitations*.1 # everything is off by a factor of 10 -- these are the precipitation values that *I* found from the LTREB sites (from 2006 to 2012)

up_to_this = end_year - 2005
vec_precips_2006_endyear_GBGLTP = table_precipitation_jenny_2005_2016$GBGLTP[1:up_to_this]*.1
vec_precips_2006_endyear_GBGPerm = table_precipitation_jenny_2005_2016$GBGPerm[1:up_to_this]*.1
vec_precips_2006_endyear_GBGPerm[7] = vector_of_precipitations[7] # since this was NOT AVAILABLE (this is why you can't pick something less than 6 years in advance of 2006)
vec_precips_2006_endyear_GBGLTP[7]  = vector_of_precipitations[7] # since this was NOT AVAILABLE
vec_precips_2006_endyear_GBGPerm[2] = vector_of_precipitations[2] # since this was NOT AVAILABLE (this is why you can't pick something less than 6 years in advance of 2006)
vec_precips_2006_endyear_GBGLTP[2]  = vector_of_precipitations[2] # since this was NOT AVAILABLE
vector_of_precipitations = (vec_precips_2006_endyear_GBGPerm + vec_precips_2006_endyear_GBGLTP)/2

table = read.csv('DataFig2Venable2007.csv')
starting_point_of_table = tail(which(table$year == beginning_year),1)
vector_of_precipitations = c(table[starting_point_of_table:dim(table)[1],]$Jmprec, vector_of_precipitations)
names(table)[6] = 'seedlings_per_m2'

table3 = read.csv('SpeciesByYear.csv')
table3 = table3[min(which(table3$year == c(beginning_year))):max(which(table3$year == c(end_year))),1:9] # table purely to create K-matrix and replicate Jenny's work
table3[155,]$lxbx = '.' # making edit so that it is consistent with Venable
names(table3)[9] = 'Ktvals'


############# END of EDITS


Kt_matrix3 = matrix(NA, nrow = num_species, ncol = (end_year-(beginning_year-1))) #  Low Seedling Density Yields
Kt_matrix3_dep = matrix(NA, nrow = num_species, ncol = (end_year-(beginning_year-1))) # Density Dependent Yields
list_indices_had_NA = list(c())
list_indices_had_NA_or_zero_lb = list(c())
add_point_this = 0.5

jenny_seedling_func = function(x){return(max(0, x-20))}
# jenny_seedling_func = function(x){return(max(0, x))} # if we don't want to incorporate this weird low-density thing.
table3$lxbx[which(table3$lxbx == '.')] = NA

for (cc in 1:num_species){
  Kt_matrix3[cc,] = as.numeric(as.character(table3[which(table3$species == new_string_of_names[cc]),]$lxbx)) + add_point_this
  Kt_matrix3_dep[cc,] = (1+a_comp_vec[cc]*sapply(as.numeric(as.character(table3[which(table3$species == new_string_of_names[cc]),]$seedlings_per_m2)), jenny_seedling_func))*(as.numeric(as.character(table3[which(table3$species == new_string_of_names[cc]),]$lxbx))+ add_point_this)
}

independent.ess.vec = array(NA, 10) # density independent (no competition term)
dependent.ess.vec = array(NA,10) # this is solving with the competition term
independent.stochastic.growth.rates.at.ess = array(NA,10) # density independent (no competition term)
dependent.stochastic.growth.rates.at.ess = array(NA,10) # growth rates with the competition term

# stopifnot(2<1)

modified_beginning = 1
modified_end = (end_year - beginning_year) # 30 should be all the years... 
vector_of_precipitations = vector_of_precipitations[modified_beginning:modified_end]
Kt_matrix3 = Kt_matrix3[,modified_beginning:modified_end]
Kt_matrix3_dep = Kt_matrix3_dep[,modified_beginning:modified_end]

take_out_these_indices_because_they_were_bad_years = c(11,23)

Kt_matrix3 = Kt_matrix3[,-take_out_these_indices_because_they_were_bad_years]
Kt_matrix3_dep = Kt_matrix3_dep[,-take_out_these_indices_because_they_were_bad_years]
vector_of_precipitations = vector_of_precipitations[-take_out_these_indices_because_they_were_bad_years] # can comment out 23... dry year? 

for (counter in 1:num_species){
  list_indices_had_NA[[counter]] = which(is.na(Kt_matrix3[counter,]))
  list_indices_had_NA_or_zero_lb[[counter]] = c(which(is.na(Kt_matrix3[counter,])), which((as.numeric(as.character(table3[which(table3$species == new_string_of_names[counter]),]$lxbx))) == 0))
}

#problem: taking out the 11th year (1994) and 23rd year (2000 whatever that Jenny mentioned) makes it so there aren't that many NA
#temp solution: manually just make them all really big values so that the code tries to subtract an index of big value -- i.e. does nothing to vector
for (index_for_nested_NA_tracker in 2:10){
  list_indices_had_NA[[index_for_nested_NA_tracker]] = 10000
}

for (ii in 1:(length(independent.ess.vec))){
  independent_ess_func = function(g){
    sample_Kts = sample(Kt_matrix3[ii,][-list_indices_had_NA[[ii]]], Tf, replace = TRUE)
    numerator_Kts = sample_Kts*s_vec_new[ii] - s_vec[ii]
    denominator_Kts = g*sample_Kts*s_vec_new[ii] + (1-g)*s_vec[ii]
    return(mean(numerator_Kts/denominator_Kts))
  }
  independent_growth_rate_func = function(g){
    sample_Kts = sample(Kt_matrix3[ii,][-list_indices_had_NA[[ii]]], Tf, replace = TRUE)
    expression = log(sample_Kts*g*s_vec_new[ii] + (1-g)*s_vec[ii])
    return(mean(expression))
  }
  independent.ess.vec[ii] =uniroot(independent_ess_func, c(0,1))$root
  independent.stochastic.growth.rates.at.ess[ii] = independent_growth_rate_func(independent.ess.vec[ii])
}

dependent_growth_rate_func = function(Ys, sval, g){ # this assumes the Ys already multiplied in s_new
  # sample_Kts_dep = sample(Kt_matrix3_dep[ii,][-list_indices_had_NA[[ii]]])
  expression_dep = log(Ys*g + (1-g)*sval)
  return(mean(expression_dep))
}


plot(independent.ess.vec, observed_ess_values_gremer_2014, ylim = c(0,1), xlim=c(0,1))
text(independent.ess.vec, observed_ess_values_gremer_2014, labels = new_string_of_names, pos = 4)
abline(c(0,0), c(1,1), lty =2)
remove_species_indices= -c(11) # I didn't want to change the rest of the code (11 is out of the species set and so removing the 11th index does nothing)

new_string_of_names_permanent = new_string_of_names # holding onto these for later use
new_string_of_names = new_string_of_names[remove_species_indices]
a_comp_vec = a_comp_vec[remove_species_indices] # make the non-working species have the same comp as before
s_vec = s_vec[remove_species_indices] # make the non-working species have the same SURVIVORSHIP as before

# TREATING THE COLORS -- I really need to automate this, cause this is ridiculous
the_colors = colorRampPalette(c('light blue', 'dark blue'))(10)
the_colors = the_colors[remove_species_indices] # make the non-working species have same COLORS as before


num_species = length(new_string_of_names) # Number of species found within the table
k = length(new_string_of_names) # Sometimes I use k; sometimes I use num_species
jen_alphas = array(NA,k)
jen_betas = array(NA,k)
alphas = array(NA, k)
betas = array(NA, num_species)# technically b, slope of log-log plot
hockey_betas = array(NA, num_species)
hockey_taus = array(NA, num_species)
logistic_hockey_b0s = array(NA, num_species)
logistic_hockey_b1s = array(NA, num_species)
conditioned_hockey_b2s= array(NA, num_species)
conditioned_hockey_taus = array(NA, num_species)
conditioned_jenny_alphas = array(NA, num_species)
conditioned_jenny_betas = array(NA, num_species)
aic_values_ls = array(NA, num_species)
r_squared_ls_vals = array(NA, num_species)
set_of_intercepts = array(NA, num_species)


# RECREATING VENABLE (2007) PLOTS 

for (array_index in 1:k){
  str = new_string_of_names[array_index]
  print(str)
  precip = vector_of_precipitations[-list_indices_had_NA[[array_index]]] #as.numeric(noquote(extended_table[index,2]))[-list_indices_had_NA[[array_index]]]
  log_precip_plus1 = log(precip + 1) # As Venable did, take the precipitation and add 1
  log_lxbx_plus_point5 = log(Kt_matrix3[array_index,][-list_indices_had_NA[[array_index]]])
  ven_seedling_densities = sapply(as.numeric(as.character(table3[which(table3$species == str),]$seedlings_per_m2)), jenny_seedling_func)[-list_indices_had_NA[[array_index]]] # read in the venable seedling densities
  log_k_vec = log(Kt_matrix3_dep[array_index,][-list_indices_had_NA[[array_index]]]) # this includes the competition
  
  separating_val = 0
  replace_or_not_vec = zeros_and_ones_func(log_k_vec, y_hyperplane = separating_val) # or put log_k_vec
  logistic_hockey = glm(replace_or_not_vec ~ log_precip_plus1, family = 'binomial')
  which_log_ks_to_regress = which(replace_or_not_vec == 1)
  
  plot(log_precip_plus1, replace_or_not_vec, ylab = 'Logistic Regression', xlab = 'log(Precip+1)', main = new_string_of_names[array_index], xlim = c(-3,3))
  domain_for_logistic = seq(-3,3, by = .1)#seq(min(log_precip_plus1), max(log_precip_plus1), by = .1)
  temp_logistic_func = function(beta0, beta1, x){return(1/(1+exp(-(beta0 + beta1*x))))}
  lines(domain_for_logistic, temp_logistic_func(coef(logistic_hockey)[[1]],coef(logistic_hockey)[[2]],domain_for_logistic), xlim = c(-3,3))
  
  set_of_intercepts[array_index] = log(.5)
  set_intercept = set_of_intercepts[array_index]
  conditioned_hockey = nls(log_k_vec[which_log_ks_to_regress] ~ hockey_stick_func(log_precip_plus1[which_log_ks_to_regress], intercept = set_intercept, m = beta, thresh = set_thresh, power = 1), start = c(beta = 1, set_thresh = 1))
  conditioned_ls = lm(log_k_vec[which_log_ks_to_regress] ~ log_precip_plus1[which_log_ks_to_regress])
  
  plot(log_precip_plus1[which_log_ks_to_regress], log_k_vec[which_log_ks_to_regress], ylim = c(-1, 8), xlim = c(-.5, 5), xlab = 'log(Precip+1)', ylab = 'log(Yield + .5)', main = paste(new_string_of_names[array_index]), cex.lab = 1.5, cex.axis = 1.5)
  domain = seq(0, 5, by = .01)
  abline(conditioned_ls, lty = 2, lwd =2)
  lines(domain, hockey_stick_func(domain, set_intercept, coef(conditioned_hockey)[1], coef(conditioned_hockey)[2], power = 1))
  abline(v = 0, h = 0, lty = 4)
  
  print(new_string_of_names[array_index])
  print(summary(conditioned_ls))
  print(summary(logistic_hockey))
  
  logistic_hockey_b0s[array_index] = coef(logistic_hockey)[[1]]
  logistic_hockey_b1s[array_index] = coef(logistic_hockey)[[2]]
  conditioned_hockey_b2s[array_index] = coef(conditioned_hockey)[[1]] 
  conditioned_hockey_taus[array_index] = coef(conditioned_hockey)[[2]]
  conditioned_jenny_alphas[array_index] = coef(conditioned_ls)[[1]] # Note that his equivalent to (-ch_b2s*ch_taus + log(.5)): look at hockey_stick formula to see why 
  conditioned_jenny_betas[array_index] = coef(conditioned_ls)[[2]]
  
  if (array_index >= (k-1)){
    plot(log_precip_plus1, log_lxbx_plus_point5, main = paste(str), ylim = c(-3,9), xlim = c(-1, 4), ylab = 'ln(lxbx + .5)', xlab = 'ln(precip + 1)', col = the_colors[array_index])
    abline(v = 0, h = 0, lty = 2)
  } else {
    plot(log_precip_plus1, log_lxbx_plus_point5, main = paste(str), ylim = c(-3,9), xlim = c(-1, 4), ylab = 'ln(lxbx+.5)', xlab = '', col = the_colors[array_index]) # Just to check if we got the same as Venable
    abline(v = 0, h = 0, lty = 2)
  }
}
dev.off()

# CREATING TABLE of all params in order of survival rates
sorted_standard_thing = sort(s_vec)
sorted_surv = sort(s_vec)
experiment_colors = colorRampPalette(c('light blue', 'dark blue'))(10)
standard_thing_and_colors = rbind(sorted_standard_thing, experiment_colors)
original_standard_thing_indices = match(sorted_surv, s_vec)

table_of_fitting_params = cbind(s_vec, s_vec_new, a_comp_vec, conditioned_jenny_alphas, conditioned_jenny_betas, conditioned_hockey_taus, logistic_hockey_b0s, logistic_hockey_b1s, jenny_baseline_ess_values, observed_ess_values_gremer_2014, set_of_intercepts)
colnames(table_of_fitting_params) = c('survs', 'survs_new', 'a_comp', 'alphas', 'betas', 'taus', 'LH_b0s', 'LH_b1s','jenESS', 'fieldESS', 'hockey_value')
rownames(table_of_fitting_params) = new_string_of_names
table_of_fitting_params_surv_sorted = table_of_fitting_params[original_standard_thing_indices,]
table_of_fitting_params_surv_sorted = as.data.frame(table_of_fitting_params_surv_sorted)
use_table = table_of_fitting_params_surv_sorted


lognormal_model_for_precip = fitdist(log(precip + 1), 'norm')#fitdist(log(precip+1), 'norm')#[-which(is.na(log(precip+1)))], 'norm')
m1 = as.numeric(lognormal_model_for_precip$estimate[1]) #mean(log_precip_plus1) # MEAN OF OUR DATA (i.e. getting the precipitation's log-mean to generate our own)
sd1 = as.numeric(lognormal_model_for_precip$estimate[2]) #sd.real(log_precip_plus1) # SD OF OUR DATA

# Creating log-log plot

matrix_of_bases = matrix(rnorm(t.mat*amt), nrow = t.mat, ncol = amt) # Non-autocorrelated normals to be transformed into autocorrelated later
sd1 = sd1
m1 = m1
P = exp(base*sd1 + m1)

Ps=seq(min(P),max(P),length=100) # For preparation of the _logprecip.pdf
Ys=matrix(0,100,k) # For preparation of the _logprecip.pdf
for(i in 1:k){
  Ys[,i]= exp(use_table$alphas[i])*Ps^(use_table$betas[i])
} # creating the Yields to be put in the log-log

rho = 0

##### CALCULATING THE BASELINE ESS VALUES

ess.vec.baseline = array(1,k) # Evolutionary Stable Strategy (empty)
P_auto = auto_cor_func(P,rho) # autocorrelate these values (if rho = 0, no autocorrelation)
##### DONE CALCULATING THE BASELINE ESS VALUES

print(paste('Creating pdfs for rho', rho))

file.name=paste('Desert_Annual_MAY7_testing','Time =', Tf,"delta =", delta, "rho =", rho) 

###############################################################
# Creating a new P (like at the beginning) -- more variance ###
###############################################################

# We need to generate new P values with increased standard deviation,
# run dynamics long enough, so that we have autocorrelated Yield values that correspond to new densities,
# and then we have to use the OLD p-matrix to keep the autocorrelation going, so that
# we have some consistency. 
sd.change.vec = seq(-1, 1, length = len_of_sdvec_bsvec) 
ext.rate.matrix.postESS.SD = array(NA, c(length(sd.change.vec), k))
ext.rate.matrix.preESS.SD = array(NA, c(length(sd.change.vec), k))
ess.vals.matrix = array(NA, c(length(sd.change.vec), k))
nhat.matrix = array(NA, c(length(sd.change.vec),k))
nsd.matrix = array(NA, c(length(sd.change.vec),k))
i.row = 1

Nhat.vec = numeric(k) # log mean density vector (empty)
Nsd.vec = numeric(k) # SD log density vector (empty)
ess.vec.newSD = array(1,k) # To hold the NEW Evolutionary Stable Strategy (empty)
ess.error.vec.newSD = array(NA, k) # Error values for candidate ESS (empty)
went.below.threshold.preESS = array(0, amt) # array of 0's to replaced with 1 whenever extinction occurs
went.below.threshold.postESS = array(0, amt)
extinction_rate_preESS = array(NA,k) # probalities of extinction
extinction_rate_postESS = array(NA,k)
threshold_factor_post_MEAN_change_mean = array(NA,k)
threshold_factor_post_SD_change_mean = array(NA,k)
threshold_factor_pre_mean = array(NA,k)
threshold_for_SD = array(NA,k)
Ydata_baseline_holder = array(NA, c(length(P_auto), k))
our_yield_sds = array(NA, k)
jenny_yield_sds = array(NA,k)
our_yield_means = array(NA,k)
jenny_yield_means = array(NA,k)

for (ii in 1:k){
  temp_func = function(x){return(one_or_zero_given_b0_and_b1_func(use_table$LH_b0s[ii], use_table$LH_b1s[ii], x))}
  Ydata_baseline = sapply(log(P_auto), temp_func)*(exp(use_table$alphas[ii])*P_auto^use_table$betas[ii])
  Ydata_baseline[which(Ydata_baseline == 0)] = 0.5
  Ydata_baseline = Ydata_baseline*(use_table$survs_new[ii])
  our_yield_sds[ii] = sd.real(Ydata_baseline)
  our_yield_means[ii] = mean(Ydata_baseline)
  ## JENNY'S YIELD -- just for double-checking things
  Ydata_baseline_jenny = sample(Kt_matrix3_dep[ii,][-list_indices_had_NA[[ii]]], Tf, replace=TRUE)*use_table$survs_new[ii]
  jenny_yield_sds[ii] = sd.real(Ydata_baseline_jenny)
  jenny_yield_means[ii] = mean(Ydata_baseline_jenny)
  # print(paste('Jenny species', ii, 'yield mean is', mean(Ydata_baseline)))
  # ess.vec.baseline[ii] = ESS_func(Ydata_baseline, N0 = 1, s_vec[ii], a_comp_vec[ii])
  # dependent.stochastic.growth.rates.at.ess[ii] = dependent_growth_rate_func(Ydata_baseline, s_vec[ii], ess.vec.baseline[ii])
  
  Ydata_baseline_holder[,ii] = Ydata_baseline
  ess.vec.baseline[ii] = ESS_func(Ydata_baseline, N0 = 1, use_table$survs[ii], use_table$a_comp[ii]) # BASELINE YIELD --> BASELINE ESS VALUES
  dependent.stochastic.growth.rates.at.ess[ii] = dependent_growth_rate_func(Ydata_baseline, use_table$survs[ii], ess.vec.baseline[ii])
  out_preESS_with_old_Ydata = Ns_func(Ydata_baseline, N0 = 1, ess.vec.baseline[ii], use_table$survs[ii], use_table$a_comp[ii])#nl_a_comp_vec[jj]) # BASE YIELD + BASE ESS --> BASELINE ADAPTED DENSITIES
  densities_preESS_with_old_Ydata = out_preESS_with_old_Ydata$pops # BASELINE DENSITIES
  threshold_for_SD[ii] = as.numeric(quantile(densities_preESS_with_old_Ydata[(Tf/2):Tf], 1/length(densities_preESS_with_old_Ydata[(Tf/2):Tf])))*delta
  print(paste('ess values', ess.vec.baseline))
  print(paste('threshold values', threshold_for_SD))
}


# SUMMARY: THIS GENERATES THE REGRESSION PLOTS BUT FOR ALL THE TRUE PARAMETERS FROM THE DATA SET

######################################
# PREDICTORS AND RESPONSE VARIABLES###
######################################

sold = use_table$survs
snew = use_table$survs_new
alpha1s = use_table$LH_b0s
beta1s = use_table$LH_b1s
CH_taus = use_table$taus
beta2s = use_table$betas
alpha2s = use_table$alphas 
field_ess_vals = use_table$fieldESS
pred_ess_vals = ess.vec.baseline
transfield_ess_vals = logit(field_ess_vals)
transpred_ess_vals = logit(pred_ess_vals)
stan_beta2s = standardize(beta2s)
stan_beta1s = standardize(beta1s)
stan_snew = standardize(snew)
stan_sold = standardize(sold)

##########################
# END OF THE PARAMETERS # 
#########################

data_table_with_params = cbind(sold, snew, CH_taus, beta2s, alpha1s, beta1s, alpha2s)
data_table_with_params = as.data.frame(data_table_with_params)

###########################
# Predicting ESS Values ## 
#########################

observedESS_which_model_to_choose = regsubsets(transfield_ess_vals ~ beta1s + beta2s + sold, data = data_table_with_params, nbest = 4, nvmax = 14, method = c("exhaustive"))
predictedESS_which_model_to_choose = regsubsets(transpred_ess_vals ~ beta1s + beta2s + sold, data = data_table_with_params, nbest = 4, nvmax = 14, method = c("exhaustive"))
plot(predictedESS_which_model_to_choose)
plot(observedESS_which_model_to_choose)

summary(lm(transpred_ess_vals ~ stan_beta2s + stan_sold + stan_beta1s))
summary(lm(transfield_ess_vals ~  stan_beta2s +  stan_sold + stan_beta1s))
summary(lm(use_table$fieldESS ~ ess.vec.baseline))
summary(lm(transfield_ess_vals ~ ess.vec.baseline)) # interesting -- this does pretty well


### In prep for plotting linear regression (transformed back) results against observed
stats_model = lm(transfield_ess_vals ~ beta2s + sold + beta1s)
predicting_logit_ESS = coef(stats_model)[[1]] + coef(stats_model)[[2]]*beta2s + coef(stats_model)[[3]]*sold + coef(stats_model)[[4]]*beta1s
predicting_ess = inverse_logit(predicting_logit_ESS)

#####################################################################################################
### Taken from earlier for loop: create figure for species of interest in file 17.10.11_prez.R    ###
#####################################################################################################

presentation_index = 8
str = new_string_of_names[presentation_index]
precip = vector_of_precipitations[-list_indices_had_NA[[presentation_index]]] #as.numeric(noquote(extended_table[index,2]))[-list_indices_had_NA[[array_index]]]
log_precip_plus1 = log(precip + 1) # As Venable did, take the precipitation and add 1
log_lxbx_plus_point5 = log(Kt_matrix3[presentation_index,][-list_indices_had_NA[[presentation_index]]])
ven_seedling_densities = sapply(as.numeric(as.character(table3[which(table3$species == str),]$seedlings_per_m2)), jenny_seedling_func)[-list_indices_had_NA[[presentation_index]]] # read in the venable seedling densities
log_k_vec = log(Kt_matrix3_dep[presentation_index,][-list_indices_had_NA[[presentation_index]]]) # this includes the competition
separating_val = 0
replace_or_not_vec = zeros_and_ones_func(log_k_vec, y_hyperplane = separating_val) # or put log_k_vec
logistic_hockey = glm(replace_or_not_vec ~ log_precip_plus1, family = 'binomial')
which_log_ks_to_regress = which(replace_or_not_vec == 1)
set_of_intercepts[array_index] = log(.5)
set_intercept = set_of_intercepts[presentation_index]
conditioned_hockey = nls(log_k_vec[which_log_ks_to_regress] ~ hockey_stick_func(log_precip_plus1[which_log_ks_to_regress], intercept = set_intercept, m = beta, thresh = set_thresh, power = 1), start = c(beta = 1, set_thresh = 1))
conditioned_ls = lm(log_k_vec[which_log_ks_to_regress] ~ log_precip_plus1[which_log_ks_to_regress])