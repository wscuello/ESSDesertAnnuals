######################################################################################
## NAME: WILLIAM CUELLO                                                             ##
## COAUTHORS: DRS. SEBASTIAN SCHREIBER, ANDY SIH, JENNIFER GREMER, PETE TRIMMER.    ##
## JOURNAL: PROC B                                                                  ##
######################################################################################

rm(list=ls())
require(MASS)
require(fitdistrplus)
require(RColorBrewer)
require(leaps)
require(fmsb)
require(lme4)
require(dunn.test)
require(lsmeans)
source('18.01.21_paper1_calledon.R')

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

Tf = 2*10^4 #6*10^4 # (Tf) precipitation values to be generated for ESS 
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

############################
# Loading in VENABLE DATA ##
############################
beginning_year = 1984
end_year = 2013

# this following table is me exploring other years from other weather stations -- not the one used by Venable/Gremer
table_precipitation_for_1922_to_2017 = read.csv('precip1922_2017.csv') # taken from https://cals.arizona.edu/SRER/data.html -- measured in hundreth of inch
desert_station_precip_1922_to_2017 = subset(table_precipitation_for_1922_to_2017, STATION == 'DESST')[,1:5]
desert_station_precip_1922_to_2017[which(desert_station_precip_1922_to_2017 == -9999,  arr.ind=TRUE)] = NA
desert_station_precip_1922_to_2017$JMprecipMM = (desert_station_precip_1922_to_2017$JAN + desert_station_precip_1922_to_2017$FEB + desert_station_precip_1922_to_2017$MAR)*.254

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

table = read.csv('DataFig2Venable2007.csv') # contains a significant portion of the precipitation data -- sent by Venable
starting_point_of_table = tail(which(table$year == beginning_year),1)
vector_of_precipitations = c(table[starting_point_of_table:dim(table)[1],]$Jmprec, vector_of_precipitations)
names(table)[6] = 'seedlings_per_m2'

table3 = read.csv('SpeciesByYear.csv')
table3 = table3[min(which(table3$year == c(beginning_year))):max(which(table3$year == c(end_year))),1:9] # table purely to create K-matrix and replicate Jenny's work
table3[155,]$lxbx = '.' # making edit so that it is consistent with Venable
names(table3)[9] = 'Ktvals'

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

modified_beginning = 1
modified_end = (end_year - beginning_year + 1) # 30 should be all the years... 
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

# - Here's something that I just didn't clean up -- it works, but it's a cruddy inefficient way of doing things:
# - took out the 11th year (1994) and 23rd year (2006) makes it so there aren't NAs, which we originally coded.
# - solution: manually just make them all really big values so that the code tries to subtract an index of big value -- i.e. does nothing to vector; otherwise R script complains that there's no value to take from

for (index_for_nested_NA_tracker in 2:10){
list_indices_had_NA[[index_for_nested_NA_tracker]] = 10000
}

# for (ii in 1:(length(independent.ess.vec))){
#   independent_ess_func = function(g){
#     sample_Kts = sample(Kt_matrix3[ii,][-list_indices_had_NA[[ii]]], Tf, replace = TRUE)
#     numerator_Kts = sample_Kts*s_vec_new[ii] - s_vec[ii]
#     denominator_Kts = g*sample_Kts*s_vec_new[ii] + (1-g)*s_vec[ii]
#     return(mean(numerator_Kts/denominator_Kts))
#   }
#   independent_growth_rate_func = function(g){
#     sample_Kts = sample(Kt_matrix3[ii,][-list_indices_had_NA[[ii]]], Tf, replace = TRUE)
#     expression = log(sample_Kts*g*s_vec_new[ii] + (1-g)*s_vec[ii])
#     return(mean(expression))
#   }
#   independent.ess.vec[ii] =uniroot(independent_ess_func, c(0,1))$root
#   independent.stochastic.growth.rates.at.ess[ii] = independent_growth_rate_func(independent.ess.vec[ii])
# }

dependent_growth_rate_func = function(Ys, sval, g){ # this assumes the Ys already multiplied in s_new
  expression_dep = log(Ys*g + (1-g)*sval)
  return(mean(expression_dep))
}


# plot(independent.ess.vec, observed_ess_values_gremer_2014, ylim = c(0,1), xlim=c(0,1))
# text(independent.ess.vec, observed_ess_values_gremer_2014, labels = new_string_of_names, pos = 4)
# abline(c(0,0), c(1,1), lty =2)
remove_species_indices= -c(11) # I didn't want to change the rest of the code (11 is out of the species set and so removing the 11th index does nothing)
new_string_of_names_permanent = new_string_of_names # holding onto these for later use
new_string_of_names = new_string_of_names[remove_species_indices]
a_comp_vec = a_comp_vec[remove_species_indices] # make the non-working species have the same comp as before
s_vec = s_vec[remove_species_indices] # make the non-working species have the same SURVIVORSHIP as before

########################
## TREATING THE COLORS##
########################
the_colors = colorRampPalette(c('light blue', 'dark blue'))(10)
the_colors = the_colors[remove_species_indices] # make the non-working species have same COLORS as before

###########################################
## CREATING VECTORS TO HOLD COEFFICIENTS ##
###########################################

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

############################################################################
### RECREATING VENABLE (2007) PLOTS AND GRABBING LINEAR MODELS SEPARATELY###
############################################################################

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

#########################################
# CREATING TABLE HOLDING ALL PARAMETERS # 
#########################################
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

#####################################################################
## CREATING MATRICES FOR FULL BINOMIAL MODEL AND POST-HURDLE MODEL ##
#####################################################################
stats_precip_yield_matrix = matrix(NA, nrow = 28*10, ncol = 3) # 28 precipitation data points;  as.data.frame(Kt_matrix3_dep)
stats_precip_yield_matrix = as.data.frame(stats_precip_yield_matrix)
colnames(stats_precip_yield_matrix) = c('species', 'yield', 'precipitation')
Kt_matrix3_dep_NAs_for_point5s = Kt_matrix3_dep 
Kt_matrix3_dep_NAs_for_point5s[which(Kt_matrix3_dep_NAs_for_point5s <= 1, arr.ind = TRUE)] = NA
stats_precip_yield_matrix_taking_Kt_directly = as.data.frame(matrix(NA, nrow = 28*10, ncol = 3))
colnames(stats_precip_yield_matrix_taking_Kt_directly) = c('species', 'yield', 'precipitation')

for (j in 1:10){
  stats_precip_yield_matrix[(1+(j-1)*28):(j*28),1] = new_string_of_names[j]
  stats_precip_yield_matrix[(1+(j-1)*28):(j*28),2] = Kt_matrix3_dep_NAs_for_point5s[j,]
  stats_precip_yield_matrix[(1+(j-1)*28):(j*28),3] = vector_of_precipitations
  stats_precip_yield_matrix_taking_Kt_directly[(1+(j-1)*28):(j*28),1] = new_string_of_names[j]
  stats_precip_yield_matrix_taking_Kt_directly[(1+(j-1)*28):(j*28),2] = Kt_matrix3_dep[j,]
  stats_precip_yield_matrix_taking_Kt_directly[(1+(j-1)*28):(j*28),3] = vector_of_precipitations
}

stats_precipplusone_yield_matrix = stats_precip_yield_matrix
stats_precipplusone_yield_matrix$precipitation = stats_precipplusone_yield_matrix$precipitation + 1
stats_precip_yieldornoyield_matrix = stats_precip_yield_matrix_taking_Kt_directly
stats_precipplusone_yieldornoyield_matrix = stats_precip_yieldornoyield_matrix
stats_precipplusone_yieldornoyield_matrix$precipitation = stats_precip_yield_matrix_taking_Kt_directly$precipitation + 1

for (j in 1:dim(stats_precipplusone_yieldornoyield_matrix)[1]){ # MODIFYING DATA FOR THE BINOMIAL PART
  if (is.na(stats_precipplusone_yieldornoyield_matrix$yield[j])){stats_precipplusone_yieldornoyield_matrix$yield[j] = stats_precipplusone_yieldornoyield_matrix$yield[j]
  }else if (stats_precipplusone_yieldornoyield_matrix$yield[j] <= 1){
    stats_precipplusone_yieldornoyield_matrix$yield[j] = 0
  }else{
    stats_precipplusone_yieldornoyield_matrix$yield[j] = 1
  }
}

#################################
# CREATING REPRODUCTIVE HURDLES #
#################################
stats_precipplusone_yieldornoyield_matrix$species = as.factor(stats_precipplusone_yieldornoyield_matrix$species)

binomial_model_full = glm(yield ~ log(precipitation)*species, family = binomial(link = 'logit'), data = subset(stats_precipplusone_yieldornoyield_matrix, yield != 'NA' ))
binomial_model_sub1 = glm(yield ~ log(precipitation) + species, family = binomial(link = 'logit'), data = subset(stats_precipplusone_yieldornoyield_matrix, yield != 'NA'))
binomial_model_sub2 = glm(yield ~ log(precipitation), family = binomial(link = 'logit'), data = subset(stats_precipplusone_yieldornoyield_matrix, yield != 'NA'))
binomial_model_sub3 = glm(yield ~ log(precipitation) + species:log(precipitation), family = binomial(link = 'logit'), data = subset(stats_precipplusone_yieldornoyield_matrix, yield !='NA'))
binomial_model_sub4 = glm(yield ~ 0 + species + species:log(precipitation), family = binomial(link = 'logit'), data = subset(stats_precipplusone_yieldornoyield_matrix, yield !='NA')) # this is the full model; significance means the coefficient is non-zero instead of significantly different from common slope. 
summary(binomial_model_sub4)

stats_precipplusone_yieldornoyield_matrix$species = relevel(stats_precipplusone_yieldornoyield_matrix$species, ref = 'stmi')
binomial_model_sub1 = glm(yield ~ log(precipitation) + species, family = binomial(link = 'logit'), data = subset(stats_precipplusone_yieldornoyield_matrix, yield != 'NA'))
summary(binomial_model_sub1)

anova(binomial_model_full, binomial_model_sub1, test = 'Chisq')
anova(binomial_model_sub3, binomial_model_full, test = 'Chisq')
anova(binomial_model_sub1, binomial_model_sub2, test = 'Chisq')

# summary(glht(binomial_model_sub1, mcp(species = "Tukey")))
lsm_binomial = lsmeans(binomial_model_sub1, pairwise ~ species, var = 'log(precipitation)')
pairs(lsm_binomial)
summary(lsm_binomial, type = 'response')
# contrast(lsm_binomial, "trt.vs.ctrl", ref = c(1))
# confint(binomial_model_full)
# dwplot(binomial_model_sub1)

logLik(binomial_model_full)
logLik(binomial_model_sub1)
logLik(binomial_model_sub3)
AIC(binomial_model_full, binomial_model_sub1, binomial_model_sub3)
BIC(binomial_model_full, binomial_model_sub1, binomial_model_sub3)

##################################################################
# GRABBING THE COEFFICIENTS FROM THE BINOMIAL MODELS WE CREATED###
##################################################################
coefficients_of_binomial_model_sub1 = coef(binomial_model_sub1)
coefficients_of_binomial_model_sub3 = coef(binomial_model_sub3)
coefficients_of_full_binomial_model = coef(binomial_model_full)
use_table$LH_b0s_multi_intercepts_same_slope = NA
use_table$LH_b1s_multi_intercepts_same_slope = as.numeric(coefficients_of_binomial_model_sub1['log(precipitation)'])
use_table$LH_b0s_precip_plus_interaction = as.numeric(coefficients_of_binomial_model_sub3['(Intercept)']) # create new column to contain the new binomial model
use_table$LH_b1s_precip_plus_interaction = NA
use_table$LH_b0s_full_binomial_model = NA
use_table$LH_b1s_full_binomial_model = NA
for (j in 1:10){
  use_table$LH_b0s_multi_intercepts_same_slope[j] = coefficients_of_binomial_model_sub1['(Intercept)'] + coefficients_of_binomial_model_sub1[paste('species', row.names(use_table)[j], sep = '')]
  use_table$LH_b1s_precip_plus_interaction[j] = coefficients_of_binomial_model_sub3['log(precipitation)'] + coefficients_of_binomial_model_sub3[paste('log(precipitation):species', row.names(use_table)[j], sep = '')]
  use_table$LH_b0s_full_binomial_model[j] = coefficients_of_full_binomial_model[['(Intercept)']] + coefficients_of_full_binomial_model[paste('species',row.names(use_table)[j], sep='')]
  use_table$LH_b1s_full_binomial_model[j] = coefficients_of_full_binomial_model[['log(precipitation)']] + coefficients_of_full_binomial_model[paste('log(precipitation):species', row.names(use_table)[j], sep = '')]
}
use_table$LH_b0s_full_binomial_model[2] = coefficients_of_full_binomial_model[['(Intercept)']]
use_table$LH_b1s_full_binomial_model[2] = coefficients_of_full_binomial_model[['log(precipitation)']]
use_table$LH_b0s_multi_intercepts_same_slope[2] = as.numeric(coefficients_of_binomial_model_sub1['(Intercept)'])
use_table$LH_b1s_precip_plus_interaction[2] = as.numeric(coefficients_of_binomial_model_sub3['log(precipitation)']) # ERCI was the reference point for the model, so we have to input this manually

#######################################
##CREATING ALL THE POST-HURDLE MODELS##
#######################################

glm_model_part_two_full = glm(log(yield) ~ species+log(precipitation)+species:log(precipitation), family = gaussian(link = 'identity'), data = subset(stats_precipplusone_yield_matrix, yield != 'NA'))
glm_model_part_two_sub1 = glm(log(yield) ~ log(precipitation) + species:log(precipitation), family = gaussian(link = 'identity'), data = subset(stats_precipplusone_yield_matrix, yield != 'NA'))
glm_model_part_two_sub2 = glm(log(yield) ~ log(precipitation) + species, family = gaussian('identity'), data = stats_precipplusone_yield_matrix)
glm_model_part_report = glm(log(yield) ~ 0 + species + log(precipitation):species, family = gaussian('identity'), data = stats_precipplusone_yield_matrix) 
# ANALYSIS OF POST-HURDLE MODELS:
stats_precipplusone_yield_matrix$species = as.factor(stats_precipplusone_yield_matrix$species)
stats_precipplusone_yield_matrix$species = relevel(stats_precipplusone_yield_matrix$species, ref = 'stmi')
glm_model_part_two_sub1 = glm(log(yield) ~ log(precipitation) + species:log(precipitation), family = gaussian(link = 'identity'), data = subset(stats_precipplusone_yield_matrix, yield != 'NA'))
summary(glm_model_part_two_sub1)

anova(glm_model_part_two_full, glm_model_part_two_sub1, test = 'Chisq')
anova(glm_model_part_two_sub2, glm_model_part_two_full, test = 'Chisq')
BIC(glm_model_part_two_full, glm_model_part_two_sub1, glm_model_part_two_sub2) #full model, multiple slopes -- common intercept, common slopes -- multi intercept
AIC(glm_model_part_two_full, glm_model_part_two_sub1, glm_model_part_two_sub2) #full model, multiple slopes -- common intercept, common slopes -- multi intercept

######################################################################################################################
#Note: We list the species that had significantly different slopes. This is seen by calling the next function. ######
#erci-evmu (p < .02), erci-scba (p<.0006), erte-evmu (p < .09), erte-scba (p < .0038), plpa-scba (p < .0851)   ######
######################################################################################################################
require(lsmeans)
posthurdle_model_slopes_significantly_different = lstrends(glm_model_part_two_sub1, pairwise ~ species, var = "log(precipitation)", adjust = 'Bonferroni') 
# EVMU: ERCI
# ERTE: SCBA
# ERCI: SCBA

###########################################
# ATTEMPT AT  LINEAR MIXED-EFFECT MODELS ##
###########################################

# # linear mixed effect models -- opted not to do this: 
# mixed_model_full = lmer(log(yield) ~ log(precipitation) + (1 + log(precipitation) | species), data = subset(stats_precipplusone_yield_matrix, species != 'erci'))
# mixed_model_random_intercept = lmer(log(yield) ~ log(precipitation) + (1|species), data = subset(stats_precipplusone_yield_matrix, species != 'erci'))
# anova(mixed_model_full, mixed_model_random_intercept)


# # SOME ANALYSIS ON THE MIXED EFFECT FOR POST-HURDLE
# random_effect_table = as.data.frame(ranef(mixed_model_full))
# sd.real(subset(random_effect_table, term == '(Intercept)')$condval)
# sd.real(subset(random_effect_table, term == 'log(precipitation)')$condval)
# cov(subset(random_effect_table, term == '(Intercept)')$condval, subset(random_effect_table, term == 'log(precipitation)')$condval)

#####################################################################
# Post-HOC Tests for DIFFERENCES in YIELD values from each species###
#####################################################################

# kruskal_info_for_posthurdle = kruskal.test(subset(stats_precipplusone_yield_matrix, yield!= 'NA')$yield, as.factor(subset(stats_precipplusone_yield_matrix, yield!='NA')$species)) # kruskal-wallis test; note we have to make the second argument as factors; this is on non-zero yield alone -- not dependent on precipitation; this shows a significant difference in yield values
# dunn_info_for_posthurdle = dunn.test(subset(stats_precipplusone_yield_matrix, yield!= 'NA')$yield, as.factor(subset(stats_precipplusone_yield_matrix, yield!='NA')$species)) # Dunn Test to see what's not from the same distribution
# wilcox_info = pairwise.wilcox.test(subset(stats_precipplusone_yield_matrix, yield!= 'NA')$yield, as.factor(subset(stats_precipplusone_yield_matrix, yield!='NA')$species)) #
# kruskal_info_for_posthurdle_slopes = kruskal.test(subset(stats_precipplusone_yield_matrix, yield!= 'NA')$individual_slopes_on_log_scale, as.factor(subset(stats_precipplusone_yield_matrix, yield!='NA')$species)) # kruskal-wallis test; note we have to make the second argument as factors; this is on non-zero yield alone -- not dependent on precipitation; this shows a significant difference in yield values
# dunn_info_for_posthurdle_slopes = dunn.test(subset(stats_precipplusone_yield_matrix, yield!= 'NA')$individual_slopes_on_log_scale, as.factor(subset(stats_precipplusone_yield_matrix, yield!='NA')$species)) # Dunn Test to see what's not from the same distribution


# ######################################################
# # TEXTBOOK ANALYSIS FOR SEEING DIFFERENCES IN SLOPES##
# ######################################################
# stats_precipplusone_yield_matrix$individual_slopes_on_log_scale = NA
# stats_precipplusone_yield_matrix$individual_squared_difference_from_mean_slope_on_log_scale = NA
# for (row in 1:(dim(stats_precipplusone_yield_matrix)[1])){ # This populates the individual slopes. 
#   stats_precipplusone_yield_matrix[row,]$individual_slopes_on_log_scale = (log(stats_precipplusone_yield_matrix[row,]$yield) - coef(glm_model_part_two_sub1)[[1]])/log(stats_precipplusone_yield_matrix[row,]$precipitation)
# }
# 
# vector_of_mean_slopes_on_log_scale = c() # to be populated later. 
# 
# for (row in 1:(dim(stats_precipplusone_yield_matrix)[1])){
#   temporary_species = stats_precipplusone_yield_matrix[row,]$species
#   temporary_subset = subset(stats_precipplusone_yield_matrix, species == temporary_species)$individual_slopes_on_log_scale
#   temporary_mean_slope_on_log_scale = mean(temporary_subset[-which(is.na(temporary_subset))])
#   vector_of_mean_slopes_on_log_scale = union(vector_of_mean_slopes_on_log_scale, temporary_mean_slope_on_log_scale)
#   stats_precipplusone_yield_matrix[row,]$individual_squared_difference_from_mean_slope_on_log_scale = (stats_precipplusone_yield_matrix[row,]$individual_slopes_on_log_scale - temporary_mean_slope_on_log_scale)^2
# }
# 
# SSE = sum(stats_precipplusone_yield_matrix$individual_squared_difference_from_mean_slope_on_log_scale[which(!is.na(stats_precipplusone_yield_matrix$individual_squared_difference_from_mean_slope_on_log_scale))])
# degrees_of_freedom_for_error = (sum(!is.na(stats_precipplusone_yield_matrix$individual_squared_difference_from_mean_slope_on_log_scale))- num_species)
# MSE = SSE/degrees_of_freedom_for_error
# qalpha = 4.55 # from table of studentized range critical values (Applied Statistics for Engineers and Scientists)
# vector_of_number_of_samples_per_species = array(NA, 10)
# running_SSTr_sum = 0
# for (index in 1:10){
#   current_species = unique(stats_precipplusone_yield_matrix$species)[index]
#   number_of_samples_within_current_species = sum(!is.na(subset(stats_precipplusone_yield_matrix, species == current_species)$individual_squared_difference_from_mean_slope_on_log_scale))
#   vector_of_number_of_samples_per_species[index] = number_of_samples_within_current_species
#   running_SSTr_sum = running_SSTr_sum + number_of_samples_within_current_species*(vector_of_mean_slopes_on_log_scale[index] - mean(stats_precipplusone_yield_matrix$individual_slopes_on_log_scale[!is.na(stats_precipplusone_yield_matrix$individual_slopes_on_log_scale)]))^2
# }
# MSTr = running_SSTr_sum/(num_species-1)
# 
# 
# statistic_value_F = MSTr/MSE
# threshold_T_value = qalpha*sqrt(MSE/min(vector_of_number_of_samples_per_species))


######################################################################
# Grabbing coefficients from the post-hurdle models created earlier###
######################################################################
coefficients_of_conditional_linear_model = coef(glm_model_part_two_sub1)
coefficients_of_simplest_conditional_linear_model = coef(glm_model_part_two_sub2)
coefficients_of_full_conditional_linear_model = coef(glm_model_part_two_full)
use_table$alpha2s_full_model = NA
use_table$beta2s_full_model = NA
use_table$alpha2s_precip_plus_interaction = as.numeric(coefficients_of_conditional_linear_model['(Intercept)']) # create new column to contain the new binomial model
use_table$beta2s_precip_plus_interaction = NA
use_table$beta2s_precip_plus_species = as.numeric(coefficients_of_simplest_conditional_linear_model['log(precipitation)'])
for (j in 1:10){
  use_table$beta2s_precip_plus_interaction[j] = coefficients_of_conditional_linear_model['log(precipitation)'] + coefficients_of_conditional_linear_model[paste('log(precipitation):species', row.names(use_table)[j], sep = '')]
  use_table$alpha2s_precip_plus_species[j] = coefficients_of_simplest_conditional_linear_model['(Intercept)']+coefficients_of_simplest_conditional_linear_model[paste('species', row.names(use_table)[j], sep = '')]
  use_table$alpha2s_full_model[j] = coefficients_of_full_conditional_linear_model[['(Intercept)']] + coefficients_of_full_conditional_linear_model[paste('species', row.names(use_table)[j], sep = '')]
  use_table$beta2s_full_model[j] = coefficients_of_full_conditional_linear_model[['log(precipitation)']] + coefficients_of_full_conditional_linear_model[paste('species', row.names(use_table)[j], ':log(precipitation)', sep = '')]
} 

use_table$alpha2s_full_model[2] = coefficients_of_full_conditional_linear_model[['(Intercept)']] 
use_table$beta2s_full_model[2] = coefficients_of_full_conditional_linear_model[['log(precipitation)']]
use_table$beta2s_precip_plus_interaction[2] = as.numeric(coefficients_of_conditional_linear_model['log(precipitation)']) # ERCI was the reference point for the model, so we have to input this manually
use_table$alpha2s_precip_plus_species[2] = coefficients_of_simplest_conditional_linear_model['(Intercept)']
# use_table$beta2s_precip_only = coef(glm_model_part_two_sub3)[[2]]
# use_table$alpha2s_precip_only = coef(glm_model_part_two_sub3)[[1]]

########################################################################################
# DONE MODIFYING THE TABLE THAT WE WILL USE FOR PREDICTING ESS... ONTO PREDICTING ESS:##
########################################################################################

#######################################################
# 1) FITTING PRECIPITATION VALUES TO A LOG-NORMAL CURVE
# 2) GENERATING PRECIPITATION VALUES FROM DRAWING FROM THIS DISTRIBUTION
# 3) SETTING RHO = 0 (PRECIPITATION VALUES ARE NOT CORRELATED ACROSS TIME)
#######################################################

lognormal_model_for_precip = fitdist(log(precip + 1), 'norm') # creating log-normal distribution for precipitation values
m1 = as.numeric(lognormal_model_for_precip$estimate[1]) # finding mean of this log-normal distribution for precip values
sd1 = as.numeric(lognormal_model_for_precip$estimate[2]) # finding standard deviation of this log-normal distribution for precip values
matrix_of_bases = matrix(rnorm(t.mat*amt), nrow = t.mat, ncol = amt) # generating a matrix distributed as N(0,1)
P = exp(base*sd1 + m1)
rho = 0
P_auto = auto_cor_func(P,rho) # autocorrelate these values (if rho = 0, no autocorrelation)

# Previously used to plot log yields against log precipitation
# Ps=seq(min(P),max(P),length=100) # For preparation of the _logprecip.pdf
# Ys=matrix(0,100,k) # For preparation of the _logprecip.pdf
# for(i in 1:k){
#   Ys[,i]= exp(use_table$alphas[i])*Ps^(use_table$betas[i])
# } # creating the Yields to be put in the log-log


print(paste('Creating pdfs for rho', rho))
file.name=paste('Desert_Annual_MAY7_testing','Time =', Tf,"delta =", delta, "rho =", rho) 

###############################################################
# Creating a new P (like at the beginning) -- more variance ###
###############################################################

# We need to generate new P values with increased standard deviation,
# run dynamics long enough, so that we have autocorrelated Yield values that correspond to new densities,
# and then we have to use the OLD p-matrix to keep the autocorrelation going, so that
# we have some consistency. 

ess.vec.baseline = array(1,k) # Evolutionary Stable Strategy vector 
sd.change.vec = seq(-1, 1, length = len_of_sdvec_bsvec) 
ext.rate.matrix.postESS.SD = array(NA, c(length(sd.change.vec), k))
ext.rate.matrix.preESS.SD = array(NA, c(length(sd.change.vec), k))
ess.vals.matrix = array(NA, c(length(sd.change.vec), k))
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

###########################################################################
# FINDING OUR EVOLUTIONARILY STABLE STRATEGIES FOR THE TEN DESERT ANNUALS #
###########################################################################

for (ii in 1:k){
  # temp_func = function(x){return(one_or_zero_given_b0_and_b1_func(use_table$LH_b0s[ii], use_table$LH_b1s[ii], x))} # regressing one at a time coefficients
  temp_func = function(x){return(one_or_zero_given_b0_and_b1_func(use_table$LH_b0s_full_binomial_model[ii], use_table$LH_b1s_full_binomial_model[ii], x))} # full model from binomial (glm)
  # temp_func = function(x){return(one_or_zero_given_b0_and_b1_func(use_table$LH_b0s_precip_plus_interaction[ii], use_table$LH_b1s_precip_plus_interaction[ii], x))}
  # temp_func = function(x){return(one_or_zero_given_b0_and_b1_func(use_table$LH_b0s_multi_intercepts_same_slope[ii], use_table$LH_b1s_multi_intercepts_same_slope[ii], x))}
  
  # Ydata_baseline = sapply(log(P_auto), temp_func)*(exp(use_table$alphas[ii])*P_auto^use_table$betas[ii]) # regressing one at a time
  Ydata_baseline = sapply(log(P_auto), temp_func)*(exp(use_table$alpha2s_full_model[ii])*P_auto^use_table$beta2s_full_model[ii])
  # Ydata_baseline = sapply(log(P_auto), temp_func)*(exp(use_table$alpha2s_precip_plus_interaction[ii])*P_auto^use_table$beta2s_precip_plus_interaction[ii])
  # Ydata_baseline = sapply(log(P_auto), temp_func)*(exp(use_table$alpha2s_precip_plus_species[ii])*P_auto^use_table$beta2s_precip_plus_species[ii])
  # Ydata_baseline = sapply(log(P_auto), temp_func)*(exp(use_table$alpha2s_precip_only[ii])*P_auto^use_table$beta2s_precip_only[ii])
  
  Ydata_baseline[which(Ydata_baseline == 0)] = 0.5 # those that are 0, we treat as 0.5 as done in literature
  Ydata_baseline = Ydata_baseline*(use_table$survs_new[ii]) # we multiply by the probability of survival to the next fall
  our_yield_sds[ii] = sd.real(Ydata_baseline) # grabbing our own standard deviations in yields
  our_yield_means[ii] = mean(Ydata_baseline) # grabbing our own mean yields
  
  ## JENNY'S YIELD -- just for double-checking things
  # Ydata_baseline_jenny = sample(Kt_matrix3_dep[ii,][-list_indices_had_NA[[ii]]], Tf, replace=TRUE)*use_table$survs_new[ii]
  # jenny_yield_sds[ii] = sd.real(Ydata_baseline_jenny)
  # jenny_yield_means[ii] = mean(Ydata_baseline_jenny)
  
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

########################################
### PREDICTORS AND RESPONSE VARIABLES###
########################################

sold = use_table$survs
snew = use_table$survs_new
alpha1s = use_table$LH_b0s
beta1s = use_table$LH_b1s
alpha1_precip_plus_interaction = use_table$LH_b0s_precip_plus_interaction
beta1_precip_plus_interaction = use_table$LH_b1s_precip_plus_interaction
alpha2_precip_plus_species = use_table$alpha2s_precip_plus_species
alpha1s_from_multi_intercept_one_slope_model = use_table$LH_b0s_multi_intercepts_same_slope
CH_taus = use_table$taus
beta2s = use_table$betas
alpha2s = use_table$alphas
alpha1s_from_full_binomial_model = use_table$LH_b0s_full_binomial_model
beta1s_from_full_binomial_model = use_table$LH_b1s_full_binomial_model
alpha2s_from_full_posthurdle_model = use_table$alpha2s_full_model
beta2s_from_full_posthurdle_model = use_table$beta2s_full_model
beta2_precip_plus_interaction = use_table$beta2s_precip_plus_interaction
field_ess_vals = use_table$fieldESS
pred_ess_vals = ess.vec.baseline
transfield_ess_vals = logit(field_ess_vals)
transpred_ess_vals = logit(pred_ess_vals)
stan_alpha1s = standardize(alpha1s)
stan_alpha2s = standardize(alpha2s)
stan_beta2s = standardize(beta2s)
stan_beta1s = standardize(beta1s)
stan_snew = standardize(snew)
stan_sold = standardize(sold)
stan_beta1_precip_plus_interaction = standardize(beta1_precip_plus_interaction)
stan_beta2_precip_plus_interaction = standardize(beta2_precip_plus_interaction)
stan_alpha2_precip_plus_species = standardize(alpha2_precip_plus_species)
stan_alpha1s_from_multi_intercept_one_slope_model = standardize(alpha1s_from_multi_intercept_one_slope_model)
stan_alpha1s_from_full_binomial_model = standardize(alpha1s_from_full_binomial_model)
stan_beta2s_from_full_posthurdle_model = standardize(beta2s_from_full_posthurdle_model)
stan_beta1s_from_full_binomial_model = standardize(beta1s_from_full_binomial_model)
stan_alpha2s_from_full_posthurdle_model = standardize(alpha2s_from_full_posthurdle_model)

############################################
## ASSIMILATING PREDICTORS INTO DATA FRAME## 
############################################

data_table_with_params = cbind(sold, snew, beta2s, alpha1s, beta1s, alpha2s, beta1_precip_plus_interaction, beta2_precip_plus_interaction, alpha2_precip_plus_species, alpha1s_from_multi_intercept_one_slope_model, alpha1s_from_full_binomial_model, beta2s_from_full_posthurdle_model, alpha2s_from_full_posthurdle_model, beta1s_from_full_binomial_model)
data_table_with_params = as.data.frame(data_table_with_params)

################################################################
### REGRESSING PREDICTED ESS VALUES AND ESTIMATED ESS VALUES ### 
################################################################
# observedESS_which_model_to_choose = regsubsets(transfield_ess_vals ~ stan_alpha1s_from_multi_intercept_one_slope_model + stan_sold + stan_alpha2_precip_plus_species, data = data_table_with_params, nbest = 4, nvmax = 14, method = c("exhaustive"))
# predictedESS_which_model_to_choose = regsubsets(transpred_ess_vals ~ stan_alpha1s_from_multi_intercept_one_slope_model + stan_sold + stan_alpha2_precip_plus_species, data = data_table_with_params, nbest = 4, nvmax = 14, method = c("exhaustive"))
# plot(predictedESS_which_model_to_choose)
# plot(observedESS_which_model_to_choose)

full_full_predicting_field_ESS_no_interactions = lm(transfield_ess_vals ~ stan_sold + stan_beta2s_from_full_posthurdle_model + stan_alpha1s_from_full_binomial_model) # these are the terms we should be using to stay consistent with our ANOVA analysis
full_full_predicting_predicted_ESS_no_interactions = lm(transpred_ess_vals ~ stan_sold + stan_beta2s_from_full_posthurdle_model + stan_alpha1s_from_full_binomial_model)
full_full_field_with_interactions = aov(transfield_ess_vals ~ stan_sold + stan_beta2s_from_full_posthurdle_model*stan_alpha1s_from_full_binomial_model)
full_full_field_without_interactions = aov(transfield_ess_vals ~ stan_sold + stan_beta2s_from_full_posthurdle_model + stan_alpha1s_from_full_binomial_model)

summary(lm(logit(use_table$fieldESS) ~ logit(ess.vec.baseline)))
summary(lm(logit(use_table$fieldESS) ~ logit(use_table$jenESS)))
summary(full_full_predicting_field_ESS_no_interactions)
summary(full_full_predicting_predicted_ESS_no_interactions)

transfield_model1 = lm(transfield_ess_vals ~ stan_alpha2_precip_plus_species + stan_sold*stan_beta1s)
transfield_model2 = lm(transfield_ess_vals ~ stan_alpha2_precip_plus_species + stan_sold + stan_beta1s)
transfield_model3 = lm(transfield_ess_vals ~ stan_alpha2_precip_plus_species + stan_sold)
transfield_model4 = lm(transfield_ess_vals ~ stan_alpha2_precip_plus_species*stan_sold)
anova(transfield_model3, transfield_model4)

#####################################################################################################
#### BACK-TRANSFORMING STATISTICAL REGRESSION RESULTS TO COMPARE WITH OBSERVED GERMINATION RATES ####
#####################################################################################################

stats_model = lm(transfield_ess_vals ~ stan_sold + stan_alpha1s_from_full_binomial_model + stan_beta2s_from_full_posthurdle_model)
predicting_logit_ESS = coef(stats_model)[[1]] + coef(stats_model)[[2]]*stan_sold + coef(stats_model)[[3]]*stan_alpha1s_from_full_binomial_model + coef(stats_model)[[4]]*stan_beta2s_from_full_posthurdle_model
predicting_ess = inverse_logit(predicting_logit_ESS)

###############################################################################################################
### Taken from earlier for loop: create figure for species of interest in file 17.10.11_prez.R              ###
### This had to be called last because it needs to be these redone variables for 18.01.21_paper1_figures.R  ###
###############################################################################################################


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

# some precipitation stuff:
first_62_years_mean_prec = mean((desert_station_precip_1922_to_2017$JMprecipMM[-which(is.na(desert_station_precip_1922_to_2017), arr.ind = TRUE)]*.1)[31:62] +1)
first_62_years_sd_prec = sd.real((desert_station_precip_1922_to_2017$JMprecipMM[-which(is.na(desert_station_precip_1922_to_2017), arr.ind = TRUE)]*.1)[31:62] +1)
our_years_mean_prec = mean((desert_station_precip_1922_to_2017$JMprecipMM[62:97] +1))
our_years_sd_prec = sd.real((desert_station_precip_1922_to_2017$JMprecipMM[62:97] +1))

percent_change_in_mean_prec = (our_years_mean_prec - first_62_years_mean_prec)/first_62_years_mean_prec
percent_change_in_sd_prec = (our_years_sd_prec - first_62_years_sd_prec)/first_62_years_sd_prec
