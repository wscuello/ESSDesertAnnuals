# Dependencies: 18.07.15_bethedging.R,  18.01.21_paper1_calledon.R

####################################
## FIGURE 1B/C AND 2A/B FOR PAPER ##
####################################

pdf('Figure1B.pdf', height = 7, width = 7)
{par(mar=c(7,7,2,1))
  plot(log_precip_plus1, replace_or_not_vec, ylab = 'Probability of Reproduction', xlab = 'log(Precip + 1)', main = 'ERLA', cex = 4, cex.axis = 2.5, cex.lab = 3, cex.main = 2, lwd = 1, col = the_colors[5])
  domain_for_logistic = seq(min(log_precip_plus1), max(log_precip_plus1), by = .1)
  temp_logistic_func = function(beta0, beta1, x){return(1/(1+exp(-(beta0 + beta1*x))))}
  lines(domain_for_logistic, temp_logistic_func(coef(logistic_hockey)[[1]],coef(logistic_hockey)[[2]],domain_for_logistic), col = the_colors[5], lwd = 5)
}
dev.off()

pdf('Figure1C.pdf', height = 7, width = 7)
{par(mar=c(7,7,2,1))
  plot(log_precip_plus1[which_log_ks_to_regress], log_k_vec[which_log_ks_to_regress], ylim = c(-1, 8), xlim = c(-.5, 5), xlab = 'log(Precip + 1)', ylab = 'log(Yield + .5)', main = 'ERLA', cex = 4, cex.main = 2, cex.lab = 3, cex.axis = 2.5, lwd = 1, col = the_colors[5])
  domain = seq(-5, 5, by = .01)
  abline(conditioned_ls, lty = 1, lwd = 5, col = the_colors[5])
  abline(v = 0, h = 0, lty = 2, lwd = 3)
}
dev.off()

pdf('Figure2A.pdf', width = 8, height = 8)
{par(mar=c(7,7,2,1))
  plot(ess.vec.baseline, use_table$fieldESS, ylim = c(0,1), xlim = c(0,1), xlab = 'Estimated ESS', ylab = 'Observed Germination Fractions', cex.lab = 2, cex.axis = 1.5, cex = 4, col = the_colors, lwd = 4)
  text(ess.vec.baseline, use_table$fieldESS, labels = row.names(use_table), pos = c(4,3,1,2,4,4,2,3,3,1), col = the_colors, cex = 1.5)
  abline(c(0,0), c(1,1), lty =2, lwd = 3, col = 'black')
  abline(lm(use_table$fieldESS ~ ess.vec.baseline), col = 'black', lwd = 2)
  # mtext(paste0('A'), side = 1, adj = .95, line = -9.0, cex = 3.5)}
}
dev.off()

pdf('Figure2B.pdf', width = 8, height = 8)
{par(mar=c(7,7,2,1))
  plot(predicting_ess, use_table$fieldESS, col = the_colors, cex = 4, cex.lab= 2, cex.axis = 1.5, ylab = 'Observed Germination Fractions', xlab = expression('logit'^-1*'(L'[obs]*')'), ylim = c(0,1), xlim = c(0,1), lwd = 4)
  text(predicting_ess, use_table$fieldESS, labels = row.names(use_table), col = the_colors, pos = c(4,3,1,2,4,4,2,3,2,1), cex = 1.5)
  abline(c(0,0),c(1,1), lty = 2, lwd = 3)
  abline(lm(use_table$fieldESS ~ predicting_ess), lwd = 2)
  # mtext(paste0('B'), side = 1, adj = .95, line = -9.0, cex = 3.5)
}
dev.off()

pdf('Figure2Aalternative.pdf', width = 8, height = 8)
{par(mar=c(7,7,2,1))
  plot(ess.vec.baseline, use_table$fieldESS, ylim = c(0,1), xlim = c(0,1), xlab = 'Estimated ESS', ylab = 'Observed Germination Fractions', cex.lab = 2, cex.axis = 1.5, cex = 4, col = the_colors, lwd = 4)
  text(ess.vec.baseline, use_table$fieldESS, labels = row.names(use_table), pos = c(4,3,1,2,4,4,2,3,3,1), col = the_colors, cex = 1.5)
  abline(c(0,0), c(1,1), lty =2, lwd = 3, col = 'black')
  # abline(lm(use_table$fieldESS ~ ess.vec.baseline), col = 'black', lwd = 2)
  # mtext(paste0('A'), side = 1, adj = .95, line = -9.0, cex = 3.5)}
}
dev.off()

pdf('Figure2Balternative.pdf', width = 8, height = 8)
{par(mar=c(7,7,2,1))
  plot(predicting_ess, use_table$fieldESS, col = the_colors, cex = 4, cex.lab= 2, cex.axis = 1.5, ylab = 'Observed Germination Fractions', xlab = expression('logit'^-1*'(L'[obs]*')'), ylim = c(0,1), xlim = c(0,1), lwd = 4)
  text(predicting_ess, use_table$fieldESS, labels = row.names(use_table), col = the_colors, pos = c(4,3,1,2,4,4,2,3,2,1), cex = 1.5)
  abline(c(0,0),c(1,1), lty = 2, lwd = 3)
  # abline(lm(use_table$fieldESS ~ predicting_ess), lwd = 2)
  # mtext(paste0('B'), side = 1, adj = .95, line = -9.0, cex = 3.5)
}
dev.off()


#############################################################
# EXPLORED SOME OTHER THINGS THAT DIDN'T GO INTO THE PAPER###
#############################################################

# FIGURE 1 -- Just hit CTRL+ENTER here
# LEFT PANEL: Observed germination fractions plotted against ESS values attained from bet-hedging model and invasion analysis
# RIGHT PANEL: Observed germination fractions plotted against germination estimates attained from statistical model based on physiological traits
# pdf('Figure2.pdf', width = 7, height = 7)
# {par(cex.lab=1.25,cex.axis=1.25,mar=c(4.5,4.5,4.5,4.5),mfrow=c(2,1))
# plot(ess.vec.baseline, use_table$fieldESS, ylim = c(0,1), xlim = c(0,1), xlab = 'Estimated ESS', ylab = 'Observed Germination', cex.lab = 1.5, cex.axis = 1.5, cex = 2, col = the_colors, lwd = 2)
# text(ess.vec.baseline, use_table$fieldESS, labels = row.names(use_table), pos = c(4,3,1,2,4,4,2,3,4,1), col = the_colors)
# abline(c(0,0), c(1,1), lty =2, lwd = 3, col = 'black')
# abline(lm(use_table$fieldESS ~ ess.vec.baseline), col = 'black', lwd = 2)
# plot(predicting_ess, use_table$fieldESS, col = the_colors, cex = 2, cex.lab=1.5, cex.axis = 1.5, ylab = 'Observed Germination', xlab = expression('logit'^-1*'(L'[obs]*')'), ylim = c(0,1), xlim = c(0,1), lwd = 2)
# text(predicting_ess, use_table$fieldESS, labels = row.names(use_table), col = the_colors, pos = c(4,3,1,2,4,4,2,3,4,1))
# abline(c(0,0),c(1,1), lty = 2, lwd = 3)
# abline(lm(use_table$fieldESS ~ predicting_ess), lwd = 2)}
# dev.off()



pdf('Figure_exploration0_paper1.pdf')
{
  morphed_1_linear_model = .87 + .23*stan_alpha1s_from_full_binomial_model - .31*stan_beta2s_from_full_posthurdle_model - .54*stan_sold
  plot(logit(predicting_ess), logit(use_table$fieldESS), cex = 2, col = the_colors, lwd = 3, ylim = c(-3,2.5), xlim = c(-2.5, 1.5), xlab='pred (circle) & ess (triangle)', ylab = 'observed', main = 'estimates on logit scale')
  text(logit(predicting_ess), logit(use_table$fieldESS), labels = row.names(use_table), pos = c(4,3,1,2,4,4,2,3,4,1), col = the_colors)
  lines(morphed_1_linear_model, logit(use_table$fieldESS), type = 'p', pch = 2, cex =2, col = rgb(0,0,0,alpha=0.5), lwd =2)
  text(morphed_1_linear_model, logit(use_table$fieldESS), labels = row.names(use_table), pos = c(4,3,1,2,4,4,2,3,4,1), col = rgb(0,0,0,alpha=0.7))
  abline(c(0,0), c(1,1), lty = 2, lwd = 2)
}
dev.off()

pdf('Figure_exploration1_paper1.pdf')
{
  morphed_1_linear_model = -.29 + .23*stan_alpha1s_from_full_binomial_model - .31*stan_beta2s_from_full_posthurdle_model - .54*stan_sold
  plot(logit(predicting_ess), logit(use_table$fieldESS), cex = 2, col = the_colors, lwd = 3, ylim = c(-3,2.5), xlim = c(-2.5, 1.5), xlab='pred (circle) & ess (triangle)', ylab = 'observed', main = 'estimates on logit scale')
  text(logit(predicting_ess), logit(use_table$fieldESS), labels = row.names(use_table), pos = c(4,3,1,2,4,4,2,3,4,1), col = the_colors)
  lines(morphed_1_linear_model, logit(use_table$fieldESS), type = 'p', pch = 2, cex =2, col = rgb(0,0,0,alpha=0.5), lwd =2)
  text(morphed_1_linear_model, logit(use_table$fieldESS), labels = row.names(use_table), pos = c(4,3,1,2,4,4,2,3,4,1), col = rgb(0,0,0,alpha=0.7))
  abline(c(0,0), c(1,1), lty = 2, lwd = 2)
}
dev.off()

pdf('Figure_exploration2_paper1.pdf')
{
  morphed_2_linear_model = -.29 + .46*stan_alpha1s_from_full_binomial_model - .31*stan_beta2s_from_full_posthurdle_model - .54*stan_sold
  plot(logit(predicting_ess), logit(use_table$fieldESS), cex = 2, col = the_colors, lwd = 3, ylim = c(-3,2.5), xlim = c(-2.5, 1.5), xlab='pred (circle) & ess (triangle)', ylab = 'observed', main = 'estimates on logit scale')
  text(logit(predicting_ess), logit(use_table$fieldESS), labels = row.names(use_table), pos = c(4,3,1,2,4,4,2,3,4,1), col = the_colors)
  lines(morphed_2_linear_model, logit(use_table$fieldESS), type = 'p', pch = 2, cex =2, col = rgb(0,0,0,alpha=0.5), lwd =2)
  text(morphed_2_linear_model, logit(use_table$fieldESS), labels = row.names(use_table), pos = c(4,3,1,2,4,4,2,3,4,1), col = rgb(0,0,0,alpha=0.7))
  abline(c(0,0), c(1,1), lty = 2, lwd = 2)
}
dev.off()

pdf('Figure_exploration3_paper1.pdf')
{
  morphed_3_linear_model = -.29 + .46*stan_alpha1s_from_full_binomial_model - .45*stan_beta2s_from_full_posthurdle_model - .54*stan_sold
  plot(logit(predicting_ess), logit(use_table$fieldESS), cex = 2, col = the_colors, lwd = 3, ylim = c(-3,2.5), xlim = c(-2.5, 1.5), xlab='pred (circle) & ess (triangle)', ylab = 'observed', main = 'estimates on logit scale')
  text(logit(predicting_ess), logit(use_table$fieldESS), labels = row.names(use_table), pos = c(4,3,1,2,4,4,2,3,4,1), col = the_colors)
  lines(morphed_3_linear_model, logit(use_table$fieldESS), type = 'p', pch = 2, cex =2, col = rgb(0,0,0,alpha=0.5), lwd =2)
  text(morphed_3_linear_model, logit(use_table$fieldESS), labels = row.names(use_table), pos = c(4,3,1,2,4,4,2,3,4,1), col = rgb(0,0,0,alpha=0.7))
  abline(c(0,0), c(1,1), lty = 2, lwd = 2)
}
dev.off()

pdf('Figure_exploration4_paper1.pdf')
{
  morphed_4_linear_model = -.29 + .46*stan_alpha1s_from_full_binomial_model - .45*stan_beta2s - .60*stan_sold
  plot(logit(predicting_ess), logit(use_table$fieldESS), cex = 2, col = the_colors, lwd = 3, ylim = c(-3,2.5), xlim = c(-2.5, 1.5), xlab='pred (circle) & ess (triangle)', ylab = 'observed', main = 'estimates on logit scale')
  text(logit(predicting_ess), logit(use_table$fieldESS), labels = row.names(use_table), pos = c(4,3,1,2,4,4,2,3,4,1), col = the_colors)
  lines(morphed_4_linear_model, logit(use_table$fieldESS), type = 'p', pch = 2, cex =2, col = rgb(0,0,0,alpha=0.5), lwd =2)
  text(morphed_4_linear_model, logit(use_table$fieldESS), labels = row.names(use_table), pos = c(4,3,1,2,4,4,2,3,4,1), col = rgb(0,0,0,alpha=0.7))
  abline(c(0,0), c(1,1), lty = 2, lwd = 2)
}
dev.off()



# FIGURE 2 -- Just hit CTRL + ENTER here
# TOP 2 PANELS: Logistic regression (reproductive success), Line of regression (Yield upon Success)
# BOTTOM 3 PANELS: beta_1 vs. alpha_1 (factors for reproductive success with more precipitation vs. on low precipitation years)
#                  beta_2 vs. alpha_2 (factors for capitalizing on bonzanzas vs. producing more upon low precipitation years)
#                  s_new vs. s_old (survivorship of fresh seeds vs. survivorship of seeds from the seed bank)
# {par(mar=c(4.5,4.5,4.5,4.5), mfrow=c(3,2))
#   plot(log_precip_plus1, replace_or_not_vec, ylab = 'Logistic', xlab = 'log(Precip+1)', main = 'ERLA', cex = 1.5, cex.axis = 2, cex.lab = 2, col = "#00008B")
#   domain_for_logistic = seq(min(log_precip_plus1), max(log_precip_plus1), by = .1)
#   temp_logistic_func = function(beta0, beta1, x){return(1/(1+exp(-(beta0 + beta1*x))))}
#   lines(domain_for_logistic, temp_logistic_func(coef(logistic_hockey)[[1]],coef(logistic_hockey)[[2]],domain_for_logistic), col = "#00008B", lwd = 2)
#   mtext(paste0("(", 'A', ")"), side = 4, adj = 0.05, line = -1.3)
#   plot(log_precip_plus1[which_log_ks_to_regress], log_k_vec[which_log_ks_to_regress], ylim = c(-1, 8), xlim = c(-.5, 5), xlab = 'log(Precip + 1)', ylab = 'log(Yield + .5)', main = 'ERLA', cex = 1.5, cex.lab = 2, cex.axis = 2, col = "#00008B")
#   domain = seq(-5, 5, by = .01)
#   abline(conditioned_ls, lty = 2, lwd =2, col = "#00008B")
#   abline(v = 0, h = 0, lty = 4)
#   mtext(paste0("(", 'B', ")"), side = 4, adj = 0.05, line = -1.3)
#   plot(use_table$LH_b0s, use_table$LH_b1s, ylim = c(min(use_table$LH_b1s),max(use_table$LH_b1s)), xlim = c(min(use_table$LH_b0s),max(use_table$LH_b0s)), cex = 2, col = the_colors, ylab = expression(beta[1]), xlab = expression(alpha[1]), cex.lab = 2, cex.axis = 2, lwd = 2)
#   abline((lm(use_table$LH_b1s ~ use_table$LH_b0s)), lty = 2, lwd = 2)
#   mtext(paste0("(", 'C', ")"), side = 2, adj = 0.05, line = -1.3)
#   plot(use_table$alphas, use_table$betas, ylim = c(min(use_table$betas),max(use_table$betas)), xlim = c(min(use_table$alphas),max(use_table$alphas)), cex = 2, col = the_colors, ylab = expression(beta[2]), xlab = expression(alpha[2]), cex.lab = 2, cex.axis = 2, lwd = 2)
#   abline((lm(use_table$betas ~ use_table$alphas)), lty = 2, lwd = 2)
#   mtext(paste0("(", 'D', ")"), side = 2, adj = 0.05, line = -1.3)
#   plot(use_table$survs, use_table$survs_new, ylim = c(min(use_table$survs_new),max(use_table$survs_new)), xlim = c(min(use_table$survs),max(use_table$survs)), cex = 2, col = the_colors, xlab = expression('s'[old]), ylab = expression('s'[new]), cex.lab = 2, cex.axis = 2, lwd = 2)
#   abline((lm(use_table$survs_new ~ use_table$survs)), lty = 2, lwd = 2)
#   mtext(paste0("(", 'E', ")"), side = 4, adj = 0.05, line = -1.3)}
# dev.off()


