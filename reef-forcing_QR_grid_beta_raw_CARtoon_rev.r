require(readr)
require(dplyr)
require(cowplot)
require(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 6)

source('get_data_CAR_rev.r')
source('plots.r')

ISLANDS <- c('KUR','MID','PHR','LIS','LAY','MAR','FFS','NEC','NIH','KAU','NII','OAH','MOL','LAN','MAI','HAW','JOH','WAK','FDP','MAU','ASC','AGR','PAG','ALA','GUG','SAR','SAI','TIN','AGU','ROT','GUA','KIN','PAL','HOW','BAK','JAR','SWA','OFU','TAU','TUT','ROS')

# #
# # ## MACROALGAE
# #
# VAR <- 'MACRO'
# model_data <- get_vars(all_data,VAR)
# 
# save(model_data,file='CAR_macro_rev.Rdata')


load('CAR_macro_rev.Rdata')
fit_MARtoon_rev <- stan('model_QR_grid_beta_raw_CARtoon_rev.stan',
                        data = model_data,
                        iter = 2000,
                        chains = 6,
                        thin = 4,
                        warmup = 1000,
                        init = function() list(rho=rep(0.95,model_data$N_islands)),
                        pars = c('island_sig',
                                 'group_sig',
                                 'grid_sig',
                                 'rho',
                                 'resid_y',
                                 'beta',
                                 't_mean',
                                 't_var',
                                 'iscale',
                                 'scale',
                                 #'log_lik',
                                 'omean',
                                 'zmean'),
                        control = list(max_treedepth = 13, adapt_delta=0.8),
                        verbose = T)

save(fit_MARtoon_rev,file='CARtoonfit_raw_macroalgae_CARfull.Rdata')
# 
# gg <- stan_trace(fit_MARtoon_rev,
#           pars=c('omean','rho[14]','grid_sig','beta[1,1]',
#                  't_var','t_mean','lp__','scale','zmean','island_sig','group_sig'),
#           #pars=c('rho_prec'),
#           include=T)
# # ggsave(gg,file='traces.pdf')
# 
# plot_betas(fit_MARtoon_rev, model_data)
# # ggsave(gg,file='betas.pdf')
# 
# plot_fit(fit_MARtoon_rev, model_data)
# 
# plot_resids(fit_MARtoon_rev, model_data, ac=T)

# #


## CCA

#' VAR <- 'CCA'
#' model_data <- get_vars(all_data,VAR)
#' save(model_data,file='CAR_cca_rev.Rdata')
#' 
#' #' 
load('CAR_cca_rev.Rdata')
#' 
fit_CCARtoon_rev <- stan('model_QR_grid_beta_raw_CARtoon_rev.stan',
                         data = model_data,
                         iter = 2000,
                         chains = 6,
                         thin = 4,
                         warmup = 1000,
                         init = function() list(rho=rep(0.95,model_data$N_islands)),
                         pars = c('island_sig',
                                  'group_sig',
                                  'grid_sig',
                                  'rho',
                                  'resid_y',
                                  'beta',
                                  't_mean',
                                  't_var',
                                  'iscale',
                                  'scale',
                                  #'log_lik',
                                  'omean',
                                  'zmean'),
                         control = list(max_treedepth = 13, adapt_delta=0.8),
                         verbose = T)

save(fit_CCARtoon_rev,file='CARtoonfit_raw_cca_CARfull.Rdata')
# 
# plot_betas(fit_CCARtoon_rev, model_data)

# #
# ggsave(gg,file='betas.pdf')


#' 
#' 
#' ## CORAL
#' 
#' VAR <- 'CORAL'
#' model_data <- get_vars(all_data,VAR)
#' 
#' save(model_data,file='CAR_coral_rev.Rdata')
#' 
#' 
load('CAR_coral_rev.Rdata')
fit_CARtoon_rev <- stan('model_QR_grid_beta_raw_CARtoon_rev.stan',
                        data = model_data,
                        iter = 2000,
                        chains = 6,
                        thin = 4,
                        warmup = 1000,
                        init = function() list(rho=rep(0.95,model_data$N_islands)),
                        pars = c('island_sig',
                                 'group_sig',
                                 'grid_sig',
                                 'rho',
                                 'resid_y',
                                 'beta',
                                 't_mean',
                                 't_var',
                                 'iscale',
                                 'scale',
                                 #'log_lik',
                                 'omean',
                                 'zmean'),
                        control = list(max_treedepth = 13, adapt_delta=0.8),
                        verbose = T)

save(fit_CARtoon_rev,file='CARtoonfit_raw_coral_CARfull.Rdata')
#
# # 
# 
# gg <- stan_trace(fit_CARtoon_rev,
#                  pars=c('omean','rho[25]','grid_sig','beta[25,1]',
#                         't_var','t_mean','lp__','scale','zmean','island_sig','group_sig'),
#                  #pars=c('rho_prec'),
#                  include=T)
# # ggsave(gg,file='traces.pdf')
# 
# plot_betas(fit_CARtoon_rev, model_data)
# get_p_beta(fit_CARtoon_rev, model_data)
