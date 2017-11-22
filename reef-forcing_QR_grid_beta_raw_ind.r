
require(readr)
require(dplyr)
require(rlang)
require(cowplot)
require(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


get_vars <- function(data,var,sample=F){
  
  adata <- data %>% filter(N_OBS_CELL!=0,
                           !is.na(TED_SUM),
                           !is.na(SURF_SUM)) %>% 
    mutate_(mm = var) %>%
    mutate(mm=mm/100)
  
 
  
  if(sample == T) adata <- adata %>% sample_n(500)

  COVS <- with(adata, as.matrix(
    data.frame(TED_SUMS   = scale(TED_SUM)[,1],
               SURF_SUMS  = scale(SURF_SUM)[,1],
               PHYS_INT   = scale(TED_SUM*SURF_SUM)[,1])))
  
  #browser()
  
  with(adata,list(
    ZEROONE   = ifelse(mm>0,1,0),
    N_notzero = sum(mm>0),
    N_groups  = length(unique(group)),
    N_islands = length(unique(ISLAND)),
    N         = nrow(COVS),
    K         = ncol(COVS),
    ii_notzero= which(mm>0),
    IMEAN     = mm[mm>0],
    COVS      = COVS,
    ISLAND    = as.numeric(as.factor(ISLAND)),
    GROUP     = as.numeric(as.factor(group)),
    GRID      = as.numeric(as.factor(paste0(ISLAND,GRID_ID))),
    N_grid  = length(unique(as.numeric(as.factor(paste0(ISLAND,GRID_ID)))))
    ))
  
}

cdata <- read_csv('./model/data/ENTIRE_CLIMATOLOGY_ALL_ISLANDS.csv')
rdata <- read_csv('./model/data/RAW_BENTHIC_ALL_ISLANDS.csv')

all_data <- inner_join(rdata,cdata)

# MACROALGAE

VAR <- 'MACRO'
model_data <- get_vars(all_data,VAR,sample = F)

fit_M <- stan('model_QR_grid_beta_raw_ind.stan',
              data = model_data,
              iter = 800,
              chains = 3,
              thin = 1,
              warmup = 500,
              #init = inits,
              pars = c('island_sig',
                       'group_sig',
                       'grid_sigs',
                       'beta',
                       #'theta',
                       'iscale',
                       'scale',
                       'zmean',
                       'resid_y',
                       'log_lik',
                       'omean'),
              control = list(max_treedepth = 10, adapt_delta=0.8),
              verbose = T)

save(fit_M,file='fit_raw_macroalgae_ind.Rdata')

ry <- get_posterior_mean(fit_M,pars='resid_y')[,3]
rys <- data.frame(g = model_data$GRID[model_data$ii_notzero],ry,i=model_data$ISLAND[model_data$ii_notzero])
ggplot(rys) + geom_line(aes(x=g,y=ry)) + facet_wrap(~as.factor(i),scales='free')
# 
require(purrr)
acc <- rys %>% split(.$i) %>% map(~acf(.$ry,plot=F,type='cov',lag.max = 300)) %>% map(~with(.,data.frame(lag,acf))) %>% bind_rows(.id='Island')
gg <- ggplot(acc) + geom_line(aes(x=lag,y=acf)) + facet_wrap(~Island,scales='free')
ggsave(gg,file='spatial.pdf')
# 
# 

# 
gg <- stan_trace(fit_M,
          pars=c('omean','grid_sig','beta[1,3]','lp__','scale','zmean','island_sig','group_sig','gscale'),
          include=T)
ggsave(gg,file='traces.pdf')

beats <- extract(fit_M,'beta')$beta

dimnames(beats)[[3]] <- dimnames(model_data$COVS)[[2]]
dimnames(beats)[[2]] <- unique(all_data$ISLAND)

beat <- tbl_df(reshape2::melt(beats))
colnames(beat)[2:3] <- c('Island','Covariate')

gg <- ggplot(beat) +
  geom_violin(aes(x=Covariate,y=value)) +
  geom_hline(aes(yintercept=0)) +
  facet_wrap(~Island) + 
  theme_cowplot() +
  geom_vline(aes(xintercept=0), linetype=2, alpha=0.4)+
  xlab('Effect size') +
  ylab('')

ggsave(gg,file='betas.pdf')

# CCA
# 
# VAR <- 'CCA'
# model_data <- get_vars(all_data,VAR,sample = F)
# 
# fit_CCA <- stan('model_QR_grid_beta_raw_h.stan',
#               data = model_data,
#               iter = 1500,
#               chains = 4,
#               thin = 2,
#               warmup = 500,
#               #init = inits,
#               pars = c('island_sig',
#                        'group_sig',
#                        'grid_sig',
#                        'beta',
#                        'gscale',
#                        'iscale',
#                        'scale',
#                        'zmean',
#                        'log_lik',
#                        'resid_y',
#                        'omean'),
#               control = list(max_treedepth = 10, adapt_delta=0.8),
#               verbose = T)
# 
# save(fit_CCA,file='fit_raw_cca_h.Rdata')


# 
# stan_trace(fit_CCA, 
#           pars=c('log-posterior','log_lik','iscale'), 
#           include=F)
# 
# stan_plot(fit_CCA, 
#           pars=c('beta'), 
#           include=T, 
#           show_density = TRUE, 
#           fill_color = "lightblue") + 
#   theme_cowplot() + 
#   geom_vline(aes(xintercept=0), linetype=2, alpha=0.4)+
#   xlab('Effect size') + 
#   ylab('')

## CORAL
# 
# VAR <- 'CORAL'
# model_data <- get_vars(all_data,VAR,sample = F)
# 
# fit_C <- stan('model_QR_grid_beta_raw_h.stan', 
#               data = model_data, 
#               iter = 1500,
#               chains = 4,
#               thin = 2,
#               warmup = 500,
#               #init = inits,
#               pars = c('island_sig',
#                        'group_sig',
#                        'grid_sig',
#                        'beta',
#                        'gscale',
#                        'iscale',
#                        'scale',
#                        'log_lik',
#                        'resid_y',
#                        'zmean',
#                        'omean'), 
#               control = list(max_treedepth = 10, adapt_delta=0.8),
#               verbose = T)
# 
# save(fit_C,file='fit_raws_coral_h.Rdata')
# 


