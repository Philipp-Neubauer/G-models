
require(readr)
require(dplyr)
require(cowplot)
require(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 2)


get_vars <- function(data,var){
  #require(INLA)
  adata <- data %>% filter(!is.na(TED_SUM),
                           !is.na(SURF_SUM)) %>% 
    mutate_(mm = var) %>%
    mutate(mm=ifelse(N_OBS_CELL>0,mm/100,NA))
  
  pmax <- 0
  island <- as.numeric(as.factor(adata$ISLAND))
  islands <- unique(island)
  nbbs <- list()
  idds <- list()
  for (i in 1:length(islands)) {
    cat(i,'\n')
    ids  <- adata$GRID_ID[island==islands[i]]
    nbbs[[i]] <- lapply(1:length(unique(ids)), function(i) {
      id <- unique(ids)[i]
      if(id<max(ids)) {
        if(!any(ids==(id+1))) NULL
        else c(i,i+1)+pmax
      }
      else {c(1,i)+pmax}
      
    })
    
    idds[[i]] <- match(ids,unique(ids))+pmax
    pmax <- pmax + length(unique(ids))
    
  }
  
  GRID <- do.call(c,idds)
  N_grid <- length(unique(GRID))
  nbs <- do.call(rbind,do.call(c,nbbs))
  N_edges <- nrow(nbs)
  #browser()
  require(Matrix)
  W = sparseMatrix(i=nbs[,1],j=nbs[,2],x=1,symmetric=TRUE,check = F)
  
  COVS <- with(adata, as.matrix(
    data.frame(
      DPOP      = scale(distant_hum_pop),
      TED_SUM   = scale(TED_SUM)[,1],
      SURF_SUM  = scale(SURF_SUM)[,1],
      TED_INT   = scale(TED_SUM*distant_hum_pop)[,1],
      SURF_INT  = scale(SURF_SUM*distant_hum_pop)[,1])))
  
  #browser()
  
  with(adata,list(
    W         = as.matrix(W),
    W_n       = N_edges,
    ZEROONE   = ifelse(mm[!is.na(mm)]>0,1,0),
    N_notzero = sum(mm>0,na.rm=T),
    N_zeroone = sum(!is.na(mm)),
    N_groups  = length(unique(group)),
    N_islands = length(unique(ISLAND)),
    N         = nrow(COVS),
    K         = ncol(COVS),
    ii_notzero= which(mm>0),
    ii_zeroone= which(!is.na(mm)),
    ii_map    = which(mm[which(!is.na(mm))]>0),
    IMEAN     = mm[which(mm>0)],
    COVS      = COVS,
    ISLAND    = as.numeric(as.factor(ISLAND)),
    GROUP     = as.numeric(as.factor(group)),
    GRID      = GRID,
    N_grid    = length(unique(GRID))
  ))
  
}



cdata <- read_csv('./model/data/ENTIRE_CLIMATOLOGY_ALL_ISLANDS.csv')
rdata <- read_csv('./model/data/RAW_BENTHIC_ALL_ISLANDS.csv')

all_data <- inner_join(rdata,cdata)

#' 
#' ## MACROALGAE
#' 
#' VAR <- 'MACRO'
#' model_data <- get_vars(all_data,VAR)
#' 
#' fit_MAR <- stan('model_QR_grid_beta_raw_AR.stan', 
#'               data = model_data, 
#'               iter = 1000,
#'               chains = 2,
#'               thin = 1,
#'               warmup = 500,
#'               #init = inits,
#'               pars = c('island_sig',
#'                        'group_sig',
#'                        'sigma',
#'                        'grid_sig',
#'                        'resid_y',
#'                        'beta',
#'                        'theta',
#'                        'gscale',
#'                        'iscale',
#'                        'scale',
#'                        'zmean',
#'                        #'log_lik',
#'                        'omean'), 
#'               control = list(max_treedepth = 10, adapt_delta=0.8),
#'               verbose = T)
#' 
#' save(fit_MAR,file='ARfit_raw_macroalgae.Rdata')
#' # 
#' stan_trace(fit_MAR,
#'           pars=c('omean','sigma','grid_sig','beta','lp__','scale','zmean','island_sig','group_sig'),
#'           include=T)
# 
# stan_plot(fit_MAR, 
#           pars=c('beta'), 
#           include=T, 
#           show_density = TRUE, 
#           fill_color = "lightblue") + 
#   theme_cowplot() + 
#   geom_vline(aes(xintercept=0), linetype=2, alpha=0.4)+
#'   xlab('Effect size') + 
#'   ylab('')
#' 
#' ## CCA
#' 
#' VAR <- 'CCA'
#' model_data <- get_vars(all_data,VAR)
#' 
#' fit_CCAR <- stan('model_QR_grid_beta_raw_AR.stan', 
#'               data = model_data, 
#'               iter = 1000,
#'               chains = 6,
#'               thin = 1,
#'               warmup = 500,
#'               #init = inits,
#'               pars = c('island_sig',
#'                        'group_sig',
#'                        'sigma',
#'                        'grid_sig',
#'                        'beta',
#'                        'gscale',
#'                        'iscale',
#'                        'scale',
#'                        #'log_lik',
#'                        'omean',
#'                        'zmean'), 
#'               control = list(max_treedepth = 10, adapt_delta=0.8),
#'               verbose = T)
#' 
#' save(fit_CCAR,file='ARfit_raw_cca.Rdata')
#' 
#' 
#' 
#' stan_trace(fit_CCAR, 
#'           pars=c('omean','sigma','grid_sig','beta','lp__','scale','zmean','island_sig','group_sig'), 
#'           include=T)
#' 
#' stan_plot(fit_CCAR, 
#'           pars=c('beta'), 
#'           include=T, 
#'           show_density = TRUE, 
#'           fill_color = "lightblue") + 
#'   theme_cowplot() + 
#'   geom_vline(aes(xintercept=0), linetype=2, alpha=0.4)+
#'   xlab('Effect size') + 
#'   ylab('')
#' 

## CORAL

VAR <- 'CORAL'
model_data <- get_vars(all_data,VAR)

fit_CAR <- stan('model_QR_grid_beta_raw_CAR.stan',
                data = model_data,
                iter = 500,
                chains = 1,
                thin = 1,
                warmup = 300,
                #init = inits,
                pars = c('island_sig',
                  'group_sig',
                  'grid_sig',
                  'rho',
                  'resid_y',
                  'beta',
                  'gscale',
                  #'iscale',
                  'scale',
                  #'log_lik',
                  'omean',
                  'zmean'),
                control = list(max_treedepth = 10, adapt_delta=0.8),
                verbose = T)

save(fit_CAR,file='ARfit_raw_coral.Rdata')
# 
# 
# stan_trace(fit_CAR,
#            pars=c('omean','sigma','beta','lp__','scale','zmean'),
#            include=T)
# 
# 
# ry <- get_posterior_mean(fit_CAR,pars='resid_y')[,3]
# rys <- data.frame(g = model_data$GRID[model_data$ii_notzero],ry,i=model_data$ISLAND[model_data$ii_notzero])
# ggplot(rys) + geom_line(aes(x=g,y=ry)) + facet_wrap(~as.factor(i),scales='free')
# 
# require(purrr)
# acc <- rys %>% split(.$i) %>% map(~acf(.$ry,plot=F,type='cov')) %>% map(~with(.,data.frame(lag,acf))) %>% bind_rows(.id='Island')
# ggplot(acc) + geom_line(aes(x=lag,y=acf)) + facet_wrap(~Island,scales='free')
# 
# 
# 
# stan_plot(fit_C,
#           pars=c('beta'),
#           include=T,
#           show_density = TRUE,
#           fill_color = "lightblue") +
#   theme_cowplot() +
#   geom_vline(aes(xintercept=0), linetype=2, alpha=0.4)+
#   xlab('Effect size') +
#   ylab('')
#' 
#' 
#' 
