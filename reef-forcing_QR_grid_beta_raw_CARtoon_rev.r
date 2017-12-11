require(readr)
require(dplyr)
require(cowplot)
require(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 3)


# get_vars <- function(data,var){
#   #require(INLA)
#   adata <- data %>% filter(!is.na(TED_SUM),
#                            !is.na(SURF_SUM)) %>%
#     mutate_(mm = var) %>%
#     mutate(mm=ifelse(N_OBS_CELL>0,mm/100,NA),
#            #DPOP      = round(scale(distant_hum_pop)[,1],2),
#            TED_SUMS   = round(scale(TED_SUM)[,1],1),
#            SURF_SUMS  = round(scale(SURF_SUM)[,1],1),
#            TED_SUM  = TED_SUMS*attr(scale(TED_SUM),'scaled:scale')+attr(scale(TED_SUM),'scaled:center'),
#            SURF_SUM  = SURF_SUMS*attr(scale(SURF_SUM),'scaled:scale')+attr(scale(SURF_SUM),'scaled:center'),
#            #TED_INT   = round(scale(TED_SUM*distant_hum_pop)[,1],1),
#            PHYS_INT  = round(scale(SURF_SUM*TED_SUM)[,1],1),
#            UCOVS = paste(TED_SUMS,SURF_SUMS,PHYS_INT))
# 
# 
# 
#   pmax <- 0
#   island <- as.numeric(as.factor(adata$ISLAND))
#   islands <- unique(island)
# 
#   for (i in 1:length(islands)) {
#     cat(i,'\n')
#     ids  <- adata$GRID_ID[island==islands[i]]
#     cvs  <- adata$UCOVS[island==islands[i]]
# 
#     li <- sum(island==islands[i])
#     grd <- vector(,li)
#     grd[1] <- 1
#     for (id in 2:li){
#       if(ids[id] == ids[id-1] | (ids[id] == ids[id-1]+1 & cvs[id] == cvs[id-1])) {
#         grd[id] <- grd[id-1]
#         } else if(ids[id] == ids[id-1]+1 & cvs[id] != cvs[id-1]) {
#           grd[id] <- grd[id-1]+1
#         } else {
#           grd[id] <- grd[id-1]+2
#         }
#     }
#     adata$GRID_ID[island==islands[i]] <- grd
#   }
# 
#   nbbs <- list()
#   idds <- list()
# 
#   D_sparses <- list()
#   lambdas <- list()
#   W_n <- vector(,length(islands))
#   grids <- vector(,length(islands))
# 
#   for (i in 1:length(islands)) {
#     cat(i,'\n')
#     ids  <- adata$GRID_ID[island==islands[i]]
#     nbbs[[i]] <- lapply(1:length(unique(ids)), function(i) {
#       id <- unique(ids)[i]
#       if(id<max(ids)) {
#         if(!any(ids==(id+1))) NULL
#         else c(i,i+1)
#       }
# 
#     })
# 
#     idds[[i]] <- match(ids,unique(ids))+pmax
#     pmax <- pmax + length(unique(ids))
#     #browser()
#     nbs <- do.call(rbind,nbbs[[i]])
#     W = as.matrix(Matrix::sparseMatrix(i=nbs[,1],j=nbs[,2],x=1,symmetric=TRUE,check = F))
#     D_sparses[[i]] <- rowSums(W)
#     # get eigenvalues of D^(-.5) * W * D^(-.5) for determinant computations
#     invsqrtD <- diag(1 / sqrt(D_sparses[[i]]))
#     quadformDAD <- drop(crossprod(crossprod(W, invsqrtD), invsqrtD))
#     lambdas[[i]] <- eigen(quadformDAD)$values
#     grids[i] <- length(unique(ids))
#     W_n[i] <- nrow(nbs)
#   }
# 
#   D_sparse <- do.call(c,D_sparses)
#   lambda <- do.call(c,lambdas)
#   GRID <- do.call(c,idds)
#   N_grid <- length(unique(GRID))
#   nbs <- do.call(rbind,do.call(c,nbbs))
#   N_edges <- nrow(nbs)
# 
#   require(Matrix)
# 
#   COVS <- with(adata, as.matrix(
#     data.frame(TED_SUMS,
#                SURF_SUMS,
#                PHYS_INT
#       )))
# 
#   #browser()
# 
#   with(adata,list(
#     W_sparse  = nbs,
#     D_sparse  = D_sparse,
#     lambda    = lambda,
#     grids     = grids,
#     W_n       = W_n,
#     W_ns      = N_edges,
#     ZEROONE   = ifelse(mm[!is.na(mm)]>0,1,0),
#     N_notzero = sum(mm>0,na.rm=T),
#     N_zeroone = sum(!is.na(mm)),
#     N_groups  = length(unique(group)),
#     N_islands = length(unique(ISLAND)),
#     N         = nrow(COVS),
#     K         = ncol(COVS),
#     ii_notzero= which(mm>0),
#     ii_zeroone= which(!is.na(mm)),
#     ii_map    = which(mm[which(!is.na(mm))]>0),
#     IMEAN     = mm[which(mm>0)],
#     COVS      = COVS,
#     ISLAND    = as.numeric(as.factor(ISLAND)),
#     GROUP     = as.numeric(as.factor(group)),
#     GRID      = GRID,
#     N_grid    = length(unique(GRID))
#   ))
# 
# }
# #
# #
# cdata <- read_csv('./model/data/ENTIRE_CLIMATOLOGY_ALL_ISLANDS.csv')
# rdata <- read_csv('./model/data/RAW_BENTHIC_ALL_ISLANDS.csv')
# #
# all_data <- inner_join(rdata,cdata)
# #
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
                               chains = 4,
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
                                        'log_lik',
                                        'omean',
                                        'zmean'),
                               control = list(max_treedepth = 13, adapt_delta=0.8),
                               verbose = T)

save(fit_MARtoon_rev,file='CARtoonfit_raw_macroalgae_hhCAR.Rdata')
# 
# gg <- stan_trace(fit_MARtoon_rev,
#           pars=c('omean','rho_prec','rho[10]','grid_sig[10]','beta[1,1]',
#                  't_var','t_mean','lp__','scale','zmean','island_sig','group_sig'),
#           #pars=c('rho_prec'),
#           include=T)
# ggsave(gg,file='traces.pdf')

# 
# beats <- extract(fit_MARtoon_rev,'beta')$beta
# 
# dimnames(beats)[[3]] <- dimnames(model_data$COVS)[[2]]
# dimnames(beats)[[2]] <- unique(all_data$ISLAND)
# 
# beat <- tbl_df(reshape2::melt(beats))
# colnames(beat)[2:3] <- c('Island','Covariate')
# 
# gg <- ggplot(beat) +
#   geom_violin(aes(x=Covariate,y=value)) +
#   geom_hline(aes(yintercept=0)) +
#   facet_wrap(~Island) + 
#   theme_cowplot() +
#   geom_vline(aes(xintercept=0), linetype=2, alpha=0.4)+
#   xlab('Effect size') +
#   ylab('')
# # 
# ggsave(gg,file='betas.pdf')


# # 
# logit <- function(p) log(p/(1-p))
# ry <- get_posterior_mean(fit_MARtoon_rev,pars='resid_y')[,1]
# rys <- data.frame(g = model_data$GRID[model_data$ii_notzero],ry,i=model_data$ISLAND[model_data$ii_notzero],imean=logit(model_data$IMEAN))
# #ggplot(rys) + geom_line(aes(x=g,y=ry)) + facet_wrap(~as.factor(i),scales='free')
# 
# require(purrr)
# acc <- rys %>% arrange(i,g) %>% split(.$i) %>% map(~acf(.$ry,plot=F,type='cov',lag.max = 300)) %>% map(~with(.,data.frame(lag,acf))) %>% bind_rows(.id='Island')
# gg <- ggplot(acc) + geom_line(aes(x=lag,y=acf)) + facet_wrap(~Island,scales='free')
# ggsave(gg,file='spatial.pdf')

## CCA

#VAR <- 'CCA'
#model_data <- get_vars(all_data,VAR)
#save(model_data,file='CAR_cca.Rdata')

#' 
#' load('CAR_cca.Rdata')
#' 
#' fit_CCARtoon <- stan('model_QR_grid_beta_raw_CARtoon_rev.stan',
#'                      data = model_data,
#'                      iter = 1000,
#'                      chains = 4,
#'                      thin = 1,
#'                      warmup = 500,
#'                      init = function() list(rho_mean=0.1),
#'                      pars = c('island_sig',
#'                               'group_sig',
#'                               'grid_mean',
#'                               'grid_sig',
#'                               'grid_var',
#'                               'rho',
#'                               'rho_mean',
#'                               'rho_prec',
#'                               'resid_y',
#'                               'beta',
#'                               'gscale',
#'                               #'iscale',
#'                               'scale',
#'                               'log_lik',
#'                               'omean',
#'                               'zmean'),
#'                      control = list(max_treedepth = 12, adapt_delta=0.8),
#'                      verbose = T)
#' 
#' save(fit_CCARtoon,file='CARtoonfit_raw_cca.Rdata')
#' 
#' 
#' 
#' stan_trace(fit_CCARtoon,
#'           pars=c('omean','beta','lp__','scale','zmean','island_sig','group_sig','rho_mean'),
#'           include=T)
#' 
#' stan_plot(fit_CCA,
#'           pars=c('beta'),
#'           include=T,
#'           show_density = TRUE,
#'           fill_color = "lightblue") +
#'   theme_cowplot() +
#'   geom_vline(aes(xintercept=0), linetype=2, alpha=0.4)+
#'   xlab('Effect size') +
#'   ylab('')
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
#' load('CAR_coral_rev.Rdata')
#' fit_CARtoon_rev <- stan('model_QR_grid_beta_raw_CARtoon_rev.stan',
#'                 data = model_data,
#'                 iter = 1000,
#'                 chains = 4,
#'                 thin = 1,
#'                 warmup = 500,
#'                 init = function() list(rho_mean=0.1),
#'                 pars = c('island_sig',
#'                   'group_sig',
#'                   'grid_mean',
#'                   'grid_sig',
#'                   'grid_var',
#'                   'rho',
#'                   'rho_mean',
#'                   'rho_prec',
#'                   'resid_y',
#'                   'beta',
#'                   'gscale',
#'                   #'iscale',
#'                   'scale',
#'                   'log_lik',
#'                   'omean',
#'                   'zmean'),
#'                 control = list(max_treedepth = 12, adapt_delta=0.8),
#'                 verbose = T)
#' 
#' save(fit_CARtoon_rev,file='CARtoonfit_raw_coral_rev.Rdata')
# 
# # 
# stan_trace(fit_CARtoons,
#             pars=c('group_sig','island_sig','omean','beta','lp__','scale','zmean','rho_mean','grid_sig'),
#             include=T)
# 
# 
# ry <- get_posterior_mean(fit_CARtoon,pars='resid_y')[,1]
# rys <- data.frame(g = model_data$GRID[model_data$ii_notzero],ry,i=model_data$ISLAND[model_data$ii_notzero])
# ggplot(rys) + geom_line(aes(x=g,y=ry)) + facet_wrap(~as.factor(i),scales='free')
# 
# require(purrr)
# acc <- rys %>% split(.$i) %>% map(~acf(.$ry,plot=F,type='cov')) %>% map(~with(.,data.frame(lag,acf))) %>% bind_rows(.id='Island')
# ggplot(acc) + geom_line(aes(x=lag,y=acf)) + facet_wrap(~Island,scales='free')

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
