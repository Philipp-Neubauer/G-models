require(INLA)
require(dplyr)
require(readr)

# 
# 
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
#     grd[1] <- pmax+1
#     for (id in 2:li){
#       if(ids[id] == ids[id-1] | (ids[id] == ids[id-1]+1 & cvs[id] == cvs[id-1])) {
#         grd[id] <- grd[id-1]
#       } else if(ids[id] == ids[id-1]+1 & cvs[id] != cvs[id-1]) {
#         grd[id] <- grd[id-1]+1
#       } else {
#         grd[id] <- grd[id-1]+2
#       }
#     }
#     pmax=max(grd)
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
#   pmax=0
#   for (i in 1:length(islands)) {
#     cat(i,'\n')
#     ids  <- adata$GRID_ID[island==islands[i]]
#     nbbs[[i]] <- lapply(1:length(unique(ids)), function(i) {
#       id <- unique(ids)[i]
#       if(id<max(ids)) {
#         if(!any(ids==(id+1))) NULL
#         else c(i,i+1)+pmax
#       }
#       else {c(i,1)+pmax}
#       
#     })
#     
#     idds[[i]] <- match(ids,unique(ids))+pmax
#     pmax <- pmax + length(unique(ids))
#     
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
#     )))
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
#     ZEROONE   = ifelse(mm>0,1,0),
#     N_notzero = sum(mm>0,na.rm=T),
#     N_zeroone = sum(!is.na(mm)),
#     N_groups  = length(unique(group)),
#     N_islands = length(unique(ISLAND)),
#     N         = nrow(COVS),
#     K         = ncol(COVS),
#     ii_notzero= which(mm>0),
#     ii_zeroone= which(!is.na(mm)),
#     ii_map    = which(mm[which(!is.na(mm))]>0),
#     IMEAN     = ifelse(mm>0,mm,NA),
#     COVS      = COVS,
#     ISLAND    = as.numeric(as.factor(ISLAND)),
#     GROUP     = as.numeric(as.factor(group)),
#     GRID      = GRID,
#     N_grid    = length(unique(GRID))
#   ))
#   
# }
# 
# 
# make_INLA_data <- function(tags){
#   
#   z <- tags$ZEROONE                   
#   y <- tags$IMEAN
#   
#   COVz <- tags$COVS 
#   colnames(COVz) <- paste0(colnames(COVz),'.z')
#   stk.z <- inla.stack(tag='est.z', 
#                       data=list(y=cbind(z, NA)), ### z at first column of y 
#                       A=list(1,1,1,1,1,1,1), 
#                       effects=list( list(grid.z=as.numeric(as.factor(paste(tags$ISLAND,tags$GRID)))), 
#                                     list(Island=tags$ISLAND), 
#                                     list(Island.z=tags$ISLAND), 
#                                     list(Island.z.c1=tags$ISLAND), 
#                                     list(Island.z.c2=tags$ISLAND),  
#                                     list(Group.z=tags$GROUP), 
#                                     list(z.b0=rep(1,length(z)),
#                                          TED_SUMS.z=COVz[,1],
#                                          SURF_SUMS.z=COVz[,2])))
#   
#   COVy <- tags$COVS 
#   colnames(COVy) <- paste0(colnames(COVy),'.y')
#   
#   nbs <- tags$W_sparse
#   g <- INLA:::inla.read.graph(as.matrix(Matrix::sparseMatrix(i=nbs[,1],j=nbs[,2],x=1,symmetric=TRUE,check = F)))
#   
#   stk.y <- inla.stack(tag='est.y', 
#                       data=list(y=cbind(NA, y)), ### at second column 
#                       A=list(1,1,1,1,1,1,1), 
#                       effects=list( list(grid.y=as.numeric(as.factor(paste(tags$ISLAND,tags$GRID)))), 
#                                     list(Island=tags$ISLAND), 
#                                     list(Island.y=tags$ISLAND),
#                                     list(Island.y.c1=tags$ISLAND), 
#                                     list(Island.y.c2=tags$ISLAND), 
#                                     list(Group.y=tags$GROUP), 
#                                     list(y.b0=rep(1,length(y)),
#                                          TED_SUMS.y=COVy[,1],
#                                          SURF_SUMS.y=COVy[,2])))
#   
#   
#   inla.stack(stk.z, stk.y) 
# }
# 

form <- y ~ 0 + z.b0 + y.b0 + 
  f(Island.y.c1,TED_SUMS.y,model='iid',
    hyper=list(theta=list(prior="loggamma",param=c(1,0.001)))) + 
  f(Island.z.c1,copy='Island.y.c1') + 
  
  f(Island.y.c2,SURF_SUMS.y,model='iid',
    hyper=list(theta=list(prior="loggamma",param=c(1,0.001)))) +  
  f(Island.z.c2,copy='Island.y.c2') +  
  
  f(Island.y,model='iid',
    hyper=list(theta=list(prior="loggamma",param=c(1,0.001)))) +  
  f(Island.z,copy='Island.y') +  
  
  f(grid.y, model='bym', graph=g,scale.model =T,control.group=list(model="exchangeable")) +
  f(grid.z, copy='grid.y')
# 
# cdata <- read_csv('./model/data/ENTIRE_CLIMATOLOGY_ALL_ISLANDS.csv')
# rdata <- read_csv('./model/data/RAW_BENTHIC_ALL_ISLANDS.csv')
# 
# all_data <- inner_join(rdata,cdata)
# 
# tags <- get_vars(all_data,'MACRO')
# 
# stk <- make_INLA_data(tags)
# save(stk,file='stk_macro.Rdata')

load('stk_macro.Rdata')
res.yz <- inla(form, 
               family=c('binomial', 'beta'), 
               data=inla.stack.data(stk),
               control.predictor=list(quantiles = c(0.025,0.25,0.5,0.75,0.975),                                               A=inla.stack.A(stk))
) 

save(res.yz,file='INLA_resyz_MACRO.Rdata')        

# 
# tags <- get_vars(all_data,'CCA')
# 
# stk <- make_INLA_data(tags)
# 
# save(stk,file='stk_CCA.Rdata')

load('stk_CCA.Rdata')

res.yz <- inla(form, 
               family=c('binomial', 'beta'), 
               data=inla.stack.data(stk),
               control.predictor=list(quantiles = c(0.025,0.25,0.5,0.75,0.975),                                               A=inla.stack.A(stk))
) 

save(res.yz,file='INLA_resyz_CCA.Rdata')          
# 
# 
# tags <- get_vars(all_data,'CORAL')
# 
# stk <- make_INLA_data(tags)
# 
# save(stk,file='stk_coral.Rdata')

load('stk_coral.Rdata')

res.yz <- inla(form, 
               family=c('binomial', 'beta'), 
               data=inla.stack.data(stk),
               control.predictor=list(quantiles = c(0.025,0.25,0.5,0.75,0.975),                                               A=inla.stack.A(stk))
) 

save(res.yz,file='INLA_resyz_CORAL.Rdata')          
