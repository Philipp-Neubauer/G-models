

get_vars <- function(data,var){
  #require(INLA)
  adata <- data %>% filter(!is.na(TED_SUM),
                           !is.na(SURF_SUM)) %>%
    mutate_(mm = var) %>%
    mutate(mm=ifelse(N_OBS_CELL>0,mm/100,NA),
           #DPOP      = round(scale(distant_hum_pop)[,1],2),
           TED_SUMS   = round(scale(TED_SUM)[,1],1),
           SURF_SUMS  = round(scale(SURF_SUM)[,1],1),
           TED_SUM  = TED_SUMS*attr(scale(TED_SUM),'scaled:scale')+attr(scale(TED_SUM),'scaled:center'),
           SURF_SUM  = SURF_SUMS*attr(scale(SURF_SUM),'scaled:scale')+attr(scale(SURF_SUM),'scaled:center'),
           #TED_INT   = round(scale(TED_SUM*distant_hum_pop)[,1],1),
           PHYS_INT  = round(scale(SURF_SUM*TED_SUM)[,1],1),
           UCOVS = paste(TED_SUMS,SURF_SUMS,PHYS_INT))



  pmax <- 0
  island <- as.numeric(as.factor(adata$ISLAND))
  islands <- unique(island)

  for (i in 1:length(islands)) {
    cat(i,'\n')
    ids  <- adata$GRID_ID[island==islands[i]]
    cvs  <- adata$UCOVS[island==islands[i]]

    li <- sum(island==islands[i])
    grd <- vector(,li)
    grd[1] <- 1
    for (id in 2:li){
      if(ids[id] == ids[id-1] | (ids[id] == ids[id-1]+1 & cvs[id] == cvs[id-1])) {
        grd[id] <- grd[id-1]
      } else if(ids[id] == ids[id-1]+1 & cvs[id] != cvs[id-1]) {
        grd[id] <- grd[id-1]+1
      } else {
        grd[id] <- grd[id-1]+2
      }
    }
    adata$GRID_ID[island==islands[i]] <- grd
  }

  nbbs <- list()
  idds <- list()

  D_sparses <- list()
  lambdas <- list()
  W_n <- vector(,length(islands))
  grids <- vector(,length(islands))

  for (i in 1:length(islands)) {
    cat(i,'\n')
    ids  <- adata$GRID_ID[island==islands[i]]
    nbbs[[i]] <- lapply(1:length(unique(ids)), function(i) {
      id <- unique(ids)[i]
      if(id<max(ids)) {
        if(!any(ids==(id+1))) NULL
        else c(i,i+1)
      }

    })

    idds[[i]] <- match(ids,unique(ids))+pmax
    pmax <- pmax + length(unique(ids))
    #browser()
    nbs <- do.call(rbind,nbbs[[i]])
    W = as.matrix(Matrix::sparseMatrix(i=nbs[,1],j=nbs[,2],x=1,symmetric=TRUE,check = F))
    D_sparses[[i]] <- rowSums(W)
    # get eigenvalues of D^(-.5) * W * D^(-.5) for determinant computations
    invsqrtD <- diag(1 / sqrt(D_sparses[[i]]))
    quadformDAD <- drop(crossprod(crossprod(W, invsqrtD), invsqrtD))
    lambdas[[i]] <- eigen(quadformDAD)$values
    grids[i] <- length(unique(ids))
    W_n[i] <- nrow(nbs)
  }

  D_sparse <- do.call(c,D_sparses)
  lambda <- do.call(c,lambdas)
  GRID <- do.call(c,idds)
  N_grid <- length(unique(GRID))
  nbs <- do.call(rbind,do.call(c,nbbs))
  N_edges <- nrow(nbs)

  require(Matrix)

  COVS <- with(adata, as.matrix(
    data.frame(TED_SUMS,
               SURF_SUMS,
               PHYS_INT
    )))

  #browser()

  with(adata,list(
    W_sparse  = nbs,
    D_sparse  = D_sparse,
    lambda    = lambda,
    grids     = grids,
    W_n       = W_n,
    W_ns      = N_edges,
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
