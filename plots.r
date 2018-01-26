plot_betas <- function(post, model_data) {
  beats <- extract(post,'beta')$beta
  
  R_ast = qr.R(qr(model_data$COVS)) / sqrt(model_data$N - 1.0);
  R_ast_inverse = solve(R_ast);
  teats <- t(apply(extract(post,'t_mean')$t_mean, 1, function(x) R_ast_inverse %*% x))
  
  dimnames(teats)[[2]] <- c('TED','SURF','INT')
  teat <- tbl_df(reshape2::melt(teats))
  colnames(teat)[c(2,3)] <- c('Covariate','Mean')
  
  dimnames(beats)[[3]] <- c('TED','SURF','INT')
  dimnames(beats)[[2]] <- unique(all_data$ISLAND)
  
  beat <- tbl_df(reshape2::melt(beats))
  colnames(beat)[2:3] <- c('Island','Covariate')
  
  bt <- inner_join(beat,teat) %>% mutate(actual=value+Mean)
  bt$Island <- factor(bt$Island,levels = ISLANDS)
  
  
  gg <- ggplot(bt) +
    geom_violin(aes(x=Covariate,y=actual)) +
    geom_hline(aes(yintercept=0)) +
    facet_wrap(~Island, scales='free',ncol = 5) +
    geom_violin(aes(x=Covariate,y=Mean),alpha=0.5,fill='skyblue',colour=NA) +
    theme_cowplot() +
    
    xlab('Effect size') +
    ylab('')
  gg
}

plot_resids <- function(posts,model_data,ac=F){
  logit <- function(p) log(p/(1-p))
  ry <- get_posterior_mean(posts,pars='resid_y')[,1]
  rys <- data.frame(g = model_data$GRID[model_data$ii_notzero],
                    ry,
                    i=unique(all_data$ISLAND)[model_data$ISLAND[model_data$ii_notzero]],imean=logit(model_data$IMEAN))
  if(ac == F) {
    
    rys$i <- factor(rys$i,levels=ISLANDS)
    
    ggplot(rys) + 
      geom_line(aes(x=g,y=ry)) + 
      facet_wrap(~as.factor(i),scales='free')+ 
      theme_cowplot()
    
  } else {
    require(purrr)
    acc <- rys %>% 
      arrange(i,g) %>% 
      split(.$i) %>% 
      map(~acf(.$ry,plot=F,type='cov',lag.max = 100)) %>% 
      map(~with(.,data.frame(lag,acf))) %>% 
      bind_rows(.id='Island')
    
    acc$Island <- factor(acc$Island,levels=ISLANDS)
    ggplot(acc) + 
      geom_line(aes(x=lag,y=acf)) + 
      facet_wrap(~Island,scales='free') + 
      theme_cowplot()
  }
  
}

HDIofMCMC = function( sampleVec , credMass=0.95 ) {
  # Computes highest density interval from a sample of representative values,
  #   estimated as shortest credible interval.
  # Arguments:
  #   sampleVec
  #     is a vector of representative values from a probability distribution.
  #   credMass
  #     is a scalar between 0 and 1, indicating the mass within the credible
  #     interval that is to be estimated.
  # Value:
  #   HDIlim is a vector containing the limits of the HDI
  sortedPts = sort( sampleVec )
  ciIdxInc = ceiling( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( HDImin , HDImax )
  return( HDIlim )
}


plot_fit <- function(posts,model_data){
  
  logit <- function(p) log(p/(1-p))
  invlogit <- function(z) 1/(1+exp(-z))
  py <- invlogit(logit(model_data$IMEAN)-extract(posts,pars='resid_y')$resid_y)
  ppy <- t(apply(py,2,HDIofMCMC))
  colnames(ppy) <- c('l','h')
  ry <- get_posterior_mean(posts,pars='resid_y')[,1]
  rys <- data.frame(g = model_data$GRID[model_data$ii_notzero],
                    pred = invlogit(logit(model_data$IMEAN)-ry),
                    data = model_data$IMEAN,
                    ppy,
                    i=unique(all_data$ISLAND)[model_data$ISLAND[model_data$ii_notzero]],
                    imean=logit(model_data$IMEAN))
  
  rys$i <- factor(rys$i,levels=ISLANDS)
  ggplot(rys) + 
    geom_ribbon(aes(x=g,y=pred,ymin=l,ymax=h),alpha=0.4,col=NA,fill='#AAFF88') + 
    geom_point(aes(x=g,y=data),col='skyblue') + 
    geom_line(aes(x=g,y=pred),col='#AAFF88',size=1.5) + 
    facet_wrap(~as.factor(i),scales='free')+
    xlab('Grid location')+
    ylab('Percent cover')+
    theme_cowplot()+ 
    theme(axis.text.x = element_blank())
}

get_p_beta <- function(post,model_data){
  R_ast = qr.R(qr(model_data$COVS)) / sqrt(model_data$N - 1.0);
  R_ast_inverse = solve(R_ast);
  teats <- t(apply(extract(post,'t_mean')$t_mean, 1, function(x) R_ast_inverse %*% x))
  
  colnames(teats) <- c('TED','SURF','INT')
  apply(teats,2,function(x) mean(x<0))
}