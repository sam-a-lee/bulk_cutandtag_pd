sva_custom <- function(dat = dat, mod = mod, mod0 = mod0, n.sv = NULL){
  mono <- function(lfdr){
    .Call("monotone", as.numeric(lfdr), PACKAGE="sva")
  }
  
  edge.lfdr <- function(p, trunc=TRUE, monotone=TRUE, transf=c("probit", "logit"), adj=1.5, eps=10^-8, lambda=0.8, ...) {
    pi0 <- mean(p >= lambda)/(1 - lambda)
    pi0 <- min(pi0, 1)
    n <- length(p)
    transf <- match.arg(transf)
    if(transf=="probit") {
      p <- pmax(p, eps)
      p <- pmin(p, 1-eps)
      x <- qnorm(p)
      myd <- density(x, adjust=adj)
      mys <- smooth.spline(x=myd$x, y=myd$y)
      y <- predict(mys, x)$y
      lfdr <- pi0*dnorm(x)/y
    }
    if(transf=="logit") {
      x <- log((p+eps)/(1-p+eps))
      myd <- density(x, adjust=adj)
      mys <- smooth.spline(x=myd$x, y=myd$y)
      y <- predict(mys, x)$y
      dx <- exp(x) / (1+exp(x))^2
      lfdr <- pi0 * dx/y
    }
    if(trunc) {
      lfdr[lfdr > 1] <- 1
    }
    if(monotone) {	
      lfdr <- lfdr[order(p)]
      lfdr <- mono(lfdr)
      lfdr <- lfdr[rank(p)]
    }
    return(lfdr)
  }
  
  
  ##################################
  # numsv function internal to sva #
  ##################################
  
  #determine the number of required SVs
  #this is numsv function code pulled out
  
  #this is what is fed into this function from the sva function
  #dat = as.matrix(wateRmelon::Beta2M(cb_betas_ordered_cc))
  dat = dat
  mod = mod
  
  #method="be" # this is the default and what is called by numsv
  #vfilter=NULL # we dont apply a filter and its not used for be method
  B=20
  set.seed(1234)
  
  # initialize some stuff
  # get the residual matrix (the remaining variance)
  n <- ncol(dat)
  m <- nrow(dat)
  H <- mod %*% solve(t(mod) %*% mod) %*% t(mod)
  res <- dat - t(H %*% t(dat))
  
  # run svd and get some eigenvalues
  uu <- svd(res)
  # get the num of "sig" eigenvalues
  ndf <- min(m, n) - ceiling(sum(diag(H)))
  # get percent variance explained by each eigenvalue
  dstat <- uu$d[1:ndf]^2/sum(uu$d[1:ndf]^2)
  # initialize null distribution matrix
  dstat0 <- matrix(0, nrow = B, ncol = ndf)
  
  # set seed here since we are calling sample
  # create the null distribution
  for (i in 1:B) {
    # this part of the code removes all relationships in the 
    res0 <- t(apply(res, 1, sample, replace = FALSE))
    res0 <- res0 - t(H %*% t(res0))
    # calculate svs on this
    uu0 <- svd(res0)
    # bind their explained variance in null matrix
    dstat0[i, ] <- uu0$d[1:ndf]^2/sum(uu0$d[1:ndf]^2)
    print(i) # print out the iteration we are on
  }
  
  # determine the number of svs required here based on pvalues
  psv <- rep(1, n)
  for (i in 1:ndf) {
    psv[i] <- mean(dstat0[, i] >= dstat[i])
  }
  for (i in 2:ndf) {
    psv[i] <- max(psv[(i - 1)], psv[i])
  }
  
  # this is the variable holding the number of SVs required 
  # number is each to all the svs with a pval less than/equal to 0.1
  nsv <- sum(psv <= 0.1)
  
  
  ####################################
  # code from the main SVA function  #
  ####################################
  n.sv = n.sv
  n.sv = ifelse(!is.null(n.sv), n.sv, nsv)
  #n.sv = nsv # from code above (run internally when sva called)
  B = 5 # number of iterations
  
  # dat =  betas # already initialized
  # mod = mod  # already initialized
  # mod0 = mod0  # already initialized 
  # controls = NULL # we dont include controls - only for supervised method
  # method = "irw" # this isnt actually used but would be called by sva
  # vfilter = NULL # we dont set this when running sva
  # numSVmethod = "be" # this is the default called by sva (code above)
  
  
  ###################################
  # irwsva function internal to sva #
  ###################################
  
  # this code is really want does all the work in sva
  # we use the irw method so this is the code for that method
  # there are also supervised and two step methods (not used by us)
  
  # irw function takes input from sva, so we have already defined below vars
  # dat = dat 
  # mod = mod
  # mod0 = mod0
  # n.sv = n.sv
  # B = B
  
  # intialize variables again
  n <- ncol(dat)
  m <- nrow(dat)
  Id <- diag(n)
  resid <- dat %*% (Id - mod %*% solve(t(mod) %*% mod) %*% t(mod))
  uu <- eigen(t(resid) %*% resid)
  vv <- uu$vectors
  ndf <- n - dim(mod)[2]
  pprob <- rep(1, m)
  one <- rep(1, n)
  Id <- diag(n)
  df1 <- dim(mod)[2] + n.sv
  df0 <- dim(mod0)[2] + n.sv
  rm(resid)
  
  cat(paste("Iteration (out of", B, "):"))
  for (i in 1:B) {
    mod.b <- cbind(mod, uu$vectors[, 1:n.sv])
    mod0.b <- cbind(mod0, uu$vectors[, 1:n.sv])
    ptmp <- f.pvalue(dat, mod.b, mod0.b)
    pprob.b <- (1 - edge.lfdr(ptmp))
    mod.gam <- cbind(mod0, uu$vectors[, 1:n.sv])
    mod0.gam <- cbind(mod0)
    ptmp <- f.pvalue(dat, mod.gam, mod0.gam)
    pprob.gam <- (1 - edge.lfdr(ptmp))
    pprob <- pprob.gam * (1 - pprob.b)
    dats <- dat * pprob
    dats <- dats - rowMeans(dats)
    uu <- eigen(t(dats) %*% dats)
    cat(paste(i, " "))
  }
  
  # this is what is normally returned by sva
  # problem is that it doesnt include the diagonal matrix we want
  sv_cb = svd(dats)$v[, 1:n.sv, drop = FALSE] 
  svobj <- list(sv = sv_cb, pprob.gam = pprob.gam, pprob.b = pprob.b, 
                n.sv = n.sv)
  
  # to get all info (including diagonal)
  svobj_all <- svd(dats)
  
  sva_list <- list("sv" = sv_cb, "svobj" = svobj, "svobj_all" = svobj_all)
  return(sva_list)
  
}

#####################
# heat scree of SVs #
#####################

# from elodie portales-casamar

# # heat scree prepartion
# vars <- svobj_all[['d']]^2
# vars <- vars/sum(vars)
# names(vars) <- paste0('SV', 1:length(vars))
# xp <- 1:svobj$n.sv
# varsdf <- as.data.frame(vars)
# varsdf$svnum <- rep(1:length(vars))
# 
# # heat scree ggplot
# no2_sv_scree <- ggplot(varsdf[1:20,], aes(x=svnum, y=vars*100)) +
#   geom_col() +
#   theme_bw() +
#   labs(x="Surrogate variables", y = "Variance (%)")
# 
# ggsave(plot = no2_sv_scree, 
#        filename = here("CANDLE_figures", 
#                        "CANDLE_ewas", 
#                        "CANDLE_cord_blood_air_pollution", 
#                        "CANDLE_cord_blood_prenatal_no2",
#                        "2023-05-15_cord_blood_no2_sv_scree_plot.png"),
#        height = 5, width = 16, units = "cm")
