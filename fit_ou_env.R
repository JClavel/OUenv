################################################################################
##                                                                            ##
##                               fit_OU_trend   (v1.1)                        ##
##                                                                            ##
##  Created by Julien Clavel - 22-05-2019                                     ##
##  (julien.clavel@hotmail.fr/ j.clavel@nhm.ac.uk)                            ##
##   require: mvMORPH                                                         ##
##                                                                            ##
################################################################################


fit_OU_trend <- function(tree, data, fun=NULL, startvalues=NULL, echo=TRUE, method="Nelder-Mead", nuisance=FALSE, mserr=0, control=list(),...){
  require(mvMORPH)
  if(is.null(startvalues)) stop("You must provide starting values for your customized function")
  
  # options
  args = list(...)
  if(is.null(args[["upper"]])) upper <- Inf else upper <- args$upper
  if(is.null(args[["lower"]])) lower <- -Inf else lower <- args$lower
  if(is.null(args[["fixedRoot"]])) fixedRoot <- TRUE else fixedRoot <- args$fixedRoot
  
  # param
  nobs <- Ntip(tree)
  C1 <- vcv.phylo(tree)
  W <- matrix(0,ncol=1,nrow=nobs) # to tweak mvmorph function
  times <- diag(C1)
  
  # reorder
  if(!is.null(names(data))) data <- data[tree$tip.label] else warning("data and tips assumed to be aligned")
  
  # Define a function to find the OU expectation
  Expectation_OU <- function(time, theta_0, alpha, fun, par){
    
    # Integral
    fun_ou_expectation <- function(x, stop_time, par, theta_0) { alpha*exp(alpha*(x-stop_time))*fun(x, par, theta_0) }
   
    # Numerical integration
    expectation_vector <- sapply(time, function(t){
       int <- try(integrate(f=fun_ou_expectation, lower=0, upper=t, stop_time=t, par=par, theta_0=theta_0, subdivisions=200L, rel.tol = .Machine$double.eps^0.05), silent=TRUE)
      if(inherits(int ,'try-error')){
        warning("An error occured during numerical integration. The integral is probably divergent or your function is maybe undefined for some values")
        integ <- NA_real_
      } else {
        integ <- int$value
      }
     #theta_0*exp(-alpha*t) + integ
      # TO CHECK but I think it's better to use as the start of the process X(0) = X(0)+beta*T(0) to take into account possible offset if the climatic curve at t=0 is not equal to 0
      exp(-alpha*t)*fun(0, par=par, theta_0) + integ
    })
    return(expectation_vector)
  }
  
  # log-likelihood function
  llik <- function(par, phy, data, fun_ou){
    
    # parameters
    sigma2 = exp(par[1])
    alpha = exp(par[2])
    theta_0 = par[3]
    if(nuisance) extra_par = par[5:length(par)] else extra_par = par[4:length(par)]
    
    # Variance-Covariance
    if(fixedRoot){
        V<-.Call("mvmorph_covar_ou_fixed", A=C1, alpha=alpha, sigma=sigma2, PACKAGE="mvMORPH")
     
        # add ME and nuisance variance
        if(nuisance){
          nuisanceME = exp(par[4])
          diag(V) <- diag(V) + mserr^2 + nuisanceME
        }else if(!is.null(mserr)){
          diag(V) <- diag(V) + mserr^2
        }
        
    }else{
        V<-.Call("mvmorph_covar_ou_random", A=C1, alpha=alpha, sigma=sigma2, PACKAGE="mvMORPH")
        
        # add ME and nuisance variance
        if(nuisance){
          nuisanceME = exp(par[4])
          diag(V) <- diag(V) + mserr^2 + nuisanceME
        }else if(!is.null(mserr)){
          diag(V) <- diag(V) + mserr^2
        }
        
    }
    
    # Expectation
    E <- Expectation_OU(times, theta_0=theta_0, alpha=alpha, fun=fun, par=extra_par)
    # check if there were an error in the integration
    if(any(is.na(E))) return(list(logl=-1e7))
    
    # ll computation
    residuals <- (data - E)
    llik <- try(mvLL(V, residuals, method="rpf", param=list(D=W, estim=FALSE, mu=1)), silent = TRUE)
    
    # Try catch if there is an error
    if(inherits(llik ,'try-error')){
        warning("An error occured during the likelihood estimation. An infinite value was returned to continue the search")
        llik <- list()
        llik$logl <- 1e6
    }
    
    return(llik)
  }
  
  # startValues
  #sig_est <- mvLL(tree, data, param=list(estim=TRUE))$sigma 
  if(nuisance){
    if(echo==TRUE) message("Finding best starting values...")
    #sig_start <- log(sig_est*c(0.2,0.5,0.8,0.9))
    #nuis_start <- log(sig_est*c(0.01, 0.05, 0.1))
    alpha_start <- log(log(2)/(max(node.depth.edgelength(tree))/c(0.1,0.5,1.5,3,8)))
    sig_start <- log(Stationary_to_scatter(exp(alpha_start), var(data))) # condition the starting values on alpha?
    nuis_start <- log(exp(mean(sig_start))*c(0.01, 0.05, 0.1))
    theta_start <- mean(data)*c(0.5,0.8,1)
    start_val_ou <- list(sigma2=sig_start, alpha=alpha_start, theta_0=theta_start, mserr=nuis_start)
    
    mat_start <- expand.grid(c(start_val_ou ,startvalues))
    start_val <- mat_start[which.min(apply(mat_start, 1, function(par) -llik(par, phy=tree, data=data, fun_ou=fun)$logl)),]#options
    names(start_val) = NULL
    
  }else{
    if(echo==TRUE) message("Finding best starting values...")
    # initialize the parameter search
    #sig_start <- log(as.numeric(sig_est)*c(0.2,0.5,0.8,0.9))
    alpha_start <- log(log(2)/(max(node.depth.edgelength(tree))/c(0.1,0.5,1.5,3,8)))
    sig_start <- log(Stationary_to_scatter(exp(alpha_start), var(data))) # condition the starting values on alpha?
    theta_start <- mean(data)*c(0.5,0.8,1)
    start_val_ou <- list(sigma2=sig_start, alpha=alpha_start, theta_0=theta_start)
    
    mat_start <- expand.grid(c(start_val_ou ,startvalues))
    start_val <- mat_start[which.min(apply(mat_start, 1, function(par) -llik(par, phy=tree, data=data, fun_ou=fun)$logl)),]#options
    names(start_val) = NULL
  } 
  
  # optimization
  if(echo==TRUE) message("Start optimization (and numerical integration). Please wait...")
  if(method=="spg"){
      require(BB)
  estimModel <- spg(unlist(start_val),
      fn = function(par) -llik(par, phy=tree, data=data, fun_ou=fun)$logl,
      method=3,
      upper=upper,
      lower=lower,
      control=control)
  }else{
  estimModel <- optim(start_val,
                      fn = function(par) -llik(par, phy=tree, data=data, fun_ou=fun)$logl,
                      method=method,
                      upper=upper,
                      lower=lower,
                      hessian=TRUE,
                      control=control)
  }
  
  # Done retrieve param
  if(echo==TRUE) message("Done. retrieve parameters and results...")
  param = estimModel$par
  LL = -estimModel$value
  nparam = length(param)
  
  # parameters
  sigma2 = exp(param[1])
  alpha = exp(param[2])
  theta_0 = param[3]
  if(nuisance) estimated_ME = exp(param[4]) else estimated_ME = NULL
  if(nuisance) custom = param[5:length(param)] else custom = param[4:length(param)]
  
  if(nuisance){
  list_par = list(sigma2=sigma2,
                        alpha=alpha,
                        theta_0=theta_0,
                        par=custom,
                        nuisance = estimated_ME,
                        root_implied=fun(0, custom, theta_0))
  }else{
      list_par = list(sigma2=sigma2,
                            alpha=alpha,
                            theta_0=theta_0,
                            par=custom,
                            root_implied=fun(0, custom, theta_0))
  }
  
  # AIC
  AIC = -2*LL+2*nparam
  # AIC corrected
  AICc = AIC+((2*nparam*(nparam+1))/(nobs-nparam-1)) 
  
  # Check the curvature of the likelihood surface from the Hessian to assess the convergence (does not guarantee that it's a local optima...)
  if(method=="spg"){
      library(numDeriv)
      hess<-eigen(hessian(function(par) -llik(par, phy=tree, data=data, fun_ou=fun)$logl, param))$values
  }else{
      hess<-eigen(estimModel$hessian)$values
  }
  if(any(hess<0)) hess.value<-1 else hess.value<-0
  
  # Return the results
  results <- list(logl=LL, AIC=AIC, AICc=AICc, list_par=list_par, nb_param=nparam, estimated_ME=estimated_ME, opt=estimModel, hess.value=hess.value)

  class(results) <- "generalized.ou"
  return(results)
  #invisible(results)
}


# Print the results
print.generalized.ou<-function(x,...){
    cat("OU model with time-dependent optimum","\n")
    cat("AIC :",x$AIC,"\n")
    cat("AICc:",x$AICc,"\n")
    cat("Log-Likelihood:",x$logl,"\n")
    cat("______________________","\n")
    cat("Parameters:","\n")
    print(unlist(x$list_par), digits=3)
    cat("______________________","\n")
    if(x$opt$convergence==0){cat("Succesful convergence","\n")}else{cat("Convergence has not been reached","\n")}
    if(x$hess.value==0 & x$opt$convergence==0){cat("A reliable solution has been reached","\n")}else{cat("Unreliable solution (Likelihood at a saddle point)","\n")}
}

Stationary_to_scatter <- function(alpha,stationary){
  return(alpha*stationary + stationary*(alpha))
}

# # Test
# n=100
# tree <- rtree(n)
# tree$edge.length <- tree$edge.length * 60/ max(nodeHeights(tree))
# 
# # define a function for the optimum theta(s)
# fun_temp <- function(x, par)  par[1]*d18(max(nodeHeights(tree))-x) # here it's max-x because the time start at present
# 
# # simulate random data
# dat <- rTraitCont(tree)
# 
# # fit the model
# fit_OU_trend(tree, dat, fun=fun_temp, startvalues=c(1))
