# likelihood function
spatial.mixed.W.likelihood <- function(param, y, X, W1, W2) {
  n_W1 <- length(W1)
  n_W2 <- length(W2)
  n_param <- length(param)
  beta <- as.matrix(param[1:(n_param - n_W1 - n_W2 - 1)])
  if(n_W1 != 0){
    rhos <- param[(n_param - n_W1 - n_W2):(n_param - n_W2 - 1)] 
  }
  if(n_W2 != 0){
    lambdas <- param[(n_param - n_W2):(n_param - 1)] 
  }
  sigma <- param[length(param)]
  if(n_W1 != 0){
    N <- nrow(W1[[1]]) 
  }else{
    N <- nrow(W2[[1]])
  }
  Lm1 <- diag(N)
  if(n_W1 != 0){
    for(i in 1:n_W1){
      Lm1 <- Lm1 - rhos[[i]]*W1[[i]]
    } 
  }
  Mm1 <- diag(N)
  if(n_W2 != 0){
    for(i in 1:n_W2){
      Mm1 <- Mm1 - lambdas[[i]]*W2[[i]]
    } 
  }
  errors_squared <- t(Lm1 %*% y - X%*%beta) %*% t(Mm1) %*% Mm1 %*% (Lm1%*%y-X%*%beta)
  logLik <- log(det(Mm1)) + log(det(Lm1)) - N / 2 * log((sigma ^ 2) * 2 * pi) - (errors_squared) / (2 * (sigma ^ 2))
  return(logLik)
}

# estimating parameters for spatial mixed-W model
spatial.mixed.W <- function(y, X, W1=list(), W2=list(), n.it = 10000, prec = 10^(-16)) {
  n_W1 <- length(W1)
  n_W2 <- length(W2)
  start.df <- data.frame(y, X)
  if(n_W1 != 0){
    WW <- W1[[1]]
  }else{
    WW <- W2[[1]]
  }
  eval(parse(text = paste("start <- sacsarlm(", colnames(start.df)[1], " ~ ", paste(colnames(start.df[2:ncol(start.df)]), collapse = "+"), ", data = start.df, listw = WW, zero.policy = TRUE)", sep = "")))
  X <- as.matrix(data.frame(rep(1, nrow(X)), X))
  y <- as.matrix(y)
  beta.0 <- as.matrix(summary(start)$coefficients)
  k <- nrow(beta.0)
  rhos.0 <- rep(summary(start)$rho/n_W1, n_W1)
  lambdas.0 <- rep(summary(start)$lambda/n_W2, n_W2)
  if(n_W1 != 0){
    for(i in 1:n_W1){
      W1[[i]] <- listw2mat(W1[[i]])
    } 
  }
  if(n_W2 != 0){
    for(i in 1:n_W2){
      W2[[i]] <- listw2mat(W2[[i]])
    } 
  }
  sigma.0 <- start$s2 ^ 0.5
  par.0 <- c(beta.0, rhos.0, lambdas.0, sigma.0)
  par.1 <- optim(par.0, fn = spatial.mixed.W.likelihood,
                 method = "BFGS", hessian = TRUE,
                 control = list(trace = 1, maxit = n.it, fnscale = -1, reltol = prec),
                 y = y, X = X, W1 = W1, W2 = W2)
  
  par.loglik <- -par.1$value
  par.coef <- as.matrix(par.1$par)
  par.se <- as.matrix(diag(-solve(par.1$hessian)) ^ 0.5)
  par.p.value <- round(2 * pnorm(abs(par.coef / par.se), lower.tail = FALSE), digits = 5)
  beta.1 <- as.matrix(par.coef[1:k])
  if(n_W1 != 0){
    rhos.1 <- par.coef[(k+1):(k+n_W1)] 
  }
  if(n_W2 != 0){
    lambdas.1 <- par.coef[(k+1+n_W1):(k+n_W1+n_W2)] 
  }
  if(n_W1 != 0){
    N <- nrow(W1[[1]]) 
  }else{
    N <- nrow(W2[[1]])
  }
  Lm1 <- diag(N)
  if(n_W1 != 0){
    for(i in 1:n_W1){
      Lm1 <- Lm1 - rhos.1[i]*W1[[i]]
    } 
  }
  Mm1 <- diag(N)
  if(n_W2 != 0){
    for(i in 1:n_W2){
      Mm1 <- Mm1 - lambdas.1[i]*W2[[i]]
    } 
  }
  resid <- Mm1%*%(Lm1%*%y-X%*%beta.1)
  if((n_W1 != 0)&(n_W2 != 0)){
    results <- list(Coef = as.vector(beta.1), p.values = par.p.value[1:k, ], coef_se = par.se[1:k, ],
                    rhos = rhos.1, rhos.p.value = par.p.value[(k+1):(k+n_W1), ], rhos_se = par.se[(k+1):(k+n_W1), ],
                    lambdas = lambdas.1, lambdas.p.value = par.p.value[(k+1+n_W1):(k+n_W1+n_W2), ], lambdas_se = par.se[(k+1+n_W1):(k+n_W1+n_W2), ],
                    residuals = resid, loglik = par.loglik, names = c('(intercept)', colnames(X)[2:ncol(X)])) 
  }else if(n_W1 == 0){
    results <- list(Coef = as.vector(beta.1), p.values = par.p.value[1:k, ], coef_se = par.se[1:k, ],
                    lambdas = lambdas.1, lambdas.p.value = par.p.value[(k+1+n_W1):(k+n_W1+n_W2), ], lambdas_se = par.se[(k+1+n_W1):(k+n_W1+n_W2), ],
                    residuals = resid, loglik = par.loglik, names = c('(intercept)', colnames(X)[2:ncol(X)])) 
  }else{
    results <- list(Coef = as.vector(beta.1), p.values = par.p.value[1:k, ], coef_se = par.se[1:k, ],
                    rhos = rhos.1, rhos.p.value = par.p.value[(k+1):(k+n_W1), ], rhos_se = par.se[(k+1):(k+n_W1), ],
                    residuals = resid, loglik = par.loglik, names = c('(intercept)', colnames(X)[2:ncol(X)])) 
  }
  return(results)
}

# summary for spatial.mixed.W
summarize <- function(model){
  formals(print.data.frame)$row.names <- FALSE
  df <- data.frame(model$names, round(model$Coef, 3), round(model$coef_se, 3), round(model$p.values, 3))
  df['signif.'] <- ifelse(model$p.values <= 0.01, '***', 
                          ifelse(model$p.values <= 0.05, '** ', 
                                 ifelse(model$p.values <= 0.1, '*  ', '')))
  colnames(df) <- c('variable', 'coeff', 'se', 'p-value', ' ')
  print(df)
  cat('\n')
  
  if(!is.null(model$rhos)){
    r <- data.frame(round(model$rhos, 3), round(model$rhos_se, 3), round(model$rhos.p.value, 3))
    r['signif.'] <- ifelse(model$rhos.p.value <= 0.01, '***', 
                           ifelse(model$rhos.p.value <= 0.05, '** ', 
                                  ifelse(model$rhos.p.value <= 0.1, '*  ','')))
    colnames(r) <- c('coeff', 'se', 'p-value', '')
    rownames(r) <- NULL
    cat('Rho: ', sep = '\n')
    print(r)
    cat('\n')
  }
  
  if(!is.null(model$lambdas)){
    l <- data.frame(round(model$lambdas, 3), round(model$lambdas_se, 3), round(model$lambdas.p.value, 3))
    l['signif.'] <- ifelse(model$lambdas.p.value <= 0.01, '***', 
                           ifelse(model$lambdas.p.value <= 0.05, '** ', 
                                  ifelse(model$lambdas.p.value <= 0.1, '*  ','')))
    colnames(l) <- c('coeff', 'se', 'p-value', '')
    rownames(l) <- NULL
    cat('Lambda: ', sep = '\n')
    print(l) 
  }
}