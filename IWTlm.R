# Interval Testing Procedure for the Linear Models (functional response, 
# coefficients and scalar covariates) framework using B-splines
IWTlm <- function(formula,B = 1000,method = 'residuals',dx=NULL,recycle=TRUE) {
  
  pval.correct <- function(pval.matrix){
    matrice_pval_2_2x <- cbind(pval.matrix,pval.matrix)
    p <- dim(pval.matrix)[2]
    matrice_pval_2_2x <- matrice_pval_2_2x[,(2*p):1]
    corrected.pval <- numeric(p)
    corrected.pval.matrix <- matrix(nrow=p,ncol=p)
    corrected.pval.matrix[p,] <- pval.matrix[p,p:1]
    for(var in 1:p){
      pval_var <- matrice_pval_2_2x[p,var]
      inizio <- var
      fine <- var #inizio fisso, fine aumenta salendo nelle righe
      for(riga in (p-1):1){
        fine <- fine + 1
        pval_cono <- matrice_pval_2_2x[riga,inizio:fine]
        pval_var <- max(pval_var,pval_cono,na.rm=TRUE)
        corrected.pval.matrix[riga,var] <- pval_var
      }
      corrected.pval[var] <- pval_var
    }
    corrected.pval <- corrected.pval[p:1]
    corrected.pval.matrix <- corrected.pval.matrix[,p:1]
    return(corrected.pval.matrix)
  }
  extract_residuals <- function(regr){
    return(regr$residuals)
  }
  extract_fitted <- function(regr){
    return(regr$fitted)
  }
  
  env <- environment(formula)
  variables = all.vars(formula)
  y.name = variables[1]
  covariates.names <- colnames(attr(terms(formula),"factors"))
  #data.all = model.frame(formula)
  cl <- match.call()
  data <- get(y.name,envir = env)
  if(is.fd(data)){ # data is a functional data object
    rangeval <- data$basis$rangeval
    if(is.null(dx)){
      dx <- (rangeval[2]-rangeval[1])*0.01
    }
    abscissa <- seq(rangeval[1],rangeval[2],by=dx)
    coeff <- t(eval.fd(fdobj=data,evalarg=abscissa))
  }else if(is.matrix(data)){
    coeff <- data
  }else{
    stop("First argument of the formula must be either a functional data object or a matrix.")
  }
  
  dummynames.all <- colnames(attr(terms(formula),"factors"))
  formula.const <- deparse(formula[[3]],width.cutoff = 500L) #extracting the part after ~ on formula. this will not work if the formula is longer than 500 char
  
  formula.discrete <- as.formula(paste('coeff ~',formula.const),env=environment())
  design_matrix = model.matrix(formula.discrete)
  mf = model.frame(formula.discrete)
  
  nvar <- dim(design_matrix)[2] - 1
  var_names <- colnames(design_matrix)
  p <- dim(coeff)[2]
  n <- dim(coeff)[1]
  # Univariate permutations
  regr0 <- lm.fit(design_matrix, coeff)
  # Test statistics
  Sigma <- chol2inv(regr0$qr$qr)
  resvar <- colSums(regr0$residuals ^ 2) / regr0$df.residual
  se <- sqrt(matrix(diag(Sigma), nrow = nvar + 1, ncol = p, byrow = FALSE) 
        * matrix(resvar, nrow = nvar + 1, ncol = p, byrow = TRUE))
  T0_part <- abs(regr0$coeff / se)^2
  if (nvar > 0) {
    T0_glob <- colSums((regr0$fitted - matrix(colMeans(regr0$fitted),
                                              nrow = n, ncol = p, 
                                              byrow = TRUE)) ^ 2) / (nvar * resvar)
  } else {
    method <- 'responses' 
    T0_glob <- numeric(p)
    T0_part <- t(as.matrix(T0_part))
  }
  # Compute residuals
  if (method == 'residuals') {
    # n residuals for each coefficient of basis expansion (1:p) 
    # and for each partial test + global test (nvar+1) 
    # Saved in array of dim (nvar+1,n,p)
    # Extracting the part after ~ on formula. 
    # This will not work if the formula 
    # is longer than 500 char
    formula_const <- deparse(formula[[3]], width.cutoff = 500L)
    design_matrix_names2 <- design_matrix
    var_names2 <- var_names
    coeffnames <- paste('coeff[,', as.character(1:p),']', sep = '')
    formula_temp <- coeff ~ design_matrix
    mf_temp <- cbind(model.frame(formula_temp)[-((p + 1):(p + nvar + 1))], 
                    as.data.frame(design_matrix[, -1]))
    if (length(grep('factor', formula_const)) > 0) {
      index_factor <- grep('factor', var_names)
      replace_names <- paste('group', (1:length(index_factor)), sep = '')
      var_names2[index_factor] <- replace_names
      colnames(design_matrix_names2) <- var_names2
    }
    residui <- array(dim=c(nvar + 1, n, p))
    fitted_part <- array(dim = c(nvar + 1, n, p)) 
    formula_coeff_part <- vector('list', nvar + 1)
    regr0_part <- vector('list',nvar + 1)
    # The first one is the intercept. Treated as special case after loop
    for (ii in 2:(nvar + 1)) { 
      var_ii <- var_names2[ii]
      variables_reduced <- var_names2[-c(1, which(var_names2 == var_ii))] 
      if (nvar > 1) {
        formula_temp <- paste(variables_reduced, collapse = ' + ')
      } else {
        # Removing the unique variable -> reduced model only has intercept ter
        formula_temp <- '1' 
      }
      formula_temp2 <- coeff ~ design_matrix_names2
      mf_temp2 <- cbind(model.frame(formula_temp2)[-((p + 1):(p + nvar + 1))], 
                       as.data.frame(design_matrix_names2[,-1]))
      formula_coeff_temp <- paste(coeffnames, '~', formula_temp) 
      formula_coeff_part[[ii]] <- sapply(formula_coeff_temp, as.formula)
      regr0_part[[ii]] <- lapply(formula_coeff_part[[ii]], lm, data = mf_temp2)
      residui[ii, , ] <- simplify2array(lapply(regr0_part[[ii]], extract_residuals))
      fitted_part[ii, , ] <- simplify2array(lapply(regr0_part[[ii]], extract_fitted))
    }
    ii <- 1 # intercept
    formula_temp <- paste(formula_const, ' -1', sep = '')
    formula_coeff_temp <- paste(coeffnames, '~', formula_temp)
    formula_coeff_part[[ii]] <- sapply(formula_coeff_temp, as.formula)
    regr0_part[[ii]] <- lapply(formula_coeff_part[[ii]], lm, data = mf_temp)
    residui[ii, , ] <- simplify2array(lapply(regr0_part[[ii]], extract_residuals))
    fitted_part[ii, , ] <- simplify2array(lapply(regr0_part[[ii]], extract_fitted))
  }
  print('Point-wise tests')
  
  # CMC algorithm
  T_glob <- matrix(ncol = p,nrow = B)
  T_part <- array(dim = c(B, nvar + 1, p))
  for (perm in 1:B) {
    # the F test is the same for both methods
    if (nvar > 0) {
      permutazioni <- sample(n)
      coeff_perm <- coeff[permutazioni, ]
    }else{ # Test on intercept permuting signs
      signs <- rbinom(n, 1, 0.5) * 2 - 1
      coeff_perm <- coeff * signs
    }
    regr_perm <- lm.fit(design_matrix, coeff_perm)
    Sigma <- chol2inv(regr_perm$qr$qr)
    resvar <- colSums(regr_perm$residuals ^ 2) / regr_perm$df.residual
    if (nvar > 0) {
      T_glob[perm, ] <- colSums((regr_perm$fitted - matrix(colMeans(regr_perm$fitted), 
                                                           nrow = n, ncol = p,
                                                           byrow=TRUE)) ^ 2)/ (nvar * resvar)
    }
    # Partial tests: differ depending on the method
    if (method == 'responses') {
      se <- sqrt(matrix(diag(Sigma), nrow = nvar + 1, ncol = p, byrow = FALSE) 
                 * matrix(resvar, nrow = nvar + 1, ncol = p, byrow = TRUE))
      T_part[perm, , ] <- abs(regr0$coeff / se)^2
    } else if (method == 'residuals'){
      residui_perm <- residui[, permutazioni, ]
      regr_perm_part <- vector('list', nvar + 1)
      for (ii in 1:(nvar + 1)) { 
        coeff_perm <- fitted_part[ii, , ] + residui_perm[ii, , ]  
        regr_perm <- lm.fit(design_matrix, coeff_perm)
        Sigma <- chol2inv(regr_perm$qr$qr)
        resvar <- colSums(regr_perm$residuals ^ 2) / regr_perm$df.residual
        se <- sqrt(matrix(diag(Sigma), nrow = nvar + 1 , ncol = p, byrow = FALSE) 
                   * matrix(resvar, nrow = nvar + 1, ncol = p, byrow = TRUE))
        T_part[perm, ii, ] <- abs(regr_perm$coeff / se)[ii, ]^2
      }
    }
  }
  pval_glob <- numeric(p)
  pval_part <- matrix(nrow = nvar + 1, ncol = p)
  for (i in 1:p) {
    pval_glob[i] <- sum(T_glob[, i] >= T0_glob[i]) / B
    pval_part[, i] <- colSums(T_part[, , i] >= matrix(T0_part[, i],nrow = B, ncol = nvar + 1,byrow = TRUE)) / B
  }
  
  print('Interval-wise tests')
  
  matrice_pval_asymm_glob <- matrix(nrow=p,ncol=p)
  matrice_pval_asymm_glob[p,] <- pval_glob[1:p]
  pval_2x_glob <- c(pval_glob,pval_glob)
  T0_2x_glob <- c(T0_glob,T0_glob)
  T_2x_glob <- cbind(T_glob,T_glob)
  
  
  matrice_pval_asymm_part <- array(dim=c(nvar+1,p,p))
  pval_2x_part <- cbind(pval_part,pval_part)
  T0_2x_part <- cbind(T0_part,T0_part)
  T_2x_part = array(dim = c(B,nvar+1, p*2))
  for(ii in 1:(nvar+1)){
    matrice_pval_asymm_part[ii,p,] <- pval_part[ii,1:p]
    T_2x_part[,ii,] <- cbind(T_part[,ii,],T_part[,ii,])
  }
  
  if(recycle==TRUE){
    for(i in (p-1):1){
      for(j in 1:p){
        inf <- j
        sup <- (p-i)+j
        T0_temp <- sum(T0_2x_glob[inf:sup])
        T_temp <- rowSums(T_2x_glob[,inf:sup])
        pval_temp <- sum(T_temp>=T0_temp)/B
        matrice_pval_asymm_glob[i,j] <- pval_temp
        for(ii in 1:(nvar + 1)){
          T0_temp <- sum(T0_2x_part[ii,inf:sup])
          T_temp <- rowSums(T_2x_part[,ii,inf:sup])
          pval_temp <- sum(T_temp>=T0_temp)/B
          matrice_pval_asymm_part[ii,i,j] <- pval_temp
        }
        
      }
      print(paste('creating the p-value matrix: end of row ',as.character(p-i+1),' out of ',as.character(p),sep=''))
    }
    
  }else{
    for(i in (p-1):1){
      for(j in 1:i){
        inf <- j
        sup <- (p-i)+j
        T0_temp <- sum(T0_2x_glob[inf:sup])
        T_temp <- rowSums(T_2x_glob[,inf:sup])
        pval_temp <- sum(T_temp>=T0_temp)/B
        matrice_pval_asymm_glob[i,j] <- pval_temp
        for(ii in 1:(nvar + 1)){
          T0_temp <- sum(T0_2x_part[ii,inf:sup])
          T_temp <- rowSums(T_2x_part[,ii,inf:sup])
          pval_temp <- sum(T_temp>=T0_temp)/B
          matrice_pval_asymm_part[ii,i,j] <- pval_temp
        }
        
      }
      print(paste('creating the p-value matrix: end of row ',as.character(p-i+1),' out of ',as.character(p),sep=''))
    }
    
  }
  
  corrected.pval.matrix_glob <- pval.correct(matrice_pval_asymm_glob)
  corrected.pval_glob <- corrected.pval.matrix_glob[1,]
  
  corrected.pval_part = matrix(nrow=nvar+1,ncol=p)  
  corrected.pval.matrix_part = array(dim=c(nvar+1,p,p))
  for(ii in 1:(nvar+1)){
    corrected.pval.matrix_part[ii,,] = pval.correct(matrice_pval_asymm_part[ii,,])
    corrected.pval_part[ii,] <- corrected.pval.matrix_part[ii,1,]
  }
  
  coeff.regr = regr0$coeff
  coeff.t <- ((coeff.regr))
  
  fitted.regr = regr0$fitted
  fitted.t <- (fitted.regr)
  
  rownames(corrected.pval_part) = var_names
  rownames(coeff.t) = var_names
  rownames(coeff.regr) = var_names
  rownames(pval_part) = var_names
  
  data.eval <- coeff
  residuals.t = data.eval - fitted.t
  ybar.t = colMeans(data.eval)
  npt <- p
  R2.t = colSums((fitted.t - matrix(data=ybar.t,nrow=n,ncol=npt,byrow=TRUE))^2)/colSums((data.eval - matrix(data=ybar.t,nrow=n,ncol=npt,byrow=TRUE))^2)
  
  print('Interval-Wise Testing completed')
  
  
  ITP_result <- list(call=cl,
                     design_matrix=design_matrix,
                     unadjusted_pval_F=pval_glob,
                     pval_matrix_F=matrice_pval_asymm_glob,
                     adjusted_pval_F=corrected.pval_glob,
                     unadjusted_pval_part=pval_part,
                     pval_matrix_part=matrice_pval_asymm_part,
                     adjusted_pval_part=corrected.pval_part,
                     data.eval=coeff,
                     coeff.regr.eval=coeff.t,
                     fitted.eval=fitted.t,
                     residuals.eval=residuals.t,
                     R2.eval=R2.t)
  class(ITP_result) = 'IWTlm'
  return(ITP_result)
}
