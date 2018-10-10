# uso:
# IWTlmpointwise(formula)
# formula: oggetto di tipo formula. 
# Esempio con due covariate: IWTlmFoF(y ~ x1 + x2)
# dove y è un dataset funzionale (matrice n*p)
# x1 e x2 possono essere funzionali (matrici n*p) o scalari (vettori n)
# il test dell'interazione tra covariate scalari e funzionali (tipo y ~ x1*x2) non funziona, per ora va fatto a mano costruendo le covariate dummy
IWTlmFoF <- function(formula,B=1000,method='residuals'){
  pval.correct <- function(pval.matrix){
    matrice_pval_2_2x <- cbind(pval.matrix,pval.matrix)
    p <- dim(pval.matrix)[2]
    matrice_pval_2_2x <- matrice_pval_2_2x[,(2*p):1]
    corrected.pval <- numeric(p)
    for(var in 1:p){
      pval_var <- matrice_pval_2_2x[p,var]
      inizio <- var
      fine <- var 
      for(riga in (p-1):1){
        fine <- fine + 1
        pval_cono <- matrice_pval_2_2x[riga,inizio:fine]
        pval_var <- max(pval_var,pval_cono)
      }
      corrected.pval[var] <- pval_var
    }
    corrected.pval <- corrected.pval[p:1]
    return(corrected.pval)
  }
  
  extract.residuals = function(regr){
    return(regr$residuals)
  }
  extract.fitted = function(regr){
    return(regr$fitted)
  }
  extract.coeff = function(regr){
    return(regr$coefficients)
  }
  env <- environment(formula)
  cl <- match.call()
  design.matrix.temp = model.matrix(formula)
  mf = model.frame(formula)
  data = model.response(mf)
  coeff <- data
  
  #nvar = dim(design.matrix.temp)[2] - 1
  variables = attr(terms(formula),"term.labels") #all.vars(formula)[-1] #colnames(design.matrix.temp)
  y.name = all.vars(formula)[1]
  n <- dim(data)[1]
  J <- dim(data)[2]
  
  assign <- attr(design.matrix.temp,'assign')
  contrast <- attr(design.matrix.temp,'contrast')
  
  length.vars <- numeric(length(variables)+1)
  for(var in 0:(length(variables))){
    length.vars[var+1] <- sum(assign==var)
  }
  
  nvar <- sum(length.vars==J) + sum(length.vars[which(length.vars!=J)]) - 1 
  var.functional <- which(length.vars==J) - 1
  var.scalar <- which(length.vars!=J) - 1
  index.scalar <- NULL
  for(ii in var.scalar){
    index.scalar <- c(index.scalar,which(assign==ii))
  }
  #index.functional <- matrix(data=(1:dim(design.matrix.temp)[2])[-index.scalar], nrow=length(var.functional),ncol=J,byrow=TRUE)
  
  if(length(index.scalar)>0){
    index.functional <- matrix(data=(1:dim(design.matrix.temp)[2])[-index.scalar], nrow=length(var.functional),ncol=J,byrow=TRUE)
  }else{
    index.functional <- matrix(data=(1:dim(design.matrix.temp)[2]), nrow=length(var.functional),ncol=J,byrow=TRUE)
  }
  
  var.names.scalar = colnames(design.matrix.temp)[index.scalar]
  var.names.functional <- variables[var.functional]
  var.names <- c(var.names.scalar,var.names.functional)
  design.matrix <- list(J)
  for(jj in 1:J){
    design.matrix.scalar <- design.matrix.temp[,index.scalar]
    design.matrix.functional <- design.matrix.temp[,index.functional[,jj]]
    design.matrix[[jj]] <- list()
    design.matrix[[jj]]$design.matrix <- cbind(design.matrix.scalar,design.matrix.functional)
    design.matrix[[jj]]$y <- coeff[,jj]
  }
  
  ############################################################
  
  print('First step: basis expansion')
  #splines coefficients:
  coeff <- eval <- data.eval
  p <- dim(coeff)[2]
  
  
  
  print('Second step: joint univariate tests')
  #univariate permutations
  #regr0old = lm.fit(design.matrix[[1]],coeff)
  lm.fit.mod <- function(X,var.keep){
    x=X$design.matrix[,var.keep,drop = F]
    return(.lm.fit(x,y=X$y))
  }
  regr0 <- lapply(design.matrix,lm.fit.mod,var.keep=1:(nvar+1))
  
  #test statistics:
  #Sigmaold <- chol2inv(regr0old$qr$qr)
  calculateTpart <- function(regr,nvar){
    Sigma <- chol2inv(regr$qr$qr)
    resvar <- sum(regr$residuals^2)/regr$df.residual
    se <- sqrt( matrix(diag(Sigma),nrow=(nvar+1),ncol=1,byrow=FALSE) * matrix(resvar,nrow=(nvar+1),ncol=1,byrow=TRUE))
    T_part_temp <- abs(regr$coeff / se)
    return(T_part_temp)
  }
  calculateTglob <- function(regr,nvar){
    Sigma <- chol2inv(regr$qr$qr)
    resvar <- sum(regr$residuals^2)/regr$df.residual
    T0_glob <- colSums((regr$fitted - matrix(mean(regr$fitted),nrow=n,ncol=1,byrow=TRUE))^2)/ ((nvar)*resvar)
  }
  T0_part <- sapply(regr0,calculateTpart,nvar=nvar)
  rownames(T0_part) <- var.names
  
  if(nvar >0){
    T0_glob <- sapply(regr0,calculateTglob,nvar=nvar)
  }else{
    method = 'responses' 
    T0_glob = numeric(p)
    T0_part = t(as.matrix(T0_part))
  }
  
  #calculate residuals
  if(method=='residuals'){
    #n residuals for each coefficient of basis expansion (1:p) 
    #and for each partial test + global test (nvar+1) 
    #saved in array of dim (nvar+1,n,p)
    #formula.const <- deparse(formula[[3]],width.cutoff = 500L) #extracting the part after ~ on formula. this will not work if the formula is longer than 500 char
    #design.matrix.names2 = design.matrix
    var.names2 = var.names
    
    residui = array(dim=c(nvar+1,n,p))
    fitted_part = array(dim=c(nvar+1,n,p)) 
    regr0_part = vector('list',nvar+1)
    #coeff.perm_part = array(dim=c(nvar+1,n,p))
    for(ii in 1:(nvar+1)){ #the first one is the intercept. treated as special case after loop
      var.ii = var.names2[ii]
      variables.reduced = var.names2[-c(1,which(var.names2==var.ii))] 
      index.keep <- (1:(nvar+1))[-ii]
      
      regr0_part[[ii]] = lapply(design.matrix,lm.fit.mod,var.keep=index.keep)
      
      residui[ii,,] = simplify2array(lapply(regr0_part[[ii]],extract.residuals))
      fitted_part[ii,,] = simplify2array(lapply(regr0_part[[ii]],extract.fitted))
    }
  }
  
  T_glob <- matrix(ncol=p,nrow=B)
  T_part = array(dim=c(B,nvar+1,p))
  
  for (perm in 1:B){
    # the F test is the same for both methods
    if(nvar >0){
      permutazioni <- sample(n)
      coeff_perm <- coeff[permutazioni,]
    }else{ # test on intercept permuting signs
      signs <- rbinom(n,1,0.5)*2 - 1
      coeff_perm <- coeff*signs
    }
    design.matrix.perm <- list(J)
    for(jj in 1:J){
      design.matrix.scalar <- design.matrix.temp[,index.scalar]
      design.matrix.functional <- design.matrix.temp[,index.functional[,jj]]
      design.matrix.perm[[jj]] <- list()
      design.matrix.perm[[jj]]$design.matrix <- cbind(design.matrix.scalar,design.matrix.functional)
      design.matrix.perm[[jj]]$y <- coeff_perm[,jj]
    }
    
    regr_perm = lapply(design.matrix.perm,lm.fit.mod,var.keep=1:(nvar+1))
    
    if(nvar > 0)
      T_glob[perm,] <- sapply(regr_perm,calculateTglob,nvar=nvar)
    
    # partial tests: differ depending on the method
    if(method=='responses'){
      T_part[perm,,] <- sapply(regr_perm,calculateTpart,nvar=nvar)
    }else if(method=='residuals'){
      residui_perm = residui[,permutazioni,]
      regr_perm_part = vector('list',nvar+1)
      for(ii in 1:(nvar+1)){ 
        coeff_perm = fitted_part[ii,,] + residui_perm[ii,,] 
        design.matrix.perm <- list(J)
        for(jj in 1:J){
          design.matrix.scalar <- design.matrix.temp[,index.scalar]
          design.matrix.functional <- design.matrix.temp[,index.functional[,jj]]
          design.matrix.perm[[jj]] <- list()
          design.matrix.perm[[jj]]$design.matrix <- cbind(design.matrix.scalar,design.matrix.functional)
          design.matrix.perm[[jj]]$y <- coeff_perm[,jj]
        }
        regr_perm = lapply(design.matrix.perm,lm.fit.mod,var.keep=1:(nvar+1))
        T_part[perm,ii,] = sapply(regr_perm,calculateTpart,nvar=nvar)[ii,]
      }
    }
  }
  
  pval_glob <- numeric(p)
  pval_part = matrix(nrow=nvar+1,ncol=p)
  for(i in 1:p){
    pval_glob[i] <- sum(T_glob[,i]>=T0_glob[i])/B
    pval_part[,i] = colSums(T_part[,,i]>=matrix(T0_part[,i],nrow=B,ncol=nvar+1,byrow=TRUE))/B
  }
  
  #combination
  print('Third step: interval-wise combination and correction')
  
  #asymmetric combination matrix:
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
  
  #######################
  
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
  
  #symmetric combination matrix
  matrice_pval_symm_glob <- matrix(nrow=p,ncol=4*p)
  matrice_pval_symm_part <- array(dim=c(nvar+1,p,4*p))
  
  for(i in 0:(p-1)){
    for(j in 1:(2*p)){
      matrice_pval_symm_glob[p-i,j+i+p] <- matrice_pval_asymm_glob[p-i,(j+1)%/%2]
      if(j+i>2*p-i){
        matrice_pval_symm_glob[p-i,j+i-p] <- matrice_pval_asymm_glob[p-i,(j+1)%/%2]
      }
      for(ii in 1:(nvar+1)){
        matrice_pval_symm_part[ii,p-i,j+i+p] <- matrice_pval_asymm_part[ii,p-i,(j+1)%/%2]
        if(j+i>2*p-i){
          matrice_pval_symm_part[ii,p-i,j+i-p] <- matrice_pval_asymm_part[ii,p-i,(j+1)%/%2]
        }
      }
    }
  }
  
  corrected.pval_glob <- pval.correct(matrice_pval_asymm_glob)
  corrected.pval_part = matrix(nrow=nvar+1,ncol=p)
  for(ii in 1:(nvar+1)){
    corrected.pval_part[ii,] = pval.correct(matrice_pval_asymm_part[ii,,])
  }
  
  coeff.regr = sapply(regr0,extract.coeff)
  coeff.t <- ((coeff.regr))
  
  fitted.regr = sapply(regr0,extract.fitted)
  fitted.t <- (fitted.regr)
  
  rownames(corrected.pval_part) = var.names
  rownames(coeff.t) = var.names
  rownames(coeff.regr) = var.names
  rownames(pval_part) = var.names
  
  residuals.t = data.eval - fitted.t
  ybar.t = colMeans(data.eval)
  npt <- J
  R2.t = colSums((fitted.t - matrix(data=ybar.t,nrow=n,ncol=npt,byrow=TRUE))^2)/colSums((data.eval - matrix(data=ybar.t,nrow=n,ncol=npt,byrow=TRUE))^2)
  
  print('Interval Testing Procedure completed')
  
  IWTresult <- list(call=cl,basis='B-spline',
                    unadjusted_pval_F=pval_glob,
                    pval_matrix_F=matrice_pval_asymm_glob,
                    adjusted_pval_F=corrected.pval_glob,
                    unadjusted_pval_part=pval_part,
                    pval_matrix_part=matrice_pval_asymm_part,
                    adjusted_pval_part=corrected.pval_part,
                    data.eval=data.eval,
                    coeff.regr.eval=coeff.t,
                    fitted.eval=fitted.t,
                    residuals.eval=residuals.t,
                    R2.eval=R2.t,
                    heatmap.matrix.F=matrice_pval_symm_glob,
                    heatmap.matrix.t=matrice_pval_symm_part)
  class(IWTresult) = 'IWTlm'
  return(IWTresult)
}


summary.IWTlm = function(object, ...){
  printresult = vector('list')
  printresult$call = object$call
  #class(printresult) <- "lm"
  printresult$ttest = matrix(data=apply(object$adjusted_pval_part,1,min),ncol=1)
  var.names = rownames(object$adjusted_pval_part)
  rownames(printresult$ttest) = var.names
  printresult$ttest = as.data.frame(printresult$ttest)
  signif = rep('',length(var.names))
  signif[which(printresult$ttest[,1] <0.001)] = '***'
  signif[which(printresult$ttest[,1] <0.01 & printresult$ttest[,1] >= 0.001)] = '**'
  signif[which(printresult$ttest[,1] <0.05 & printresult$ttest[,1] >= 0.01)] = '*'
  signif[which(printresult$ttest[,1] <0.1 & printresult$ttest[,1] >= 0.05)] = '.'
  printresult$ttest[,2] = signif
  colnames(printresult$ttest) = c('Minimum p-value','')
  
  printresult$R2 = as.matrix(range(object$R2.eval))
  colnames(printresult$R2) = 'Range of functional R-squared'
  rownames(printresult$R2) = c('Min R-squared', 'Max R-squared')
  printresult$ftest = as.matrix(min(object$adjusted_pval_F))
  printresult$ftest = as.data.frame(printresult$ftest)
  signif.f = ''
  signif.f[which(printresult$ftest[,1] <0.001)] = '***'
  signif.f[which(printresult$ftest[,1] <0.01 & printresult$ftest[,1] >= 0.001)] = '**'
  signif.f[which(printresult$ftest[,1] <0.05 & printresult$ftest[,1] >= 0.01)] = '*'
  signif.f[which(printresult$ftest[,1] <0.1 & printresult$ftest[,1] >= 0.05)] = '.'
  printresult$ftest[,2] = signif.f
  colnames(printresult$ftest) = c('Minimum p-value','')
  printresult
  
}



