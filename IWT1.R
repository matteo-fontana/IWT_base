IWT1 <- function(data,mu=0,B=1000,dx=NULL,recycle=TRUE){
  pval.correct <- function(pval.matrix){
    matrice_pval_2_2x <- cbind(pval.matrix,pval.matrix)
    p <- dim(pval.matrix)[2]
    matrice_pval_2_2x <- matrice_pval_2_2x[,(2*p):1]
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
    }
    corrected.pval.matrix <- corrected.pval.matrix[,p:1]
    return(corrected.pval.matrix)
  }
  
  # data preprocessing
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
    stop("First argument must be either a functional data object or a matrix.")
  }
  
  if (is.fd(mu)){ # mu is a functional data
    rangeval.mu <- mu$basis$rangeval
    if(sum(rangeval.mu == rangeval)!=2){
      stop("rangeval of mu must be the same as rangeval of data.")
    }
    if(is.null(dx)){
      dx <- (rangeval.mu[2]-rangeval.mu[1])*0.01
    }
    abscissa <- seq(rangeval.mu[1],rangeval.mu[2],by=dx)
    mu.eval <- t(eval.fd(fdobj=mu,evalarg=abscissa))
  }else if(is.vector(mu)){
    mu.eval <- mu
  }else{
    stop("Second argument must be either a functional data object or a numeric vector.")
  }
  
  n <- dim(coeff)[1]
  p <- dim(coeff)[2]
  data.eval <- coeff <- coeff - matrix(data=mu.eval,nrow=n,ncol=p,byrow=TRUE)
  
  #univariate permutations
  print('Point-wise tests')
  T0 <- abs(colMeans(coeff))^2  #sample mean
  T_coeff <- matrix(ncol=p,nrow=B)
  for (perm in 1:B){
    signs <- rbinom(n,1,0.5)*2 - 1
    coeff_perm <- coeff*signs
    T_coeff[perm,] <- abs(colMeans(coeff_perm))^2
  }
  pval <- numeric(p)
  for(i in 1:p){
    pval[i] <- sum(T_coeff[,i]>=T0[i])/B
  }
  
  #combination
  print('Interval-wise tests')
  
  #asymmetric combination matrix:
  matrice_pval_asymm <- matrix(nrow=p,ncol=p)
  matrice_pval_asymm[p,] <- pval[1:p]
  T0_2x <- c(T0,T0)
  T_coeff_2x <- cbind(T_coeff,T_coeff)
  
  maxrow <- 1
  # con parametro scale
  #maxrow <- p-scale+1
  
  if(recycle==TRUE){
    for(i in (p-1):maxrow){ # rows
      for(j in 1:p){ # columns
        inf <- j
        sup <- (p-i)+j
        T0_temp <- sum(T0_2x[inf:sup])
        T_temp <- rowSums(T_coeff_2x[,inf:sup])
        pval_temp <- sum(T_temp>=T0_temp)/B
        matrice_pval_asymm[i,j] <- pval_temp
      }
      print(paste('creating the p-value matrix: end of row ',as.character(p-i+1),' out of ',as.character(p),sep=''))
    }
  }else{ # without recycling
    for(i in (p-1):maxrow){ # rows
      for(j in 1:i){ # columns
        inf <- j
        sup <- (p-i)+j
        T0_temp <- sum(T0_2x[inf:sup])
        T_temp <- rowSums(T_coeff_2x[,inf:sup])
        pval_temp <- sum(T_temp>=T0_temp)/B
        matrice_pval_asymm[i,j] <- pval_temp
      }
      print(paste('creating the p-value matrix: end of row ',as.character(p-i+1),' out of ',as.character(p),sep=''))
    }
  }
  
  corrected.pval.matrix <- pval.correct(matrice_pval_asymm)
  corrected.pval <- corrected.pval.matrix[1,]
  
  print('Interval-Wise Testing completed')
  IWT.result <- list(
    test = '1pop', mu = mu.eval,
    adjusted_pval = corrected.pval,
    unadjusted_pval = pval,
    pval_matrix = matrice_pval_asymm,
    data.eval=data.eval)
  class(IWT.result) = 'IWT1'
  return(IWT.result)
}
