summary.IWTlm <-
function(object, ...){
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
