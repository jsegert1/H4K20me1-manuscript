#' Compute ROC
#'
#' @param elements: GRanges object of silencer screen results
#' @param filename: filename of histone ChIP-seq peaks 
#'
#' @return
#' @export
#'
#' @examples
compute_roc = function(elements, filename){
  mark = readPeakFile(filename)
  #make table
  roc_table = data.frame()
  
  
  #intersect with mark
  elements$mark = !is.na(findOverlaps(elements, mark, select = 'first'))
  quantiles = quantile(elements$Fold_change, seq(0, 1,.01))
  
  
  for (i in 1:length(quantiles)){
    #TP = sum(mark+ and silencer)
    TP = sum(elements$mark & elements$Fold_change >= quantiles[i])
    
    #FP= sum(mark- and silencer)
    FP = sum(!elements$mark & elements$Fold_change >= quantiles[i])
    
    #FN = sum(mark+ and nonsilencer)
    FN = sum(elements$mark & elements$Fold_change < quantiles[i])
    
    #TN = sum(mark- and non-silencer)
    TN = sum(!elements$mark & elements$Fold_change < quantiles[i])
    
    TPR = TP/(TP+FN)
    FPR = FP/(TN+FP)
    roc_table = rbind(roc_table, c(names(quantiles)[i], TPR, FPR, TP, FP, FN, TN))
  }
  colnames(roc_table) = c('quantile', 'TPR', 'FPR', 'TP', 'FP', 'FN', 'TN')
  roc_table$TPR = as.numeric(roc_table$TPR)
  roc_table$FPR = as.numeric(roc_table$FPR)
  return(roc_table)
}