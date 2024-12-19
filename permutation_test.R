library(ChIPseeker)
library(GenomicRanges)
library(regioneR)
#' ReSE permutation test
#'
#' @param fg_elements: GRanges object of silencer elements
#' @param bg_elements: GRanges object of background elements
#' @param annotationFile: path to bed file for peaks of interest or GRanges of bed file
#' @param numPermutations: number of permutations to run. Default: 20,000
#'
#' @return list of:
#' 
#' p_val: empirical p value of enrichment
#' 
#' cover_ratio: proportion of silencer elements overlapping the given annotation
#' 
#' fold_change: the cover ratio of the foreground over the cover ratio of background elements
#' 
#' background: vector of overlaps for the background sets
histone_permutation_test = function(fg_elements, bg_elements, annotationFile, numPermutations = 20000, blacklist = NA){
  if (is(annotationFile, 'GRanges')){
    mark = annotationFile
  }else{
  mark = readPeakFile(annotationFile)
  }
  
  fg_elements$intersect = !is.na(findOverlaps(fg_elements, mark, select = 'first'))
  bg_elements$intersect = !is.na(findOverlaps(bg_elements, mark, select = 'first'))
  
  fg_cover_ratio = sum(fg_elements$intersect)/ length(fg_elements)
  
  bg_cover_ratios = replicate(numPermutations, sum(sample(bg_elements$intersect, length(fg_elements))))/ length(fg_elements)
  
  p_val = sum(bg_cover_ratios > fg_cover_ratio)/ numPermutations
  
  fold_change = fg_cover_ratio/ mean(bg_cover_ratios)
  
  results = list('p_val' = p_val, 'cover_ratio' = fg_cover_ratio, 'fold_change' = fold_change, 'background' = bg_cover_ratios)
  return(results)
}


#' Histone Permutation Test for STARR-seq
#'
#' @param elements: GRanges object of silencer elements
#' @param annotationFile: path to bed file for peaks of interest or GRanges of bed file
#' @param numPermutations: number of permutations to run. Default: 20,000
#' @param blacklist: GRanges object of blacklisted regions
#'
#' @return
#' @export
#'
#' @examples
histone_permutation_test_starr = function(elements, annotationFile, numPermutations = 20000, blacklist = NA){
  if (is(annotationFile, 'GRanges')){
    mark = annotationFile
  }else{
    mark = readPeakFile(annotationFile)
  }
  
  intersects = !is.na(findOverlaps(elements, mark, select = 'first'))
  
  intersected_elements = elements
  intersected_elements$intersect = intersects
  
  bg_cover_ratios = c()
  #draw background elements
  for (i in 1:numPermutations){
    bg = createRandomRegions(nregions = length(elements), length.mean = 200, length.sd = 0, genome = 'dm3', mask = blacklist)
    bg$intersects = !is.na(findOverlaps(bg, mark, select = 'first'))
    bg_cover_ratios = c(bg_cover_ratios, (sum(bg$intersects)/ length(bg)))
  }
  
  fg_cover_ratio = sum(intersected_elements$intersect)/ length(elements)
  
  p_val = sum(bg_cover_ratios > fg_cover_ratio)/ numPermutations
  
  fold_change = fg_cover_ratio/ mean(bg_cover_ratios)
  
  results = list('p_val' = p_val, 'cover_ratio' = fg_cover_ratio, 'fold_change' = fold_change, 'backround' = bg_cover_ratios)
  return(results)
}