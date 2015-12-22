compare_overlap_between_two_genesets = function ( geneset1, geneset2, background_list, granularity=10,
              color = colorRampPalette(c("#3362A5", "dodgerblue1", "deepskyblue", "gray99", "gold", "firebrick2","#A31D1D"), space="Lab")) { 
  # This function finds and plots the degree of overlap between two sets of genes
  # This function outputs a filled contour graph of the relationship
  # The return value is a list that contains the odds ratio and the associated p-values
  # There two optional parameters: granularity controls the number of genes to consider at each step
  # Additionally, a color palette can be specified which will be passed on to the filled.contour
  # The input is two indices of the enriched regions table
  # The degree of overlap is quantified using Fisher's exact test
  stopifnot(is.character(geneset1), is.character(geneset2), is.character(background_list), is.numeric(granularity))
 
  ## TEST THE REGION BELOW
  total_overlaps = length(background_list )
  fmat = matrix ( nrow = 2, ncol = 2)
  enrichment_matrix_odds = matrix(nrow = floor( length(geneset1) / granularity ), 
                                  ncol = floor( length(geneset2) / granularity ))
  enrichment_matrix_pval = matrix(nrow = floor( length(geneset1) / granularity ),
                                  ncol = floor(length(geneset2) / granularity ))
  
  
  for (k in 1: floor(length(geneset1) / granularity )) { 
    for (l in 1: floor(length(geneset2) / granularity )) {
      in_both = length ( intersect (geneset1[1:(k*granularity)], geneset2[1:(l*granularity)] ) )
      in_one = length ( setdiff(geneset1[1:(k*granularity)], geneset2[1:(l*granularity)]) )
      in_two = length ( setdiff(geneset2[1:(l*granularity)], geneset1[1:(k*granularity)]) )
      neither = total_overlaps - length(union(geneset1[1:(k*granularity)], geneset2[1:(l*granularity)]))
      fmat[1,]  = c(in_both, in_one)
      fmat[2,] = c(in_two, neither)
      fisher_test = fisher.test(fmat)
      enrichment_matrix_odds[k,l] = fisher_test$estimate
      enrichment_matrix_pval[k,l] = fisher_test$p.value
      
    }
  }
  filled.contour(enrichment_matrix_odds, color.palette = color, 
                 zlim = c(0, max(enrichment_matrix_odds)), 
                 main = "Odds Ratio of Enrichment" )
  return ( list(ODDS = enrichment_matrix_odds, PVAL= enrichment_matrix_pval) )  
}



