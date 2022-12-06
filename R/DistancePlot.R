#' Jaccard Similarity C++ Function
#' 
#' @details Returns the Jaccard Similarity index of a selected gene with each of the genes in w
#'
#' @param w         matrix giving genes as rows and samples as columns
#' @param YFG_idx   index of selected gene
#' 
#' @return res      jaccard similarity
#' 
#' @export
#'  
Rcpp::cppFunction(
     'Rcpp::NumericVector c_jaccard_similarity(const Rcpp::NumericMatrix w, size_t YFG_idx) {
       --YFG_idx; // C++ indexing is 1 less than R indexing
       Rcpp::NumericVector res(w.rows());
       Rcpp::NumericVector p = w.row(YFG_idx);
       for(size_t i = 0; i < w.rows(); ++i){
         Rcpp::NumericVector q = w.row(i);
         double pq  = 0, p2 = 0, q2 = 0;
         for(size_t j = 0; j < w.cols(); ++j){
           pq += p[j] * q[j];
           p2 += p[j] * p[j];
           q2 += q[j] * q[j];
         }
         // divide intersection by union (p^2 + q^2 - intersection)
         res(i) = pq / (p2 + q2 - pq);     
       }
       return res;
     }')
     
# Creating weighted matrix w_0 based on celltype - should this part be in this file or just the function?
w <- data@reductions$nmf@feature.loadings
h <- data@reductions$nmf@cell.embeddings
celltype <- data@meta.data$cell_type

h_0 <- h[which(celltype == "Naive CD4 T"),]
column_mean <- colMeans(h_0)
library(Matrix)
w_0 <- as.matrix(w %*% Diagonal(x = column_mean))

#' Plot overall distance vs cell type distance
#'
#' @details returns a ggplot2 with the overall distance vs cell type distance 
#'  
#' @param w     matrix giving genes as rows and samples as columns
#' @param w_0   weighted matrix giving genes as rows and samples as columns
#' @param YFG   gene to compare
#' 
#' @return      a ggplot2 object
#' 
#' @export 
#' 
plot_YFG <- function(w, w_0, YFG = NULL){
   
  dist_YFG <- rep(0, nrow(w)) 
  sapply(YFG, function(gene) {
    idx_YFG <<- which(rownames(w) == gene)
    dist_YFG <<- dist_YFG + c_jaccard_similarity(w, idx_YFG)
  })
  
  dist_YFG <- dist_YFG/nrow(w)
   
  dist_w0 <- rep(0, nrow(w_0))
  #for(gene in YFG){
  #  idx_w0 <- which(rownames(w_0) == gene)
  #  dist_w0 <- dist_w0 + c_jaccard_similarity(w_0, idx_w0)
  #}
  sapply(YFG, function(gene) {
    idx_w0 <<- which(rownames(w_0) == gene)
    dist_w0 <<- dist_w0 + c_jaccard_similarity(w_0, idx_w0)
  })
  dist_w0 <- dist_w0/nrow(w_0)
  
  library(ggplot2)
  df <- data.frame("w" = dist_YFG[-idx_YFG], "w_0" = dist_w0[-idx_w0])
  df <- df[order(df$w, decreasing = TRUE), ]
  df <- df[1:10, ]
  ggplot(df, aes(x = w, y = w_0)) + geom_point()+ 
   labs(title = "Jaccard Distance Between Weighted and Unweighted Genes", x = "Unweighted (w)", y = "Weighted (w_0)") +
   coord_cartesian(xlim = NULL, ylim = NULL) + 
   theme(axis.text.x = element_text(angle = 0, vjust = 0.9, hjust = 0.5))
}
