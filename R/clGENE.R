#' This function internally for screening key genes.
#' @keywords internal
#' @param scores, n_top
#' @return Key similar genes
select_top_n<-function(scores,n_top){
  d <- data.frame(
    x   = data.table::copy(scores),
    indice=seq(1,length(scores)))
  
  data.table::setDT(d)
  data.table::setorder(d,-x)
  n_top_indice<-d$indice[1:n_top]
  return(n_top_indice)
}
###########################################
#' The inner function is used for L1 regularization with cosine similarity algorithm, screening of key genes.
#' @keywords internal
#' @param  scores,penalty_factor
#' @return Key similar genes
# L1 penalty function
penalty_function_L1 <- function(scores, penalty_factor) {
  penalized_scores <- sapply(scores, function(x) {
    sign(x) * max(0, abs(x) - penalty_factor)
  })
  return(penalized_scores)
}
#########################################
#' This function is used to for the normalized matrix
#' @param  X
#' @return normalize_matrix
normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}
##########################################
#' Calculate Cosine Similarity for Gene Expression Data
#'
#' This function calculates the cosine similarity for gene expression data
#' and performs L1 regularization.
#' @param data A matrix of gene expression data.
#' @param group_info A vector or factor of group classifications.
#' @return A data frame of cosine similarity scores.
#' @importFrom dplyr select filter
#' @examples
#' data <- matrix(rnorm(100), 10, 10)
#' group_info <- rep(1:2, each=5)
#' results <- ssGED(data, group_info)
ssGED <- function(data, group_info, penalty_factor = 0.1) {
  normalized_data <- t(apply(data, 1, normalize))
  genexcell <- normalized_data
  n_sample <- ncol(normalized_data)
  n_gene <- nrow(normalized_data)
  gene_name <- rownames(genexcell)
  
  groups_order <- sort(unique(group_info$group))
  n_cluster <- length(groups_order)
  
  cluster_mat <- matrix(0, nrow = n_cluster, ncol = n_sample)
  order_i <- 1
  for (group_i in groups_order) {
    idx_i <- group_info == group_i
    cluster_mat[order_i, idx_i] <- 1
    order_i <- order_i + 1
  }
  
  genexcell <- Matrix(as.matrix(genexcell), sparse = TRUE)
  cosine_sim <- proxyC::simil(genexcell, cluster_mat, method = "cosine", drop0 = TRUE)
  
  # Apply L1 penalty to cosine similarity
  penalized_cosine_sim <- apply(cosine_sim, 2, function(column) penalty_function_L1(column, penalty_factor))
  
  # Initializes the resulting data frame
  rank_stats_names <- data.frame(matrix(NA, n_gene, length(groups_order),
                                        dimnames = list(seq(1, n_gene), groups_order)),
                                 stringsAsFactors = FALSE)
  rank_stats_scores <- data.frame(matrix(NA, n_gene, length(groups_order),
                                         dimnames = list(seq(1, n_gene), groups_order)),
                                  stringsAsFactors = FALSE)
  
  order_i <- 1
  for (group_i in groups_order) {
    idx_i <- group_info == group_i
    scores <- penalized_cosine_sim[, order_i]
    
    # The select_top_n function was used for scoring
    global_indices <- select_top_n(scores, n_gene)
    rank_stats_names[, order_i] <- gene_name[global_indices]
    rank_stats_scores[, order_i] <- scores[global_indices]
    
    order_i <- order_i + 1
  }
  
  colnames(rank_stats_names) <- groups_order
  colnames(rank_stats_scores) <- groups_order
  
  ranks_stats <- list(names = rank_stats_names,
                      scores = rank_stats_scores)
  
  
  first_5 <- apply(ranks_stats$name, 2, head, n = 20)
  merged_column <- c(first_5)
  TOPexp<- ncol(normalized_data)
  TOPgene <-normalized_data[merged_column, (1:TOPexp)]
  library(pheatmap)
  
  # Standardize the top gene data using Z-scores
  z_scores <- t(scale(t(TOPgene)))
  
  # Limit the Z-scores to a maximum and minimum of 4 and -4
  z_scores_limited <- pmax(pmin(z_scores, 4), -4)
  
  # Create a heatmap of the limited Z-scores without clustering
  pheatmap(z_scores_limited,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           scale = "none",  # No scaling is needed as manual adjustment was done
           show_rownames = TRUE,
           show_colnames = FALSE,  # Column names are not shown based on the function input
           legend = TRUE,
           color = colorRampPalette(c("blue", "white", "red"))(255))  # Use a blue-white-red color gradient
  
  return(ranks_stats)
}
####################################
#' TOPGENE Screening
#'
#' This function uses cosine similarity to screen key genes.
#' @param matrix A numeric matrix representing gene expressions.
#' @param threshold A numeric value setting the threshold for gene selection.
#' @param genenumber An integer indicating the number of top genes to select.
#' @return A matrix containing the selected top genes.
#' @examples
#' data <- matrix(rnorm(100), 10, 10)
#' selected_genes <- TOPGENE(data, threshold = 0.5, genenumber = 5)
#' @export
TOPGENE <- function(matrix, threshold = 0.5, genenumber = 200) {
  if (genenumber <= 0) {
    stop("genenumber must be a positive integer")
  }
  # Assuming that matrices function is defined, used for data preprocessing
  normalized_data <- normalize(matrix)  # Suppose the normalize function is defined
  
  # The first row is extracted and removed from the matrix
  first_row <- normalized_data[1, , drop = FALSE]
  normalized_data <- normalized_data[-1, , drop = FALSE]
  
  # The first row was converted to a data frame for subsequent analysis
  column_df <- as.data.frame(t(first_row))
  
  # Assuming ssGED functions have been defined, return a list of all contain gene name
  ranks_stats <- ssGED(normalized_data, column_df)
  
  # If ranks_stats $name is a list of lists, appropriate to modify this line
  merged_column <- apply(ranks_stats$name, 2, head, n = genenumber)
  # The number of columns using normalized_data
  TOPexp <- ncol(normalized_data)
  TOPgene <- normalized_data[merged_column, 1:TOPexp, drop = FALSE]
  n_rows <- nrow(TOPgene)
  row_indices <- seq(1, n_rows, genenumber)
  list_of_matrices <- lapply(row_indices, function(idx) {
    end_row <- min(idx + genenumber - 1, n_rows)
    TOPgene[idx:end_row, ]
  })
  
  return(list_of_matrices)
}
####################################
#' The matrix was debatched using PCA
#' @param   matrix
#' @return Debatching the matrix
#' @export
PCA_adjust <- function(matrix) {
  matrix <- t(matrix)
  
  # Detect and remove often zero sequence or all of the columns
  constant_or_zero_cols <- apply(matrix, 2, function(col) {
    sd(col) == 0
  })
  matrix <- matrix[, !constant_or_zero_cols]
  
  if (ncol(matrix) %% 2 != 0) {
    matrix <- matrix[, -ncol(matrix)]
  }
  
  split_point <- ceiling(ncol(matrix) / 2)
  m1 <- matrix[, 1:split_point]  # Before the first matrix contains split_point columns
  m2 <- matrix[, (split_point + 1):ncol(matrix)]  # The second matrix contains the remaining columns
  
  # Perform PCA analysis
  pca_results  <- prcomp(matrix)
  pca_results1 <- prcomp(m1, scale. = TRUE)
  pca_results2 <- prcomp(m2, scale. = TRUE)
  
  result1 <- pca_results1$rotation
  result2 <- pca_results2$rotation
  cosine_similarity <- sapply(1:ncol(result1), function(i) {
    dot_product <- sum(result1[,i] * result2[,i])
    norm_m1 <- sqrt(sum(result1[,i]^2))
    norm_m2 <- sqrt(sum(result2[,i]^2))
    dot_product / (norm_m1 * norm_m2)
  })
  min_index <- which.min(cosine_similarity)
  
  # The maximum variance of difference was obtained
  PC1 <- pca_results$rotation[, min_index]
  
  # Adjust each column to match the first principal component of the standard deviation
  adjusted_matrix <- t(apply(matrix, 1, function(row) {
    standardized_row <- row / sd(row)
    adjusted_row <- standardized_row - mean(standardized_row)
    return(standardized_row)
  }))
  adjusted_matrix <- t(adjusted_matrix)
  return(adjusted_matrix)
}
########################
#' Orthogonal Transformation
#'
#' This function applies an orthogonal transformation to the given matrix to remove dependency.
#' @param A A numeric matrix.
#' @return A matrix that has been orthogonally transformed.
#' @examples
#' A <- matrix(rnorm(25), 5, 5)
#' ortho_A <- Orthogonal(A)
#' @export
Orthogonal <- function(A) {
  # Saving metadata
  row_names <- rownames(A)
  col_names <- colnames(A)
  
  # Gets the dimensions of the matrix
  n <- nrow(A)
  m <- ncol(A)
  Q <- matrix(0, n, m)
  
  non_zero_columns <- logical(m)  # Create a logical vector to record nonzero columns
  
  for (j in 1:m) {
    v <- A[, j]
    
    for (i in 1:(j - 1)) {
      v <- v - sum(Q[, i] * v) * Q[, i]
    }
    
    norm_v <- sqrt(sum(v^2))
    if (norm_v < 1e-10) {
      warning(paste("向量", j, "近似为零向量，被忽略"))
      next
    }
    
    Q[, j] <- v / norm_v
    non_zero_columns[j] <- TRUE  # Marked as a nonzero column
  }
  
  Q <- Q[, non_zero_columns]  # Keep only the non-zero columns
  
  # Set dimnames only for existing columns
  rownames(Q) <- row_names
  colnames(Q) <- col_names[non_zero_columns]  # Note changes here
  # Select the columns shared by A and Q
  common_columns <- intersect(colnames(A), colnames(Q))
  matrix1 <- Q[, colnames(Q) %in% common_columns]
  matrix2 <- A[, colnames(A) %in% common_columns]
  
  # Calculate and return the difference before and after the orthogonalization
  result_matrix <- matrix2 - matrix1
  
  return(result_matrix)
}
####################################
#' Key Program Identification
#'
#' This function identifies key molecular programs using NMF.
#' @param data A numeric matrix of molecular expression profiles.
#' @param Clustering An integer specifying the number of clusters to use.
#' @return A data frame with clustering results.
#' @examples
#' data <- matrix(rnorm(100), 10, 10)
#' key_programs <- key_Program(data, Clustering = 5)
#' @export
key_Program <- function(data, Clustering = 10){
  data <- data.frame(data)
  nmf_result <- nmf(data, rank = Clustering, nrun = 10, seed = 123)
  ##The clustering results were extracted
  clusters <- predict(nmf_result)
  annotation <- data.frame(Cluster = factor(clusters))
  annotation <- annotation[order(as.numeric(as.character(annotation$Cluster))), , drop=FALSE]
  #The original matrix rank
  sort_index <- match(row.names(annotation), colnames(data))
  # Sort the columns of a matrix by index
  sort_index <- sort_index[!is.na(sort_index)] 
  sort_index <- sort_index[sort_index <= ncol(data)]
  m_sorted <- data[, sort_index]
  #Convert a matrix to a data frame
  data_df <- as.data.frame(data)
  #Grouping information
  group_info <- data.frame(
    Column = colnames(data_df),
    Group = (annotation$Cluster)
  )
  #The mean value of each row was calculated by grouping
  group_means <- sapply(unique(group_info$Group), function(g) {
    cols <- group_info$Column[group_info$Group == g]
    rowMeans(data_df[, cols])
  })
  #Transpose the result to a matrix  
  group_means_matrix <- as.matrix(group_means)
  z_scores <- t(scale(t(group_means_matrix)))
  z_scores_limited <- pmax(pmin(z_scores, 4), -4)
  pheatmap(z_scores_limited,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           scale = "none",  # No scaling is needed as manual adjustment was done
           show_rownames = TRUE,
           show_colnames = FALSE,  # Column names are not shown based on the function input
           legend = TRUE,
           color = colorRampPalette(c("#2878B5", "white", "#C82423"))(255))  # Use a blue-white-red color gradient
  return(z_scores_limited)
}

