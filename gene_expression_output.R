
create_expression_summary <- function(
    seurat_obj,
    genes,
    first_meta,
    second_meta = "",
    assay = "RNA",
    layer = "counts"
) {
  
  library(Seurat)
  library(Matrix)
  
  expr <- GetAssayData(
    seurat_obj,
    assay = assay,
    layer = layer
  )
  
  genes <- intersect(genes, rownames(expr))
  expr <- expr[genes, , drop = FALSE]
  
  meta <- seurat_obj@meta.data
  
  use_second_meta <- !is.null(second_meta) && second_meta != ""
  
  # ----------------- preserve factor order -----------------
  
  if (!is.factor(meta[[first_meta]])) {
    meta[[first_meta]] <- factor(
      meta[[first_meta]],
      levels = unique(as.character(meta[[first_meta]]))
    )
  }
  
  if (use_second_meta) {
    
    if (!is.factor(meta[[second_meta]])) {
      meta[[second_meta]] <- factor(
        meta[[second_meta]],
        levels = unique(as.character(meta[[second_meta]]))
      )
    }
  }
  
  # ----------------- group ids -----------------
  
  if (use_second_meta) {
    
    group_ids <- paste(
      as.character(meta[[first_meta]]),
      as.character(meta[[second_meta]]),
      sep = "||"
    )
    
  } else {
    
    group_ids <- as.character(meta[[first_meta]])
  }
  
  # ----------------- ordered group dataframe -----------------
  
  if (use_second_meta) {
    
    group_df <- data.frame(
      group  = group_ids,
      first  = meta[[first_meta]],
      second = meta[[second_meta]],
      stringsAsFactors = FALSE
    )
    
    group_df <- unique(group_df)
    
    group_df <- group_df[
      order(group_df$first, group_df$second),
      ,
      drop = FALSE
    ]
    
  } else {
    
    group_df <- data.frame(
      group = group_ids,
      first = meta[[first_meta]],
      stringsAsFactors = FALSE
    )
    
    group_df <- unique(group_df)
    
    group_df <- group_df[
      order(group_df$first),
      ,
      drop = FALSE
    ]
  }
  
  unique_groups <- group_df$group
  
  # ----------------- result containers -----------------
  
  result_cols <- list()
  
  header_first  <- c()
  header_second <- c()
  header_value  <- c()
  
  # ----------------- calculate expression -----------------
  
  for (group in unique_groups) {
    
    cells_use <- rownames(meta)[group_ids == group]
    
    sub_expr <- expr[, cells_use, drop = FALSE]
    
    mean_vec <- Matrix::rowMeans(sub_expr)
    
    prop_vec <- Matrix::rowMeans(sub_expr != 0) * 100
    
    if (use_second_meta) {
      
      split_group <- strsplit(group, "\\|\\|")[[1]]
      
      current_first  <- split_group[1]
      current_second <- split_group[2]
      
    } else {
      
      current_first  <- group
      current_second <- NULL
    }
    
    header_first <- c(
      header_first,
      current_first,
      current_first,
      current_first
    )
    
    if (use_second_meta) {
      
      header_second <- c(
        header_second,
        current_second,
        current_second,
        current_second
      )
    }
    
    header_value <- c(
      header_value,
      "mean",
      "relative expression",
      "non-zero proportion (%)"
    )
    
    result_cols[[length(result_cols) + 1]] <- mean_vec
    
    result_cols[[length(result_cols) + 1]] <- mean_vec
    # placeholder for relative expression
    
    result_cols[[length(result_cols) + 1]] <- prop_vec
  }
  
  # ----------------- combine matrix -----------------
  
  mat <- do.call(cbind, result_cols)
  
  rownames(mat) <- genes
  
  # ----------------- relative expression -----------------
  
  n_groups <- length(unique_groups)
  
  mean_mat <- mat[, seq(1, n_groups * 3, by = 3), drop = FALSE]
  
  max_vals <- apply(mean_mat, 1, max)
  
  relative_mat <- mean_mat
  
  nonzero_idx <- max_vals > 0
  
  relative_mat[nonzero_idx, ] <-
    mean_mat[nonzero_idx, , drop = FALSE] /
    max_vals[nonzero_idx]
  
  relative_mat[!nonzero_idx, ] <- 0
  
  mat[, seq(2, ncol(mat), by = 3)] <- relative_mat
  
  # ----------------- column names -----------------
  
  if (use_second_meta) {
    
    colnames(mat) <- paste(
      header_first,
      header_second,
      header_value,
      sep = "|"
    )
    
  } else {
    
    colnames(mat) <- paste(
      header_first,
      header_value,
      sep = "|"
    )
  }
  
  # ----------------- export dataframe -----------------
  
  data_df <- as.data.frame(mat)
  
  if (use_second_meta) {
    
    header_df <- rbind(
      c(first_meta, header_first),
      c(second_meta, header_second),
      c("value", header_value)
    )
    
  } else {
    
    header_df <- rbind(
      c(first_meta, header_first),
      c("value", header_value)
    )
  }
  
  header_df <- as.data.frame(
    header_df,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  export_df <- data.frame(
    gene = genes,
    data_df,
    check.names = FALSE
  )
  
  colnames(header_df) <- colnames(export_df)
  
  export_df <- rbind(header_df, export_df)
  
  rownames(export_df) <- NULL
  
  colnames(export_df) <- export_df[1, ]
  
  export_df <- export_df[-1, ]
  
  return(list(
    data = data_df,
    header = header_df,
    export = export_df
  ))
}
# 
# res <- create_expression_summary(
#   seurat_obj = seurat_obj,
#   genes = rownames(seurat_obj),
#   first_meta = "sample_id",
#   second_meta = ""
# )
