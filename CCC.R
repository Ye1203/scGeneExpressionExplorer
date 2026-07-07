run_ccc <- function(seurat_file_path,
                    sample_info,
                    cluster_column,
                    species,
                    pathway_type,
                    prob_type,
                    min_cells,
                    save_path,
                    project_name,
                    runtime = 12,
                    cores = 8,
                    email = NA) {
  
  if (!dir.exists(save_path)) {
    dir.create(save_path, showWarnings = FALSE, recursive = TRUE)
  }
  
  sample_info_file <- file.path(save_path, "cellchat_sample_info.rds")
  saveRDS(sample_info, sample_info_file)
  
  params_file <- file.path(save_path, "cellchat_parameters.txt")
  
  sample_text <- paste(capture.output(str(sample_info)), collapse = "\n")
  
  params_text <- paste0(
    "=== CellChat Analysis Parameters ===\n",
    "Project Name: ", project_name, "\n",
    "Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
    "\n--- Data ---\n",
    "Seurat File: ", seurat_file_path, "\n",
    "Sample Info File: ", sample_info_file, "\n",
    "Sample Setting:\n", sample_text, "\n",
    "Cluster Column: ", cluster_column, "\n",
    "\n--- CellChat Parameters ---\n",
    "Species: ", species, "\n",
    "Pathway Types: ", paste(pathway_type, collapse = "; "), "\n",
    "Probability Method: ", prob_type, "\n",
    "Min Cells Filter: ", min_cells, "\n",
    "\n--- HPC Parameters ---\n",
    "CPU Cores: ", cores, "\n",
    "Runtime (hours): ", runtime, "\n",
    "Email: ", if (is.na(email) || !nzchar(email)) "N/A" else email, "\n"
  )
  
  writeLines(params_text, params_file)
  
  seurat_filename <- basename(seurat_file_path)
  seurat_copy_path <- file.path(save_path, seurat_filename)
  file.copy(seurat_file_path, seurat_copy_path, overwrite = TRUE)
  
  qsub_file <- file.path(save_path, "cellchat_analysis.qsub")
  
  qsub_content <- generate_qsub_script(
    project_name = project_name,
    cores = cores,
    runtime = runtime,
    email = email,
    save_path = save_path,
    seurat_file = seurat_copy_path,
    sample_info_file = sample_info_file,
    cluster_column = cluster_column,
    species = species,
    pathway_type = pathway_type,
    prob_type = prob_type,
    min_cells = min_cells
  )
  
  writeLines(qsub_content, qsub_file)
  
  mission_id <- submit_qsub_job(qsub_file)
  
  return(mission_id)
}


generate_qsub_script <- function(project_name,
                                 cores,
                                 runtime,
                                 email,
                                 save_path,
                                 seurat_file,
                                 sample_info_file,
                                 cluster_column,
                                 species,
                                 pathway_type,
                                 prob_type,
                                 min_cells) {
  
  hours <- as.integer(runtime)
  minutes <- as.integer((runtime - hours) * 60)
  time_format <- sprintf("%02d:%02d:00", hours, minutes)
  
  pathway_str <- paste(pathway_type, collapse = "|")
  
  script_dir <- getwd()
  cellchat_runner <- file.path(script_dir, "cellchat_runner.R")
  log_path <- file.path(save_path, "ccc.log")
  
  email_str <- if (is.na(email) || !nzchar(email)) "NA" else email
  
  qsub_script <- paste0(
    "#!/bin/bash\n",
    "#$ -N ccc\n",
    "#$ -cwd\n",
    "#$ -j y\n",
    "#$ -o ", shQuote(log_path), "\n",
    "#$ -pe omp ", cores, "\n",
    "#$ -l h_rt=", time_format, "\n",
    "#$ -l mem_per_core=16G\n",
    "#$ -V\n",
    "#$ -P ", project_name, "\n",
    "\n",
    "echo \"==========================================================\"\n",
    "echo \"Starting on       : $(date)\"\n",
    "echo \"Running on node   : $(hostname)\"\n",
    "echo \"Current job ID    : $JOB_ID\"\n",
    "echo \"Current job name  : $JOB_NAME\"\n",
    "echo \"==========================================================\"\n",
    "\n",
    "module load R/4.4.3\n",
    "\n",
    "Rscript ", shQuote(cellchat_runner), " \\\n",
    "  --seurat_file=", shQuote(seurat_file), " \\\n",
    "  --sample_info_file=", shQuote(sample_info_file), " \\\n",
    "  --cluster_column=", shQuote(cluster_column), " \\\n",
    "  --species=", shQuote(species), " \\\n",
    "  --pathway_type=", shQuote(pathway_str), " \\\n",
    "  --prob_type=", shQuote(prob_type), " \\\n",
    "  --min_cells=", min_cells, " \\\n",
    "  --cores=", cores, " \\\n",
    "  --save_path=", shQuote(save_path), " \\\n",
    "  --email=", shQuote(email_str), "\n"
  )
  
  return(qsub_script)
}


submit_qsub_job <- function(qsub_file) {
  
  Sys.chmod(qsub_file, mode = "0755")
  
  result <- system2("qsub", args = qsub_file, stdout = TRUE, stderr = TRUE)
  
  mission_id <- sub(".*?([0-9]+).*", "\\1", result)
  
  if (is.na(mission_id) || mission_id == "" || identical(mission_id, result)) {
    stop("Failed to submit qsub job. Result: ", paste(result, collapse = "\n"))
  }
  
  return(mission_id)
}

diagnose_cellchat_summary_method <- function(
    seurat_obj,
    group.by,
    assay = "RNA",
    layer = "data",
    cellchat_db = NULL,
    trim_values = c(0.1, 0.05),
    min_positive_value = 0,
    min_pct_cutoffs = c(0.05, 0.1, 0.25)
) {
  
  library(Seurat)
  library(Matrix)
  library(dplyr)
  
  expr <- GetAssayData(
    seurat_obj,
    assay = assay,
    layer = layer
  )
  
  meta <- seurat_obj@meta.data
  
  if (!group.by %in% colnames(meta)) {
    stop("group.by not found in seurat_obj@meta.data")
  }
  
  groups <- as.factor(meta[[group.by]])
  names(groups) <- rownames(meta)
  
  common_cells <- intersect(colnames(expr), names(groups))
  expr <- expr[, common_cells, drop = FALSE]
  groups <- groups[common_cells]
  
  # ------------------------------------------------------------
  # 1. Get CellChat ligand/receptor genes
  # ------------------------------------------------------------
  if (is.character(cellchat_db)) {
    cellchat_db <- if (cellchat_db == "human") {
      CellChatDB.human
    } else {
      CellChatDB.mouse
    }
  }
  
  if (!is.null(cellchat_db)) {
    
    if (!"interaction" %in% names(cellchat_db)) {
      stop("cellchat_db should be like CellChatDB.mouse or CellChatDB.human")
    }
    
    db <- cellchat_db$interaction
    
    lr_genes <- unique(c(
      db$ligand,
      db$receptor
    ))
    
    lr_genes <- unique(unlist(strsplit(lr_genes, "_")))
    lr_genes <- lr_genes[lr_genes %in% rownames(expr)]
    
    gene_pathway_type <- data.frame(
      gene = character(),
      annotation = character(),
      stringsAsFactors = FALSE
    )
    
    for (i in seq_len(nrow(db))) {
      
      genes_i <- unique(unlist(strsplit(
        c(db$ligand[i], db$receptor[i]),
        "_"
      )))
      
      genes_i <- genes_i[genes_i %in% lr_genes]
      
      if (length(genes_i) > 0) {
        gene_pathway_type <- rbind(
          gene_pathway_type,
          data.frame(
            gene = genes_i,
            annotation = db$annotation[i],
            stringsAsFactors = FALSE
          )
        )
      }
    }
    
    gene_pathway_type <- unique(gene_pathway_type)
    
  } else {
    
    lr_genes <- rownames(expr)
    
    gene_pathway_type <- data.frame(
      gene = lr_genes,
      annotation = "All genes",
      stringsAsFactors = FALSE
    )
  }
  
  expr <- expr[lr_genes, , drop = FALSE]
  
  # ------------------------------------------------------------
  # 2. Summary functions
  # ------------------------------------------------------------
  tri_mean <- function(x) {
    qs <- quantile(x, probs = c(0.25, 0.5, 0.75), names = FALSE)
    (qs[1] + 2 * qs[2] + qs[3]) / 4
  }
  
  truncated_mean <- function(x, trim) {
    mean(x, trim = trim)
  }
  
  # ------------------------------------------------------------
  # 3. Compute gene-cluster summaries
  # ------------------------------------------------------------
  result_list <- list()
  
  for (cl in levels(groups)) {
    
    cells_use <- names(groups)[groups == cl]
    mat <- expr[, cells_use, drop = FALSE]
    
    pct_exp <- Matrix::rowMeans(mat > min_positive_value)
    mean_exp <- Matrix::rowMeans(mat)
    
    mat_dense <- as.matrix(mat)
    
    median_exp <- apply(mat_dense, 1, median)
    trimean_exp <- apply(mat_dense, 1, tri_mean)
    
    df <- data.frame(
      gene = rownames(mat),
      cluster = cl,
      n_cells = length(cells_use),
      pct_exp = pct_exp,
      mean = mean_exp,
      median = median_exp,
      triMean = trimean_exp,
      stringsAsFactors = FALSE
    )
    
    for (tr in trim_values) {
      colname <- paste0("truncatedMean_trim_", tr)
      df[[colname]] <- apply(
        mat_dense,
        1,
        truncated_mean,
        trim = tr
      )
    }
    
    result_list[[cl]] <- df
  }
  
  gene_cluster_summary <- bind_rows(result_list)
  
  # ------------------------------------------------------------
  # 4. Method-level signal retention summary
  # ------------------------------------------------------------
  
  method_cols <- c(
    "median",
    "triMean",
    paste0("truncatedMean_trim_", trim_values)
  )
  
  display_method_name <- function(m) {
    if (m == "median") return("median")
    if (m == "triMean") return("triMean")
    
    if (grepl("^truncatedMean_trim_", m)) {
      trim_value <- sub("^truncatedMean_trim_", "", m)
      return(paste0("truncatedMean(", trim_value, ")"))
    }
    
    m
  }
  
  gene_cluster_summary_with_type <- gene_cluster_summary %>%
    left_join(
      gene_pathway_type,
      by = "gene"
    )
  
  pathway_types <- c("Secreted Signaling", "Cell-Cell Contact", "ECM-Receptor", "Non-protein Signaling")
  
  method_summary <- lapply(method_cols, function(m) {
    
    total_n <- sum(gene_cluster_summary[[m]] > 0)
    
    row_out <- data.frame(
      `Probability Method` = display_method_name(m),
      Total = total_n,
      check.names = FALSE
    )
    
    for (pt in pathway_types) {
      
      n_pt <- gene_cluster_summary_with_type %>%
        filter(annotation == pt) %>%
        summarise(n = sum(.data[[m]] > 0)) %>%
        pull(n)
      
      row_out[[pt]] <- n_pt
    }
    
    row_out
    
  }) %>%
    bind_rows()
  # ------------------------------------------------------------
  # 5. Expression sparsity summary
  # ------------------------------------------------------------
  pct_summary <- data.frame(
    cutoff = min_pct_cutoffs,
    gene_cluster_pct_above_cutoff = sapply(
      min_pct_cutoffs,
      function(cut) mean(gene_cluster_summary$pct_exp >= cut)
    )
  )
  
  median_pct <- median(gene_cluster_summary$pct_exp)
  
  pct_10_25 <- mean(
    gene_cluster_summary$pct_exp >= 0.1 &
      gene_cluster_summary$pct_exp < 0.25
  )
  
  pct_above_25 <- mean(gene_cluster_summary$pct_exp >= 0.25)
  pct_below_10 <- mean(gene_cluster_summary$pct_exp < 0.1)
  
  # ------------------------------------------------------------
  # 6. Simplified recommendation logic
  # ------------------------------------------------------------
  get_pos_n <- function(method_name) {
    method_summary$Total[
      method_summary$`Probability Method` == method_name
    ]
  }
  safe_ratio <- function(a, b) {
    if (
      length(a) == 0 ||
      length(b) == 0 ||
      is.na(a) ||
      is.na(b) ||
      b == 0
    ) {
      return(NA_real_)
    }
    a / b
  }
  
  trimean_n <- get_pos_n("triMean")
  trim005_n <- get_pos_n("truncatedMean(0.05)")
  trim010_n <- get_pos_n("truncatedMean(0.1)")
  
  ratio_trim010_vs_trimean <- safe_ratio(trim010_n, trimean_n)
  ratio_trim005_vs_trim010 <- safe_ratio(trim005_n, trim010_n)
  
  if (is.na(ratio_trim010_vs_trimean)) {
    
    final_recommendation <- paste0(
      "Unable to compare triMean and truncatedMean(0.1). ",
      "Please check whether the method summary contains valid nonzero signals."
    )
    
  } else if (ratio_trim010_vs_trimean < 1.5) {
    
    final_recommendation <- paste0(
      "Use triMean as the main CellChat setting. ",
      "truncatedMean(0.1) does not add many extra gene-cluster signals."
    )
    
  } else if (ratio_trim010_vs_trimean < 3) {
    
    final_recommendation <- paste0(
      "Use triMean as the conservative main result and truncatedMean(0.1) ",
      "as sensitivity analysis. truncatedMean(0.1) adds a substantial number ",
      "of extra signals."
    )
    
  } else {
    
    final_recommendation <- paste0(
      "The result is highly sensitive to the summary method. ",
      "Use triMean as a conservative result, and inspect cluster heterogeneity ",
      "before interpreting truncatedMean(0.1) results."
    )
  }
  
  if (!is.na(ratio_trim005_vs_trim010)) {
    
    if (ratio_trim005_vs_trim010 < 1.3) {
      trim005_comment <- "truncatedMean(0.05) is similar to truncatedMean(0.1)."
    } else if (ratio_trim005_vs_trim010 < 1.8) {
      trim005_comment <- paste0(
        "truncatedMean(0.05) adds additional rare-cell signals compared with ",
        "truncatedMean(0.1); use it cautiously."
      )
    } else {
      trim005_comment <- paste0(
        "Avoid using truncatedMean(0.05) as the main result because it adds many ",
        "extra rare-cell signals compared with truncatedMean(0.1)."
      )
    }
    
  } else {
    trim005_comment <- "Unable to evaluate truncatedMean(0.05)."
  }
  
  if (
    !is.na(ratio_trim010_vs_trimean) &&
    ratio_trim010_vs_trimean < 1.5
  ) {
    
    recommended_method <- "triMean"
    
  } else if (
    !is.na(ratio_trim010_vs_trimean) &&
    ratio_trim010_vs_trimean >= 3 &&
    !is.na(ratio_trim005_vs_trim010) &&
    ratio_trim005_vs_trim010 < 1.3
  ) {
    
    recommended_method <- "truncatedMean(0.05)"
    
  } else {
    
    recommended_method <- "truncatedMean(0.1)"
  }
  
  recommendation <- list(
    recommended_method = recommended_method,
    final_recommendation = final_recommendation,
    trim005_comment = trim005_comment,
    diagnostic_ratios = data.frame(
      `truncatedMean(0.1) / triMean` = round(ratio_trim010_vs_trimean, 3),
      `truncatedMean(0.05) / truncatedMean(0.1)` = round(ratio_trim005_vs_trim010, 3),
      check.names = FALSE
    )
  )
  
  return(list(
    gene_cluster_summary = gene_cluster_summary,
    method_summary = method_summary,
    pct_summary = pct_summary,
    sparsity_summary = data.frame(
      median_pct_exp = median_pct,
      pct_gene_cluster_below_10pct = pct_below_10,
      pct_gene_cluster_10_to_25pct = pct_10_25,
      pct_gene_cluster_above_25pct = pct_above_25
    ),
    recommendation = recommendation
  ))
}

