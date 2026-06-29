#!/usr/bin/env Rscript

project_root <- "/projectnb/wax-es/00_shinyapp/scGeneExpressionExplorer"

send_email <- function(to, subject, body) {
  if (!is.na(to) && nzchar(to)) {
    f <- tempfile(fileext = ".txt")
    mail_content <- c(
      paste0("To: ", to),
      paste0("Subject: ", subject),
      "",
      body
    )
    writeLines(mail_content, f)
    cmd <- paste0("/usr/sbin/sendmail -t < ", shQuote(f))
    ret <- system(cmd)
    unlink(f)
    if (ret != 0) {
      warning("sendmail command failed with exit code ", ret)
    }
  }
}


save_cellchat_tables_to_excel <- function(cellChat_obj, save_path, sample_name = NULL) {
  tryCatch({
    if (!require("openxlsx", quietly = TRUE)) {
      install.packages("openxlsx", repos = "https://cloud.r-project.org")
      library(openxlsx)
    }
    
    net_data <- tryCatch({
      df <- subsetCommunication(cellChat_obj, slot.name = "net")
      if (nrow(df) > 0) {
        df[order(df$prob, decreasing = TRUE), ]
      } else {
        data.frame(Message = "No net data available")
      }
    }, error = function(e) {
      data.frame(Error = paste("Failed to extract net data:", conditionMessage(e)))
    })
    
    netP_data <- tryCatch({
      df <- subsetCommunication(cellChat_obj, slot.name = "netP")
      if (nrow(df) > 0) {
        df[order(df$prob, decreasing = TRUE), ]
      } else {
        data.frame(Message = "No netP data available")
      }
    }, error = function(e) {
      data.frame(Error = paste("Failed to extract netP data:", conditionMessage(e)))
    })
    
    return(list(net = net_data, netP = netP_data))
    
  }, error = function(e) {
    return(list(
      net = data.frame(Error = "Failed to extract net data"),
      netP = data.frame(Error = "Failed to extract netP data")
    ))
  })
}


create_cellchat_excel <- function(cellChat_list, save_path, single_analysis = FALSE) {
  tryCatch({
    if (!require("openxlsx", quietly = TRUE)) {
      install.packages("openxlsx", repos = "https://cloud.r-project.org")
      library(openxlsx)
    }
    
    wb <- createWorkbook()
    
    if (single_analysis) {
      tables <- save_cellchat_tables_to_excel(cellChat_list, save_path)
      
      addWorksheet(wb, "net")
      writeData(wb, "net", tables$net, rowNames = FALSE)
      
      addWorksheet(wb, "netP")
      writeData(wb, "netP", tables$netP, rowNames = FALSE)
      
    } else {
      for (sample_name in names(cellChat_list)) {
        cellChat_obj <- cellChat_list[[sample_name]]
        tables <- save_cellchat_tables_to_excel(cellChat_obj, save_path, sample_name)
        
        net_sheet_name <- paste0(substr(sample_name, 1, 18), "_net")
        netP_sheet_name <- paste0(substr(sample_name, 1, 17), "_netP")
        
        addWorksheet(wb, net_sheet_name)
        writeData(wb, net_sheet_name, tables$net, rowNames = FALSE)
        
        addWorksheet(wb, netP_sheet_name)
        writeData(wb, netP_sheet_name, tables$netP, rowNames = FALSE)
      }
    }
    
    excel_file <- file.path(save_path, "cellchat_tables.xlsx")
    saveWorkbook(wb, excel_file, overwrite = TRUE)
    return(excel_file)
    
  }, error = function(e) {
    cat("[CellChat Runner] ERROR: Failed to create Excel file -", conditionMessage(e), "\n")
    return(NULL)
  })
}


parse_args <- function(args) {
  result <- list()
  
  for (arg in args) {
    if (grepl("^--", arg)) {
      key <- sub("^--", "", arg)
      key <- sub("=.*$", "", key)
      value <- sub("^--[^=]+=", "", arg)
      result[[key]] <- value
    }
  }
  
  return(result)
}


filter_cellchat_db <- function(CellChatDB, pathway_types) {
  db_list <- lapply(pathway_types, function(ptype) {
    tryCatch({
      subsetDB(CellChatDB, search = ptype)
    }, error = function(e) {
      cat("[CellChat Runner] Warning: Could not filter pathway:", ptype, "\n")
      return(CellChatDB)
    })
  })
  
  interactions_combined <- do.call(rbind, lapply(db_list, function(db) db$interaction))
  interactions_combined <- interactions_combined[!duplicated(interactions_combined), ]
  
  CellChatDB$interaction <- interactions_combined
  
  return(CellChatDB)
}


get_prob_setting <- function(prob_type) {
  if (tolower(prob_type) == "trimean") {
    return(list(method = "triMean", trim = 0.1))
  }
  
  if (tolower(prob_type) == "truncatedmean_0.1") {
    return(list(method = "truncatedMean", trim = 0.1))
  }
  
  if (tolower(prob_type) == "truncatedmean_0.05") {
    return(list(method = "truncatedMean", trim = 0.05))
  }
  
  warning("Unknown prob_type: ", prob_type, ". Falling back to triMean.")
  return(list(method = "triMean", trim = 0.1))
}


run_cellchat_core <- function(seu,
                              sample_name,
                              cluster_col,
                              species,
                              pathway_types,
                              prob_method,
                              prob_trim,
                              min_cells,
                              seurat_full = NULL) {
  
  cat("[CellChat Runner] Processing:", sample_name, "\n")
  cat("[CellChat Runner] Cells:", ncol(seu), "\n")
  cat("[CellChat Runner] Clusters:", paste(unique(seu@meta.data[[cluster_col]]), collapse = ", "), "\n")
  
  cellChat <- createCellChat(object = seu, group.by = cluster_col, assay = "RNA")
  
  if (species == "human") {
    CellChatDB <- CellChatDB.human
  } else {
    CellChatDB <- CellChatDB.mouse
  }
  
  CellChatDB <- filter_cellchat_db(CellChatDB, pathway_types)
  cellChat@DB <- CellChatDB
  
  cellChat <- subsetData(cellChat)
  
  cat("[CellChat Runner] Identifying over-expressed genes\n")
  cellChat <- identifyOverExpressedGenes(cellChat, do.fast = FALSE)
  
  cat("[CellChat Runner] Identifying over-expressed interactions\n")
  cellChat <- identifyOverExpressedInteractions(cellChat)
  
  cat("[CellChat Runner] Computing communication probability\n")
  cellChat <- computeCommunProb(cellChat, type = prob_method, trim = prob_trim)
  
  cat("[CellChat Runner] Filtering communication\n")
  cellChat <- filterCommunication(cellChat, min.cells = min_cells)
  
  cat("[CellChat Runner] Computing pathway-level probability\n")
  cellChat <- computeCommunProbPathway(cellChat)
  
  cat("[CellChat Runner] Aggregating network\n")
  cellChat <- aggregateNet(cellChat)
  
  cat("[CellChat Runner] Computing centrality\n")
  cellChat <- netAnalysis_computeCentrality(cellChat, slot.name = "netP")
  
  if (!is.null(seurat_full)) {
    cat("[CellChat Runner] Adding dimensionality reduction\n")
    cellChat <- addReduction(cellChat, seu.obj = seurat_full)
  }
  
  return(cellChat)
}


args <- commandArgs(trailingOnly = TRUE)
params <- parse_args(args)

SEURAT_FILE <- params$seurat_file
SAMPLE_INFO_FILE <- params$sample_info_file
CLUSTER_COLUMN <- params$cluster_column
SPECIES <- params$species
PATHWAY_TYPES <- strsplit(params$pathway_type, "\\|")[[1]]
PROB_TYPE <- params$prob_type
MIN_CELLS <- as.numeric(params$min_cells)
CORES <- as.numeric(params$cores)
SAVE_PATH <- params$save_path
EMAIL <- if (is.null(params$email) || params$email == "NA") NA else params$email

if (is.null(SAMPLE_INFO_FILE) || !file.exists(SAMPLE_INFO_FILE)) {
  stop("sample_info_file is missing or does not exist.")
}

SAMPLE_INFO <- readRDS(SAMPLE_INFO_FILE)
SAMPLE_COLUMN <- SAMPLE_INFO$sample_column

prob_setting <- get_prob_setting(PROB_TYPE)
prob_method <- prob_setting$method
prob_trim <- prob_setting$trim


tryCatch({
  
  if (file.exists(file.path(project_root, "renv/activate.R"))) {
    source(file.path(project_root, "renv/activate.R"))
    cat("[CellChat Runner] renv environment activated\n")
  } else {
    cat("[CellChat Runner] Warning: renv/activate.R not found\n")
  }
  
  cat("[CellChat Runner] Starting analysis...\n")
  cat("[CellChat Runner] Seurat file:", SEURAT_FILE, "\n")
  cat("[CellChat Runner] Sample column:", SAMPLE_COLUMN, "\n")
  cat("[CellChat Runner] Cluster column:", CLUSTER_COLUMN, "\n")
  cat("[CellChat Runner] Species:", SPECIES, "\n")
  cat("[CellChat Runner] Cores:", CORES, "\n")
  cat("[CellChat Runner] Email:", ifelse(is.na(EMAIL), "N/A", EMAIL), "\n")
  
  library(CellChat)
  library(Seurat)
  library(future)
  
  options(parallelly.maxWorkers.localhost = Inf)
  options(future.globals.maxSize = 10 * 1024^3)
  
  if (CORES > 1) {
    future::plan("multisession", workers = CORES)
    cat("[CellChat Runner] Using multisession with", CORES, "workers\n")
  } else {
    future::plan("sequential")
    cat("[CellChat Runner] Using sequential execution\n")
  }
  
  cat("[CellChat Runner] Loading Seurat object...\n")
  seurat_obj <- readRDS(SEURAT_FILE)
  
  if (!CLUSTER_COLUMN %in% colnames(seurat_obj@meta.data)) {
    stop("Column '", CLUSTER_COLUMN, "' not found in Seurat object metadata")
  }
  
  if (any(as.character(seurat_obj@meta.data[[CLUSTER_COLUMN]]) == "0", na.rm = TRUE)) {
    cat("[CellChat Runner] Adding 'C' prefix to cluster values\n")
    seurat_obj@meta.data$cellchat_cluster <- paste0("C", seurat_obj@meta.data[[CLUSTER_COLUMN]])
    CLUSTER_COLUMN <- "cellchat_cluster"
  }
  
  cat("[CellChat Runner] Final cluster column:", CLUSTER_COLUMN, "\n")
  cat("[CellChat Runner] Unique clusters:", paste(unique(seurat_obj@meta.data[[CLUSTER_COLUMN]]), collapse = ", "), "\n")
  
  RUN_BY_GROUP <- tolower(SAMPLE_COLUMN) != "all together"
  
  if (RUN_BY_GROUP) {
    
    if (!SAMPLE_COLUMN %in% colnames(seurat_obj@meta.data)) {
      stop("Column '", SAMPLE_COLUMN, "' not found in Seurat object metadata")
    }
    
    control_samples <- SAMPLE_INFO$control
    treatment_samples <- SAMPLE_INFO$treatment
    
    if (length(control_samples) == 0 || length(treatment_samples) == 0) {
      stop("CONTROL and TREATMENT groups must both contain at least one sample.")
    }
    
    cat("[CellChat Runner] Running CONTROL vs TREATMENT analysis\n")
    cat("[CellChat Runner] CONTROL samples:", paste(control_samples, collapse = ", "), "\n")
    cat("[CellChat Runner] TREATMENT samples:", paste(treatment_samples, collapse = ", "), "\n")
    
    meta_sample <- as.character(seurat_obj@meta.data[[SAMPLE_COLUMN]])
    names(meta_sample) <- rownames(seurat_obj@meta.data)
    
    control_cells <- names(meta_sample)[meta_sample %in% control_samples]
    treatment_cells <- names(meta_sample)[meta_sample %in% treatment_samples]
    
    if (length(control_cells) == 0) {
      stop("No cells found for CONTROL samples.")
    }
    
    if (length(treatment_cells) == 0) {
      stop("No cells found for TREATMENT samples.")
    }
    
    seu_control <- subset(seurat_obj, cells = control_cells)
    seu_treatment <- subset(seurat_obj, cells = treatment_cells)
    
    cellChat_list <- list()
    
    cellChat_list[["CONTROL"]] <- run_cellchat_core(
      seu = seu_control,
      sample_name = "CONTROL",
      cluster_col = CLUSTER_COLUMN,
      species = SPECIES,
      pathway_types = PATHWAY_TYPES,
      prob_method = prob_method,
      prob_trim = prob_trim,
      min_cells = MIN_CELLS,
      seurat_full = seurat_obj
    )
    
    cellChat_list[["TREATMENT"]] <- run_cellchat_core(
      seu = seu_treatment,
      sample_name = "TREATMENT",
      cluster_col = CLUSTER_COLUMN,
      species = SPECIES,
      pathway_types = PATHWAY_TYPES,
      prob_method = prob_method,
      prob_trim = prob_trim,
      min_cells = MIN_CELLS,
      seurat_full = seurat_obj
    )
    
    cat("[CellChat Runner] Saving CellChat list\n")
    cellchat_list_file <- file.path(SAVE_PATH, "cellchat_compare.rds")
    saveRDS(cellChat_list, cellchat_list_file)
    
    cat("[CellChat Runner] Generating Excel file\n")
    create_cellchat_excel(cellChat_list, SAVE_PATH, single_analysis = FALSE)
    
    summary_text <- paste0(
      "=== CellChat Analysis Summary (CONTROL vs TREATMENT) ===\n",
      "Analysis Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
      "Sample Column: ", SAMPLE_COLUMN, "\n",
      "Cluster Column: ", CLUSTER_COLUMN, "\n",
      "Species: ", SPECIES, "\n",
      "Pathway Types: ", paste(PATHWAY_TYPES, collapse = ", "), "\n",
      "Probability Method: ", PROB_TYPE, "\n",
      "Min Cells: ", MIN_CELLS, "\n",
      "\n--- Group Information ---\n",
      "CONTROL samples: ", paste(control_samples, collapse = ", "), "\n",
      "TREATMENT samples: ", paste(treatment_samples, collapse = ", "), "\n",
      "CONTROL cells: ", ncol(seu_control), "\n",
      "TREATMENT cells: ", ncol(seu_treatment), "\n",
      "\n--- Output Files ---\n",
      "CellChat list: ", cellchat_list_file, "\n"
    )
    
    summary_file <- file.path(SAVE_PATH, "analysis_summary.txt")
    writeLines(summary_text, summary_file)
    
  } else {
    
    cat("[CellChat Runner] Running CellChat on all samples together\n")
    
    cellChat <- run_cellchat_core(
      seu = seurat_obj,
      sample_name = "All Together",
      cluster_col = CLUSTER_COLUMN,
      species = SPECIES,
      pathway_types = PATHWAY_TYPES,
      prob_method = prob_method,
      prob_trim = prob_trim,
      min_cells = MIN_CELLS,
      seurat_full = seurat_obj
    )
    
    cat("[CellChat Runner] Saving CellChat object\n")
    cellchat_file <- file.path(SAVE_PATH, "cellchat_object.rds")
    saveRDS(cellChat, cellchat_file)
    
    cat("[CellChat Runner] Generating Excel file\n")
    create_cellchat_excel(cellChat, SAVE_PATH, single_analysis = TRUE)
    
    summary_text <- paste0(
      "=== CellChat Analysis Summary (All Together) ===\n",
      "Analysis Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
      "Sample Column: All Together\n",
      "Cluster Column: ", CLUSTER_COLUMN, "\n",
      "Species: ", SPECIES, "\n",
      "Pathway Types: ", paste(PATHWAY_TYPES, collapse = ", "), "\n",
      "Probability Method: ", PROB_TYPE, "\n",
      "Min Cells: ", MIN_CELLS, "\n",
      "Cell Types: ", length(unique(cellChat@idents)), "\n",
      "Significant Pathways: ", length(cellChat@netP$pathways), "\n",
      "\n--- Output Files ---\n",
      "CellChat object: ", cellchat_file, "\n"
    )
    
    summary_file <- file.path(SAVE_PATH, "analysis_summary.txt")
    writeLines(summary_text, summary_file)
  }
  
  future::plan("sequential")
  
  cat("[CellChat Runner] Analysis completed successfully!\n")
  
  if (!is.na(EMAIL) && nzchar(EMAIL)) {
    body <- c(
      "Hi,",
      "",
      "The CellChat analysis has been completed successfully.",
      "",
      paste0("Results are saved in: ", SAVE_PATH),
      "",
      "Best,",
      "Bingtian"
    )
    send_email(to = EMAIL, subject = "CellChat Analysis - COMPLETED", body = body)
    cat("[CellChat Runner] Success email sent to:", EMAIL, "\n")
  }
  
}, error = function(e) {
  
  err_msg <- conditionMessage(e)
  cat("[CellChat Runner] ERROR:", err_msg, "\n")
  cat("[CellChat Runner] Stack trace:\n")
  print(traceback())
  
  tryCatch(future::plan("sequential"), error = function(e) {})
  
  if (!is.na(EMAIL) && nzchar(EMAIL)) {
    body <- c(
      "Hi,",
      "",
      "There was an error in the CellChat analysis.",
      "",
      paste0("Error message: ", err_msg),
      "",
      paste0("Output folder: ", SAVE_PATH),
      "",
      "Please check the log file for more details.",
      "",
      "Best,",
      "Bingtian"
    )
    send_email(to = EMAIL, subject = "CellChat Analysis - ERROR", body = body)
    cat("[CellChat Runner] Error email sent to:", EMAIL, "\n")
  }
  
  quit(status = 1)
})