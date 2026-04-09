# Since the de_analysis uses a pseudo-bulk approach, low expression genes in 
# scRNA-seq data might not be detected due to the inherently low signal. However, 
# to maintain consistency with the SEGEX dataset, DW would likely not allow modifying
# this method. Therefore, if someone intends to perform differential expression
# analysis on scRNA-seq or snRNA-seq data, I would recommend using Seurat’s FindMarkers function instead.
# --Bingtian write on Apr. 6th, 2026
library(stringr)
library(tidyverse)
library(NotationConverter)
library(FindMarkersLoupe)
library(patchwork)
library(gridExtra)
library(cowplot)
# devtools::install_local("/projectnb/wax-es/00_shinyapp/DEG/SOFT/FindMarkersLoupe", force =TRUE)
## TMP: original location FindMarkersLoupe, mm10 converter and Utils.R: /projectnb/wax-dk/max/RSRC/G190

DEBUG <- FALSE

if(DEBUG){
  seurat_obj <- readRDS("/projectnb/wax-es/00_shinyapp/Clustering/file/G193_Male.rds")
  de_config_validated <- read.csv("test/DEG_config.csv", sep = "\t")
  sample_column = "sample_id"
  cluster_column = "seurat_clusters"
  # s1 = c("G193M1", "G193M2")
  # s2 = c("G193M1", "G193M2")
  # c1 = 2
  # c2 = 3
  # de_config_validated <- tibble(
  #   SAMPLE_ID_1 = list(c("G193M1", "G193M2")),
  #   SAMPLE_ID_2 = list(c("G193M1", "G193M2")),
  #   CLUSTER_ID_1 = 2,
  #   CLUSTER_ID_2 = 3,
  #   ix = 1
  # )
}

# de_config <- read_csv(argv$de_config, col_names = T, show_col_types = FALSE) %>% 
#   distinct()

get_condition_by_sample_id <- function(sample_column, sid, sobj) {
  sobj@meta.data %>% 
    dplyr::select(all_of(sample_column), condition) %>% 
    dplyr::filter(.data[[sample_column]] %in% sid) %>% 
    dplyr::distinct(condition) %>% 
    dplyr::pull(condition) %>% 
    paste(collapse = "-")
}

mp_extract_n_top_bottom <- function(t,n){
  rbind(slice_head(t,n = n), slice_tail(t, n = n)) %>% 
    distinct()
}

compute_de <- function(seurat_obj, s1,c1,s2,c2, ix, sample_column, cluster_column) {
  #' Calculate pairwise clusters comparisons
  #'
  #' @param seurat_obj 
  #' @param s1 sample id 1 
  #' @param s2 sample id 2
  #' @param c1 cluster id 1
  #' @param c2 cluster id 2
  #'
  #' @return list of the following objects:
  #' - segex_output table with differential expression of clusters
  #' - segex_filename output file for segex filename
  #' - pdf with dotplot
  #' 
  #' @export
  #'
  #' @examples
  ## create new Idents and use them 
    ## processing of sample/sample
    ## using "sampleid_clusterid"
  split_vec <- function(x) {
    strsplit(as.character(x), ",")[[1]]
  }
  
  s1 <- split_vec(s1)
  s2 <- split_vec(s2)
  c1 <- split_vec(c1)
  c2 <- split_vec(c2)
  meta <- seurat_obj@meta.data
  new_idents_df <- meta %>%
    mutate(new_idents = case_when(
      .data[[sample_column]] %in% s1 & .data[[cluster_column]] %in% c1 ~ "GROUP1",
      .data[[sample_column]] %in% s2 & .data[[cluster_column]] %in% c2 ~ "GROUP2",
      TRUE ~ "Other"
    )) %>%
    select(CB, new_idents)
  
    id.1 <- "GROUP1"
    id.2 <- "GROUP2"
  
  print("Number of cellls in each group")
  print(table(new_idents_df$new_idents))
  
  ## add new Idents to meta.data
  seurat_obj <- AddMetaData(seurat_obj, new_idents_df)
  
  ## activate new Idents
  Idents(seurat_obj) <- "new_idents"
  
  ## PROCESSING
  
  ## findMarkersLoupe (we cannot use the "short" Seurat object, because the entire 
  ## matrix must be used to calculate the intensities)
  markers_short <- FindMarkersLoupe(seurat_obj, id.1 = id.1, id.2 = id.2, formatted = "short")

  ## prepare Segex output data.frame
  segex_output <- exportToSegex(input_df = markers_short, from = "mm10")
  
  s1_label <- get_condition_by_sample_id(sample_column, s1, seurat_obj)
  s2_label <- get_condition_by_sample_id(sample_column, s2, seurat_obj)
  s1_str <- paste(s1, collapse = "-")
  s2_str <- paste(s2, collapse = "-")
  
  c1_str <- paste(c1, collapse = "-")
  c2_str <- paste(c2, collapse = "-")
  segex_fn <- stringr::str_glue(
    "{ix}_scLoupe_{s1_str}_{s1_label}_{c1_str}_vs_{s2_str}_{s2_label}_{c2_str}_DiffExp_IntronicMonoExonic.tsv"
  )
  

  
  # # DOTPLOT
  # # subset of only required clusters (need to create 2 lines DotPlot)
  # only_clusters_CB <- WhichCells(seurat_obj, idents = c(id.1, id.2))
  # # short 'Seurat' object only CB related to clusters
  # only_clusters_seurat <- subset(seurat_obj, cells = only_clusters_CB)
  # 
  # ## get top 30 genes from both sides
  # ## TODO: double plot top +/-30 genes and +/-30 lncRNA
  # ## TODO: need to create histograms which show shift of pvalue if we compare two
  # ## different is size clusters (small clusters will not have any significant 
  # ## genes by pvalue, but some genes will have good Log2FC)
  # 
  # top60genes <- markers_short %>% 
  #   #filter( 0.5*((id.2.intensity+1)+ (id.2.intensity+1)) > 1, !grepl("lnc", gname))
  #   filter(!grepl("lnc", gname)) %>% 
  #   arrange(desc(log2_fold_change)) %>% 
  #   mp_extract_n_top_bottom(., n = 30) %>% 
  #   pull(gname)
  # 
  # top60lncrna <- markers_short %>% 
  #   filter(grepl("lnc",gname)) %>% 
  #   arrange(desc(log2_fold_change)) %>% 
  #   mp_extract_n_top_bottom(.,  n = 30) %>% 
  #   pull(gname)
  # 
  # ncells <- table(new_idents_df$new_idents)
  # 
  # id.1.ncells <- ncells[["GROUP1"]]
  # id.2.ncells <- ncells[["GROUP2"]]
  # id.1p <- str_glue("{s1_str}_{c1_str}")
  # id.2p <- str_glue("{s2_str}_{c2_str}")
  # 
  # cols <- rev(c("#225ea8","#6baed6","#eff3ff","#fbb4b9","#fbb4b9"))
  # suppressWarnings({
  # top60genes_dotplot <- wrap_elements(mp_dotplot(only_clusters_seurat, 
  #                                                top60genes, 
  #                                                title = str_glue("{id.1p}({id.1.ncells} cells) vs {id.2p}({id.2.ncells} cells). Top 30 genes with average intensity > 1")) + 
  #                                       scale_colour_gradientn(colors = cols))
  # 
  # top60lncrna_dotplot <- wrap_elements(mp_dotplot(only_clusters_seurat, 
  #                                                 top60lncrna, 
  #                                                 title = str_glue("{id.1p}({id.1.ncells} cells) vs {id.2p}({id.2.ncells} cells). Top 30 lncRNA")) + 
  #                                        scale_colour_gradientn(colors = cols))
  # })
  # top60_dotplot <- top60genes_dotplot/top60lncrna_dotplot+plot_layout(heights = c(5,4))
  # 
  ## back to usual Idents and clean meta.data
  Idents(object = seurat_obj) <- seurat_obj@meta.data$seurat_clusters
  seurat_obj@meta.data$new_idents <- NULL
  
  # if (DEBUG) {
  #   list(segex_output = segex_output$segex,
  #        lp = markers_short,
  #        sr = markers_short_seurat,
  #        segex_filename = segex_fn,
  #        pdf = top60_dotplot)
  # } else {
  list(segex_output = segex_output$segex,
       segex_filename = segex_fn)
    # list(segex_output = segex_output$segex,
    #      segex_filename = segex_fn,
    #      pdf = top60_dotplot)
    
  # }
}

# mp_dotplot <- function(sobj, genes, title = "", ...) {
#   suppressWarnings(
#   DotPlot(sobj, features = genes, assay = "RNA", ...))+
#     RotatedAxis()+ 
#     {if(title !="") ggtitle(title)} +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#           axis.title.x = element_blank(),
#           axis.title.y = element_blank(),
#           legend.direction = "horizontal", 
#           legend.position = "bottom",
#           legend.box = "horizontal", 
#           legend.justification = "center",
#           legend.text=element_text(size=8),
#           legend.title = element_text(size = 10))
# }

# TODO: average intensity does not work need replacing for something else
# 

# z.t3 <- compute_de(input_seurat, "G190M2", 1, "G190M2", 2, 1)
# res <- left_join(z.t3$lp,z.t3$sr)
# z.t3$pdf

# if (DEBUG){
#   z.t1 <- compute_de(input_seurat, "AGGR", 1, "AGGR", 2, 1)
#   z.t2 <- compute_de(input_seurat, "AGGR", 1, "G190M2", 0, 1)
#   z.t3 <- compute_de(input_seurat, "G190M2", 1, "G190M2", 2, 1)
#   z.t4 <- compute_de(input_seurat, "G190M2", 1, "G183M1", 1, 1)
# }
# 
# View(z.t1$tmp)
# summary(sd(z.t2$tmp$adjusted_p_value))
# #zscore <- (z.t1$tmp$adjusted_p_value-median(z.t1$tmp$adjusted_p_value))/sd(z.t1$tmp$adjusted_p_value)
# zscore <- z.t2$tmp$adjusted_p_value
# zscore <- zscore[zscore != 1.0]
# hist(zscore, breaks = 100)

# z.t2 <- compute_de(input_seurat, "AGGR", 1, "G190M2", 0, 1)
# z.t1$pdf
# z.tmp <- z.t2$tmp %>% 
#   filter(!grepl("lnc",gname) & adjusted_p_value < 0.05)
# View(z.tmp)  
#
#
# z.t4 <- compute_de(input_seurat, "G190M2", 1, "G183M1", 1, 1)
# z.tmp <- z.t4$tmp %>% 
#   filter(!grepl("lnc",gname) & adjusted_p_value < 0.05)
# 
# View(z.t4$tmp)
#
# 
# View(z.t2$segex_output)
# z.t1 <- compute_de(input_seurat, "AGGR", 0, "AGGR", 2, 1)
# z.tmp <- z.t1$tmp %>% 
#   filter(!grepl("lnc",gname) & adjusted_p_value < 0.05)
# View(z.t1$tmp)
# View(z.tmp)

# HERE

# print("Saving Segex TSVs and prepare list of pdfs to save")
# pdf_list <- pmap(de_config_validated, function(SAMPLE_ID_1, SAMPLE_ID_2, CLUSTER_ID_1, CLUSTER_ID_2, ix) {
# 
#   res_list <- compute_de(seurat_obj, SAMPLE_ID_1, CLUSTER_ID_1, SAMPLE_ID_2, CLUSTER_ID_2, ix, sample_column, cluster_column)
#   if (!DEBUG)  write_tsv(res_list$segex_output, res_list$segex_filename, col_names = T)
# 
#   res_list$pdf
# })
# start <- Sys.time()
# res_list_all <- pmap(de_config_validated, function(SAMPLE_ID_1, SAMPLE_ID_2, CLUSTER_ID_1, CLUSTER_ID_2, ix) {
#   res_list <- compute_de(seurat_obj, SAMPLE_ID_1, CLUSTER_ID_1, SAMPLE_ID_2, CLUSTER_ID_2, ix, sample_column, cluster_column)
#   write_tsv(res_list$segex_output, file.path(path, res_list$segex_filename), col_names = T)
# })
# end <- Sys.time()
# print(start-end)
# print("Saving PDFs")
# pdf_list %>%
#   map(function(p){plot_grid(wrap_plots(p))}) %>%
#   marrangeGrob(nrow = 3, ncol = 1) %>%
#   ggsave(filename = "dotplots_by_comparision.pdf", width = 15.50, height = 20)
