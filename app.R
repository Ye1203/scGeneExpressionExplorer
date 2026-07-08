library(shiny)
library(shinyjs)
library(Seurat)
library(openxlsx)
library(ggplot2)
library(sortable)
library(colourpicker)
library(cowplot)
library(dplyr)
library(stringr)
library(tidyr)
library(ggrepel)
library(sortable)
library(CellChat)
# devtools::install_github("mpyatkov/NotationConverter")
library(NotationConverter)
# source(file.path(getwd(), "DE_analysis.R"))
# ============================================================
# -------------------- Utility Functions ---------------------
# ============================================================

# Load differnt type of data (.rds, .tsv, .xlsx)
load_data_object <- function(path, type = c("seurat", "deg", "gsea", "cccd")) {
  type <- match.arg(type)
  validate(need(file.exists(path), paste(type, "file does not exist")))
  
  if (type == "seurat") {
    readRDS(path)
  } else if (type == "deg") {
    data <- read.delim(path, header = TRUE, stringsAsFactors = FALSE)
    data <- NotationConverter::notationConverter(data, from = "segex", to = "mm10", column = "segex", replace_column = TRUE)
    data$ratio <- log2(data$ratio)
    colnames(data)[colnames(data) == "segex"] <- "gene"
    colnames(data)[colnames(data) == "ratio"] <- "log2FC"
    colnames(data)[colnames(data) == "fc"] <- "FC (linear)"
    return(data)
  } else if (type == "gsea") {
    data <- openxlsx::read.xlsx(path, sheet = 1)
    return(data)
  } else if (type == "cccd"){
    data <- readRDS(path)
    if(is.list(data)){
      data_merge <- mergeCellChat(data, add.names = names(data))
    return(list(data = data,
                data_merge = data_merge
    ))}else{
             return(list(data = data,
                        data_merge = NULL))
           }
}}

extract_meta_columns <- function(data) {
  if (!"meta.data" %in% slotNames(data)) return(NULL)
  meta_cols <- colnames(data@meta.data)
  
  meta_df <- data@meta.data
  valid_cols <- meta_cols[
    sapply(meta_df[, meta_cols, drop = FALSE], function(x) {
      !is.numeric(x)
    })
  ]
  
  c("", setdiff(valid_cols, c("CB", "CB_original")))
}

extract_unique_values <- function(data, column_name) {
  
  if (is.null(data) || is.null(column_name)) {
    return(NULL)
  }
  
  x <- data@meta.data[[column_name]]
  
  if (is.factor(x)) {
    levels(x)
  } else {
    unique(as.character(x))
  }
}

generate_combinations <- function(values) {
  if (length(values) < 2) return(data.frame())
  comb <- expand.grid(numerator = values, denominator = values, stringsAsFactors = FALSE)
  comb <- comb[comb$numerator != comb$denominator, ]
  rownames(comb) <- NULL
  comb
}

# ============================================================
# -------------------- Comparison Module ---------------------
# ============================================================

comparisonTableUI <- function(id, title_prefix) {
  ns <- NS(id)
  uiOutput(ns("table_ui"))
}

# Step 1 UI
comparisonTableServer <- function(id, column_values, prefix) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    rv <- reactiveValues(rows = list())
    row_count <- reactiveVal(0)
    
    output$table_ui <- renderUI({
      req(column_values())
      tagList(
        fluidRow(
          div(
            style = "display: inline-block; vertical-align: middle; margin-left:20px;",
            actionButton(ns("add_all"), "Add All Comparison", class = "btn-primary", width = "100%")
          ),
          div(
            style = "display: inline-block; vertical-align: middle; margin-left:5px;",
            actionButton(ns("remove_all"), "Remove All", class = "btn-warning", width = "100%")
          ),
          div(
            style = "display: inline-block; vertical-align: middle; margin-left:5px;",
            actionButton(ns("add_row"), "Add One Row", class = "btn-success", width = "100%")
          )
        ),
        br(),
        fluidRow(
          column(3, strong("Treatment (numerator)"), style = "text-align:center;"),
          column(3, strong("Control (denominator)"), style = "text-align:center;"),
          column(2, strong("Swap N/E")),
          column(2, strong("Action"))
        ),
        br(),
        div(id = ns("rows_container"))
      )
    })
    
    add_row_ui <- function(num = NULL, den = NULL) {
      row_idx <- as.character(Sys.time())
      row_idx <- gsub("[^0-9]", "", row_idx)  
      row_id <- paste0("row_", row_idx)
      
      insertUI(
        selector = paste0("#", ns("rows_container")),
        where = "beforeEnd",
        ui = div(
          id = ns(row_id),
          fluidRow(
            style = "margin-top: 5px;",
            column(3, selectInput(ns(paste0("num_", row_id)), NULL,
                                  choices = column_values(), selected = num,
                                  width = "100%")),
            column(3, selectInput(ns(paste0("den_", row_id)), NULL,
                                  choices = column_values(), selected = den,
                                  width = "100%")),
            column(2, actionButton(ns(paste0("exchange_", row_id)), "⇄ Swap",
                                   style = "color:#fff;background-color:#337ab7;border-color:#2e6da4;",
                                   width = "100%")),
            column(2, actionButton(ns(paste0("delete_", row_id)), "✕ Delete",
                                   style = "color:#fff;background-color:#d9534f;border-color:#d43f3a;",
                                   width = "100%"))
          )
        )
      )
      
      rv$rows[[row_id]] <- list(numerator = num, denominator = den)
      
      # DELETE
      observeEvent(input[[paste0("delete_", row_id)]], {
        removeUI(selector = paste0("#", ns(row_id)))
        rv$rows[[row_id]] <- NULL
      }, ignoreInit = TRUE)
      
      # SWAP
      observeEvent(input[[paste0("exchange_", row_id)]], {
        row <- isolate(rv$rows[[row_id]])
        if (is.null(row)) return()
        
        updateSelectInput(session, paste0("num_", row_id), selected = row$denominator)
        updateSelectInput(session, paste0("den_", row_id), selected = row$numerator)
        
        rv$rows[[row_id]]$numerator <- row$denominator
        rv$rows[[row_id]]$denominator <- row$numerator
      }, ignoreInit = TRUE)
      
      observeEvent(input[[paste0("num_", row_id)]], {
        row <- isolate(rv$rows[[row_id]])
        if (is.null(row)) return()
        rv$rows[[row_id]]$numerator <- input[[paste0("num_", row_id)]]
      }, ignoreInit = TRUE)
      
      observeEvent(input[[paste0("den_", row_id)]], {
        row <- isolate(rv$rows[[row_id]])
        if (is.null(row)) return()
        rv$rows[[row_id]]$denominator <- input[[paste0("den_", row_id)]]
      }, ignoreInit = TRUE)
      
      row_count(length(rv$rows))
    }
    
    clear_all_rows <- function() {
      lapply(names(rv$rows), function(rid) {
        removeUI(selector = paste0("#", ns(rid)), immediate = TRUE)
      })
      rv$rows <- list()
      row_count(0)
    }
    
    observeEvent(input$add_row, { add_row_ui() })
    
    observeEvent(input$add_all, {
      vals <- column_values()
      if (length(vals) > 1) {
        clear_all_rows()
        comb <- t(combn(vals, 2))
        apply(comb, 1, function(row) add_row_ui(num = row[1], den = row[2]))
      }
    })
    
    observeEvent(input$remove_all, { clear_all_rows() })
    
    return(reactive(rv$rows))
  })
}

# Step 2 UI
step2UI <- function(id) {
  ns <- NS(id)
  
  tagList(
    h4("Sample Comparisons:"),
    uiOutput(ns("add_sample_comp_ui")),
    br(),
    div(id = ns("sample_cards_container"), class = "row"),
    
    fluidRow(column(12, hr())),
    
    h4("Cluster Comparisons:"),
    uiOutput(ns("add_cluster_comp_ui")),
    br(),
    div(id = ns("cluster_cards_container"), class = "row"),
    
    fluidRow(
      column(12,hr()),
      column(
        12,
        actionButton(
          "back_step",
          "Back to Step 1",
          class = "btn-warning",
          style = "width:200px;"
        ),
        actionButton(
          "analysis_summary",
          "Review and Analysis",
          class = "btn-success",
          stype = "width:200px;"
        )))
  )
}

color_selected <- function(color_length) {
  color_total <- c(
    "#e6194b","#ffe119","#46f0f0","#f58231","#bcf60c", 
    "#ff00ff","#9a6324","#fffac8","#e6beff","#00bfff", 
    "#ffd8b1","#00ff7f","#f5a9bc","#1e90ff","#ffa500",
    "#98fb98","#911eb4","#afeeee","#fa8072","#9acd32",
    "#3cb44b","#000075","#808000","#cd5c5c","#dda0dd",
    "#40e0d0","#ff69b4","#8a2be2","#c71585","#5f9ea0",
    "#dc143c","#87cefa","#ff6347","#9932cc","#00ced1",
    "#ff4500","#6a5acd","#b0e0e6","#d2691e","#a9a9f5",
    "#adff2f","#8b0000","#7fffd4","#00fa9a","#ba55d3",
    "#2e8b57","#ffdab9","#b22222","#ffe4e1","#7b68ee"
  )
  
  if (color_length <= length(color_total)) {
    return(color_total[1:color_length])
  }
  
  warning("groups is larger than 50, color will randomly select")
  palette_fn <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))
  extra_colors <- palette_fn(color_length - length(color_total))
  colors <- c(color_total, extra_colors)
  return(colors)
}

# ----------------- Step2 Server Module -----------------
step2Server <- function(id, sample_rows, cluster_rows, sample_meta, cluster_meta) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Store observers for cleanup
    card_observers <- reactiveValues(
      sample = list(),
      cluster = list()
    )
    
    output$add_sample_comp_ui <- renderUI({
      rows <- sample_rows()
      if (is.null(rows) || length(rows) == 0){
        div(
          style = "color:#555; font-style:italic;",
          "No sample comparisons found. You can add comparisons in Step 1 using the 'Back to Step 1' button."
        )
      }else{
        div(
          style="display:flex; gap:10px;",
          actionButton(
            ns("add_sample_comp"),
            "Add One Card",
            class = "btn-primary"
          ),
          actionButton(
            ns("add_all_sample_cards"),
            "Add All",
            class = "btn-success"
          ),
          actionButton(
            ns("delete_all_sample_cards"),
            "Delete All",
            class = "btn-danger"
          ),
          p("If Step 1 comparisons change, recreate the cards.",
            style = "color:red;")
        )
      }
    })
    
    output$add_cluster_comp_ui <- renderUI({
      rows <- cluster_rows()
      if (is.null(rows) || length(rows) == 0){
        div(
          style = "color:#555; font-style:italic;",
          "No celltype cluster comparisons found."
        )
      }else{
        div(
          style="display:flex; gap:10px;",
          actionButton(
            ns("add_cluster_comp"),
            "Add One Card",
            class = "btn-primary"
          ),
          actionButton(
            ns("add_all_cluster_cards"),
            "Add All",
            class = "btn-success"
          ),
          actionButton(
            ns("delete_all_cluster_cards"),
            "Delete All",
            class = "btn-danger"
          ),
          p("If Step 1 comparisons change, recreate the cards.",
            style = "color:red;")
        )
      }
    })
    
    # Helper: Create card UI
    createCard <- function(card_title, available_comps, meta_choices, card_id, is_sample = TRUE) {
      # Create checkbox choices from all available comparisons
      comp_choices <- list()
      for (comp in available_comps) {
        comp_choices[[comp$name]] <- comp$id
      }
      
      div(
        id = ns(card_id),  # Add ID to the main div for easy removal
        class = "col-sm-4",  
        style = "margin-bottom: 15px;",
        wellPanel(
          h5(card_title, style = "margin-top: 0; font-weight: bold;"),
          
          selectInput(
            ns(paste0(card_id, "_meta_selector")),
            label = ifelse(is_sample, "Cluster to use in comparison:", "Sample to use in comparison:"),
            choices = meta_choices,
            selected = NULL,
            multiple = TRUE,
            width = "100%"
          ),
          
          checkboxGroupInput(
            ns(paste0(card_id, "_comparison_checkboxes")),
            label = "Select comparisons to perform:",
            choices = comp_choices,
            selected = NULL,
            width = "100%"
          ),
          
          fluidRow(
            column(6, actionButton(ns(paste0(card_id, "_select_all")), "Select All", 
                                   class = "btn-sm btn-primary", width = "100%")),
            column(6, actionButton(ns(paste0(card_id, "_clear_all")), "Clear", 
                                   class = "btn-sm btn-warning", width = "100%"))
          ),
          br(),
          actionButton(ns(paste0(card_id, "_delete")), "Delete", 
                       class = "btn-danger", style = "width:100%; margin-top:5px;")
        )
      )
    }
    
    # Get all valid sample comparisons from Step1
    all_sample_comps <- reactive({
      comps <- sample_rows()
      
      if (is.null(comps) || length(comps) == 0) {
        return(list())
      }
      
      comps <- comps[unlist(lapply(comps, function(x) {
        is.list(x) && !is.null(x$numerator) && !is.null(x$denominator)
      }))]
      
      comp_list <- list()
      for (i in seq_along(comps)) {
        comp <- comps[[i]]
        comp_id <- paste0("S", i, "_", comp$numerator, "_vs_", comp$denominator)
        comp_name <- paste0("S", i, ": ", comp$numerator, " vs ", comp$denominator)
        comp_list[[i]] <- list(id = comp_id, name = comp_name, 
                               numerator = comp$numerator, denominator = comp$denominator)
      }
      comp_list
    })
    
    # Get all valid cluster comparisons from Step1
    all_cluster_comps <- reactive({
      comps <- cluster_rows()
      
      if (is.null(comps) || length(comps) == 0) {
        return(list())
      }
      
      comps <- comps[unlist(lapply(comps, function(x) {
        is.list(x) && !is.null(x$numerator) && !is.null(x$denominator)
      }))]
      
      comp_list <- list()
      for (i in seq_along(comps)) {
        comp <- comps[[i]]
        comp_id <- paste0("C", i, "_", comp$numerator, "_vs_", comp$denominator)
        comp_name <- paste0("C", i, ": ", comp$numerator, " vs ", comp$denominator)
        comp_list[[i]] <- list(id = comp_id, name = comp_name, 
                               numerator = comp$numerator, denominator = comp$denominator)
      }
      comp_list
    })
    
    # Function to destroy observers for a card
    destroyCardObservers <- function(card_id, card_type) {
      if (card_type == "sample") {
        if (!is.null(card_observers$sample[[card_id]])) {
          for (obs in card_observers$sample[[card_id]]) {
            obs$destroy()
          }
          card_observers$sample[[card_id]] <- NULL
        }
      } else {
        if (!is.null(card_observers$cluster[[card_id]])) {
          for (obs in card_observers$cluster[[card_id]]) {
            obs$destroy()
          }
          card_observers$cluster[[card_id]] <- NULL
        }
      }
    }
    
    # ----------------- Sample Cards Management -----------------
    sample_cards <- reactiveValues(ids = character(0), data = list())
    
    # Add new sample comparison card
    observeEvent(input$add_sample_comp, {
      sample_comps <- all_sample_comps()
      if (length(sample_comps) == 0) {
        showModal(modalDialog(
          title = "No Sample comparisons", 
          "Please add at least one sample comparison in Step1."
        ))
        return()
      }
      
      # Generate a truly unique card ID using timestamp or random number
      card_id <- paste0("sample_card_", format(Sys.time(), "%Y%m%d%H%M%S"), "_", 
                        paste(sample(letters, 4, replace = TRUE), collapse = ""))
      
      sample_cards$ids <- c(sample_cards$ids, card_id)
      
      sample_cards$data[[card_id]] <- list(
        type = "sample",
        selected_comps = character(0),
        selected_meta = character(0)
      )
      
      # Insert the card UI
      insertUI(
        selector = paste0("#", ns("sample_cards_container")),
        ui = createCard(
          card_title = "Sample Comparisons",
          available_comps = sample_comps,
          meta_choices = sample_meta(),
          card_id = card_id,
          is_sample = TRUE
        ),
        immediate = TRUE
      )
      
      # Set up observers for this card and store them
      card_observers$sample[[card_id]] <- setupCardObservers(card_id, "sample", sample_cards, all_sample_comps)
    })
    
    # ----------------- Cluster Cards Management -----------------
    cluster_cards <- reactiveValues(ids = character(0), data = list())
    
    # Add new cluster comparison card
    observeEvent(input$add_cluster_comp, {
      cluster_comps <- all_cluster_comps()
      if (length(cluster_comps) == 0) {
        showModal(modalDialog(
          title = "No Cluster comparisons", 
          "Please add at least one cluster comparison in Step1."
        ))
        return()
      }
      
      # Generate a truly unique card ID using timestamp or random number
      card_id <- paste0("cluster_card_", format(Sys.time(), "%Y%m%d%H%M%S"), "_", 
                        paste(sample(letters, 4, replace = TRUE), collapse = ""))
      
      cluster_cards$ids <- c(cluster_cards$ids, card_id)
      
      cluster_cards$data[[card_id]] <- list(
        type = "cluster",
        selected_comps = character(0),
        selected_meta = character(0)
      )
      
      # Insert the card UI
      insertUI(
        selector = paste0("#", ns("cluster_cards_container")),
        ui = createCard(
          card_title = "Cluster Comparisons",
          available_comps = cluster_comps,
          meta_choices = cluster_meta(),
          card_id = card_id,
          is_sample = FALSE
        ),
        immediate = TRUE
      )
      
      # Set up observers for this card and store them
      card_observers$cluster[[card_id]] <- setupCardObservers(card_id, "cluster", cluster_cards, all_cluster_comps)
    })
    
    # Helper function to set up card observers
    setupCardObservers <- function(card_id, card_type, cards_reactive, all_comps_func) {
      
      # Create input ID patterns
      delete_id <- paste0(card_id, "_delete")
      select_id <- paste0(card_id, "_select_all")
      clear_id <- paste0(card_id, "_clear_all")
      checkbox_id <- paste0(card_id, "_comparison_checkboxes")
      meta_id <- paste0(card_id, "_meta_selector")
      
      # Store observers in a list to return
      observers <- list()
      
      # Delete observer
      observers$delete <- observeEvent(input[[delete_id]], {
        # First destroy all observers for this card
        destroyCardObservers(card_id, card_type)
        
        # Then remove the UI
        removeUI(selector = paste0("#", ns(card_id)))
        
        # Update reactive values
        cards_reactive$ids <- setdiff(cards_reactive$ids, card_id)
        cards_reactive$data[[card_id]] <- NULL
      }, ignoreInit = TRUE, autoDestroy = TRUE)
      
      # Select All observer
      observers$select <- observeEvent(input[[select_id]], {
        all_comp_choices <- sapply(all_comps_func(), function(x) x$id)
        
        updateCheckboxGroupInput(
          session,
          checkbox_id,
          selected = all_comp_choices
        )
        
        # Only update if this card still exists
        if (!is.null(cards_reactive$data[[card_id]])) {
          cards_reactive$data[[card_id]]$selected_comps <- all_comp_choices
        }
      }, ignoreInit = TRUE, autoDestroy = TRUE)
      
      # Clear All observer
      observers$clear <- observeEvent(input[[clear_id]], {
        updateCheckboxGroupInput(
          session,
          checkbox_id,
          selected = character(0)
        )
        
        # Only update if this card still exists
        if (!is.null(cards_reactive$data[[card_id]])) {
          cards_reactive$data[[card_id]]$selected_comps <- character(0)
        }
      }, ignoreInit = TRUE, autoDestroy = TRUE)
      
      # Monitor checkbox changes
      observers$checkbox <- observeEvent(input[[checkbox_id]], {
        selected <- input[[checkbox_id]]
        # Only update if this card still exists
        if (!is.null(cards_reactive$data[[card_id]])) {
          cards_reactive$data[[card_id]]$selected_comps <- selected
        }
      }, ignoreNULL = FALSE, autoDestroy = TRUE)
      
      # Monitor meta selector changes
      observers$meta <- observeEvent(input[[meta_id]], {
        selected <- input[[meta_id]]
        # Only update if this card still exists
        if (!is.null(cards_reactive$data[[card_id]])) {
          cards_reactive$data[[card_id]]$selected_meta <- selected
        }
      }, ignoreNULL = FALSE, autoDestroy = TRUE)
      
      return(observers)
    }
    
    observeEvent(input$add_all_sample_cards, {
      
      metas <- sample_meta()
      sample_comps <- all_sample_comps()
      
      if (length(metas) == 0 || length(sample_comps) == 0) return()
      
      all_comp_ids <- sapply(sample_comps, function(x) x$id)
      
      for (meta in metas){
        
        card_id <- paste0(
          "sample_card_",
          format(Sys.time(), "%Y%m%d%H%M%S"),
          "_",
          paste(sample(letters, 4, replace = TRUE), collapse = "")
        )
        
        sample_cards$ids <- c(sample_cards$ids, card_id)
        
        sample_cards$data[[card_id]] <- list(
          type = "sample",
          selected_comps = all_comp_ids,
          selected_meta = meta
        )
        
        insertUI(
          selector = paste0("#", ns("sample_cards_container")),
          ui = createCard(
            card_title = "Sample Comparisons",
            available_comps = sample_comps,
            meta_choices = metas,
            card_id = card_id,
            is_sample = TRUE
          ),
          immediate = TRUE
        )
        
        card_observers$sample[[card_id]] <-
          setupCardObservers(card_id, "sample", sample_cards, all_sample_comps)
        
        updateSelectInput(
          session,
          paste0(card_id, "_meta_selector"),
          selected = meta
        )
        
        updateCheckboxGroupInput(
          session,
          paste0(card_id, "_comparison_checkboxes"),
          selected = all_comp_ids
        )
      }
      
    })
    
    observeEvent(input$add_all_cluster_cards, {
      
      metas <- cluster_meta()
      cluster_comps <- all_cluster_comps()
      
      if (length(metas) == 0 || length(cluster_comps) == 0) return()
      
      all_comp_ids <- sapply(cluster_comps, function(x) x$id)
      
      for (meta in metas){
        
        card_id <- paste0(
          "cluster_card_",
          format(Sys.time(), "%Y%m%d%H%M%S"),
          "_",
          paste(sample(letters, 4, replace = TRUE), collapse = "")
        )
        
        cluster_cards$ids <- c(cluster_cards$ids, card_id)
        
        cluster_cards$data[[card_id]] <- list(
          type = "cluster",
          selected_comps = all_comp_ids,
          selected_meta = meta
        )
        
        insertUI(
          selector = paste0("#", ns("cluster_cards_container")),
          ui = createCard(
            card_title = "Cluster Comparisons",
            available_comps = cluster_comps,
            meta_choices = metas,
            card_id = card_id,
            is_sample = FALSE
          ),
          immediate = TRUE
        )
        
        card_observers$cluster[[card_id]] <-
          setupCardObservers(card_id, "cluster", cluster_cards, all_cluster_comps)
        
        updateSelectInput(
          session,
          paste0(card_id, "_meta_selector"),
          selected = meta
        )
        
        updateCheckboxGroupInput(
          session,
          paste0(card_id, "_comparison_checkboxes"),
          selected = all_comp_ids
        )
      }
      
    })
    
    observeEvent(input$delete_all_sample_cards, {
      
      ids <- sample_cards$ids
      
      for(id in ids){
        
        destroyCardObservers(id, "sample")
        
        removeUI(
          selector = paste0("#", ns(id))
        )
        
      }
      
      sample_cards$ids <- character(0)
      sample_cards$data <- list()
      
    })
    
    observeEvent(input$delete_all_cluster_cards, {
      
      ids <- cluster_cards$ids
      
      for(id in ids){
        
        destroyCardObservers(id, "cluster")
        
        removeUI(
          selector = paste0("#", ns(id))
        )
        
      }
      
      cluster_cards$ids <- character(0)
      cluster_cards$data <- list()
      
    })
    # Return all selected data
    return(
      list(
        sample_cards_data = reactive(sample_cards$data),
        cluster_cards_data = reactive(cluster_cards$data),
        all_sample_comps = all_sample_comps,
        all_cluster_comps = all_cluster_comps
      )
    )
  })
}
# ============================================================
# ------------------------- UI -------------------------------
# ============================================================

ui <- fluidPage(
  htmltools::findDependencies(selectizeInput("dummy", NULL, choices = NULL)),
  useShinyjs(),
  titlePanel("scGene Expression Explorer"),
  
  tabsetPanel(
    id = "main_tabs",
    
    tabPanel("DEG", uiOutput("deg_ui")),
    tabPanel("GSEA", uiOutput("gsea_ui")),
    tabPanel("Cell-Cell Communication Analysis", uiOutput("cellchat_ui")),
    tabPanel(
      "Visualization",
      fluidRow(
        
        column(
          width = 2,
          tags$div(
            style = "
              background-color: #f8f9fa;
              padding: 10px;
              border-radius: 10px;
              height: 100%;
            ",
            
            tabsetPanel(
              id = "viz_tabs",
              type = "pills",
              tabPanel("Violin Plot", value = "violin"),
              tabPanel("Dot Plot", value = "dotplot"),
              tabPanel("Feature Plot", value = "feature"),
              tabPanel("Volcano Plot (Genes)", value = "volcano"),
              tabPanel("Volcano Plot (Genes in Pathways)", value = "pathway_volcano"),
              tabPanel("GSEA NES", value = "gsea_nes"),
              tabPanel("Cell-Cell-Communication", value = "cccv"),
              tabPanel("Cell-Cell-Communication (condition comparison)", value = "cccv_condition")
            )
          )
        ),
        column(
          width = 10,
          uiOutput("viz_panel")
        )
      )
    ),
    tabPanel(
      "Data Output", 
      fluidRow(
        
        column(
          width = 2,
          tags$div(
            style = "
              background-color: #f8f9fa;
              padding: 10px;
              border-radius: 10px;
              height: 100%;
            ",
            
            tabsetPanel(
              id = "do_tabs",
              type = "pills",
              tabPanel("Gene Expression", value = "gene_expression")
            )
          )
        ),
        column(
          width = 10,
          uiOutput("data_output_ui")
        )
      )
    )
  ),
  
  tags$style(HTML("
    .nav-pills {
      width: 100%;
    }
    .nav-pills > li {
      float: none;
      width: 100%;
      margin-bottom: 6px;
    }
    .nav-pills > li > a {
      border-radius: 8px;
      width: 100%;
      text-align: left;
      padding: 10px 12px;
    }
    
    .token-pool {
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
      padding: 10px;
      border: 1px solid #ddd;
      border-radius: 10px;
      background: #fafafa;
      margin-bottom: 15px;
    }
    
    .token {
      padding: 6px 12px;
      border-radius: 999px;
      background: #2c7be5;
      color: white;
      font-size: 13px;
      cursor: grab;
      user-select: none;
      box-shadow: 0 1px 3px rgba(0,0,0,0.15);
      white-space: nowrap;
    }
    
    .token:active {
      cursor: grabbing;
    }
    
    .sheet {
      height: 150px;
      border: 1px solid #ddd;
      border-radius: 10px;
      background: white;
      padding: 10px;
      overflow-y: auto;
      display: flex;
      flex-wrap: wrap;
      align-content: flex-start;
      gap: 8px;
    }
    
    .sheet-title {
      font-weight: bold;
      margin-bottom: 8px;
    }
    
    .sortable-ghost {
      opacity: 0.4;
    }
  "))
)

# ============================================================
# ------------------------ SERVER ----------------------------
# ============================================================

server <- function(input, output, session) {
  sortable::enable_modules()
  # Reactive value
  rv <- reactiveValues(
    data_obj = NULL,
    data_obj_path = NULL,
    
    deg_obj = NULL,
    deg_obj_path = NULL,
    
    gsea_obj = NULL,
    gsea_obj_path = NULL,
    
    cccd_obj = NULL,
    cccd_obj_path = NULL,
    
    visualization_marker = "",
    visualization_pathway = "",
    pathway_rename = NULL,
    
    visualization_sample_column = "",
    
    save_path = "",
    ui_freeze = FALSE
    
  )
  
  # ----------------- Load Data -----------------
  
  # Enhanced path_input_ui with automatic dependency tracking
  path_input_ui <- function(id, label, placeholder, value, load_btn_class = "btn-success") {
    
    input_path_id <- paste0(id, "_path")
    
    if (!is.null(value)) {
      fluidRow(
        column(
          10,
          tags$div(
            style = "padding: 8px; background-color: #e8f5e8;
               border-radius: 4px; border: 1px solid #c8e6c9;",
            "📁 ",
            tags$strong("Loaded: "),
            tags$code(value)
          )
        ),
        
        column(
          2,
          actionButton(
            inputId = paste0("reset_", id, "_obj_btn"),
            label = "Reset",
            class = "btn-warning",
            style = "width: 100%; height: 35px;"
          )
        )
      )%>% 
        tagAppendAttributes(
          style = "margin-bottom: 15px;"
        )
    } else {
      fluidRow(
        column(
          10,
          textInput(
            inputId = input_path_id,
            value = "",
            placeholder = placeholder,
            label = NULL,
            width = "100%"
          )
        ),
        column(
          2,
          actionButton(
            inputId = paste0("load_", id),
            label = "Load",
            class = load_btn_class,
            style = "width: 100%; height: 35px;"
          )
        )
      )
    }
  }
  
  # Seurat data loading
  observeEvent(input$load_data, {
    showModal(modalDialog(
      title = "Loading Data",
      "Please wait...",
      footer = NULL,
      easyClose = FALSE
    ))
    
    tryCatch({
      obj <- load_data_object(input$data_path, type = "seurat")
      rv$visualization_sample_column <- ""
      rv$data_obj <- obj
      rv$data_obj_path <- input$data_path
      rv$save_path <- paste0(gsub("\\.rds$", "", rv$data_obj_path), "_visualization")
      removeModal()
      
      showNotification("Seurat data loaded successfully!", type = "message")
      
    }, error = function(e) {
      removeModal()
      
      showModal(modalDialog(
        title = "Error",
        paste("Error:", e$message),
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
    })
  })
  
  # DEG data loading
  observeEvent(input$load_deg, {
    showModal(modalDialog(title = "Loading DEG Data", "Please wait...", footer = NULL, easyClose = FALSE))
    tryCatch({
      obj <- load_data_object(input$deg_path, type = "deg")
      rv$deg_obj <- obj
      rv$deg_obj_path <- input$deg_path
      if(is.null(rv$data_obj)){rv$save_path <- paste0(gsub("\\.tsv$", "", rv$deg_obj_path), "_visualization")}
      removeModal()
      showNotification("DEG data loaded successfully!", type = "message")
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # GSEA data loading
  observeEvent(input$load_gsea, {
    showModal(modalDialog(title = "Loading GSEA Data", "Please wait...", footer = NULL, easyClose = FALSE))
    tryCatch({
      obj <- load_data_object(input$gsea_path, type = "gsea")
      rv$gsea_obj <- obj
      rv$gsea_obj_path <- input$gsea_path
      if(is.null(rv$data_obj) & is.null(rv$deg_obj)){rv$save_path <- paste0(gsub("\\.xlsx$", "", rv$gsea_obj_path), "_visualization")}
      removeModal()
      showNotification("GSEA data loaded successfully!", type = "message")
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # GSEA data loading
  observeEvent(input$load_cccd, {
    showModal(modalDialog(title = "Loading cell cell communication data.", "Please wait...", footer = NULL, easyClose = FALSE))
    tryCatch({
      obj <- load_data_object(input$cccd_path, type = "cccd")
      rv$cccd_obj <- obj
      rv$cccd_obj_path <- input$cccd_path
      if(is.null(rv$data_obj)){rv$save_path <- paste0(gsub("\\.rds$", "", rv$cccd_obj_path), "_visualization")}
      removeModal()
      showNotification("Cell cell communication data loaded successfully!", type = "message")
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # ----------------- Reactive unique values sample value, cluster value -----------------
  sample_values <- reactive({
    req(rv$data_obj, input$sample_column)
    extract_unique_values(rv$data_obj, input$sample_column)
  })
  
  cluster_values <- reactive({
    req(rv$data_obj, input$cluster_column)
    extract_unique_values(rv$data_obj, input$cluster_column)
  })
  
  # ----------------- Modules -----------------
  sample_rows <- comparisonTableServer("sample_table", sample_values, prefix = "S")
  cluster_rows <- comparisonTableServer("cluster_table", cluster_values, prefix = "C")
  
  # ----------------- Step navigation -----------------
  observeEvent(input$next_step, {
    sample_comparisons <- sample_rows()
    cluster_comparisons <- cluster_rows()
    
    # Check if at least one comparison is selected (not NULL)
    has_sample <- any(sapply(sample_comparisons, function(x) !is.null(x$numerator) && !is.null(x$denominator)))
    has_cluster <- any(sapply(cluster_comparisons, function(x) !is.null(x$numerator) && !is.null(x$denominator)))
    if (!has_sample && !has_cluster) {
      showModal(modalDialog(
        title = "No Comparisons Selected",
        "Please add at least one comparison for Sample or Cluster analysis before proceeding.",
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
    }else if(input$sample_column == "" | input$cluster_column == ""){
      showModal(modalDialog(
        title = "Sample meta.data or Cluster meta.data was selected.",
        "Please make sure you have select meta.data for both sample and cluster.",
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
    }
    else {
      shinyjs::hide("step1_ui")
      shinyjs::show("step2_ui")
    }
  })
  
  observeEvent(input$back_step, {
    shinyjs::show("step1_ui")
    shinyjs::hide("step2_ui")
  })
  
  # ----------------- DEG UI -----------------
  output$deg_ui <- renderUI({
    
    if (is.null(rv$data_obj)) {
      tagList(
        h4("Input Data Path"),
        path_input_ui(id = "data", label = "Seurat File Path", 
                      placeholder = "/path/to/seurat_file.rds", 
                      value = rv$data_obj_path),
        h4("Please load data first.")
      )
    } else {
      tagList(
        # Step 1 (original content wrapped in div)
        br(),
        tags$div(
          style = "padding: 8px; background-color: #e8f5e8; border-radius: 4px; border: 1px solid #c8e6c9; display: flex; align-items: center; justify-content: space-between;",
          tags$p(
            style = "margin: 0;",
            "📁 ",
            tags$strong("Loaded: "),
            tags$code(input$data_path)
          ),
          actionButton(
            "reset_data_obj_btn", 
            "Reset (change RDS file)", 
            class = "btn-warning", 
            style = "min-width: 180px; height: 35px; padding: 5px 10px;"
          )
        ),
        br(),
        div(id = "step1_ui",
            fluidRow(
              column(6, selectInput("sample_column", "meta.data sample identification", 
                                    choices = extract_meta_columns(rv$data_obj), width = "100%")),
              column(6, selectInput("cluster_column", "meta.data celltype cluster identification", 
                                    choices = extract_meta_columns(rv$data_obj), width = "100%"))
            ),
            hr(),
            h3("Sample Analysis"),
            comparisonTableUI("sample_table", "S"),
            hr(),
            h3("Cluster Analysis"),
            comparisonTableUI("cluster_table", "C"),
            hr(),
            actionButton("next_step", "NEXT STEP", class = "btn-primary", style = "width:200px;")
        ),
        
        # Step 2 (hidden initially)
        div(id = "step2_ui", style="display:none;", step2UI("step2_module"))
      )
    }
  })
  
  step2_data <- step2Server(
    "step2_module",
    sample_rows = sample_rows,
    cluster_rows = cluster_rows,
    sample_meta = cluster_values, 
    cluster_meta = sample_values  
  )
  
  # ----------------- GSEA UI -----------------
  output$gsea_ui <- renderUI({
    fluidPage(
      useShinyjs(),  
      tags$style(HTML("
      .wrap-selected .selectize-control.multi .selectize-input > div {
        white-space: normal !important;
        display: block !important;
        width: 95%;
        margin: 4px 2px;
        padding: 4px 8px;
        background-color: #e6f3ff;
        border: 1px solid #b8d4f0;
        border-radius: 3px;
      }
      
      .wrap-selected .selectize-control.multi .selectize-input {
        min-height: 100px;
        max-height: 300px;
        overflow-y: auto;
        padding: 8px;
      }
    ")),
      h3("GSEA Analysis"),
      
      # GSEA file - full row
      fluidRow(
        column(12,
               textInput("gsea_file", "Input file path (.tsv or folder):", "", width = "100%")
        )
      ),
      hr(),
      # Method + Species - same row, 1/4 each
      fluidRow(
        column(3,
               selectInput(
                 "gsea_method",
                 "Ranking method:",
                 choices = c("fc-pvalue","log2fc-pvalue","log2fc","fc"),
                 selected = "fc",
                 width = "100%"
               )
        ),
        column(6,
               radioButtons(
                 "species",
                 "Species (select one):",
                 choices = c("Mus musculus","Homo sapiens"),
                 inline = TRUE
               )
        )
      ),
      hr(),
      # DB selection - full row
      div(
        id = "db_container",
        class = "wrap-selected",  
        fluidRow(
          column(12,
                 selectizeInput(
                   "gsea_db",
                   "Select gene set databases:",
                   choices = NULL,  
                   multiple = TRUE,
                   width = "100%"
                 )
          )
        )
      ),
      
      br(),
      actionButton("preview_gsea", "Preview GSEA", class = "btn-success"),
      br(),br(),br(),br(),br()
    )
  })
  
  observeEvent(input$species, {
    req(input$species)
    species <- input$species
    if(species == "Homo sapiens"){
      db_choices <- c(
        "H - Hallmark gene sets — Coherent biological state or process signatures." = "H(Hallmark)",  
        "C1 - Positional gene sets — Gene sets grouped by chromosomal location." = "C1(Positional)", 
        "C2 (CGP) - Chemical and Genetic Perturbation gene sets." = "C2(CGP)", 
        "C2 (CP) - Canonical pathway gene sets." = "C2(CP)", 
        "C2 (CP:BIOCARTA) - BioCarta curated signaling pathways." = "C2(CP_BioCarta)",  
        "C2 (CP:KEGG_LEGACY) - KEGG legacy pathway gene sets." = "C2(CP_KEGG_LEGACY)",  
        "C2 (CP:KEGG_MEDICUS) - KEGG Medicus pathway gene sets." = "C2(CP_KEGG_MEDICUS)",  
        "C2 (CP:PID) - Pathway Interaction Database gene sets." = "C2(CP_PID)",  
        "C2 (CP:REACTOME) - Reactome curated pathway gene sets." = "C2(CP_Reactome)", 
        "C2 (CP:WIKIPATHWAYS) - WikiPathways curated pathway gene sets." = "C2(CP_WikiPathways)", 
        "C3 (MIR) - microRNA target gene sets." = "C3(MIR)",  
        "C3 (MIR:MIRDB) - miRDB microRNA target predictions." = "C3(MIR_MIRDB)", 
        "C3 (MIR:MIR_LEGACY) - Legacy microRNA target gene sets." = "C3(MIR_MIR_LEGACY)", 
        "C3 (TFT) - Transcription factor target gene sets." = "C3(TFT)",  
        "C3 (TFT:GTRD) - GTRD transcription factor target gene sets." = "C3(TFT_GTRD)", 
        "C3 (TFT:TFT_LEGACY) - Legacy transcription factor target gene sets." = "C3(TFT_TFT_LEGACY)", 
        "C4 - Computational gene sets." = "C4(Computational)", 
        "C4 (3CA) - Curated Cancer Cell Atlas gene sets." = "C4(3CA)", 
        "C4 (CGN) - Cancer Gene Neighborhood gene sets." = "C4(CGN)", 
        "C4 (CM) - Cancer Module gene sets." = "C4(CM)", 
        "C5 (GO:BP) - GO Biological Process ontology gene sets." = "C5(GO_BP)",  
        "C5 (GO:CC) - GO Cellular Component ontology gene sets." = "C5(GO_CC)",  
        "C5 (GO:MF) - GO Molecular Function ontology gene sets." = "C5(GO_MF)", 
        "C5 (HPO) - Human Phenotype Ontology gene sets." = "C5(HPO)", 
        "C6 - Oncogenic signature gene sets." = "C6(Oncogenic)", 
        "C7 (IMMUNESIGDB) - Immunologic signature gene sets." = "C7(ImmuneSigDB)",
        "C7 (VAX) - Vaccine response gene sets." = "C7(VAX)",  
        "C8 - Cell type signature gene sets." = "C8(CellType)",  
        "C9 - Computational perturbation signature gene sets." = "C9(CompPerturb)"
      )
      
    } 
    else {db_choices <- c(
      "M1 - Positional gene sets — Gene sets corresponding to mouse chromosome cytogenetic bands." = "M1(Positional)",
      "M2 (CGP) - Chemical and genetic perturbations — Gene sets represent expression signatures of genetic and chemical perturbations. A number of these gene sets come in pairs: xxx_UP (and xxx_DN) gene set representing genes induced (and repressed) by the perturbation." = "M2(CGP)",
      "M2 (CP:BIOCARTA) - BioCarta subset of Canonical pathways — Classical signaling pathways." = "M2(CP_BioCarta)",
      "M2 (CP:REACTOME) - Reactome subset of Canonical pathways — Curated pathway database." = "M2(CP_Reactome)",
      "M2 (CP:WIKIPATHWAYS) - WikiPathways subset of Cononical pathways — Curated signaling pathways." = "M2(CP_WikiPathways)",
      "M3 (MIRDB) - Motif_miR gene sets — microRNA target predictions." = "M3(MIRDB)",
      "M3 (GTRD) - Motif_TF gene sets — Transcription factor target gene sets." = "M3(GTRD)",
      "M5 (GO:BP) - GO_BP subset of Gene Ontology gene set — Biological Process ontology." = "M5(GO_BP)",
      "M5 (GO:CC) - GO_CC subset of Gene Ontology gene sets — Cellular Component ontology." = "M5(GO_CC)",
      "M5 (GO:MF) - GO_MF subset of Gene Ontology gene sets — Molecular Function ontology." = "M5(GO_MF)",
      "M5 (MPT) - Tumor phenotype ontology — Mouse phenotype tumor gene sets." = "M5(MPT)",
      "M7 - Immunologic signature gene sets — Gene sets that represent cell states and perturbations within the immune system." = "M7(Immune)",
      "M8 - Cell type signature gene sets — Gene sets that contain curated cluster markers for cell types identified in single-cell sequencing studies of mouse tissue." = "M8(Celltype)",
      "MH - Hallmark gene sets— Hallmark gene sets summarize and represent specific well-defined biological states or processes and display coherent expression. These gene sets were generated by a computational methodology based on identifying overlaps between gene sets in other MSigDB collections and retaining genes that display coordinate expression." = "MH(Hallmark)"
    )
    }
    
    updateSelectInput(
      session,
      "gsea_db",
      choices = db_choices,
      selected = unname(db_choices)
    )
    
  })
  
  # ----------------- Cell Cell Communication UI -----------------
  output$cellchat_ui <- renderUI({
    
    if (is.null(rv$data_obj)) {
      tagList(
        h4("Input Data Path"),
        path_input_ui(id = "data", label = "Seurat File Path", 
                      placeholder = "/path/to/seurat_file.rds", 
                      value = rv$data_obj_path),
        h4("Please load data first.")
      )
    } else {
      tagList(
        br(),
        tags$div(
          style = "padding: 8px; background-color: #e8f5e8; border-radius: 4px; border: 1px solid #c8e6c9; display: flex; align-items: center; justify-content: space-between;",
          tags$p(
            style = "margin: 0;",
            "📁 ",
            tags$strong("Loaded: "),
            tags$code(input$data_path)
          ),
          actionButton(
            "reset_data_obj_btn", 
            "Reset (change RDS file)", 
            class = "btn-warning", 
            style = "min-width: 180px; height: 35px; padding: 5px 10px;"
          )
        ),
        br(),
        fluidRow(
          column(
            6,
            actionButton(
              "open_sample_group_modal",
              HTML("Set <b style='color: red;'>sample</b> groups"),
              class = "btn-primary",
              width = "100%",
              style = "margin-top:20px;"
            )
          ),
          column(
            6,
            selectInput(
              "cluster_column",
              HTML("meta.data <b style='color: red;'>celltype</b> cluster identification"),
              choices = extract_meta_columns(rv$data_obj),
              width = "100%"
            )
          )
        ),
        hr(),
        radioButtons(
          "cc_species",
          "Species",
          choices = c("Human" = "human", "Mouse" = "mouse"),
          selected = "mouse",
          inline = TRUE
        ),
        
        checkboxGroupInput(
          "cc_pathway_type",
          "Pathway Type",
          choices = c(
            "Secreted Signaling" = "Secreted Signaling",
            "Cell-Cell Contact" = "Cell-Cell Contact",
            "ECM-Receptor" = "ECM-Receptor",
            "Non-protein Signaling" = "Non-protein Signaling"
          ),
          selected = c("Secreted Signaling", "Cell-Cell Contact", "ECM-Receptor", "Non-protein Signaling"),
          inline = TRUE
        ),
        div(
          style = "display:flex; align-items:center;",
          
          radioButtons(
            "cc_prob_type",
            "Probability Method",
            choices = c(
              "triMean" = "triMean",
              "truncatedMean(0.1)" = "truncatedmean_0.1",
              "truncatedMean(0.05)" = "truncatedmean_0.05"
            ),
            selected = "triMean",
            inline = TRUE
          ),
          
          div(
            style = "margin-left:15px;",
            actionButton(
              "recommend_cc_prob_type",
              "Calculate recommended method",
              class = "btn-primary"
            )
          )
        ),
        numericInput("cc_min_cells", "Min Cells (filter)", value = 10, min = 1),
        actionButton("cc_run_analysis", "Run CellChat Analysis",
                     class = "btn-success",
                     style = "margin-top: 20px;")
      )
    }
  })
  # ----------------- Visualization UI (VIOLIN, DOTPLOT, VOLCANO, NES...) -----------------
  output$viz_panel <- renderUI({
    req(input$viz_tabs)
    tagList(
      switch(input$viz_tabs,
             "violin" = uiOutput("violin_ui"),
             "dotplot" = uiOutput("dotplot_ui"),
             "feature" = uiOutput("feature_ui"),
             "volcano" = uiOutput("volcano_ui"),
             "pathway_volcano" = uiOutput("pathway_volcano_ui"),
             "gsea_nes" = uiOutput("gsea_nes_ui"),
             "cccv" = uiOutput("cccv_ui"),
             "cccv_condition" = uiOutput("cccv_condition_ui")
      ),
      br(),br(),br(),br(),br()
    )
  })
  
  # VIOLIN UI
  output$violin_ui <- renderUI({
    tagList(
      h4("Seurat RDS Obj:"),
      path_input_ui(
        id = "data",
        label = "Seurat File Path",
        placeholder = "/path/to/seurat_file.rds",
        value = rv$data_obj_path),
      h4("DEG result (Optional):"),
      path_input_ui(
        id = "deg",
        label = "DEG File Path (.tsv)",
        placeholder = "/path/to/deg_result.tsv",
        value = rv$deg_obj_path),
      fluidRow(
        column(6, actionButton("marker_select_btn", "Select Genes", class = "btn-success", width = "100%")),
        column(6, actionButton("sample_ident_select_btn", "Select meta.data field", class = "btn-primary", width = "100%"))
      ),
      br(),
      uiOutput("violin_plot_tabs_ui"))
    
  })
  
  # DOTPLOT UI
  output$dotplot_ui <- renderUI({
    tagList(
      h4("Seurat RDS Obj:"),
      path_input_ui(
        id = "data",
        label = "Seurat File Path",
        placeholder = "/path/to/seurat_file.rds",
        value = rv$data_obj_path),
      h4("DEG result (Optional):"),
      path_input_ui(
        id = "deg",
        label = "DEG File Path (.tsv)",
        placeholder = "/path/to/deg_result.tsv",
        value = rv$deg_obj_path),
      fluidRow(
        column(6, actionButton("marker_select_btn", "Select Genes", class = "btn-success", width = "100%")),
        column(6, actionButton("sample_ident_select_btn", "Select meta.data field", class = "btn-primary", width = "100%"))
      ),
      br(), 
      uiOutput("dot_plot_ui")
    )
  })
  
  #feature Plot UI
  output$feature_ui <- renderUI({
    tagList(
      h4("Seurat RDS Obj:"),
      path_input_ui(
        id = "data",
        label = "Seurat File Path",
        placeholder = "/path/to/seurat_file.rds",
        value = rv$data_obj_path),
      h4("DEG result (Optional):"),
      path_input_ui(
        id = "deg",
        label = "DEG File Path (.tsv)",
        placeholder = "/path/to/deg_result.tsv",
        value = rv$deg_obj_path),
      fluidRow(
        column(6, actionButton("marker_select_btn", "Select Genes", class = "btn-success", width = "100%")),
        column(6, actionButton("sample_ident_select_btn", "Select meta.data field", class = "btn-primary", width = "100%"))
      ),
      br(), 
      uiOutput("featureplot_ui")
    )
  })
  
  # ORIGINAL VOLCANO UI
  output$volcano_ui <- renderUI({
    tagList(
      h4("DEG result:"),
      path_input_ui(
        id = "deg",
        label = "DEG File Path (.tsv)",
        placeholder = "/path/to/deg_result.tsv",
        value = rv$deg_obj_path),
      actionButton("marker_select_btn", "Select genes to be labeled (Optional)", class = "btn-success", width = "100%"),
      br(), 
      uiOutput("volcano_plot_ui")
    )
  })
  
  # PATHWAY VOLCANO UI
  output$pathway_volcano_ui <- renderUI({
    tagList(
      h4("DEG result:"),
      path_input_ui(
        id = "deg",
        label = "DEG File Path (.tsv)",
        placeholder = "/path/to/deg_result.tsv",
        value = rv$deg_obj_path),
      h4("GSEA result:"),
      path_input_ui(
        id = "gsea",
        label = "GSEA File Path (.xlsx)",
        placeholder = "/path/to/gsea_file.xlsx",
        value = rv$gsea_obj_path),
      fluidRow(
        column(6, actionButton("pathway_select_btn", "Select Pathway", class = "btn-warning", width = "100%")),
        column(6, 
               actionButton("pathway_rename_btn", "Rename Pathway", width = "100%"),
               tags$style("#pathway_rename_btn { background-color: #967969; color: #ffffff; border-color: #ccc; }"))
      ),
      br(),
      uiOutput("pathway_volcano_plot_ui"))
    
  })
  
  # GSEA NES UI
  output$gsea_nes_ui <- renderUI({
    tagList(
      h4("GSEA result:"),
      path_input_ui(
        id = "gsea",
        label = "GSEA File Path (.xlsx)",
        placeholder = "/path/to/gsea_file.xlsx",
        value = rv$gsea_obj_path),
      fluidRow(
        column(6, actionButton("pathway_select_btn", "Select Pathway", class = "btn-warning", width = "100%")),
        column(6, 
               actionButton("pathway_rename_btn", "Rename Pathway", width = "100%"),
               tags$style("#pathway_rename_btn { background-color: #967969; color: #ffffff; border-color: #ccc; }"))
      ),
      br(), 
      uiOutput("nes_barplot_ui")
    )
  })
  
  # CELL CELL COMMUNICATION UI
  output$cccv_ui <- renderUI({
    tagList(
    h4("Cell Cell Communication result:"),
    path_input_ui(
      id = "cccd",
      label = "Cell Cell Communication File Path (.rds)",
      placeholder = "/path/to/cellchat_compare.rds or /path/to/cellchat_object.rds",
      value = rv$cccd_obj_path),
    br(),
    uiOutput("cccv_control_ui")
    )
  })
  
  # CELL CELL COMMUNICATION - CONDITON UI
  output$cccv_condition_ui <- renderUI({
    tagList(
      h4("Cell Cell Communication result:"),
      path_input_ui(
        id = "cccd",
        label = "Cell Cell Communication File Path (.rds), must be included control and treatment information",
        placeholder = "/path/to/cellchat_compare.rds",
        value = rv$cccd_obj_path),
      br(),
      uiOutput("cccv_condition_control_ui")
    )
  })
  
  # ----------------- Violin plot UI, Check whether data_obj exist -----------------
  output$violin_plot_tabs_ui <- renderUI({
    
    if (is.null(rv$data_obj) | rv$visualization_marker == "" | rv$visualization_sample_column == "") {
      
      tagList(
        div(
          style = "padding: 15px; background-color: #fff3cd; border: 1px solid #ffeeba; border-radius: 5px;",
          "Violin plot requires a loaded Seurat RDS object, gene list (`Select Gene`) and meta.data field (Select meta.data field). Please load data to process"
        )
      )
      
    } else {
      
      tabsetPanel(
        tabPanel("Original Violin Plot", uiOutput("original_violin_ui")),
        tabPanel("Stack Violin Plot", uiOutput("stack_violin_ui"))
      )
      
    }
  })
  
  # ----------------- ORIGINAL VIOLIN UI - VIOLIN UI - violin_plot_tabs_ui -----------------
  output$original_violin_ui <- renderUI({
    req(rv$data_obj)
    req(rv$visualization_marker)
    req(rv$visualization_sample_column != "")
    req(!rv$ui_freeze)
    tagList(
      
      h3("Original Violin Plot"),
      
      hr(),
      div(
        style = "display: flex; justify-content: center;",
        uiOutput("original_violin_plot_ui")),
      
      hr(),
      
      fluidRow(
        column(3, numericInput("v_x_lab_size", "X label size", value = 12)),
        column(3, numericInput("v_y_lab_size", "Y label size", value = 12)),
        column(3, numericInput("v_x_text_size", "X text size", value = 10)),
        column(3, numericInput("v_y_text_size", "Y text size", value = 10))
      ),
      fluidRow(
        column(3, numericInput("v_title_size", "Title size", value = 14)),
        column(3, numericInput("v_angle", "X rotation angle", value = 0)),
        column(3, numericInput("v_width", "Graph width (inch)", value = 4)),
        column(3, numericInput("v_height", "Graph height (inch)", value = 3))
      ),
      fluidRow(
        column(4, checkboxInput("v_show_points", "Show points", value = TRUE)),
        column(4, checkboxInput("v_show_non_zero", "Show non-zero proportion", value = TRUE))
      ),
      hr(),
      textInput(
        "violin_save_path",
        "Violin output directory",
        value = file.path(rv$save_path, "Original_Violin"),
        width ="100%"
      ),
      hr(),
      fluidRow(
        column(6,
               actionButton(
                 "generate_all_violin",
                 "Generate all genes' violin plots and save",
                 style = "background-color: #87cefa; color: white;",
                 width = "100%"
               )
        ),
        column(6,
               downloadButton(
                 "download_violin",
                 "Download Violin Plots to Local",
                 style = "background-color: #4CAF50; color: white; width: 100%;"
               )
        )
      )
    )
  })
  
  # ----------------- STACK VIOLIN UI - VIOLIN UI - violin_plot_tabs_ui -----------------
  output$stack_violin_ui <- renderUI({
    req(rv$data_obj)
    req(rv$visualization_marker)
    req(rv$visualization_sample_column != "")
    req(!rv$ui_freeze)
    tagList(
      h3("Stack Violin Plot"),
      hr(),
      actionButton("set_group_gene", "Set the grouping and color of genes", class = "btn-info", width = "100%"),
      br(),br(),
      div(
        style = "display: flex; justify-content: center;",
        plotOutput("stack_violin_plot", height = "768px", width = "864px")),
      fluidRow(
        # LEFT
        column(3, numericInput("stack_violin_y_left_sample_size", "Left sample text size", value = 2.5)),
        column(3, numericInput("stack_violin_y_left_gene_size", "Left gene size", value = 8)),
        column(3, numericInput("stack_violin_y_left_label_size", "Left axis label size", value = 8)),
        # RIGHT
        column(3, numericInput("stack_violin_x_size", "X axis sample size", value = 8)),
        column(3, numericInput("stack_violin_y_right_value_size", "Right Y-axis value size", value = 8)),
        column(3, numericInput("stack_violin_y_right_label_size", "Right Y-axis label size", value = 8)),
        # LAYOUT
        column(6, sliderInput("stack_violin_left_proportion", "Left panel proportion", min = 0.001, max = 0.999, value = 0.3, width = "100%")),
        column(3, numericInput("stack_width", "Graph width (inch)", value = 9)),
        column(3, numericInput("stack_height", "Graph height (inch)", value = 8))
      ),
      
      fluidRow(
        column(12,
               textInput(
                 "stack_save_path",
                 "Violin output directory",
                 value = file.path(rv$save_path, "Stack_Violin"),
                 width ="100%"
               ),
        ),
        column(6,
               actionButton("save_stack_plot", "Save Graph on SCC",
                            style = "background-color: #87cefa; color: white; width: 100%;")
        ),
        column(6,
               downloadButton("download_stack_plot", "Download Violin Plot to Local",
                              style = "background-color: #4CAF50; color: white; width: 100%;")
        )
      )
      
    )
  })
  
  # ----------------- DOTPLOT UI -----------------
  output$dot_plot_ui <- renderUI({
    if (is.null(rv$data_obj) | rv$visualization_marker == "" | rv$visualization_sample_column == "") {
      tagList(
        div(
          style = "padding: 15px; background-color: #fff3cd; border: 1px solid #ffeeba; border-radius: 5px;",
          "Dotplot requires a loaded Seurat RDS object, gene list (`Select Gene`) and meta.data field (Select meta.data field). Please load data to process."
        )
      )
    }else{
      tagList(
        div(
          style = "display: flex; justify-content: center;",
          uiOutput("dotplot_plot_ui")),
        fluidRow(
          column(3, numericInput("dotplot_x_text_size", "X text size", value = 10)),
          column(3, numericInput("dotplot_y_text_size", "Y text size", value = 10)),
          column(3, numericInput("dotplot_dot_max_size", "Max dot size", value = 6)),
          column(3, numericInput("dotplot_x_rotate", "X test rotate degree", value = 45))
        ),
        fluidRow(
          column(3,
                 selectInput(
                   "dotplot_legend_position",
                   "Legend position",
                   choices = c("right", "left", "top", "bottom", "none"),
                   selected = "right"
                 )
          ),
          column(3, numericInput("dotplot_legend_text_size", "Legend text size", value = 10)),
          column(3, numericInput("dotplot_legend_title_size", "Legend title size", value = 11)),
          column(3, numericInput("dotplot_legend_key_size", "Key size", value = 0.5))
        ),
        fluidRow(
          column(3, numericInput("dotplot_width", "Graph width (inch)", value = 5)),
          column(3, numericInput("dotplot_height", "Graph height (inch)", value = 3)),
          column(6, div(style = "padding-top: 22px;", checkboxInput("dotplot_scale", "Scale the color of dotplot (average expression)", value = TRUE, width = "100%")))
        ),
        fluidRow(
          column(12,
                 textInput(
                   "dotplot_save_path",
                   "Dotplot output directory",
                   value = file.path(rv$save_path, "Dotplot"),
                   width ="100%"
                 ),
          ),
          column(6,
                 actionButton("save_dotplot", "Save Graph on SCC",
                              style = "background-color: #87cefa; color: white; width: 100%;")
          ),
          column(6,
                 downloadButton("download_dotplot", "Download Dot Plot to Local",
                                style = "background-color: #4CAF50; color: white; width: 100%;")
          )
        )
      )
    }
  })
  # ----------------- FEATURE PLOT UI -----------------
  output$featureplot_ui <- renderUI({
    if (is.null(rv$data_obj) | rv$visualization_marker == "" | rv$visualization_sample_column == "") {
      tagList(
        div(
          style = "padding: 15px; background-color: #fff3cd; border: 1px solid #ffeeba; border-radius: 5px;",
          "Feature plot requires a loaded Seurat RDS object, gene list (`Select Gene`) and meta.data field (Select meta.data field). Please load data to process."
        )
      )
    }else{
      tagList(
        div(
          style = "display: flex; justify-content: center;",
          plotOutput("feature_plot", height = "360px")),
        fluidRow(
          column(3, numericInput("feature_plot_point_size", "Point size", value = 0.5)),
          column(3, numericInput("feature_plot_title_size", "Title size", value = 15)),
          column(3, numericInput("feature_plot_legend_text_size", "Legend text size", value = 12)),
          column(
            3,
            selectInput(
              "feature_plot_legend_position",
              "Legend position",
              choices = c("right", "left", "top", "bottom", "none"),
              selected = "right"
            )
          )
        ),
        fluidRow(
          column(3, numericInput("feature_plot_title_size", "Axis title size", value = 12)),
          column(3, numericInput("feature_plot_text_size", "Axis text size", value = 10)),
          column(3, numericInput("feature_plot_width", "Graph width (inch)", value = 10)),
          column(3, numericInput("feature_plot_height", "Graph height (inch)", value = 3))
        ),
        hr(),
        textInput(
          "feature_plot_save_path",
          "Feature Plot output directory",
          value = file.path(rv$save_path, "Feature_Plot"),
          width ="100%"
        ),
        hr(),
        fluidRow(
          column(6,
                 actionButton(
                   "generate_all_feature_plot",
                   "Generate all genes' feature plots and save",
                   style = "background-color: #87cefa; color: white;",
                   width = "100%"
                 )
          ),
          column(6,
                 downloadButton(
                   "download_feature_plot",
                   "Download Feature Plot to Local",
                   style = "background-color: #4CAF50; color: white; width: 100%;"
                 )
          )
        )
        
        )
    }
  })
  # ----------------- Original Volcano UI -----------------
  output$volcano_plot_ui <- renderUI({
    
    if (is.null(rv$deg_obj)) {
      return(
        tagList(
          div(
            style = "padding: 15px; background-color: #fff3cd; border: 1px solid #ffeeba; border-radius: 5px;",
            "Volcano Plot (Genes) requires a loaded DEG dataset. Please load the data to proceed."
          )
        )
      )
    }else{
      tagList(
        div(
          style = "display: flex; justify-content: center;",
          plotOutput("volcano_plot", 
                     height = "307.2px",   
                     width = "360px")    
        ),
        fluidRow(
          column(3, numericInput("original_volcano_ylimit", "Y limit", value = 150)),
          column(3, numericInput("original_volcano_xlimit", "X limit", value = 10)),
          column(3, numericInput("original_volcano_pvalue_limit", "P-value cutoff", value = 0.05)),
          column(3, numericInput("original_volcano_log2FC_limit", "log2FC cutoff", value = 0.25))
        ),
        fluidRow(
          column(3, numericInput("original_volcano_size_point", "Point size", value = 0.8)),
          column(3, numericInput("original_volcano_size_label", "Label size", value = 3)),
          column(3, numericInput("original_volcano_text_size", "Axis text size", value = 10)),
          column(3, numericInput("original_volcano_title_size", "Axis title size", value = 10))
        ),
        fluidRow(
          column(3, numericInput("original_volcano_width", "Graph width (inch)", value = 3.6)),
          column(3, numericInput("original_volcano_height", "Graph height (inch)", value = 3.2)),
          column(6, div(style = "padding-top: 22px;", checkboxInput("original_volcano_auto_label", "Automatically label gene based on volcano graph", value = FALSE, width = "100%")))
        ),
        fluidRow(
          column(12,
                 textInput(
                   "original_volcano_save_path",
                   "Original Volcano Plot output directory",
                   value = file.path(rv$save_path, "Original_Volcano"),
                   width ="100%"
                 ),
          ),
          column(6,
                 actionButton("save_original_volcano_plot", "Save Graph on SCC",
                              style = "background-color: #87cefa; color: white; width: 100%;")
          ),
          column(6,
                 downloadButton("download_original_volcano_plot", "Download Volcano Plot to local",
                                style = "background-color: #4CAF50; color: white; width: 100%;")
          )
        )
      )}
  })
  # ----------------- Pathway Volcano UI -----------------
  output$pathway_volcano_plot_ui <- renderUI({
    
    if (is.null(rv$deg_obj) | is.null(rv$gsea_obj) | rv$visualization_pathway == "") {
      
      tagList(
        div(
          style = "padding: 15px; background-color: #fff3cd; border: 1px solid #ffeeba; border-radius: 5px;",
          "Volcano Plot (Pathways) requires a loaded DEG dataset, GSEA dataset and selected pathways. Please load the data and select pathway to proceed."
        )
      )
      
    } else {
      
      tabsetPanel(
        tabPanel("Generate multi pathways on separate plot", uiOutput("pathway_volcano_multi_ui")),
        tabPanel("Generate multi pathways on one plot", uiOutput("pathway_volcano_single_ui"))
      )
    }
  })
  # ----------------- Pathway Volcano UI (Multi)-----------------
  output$pathway_volcano_multi_ui <- renderUI({
    tagList(
      h5("Example:"),
      div(
        style = "display: flex; justify-content: center;",
        plotOutput(
          "pathway_volcano_multi_plot",
          height = "500px",
          width = "720px"
        )
      ),
      
      fluidRow(
        column(3, numericInput("pathway_volcano_pv_limit", "P value limit", value = 0.05)),
        column(3, numericInput("pathway_volcano_log2fc_limit", "log2FC limit", value = 0.25)),
        column(3, numericInput("pathway_volcano_point_size", "Point size", value = 5)),
        column(3,
               selectInput(
                 "pathway_volcano_legend_position",
                 "Legend position",
                 choices = c("right", "left", "top", "bottom", "none"),
                 selected = "right"
               )
        )
      ),
      fluidRow(
        column(3, numericInput("pathway_volcano_label_size", "Label size", value = 4)),
        column(3, numericInput("pathway_volcano_axis_title_size", "Axis title size", value = 10)),
        column(3, numericInput("pathway_volcano_axis_text_size", "Axis text size", value = 10)),
        column(3, numericInput("pathway_volcano_legend_text_size", "Legend text size", value = 10))
      ),
      
      fluidRow(
        column(3, numericInput("pathway_volcano_label_n", "The number of labels (based on top log2FC)", value = 10)),
        column(3, numericInput("pathway_volcano_pdf_width", "Graph width (inch)", value = 4)),
        column(3, numericInput("pathway_volcano_pdf_height", "Graph height (inch)", value = 4))
      ),
      
      fluidRow(
        column(
          12,
          textInput(
            "pathway_volcano_save_path",
            "Output directory",
            value = file.path(rv$save_path, "Pathway_Volcano"),
            width = "100%"
          )
        )
      ),
      
      fluidRow(
        column(
          6,
          actionButton(
            "save_pathway_multi_volcano",
            "Save Graph on SCC",
            style = "background-color: #87cefa; color: white; width: 100%;"
          )
        ),
        column(
          6,
          downloadButton(
            "download_pathway_multi_volcano",
            "Download Plot to Local",
            style = "background-color: #4CAF50; color: white; width: 100%;"
          )
        )
      )
    )
  })
  
  # ----------------- Pathway Volcano UI (Single) -----------------
  output$pathway_volcano_single_ui <- renderUI({
    tagList(
      
      div(
        style = "display: flex; justify-content: center;",
        plotOutput(
          "pathway_volcano_plot",
          height = "500px",
          width = "720px"
        )
      ),
      
      fluidRow(
        column(3, numericInput("pathway_volcano_pv_limit", "P value limit", value = 0.05)),
        column(3, numericInput("pathway_volcano_log2fc_limit", "log2FC limit", value = 0.25)),
        column(3, numericInput("pathway_volcano_point_size", "Point size", value = 5)),
        column(3,
               selectInput(
                 "pathway_volcano_legend_position",
                 "Legend position",
                 choices = c("right", "left", "top", "bottom", "none"),
                 selected = "right"
               )
        )
      ),
      fluidRow(
        column(3, numericInput("pathway_volcano_label_size", "Label size", value = 4)),
        column(3, numericInput("pathway_volcano_axis_title_size", "Axis title size", value = 10)),
        column(3, numericInput("pathway_volcano_axis_text_size", "Axis text size", value = 10)),
        column(3, numericInput("pathway_volcano_legend_text_size", "Legend text size", value = 10))
      ),
      
      fluidRow(
        column(3, numericInput("pathway_volcano_label_n", "The number of labels (based on top log2FC)", value = 10)),
        column(3, numericInput("pathway_volcano_pdf_width", "Graph width (inch)", value = 6)),
        column(3, numericInput("pathway_volcano_pdf_height", "Graph height (inch)", value = 8))
      ),
      
      fluidRow(
        column(
          12,
          textInput(
            "pathway_volcano_save_path",
            "Output directory",
            value = file.path(rv$save_path, "Pathway_Volcano"),
            width = "100%"
          )
        )
      ),
      
      fluidRow(
        column(
          6,
          actionButton(
            "save_pathway_volcano",
            "Save Graph on SCC",
            style = "background-color: #87cefa; color: white; width: 100%;"
          )
        ),
        column(
          6,
          downloadButton(
            "download_pathway_volcano",
            "Download Plot to Local",
            style = "background-color: #4CAF50; color: white; width: 100%;"
          )
        )
      )
    )
    })
  # ----------------- NES Barplot UI -----------------
  output$nes_barplot_ui <- renderUI({
    
    if (is.null(rv$gsea_obj) | rv$visualization_pathway == "") {
      
      tagList(
        div(
          style = "padding: 15px; background-color: #fff3cd; border: 1px solid #ffeeba; border-radius: 5px;",
          "NES Barplot requires a loaded GSEA dataset and selected pathways. Please load the data and select pathway to proceed."
        )
      )
      
    } else {
      tagList(
        
        div(
          style = "display: flex; justify-content: center;",
          plotOutput(
            "nes_barplot",
            height = "600px",
            width = "600px"
          )
        ),
        fluidRow(
          column(3, numericInput("nes_barplot_pathway_size", "y axis text (pathway) size", value = 12)),
          column(3, numericInput("nes_barplot_y_title_size", "y axis title size", value = 10)),
          column(3, numericInput("nes_barplot_nes_value_size", "x axis text (NES value) size", value = 12)),
          column(3, numericInput("nes_barplot_x_title_size", "x axis title size", value = 10)),
        ),
        fluidRow(
          column(3, numericInput("nes_barplot_pathway_length", "Pathway line break length", value = 40)),
          column(3, numericInput("nes_barplot_bar_width", "Bar width", value = 0.7)),
          column(3, numericInput("nes_barplot_pdf_width", "Graph width (inch)", value = 6)),
          column(3, numericInput("nes_barplot_pdf_height", "Graph height (inch)", value = 6)),
        ),
        fluidRow(
          column(
            12,
            textInput(
              "nes_barplot_save_path",
              "Output directory",
              value = file.path(rv$save_path, "NES_barplot"),
              width = "100%"
            )
          )
        ),
        
        fluidRow(
          column(
            6,
            actionButton(
              "save_nes_barplot",
              "Save Graph on SCC",
              style = "background-color: #87cefa; color: white; width: 100%;"
            )
          ),
          column(
            6,
            downloadButton(
              "download_nes_barplot",
              "Download Plot to Local",
              style = "background-color: #4CAF50; color: white; width: 100%;"
            )
          )
        )
        )
      }})
  
  # ----------------- Cell Cell Communication Control UI -----------------
  output$cccv_control_ui <- renderUI({
    if (is.null(rv$cccd_obj)) {
        div(
          style = "padding: 15px; background-color: #fff3cd; border: 1px solid #ffeeba; border-radius: 5px;",
          "Please upload the cell cell communication result from `Cell-Cell Communication Analysis` tab."
        )
      
    } else {
      tabsetPanel(
        tabPanel("Cluster Level", uiOutput("ccc_cluster_control_ui")),
        tabPanel("Pathway Level", uiOutput("ccc_pathway_control_ui")),
        tabPanel("L-R Pair Level", uiOutput("ccc_pair_control_ui"))
      )
    }
})
  
  output$cccv_condition_control_ui <- renderUI({
    if (is.null(rv$cccd_obj) | is.null(rv$cccd_obj$data_merge)){
      div(
        style = "padding: 15px; background-color: #fff3cd; border: 1px solid #ffeeba; border-radius: 5px;",
        "Please upload the cell cell communication result from `Cell-Cell Communication Analysis` tab. Also make sure the upload file includes control treatment information (cellchat_compare.rds)"
      )
      
    } else {
      tabsetPanel(
        tabPanel("Overall Level", uiOutput("cccc_overall_control_ui")),
        tabPanel("Cluster Level", uiOutput("cccc_cluster_control_ui")),
        tabPanel("Pathway Level", uiOutput("cccc_pathway_control_ui"))
      )
    }
  })
  # ----------------- Cell Cell Communication UI - cluster-----------------
  output$ccc_cluster_control_ui <- renderUI({
    tabsetPanel(
      
      tabPanel(
        "Communication Network",
        uiOutput("netVisual_circle_ui")# netVisual_circle
      ),
      
      tabPanel(
        "Sender vs Receiver Roles",
        uiOutput("netAnalysis_signalingRole_scatter_ui")# netAnalysis_signalingRole_scatter
      ),
      
      tabPanel(
        "Pathway Role Heatmap",
        uiOutput("netAnalysis_signalingRole_heatmap_ui")# netAnalysis_signalingRole_heatmap
      )
      
    )
  })
  
  output$netVisual_circle_ui <- renderUI({
    tagList(
      div(
        style = "display: flex; justify-content: center;",
        uiOutput("netVisual_circle_control_ui")
      ),
      
      fluidRow(
        column(
          3,
          selectInput(
            "circle_measure",
            "Network measure",
            choices = c(
              "Interaction count" = "count",
              "Interaction weight" = "weight"
            ),
            selected = "weight"
          )),
        column(3, numericInput("circle_vertex_label_size", "Cell group label size", value = 1.0)),
        column(3, numericInput("circle_vertex_size", "Node size", value = 3)),
        column(3, numericInput("circle_edge_width_max", "Max line width", value = 10))
      ),
      
      fluidRow(
        column(3, numericInput("circle_arrow_size", "Arrow size", value = 0.5)),
        column(3, numericInput("circle_title_size", "Title size", value = 1.2)),
        column(3, div(style = "padding-top: 22px;", checkboxInput("circle_show_edge_label", "Show edge label", value = FALSE)))
      ),
      
      fluidRow(
        column(3, numericInput("circle_width", "PDF width (inch)", value = if(is.null(rv$cccd_obj$data_merge)){4}else{8})),
        column(3, numericInput("circle_height", "PDF height (inch)", value = 4))
    ),
    fluidRow(
      column(
        12,
        textInput(
          "circle_save_path",
          "Output directory",
          value = file.path(rv$save_path, "CellCellCommunication"),
          width = "100%"
        )
      )
    ),
    
    fluidRow(
      column(
        6,
        actionButton(
          "circle_barplot",
          "Save Graph on SCC",
          style = "background-color: #87cefa; color: white; width: 100%;"
        )
      ),
      column(
        6,
        downloadButton(
          "download_circle",
          "Download Plot to Local",
          style = "background-color: #4CAF50; color: white; width: 100%;"
        )
      )
    )
    )
  })
  
  output$netAnalysis_signalingRole_scatter_ui <- renderUI({
    tagList(
      div(
        style = "display: flex; justify-content: center;",
        uiOutput("netAnalysis_signalingRole_scatter_control_ui")
      ),
      
      fluidRow(
        column(3, numericInput("role_scatter_label_size", "Label size", value = 8)),
        column(3, numericInput("role_scatter_point_size", "Point size", value = 5)),
        column(3, numericInput("role_scatter_title_size", "Title size", value = 14)),
        column(3, numericInput("role_scatter_axis_text_size", "Axis text size", value = 10))
      ),
      
      fluidRow(
        column(3, numericInput("role_scatter_axis_title_size", "Axis title size", value = 12)),
        column(3, numericInput("role_scatter_legend_text_size", "Legend text size", value = 10)),
        column(3, numericInput("role_scatter_legend_title_size", "Legend title size", value = 11))
      ),
      
      fluidRow(
        column(3, numericInput("role_scatter_width", "PDF width (inch)", value = if(is.null(rv$cccd_obj$data_merge)){4}else{8})),
        column(3, numericInput("role_scatter_height", "PDF height (inch)", value = 3))
      ),
      
      fluidRow(
        column(
          12,
          textInput(
            "role_scatter_save_path",
            "Output directory",
            value = file.path(rv$save_path, "CellCellCommunication"),
            width = "100%"
          )
        )
      ),
      
      fluidRow(
        column(
          6,
          actionButton(
            "save_role_scatter",
            "Save Graph on SCC",
            style = "background-color: #87cefa; color: white; width: 100%;"
          )
        ),
        column(
          6,
          downloadButton(
            "download_role_scatter",
            "Download Plot to Local",
            style = "background-color: #4CAF50; color: white; width: 100%;"
          )
        )
      )
    )
  })
  
  output$netAnalysis_signalingRole_heatmap_ui <- renderUI({
    
    tagList(
      
      div(
        style = "display: flex; justify-content: center;",
        uiOutput("netAnalysis_signalingRole_heatmap_control")
      ),
      
      fluidRow(
        column(3,
               selectInput(
                 "netAnalysis_signalingRole_heatmap_pattern",
                 "Role pattern",
                 choices = c(
                   "Outgoing" = "outgoing",
                   "Incoming" = "incoming",
                   "Both" = "both"
                 ),
                 selected = "both"
               )
        ),
        column(3,
               selectInput(
                 "netAnalysis_signalingRole_heatmap_color",
                 "Heatmap color",
                 choices = c(
                   "Reds" = "Reds",
                   "Blues" = "Blues",
                   "Greens" = "Greens",
                   "Purples" = "Purples",
                   "Oranges" = "Oranges"
                 ),
                 selected = "Greens"
               )
        ),
        column(3, numericInput("netAnalysis_signalingRole_heatmap_axis_title_size", "Axis title size", value = 10)),
        column(3, numericInput("netAnalysis_signalingRole_heatmap_axis_text_size", "Axis text size", value = 8))
      ),
      
      fluidRow(
        column(3, numericInput("netAnalysis_signalingRole_heatmap_legend_size", "Legend size", value = 10)),
        column(3, numericInput("netAnalysis_signalingRole_heatmap_width", "PDF width (inch)", value = if (is.null(rv$cccd_obj$data_merge)) {6} else {12})),
        column(3, numericInput("netAnalysis_signalingRole_heatmap_height", "PDF height (inch)", value = 10))
      ),
      
      fluidRow(
        column(
          12,
          textInput(
            "netAnalysis_signalingRole_heatmap_save_path",
            "Output directory",
            value = file.path(rv$save_path, "CellCellCommunication"),
            width = "100%"
          )
        )
      ),
      
      fluidRow(
        column(
          6,
          actionButton(
            "save_netAnalysis_signalingRole_heatmap",
            "Save Graph on SCC",
            style = "background-color: #87cefa; color: white; width: 100%;"
          )
        ),
        column(
          6,
          downloadButton(
            "download_netAnalysis_signalingRole_heatmap",
            "Download Plot to Local",
            style = "background-color: #4CAF50; color: white; width: 100%;"
          )
        )
      )
    )
  })
  # ----------------- Cell Cell Communication UI - pathway-----------------
  output$ccc_pathway_control_ui <- renderUI({
    tagList(
      actionButton("ccc_pathway_select_btn", "Select Cell Cell Communication pathway to display", class = "btn-success", width = "100%"),
      uiOutput("ccc_pathway_result_ui")
      )
  })
  
  output$ccc_pathway_result_ui <- renderUI({
    
    tabsetPanel(
      tabPanel(
        "Pathway Network",
        uiOutput("netVisual_aggregate_ui")
      ),
      tabPanel(
        "Pathway Heatmap",
        uiOutput("netVisual_heatmap_ui")
      ),
      tabPanel(
        "Pathway Role Network",
        uiOutput("netAnalysis_signalingRole_network_ui")
      ),
      tabPanel(
        "L-R Contribution",
        uiOutput("netAnalysis_contribution_ui")
      )
    )
  })
  
  output$netVisual_aggregate_ui <- renderUI({
    pathway <- rv$ccc_pathway_select
    
    if (is.null(pathway) || pathway == "") {
      return(
        div(
          style = "padding:15px;background-color:#fff3cd;border:1px solid #ffeeba;border-radius:5px;",
          strong("Pathway required"),
          br(),
          "Please select and confirm a pathway first."
        )
      )
    }
    
    tagList(
      div(
        style = "display: flex; justify-content: center;",
        uiOutput("netVisual_aggregate_control")
      ),
      
      fluidRow(
        column(
          3,
          selectInput(
            "aggregate_layout",
            "Layout",
            choices = c(
              "Circle" = "circle",
              "Hierarchy" = "hierarchy",
              "Chord" = "chord"
            ),
            selected = "circle"
          )
        ),
        column(3, numericInput("aggregate_vertex_label_size", "Label size", value = 1)),
        column(3, numericInput("aggregate_edge_width_max", "Max line width", value = 8)),
        column(3, numericInput("aggregate_arrow_size", "Arrow size", value = 0.5))
      ),
      
      conditionalPanel(
        condition = "input.aggregate_layout == 'hierarchy'",
        fluidRow(
          column(
            6,
            textInput(
              "aggregate_vertex_receiver",
              "Receiver cell index for hierarchy",
              value = "1,2",
              placeholder = "e.g. 1 or 1,2,3, at least two value"
            )
          ),
          column(
            6,
            helpText("Hierarchy layout requires receiver cell indices, e.g. 1,2,3.")
          )
        )
      ),
      
      fluidRow(
        column(3, numericInput("aggregate_width", "PDF width (inch)", value = if(is.null(rv$cccd_obj$data_merge)){6}else{12})),
        column(3, numericInput("aggregate_height", "PDF height (inch)", value = 6))
      ),
      
      fluidRow(
        column(
          12,
          textInput(
            "aggregate_save_path",
            "Output directory",
            value = file.path(rv$save_path, "CellCellCommunication"),
            width = "100%"
          )
        )
      ),
      
      fluidRow(
        column(
          6,
          actionButton(
            "save_netVisual_aggregate",
            "Save Graph on SCC",
            style = "background-color: #87cefa; color: white; width: 100%;"
          )
        ),
        column(
          6,
          downloadButton(
            "download_netVisual_aggregate",
            "Download Plot to Local",
            style = "background-color: #4CAF50; color: white; width: 100%;"
          )
        )
      )
    )
  })
  
  output$netVisual_heatmap_ui <- renderUI({
    pathway <- rv$ccc_pathway_select
    
    if (is.null(pathway) || pathway == "") {
      return(
        div(
          style = "padding:15px;background-color:#fff3cd;border:1px solid #ffeeba;border-radius:5px;",
          strong("Pathway required"),
          br(),
          "Please select and confirm a pathway first."
        )
      )
    }
    
    tagList(
      div(
        style = "display: flex; justify-content: center;",
        uiOutput("netVisual_heatmap_control")
      ),
      
      fluidRow(
        column(3,
          selectInput(
            "netVisual_heatmap_color",
            "Heatmap color",
            choices = c(
              "Reds" = "Reds",
              "Blues" = "Blues",
              "Greens" = "Greens",
              "Purples" = "Purples",
              "Oranges" = "Oranges"
            ), selected = "Reds"
          )),
        column(3, numericInput("netVisual_heatmap_axis_title_size", "Axis title size", value = 10)),
        column(3, numericInput("netVisual_heatmap_axis_text_size", "Axis text size", value = 8))
        ),
      
      fluidRow(
        column(3, numericInput("netVisual_heatmap_legend_size", "Legend size", value = 10)),
        column(3, numericInput("netVisual_heatmap_width", "PDF width (inch)", value = if(is.null(rv$cccd_obj$data_merge)){5}else{10})),
        column(3, numericInput("netVisual_heatmap_height", "PDF height (inch)", value = 4))
      ),
      
      fluidRow(
        column(
          12,
          textInput(
            "netVisual_heatmap_save_path",
            "Output directory",
            value = file.path(rv$save_path, "CellCellCommunication"),
            width = "100%"
          )
        )
      ),
      
      fluidRow(
        column(
          6,
          actionButton(
            "save_netVisual_heatmap",
            "Save Graph on SCC",
            style = "background-color: #87cefa; color: white; width: 100%;"
          )
        ),
        column(
          6,
          downloadButton(
            "download_netVisual_heatmap",
            "Download Plot to Local",
            style = "background-color: #4CAF50; color: white; width: 100%;"
          )
        )
      )
    )
  })
  
  output$netAnalysis_signalingRole_network_ui <- renderUI({
    pathway <- rv$ccc_pathway_select
    
    if (is.null(pathway) || pathway == "") {
      return(
        div(
          style = "padding:15px;background-color:#fff3cd;border:1px solid #ffeeba;border-radius:5px;",
          strong("Pathway required"),
          br(),
          "Please select and confirm a pathway first."
        )
      )
    }
    
    tagList(
      
      div(
        style = "display: flex; justify-content: center;",
        uiOutput("netAnalysis_signalingRole_network_control")
      ),
      
      fluidRow(
        column(3,
          selectInput(
            "netAnalysis_signalingRole_network_color",
            "Heatmap color",
            choices = c(
              "Reds" = "Reds",
              "Blues" = "Blues",
              "Greens" = "Greens",
              "Purples" = "Purples",
              "Oranges" = "Oranges",
              "BuGn" = "BuGn"
            ),
            selected = "BuGn"
          )
        ),
        column(3, numericInput("netAnalysis_signalingRole_network_font_size", "Font size", value = 12)),
        column(3, numericInput("netAnalysis_signalingRole_network_title_size", "Title size", value = 12))
      ),
      fluidRow(
        column(3, numericInput("netAnalysis_signalingRole_network_width", "PDF width (inch)", value = if (is.null(rv$cccd_obj$data_merge)) {8} else {16})),
        column(3, numericInput("netAnalysis_signalingRole_network_height", "PDF height (inch)", value = 4))
      ),
      
      fluidRow(
        column(
          12,
          textInput(
            "netAnalysis_signalingRole_network_save_path",
            "Output directory",
            value = file.path(rv$save_path, "CellCellCommunication"),
            width = "100%"
          )
        )
      ),
      
      fluidRow(
        column(
          6,
          actionButton(
            "save_netAnalysis_signalingRole_network",
            "Save Graph on SCC",
            style = "background-color: #87cefa; color: white; width: 100%;"
          )
        ),
        column(
          6,
          downloadButton(
            "download_netAnalysis_signalingRole_network",
            "Download Plot to Local",
            style = "background-color: #4CAF50; color: white; width: 100%;"
          )
        )
      )
    )
  })
  
  output$netAnalysis_contribution_ui <- renderUI({
    pathway <- rv$ccc_pathway_select
    
    if (is.null(pathway) || pathway == "") {
      return(
        div(
          style = "padding:15px;background-color:#fff3cd;border:1px solid #ffeeba;border-radius:5px;",
          strong("Pathway required"),
          br(),
          "Please select and confirm a pathway first."
        )
      )
    }
    tagList(
      
      div(
        style = "display: flex; justify-content: center;",
        uiOutput("netAnalysis_contribution_control")
      ),
      fluidRow(
        column(3, numericInput("netAnalysis_contribution_x_text_size", "X-axis text size", value = 12)),
        column(3, numericInput("netAnalysis_contribution_lr_pair_size", "L-R pair text size", value = 12)),
        column(3, numericInput("netAnalysis_contribution_title_size", "Title size", value = 10)),
        column(3, div(style = "padding-top: 22px;", checkboxInput("netAnalysis_contribution_show_value", "Show bar values", value = FALSE)))
      ),
      fluidRow(
        column(3, numericInput("netAnalysis_contribution_width", "PDF width (inch)", value = if (is.null(rv$cccd_obj$data_merge)) {6} else {12})),
        column(3, numericInput("netAnalysis_contribution_height", "PDF height (inch)", value = 5))
      ),
      
      fluidRow(
        column(
          12,
          textInput(
            "netAnalysis_contribution_save_path",
            "Output directory",
            value = file.path(rv$save_path, "CellCellCommunication"),
            width = "100%"
          )
        )
      ),
      
      fluidRow(
        column(
          6,
          actionButton(
            "save_netAnalysis_contribution",
            "Save Graph on SCC",
            style = "background-color: #87cefa; color: white; width: 100%;"
          )
        ),
        column(
          6,
          downloadButton(
            "download_netAnalysis_contribution",
            "Download Plot to Local",
            style = "background-color: #4CAF50; color: white; width: 100%;"
          )
        )
      )
    )
  })
  # ----------------- Cell Cell Communication UI - LRpair-----------------
  output$ccc_pair_control_ui <- renderUI({
    
    tagList(
      actionButton(
        "select_lr_pair",
        "Select L-R pair for analysis",
        class = "btn-success",
        style = "width: 100%;"
      ),
      
      if (!is.null(rv$ccc_lr_pair_select) &&
          length(rv$ccc_lr_pair_select) > 0 &&
          any(trimws(rv$ccc_lr_pair_select) != "")) {
        uiOutput("netVisual_individual_ui")
      } else {
        
        div(
          style = "padding: 15px; background-color: #fff3cd; border: 1px solid #ffeeba; border-radius: 5px;",
          "Please select one or more L-R pairs first."
        )
        
      }
    )
  })
  
  output$netVisual_individual_ui <- renderUI({
    
    tagList(
      
      div(
        style = "display: flex; justify-content: center;",
        uiOutput("netVisual_individual_control")
      ),
      
      fluidRow(
        column(3, numericInput("individual_vertex_label_size", "Vertex label size", value = 1)),
        column(3, numericInput("individual_edge_width_max", "Max edge width", value = 8)),
        column(3, numericInput("individual_arrow_size", "Arrow size", value = 0.5))
      ),
      
      fluidRow(
        column(
          3,
          numericInput(
            "individual_width",
            "PDF width (inch)",
            value = if (is.null(rv$cccd_obj$data_merge)) 6 else 12
          )
        ),
        column(3, numericInput("individual_height", "PDF height (inch)", value = 6))
      ),
      
      fluidRow(
        column(
          12,
          textInput(
            "individual_save_path",
            "Output directory",
            value = file.path(rv$save_path, "CellCellCommunication"),
            width = "100%"
          )
        )
      ),
      
      fluidRow(
        column(
          6,
          actionButton(
            "generate_netVisual_individual",
            "Generate Graph to Save Folder",
            style = "background-color: #87cefa; color: white; width: 100%;"
          )
        ),
        column(
          6,
          downloadButton(
            "download_netVisual_individual",
            "Download Generated Graphs",
            style = "background-color: #4CAF50; color: white; width: 100%;"
          )
        )
      )
    )
  })
  # ----------------- Cell Cell Communication Comparison UI - overall-----------------
  output$cccc_overall_control_ui <- renderUI({
    
    tagList(
      
      h4("Overall Communication Comparison"),
      
      div(
        style = "display: flex; justify-content: center;",
        uiOutput("compareInteractions_control")
      ),
      
      br(),
      
      fluidRow(
        column(3, numericInput("compareInteractions_font_size", "Font size", value = 12)),
        column(3, numericInput("compareInteractions_title_size", "Title size", value = 12)),
        column(3, numericInput("compareInteractions_width", "PDF width (inch)", value = 6)),
        column(3, numericInput("compareInteractions_height", "PDF height (inch)", value = 4))
      ),
      
      fluidRow(
        column(
          12,
          textInput(
            "compareInteractions_save_path",
            "Output directory",
            value = file.path(rv$save_path, "CellCellCommunication"),
            width = "100%"
          )
        )
      ),
      
      fluidRow(
        column(
          6,
          actionButton(
            "save_compareInteractions",
            "Generate Graph to Save Folder",
            style = "background-color:#87CEFA;color:white;width:100%;"
          )
        ),
        column(
          6,
          downloadButton(
            "download_compareInteractions",
            "Download Plot to Local",
            style = "background-color:#4CAF50;color:white;width:100%;"
          )
        )
      )
    )
  })
  # ----------------- Cell Cell Communication Comparison UI - cluster-----------------
  output$cccc_cluster_control_ui <- renderUI({
    tabsetPanel(
      tabPanel(
        "Communication Comparison",
        uiOutput("netVisual_diffInteraction_ui")
      ),
        tabPanel(
      "Communication Comparison Heatmap",
      uiOutput("netVisual_relative_heatmap_ui")
    ))
  })
   output$netVisual_diffInteraction_ui <- renderUI({

    tagList(
      
      div(
        style = "display: flex; justify-content: center;",
        uiOutput("cccc_cluster_compare_control")
      ),
      div(
            style = "
        display: flex;
        justify-content: center;
        align-items: center;
        gap: 30px;
        margin-top: 8px;
        margin-bottom: 10px;
        font-size: 16px;
        font-weight: bold;
      ",
        tags$span(style = "color: #2166ac;", "Control higher"),
        tags$span(style = "color: #b2182b;", "Treatment higher")
      ),
      br(),
      
      fluidRow(
        column(3, numericInput("cccc_cluster_vertex_label_size", "Vertex label size", value = 1)),
        column(3, numericInput("cccc_cluster_edge_width_max", "Max edge width", value = 8)),
        column(3, numericInput("cccc_cluster_arrow_size", "Arrow size", value = 0.5)),
        column(3, numericInput("cccc_cluster_title_size", "Title size", value = 12))
      ),
      
      fluidRow(
        column(3, numericInput("cccc_cluster_width", "PDF width (inch)", value = 12)),
        column(3, numericInput("cccc_cluster_height", "PDF height (inch)", value = 6))
      ),
      
      fluidRow(
        column(
          12,
          textInput(
            "cccc_cluster_save_path",
            "Output directory",
            value = file.path(rv$save_path, "CellCellCommunication"),
            width = "100%"
          )
        )
      ),
      
      fluidRow(
        column(
          6,
          actionButton(
            "save_cccc_cluster_compare",
            "Generate Graph to Save Folder",
            style = "background-color:#87CEFA;color:white;width:100%;"
          )
        ),
        column(
          6,
          downloadButton(
            "download_cccc_cluster_compare",
            "Download Plot to Local",
            style = "background-color:#4CAF50;color:white;width:100%;"
          )
        )
      )
    )
  })
   
   output$netVisual_relative_heatmap_ui <- renderUI({
     
     if (is.null(rv$cccd_obj$data_merge)) {
       return(
         div(
           style = "padding: 15px; background-color: #fff3cd; border: 1px solid #ffeeba; border-radius: 5px;",
           strong("Merged CellChat object required"),
           br(),
           "Relative communication heatmap is only available for merged CellChat objects."
         )
       )
     }
     
     tagList(
       
       h4("Relative Communication Heatmap"),
       
       div(
         style = "display: flex; justify-content: center;",
         uiOutput("netVisual_relative_heatmap_control")
       ),
       div(
         style = "
        display: flex;
        justify-content: center;
        align-items: center;
        gap: 30px;
        margin-top: 8px;
        margin-bottom: 10px;
        font-size: 16px;
        font-weight: bold;
      ",
         tags$span(style = "color: #2166ac;", "Control higher"),
         tags$span(style = "color: #b2182b;", "Treatment higher")
       ),
       br(),
       
       fluidRow(

         column(
           3,
           numericInput(
             "netVisual_relative_heatmap_axis_title_size",
             "Axis title size",
             value = 10,
             min = 1
           )
         ),
         column(
           3,
           numericInput(
             "netVisual_relative_heatmap_axis_text_size",
             "Axis text size",
             value = 8,
             min = 1
           )
         ),
         column(
           3,
           numericInput(
             "netVisual_relative_heatmap_legend_size",
             "Legend size",
             value = 10,
             min = 1
           )
         )
       ),
       
       fluidRow(
         column(
           3,
           numericInput(
             "netVisual_relative_heatmap_width",
             "PDF width (inch)",
             value = 12,
             min = 1
           )
         ),
         column(
           3,
           numericInput(
             "netVisual_relative_heatmap_height",
             "PDF height (inch)",
             value = 6,
             min = 1
           )
         )
       ),
       
       fluidRow(
         column(
           12,
           textInput(
             "netVisual_relative_heatmap_save_path",
             "Output directory",
             value = file.path(rv$save_path, "CellCellCommunication"),
             width = "100%"
           )
         )
       ),
       
       fluidRow(
         column(
           6,
           actionButton(
             "save_netVisual_relative_heatmap",
             "Generate Graph to Save Folder",
             style = "background-color:#87CEFA;color:white;width:100%;"
           )
         ),
         column(
           6,
           downloadButton(
             "download_netVisual_relative_heatmap",
             "Download Plot to Local",
             style = "background-color:#4CAF50;color:white;width:100%;"
           )
         )
       )
     )
   })
  # ----------------- Cell Cell Communication Comparison UI - pathway-----------------
   output$cccc_pathway_control_ui <- renderUI({
     tabsetPanel(
       tabPanel(
         "Differential Interaction Network",
         tagList(actionButton("ccc_pathway_select_btn", "Select Cell Cell Communication pathway to display", class = "btn-success", width = "100%"),
         uiOutput("diffInteraction_pathway_ui"))
       ),
       tabPanel(
         "Signaling Pathway Information Flow Ranking",
         uiOutput("rankNet_ui")
       )
     )
   })
   output$diffInteraction_pathway_ui <- renderUI({
     
     tagList(
       
       div(
         style = "display: flex; justify-content: center;",
         uiOutput("diffInteraction_pathway_control")
       ),
       div(
         style = "
        display: flex;
        justify-content: center;
        align-items: center;
        gap: 30px;
        margin-top: 8px;
        margin-bottom: 10px;
        font-size: 16px;
        font-weight: bold;
      ",
         tags$span(style = "color: #2166ac;", "Control higher"),
         tags$span(style = "color: #b2182b;", "Treatment higher")
       ),
       br(),
       
       fluidRow(
         column(3, numericInput("diffInteraction_pathway_vertex_label_size", "Vertex label size", value = 1)),
         column(3, numericInput("diffInteraction_pathway_edge_width_max", "Max edge width", value = 3)),
         column(3, numericInput("diffInteraction_pathway_arrow_size", "Arrow size", value = 0.5)),
         column(3, numericInput("diffInteraction_pathway_title_size", "Title size", value = 12))
       ),
       
       fluidRow(
         column(3, numericInput("diffInteraction_pathway_width", "PDF width (inch)", value = 12)),
         column(3, numericInput("diffInteraction_pathway_height", "PDF height (inch)", value = 6))
       ),
       
       fluidRow(
         column(
           12,
           textInput(
             "diffInteraction_pathway_save_path",
             "Output directory",
             value = file.path(rv$save_path, "CellCellCommunication"),
             width = "100%"
           )
         )
       ),
       
       fluidRow(
         column(
           6,
           actionButton(
             "save_diffInteraction_pathway",
             "Generate Graph to Save Folder",
             style = "background-color:#87CEFA;color:white;width:100%;"
           )
         ),
         column(
           6,
           downloadButton(
             "download_diffInteraction_pathway",
             "Download Plot to Local",
             style = "background-color:#4CAF50;color:white;width:100%;"
           )
         )
       )
     )
   })
    output$rankNet_ui <- renderUI({
      pathway <- rv$ccc_pathway_select
    
    tagList(
      
      div(
        style = "display: flex; justify-content: center;",
        uiOutput("rankNet_control")
      ),
      
      br(),
      
      fluidRow(
        column(3, numericInput("rankNet_font_size", "Font size", value = 10, min = 1)),
        column(3, numericInput("rankNet_title_size", "Title size", value = 12, min = 1)),
        column(3, numericInput("rankNet_width", "PDF width (inch)", value = 10, min = 1)),
        column(3, numericInput("rankNet_height", "PDF height (inch)", value = 8, min = 1))
      ),
      
      fluidRow(
        column(
          12,
          textInput(
            "rankNet_save_path",
            "Output directory",
            value = file.path(rv$save_path, "CellCellCommunication"),
            width = "100%"
          )
        )
      ),
      
      fluidRow(
        column(
          6,
          actionButton(
            "save_rankNet",
            "Generate Graph to Save Folder",
            style = "background-color:#87CEFA;color:white;width:100%;"
          )
        ),
        column(
          6,
          downloadButton(
            "download_rankNet",
            "Download Plot to Local",
            style = "background-color:#4CAF50;color:white;width:100%;"
          )
        )
      )
    )
  })

  # ----------------- Data Output UI-----------------
  output$data_output_ui <- renderUI({
    req(input$do_tabs)
    tagList(
      switch(input$do_tabs,
             "gene_expression" = uiOutput("gene_expression_ui")
      ),
      br(),br(),br(),br(),br()
    )
  })
  # ----------------- Data Output - Gene Expression UI-----------------
  output$gene_expression_ui <- renderUI({
    tagList(
      h4("Seurat RDS Obj:"),
      path_input_ui(
        id = "data",
        label = "Seurat File Path",
        placeholder = "/path/to/seurat_file.rds",
        value = rv$data_obj_path),
      h4("DEG result (Optional):"),
      path_input_ui(
        id = "deg",
        label = "DEG File Path (.tsv)",
        placeholder = "/path/to/deg_result.tsv",
        value = rv$deg_obj_path),
      uiOutput("gene_expression_sub_ui")
    )
  })
  
  # ----------------- DEG FUNCTION -----------------
  
  # Reset handlers
  observeEvent(input$reset_data_obj_btn, {
    rv$data_obj <- NULL
    rv$data_obj_path <- NULL
  })
  
  observeEvent(input$reset_deg_obj_btn, {
    rv$deg_obj <- NULL
    rv$deg_obj_path <- NULL
  })
  
  observeEvent(input$reset_gsea_obj_btn, {
    rv$gsea_obj <- NULL
    rv$gsea_obj_path <- NULL
  })
  
  observeEvent(input$reset_cccd_obj_btn, {
    rv$cccd_obj <- NULL
    rv$cccd_obj_path <- NULL
  })
  
  # DEG window 1 - Analysis Summary
  observeEvent(input$analysis_summary, {
    sample_cards_data <- step2_data$sample_cards_data()
    cluster_cards_data <- step2_data$cluster_cards_data()
    make_table_rows <- function(cards_data, all_comps_func, type_label) {
      rows <- lapply(seq_along(cards_data), function(i) {
        card <- cards_data[[i]]
        if (is.null(card)) return(NULL)
        
        meta_col <- paste(card$selected_meta, collapse = ", ")
        meta_col <- gsub("([^,]+)", "<span style='color:darkred;'>\\1</span>", meta_col)
        
        comp_names <- sapply(card$selected_comps, function(comp_id) {
          comp_name <- sapply(all_comps_func(), function(x) if (x$id==comp_id) x$name else NULL)
          comp_name <- comp_name[!sapply(comp_name, is.null)]
          comp_name <- sapply(comp_name, function(name) {
            parts <- strsplit(name, ":")[[1]]
            prefix <- parts[1]        
            rest <- trimws(parts[2])  
            vs_parts <- strsplit(rest, " vs ")[[1]]
            if(length(vs_parts) == 2){
              paste0(
                prefix, ": ",
                "<span style='color:darkred;'>", vs_parts[1], "</span> vs <span style='color:darkred;'>", vs_parts[2], "</span>"
              )
            } else name
          })
          paste(comp_name, collapse = "<br>")  
        })
        comp_names <- paste(comp_names, collapse = "<br>")
        
        tags$tr(
          tags$td(paste0(type_label, " Card ", i), style="word-break: break-word;"),
          tags$td(HTML(meta_col), style="word-break: break-word;"),
          tags$td(HTML(comp_names), style="word-break: break-word;")
        )
      })
      do.call(tagList, rows)
    }
    showModal(
      modalDialog(
        title = "Analysis Summary",
        size = "l",
        easyClose = FALSE,
        tags$div(style="overflow-y:auto;",
                 tags$table(
                   class = "table table-bordered table-striped",
                   tags$thead(
                     tags$tr(
                       tags$th(style="width:33%", "Card"),
                       tags$th(style="width:33%", "Meta Column"),
                       tags$th(style="width:33%", "Selected Comparisons")
                     )
                   ),
                   tags$tbody(
                     tagList(
                       make_table_rows(step2_data$sample_cards_data(), step2_data$all_sample_comps, "Sample"),
                       make_table_rows(step2_data$cluster_cards_data(), step2_data$all_cluster_comps, "Cluster")
                     )
                   )
                 )
        ),
        br(),
        
        footer = tagList(
          modalButton("Close"),
          actionButton("confirm_summary", "Next Step", class = "btn-success")
        )
      )
    )
  })
  
  # DEG window 2 - Save Directory
  observeEvent(input$confirm_summary, {
    removeModal()  
    
    showModal(
      modalDialog(
        title = "Save Directory and Start Analysis",
        size = "l",
        easyClose = FALSE,
        
        textInput("save_path", "Save Directory:", value = "", width = "100%"),
        uiOutput("path_status"),
        br(),
        
        footer = tagList(
          modalButton("Close"),
          actionButton("start_deg", "Start DEG Analysis", class = "btn-success")
        )
      )
    )
    
    observe({
      path <- input$save_path
      
      if (is.null(path) || path == "") {
        status <- span("Please enter a path.", style="color:red;")
      } else if (!dir.exists(path)) {
        status <- span("Path does not exist, it will be created.", style="color:red;")
      } else {
        files <- list.files(path)
        if (length(files) == 0) {
          status <- span("Path exists and is empty.", style="color:green;")
        } else {
          status <- span("Path exists and contains files. Files may be overwritten.", style="color:orange;")
        }
      }
      
      output$path_status <- renderUI({ status })
    })
  })
  
  # DEG analysis execution
  observeEvent(input$start_deg, {
    source("DE_analysis.R")
    save_path <- input$save_path
    if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
    sample_cards_data <- step2_data$sample_cards_data()
    cluster_cards_data <- step2_data$cluster_cards_data()
    
    all_sample_comps <- step2_data$all_sample_comps()
    all_cluster_comps <- step2_data$all_cluster_comps()
    
    result_list <- list()
    ix <- 1
    # --- Sample Cards ---
    for (card in sample_cards_data) {
      if (is.null(card)) next
      metas <- card$selected_meta
      meta_str <- paste(metas, collapse = ",")
      comps <- card$selected_comps
      for (meta in metas) {
        for (comp_id in comps) {
          comp <- all_sample_comps[[which(sapply(all_sample_comps, function(x) x$id == comp_id))]]
          result_list[[ix]] <- data.frame(
            SAMPLE_ID_1 = comp$numerator,
            SAMPLE_ID_2 = comp$denominator,
            CLUSTER_ID_1 = meta_str,
            CLUSTER_ID_2 = meta_str,
            ix = ix,
            stringsAsFactors = FALSE
          )
          ix <- ix + 1
        }
      }
    }
    for (card in cluster_cards_data) {
      if (is.null(card)) next
      metas <- card$selected_meta
      meta_str <- paste(metas, collapse = ",")
      comps <- card$selected_comps
      for (meta in metas) {
        for (comp_id in comps) {
          comp <- all_cluster_comps[[which(sapply(all_cluster_comps, function(x) x$id == comp_id))]]
          result_list[[ix]] <- data.frame(
            SAMPLE_ID_1 = meta_str,
            SAMPLE_ID_2 = meta_str,
            CLUSTER_ID_1 = comp$numerator,
            CLUSTER_ID_2 = comp$denominator,
            ix = ix,
            stringsAsFactors = FALSE
          )
          ix <- ix + 1
        }
      }
    }
    de_config_validated <- do.call(rbind, result_list)
    file_path <- file.path(save_path, "DEG_config.csv")
    write.table(de_config_validated, file = file_path, sep = "\t", row.names = FALSE, quote = FALSE)
    sample_column <- input$sample_column
    cluster_column <- input$cluster_column
    seurat_obj <- rv$data_obj
    n <- nrow(de_config_validated)
    showModal(modalDialog(
      title = "DE Analysis in Progress",
      "Please do not close the Shiny app until the process completes.",
      footer = NULL,
      easyClose = FALSE
    ))
    
    withProgress(message = "Running DE analysis...", value = 0, {
      for (i in seq_len(n)) {
        row <- de_config_validated[i, ]
        res_list <- compute_de(
          seurat_obj,
          row$SAMPLE_ID_1,
          row$CLUSTER_ID_1,
          row$SAMPLE_ID_2,
          row$CLUSTER_ID_2,
          row$ix,
          sample_column,
          cluster_column
        )
        
        write_tsv(res_list$segex_output, file.path(save_path, res_list$segex_filename), col_names = TRUE)
        incProgress(1/n, detail = paste("Processing comparison", i, "of", n))
      }
    })
    
    showModal(modalDialog(
      title = "DEG Analysis Completed",
      tagList(
        tags$p("✅ DEG analysis has been completed successfully!"),
        tags$p("You can check the output files here:"),
        tags$pre(style = "color: red; user-select: text;", input$save_path)
      ),
      easyClose = TRUE,            
      size = "m"
    ))
  })
  
  # ----------------- GSEA FUNCTION -----------------
  
  # GSEA window 1 - Preview
  observeEvent(input$preview_gsea, {
    
    file_path <- input$gsea_file
    status_msg <- NULL
    files <- NULL
    
    if (is.null(file_path) || file_path == "") {
      status_msg <- HTML("<span style='color:red; font-weight:bold;'>Error: Please enter a file path.</span>")
    } else if (!file.exists(file_path)) {
      status_msg <- HTML("<span style='color:red; font-weight:bold;'>Error: File or folder does not exist.</span>")
    } else if (dir.exists(file_path)) {
      files <- list.files(file_path, pattern = "\\.tsv$", full.names = TRUE)
      if (length(files) == 0) {
        status_msg <- HTML("<span style='color:red; font-weight:bold;'>Error: No .tsv files found in the folder.</span>")
      }
    } else {
      if (!grepl("\\.tsv$", file_path, ignore.case = TRUE)) {
        status_msg <- HTML("<span style='color:red; font-weight:bold;'>Error: File is not a .tsv file.</span>")
      } else {
        files <- file_path
      }
    }
    
    showModal(modalDialog(
      title = "Confirm GSEA Run",
      size = "l",
      if (!is.null(status_msg)) {
        status_msg
      } else {
        tagList(
          tags$b("Files to be processed:"),
          tags$ul(lapply(files, tags$li)),
          
          tags$hr(),
          
          tags$b("Parameters:"),
          tags$ul(
            tags$li(paste("Species:", input$species)),
            tags$li(paste("Method:", input$gsea_method)),
            tags$li(paste("DB:", paste(input$gsea_db, collapse = ", ")))
          )
        )
      },
      
      footer = tagList(
        modalButton("Close"),
        if (is.null(status_msg)) actionButton("confirm_gsea", "Confirm", class = "btn-success")
      )
    ))
  })
  
  # GSEA window 2 - Output path
  observeEvent(input$confirm_gsea, {
    
    showModal(modalDialog(
      title = "Select Output Path",
      size = "l",
      textInput("output_path_gsea", "Output directory:", "", width = "100%"),
      
      uiOutput("path_status_gsea"),
      
      footer = tagList(
        modalButton("Close"),
        actionButton("start_gsea", "Start GSEA analysis", class = "btn-success")
      )
    ))
  })
  
  output$path_status_gsea <- renderUI({
    
    path <- input$output_path_gsea
    
    if (is.null(path) || path == "") {
      status <- span("Please enter a path.", style="color:red;")
      
    } else if (!dir.exists(path)) {
      status <- span("Path does not exist, it will be created.", style="color:red;")
      
    } else {
      files <- list.files(path)
      
      if (length(files) == 0) {
        status <- span("Path exists and is empty.", style="color:green;")
      } else {
        status <- span("Path exists and contains files. Files may be overwritten.", style="color:orange;")
      }
    }
    
  })
  
  # GSEA analysis execution
  observeEvent(input$start_gsea, {
    
    req(input$output_path_gsea)
    req(input$gsea_file)
    
    out_dir <- input$output_path_gsea
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    
    files <- if (dir.exists(input$gsea_file)) {
      list.files(input$gsea_file, pattern = "\\.tsv$", full.names = TRUE)
    } else {
      input$gsea_file
    }
    
    if (length(files) == 0) {
      showNotification("No TSV files to process.", type = "error", duration = 5, closeButton = TRUE)
      return()
    }
    
    source("GSEA_analysis.R")  
    
    showModal(modalDialog(
      title = "GSEA Analysis Running",
      size = "l",
      tags$div(id = "gsea_progress_detail", "Starting..."),
      footer = NULL,
      easyClose = FALSE
    ))
    
    n_files <- length(files)
    
    for (i in seq_along(files)) {
      file <- files[i]
      done <- i - 1
      remaining <- length(files) - i + 1
      shinyjs::html(
        "gsea_progress_detail",
        paste0(
          "Processing file: ", basename(file), "<br>",
          "Files completed: ", done, "<br>",
          "Files remaining: ", remaining
        )
      )
      
      gsea_res <- gsea_analysis_from_tsv(
        file_path = file,
        species = input$species[1],
        method = input$gsea_method,
        db = input$gsea_db
      )
      
      out_file <- file.path(out_dir, paste0("GSEA_", tools::file_path_sans_ext(basename(file)), ".xlsx"))
      wb <- createWorkbook()
      for (db in names(gsea_res)) {
        df <- as.data.frame(gsea_res[[db]])
        addWorksheet(wb, sheetName = db)
        writeData(wb, sheet = db, x = df)
      }
      saveWorkbook(wb, out_file, overwrite = TRUE)
    }
    
    showModal(modalDialog(
      title = "GSEA Analysis Completed",
      tagList(
        tags$p("✅ GSEA analysis has been completed successfully!"),
        tags$p("You can check the output files here:"),
        tags$pre(style = "color: red; user-select: text;", input$output_path_gsea)
      ),
      easyClose = TRUE,             
      size = "m"
    ))
    
    showNotification(
      HTML("<div style='font-size:20px; font-weight:bold;'>✅ GSEA analysis finished!</div>"),
      type = "message", duration = 10, closeButton = TRUE
    )
  })
  
  # ----------------- Cell Cell Communication FUNCTION -----------------
  observeEvent(input$open_sample_group_modal, {
    req(rv$data_obj)
    rv$ui_freeze = TRUE
    showModal(
      modalDialog(
        title = "Set sample groups",
        size = "l",
        easyClose = FALSE,
        
        selectInput(
          "sample_column_modal",
          HTML("Select meta.data <b style='color: red;'>sample</b> column"),
          choices = c("All Together", extract_meta_columns(rv$data_obj)),
          selected = if (!is.null(rv$sample_infor$sample_column)) rv$sample_infor$sample_column else NULL,
          width = "100%"
        ),
        
        uiOutput("sample_group_bucket_ui"),
        
        footer = tagList(
          actionButton("confirm_sample_groups", "Confirm", class = "btn-primary")
        )
      )
    )
  })
  
  output$sample_group_bucket_ui <- renderUI({
    req(rv$data_obj)
    req(input$sample_column_modal)
    
    if (input$sample_column_modal == "All Together") {
      sample_levels <- "All Together"
    } else {
      sample_levels <- unique(as.character(
        rv$data_obj@meta.data[[input$sample_column_modal]]
      ))
      sample_levels <- sample_levels[!is.na(sample_levels)]
      sample_levels <- sort(sample_levels)
    
    tagList(
      
      rank_list(
        text = "Initial Pool",
        labels = sample_levels,
        input_id = "initial_pool",
        options = sortable_options(group = "sample_groups")
      ),
      
      br(),
      
      fluidRow(
        
        column(
          6,
          rank_list(
            text = "CONTROL",
            labels = character(0),
            input_id = "control_pool",
            options = sortable_options(group = "sample_groups")
          )
        ),
        
        column(
          6,
          rank_list(
            text = "TREATMENT",
            labels = character(0),
            input_id = "treatment_pool",
            options = sortable_options(group = "sample_groups")
          )
        )
        
      )
    )}
  })
  
  observeEvent(input$confirm_sample_groups, {
    
    req(input$sample_column_modal)
    
    rv$ui_freeze <- FALSE
    
    if (input$sample_column_modal == "All Together") {
      
      rv$sample_infor <- list(
        sample_column = "All Together",
        control = NULL,
        treatment = NULL,
        unused = NULL
      )
      
    } else {
      
      validate(
        need(length(input$control_pool) > 0,
             "Please assign at least one sample to CONTROL"),
        need(length(input$treatment_pool) > 0,
             "Please assign at least one sample to TREATMENT")
      )
      
      rv$sample_infor <- list(
        sample_column = input$sample_column_modal,
        control = input$control_pool,
        treatment = input$treatment_pool,
        unused = input$initial_pool
      )
      
    }
    
    removeModal()
  })
  
  observeEvent(input$recommend_cc_prob_type, {
    if (
      is.null(input$cluster_column) ||
      input$cluster_column == "" ||
      input$cluster_column == "None"
    ) {
      
      showModal(
        modalDialog(
          title = "Missing Cell Type Annotation",
          
          tags$div(
            style = "text-align:center;",
            
            tags$h4(
              "Please add meta.data celltype cluster identification first.",
              style = "color:red;"
            ),
            
            tags$p(
              "Select `meta.data celltype cluster identification` before calculating the recommended CellChat probability method."
            )
          ),
          
          easyClose = TRUE,
          footer = modalButton("Close")
        )
      )
      
      return(NULL)
    }
    req(rv$data_obj)
    
    showModal(
      modalDialog(
        title = "Calculating recommended method",
        div(
          style = "text-align:center;",
          tags$h4("Calculating recommended CellChat probability method..."),
          tags$p("Please wait.")
        ),
        footer = NULL,
        easyClose = FALSE
      )
    )
    source("CCC.R")
    res <- diagnose_cellchat_summary_method(
      seurat_obj = rv$data_obj,
      group.by = input$cluster_column,
      assay = "RNA",
      layer = "data",
      cellchat_db = input$cc_species,
      trim_values = c(0.1, 0.05)
    )
    
    recommended_method <- res$recommendation$recommended_method
    
    recommended_method_display <- switch(
      recommended_method,
      "triMean" = "triMean",
      "truncatedMean(0.1)" = "truncatedmean_0.1",
      "truncatedMean(0.05)" = "truncatedmean_0.05",
      "triMean"
    )
    
    updateRadioButtons(
      session,
      "cc_prob_type",
      selected = recommended_method_display
    )
    
    showModal(
      modalDialog(
        title = "Recommended CellChat Probability Method",
        
        tags$div(
          tags$h4("Recommended method:"),
          tags$h3(
            recommended_method,
            style = "color:red; font-weight:bold;"
          ),
          tags$p(
            res$recommendation$final_recommendation,
            style = "font-size:14px;"
          )
        ),
        
        tags$hr(),
        
        tags$h4("Method summary"),
        tags$p("Using the list of genes from all LR (ligand-receptor) pairs, identify interactions between these genes and the gene list in our dataset; perform calculations using all four methods and count the number of genes with values greater than zero."),
        tags$div(
          style = "overflow-x:auto; width:100%;",
          tableOutput("cc_method_summary_table")
        ),
        
        easyClose = TRUE,
        footer = modalButton("Close"),
        size = "l"
      )
    )
    
    output$cc_method_summary_table <- renderTable({
      res$method_summary
    })
  })
  
  observeEvent(input$cc_run_analysis, {
    
    req(rv$data_obj)
    req(input$cluster_column)
    showModal(
      modalDialog(
        title = "Review CellChat Parameters",
        size = "l",
        
        h4("📋 Review Parameters"),
        
        tags$div(
          style = "background-color: #f5f5f5; padding: 15px; border-radius: 5px; border-left: 4px solid #2196F3;",
          
          tags$table(
            style = "width: 100%;",
            tags$tr(
              tags$td(style = "padding: 8px; font-weight: bold; width: 150px;", "Data File:"),
              tags$td(style = "padding: 8px;", tags$code(input$data_path))
            ),
            tags$tr(
              tags$td(style = "padding: 8px; font-weight: bold;", "Sample Setting:"),
              tags$td(
                style = "padding: 8px;",
                if (is.null(rv$sample_infor)) {
                  "Not set"
                } else if (rv$sample_infor$sample_column == "All Together") {
                  "All Together"
                } else {
                  paste0(
                    rv$sample_infor$sample_column,
                    " | CONTROL: ", paste(rv$sample_infor$control, collapse = ", "),
                    " | TREATMENT: ", paste(rv$sample_infor$treatment, collapse = ", ")
                  )
                }
              )
            ),
            tags$tr(
              tags$td(style = "padding: 8px; font-weight: bold;", "Cell Type Column:"),
              tags$td(style = "padding: 8px;", input$cluster_column)
            ),
            tags$tr(
              tags$td(style = "padding: 8px; font-weight: bold;", "Species:"),
              tags$td(style = "padding: 8px;", input$cc_species)
            ),
            tags$tr(
              tags$td(style = "padding: 8px; font-weight: bold;", "Pathways:"),
              tags$td(style = "padding: 8px;", paste(input$cc_pathway_type, collapse = ", "))
            ),
            tags$tr(
              tags$td(style = "padding: 8px; font-weight: bold;", "Probability Method:"),
              tags$td(style = "padding: 8px;", input$cc_prob_type)
            ),
            tags$tr(
              tags$td(style = "padding: 8px; font-weight: bold;", "Min Cells:"),
              tags$td(style = "padding: 8px;", input$cc_min_cells)
            )
          )
        ),
        
        hr(),
        
        h4("📝 Analysis Configuration"),
        
        fluidRow(
          column(6, textInput("cc_project_name", "Project Name on SCC", 
                              placeholder = "e.g., wax-es",
                              width = "100%")),
          column(6, numericInput("cc_core_number", "CPU Cores", 
                                 value = 8, min = 1, max = 32, width = "100%"))
        ),
        
        fluidRow(
          column(6, numericInput("cc_task_hours", "Task Time (hours)", 
                                 value = 12, min = 0.5, max = 24, step = 0.5, width = "100%")),
          
          column(6, textInput("cc_email", "Email (for notification)", 
                               placeholder = "user@example.com",
                               width = "100%"))

        ),
        textInput("cc_save_path", "Result Folder", 
                            placeholder = "/path/to/save/folder",
                            width = "100%"),
        
        
        footer = tagList(
          modalButton("Close"),
          actionButton("cc_start_analysis", "Start Analysis on Background", 
                       class = "btn-success btn-lg",
                       style = "margin-left: 10px;")
        )
      )
    )
  })
  
  observeEvent(input$cc_start_analysis,{
    source("CCC.R")
    
    if (is.null(input$cc_project_name) || input$cc_project_name == "") {
      showNotification("Please enter project name on SCC", type = "error")
      return()
    }
    
    if (is.null(input$cc_save_path) || input$cc_save_path == "") {
      showNotification("Please enter result folder path", type = "error")
      return()
    }
    req(rv$sample_infor)
    
    sample_info <- rv$sample_infor
    
    if (sample_info$sample_column != "All Together") {
      if (length(sample_info$control) == 0 || length(sample_info$treatment) == 0) {
        showNotification("Please assign samples to both CONTROL and TREATMENT.", type = "error")
        return()
      }
    }
    cluster_column <- input$cluster_column
    cc_species <- input$cc_species
    cc_pathway_type <- input$cc_pathway_type
    cc_prob_type <- input$cc_prob_type
    cc_min_cells <- input$cc_min_cells
    cc_save_path <- input$cc_save_path
    cc_project_name <- input$cc_project_name
    cc_task_hours <- input$cc_task_hours
    cc_core_number <- input$cc_core_number
    cc_email <- input$cc_email
    
    showModal(modalDialog(
      title = "Preparing Cell Cell Communication mission file",
      "Please wait, after data preparing completed you will get the qsub mission id.",
      easyClose = FALSE,
      footer = NULL
    ))
    
    tryCatch({
      mission_id <- run_ccc(
        seurat_file_path = input$data_path,
        sample_info = sample_info,
        cluster_column = cluster_column,
        species = cc_species,
        pathway_type = cc_pathway_type,
        prob_type = cc_prob_type,
        min_cells = cc_min_cells,
        save_path = cc_save_path,
        project_name = cc_project_name,
        runtime = cc_task_hours,
        cores = cc_core_number,
        email = cc_email
      )
      
      rv$cc_mission_id <- mission_id
      
      showModal(
        modalDialog(
          title = "Clustering Visualization Job Submitted",
          tagList(
            p("The Clustering analysis is now running in the background."),
            p("Your job ID is:"),
            tags$pre(style = "color: red; user-select: text;", mission_id),
            p("You can copy the following command in the SCC terminal to check the job status:"),
            tags$pre(style = "color: red; user-select: text;", paste0("qstat -j ", mission_id)),
            p("You can now close the shinyapp and wait for the results."),
            p("For subsequent analysis, it is recommended to use the following path of seurat rds file as input:"),
            tags$pre(style = "color: blue; user-select: text;", 
                     file.path(input$cc_save_path, "input_seurat_file.rds")),
            if (nchar(input$cc_email) > 0) {
              p("The results will be sent to your email address: ", input$cc_email)
            } else {
              NULL
            }
          ),
          footer = modalButton("Close"),
          easyClose = TRUE
        )
      )
      
    }, error = function(e) {
      removeModal()
      showModal(modalDialog(title = "Error", e$message, easyClose = TRUE))
    })
  })
  # ----------------- VISUALIZATION FUNCTION Marker Select -----------------
  
  observeEvent(input$marker_select_btn, {
    rv$ui_freeze = TRUE
    showModal(
      modalDialog(
        title = "Selection of genes for analysis",
        size = "l",
        easyClose = FALSE,
        
        fluidRow(
          column(
            width = 6,
            textAreaInput(
              inputId = "marker_text",
              label = "Gene List:",
              value = rv$visualization_marker,
              rows = 20,
              placeholder = "Paste genes here:\nGene_A\nGene_B\nGene_C",
              width = "100%"
            ),
            actionButton(
              "clear_marker",
              "Clear All Genes",
              class = "btn-danger",
              width = "100%"
            )
          ),
          column(
            width = 6,
            tagList(
              # --- Marker count ---
              tags$div(
                style = "margin-bottom: 10px; display: flex; align-items: center; gap: 6px;",
                tags$span("Current Number of Genes:"),
                tags$span(
                  style = "color: red; font-weight: bold;",
                  textOutput("marker_count")
                )
              ),
              # =========================
              # CASE 1: DEG NOT AVAILABLE
              # =========================
              conditionalPanel(
                condition = "output.deg_available == false",
                
                tags$div(
                  style = "padding: 15px; color: #777;",
                  tags$h5("DEG not loaded"),
                  tags$p("Upload DEG file to enable advanced marker selection.")
                )
              ),
              
              # =========================
              # CASE 2: DEG AVAILABLE
              # =========================
              conditionalPanel(
                condition = "output.deg_available == true",
                
                # --- DEG button ---
                actionButton(
                  "use_deg_marker_btn",
                  "Use DEG file to search for genes",
                  class = "btn btn-primary",
                  width = "100%"
                ),
                
                tags$hr(),
                
                # =========================
                # 1. Top N markers
                # =========================
                fluidRow(
                  column(
                    width = 8,
                    numericInput(
                      "top_n_pvalue",
                      "Top N (smallest p-value):",
                      value = 50,
                      min = 1
                    )
                  ),
                  column(
                    width = 4,
                    br(),
                    actionButton(
                      "apply_top_n",
                      "Apply",
                      class = "btn btn-success",
                      width = "100%"
                    )
                  )
                ),
                
                # =========================
                # 2. Top % markers
                # =========================
                fluidRow(
                  column(
                    width = 8,
                    numericInput(
                      "top_pct_pvalue",
                      "Top % (smallest p-value):",
                      value = 10,
                      min = 1,
                      max = 100
                    )
                  ),
                  column(
                    width = 4,
                    br(),
                    actionButton(
                      "apply_top_pct",
                      "Apply",
                      class = "btn btn-success",
                      width = "100%"
                    )
                  )
                ),
                
                # =========================
                # 3. pvalue threshold
                # =========================
                fluidRow(
                  column(
                    width = 8,
                    numericInput(
                      "pvalue_threshold",
                      "p-value < threshold:",
                      value = 0.05,
                      min = 0,
                      step = 0.001
                    )
                  ),
                  column(
                    width = 4,
                    br(),
                    actionButton(
                      "apply_pvalue",
                      "Apply",
                      class = "btn btn-success",
                      width = "100%"
                    ))),
                # =========================
                # 4. log2fc
                # =========================
                fluidRow(
                  column(
                    width = 8,
                    numericInput(
                      "log2fc_updown_number",
                      "log2FC upregulated and downregulated number:",
                      value = 10,
                      min = 1,
                      step = 1
                    )
                  ),
                  column(
                    width = 4,
                    br(),
                    actionButton(
                      "apply_log2FC",
                      "Apply",
                      class = "btn btn-success",
                      width = "100%"
                    )))
              )))),
        footer = actionButton("close_modal", "Close")
      ))
  })
  
  observeEvent(input$marker_text, {
    
    lines <- unlist(strsplit(input$marker_text, "\n"))
    
    lines <- lines[trimws(lines) != ""]
    
    rv$visualization_marker <- paste(lines, collapse = "\n")
    
  }, ignoreInit = TRUE)
  
  observeEvent(input$clear_marker, {
    rv$visualization_marker <- ""
    
    updateTextAreaInput(
      session,
      "marker_text",
      value = ""
    )
  })
  
  output$deg_available <- reactive({
    !is.null(rv$deg_obj)
  })
  outputOptions(output, "deg_available", suspendWhenHidden = FALSE)
  
  output$marker_count <- renderText({
    if (is.null(rv$visualization_marker) || rv$visualization_marker == "") {
      return(0)
    }
    length(strsplit(rv$visualization_marker, "\n")[[1]])
  })
  # p value select method
  # top p vlaue select
  observeEvent(input$apply_top_n, {
    req(rv$deg_obj)
    
    df <- rv$deg_obj
    
    df2 <- df[order(df$pvalue_1), ]
    df2 <- head(df2, input$top_n_pvalue)
    
    rv$visualization_marker <- paste(df2$gene, collapse = "\n")
    
    updateTextAreaInput(session, "marker_text", value = rv$visualization_marker)
  })
  # pct select
  observeEvent(input$apply_top_pct, {
    req(rv$deg_obj)
    df <- rv$deg_obj
    df <- df[order(df$pvalue_1), ]
    n <- ceiling(nrow(df) * input$top_pct_pvalue / 100)
    df2 <- head(df, n)
    rv$visualization_marker <- paste(df2$gene, collapse = "\n")
    updateTextAreaInput(session, "marker_text", value = rv$visualization_marker)
  })
  # less than p value selct
  observeEvent(input$apply_pvalue, {
    req(rv$deg_obj)
    
    df <- rv$deg_obj
    
    df2 <- df[df$pvalue_1 < input$pvalue_threshold, ]
    
    rv$visualization_marker <- paste(df2$gene, collapse = "\n")
    
    updateTextAreaInput(session, "marker_text", value = rv$visualization_marker)
  })
  # log2fc upregulated and downregulated
  observeEvent(input$apply_log2FC, {
    req(rv$deg_obj)
    req(input$log2fc_updown_number)
    df <- rv$deg_obj
    n <- input$log2fc_updown_number
    ord <- order(df$log2FC)
    df_down <- df[ord[1:n], ]
    df_up   <- df[rev(ord)[1:n], ]
    df2 <- unique(rbind(df_up, df_down))
    rv$visualization_marker <- paste(df2$gene, collapse = "\n")
    updateTextAreaInput(session, "marker_text", value = rv$visualization_marker)
  })
  # select marker on deg file
  observeEvent(input$use_deg_marker_btn, {
    library(DT)
    req(rv$deg_obj)
    
    showModal(
      modalDialog(
        title = "Select Genes from DEG",
        size = "l",
        easyClose = FALSE,
        
        DT::DTOutput("deg_table"),
        
        footer = tagList(
          actionButton(
            "confirm_deg_marker",
            "Confirm Selection",
            class = "btn btn-primary"
          ),
          actionButton("close_modal", "Close")
        )
      )
    )
  })
  
  output$deg_table <- DT::renderDT({
    
    req(rv$deg_obj)
    
    df <- rv$deg_obj
    
    # ---- Get current marker genes safely ----
    marker_genes <- character(0)
    
    if (!is.null(rv$visualization_marker) && rv$visualization_marker != "") {
      marker_genes <- unlist(strsplit(rv$visualization_marker, "\n"))
      marker_genes <- trimws(marker_genes)
    }
    
    # ---- Only keep genes that exist in DEG table ----
    valid_idx <- which(df$gene %in% marker_genes)
    
    DT::datatable(
      df,
      class = "compact",
      # ---- Preselect valid rows ----
      selection = list(
        mode = "multiple",
        selected = valid_idx
      ),
      
      filter = "top",
      
      options = list(
        pageLength = 15,
        
        lengthMenu = list(
          c(10, 15, 20, 30, 50),   
          c("10", "15", "20", "30", "50")  
        ),
        
        scrollX = TRUE
      )
    )%>%
      DT::formatRound(c("log2FC", "FC (linear)", "intensity_1", "intensity_2", "pvalue_1"), 2)
  })
  
  get_selected_genes <- function() {
    
    req(rv$deg_obj)
    
    sel_idx <- input$deg_table_rows_selected
    
    if (is.null(sel_idx) || length(sel_idx) == 0) {
      return(character(0))
    }
    
    rv$deg_obj$gene[sel_idx]
  }
  
  observeEvent(input$confirm_deg_marker, {
    
    genes <- get_selected_genes()
    
    # ---- Update reactive value ----
    rv$visualization_marker <- paste(genes, collapse = "\n")
    
    # ---- Sync back to textarea ----
    updateTextAreaInput(
      session,
      "marker_text",
      value = rv$visualization_marker
    )
    rv$ui_freeze = TRUE
    showModal(
      modalDialog(
        title = "Marker Selection",
        size = "l",
        easyClose = FALSE,
        
        fluidRow(
          column(
            width = 6,
            textAreaInput(
              inputId = "marker_text",
              label = "Current Marker:",
              value = rv$visualization_marker,
              rows = 20,
              placeholder = "Paste genes here:\nGene_A\nGene_B\nGene_C",
              width = "100%"
            ),
            actionButton(
              "clear_marker",
              "Clear All Markers",
              class = "btn-danger",
              width = "100%"
            )
          ),
          column(
            width = 6,
            tagList(
              # --- Marker count ---
              tags$div(
                style = "margin-bottom: 10px; display: flex; align-items: center; gap: 6px;",
                tags$span("Current Number of Marker:"),
                tags$span(
                  style = "color: red; font-weight: bold;",
                  textOutput("marker_count")
                )
              ),
              # =========================
              # CASE 1: DEG NOT AVAILABLE
              # =========================
              conditionalPanel(
                condition = "output.deg_available == false",
                
                tags$div(
                  style = "padding: 15px; color: #777;",
                  tags$h5("DEG not loaded"),
                  tags$p("Upload DEG file to enable advanced marker selection.")
                )
              ),
              
              # =========================
              # CASE 2: DEG AVAILABLE
              # =========================
              conditionalPanel(
                condition = "output.deg_available == true",
                
                # --- DEG button ---
                actionButton(
                  "use_deg_marker_btn",
                  "Use DEG file to search marker",
                  class = "btn btn-primary",
                  width = "100%"
                ),
                
                tags$hr(),
                
                # =========================
                # 1. Top N markers
                # =========================
                fluidRow(
                  column(
                    width = 8,
                    numericInput(
                      "top_n_pvalue",
                      "Top N (smallest p-value):",
                      value = 50,
                      min = 1
                    )
                  ),
                  column(
                    width = 4,
                    br(),
                    actionButton(
                      "apply_top_n",
                      "Apply",
                      class = "btn btn-success",
                      width = "100%"
                    )
                  )
                ),
                
                # =========================
                # 2. Top % markers
                # =========================
                fluidRow(
                  column(
                    width = 8,
                    numericInput(
                      "top_pct_pvalue",
                      "Top % (smallest p-value):",
                      value = 10,
                      min = 1,
                      max = 100
                    )
                  ),
                  column(
                    width = 4,
                    br(),
                    actionButton(
                      "apply_top_pct",
                      "Apply",
                      class = "btn btn-success",
                      width = "100%"
                    )
                  )
                ),
                
                # =========================
                # 3. pvalue threshold
                # =========================
                fluidRow(
                  column(
                    width = 8,
                    numericInput(
                      "pvalue_threshold",
                      "p-value < threshold:",
                      value = 0.05,
                      min = 0,
                      step = 0.001
                    )
                  ),
                  column(
                    width = 4,
                    br(),
                    actionButton(
                      "apply_pvalue",
                      "Apply",
                      class = "btn btn-success",
                      width = "100%"
                    ))),
                # =========================
                # 4. log2fc
                # =========================
                fluidRow(
                  column(
                    width = 8,
                    numericInput(
                      "log2fc_updown_number",
                      "log2FC upregulated and downregulated number:",
                      value = 10,
                      min = 1,
                      step = 1
                    )
                  ),
                  column(
                    width = 4,
                    br(),
                    actionButton(
                      "apply_log2FC",
                      "Apply",
                      class = "btn btn-success",
                      width = "100%"
                    ))))))),
        footer = actionButton("close_modal", "Close")
      ))
  })
  observeEvent(input$close_modal, {
    
    removeModal()
    rv$ui_freeze <- FALSE
  })
  # ----------------- VISUALIZATION FUNCTION Sample Identification Select -----------------
  observeEvent(input$sample_ident_select_btn, {
    rv$ui_freeze <- TRUE
    if (is.null(rv$data_obj)) {
      
      showModal(modalDialog(
        title = "sample identification select",
        "Please load a Seurat RDS obj to select identification",
        easyClose = FALSE,
        footer = actionButton("close_modal", "Close")
      ))
      
      return()
    }
    
    meta_cols <- colnames(rv$data_obj@meta.data)
    
    meta_df <- rv$data_obj@meta.data
    valid_cols <- meta_cols[
      sapply(meta_df[, meta_cols, drop = FALSE], function(x) {
        !is.numeric(x)
      })
    ]
    
    sample_cols <- c("", setdiff(valid_cols, c("CB", "CB_original")))
    
    if (is.null(rv$visualization_sample_column)) {
      rv$visualization_sample_column <- ""
    }
    
    if (is.null(rv$visualization_split_column)) {
      rv$visualization_split_column <- ""
    }
    
    if (is.null(rv$visualization_use_split)) {
      rv$visualization_use_split <- FALSE
    }
    
    showModal(modalDialog(
      title = "sample/celltype identification select",
      
      tagList(
        
        # ====== select sample column (in graph) ======
        selectInput(
          inputId = "sample_ident_col",
          label   = "Sample/celltype identification column in graph",
          choices = sample_cols,
          selected = rv$visualization_sample_column,
          width = "100%"
        ),
        
        # ====== UMAP plot（conditional UI）======
        uiOutput("umap_in_modal_ui"),
        
        # ====== checkbox for split ======
        checkboxInput(
          inputId = "enable_split_by_metadata",
          label = "Split by additional metadata",
          value = rv$visualization_use_split
        ),
        
        # ====== select split column (conditional) ======
        uiOutput("split_meta_select_ui"),
        
        # ====== UMAP plot for split (conditional) ======
        uiOutput("umap_split_in_modal_ui")
      ),
      
      footer = tagList(
        actionButton("close_modal", "Close"),
        actionButton("confirm_sample_ident", "Confirm", class = "btn-primary")
      ),
      easyClose = FALSE
    ))
  })
  
  output$umap_in_modal_ui <- renderUI({
    
    req(input$sample_ident_col)
    req(input$sample_ident_col != "")
    
    has_umap <- "umap" %in% names(rv$data_obj@reductions)
    
    if (!has_umap) {
      return(tags$div("No UMAP found in this Seurat object"))
    }
    
    tagList(
      hr(),
      h4("Main UMAP Preview"),
      plotOutput("umap_modal_plot", height = "400px")
    )
  })
  
  output$umap_modal_plot <- renderPlot({
    
    req(input$sample_ident_col)
    req(input$sample_ident_col != "")
    
    req("umap" %in% names(rv$data_obj@reductions))
    
    DimPlot(
      rv$data_obj,
      reduction = "umap",
      group.by  = input$sample_ident_col
    )
  })
  
  output$split_meta_select_ui <- renderUI({
    
    req(input$enable_split_by_metadata)
    
    meta_cols <- colnames(rv$data_obj@meta.data)
    
    meta_df <- rv$data_obj@meta.data
    valid_cols <- meta_cols[
      sapply(meta_df[, meta_cols, drop = FALSE], function(x) {
        !is.numeric(x)
      })
    ]
    
    split_cols <- c("", setdiff(valid_cols, c("CB", "CB_original", input$sample_ident_col)))
    
    selectInput(
      inputId = "split_ident_col",
      label   = "Metadata column to split by",
      choices = split_cols,
      selected = if (!is.null(rv$visualization_split_column)) rv$visualization_split_column else "",
      width = "100%"
    )
  })
  
  output$umap_split_in_modal_ui <- renderUI({
    
    req(input$enable_split_by_metadata)
    req(input$split_ident_col)
    req(input$split_ident_col != "")
    
    has_umap <- "umap" %in% names(rv$data_obj@reductions)
    
    if (!has_umap) {
      return(tags$div("No UMAP found in this Seurat object"))
    }
    
    tagList(
      hr(),
      h4("Split UMAP Preview"),
      plotOutput("umap_split_modal_plot", height = "500px")
    )
  })
  
  output$umap_split_modal_plot <- renderPlot({
    
    req(input$enable_split_by_metadata)
    req(input$split_ident_col)
    req(input$split_ident_col != "")
    
    req("umap" %in% names(rv$data_obj@reductions))
    
    DimPlot(
      rv$data_obj,
      reduction = "umap",
      group.by  = input$sample_ident_col,
      split.by  = input$split_ident_col,
      ncol      = 3
    )
  })
  
  observeEvent(input$confirm_sample_ident, {
    
    req(rv$data_obj)
    
    rv$visualization_sample_column <- input$sample_ident_col
    rv$visualization_use_split <- input$enable_split_by_metadata
    
    if (input$enable_split_by_metadata && !is.null(input$split_ident_col)) {
      rv$visualization_split_column <- input$split_ident_col
    } else {
      rv$visualization_split_column <- ""
    }
    
    rv$ui_freeze = FALSE
    removeModal()
  })
  
  # ----------------- VISUALIZATION FUNCTION Pathway Select ----------------- 
  observeEvent(input$pathway_select_btn, {
    rv$ui_freeze = TRUE
    showModal(
      modalDialog(
        title = "Selection of pathways for analysis",
        size = "l",
        easyClose = FALSE,
        
        fluidRow(
          column(
            width = 6,
            textAreaInput(
              inputId = "pathway_text",
              label = "Pathway List:",
              value = rv$visualization_pathway,
              rows = 20,
              placeholder = "Pathway_ID_1\nPathway_ID_2",
              width = "100%"
            ),
            actionButton(
              "clear_pathway",
              "Clear All Pathways",
              class = "btn-danger",
              width = "100%"
            )
          ),
          column(
            width = 6,
            tagList(
              # --- Marker count ---
              tags$div(
                style = "margin-bottom: 10px; display: flex; align-items: center; gap: 6px;",
                tags$span("Current Number of Pathways:"),
                tags$span(
                  style = "color: red; font-weight: bold;",
                  textOutput("pathway_count")
                )
              ),
              # =========================
              # CASE 1: GSEA NOT AVAILABLE
              # =========================
              conditionalPanel(
                condition = "output.gsea_available == false",
                
                tags$div(
                  style = "padding: 15px; color: #777;",
                  tags$h5("GSEA not loaded"),
                  tags$p("Upload GSEA file to enable advanced pathway selection.")
                )
              ),
              
              # =========================
              # CASE 2: GSEA AVAILABLE
              # =========================
              conditionalPanel(
                condition = "output.gsea_available == true",
                
                # --- DEG button ---
                actionButton(
                  "use_gsea_table_btn",
                  "Use GSEA file to search for pathways",
                  class = "btn btn-primary",
                  width = "100%"
                ),
                
                tags$hr(),
                
                # ---- Top N NES ----
                fluidRow(
                  column(8,
                         numericInput("top_n_nes", "Top N |NES|:", 20, min = 1)
                  ),
                  column(4,
                         br(),
                         actionButton("apply_top_n_nes", "Apply", class = "btn-success")
                  )
                ),
                
                # ---- NES threshold ----
                fluidRow(
                  column(8,
                         numericInput("nes_threshold", "|NES| >", 1.5)
                  ),
                  column(4,
                         br(),
                         actionButton("apply_nes", "Apply", class = "btn-success")
                  )
                ),
                
                # ---- pvalue ----
                fluidRow(
                  column(8,
                         numericInput("gsea_pval", "pvalue <", 0.05)
                  ),
                  column(4,
                         br(),
                         actionButton("apply_gsea_pval", "Apply", class = "btn-success")
                  )
                ),
                # =========================
                # 4. log2fc
                # =========================
                # ---- padj ----
                fluidRow(
                  column(8,
                         numericInput("gsea_padj", "p.adjust <", 0.05)
                  ),
                  column(4,
                         br(),
                         actionButton("apply_gsea_padj", "Apply", class = "btn-success")
                  )))))),
        footer = actionButton("close_modal", "Close")
      ))
  })
  
  observeEvent(input$pathway_text, {
    
    lines <- unlist(strsplit(input$pathway_text, "\n"))
    
    lines <- lines[trimws(lines) != ""]
    
    rv$visualization_pathway <- paste(lines, collapse = "\n")
    rv$pathway_rename <- data.frame(
      original_pathway = lines,
      renamed_pathway = lines |>
        gsub("^[^_]+_", "", x = _) |>
        gsub("_", " ", x = _) |>
        str_to_title()
    )
  }, ignoreInit = TRUE)
  
  observeEvent(input$clear_pathway, {
    rv$visualization_pathway <- ""
    
    updateTextAreaInput(
      session,
      "pathway_text",
      value = ""
    )
  })
  
  output$gsea_available <- reactive({
    !is.null(rv$gsea_obj)
  })
  
  outputOptions(output, "gsea_available", suspendWhenHidden = FALSE)
  
  output$pathway_count <- renderText({
    if (is.null(rv$visualization_pathway) || rv$visualization_pathway == "") {
      return(0)
    }
    length(strsplit(rv$visualization_pathway, "\n")[[1]])
  })
  #top n NES
  observeEvent(input$apply_top_n_nes, {
    req(rv$gsea_obj)
    df <- rv$gsea_obj
    
    df <- df[order(abs(df$NES), decreasing = TRUE), ]
    df2 <- head(df, input$top_n_nes)
    
    rv$visualization_pathway <- paste(df2$ID, collapse = "\n")
    
    updateTextAreaInput(session, "pathway_text", value = rv$visualization_pathway)
  })
  
  # abs NES >
  observeEvent(input$apply_nes, {
    req(rv$gsea_obj)
    df <- rv$gsea_obj
    
    df2 <- df[abs(df$NES) > input$nes_threshold, ]
    
    rv$visualization_pathway <- paste(df2$ID, collapse = "\n")
    
    updateTextAreaInput(session, "pathway_text", value = rv$visualization_pathway)
  })
  
  #pvalue <
  observeEvent(input$apply_gsea_pval, {
    req(rv$gsea_obj)
    df <- rv$gsea_obj
    
    df2 <- df[df$pvalue < input$gsea_pval, ]
    
    rv$visualization_pathway <- paste(df2$ID, collapse = "\n")
    
    updateTextAreaInput(session, "pathway_text", value = rv$visualization_pathway)
  })
  
  #p.adjust <
  observeEvent(input$apply_gsea_padj, {
    req(rv$gsea_obj)
    df <- rv$gsea_obj
    
    df2 <- df[df$p.adjust < input$gsea_padj, ]
    
    rv$visualization_pathway <- paste(df2$ID, collapse = "\n")
    
    updateTextAreaInput(session, "pathway_text", value = rv$visualization_pathway)
  })
  
  # select marker on gsea file
  observeEvent(input$use_gsea_table_btn, {
    library(DT)
    req(rv$gsea_obj)
    
    showModal(
      modalDialog(
        title = "Select Genes from GSEA",
        size = "l",
        easyClose = FALSE,
        
        DT::DTOutput("gsea_table"),
        
        footer = tagList(
          actionButton(
            "confirm_gsea_pathway",
            "Confirm Selection",
            class = "btn btn-primary"
          ),
          actionButton("close_modal", "Close")
        )
      )
    )
  })
  
  output$gsea_table <- DT::renderDT({
    
    req(rv$gsea_obj)
    
    df <- rv$gsea_obj
    
    selected_pathway_name <- character(0)
    
    if (!is.null(rv$visualization_pathway) && rv$visualization_pathway != "") {
      selected_pathway_name <- unlist(strsplit(rv$visualization_pathway, "\n"))
      selected_pathway_name <- trimws(selected_pathway_name)
    }
    
    valid_idx <- which(df$ID %in% selected_pathway_name)
    
    df_show <- df[, !(colnames(df) %in% c(
      "Description",
      "enrichmentScore",
      "qvalue",
      "rank",
      "leading_edge"
    )), drop = FALSE]
    
    df_show$core_enrichment <- ifelse(
      nchar(df_show$core_enrichment) > 100,
      paste0(substr(df_show$core_enrichment, 1, 100), "..."),
      df_show$core_enrichment
    )
    
    DT::datatable(
      df_show,
      class = "compact",
      
      selection = list(
        mode = "multiple",
        selected = valid_idx
      ),
      
      filter = "top",
      
      options = list(
        pageLength = 15,
        
        lengthMenu = list(
          c(10, 15, 20, 30, 50),
          c("10", "15", "20", "30", "50")
        ),
        
        scrollX = TRUE
      )
    ) %>%
      DT::formatRound(
        c("NES", "pvalue", "p.adjust"),
        2
      )
  })
  
  get_selected_pathways <- function() {
    
    req(rv$gsea_obj)
    
    sel_idx <- input$gsea_table_rows_selected
    
    if (is.null(sel_idx) || length(sel_idx) == 0) {
      return(character(0))
    }
    
    rv$gsea_obj$ID[sel_idx]
  }
  
  observeEvent(input$confirm_gsea_pathway, {
    pathways <- get_selected_pathways()
    # ---- Update reactive value ----
    rv$visualization_pathway <- paste(pathways, collapse = "\n")
    # ---- Sync back to textarea ----
    updateTextAreaInput(
      session,
      "pathway_text",
      value = rv$visualization_pathway
    )
    rv$ui_freeze = TRUE
    showModal(
      modalDialog(
        title = "Selection of pathways for analysis",
        size = "l",
        easyClose = FALSE,
        
        fluidRow(
          column(
            width = 6,
            textAreaInput(
              inputId = "pathway_text",
              label = "Pathway List:",
              value = rv$visualization_pathway,
              rows = 20,
              placeholder = "Pathway_ID_1\nPathway_ID_2",
              width = "100%"
            ),
            actionButton(
              "clear_pathway",
              "Clear All Pathways",
              class = "btn-danger",
              width = "100%"
            )
          ),
          column(
            width = 6,
            tagList(
              # --- Marker count ---
              tags$div(
                style = "margin-bottom: 10px; display: flex; align-items: center; gap: 6px;",
                tags$span("Current Number of Pathways:"),
                tags$span(
                  style = "color: red; font-weight: bold;",
                  textOutput("pathway_count")
                )
              ),
              # =========================
              # CASE 1: GSEA NOT AVAILABLE
              # =========================
              conditionalPanel(
                condition = "output.gsea_available == false",
                
                tags$div(
                  style = "padding: 15px; color: #777;",
                  tags$h5("GSEA not loaded"),
                  tags$p("Upload GSEA file to enable advanced pathway selection.")
                )
              ),
              
              # =========================
              # CASE 2: GSEA AVAILABLE
              # =========================
              conditionalPanel(
                condition = "output.gsea_available == true",
                
                # --- DEG button ---
                actionButton(
                  "use_gsea_table_btn",
                  "Use GSEA file to search for pathways",
                  class = "btn btn-primary",
                  width = "100%"
                ),
                
                tags$hr(),
                
                # ---- Top N NES ----
                fluidRow(
                  column(8,
                         numericInput("top_n_nes", "Top N |NES|:", 20, min = 1)
                  ),
                  column(4,
                         br(),
                         actionButton("apply_top_n_nes", "Apply", class = "btn-success")
                  )
                ),
                
                # ---- NES threshold ----
                fluidRow(
                  column(8,
                         numericInput("nes_threshold", "|NES| >", 1.5)
                  ),
                  column(4,
                         br(),
                         actionButton("apply_nes", "Apply", class = "btn-success")
                  )
                ),
                
                # ---- pvalue ----
                fluidRow(
                  column(8,
                         numericInput("gsea_pval", "pvalue <", 0.05)
                  ),
                  column(4,
                         br(),
                         actionButton("apply_gsea_pval", "Apply", class = "btn-success")
                  )
                ),
                # =========================
                # 4. log2fc
                # =========================
                # ---- padj ----
                fluidRow(
                  column(8,
                         numericInput("gsea_padj", "p.adjust <", 0.05)
                  ),
                  column(4,
                         br(),
                         actionButton("apply_gsea_padj", "Apply", class = "btn-success")
                  )))))),
        footer = actionButton("close_modal", "Close")
      ))
  })
  
  # ----------------- VISUALIZATION FUNCTION Pathway Rename -----------------
  observeEvent(input$pathway_rename_btn, {
    req(rv$visualization_pathway != "")
    req(rv$pathway_rename)
    
    rv$ui_freeze <- TRUE
    
    showModal(
      modalDialog(
        title = "Rename Pathways",
        
        div(
          style = "
        max-height: 70vh;
        overflow-y: auto;
        overflow-x: auto;
        white-space: nowrap;
      ",
          uiOutput("pathway_rename_ui")
        ),
        
        footer = tagList(
          actionButton("apply_pathway_rename", "Apply", class = "btn-success", width = "200px"),
          actionButton("close_modal", "Close")
        ),
        size = "l",
        easyClose = FALSE
      )
    )
  })
  
  output$pathway_rename_ui <- renderUI({
    
    req(rv$pathway_rename)
    
    df <- rv$pathway_rename
    
    tagList(
      
      tags$table(
        id = "rename_table",  
        style = "width: auto; border-collapse: collapse; table-layout: fixed; border-spacing: 0;",
        
        # ===== header =====
        tags$thead(
          tags$tr(
            style = "border-bottom: 1px solid #ddd; font-weight: 700;",
            tags$th(
              style = "padding: 8px 20px 8px 0; white-space: nowrap; text-align: left;",
              "Original Pathway"
            ),
            tags$th(
              style = "padding: 8px 0; white-space: nowrap; text-align: left;",
              "Renamed Pathway"
            )
          )
        ),
        
        # ===== body =====
        tags$tbody(
          lapply(seq_len(nrow(df)), function(i) {
            tags$tr(
              style = "margin: 0; padding: 0;", 
              tags$td(
                style = "margin: 0; padding: 0; white-space: nowrap; color: #666; font-size: 13px;",
                df$original_pathway[i]
              ),
              tags$td(
                style = "margin: 0; padding-top: 10px;",  
                textInput(
                  inputId = paste0("rename_", i),
                  label = NULL,
                  value = df$renamed_pathway[i],
                  width = "100%"
                )
              )
            )
          })
        )
      ),
      
      tags$script(HTML("
      (function() {
        setTimeout(function() {
          var table = document.getElementById('rename_table');
          if (!table) {
            console.warn('Table not found');
            return;
          }
          
          var firstColCells = table.querySelectorAll('tr th:first-child, tr td:first-child');
          if (firstColCells.length === 0) return;
          
          var maxWidth = 0;
          firstColCells.forEach(function(cell) {
            var width = cell.offsetWidth;
            if (width > maxWidth) maxWidth = width;
          });
          
          maxWidth = maxWidth + 10;
          table.style.tableLayout = 'fixed';
          
          var rows = table.querySelectorAll('tr');
          rows.forEach(function(row) {
            var cells = row.querySelectorAll('th, td');
            if (cells.length >= 2) {
              cells[0].style.width = maxWidth + 'px';
              cells[0].style.minWidth = maxWidth + 'px';
              cells[0].style.maxWidth = maxWidth + 'px';
              cells[0].style.whiteSpace = 'nowrap';
              cells[1].style.width = maxWidth + 'px';
              cells[1].style.minWidth = maxWidth + 'px';
              cells[1].style.maxWidth = maxWidth + 'px';
            }
          });
          
          var inputs = table.querySelectorAll('input');
          inputs.forEach(function(input) {
            input.style.width = '100%';
            input.style.boxSizing = 'border-box';
          });
          
        }, 150);
      })();
    "))
    )
  })
  
  observeEvent(input$apply_pathway_rename, {
    
    req(rv$pathway_rename)
    
    df <- rv$pathway_rename
    
    new_names <- sapply(seq_len(nrow(df)), function(i) {
      input[[paste0("rename_", i)]]
    })
    
    rv$pathway_rename$renamed_pathway <- new_names
    
    rv$pathway_map <- setNames(
      rv$pathway_rename$renamed_pathway,
      rv$pathway_rename$original_pathway
    )
    
    removeModal()
    rv$ui_freeze <- FALSE
  })
  # ----------------- VISUALIZATION FUNCTION CCC Pathway Select -----------------
  observeEvent(input$ccc_pathway_select_btn, {
    req(rv$cccd_obj$data)
    rv$ui_freeze = TRUE
    is_compare <- !is.null(rv$cccd_obj$data_merge)
    
    get_pathway_summary <- function(cellchat_obj, group_name) {
      subsetCommunication(
        cellchat_obj,
        slot.name = "netP"
      ) %>%
        dplyr::group_by(pathway_name) %>%
        dplyr::summarise(
          total_prob = sum(prob),
          n_interactions = dplyr::n(),
          .groups = "drop"
        ) %>%
        dplyr::rename(
          !!paste0(group_name, "_total_prob") := total_prob,
          !!paste0(group_name, "_count") := n_interactions
        )
    }
    
    get_netP_matrix <- function(cellchat_obj, pathway) {
      prob_obj <- cellchat_obj@netP$prob
      
      if (is.array(prob_obj) && length(dim(prob_obj)) == 3) {
        pathway_names <- dimnames(prob_obj)[[3]]
        
        if (!(pathway %in% pathway_names)) {
          return(data.frame(Message = "Current pathway does not exist"))
        }
        
        mat <- prob_obj[, , pathway, drop = FALSE]
        mat <- mat[, , 1]
        
      } else if (is.list(prob_obj)) {
        if (!(pathway %in% names(prob_obj))) {
          return(data.frame(Message = "Current pathway does not exist"))
        }
        
        mat <- prob_obj[[pathway]]
        
      } else {
        return(data.frame(Message = "Unsupported netP$prob structure"))
      }
      
      as.data.frame.matrix(round(mat, 2))
    }
    
    if (is_compare) {
      cellchat_list <- rv$cccd_obj$data
      
      control_obj <- cellchat_list[["CONTROL"]]
      treatment_obj <- cellchat_list[["TREATMENT"]]
      
      req(control_obj)
      req(treatment_obj)
      
      control_summary <- get_pathway_summary(control_obj, "C")
      treatment_summary <- get_pathway_summary(treatment_obj, "T")
      
      pathway_summary <- dplyr::full_join(
        control_summary,
        treatment_summary,
        by = "pathway_name"
      ) %>%
        dplyr::mutate(
          C_total_prob = tidyr::replace_na(C_total_prob, 0),
          T_total_prob = tidyr::replace_na(T_total_prob, 0),
          C_count = tidyr::replace_na(C_count, 0),
          T_count = tidyr::replace_na(T_count, 0),
          
          total_prob = C_total_prob + T_total_prob,
          prob_diff = T_total_prob - C_total_prob,
          count_diff = T_count - C_count,
          
          C_total_prob = round(C_total_prob, 2),
          T_total_prob = round(T_total_prob, 2),
          total_prob = round(total_prob, 2),
          prob_diff = round(prob_diff, 2)
        ) %>%
        dplyr::arrange(dplyr::desc(total_prob))
    } else {
      
      cellchat_obj <- rv$cccd_obj$data
      
      pathway_summary <- subsetCommunication(
        cellchat_obj,
        slot.name = "netP"
      ) %>%
        dplyr::group_by(pathway_name) %>%
        dplyr::summarise(
          total_prob = sum(prob),
          n_interactions = dplyr::n(),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          total_prob = round(total_prob, 2)
        ) %>%
        dplyr::arrange(dplyr::desc(total_prob))
    }
    
    showModal(
      modalDialog(
        title = "netP Pathway Probability",
        size = "l",
        easyClose = FALSE,
        footer = tagList(
          actionButton("confirm_netP_pathway", "Confirm", class = "btn-success"),
          actionButton("close_modal", "Close")
        ),
        fluidRow(
          column(
            width = 5,
            
            textInput(
              "netP_pathway_input",
              "Pathway name",
              value = ifelse(
                nrow(pathway_summary) > 0,
                pathway_summary$pathway_name[1],
                ""
              ),
              width = "100%"
            ),
            
            if (is_compare) {
              tagList(
                h4("CONTROL"),
                div(
                  style = "
                  max-height: 300px;
                  max-width: 100%;
                  overflow-x: auto;
                  overflow-y: auto;
                  border: 1px solid #ddd;
                  padding: 8px;
                ",
                  tableOutput("netP_control_matrix")
                ),
                
                br(),
                
                h4("TREATMENT"),
                div(
                  style = "
                  max-height: 300px;
                  max-width: 100%;
                  overflow-x: auto;
                  overflow-y: auto;
                  border: 1px solid #ddd;
                  padding: 8px;
                ",
                  tableOutput("netP_treatment_matrix")
                )
              )
            } else {
              tagList(
                div(
                  style = "
                            max-height: 500px;
                            max-width: 100%;
                            overflow-x: auto;
                            overflow-y: auto;
                            border: 1px solid #ddd;
                            padding: 8px;
                          ",
                  tableOutput("netP_pathway_matrix")
                ),
                
                br(),
                
                tags$div(
                  style = "
                            background-color:#f8f9fa;
                            border-left:4px solid #007bff;
                            padding:10px;
                            font-size:13px;
                          ",
                  HTML(
                          "<b>Interpretation:</b><br/>
                          Rows = Sender cell types<br/>
                          Columns = Receiver cell types<br/>
                          Example: value at row <b>C2</b> and column <b>C0</b> means
                          <b>C2 → C0</b> communication probability."
                  )
                )
              )
            }
          ),
          
          column(
            width = 7,
            
            div(
              style = "
              max-height: 650px;
              max-width: 100%;
              overflow-x: auto;
              overflow-y: auto;
            ",
              DT::DTOutput("netP_pathway_summary_table")
            )
          )
        )
      )
    )
    
    output$netP_pathway_summary_table <- DT::renderDT({
      DT::datatable(
        pathway_summary,
        selection = list(mode = "single", target = "row"),
        rownames = FALSE,
        extensions = "FixedColumns",
        options = list(
          dom = "t",
          paging = FALSE,
          searching = FALSE,
          info = FALSE,
          ordering = TRUE,
          order = list(list(
            which(colnames(pathway_summary) == "total_prob") - 1,
            "desc"
          )),
          scrollX = TRUE,
          scrollY = "600px",
          fixedColumns = list(leftColumns = 1)
        )
      )
    })
    
    observeEvent(input$netP_pathway_summary_table_rows_selected, {
      selected_row <- input$netP_pathway_summary_table_rows_selected
      
      if (length(selected_row) == 1) {
        updateTextInput(
          session,
          "netP_pathway_input",
          value = pathway_summary$pathway_name[selected_row]
        )
      }
    })
    
    if (is_compare) {
      
      output$netP_control_matrix <- renderTable({
        req(input$netP_pathway_input)
        
        pathway <- trimws(input$netP_pathway_input)
        get_netP_matrix(control_obj, pathway)
        
      }, rownames = TRUE)
      
      output$netP_treatment_matrix <- renderTable({
        req(input$netP_pathway_input)
        
        pathway <- trimws(input$netP_pathway_input)
        get_netP_matrix(treatment_obj, pathway)
        
      }, rownames = TRUE)
      
    } else {
      
      output$netP_pathway_matrix <- renderTable({
        req(input$netP_pathway_input)
        
        pathway <- trimws(input$netP_pathway_input)
        get_netP_matrix(cellchat_obj, pathway)
        
      }, rownames = TRUE)
    }
  })
  
  observeEvent(input$confirm_netP_pathway, {
    req(input$netP_pathway_input)
    
    pathway <- trimws(input$netP_pathway_input)
    
    validate(
      need(pathway != "", "Pathway cannot be empty.")
    )
    
    rv$ccc_pathway_select <- pathway
    rv$ui_freeze = FALSE
    removeModal()
  })
  
  # ----------------- VISUALIZATION FUNCTION CCC l-R pair Select -----------------
  observeEvent(input$select_lr_pair, {
    open_lr_pair_modal(input$ccc_pair_active_tab)
  })
  
  open_lr_pair_modal <- function(target = NULL) {
    
    rv$ui_freeze <- TRUE
    
    target <- input$ccc_pair_active_tab
    
    if (is.null(target) || target == "") {
      target <- "netVisual_individual"
    }
    
    rv$ccc_lr_pair_select_target <- target
    
    if (is.null(rv$ccc_lr_pair_select)) {
      rv$ccc_lr_pair_select <- character(0)
    }
    
    
    showModal(
      modalDialog(
        title = "Select L-R pair",
        size = "l",
        easyClose = FALSE,
        
        fluidRow(
          column(
            width = 6,
            
            textAreaInput(
              inputId = "ccc_lr_pair_text",
              label = "Selected L-R pairs:",
              value = paste(unique(rv$ccc_lr_pair_select), collapse = "\n"),
              rows = 20,
              placeholder = "Selected L-R pairs will appear here...",
              width = "100%"
            ),
            
            actionButton(
              "clear_ccc_lr_pair",
              "Clear All L-R pairs",
              class = "btn-danger",
              width = "100%"
            )
          ),
          
          column(
            width = 6,
            
            tags$div(
              style = "margin-bottom: 10px; display: flex; align-items: center; gap: 6px;",
              tags$span("Current Number of Unique L-R pairs:"),
              tags$span(
                style = "color: red; font-weight: bold;",
                textOutput("ccc_lr_pair_count")
              )
            ),
            
            tags$hr(),
            actionButton(
              "use_pathway_lr_pair_btn",
              "Use pathway to select L-R pairs",
              class = "btn btn-primary",
              width = "100%"
            ),
            
            tags$hr(),
            
            actionButton(
              "use_net_lr_pair_btn",
              "Use communication table to select L-R pairs",
              class = "btn btn-primary",
              width = "100%"
            ),
            
            tags$hr(),
            
            uiOutput("ccc_lr_annotation_filter_ui"),
            
            fluidRow(
              column(
                width = 8,
                numericInput(
                  "ccc_lr_pvalue_threshold",
                  "L-R p-value < threshold:",
                  value = 0.01
                )
              ),
              column(
                width = 4,
                br(),
                actionButton(
                  "apply_ccc_lr_pvalue",
                  "Apply",
                  class = "btn btn-success",
                  width = "100%"
                )
              )
            ),
            
            fluidRow(
              column(
                width = 8,
                numericInput(
                  "ccc_lr_prob_threshold",
                  "L-R probability > threshold:",
                  value = 0.1
                )
              ),
              column(
                width = 4,
                br(),
                actionButton(
                  "apply_ccc_lr_prob",
                  "Apply",
                  class = "btn btn-success",
                  width = "100%"
                )
              )
            )
          )
        ),
        
        footer = tagList(
          actionButton("close_modal", "Close"),
          actionButton(
            "confirm_lr_pair_select",
            "Confirm",
            class = "btn-success"
          )
        )
      )
    )
  }
  
  get_current_lr_pairs <- function() {
    
    if (is.null(rv$ccc_lr_pair_select)) {
      return(character(0))
    }
    
    x <- rv$ccc_lr_pair_select
    x <- trimws(x)
    x <- x[x != ""]
    unique(x)
  }
  
  update_lr_pair_text <- function(session) {
    
    updateTextAreaInput(
      session,
      "ccc_lr_pair_text",
      value = paste(get_current_lr_pairs(), collapse = "\n")
    )
  }
  
  extract_lr_pair_name <- function(df) {
    if ("interaction_name" %in% colnames(df)) {
      return(as.character(df$interaction_name))
    }
    
    if ("interaction_name_2" %in% colnames(df)) {
      return(as.character(df$interaction_name_2))
    }
    stop("Cannot find L-R pair name column.")
  }
  
  get_net_lr_table <- function() {
    
    req(rv$cccd_obj)
    
    obj <- rv$cccd_obj$data
    
    clean_one_table <- function(df, dataset_name = NULL) {
      
      df$lr_pair <- extract_lr_pair_name(df)
      
      if ("pathway_name" %in% colnames(df)) {
        df$pathway <- df$pathway_name
      } else if (!("pathway" %in% colnames(df))) {
        df$pathway <- NA
      }
      
      if (!is.null(dataset_name)) {
        df$dataset <- dataset_name
      }
      
      remove_cols <- intersect(
        c("interaction_name", "interaction_name_2", "pathway_name"),
        colnames(df)
      )
      
      df <- df[, !(colnames(df) %in% remove_cols), drop = FALSE]
      
      first_cols <- if (!is.null(dataset_name)) {
        c("lr_pair", "pathway", "dataset")
      } else {
        c("lr_pair", "pathway")
      }
      
      other_cols <- setdiff(colnames(df), first_cols)
      
      df <- df[, c(first_cols, other_cols), drop = FALSE]
      
      df
    }
    
    if (!is.null(rv$cccd_obj$data_merge)) {
      
      object.list <- obj
      object.list <- object.list[
        sapply(object.list, function(x) methods::is(x, "CellChat"))
      ]
      
      dataset_use <- names(object.list)
      
      df_list <- lapply(dataset_use, function(dataset_name) {
        df <- subsetCommunication(object.list[[dataset_name]])
        clean_one_table(df, dataset_name = dataset_name)
      })
      
      df <- do.call(rbind, df_list)
      
    } else {
      
      df <- subsetCommunication(obj)
      df <- clean_one_table(df, dataset_name = NULL)
    }
    
    rownames(df) <- NULL
    df
  }
  
  output$ccc_lr_pair_count <- renderText({
    length(get_current_lr_pairs())
  })
  
  observeEvent(input$ccc_lr_pair_text, {
    
    x <- unlist(strsplit(input$ccc_lr_pair_text, "\n"))
    x <- trimws(x)
    x <- x[x != ""]
    
    rv$ccc_lr_pair_select <- unique(x)
    
  }, ignoreInit = TRUE)
  
  observeEvent(input$clear_ccc_lr_pair, {
    
    rv$ccc_lr_pair_select <- character(0)
    
    updateTextAreaInput(
      session,
      "ccc_lr_pair_text",
      value = ""
    )
  })
  
  observeEvent(input$use_net_lr_pair_btn, {
    
    showModal(
      modalDialog(
        title = "Select L-R pairs from communication table",
        size = "l",
        easyClose = FALSE,
        
        DT::DTOutput("ccc_net_lr_pair_table"),
        
        footer = tagList(
          actionButton(
            "confirm_net_lr_pair_select",
            "Confirm Selection",
            class = "btn btn-primary"
          ),
          actionButton("back_to_lr_pair_modal_from_net", "Back")
        )
      )
    )
  })
  
  output$ccc_net_lr_pair_table <- DT::renderDT({
    
    df <- get_net_lr_table()
    
    selected_pairs <- get_current_lr_pairs()
    
    selected_idx <- which(df$lr_pair %in% selected_pairs)
    
    dt <- DT::datatable(
      df,
      class = "compact stripe hover nowrap",
      selection = list(
        mode = "multiple",
        selected = selected_idx
      ),
      filter = "top",
      options = list(
        pageLength = 15,
        lengthMenu = list(
          c(10, 15, 20, 30, 50),
          c("10", "15", "20", "30", "50")
        ),
        scrollX = TRUE,
        autoWidth = TRUE
      )
    )
    
    if ("prob" %in% colnames(df)) {
      dt <- DT::formatRound(dt, "prob", 3)
    }
    
    dt
  })
  
  
  observeEvent(input$confirm_net_lr_pair_select, {
    
    df <- get_net_lr_table()
    sel_idx <- input$ccc_net_lr_pair_table_rows_selected
    
    if (!is.null(sel_idx) && length(sel_idx) > 0) {
      selected_pairs <- unique(df$lr_pair[sel_idx])
      
      rv$ccc_lr_pair_select <- unique(c(
        get_current_lr_pairs(),
        selected_pairs
      ))
    }
    
    removeModal()
    open_lr_pair_modal(rv$ccc_lr_pair_select_target)
  })
  
  
  observeEvent(input$back_to_lr_pair_modal_from_net, {
    removeModal()
    open_lr_pair_modal(rv$ccc_lr_pair_select_target)
  })
  
  observeEvent(input$use_pathway_lr_pair_btn, {
    
    showModal(
      modalDialog(
        title = "Select L-R pairs by pathway",
        size = "l",
        easyClose = FALSE,
        
        selectInput(
          "ccc_lr_pathway_select",
          "Select pathway:",
          choices = get_available_pathways_for_lr(),
          selected = rv$ccc_pathway_select
        ),
        
        tags$hr(),
        
        uiOutput("ccc_pathway_lr_pair_checkbox_ui"),
        
        footer = tagList(
          actionButton(
            "confirm_pathway_lr_pair_select",
            "Confirm Selection",
            class = "btn btn-primary"
          ),
          actionButton("back_to_lr_pair_modal_from_pathway", "Back")
        )
      )
    )
  })
  
  get_available_pathways_for_lr <- function() {
    
    req(rv$cccd_obj)
    
    obj <- rv$cccd_obj$data
    
    if (!is.null(rv$cccd_obj$data_merge)) {
      object.list <- obj
      object.list <- object.list[
        sapply(object.list, function(x) methods::is(x, "CellChat"))
      ]
      
      pathways <- unique(unlist(lapply(object.list, function(x) {
        if (!is.null(x@netP$pathways)) x@netP$pathways else character(0)
      })))
      
    } else {
      pathways <- obj@netP$pathways
    }
    
    sort(unique(pathways))
  }
  
  
  get_lr_pairs_by_pathway <- function(pathway_use) {
    
    req(rv$cccd_obj)
    
    obj <- rv$cccd_obj$data
    
    if (!is.null(rv$cccd_obj$data_merge)) {
      
      object.list <- obj
      object.list <- object.list[
        sapply(object.list, function(x) methods::is(x, "CellChat"))
      ]
      
      df_list <- lapply(object.list, function(x) {
        if (!(pathway_use %in% x@netP$pathways)) {
          return(NULL)
        }
        
        df <- subsetCommunication(
          x,
          signaling = pathway_use
        )
        
        if (is.null(df) || nrow(df) == 0) {
          return(NULL)
        }
        
        df$lr_pair <- extract_lr_pair_name(df)
        df
      })
      
      df <- do.call(rbind, df_list[!sapply(df_list, is.null)])
      
    } else {
      
      if (!(pathway_use %in% obj@netP$pathways)) {
        return(character(0))
      }
      
      df <- subsetCommunication(
        obj,
        signaling = pathway_use
      )
      
      if (is.null(df) || nrow(df) == 0) {
        return(character(0))
      }
      
      df$lr_pair <- extract_lr_pair_name(df)
    }
    
    unique(df$lr_pair)
  }
  
  
  output$ccc_pathway_lr_pair_checkbox_ui <- renderUI({
    
    req(input$ccc_lr_pathway_select)
    
    lr_pairs <- get_lr_pairs_by_pathway(input$ccc_lr_pathway_select)
    
    if (length(lr_pairs) == 0) {
      return(
        div(
          style = "padding: 15px; color: red;",
          paste0("No L-R pairs found for pathway: ", input$ccc_lr_pathway_select)
        )
      )
    }
    
    checkboxGroupInput(
      "ccc_pathway_lr_pair_select",
      "Select L-R pairs:",
      choices = lr_pairs,
      selected = lr_pairs
    )
  })
  
  
  observeEvent(input$confirm_pathway_lr_pair_select, {
    
    selected_pairs <- input$ccc_pathway_lr_pair_select
    
    if (!is.null(selected_pairs) && length(selected_pairs) > 0) {
      rv$ccc_lr_pair_select <- unique(c(
        get_current_lr_pairs(),
        selected_pairs
      ))
    }
    
    removeModal()
    open_lr_pair_modal(rv$ccc_lr_pair_select_target)
  })
  
  
  observeEvent(input$back_to_lr_pair_modal_from_pathway, {
    removeModal()
    open_lr_pair_modal(rv$ccc_lr_pair_select_target)
  })
  
  observeEvent(input$apply_ccc_lr_pvalue, {
    
    df <- get_net_lr_table()
    
    if (!("pval" %in% colnames(df))) {
      showNotification("Column 'pval' was not found in communication table.", type = "error")
      return()
    }
    
    selected_pairs <- unique(df$lr_pair[df$pval < input$ccc_lr_pvalue_threshold])
    
    rv$ccc_lr_pair_select <- unique(c(
      get_current_lr_pairs(),
      selected_pairs
    ))
    
    update_lr_pair_text(session)
  })
  
  
  observeEvent(input$apply_ccc_lr_prob, {
    
    df <- get_net_lr_table()
    
    if (!("prob" %in% colnames(df))) {
      showNotification("Column 'prob' was not found in communication table.", type = "error")
      return()
    }
    
    selected_pairs <- unique(df$lr_pair[df$prob > input$ccc_lr_prob_threshold])
    
    rv$ccc_lr_pair_select <- unique(c(
      get_current_lr_pairs(),
      selected_pairs
    ))
    
    update_lr_pair_text(session)
  })
  
  observeEvent(input$confirm_lr_pair_select, {
    
    rv$ccc_lr_pair_select <- get_current_lr_pairs()
    rv$ui_freeze <- FALSE
    removeModal()
  })
  
  # ----------------- VISUALIZATION FUNCTION Original Violin Plot -----------------

  visualization_data_prepared <- reactive({
    
    req(rv$data_obj)
    req(rv$visualization_marker)
    req(rv$visualization_sample_column != "")
    req(!rv$ui_freeze)
    showModal(modalDialog(
      title = "Processing Data for Visualization",
      tags$div("Preparing ..."),
      easyClose = FALSE,
      footer = NULL
    ))
    sample_col <- rv$visualization_sample_column
    use_split   <- rv$visualization_use_split
    split_col   <- rv$visualization_split_column
    
    markers <- unlist(strsplit(rv$visualization_marker, "\n"))
    markers <- trimws(markers)
    markers <- markers[markers != "" & !is.na(markers)]
    valid_genes <- rownames(rv$data_obj)
    genes <- markers[markers %in% valid_genes][1]
    
    if (isTRUE(use_split) && !is.null(split_col) && split_col != "") {
      
      meta_vec <- rv$data_obj@meta.data[[split_col]]
      
      if (is.factor(meta_vec)) {
        split_vals <- levels(meta_vec)
        split_vals <- split_vals[split_vals %in% unique(meta_vec)]
      } else {
        split_vals <- sort(unique(meta_vec))
        split_vals <- split_vals[!is.na(split_vals)]
      }
      
      split_data <- lapply(split_vals, function(split_val) {
        subset(rv$data_obj, subset = !!as.name(split_col) == split_val)
      })
      removeModal()
      list(
        gene = genes,
        split_data = split_data,
        split_vals = split_vals,
        has_split = TRUE,
        sample_col = sample_col,
        split_col = split_col
      )
      
    } else {
      removeModal()
      list(
        gene = genes,
        data = rv$data_obj,
        has_split = FALSE,
        sample_col = sample_col,
        split_col = NULL
      )
    }
  })
  
  output$original_violin_plot_ui <- renderUI({
    
    if (!is.null(rv$visualization_split_column) && rv$visualization_split_column != "") {
      n_splits <- length(unique(rv$data_obj@meta.data[[rv$visualization_split_column]]))
      width_val <- paste0(n_splits * 350, "px")  
    } else {
      width_val <- "384px"
    }
    
    div(
      style = "max-width: 800px; overflow-x: auto; width: 100%;",
      plotOutput("original_violin_plot", height = "288px", width = width_val)
    )
  })
  
  output$original_violin_plot <- renderPlot({
    
    data_prep <- visualization_data_prepared()
    genes <- data_prep$gene
    
    if (data_prep$has_split) {
      plot_list <- lapply(seq_along(data_prep$split_data), function(i) {
        subset_obj <- data_prep$split_data[[i]]
        split_val <- data_prep$split_vals[i]
        
        x_labels <- names(table(subset_obj@meta.data[[rv$visualization_sample_column]]))
        
        if (input$v_show_non_zero) {
          percentages <- get_non_zero_percentage(
            subset_obj, 
            genes, 
            rv$visualization_sample_column
          )
          
          if (!is.null(percentages)) {
            x_labels <- sapply(x_labels, function(label) {
              pct <- percentages[label]
              if (!is.na(pct)) {
                paste0(label, "\n(", sprintf("%.1f", pct), "%)")
              } else {
                label
              }
            })
          }
        }
        
        VlnPlot(
          subset_obj,
          features = genes,
          group.by = rv$visualization_sample_column,
          pt.size = if (input$v_show_points) 0.5 else 0
        ) +
          ggplot2::scale_x_discrete(labels = x_labels) +
          ggplot2::labs(
            title = paste0(genes, " (", split_val, ")"),
            x = rv$visualization_sample_column,
            y = "Expression"
          ) +
          ggplot2::theme(
            axis.title.x = ggplot2::element_text(size = input$v_x_lab_size),
            axis.title.y = ggplot2::element_text(size = input$v_y_lab_size),
            axis.text.x = ggplot2::element_text(
              size = input$v_x_text_size,
              angle = input$v_angle, 
              hjust = 0.5, vjust = 0.5
            ),
            axis.text.y = ggplot2::element_text(
              size = input$v_y_text_size
            ),
            plot.title = ggplot2::element_text(
              size = input$v_title_size,
              hjust = 0.5
            ),
            legend.position = "none"
          )
      })
      
      patchwork::wrap_plots(plot_list, ncol = length(data_prep$split_vals))
      
    } else {
      x_labels <- names(table(data_prep$data@meta.data[[rv$visualization_sample_column]]))
      
      if (input$v_show_non_zero) {
        percentages <- get_non_zero_percentage(
          data_prep$data, 
          genes, 
          rv$visualization_sample_column
        )
        
        if (!is.null(percentages)) {
          x_labels <- sapply(x_labels, function(label) {
            pct <- percentages[label]
            if (!is.na(pct)) {
              paste0(label, "\n(", sprintf("%.1f", pct), "%)")
            } else {
              label
            }
          })
        }
      }
      
      VlnPlot(
        data_prep$data,
        features = genes,
        group.by = rv$visualization_sample_column,
        pt.size = if (input$v_show_points) 0.5 else 0
      ) +
        ggplot2::scale_x_discrete(labels = x_labels) +
        ggplot2::labs(
          title = paste0(genes,"_example"),
          x = rv$visualization_sample_column,
          y = "Expression"
        ) +
        ggplot2::theme(
          axis.title.x = ggplot2::element_text(size = input$v_x_lab_size),
          axis.title.y = ggplot2::element_text(size = input$v_y_lab_size),
          axis.text.x = ggplot2::element_text(
            size = input$v_x_text_size,
            angle = input$v_angle, 
            hjust = 0.5, vjust = 0.5
          ),
          axis.text.y = ggplot2::element_text(
            size = input$v_y_text_size
          ),
          plot.title = ggplot2::element_text(
            size = input$v_title_size,
            hjust = 0.5
          ),
          legend.position = "none"
        )
    }
  }, height = function() 300)
  
  get_non_zero_percentage <- function(data_obj, feature, group_column) {
    tryCatch({
      expr_matrix <- GetAssayData(data_obj, layer = "data")
      
      if (!feature %in% rownames(expr_matrix)) {
        return(NULL)
      }
      
      gene_expr <- expr_matrix[feature, ]
      groups <- data_obj@meta.data[[group_column]]
      unique_groups <- unique(groups[!is.na(groups)])
      
      percentages <- sapply(unique_groups, function(g) {
        group_mask <- groups == g
        group_expr <- gene_expr[group_mask]
        non_zero_count <- sum(group_expr > 0)
        total_count <- sum(group_mask)
        percentage <- (non_zero_count / total_count) * 100
        percentage
      })
      
      names(percentages) <- unique_groups
      percentages
    }, error = function(e) {
      warning(paste("Error calculating non-zero percentage:", e$message))
      return(NULL)
    })
  }
  
  combine_pngs_to_pdf_violin <- function(png_files, output_path,
                                         has_split = FALSE,
                                         n_splits = 1,
                                         genes_per_page = 8,
                                         genes_per_row = 2,
                                         final_width = 8.5,
                                         final_height = 11,
                                         update_progress_fn = NULL) {
    
    library(magick)
    
    if (length(png_files) == 0) {
      stop("No PNG files to combine")
    }
    
    if (!is.null(update_progress_fn)) {
      update_progress_fn("<b>Combining images into PDF...</b><br>Reading PNG files...")
    }
    
    tryCatch({
      
      all_images <- list()
      for (i in seq_along(png_files)) {
        if (!is.null(update_progress_fn) && i %% 5 == 0) {
          update_progress_fn(paste0(
            "<b>Combining images into PDF...</b><br>",
            "Loading images: ", i, " / ", length(png_files)
          ))
        }
        
        tryCatch({
          all_images[[i]] <- magick::image_read(png_files[i])
        }, error = function(e) {
          warning(paste("Failed to read", png_files[i]))
        })
      }
      
      if (!is.null(update_progress_fn)) {
        update_progress_fn("<b>Combining images into PDF...</b><br>Creating pages...")
      }
      
      page_images <- list()
      num_pages <- ceiling(length(all_images) / genes_per_page)
      
      for (page in seq_len(num_pages)) {
        start_idx <- (page - 1) * genes_per_page + 1
        end_idx <- min(page * genes_per_page, length(all_images))
        
        page_imgs <- all_images[start_idx:end_idx]
        page_imgs <- page_imgs[!sapply(page_imgs, is.null)]
        
        if (length(page_imgs) == 0) next
        
        if (has_split) {
          combined <- page_imgs[[1]]
          if (length(page_imgs) > 1) {
            for (k in 2:length(page_imgs)) {
              combined <- magick::image_append(c(combined, page_imgs[[k]]), stack = TRUE)
            }
          }
        } else {
          rows <- list()
          n_rows <- ceiling(length(page_imgs) / genes_per_row)
          
          for (r in seq_len(n_rows)) {
            row_start <- (r - 1) * genes_per_row + 1
            row_end <- min(r * genes_per_row, length(page_imgs))
            row_imgs <- page_imgs[row_start:row_end]
            
            row_combined <- row_imgs[[1]]
            if (length(row_imgs) > 1) {
              for (k in 2:length(row_imgs)) {
                row_combined <- magick::image_append(c(row_combined, row_imgs[[k]]), stack = FALSE)
              }
            }
            
            if (length(row_imgs) == 1) {
              info <- magick::image_info(row_combined)
              white_space <- magick::image_blank(
                width = info$width,
                height = info$height,
                color = "white"
              )
              row_combined <- magick::image_append(c(row_combined, white_space), stack = FALSE)
            }
            
            rows[[r]] <- row_combined
          }
          
          combined <- rows[[1]]
          if (length(rows) > 1) {
            for (r in 2:length(rows)) {
              combined <- magick::image_append(c(combined, rows[[r]]), stack = TRUE)
            }
          }
        }
        
        combined <- magick::image_scale(
          combined,
          geometry = paste0(
            as.integer(final_width * 300), "x",
            as.integer(final_height * 300)
          )
        )
        
        page_images[[page]] <- combined
      }
      
      if (!is.null(update_progress_fn)) {
        update_progress_fn("<b>Combining images into PDF...</b><br>Writing PDF...")
      }
      
      final_images <- page_images[[1]]
      if (length(page_images) > 1) {
        for (p in 2:length(page_images)) {
          final_images <- c(final_images, page_images[[p]])
        }
      }
      
      magick::image_write(final_images, output_path, format = "pdf")
      
      if (!is.null(update_progress_fn)) {
        update_progress_fn("<b>PDF created successfully!</b>")
      }
      
      output_path
      
    }, error = function(e) {
      stop(paste("Error combining PNGs:", e$message))
    })
  }
  
  observeEvent(input$generate_all_violin, {
    
    req(rv$data_obj)
    
    showModal(modalDialog(
      title = "Processing Genes",
      tags$div(id = "violin_progress_detail", "Starting..."),
      easyClose = FALSE,
      footer = NULL
    ))
    
    markers <- unlist(strsplit(rv$visualization_marker, "\n"))
    markers <- trimws(markers)
    markers <- markers[markers != "" & !is.na(markers)]
    valid_genes <- rownames(rv$data_obj)
    
    markers <- markers[markers %in% valid_genes]
    if (length(markers) == 0) {
      removeModal()
      showNotification("No valid genes found in Seurat object", type = "error")
      return()
    }
    
    has_split <- !is.null(rv$visualization_split_column) && rv$visualization_split_column != ""
    
    if (has_split) {
      split_vals <- unique(rv$data_obj@meta.data[[rv$visualization_split_column]])
      split_vals <- split_vals[!is.na(split_vals)]
      n_splits <- length(split_vals)
      
      split_data_list <- lapply(split_vals, function(split_val) {
        subset(rv$data_obj, subset = !!as.name(rv$visualization_split_column) == split_val)
      })
    } else {
      n_splits <- 1
      split_data_list <- NULL
      split_vals <- NULL
    }
    
    total <- length(markers)
    
    out_dir <- input$violin_save_path
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    
    generated_files <- c()
    png_dir <- file.path(out_dir, ".temp_png")
    if (!dir.exists(png_dir)) dir.create(png_dir, recursive = TRUE)
    png_files <- c()
    
    for (i in seq_along(markers)) {
      genes <- markers[i]
      done <- i - 1
      remaining <- total - i + 1
      
      shinyjs::html(
        "violin_progress_detail",
        paste0(
          "<b>Generating violin plots:</b> ", genes, "<br>",
          "Completed: ", done, "<br>",
          "Remaining: ", remaining
        )
      )
      
      Sys.sleep(0.01)
      
      png_file <- file.path(png_dir, paste0(genes, ".png"))
      pdf_file <- file.path(out_dir, paste0(genes, ".pdf"))
      
      tryCatch({
        
        # ---------- PNG ----------
        if (has_split) {
          png_width <- input$v_width * n_splits * 300
          png_height <- input$v_height * 300
          
          png(png_file, width = png_width, height = png_height, res = 300)
          
          plot_list <- lapply(seq_along(split_data_list), function(j) {
            subset_obj <- split_data_list[[j]]
            split_val <- split_vals[j]
            
            x_labels <- names(table(subset_obj@meta.data[[rv$visualization_sample_column]]))
            
            if (input$v_show_non_zero) {
              percentages <- get_non_zero_percentage(
                subset_obj, 
                genes, 
                rv$visualization_sample_column
              )
              
              if (!is.null(percentages)) {
                x_labels <- sapply(x_labels, function(label) {
                  pct <- percentages[label]
                  if (!is.na(pct)) {
                    paste0(label, "\n(", sprintf("%.1f", pct), "%)")
                  } else {
                    label
                  }
                })
              }
            }
            
            VlnPlot(
              subset_obj,
              features = genes,
              group.by = rv$visualization_sample_column,
              pt.size = if (input$v_show_points) 0.5 else 0
            ) +
              ggplot2::scale_x_discrete(labels = x_labels) +
              ggplot2::labs(
                title = paste0(genes, " (", split_val, ")"),
                x = rv$visualization_sample_column,
                y = "Expression"
              ) +
              ggplot2::theme(
                axis.title.x = ggplot2::element_text(size = input$v_x_lab_size),
                axis.title.y = ggplot2::element_text(size = input$v_y_lab_size),
                axis.text.x = ggplot2::element_text(
                  size = input$v_x_text_size,
                  angle = input$v_angle,
                  hjust = 0.5, vjust = 0.5
                ),
                axis.text.y = ggplot2::element_text(size = input$v_y_text_size),
                plot.title = ggplot2::element_text(size = input$v_title_size, hjust = 0.5),
                legend.position = "none"
              )
          })
          
          print(patchwork::wrap_plots(plot_list, ncol = n_splits))
          
        } else {
          png(png_file, width = input$v_width * 300, height = input$v_height * 300, res = 300)
          
          x_labels <- names(table(rv$data_obj@meta.data[[rv$visualization_sample_column]]))
          
          if (input$v_show_non_zero) {
            percentages <- get_non_zero_percentage(
              rv$data_obj, 
              genes, 
              rv$visualization_sample_column
            )
            
            if (!is.null(percentages)) {
              x_labels <- sapply(x_labels, function(label) {
                pct <- percentages[label]
                if (!is.na(pct)) {
                  paste0(label, "\n(", sprintf("%.1f", pct), "%)")
                } else {
                  label
                }
              })
            }
          }
          
          print(
            VlnPlot(
              rv$data_obj,
              features = genes,
              group.by = rv$visualization_sample_column,
              pt.size = if (input$v_show_points) 0.5 else 0
            ) +
              ggplot2::scale_x_discrete(labels = x_labels) +
              ggplot2::labs(
                title = genes,
                x = rv$visualization_sample_column,
                y = "Expression"
              ) +
              ggplot2::theme(
                axis.title.x = ggplot2::element_text(size = input$v_x_lab_size),
                axis.title.y = ggplot2::element_text(size = input$v_y_lab_size),
                axis.text.x = ggplot2::element_text(
                  size = input$v_x_text_size,
                  angle = input$v_angle,
                  hjust = 0.5, vjust = 0.5
                ),
                axis.text.y = ggplot2::element_text(size = input$v_y_text_size),
                plot.title = ggplot2::element_text(size = input$v_title_size, hjust = 0.5),
                legend.position = "none"
              )
          )
        }
        
        dev.off()
        
        pdf(pdf_file,
            width = if (has_split) input$v_width * n_splits else input$v_width,
            height = input$v_height)
        
        if (has_split) {
          plot_list <- lapply(seq_along(split_data_list), function(j) {
            subset_obj <- split_data_list[[j]]
            split_val <- split_vals[j]
            
            x_labels <- names(table(subset_obj@meta.data[[rv$visualization_sample_column]]))
            
            if (input$v_show_non_zero) {
              percentages <- get_non_zero_percentage(
                subset_obj, 
                genes, 
                rv$visualization_sample_column
              )
              
              if (!is.null(percentages)) {
                x_labels <- sapply(x_labels, function(label) {
                  pct <- percentages[label]
                  if (!is.na(pct)) {
                    paste0(label, "\n(", sprintf("%.1f", pct), "%)")
                  } else {
                    label
                  }
                })
              }
            }
            
            VlnPlot(
              subset_obj,
              features = genes,
              group.by = rv$visualization_sample_column,
              pt.size = if (input$v_show_points) 0.5 else 0
            ) +
              ggplot2::scale_x_discrete(labels = x_labels) +
              ggplot2::labs(
                title = paste0(genes, " (", split_val, ")"),
                x = rv$visualization_sample_column,
                y = "Expression"
              ) +
              ggplot2::theme(
                axis.title.x = ggplot2::element_text(size = input$v_x_lab_size),
                axis.title.y = ggplot2::element_text(size = input$v_y_lab_size),
                axis.text.x = ggplot2::element_text(
                  size = input$v_x_text_size,
                  angle = input$v_angle,
                  hjust = 0.5, vjust = 0.5
                ),
                axis.text.y = ggplot2::element_text(size = input$v_y_text_size),
                plot.title = ggplot2::element_text(size = input$v_title_size, hjust = 0.5),
                legend.position = "none"
              )
          })
          
          print(patchwork::wrap_plots(plot_list, ncol = n_splits))
          
        } else {
          x_labels <- names(table(rv$data_obj@meta.data[[rv$visualization_sample_column]]))
          
          if (input$v_show_non_zero) {
            percentages <- get_non_zero_percentage(
              rv$data_obj, 
              genes, 
              rv$visualization_sample_column
            )
            
            if (!is.null(percentages)) {
              x_labels <- sapply(x_labels, function(label) {
                pct <- percentages[label]
                if (!is.na(pct)) {
                  paste0(label, "\n(", sprintf("%.1f", pct), "%)")
                } else {
                  label
                }
              })
            }
          }
          
          print(
            VlnPlot(
              rv$data_obj,
              features = genes,
              group.by = rv$visualization_sample_column,
              pt.size = if (input$v_show_points) 0.5 else 0
            ) +
              ggplot2::scale_x_discrete(labels = x_labels) +
              ggplot2::labs(
                title = genes,
                x = rv$visualization_sample_column,
                y = "Expression"
              ) +
              ggplot2::theme(
                axis.title.x = ggplot2::element_text(size = input$v_x_lab_size),
                axis.title.y = ggplot2::element_text(size = input$v_y_lab_size),
                axis.text.x = ggplot2::element_text(
                  size = input$v_x_text_size,
                  angle = input$v_angle,
                  hjust = 0.5, vjust = 0.5
                ),
                axis.text.y = ggplot2::element_text(size = input$v_y_text_size),
                plot.title = ggplot2::element_text(size = input$v_title_size, hjust = 0.5),
                legend.position = "none"
              )
          )
        }
        dev.off()
        
        png_files <- c(png_files, png_file)
        generated_files <- c(generated_files, pdf_file)
        
      }, error = function(e) {
        warning(paste("Failed to generate plot for", genes, ":", e$message))
      })
    }
    
    rv$violin_files <- generated_files
    
    tryCatch({
      shinyjs::html(
        "violin_progress_detail",
        "<b>Combined file</b><br>Creating multi-page PDF..."
      )
      
      combined_output <- file.path(
        out_dir,
        paste0("combined_violin_plots_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
      )
      
      combine_pngs_to_pdf_violin(
        png_files = png_files,
        output_path = combined_output,
        has_split = has_split,
        n_splits = n_splits,
        genes_per_page = if (has_split) 4 else 8,
        genes_per_row = 2,
        final_width = 8.5,
        final_height = if (has_split) 8.5 else 11,
        update_progress_fn = function(msg) {
          shinyjs::html("violin_progress_detail", paste0("<b>Combined file</b><br>", msg))
        }
      )
      
      rv$violin_files <- c(rv$violin_files, combined_output)
      
      unlink(png_dir, recursive = TRUE)
      removeModal()
      
      showNotification(
        paste0(
          "All violin plots generated and combined!",
          "\nSplit groups: ", n_splits,
          "\nOutput directory: ", out_dir
        ),
        type = "message",
        duration = 10
      )
      
    }, error = function(e) {
      removeModal()
      showNotification(
        paste("Error combining PDFs:", e$message),
        type = "error"
      )
    })
  })
  
  output$download_violin <- downloadHandler(
    
    filename = function() {
      "violin_plots.zip"
    },
    
    content = function(file) {
      
      files <- rv$violin_files
      
      if (is.null(files) || length(files) == 0) {
        showNotification(
          "No violin plots generated yet.",
          type = "error",
          duration = 5
        )
        return(NULL)
      }
      
      old <- getwd()
      on.exit(setwd(old), add = TRUE)
      
      setwd(dirname(files[1]))
      
      zip::zip(
        zipfile = file,
        files = basename(files)
      )
      
      showNotification(
        "Violin plots downloaded successfully.",
        type = "message",
        duration = 3
      )
    }
  )
  
  # ----------------- VISUALIZATION FUNCTION Stack Violin Plot -----------------
  observeEvent(input$set_group_gene, {
    
    req(rv$data_obj)
    req(rv$visualization_sample_column)
    rv$ui_freeze = TRUE
    sample_group <- rv$data_obj@meta.data[[rv$visualization_sample_column]]
    
    groups <- unique(as.character(sample_group))
    
    if (is.factor(sample_group)) {
      groups <- levels(sample_group)
    }
    
    if (!is.null(rv$gene_grouping)) {
      colors <- sapply(seq_along(groups), function(i) {
        g <- groups[i]
        idx <- which(sapply(rv$gene_grouping, function(x) x$group == g))
        if (length(idx) > 0) rv$gene_grouping[[idx]]$color else color_selected(length(groups))[i]
      })
    } else {
      colors <- color_selected(length(groups))
    }
    
    all_genes <- unique(trimws(unlist(strsplit(rv$visualization_marker, "\n"))))
    all_genes <- all_genes[all_genes != ""]
    
    assigned_genes <- c()
    
    group_gene_map <- lapply(seq_along(groups), function(i) {
      g <- groups[i]
      
      if (!is.null(rv$gene_grouping)) {
        idx <- which(sapply(rv$gene_grouping, function(x) x$group == g))
        if (length(idx) > 0) {
          genes <- rv$gene_grouping[[idx]]$genes
          assigned_genes <<- c(assigned_genes, genes)
          return(genes)
        }
      }
      
      return(character(0))
    })
    
    gene_pool <- setdiff(all_genes, assigned_genes)
    showModal(modalDialog(
      title = "Gene Grouping & Color Setting",
      size = "l",
      easyClose = FALSE,
      
      div(
        style = "overflow-x: auto; width: 100%; padding-bottom: 10px;",
        
        div(
          style = "display: flex; gap: 20px; margin-bottom: 10px;",
          lapply(seq_along(groups), function(i) {
            div(
              style = "flex: 0 0 200px;",
              tags$b(groups[i])
            )
          })
        ),
        
        div(
          style = "display: flex; gap: 20px; margin-bottom: 10px;",
          lapply(seq_along(groups), function(i) {
            div(
              style = "width: 200px;",
              colourInput(
                inputId = paste0("group_color_", i),
                label = NULL,
                value = colors[i]
              )
            )
          })
        ),
        div(
          style = "display: flex; gap: 20px;",
          lapply(seq_along(groups), function(i) {
            div(
              style = "flex: 0 0 200px; border: 1px solid #ddd; padding: 5px;",
              rank_list(
                text = NULL,
                labels = group_gene_map[[i]],
                input_id = paste0("group_genes_", i),
                options = sortable_options(
                  group = "genes_shared",
                  animation = 150
                )
              )
            )
          })
        )
      ),
      
      tags$hr(),
      
      tags$h5("Gene Pool (Drag genes into groups)"),
      
      div(
        style = "
      height: 200px;
      overflow-y: auto;
      border: 1px solid #ddd;
      padding: 5px;
      background: white;
    ",
        
        rank_list(
          text = NULL,
          labels = gene_pool,
          input_id = "gene_pool",
          options = sortable_options(
            group = "genes_shared",
            animation = 150
          )
        )
      ),
      
      footer = tagList(
        fluidRow(
          column(12,
                 actionButton(
                   "save_group_gene",
                   "Confirm and Close",
                   class = "btn-success"
                 )
          )
        )
      )
    ))
  })
  
  observeEvent(input$save_group_gene, {
    
    req(rv$data_obj)
    groups <- unique(rv$data_obj@meta.data[[rv$visualization_sample_column]])
    
    res <- lapply(seq_along(groups), function(i) {
      
      list(
        group = groups[i],
        genes = input[[paste0("group_genes_", i)]],
        color = input[[paste0("group_color_", i)]]
      )
    })
    
    rv$gene_grouping <- res
    rv$ui_freeze = FALSE
    removeModal()
    
    showNotification("Gene grouping saved!", type = "message")
  })
  
  stack_violin_data <- reactive({
    
    req(rv$data_obj)
    req(rv$visualization_marker)
    req(rv$visualization_sample_column)
    req(rv$gene_grouping)
    req(!rv$ui_freeze)
    genes <- unique(unlist(lapply(rv$gene_grouping, function(x) x$genes)))
    
    genes <- trimws(genes)
    genes <- genes[genes != "" & !is.na(genes)]
    
    genes <- genes[genes %in% rownames(rv$data_obj)]
    
    validate(
      need(length(genes) > 0, "No valid genes in gene_grouping")
    )
    expr <- GetAssayData(
      rv$data_obj,
      assay = "RNA",
      layer = "data"
    )[genes, , drop = FALSE]
    
    expr <- as.data.frame(t(expr))
    
    cluster_vec <- rv$data_obj@meta.data[[rv$visualization_sample_column]]
    
    if (is.factor(cluster_vec)) {
      cluster_levels <- levels(cluster_vec)
      cluster_vec <- factor(cluster_vec, levels = cluster_levels)
    } else {
      cluster_levels <- sort(unique(cluster_vec))
      cluster_vec <- factor(cluster_vec, levels = cluster_levels)
    }
    
    expr$cluster <- cluster_vec
    
    df <- expr %>%
      tibble::as_tibble() %>%
      tidyr::pivot_longer(
        cols = -cluster,
        names_to = "Feat",
        values_to = "Expr"
      )
    
    df$Expr <- as.numeric(df$Expr)
    
    noise <- rnorm(n = nrow(df), mean = 0, sd = 1e-5)
    
    # only add noise if not constant
    if (!all(df$Expr == df$Expr[1])) {
      df$Expr <- df$Expr + noise
    } else {
      warning("All cells have identical expression for plotted features.")
    }
    
    gene_grouping <- rv$gene_grouping
    
    group_map <- do.call(rbind, lapply(gene_grouping, function(g) {
      
      if (is.null(g$genes) || length(g$genes) == 0) return(NULL)
      
      data.frame(
        gene = g$genes,
        group = g$group,
        color = g$color,
        stringsAsFactors = FALSE
      )
    }))
    
    group_map <- group_map[!is.na(group_map$gene), ]
    group_levels <- unique(sapply(rv$gene_grouping, function(g) as.character(g$group)))
    group_map$group <- factor(as.character(group_map$group), levels = rev(group_levels))
    
    df <- dplyr::left_join(df, group_map, by = c("Feat" = "gene"))
    
    df$group[is.na(df$group)] <- "Unknown"
    df$color[is.na(df$color)] <- "#cccccc"
    df$Feat <- factor(df$Feat, levels = genes)
    df$cluster <- factor(df$cluster, levels = cluster_levels)
    group_levels <- unique(sapply(rv$gene_grouping, function(g) g$group))
    
    gene_levels <- unique(unlist(lapply(rv$gene_grouping, function(g) {
      g$genes
    })))
    
    gene_levels <- trimws(gene_levels)
    gene_levels <- gene_levels[gene_levels != "" & !is.na(gene_levels)]
    gene_levels <- gene_levels[gene_levels %in% genes]
    
    df_com <- group_map
    df_com$group <- factor(df_com$group, levels = rev(levels(group_levels)))
    
    df_com$gene <- factor(df_com$gene, levels = gene_levels)
    df_com <- df_com[order(df_com$group, decreasing = TRUE), ]
    df_com$gene <- factor(df_com$gene, levels = df_com$gene)
    return(list(df = df,
           df_com = df_com))
  })
    
  stack_violin_plot <- reactive({
      req(stack_violin_data())
    df <- stack_violin_data()$df
    df_com <- stack_violin_data()$df_com
    right <- ggplot(df, aes(factor(cluster), Expr, fill = Feat)) + 
      geom_violin(scale = "width", adjust = 1, trim = TRUE) +
      scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
        c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
      scale_fill_manual(
        values = setNames(df$color[match(levels(df$Feat), df$Feat)], levels(df$Feat))
      ) +
      facet_grid(rows = vars(Feat), scales = "free_y", switch = "y") +
      theme_cowplot(font_size = input$stack_violin_y_right_label_size) +# Expression level size 
      theme(legend.position = "none", panel.spacing = unit(0, "lines"),
            plot.title = element_blank(),
            panel.background = element_rect(fill = NA, color = "black"),
            strip.background = element_blank(),
            strip.text = element_text(face = "bold"),
            strip.text.y.left = element_blank(),
            axis.text.y.right = element_text(size = input$stack_violin_y_right_value_size), # Expression level value size 
            axis.text.x = element_text(size = input$stack_violin_x_size), # x sample id size
            axis.title.x =  element_blank(),
            plot.margin = margin(0, 0, 0, 0, "pt"))  + 
      xlab(NULL) + ylab("Expression Level")
    left <- ggplot(df_com, aes(x = 1, y = gene, fill = group, label = group)) + geom_tile(NULL) +
      geom_text(fontface = "bold", size = input$stack_violin_y_left_sample_size) + # y left sample id size
      theme_bw(base_size = 12) + 
      scale_y_discrete(limits = rev, expand = expansion(add = c(0, 0))) +
      scale_fill_manual(values = setNames(df_com$color, df_com$group)) + 
      scale_x_continuous(expand = c(0, 0)) +
      theme(legend.position = "none", panel.spacing = unit(0, "lines"),
            panel.background = element_blank(), 
            panel.border = element_blank(),
            plot.background = element_blank(), 
            plot.margin = margin(0, 4, 0, 0, "pt"),
            axis.text.y = element_text(size = input$stack_violin_y_left_gene_size, angle = 0, hjust = 1, vjust = 0.5, color = "black"), # y left gene size
            axis.title.y = element_text(size = input$stack_violin_y_left_label_size), # y left label size
            axis.text.x = element_text(size = input$stack_violin_x_size, angle = 0, hjust = 1, vjust = 0.5, color = "white"),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid = element_blank()) + ylab("Feature") + xlab(NULL) 
    
    cowplot::plot_grid(
      left, right,
      rel_widths = c(input$stack_violin_left_proportion, 1-input$stack_violin_left_proportion), # left proportion
      align = "h",
      axis = "lc"
    )
  })
  
  output$stack_violin_plot <- renderPlot({
    req(stack_violin_plot())
    stack_violin_plot()
  })
  
  observeEvent(input$save_stack_plot, {
    
    req(stack_violin_plot())
    out_dir <- input$stack_save_path
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    out_file <- file.path(
      out_dir,
      paste0("stack_violin_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
    )
    
    pdf(out_file,
        width = input$stack_width,
        height = input$stack_height)
    
    print(stack_violin_plot())
    
    dev.off()
    
    showNotification("Saved stack violin plot", type = "message")
  })
  
  output$download_stack_plot <- downloadHandler(
    
    filename = function() {
      paste0("stack_violin_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
    },
    
    content = function(file) {
      
      req(stack_violin_plot())
      pdf(file,
          width = input$stack_width,
          height = input$stack_height)
      
      print(stack_violin_plot())
      
      dev.off()
    }
  )
  # ----------------- VISUALIZATION FUNCTION Dot Plot -----------------

  dotplot_plot_reactive <- reactive({
    
    data_prep <- visualization_data_prepared()
    
    req(data_prep)
    
    genes <- unique(trimws(unlist(strsplit(rv$visualization_marker, "\n"))))
    genes <- genes[genes != ""]
    genes <- genes[genes %in% rownames(rv$data_obj)]
    
    validate(need(length(genes) > 0, "No valid genes found."))
    
    # Check if split is enabled
    if (data_prep$has_split) {
      # Calculate global min/max for both percent.exp and avg.exp
      expr_stats <- lapply(data_prep$split_data, function(subset_obj) {
        # Get DotPlot data
        dot_data <- Seurat::DotPlot(
          object = subset_obj,
          features = genes,
          group.by = rv$visualization_sample_column,
          scale = input$dotplot_scale
        )$data
        
        list(
          pct_min = min(dot_data$pct.exp, na.rm = TRUE),
          pct_max = max(dot_data$pct.exp, na.rm = TRUE),
          avg_min = min(dot_data$avg.exp.scaled, na.rm = TRUE),
          avg_max = max(dot_data$avg.exp.scaled, na.rm = TRUE)
        )
      })
      
      global_pct_min <- min(sapply(expr_stats, function(x) x$pct_min))
      global_pct_max <- max(sapply(expr_stats, function(x) x$pct_max))
      global_avg_min <- min(sapply(expr_stats, function(x) x$avg_min))
      global_avg_max <- max(sapply(expr_stats, function(x) x$avg_max))
      
      # Create dotplot for each split
      plot_list <- lapply(seq_along(data_prep$split_data), function(i) {
        subset_obj <- data_prep$split_data[[i]]
        split_val <- data_prep$split_vals[i]
        
        p <- suppressWarnings(Seurat::DotPlot(
          object = subset_obj,
          features = genes,
          group.by = rv$visualization_sample_column,
          dot.scale = input$dotplot_dot_max_size,
          scale = input$dotplot_scale
        ))
        
        # Apply global scale limits
        p <- p +
          ggplot2::scale_color_gradient(
            low = "lightgrey", 
            high = "blue",
            limits = c(global_avg_min, global_avg_max)
          ) +
          ggplot2::scale_size_continuous(
            limits = c(global_pct_min, global_pct_max)
          ) +
          labs(title = split_val) +
          theme(
            axis.text.x = element_text(
              angle = input$dotplot_x_rotate,
              hjust = 1,
              size = input$dotplot_x_text_size
            ),
            axis.text.y = element_text(
              size = input$dotplot_y_text_size
            ),
            legend.position = "none",
            plot.title = element_text(hjust = 0.5, face = "bold")
          ) +
          labs(
            x = NULL,
            y = NULL
          )
        
        p
      })
      
      # Combine plots with shared legend
      combined_plot <- patchwork::wrap_plots(plot_list, ncol = length(data_prep$split_vals)) &
        theme(
          legend.position = input$dotplot_legend_position,
          legend.text = element_text(size = input$dotplot_legend_text_size),
          legend.title = element_text(size = input$dotplot_legend_title_size),
          legend.key.size = unit(input$dotplot_legend_key_size, "cm")
        ) &
        patchwork::plot_layout(guides = "collect")
      
      combined_plot
      
    } else {
      suppressWarnings(Seurat::DotPlot(
        object = data_prep$data,
        features = genes,
        group.by = rv$visualization_sample_column,
        dot.scale = input$dotplot_dot_max_size,
        scale = input$dotplot_scale
      ) +
        theme(
          axis.text.x = element_text(
            angle = input$dotplot_x_rotate,
            hjust = 1,
            size = input$dotplot_x_text_size
          ),
          axis.text.y = element_text(
            size = input$dotplot_y_text_size
          ),
          legend.position = input$dotplot_legend_position,
          legend.text = element_text(size = input$dotplot_legend_text_size),
          legend.title = element_text(size = input$dotplot_legend_title_size),
          legend.key.size = unit(input$dotplot_legend_key_size, "cm")
        ) +
        labs(
          x = NULL,
          y = NULL,
          title = NULL
        ))
    }
  })
  
  output$dotplot_plot_ui <- renderUI({
    
    if (!is.null(rv$visualization_split_column) && rv$visualization_split_column != "") {
      n_splits <- length(unique(rv$data_obj@meta.data[[rv$visualization_split_column]]))
      width_val <- paste0(n_splits * 480, "px")  
    } else {
      width_val <- "480px"
    }
    
    div(
      style = "max-width: 800px; overflow-x: auto; width: 100%;",
      plotOutput("dotplot_plot", height = "288px", width = width_val)
    )
  })
  
  output$dotplot_plot <- renderPlot({
    req(dotplot_plot_reactive())
    req(!rv$ui_freeze)
    dotplot_plot_reactive()
  })
  
  observeEvent(input$save_dotplot, {
    
    req(dotplot_plot_reactive())
    req(input$dotplot_save_path)
    
    data_prep <- visualization_data_prepared()
    
    # Check if split is enabled
    has_split <- data_prep$has_split
    
    if (has_split) {
      n_splits <- length(data_prep$split_vals)
      pdf_width <- input$dotplot_width * n_splits
    } else {
      pdf_width <- input$dotplot_width
    }
    
    out_dir <- input$dotplot_save_path
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    
    out_file <- file.path(out_dir,
                          paste0("dotplot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"))
    
    pdf(out_file,
        width = pdf_width,
        height = input$dotplot_height)
    
    print(dotplot_plot_reactive())
    
    dev.off()
    
    showNotification("Saved dotplot", type = "message")
  })
  
  output$download_dotplot <- downloadHandler(
    
    filename = function() {
      paste0("dotplot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
    },
    
    content = function(file) {
      
      req(dotplot_plot_reactive())
      
      data_prep <- visualization_data_prepared()
      
      # Check if split is enabled
      has_split <- data_prep$has_split
      
      if (has_split) {
        n_splits <- length(data_prep$split_vals)
        pdf_width <- input$dotplot_width * n_splits
      } else {
        pdf_width <- input$dotplot_width
      }
      
      pdf(file,
          width = pdf_width,
          height = input$dotplot_height)
      
      print(dotplot_plot_reactive())
      
      dev.off()
    }
  )
  # ----------------- VISUALIZATION FUNCTION Feature Plot -----------------
  create_featureplot <- function(
    seurat_obj,
    gene,
    split.by = NULL,
    assay = "RNA",
    reduction = "umap",
    pt.size = 0.5,
    alpha = 1,
    legend.position = "right",
    title.size = 12,
    legend.text.size = 10,
    axis.title.size = 12,
    axis.text.size = 10
  ) {
    library(patchwork)
    
    all_genes <- rownames(seurat_obj)
    
    if (!gene %in% all_genes) {
      stop(paste0("Gene '", gene, "' not found in Seurat object."))
    }
    
    if (!is.null(assay)) {
      DefaultAssay(seurat_obj) <- assay
    }
    
    p <- FeaturePlot(
      object = seurat_obj,
      features = gene,
      split.by = split.by,
      reduction = reduction,
      pt.size = pt.size,
      order = TRUE,
      combine = TRUE
    )
    
    if (inherits(p, "patchwork")) {
      for (i in seq_along(p$patches$plots)) {
        p$patches$plots[[i]]$layers[[1]]$aes_params$alpha <- alpha
      }
    } else {
      p$layers[[1]]$aes_params$alpha <- alpha
    }
    
    if (!is.null(split.by)) {
      p <- p & theme(
        axis.title.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.line.y.right = element_blank()
      )
    }
    
    p <- p +
      plot_layout(guides = "collect") +
      plot_annotation(
        title = gene,
        theme = theme(
          plot.title = element_text(
            size = title.size, 
            hjust = 0.5, 
            face = "bold",
            margin = margin(b = 10) 
          )
        )
      ) &
      theme(
        axis.title = element_text(size = axis.title.size),
        axis.text = element_text(size = axis.text.size),
        
        legend.position = legend.position,
        legend.title = element_text(size = title.size),
        legend.text = element_text(size = legend.text.size)
      )
    
    return(p)
  }
  
  output$feature_plot <- renderPlot({
    req(rv$data_obj)
    req(rv$visualization_marker)
    req(rv$visualization_sample_column)
    req(!rv$ui_freeze)
    
    genes <- unique(trimws(unlist(strsplit(rv$visualization_marker, "\n"))))
    genes <- genes[genes != ""]
    gene <- genes[genes %in% rownames(rv$data_obj)][1]
    create_featureplot(seurat_obj = rv$data_obj,
                       gene = gene,
                       split.by = rv$visualization_sample_column,
                       pt.size = input$feature_plot_point_size,
                       legend.position = input$feature_plot_legend_position,
                       title.size = input$feature_plot_title_size,
                       legend.text.size = input$feature_plot_legend_text_size,
                       axis.title.size = input$feature_plot_axis_title_size,
                       axis.text.size = input$feature_plot_text_size)
  })
  
  # Helper function to generate feature plots directly as PNG
  generate_feature_plot_png <- function(seurat_obj, gene, output_file,
                                        split.by = NULL,
                                        pt.size = 1,
                                        legend.position = "right",
                                        title.size = 12,
                                        legend.text.size = 10,
                                        axis.title.size = 10,
                                        axis.text.size = 8,
                                        dpi = 300,
                                        width = 7,
                                        height = 7) {
    
    library(png)
    
    # Create PNG directly (much smaller than PDF)
    png(output_file, width = width * dpi, height = height * dpi, res = dpi)
    
    print(
      create_featureplot(seurat_obj = seurat_obj,
                         gene = gene,
                         split.by = split.by,
                         pt.size = pt.size,
                         legend.position = legend.position,
                         title.size = title.size,
                         legend.text.size = legend.text.size,
                         axis.title.size = axis.title.size,
                         axis.text.size = axis.text.size)
    )
    
    dev.off()
  }
  
  # Helper function to combine PNGs into multi-page PDF
  combine_pngs_to_pdf_feature <- function(png_files, output_path, 
                                          genes_per_page = 4,
                                          final_width = 8.5,
                                          final_height = 11,
                                          update_progress_fn = NULL) {
    
    library(magick)
    
    if (length(png_files) == 0) {
      stop("No PNG files to combine")
    }
    
    if (!is.null(update_progress_fn)) {
      update_progress_fn("<b>Combining images into PDF...</b><br>Reading PNG files...")
    }
    
    num_pages <- ceiling(length(png_files) / genes_per_page)
    
    tryCatch({
      # Read all PNG files
      all_images <- list()
      for (i in seq_along(png_files)) {
        if (!is.null(update_progress_fn) && i %% 5 == 0) {
          update_progress_fn(paste0(
            "<b>Combining images into PDF...</b><br>",
            "Loading images: ", i, " / ", length(png_files)
          ))
        }
        
        tryCatch({
          img <- magick::image_read(png_files[i])
          all_images[[i]] <- img
        }, error = function(e) {
          warning(paste("Failed to read", png_files[i]))
        })
      }
      
      # Create multi-page PDF
      if (!is.null(update_progress_fn)) {
        update_progress_fn("<b>Combining images into PDF...</b><br>Creating pages...")
      }
      
      page_images <- list()
      
      for (page in 1:num_pages) {
        start_idx <- (page - 1) * genes_per_page + 1
        end_idx <- min(page * genes_per_page, length(all_images))
        
        # Get images for this page
        page_imgs <- all_images[start_idx:end_idx]
        page_imgs <- page_imgs[!sapply(page_imgs, is.null)]
        
        if (length(page_imgs) > 0) {
          # Stack images vertically
          combined <- page_imgs[[1]]
          
          if (length(page_imgs) > 1) {
            for (j in 2:length(page_imgs)) {
              combined <- magick::image_append(c(combined, page_imgs[[j]]), 
                                               stack = TRUE)
            }
          }
          
          # Add white space at bottom if less than 4 images
          if (length(page_imgs) < genes_per_page) {
            img_height <- magick::image_info(combined)$height
            padding_height <- (genes_per_page - length(page_imgs)) * 
              (img_height / length(page_imgs))
            
            white_space <- magick::image_blank(width = magick::image_info(combined)$width,
                                               height = padding_height, 
                                               color = "white")
            combined <- magick::image_append(c(combined, white_space), 
                                             stack = TRUE)
          }
          
          # Resize to fit page
          info <- magick::image_info(combined)
          combined <- magick::image_scale(combined, 
                                          geometry = paste0(
                                            as.integer(final_width * 300), "x",
                                            as.integer(final_height * 300)
                                          ))
          
          page_images[[page]] <- combined
        }
      }
      
      # Combine all pages into single PDF
      if (!is.null(update_progress_fn)) {
        update_progress_fn("<b>Combining images into PDF...</b><br>Writing PDF...")
      }
      
      final_images <- page_images[[1]]
      if (length(page_images) > 1) {
        for (p in 2:length(page_images)) {
          final_images <- c(final_images, page_images[[p]])
        }
      }
      
      magick::image_write(final_images, output_path, format = "pdf")
      
      if (!is.null(update_progress_fn)) {
        update_progress_fn("<b>PDF created successfully!</b>")
      }
      
      return(output_path)
      
    }, error = function(e) {
      stop(paste("Error combining PNGs:", e$message))
    })
  }
  
  # Modified observeEvent
  observeEvent(input$generate_all_feature_plot, {
    
    req(rv$data_obj)
    
    showModal(modalDialog(
      title = "Processing Genes",
      
      tags$div(
        id = "feature_plot_progress_detail",
        "Starting..."
      ),
      
      easyClose = FALSE,
      footer = NULL
    ))
    
    markers <- unlist(strsplit(rv$visualization_marker, "\n"))
    markers <- trimws(markers)
    markers <- markers[markers != "" & !is.na(markers)]
    valid_genes <- rownames(rv$data_obj)
    
    markers <- markers[markers %in% valid_genes]
    if (length(markers) == 0) {
      removeModal()
      showNotification("No valid genes found in Seurat object", type = "error")
      return()
    }
    
    total <- length(markers)
    
    out_dir <- input$feature_plot_save_path
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    generated_files <- c()
    png_dir <- file.path(out_dir, ".temp_png")
    if (!dir.exists(png_dir)) dir.create(png_dir, recursive = TRUE)
    png_files <- c()
    
    # Step 1: Generate feature plots as PNG files
    for (i in seq_along(markers)) {
      genes <- markers[i]
      done <- i - 1
      remaining <- total - i + 1
      
      shinyjs::html(
        "feature_plot_progress_detail",
        paste0(
          "<b>Generating feature plots:</b> ", genes, "<br>",
          "Completed: ", done, "<br>",
          "Remaining: ", remaining
        )
      )
      
      Sys.sleep(0.01)
      
      # Save as PNG directly (much smaller than PDF)
      png_file <- file.path(png_dir, paste0(genes, ".png"))
      pdf_file <- file.path(out_dir, paste0(genes, ".pdf"))
      
      tryCatch({
        # Generate PNG
        generate_feature_plot_png(
          seurat_obj = rv$data_obj,
          gene = genes,
          output_file = png_file,
          split.by = rv$visualization_sample_column,
          pt.size = input$feature_plot_point_size,
          legend.position = input$feature_plot_legend_position,
          title.size = input$feature_plot_title_size,
          legend.text.size = input$feature_plot_legend_text_size,
          axis.title.size = input$feature_plot_axis_title_size,
          axis.text.size = input$feature_plot_text_size,
          dpi = 100,  # Adjust DPI for quality/size tradeoff
          width = input$feature_plot_width,
          height = input$feature_plot_height
        )
        
        # Also create PDF for individual access
        pdf(pdf_file,
            width = input$feature_plot_width,
            height = input$feature_plot_height)
        
        print(
          create_featureplot(seurat_obj = rv$data_obj,
                             gene = genes,
                             split.by = rv$visualization_sample_column,
                             pt.size = input$feature_plot_point_size,
                             legend.position = input$feature_plot_legend_position,
                             title.size = input$feature_plot_title_size,
                             legend.text.size = input$feature_plot_legend_text_size,
                             axis.title.size = input$feature_plot_axis_title_size,
                             axis.text.size = input$feature_plot_text_size)
        )
        
        dev.off()
        
        png_files <- c(png_files, png_file)
        generated_files <- c(generated_files, pdf_file)
        
      }, error = function(e) {
        warning(paste("Failed to generate plot for", genes, ":", e$message))
      })
    }
    
    rv$featureplot_files <- generated_files
    
    # Step 2: Combine PNG files into single multi-page PDF
    tryCatch({
      shinyjs::html(
        "feature_plot_progress_detail",
        "<b>Combined file</b><br>Creating multi-page PDF..."
      )
      
      combined_output <- file.path(out_dir, paste0("combined_feature_plots_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"))
      
      combine_pngs_to_pdf_feature(
        png_files = png_files,
        output_path = combined_output,
        genes_per_page = 4,
        final_width = input$feature_plot_width,
        final_height = input$feature_plot_height * 4 + 1,  # 4 plots per page
        update_progress_fn = function(msg) {
          shinyjs::html("feature_plot_progress_detail", paste0("<b>Combined file</b><br>", msg))
        }
      )
      
      rv$featureplot_files <- c(rv$featureplot_files, combined_output)
      
      # Step 3: Clean up temporary PNG files
      shinyjs::html(
        "feature_plot_progress_detail",
        "<b>Combined file</b><br>Cleaning up temporary files..."
      )
      
      unlink(png_dir, recursive = TRUE)
      
      removeModal()
      
      showNotification(
        paste("All feature plots generated and combined!",
              "\nIndividual PDFs: ", out_dir,
              "\nCombined PDF: combined_feature_plots.pdf"),
        type = "message",
        duration = 10
      )
      
    }, error = function(e) {
      removeModal()
      showNotification(
        paste("Error combining PDFs:", e$message),
        type = "error"
      )
    })
  })
  
  output$download_feature_plot <- downloadHandler(
    
    filename = function() {
      "feature_plots.zip"
    },
    
    content = function(file) {
      
      files <- rv$featureplot_files
      
      if (is.null(files) || length(files) == 0) {
        
        showNotification(
          "No feature plots generated yet.",
          type = "error",
          duration = 5
        )
        
        return(NULL)
      }
      
      old <- getwd()
      on.exit(setwd(old), add = TRUE)
      
      setwd(dirname(files[1]))
      
      zip::zip(
        zipfile = file,
        files = basename(files)
      )
      
      showNotification(
        "Feature plots downloaded successfully.",
        type = "message",
        duration = 3
      )
    }
  )
  # ----------------- VISUALIZATION FUNCTION Original Volcano Plot -----------------
  
  volcano_plot_reactive <- reactive({
    
    req(rv$deg_obj)
    req(!rv$ui_freeze)
    
    volcanoplot <- function(de, labeled_genes, ylimit = 150, xlimit = 5,  
                            pvalue_limit = 0.05, log2FC_limit = 0.25, 
                            size_point = 0.8, size_labeled_gene = 3,
                            text_size = 10, title_size = 10) {
      
      de$p_val_adj <- pmax(de$pvalue_1, 1e-300)  
      de <- de[!(de$intensity_1 == 0 & de$intensity_2 == 0), ]
      de$neg_log_pval <- -log10(de$p_val_adj)
      
      de$diffexpressed <- "NO"
      
      de$diffexpressed[de$log2FC < -log2FC_limit & de$pvalue_1 < pvalue_limit] <- "DOWN"
      de$diffexpressed[de$log2FC > log2FC_limit & de$pvalue_1 < pvalue_limit] <- "UP"
      
      de$gene <- trimws(de$gene) 
      labeled_genes <- trimws(labeled_genes)  
      
      de$delabel <- if (input$original_volcano_auto_label) {
        ifelse(de$diffexpressed != "NO", de$gene, NA)
      } else {
        ifelse(de$gene %in% labeled_genes & de$diffexpressed != "NO",
               de$gene, NA)
      }
      ggplot(data = de, aes(x = log2FC, y = neg_log_pval, col = diffexpressed, label = delabel)) +
        geom_vline(xintercept = c(-log2FC_limit, log2FC_limit), col = "gray", linetype = 'dashed') +
        geom_hline(yintercept = -log10(pvalue_limit), col = "gray", linetype = 'dashed') +
        geom_point(data = de %>% dplyr::filter(diffexpressed == "DOWN"), color = "#00AFBB", size = size_point, alpha = 0.5) +
        geom_point(data = de %>% dplyr::filter(diffexpressed == "UP"), color = "#bb0c00", size = size_point, alpha = 0.5) +
        geom_point(data = de %>% dplyr::filter(diffexpressed == "NO"), color = "gray", size = size_point, alpha = 0.5) +
        geom_text_repel(show.legend  = FALSE, 
                        box.padding = 0.3, 
                        min.segment.length = 0.1, 
                        direction = "both", 
                        segment.color = 'darkgray', 
                        max.overlaps = 20, 
                        size = size_labeled_gene, 
                        fontface = "bold.italic",
                        na.rm = TRUE) +
        scale_color_manual(values = c("DOWN" = "black", "UP" = "black")) +
        coord_cartesian(xlim = c(-xlimit, xlimit), ylim = c(0, ylimit)) +
        labs(x = expression("CONTROL <- log"[2]*"FC -> TREATMENT"), 
             y = expression("-log"[10]*"p-value")) +
        theme_bw() +
        theme(
          axis.text = element_text(size = text_size),
          axis.title = element_text(size = title_size, face = "bold"),
          legend.position = "none"
        )
    }
    
    labeled_genes <- unique(trimws(unlist(strsplit(rv$visualization_marker, "\n"))))
    labeled_genes <- labeled_genes[labeled_genes != ""]
    
    volcanoplot(
      de = rv$deg_obj, 
      labeled_genes = labeled_genes, 
      ylimit = input$original_volcano_ylimit,
      xlimit = input$original_volcano_xlimit,
      pvalue_limit = input$original_volcano_pvalue_limit,
      log2FC_limit = input$original_volcano_log2FC_limit,
      size_point = input$original_volcano_size_point,
      size_labeled_gene = input$original_volcano_size_label,
      text_size = input$original_volcano_text_size,
      title_size = input$original_volcano_title_size
    )
  })
  
  output$volcano_plot <- renderPlot({
    req(volcano_plot_reactive())
    req(!rv$ui_freeze)
    volcano_plot_reactive()
  })
  
  observeEvent(input$save_original_volcano_plot, {
    
    req(volcano_plot_reactive())
    req(input$original_volcano_save_path)
    
    out_dir <- input$original_volcano_save_path
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    
    out_file <- file.path(out_dir,
                          paste0("volcano(gene)_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"))
    
    pdf(out_file,
        width = input$original_volcano_width,
        height = input$original_volcano_height)
    
    print(volcano_plot_reactive())
    
    dev.off()
    
    showNotification("Saved volcano(gene) plot", type = "message")
  })
  
  output$download_original_volcano_plot <- downloadHandler(
    
    filename = function() {
      paste0("volcano(gene)_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
    },
    
    content = function(file) {
      
      req(volcano_plot_reactive())
      
      pdf(file,
          width = input$original_volcano_width,
          height = input$original_volcano_height)
      
      print(volcano_plot_reactive())
      
      dev.off()
    }
  )
  # ----------------- VISUALIZATION FUNCTION Pathway Volcano Plot (Single) -----------------

  pathway_volcano_plot_reactive <- reactive({

    req(rv$deg_obj)
    req(rv$gsea_obj)
    req(rv$visualization_pathway != "")
    req(rv$pathway_rename)
    req(!rv$ui_freeze)

    selected_pathway <- strsplit(
      rv$visualization_pathway,
      "\n"
    )[[1]]

    selected_pathway <- selected_pathway[
      selected_pathway != ""
    ]
    if (length(selected_pathway) > 30) {

      return(
        list(
          type = "message",
          content = paste0(
            "Error: Too many pathways selected (>30). ",
            "Please reduce selection."
          )
        )
      )
    }
    
    if(sum(rv$gsea_obj$ID %in% selected_pathway) == 0){
      return(
        list(
          type = "message",
          content = paste0(
            "Error: No matched pathway. ",
            "Please check the pathway name input."
          )
        )
      )
    }

    selected_pathway <- intersect(
      selected_pathway,
      rv$gsea_obj$ID
    )

    draw_pathway_volcano <- function(
    deg,
    gsea,
    selected_pathway,
    pathway_rename,
    pvalue_limit = 0.05,
    log2FC_limit = 0.25,
    top_n_label = 10,
    point_size = 2.5,
    label_size = 3,
    axis_title_size = 14,
    axis_text_size = 12,
    legend_position = "right",
    legend_text_size = 11
    ) {
      
      library(dplyr)
      library(tidyr)
      library(ggplot2)
      library(ggrepel)
      library(scales)
      library(scatterpie)
      
      deg <- deg %>%
        mutate(
          p_val_adj = pmax(pvalue_1, 1e-300),
          neg_log_pval = -log10(p_val_adj)
        )
      
      marker_df <- gsea %>%
        filter(ID %in% selected_pathway) %>%
        select(ID, core_enrichment) %>%
        mutate(gene = strsplit(core_enrichment, "/")) %>%
        unnest(gene) %>%
        rename(pathway = ID)
      
      rename_map <- setNames(
        pathway_rename$renamed_pathway,
        pathway_rename$original_pathway
      )
      
      marker_df <- marker_df %>%
        mutate(pathway = unname(rename_map[pathway]))
      
      pathway_cols <- unique(marker_df$pathway)
      
      pie_df <- marker_df %>%
        distinct(gene, pathway) %>%
        mutate(value = 1) %>%
        group_by(gene) %>%
        mutate(value = value / n()) %>%
        ungroup() %>%
        pivot_wider(
          names_from = pathway,
          values_from = value,
          values_fill = 0
        )
      
      plot_df <- pie_df %>%
        inner_join(deg, by = "gene") %>%
        mutate(
          n_pathway = rowSums(across(all_of(pathway_cols)) > 0),
          significant = (p_val_adj < pvalue_limit &
                           abs(log2FC) > log2FC_limit)
        )
      
      label_df <- plot_df %>%
        filter(significant) %>%
        arrange(desc(neg_log_pval)) %>%
        slice_head(n = top_n_label)
      
      label_df$label_color <- "black"
      
      pathway_colors <- setNames(
        hue_pal()(length(pathway_cols)),
        pathway_cols
      )
      
      xlim <- max(abs(plot_df$log2FC), na.rm = TRUE) + 2
      ylim <- max(plot_df$neg_log_pval, na.rm = TRUE) + 1
      
      p <- ggplot() +
        
      scatterpie::geom_scatterpie(
        data = plot_df,
        aes(x = log2FC, y = neg_log_pval),
        cols = pathway_cols,
        pie_scale = point_size,
        color = NA,
        alpha = 0.9
      ) +
        
      geom_point(
        data = plot_df,
        aes(x = log2FC, y = neg_log_pval),
        alpha = 0,   # invisible but REAL geometry
        size = 0.5
      ) +
        
      geom_text_repel(
        data = label_df,
        aes(x = log2FC, y = neg_log_pval, label = gene),
        size = label_size,
        fontface = "bold.italic",
        segment.color = "gray40",
        segment.size = 0.4,
        min.segment.length = 0,   
        box.padding = 1,
        point.padding = 0,
        force = 2,               
        force_pull = 1,
        max.overlaps = Inf,
        seed = 123,
        show.legend = FALSE
      ) +
        
        # thresholds
        geom_hline(yintercept = -log10(pvalue_limit),
                   linetype = "dashed",
                   color = "gray") +
        geom_vline(xintercept = c(-log2FC_limit, log2FC_limit),
                   linetype = "dashed",
                   color = "gray") +
        
        coord_cartesian(      
          xlim = c(-xlim, xlim),
          ylim = c(0, ylim)) +
        
        scale_fill_manual(values = pathway_colors) +
        
        labs(
          x = expression("log"[2] * "FC"),
          y = expression("-log"[10] * "p-value"),
          fill = "Pathway"
        ) +
        
        theme_minimal(base_size = 14) +
        theme(
          axis.title = element_text(size = axis_title_size),
          axis.text = element_text(size = axis_text_size),
          legend.position = legend_position,
          legend.text = element_text(size = legend_text_size),
          legend.title = element_text(size = legend_text_size + 1)
        )
      
      return(p)
    }

    draw_pathway_volcano(
      deg = rv$deg_obj,
      gsea = rv$gsea_obj,
      selected_pathway = selected_pathway,
      pathway_rename = rv$pathway_rename,
      pvalue_limit = input$pathway_volcano_pv_limit,
      log2FC_limit = input$pathway_volcano_log2fc_limit,
      top_n_label = input$pathway_volcano_label_n,
      point_size = input$pathway_volcano_point_size,
      label_size = input$pathway_volcano_label_size,
      axis_title_size = input$pathway_volcano_axis_title_size,
      axis_text_size = input$pathway_volcano_axis_text_size,
      legend_position = input$pathway_volcano_legend_position,
      legend_text_size = input$pathway_volcano_legend_text_size
    )

  })

  output$pathway_volcano_plot <- renderPlot({
    
    req(pathway_volcano_plot_reactive())
    req(!rv$ui_freeze)
    res <- pathway_volcano_plot_reactive()
    
    if (is.list(res) && res$type == "message") {
      
      plot.new()
      text(0.5, 0.5, res$content, cex = 1.2, col = "red")
      
      return()
    }
    res
  })
  
  observeEvent(input$save_pathway_volcano, {
    
    req(pathway_volcano_plot_reactive())
    req(input$pathway_volcano_save_path)
    
    out_dir <- input$pathway_volcano_save_path
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    
    out_file <- file.path(
      out_dir,
      paste0("volcano(pathway)_multi_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
    )
    
    pdf(out_file,
        width = input$pathway_volcano_pdf_width,
        height = input$pathway_volcano_pdf_height)
    
    print(pathway_volcano_plot_reactive())
    
    dev.off()
    
    showNotification("Saved volcano(pathway) plot", type = "message")
  })
  
  output$download_pathway_volcano <- downloadHandler(
    
    filename = function() {
      paste0("volcano(pathway)_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
    },
    
    content = function(file) {
      
      req(pathway_volcano_plot_reactive())
      
      pdf(file,
          width = input$pathway_volcano_pdf_width,
          height = input$pathway_volcano_pdf_height)
      
      print(pathway_volcano_plot_reactive())
      
      dev.off()
    }
  )
  # ----------------- VISUALIZATION FUNCTION Pathway Volcano Plot (Multi) -----------------
  draw_pathway_volcano_single <- function(deg,
                                          gsea,
                                          selected_pathway,
                                          pathway_rename,
                                          pathway_color = "#E64B35",
                                          pvalue_limit = 0.05,
                                          log2FC_limit = 0.25,
                                          top_n_label = 10,
                                          point_size = 2.5,
                                          label_size = 3,
                                          axis_title_size = 14,
                                          axis_text_size = 12,
                                          legend_position = "none",
                                          legend_text_size = 11,
                                          xlim = NULL,
                                          ylim = NULL) {
    
    deg <- deg %>%
      mutate(
        p_val_adj = pmax(pvalue_1, 1e-300),
        neg_log_pval = -log10(p_val_adj)
      )
    
    pathway_name <- pathway_rename %>%
      filter(original_pathway == selected_pathway) %>%
      pull(renamed_pathway)
    
    if (length(pathway_name) == 0) {
      pathway_name <- selected_pathway
    }
    
    marker_df <- gsea %>%
      filter(ID == selected_pathway) %>%
      select(ID, NES, core_enrichment) %>%
      mutate(gene = strsplit(core_enrichment, "/")) %>%
      tidyr::unnest(gene) %>%
      rename(pathway = ID)
    
    plot_df <- marker_df %>%
      inner_join(deg, by = "gene") %>%
      mutate(
        significant = (p_val_adj < pvalue_limit) &
          (abs(log2FC) > log2FC_limit)
      )
    
    label_candidate_df <- plot_df %>%
      filter(
        p_val_adj < pvalue_limit,
        abs(log2FC) > log2FC_limit
      )
    
    top_label_df <- label_candidate_df %>%
      arrange(desc(neg_log_pval)) %>%
      slice_head(n = top_n_label)
    
    # =========================
    # AUTO LIMITS IF NOT GIVEN
    # =========================
    
    if (is.null(xlim)) {
      xlim <- max(abs(deg$log2FC), na.rm = TRUE) + 0.5
    }
    
    if (is.null(ylim)) {
      ylim <- max(deg$neg_log_pval, na.rm = TRUE) + 1
    }
    
    p <- ggplot(
      plot_df,
      aes(x = log2FC, y = neg_log_pval)
    ) +
      
      geom_point(
        aes(alpha = significant),
        color = pathway_color,
        size = point_size
      ) +
      
      scale_alpha_manual(
        values = c("TRUE" = 1, "FALSE" = 0.2),
        guide = "none"
      ) +
      
      coord_cartesian(
        xlim = c(-xlim, xlim),
        ylim = c(0, ylim)
      ) +
      
      geom_hline(
        yintercept = -log10(pvalue_limit),
        color = "gray",
        linetype = "dashed"
      ) +
      
      geom_vline(
        xintercept = c(-log2FC_limit, log2FC_limit),
        color = "gray",
        linetype = "dashed"
      ) +
      
      geom_text_repel(
        data = top_label_df,
        aes(label = gene),
        show.legend = FALSE,
        box.padding = 0.3,
        min.segment.length = 0.1,
        direction = "both",
        segment.color = "darkgray",
        max.overlaps = 20,
        size = label_size,
        na.rm = TRUE,
        fontface = "bold.italic"
      ) +
      
      labs(
        title = pathway_name,
        x = expression("CONTROL <- log"[2] * "FC -> TREATMENT"),
        y = expression("-log"[10] * "p-value")
      ) +
      
      theme_minimal(base_size = 14) +
      
      theme(
        axis.title = element_text(size = axis_title_size),
        axis.text = element_text(size = axis_text_size),
        legend.position = legend_position,
        plot.title = element_text(
          hjust = 0.5,
          face = "bold"
        )
      )
    
    return(p)
  }
  
  output$pathway_volcano_multi_plot <- renderPlot({
    req(rv$deg_obj)
    req(rv$gsea_obj)
    req(rv$visualization_pathway!="")
    req(rv$pathway_rename)
    req(!rv$ui_freeze)

    selected_pathway <- strsplit(rv$visualization_pathway, "\n")[[1]]
    selected_pathway <- selected_pathway[selected_pathway != ""]

    selected_pathway <- intersect(selected_pathway, rv$gsea_obj$ID)
    selected_pathway <- selected_pathway[1]
    draw_pathway_volcano_single(
      deg = rv$deg_obj,
      gsea = rv$gsea_obj,
      selected_pathway = selected_pathway,
      pathway_rename = rv$pathway_rename,
      pvalue_limit = input$pathway_volcano_pv_limit,
      log2FC_limit = input$pathway_volcano_log2fc_limit,
      top_n_label = input$pathway_volcano_label_n,
      point_size = input$pathway_volcano_point_size,
      label_size = input$pathway_volcano_label_size,
      axis_title_size = input$pathway_volcano_axis_title_size,
      axis_text_size = input$pathway_volcano_axis_text_size,
      legend_position = input$pathway_volcano_legend_position,
      legend_text_size = input$pathway_volcano_legend_text_size
    )
  })
  
  observeEvent(input$save_pathway_multi_volcano, {
    
    req(rv$deg_obj)
    req(rv$gsea_obj)
    req(rv$visualization_pathway != "")
    req(rv$pathway_rename)
    req(!rv$ui_freeze)
    
    showModal(
      modalDialog(
        title = "Generating Pathway Volcano Plots",
        
        tags$div(
          id = "pathway_volcano_progress_detail",
          "Starting..."
        ),
        
        easyClose = FALSE,
        footer = NULL
      )
    )
    
    color_total <- c(
      "#e6194b","#ffe119","#46f0f0","#f58231","#bcf60c", 
      "#ff00ff","#9a6324","#fffac8","#e6beff","#00bfff", 
      "#ffd8b1","#00ff7f","#f5a9bc","#1e90ff","#ffa500",
      "#98fb98","#911eb4","#afeeee","#fa8072","#9acd32",
      "#3cb44b","#000075","#808000","#cd5c5c","#dda0dd",
      "#40e0d0","#ff69b4","#8a2be2","#c71585","#5f9ea0",
      "#dc143c","#87cefa","#ff6347","#9932cc","#00ced1",
      "#ff4500","#6a5acd","#b0e0e6","#d2691e","#a9a9f5",
      "#adff2f","#8b0000","#7fffd4","#00fa9a","#ba55d3",
      "#2e8b57","#ffdab9","#b22222","#ffe4e1","#7b68ee"
    )
    
    selected_pathway <- strsplit(rv$visualization_pathway, "\n")[[1]]
    selected_pathway <- trimws(selected_pathway)
    selected_pathway <- selected_pathway[selected_pathway != ""]
    
    selected_pathway <- intersect(
      selected_pathway,
      rv$gsea_obj$ID
    )
    
    selected_pathway <- strsplit(rv$visualization_pathway, "\n")[[1]]
    selected_pathway <- trimws(selected_pathway)
    selected_pathway <- selected_pathway[selected_pathway != ""]
    
    selected_pathway <- intersect(
      selected_pathway,
      rv$gsea_obj$ID
    )
    all_pathway_genes <- rv$gsea_obj %>%
      dplyr::filter(ID %in% selected_pathway) %>%
      dplyr::pull(core_enrichment) %>%
      strsplit("/") %>%
      unlist() %>%
      unique()
    
    global_deg_df <- rv$deg_obj %>%
      dplyr::filter(gene %in% all_pathway_genes) %>%
      dplyr::mutate(
        neg_log_pval = -log10(
          pmax(pvalue_1, 1e-300)
        )
      )
    global_xlim <- max(
      abs(global_deg_df$log2FC),
      na.rm = TRUE
    ) + 0.5
    
    global_ylim <- max(
      global_deg_df$neg_log_pval,
      na.rm = TRUE
    ) + 1
    
    if (length(selected_pathway) == 0) {
      
      removeModal()
      
      showNotification(
        "No valid pathways found in GSEA object",
        type = "error"
      )
      
      return()
    }
    
    total <- length(selected_pathway)
    
    out_dir <- input$pathway_volcano_save_path
    
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }
    
    generated_files <- c()
    
    for (i in seq_along(selected_pathway)) {
      
      pathway <- selected_pathway[i]
      
      done <- i - 1
      remaining <- total - i + 1
      
      shinyjs::html(
        "pathway_volcano_progress_detail",
        paste0(
          "<b>Processing pathway:</b> ", pathway, "<br>",
          "Completed: ", done, "<br>",
          "Remaining: ", remaining
        )
      )
      
      Sys.sleep(0.01)
      
      color_index <- ((i - 1) %% length(color_total)) + 1
      pathway_color <- color_total[color_index]
      
      pathway_title <- rv$pathway_rename %>%
        dplyr::filter(original_pathway == pathway) %>%
        dplyr::pull(renamed_pathway)
      
      if (length(pathway_title) == 0) {
        pathway_title <- pathway
      }
      
      safe_name <- gsub("[^A-Za-z0-9_\\-]", "_", pathway_title)
      
      file_path <- file.path(
        out_dir,
        paste0(safe_name, ".pdf")
      )

      p <- draw_pathway_volcano_single(
        deg = rv$deg_obj,
        gsea = rv$gsea_obj,
        selected_pathway = pathway,
        pathway_rename = rv$pathway_rename,
        pathway_color = pathway_color,
        pvalue_limit = input$pathway_volcano_pv_limit,
        log2FC_limit = input$pathway_volcano_log2fc_limit,
        top_n_label = input$pathway_volcano_label_n,
        point_size = input$pathway_volcano_point_size,
        label_size = input$pathway_volcano_label_size,
        axis_title_size = input$pathway_volcano_axis_title_size,
        axis_text_size = input$pathway_volcano_axis_text_size,
        legend_position = input$pathway_volcano_legend_position,
        legend_text_size = input$pathway_volcano_legend_text_size,
        xlim = global_xlim,
        ylim = global_ylim
      )
      pdf(
        file_path,
        width = input$pathway_volcano_pdf_width,
        height = input$pathway_volcano_pdf_height
      )
      
      print(p)
      
      dev.off()
      
      generated_files <- c(
        generated_files,
        file_path
      )
    }
    
    rv$pathway_multi_volcano_files <- generated_files
    
    removeModal()
    
    showNotification(
      "All pathway volcano plots generated!",
      type = "message"
    )
  })
  
  output$download_pathway_multi_volcano <- downloadHandler(
    
    filename = function() {
      "pathway_multi_volcano_plots.zip"
    },
    
    content = function(file) {
      
      files <- rv$pathway_multi_volcano_files
      
      if (is.null(files) || length(files) == 0) {
        
        showNotification(
          "No pathway volcano plots generated yet.",
          type = "error",
          duration = 5
        )
        
        return(NULL)
      }
      
      old <- getwd()
      on.exit(setwd(old), add = TRUE)
      
      setwd(dirname(files[1]))
      
      zip::zip(
        zipfile = file,
        files = basename(files)
      )
      
      showNotification(
        "Pathway volcano plots downloaded successfully.",
        type = "message",
        duration = 3
      )
    }
  )
  
  # ----------------- VISUALIZATION FUNCTION NES Bar Plot -----------------
  nes_barplot_reactive <- reactive({

    req(rv$gsea_obj)
    req(rv$visualization_pathway != "")
    req(rv$pathway_rename)
    req(!rv$ui_freeze)
    
    
    selected_pathway <- strsplit(
      rv$visualization_pathway,
      "\n"
    )[[1]]
    
    selected_pathway <- selected_pathway[
      selected_pathway != ""
    ]
    if (length(selected_pathway) > 40) {
      
      return(
        list(
          type = "message",
          content = paste0(
            "Error: Too many pathways selected (>40). ",
            "Please reduce selection."
          )
        )
      )
    }
    
    if(sum(rv$gsea_obj$ID %in% selected_pathway) == 0){
      return(
        list(
          type = "message",
          content = paste0(
            "Error: No matching pathways found in gsea_obj. ",
            "Please check the pathway name input."
          )
        )
      )
    }
    
    plot_gsea_nes_barplot <- function(gsea_obj,
                                      pathway_rename,
                                      wrap_width = 40,
                                      pathway_size = 12,
                                      nes_value_size = 12,
                                      y_title_size = 10,
                                      x_title_size = 10,
                                      bar_width = 0.7) {
      rename_map <- setNames(
        pathway_rename$renamed_pathway,
        pathway_rename$original_pathway
      )
      valid_ids <- pathway_rename$original_pathway
      gsea_obj <- gsea_obj[gsea_obj$ID %in% valid_ids, , drop = FALSE]
      gsea_obj$ID <- ifelse(
        gsea_obj$ID %in% names(rename_map),
        rename_map[gsea_obj$ID],
        gsea_obj$ID
      )
      gsea_obj$ID <- stringr::str_wrap(
        gsea_obj$ID,
        width = wrap_width
      )
      gsea_obj$ID <- reorder(
        gsea_obj$ID,
        gsea_obj$NES
      )
      p <- ggplot2::ggplot(gsea_obj, ggplot2::aes(x = ID, y = NES, fill = NES > 0)) +
        ggplot2::geom_bar(stat = "identity", width = bar_width) +
        ggplot2::coord_flip() +
        ggplot2::scale_fill_manual(values = c("lightblue", "darkred")) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(size = nes_value_size),
          axis.text.y = ggplot2::element_text(size = pathway_size, face = "bold"),
          axis.title.x = ggplot2::element_text(size = x_title_size, face = "bold"),
          axis.title.y = ggplot2::element_text(size = y_title_size, face = "bold")) +
        ggplot2::labs(x = "Pathway", y = "NES", fill = "Direction") +
        ggplot2::guides(fill = "none") +
        ggplot2::ylim(
          min(gsea_obj$NES, na.rm = TRUE) - 0.2,
          max(gsea_obj$NES, na.rm = TRUE) + 0.2
        )
      
      return(p)
    }
    plot_gsea_nes_barplot(gsea_obj = rv$gsea_obj, 
                          pathway_rename = rv$pathway_rename,
                          wrap_width = input$nes_barplot_pathway_length,
                          pathway_size = input$nes_barplot_pathway_size,
                          nes_value_size = input$nes_barplot_nes_value_size,
                          y_title_size = input$nes_barplot_y_title_size,
                          x_title_size = input$nes_barplot_x_title_size,
                          bar_width = input$nes_barplot_bar_width)
  })
  
  output$nes_barplot <- renderPlot({
    
    req(nes_barplot_reactive())
    req(!rv$ui_freeze)
    res <- nes_barplot_reactive()
    
    if (is.list(res) && res$type == "message") {
      
      plot.new()
      text(0.5, 0.5, res$content, cex = 1.2, col = "red")
      
      return()
    }
    res
  })
  
  observeEvent(input$save_nes_barplot, {
    
    req(nes_barplot_reactive())
    req(input$nes_barplot_save_path)
    
    out_dir <- input$nes_barplot_save_path
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    
    out_file <- file.path(
      out_dir,
      paste0("NES_barplot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
    )
    
    pdf(out_file,
        width = input$nes_barplot_pdf_width,
        height = input$nes_barplot_pdf_height)
    
    print(nes_barplot_reactive())
    
    dev.off()
    
    showNotification("Saved NES barplot", type = "message")
  })
  
  output$download_nes_barplot <- downloadHandler(
    
    filename = function() {
      paste0("NES_barplot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
    },
    
    content = function(file) {
      
      req(nes_barplot_reactive())
      
      pdf(file,
          width = input$nes_barplot_pdf_width,
          height = input$nes_barplot_pdf_height)
      
      print(nes_barplot_reactive())
      
      dev.off()
      showNotification(
        "NES barplots downloaded successfully.",
        type = "message",
        duration = 3
      )
    }
  )
  # ----------------- VISUALIZATION FUNCTION Cell Cell Communication netVisual_circle_plot-----------------
  netVisual_circle_reactive <- reactive({
    
    req(rv$cccd_obj)
    req(input$circle_measure)
    req(!rv$ui_freeze)
    
    cellchat_obj <- rv$cccd_obj$data
    cellchat_merge <- rv$cccd_obj$data_merge
    
    measure_title <- if (input$circle_measure == "count") {
      "Interaction Count"
    } else {
      "Interaction Weight"
    }
    
    plot_circle_one <- function(mat, title_name = NULL) {
      
      if (is.null(mat)) {
        plot.new()
        text(0.5, 0.5, "Selected network measure is not available.", col = "red", cex = 1.2)
        return(invisible(NULL))
      }
      
      if (nrow(mat) == 0 || ncol(mat) == 0) {
        plot.new()
        text(0.5, 0.5, "Network matrix is empty.", col = "red", cex = 1.2)
        return(invisible(NULL))
      }
      
      if (nrow(mat) == 0 || ncol(mat) == 0) {
        plot.new()
        text(0.5, 0.5, "No non-isolated cell groups remain.", col = "red", cex = 1.2)
        return(invisible(NULL))
      }
      
      final_title <- if (is.null(title_name)) {
        measure_title
      } else {
        paste0(title_name, " - ", measure_title)
      }
      
      group_size <- rep(1, nrow(mat))
      names(group_size) <- rownames(mat)
      
      netVisual_circle(
        net = mat,
        vertex.weight = group_size,
        weight.scale = TRUE,
        label.edge = isTRUE(input$circle_show_edge_label),
        edge.weight.max = max(mat, na.rm = TRUE),
        edge.width.max = input$circle_edge_width_max,
        vertex.size.max = input$circle_vertex_size,
        vertex.label.cex = input$circle_vertex_label_size,
        arrow.size = input$circle_arrow_size,
        title.name = final_title
      )
      
      invisible(NULL)
    }
    
    if (!is.null(cellchat_merge)) {
      
      dataset_names <- names(cellchat_merge@net)
      
      if (length(dataset_names) < 2) {
        return(list(
          type = "message",
          content = "Merged CellChat object should contain at least two datasets."
        ))
      }
      
      plot_fun <- function() {
        old_par <- par(no.readonly = TRUE)
        on.exit(par(old_par), add = TRUE)
        
        dataset_use <- dataset_names[1:min(2, length(dataset_names))]
        
        par(mfrow = c(1, length(dataset_use)), xpd = TRUE)
        
        for (dataset_name in dataset_use) {
          mat <- cellchat_merge@net[[dataset_name]][[input$circle_measure]]
          plot_circle_one(mat, title_name = dataset_name)
        }
        
        invisible(NULL)
      }
      
    } else {
      
      plot_fun <- function() {
        old_par <- par(no.readonly = TRUE)
        on.exit(par(old_par), add = TRUE)
        
        par(mfrow = c(1, 1), xpd = TRUE)
        
        mat <- cellchat_obj@net[[input$circle_measure]]
        plot_circle_one(mat, title_name = "CellChat Network")
        
        invisible(NULL)
      }
    }
    
    list(
      type = "plot",
      plot_fun = plot_fun
    )
  })
  
  output$netVisual_circle_control_ui <- renderUI({
    
    if (!is.null(rv$cccd_obj$data_merge)) {
      width_val <- "800px"
    } else {
      width_val <- "400px"
    }
    plotOutput("netVisual_circle_plot_ui", height = "400px", width = width_val)
  })
  
  output$netVisual_circle_plot_ui <- renderPlot({
    
    req(netVisual_circle_reactive())
    req(!rv$ui_freeze)
    
    res <- netVisual_circle_reactive()
    
    if (is.list(res) && res$type == "message") {
      plot.new()
      text(0.5, 0.5, res$content, cex = 1.2, col = "red")
      return()
    }
    
    res$plot_fun()
  })
  
  observeEvent(input$circle_barplot, {
    
    req(netVisual_circle_reactive())
    req(input$circle_save_path)
    
    res <- netVisual_circle_reactive()
    
    if (is.list(res) && res$type == "message") {
      showNotification(res$content, type = "error")
      return()
    }
    
    out_dir <- input$circle_save_path
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    
    out_file <- file.path(
      out_dir,
      paste0("CommunicationNetwork_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
    )
    
    pdf(
      out_file,
      width = input$circle_width,
      height = input$circle_height
    )
    
    res$plot_fun()
    
    dev.off()
    
    showNotification(
      paste0("Saved CellChat circle plot: ", out_file),
      type = "message"
    )
  })
  
  output$download_circle <- downloadHandler(
    
    filename = function() {
      paste0("CommunicationNetwork_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
    },
    
    content = function(file) {
      
      req(netVisual_circle_reactive())
      
      res <- netVisual_circle_reactive()
      
      pdf(
        file,
        width = input$circle_width,
        height = input$circle_height
      )
      
      if (is.list(res) && res$type == "message") {
        plot.new()
        text(0.5, 0.5, res$content, cex = 1.2, col = "red")
      } else {
        res$plot_fun()
      }
      
      dev.off()
      
      showNotification(
        "CellChat circle plot downloaded successfully.",
        type = "message",
        duration = 3
      )
    }
  )
  # ----------------- VISUALIZATION FUNCTION Cell Cell Communication netAnalysis_signalingRole_scatter-----------------
  netAnalysis_signalingRole_scatter_reactive <- reactive({
    
    req(rv$cccd_obj)
    req(!rv$ui_freeze)
    
    cellchat_obj <- rv$cccd_obj$data
    cellchat_merge <- rv$cccd_obj$data_merge
    
    make_role_scatter <- function(obj, title_name = NULL, weight.MinMax = NULL) {
      
      p <- netAnalysis_signalingRole_scatter(
        obj,
        title = title_name,
        weight.MinMax = weight.MinMax,
        label.size = input$role_scatter_label_size,
        dot.size = input$role_scatter_point_size
      )
      
      p <- p +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            size = input$role_scatter_title_size,
            face = "bold",
            hjust = 0.5
          ),
          axis.text = ggplot2::element_text(size = input$role_scatter_axis_text_size),
          axis.title = ggplot2::element_text(size = input$role_scatter_axis_title_size, face = "bold"),
          legend.text = ggplot2::element_text(size = input$role_scatter_legend_text_size),
          legend.title = ggplot2::element_text(size = input$role_scatter_legend_title_size),
          legend.position = "none"
        )
      
      return(p)
    }
    
    if (!is.null(cellchat_merge)) {
      
      object.list <- cellchat_obj
      
      if (!is.list(object.list) || length(object.list) < 2) {
        return(list(
          type = "message",
          content = "The uploaded file has a merged CellChat object, but the original object list is not available."
        ))
      }
      
      num.link <- sapply(object.list, function(x) {
        rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
      })
      
      weight.MinMax <- c(
        min(num.link, na.rm = TRUE),
        max(num.link, na.rm = TRUE)
      )
      
      gg <- list()
      
      for (i in seq_along(object.list)) {
        gg[[i]] <- make_role_scatter(
          obj = object.list[[i]],
          title_name = names(object.list)[i],
          weight.MinMax = weight.MinMax
        )
      }
      
      p <- patchwork::wrap_plots(plots = gg, nrow = 1)
      
    } else {
      
      p <- make_role_scatter(
        obj = cellchat_obj,
        title_name = "CellChat Signaling Roles",
        weight.MinMax = NULL
      )
    }
    
    return(p)
  })
  output$netAnalysis_signalingRole_scatter_control_ui <- renderUI({
    
    if (!is.null(rv$cccd_obj$data_merge)) {
      width_val <- "800px"
    } else {
      width_val <- "400px"
    }
      plotOutput("netAnalysis_signalingRole_scatter_plot", height = "300px", width = width_val)
  })
  
  output$netAnalysis_signalingRole_scatter_plot <- renderPlot({
    
    req(netAnalysis_signalingRole_scatter_reactive())
    req(!rv$ui_freeze)
    
    res <- netAnalysis_signalingRole_scatter_reactive()
    
    if (is.list(res) && !inherits(res, "ggplot") && res$type == "message") {
      plot.new()
      text(0.5, 0.5, res$content, cex = 1.2, col = "red")
      return()
    }
    
    print(res)
  })
  
  observeEvent(input$save_role_scatter, {
    
    req(netAnalysis_signalingRole_scatter_reactive())
    req(input$role_scatter_save_path)
    
    res <- netAnalysis_signalingRole_scatter_reactive()
    
    if (is.list(res) && !inherits(res, "ggplot") && res$type == "message") {
      showNotification(res$content, type = "error")
      return()
    }
    
    out_dir <- input$role_scatter_save_path
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    
    out_file <- file.path(
      out_dir,
      paste0("SendervsReceiverRoles", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
    )
    
    ggplot2::ggsave(
      filename = out_file,
      plot = res,
      width = input$role_scatter_width,
      height = input$role_scatter_height
    )
    
    showNotification(
      paste0("Saved signaling role scatter: ", out_file),
      type = "message"
    )
  })
  
  output$download_role_scatter <- downloadHandler(
    
    filename = function() {
      paste0("SendervsReceiverRoles", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
    },
    
    content = function(file) {
      
      req(netAnalysis_signalingRole_scatter_reactive())
      
      res <- netAnalysis_signalingRole_scatter_reactive()
      
      if (is.list(res) && !inherits(res, "ggplot") && res$type == "message") {
        pdf(file, width = input$role_scatter_width, height = input$role_scatter_height)
        plot.new()
        text(0.5, 0.5, res$content, cex = 1.2, col = "red")
        dev.off()
      } else {
        ggplot2::ggsave(
          filename = file,
          plot = res,
          width = input$role_scatter_width,
          height = input$role_scatter_height
        )
      }
      
      showNotification(
        "Signaling role scatter downloaded successfully.",
        type = "message",
        duration = 3
      )
    }
  )
  
  # ----------------- VISUALIZATION FUNCTION Cell Cell Communication netAnalysis_signalingRole_heatmap -----------------
  
  netAnalysis_signalingRole_heatmap_reactive <- reactive({
    
    req(rv$cccd_obj)
    req(input$netAnalysis_signalingRole_heatmap_pattern)
    req(input$netAnalysis_signalingRole_heatmap_color)
    req(input$netAnalysis_signalingRole_heatmap_axis_title_size)
    req(input$netAnalysis_signalingRole_heatmap_axis_text_size)
    req(input$netAnalysis_signalingRole_heatmap_legend_size)
    req(!rv$ui_freeze)
    
    cellchat_obj <- rv$cccd_obj$data
    cellchat_merge <- rv$cccd_obj$data_merge
    pattern_use <- input$netAnalysis_signalingRole_heatmap_pattern
    
    is_cellchat_obj <- function(x) {
      methods::is(x, "CellChat")
    }
    
    pattern_list <- if (pattern_use == "both") {
      c("incoming", "outgoing")
    } else {
      pattern_use
    }
    
    get_cell_levels <- function(obj) {
      if (!is.null(obj@idents)) {
        lev <- levels(obj@idents)
        if (!is.null(lev) && length(lev) > 0) {
          return(lev)
        }
        return(unique(as.character(obj@idents)))
      }
      
      if (!is.null(obj@netP$prob)) {
        return(dimnames(obj@netP$prob)[[1]])
      }
      
      NULL
    }
    
    get_scPalette <- function(n) {
      if (exists("scPalette", envir = asNamespace("CellChat"), inherits = FALSE)) {
        get("scPalette", envir = asNamespace("CellChat"))(n)
      } else if (exists("ggPalette", envir = asNamespace("CellChat"), inherits = FALSE)) {
        get("ggPalette", envir = asNamespace("CellChat"))(n)
      } else {
        grDevices::rainbow(n)
      }
    }
    
    get_role_matrix_raw <- function(obj, pattern_single, slot.name = "netP") {
      
      if (is.null(obj) || !methods::is(obj, "CellChat")) {
        return(NULL)
      }
      
      if (is.null(methods::slot(obj, slot.name)$centr) ||
          length(methods::slot(obj, slot.name)$centr) == 0) {
        obj <- netAnalysis_computeCentrality(obj, slot.name = slot.name)
      }
      
      centr <- methods::slot(obj, slot.name)$centr
      
      if (is.null(centr) || length(centr) == 0) {
        return(NULL)
      }
      
      cell_names <- get_cell_levels(obj)
      
      if (is.null(cell_names) || length(cell_names) == 0) {
        first_centr <- centr[[1]]
        if (!is.null(first_centr$outdeg)) {
          cell_names <- names(first_centr$outdeg)
        } else if (!is.null(first_centr$indeg)) {
          cell_names <- names(first_centr$indeg)
        }
      }
      
      if (is.null(cell_names) || length(cell_names) == 0) {
        return(NULL)
      }
      
      outgoing <- matrix(
        0,
        nrow = length(cell_names),
        ncol = length(centr),
        dimnames = list(cell_names, names(centr))
      )
      
      incoming <- matrix(
        0,
        nrow = length(cell_names),
        ncol = length(centr),
        dimnames = list(cell_names, names(centr))
      )
      
      for (i in seq_along(centr)) {
        
        out_use <- centr[[i]]$outdeg
        in_use <- centr[[i]]$indeg
        
        if (!is.null(names(out_use))) {
          outgoing[names(out_use), i] <- as.numeric(out_use)
        } else {
          outgoing[, i] <- as.numeric(out_use)
        }
        
        if (!is.null(names(in_use))) {
          incoming[names(in_use), i] <- as.numeric(in_use)
        } else {
          incoming[, i] <- as.numeric(in_use)
        }
      }
      
      if (pattern_single == "outgoing") {
        mat <- t(outgoing)
        legend.name <- "Outgoing"
      } else if (pattern_single == "incoming") {
        mat <- t(incoming)
        legend.name <- "Incoming"
      } else {
        mat <- t(outgoing + incoming)
        legend.name <- "Overall"
      }
      
      mat[is.nan(mat)] <- NA
      mat[is.infinite(mat)] <- NA
      
      list(
        mat = mat,
        legend.name = legend.name
      )
    }
    
    build_cellchat_style_role_heatmap_from_mat <- function(
    mat.ori,
    legend.name,
    title = NULL,
    color.use = NULL,
    color.heatmap = "BuGn",
    width = 10,
    height = 8,
    font.size = 8,
    font.size.title = 10,
    legend.size = 8
    ) {
      
      title_use <- if (is.null(title)) {
        paste0(legend.name, " signaling patterns")
      } else {
        paste0(legend.name, " signaling patterns - ", title)
      }
      
      mat.ori <- as.matrix(mat.ori)
      
      rownames(mat.ori) <- as.character(rownames(mat.ori))
      colnames(mat.ori) <- as.character(colnames(mat.ori))
      
      mat.ori <- mat.ori[
        !is.na(rownames(mat.ori)) & rownames(mat.ori) != "",
        !is.na(colnames(mat.ori)) & colnames(mat.ori) != "",
        drop = FALSE
      ]
      
      if (nrow(mat.ori) == 0 || ncol(mat.ori) == 0) {
        return(grid::textGrob(
          "No valid role matrix.",
          gp = grid::gpar(col = "red", fontsize = 14)
        ))
      }
      
      row_max <- apply(mat.ori, 1, function(x) {
        x <- x[is.finite(x)]
        if (length(x) == 0) return(NA_real_)
        max(x, na.rm = TRUE)
      })
      
      row_max[!is.finite(row_max) | row_max <= 0] <- NA_real_
      
      mat <- sweep(mat.ori, 1L, row_max, "/", check.margin = FALSE)
      mat[is.nan(mat)] <- NA
      mat[is.infinite(mat)] <- NA
      mat[mat == 0] <- NA
      
      color.heatmap.use <- circlize::colorRamp2(
        c(0, 0.5, 1),
        RColorBrewer::brewer.pal(
          n = 9,
          name = color.heatmap
        )[c(1, 5, 9)]
      )
      
      if (is.null(color.use)) {
        color.use <- get_scPalette(ncol(mat))
      }
      
      color.use <- as.character(color.use)
      
      if (length(color.use) < ncol(mat)) {
        extra_cols <- get_scPalette(ncol(mat) - length(color.use))
        color.use <- c(color.use, extra_cols)
      }
      
      color.use <- color.use[seq_len(ncol(mat))]
      names(color.use) <- colnames(mat)
      
      df <- data.frame(
        group = as.character(colnames(mat)),
        stringsAsFactors = FALSE
      )
      rownames(df) <- colnames(mat)
      
      col_annotation <- ComplexHeatmap::HeatmapAnnotation(
        df = df,
        col = list(group = color.use),
        which = "column",
        show_legend = FALSE,
        show_annotation_name = FALSE,
        simple_anno_size = grid::unit(0.2, "cm")
      )
      
      ha_top <- ComplexHeatmap::HeatmapAnnotation(
        Strength = ComplexHeatmap::anno_barplot(
          colSums(mat.ori, na.rm = TRUE),
          border = FALSE,
          gp = grid::gpar(
            fill = color.use,
            col = color.use
          )
        ),
        show_annotation_name = FALSE
      )
      
      pSum.original <- rowSums(mat.ori, na.rm = TRUE)
      pSum <- pSum.original
      
      pSum[pSum <= 0 | !is.finite(pSum)] <- NA_real_
      pSum <- -1 / log(pSum)
      pSum[is.na(pSum)] <- 0
      
      idx1 <- which(is.infinite(pSum) | pSum < 0)
      
      if (length(idx1) > 0) {
        max_p <- max(pSum[is.finite(pSum) & pSum >= 0], na.rm = TRUE)
        if (!is.finite(max_p)) max_p <- 1
        
        values.assign <- seq(
          max_p * 1.1,
          max_p * 1.5,
          length.out = length(idx1)
        )
        
        position <- sort(
          pSum.original[idx1],
          index.return = TRUE
        )$ix
        
        pSum[idx1] <- values.assign[
          match(seq_len(length(idx1)), position)
        ]
      }
      
      ha_right <- ComplexHeatmap::rowAnnotation(
        Strength = ComplexHeatmap::anno_barplot(
          pSum,
          border = FALSE
        ),
        show_annotation_name = FALSE
      )
      
      ComplexHeatmap::Heatmap(
        mat,
        col = color.heatmap.use,
        na_col = "white",
        name = "Relative strength",
        bottom_annotation = col_annotation,
        top_annotation = ha_top,
        right_annotation = ha_right,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_side = "left",
        row_names_rot = 0,
        row_names_gp = grid::gpar(
          fontsize = font.size
        ),
        column_names_gp = grid::gpar(
          fontsize = font.size
        ),
        width = grid::unit(width, "cm"),
        height = grid::unit(height, "cm"),
        column_title = title_use,
        column_title_gp = grid::gpar(
          fontsize = font.size.title
        ),
        column_names_rot = 90,
        heatmap_legend_param = list(
          title_gp = grid::gpar(
            fontsize = legend.size,
            fontface = "plain"
          ),
          title_position = "leftcenter-rot",
          border = NA,
          at = c(0, 0.5, 1),
          labels = c("0", "0.5", "1"),
          legend_height = grid::unit(20, "mm"),
          labels_gp = grid::gpar(
            fontsize = legend.size
          ),
          grid_width = grid::unit(2, "mm")
        )
      )
    }
    
    align_merged_role_matrices <- function(object.list, dataset_use, pattern_single) {
      
      raw_list <- lapply(dataset_use, function(dataset_name) {
        get_role_matrix_raw(
          obj = object.list[[dataset_name]],
          pattern_single = pattern_single,
          slot.name = "netP"
        )
      })
      
      names(raw_list) <- dataset_use
      
      if (any(sapply(raw_list, is.null))) {
        return(NULL)
      }
      
      mat_list <- lapply(raw_list, function(x) {
        x$mat
      })
      
      all_rows <- unique(unlist(lapply(mat_list, rownames)))
      all_cols <- unique(unlist(lapply(mat_list, colnames)))
      
      all_rows <- all_rows[!is.na(all_rows) & all_rows != ""]
      all_cols <- all_cols[!is.na(all_cols) & all_cols != ""]
      
      if (length(all_rows) == 0 || length(all_cols) == 0) {
        return(NULL)
      }
      
      mat_full_list <- lapply(mat_list, function(mat) {
        
        mat_full <- matrix(
          NA_real_,
          nrow = length(all_rows),
          ncol = length(all_cols),
          dimnames = list(all_rows, all_cols)
        )
        
        mat_full[rownames(mat), colnames(mat)] <- mat
        mat_full
      })
      
      row_score <- Reduce(
        "+",
        lapply(mat_full_list, function(mat) {
          rowSums(mat, na.rm = TRUE)
        })
      )
      
      row_keep <- Reduce(
        "|",
        lapply(mat_full_list, function(mat) {
          rowSums(!is.na(mat) & mat > 0, na.rm = TRUE) > 0
        })
      )
      
      row_use <- names(row_keep)[row_keep]
      
      if (length(row_use) == 0) {
        return(NULL)
      }
      
      row_use <- row_use[
        order(row_score[row_use], decreasing = TRUE)
      ]
      
      mat_full_list <- lapply(mat_full_list, function(mat) {
        mat[row_use, all_cols, drop = FALSE]
      })
      
      list(
        mat_list = mat_full_list,
        legend.name = raw_list[[1]]$legend.name
      )
    }
    
    draw_heatmap_to_grob <- function(ht) {
      grid::grid.grabExpr({
        if (inherits(ht, "grob") || inherits(ht, "gTree")) {
          grid::grid.draw(ht)
        } else {
          ComplexHeatmap::draw(ht)
        }
      })
    }
    
    if (!is.null(cellchat_merge)) {
      
      object.list <- cellchat_obj
      
      if (!is.list(object.list) || length(object.list) < 2) {
        return(list(
          type = "message",
          content = "The uploaded file has a merged CellChat object, but the original object list is not available."
        ))
      }
      
      object.list <- object.list[
        sapply(object.list, is_cellchat_obj)
      ]
      
      if (length(object.list) == 0) {
        return(list(
          type = "message",
          content = "No valid CellChat object was found in the object list."
        ))
      }
      
      dataset_use <- names(object.list)[1:min(2, length(object.list))]
      
      grob_matrix <- lapply(pattern_list, function(pattern_single) {
        
        aligned <- align_merged_role_matrices(
          object.list = object.list,
          dataset_use = dataset_use,
          pattern_single = pattern_single
        )
        
        if (is.null(aligned)) {
          return(list(
            grid::textGrob(
              paste0("No role matrix found for pattern: ", pattern_single),
              gp = grid::gpar(col = "red", fontsize = 14)
            )
          ))
        }
        
        ht_list <- lapply(seq_along(dataset_use), function(i) {
          
          dataset_name <- dataset_use[i]
          mat_use <- aligned$mat_list[[dataset_name]]
          
          n_pathway <- nrow(mat_use)
          n_celltype <- ncol(mat_use)
          
          build_cellchat_style_role_heatmap_from_mat(
            mat.ori = mat_use,
            legend.name = aligned$legend.name,
            title = dataset_name,
            color.use = NULL,
            color.heatmap = input$netAnalysis_signalingRole_heatmap_color,
            width = max(8, n_celltype * 0.8),
            height = max(10, n_pathway * 0.3),
            font.size = input$netAnalysis_signalingRole_heatmap_axis_text_size,
            font.size.title = input$netAnalysis_signalingRole_heatmap_axis_title_size,
            legend.size = input$netAnalysis_signalingRole_heatmap_legend_size
          )
        })
        
        lapply(ht_list, draw_heatmap_to_grob)
      })
      
      grob_list <- unlist(grob_matrix, recursive = FALSE)
      
      plot_fun <- function() {
        
        if (pattern_use == "both") {
          combined_grob <- gridExtra::arrangeGrob(
            grobs = grob_list,
            nrow = 2,
            ncol = 2,
            byrow = TRUE
          )
        } else {
          combined_grob <- gridExtra::arrangeGrob(
            grobs = grob_list,
            nrow = 1,
            ncol = 2
          )
        }
        
        grid::grid.newpage()
        grid::grid.draw(combined_grob)
        invisible(NULL)
      }
      
    } else {
      
      grob_list <- lapply(seq_along(pattern_list), function(i) {
        
        pattern_single <- pattern_list[i]
        
        raw <- get_role_matrix_raw(
          obj = cellchat_obj,
          pattern_single = pattern_single,
          slot.name = "netP"
        )
        
        if (is.null(raw)) {
          return(
            grid::textGrob(
              paste0("No role matrix found for pattern: ", pattern_single),
              gp = grid::gpar(col = "red", fontsize = 14)
            )
          )
        }
        
        mat_use <- raw$mat
        
        row_score <- rowSums(mat_use, na.rm = TRUE)
        row_use <- names(row_score)[row_score > 0]
        
        if (length(row_use) == 0) {
          return(
            grid::textGrob(
              paste0("No signaling role detected for pattern: ", pattern_single),
              gp = grid::gpar(col = "red", fontsize = 14)
            )
          )
        }
        
        row_use <- row_use[
          order(row_score[row_use], decreasing = TRUE)
        ]
        
        mat_use <- mat_use[row_use, , drop = FALSE]
        
        n_pathway <- nrow(mat_use)
        n_celltype <- ncol(mat_use)
        
        ht <- build_cellchat_style_role_heatmap_from_mat(
          mat.ori = mat_use,
          legend.name = raw$legend.name,
          title = NULL,
          color.use = NULL,
          color.heatmap = input$netAnalysis_signalingRole_heatmap_color,
          width = max(8, n_celltype * 0.8),
          height = max(10, n_pathway * 0.3),
          font.size = input$netAnalysis_signalingRole_heatmap_axis_text_size,
          font.size.title = input$netAnalysis_signalingRole_heatmap_axis_title_size,
          legend.size = input$netAnalysis_signalingRole_heatmap_legend_size
        )
        
        draw_heatmap_to_grob(ht)
      })
      
      plot_fun <- function() {
        
        if (pattern_use == "both") {
          combined_grob <- gridExtra::arrangeGrob(
            grobs = grob_list,
            nrow = 1,
            ncol = 2
          )
        } else {
          combined_grob <- gridExtra::arrangeGrob(
            grobs = grob_list,
            nrow = 1,
            ncol = 1
          )
        }
        
        grid::grid.newpage()
        grid::grid.draw(combined_grob)
        invisible(NULL)
      }
    }
    
    list(
      type = "plot",
      plot_fun = plot_fun
    )
  })
  output$netAnalysis_signalingRole_heatmap_control <- renderUI({
    
    has_merge <- !is.null(rv$cccd_obj$data_merge)
    
    is_both <- !is.null(input$netAnalysis_signalingRole_heatmap_pattern) &&
      input$netAnalysis_signalingRole_heatmap_pattern == "both"
    
    get_pathway_n <- function(obj) {
      if (is.null(obj)) return(30)
      if (methods::is(obj, "CellChat")) {
        if (!is.null(obj@netP$pathways)) {
          return(length(obj@netP$pathways))
        }
      }
      return(30)
    }
    
    get_cluster_n <- function(obj) {
      if (is.null(obj)) return(10)
      if (methods::is(obj, "CellChat")) {
        if (!is.null(obj@idents)) {
          return(length(levels(obj@idents)))
        }
        if (!is.null(obj@net$weight)) {
          return(ncol(obj@net$weight))
        }
      }
      return(10)
    }
    
    if (has_merge) {
      object.list <- rv$cccd_obj$data
      object.list <- object.list[
        sapply(object.list, function(x) methods::is(x, "CellChat"))
      ]
      
      pathway_n <- max(sapply(object.list, get_pathway_n))
      cluster_n <- max(sapply(object.list, get_cluster_n))
    } else {
      pathway_n <- get_pathway_n(rv$cccd_obj$data)
      cluster_n <- get_cluster_n(rv$cccd_obj$data)
    }
    
    single_panel_height <- max(100, pathway_n * 3) * 3.78
    single_panel_width <- max(80, cluster_n * 8) * 3.78
    
    if (is_both && !has_merge) {
      width_val <- paste0(single_panel_width * 2, "px")
      height_val <- paste0(single_panel_height, "px")
      
    } else if (is_both && has_merge) {
      width_val <- paste0(single_panel_width * 2, "px")
      height_val <- paste0(single_panel_height * 2, "px")
      
    } else if (!is_both && !has_merge) {
      width_val <- paste0(single_panel_width, "px")
      height_val <- paste0(single_panel_height, "px")
      
    } else {
      width_val <- paste0(single_panel_width * 2, "px")
      height_val <- paste0(single_panel_height, "px")
    }
    
    div(
      style = "max-width: 800px; max-height: 1600px; overflow-x: auto; overflow-y: auto; width: 100%;",
      plotOutput(
        "netAnalysis_signalingRole_heatmap_plot",
        height = height_val,
        width = width_val
      )
    )
  })
  output$netAnalysis_signalingRole_heatmap_plot <- renderPlot({
    
    req(netAnalysis_signalingRole_heatmap_reactive())
    req(!rv$ui_freeze)
    
    res <- netAnalysis_signalingRole_heatmap_reactive()
    
    if (is.list(res) && res$type == "message") {
      plot.new()
      text(0.5, 0.5, res$content, cex = 1.2, col = "red")
      return()
    }
    
    res$plot_fun()
  })
  
  observeEvent(input$save_netAnalysis_signalingRole_heatmap, {
    
    req(netAnalysis_signalingRole_heatmap_reactive())
    req(input$netAnalysis_signalingRole_heatmap_save_path)
    
    res <- netAnalysis_signalingRole_heatmap_reactive()
    
    if (is.list(res) && res$type == "message") {
      showNotification(res$content, type = "error")
      return()
    }
    
    out_dir <- input$netAnalysis_signalingRole_heatmap_save_path
    
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }
    
    out_file <- file.path(
      out_dir,
      paste0(
        "PathwayRoleHeatmap_",
        input$netAnalysis_signalingRole_heatmap_pattern,
        "_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".pdf"
      )
    )
    
    pdf(
      out_file,
      width = input$netAnalysis_signalingRole_heatmap_width,
      height = input$netAnalysis_signalingRole_heatmap_height
    )
    
    res$plot_fun()
    
    dev.off()
    
    showNotification(
      paste0("Saved pathway role heatmap: ", out_file),
      type = "message"
    )
  })
  output$download_netAnalysis_signalingRole_heatmap <- downloadHandler(
    
    filename = function() {
      paste0(
        "PathwayRoleHeatmap_",
        input$netAnalysis_signalingRole_heatmap_pattern,
        "_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".pdf"
      )
    },
    
    content = function(file) {
      
      req(netAnalysis_signalingRole_heatmap_reactive())
      
      res <- netAnalysis_signalingRole_heatmap_reactive()
      
      pdf(
        file,
        width = input$netAnalysis_signalingRole_heatmap_width,
        height = input$netAnalysis_signalingRole_heatmap_height
      )
      
      if (is.list(res) && res$type == "message") {
        plot.new()
        text(0.5, 0.5, res$content, cex = 1.2, col = "red")
      } else {
        res$plot_fun()
      }
      
      dev.off()
      
      showNotification(
        "Pathway role heatmap downloaded successfully.",
        type = "message",
        duration = 3
      )
    }
  )
  # ----------------- VISUALIZATION FUNCTION Cell Cell Communication netVisual_aggregate -----------------
  netVisual_aggregate_reactive <- reactive({
    
    req(rv$cccd_obj)
    req(rv$ccc_pathway_select)
    req(input$aggregate_layout)
    req(!rv$ui_freeze)
    
    cellchat_obj <- rv$cccd_obj$data
    cellchat_merge <- rv$cccd_obj$data_merge
    signaling_use <- rv$ccc_pathway_select
    
    if (is.null(signaling_use) || signaling_use == "") {
      return(list(type = "message", content = "Please select a pathway first."))
    }
    
    parse_vertex_receiver <- function(x, n_cell_groups) {
      x <- gsub(" ", "", x)
      x <- unlist(strsplit(x, ","))
      x <- suppressWarnings(as.integer(x))
      x <- unique(x[!is.na(x)])
      x <- x[x >= 1 & x <= n_cell_groups]
      if (length(x) == 0) return(NULL)
      x
    }
    
    plot_aggregate_one <- function(obj, dataset_name = NULL) {
      
      if (is.null(obj) || !methods::is(obj, "CellChat")) {
        plot.new()
        text(0.5, 0.5, "Invalid CellChat object.", col = "red", cex = 1.2)
        return(invisible(NULL))
      }
      
      if (!(signaling_use %in% obj@netP$pathways)) {
        plot.new()
        text(
          0.5, 0.5,
          paste0(
            "Pathway '", signaling_use, "' is not found",
            if (!is.null(dataset_name)) paste0(" in ", dataset_name) else "",
            "."
          ),
          col = "red",
          cex = 1.2
        )
        return(invisible(NULL))
      }
      
      title_use <- if (is.null(dataset_name)) {
        signaling_use
      } else {
        paste0(signaling_use, "_", dataset_name)
      }
      
      if (input$aggregate_layout == "circle") {
        
        netVisual_aggregate(
          obj,
          signaling = signaling_use,
          signaling.name = title_use,
          layout = "circle",
          vertex.label.cex = input$aggregate_vertex_label_size,
          edge.width.max = input$aggregate_edge_width_max,
          arrow.size = input$aggregate_arrow_size,
          title.space = 2
        )
        
        return(invisible(NULL))
      }
      
      if (input$aggregate_layout == "hierarchy") {
        
        n_cell_groups <- length(levels(obj@idents))
        vertex_receiver_use <- parse_vertex_receiver(
          input$aggregate_vertex_receiver,
          n_cell_groups
        )
        
        if (is.null(vertex_receiver_use)) {
          plot.new()
          text(
            0.5, 0.5,
            paste0(
              "Please enter valid receiver cell indices.\n",
              "Valid range: 1 to ", n_cell_groups
            ),
            col = "red",
            cex = 1.2
          )
          return(invisible(NULL))
        }
        
        netVisual_aggregate(
          obj,
          signaling = signaling_use,
          signaling.name = title_use,
          layout = "hierarchy",
          vertex.receiver = vertex_receiver_use,
          vertex.label.cex = input$aggregate_vertex_label_size,
          edge.width.max = input$aggregate_edge_width_max,
          arrow.size = input$aggregate_arrow_size,
          title.space = 2
        )
        
        return(invisible(NULL))
      }
      
      if (input$aggregate_layout == "chord") {
        
        netVisual_aggregate(
          obj,
          signaling = signaling_use,
          signaling.name = title_use,
          layout = "chord",
          lab.cex = input$aggregate_vertex_label_size,
          edge.width.max = input$aggregate_edge_width_max
        )
        
        return(invisible(NULL))
      }
      
      invisible(NULL)
    }
    
    if (!is.null(cellchat_merge)) {
      
      object.list <- cellchat_obj
      
      if (!is.list(object.list) || length(object.list) < 2) {
        return(list(
          type = "message",
          content = "The uploaded file has a merged CellChat object, but the original object list is not available."
        ))
      }
      
      object.list <- object.list[
        sapply(object.list, function(x) methods::is(x, "CellChat"))
      ]
      
      if (length(object.list) == 0) {
        return(list(
          type = "message",
          content = "No valid CellChat object was found in the object list."
        ))
      }
      
      plot_fun <- function() {
        
        old_par <- par(no.readonly = TRUE)
        on.exit(par(old_par), add = TRUE)
        
        dataset_use <- names(object.list)[1:min(2, length(object.list))]
        
        par(mfrow = c(1, length(dataset_use)), xpd = TRUE)
        
        for (dataset_name in dataset_use) {
          plot_aggregate_one(
            obj = object.list[[dataset_name]],
            dataset_name = dataset_name
          )
        }
        
        invisible(NULL)
      }
      
    } else {
      
      plot_fun <- function() {
        
        old_par <- par(no.readonly = TRUE)
        on.exit(par(old_par), add = TRUE)
        
        par(mfrow = c(1, 1), xpd = TRUE)
        
        plot_aggregate_one(
          obj = cellchat_obj,
          dataset_name = NULL
        )
        
        invisible(NULL)
      }
    }
    
    list(type = "plot", plot_fun = plot_fun)
  })
  
  output$netVisual_aggregate_control <- renderUI({
    
    if (!is.null(rv$cccd_obj$data_merge)) {
      width_val <- "800px"
    } else {
      width_val <- "400px"
    }
    plotOutput("netVisual_aggregate_plot", height = "400px", width = width_val)
  })
  
  output$netVisual_aggregate_plot <- renderPlot({
    
    req(netVisual_aggregate_reactive())
    req(!rv$ui_freeze)
    
    res <- netVisual_aggregate_reactive()
    
    if (is.list(res) && res$type == "message") {
      plot.new()
      text(0.5, 0.5, res$content, cex = 1.2, col = "red")
      return()
    }
    
    res$plot_fun()
  })
  
  observeEvent(input$save_netVisual_aggregate, {
    
    req(netVisual_aggregate_reactive())
    req(input$aggregate_save_path)
    
    res <- netVisual_aggregate_reactive()
    
    if (is.list(res) && res$type == "message") {
      showNotification(res$content, type = "error")
      return()
    }
    
    out_dir <- input$aggregate_save_path
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    
    out_file <- file.path(
      out_dir,
      paste0(
        "PathwayNetwork_",
        rv$ccc_pathway_select,
        "_",
        input$aggregate_layout,
        "_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".pdf"
      )
    )
    
    pdf(
      out_file,
      width = input$aggregate_width,
      height = input$aggregate_height
    )
    
    res$plot_fun()
    
    dev.off()
    
    showNotification(
      paste0("Saved pathway network: ", out_file),
      type = "message"
    )
  })
  
  output$download_netVisual_aggregate <- downloadHandler(
    
    filename = function() {
      paste0(
        "PathwayNetwork_",
        rv$ccc_pathway_select,
        "_",
        input$aggregate_layout,
        "_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".pdf"
      )
    },
    
    content = function(file) {
      
      req(netVisual_aggregate_reactive())
      
      res <- netVisual_aggregate_reactive()
      
      pdf(
        file,
        width = input$aggregate_width,
        height = input$aggregate_height
      )
      
      if (is.list(res) && res$type == "message") {
        plot.new()
        text(0.5, 0.5, res$content, cex = 1.2, col = "red")
      } else {
        res$plot_fun()
      }
      
      dev.off()
      
      showNotification(
        "Pathway network downloaded successfully.",
        type = "message",
        duration = 3
      )
    }
  )
  # ----------------- VISUALIZATION FUNCTION Cell Cell Communication netVisual_heatmap -----------------
  netVisual_heatmap_reactive <- reactive({
    
    req(rv$cccd_obj)
    req(rv$ccc_pathway_select)
    req(input$netVisual_heatmap_color)
    req(input$netVisual_heatmap_axis_title_size)
    req(input$netVisual_heatmap_axis_text_size)
    req(input$netVisual_heatmap_legend_size)
    req(!rv$ui_freeze)
    
    cellchat_obj <- rv$cccd_obj$data
    cellchat_merge <- rv$cccd_obj$data_merge
    signaling_use <- trimws(rv$ccc_pathway_select)
    
    if (is.null(signaling_use) || signaling_use == "") {
      return(list(type = "message", content = "Please select a pathway first."))
    }
    
    is_cellchat_obj <- function(x) {
      methods::is(x, "CellChat")
    }
    
    pathway_exists <- function(obj) {
      if (is.null(obj)) return(FALSE)
      if (!is_cellchat_obj(obj)) return(FALSE)
      if (is.null(obj@netP$pathways)) return(FALSE)
      signaling_use %in% trimws(obj@netP$pathways)
    }
    
    make_message_grob <- function(msg) {
      grid::grobTree(
        grid::rectGrob(gp = grid::gpar(col = NA, fill = "white")),
        grid::textGrob(
          msg,
          x = 0.5,
          y = 0.5,
          gp = grid::gpar(
            col = "red",
            fontsize = input$netVisual_heatmap_axis_text_size + 4
          )
        )
      )
    }
    
    make_pathway_heatmap <- function(obj, dataset_name = NULL) {
      
      if (!pathway_exists(obj)) {
        msg <- paste0(
          "Pathway '", signaling_use, "' is not found",
          if (!is.null(dataset_name)) paste0(" in ", dataset_name) else "",
          "."
        )
        return(make_message_grob(msg))
      }
      
      title_use <- if (is.null(dataset_name)) {
        signaling_use
      } else {
        paste0(signaling_use, "_", dataset_name)
      }
      
      ht <- netVisual_heatmap(
        obj,
        signaling = signaling_use,
        slot.name = "netP",
        color.heatmap = input$netVisual_heatmap_color,
        title.name = title_use,
        font.size = input$netVisual_heatmap_axis_text_size
      )
      
      ht@row_names_param$gp <- grid::gpar(fontsize = input$netVisual_heatmap_axis_text_size)
      ht@column_names_param$gp <- grid::gpar(fontsize = input$netVisual_heatmap_axis_text_size)
      ht@row_title_param$gp <- grid::gpar(fontsize = input$netVisual_heatmap_axis_title_size)
      ht@column_title_param$gp <- grid::gpar(fontsize = input$netVisual_heatmap_axis_title_size)
      ht@matrix_legend_param$title_gp <- grid::gpar(fontsize = input$netVisual_heatmap_legend_size)
      ht@matrix_legend_param$labels_gp <- grid::gpar(fontsize = input$netVisual_heatmap_legend_size)
      
      grid::grid.grabExpr({
        ComplexHeatmap::draw(ht)
      })
    }
    
    if (!is.null(cellchat_merge)) {
      
      object.list <- cellchat_obj
      
      if (!is.list(object.list) || length(object.list) < 2) {
        return(list(
          type = "message",
          content = "The uploaded file has a merged CellChat object, but the original object list is not available."
        ))
      }
      
      object.list <- object.list[sapply(object.list, is_cellchat_obj)]
      
      if (length(object.list) == 0) {
        return(list(type = "message", content = "No valid CellChat object was found in the object list."))
      }
      
      dataset_use <- names(object.list)[1:min(2, length(object.list))]
      
      grob_list <- lapply(dataset_use, function(dataset_name) {
        make_pathway_heatmap(
          obj = object.list[[dataset_name]],
          dataset_name = dataset_name
        )
      })
      
      plot_fun <- function() {
        combined_grob <- gridExtra::arrangeGrob(
          grobs = grob_list,
          nrow = 1,
          ncol = length(dataset_use)
        )
        
        grid::grid.newpage()
        grid::grid.draw(combined_grob)
        invisible(NULL)
      }
      
    } else {
      
      grob_use <- make_pathway_heatmap(
        obj = cellchat_obj,
        dataset_name = NULL
      )
      
      plot_fun <- function() {
        grid::grid.newpage()
        grid::grid.draw(grob_use)
        invisible(NULL)
      }
    }
    
    list(type = "plot", plot_fun = plot_fun)
  })
  
  output$netVisual_heatmap_control <- renderUI({
    
    if (!is.null(rv$cccd_obj$data_merge)) {
      width_val <- "800px"
    } else {
      width_val <- "400px"
    }
    div(
      style = "max-width: 800px; overflow-x: auto; width: 100%;",
    plotOutput(
      "netVisual_heatmap_plot",
      height = "350px",
      width = width_val
    )
    )
    
  })
  
  output$netVisual_heatmap_plot <- renderPlot({
    
    req(netVisual_heatmap_reactive())
    req(!rv$ui_freeze)
    
    res <- netVisual_heatmap_reactive()
    
    if (is.list(res) && res$type == "message") {
      
      plot.new()
      
      text(
        0.5,
        0.5,
        res$content,
        cex = 1.2,
        col = "red"
      )
      
      return()
    }
    
    res$plot_fun()
    
  })
  
  observeEvent(input$save_netVisual_heatmap, {
    
    req(netVisual_heatmap_reactive())
    req(input$netVisual_heatmap_save_path)
    
    res <- netVisual_heatmap_reactive()
    
    if (is.list(res) && res$type == "message") {
      
      showNotification(
        res$content,
        type = "error"
      )
      
      return()
    }
    
    out_dir <- input$netVisual_heatmap_save_path
    
    if (!dir.exists(out_dir)) {
      dir.create(
        out_dir,
        recursive = TRUE
      )
    }
    
    out_file <- file.path(
      out_dir,
      paste0(
        "PathwayHeatmap_",
        rv$ccc_pathway_select,
        "_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".pdf"
      )
    )
    
    pdf(
      out_file,
      width = input$netVisual_heatmap_width,
      height = input$netVisual_heatmap_height
    )
    
    res$plot_fun()
    
    dev.off()
    
    showNotification(
      paste0(
        "Saved pathway heatmap: ",
        out_file
      ),
      type = "message"
    )
    
  })
  
  output$download_netVisual_heatmap <- downloadHandler(
    
    filename = function() {
      
      paste0(
        "PathwayHeatmap_",
        rv$ccc_pathway_select,
        "_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".pdf"
      )
      
    },
    
    content = function(file) {
      
      req(netVisual_heatmap_reactive())
      
      res <- netVisual_heatmap_reactive()
      
      pdf(
        file,
        width = input$netVisual_heatmap_width,
        height = input$netVisual_heatmap_height
      )
      
      if (is.list(res) && res$type == "message") {
        
        plot.new()
        
        text(
          0.5,
          0.5,
          res$content,
          cex = 1.2,
          col = "red"
        )
        
      } else {
        
        res$plot_fun()
        
      }
      
      dev.off()
      
      showNotification(
        "Pathway heatmap downloaded successfully.",
        type = "message",
        duration = 3
      )
      
    }
    
  )
  

  # ----------------- VISUALIZATION FUNCTION Cell Cell Communication netAnalysis_signalingRole_network -----------------
  netAnalysis_signalingRole_network_reactive <- reactive({
    
    req(rv$cccd_obj)
    req(rv$ccc_pathway_select)
    req(input$netAnalysis_signalingRole_network_color)
    req(input$netAnalysis_signalingRole_network_font_size)
    req(input$netAnalysis_signalingRole_network_title_size)
    req(!rv$ui_freeze)
    
    cellchat_obj <- rv$cccd_obj$data
    cellchat_merge <- rv$cccd_obj$data_merge
    signaling_use <- trimws(rv$ccc_pathway_select)
    
    if (is.null(signaling_use) || signaling_use == "") {
      return(list(type = "message", content = "Please select and confirm a pathway first."))
    }
    
    is_cellchat_obj <- function(x) {
      methods::is(x, "CellChat")
    }
    
    pathway_exists <- function(obj) {
      if (is.null(obj)) return(FALSE)
      if (!is_cellchat_obj(obj)) return(FALSE)
      if (is.null(obj@netP$pathways)) return(FALSE)
      signaling_use %in% trimws(obj@netP$pathways)
    }
    
    make_message_grob <- function(msg) {
      grid::grobTree(
        grid::rectGrob(gp = grid::gpar(col = NA, fill = "white")),
        grid::textGrob(
          msg,
          x = 0.5,
          y = 0.5,
          gp = grid::gpar(
            col = "red",
            fontsize = input$netAnalysis_signalingRole_network_font_size + 4
          )
        )
      )
    }
    
    draw_single_network <- function(obj, dataset_name = NULL) {
      
      if (!pathway_exists(obj)) {
        msg <- paste0(
          "Pathway '", signaling_use, "' is not found",
          if (!is.null(dataset_name)) paste0(" in ", dataset_name) else "",
          "."
        )
        return(make_message_grob(msg))
      }
      
      grid::grid.grabExpr({
        netAnalysis_signalingRole_network(
          obj,
          signaling = signaling_use,
          slot.name = "netP",
          color.heatmap = input$netAnalysis_signalingRole_network_color,
          width = if (!is.null(rv$cccd_obj$data_merge)) {
            input$netAnalysis_signalingRole_network_width / 2
          } else {
            input$netAnalysis_signalingRole_network_width
          },
          height = input$netAnalysis_signalingRole_network_height - 1,
          font.size = input$netAnalysis_signalingRole_network_font_size,
          font.size.title = input$netAnalysis_signalingRole_network_title_size
        )
        
        if (!is.null(dataset_name)) {
          grid::grid.text(
            dataset_name,
            x = 0.5,
            y = 0.98,
            gp = grid::gpar(
              fontsize = input$netAnalysis_signalingRole_network_title_size + 2,
              fontface = "bold"
            )
          )
        }
      })
    }
    
    if (!is.null(cellchat_merge)) {
      
      object.list <- cellchat_obj
      
      if (!is.list(object.list) || length(object.list) < 2) {
        return(list(
          type = "message",
          content = "The uploaded file has a merged CellChat object, but the original object list is not available."
        ))
      }
      
      object.list <- object.list[sapply(object.list, is_cellchat_obj)]
      
      if (length(object.list) == 0) {
        return(list(type = "message", content = "No valid CellChat object was found in the object list."))
      }
      
      dataset_use <- names(object.list)[1:min(2, length(object.list))]
      
      grob_list <- lapply(dataset_use, function(dataset_name) {
        draw_single_network(
          obj = object.list[[dataset_name]],
          dataset_name = dataset_name
        )
      })
      
      plot_fun <- function() {
        combined_grob <- gridExtra::arrangeGrob(
          grobs = grob_list,
          nrow = 1,
          ncol = length(dataset_use)
        )
        
        grid::grid.newpage()
        grid::grid.draw(combined_grob)
        invisible(NULL)
      }
      
    } else {
      
      grob_use <- draw_single_network(
        obj = cellchat_obj,
        dataset_name = NULL
      )
      
      plot_fun <- function() {
        grid::grid.newpage()
        grid::grid.draw(grob_use)
        invisible(NULL)
      }
    }
    
    list(type = "plot", plot_fun = plot_fun)
  })
  
  output$netAnalysis_signalingRole_network_control <- renderUI({
    
    has_merge <- !is.null(rv$cccd_obj$data_merge)
    
    if (has_merge) {
      width_val <- "1080px"
    } else {
      width_val <- "540px"
    }
    
    height_val <- "300px"
    
    div(
      style = "max-width: 800px; overflow-x: auto; width: 100%;",
      plotOutput(
        "netAnalysis_signalingRole_network_plot",
        height = height_val,
        width = width_val
      )
    )
  })
  output$netAnalysis_signalingRole_network_plot <- renderPlot({
    
    req(netAnalysis_signalingRole_network_reactive())
    req(!rv$ui_freeze)
    
    res <- netAnalysis_signalingRole_network_reactive()
    
    if (is.list(res) && res$type == "message") {
      plot.new()
      text(0.5, 0.5, res$content, cex = 1.2, col = "red")
      return()
    }
    
    res$plot_fun()
  })
  observeEvent(input$save_netAnalysis_signalingRole_network, {
    
    req(netAnalysis_signalingRole_network_reactive())
    req(input$netAnalysis_signalingRole_network_save_path)
    
    res <- netAnalysis_signalingRole_network_reactive()
    
    if (is.list(res) && res$type == "message") {
      showNotification(res$content, type = "error")
      return()
    }
    
    out_dir <- input$netAnalysis_signalingRole_network_save_path
    
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }
    
    out_file <- file.path(
      out_dir,
      paste0(
        "PathwayRoleNetwork_",
        rv$ccc_pathway_select,
        "_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".pdf"
      )
    )
    
    pdf(
      out_file,
      width = input$netAnalysis_signalingRole_network_width,
      height = input$netAnalysis_signalingRole_network_height
    )
    
    res$plot_fun()
    
    dev.off()
    
    showNotification(
      paste0("Saved pathway role network: ", out_file),
      type = "message"
    )
  })
  output$download_netAnalysis_signalingRole_network <- downloadHandler(
    
    filename = function() {
      paste0(
        "PathwayRoleNetwork_",
        rv$ccc_pathway_select,
        "_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".pdf"
      )
    },
    
    content = function(file) {
      
      req(netAnalysis_signalingRole_network_reactive())
      
      res <- netAnalysis_signalingRole_network_reactive()
      
      pdf(
        file,
        width = input$netAnalysis_signalingRole_network_width,
        height = input$netAnalysis_signalingRole_network_height
      )
      
      if (is.list(res) && res$type == "message") {
        plot.new()
        text(0.5, 0.5, res$content, cex = 1.2, col = "red")
      } else {
        res$plot_fun()
      }
      
      dev.off()
      
      showNotification(
        "Pathway role network downloaded successfully.",
        type = "message",
        duration = 3
      )
    }
  )
  
  # ----------------- VISUALIZATION FUNCTION Cell Cell Communication netAnalysis_contribution -----------------
  netAnalysis_contribution_reactive <- reactive({
    
    req(rv$cccd_obj)
    req(rv$ccc_pathway_select)
    req(input$netAnalysis_contribution_x_text_size)
    req(input$netAnalysis_contribution_lr_pair_size)
    req(input$netAnalysis_contribution_title_size)
    req(!rv$ui_freeze)
    
    cellchat_obj <- rv$cccd_obj$data
    cellchat_merge <- rv$cccd_obj$data_merge
    signaling_use <- trimws(rv$ccc_pathway_select)
    
    if (is.null(signaling_use) || signaling_use == "") {
      return(list(
        type = "message",
        content = "Please select and confirm a pathway first."
      ))
    }
    
    is_cellchat_obj <- function(x) {
      methods::is(x, "CellChat")
    }
    
    pathway_exists <- function(obj) {
      if (is.null(obj)) return(FALSE)
      if (!is_cellchat_obj(obj)) return(FALSE)
      if (is.null(obj@netP$pathways)) return(FALSE)
      signaling_use %in% trimws(obj@netP$pathways)
    }
    
    make_message_grob <- function(msg) {
      grid::grobTree(
        grid::rectGrob(
          gp = grid::gpar(col = NA, fill = "white")
        ),
        grid::textGrob(
          msg,
          x = 0.5,
          y = 0.5,
          gp = grid::gpar(
            col = "red",
            fontsize = input$netAnalysis_contribution_font_size + 4
          )
        )
      )
    }
    
    make_contribution_grob <- function(obj, dataset_name = NULL) {
      
      if (!pathway_exists(obj)) {
        msg <- paste0(
          "Pathway '", signaling_use, "' is not found",
          if (!is.null(dataset_name)) paste0(" in ", dataset_name) else "",
          "."
        )
        return(make_message_grob(msg))
      }
      
      title_use <- if (is.null(dataset_name)) {
        paste0(signaling_use, " signaling contribution")
      } else {
        paste0(signaling_use, " signaling contribution - ", dataset_name)
      }
      
      p0 <- netAnalysis_contribution(
        obj,
        signaling = signaling_use
      )
      
      plot_df <- p0$data
      
      plot_df <- plot_df[
        plot_df$contribution > 1e-10 &
          !is.na(plot_df$name) &
          !grepl("^[0-9]+$", as.character(plot_df$name)),
      ]
      
      plot_df$name <- factor(
        plot_df$name,
        levels = rev(as.character(plot_df$name))
      )
      
      p <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(
          x = contribution,
          y = name
        )
      ) +
        ggplot2::geom_col(
          fill = "gray35",
          width = 0.7
        ) +
        ggplot2::labs(
          title = title_use,
          x = "Relative contribution",
          y = NULL
        ) +
        ggplot2::theme_classic() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            size = input$netAnalysis_contribution_title_size,
            face = "bold",
            hjust = 0.5
          ),
          axis.title.x = ggplot2::element_text(
            size = input$netAnalysis_contribution_x_text_size
          ),
          axis.text.x = ggplot2::element_text(
            size = input$netAnalysis_contribution_x_text_size
          ),
          axis.text.y = ggplot2::element_text(
            size = input$netAnalysis_contribution_lr_pair_size
          )
        )
      
      if (isTRUE(input$netAnalysis_contribution_show_value)) {
        p <- p +
          ggplot2::geom_text(
            ggplot2::aes(
              x = contribution + 0.01,
              label = scales::number(contribution, accuracy = 0.001)
            ),
            hjust = 0
          ) +
          ggplot2::coord_cartesian(
            xlim = c(0, max(plot_df$contribution, na.rm = TRUE) * 1.15)
          )
      }
      
      ggplot2::ggplotGrob(p)
      
      ggplot2::ggplotGrob(p)
    }
    
    if (!is.null(cellchat_merge)) {
      
      object.list <- cellchat_obj
      
      if (!is.list(object.list) || length(object.list) < 2) {
        return(list(
          type = "message",
          content = "The uploaded file has a merged CellChat object, but the original object list is not available."
        ))
      }
      
      object.list <- object.list[
        sapply(object.list, is_cellchat_obj)
      ]
      
      if (length(object.list) == 0) {
        return(list(
          type = "message",
          content = "No valid CellChat object was found in the object list."
        ))
      }
      
      dataset_use <- names(object.list)[1:min(2, length(object.list))]
      
      grob_list <- lapply(dataset_use, function(dataset_name) {
        make_contribution_grob(
          obj = object.list[[dataset_name]],
          dataset_name = dataset_name
        )
      })
      
      plot_fun <- function() {
        combined_grob <- gridExtra::arrangeGrob(
          grobs = grob_list,
          nrow = 1,
          ncol = length(grob_list)
        )
        
        grid::grid.newpage()
        grid::grid.draw(combined_grob)
        invisible(NULL)
      }
      
    } else {
      
      grob_use <- make_contribution_grob(
        obj = cellchat_obj,
        dataset_name = NULL
      )
      
      plot_fun <- function() {
        grid::grid.newpage()
        grid::grid.draw(grob_use)
        invisible(NULL)
      }
    }
    
    list(
      type = "plot",
      plot_fun = plot_fun
    )
  })
  
  output$netAnalysis_contribution_control <- renderUI({
    
    has_merge <- !is.null(rv$cccd_obj$data_merge)
    
    width_val <- if (has_merge) {
      "900px"
    } else {
      "450px"
    }
    
    div(
      style = "max-width: 800px; overflow-x: auto; width: 100%;",
      plotOutput(
        "netAnalysis_contribution_plot",
        height = "400px",
        width = width_val
      )
    )
  })
  
  
  output$netAnalysis_contribution_plot <- renderPlot({
    
    req(netAnalysis_contribution_reactive())
    req(!rv$ui_freeze)
    
    res <- netAnalysis_contribution_reactive()
    
    if (is.list(res) && res$type == "message") {
      plot.new()
      text(0.5, 0.5, res$content, cex = 1.2, col = "red")
      return()
    }
    
    res$plot_fun()
  })
  
  observeEvent(input$save_netAnalysis_contribution, {
    
    req(netAnalysis_contribution_reactive())
    req(input$netAnalysis_contribution_save_path)
    
    res <- netAnalysis_contribution_reactive()
    
    if (is.list(res) && res$type == "message") {
      showNotification(res$content, type = "error")
      return()
    }
    
    out_dir <- input$netAnalysis_contribution_save_path
    
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }
    
    out_file <- file.path(
      out_dir,
      paste0(
        "LRContribution_",
        rv$ccc_pathway_select,
        "_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".pdf"
      )
    )
    
    pdf(
      out_file,
      width = input$netAnalysis_contribution_width,
      height = input$netAnalysis_contribution_height
    )
    
    res$plot_fun()
    
    dev.off()
    
    showNotification(
      paste0("Saved L-R contribution plot: ", out_file),
      type = "message"
    )
  })
  
  output$download_netAnalysis_contribution <- downloadHandler(
    
    filename = function() {
      paste0(
        "LRContribution_",
        rv$ccc_pathway_select,
        "_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".pdf"
      )
    },
    
    content = function(file) {
      
      req(netAnalysis_contribution_reactive())
      
      res <- netAnalysis_contribution_reactive()
      
      pdf(
        file,
        width = input$netAnalysis_contribution_width,
        height = input$netAnalysis_contribution_height
      )
      
      if (is.list(res) && res$type == "message") {
        plot.new()
        text(0.5, 0.5, res$content, cex = 1.2, col = "red")
      } else {
        res$plot_fun()
      }
      
      dev.off()
      
      showNotification(
        "L-R contribution plot downloaded successfully.",
        type = "message",
        duration = 3
      )
    }
  )
  #----------------- VISUALIZATION FUNCTION Cell Cell Communication netVisual_individual -----------------
  
  plot_individual_message <- function(msg) {
    plot.new()
    text(0.5, 0.5, msg, col = "red", cex = 1.2)
    invisible(NULL)
  }
  
  plot_individual_one <- function(obj, lr_pair, dataset_name = NULL) {
    
    if (is.null(obj) || !methods::is(obj, "CellChat")) {
      plot_individual_message("Invalid CellChat object.")
      return(invisible(NULL))
    }
    
    if (is.null(obj@LR$LRsig) || !("interaction_name" %in% colnames(obj@LR$LRsig))) {
      plot_individual_message("L-R pair information is not available.")
      return(invisible(NULL))
    }
    
    lr_index <- which(obj@LR$LRsig$interaction_name == lr_pair)[1]
    
    if (is.na(lr_index)) {
      plot_individual_message(
        paste0(
          "L-R pair '", lr_pair, "' is not found",
          if (!is.null(dataset_name)) paste0(" in ", dataset_name) else "",
          "."
        )
      )
      return(invisible(NULL))
    }
    
    if (is.null(obj@net$prob)) {
      plot_individual_message("obj@net$prob is not available.")
      return(invisible(NULL))
    }
    
    mat <- obj@net$prob[, , lr_index]
    
    if (!is.null(dimnames(obj@net$weight))) {
      dimnames(mat) <- dimnames(obj@net$weight)
    }
    
    if (all(mat == 0, na.rm = TRUE)) {
      plot_individual_message(
        paste0(
          "L-R pair '", lr_pair, "' has no detected communication",
          if (!is.null(dataset_name)) paste0(" in ", dataset_name) else "",
          "."
        )
      )
      return(invisible(NULL))
    }
    
    title_use <- if (is.null(dataset_name)) {
      lr_pair
    } else {
      paste0(lr_pair, "_", dataset_name)
    }
    
    g <- igraph::graph_from_adjacency_matrix(
      mat,
      mode = "directed",
      weighted = TRUE,
      diag = TRUE
    )
    
    n <- nrow(mat)
    theta <- seq(0, 2 * pi, length.out = n + 1)[1:n]
    layout_circle <- cbind(cos(theta), sin(theta))
    rownames(layout_circle) <- rownames(mat)
    
    vertex_colors <- tryCatch(
      CellChat::scPalette(n),
      error = function(e) grDevices::rainbow(n)
    )
    names(vertex_colors) <- rownames(mat)
    
    edge_ends <- igraph::ends(g, igraph::E(g), names = TRUE)
    edge_from <- edge_ends[, 1]
    edge_to <- edge_ends[, 2]
    edge_weights <- igraph::E(g)$weight
    
    max_weight <- max(edge_weights, na.rm = TRUE)
    
    if (!is.finite(max_weight) || max_weight <= 0) {
      plot_individual_message(
        paste0(
          "L-R pair '", lr_pair, "' has no positive edge weight",
          if (!is.null(dataset_name)) paste0(" in ", dataset_name) else "",
          "."
        )
      )
      return(invisible(NULL))
    }
    
    edge_width <- edge_weights / max_weight * input$individual_edge_width_max
    edge_colors <- grDevices::adjustcolor(vertex_colors[edge_from], alpha.f = 0.6)
    
    loop_angles <- rep(0, igraph::ecount(g))
    
    for (i in seq_len(igraph::ecount(g))) {
      if (edge_from[i] == edge_to[i]) {
        vertex_index <- match(edge_from[i], rownames(mat))
        loop_angles[i] <- -theta[vertex_index]
      }
    }
    
    plot(
      g,
      layout = layout_circle,
      vertex.label = NA,
      vertex.size = 15,
      vertex.color = vertex_colors,
      vertex.frame.color = NA,
      edge.width = edge_width,
      edge.color = edge_colors,
      edge.arrow.size = input$individual_arrow_size,
      edge.curved = 0.25,
      loop.angle = loop_angles,
      loop.size = 1.2,
      main = "",
      margin = 0.2
    )
    
    label_radius <- 1.15
    label_adj <- ifelse(cos(theta) >= 0, 0, 1)
    
    text(
      x = cos(theta) * label_radius,
      y = sin(theta) * label_radius,
      labels = rownames(mat),
      cex = input$individual_vertex_label_size,
      col = "black",
      font = 1,
      family = "sans",
      adj = label_adj
    )
    
    title(main = title_use)
    invisible(NULL)
  }
  
  netVisual_individual_reactive <- reactive({
    
    req(rv$cccd_obj)
    req(rv$ccc_lr_pair_select)
    req(!rv$ui_freeze)
    
    cellchat_obj <- rv$cccd_obj$data
    cellchat_merge <- rv$cccd_obj$data_merge
    
    lr_pairs <- rv$ccc_lr_pair_select
    
    if (length(lr_pairs) == 0) {
      return(list(
        type = "message",
        content = "Please select at least one L-R pair first."
      ))
    }
    
    lr_pair_use <- lr_pairs[1]
    
    if (!is.null(cellchat_merge)) {
      
      object.list <- cellchat_obj
      
      if (!is.list(object.list) || length(object.list) < 2) {
        return(list(
          type = "message",
          content = "The uploaded file has a merged CellChat object, but the original object list is not available."
        ))
      }
      
      dataset_use <- names(object.list)[1:min(2, length(object.list))]
      
      plot_fun <- function() {
        
        old_par <- par(no.readonly = TRUE)
        on.exit(par(old_par), add = TRUE)
        
        layout(matrix(seq_along(dataset_use), nrow = 1))
        par(xpd = TRUE)
        
        for (dataset_name in dataset_use) {
          plot_individual_one(
            obj = object.list[[dataset_name]],
            lr_pair = lr_pair_use,
            dataset_name = dataset_name
          )
        }
        
        invisible(NULL)
      }
      
    } else {
      
      plot_fun <- function() {
        old_par <- par(no.readonly = TRUE)
        on.exit(par(old_par), add = TRUE)
        
        layout(matrix(1, nrow = 1))
        par(xpd = TRUE)
        
        plot_individual_one(
          obj = cellchat_obj,
          lr_pair = lr_pair_use,
          dataset_name = NULL
        )
        
        invisible(NULL)
      }
    }
    
    list(
      type = "plot",
      plot_fun = plot_fun,
      lr_pair = lr_pair_use
    )
  })
  output$netVisual_individual_control <- renderUI({
    
    width_val <- if (!is.null(rv$cccd_obj$data_merge)) {
      "800px"
    } else {
      "400px"
    }
    
    div(
      style = "max-width: 800px; overflow-x: auto; width: 100%;",
      plotOutput(
        "netVisual_individual_plot",
        height = "400px",
        width = width_val
      )
    )
  })
  
  output$netVisual_individual_plot <- renderPlot({
    
    req(netVisual_individual_reactive())
    req(!rv$ui_freeze)
    
    res <- netVisual_individual_reactive()
    
    if (is.list(res) && res$type == "message") {
      plot.new()
      text(0.5, 0.5, res$content, cex = 1.2, col = "red")
      return()
    }
    
    res$plot_fun()
  })
  
  observeEvent(input$generate_netVisual_individual, {
    
    req(rv$cccd_obj)
    req(input$individual_save_path)
    
    lr_pairs <- rv$ccc_lr_pair_select
    
    if (is.null(lr_pairs) || length(lr_pairs) == 0) {
      showNotification(
        "Please select at least one L-R pair first.",
        type = "error"
      )
      return()
    }
    
    out_dir <- input$individual_save_path
    
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }
    
    cellchat_obj <- rv$cccd_obj$data
    cellchat_merge <- rv$cccd_obj$data_merge
    has_merge <- !is.null(cellchat_merge)
    
    if (has_merge) {
      object.list <- cellchat_obj
      
      if (!is.list(object.list) || length(object.list) < 2) {
        showNotification(
          "The uploaded file has a merged CellChat object, but the original object list is not available.",
          type = "error"
        )
        return()
      }
      
      dataset_use <- names(object.list)[1:min(2, length(object.list))]
    }
    
    showModal(
      modalDialog(
        title = "Generating L-R Pair Network",
        tags$div(id = "netVisual_individual_progress_detail", "Starting..."),
        easyClose = FALSE,
        footer = NULL
      )
    )
    
    generated_files <- character(0)
    
    for (i in seq_along(lr_pairs)) {
      
      lr_pair <- lr_pairs[i]
      
      shinyjs::html(
        "netVisual_individual_progress_detail",
        paste0(
          "<b>Generating L-R pair network:</b> ", lr_pair, "<br>",
          "Completed: ", i - 1, "<br>",
          "Remaining: ", length(lr_pairs) - i + 1
        )
      )
      
      Sys.sleep(0.01)
      
      out_file <- file.path(
        out_dir,
        paste0(
          "LRPairNetwork_",
          lr_pair,
          "_",
          format(Sys.time(), "%Y%m%d_%H%M%S"),
          ".pdf"
        )
      )
      
      tryCatch({
        
        pdf(
          out_file,
          width = input$individual_width,
          height = input$individual_height
        )
        
        if (has_merge) {
          
          old_par <- par(no.readonly = TRUE)
          on.exit(par(old_par), add = TRUE)
          
          layout(matrix(seq_along(dataset_use), nrow = 1))
          par(xpd = TRUE)
          
          for (dataset_name in dataset_use) {
            plot_individual_one(
              obj = object.list[[dataset_name]],
              lr_pair = lr_pair,
              dataset_name = dataset_name
            )
          }
          
        } else {
          
          old_par <- par(no.readonly = TRUE)
          on.exit(par(old_par), add = TRUE)
          
          layout(matrix(1, nrow = 1))
          par(xpd = TRUE)
          
          plot_individual_one(
            obj = cellchat_obj,
            lr_pair = lr_pair,
            dataset_name = NULL
          )
        }
        
        dev.off()
        
        generated_files <- c(generated_files, out_file)
        
      }, error = function(e) {
        
        try(dev.off(), silent = TRUE)
        
        warning(
          paste0(
            "Failed to generate plot for ",
            lr_pair,
            ": ",
            e$message
          )
        )
      })
    }
    
    rv$netVisual_individual_files <- generated_files
    
    removeModal()
    
    showNotification(
      paste0(
        "L-R pair network plots generated successfully. ",
        "Generated files: ", length(generated_files),
        ". Output directory: ", out_dir
      ),
      type = "message",
      duration = 8
    )
  })
  
  output$download_netVisual_individual <- downloadHandler(
    
    filename = function() {
      paste0(
        "LRPairNetwork_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".zip"
      )
    },
    
    content = function(file) {
      
      files <- rv$netVisual_individual_files
      
      if (is.null(files) || length(files) == 0) {
        showNotification(
          "No L-R pair network plots generated yet.",
          type = "error",
          duration = 5
        )
        return(NULL)
      }
      
      files <- files[file.exists(files)]
      
      if (length(files) == 0) {
        showNotification(
          "Generated files were not found.",
          type = "error",
          duration = 5
        )
        return(NULL)
      }
      
      old <- getwd()
      on.exit(setwd(old), add = TRUE)
      
      setwd(dirname(files[1]))
      
      zip::zip(
        zipfile = file,
        files = basename(files)
      )
      
      showNotification(
        "L-R pair network plots downloaded successfully.",
        type = "message",
        duration = 3
      )
    }
  )
  # ----------------- VISUALIZATION FUNCTION Cell Cell Communication compareInteractions -----------------
  compareInteractions_reactive <- reactive({
    
    req(rv$cccd_obj)
    req(rv$cccd_obj$data_merge)
    req(input$compareInteractions_font_size)
    req(input$compareInteractions_title_size)
    req(!rv$ui_freeze)
    
    cellchat_merge <- rv$cccd_obj$data_merge
    
    make_compare_plot <- function(measure_use, title_use) {
      
      p <- compareInteractions(
        cellchat_merge,
        show.legend = FALSE,
        group = c(1, 2),
        measure = measure_use
      )
      
      p +
        ggplot2::ggtitle(title_use) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            size = input$compareInteractions_title_size,
            face = "bold",
            hjust = 0.5
          ),
          axis.title = ggplot2::element_text(
            size = input$compareInteractions_font_size
          ),
          axis.text = ggplot2::element_text(
            size = input$compareInteractions_font_size
          ),
          legend.text = ggplot2::element_text(
            size = input$compareInteractions_font_size
          ),
          legend.title = ggplot2::element_text(
            size = input$compareInteractions_font_size
          )
        )
    }
    
    plot_count <- make_compare_plot(
      measure_use = "count",
      title_use = "Number of interactions"
    )
    
    plot_weight <- make_compare_plot(
      measure_use = "weight",
      title_use = "Interaction strength"
    )
    
    plot_fun <- function() {
      print(
        patchwork::wrap_plots(
          plot_weight,         
           plot_count,
          nrow = 1
        )
      )
      invisible(NULL)
    }
    
    list(
      type = "plot",
      plot_fun = plot_fun
    )
  })
  
  output$compareInteractions_control <- renderUI({
    
    div(
      style = "max-width: 800px; overflow-x: auto; width: 100%;",
      plotOutput(
        "compareInteractions_plot",
        height = "300px",
        width = "600px"
      )
    )
  })
  
  
  output$compareInteractions_plot <- renderPlot({
    
    req(compareInteractions_reactive())
    req(!rv$ui_freeze)
    
    res <- compareInteractions_reactive()
    
    if (is.list(res) && res$type == "message") {
      plot.new()
      text(0.5, 0.5, res$content, cex = 1.2, col = "red")
      return()
    }
    
    res$plot_fun()
  })
  observeEvent(input$save_compareInteractions, {
    
    req(compareInteractions_reactive())
    req(input$compareInteractions_save_path)
    
    res <- compareInteractions_reactive()
    
    if (is.list(res) && res$type == "message") {
      showNotification(res$content, type = "error")
      return()
    }
    
    out_dir <- input$compareInteractions_save_path
    
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }
    
    out_file <- file.path(
      out_dir,
      paste0(
        "OverallCommunicationComparison_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".pdf"
      )
    )
    
    pdf(
      out_file,
      width = input$compareInteractions_width,
      height = input$compareInteractions_height
    )
    
    res$plot_fun()
    
    dev.off()
    
    showNotification(
      paste0("Saved overall communication comparison: ", out_file),
      type = "message"
    )
  })
  output$download_compareInteractions <- downloadHandler(
    
    filename = function() {
      paste0(
        "OverallCommunicationComparison_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".pdf"
      )
    },
    
    content = function(file) {
      
      req(compareInteractions_reactive())
      
      res <- compareInteractions_reactive()
      
      pdf(
        file,
        width = input$compareInteractions_width,
        height = input$compareInteractions_height
      )
      
      if (is.list(res) && res$type == "message") {
        plot.new()
        text(0.5, 0.5, res$content, cex = 1.2, col = "red")
      } else {
        res$plot_fun()
      }
      
      dev.off()
      
      showNotification(
        "Overall communication comparison downloaded successfully.",
        type = "message",
        duration = 3
      )
    }
  )
  # ----------------- VISUALIZATION FUNCTION Cell Cell Communication netVisual_diffInteraction -----------------
  cccc_cluster_compare_reactive <- reactive({
    
    req(rv$cccd_obj)
    req(rv$cccd_obj$data_merge)
    req(input$cccc_cluster_vertex_label_size)
    req(input$cccc_cluster_edge_width_max)
    req(input$cccc_cluster_arrow_size)
    req(input$cccc_cluster_title_size)
    req(!rv$ui_freeze)
    
    cellchat_merge <- rv$cccd_obj$data_merge
    
    plot_fun <- function() {
      
      old_par <- par(no.readonly = TRUE)
      on.exit(par(old_par), add = TRUE)
      
      par(mfrow = c(1, 2), xpd = TRUE)
      
      netVisual_diffInteraction(
        cellchat_merge,
        weight.scale = TRUE,
        measure = "weight",
        vertex.label.cex = input$cccc_cluster_vertex_label_size,
        edge.width.max = input$cccc_cluster_edge_width_max,
        arrow.size = input$cccc_cluster_arrow_size,
        title.name = ""
      )
      title(
        main = "Differential interaction strength",
        cex.main = input$cccc_cluster_title_size / 12
      )
      
      netVisual_diffInteraction(
        cellchat_merge,
        weight.scale = TRUE,
        measure = "count",
        vertex.label.cex = input$cccc_cluster_vertex_label_size,
        edge.width.max = input$cccc_cluster_edge_width_max,
        arrow.size = input$cccc_cluster_arrow_size,
        title.name = ""
      )
      title(
        main = "Differential number of interactions",
        cex.main = input$cccc_cluster_title_size / 12
      )
      

      
      invisible(NULL)
    }
    
    list(
      type = "plot",
      plot_fun = plot_fun
    )
  })
  output$cccc_cluster_compare_control <- renderUI({
    
    div(
      style = "max-width: 800px; overflow-x: auto; width: 100%;",
      plotOutput(
        "cccc_cluster_compare_plot",
        height = "450px",
        width = "900px"
      )
    )
  })
  output$cccc_cluster_compare_plot <- renderPlot({
    
    req(cccc_cluster_compare_reactive())
    req(!rv$ui_freeze)
    
    res <- cccc_cluster_compare_reactive()
    
    if (is.list(res) && res$type == "message") {
      plot.new()
      text(0.5, 0.5, res$content, cex = 1.2, col = "red")
      return()
    }
    
    res$plot_fun()
  })
  observeEvent(input$save_cccc_cluster_compare, {
    
    req(cccc_cluster_compare_reactive())
    req(input$cccc_cluster_save_path)
    
    res <- cccc_cluster_compare_reactive()
    
    out_dir <- input$cccc_cluster_save_path
    
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }
    
    out_file <- file.path(
      out_dir,
      paste0(
        "ClusterCommunicationComparison_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".pdf"
      )
    )
    
    pdf(
      out_file,
      width = input$cccc_cluster_width,
      height = input$cccc_cluster_height
    )
    
    res$plot_fun()
    
    dev.off()
    
    showNotification(
      paste0("Saved cluster communication comparison: ", out_file),
      type = "message"
    )
  })
  output$download_cccc_cluster_compare <- downloadHandler(
    
    filename = function() {
      paste0(
        "ClusterCommunicationComparison_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".pdf"
      )
    },
    
    content = function(file) {
      
      req(cccc_cluster_compare_reactive())
      
      res <- cccc_cluster_compare_reactive()
      
      pdf(
        file,
        width = input$cccc_cluster_width,
        height = input$cccc_cluster_height
      )
      
      res$plot_fun()
      
      dev.off()
      
      showNotification(
        "Cluster communication comparison downloaded successfully.",
        type = "message",
        duration = 3
      )
    }
  )
  # ----------------- VISUALIZATION FUNCTION Cell Cell Communication netVisual_heatmap -----------------
  netVisual_relative_heatmap_reactive <- reactive({
    
    req(rv$cccd_obj)
    req(rv$cccd_obj$data_merge)
    req(input$netVisual_relative_heatmap_axis_title_size)
    req(input$netVisual_relative_heatmap_axis_text_size)
    req(input$netVisual_relative_heatmap_legend_size)
    req(!rv$ui_freeze)
    
    cellchat_merge <- rv$cccd_obj$data_merge
    
    make_relative_heatmap <- function(measure_use) {
      
      title_use <- if (measure_use == "count") {
        "Differential number of interactions"
      } else {
        "Differential interaction strength"
      }
      
      ht <- netVisual_heatmap(
        cellchat_merge,
        comparison = c(1, 2),
        measure = measure_use,
        slot.name = "net",
        title.name = title_use,
        font.size = input$netVisual_relative_heatmap_axis_text_size
      )
      
      ht@row_names_param$gp <- grid::gpar(
        fontsize = input$netVisual_relative_heatmap_axis_text_size
      )
      ht@column_names_param$gp <- grid::gpar(
        fontsize = input$netVisual_relative_heatmap_axis_text_size
      )
      ht@row_title_param$gp <- grid::gpar(
        fontsize = input$netVisual_relative_heatmap_axis_title_size
      )
      ht@column_title_param$gp <- grid::gpar(
        fontsize = input$netVisual_relative_heatmap_axis_title_size
      )
      ht@matrix_legend_param$title_gp <- grid::gpar(
        fontsize = input$netVisual_relative_heatmap_legend_size
      )
      ht@matrix_legend_param$labels_gp <- grid::gpar(
        fontsize = input$netVisual_relative_heatmap_legend_size
      )
      
      ht
    }
    
    draw_heatmap_to_grob <- function(ht) {
      grid::grid.grabExpr({
        ComplexHeatmap::draw(ht)
      })
    }
    
    grob_list <- list(
      draw_heatmap_to_grob(make_relative_heatmap("weight")),
      draw_heatmap_to_grob(make_relative_heatmap("count"))

    )
    
    plot_fun <- function() {
      combined_grob <- gridExtra::arrangeGrob(
        grobs = grob_list,
        nrow = 1,
        ncol = 2
      )
      grid::grid.newpage()
      grid::grid.draw(combined_grob)
      invisible(NULL)
    }
    
    list(type = "plot", plot_fun = plot_fun)
  })
  output$netVisual_relative_heatmap_control <- renderUI({
    
    get_cluster_n <- function(obj) {
      if (is.null(obj)) return(10)
      if (methods::is(obj, "CellChat")) {
        if (!is.null(obj@idents)) return(length(levels(obj@idents)))
        if (!is.null(obj@net$weight)) return(ncol(obj@net$weight))
      }
      return(10)
    }
    
    cluster_n <- get_cluster_n(rv$cccd_obj$data_merge)
    single_panel_size <- max(350, cluster_n * 45)
    
    div(
      style = "max-width: 800px; max-height: 1600px; overflow-x: auto; overflow-y: auto; width: 100%;",
      plotOutput(
        "netVisual_relative_heatmap_plot",
        height = paste0(single_panel_size, "px"),
        width = paste0(single_panel_size * 2, "px")
      )
    )
  })
  output$netVisual_relative_heatmap_plot <- renderPlot({
    
    req(netVisual_relative_heatmap_reactive())
    req(!rv$ui_freeze)
    
    res <- netVisual_relative_heatmap_reactive()
    
    if (is.list(res) && res$type == "message") {
      plot.new()
      text(0.5, 0.5, res$content, cex = 1.2, col = "red")
      return()
    }
    
    res$plot_fun()
  })
  observeEvent(input$save_netVisual_relative_heatmap, {
    
    req(netVisual_relative_heatmap_reactive())
    req(input$netVisual_relative_heatmap_save_path)
    
    res <- netVisual_relative_heatmap_reactive()
    
    out_dir <- input$netVisual_relative_heatmap_save_path
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    
    out_file <- file.path(
      out_dir,
      paste0(
        "RelativeCommunicationHeatmap_Count_Weight_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".pdf"
      )
    )
    
    pdf(
      out_file,
      width = input$netVisual_relative_heatmap_width,
      height = input$netVisual_relative_heatmap_height
    )
    
    res$plot_fun()
    dev.off()
    
    showNotification(
      paste0("Saved relative communication heatmap: ", out_file),
      type = "message"
    )
  })
  output$download_netVisual_relative_heatmap <- downloadHandler(
    
    filename = function() {
      paste0(
        "RelativeCommunicationHeatmap_Count_Weight_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".pdf"
      )
    },
    
    content = function(file) {
      
      req(netVisual_relative_heatmap_reactive())
      
      res <- netVisual_relative_heatmap_reactive()
      
      pdf(
        file,
        width = input$netVisual_relative_heatmap_width,
        height = input$netVisual_relative_heatmap_height
      )
      
      res$plot_fun()
      dev.off()
      
      showNotification(
        "Relative communication heatmap downloaded successfully.",
        type = "message",
        duration = 3
      )
    }
  )
  # ----------------- VISUALIZATION FUNCTION Cell Cell Communication relatived netVisual_diffInteraction -----------------
  plot_diffInteraction_pathway_message <- function(msg) {
    plot.new()
    text(0.5, 0.5, msg, col = "red", cex = 1.2)
    invisible(NULL)
  }
  netVisual_diffInteraction_pathway <- function(
    control_obj,
    treatment_obj,
    signaling,
    measure = c("weight", "count"),
    thresh = 0.05,
    color.use = NULL,
    color.edge = c("#b2182b", "#2166ac"),
    title.name = NULL,
    sources.use = NULL,
    targets.use = NULL,
    remove.isolate = FALSE,
    top = 1,
    weight.scale = TRUE,
    vertex.weight = NULL,
    vertex.weight.max = NULL,
    vertex.size.max = 15,
    vertex.label.cex = 1,
    edge.weight.max = NULL,
    edge.width.max = 8,
    arrow.size = 0.5,
    title.size = 12
  ) {
    
    measure <- match.arg(measure)
    
    get_cellchat_palette <- function(n) {
      if (exists("ggPalette", envir = asNamespace("CellChat"), inherits = FALSE)) {
        get("ggPalette", envir = asNamespace("CellChat"))(n)
      } else {
        grDevices::rainbow(n)
      }
    }
    
    get_group_colors <- function(obj, cells, color.use = NULL) {
      
      if (!is.null(color.use)) {
        if (is.null(names(color.use))) {
          color.use <- color.use[seq_along(cells)]
          names(color.use) <- cells
        }
        return(color.use[cells])
      }
      
      if (!is.null(obj@idents)) {
        cell.levels <- levels(obj@idents)
        color.use <- get_cellchat_palette(length(cell.levels))
        names(color.use) <- cell.levels
        return(color.use[cells])
      }
      
      color.use <- get_cellchat_palette(length(cells))
      names(color.use) <- cells
      color.use[cells]
    }
    
    get_weight_mat <- function(obj, signaling) {
      
      if (is.null(obj@netP$prob) || is.null(obj@netP$pathways)) {
        return(NULL)
      }
      
      pathway_idx <- match(signaling, obj@netP$pathways)
      
      if (is.na(pathway_idx)) {
        return(NULL)
      }
      
      mat <- obj@netP$prob[, , pathway_idx, drop = FALSE]
      mat <- mat[, , 1]
      mat[is.na(mat)] <- 0
      
      mat
    }
    
    get_count_mat <- function(obj, signaling, thresh = 0.05) {
      
      if (is.null(obj@net$prob) || is.null(obj@net$pval)) {
        return(NULL)
      }
      
      if (is.null(obj@LR$LRsig)) {
        return(NULL)
      }
      
      lr_info <- obj@LR$LRsig
      
      if (!("pathway_name" %in% colnames(lr_info))) {
        return(NULL)
      }
      
      lr_use <- lr_info[lr_info$pathway_name == signaling, , drop = FALSE]
      
      if (nrow(lr_use) == 0) {
        return(NULL)
      }
      
      if ("interaction_name" %in% colnames(lr_use) &&
          !is.null(dimnames(obj@net$prob)[[3]])) {
        
        lr_idx <- match(lr_use$interaction_name, dimnames(obj@net$prob)[[3]])
        lr_idx <- lr_idx[!is.na(lr_idx)]
        
      } else {
        
        lr_idx <- which(lr_info$pathway_name == signaling)
      }
      
      if (length(lr_idx) == 0) {
        return(NULL)
      }
      
      prob_array <- obj@net$prob[, , lr_idx, drop = FALSE]
      pval_array <- obj@net$pval[, , lr_idx, drop = FALSE]
      
      count_mat <- apply(
        (prob_array > 0) & (pval_array < thresh),
        c(1, 2),
        sum,
        na.rm = TRUE
      )
      
      count_mat[is.na(count_mat)] <- 0
      count_mat
    }
    
    resolve_cells <- function(x, cells) {
      if (is.null(x)) return(cells)
      if (is.numeric(x)) return(cells[x])
      intersect(x, cells)
    }
    
    if (measure == "weight") {
      control_mat <- get_weight_mat(control_obj, signaling)
      treatment_mat <- get_weight_mat(treatment_obj, signaling)
    } else {
      control_mat <- get_count_mat(control_obj, signaling, thresh)
      treatment_mat <- get_count_mat(treatment_obj, signaling, thresh)
    }
    if (is.null(control_mat) || is.null(treatment_mat)) {
      plot.new()
      text(
        0.5, 0.5,
        paste0("Pathway '", signaling, "' is not found in CONTROL or TREATMENT."),
        col = "red",
        cex = 1.2
      )
      return(invisible(NULL))
    }
    
    common_cells <- intersect(rownames(control_mat), rownames(treatment_mat))
    
    control_mat <- control_mat[common_cells, common_cells, drop = FALSE]
    treatment_mat <- treatment_mat[common_cells, common_cells, drop = FALSE]
    
    net.diff <- treatment_mat - control_mat
    net.diff[is.na(net.diff)] <- 0
    
    if (!is.null(sources.use)) {
      source.cells <- resolve_cells(sources.use, common_cells)
      net.diff[setdiff(common_cells, source.cells), ] <- 0
    }
    
    if (!is.null(targets.use)) {
      target.cells <- resolve_cells(targets.use, common_cells)
      net.diff[, setdiff(common_cells, target.cells)] <- 0
    }
    
    if (top < 1) {
      edge_values <- abs(as.vector(net.diff))
      edge_values <- edge_values[edge_values > 0]
      
      if (length(edge_values) > 0) {
        cutoff <- stats::quantile(edge_values, probs = 1 - top)
        net.diff[abs(net.diff) < cutoff] <- 0
      }
    }
    
    if (remove.isolate) {
      keep <- rowSums(abs(net.diff)) + colSums(abs(net.diff)) > 0
      net.diff <- net.diff[keep, keep, drop = FALSE]
      common_cells <- rownames(net.diff)
    }
    
    if (length(common_cells) == 0 || all(net.diff == 0)) {
      plot.new()
      text(
        0.5, 0.5,
        paste0("No differential interaction detected for pathway: ", signaling),
        col = "red",
        cex = 1.2
      )
      return(invisible(NULL))
    }
    
    net.abs <- abs(net.diff)
    
    g <- igraph::graph_from_adjacency_matrix(
      net.abs,
      mode = "directed",
      weighted = TRUE,
      diag = TRUE
    )
    
    edge_df <- igraph::as_data_frame(g, what = "edges")
    
    edge.sign <- mapply(
      function(from, to) {
        net.diff[from, to]
      },
      edge_df$from,
      edge_df$to
    )
    
    igraph::E(g)$color <- ifelse(
      edge.sign > 0,
      color.edge[1],
      color.edge[2]
    )
    
    if (is.null(edge.weight.max)) {
      edge.weight.max <- max(igraph::E(g)$weight)
    }
    
    if (weight.scale) {
      igraph::E(g)$width <- igraph::E(g)$weight / edge.weight.max * edge.width.max
    } else {
      igraph::E(g)$width <- igraph::E(g)$weight
    }
    
    if (is.null(vertex.weight)) {
      
      group.size <- as.numeric(table(control_obj@idents))
      names(group.size) <- names(table(control_obj@idents))
      
      vertex.weight <- group.size[common_cells]
      vertex.weight[is.na(vertex.weight)] <- 1
    }
    
    if (is.null(vertex.weight.max)) {
      vertex.weight.max <- max(vertex.weight)
    }
    
    vertex.size <- vertex.weight / vertex.weight.max * vertex.size.max
    vertex.size[vertex.size < 3] <- 3
    
    
    vertex.color <- get_group_colors(
      control_obj,
      common_cells,
      color.use = color.use
    )
    
    layout.circle <- igraph::layout_in_circle(g)
    
    plot(
      g,
      layout = layout.circle,
      vertex.label = common_cells,
      vertex.label.cex = vertex.label.cex,
      vertex.label.color = "black",
      vertex.size = vertex.size,
      vertex.color = vertex.color,
      vertex.frame.color = NA,
      edge.color = igraph::E(g)$color,
      edge.width = igraph::E(g)$width,
      edge.arrow.size = arrow.size * 0.8,
      edge.curved = 0.2,
      margin = 0.15,
      main = title.name,
      cex.main = title.size / 12
    )
    
    invisible(NULL)
  }
  diffInteraction_pathway_reactive <- reactive({
    
    req(rv$cccd_obj)
    req(rv$cccd_obj$data)
    req(rv$ccc_pathway_select)
    req(input$diffInteraction_pathway_vertex_label_size)
    req(input$diffInteraction_pathway_edge_width_max)
    req(input$diffInteraction_pathway_arrow_size)
    req(input$diffInteraction_pathway_title_size)
    req(!rv$ui_freeze)
    
    cellchat_list <- rv$cccd_obj$data
    
    control_obj <- cellchat_list[["CONTROL"]]
    treatment_obj <- cellchat_list[["TREATMENT"]]
    
    signaling_use <- trimws(rv$ccc_pathway_select)
    
    plot_fun <- function() {
      
      old_par <- par(no.readonly = TRUE)
      on.exit(par(old_par), add = TRUE)
      
      par(mfrow = c(1, 2), xpd = TRUE)
      
      netVisual_diffInteraction_pathway(
        control_obj = control_obj,
        treatment_obj = treatment_obj,
        signaling = signaling_use,
        measure = "weight",
        weight.scale = TRUE,
        vertex.label.cex = input$diffInteraction_pathway_vertex_label_size,
        edge.width.max = input$diffInteraction_pathway_edge_width_max,
        arrow.size = input$diffInteraction_pathway_arrow_size,
        title.size = input$diffInteraction_pathway_title_size,
        title.name = "Differential interaction strength"
      )
      
      netVisual_diffInteraction_pathway(
        control_obj = control_obj,
        treatment_obj = treatment_obj,
        signaling = signaling_use,
        measure = "count",
        weight.scale = TRUE,
        vertex.label.cex = input$diffInteraction_pathway_vertex_label_size,
        edge.width.max = input$diffInteraction_pathway_edge_width_max,
        arrow.size = input$diffInteraction_pathway_arrow_size,
        title.size = input$diffInteraction_pathway_title_size,
        title.name = "Differential number of interactions"
      )
      
      invisible(NULL)
    }
    
    list(
      type = "plot",
      plot_fun = plot_fun
    )
  })
  
  output$diffInteraction_pathway_control <- renderUI({
    
    div(
      style = "max-width: 800px; overflow-x: auto; width: 100%;",
      plotOutput(
        "diffInteraction_pathway_plot",
        height = "400px",
        width = "800px"
      )
    )
  })
  
  output$diffInteraction_pathway_plot <- renderPlot({
    
    req(diffInteraction_pathway_reactive())
    req(!rv$ui_freeze)
    
    res <- diffInteraction_pathway_reactive()
    
    if (is.list(res) && res$type == "message") {
      plot_diffInteraction_pathway_message(res$content)
      return()
    }
    
    res$plot_fun()
  })
  
  observeEvent(input$save_diffInteraction_pathway, {
    
    req(diffInteraction_pathway_reactive())
    req(input$diffInteraction_pathway_save_path)
    
    res <- diffInteraction_pathway_reactive()
    out_dir <- input$diffInteraction_pathway_save_path
    
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }
    
    pathway_safe <- gsub("[^A-Za-z0-9_\\-]", "_", rv$ccc_pathway_select)
    
    out_file <- file.path(
      out_dir,
      paste0(
        "DifferentialInteractionPathway_",
        pathway_safe,
        "_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".pdf"
      )
    )
    
    pdf(
      out_file,
      width = input$diffInteraction_pathway_width,
      height = input$diffInteraction_pathway_height
    )
    
    res$plot_fun()
    
    dev.off()
    
    showNotification(
      paste0("Saved differential pathway interaction network: ", out_file),
      type = "message"
    )
  })
  
  output$download_diffInteraction_pathway <- downloadHandler(
    
    filename = function() {
      
      pathway_safe <- gsub("[^A-Za-z0-9_\\-]", "_", rv$ccc_pathway_select)
      
      paste0(
        "DifferentialInteractionPathway_",
        pathway_safe,
        "_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".pdf"
      )
    },
    
    content = function(file) {
      
      req(diffInteraction_pathway_reactive())
      
      res <- diffInteraction_pathway_reactive()
      
      pdf(
        file,
        width = input$diffInteraction_pathway_width,
        height = input$diffInteraction_pathway_height
      )
      
      res$plot_fun()
      
      dev.off()
      
      showNotification(
        "Differential pathway interaction network downloaded successfully.",
        type = "message",
        duration = 3
      )
    }
  )
  
  # ----------------- VISUALIZATION FUNCTION Cell Cell Communication rankNet -----------------
  rankNet_reactive <- reactive({
    
    req(rv$cccd_obj)
    req(rv$cccd_obj$data_merge)
    req(input$rankNet_font_size)
    req(input$rankNet_title_size)
    req(!rv$ui_freeze)
    
    cellchat_merge <- rv$cccd_obj$data_merge
    
    plot_fun <- function() {
      
      p1 <- rankNet(
        cellchat_merge,
        mode = "comparison",
        show.raw = TRUE,
        stacked = TRUE,
        do.stat = TRUE
      ) +
        ggplot2::ggtitle("Stacked") +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            size = input$rankNet_title_size,
            face = "bold",
            hjust = 0.5
          ),
          axis.title = ggplot2::element_text(
            size = input$rankNet_font_size
          ),
          axis.text = ggplot2::element_text(
            size = input$rankNet_font_size
          ),
          legend.title = ggplot2::element_text(
            size = input$rankNet_font_size
          ),
          legend.text = ggplot2::element_text(
            size = input$rankNet_font_size
          )
        )
      
      p2 <- rankNet(
        cellchat_merge,
        mode = "comparison",
        show.raw = TRUE,
        stacked = FALSE,
        do.stat = TRUE
      ) +
        ggplot2::ggtitle("Unstacked") +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            size = input$rankNet_title_size,
            face = "bold",
            hjust = 0.5
          ),
          axis.title = ggplot2::element_text(
            size = input$rankNet_font_size
          ),
          axis.text = ggplot2::element_text(
            size = input$rankNet_font_size
          ),
          legend.title = ggplot2::element_text(
            size = input$rankNet_font_size
          ),
          legend.text = ggplot2::element_text(
            size = input$rankNet_font_size
          )
        )
      
      print(
        p1 + p2 +
          patchwork::plot_layout(ncol = 2)
      )
      
      invisible(NULL)
    }
    
    list(
      type = "plot",
      plot_fun = plot_fun
    )
  })
  
  output$rankNet_control <- renderUI({
    
    pathway_n <- length(unique(rv$cccd_obj$data_merge@netP$pathways))
    
    height_val <- max(600, pathway_n * 18)
    
    div(
      style = "max-width: 800px; max-width: 1000px; overflow-x: auto; overflow-y: auto; width: 100%;",
      plotOutput(
        "rankNet_plot",
        height = paste0(height_val, "px"),
        width = "800px"
      )
    )
  })

  
  
  output$rankNet_plot <- renderPlot({
    
    req(rankNet_reactive())
    req(!rv$ui_freeze)
    
    res <- rankNet_reactive()
    
    if (is.list(res) && res$type == "message") {
      plot.new()
      text(0.5, 0.5, res$content, cex = 1.2, col = "red")
      return()
    }
    
    res$plot_fun()
  })
  observeEvent(input$save_rankNet, {
    
    req(rankNet_reactive())
    req(input$rankNet_save_path)
    
    res <- rankNet_reactive()
    
    out_dir <- input$rankNet_save_path
    
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }
    
    out_file <- file.path(
      out_dir,
      paste0(
        "RankNet_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".pdf"
      )
    )
    
    pdf(
      out_file,
      width = input$rankNet_width,
      height = input$rankNet_height
    )
    
    res$plot_fun()
    
    dev.off()
    
    showNotification(
      paste0("Saved rankNet plot: ", out_file),
      type = "message"
    )
  })
  output$download_rankNet <- downloadHandler(
    
    filename = function() {
      paste0(
        "RankNet_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".pdf"
      )
    },
    
    content = function(file) {
      
      req(rankNet_reactive())
      
      res <- rankNet_reactive()
      
      pdf(
        file,
        width = input$rankNet_width,
        height = input$rankNet_height
      )
      
      res$plot_fun()
      
      dev.off()
      
      showNotification(
        "rankNet plot downloaded successfully.",
        type = "message",
        duration = 3
      )
    }
  )
  # ----------------- Data Output Gene expression -----------------
  output$gene_expression_sub_ui <- renderUI({
    
    if (is.null(rv$data_obj)) {
      
      tagList(
        div(
          style = "padding: 15px; background-color: #fff3cd; border: 1px solid #ffeeba; border-radius: 5px;",
          "Gene expression data output requires a loaded Seurat Obj data. Please load the data to proceed."
        )
      )
      
    } else {
      meta_cols <- colnames(rv$data_obj@meta.data)
      
      meta_df <- rv$data_obj@meta.data
      valid_cols <- meta_cols[
        sapply(meta_df[, meta_cols, drop = FALSE], function(x) {
          !is.numeric(x)
        })
      ]
      
      valid_cols <- c("", setdiff(valid_cols, c("CB", "CB_original")))
      tagList(
        fluidRow(
          column(6, selectInput("gene_expression_first_meta", "First layer meta.data field", choices = valid_cols, width = "100%")),
          column(6, selectInput("gene_expression_second_meta", "Second layer meta.data field", choices = valid_cols, width = "100%"))
        ),
        fluidRow(
          column(6, actionButton("reorder_first_meta", "Reorder First meta.data field", class = "btn-warning", width = "100%")),
          column(6, actionButton("reorder_second_meta", "Reorder Second meta.data field", class = "btn-warning", width = "100%"))
        ),
        br(),
        radioButtons("gene_expression_data_slot", "Data slot for visualization",
                     choices = c("Normalized UMI (Recommended)" = "data",
                                 "Raw UMI per cell (counts)" = "counts"),
                     selected = "data", width = "100%", inline = TRUE),
        actionButton("marker_select_btn", "Select genes to be labeled (All genes will be analyzed by default)", class = "btn-info", width = "100%"),
        actionButton("generate_gene_expression", "Generate and Preview Data", class = "btn-success", width = "100%"),
      br(),
      uiOutput("preview_gene_expression")
      )
    }
  })
  
  observeEvent(input$reorder_first_meta, {
    
    req(input$gene_expression_first_meta)
    
    meta_name <- input$gene_expression_first_meta
    
    req(meta_name != "")
    
    x <- rv$data_obj@meta.data[[meta_name]]
    
    current_levels <- if (is.factor(x)) {
      levels(x)
    } else {
      unique(as.character(x))
    }
    
    showModal(
      modalDialog(
        
        title = paste("Reorder:", meta_name),
        
        rank_list(
          text = "Drag to reorder",
          labels = current_levels,
          input_id = "first_meta_order"
        ),
        
        footer = tagList(
          modalButton("Cancel"),
          actionButton(
            "apply_first_meta_order",
            "Apply",
            class = "btn-success"
          )
        ),
        
        easyClose = TRUE,
        size = "m"
      )
    )
  })
  
  observeEvent(input$apply_first_meta_order, {
    
    req(input$gene_expression_first_meta)
    req(input$first_meta_order)
    
    meta_name <- input$gene_expression_first_meta
    
    new_order <- input$first_meta_order
    
    rv$data_obj@meta.data[[meta_name]] <- factor(
      rv$data_obj@meta.data[[meta_name]],
      levels = new_order
    )
    
    updateSelectInput(
      session,
      "gene_expression_first_meta",
      selected = meta_name
    )
    meta.name2 = input$gene_expression_second_meta
    updateSelectInput(
      session,
      "gene_expression_second_meta",
      selected = meta.name2
    )
    
    removeModal()
    
    showNotification(
      paste("Updated order for", meta_name),
      type = "message",
      duration = 3
    )
  })
  
  observeEvent(input$reorder_second_meta, {
    
    req(input$gene_expression_second_meta)
    
    meta_name <- input$gene_expression_second_meta
    
    req(meta_name != "")
    
    x <- rv$data_obj@meta.data[[meta_name]]
    
    current_levels <- if (is.factor(x)) {
      levels(x)
    } else {
      unique(as.character(x))
    }
    
    showModal(
      modalDialog(
        
        title = paste("Reorder:", meta_name),
        
        rank_list(
          text = "Drag to reorder",
          labels = current_levels,
          input_id = "second_meta_order"
        ),
        
        footer = tagList(
          modalButton("Cancel"),
          actionButton(
            "apply_second_meta_order",
            "Apply",
            class = "btn-success"
          )
        ),
        
        easyClose = TRUE,
        size = "m"
      )
    )
  })
  
  observeEvent(input$apply_second_meta_order, {
    
    req(input$gene_expression_second_meta)
    req(input$second_meta_order)
    
    meta_name <- input$gene_expression_second_meta
    
    new_order <- input$second_meta_order
    
    rv$data_obj@meta.data[[meta_name]] <- factor(
      rv$data_obj@meta.data[[meta_name]],
      levels = new_order
    )
    meta_name1 <- input$gene_expression_first_meta
    updateSelectInput(
      session,
      "gene_expression_first_meta",
      selected = meta_name1
    )
    
    updateSelectInput(
      session,
      "gene_expression_second_meta",
      selected = meta_name
    )
    
    removeModal()
    
    showNotification(
      paste("Updated order for", meta_name),
      type = "message",
      duration = 3
    )
  })
  
  observeEvent(input$generate_gene_expression, {
    if(input$gene_expression_first_meta == ""){
      showModal(
        modalDialog(
          title = "Error",
          p("You must select the First layer meta.data field to analysis."),
          easyClose = TRUE,
          footer = actionButton("close_modal", "Close")
        )
      )
      return(NULL)
    }
    source("gene_expression_output.R")
    
    showModal(
      modalDialog(
        title = "Generate Data",
        p("Generate the gene expression result"),
        easyClose = FALSE,
        footer = NULL
      )
    )
    
    ana_genes <- unique(trimws(unlist(strsplit(rv$visualization_marker, "\n"))))
    ana_genes <- ana_genes[ana_genes != ""]
    if (is.null(ana_genes) || length(ana_genes) == 0 || all(is.na(ana_genes))) {
      ana_genes <- rownames(rv$data_obj)
    }
    res <- create_expression_summary(
      seurat_obj = rv$data_obj,
      genes = ana_genes,
      first_meta = input$gene_expression_first_meta,
      second_meta = input$gene_expression_second_meta,
      layer = input$gene_expression_data_slot
    )
    rv$gene_expression_result <- res$export
    removeModal()
  })
  
  output$preview_gene_expression <- renderUI({
    req(rv$gene_expression_result)
    
    tagList(
      br(),
      h4("Preview:"),
      div(
        style = "overflow-x: auto; width: 100%;",
        tableOutput("gene_expression_table")
      ),
      br(),
      downloadButton("download_gene_expression", "Download result", style = "background-color: #4CAF50; color: white; width: 100%;")
    )
  })
  
  output$gene_expression_table <- renderTable({
    req(rv$gene_expression_result)
    head(rv$gene_expression_result, 10)
  }, rownames = FALSE)
  
  output$download_gene_expression <- downloadHandler(
    
    filename = function() {
      paste0(
        "Gene_expression_",
        tools::file_path_sans_ext(basename(rv$data_obj_path)),
        "_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".xlsx"
      )
    },
    
    content = function(file) {
      
      res <- rv$gene_expression_result
      
      if (is.null(res) || nrow(res) == 0) {
        stop("gene_expression_result is empty")
      }
      
      df <- as.data.frame(res, stringsAsFactors = FALSE)
      
      x_row <- which(df[[1]] == "value")[1] + 1
      
      if (is.na(x_row)) {
        stop("Cannot find 'value' in first column")
      }
      
      if (x_row < nrow(df)) {
        rows_to_convert <- (x_row + 1):nrow(df)
        cols_to_convert <- 2:ncol(df)
        
        for (j in cols_to_convert) {
          df[rows_to_convert, j] <- as.numeric(as.character(df[rows_to_convert, j]))
        }
      }
      
      wb <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb, "Gene_expression")
      header_rows <- 1:x_row
      openxlsx::writeData(wb, "Gene_expression", df[header_rows, , drop = FALSE])
      
      if (x_row < nrow(df)) {
        data_rows <- (x_row + 1):nrow(df)
        data_df <- df[data_rows, , drop = FALSE]
        
        numeric_cols <- 2:ncol(data_df)
        for (j in numeric_cols) {
          data_df[[j]] <- as.numeric(as.character(data_df[[j]]))
        }
        
        openxlsx::writeData(
          wb, 
          "Gene_expression", 
          data_df,
          startRow = x_row + 1,
          startCol = 1,
          colNames = FALSE
        )
      }
      
      num_style <- openxlsx::createStyle(numFmt = "0.000")
      
      if (x_row < nrow(df)) {
        openxlsx::addStyle(
          wb,
          sheet = "Gene_expression",
          style = num_style,
          rows = (x_row + 1):nrow(df),
          cols = 2:ncol(df),
          gridExpand = TRUE,
          stack = TRUE
        )
      }
      
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
      
      showNotification(
        "Gene expression data downloaded successfully.",
        type = "message",
        duration = 3
      )
    }
  )
  
  
}

shinyApp(ui, server)