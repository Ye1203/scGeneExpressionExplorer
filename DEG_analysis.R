library(shiny)
library(shinyjs)
library(Seurat)
library(openxlsx)
# source(file.path(getwd(), "DE_analysis.R"))
# ============================================================
# -------------------- Utility Functions ---------------------
# ============================================================

load_data_object <- function(path, type = c("seurat", "deg", "gsea")) {
  type <- match.arg(type)
  validate(need(file.exists(path), paste(type, "file does not exist")))
  
  if (type == "seurat") {
    readRDS(path)
  } else if (type == "deg") {
    read.delim(path, header = TRUE, stringsAsFactors = FALSE)
  } else if (type == "gsea") {
    openxlsx::read.xlsx(path)
  }
}

extract_meta_columns <- function(data) {
  if (!"meta.data" %in% slotNames(data)) return(NULL)
  c("",colnames(data@meta.data))
}

extract_unique_values <- function(data, column_name) {
  if (is.null(data) || is.null(column_name)) return(NULL)
  unique(as.character(data@meta.data[[column_name]]))
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
      comps <- comps[sapply(comps, function(x) !is.null(x$numerator) && !is.null(x$denominator))]
      
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
      comps <- comps[sapply(comps, function(x) !is.null(x$numerator) && !is.null(x$denominator))]
      
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
  useShinyjs(),
  titlePanel("DEG/GSEA prealpha"),
  
  tabsetPanel(
    id = "main_tabs",
    
    tabPanel("DEG", uiOutput("deg_ui")),
    tabPanel("GSEA", uiOutput("gsea_ui")),
    
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
              tabPanel("Original Volcano", value = "volcano"),
              tabPanel("Pathway Volcano", value = "pathway_volcano"),
              tabPanel("GSEA NES", value = "gsea_nes")
            )
          )
        ),
        column(
          width = 10,
          uiOutput("viz_panel")
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
  "))
)

# ============================================================
# ------------------------ SERVER ----------------------------
# ============================================================

server <- function(input, output, session) {
  
  # Reactive value
  rv <- reactiveValues(
    data_obj = NULL,
    data_path = NULL,
    
    deg_obj = NULL,
    deg_path = NULL,
    
    gsea_obj = NULL,
    gsea_path = NULL,
    
    visualization_marker = ""
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
      )
    } else {
      fluidRow(
        column(
          9,
          textInput(
            inputId = input_path_id,
            value = "",
            placeholder = placeholder,
            label = NULL,
            width = "100%"
          )
        ),
        column(
          3,
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
    showModal(modalDialog(title = "Loading Data", "Please wait...", footer = NULL, easyClose = FALSE))
    tryCatch({
      obj <- load_data_object(input$data_path, type = "seurat")
      rv$data_obj <- obj
      rv$data_obj_path <- input$data_path
      removeModal()
      showNotification("Seurat data loaded successfully!", type = "message")
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # DEG data loading
  observeEvent(input$load_deg, {
    showModal(modalDialog(title = "Loading DEG Data", "Please wait...", footer = NULL, easyClose = FALSE))
    tryCatch({
      obj <- load_data_object(input$deg_path, type = "deg")
      rv$deg_obj <- obj
      rv$deg_obj_path <- input$deg_path
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
      removeModal()
      showNotification("GSEA data loaded successfully!", type = "message")
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # ----------------- Reactive unique values -----------------
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
              column(6, selectInput("sample_column", "meta.data sample column", 
                                    choices = extract_meta_columns(rv$data_obj))),
              column(6, selectInput("cluster_column", "meta.data cluster column", 
                                    choices = extract_meta_columns(rv$data_obj)))
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
        "Positional (C1) — Gene sets grouped by chromosomal location." = "Positional",
        "CGP (C2, CGP) — Chemical and Genetic Perturbation gene sets." = "CGP",
        "CP (C2, CP) — Canonical pathways." = "CP",
        "BioCarta (C2, CP:BIOCARTA) — Curated signaling pathways." = "BioCarta",
        "KEGG_LEGACY (C2, CP:KEGG_LEGACY) — KEGG legacy pathways." = "KEGG_LEGACY",
        "KEGG (C2, CP:KEGG_MEDICUS) — Curated metabolic and signaling pathways." = "KEGG",
        "PID (C2, CP:PID) — Pathway Interaction Database." = "PID",
        "Reactome (C2, CP:REACTOME) — Expert-curated pathways." = "Reactome",
        "WikiPathways (C2, CP:WIKIPATHWAYS) — Curated signaling pathways." = "WikiPathways",
        "CP_CUSTOM (C2) — Custom curated canonical pathways." = "CP_CUSTOM",
        "Motif_miR (C3, MIR:MIRDB) — microRNA target gene sets." = "Motif_miR",
        "Motif_miR_Legacy (C3, MIR:MIR_LEGACY) — Legacy microRNA target gene sets." = "Motif_miR_Legacy",
        "Motif_TF (C3, TFT:GTRD) — Transcription factor target gene sets." = "Motif_TF",
        "Motif_TF_Legacy (C3, TFT:TFT_LEGACY) — Legacy transcription factor target sets." = "Motif_TF_Legacy",
        "Computational (C4) — Computational gene sets." = "Computational",
        "GO_BP (C5, GO:BP) — Biological Process ontology." = "GO_BP",
        "GO_CC (C5, GO:CC) — Cellular Component ontology." = "GO_CC",
        "GO_MF (C5, GO:MF) — Molecular Function ontology." = "GO_MF",
        "HPO (C5, HPO) — Human Phenotype Ontology." = "HPO",
        "Oncogenic (C6) — Oncogenic signatures." = "Oncogenic",
        "Immune (C7, IMMUNESIGDB) — Immunologic gene sets." = "Immune",
        "VAX (C7, VAX) — Vaccine response gene sets." = "VAX",
        "Hallmark (H) — Hallmark gene sets." = "Hallmark"
      )
    } 
    else {db_choices <- c(
      "Positional (M1) — Gene sets grouped by genomic location." = "Positional",
      "CGP (M2, CGP) — Chemical and Genetic Perturbation gene sets." = "CGP",
      "BioCarta (M2, CP:BIOCARTA) — Classical signaling pathways." = "BioCarta",
      "Reactome (M2, CP:REACTOME) — Curated pathway database." = "Reactome",
      "WikiPathways (M2, CP:WIKIPATHWAYS) — Curated signaling pathways." = "WikiPathways",
      "Motif_TF (M3, GTRD) — Transcription factor target gene sets." = "Motif_TF",
      "Motif_miR (M3, MIRDB) — microRNA target predictions." = "Motif_miR",
      "GO_BP (M5, GO:BP) — Biological Process ontology." = "GO_BP",
      "GO_CC (M5, GO:CC) — Cellular Component ontology." = "GO_CC",
      "GO_MF (M5, GO:MF) — Molecular Function ontology." = "GO_MF",
      "MP Tumor (M5, MPT) — Mouse phenotype tumor gene sets." = "MP_Tumor",
      "Immune (M7) — Immunologic gene sets." = "Immune",
      "Hallmark (MH) — Mouse hallmark gene sets." = "Hallmark"
    )
    }
    
    updateSelectInput(session, "gsea_db", choices = db_choices, selected = c("GO_BP", "Reactome", "Hallmark"))
    
  })
  
  # ----------------- Visualization UI -----------------
  output$viz_panel <- renderUI({
    req(input$viz_tabs)
    
    switch(input$viz_tabs,
           "violin" = uiOutput("violin_ui"),
           "dotplot" = uiOutput("dotplot_ui"),
           "volcano" = uiOutput("volcano_ui"),
           "pathway_volcano" = uiOutput("pathway_volcano_ui"),
           "gsea_nes" = uiOutput("gsea_nes_ui")
    )
  })
  
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
        column(6, actionButton("marker_select_btn", "Select Marker", class = "btn-success", width = "100%")),
        column(6, actionButton("sample_ident_select_btn", "Select Sample Identification", class = "btn-primary", width = "100%"))
      ))
    
  })
  
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
        column(6, actionButton("marker_select_btn", "Select Marker", class = "btn-success", width = "100%")),
        column(6, actionButton("sample_ident_select_btn", "Select Sample Identification", class = "btn-primary", width = "100%"))
      ))
  })
  
  output$volcano_ui <- renderUI({
    tagList(
      h4("DEG result:"),
      path_input_ui(
        id = "deg",
        label = "DEG File Path (.tsv)",
        placeholder = "/path/to/deg_result.tsv",
        value = rv$deg_obj_path),
      actionButton("marker_select_btn", "Select Marker", class = "btn-success", width = "100%")
    )
  })
  
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
      actionButton("pathway_select_btn", "Select Pathway", class = "btn-warning", width = "100%"))
  })
  
  output$gsea_nes_ui <- renderUI({
    tagList(
      h4("GSEA result:"),
      path_input_ui(
        id = "gsea",
        label = "GSEA File Path (.xlsx)",
        placeholder = "/path/to/gsea_file.xlsx",
        value = rv$gsea_obj_path),
      actionButton("pathway_select_btn", "Select Pathway", class = "btn-warning", width = "100%"))
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
        easyClose = TRUE,
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
    source("/projectnb/wax-es/00_shinyapp/DEG/DEG/DE_analysis.R")
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
    
    source("/projectnb/wax-es/00_shinyapp/DEG/DEG/GSEA_analysis.R")  
    
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
  
  # ----------------- VISUALIZATION FUNCTION Marker Select -----------------
  
  observeEvent(input$marker_select_btn, {
    showModal(
      modalDialog(
        title = "Marker Selection",
        size = "l",
        easyClose = TRUE,
        
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
                  ))))))),
      footer = modalButton("Close"))
    )
  })
  
  observeEvent(input$marker_text, {
    rv$visualization_marker <- input$marker_text
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
    
    rv$visualization_marker <- paste(df2$segex, collapse = "\n")
    
    updateTextAreaInput(session, "marker_text", value = rv$visualization_marker)
  })
  # pct select
  observeEvent(input$apply_top_pct, {
    req(rv$deg_obj)
    df <- rv$deg_obj
    df <- df[order(df$pvalue_1), ]
    n <- ceiling(nrow(df) * input$top_pct_pvalue / 100)
    df2 <- head(df, n)
    rv$visualization_marker <- paste(df2$segex, collapse = "\n")
    updateTextAreaInput(session, "marker_text", value = rv$visualization_marker)
  })
  # less than p value selct
  observeEvent(input$apply_pvalue, {
    req(rv$deg_obj)
    
    df <- rv$deg_obj
    
    df2 <- df[df$pvalue_1 < input$pvalue_threshold, ]
    
    rv$visualization_marker <- paste(df2$segex, collapse = "\n")
    
    updateTextAreaInput(session, "marker_text", value = rv$visualization_marker)
  })
  
  # select marker on deg file
  observeEvent(input$use_deg_marker_btn, {
    library(DT)
    req(rv$deg_obj)
    
    showModal(
      modalDialog(
        title = "Select Markers from DEG",
        size = "l",
        easyClose = FALSE,
        
        DT::DTOutput("deg_table"),
        
        footer = tagList(
          actionButton(
            "confirm_deg_marker",
            "Confirm Selection",
            class = "btn btn-primary"
          ),
          modalButton("Cancel")
        )
      )
    )
  })
  
  output$deg_table <- DT::renderDT({
    
    req(rv$deg_obj)
    
    df <- rv$deg_obj
    df_display <- df
    
    num_cols <- sapply(df, is.numeric)
    
    df_display[num_cols] <- lapply(df_display[num_cols], function(x) {
      round(x, 2)
    })
    
    # ---- Get current marker genes safely ----
    marker_genes <- character(0)
    
    if (!is.null(rv$visualization_marker) && rv$visualization_marker != "") {
      marker_genes <- unlist(strsplit(rv$visualization_marker, "\n"))
      marker_genes <- trimws(marker_genes)
    }
    
    # ---- Only keep genes that exist in DEG table ----
    valid_idx <- which(df$segex %in% marker_genes)
    
    DT::datatable(
      df_display,
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
    )
  })
  
  get_selected_genes <- function() {
    
    req(rv$deg_obj)
    
    sel_idx <- input$deg_table_rows_selected
    
    if (is.null(sel_idx) || length(sel_idx) == 0) {
      return(character(0))
    }
    
    rv$deg_obj$segex[sel_idx]
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
    
    showModal(
      modalDialog(
        title = "Marker Selection",
        size = "l",
        easyClose = TRUE,
        
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
                    ))))))),
        footer = modalButton("Close"))
    )
  })
  
  # ----------------- VISUALIZATION FUNCTION Sample Identification select -----------------
}

shinyApp(ui, server)