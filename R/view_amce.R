#' @import shiny 
#' @export view

## view(amce): R Shiny exploratory visualization tool that ships with `amce`; useed 
## to quickly perform `amce` on a pre-loaded dataset/design.
## June, 2018
## Soubhik Barari

PKG.ENV <- new.env(parent = emptyenv())

# ==================================================================================
# ==================================== USER API ==================================== 
# ==================================================================================

view <- function(object) {
    ## Top-level function for the user to call in order to view 
    ## the `amce` Shiny dashboard
    
    data("japan2014conjoint", envir = environment())
    data("immigrationconjoint", envir = environment())
    data("immigrationdesign", envir = environment())

    if (typeof(object) == "closure" & deparse(substitute(object)) == "amce") {
        
        if (!requireNamespace("shiny", quietly=TRUE)) {
            stop("Package \"shiny\" needed for `view(amce)` to work. Please install it")
        }

        if (!requireNamespace("shinyjs", quietly=TRUE)) {
            stop("Package \"shinyjs\" needed for `view(amce)` to work. Please install it")
        }
        
        if (!requireNamespace("shinyBS", quietly=TRUE)) {
            stop("Package \"shinyBS\" needed for `view(amce)` to work. Please install it")
        }
        
        if (!requireNamespace("DT", quietly=TRUE)) {
            stop("Package \"DT\" needed for `view(amce)` to work. Please install it")
        }
        
        #--------------------------------------------
        # Set DESIGN/DATA globals
        #--------------------------------------------
        
        ## Global collection of design + data objects
        ## (update below)
        PKG.ENV$DATA_OPTIONS    <- NULL
        PKG.ENV$DESIGN_OPTIONS  <- NULL
        
        ## Default design + data to load on start up
        PKG.ENV$DEFAULT_DATA_NAME   <- "immigrationconjoint"
        PKG.ENV$DEFAULT_DESIGN_NAME <- "immigrationdesign"
        
        ## Names of current design + data
        PKG.ENV$SELECTED_DATA_NAME   <- PKG.ENV$DEFAULT_DATA_NAME
        PKG.ENV$SELECTED_DESIGN_NAME <- PKG.ENV$DEFAULT_DESIGN_NAME
        
        ## Objects for current design + data
        ## (update below)
        PKG.ENV$SELECTED_DATA   <- NULL
        PKG.ENV$SELECTED_DESIGN <- NULL
        
        ## Last set of results
        PKG.ENV$LAST_RESULTS    <- NULL
        
        #--------------------------------------------
        # Set UI globals
        #--------------------------------------------
        
        PKG.ENV$TITLE_TEXT        <- "cjoint"
        PKG.ENV$WINDOW_TITLE_TEXT <- "cjoint AMCE dashboard"
        
        ## Name of UI elements shown/hidden
        ## when AMCE is run
        PKG.ENV$RESULTS_ID               <- "analysisResultsPanel"
        PKG.ENV$RESULTS_SUMMARY_ID       <- "resultsSummary"
        PKG.ENV$RESULTS_PLOT_ID          <- "resultsPlot"
        PKG.ENV$RESULTS_UNCOND_PLOT_ID   <- "resultsUncondPlot"
        PKG.ENV$RESULTS_INTERACT_PLOT_ID <- "resultsInteractPlot"
        PKG.ENV$DEFAULT_RESULTS_HEIGHT   <- 1000
        
        ## Temporary "baselines" UI elements
        ## that appear/disappear based on user
        ## selection
        PKG.ENV$RESP_BASELINE_ID   <- c()
        PKG.ENV$ATTR_BASELINE_IDS  <- c()
        
        ## Plot UI elements
        PKG.ENV$PLOT_SETTINGS_ID        <- "plotSettings"
        
        PKG.ENV$PLOT_TEXT_SIZE_DEFAULT_MAX  <- 25
        PKG.ENV$PLOT_TEXT_SIZE_DEFAULT_MIN  <- 5
        PKG.ENV$PLOT_TEXT_SIZE_DEFAULT_VAL  <- 11
        PKG.ENV$PLOT_TEXT_SIZE_SELECTED            <- PKG.ENV$PLOT_TEXT_SIZE_DEFAULT_VAL
        PKG.ENV$PLOT_UNCOND_TEXT_SIZE_SELECTED     <- PKG.ENV$PLOT_TEXT_SIZE_DEFAULT_VAL
        PKG.ENV$PLOT_INTERACT_TEXT_SIZE_SELECTED   <- PKG.ENV$PLOT_TEXT_SIZE_DEFAULT_VAL
        
        PKG.ENV$PLOT_DEFAULT_CI  <- 0.95
        PKG.ENV$PLOT_CI          <- 0.95
        PKG.ENV$PLOT_UNCOND_CI   <- 0.95
        PKG.ENV$PLOT_INTERACT_CI <- 0.95
        
        ## Number of significant digits to display
        PKG.ENV$SIG_DIGITS <- 2
        
        PKG.ENV$PLOT_FACET_NAMES         <- c()
        PKG.ENV$PLOT_UNCOND_FACET_NAMES  <- c()
        PKG.ENV$PLOT_INTERACT_FACE_NAMES <- c()
        
        #--------------------------------------------
        # Set SERVER globals
        #--------------------------------------------
        
        ## Keep track of baseline selectors that need 
        ## to be rendered
        PKG.ENV$NUM_ATTR_BASELINE_IDS <- 0
        
        ## Keep track of whether we can
        ## actually plot interaction effects
        PKG.ENV$NO_RESP_VARYING <- TRUE
        PKG.ENV$NO_FACET_NAMES  <- TRUE
        PKG.ENV$NO_INTERACTIONS <- TRUE
        
        #--------------------------------------------
        # Initialize + run app
        #--------------------------------------------
        
        ## Update available data + design options
        PKG.ENV$DATA_OPTIONS          <- list(environment()$immigrationconjoint)
        names(PKG.ENV$DATA_OPTIONS)   <- c(PKG.ENV$DEFAULT_DATA_NAME)
        
        PKG.ENV$DESIGN_OPTIONS        <- list(environment()$immigrationdesign, "uniform")
        names(PKG.ENV$DESIGN_OPTIONS) <- c(PKG.ENV$DEFAULT_DESIGN_NAME, "'uniform'")
        
        .LOAD_POSSIBLE_USER_DATA()
        .LOAD_POSSIBLE_USER_DESIGNS()
        
        ## Update selected data + design (default)
        PKG.ENV$SELECTED_DATA   <- PKG.ENV$DATA_OPTIONS[[PKG.ENV$SELECTED_DATA_NAME]]
        PKG.ENV$SELECTED_DESIGN <- PKG.ENV$DESIGN_OPTIONS[[PKG.ENV$SELECTED_DESIGN_NAME]]
        
        ## Create app components
        UI     <- .buildAmceAppUI()
        SERVER <- .buildAmceAppSERVER()
        
        ## Start up
        cat("Launching `amce` browser dashboard")
        runApp(launch.browser=TRUE, list(ui = UI, server = SERVER))
    } else {
        View(object)
    }
}

# =================================================================================
# ==================================== HELPERS ==================================== 
# =================================================================================


.LOAD_POSSIBLE_USER_DATA <- function() {
    ## Scan workspace and load in possible user-generated
    ## results dataframes. Always do this after loading in 
    ## the default results options.
    
    ### get names of applicable data frames
    data_frame_names <- names(which(sapply(.GlobalEnv, is.data.frame)))
    
    data_frame_names <- data_frame_names[!(data_frame_names %in% names(PKG.ENV$DATA_OPTIONS))]
    data_frame_names <- data_frame_names[!(data_frame_names %in% names(PKG.ENV$DESIGN_OPTIONS))]
    data_frame_names <- data_frame_names[(data_frame_names != "SELECTED_DATA") & (data_frame_names != "SELECTED_DESIGN")]
    
    ### get matching data frames themselves
    data_frames <- lapply(data_frame_names, get)
    names(data_frames) <- data_frame_names
    
    ### add to globals
    PKG.ENV$DATA_OPTIONS <- c(PKG.ENV$DATA_OPTIONS, data_frames)
    
}

.LOAD_POSSIBLE_USER_DESIGNS <- function() {
    ## Scan workspace and load in possible user-generated 
    ## conjoint designs. Always do this after loading in 
    ## the default design options.
    
    ### get names of applicable designs
    design_names <- names(which(sapply(.GlobalEnv, function(x){ typeof(x) == "conjointDesign" })))
    design_names <- design_names[!(design_names %in% names(PKG.ENV$DESIGN_OPTIONS))]
    
    ### get matching designs themselves
    designs <- lapply(design_names, get)
    names(designs) <- design_names
    
    ### add to globals
    PKG.ENV$DESIGN_OPTIONS <- c(PKG.ENV$DESIGN_OPTIONS, design_names)
    
}

# =====================================================================================
# ==================================== UI BACK-END ==================================== 
# =====================================================================================

.buildAmceAppUI <- function() {
    # Function that builds the AMCE dashboard RShiny app's UI object
    # in a componentwise manner.
    
#######################################
    analysisInputPanel <- sidebarPanel(
#######################################
        # ------------------------------------------------------       
        h4("INPUT"),
        # ------------------------------------------------------       
        selectInput("designName", HTML("Choose a conjoint design (<code>design</code> object ):"), 
                    choices = c()),
        selectInput("dataName", HTML("Choose a results dataset (<code>data.frame</code>):"),
                    choices = c()),
        # ------------------------------------------------------       
        h4("RESPONSES"),
        # ------------------------------------------------------       
        selectInput("respondent.id", "Respondent ID:",
                    choices = c()),
        strong("Responses to exclude:"),
        verbatimTextOutput('selectedTableRows'),
        # ------------------------------------------------------       
        h4("ESTIMATION"),
        # ------------------------------------------------------       
        selectInput("outcome", "Outcome variable:",
                    choices=c()),
        selectInput("attributes", "Attribute variables:",
                    multiple=TRUE,
                    selectize=TRUE,
                    choices=c()),
        div(id="baselines"),
        selectInput("interactions", "Interaction terms:",
                    multiple=TRUE,
                    selectize=TRUE,
                    choices=c()),
        radioButtons("cluster", "Cluster:",
                     choices = c("Yes", "No")),
        radioButtons("na.ignore", "Ignore missing data:",
                     choices = c("Yes", "No")),
        selectInput("weights", "Weights column:",
                    choices = c()),
        selectInput("respondent.varying", "Respondent varying:",
                    choices = c()),
        
        div(id="baselinesResp"),
        # ------------------------------------------------------       
        actionButton("run", "Estimate AMCE", class="btn-primary")
        # ------------------------------------------------------       
    )

#######################################
    analysisTablePanel <- wellPanel(
#######################################
        height="500px",
        # ------------------------------------------------------
        h4("DATA"),
        # ------------------------------------------------------
        p("(Selected rows will be excluded from analysis)", style="font-size:75%"),
        div(DT::dataTableOutput('selectedDataTable'), style = "font-size: 75%;"),
        # actionButton('selectAllRows',label="Select all entries", class="btn-warning btn-sm"),
        actionButton('resetSelectedRows',label="Deselect all entries", class="btn-warning btn-sm")
    )
    
#######################################
    analysisResultsPanel <- wellPanel(
#######################################
        # ------------------------------------------------------  
        h4("RESULTS"),
        # ------------------------------------------------------
        tabsetPanel(
            tabPanel(
                # ------------------------------------------------------
                h5("Summary"),
                # ------------------------------------------------------
                br(),
                div(id=PKG.ENV$RESULTS_SUMMARY_ID,
                    verbatimTextOutput('results_summary'),
                    tags$style(type='text/css', '#results_summary {max-height: 1000px; overflow-y:scroll; }'),
                    sliderInput("sigdigits",
                                HTML("<small>Significant digits:</small>"),
                                min=1,
                                max=4,
                                value=2),
                    downloadButton("downloadSummary", "Download report")
                )
            ),
            tabPanel(
                h5("Plot"),
                br(),
                div(id=PKG.ENV$RESULTS_PLOT_ID, 
                    plotOutput("results_plot", height="800"),
                    br(),
                    selectInput("plotFacetNames",
                                HTML("<small>Facet names (optional):</small>"),
                                multiple=TRUE,
                                selectize=TRUE,
                                choices=c()),
                    textInput("plotFontFamily",
                              HTML("<small>Font family (specify if non-latin characters do not render, else leave blank):</small>")),
                    sliderInput("plotTextSize",
                                HTML("<small>Text size:</small>"),
                                min=PKG.ENV$PLOT_TEXT_SIZE_DEFAULT_MIN,
                                max=PKG.ENV$PLOT_TEXT_SIZE_DEFAULT_MAX,
                                value=PKG.ENV$PLOT_TEXT_SIZE_DEFAULT_VAL),
                    sliderInput("plotCI",
                                HTML("<small>Confidence intervals:</small>"),
                                min=0.0,
                                max=1.0,
                                value=PKG.ENV$PLOT_DEFAULT_CI,
                                step=0.01),
                    downloadButton("downloadPlot", "Download PDF")
                )
            ),
            tabPanel(
                h5("Plot (unconditional)"),
                br(),
                div(id=PKG.ENV$RESULTS_UNCOND_PLOT_ID, 
                    plotOutput("results_uncond_plot", height="800"),
                    shinyjs::hidden(
                        div(id="errorUncondMessage",
                            HTML("<br><p><i>Could not render plot. Please specify different settings below.</i></p><br>")
                        )
                    ),
                    br(),
                    textInput("plotUncondFontFamily",
                              HTML("<small>Font family (specify if non-latin characters do not render, else leave blank):</small>")),
                    # selectInput("plotUncondFacetNames",
                    #             HTML("<small>Facet names (optional):</small>"),
                    #             multiple=TRUE,
                    #             selectize=TRUE,
                    #             choices=c()),
                    sliderInput("plotUncondTextSize",
                                HTML("<small>Text size:</small>"),
                                min=PKG.ENV$PLOT_TEXT_SIZE_DEFAULT_MIN,
                                max=PKG.ENV$PLOT_TEXT_SIZE_DEFAULT_MAX,
                                value=PKG.ENV$PLOT_TEXT_SIZE_DEFAULT_VAL),
                    sliderInput("plotUncondCI",
                                HTML("<small>Confidence intervals:</small>"),
                                min=0.0,
                                max=1.0,
                                value=PKG.ENV$PLOT_DEFAULT_CI,
                                step=0.01),
                    downloadButton("downloadPlotUncond", "Download PDF")
                )
            ),
            tabPanel(
                h5("Plot (interaction)"),
                br(),
                div(id=PKG.ENV$RESULTS_INTERACT_PLOT_ID, 
                    shinyjs::hidden(
                        div(id="noInteractMessage",
                            HTML("<p><i>No interactions specified in model <i>or</i> no <code>facet.names</code>
                                 provided. Please re-run model with <code>respondent.varying</code>, and/or
                                 interaction terms specified, and (optionally) provide <code>facet.names</code> 
                                 below.</i></p><br>")
                            )
                    ),
                    plotOutput("results_interact_plot", height="800"),
                    shinyjs::hidden(
                        div(id="errorInteractMessage",
                            HTML("<br><p><i>Could not render plot. Please specify different settings below.</i></p><br>")
                        )
                    ),
                    br(),
                    selectInput("plotInteractFacetNames",
                                multiple=TRUE,
                                selectize=TRUE,
                                HTML("<small>Facet names (optional):</small>"),
                                choices=c()),
                    textInput("plotInteractFontFamily",
                              HTML("<small>Font family (specify if non-latin characters do not render, else leave blank):</small>")),
                    sliderInput("plotInteractTextSize",
                                HTML("<small>Text size:</small>"),
                                min=PKG.ENV$PLOT_TEXT_SIZE_DEFAULT_MIN,
                                max=PKG.ENV$PLOT_TEXT_SIZE_DEFAULT_MAX,
                                value=PKG.ENV$PLOT_TEXT_SIZE_DEFAULT_VAL),
                    sliderInput("plotInteractCI",
                                HTML("<small>Confidence intervals:</small>"),
                                min=0.0,
                                max=0.99,
                                value=PKG.ENV$PLOT_DEFAULT_CI,
                                step=0.01),
                    downloadButton("downloadPlotInteract", "Download PDF")
                            )
                    )
            )
        )
#######################################
    analysisTabPanel <- tagList("",
#######################################
                                analysisInputPanel,
                                mainPanel(
                                    analysisTablePanel,
                                    shinyjs::hidden(
                                        div(id=PKG.ENV$RESULTS_ID,
                                            analysisResultsPanel
                                        )
                                    )
                                )
    )
    
#######################################
    UI <- shinyUI(fluidPage(shinyjs::useShinyjs(),
#######################################
                            shinyjs::extendShinyjs(text = "shinyjs.closeWindow = function() { window.close(); }", functions = c("closeWindow")),
                            br(),
                            ## SITE TITLE
                            titlePanel(HTML("<code>cjoint</code> AMCE dashboard")),
                            ## SITE SUBTITLE
                            p(("For any loaded conjoint design and a corresponding dataset of conjoint experiment results, 
                               quickly estimate and analyze the AMCE (Average Marginal Component Effects)."),
                              ## EXIT BUTTON
                              actionButton("close", "Close app", icon = icon("times"), style='padding:4px; font-size:80%; position:absolute;right:3em; top:3em;')
                            ),
                            br(),
                            ## SITE CONTENT
                            analysisTabPanel, 
                            ## WINDOW TITLE
                            title=PKG.ENV$WINDOW_TITLE_TEXT))
    
    return(UI)
    
}


# =========================================================================================
# ==================================== SERVER BACK-END ==================================== 
# =========================================================================================

.buildAmceAppSERVER <- function() {
    # Function that builds the AMCE dashboard RShiny app's SERVER function
    # in a componentwise manner. 
    SERVER <- function(input, output, clientData, session) {
        
        observe({if (input$close > 0) {
            shinyjs::js$closeWindow()
            stopApp()
        }})
        
        .amceApp.initSession(input, output, clientData, session)
    
        .amceApp.waitForUIUpdates(input, output, clientData, session)
        .amceApp.waitForRunAction(input, output, clientData, session)
        .amceApp.waitForResultsSettingsUpdates(input, output, clientData, session)
        .amceApp.waitForPlotSettingsUpdates(input, output, clientData, session)
        .amceApp.waitForPlotUncondSettingsUpdates(input, output, clientData, session)
        .amceApp.waitForPlotInteractSettingsUpdates(input, output, clientData, session)
    }
    return(SERVER)
}

#--------------------------------------------
# INIT fxn
#--------------------------------------------

.amceApp.initSession <- function(input, output, clientData, session){
    ######################################################
    ## Initialize the current server session            ##
    ## (first thing to call upon running the Shiny app) ##
    ######################################################
    withProgress(message="Initializing dashboard session", value=0, {
        
        ### Hide elements that must conditionally appear
        # shinyjs::hide(RESULTS_SUMMARY_ID)
        # shinyjs::hide(RESULTS_PLOT_ID)
        # shinyjs::hide(PLOT_SETTINGS_ID)
        shinyjs::hide(PKG.ENV$RESULTS_ID)
        shinyjs::hide("baselines")
        shinyjs::hide("baselinesResp")
        
        # ### Initialize all UI elements with default values 
        # .amceApp.updateUI() 
        
    })
}


#--------------------------------------------
# ACTION fxns
#--------------------------------------------


.amceApp.updateUI <- function(session) {
    ## Update user interface elements based on user selected data/design
    ## which is stored in global variables
    updateSelectInput(session, 
                      "designName",
                      choices=names(PKG.ENV$DESIGN_OPTIONS),
                      selected=PKG.ENV$SELECTED_DESIGN_NAME)
    
    updateSelectInput(session, 
                      "dataName",
                      choices=names(PKG.ENV$DATA_OPTIONS),
                      selected=PKG.ENV$SELECTED_DATA_NAME)
    
    updateSelectInput(session,
                      "outcome",
                      choices=names(PKG.ENV$SELECTED_DATA),
                      selected=NULL)
    
    updateSelectInput(session,
                      "attributes",
                      choices=names(PKG.ENV$SELECTED_DATA),
                      selected=NULL)
    
    updateSelectInput(session,
                      "interactions",
                      choices=NULL,
                      selected=NULL)
    
    updateSelectInput(session, 
                      "respondent.id",
                      choices=names(PKG.ENV$SELECTED_DATA),
                      selected=names(PKG.ENV$SELECTED_DATA)[1])
    
    updateSelectInput(session, 
                      "respondent.varying",
                      choices=c("None", names(PKG.ENV$SELECTED_DATA)),
                      selected="None")
    
    updateSelectInput(session, 
                      "weights",
                      choices=c("None", names(PKG.ENV$SELECTED_DATA)),
                      selected="None") 
}

.amceApp.getInteractionTerms <- function(varNames) {
    ## From a list of variables, e.g. [a, b], output
    ## all possible interaction terms, e.g.
    ## [a:b].
    cartesProd <- expand.grid(varNames, varNames)
    interTerms <- sapply(1:nrow(cartesProd), 
                         function(s) { 
                             paste(cartesProd[s, "Var1"], cartesProd[s, "Var2"], sep=":") 
                         })
    return(interTerms)
}

.amceApp.getResults <- function(input, excl) {
    ########################################
    ## Perform an AMCE run given UI input ##
    ########################################
    encFxn <- function(x){return(paste0("`",x,"`"))}
    
    ### Create formula
    interactions <- sapply(input$interactions, function(x){paste(sapply(strsplit(x, ":"), encFxn), collapse=":")}, USE.NAMES=FALSE)
    all_terms    <- c(sapply(c(input$attributes), encFxn, USE.NAMES=FALSE), interactions)
    fmula        <- as.formula(paste(encFxn(input$outcome), paste(all_terms, collapse=" + "), sep=" ~ "))
    
    ### Subset on selected rows
    if (!is.null(input$selectedDataTable_rows_selected)) {
        rowSubset <- sapply(1:nrow(PKG.ENV$SELECTED_DATA), function(x){ return(ifelse(x %in% input$selectedDataTable_rows_selected, FALSE, TRUE)) })
    } else {
        rowSubset <- NULL
    }
    
    if (input$respondent.id == "None"){
        respondentID <-  NULL
    } else {
        respondentID <-  input$respondent.id
    }
    
    if (input$respondent.varying == "None"){
        respondentVarying <-  NULL
        PKG.ENV$NO_RESP_VARYING <- TRUE
    } else {
        respondentVarying <- input$respondent.varying
        PKG.ENV$NO_RESP_VARYING <- FALSE
    }
    
    if (length(input$interactions) == 0) {
        PKG.ENV$NO_INTERACTIONS <- TRUE
    } else {
        PKG.ENV$NO_INTERACTIONS <- FALSE
    }
    
    ### Collect user inputs for attr baselines if any 
    baselines <- c()
    for (name in colnames(PKG.ENV$SELECTED_DATA)){
        baseline <- input[[paste0("btn", name)]]
        if (!(is.null(baseline))){
            baselines[[name]] <- baseline
        }
    }
    ### Collect user inputs for resp. varying baselines if any 
    for (name in colnames(PKG.ENV$SELECTED_DATA)){
        baseline <- input[[paste0("btnR", name)]]
        baseline <- input[[paste0("btn", name)]]
        if (!(is.null(baseline))){
            baselines[[name]] <- baseline
        }
    }

    ### Run AMCE fxn with all other user inputs as options
    results <- amce(fmula,
                    data=PKG.ENV$SELECTED_DATA,
                    subset=rowSubset,
                    cluster=ifelse(input$cluster == "Yes", TRUE, FALSE),
                    respondent.id=respondentID,
                    respondent.varying=respondentVarying, 
                    na.ignore=ifelse(input$na.ignore == "Yes", TRUE, FALSE),
                    baselines=baselines,
                    design=PKG.ENV$SELECTED_DESIGN)
    return(results)
}


.amceApp.plotResults <- function(results, ylim, textsize, ci, facetnames, font.family=NULL){
    #####################################
    ## Plot the results of an AMCE run ##
    #####################################
    graphics::plot(results, 
         height=PKG.ENV$DEFAULT_RESULTS_HEIGHT,
         ylim=ylim,
         ci=ci,
         # breaks=c(-.2, 0, .2),
         # labels=c("-.2","0",".2"), 
         text.size=textsize,
         plot.theme=theme(text=element_text(family=font.family, size=textsize)),
         plot.display="all",
         facet.names=facetnames)
}


.amceApp.plotUncondResults <- function(results, ylim, textsize, ci, facetnames, font.family=NULL){
    ############################################
    ## Plot the UNCOND results of an AMCE run ##
    ############################################
    graphics::plot(results, 
         height=PKG.ENV$DEFAULT_RESULTS_HEIGHT,
         ylim=ylim,
         ci=ci,
         # breaks=c(-.2, 0, .2),
         # labels=c("-.2","0",".2"), 
         text.size=textsize,
         plot.theme=theme(text=element_text(family=font.family, size=textsize)),
         plot.display="unconditional",
         facet.names=facetnames)
}


.amceApp.plotInteractResults <- function(results, ylim, textsize, ci, facetnames, font.family=NULL){
    
    if (PKG.ENV$NO_RESP_VARYING & PKG.ENV$NO_FACET_NAMES) {
        shinyjs::show("noInteractMessage")
    } else {
        shinyjs::hide("noInteractMessage")
    }
    
    ##############################################
    ## Plot the INTERACT results of an AMCE run ##
    ##############################################
    graphics::plot(results, 
         height=PKG.ENV$DEFAULT_RESULTS_HEIGHT,
         ylim=ylim,
         ci=ci,
         # breaks=c(-.2, 0, .2),
         # labels=c("-.2","0",".2"), 
         text.size=textsize,
         plot.display="interaction",
         plot.theme=theme(text=element_text(family=font.family, size=textsize)),
         facet.names=facetnames)
}

#--------------------------------------------
# LISTENER fxns
#--------------------------------------------

.amceApp.waitForUIUpdates <- function(input, output, clientData, session) {
    ###############################################
    ## Update the dashboard state based on input ## 
    ## selected by the user                      ##
    ###############################################
    
    observeEvent(input$designName, {
        ## (Design is changed)
        if(input$designName != PKG.ENV$SELECTED_DESIGN_NAME) {
            PKG.ENV$SELECTED_DESIGN_NAME <- input$designName
            PKG.ENV$SELECTED_DESIGN <- PKG.ENV$DESIGN_OPTIONS[[PKG.ENV$SELECTED_DESIGN_NAME]]
        }
    })
    
    observeEvent(input$dataName, {
        ## (Dataset is changed)
        
        if(input$dataName != PKG.ENV$SELECTED_DATA_NAME) {
            
            ### Remove baseline selectors for existing dataset
            for(i in 1:length(colnames(PKG.ENV$SELECTED_DATA))) {
                attrBaselineID <- gsub(" ", "_", colnames(PKG.ENV$SELECTED_DATA)[i])
                removeUI(selector=paste0("div#", attrBaselineID), immediate=TRUE)
            }
            
            ### Extract chosen input or choose default input
            if(input$designName == "" | input$dataName == "") {
                PKG.ENV$SELECTED_DESIGN_NAME <- PKG.ENV$DEFAULT_DESIGN_NAME
                PKG.ENV$SELECTED_DATA_NAME   <- PKG.ENV$DEFAULT_DATA_NAME
            } else {
                PKG.ENV$SELECTED_DESIGN_NAME <- input$designName
                PKG.ENV$SELECTED_DATA_NAME   <- input$dataName
            }
            
            PKG.ENV$SELECTED_DESIGN <-  PKG.ENV$DESIGN_OPTIONS[[PKG.ENV$SELECTED_DESIGN_NAME]]
            PKG.ENV$SELECTED_DATA   <-  PKG.ENV$DATA_OPTIONS[[PKG.ENV$SELECTED_DATA_NAME]]
            
            ### Update data/design specific elements dynamically
            withProgress(message="Updating UI", value=0, {
                
                .amceApp.updateUI(session)
                
                # updateSelectInput(session, 
                #                   "respondent.id",
                #                   choices=names(SELECTED_DATA),
                #                   selected=names(SELECTED_DATA)[1])
                # 
                # updateSelectInput(session, 
                #                   "weights",
                #                   choices=c("None", names(SELECTED_DATA)),
                #                   selected="None")
                
            })
            
            ### Update data table
            selectedRows <- reactiveValues(sel=c())
            
            output$selectedDataTable <- DT::renderDataTable(
                DT::datatable(
                    input$resetSelectedRows,
                    data=PKG.ENV$SELECTED_DATA,
                    options=list(
                        pageLength=5,
                        scrollX=TRUE,
                        scrollY=TRUE,
                        lengthChange=FALSE,
                        processing=FALSE,
                        rowCallback=DT::JS("function(r,d) {$(r).attr('min-height', '20px')}")
                    ),
                    selection=list(mode="multiple", selected=selectedRows$sel),
                    rownames = TRUE,
                    width="100%"
                ))
            
            output$selectedTableRows <- renderPrint({
                s = input$selectedDataTable_rows_selected
                if (length(s)) {
                    cat(s, sep=",\n")
                } else {
                    cat("None")
                }
            })
            
            
            
        }
    }, ignoreInit=FALSE)
    
    observeEvent(input$attributes, {
        ## (Input attributes are changed)
        
        updateSelectInput(session,
                          "interactions",
                          choices=.amceApp.getInteractionTerms(input$attributes),
                          selected=NULL)
        
        if(!is.null(input$attributes)) {
            shinyjs::show("baselines")
        }
        
        if (length(input$attributes) > PKG.ENV$NUM_ATTR_BASELINE_IDS) {
            ### ADD baseline factors if needed
            PKG.ENV$NUM_ATTR_BASELINE_IDS <- length(input$attributes)
            
            for(i in 1:length(input$attributes)) {
                attrBaselineID <- gsub(" ", "_", input$attributes[i])
                if (is.factor(PKG.ENV$SELECTED_DATA[[input$attributes[i]]]) & !(attrBaselineID %in% PKG.ENV$ATTR_BASELINE_IDS) ) {
                    PKG.ENV$ATTR_BASELINE_IDS <- unique(c(PKG.ENV$ATTR_BASELINE_IDS, attrBaselineID))
                    insertUI(selector="#baselines",
                             ui=tags$div(id=attrBaselineID)
                    )
                    
                    insertUI(selector=paste0("#",attrBaselineID),
                             ui=selectInput(paste0("btn",attrBaselineID), 
                                            HTML(paste0("<small>&emsp;&#9830; ",input$attributes[i]," baseline</small>:")), 
                                            levels(PKG.ENV$SELECTED_DATA[[input$attributes[i]]]))
                    )
                }    
            }
            
        } else {
            ### REMOVE baseline factors if needed
            PKG.ENV$NUM_ATTR_BASELINE_IDS <- length(input$attributes)
            
            for(i in 1:length(colnames(PKG.ENV$SELECTED_DATA))) {
                attrBaselineID <- gsub(" ", "_", colnames(PKG.ENV$SELECTED_DATA)[i])
                if ((attrBaselineID %in% PKG.ENV$ATTR_BASELINE_IDS) & !(colnames(PKG.ENV$SELECTED_DATA)[i] %in% input$attributes)){
                    removeUI(selector=paste0("div#", attrBaselineID), immediate=TRUE)
                }
            }
        }
        
    })
    
    observeEvent(input$respondent.varying, {
        ### (Respondent varying variable is changed)
        
        if(input$respondent.varying != "None") {
            shinyjs::show("baselinesResp")
            if (is.factor(PKG.ENV$SELECTED_DATA[[input$respondent.varying]])) {
                ### REMOVE previous baseline factor (if exists)
                removeUI(selector=paste0("div#", PKG.ENV$RESP_BASELINE_ID), immediate=TRUE)
                
                ### ADD baseline factor
                PKG.ENV$RESP_BASELINE_ID <- gsub(" ", "_", input$respondent.varying)
                PKG.ENV$RESP_BASELINE_ID <- paste0("R", PKG.ENV$RESP_BASELINE_ID)
                
                insertUI(selector="#baselinesResp",
                         ui=tags$div(id=PKG.ENV$RESP_BASELINE_ID)
                )
                
                insertUI(selector=paste0("#", PKG.ENV$RESP_BASELINE_ID),
                         ui=selectInput(paste0("btn", PKG.ENV$RESP_BASELINE_ID), 
                                        HTML(paste0("<small>&emsp;&emsp;&emsp;&emsp;&#9830; ",input$respondent.varying," baseline</small>:")), 
                                        levels(PKG.ENV$SELECTED_DATA[[input$respondent.varying]]))
                ) 
            } else {
                removeUI(selector=paste0("div#", PKG.ENV$RESP_BASELINE_ID), immediate=TRUE)
                shinyjs::hide("baselinesResp")
            }
        } else {
            shinyjs::hide("baselinesResp")
        }
        
    })
    
}

.amceApp.waitForRunAction <- function(input, output, clientData, session) {
    #############################################################
    ## Update the dashboard when the user chooses to hit "run" ##
    #############################################################
    observeEvent(input$run, {
        withProgress(message="Computing results", value=0, {
            incProgress(1/4, detail = "Running AMCE estimator")
            
            tryCatch({
                ## Hide `facet.names` selectors in case there's an
                ## existing plot already shown (unsafe to change b/t plots!)
                shinyjs::hide("plotFacetNames")
                shinyjs::hide("plotUncondFacetNames")
                shinyjs::hide("plotInteractFacetNames")
                
                PKG.ENV$LAST_RESULTS <- .amceApp.getResults(input)
                PKG.ENV$SIG_DIGITS   <- input$sigdigits
                
                incProgress(2/4, detail="Getting summary")
                output$results_summary <- renderPrint( 
                    ## Show according to significant digits specified by user
                    print(summary(PKG.ENV$LAST_RESULTS), digits=PKG.ENV$SIG_DIGITS)
                )
                
                incProgress(3/4, detail="Creating plot")
                
                ## Update possible `facet.names` the user can select for each plot
                interacted_terms <- unique(Reduce(c,sapply(input$interactions, function(x) { strsplit(x, ":") }) ))
                if (input$respondent.varying == "None") {
                    possible_facet_names <- interacted_terms
                } else {
                    possible_facet_names <- c(input$respondent.varying, interacted_terms)
                }
                
                PKG.ENV$PLOT_FONT_FAMILY          <- NULL
                PKG.ENV$PLOT_UNCOND_FONT_FAMILY   <- NULL
                PKG.ENV$PLOT_INTERACT_FONT_FAMILY <- NULL

                PKG.ENV$PLOT_FACET_NAMES <- NULL
                updateSelectInput(session, "plotFacetNames",
                                  choices=possible_facet_names,
                                  selected=NULL)
                
                PKG.ENV$PLOT_UNCOND_FACET_NAMES <- NULL
                updateSelectInput(session, "plotUncondFacetNames",
                                  choices=possible_facet_names,
                                  selected=NULL)
                
                PKG.ENV$PLOT_INTERACT_FACET_NAMES <- NULL
                updateSelectInput(session, "plotInteractFacetNames",
                                  choices=possible_facet_names,
                                  selected=NULL)
                
                ## Show `facet.names` selectors when it's safe again
                shinyjs::show("plotFacetNames")
                shinyjs::show("plotUncondFacetNames")
                shinyjs::show("plotInteractFacetNames")
                
                output$results_plot <- renderPlot({
                    withProgress(message="Rendering plot (all estimates) in UI", {
                        .amceApp.plotResults(PKG.ENV$LAST_RESULTS, 
                                    ci=PKG.ENV$PLOT_DEFAULT_CI,
                                    ylim=c(PKG.ENV$PLOT_YLIM_DEFAULT_MIN, PKG.ENV$PLOT_YLIM_DEFAULT_MAX),
                                    textsize=PKG.ENV$PLOT_TEXT_SIZE_DEFAULT_VAL,
                                    facetnames=PKG.ENV$PLOT_FACET_NAMES)       
                    })
                })
                
                output$results_uncond_plot <- renderPlot({
                    withProgress(message="Rendering plot (unconditional estimates) in UI", {
                        .amceApp.plotUncondResults(PKG.ENV$LAST_RESULTS, 
                                          ci=PKG.ENV$PLOT_DEFAULT_CI,
                                          ylim=c(PKG.ENV$PLOT_YLIM_DEFAULT_MIN, PKG.ENV$PLOT_YLIM_DEFAULT_MAX),
                                          textsize=PKG.ENV$PLOT_TEXT_SIZE_DEFAULT_VAL,
                                          facetnames=PKG.ENV$PLOT_UNCOND_FACET_NAMES)      
                    })
                })
                
                output$results_interact_plot <- renderPlot({
                    withProgress(message="Rendering plot (interaction estimates) in UI", {
                        .amceApp.plotInteractResults(PKG.ENV$LAST_RESULTS, 
                                            ci=PKG.ENV$PLOT_DEFAULT_CI,
                                            ylim=c(PKG.ENV$PLOT_YLIM_DEFAULT_MIN, PKG.ENV$PLOT_YLIM_DEFAULT_MAX),
                                            textsize=PKG.ENV$PLOT_TEXT_SIZE_DEFAULT_VAL,
                                            facetnames=PKG.ENV$PLOT_INTERACT_FACET_NAMES)      
                    })
                })
                
                
                incProgress(4/4, detail="Generating downloads")
                output$downloadSummary <- downloadHandler(
                    filename =  paste(PKG.ENV$SELECTED_DATA_NAME, "summary", Sys.Date(), ".txt", sep="_"),
                    content = function(file){
                        sink(file=file)
                        print(summary(PKG.ENV$LAST_RESULTS), digits=PKG.ENV$SIG_DIGITS)
                        sink()
                    })
                output$downloadPlot <- downloadHandler(
                    filename =  paste(PKG.ENV$SELECTED_DATA_NAME, "plot", Sys.Date(), ".pdf"),
                    content = function(file){
                        pdf(file, width=8.5, height=15, paper="special")
                        .amceApp.plotResults(PKG.ENV$LAST_RESULTS, 
                                    ylim=PKG.ENV$PLOT_YLIM_SELECTED,
                                    ci=PKG.ENV$PLOT_CI,
                                    textsize=PKG.ENV$PLOT_TEXT_SIZE_SELECTED,
                                    facetnames=PKG.ENV$PLOT_FACET_NAMES)                          
                        dev.off()
                    })
                output$downloadPlotUncond <- downloadHandler(
                    filename =  paste(PKG.ENV$SELECTED_DATA_NAME, "plotUnconditional", Sys.Date(), ".pdf", sep="_"),
                    content = function(file){
                        pdf(file, width=8.5, height=15, paper="special")
                        .amceApp.plotUncondResults(PKG.ENV$LAST_RESULTS, 
                                          ylim=PKG.ENV$PLOT_UNCOND_YLIM_SELECTED,
                                          ci=PKG.ENV$PLOT_UNCOND_CI,
                                          textsize=PKG.ENV$PLOT_UNCOND_TEXT_SIZE_SELECTED,
                                          facetnames=PKG.ENV$PLOT_UNCOND_FACET_NAMES)                          
                        dev.off()
                    })
                output$downloadPlotInteract <- downloadHandler(
                    filename =  paste(PKG.ENV$SELECTED_DATA_NAME, "plotInteraction", Sys.Date(), ".pdf", sep="_"),
                    content = function(file){
                        pdf(file, width=8.5, height=15, paper="special")
                        .amceApp.plotInteractResults(PKG.ENV$LAST_RESULTS, 
                                            ylim=PKG.ENV$PLOT_INTERACT_YLIM_SELECTED,
                                            ci=PKG.ENV$PLOT_INTERACT_CI,
                                            textsize=PKG.ENV$PLOT_INTERACT_TEXT_SIZE_SELECTED,
                                            facetnames=PKG.ENV$PLOT_INTERACT_FACET_NAMES)                          
                        dev.off()
                    })
                
                
                shinyjs::show(PKG.ENV$RESULTS_ID)
            },
            error=function(e){
                
                # shinyjs::hide(RESULTS_ID)
                
                ### Check for obvious user input errors first
                if (is.null(input$attributes)) {
                    e <- "Please enter a valid set of attributes."
                }
                
                showModal(modalDialog(e, title="Error", size="m", easyClose=TRUE))
                # print(e)
            })
        })
    })
}


.amceApp.waitForResultsSettingsUpdates <- function(input, output, clientData, session){
    #################################################################
    ## Functions to run when user adjusts results summary settings ##
    #################################################################
    observeEvent(input$sigdigits, {
        
        PKG.ENV$SIG_DIGITS <- input$sigdigits
        
        output$results_summary <- renderPrint( 
            ## Show according to significant digits specified by user
            print(summary(PKG.ENV$LAST_RESULTS), digits=PKG.ENV$SIG_DIGITS)
        )
    })
    
}

.amceApp.waitForPlotSettingsUpdates <- function(input, output, clientData, session){
    ######################################################
    ## Functions to run when user adjusts plot settings ##
    ######################################################
    observeEvent(input$plotTextSize, {
        PKG.ENV$PLOT_TEXT_SIZE_SELECTED <- input$plotTextSize
        
        output$results_plot <- renderPlot({
            withProgress(message="Rendering plot (all estimates) in UI", {
                .amceApp.plotResults(PKG.ENV$LAST_RESULTS, 
                            ci=PKG.ENV$PLOT_CI,
                            ylim=PKG.ENV$PLOT_YLIM_SELECTED,
                            textsize=PKG.ENV$PLOT_TEXT_SIZE_SELECTED,
                            facetnames=PKG.ENV$PLOT_FACET_NAMES,
                            font.family=PKG.ENV$PLOT_FONT_FAMILY)       
            })
        })
    })
    
    observeEvent(input$plotFontFamily, {
        PKG.ENV$PLOT_FONT_FAMILY <- input$plotFontFamily
        
        output$results_plot <- renderPlot({
            withProgress(message="Rendering plot (all estimates) in UI", {
                .amceApp.plotResults(PKG.ENV$LAST_RESULTS, 
                                    ci=PKG.ENV$PLOT_CI,
                                    ylim=PKG.ENV$PLOT_YLIM_SELECTED,
                                    textsize=PKG.ENV$PLOT_TEXT_SIZE_SELECTED,
                                    facetnames=PKG.ENV$PLOT_FACET_NAMES,
                                    font.family=PKG.ENV$PLOT_FONT_FAMILY)       
            })
        })
    }) 
    
    observeEvent(input$plotCI, {
        PKG.ENV$PLOT_CI <- input$plotCI
        
        output$results_plot <- renderPlot({
            withProgress(message="Rendering plot (all estimates) in UI", {
                .amceApp.plotResults(PKG.ENV$LAST_RESULTS, 
                            ci=PKG.ENV$PLOT_CI,
                            ylim=PKG.ENV$PLOT_YLIM_SELECTED,
                            textsize=PKG.ENV$PLOT_TEXT_SIZE_SELECTED,
                            facetnames=PKG.ENV$PLOT_FACET_NAMES,
                            font.family=PKG.ENV$PLOT_FONT_FAMILY)       
            })
        })
    })
    
    observeEvent(input$plotFacetNames, {
        
        PKG.ENV$PLOT_FACET_NAMES <- input$plotFacetNames
        
        output$results_plot <- renderPlot({
            withProgress(message="Rendering plot (all estimates) in UI", {
                .amceApp.plotResults(PKG.ENV$LAST_RESULTS, 
                            ci=PKG.ENV$PLOT_CI,
                            ylim=PKG.ENV$PLOT_YLIM_SELECTED,
                            textsize=PKG.ENV$PLOT_TEXT_SIZE_SELECTED,
                            facetnames=PKG.ENV$PLOT_FACET_NAMES,
                            font.family=PKG.ENV$PLOT_FONT_FAMILY)       
            })
        })  
    }, ignoreNULL=FALSE)
}


.amceApp.waitForPlotUncondSettingsUpdates <- function(input, output, clientData, session){
    #############################################################
    ## Functions to run when user adjusts UNCOND plot settings ##
    #############################################################
    observeEvent(input$plotUncondTextSize, {
        
        PKG.ENV$PLOT_UNCOND_TEXT_SIZE_SELECTED <- input$plotUncondTextSize
        
        output$results_uncond_plot <- renderPlot({
            withProgress(message="Rendering plot (unconditional estimates) in UI", {
                .amceApp.plotUncondResults(PKG.ENV$LAST_RESULTS, 
                                  ci=PKG.ENV$PLOT_UNCOND_CI,
                                  ylim=PKG.ENV$PLOT_UNCOND_YLIM_SELECTED,
                                  textsize=PKG.ENV$PLOT_UNCOND_TEXT_SIZE_SELECTED,
                                  facetnames=PKG.ENV$PLOT_UNCOND_FACET_NAMES,
                                  font.family=PKG.ENV$PLOT_UNCOND_FONT_FAMILY)       
            })
        })
    })
    
    observeEvent(input$plotUncondFontFamily, {
        PKG.ENV$PLOT_UNCOND_FONT_FAMILY <- input$plotUncondFontFamily
        
        output$results_uncond_plot <- renderPlot({
            withProgress(message="Rendering plot (all estimates) in UI", {
                .amceApp.plotUncondResults(PKG.ENV$LAST_RESULTS, 
                                     ci=PKG.ENV$PLOT_UNCOND_CI,
                                     ylim=PKG.ENV$PLOT_UNCOND_YLIM_SELECTED,
                                     textsize=PKG.ENV$PLOT_UNCOND_TEXT_SIZE_SELECTED,
                                     facetnames=PKG.ENV$PLOT_UNCOND_FACET_NAMES,
                                     font.family=PKG.ENV$PLOT_UNCOND_FONT_FAMILY)       
            })
        })        
    }) 
    
    observeEvent(input$plotUncondCI, {
        PKG.ENV$PLOT_UNCOND_CI <- input$plotUncondCI
        
        output$results_uncond_plot <- renderPlot({
            withProgress(message="Rendering plot (unconditional estimates) in UI", {
                .amceApp.plotUncondResults(PKG.ENV$LAST_RESULTS, 
                                  ci=PKG.ENV$PLOT_UNCOND_CI,
                                  ylim=PKG.ENV$PLOT_UNCOND_YLIM_SELECTED,
                                  textsize=PKG.ENV$PLOT_UNCOND_TEXT_SIZE_SELECTED,
                                  facetnames=PKG.ENV$PLOT_UNCOND_FACET_NAMES,
                                  font.family=PKG.ENV$PLOT_UNCOND_FONT_FAMILY)       
            })
        })
    })
    
    observeEvent(input$plotUncondFacetNames, {
        
        PKG.ENV$PLOT_UNCOND_FACET_NAMES <- input$plotUncondFacetNames
        
        output$results_uncond_plot <- renderPlot({
            withProgress(message="Rendering plot (unconditional estimates) in UI", {
                .amceApp.plotUncondResults(PKG.ENV$LAST_RESULTS, 
                                  ci=PKG.ENV$PLOT_UNCOND_CI,
                                  ylim=PKG.ENV$PLOT_UNCOND_YLIM_SELECTED,
                                  textsize=PKG.ENV$PLOT_UNCOND_TEXT_SIZE_SELECTED,
                                  facetnames=PKG.ENV$PLOT_UNCOND_FACET_NAMES,
                                  font.family=PKG.ENV$PLOT_UNCOND_FONT_FAMILY)       
            })
        })  
    }, ignoreNULL=FALSE)
}


.amceApp.waitForPlotInteractSettingsUpdates <- function(input, output, clientData, session){
    #############################################################
    ## Functions to run when user adjusts UNCOND plot settings ##
    #############################################################
    observeEvent(input$plotInteractTextSize, {
        
        PKG.ENV$PLOT_INTERACT_TEXT_SIZE_SELECTED <- input$plotInteractTextSize
        
        output$results_interact_plot <- renderPlot({
            withProgress(message="Rendering plot (interaction estimates) in UI", {
                .amceApp.plotInteractResults(PKG.ENV$LAST_RESULTS, 
                                             ci=PKG.ENV$PLOT_INTERACT_CI,
                                             ylim=PKG.ENV$PLOT_INTERACT_YLIM_SELECTED,
                                             textsize=PKG.ENV$PLOT_INTERACT_TEXT_SIZE_SELECTED,
                                             facetnames=PKG.ENV$PLOT_INTERACT_FACET_NAMES,
                                             font.family=PKG.ENV$PLOT_INTERACT_FONT_FAMILY)       
            })
        })
        
    })
    
    
    observeEvent(input$plotInteractFontFamily, {
        PKG.ENV$PLOT_INTERACT_FONT_FAMILY <- input$plotInteractFontFamily
        
        output$results_interact_plot <- renderPlot({
            withProgress(message="Rendering plot (all estimates) in UI", {
                .amceApp.plotInteractResults(PKG.ENV$LAST_RESULTS, 
                                             ci=PKG.ENV$PLOT_INTERACT_CI,
                                             ylim=PKG.ENV$PLOT_INTERACT_YLIM_SELECTED,
                                             textsize=PKG.ENV$PLOT_INTERACT_TEXT_SIZE_SELECTED,
                                             facetnames=PKG.ENV$PLOT_INTERACT_FACET_NAMES,
                                             font.family=PKG.ENV$PLOT_INTERACT_FONT_FAMILY)       
            })
        })        
    }) 
    
    
    observeEvent(input$plotInteractCI, {
        
        PLOT_INTERACT_CI <- input$plotInteractCI
        
        output$results_interact_plot <- renderPlot({
            withProgress(message="Rendering plot (interaction estimates) in UI", {
                .amceApp.plotInteractResults(PKG.ENV$LAST_RESULTS, 
                                             ci=PKG.ENV$PLOT_INTERACT_CI,
                                             ylim=PKG.ENV$PLOT_INTERACT_YLIM_SELECTED,
                                             textsize=PKG.ENV$PLOT_INTERACT_TEXT_SIZE_SELECTED,
                                             facetnames=PKG.ENV$PLOT_INTERACT_FACET_NAMES,
                                             font.family=PKG.ENV$PLOT_INTERACT_FONT_FAMILY)       
            })
        })
    })
    
    observeEvent(input$plotInteractFacetNames, {
        
        PKG.ENV$PLOT_INTERACT_FACET_NAMES <- input$plotInteractFacetNames
        
        if (length(PKG.ENV$PLOT_INTERACT_FACET_NAMES) == 0) {
            PKG.ENV$NO_FACET_NAMES <- TRUE
        } else {
            PKG.ENV$NO_FACET_NAMES <- FALSE
        }
        
        output$results_interact_plot <- renderPlot({
            withProgress(message="Rendering plot (interaction estimates) in UI", {
                .amceApp.plotInteractResults(PKG.ENV$LAST_RESULTS, 
                                             ci=PKG.ENV$PLOT_INTERACT_CI,
                                             ylim=PKG.ENV$PLOT_INTERACT_YLIM_SELECTED,
                                             textsize=PKG.ENV$PLOT_INTERACT_TEXT_SIZE_SELECTED,
                                             facetnames=PKG.ENV$PLOT_INTERACT_FACET_NAMES,
                                             font.family=PKG.ENV$PLOT_INTERACT_FONT_FAMILY)       
            })
        })  
    }, ignoreNULL=FALSE)
    
}

