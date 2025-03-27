library(shiny)
library(bs4Dash)
library(metafor)
library(esc)
library(DT)
library(zip)

### Custom Function: Compute standard error from p-value
se.from.p <- function(effect.size, p, N, effect.size.type) {
  df <- N - 2
  t_val <- qt(1 - p/2, df = df)
  se <- effect.size / t_val
  return(list(EffectSize = effect.size, StandardError = se))
}

### Helper: Convert a value to numeric if possible
maybe_numeric <- function(x) {
  num <- suppressWarnings(as.numeric(x))
  if (is.na(num) && !is.na(x) && x != "") return(x) else return(num)
}

### Helper: Compute effect size and variance from raw data with extra debugging
computeES <- function(conv, params) {
  message("Inside computeES for conversion: ", conv)
  message("Parameters: ", paste(names(params), params, collapse = ", "))
  res <- switch(conv,
                "Mean & SE" = {
                  out <- esc_mean_se(grp1m = as.numeric(params$grp1m),
                                     grp1se = as.numeric(params$grp1se),
                                     grp1n = as.numeric(params$grp1n),
                                     grp2m = as.numeric(params$grp2m),
                                     grp2se = as.numeric(params$grp2se),
                                     grp2n = as.numeric(params$grp2n),
                                     es.type = params$esType)
                  message("esc_mean_se raw output:")
                  message(paste(capture.output(str(out)), collapse = "\n"))
                  list(yi = as.numeric(out$es), vi = (as.numeric(out$se))^2)
                },
                "Reg_B" = {
                  out <- esc_B(b = as.numeric(params$b),
                               sdy = as.numeric(params$sdy),
                               grp1n = as.numeric(params$grp1n),
                               grp2n = as.numeric(params$grp2n),
                               es.type = params$esType)
                  message("esc_B raw output:")
                  message(paste(capture.output(str(out)), collapse = "\n"))
                  list(yi = as.numeric(out$es), vi = (as.numeric(out$se))^2)
                },
                "Reg_beta" = {
                  out <- esc_beta(beta = as.numeric(params$beta),
                                  sdy = as.numeric(params$sdy),
                                  grp1n = as.numeric(params$grp1n),
                                  grp2n = as.numeric(params$grp2n),
                                  es.type = params$esType)
                  message("esc_beta raw output:")
                  message(paste(capture.output(str(out)), collapse = "\n"))
                  list(yi = as.numeric(out$es), vi = (as.numeric(out$se))^2)
                },
                "PointBiserial" = {
                  out <- esc_rpb(r = as.numeric(params$rpb),
                                 grp1n = as.numeric(params$grp1n),
                                 grp2n = as.numeric(params$grp2n),
                                 es.type = params$esType)
                  message("esc_rpb raw output:")
                  message(paste(capture.output(str(out)), collapse = "\n"))
                  list(yi = as.numeric(out$es), vi = (as.numeric(out$se))^2)
                },
                "OneWayANOVA" = {
                  out <- esc_f(f = as.numeric(params$f),
                               grp1n = as.numeric(params$grp1n),
                               grp2n = as.numeric(params$grp2n),
                               es.type = params$esType)
                  message("esc_f raw output:")
                  message(paste(capture.output(str(out)), collapse = "\n"))
                  list(yi = as.numeric(out$es), vi = (as.numeric(out$se))^2)
                },
                "TwoSample_t" = {
                  out <- esc_t(t = as.numeric(params$t),
                               grp1n = as.numeric(params$grp1n),
                               grp2n = as.numeric(params$grp2n),
                               es.type = params$esType)
                  message("esc_t raw output:")
                  message(paste(capture.output(str(out)), collapse = "\n"))
                  list(yi = as.numeric(out$es), vi = (as.numeric(out$se))^2)
                },
                "p_to_SE" = {
                  out <- se.from.p(effect.size = as.numeric(params$effect_size),
                                   p = as.numeric(params$p),
                                   N = as.numeric(params$N),
                                   effect.size.type = params$effectSizeType)
                  message("se.from.p raw output:")
                  message(paste(capture.output(str(out)), collapse = "\n"))
                  list(yi = as.numeric(out$EffectSize), vi = (as.numeric(out$StandardError))^2)
                },
                "ChiSquared" = {
                  out <- esc_chisq(chisq = as.numeric(params$chisq),
                                   totaln = as.numeric(params$totaln),
                                   es.type = params$esType)
                  message("esc_chisq raw output:")
                  message(paste(capture.output(str(out)), collapse = "\n"))
                  list(yi = as.numeric(out$es), vi = (as.numeric(out$se))^2)
                },
                { stop("Unsupported conversion type") }
  )
  if(length(res$yi) == 0) res$yi <- NA
  if(length(res$vi) == 0) res$vi <- NA
  message("Returning: yi = ", res$yi, "; vi = ", res$vi)
  return(res)
}

### Sample CSV Generation Helpers
generateSampleCSV <- function(conv) {
  switch(conv,
         "Mean & SE" = data.frame(
           Study = "Study1",
           ConversionType = "Mean & SE",
           grp1m = 8.5,
           grp1se = 1.5,
           grp1n = 50,
           grp2m = 11,
           grp2se = 1.8,
           grp2n = 60,
           b = NA,
           beta = NA,
           sdy = NA,
           rpb = NA,
           f = NA,
           t = NA,
           effect_size = NA,
           p = NA,
           N = NA,
           effectSizeType = NA,
           chisq = NA,
           totaln = NA,
           esType = "d",
           stringsAsFactors = FALSE
         ),
         "Reg_B" = data.frame(
           Study = "Study2",
           ConversionType = "Reg_B",
           grp1m = NA,
           grp1se = NA,
           grp1n = 100,
           grp2m = NA,
           grp2se = NA,
           grp2n = 150,
           b = 3.3,
           beta = NA,
           sdy = 5,
           rpb = NA,
           f = NA,
           t = NA,
           effect_size = NA,
           p = NA,
           N = NA,
           effectSizeType = NA,
           chisq = NA,
           totaln = NA,
           esType = "d",
           stringsAsFactors = FALSE
         ),
         "Reg_beta" = data.frame(
           Study = "Study3",
           ConversionType = "Reg_beta",
           grp1m = NA,
           grp1se = NA,
           grp1n = 100,
           grp2m = NA,
           grp2se = NA,
           grp2n = 150,
           b = NA,
           beta = 0.32,
           sdy = 5,
           rpb = NA,
           f = NA,
           t = NA,
           effect_size = NA,
           p = NA,
           N = NA,
           effectSizeType = NA,
           chisq = NA,
           totaln = NA,
           esType = "d",
           stringsAsFactors = FALSE
         ),
         "PointBiserial" = data.frame(
           Study = "Study4",
           ConversionType = "PointBiserial",
           grp1m = NA,
           grp1se = NA,
           grp1n = 99,
           grp2m = NA,
           grp2se = NA,
           grp2n = 120,
           b = NA,
           beta = NA,
           sdy = NA,
           rpb = 0.25,
           f = NA,
           t = NA,
           effect_size = NA,
           p = NA,
           N = NA,
           effectSizeType = NA,
           chisq = NA,
           totaln = NA,
           esType = "d",
           stringsAsFactors = FALSE
         ),
         "OneWayANOVA" = data.frame(
           Study = "Study5",
           ConversionType = "OneWayANOVA",
           grp1m = NA,
           grp1se = NA,
           grp1n = 519,
           grp2m = NA,
           grp2se = NA,
           grp2n = 528,
           b = NA,
           beta = NA,
           sdy = NA,
           rpb = NA,
           f = 5.04,
           t = NA,
           effect_size = NA,
           p = NA,
           N = NA,
           effectSizeType = NA,
           chisq = NA,
           totaln = NA,
           esType = "g",
           stringsAsFactors = FALSE
         ),
         "TwoSample_t" = data.frame(
           Study = "Study6",
           ConversionType = "TwoSample_t",
           grp1m = NA,
           grp1se = NA,
           grp1n = 100,
           grp2m = NA,
           grp2se = NA,
           grp2n = 150,
           b = NA,
           beta = NA,
           sdy = NA,
           rpb = NA,
           f = NA,
           t = 3.3,
           effect_size = NA,
           p = NA,
           N = NA,
           effectSizeType = NA,
           chisq = NA,
           totaln = NA,
           esType = "d",
           stringsAsFactors = FALSE
         ),
         "p_to_SE" = data.frame(
           Study = "Study7",
           ConversionType = "p_to_SE",
           grp1m = NA,
           grp1se = NA,
           grp1n = NA,
           grp2m = NA,
           grp2se = NA,
           grp2n = NA,
           b = NA,
           beta = NA,
           sdy = NA,
           rpb = NA,
           f = NA,
           t = NA,
           effect_size = 0.71,
           p = 0.013,
           N = 71,
           effectSizeType = "difference",
           chisq = NA,
           totaln = NA,
           esType = NA,
           stringsAsFactors = FALSE
         ),
         "ChiSquared" = data.frame(
           Study = "Study8",
           ConversionType = "ChiSquared",
           grp1m = NA,
           grp1se = NA,
           grp1n = NA,
           grp2m = NA,
           grp2se = NA,
           grp2n = NA,
           b = NA,
           beta = NA,
           sdy = NA,
           rpb = NA,
           f = NA,
           t = NA,
           effect_size = NA,
           p = NA,
           N = NA,
           effectSizeType = NA,
           chisq = 7.9,
           totaln = 100,
           esType = "cox.or",
           stringsAsFactors = FALSE
         )
  )
}

allConvTypes <- c("Mean & SE", "Reg_B", "Reg_beta", "PointBiserial", 
                  "OneWayANOVA", "TwoSample_t", "p_to_SE", "ChiSquared")

### Reactive storage for combined study data ###
metaData <- reactiveVal(data.frame(
  Study = character(),
  ConversionType = character(),
  yi = numeric(),
  vi = numeric(),
  stringsAsFactors = FALSE
))

### UI ###
ui <- bs4DashPage(
  title = "Meta-Analysis Calculator",
  header = bs4DashNavbar(title = "Meta-Analysis Calculator"),
  sidebar = bs4DashSidebar(
    sidebarMenu(
      menuItem("Upload Data", tabName = "upload", icon = icon("upload")),
      menuItem("Meta-Analysis", tabName = "meta", icon = icon("chart-line")),
      menuItem("Forest Plot", tabName = "plot", icon = icon("image")),
      menuItem("Funnel Plot", tabName = "funnel", icon = icon("chart-area")),
      menuItem("Baujat Plot", tabName = "baujat", icon = icon("chart-area")),
      menuItem("Influence Diagnostics", tabName = "influence", icon = icon("chart-area")),
      menuItem("Cumulative Analysis", tabName = "cumulative", icon = icon("chart-bar")),
      menuItem("Moderator Analysis", tabName = "moderator", icon = icon("sliders-h")),
      menuItem("Radial Plot", tabName = "radial", icon = icon("chart-pie")),
      menuItem("Leave-One-Out", tabName = "leaveone", icon = icon("chart-line")),
      menuItem("Cook's Distance", tabName = "cooks", icon = icon("exclamation-triangle")),
      menuItem("Data Summary", tabName = "datasummary", icon = icon("table")),
      menuItem("Additional Results", tabName = "additional", icon = icon("info-circle")),
      menuItem("Download Samples", tabName = "downloads", icon = icon("download")),
      menuItem("Instructions", tabName = "instructions", icon = icon("info-circle"))
    )
  ),
  body = bs4DashBody(
    tabItems(
      # Upload Data Tab
      tabItem(tabName = "upload",
              fluidRow(
                bs4Card(
                  title = "Upload CSV Data",
                  width = 12,
                  status = "primary",
                  solidHeader = TRUE,
                  fileInput("datafile", "Choose CSV File", accept = ".csv"),
                  actionButton("runMeta", "Run Meta-Analysis", icon = icon("play"))
                )
              ),
              fluidRow(
                bs4Card(
                  title = "Uploaded Data Preview",
                  width = 12,
                  status = "info",
                  solidHeader = TRUE,
                  DTOutput("dataPreview")
                )
              )
      ),
      # Meta-Analysis Tab
      tabItem(tabName = "meta",
              fluidRow(
                bs4Card(
                  title = "Meta-Analysis Summary",
                  width = 12,
                  status = "success",
                  solidHeader = TRUE,
                  verbatimTextOutput("metaSummary")
                )
              )
      ),
      # Forest Plot Tab (Full width, options below)
      tabItem(tabName = "plot",
              fluidRow(
                bs4Card(
                  title = "Forest Plot",
                  width = 12,
                  status = "warning",
                  solidHeader = TRUE,
                  plotOutput("forestPlot", height = "600px")
                )
              ),
              fluidRow(
                bs4Card(
                  title = "Forest Plot Options",
                  width = 12,
                  status = "info",
                  solidHeader = TRUE,
                  sliderInput("forestXRange", "X-axis Range", min = -5, max = 5, value = c(-2, 2), step = 0.1),
                  numericInput("forestCex", "Study Label Size", value = 1, min = 0.5, max = 3, step = 0.1),
                  checkboxInput("showLabels", "Show Study Names", value = TRUE)
                )
              )
      ),
      # Funnel Plot Tab with Options
      tabItem(tabName = "funnel",
              fluidRow(
                bs4Card(
                  title = "Funnel Plot",
                  width = 12,
                  status = "warning",
                  solidHeader = TRUE,
                  plotOutput("funnelPlot", height = "600px")
                )
              ),
              fluidRow(
                bs4Card(
                  title = "Funnel Plot Options",
                  width = 12,
                  status = "info",
                  solidHeader = TRUE,
                  checkboxInput("funnelContours", "Add Contour Lines", value = FALSE),
                  sliderInput("funnelLevel", "Contour Level (%)", min = 80, max = 99, value = 95, step = 1),
                  numericInput("funnelCex", "Funnel Plot Text Size", value = 1, min = 0.5, max = 3, step = 0.1)
                )
              )
      ),
      # Baujat Plot Tab
      tabItem(tabName = "baujat",
              fluidRow(
                bs4Card(
                  title = "Baujat Plot",
                  width = 12,
                  status = "warning",
                  solidHeader = TRUE,
                  plotOutput("baujatPlot", height = "600px")
                )
              )
      ),
      # Influence Diagnostics Tab
      tabItem(tabName = "influence",
              fluidRow(
                bs4Card(
                  title = "Influence Diagnostics",
                  width = 12,
                  status = "warning",
                  solidHeader = TRUE,
                  plotOutput("influencePlot", height = "600px")
                )
              )
      ),
      # Cumulative Analysis Tab
      tabItem(tabName = "cumulative",
              fluidRow(
                bs4Card(
                  title = "Cumulative Meta-Analysis",
                  width = 12,
                  status = "primary",
                  solidHeader = TRUE,
                  plotOutput("cumulativePlot", height = "600px")
                )
              )
      ),
      # Moderator Analysis Tab (Text output and scatter plot)
      tabItem(tabName = "moderator",
              fluidRow(
                bs4Card(
                  title = "Moderator Analysis",
                  width = 12,
                  status = "primary",
                  solidHeader = TRUE,
                  uiOutput("moderatorSelectUI"),
                  actionButton("runMod", "Run Meta-Regression", icon = icon("play")),
                  br(), br(),
                  verbatimTextOutput("moderatorOutput"),
                  plotOutput("moderatorScatter", height = "600px")
                )
              )
      ),
      # Radial Plot Tab
      tabItem(tabName = "radial",
              fluidRow(
                bs4Card(
                  title = "Radial Plot",
                  width = 12,
                  status = "primary",
                  solidHeader = TRUE,
                  plotOutput("radialPlot", height = "600px")
                )
              )
      ),
      # Leave-One-Out Analysis Tab
      tabItem(tabName = "leaveone",
              fluidRow(
                bs4Card(
                  title = "Leave-One-Out Analysis",
                  width = 12,
                  status = "primary",
                  solidHeader = TRUE,
                  verbatimTextOutput("leaveOneOutText"),
                  plotOutput("leaveOneOutPlot", height = "600px")
                )
              )
      ),
      # Cook's Distance Tab
      tabItem(tabName = "cooks",
              fluidRow(
                bs4Card(
                  title = "Cook's Distance Plot",
                  width = 12,
                  status = "primary",
                  solidHeader = TRUE,
                  plotOutput("cooksPlot", height = "600px")
                )
              )
      ),
      # Data Summary Tab
      tabItem(tabName = "datasummary",
              fluidRow(
                bs4Card(
                  title = "Data Summary",
                  width = 12,
                  status = "info",
                  solidHeader = TRUE,
                  verbatimTextOutput("dataSummary")
                )
              )
      ),
      # Additional Results Tab
      tabItem(tabName = "additional",
              fluidRow(
                bs4Card(
                  title = "Additional Meta-Analysis Results",
                  width = 12,
                  status = "info",
                  solidHeader = TRUE,
                  verbatimTextOutput("additionalResults")
                )
              )
      ),
      # Download Samples Tab
      tabItem(tabName = "downloads",
              fluidRow(
                bs4Card(
                  title = "Download Sample CSV Files",
                  width = 12,
                  status = "info",
                  solidHeader = TRUE,
                  selectInput("sampleType", "Select Conversion Type:",
                              choices = allConvTypes, selected = "Mean & SE"),
                  downloadButton("downloadSample", "Download Sample CSV"),
                  br(), br(),
                  downloadButton("downloadAllSamples", "Download Combined Sample CSV"),
                  br(), br(),
                  h4("Data Explanation:"),
                  HTML("<p><strong>Study:</strong> Identifier for the study.</p>
                        <p><strong>ConversionType:</strong> Format of the reported data (e.g., 'Mean & SE', 'Reg_B', etc.).</p>
                        <p><strong>grp1m, grp1se, grp1n:</strong> Mean, standard error, and sample size for group 1.</p>
                        <p><strong>grp2m, grp2se, grp2n:</strong> Mean, standard error, and sample size for group 2.</p>
                        <p><strong>b, beta, sdy:</strong> Regression parameters and standard deviation for converting regression coefficients.</p>
                        <p><strong>rpb:</strong> Point biserial correlation value.</p>
                        <p><strong>f, t:</strong> F-statistic and t-statistic for one-way ANOVA and two-sample t-tests.</p>
                        <p><strong>effect_size, p, N:</strong> Effect size value, p-value, and sample size used in conversion from p-value.</p>
                        <p><strong>effectSizeType:</strong> Type of effect size (e.g., 'd', 'g', 'difference').</p>
                        <p><strong>chisq, totaln:</strong> Chi-squared statistic and total sample size for chi-squared conversion.</p>
                        <p><strong>esType:</strong> The effect size metric to be computed (e.g., 'd', 'g', 'cox.or').</p>")
                )
              )
      ),
      # Instructions Tab
      tabItem(tabName = "instructions",
              fluidRow(
                bs4Card(
                  title = "Instructions & Sample CSV",
                  width = 12,
                  status = "info",
                  solidHeader = TRUE,
                  HTML("
              <h4>Meta-Analysis Calculator Instructions</h4>
              <p>This app lets you upload a CSV file containing raw study data from various reporting formats.
              The app converts each study’s raw data into a common effect size (yi) and its variance (vi) using functions from the <code>esc</code> package, and then performs a random‐effects meta‐analysis using <code>rma()</code> from <code>metafor</code>.</p>
              <p><strong>Required CSV Format:</strong></p>
              <p>Your CSV file must have the following header (exactly):</p>
              <pre>
Study,ConversionType,grp1m,grp1se,grp1n,grp2m,grp2se,grp2n,b,beta,sdy,rpb,f,t,effect_size,p,N,effectSizeType,chisq,totaln,esType
              </pre>
              <p>Fill in only the columns needed for the given conversion type; leave the others blank or as NA.</p>
              <p><strong>Sample CSV File:</strong></p>
              <pre>
Study,ConversionType,grp1m,grp1se,grp1n,grp2m,grp2se,grp2n,b,beta,sdy,rpb,f,t,effect_size,p,N,effectSizeType,chisq,totaln,esType
Study1,Mean & SE,8.5,1.5,50,11,1.8,60,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,d
Study2,Reg_B,NA,NA,100,NA,NA,150,3.3,NA,5,NA,NA,NA,NA,NA,NA,NA,NA,NA,d
Study3,Reg_beta,NA,NA,100,NA,NA,150,NA,0.32,5,NA,NA,NA,NA,NA,NA,NA,NA,NA,d
Study4,PointBiserial,NA,NA,99,NA,NA,120,NA,NA,NA,0.25,NA,NA,NA,NA,NA,NA,NA,NA,d
Study5,OneWayANOVA,NA,NA,519,NA,NA,528,NA,NA,NA,NA,5.04,NA,NA,NA,NA,NA,NA,NA,g
Study6,TwoSample_t,NA,NA,100,NA,NA,150,NA,NA,NA,NA,NA,3.3,NA,NA,NA,NA,NA,NA,d
Study7,p_to_SE,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0.71,0.013,71,NA,NA,NA,NA,difference
Study8,ChiSquared,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,7.9,100,cox.or
              </pre>
              <p>Save the text above as <em>sample_data.csv</em> and upload it to test the app.</p>
            ")
                )
              )
      )
    )
  ),
  footer = bs4DashFooter("Meta-Analysis Calculator © 2025")
)

### Server ###
server <- function(input, output, session) {
  
  # Reactive: Read uploaded CSV data (with proper NA conversion)
  rawData <- reactive({
    req(input$datafile)
    data <- read.csv(input$datafile$datapath, stringsAsFactors = FALSE, na.strings = c("NA", ""))
    message("Raw CSV data loaded:")
    print(head(data))
    return(data)
  })
  
  # Display a preview of the uploaded data
  output$dataPreview <- renderDT({
    req(rawData())
    datatable(rawData(), options = list(pageLength = 5, autoWidth = TRUE))
  })
  
  # Process each row: compute yi and vi
  computedData <- reactive({
    req(rawData())
    df <- rawData()
    rows <- lapply(seq_len(nrow(df)), function(i) {
      row <- df[i, ]
      if (is.null(row$ConversionType) || is.na(row$ConversionType) || row$ConversionType == "") {
        message("Row ", i, ": No ConversionType provided, skipping.")
        return(NULL)
      }
      conv <- as.character(row$ConversionType)
      study <- as.character(row$Study)
      if (is.na(study) || study == "") study <- paste0("Row_", i)
      message("Processing Study: ", study, " with ConversionType: ", conv)
      params <- as.list(row[!(names(row) %in% c("Study", "ConversionType"))])
      params <- lapply(params, maybe_numeric)
      es <- tryCatch({
        computeES(conv, params)
      }, error = function(e) {
        message("Error computing effect size for Study: ", study, " : ", e)
        return(list(yi = NA, vi = NA))
      })
      message("Computed for Study: ", study, " -> yi: ", es$yi, " ; vi: ", es$vi)
      data.frame(Study = study,
                 ConversionType = conv,
                 yi = ifelse(length(es$yi)==0, NA, es$yi),
                 vi = ifelse(length(es$vi)==0, NA, es$vi),
                 stringsAsFactors = FALSE)
    })
    rows <- Filter(Negate(is.null), rows)
    if (length(rows) == 0) {
      message("No valid rows found in computedData.")
      return(data.frame(Study = character(), ConversionType = character(), yi = numeric(), vi = numeric(), stringsAsFactors = FALSE))
    }
    dat <- do.call(rbind, rows)
    dat$yi <- as.numeric(dat$yi)
    dat$vi <- as.numeric(dat$vi)
    message("Final computedData:")
    print(dat)
    return(dat)
  })
  
  # Run meta-analysis using computed effect sizes
  metaResult <- reactive({
    dat <- computedData()
    message("Inside metaResult; computedData has ", nrow(dat), " rows.")
    req(nrow(dat) > 0)
    # Exclude rows with NA in yi or vi
    validData <- dat[!is.na(dat$yi) & !is.na(dat$vi), ]
    message("Data used for meta-analysis (n = ", nrow(validData), "):")
    print(validData)
    req(nrow(validData) > 0)
    res <- rma(yi, vi, data = validData, method = "REML")
    message("Meta-analysis result computed:")
    print(res)
    return(list(result = res, validData = validData))
  })
  
  output$metaSummary <- renderPrint({
    message("Rendering metaSummary output")
    req(metaResult())
    res <- metaResult()$result
    summary(res)
  })
  
  output$forestPlot <- renderPlot({
    message("Rendering forest plot")
    req(metaResult())
    validData <- metaResult()$validData
    xRange <- input$forestXRange
    cexVal <- input$forestCex
    slabText <- if(input$showLabels) validData$Study else rep("", nrow(validData))
    forest(metaResult()$result, slab = slabText, xlim = xRange, cex = cexVal)
  })
  
  output$funnelPlot <- renderPlot({
    message("Rendering funnel plot")
    req(metaResult())
    res <- metaResult()$result
    if(input$funnelContours){
      funnel(res, cex = input$funnelCex, contour = TRUE, level = input$funnelLevel/100, legend = TRUE)
    } else {
      funnel(res, cex = input$funnelCex)
    }
  })
  
  output$baujatPlot <- renderPlot({
    message("Rendering Baujat plot")
    req(metaResult())
    baujat(metaResult()$result)
  })
  
  output$influencePlot <- renderPlot({
    message("Rendering influence diagnostics plot")
    req(metaResult())
    infl <- influence(metaResult()$result)
    plot(infl)
  })
  
  output$cumulativePlot <- renderPlot({
    message("Rendering cumulative meta-analysis plot")
    req(metaResult())
    cumulRes <- cumul(metaResult()$result)
    forest(cumulRes, main="Cumulative Meta-Analysis")
  })
  
  # Moderator Analysis: Text output and scatter plot (no bubble plot)
  output$moderatorSelectUI <- renderUI({
    req(rawData())
    numCols <- names(rawData())[sapply(rawData(), is.numeric)]
    if (length(numCols) == 0) {
      h4("No numeric moderator variables available in the data.")
    } else {
      selectInput("moderatorVar", "Select Moderator Variable:", choices = numCols)
    }
  })
  
  moderatorRes <- eventReactive(input$runMod, {
    req(rawData(), input$moderatorVar)
    dat <- rawData()
    if (!(input$moderatorVar %in% names(dat))) {
      return(NULL)
    }
    dat[[input$moderatorVar]] <- as.numeric(dat[[input$moderatorVar]])
    compData <- computedData()
    compData[[input$moderatorVar]] <- dat[[input$moderatorVar]][1:nrow(compData)]
    validData <- compData[!is.na(compData$yi) & !is.na(compData$vi), ]
    if(nrow(validData) < 2) return(NULL)
    res <- rma(yi, vi, mods = ~ get(input$moderatorVar), data = validData, method = "REML")
    return(list(result = res, validData = validData))
  })
  
  output$moderatorOutput <- renderPrint({
    req(moderatorRes())
    cat("Meta-Regression Output:\n")
    print(summary(moderatorRes()$result))
  })
  
  output$moderatorScatter <- renderPlot({
    req(moderatorRes())
    dat <- moderatorRes()$validData
    modVar <- input$moderatorVar
    validRows <- !is.na(dat[[modVar]]) & is.finite(dat[[modVar]]) & !is.na(dat$yi) & is.finite(dat$yi)
    if(sum(validRows) < 2){
      plot.new()
      text(0.5, 0.5, "Not enough valid data to compute regression line.", cex = 1.5)
    } else {
      lm_fit <- lm(yi ~ dat[[modVar]], data = dat[validRows,])
      plot(dat[[modVar]][validRows], dat$yi[validRows], xlab = modVar, ylab = "Effect Size",
           main = paste("Scatter Plot:", modVar, "vs Effect Size"))
      abline(lm_fit, col = "red")
    }
  })
  
  output$radialPlot <- renderPlot({
    message("Rendering radial plot")
    req(metaResult())
    radial(metaResult()$result)
  })
  
  output$leaveOneOutText <- renderPrint({
    message("Rendering leave-one-out text output")
    req(metaResult())
    loo <- leave1out(metaResult()$result)
    print(loo)
  })
  
  output$leaveOneOutPlot <- renderPlot({
    message("Rendering leave-one-out plot")
    req(metaResult())
    loo <- leave1out(metaResult()$result)
    plot(loo$estimate, type = "b", xlab = "Study Removed", 
         ylab = "Effect Size Estimate", main = "Leave-One-Out Analysis")
    abline(h = metaResult()$result$b, col = "red", lty = 2)
  })
  
  output$cooksPlot <- renderPlot({
    message("Rendering Cook's Distance plot")
    req(metaResult())
    inf <- influence(metaResult()$result)
    cooks <- inf$cooks.distance
    if(length(cooks) == 0){
      plot.new()
      text(0.5, 0.5, "No Cook's distance values available.", cex = 1.5)
    } else {
      plot(cooks, type = "h", main = "Cook's Distance", xlab = "Study", ylab = "Cook's Distance")
      abline(h = 4/length(cooks), col = "red", lty = 2)
    }
  })
  
  output$dataSummary <- renderPrint({
    req(computedData())
    cat("Summary of Computed Data:\n")
    print(summary(computedData()))
    cat("\nHead of Computed Data:\n")
    print(head(computedData()))
  })
  
  output$additionalResults <- renderPrint({
    message("Rendering additional meta-analysis results")
    req(metaResult())
    res <- metaResult()$result
    cat("Detailed Meta-Analysis Diagnostics:\n")
    cat("QE (heterogeneity Q): ", round(res$QE, 3), "\n")
    cat("df: ", res$k - 1, "\n")
    cat("p-value for Q: ", round(res$QEp, 3), "\n")
    cat("I-squared: ", round(res$I2, 1), "%\n")
    cat("Tau-squared: ", round(res$tau2, 3), "\n")
    predInt <- predict(res, transf = FALSE)
    cat("\nPrediction Interval:\n")
    print(predInt)
    egger <- tryCatch({
      regtest(res, model = "lm")
    }, error = function(e) {
      "Egger's test could not be computed."
    })
    cat("\nEgger's Test:\n")
    print(egger)
    cat("\nFull meta-analysis result:\n")
    print(res)
  })
  
  output$downloadSample <- downloadHandler(
    filename = function() {
      paste0(gsub(" ", "_", input$sampleType), "_sample.csv")
    },
    content = function(file) {
      sampleData <- generateSampleCSV(input$sampleType)
      write.csv(sampleData, file, row.names = FALSE)
    }
  )
  
  output$downloadAllSamples <- downloadHandler(
    filename = function() {
      "all_samples.csv"
    },
    content = function(file) {
      allSamples <- do.call(rbind, lapply(allConvTypes, generateSampleCSV))
      write.csv(allSamples, file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
