library(dplyr)
library(DT)
library(data.table)
library(plotly)
library(reshape2)
options(shiny.maxRequestSize=100*1024^2) 

server <- function(input, output, session) {
  dsnames <- c()
  
  dat <- reactive({
    req(input$file1)
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    dat <- read.table(inFile$datapath, sep="\t", header=T)
  })
  
  cent_dat <- reactive({
    req(input$file2)
    inFile <- input$file2
    if (is.null(inFile))
      return(NULL)
    cent_dat <- read.table(inFile$datapath, sep="\t", header=F)
  })
  
 output$contents <- renderDT({
    # If missing input, return to avoid error later in function
    if(is.null(dat()))
      return(NULL)

    datatable(data.frame(dat()), extensions = 'FixedColumns', filter = 'bottom', options = list(
      dom = 'RMDCT<"clear">lfrtip',
      searchHighlight = TRUE,
      pageLength = 5,
      lengthMenu = list(c(5, 10, 15, 20, 25, -1), c('5', '10', '15', '20', '25', 'All')),
      initComplete = JS("function(settings, json) {",
                        "$(this.api().table().header()).css({'background-color': '#005ab3', 'color': '#fff'});",
                        "}"),
      scrollX = TRUE,
      scrollY = TRUE,
      fixedColumns = list(leftColumns = 2)
    ), rownames=F)
  })

  
  output$choose_xaxis <- renderUI({
    # If missing input, return to avoid error later in function
    if(is.null(dat()))
      return(NULL)
    
    # Get the data set with the appropriate name
    nums <- sapply(dat(), is.numeric)
    df <- dat()[ , nums]
    colnames <- colnames(df)
    
    # Create the checkboxes and select them all by default
    selectInput("xaxis", "Choose X-axis", 
                       choices  = colnames)
  })

  output$choose_yaxis <- renderUI({
    # If missing input, return to avoid error later in function
    if(is.null(dat()))
      return(NULL)
    
    # Get the data set with the appropriate name
    nums <- sapply(dat(), is.numeric)
    df <- dat()[ , nums]
    colnames <- colnames(df)
    
    # Create the checkboxes and select them all by default
    selectInput("yaxis", "Choose Y-axis", 
                choices  = colnames)
  })
  
  output$choose_variable <- renderUI({
    # If missing input, return to avoid error later in function
    if(is.null(dat()))
      return(NULL)
    
    # Get the data set with the appropriate name
    nums <- sapply(dat(), is.numeric)
    df <- dat()[ , nums]
    colnames <- colnames(df)
    
    # Create the checkboxes and select them all by default
    selectInput("variable", "Choose Variable", 
                choices  = colnames)
  })
  
  output$choose_color <- renderUI({
    # If missing input, return to avoid error later in function
    if(is.null(dat()))
      return()
    
    # Get the data set with the appropriate name
    
    colnames <- colnames(dat())
    
    # Create the checkboxes and select them all by default
    selectInput("color", "Choose Coloring", 
                choices  = colnames)
  })
  
  output$choose_group <- renderUI({
    # If missing input, return to avoid error later in function
    if(is.null(dat()))
      return(NULL)
    
    # Get the data set with the appropriate name
    nums <- sapply(dat(), is.numeric)
    df <- dat()[ , !(nums)]
    colnames <- colnames(df)
    
    # Create the checkboxes and select them all by default
    selectInput("group", "Choose Grouping", 
                choices  = colnames)
  })
  
  output$choose_plottype <- renderUI({
    # If missing input, return to avoid error later in function
    if(is.null(dat()))
      return()
    
    # Get the data set with the appropriate name
    
    colnames <- c("boxplot", "histogram")
    
    # Create the checkboxes and select them all by default
    selectInput("plottype", "Choose Plot Type", 
                choices  = colnames)
  })
  
  output$choose_tax_class <- renderUI({
    # If missing input, return to avoid error later in function
    if(is.null(cent_dat()))
      return()
    
    # Get the data set with the appropriate name
    tax_classes = unique(cent_dat()$V2)
    
    # Create the checkboxes and select them all by default
    selectInput("tax_class", "Choose Taxonomic Class", 
                choices  = tax_classes)
  })    
  
  output$mainPlot <- renderPlotly({
    req(input$xaxis)
    req(input$yaxis)
    req(input$color)
    df <- dat()
    xvalue <- df[[input$xaxis]]
    yvalue <- df[[input$yaxis]]
    colorvalue <- df[[input$color]]
    if (is.null(df) | is.null(xvalue) | is.null(yvalue) | is.null(colorvalue)) {
      return(NULL)
    }
    if (input$overview_mode) {
      # Get the data set with the appropriate name
      nums <- sapply(dat(), is.numeric)
      df <- dat()[ , nums]
      colnames <- colnames(df)
      df$sequencing_file <- dat()$sample_id
      melted.df <- melt(df, measure.vars=colnames, variable.name='quantitative_metric')
      p <- ggplot(melted.df, aes(x=value)) + geom_histogram(alpha=0.7) + theme_bw() + facet_wrap(~quantitative_metric, ncol=3, scales = "free") + scale_x_log10()
      p <- ggplotly(p) %>% layout(autosize=T)
    } else {
      xvalue_nona <- xvalue[!is.na(xvalue)]
      yvalue_nona <- yvalue[!is.na(yvalue)]
      if (input$fit_lm) {
        fit <- lm(yvalue ~ xvalue)
        rsq <- as.character(summary(fit)$r.squared)
        maxx = max(xvalue)
        maxy = max(yvalue)
        p <- plot_ly(x = ~xvalue) %>% add_markers(y = ~yvalue, color = ~colorvalue, hoverinfo = 'text', text = ~paste('Sequencing File: ', df$sample_id))
        p <- p %>% add_lines(x = ~xvalue, y = fitted(fit)) %>% add_annotations(x= maxx, y= maxy, text = ~paste("R-Squared = ", rsq, sep=""), showarrow=F)
      } else {
        p <- plot_ly(x = ~xvalue) %>% add_markers(y = ~yvalue, color = ~colorvalue, hoverinfo = 'text', text = ~paste('Sequencing File: ', df$sample_id))
      }
    }
  })
  
  output$projPlot <- renderPlotly({
    req(input$group)
    req(input$variable)
    req(input$plottype)
    df <- dat()
    variable.vec <- df[[input$variable]]
    group.vec <- df[[input$group]]
    if (is.null(df) | is.null(variable.vec) | is.null(group.vec)) {
      return(NULL)
    }
    if (length(unique(group.vec)) > 10) {
      print("HELLO")
      blank.df <- data.frame()
      p <- ggplot(blank.df) + geom_point() + xlim(-1, 1) + ylim(-1, 1) + annotate('text', label='Grouping variable contains more than 10 groups! Unable to process.', x=0, y=0) + theme_void()
      p <- ggplotly(p) %>% layout(autosize=T)
    }
    else {
      new.df <- data.frame(variable=variable.vec, grouping=group.vec)
      if (input$facet) {
          ggplot(new.df, aes(x=variable)) + geom_histogram(alpha=0.7) + theme_bw() + facet_wrap(~grouping, ncol=3, scales = "free")
      } else {
        if (input$plottype == 'boxplot') {
          p <- ggplot(new.df, aes(x=grouping, y=variable, fill=grouping)) + geom_boxplot(alpha=0.7) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
          p <- ggplotly(p) %>% layout(autosize=T)
        } else if (input$plottype == 'histogram') {
          p <- ggplot(new.df, aes(x=variable, fill=grouping)) + geom_histogram(color="black") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
          p <- ggplotly(p) %>% layout(autosize=T)
        }
      }
    }
  })
  
  output$centPlot <- renderPlotly({
    df <- cent_dat()
    filt.df <- df[df$V2 == input$tax_class,]
    filt.mat <- acast(filt.df, V1~V3, value.var="V4")
    if (input$quantile_mode) {
      quantiles = unique(quantile(filt.df$V4, seq(0,1,0.1)))
      filt.df$V5 <- cut(filt.df$V4, breaks = c(quantiles[1]-1, quantiles[-1]), right=T, labels=F)
      filt.mat <- acast(filt.df, V1~V3, value.var="V5")
    }
    p <- plot_ly(x = colnames(filt.mat), y = rownames(filt.mat), z = filt.mat, type = "heatmap") %>% layout(xaxis=list(title="Taxonomic Clades", showticklabels=F), yaxis=list(title="Samples", showticklabels=F))
  })
}