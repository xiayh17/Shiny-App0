library(shiny)
library(shinycssloaders)
library(shinydashboard)
library(DT)
#library(collapsibleTree)
library(shinyBS)
library(AnnoProbe)
library(GEOquery)
library(stringr)
library(limma)
library(markdown)
library(shinyjs)
library(plotly)
library(tidyr)
createLink <- function(val) {
  sprintf('<a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s" target="_blank" class="btn btn-primary">GeneCards</a>',val)
}
#source("helper.R")
ui <- fluidPage(
  useShinyjs(),
  # load custom stylesheet
  includeCSS("www/style.css"),
  includeCSS("www/iconfont.css"),
  # load google analytics script
  # tags$head(includeScript("www/google-analytics-bioNPS.js")),
  # custom title and ico of web tab
  tags$head(
    tags$title('CP2S'),
    tags$link(rel = 'shortcut icon',
              href='biotree.ico',
              type='image/x-icon')
  ),
  # remove shiny "red" warning messages on GUI
  tags$style(type="text/css",
             ".shiny-output-error { visibility: hidden; }",
             ".shiny-output-error:before { visibility: hidden; }"
  ),

  # load page layout
  dashboardPage(
    # color of dashboard
    skin = "green",
    # title and size of dashboard head
    dashboardHeader(title = "Converting Probes into Symbols", titleWidth = 300),
    # sidebar of web
    dashboardSidebar(width = 300,
                     sidebarMenu(
                       # add link, picture and text
                       HTML(paste0(
                         "<br>",
                         "<a href='https://mp.weixin.qq.com/s/OPhXbBJQC-3gQ3dr5tLZ7Q' target='_blank'><img style = 'display: block; margin-left: auto; margin-right: auto;' src='log.png' width = '186'></a>",
                         "<br>",
                         "<p style = 'text-align: center;'><small><a href='https://www.github.com//xiayh17' target='_blank'>Probes2Symbol logo designer</a></small></p>",
                         "<br>"
                       )),
                       # add menu of sidebar
                       menuItem("Home", tabName = "home", icon = icon("robot")),
                       menuItem("Download Data", tabName = "download", icon = icon("cloud-download-alt")),
                       menuItem("Clinic infomation", tabName = "pdata", icon = icon("diagnoses")),
                       menuItem("Probe Annotation", tabName = "probe", icon = icon("random", lib = "glyphicon")),
                       menuItem("Filter expression", tabName = "filter", icon = icon("stats", lib = "glyphicon")),
                       menuItem("Normalization", tabName = "normal", icon = icon("align-right")),
                       menuItem("Heatmap", tabName = "heatmap", icon = icon("map marked alt")),
                       menuItem("Release", tabName = "release", icon = icon("code-branch")),
                       HTML(paste0(
                         "<br><br><br><br><br><br><br><br><br>",
                         "<table style='margin-left:auto; margin-right:auto;'>",
                         "<tr>",
                         "<td style='padding: 5px;'><a href='https://space.bilibili.com/338686099' target='_blank'><i class='iconfont iconbilibili'></i></a></td>",
                         "<td style='padding: 5px;'><a href='https://www.youtube.com/channel/UC67sImqK7V8tSWHMG8azIVA' target='_blank'><i class='fab fa-youtube fa-lg'></i></a></td>",
                         "<td style='padding: 5px;'><a href='weixin.png' target='_blank'><i class='iconfont iconweixin'></i></a></td>",
                         "</tr>",
                         "</table>",
                         "<br>"),
                         HTML(paste0(
                           "<script>",
                           "var today = new Date();",
                           "var yyyy = today.getFullYear();",
                           "</script>",
                           "<p style = 'text-align: center;'><small>&copy; - <a href='http://www.biotrainee.com' target='_blank'>biotrainee.com</a> - <script>document.write(yyyy);</script></small></p>")
                         ))
                     )

    ), # end dashboardSidebar

    dashboardBody(
      tabItems(

        # home section
        tabItem(tabName = "home",includeMarkdown("www/home.md")),

        # download section
        tabItem(
          tabName = "download",
          # a panel to input GEO accession and start download
          h3(em(strong("Input GEO Accession"))),
          br(),
          helpText("Input a GEO Accession.
                    Data will be collected with geoChina"),
          textInput("geoacc", "", "GSE1009"),
          actionButton("applyDownload", "Start Download"),
          h3(em(strong("A preview of data after download success"))),
          br(),
          helpText("Preview expression data just downloaded"),
          htmlOutput("dim"),
          br(),
          #dataTableOutput("preview1") %>% withSpinner(type = 6) %>% hidden(),
          hidden(div(id = 'show1', withSpinner((dataTableOutput("preview1")),type = 6))),
          br(),
          downloadLink("downloadpreview1", "Download Table"),
          h3(em(strong("A boxplot check"))),
          br(),
          helpText("Quick check expression data just downloaded with boxplot"),
          actionButton("applyBoxpot1", "Box plot"),
          checkboxInput("log", "log gene expression",
                    value = FALSE),
          #plotOutput("boxplot1") %>% withSpinner(type = 6) %>% hidden()
          hidden(div(id = 'show2', withSpinner((plotlyOutput("boxplot1")),type = 6))),
          br(),
          #downloadLink("downloadboxplot1", "Download Plot"),
          br(),
          br()
        ),

        tabItem(
          # clinic data section
          tabName = "pdata",
          #dataTableOutput("clinicDataTable") %>% withSpinner(color = "green")
          h3(em(strong("Group of Clinic infomation"))),
          textInput("gc","Which colunm to Group","title"),
          textInput("g1","Group1 keywords","Control"),
          textInput("g2","Group2 Keywords","Diabetes"),
          actionButton("applyGroup","Access Group"),
          helpText("The robot point out Group of Clinic infomation below"),
          #tableOutput("group") %>% withSpinner(type = 6) %>% hidden(),
          hidden(div(id = 'show3', withSpinner((tableOutput("group")),type = 6))),
          h3(em(strong("Preview Clinic infomation"))),
          helpText("Preview phenotypic data after download"),
          dataTableOutput("preview2") %>% withSpinner(type = 6)
        ),

        tabItem(
          # probe Annotation section
          tabName = "probe",
          h3(em(strong("Select your type"))),
          helpText("source of probe anntation stored, one of 'pipe', 'bioc', 'soft'"),
          #textInput("type","","bioc"),
          selectInput(inputId = "type",
                      label = "Select a type:",
                      choices = c('pipe', 'bioc', 'soft'),
                      selected = 'bioc'),
          helpText("choose human or mouse, or rat, default: human"),
          selectInput(inputId = "species",
                      label = "Choose species:",
                      choices = c('human', 'mouse', 'rat'),
                      selected = 'human'),
          actionButton("anno","Start Probe Annotation"),
          h3(em(strong("Preview Probe Anntation"))),
          #dataTableOutput("preview3") %>% withSpinner(type = 6) %>% hidden()
          hidden(div(id = 'show4', withSpinner((dataTableOutput("preview3")),type = 6)))
        ),

        tabItem(
          # Filter expression data
          tabName = "filter",
          h3(em(strong("Filter expression matrix based on annotation"))),
          #actionButton("filter","Start Filter"),
          br(),
          dataTableOutput("preview4") %>% withSpinner(type = 6)
        ),

        tabItem(
          # Normalization section
          tabName = "normal",
          h3(em(strong("Normalization with limma"))),
          actionButton("norm","Start Normalization"),
          #dataTableOutput("preview5")  %>% withSpinner(type = 6) %>% hidden()
          hidden(div(id = 'show5', withSpinner((dataTableOutput("preview5")),type = 6))),
          downloadButton("downloadpreview5", "Download Table")
        ),

        tabItem(
          # heatmap section
          tabName = "heatmap",
          h3(em(strong("Plot Heatmap"))),
          actionButton("ph","Plot Heatmap"),
          #plotOutput("hplot")  %>% withSpinner(type = 6),
          hidden(div(id = 'show6', withSpinner((plotOutput("hplot")),type = 6))),
          h3(em(strong("Plot Volcano"))),
          br(),
          textInput("style","Style of Volcano",1),
          helpText("you can try 1 or 2"),
          br(),
          textInput("p_thred","P Thred",0.05),
          textInput("logFC_thred","logFC Thred",1),
          actionButton("pv","Plot Volcano"),
          #plotOutput("vplot")  %>% withSpinner(type = 6) %>% hidden()
          hidden(div(id = 'show7', withSpinner((plotOutput("vplot")),type = 6)))
        ),

        tabItem(tabName = "release", includeMarkdown("www/releases.md"))

      )

    ) # end dashboardBody

  )# end dashboardPage

)

server <- function(input, output){
  ## 1. download section
  # download data
  eSet <- eventReactive(input$applyDownload, {
    show("show1")
    gse=geoChina(input$geoacc)
    gse[[1]]
  })

  # access the expression of assay data
  probes_expr <- reactive({
    probes_expr <- exprs(eSet())
    if(input$log) {
      log(probes_expr)
    } else {
      probes_expr
    }
  })

  # preview data table
  output$preview1 <- renderDataTable({
      probes_expr()
  })

  # describe data dimension
  output$dim <- renderText({
    paste("Your data have",
          "<font color=\"#0275D8\"><b>",dim(probes_expr())[1], "</b></font>",
          "rows and",
          "<font color=\"#0275D8\"><b>",dim(probes_expr())[2], "</b></font>",
          "column")
  })

  #download the expression of assay data
  output$downloadpreview1 <- downloadHandler(
    filename = function() {
      paste("exp-data-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(probes_expr(), file)
    }
  )
  # boxplot1
  observeEvent(input$applyBoxpot1, {
    show("show2")
    data <- as.data.frame(probes_expr())
    res1 <- gather(data) #宽转长
    output$boxplot1 <- renderPlotly({
      #boxplot(probes_expr(),las=2)
      plot_ly(res1,y=~value,x=~key,color=~key,type = "box") %>%
      layout(
         xaxis = list(title = 'Samples'),
         yaxis = list(title = 'Expression'))
    })
  })

  # 2. clinic section
  # access the phenotypic data
  phenoDat <- reactive({
    pData(eSet())
  })

  # preview phenotypic
  output$preview2 <- renderDataTable({
    phenoDat()
  })

  # group after press button
  group_list <- reactive({
    g1 <- input$g1
    g2 <- input$g2
    gc <- input$gc
    ifelse(grepl(g1,phenoDat()[,gc]),g1,g2)
  })
  observeEvent(input$applyGroup, {
    show("show3")
    # group table
    output$group <- renderTable({
      test <- table(group_list())
      data <- as.data.frame(test)
      colnames(data) <- c("Group","Freq")
      data
    })

    #3. Probe Annotation section

    # annotate
    probes_anno <- reactive({
      GPL=eSet()@annotation
      idmap(GPL,type = input$type)
    })

    observeEvent(input$anno, {
      show("show4")
      # merge more info
      t2 <- reactive({
        load('GPL8300_bioc.rda')
        tmp=annoGene(GPL8300_bioc$symbol,'SYMBOL',input$species)
        t2=merge(tmp,GPL8300_bioc,by.y='symbol',by.x='SYMBOL')
        t2$links <- createLink(t2$SYMBOL)
        t2
      })
      # symbol  链接到genecard数据库
      # preview data
      output$preview3 <- renderDataTable({
        t2()
      }, escape = FALSE)
    })

    #4. Filter Expression
    # annotate
    genes_expr <- reactive({
      dat1 <- probes_expr()
      dat2 <- probes_anno()
      filterEM(dat1,dat2)
    })

    #observeEvent(input$filter, {

      # preview data
      output$preview4 <- renderDataTable({
        genes_expr()
      })
    #})

    #5. Normalization
    deg <- reactive({
      group_list <- group_list()
      genes_expr <- genes_expr()
      design=model.matrix(~factor(group_list))
      fit=lmFit(genes_expr,design)
      fit=eBayes(fit)
      topTable(fit,coef=2,n=Inf)
    })
    # preview after press button
    observeEvent(input$norm, {
      show("show5")
      output$preview5 <- renderDataTable({
        deg()
      })
    })
    #download the expression of assay data
    output$downloadpreview5 <- downloadHandler(
      filename = function() {
        paste("normalization-data-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.csv(deg(), file)
      }
    )

    # heatmap
    observeEvent(input$ph,{
      show("show6")
      output$hplot <- renderPlot({
        DEG <- deg()
        genes_expr <- genes_expr()
        group_list <- group_list()
        deg_heatmap(DEG, genes_expr, group_list, topn = 20)
      })
    })


    need_deg <- reactive({
      DEG <- deg()
      data.frame(symbols=rownames(DEG), logFC=DEG$logFC, p=DEG$P.Value)
    })

    observeEvent(input$pv,{
      show("show7")
      output$vplot <- renderPlot({
        need_deg <- need_deg()
        st <- as.numeric(input$style)
        p <- as.numeric(input$p_thred)
        logFC <- as.numeric(input$logFC_thred)
        deg_volcano(need_deg, style = st, p_thred = p, logFC_thred = logFC)
      })
    })
  })



}

shinyApp(ui, server)
