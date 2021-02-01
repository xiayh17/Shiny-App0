library(shiny)
library(shinycssloaders)
library(shinydashboard)
library(DT)
library(collapsibleTree)
library(shinyBS)
library(AnnoProbe)
library(GEOquery)
library(stringr)
library(limma)
library(markdown)
#source("helper.R")
ui <- fluidPage(

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
          dataTableOutput("preview1") %>% withSpinner(type = 6),
          h3(em(strong("A boxplot check"))),
          br(),
          helpText("Quick check expression data just downloaded with boxplot"),
          actionButton("applyBoxpot1", "Box plot"),
          plotOutput("boxplot1") %>% withSpinner(type = 6)
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
          tableOutput("group") %>% withSpinner(type = 6),
          h3(em(strong("Preview Clinic infomation"))),
          helpText("Preview phenotypic data after download"),
          dataTableOutput("preview2") %>% withSpinner(type = 6)
        ),

        tabItem(
          # probe Annotation section
          tabName = "probe",
          h3(em(strong("Input your type"))),
          helpText("source of probe anntation stored, one of 'pipe', 'bioc', 'soft', default:'pipe'"),
          textInput("type","","bioc"),
          actionButton("anno","Start Probe Annotation"),
          h3(em(strong("Preview Probe Anntation"))),
          dataTableOutput("preview3") %>% withSpinner(type = 6)
        ),

        tabItem(
          # Filter expression data
          tabName = "filter",
          h3(em(strong("Filter expression matrix based on annotation"))),
          actionButton("filter","Start Filter"),
          br(),
          dataTableOutput("preview4") %>% withSpinner(type = 6)
        ),

        tabItem(
          # Normalization section
          tabName = "normal",
          h3(em(strong("Normalization with limma"))),
          actionButton("norm","Start Normalization"),
          dataTableOutput("preview5")  %>% withSpinner(type = 6)
        ),

        tabItem(
          # heatmap section
          tabName = "heatmap",
          h3(em(strong("Plot Heatmap"))),
          actionButton("ph","Plot Heatmap"),
          plotOutput("hplot")  %>% withSpinner(type = 6),
          h3(em(strong("Plot Volcano"))),
          br(),
          textInput("style","Style of Volcano",1),
          helpText("you can try 1 or 2"),
          br(),
          textInput("p_thred","P Thred",0.05),
          textInput("logFC_thred","logFC Thred",1),
          actionButton("pv","Plot Volcano"),
          plotOutput("vplot")  %>% withSpinner(type = 6)
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
    gse=geoChina(input$geoacc)
    gse[[1]]
  })

  # access the expression of assay data
  probes_expr <- reactive({
    exprs(eSet())
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

  # boxplot1
  observeEvent(input$applyBoxpot1, {
    output$boxplot1 <- renderPlot({
      boxplot(probes_expr(),las=2)
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
      # preview data
      output$preview3 <- renderDataTable({
        probes_anno()
      })
    })

    #4. Filter Expression
    # annotate
    genes_expr <- reactive({
      dat1 <- probes_expr()
      dat2 <- probes_anno()
      filterEM(dat1,dat2)
    })

    observeEvent(input$filter, {
      # preview data
      output$preview4 <- renderDataTable({
        genes_expr()
      })
    })

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
      output$preview5 <- renderDataTable({
        deg()
      })
    })

    # heatmap
    observeEvent(input$ph,{
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
