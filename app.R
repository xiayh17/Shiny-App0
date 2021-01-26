library(shiny)
library(shinycssloaders)
library(shinydashboard)
library(DT)
library(collapsibleTree)
# sudo su - -c "R -e \"install.packages(c('shinydashboard','shinycssloaders','DT','collapsibleTree'), repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')\""
ui <- fluidPage(

  # load custom stylesheet
  includeCSS("www/style.css"),
  includeCSS("www/iconfont.css"),
  # load google analytics script
  # tags$head(includeScript("www/google-analytics-bioNPS.js")),
  tags$head(
    tags$title('CP2S'),
    tags$link(rel = 'shortcut icon',
              href='biotree.ico',
              type='image/x-icon')
  ),
  # remove shiny "red" warning messages on GUI
  # tags$style(type="text/css",
  #            ".shiny-output-error { visibility: hidden; }",
  #            ".shiny-output-error:before { visibility: hidden; }"
  # ),

  # load page layout
  dashboardPage(

    skin = "green",

    dashboardHeader(title = "Converting Probes into Symbols", titleWidth = 300),

    dashboardSidebar(width = 300,
                     sidebarMenu(
                       HTML(paste0(
                         "<br>",
                         "<a href='https://mp.weixin.qq.com/s/OPhXbBJQC-3gQ3dr5tLZ7Q' target='_blank'><img style = 'display: block; margin-left: auto; margin-right: auto;' src='log.png' width = '186'></a>",
                         "<br>",
                         "<p style = 'text-align: center;'><small><a href='https://www.github.com//xiayh17' target='_blank'>Probes2Symbol logo designer</a></small></p>",
                         "<br>"
                       )),
                       menuItem("Home", tabName = "home", icon = icon("home")),
                       menuItem("Download Data", tabName = "download", icon = icon("thumbtack")),
                       menuItem("Previous Tables", tabName = "table", icon = icon("table")),
                       menuItem("Species Tree", tabName = "tree", icon = icon("random", lib = "glyphicon")),
                       menuItem("Species Charts", tabName = "charts", icon = icon("stats", lib = "glyphicon")),
                       menuItem("Species Choropleth Map", tabName = "choropleth", icon = icon("map marked alt")),
                       menuItem("Releases", tabName = "releases", icon = icon("tasks")),
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

        tabItem(tabName = "home",

                # home section
                includeMarkdown("www/home.md")

        ),

        tabItem(tabName = "map",

                # parks map section
                # leafletOutput("parksMap") %>% withSpinner(color = "green")

        ),

        tabItem(
          # species data section
          tabName = "table", dataTableOutput("speciesDataTable") %>% withSpinner(color = "green")

        ),

        tabItem(tabName = "tree",

                # collapsible species tree section
                # includeMarkdown("www/tree.md"),
                column(3, uiOutput("parkSelectComboTree")),
                column(3, uiOutput("categorySelectComboTree")),
                collapsibleTreeOutput('tree', height='700px') %>% withSpinner(color = "green")

        ),

        tabItem(tabName = "charts",

                # ggplot2 species charts section
                # includeMarkdown("www/charts.md"),
                fluidRow(column(3, uiOutput("categorySelectComboChart"))),
                column(6, plotOutput("ggplot2Group1") %>% withSpinner(color = "green")),
                column(6, plotOutput("ggplot2Group2") %>% withSpinner(color = "green"))

        ),

        tabItem(tabName = "choropleth",

                # choropleth species map section
                # includeMarkdown("www/choropleth.md"),
                fluidRow(
                  column(3, uiOutput("statesSelectCombo")),
                  column(3, uiOutput("categorySelectComboChoro"))
                ),
                fluidRow(
                  column(3,tableOutput('stateCategoryList') %>% withSpinner(color = "green")),
                  #column(9,leafletOutput("choroplethCategoriesPerState") %>% withSpinner(color = "green"))
                )

        )#,

        #tabItem(tabName = "releases", includeMarkdown("www/releases.md"))

      )

    ) # end dashboardBody

  )# end dashboardPage

)

server <- function(input, output, session) {

}

shinyApp(ui, server)
