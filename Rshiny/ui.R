suppressPackageStartupMessages(library(shinydashboard))
suppressPackageStartupMessages(library(shinyWidgets))
#https://rstudio.github.io/shinydashboard/get_started.html
# Put them together into a dashboardPage
dashboardPage(
  dashboardHeader(
    title = "Nanome - Nanopore methylation database",
    titleWidth = 1200
  ),
  #Sidebar content
  dashboardSidebar(
    sidebarMenu(# dashboard should display tutorial or instructions of how to use the app
      # icon list: https://rstudio.github.io/shinydashboard/appearance.html#icons
      menuItem("About", tabName = "about", icon = icon("dashboard")),
      menuItem("Tool comparison", tabName = "table", icon = icon("table")),
      menuItem("Correlation plotting", tabName = "plot", icon = icon("chart-line"))
    )
  ),
  
  # Body content
  dashboardBody(
    tabItems(
      tabItem(
        ## 1st tab content
        tabName = "about", 
        fluidRow(
          column(12,
                 titlePanel("About"),
                 htmlOutput('about'))
        )
      ),
      
      ## 2nd tab content
      tabItem(
        tabName = "table", 
        h2("Methylation Calling Results at site level"),
        fluidRow(
          column(width = 12,
                 "*All fields are required to display the result."),
          column(width = 4,
                 pickerInput(inputId = "Dataset", label = "Dataset", 
                             choices = levels(factor(dataset_order)),
                             options = list(`actions-box` = TRUE), 
                             multiple = TRUE)
          ),
          column(width = 4,
                 pickerInput(inputId = "Chrom", label = "Chromosome", 
                             choices = levels(factor(chr_order)),
                             options = list(`actions-box` = TRUE), 
                             multiple = TRUE)
          ),
          column(width = 4,
                 pickerInput(inputId = "Strand", label = "Strand", 
                             choices = levels(factor(strand_order)),
                             options = list(`actions-box` = TRUE), 
                             multiple = TRUE)
          ),
          column(width = 2,
                 pickerInput(inputId = 'Singleton',label = "Singleton/Non-singleton",
                             choices = levels(factor(singleton_order)),
                             options = list(`actions-box` = TRUE), 
                             multiple = TRUE)
          ),
          column(width = 2,
                 pickerInput(inputId = 'Genomic_location',label = "Genomic location",
#                             choices= c(sort(unique(as.character(df$Genomic_location)))), 
                             choices = levels(factor(df$Genomic_location)),
                             options = list(`actions-box` = TRUE), 
                             multiple = TRUE)
          ),
          column(width = 2,
                 pickerInput(inputId = 'CpG_location',label = "CpG location",
                             choices= c(sort(unique(as.character(df$CpG_location)))), 
                             options = list(`actions-box` = TRUE), 
                             multiple = TRUE)
          ),
          column(width = 2,
                 pickerInput(inputId = 'CG_density',label = "CG density",
                             choices= levels(factor(CGdensity_order)), 
                             options = list(`actions-box` = TRUE), 
                             multiple = TRUE)
          ),
          column(width = 2,
                 pickerInput(inputId = 'Repetitive_regions',label = "Repetitive regions",
                             choices= c(sort(unique(as.character(df$Repetitive_regions)))), 
                             options = list(`actions-box` = TRUE), 
                             multiple = TRUE)
          ),
          column(width = 12,
                 div(style = "font-size:medium; font-weight:bold",'Methylation calling results'), 
                 DT::dataTableOutput(outputId = "table1")
          ),
          column(width = 12,
                 "All start coordinates shown here are 1-based, not 0-based.")
        )
      ),
      #3rd tab content
      tabItem(
        tabName = "plot", 
        h2("Methylation correlation plotting for each methylation-calling tool compared to BS-seq"),
        fluidRow(
          tabBox(width = 7,
                 tabPanel('NA19240',
                          img(src = "NA19240-baseFormat1.jpg", height = 520, width = 520)),
                 tabPanel('NA12878',
                          img(src = "NA12878-baseFormat1.jpg", height = 520, width = 520)),
                 tabPanel('APL',                                         
                          img(src = "APL-baseFormat1.jpg", height = 520, width = 520)),
                 tabPanel('K562',                                        
                          img(src = "K562-baseFormat1.jpg", height = 520, width = 520))
          )
        )
      )
    )
  ))