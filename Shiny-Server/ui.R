library(shiny)
library(shinyjs)
library(rmarkdown)
#library(gmailr)
library(NGLVieweR)
library(dplyr)
library(shinycssloaders)
library(htmltools)
library(DT)
library(stringr)
library(shinyFeedback)
### http://shiny.russelllab.org/proteorizer/

ui <- fluidPage(
  useShinyFeedback(),
  navbarPage(
    "Proteorizer",
    
    tags$style(HTML("
              .navbar-default .navbar-nav > li > a[data-value='Submit'] {background-color: firebrick;   color:white}
              .navbar-default .navbar-nav > li > a[data-value='Introduction'] {background-color: #F0F8FF;   color:black}
              .navbar-default .navbar-nav > li > a[data-value='Explore the results'] {background-color: lightcoral;   color:black}
              .navbar-default .navbar-nav > li > a[data-value='Discover'] {background-color: lightsalmon;   color:black}
              .navbar-default .navbar-nav > li > a[data-value='Contact'] {background-color: aliceblue;   color:black}
            ")),
    
    tabPanel("Introduction",
             icon = icon("flask"),lib="font-awesome",
             includeMarkdown("instructions.Rmd")
    ),
    # Second tab with input form
    tabPanel("Submit",
             fluidRow(
               column(3,textInput("request", "Enter your request (single line):", 
                                  placeholder = "RHOA/Y34C"))),
             fluidRow(
               column(2,selectInput("option_input","Select a clustering method:",
                                    choices = c("RandomWalk","HClust"))),
               column(2,selectInput("dataset_input","Select dataset to use:",
                                    choices = c("Uniprot","Humsavar","Both")))),
             fluidRow(
               column(3,fileInput("file_input","OPTIONAL: Select a custom alignment:"))),
             fluidRow(
               column(2,actionButton("example","Example")),
               column(2,actionButton("submit", "Submit"))),
             icon = icon("envelope"),lib="font-awesome"
             
    ),
    tabPanel("Explore the results",
             icon = icon("magnifying-glass"),lib="font-awesome",
             tabsetPanel(
               tabPanel("Data Explanation",
                        includeMarkdown("explanation.Rmd"),
                        icon = icon("eye"), lib = "font-awesome",
                        id = "dataexplana"),
               tabPanel("Data Table", 
                        downloadButton("download_table"),
                        tags$style(type='text/css', "#download_table { margin-top: 25px;}"),
                        shinycssloaders::withSpinner(DTOutput("mytable1"), id="spinner1"),
                        icon = icon("keyboard"),lib="font-awesome",
                        id = "datatab"),
               tabPanel("Clusters",
                        shinycssloaders::withSpinner(imageOutput("image_1"), id="spinner2"),
                        icon = icon("globe"),lib="font-awesome",
                        id = "clust"),
               tabPanel("AnnotatedAlignment",
                        fluidRow(
                          column(2, numericInput("windo", "Windowsize:", 30, min = 1, max = 200)),
                          column(2, numericInput("topseqs", "Sequences:", 15, min = 1, max = 30)),
                          column(2, actionButton("filter", "Apply settings")),
                          column(3, downloadButton("download_alignment", "Download Alignment"))
                        ),
                        tags$style(type='text/css', "#download_alignment { width:100%; margin-top: 25px;}"),
                        tags$style(type='text/css', "#filter { width:100%; margin-top: 25px;}"),
                        shinycssloaders::withSpinner(uiOutput("image_2"), id="spinner3"),
                        icon = icon("grip-lines"),lib="font-awesome",
                        id = "align"),
               tabPanel("Structure ",
                        icon = icon("dice-d20"),lib="font-awesome",
                        id = "struct",
                        fluidRow(
                          column(10,
                                 shinycssloaders::withSpinner(NGLVieweROutput("structure")),
                                 fluidRow(
                                   column(2,actionButton("add", "Add Clusterinfo")),
                                   column(2,actionButton("reset", "Reset Clusterinfo")),
                                   column(2,selectInput("color","Background",c('black','white'))),
                                   column(2,actionButton("snapshot", "Snapshot"))
                                 ),
                                 tags$style(type='text/css', "#add { margin-top: 25px;}"),
                                 tags$style(type='text/css', "#reset { margin-top: 25px;}"),
                                 tags$style(type='text/css', "#snapshot { margin-top: 25px;}"),
                          )
                        )
               )
               
             )
    ),
    tabPanel("Discover",
             icon = icon("compass"),lib="font-awesome",
             tabsetPanel(
               tabPanel("Data Selection",
                        selectInput("dataset","Select Dataset",choices = c("",
                                                                           "Proteorizer_RW_VUS_Humsavar",
                                                                           "Proteorizer_HClust_VUS_Humsavar"
                        )),
                        selectInput("proteincase",
                                    "Waiting for dataset selection. Loading...",
                                    choices = NULL),
                        includeMarkdown("discover.Rmd"),
                        icon = icon("binoculars"), lib = "font-awesome"),
               tabPanel("Data Table", 
                        downloadButton("download_table_discover"),
                        tags$style(type='text/css', "#download_table_discover { margin-top: 25px;}"),
                        DTOutput("mytable1_discover"),
                        icon = icon("keyboard"),lib="font-awesome"),
               tabPanel("Clusters",
                        imageOutput("image_1_discover"),
                        icon = icon("globe"),lib="font-awesome"),
               
               tabPanel("AnnotatedAlignment",
                        p("Note: If the alignment is completely shown (with multiple Input positions) and appears too small, try downloading (right click) and viewing it externally, for example via InkScape."),
                        imageOutput("image_2_discover"),
                        icon = icon("grip-lines"),lib="font-awesome"),
               
               tabPanel("Structure ",
                        icon = icon("dice-d20"),lib="font-awesome",
                        fluidRow(
                          column(10,
                                 NGLVieweROutput("structure_discover"),
                                 fluidRow(
                                   column(2,actionButton("add_discover", "Add Clusterinfo")),
                                   column(2,actionButton("reset_discover", "Reset Clusterinfo")),
                                   column(2,selectInput("color_discover","Background",c('black','white'))),
                                   column(2,actionButton("snapshot_discover", "Snapshot"))
                                 ),
                                 tags$style(type='text/css', "#add_discover { margin-top: 25px;}"),
                                 tags$style(type='text/css', "#reset_discover { margin-top: 25px;}"),
                                 tags$style(type='text/css', "#snapshot_discover { margin-top: 25px;}"),
                          )
                        )
               )
               
             )
    ),
    tabPanel("Contact",
             icon = icon("address-card"),lib="font-awesome",
             includeMarkdown("contact.Rmd")
    )
  )
)
#######################################################################################################################################