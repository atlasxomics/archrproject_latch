# this version produced to use when the project has only one Sample!

library(shiny) 
library(shinyhelper) 
library(data.table) 
library(Matrix) 
library(DT) 
library(magrittr) 
library(plotly)

sc1conf = readRDS("sc1conf.rds")
sc1def  = readRDS("sc1def.rds")
sc2conf = readRDS("sc2conf.rds")
sc2def  = readRDS("sc2def.rds")


### Start server code 
shinyUI(fluidPage(  
### HTML formatting of error messages 
tags$head(tags$style(HTML(".shiny-output-error-validation {color: red; font-weight: bold;}"))), 
list(tags$style(HTML(".navbar-default .navbar-nav { font-weight: bold; font-size: 16px; background-color:#D3D3D3; }"))), 
 
   
### Page title 
# titlePanel("Data Analysis"),  


navbarPage( 
  img(src='AtlasXomicslogo.png', height = '30px', width = '200px'),
  NULL,   
 
 # navbarMenu("Gene Accessibility",### Tab1.a1: Metadata vs geneAccs on dimRed 
  tabPanel( 
    HTML("Metadata vs GeneAccs"), 
    h4("Cell information vs gene accessibility on reduced dimensions"), 
    "In this tab, users can visualise both cell information and gene ",  
    "accessibility side-by-side on low-dimensional representions.", 
    br(),br(), 
    
    
    
    
    fluidRow( 
      column( 
        6, 
        # style="border-right: 2px solid black", h4(""), 

        fluidPage(
          sidebarLayout(
            sidebarPanel( width = 6,
                           
                          selectInput("sc1a1drX","Dimension Reduction", 
                                      choices = sc1conf[-c(which(sc1conf$dimred==TRUE & grepl("_2",sc1conf$UI))),][dimred == TRUE]$UI, 
                                      selected = sc1def$dimred[1]), 
                          # selectInput("sc1a1drY", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI,
                          #             selected = sc1def$dimred[2]),
                          selectInput("sc1a1inp1", "Cell information:", 
                                      choices = sc1conf$UI, 
                                      selected = sc1def$meta1) %>%  
                            helper(type = "inline", size = "m", fade = TRUE, 
                                   title = "Cell information to colour cells by", 
                                   content = c("Select cell information to colour cells", 
                                               "- Categorical covariates have a fixed colour palette", 
                                               paste0("- Continuous covariates are coloured in a ",  
                                                      "Blue-Yellow-Red colour scheme, which can be ", 
                                                      "changed in the plot controls"))) , 
                          
                          actionButton("sc1a1tog1", "Toggle plot controls"), 
                          conditionalPanel( 
                            condition = "input.sc1a1tog1 % 2 == 1", 
                            radioButtons("sc1a1col1", "Colour (Continuous data):", 
                                         choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                                         selected = "Blue-Yellow-Red"), 
                            radioButtons("sc1a1ord1", "Plot order:", 
                                         choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                                         selected = "Original", inline = TRUE), 
                            checkboxInput("sc1a1lab1", "Show cell info labels", value = TRUE) 
                          ),br(),br(),
                          actionButton("sc1a1togL", "Toggle to subset cells"), 
                          conditionalPanel( 
                            condition = "input.sc1a1togL % 2 == 1", 
                            selectInput("sc1a1sub1", "Cell information to subset:", 
                                        choices = sc1conf[grp == TRUE]$UI, 
                                        selected = sc1def$grp1), 
                            uiOutput("sc1a1sub1.ui"), 
                            actionButton("sc1a1sub1all", "Select all groups", class = "btn btn-primary"), 
                            actionButton("sc1a1sub1non", "Deselect all groups", class = "btn btn-primary") 
                          ) , br(), br(),
                          downloadButton("sc1a1oup1.pdf", "Download PDF", style = "width:123px;font-size: 13px;"), 
                          downloadButton("sc1a1oup1.png", "Download PNG", style = "width:123px;font-size: 13px;"), br(), 
                          div(style="display:inline-block", 
                              numericInput("sc1a1oup1.h", "PDF / PNG height:", width = "123px", 
                                           min = 4, max = 20, value = 6, step = 0.5)), 
                          div(style="display:inline-block", 
                              numericInput("sc1a1oup1.w", "PDF / PNG width:", width = "123px", 
                                           min = 4, max = 20, value = 8, step = 0.5)), br(), 
                          actionButton("sc1a1tog9", "Toggle to show cell numbers / statistics"), 
                          conditionalPanel( 
                            condition = "input.sc1a1tog9 % 2 == 1", 
                            h4("Cell numbers / statistics"), 
                            radioButtons("sc1a1splt", "Split continuous cell info into:", 
                                         choices = c("Quartile", "Decile"), 
                                         selected = "Decile", inline = TRUE), 
                            dataTableOutput("sc1a1.dt") 
                          )          
                          
                          
                          
                          ),
            mainPanel = mainPanel(width = 6,uiOutput("sc1a1oup1.ui"))
            
            
          ))
        
        
      ), # End of column (6 space) 
      
      
      
      column(
        6, h4(""),
        
        fluidPage(
          sidebarLayout(
            sidebarPanel( width = 6,
                          
                          selectInput("sc1a1drX2", "Dimension Reduction", 
                                      choices = sc1conf[-c(which(sc1conf$dimred==TRUE & grepl("_2",sc1conf$UI))),][dimred == TRUE]$UI, 
                                      selected = sc1def$dimred[3]), 
                          # selectInput("sc1a1drY2", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI,
                          #             selected = sc1def$dimred[4]),
              selectInput("sc1a1inp2", "Gene name:", choices=NULL) %>%  
                helper(type = "inline", size = "m", fade = TRUE, 
                       title = "Gene accessibility to colour cells by", 
                       content = c("Select gene to colour cells by gene accessibility", 
                                   paste0("- Gene accessibility are coloured in a ", 
                                          "White-Red colour scheme which can be ", 
                                          "changed in the plot controls"))) ,
              actionButton("sc1a1tog2", "Toggle plot controls"), 
              conditionalPanel( 
                condition = "input.sc1a1tog2 % 2 == 1", 
                radioButtons("sc1a1col2", "Colour:", 
                             choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                             selected = "Blue-Yellow-Red"), 
                radioButtons("sc1a1ord2", "Plot order:", 
                             choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                             selected = "Max-1st", inline = TRUE) 
              ) , br(), br(),
              actionButton("sc1a1tog0", "Toggle graphics controls"), 
              conditionalPanel( 
                condition = "input.sc1a1tog0 % 2 == 1", 
                fluidRow( 
                  column( 
                    6, sliderInput("sc1a1siz", "Point size:", 
                                   min = 0, max = 4, value = 2.00, step = 0.25), 
                    radioButtons("sc1a1psz", "Plot size:", 
                                 choices = c("Small", "Medium", "Large"), 
                                 selected = "Small", inline = TRUE), 
                    radioButtons("sc1a1fsz", "Font size:", 
                                 choices = c("Small", "Medium", "Large"), 
                                 selected = "Small", inline = TRUE) 
                  ), 
                  column( 
                    6, radioButtons("sc1a1asp", "Aspect ratio:", 
                                    choices = c("Square", "Fixed", "Free"), 
                                    selected = "Square", inline = TRUE), 
                    checkboxInput("sc1a1txt", "Show axis text", value = FALSE) 
                  ) 
                ) 
              ) , br(), br(),
              downloadButton("sc1a1oup2.pdf", "Download PDF", style = "width:123px;font-size: 13px;"),
              downloadButton("sc1a1oup2.png", "Download PNG", style = "width:123px;font-size: 13px;"), br(),
              div(style="display:inline-block",
                  numericInput("sc1a1oup2.h", "PDF / PNG height:", width = "123px",
                               min = 4, max = 20, value = 6, step = 0.5)),
              div(style="display:inline-block",
                  numericInput("sc1a1oup2.w", "PDF / PNG width:", width = "123px",
                               min = 4, max = 20, value = 8, step = 0.5))
              
              
              
            ),
            mainPanel = mainPanel( width = 6,plotOutput("sc1a1oup2")
                                  
                                  ),
            
         
         
          ))
        
        
        
        
        
      )  # End of column (6 space)
    )    # End of fluidRow (4 space) 
    
 
 
  ),# End of tabPanel
 
 
 
 
 
 
 
 
  ### Tab1.a2: Metadata vs Metadata on dimRed 
  tabPanel( 
    HTML("Metadata vs Metadata"), 
    h4("Cell information vs cell information on dimension reduction"), 
    "In this tab, users can visualise two cell informations side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
  
    fluidRow( 
      column( 
        6, h4(""), 
     
        
        
        fluidPage(
          sidebarLayout(
            sidebarPanel(width = 6,
                         selectInput("sc1a2drX", "Dimension Reduction", 
                                     choices = sc1conf[-c(which(sc1conf$dimred==TRUE & grepl("_2",sc1conf$UI))),][dimred == TRUE]$UI, 
                                     selected = sc1def$dimred[5]), 
                         # selectInput("sc1a2drY", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI, 
                         #             selected = sc1def$dimred[6]),
                         selectInput("sc1a2inp1", "Cell information:", 
                                     choices = sc1conf$UI, 
                                     selected = sc1def$meta1) %>%  
                           helper(type = "inline", size = "m", fade = TRUE, 
                                  title = "Cell information to colour cells by", 
                                  content = c("Select cell information to colour cells", 
                                              "- Categorical covariates have a fixed colour palette", 
                                              paste0("- Continuous covariates are coloured in a ",  
                                                     "Blue-Yellow-Red colour scheme, which can be ", 
                                                     "changed in the plot controls"))),
                         actionButton("sc1a2tog1", "Toggle plot controls"), 
                         conditionalPanel( 
                           condition = "input.sc1a2tog1 % 2 == 1", 
                           radioButtons("sc1a2col1", "Colour (Continuous data):", 
                                        choices = c("White-Red", "Blue-Yellow-Red", 
                                                    "Yellow-Green-Purple"), 
                                        selected = "Blue-Yellow-Red"), 
                           radioButtons("sc1a2ord1", "Plot order:", 
                                        choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                                        selected = "Original", inline = TRUE), 
                           checkboxInput("sc1a2lab1", "Show cell info labels", value = TRUE) 
                         ) ,br(),br(),
                         
                         actionButton("sc1a2togL", "Toggle to subset cells"), 
                         conditionalPanel( 
                           condition = "input.sc1a2togL % 2 == 1", 
                           selectInput("sc1a2sub1", "Cell information to subset:", 
                                       choices = sc1conf[grp == TRUE]$UI, 
                                       selected = sc1def$grp1), 
                           uiOutput("sc1a2sub1.ui"), 
                           actionButton("sc1a2sub1all", "Select all groups", class = "btn btn-primary"), 
                           actionButton("sc1a2sub1non", "Deselect all groups", class = "btn btn-primary") 
                         ),br(),br(),
                         
                         downloadButton("sc1a2oup1.pdf", "Download PDF", style = "width:123px;font-size: 13px;"), 
                         downloadButton("sc1a2oup1.png", "Download PNG", style = "width:123px;font-size: 13px;"), br(), 
                         div(style="display:inline-block", 
                             numericInput("sc1a2oup1.h", "PDF / PNG height:", width = "123px", 
                                          min = 4, max = 20, value = 6, step = 0.5)), 
                         div(style="display:inline-block", 
                             numericInput("sc1a2oup1.w", "PDF / PNG width:", width = "123px", 
                                          min = 4, max = 20, value = 8, step = 0.5))             
                         
                         ),
            mainPanel = mainPanel(width = 6,uiOutput("sc1a2oup1.ui")
                                  
            )
          ))
        
        
      ), # End of column (6 space) 
      column( 
        6, h4(""), 
      
        
        
        fluidPage(
          sidebarLayout(
            sidebarPanel(width = 6,
                         selectInput("sc1a2drX2", "Dimension Reduction", 
                                     choices = sc1conf[-c(which(sc1conf$dimred==TRUE & grepl("_2",sc1conf$UI))),][dimred == TRUE]$UI, 
                                     selected = sc1def$dimred[3]), 
                         # selectInput("sc1a2drY2", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI, 
                         #             selected = sc1def$dimred[4]),
                         selectInput("sc1a2inp2", "Cell information:", 
                                     choices = sc1conf$UI, 
                                     selected = sc1def$meta1) %>%  
                           helper(type = "inline", size = "m", fade = TRUE, 
                                  title = "Cell information to colour cells by", 
                                  content = c("Select cell information to colour cells", 
                                              "- Categorical covariates have a fixed colour palette", 
                                              paste0("- Continuous covariates are coloured in a ",  
                                                     "Blue-Yellow-Red colour scheme, which can be ", 
                                                     "changed in the plot controls"))) ,
                         
                         actionButton("sc1a2tog2", "Toggle plot controls"), 
                         conditionalPanel( 
                           condition = "input.sc1a2tog2 % 2 == 1", 
                           radioButtons("sc1a2col2", "Colour (Continuous data):", 
                                        choices = c("White-Red", "Blue-Yellow-Red", 
                                                    "Yellow-Green-Purple"), 
                                        selected = "Blue-Yellow-Red"), 
                           radioButtons("sc1a2ord2", "Plot order:", 
                                        choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                                        selected = "Original", inline = TRUE), 
                           checkboxInput("sc1a2lab2", "Show cell info labels", value = TRUE) 
                         ),br(),br(),
                         actionButton("sc1a2tog0", "Toggle graphics controls"), 
                         conditionalPanel( 
                           condition = "input.sc1a2tog0 % 2 == 1", 
                           fluidRow( 
                             column( 
                               6, sliderInput("sc1a2siz", "Point size:", 
                                              min = 0, max = 4, value = 2.00, step = 0.25), 
                               radioButtons("sc1a2psz", "Plot size:", 
                                            choices = c("Small", "Medium", "Large"), 
                                            selected = "Small", inline = TRUE), 
                               radioButtons("sc1a2fsz", "Font size:", 
                                            choices = c("Small", "Medium", "Large"), 
                                            selected = "Small", inline = TRUE) 
                             ), 
                             column( 
                               6, radioButtons("sc1a2asp", "Aspect ratio:", 
                                               choices = c("Square", "Fixed", "Free"), 
                                               selected = "Square", inline = TRUE), 
                               checkboxInput("sc1a2txt", "Show axis text", value = FALSE) 
                             ) 
                           ) 
                         ) ,br(),br(),
                         downloadButton("sc1a2oup2.pdf", "Download PDF", style = "width:123px;font-size: 13px;"), 
                         downloadButton("sc1a2oup2.png", "Download PNG", style = "width:123px;font-size: 13px;"), br(), 
                         div(style="display:inline-block", 
                             numericInput("sc1a2oup2.h", "PDF / PNG height:", width = "123px", 
                                          min = 4, max = 20, value = 6, step = 0.5)), 
                         div(style="display:inline-block", 
                             numericInput("sc1a2oup2.w", "PDF / PNG width:", width = "123px", 
                                          min = 4, max = 20, value = 8, step = 0.5)) ),
            mainPanel = mainPanel(width = 6,uiOutput("sc1a2oup2.ui")
                                  
                                  )
          ))
        
        
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
    
    
  ),     # End of tab (2 space) 
   
  ### Tab1.a3: geneAccs vs geneAccs on dimRed 
  tabPanel( 
    HTML("GeneAccs vs GeneAccs"), 
    h4("Gene accessibility vs gene accessibility on dimension reduction"), 
    "In this tab, users can visualise two gene accessibilities side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
    
    
    
    
    fluidRow( 
      column( 
        6, 
         
        fluidPage(
          sidebarLayout(
            sidebarPanel(width = 6,
                         selectInput("sc1a3drX", "Dimension Reduction", 
                                     choices = sc1conf[-c(which(sc1conf$dimred==TRUE & grepl("_2",sc1conf$UI))),][dimred == TRUE]$UI, 
                                     selected = sc1def$dimred[5]), 
                         # selectInput("sc1a3drY", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI, 
                         #             selected = sc1def$dimred[6]),
                         selectInput("sc1a3inp1", "Gene name:", choices=NULL) %>%  
                           helper(type = "inline", size = "m", fade = TRUE, 
                                  title = "Gene accessibility to colour cells by", 
                                  content = c("Select gene to colour cells by gene accessibility", 
                                              paste0("- Gene accessibility are coloured in a ", 
                                                     "White-Red colour scheme which can be ", 
                                                     "changed in the plot controls"))) ,
                         actionButton("sc1a3tog1", "Toggle plot controls"), 
                         conditionalPanel( 
                           condition = "input.sc1a3tog1 % 2 == 1", 
                           radioButtons("sc1a3col1", "Colour:", 
                                        choices = c("White-Red", "Blue-Yellow-Red", 
                                                    "Yellow-Green-Purple"), 
                                        selected = "Blue-Yellow-Red"), 
                           radioButtons("sc1a3ord1", "Plot order:", 
                                        choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                                        selected = "Max-1st", inline = TRUE) 
                         ),br(),br(),
                         actionButton("sc1a3togL", "Toggle to subset cells"), 
                         conditionalPanel( 
                           condition = "input.sc1a3togL % 2 == 1", 
                           selectInput("sc1a3sub1", "Cell information to subset:", 
                                       choices = sc1conf[grp == TRUE]$UI, 
                                       selected = sc1def$grp1), 
                           uiOutput("sc1a3sub1.ui"), 
                           actionButton("sc1a3sub1all", "Select all groups", class = "btn btn-primary"), 
                           actionButton("sc1a3sub1non", "Deselect all groups", class = "btn btn-primary") 
                         ),br(),br(),
                         downloadButton("sc1a3oup1.pdf", "Download PDF", style = "width:123px;font-size: 13px;"), 
                         downloadButton("sc1a3oup1.png", "Download PNG", style = "width:123px;font-size: 13px;"), br(), 
                         div(style="display:inline-block", 
                             numericInput("sc1a3oup1.h", "PDF / PNG height:", width = "123px", 
                                          min = 4, max = 20, value = 6, step = 0.5)), 
                         div(style="display:inline-block", 
                             numericInput("sc1a3oup1.w", "PDF / PNG width:", width = "123px", 
                                          min = 4, max = 20, value = 8, step = 0.5))
                         
                         ),
            mainPanel = mainPanel(width = 6,uiOutput("sc1a3oup1.ui"))
          ))
        
        
      ), # End of column (6 space) 
      column( 
        6, h4(""), 

        fluidPage(
          sidebarLayout(
            sidebarPanel(width = 6,
                         
                         selectInput("sc1a3drX2", "Dimension Reduction", 
                                     choices = sc1conf[-c(which(sc1conf$dimred==TRUE & grepl("_2",sc1conf$UI))),][dimred == TRUE]$UI, 
                                     selected = sc1def$dimred[3]), 
                         # selectInput("sc1a3drY2", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI, 
                         #             selected = sc1def$dimred[4]),
                         selectInput("sc1a3inp2", "Gene name:", choices=NULL) %>%  
                           helper(type = "inline", size = "m", fade = TRUE, 
                                  title = "Gene accessibility to colour cells by", 
                                  content = c("Select gene to colour cells by gene accessibility", 
                                              paste0("- Gene accessibility are coloured in a ", 
                                                     "White-Red colour scheme which can be ", 
                                                     "changed in the plot controls"))) ,
                         actionButton("sc1a3tog2", "Toggle plot controls"), 
                         conditionalPanel( 
                           condition = "input.sc1a3tog2 % 2 == 1", 
                           radioButtons("sc1a3col2", "Colour:", 
                                        choices = c("White-Red", "Blue-Yellow-Red", 
                                                    "Yellow-Green-Purple"), 
                                        selected = "Blue-Yellow-Red"), 
                           radioButtons("sc1a3ord2", "Plot order:", 
                                        choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                                        selected = "Max-1st", inline = TRUE) 
                         ) , br(),br(),
                         actionButton("sc1a3tog0", "Toggle graphics controls"), 
                         conditionalPanel( 
                           condition = "input.sc1a3tog0 % 2 == 1", 
                           fluidRow( 
                             column( 
                               6, sliderInput("sc1a3siz", "Point size:", 
                                              min = 0, max = 4, value = 2.00, step = 0.25), 
                               radioButtons("sc1a3psz", "Plot size:", 
                                            choices = c("Small", "Medium", "Large"), 
                                            selected = "Small", inline = TRUE), 
                               radioButtons("sc1a3fsz", "Font size:", 
                                            choices = c("Small", "Medium", "Large"), 
                                            selected = "Small", inline = TRUE) 
                             ), 
                             column( 
                               6, radioButtons("sc1a3asp", "Aspect ratio:", 
                                               choices = c("Square", "Fixed", "Free"), 
                                               selected = "Square", inline = TRUE), 
                               checkboxInput("sc1a3txt", "Show axis text", value = FALSE) 
                             ) 
                           ) 
                         ) , br(),br(),
                         downloadButton("sc1a3oup2.pdf", "Download PDF", style = "width:123px;font-size: 13px;"), 
                         downloadButton("sc1a3oup2.png", "Download PNG", style = "width:123px;font-size: 13px;"), br(), 
                         div(style="display:inline-block", 
                             numericInput("sc1a3oup2.h", "PDF / PNG height:", width = "123px", 
                                          min = 4, max = 20, value = 6, step = 0.5)), 
                         div(style="display:inline-block", 
                             numericInput("sc1a3oup2.w", "PDF / PNG width:", width = "123px", 
                                          min = 4, max = 20, value = 8, step = 0.5))   
                         ),
            mainPanel = mainPanel(width = 6,uiOutput("sc1a3oup2.ui"))
          ))
        
        
      )  # End of column (6 space) 

    )    # End of fluidRow (4 space) 
    
    
    
    
    
  ),     # End of tab (2 space) 
 
 ### Tab1.b2: Gene coAccsession plot 
 ### Tab1.b2: Gene coAccsession plot ###########################################
 tabPanel(
   HTML("Gene coAccsession"),
   h4("CoAccsession of multiple genes on reduced dimensions"),
   "In this tab, users can visualise the average coAccsessibility of programs (gene sets)",
   "on low-dimensional representions.",
   br(),br(),
   fluidPage(
     sidebarLayout(
       sidebarPanel(
         
         selectInput("sc1b2drX", "Dimension Reduction", 
                     choices = sc1conf[-c(which(sc1conf$dimred==TRUE & grepl("_2",sc1conf$UI))),][dimred == TRUE]$UI,
                     selected = sc1def$dimred[1]),
         # selectInput("sc1b2drY", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI,
         #             selected = sc1def$dimred[2]),
         
         actionButton("sc1b2tog0", "Toggle graphics controls"),
         conditionalPanel(
           condition = "input.sc1b2tog0 % 2 == 1"
           
           
           , sliderInput("sc1b2siz", "Point size:",
                         min = 0, max = 4, value = 2.0, step = 0.25),
           radioButtons("sc1b2psz", "Plot size:",
                        choices = c("Small", "Medium", "Large"),
                        selected = "Small", inline = TRUE),
           radioButtons("sc1b2fsz", "Font size:",
                        choices = c("Small", "Medium", "Large"),
                        selected = "Small", inline = TRUE),
           
           
           radioButtons("sc1b2asp", "Aspect ratio:",
                        choices = c("Square", "Fixed", "Free"),
                        selected = "Square", inline = TRUE),
           checkboxInput("sc1b2txt", "Show axis text", value = FALSE)
           
           
         ),br(), br(),
         
         downloadButton("sc1b2oup1.pdf", "Download PDF"),
         downloadButton("sc1b2oup1.png", "Download PNG"), br(),
         div(style="display:inline-block",
             numericInput("sc1b2oup1.h", "PDF / PNG height:", width = "138px",
                          min = 4, max = 20, value = 8, step = 0.5)),
         div(style="display:inline-block",
             numericInput("sc1b2oup1.w", "PDF / PNG width:", width = "138px",
                          min = 4, max = 20, value = 10, step = 0.5))
         
         
       ), mainPanel = mainPanel(
         
         
         
         h4(htmlOutput("sc1b2oupTxt")),
         uiOutput("sc1b2oup1.ui"),
         
       )# End of column (6 space)
       
 
     ))
 
 ),     # End of tab (2 space)
 
 ### Tab1.c1: violinplot / boxplot 
 tabPanel( 
    HTML("Violinplot / Boxplot"),  
   h4("Cell information / gene accessibility violin plot / box plot"), 
   "In this tab, users can visualise the gene accessibility or continuous cell information ",  
   "(e.g. Number of Fragment / TSS) across groups of cells (e.g. libary / clusters).", 
   br(),br(), 
   fluidPage(
     sidebarLayout(
       sidebarPanel(
       selectInput("sc1c1inp1", "Cell information (X-axis):", 
                   choices = sc1conf[grp == TRUE]$UI, 
                   selected = sc1def$grp1) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell information to group cells by",  
                content = c("Select categorical cell information to group cells by",  
                            "- Single cells are grouped by this categorical covariate",  
                            "- Plotted as the X-axis of the violin plot / box plot")),  
       selectInput("sc1c1inp2", "Cell Info / Gene name (Y-axis):", choices=NULL) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell Info / Gene to plot", 
                content = c("Select cell info / gene to plot on Y-axis", 
                            "- Can be continuous cell information (e.g. nFragments / TSS)", 
                            "- Can also be gene accessibility")), 
       radioButtons("sc1c1typ", "Plot type:", 
                    choices = c("violin", "boxplot"), 
                    selected = "violin", inline = TRUE), 
       checkboxInput("sc1c1pts", "Show data points", value = FALSE), 
       actionButton("sc1c1togL", "Toggle to subset cells"), 
       conditionalPanel( 
         condition = "input.sc1c1togL % 2 == 1", 
         selectInput("sc1c1sub1", "Cell information to subset:", 
                     choices = sc1conf[grp == TRUE]$UI, 
                     selected = sc1def$grp1), 
         uiOutput("sc1c1sub1.ui"), 
         actionButton("sc1c1sub1all", "Select all groups", class = "btn btn-primary"), 
         actionButton("sc1c1sub1non", "Deselect all groups", class = "btn btn-primary") 
       ), br(), br(), 
       actionButton("sc1c1tog", "Toggle graphics controls"), 
       conditionalPanel( 
         condition = "input.sc1c1tog % 2 == 1", 
         sliderInput("sc1c1siz", "Data point size:",  
                     min = 0, max = 4, value = 0.25, step = 0.25),  
         radioButtons("sc1c1psz", "Plot size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Small", inline = TRUE), 
         radioButtons("sc1c1fsz", "Font size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Small", inline = TRUE)) ,br(), br(), 
       downloadButton("sc1c1oup.pdf", "Download PDF"),  
       downloadButton("sc1c1oup.png", "Download PNG"), br(), 
       div(style="display:inline-block", 
           numericInput("sc1c1oup.h", "PDF / PNG height:", width = "138px", 
                        min = 4, max = 20, value = 8, step = 0.5)), 
       div(style="display:inline-block", 
           numericInput("sc1c1oup.w", "PDF / PNG width:", width = "138px", 
                        min = 4, max = 20, value = 10, step = 0.5)) 
       
       
     ),mainPanel = mainPanel(
       plotOutput('sc1c1oup')
     ) # End of column (6 space) 
      

     )  # End of column (6 space) 
   )    # End of fluidRow (4 space) 
 ),     # End of tab (2 space) 
 
### Tab1.c2: Proportion plot 
tabPanel( 
  HTML("Proportion plot"), 
  h4("Proportion / cell numbers across different cell information"), 
  "In this tab, users can visualise the composition of single cells based on one discrete ", 
  "cell information across another discrete cell information. ",  
  "Usage examples include the library or cellcycle composition across clusters.", 
  br(),br(), 
  fluidPage(
    sidebarLayout(
      sidebarPanel(
      selectInput("sc1c2inp1", "Cell information to plot (X-axis):", 
                  choices = sc1conf[grp == TRUE]$UI, 
                  selected = sc1def$grp2) %>%  
        helper(type = "inline", size = "m", fade = TRUE, 
               title = "Cell information to plot cells by",  
               content = c("Select categorical cell information to plot cells by", 
                           "- Plotted as the X-axis of the proportion plot")), 
      selectInput("sc1c2inp2", "Cell information to group / colour by:", 
                  choices = sc1conf[grp == TRUE]$UI, 
                  selected = sc1def$grp1) %>%  
        helper(type = "inline", size = "m", fade = TRUE, 
               title = "Cell information to group / colour cells by", 
               content = c("Select categorical cell information to group / colour cells by", 
                           "- Proportion / cell numbers are shown in different colours")), 
      radioButtons("sc1c2typ", "Plot value:", 
                   choices = c("Proportion", "CellNumbers"), 
                   selected = "Proportion", inline = TRUE), 
      checkboxInput("sc1c2flp", "Flip X/Y", value = FALSE), 
      actionButton("sc1c2togL", "Toggle to subset cells"), 
      conditionalPanel( 
        condition = "input.sc1c2togL % 2 == 1", 
        selectInput("sc1c2sub1", "Cell information to subset:", 
                    choices = sc1conf[grp == TRUE]$UI, 
                    selected = sc1def$grp1), 
        uiOutput("sc1c2sub1.ui"), 
        actionButton("sc1c2sub1all", "Select all groups", class = "btn btn-primary"), 
        actionButton("sc1c2sub1non", "Deselect all groups", class = "btn btn-primary") 
      ), br(), br(), 
      actionButton("sc1c2tog", "Toggle graphics controls"), 
      conditionalPanel( 
        condition = "input.sc1c2tog % 2 == 1", 
        radioButtons("sc1c2psz", "Plot size:", 
                     choices = c("Small", "Medium", "Large"), 
                     selected = "Small", inline = TRUE), 
        radioButtons("sc1c2fsz", "Font size:", 
                     choices = c("Small", "Medium", "Large"), 
                     selected = "Small", inline = TRUE)),br(), br(),   
      downloadButton("sc1c2oup.pdf", "Download PDF"),  
      downloadButton("sc1c2oup.png", "Download PNG"), br(), 
      div(style="display:inline-block", 
          numericInput("sc1c2oup.h", "PDF / PNG height:", width = "138px", 
                       min = 4, max = 20, value = 8, step = 0.5)), 
      div(style="display:inline-block", 
          numericInput("sc1c2oup.w", "PDF / PNG width:", width = "138px", 
                       min = 4, max = 20, value = 10, step = 0.5))  
    ), mainPanel = mainPanel(
      plotOutput('sc1c2oup')
    )# End of column (6 space) 
    )  # End of column (6 space) 
  )    # End of fluidRow (4 space) 
),     # End of tab (2 space) 
 
##### START OF NEW CODE TO ADD #####

##### END OF NEW CODE TO ADD #####


  ### Tab1.d1: Multiple gene Accs 
  tabPanel( 
    HTML(" Gene Access Heatmap"),
    h4("Gene accessibility heatmap"),
    "In this tab, users can visualise the gene accessibility patterns of ",
    "multiple genes grouped by categorical cell information (e.g. library / cluster).", br(),
    "The normalised accessibility are averaged, log-transformed and then plotted.",
    br(),br(),
    fluidPage(
      sidebarLayout(
        sidebarPanel(

        selectInput("sc1d1grp", "Group by:",
                    choices = 'Clusters',
                    selected = 'Clusters') %>%
          helper(type = "inline", size = "m", fade = TRUE,
                 title = "Cell information to group cells by",
                 content = c("Select categorical cell information to group cells by",
                             "- Single cells are grouped by this categorical covariate",
                             "- Plotted as the X-aqxis of the heatmap")),

    #     
        # textAreaInput("sc1d1inp", HTML("List of gene names <br />
        #                                   (Max 50 genes, separated <br />
        #                                    by , or ; or newline):"),
        #               height = "200px",
        #               value = paste0(sc1def$genes, collapse = ", ")
        #             ) %>%
        # helper(type = "inline", size = "m", fade = TRUE,
        #        title = "List of genes to plot on bubbleplot / heatmap",
        #        content = c("Input genes to plot",
        #                    "- Maximum 50 genes (due to ploting space limitations)",
        #                    "- Genes should be separated by comma, semicolon or newline")),
    # 
        radioButtons("sc1d1plt", "Plot type:",
                     choices = c(
                       # "Bubbleplot", 
                                 "Heatmap"),
                     selected = "Heatmap", inline = TRUE),
 #       checkboxInput("sc1d1scl", "Scale gene accessibility", value = TRUE),
        checkboxInput("sc1d1row", "Cluster rows (genes)", value = FALSE),
        checkboxInput("sc1d1col", "Cluster columns (samples)", value = FALSE),
        br(),
        actionButton("sc1d1togL", "Toggle to subset cells"),
        conditionalPanel(
          condition = "input.sc1d1togL % 2 == 1",
          selectInput("sc1d1sub1", "Cell information to subset:",
                      choices = sc1conf[grp == TRUE]$UI[1:4],
                      selected = sc1def$grp1),
          uiOutput("sc1d1sub1.ui"),
          actionButton("sc1d1sub1all", "Select all groups", class = "btn btn-primary"),
          actionButton("sc1d1sub1non", "Deselect all groups", class = "btn btn-primary")
        ), br(), br(),
        actionButton("sc1d1tog", "Toggle graphics controls"),
        conditionalPanel(
          condition = "input.sc1d1tog % 2 == 1",
          radioButtons("sc1d1cols", "Colour scheme:",
                       choices = c("White-Red", "Blue-Yellow-Red",
                                   "Yellow-Green-Purple","comet","blueYellow"),
                       selected = "blueYellow"),
          radioButtons("sc1d1psz", "Plot size:",
                       choices = c("Small", "Medium", "Large"),
                       selected = "Small", inline = TRUE),
          radioButtons("sc1d1fsz", "Font size:",
                       choices = c("Small", "Medium", "Large"),
                       selected = "Small", inline = TRUE)), br(), br(),
 downloadButton("sc1d1oup.pdf", "Download PDF"),
 downloadButton("sc1d1oup.png", "Download PNG"), br(),
 div(style="display:inline-block",
     numericInput("sc1d1oup.h", "PDF / PNG height:", width = "138px",
                  min = 4, max = 20, value = 10, step = 0.5)),
 div(style="display:inline-block",
     numericInput("sc1d1oup.w", "PDF / PNG width:", width = "138px",
                  min = 4, max = 20, value = 10, step = 0.5))
 
      ), mainPanel = mainPanel(
        h4(htmlOutput("sc1d1oupTxt")),
        plotOutput('sc1d1oup')
      )# End of column (6 space)


      )  # End of column (6 space)
    )  # End of fluidRow (4 space) 
  ) ,     # End of tab (2 space) 


### Tab1.x2: New tab for tracks
# tabPanel(
#   HTML("Genome Tracks"),
#   h4("Genome Tracks"),
#   "In this tab, users can visualise the genome browser for gene accessibility",
#   br(),br(),
#   
#     
#   fluidRow(
#    
#     
#     column(width = 3, offset = 0, style="background-color: #D3D3D3",
#       selectInput("sc1trackinp", "Gene Name", choices=NULL),
#       
#       br(),
#       
#       selectInput("sc1trackgrp", "Group Name", choices=NULL),
#       
#       br(),
#       
#       
#       sliderInput("range_1", "Distance From Center (kb):", min = -50, max = 50, value = dist_range),
#       # splitLayout(cellWidths = c("50%","50%"),
#       #             numericInput("range_min_1", "Distance (-kb):", min = dist_range[1], max = dist_range[2], value = dist_range[1]),
#       #             numericInput("range_max_1", "Distance (+kb):", min = dist_range[1], max = dist_range[2], value = dist_range[2])
#       # ),
#       br(),
#       downloadButton("sc1trackoup.pdf", "Download PDF"),
#       downloadButton("sc1trackoup.png", "Download PNG"), br(),
#       div(style="display:inline-block",
#           numericInput("sc1trackoup.h", "PDF / PNG height:", width = "138px",
#                        min = 4, max = 20, value = 8, step = 0.5)),
#       div(style="display:inline-block",
#           numericInput("sc1trackoup.w", "PDF / PNG width:", width = "138px",
#                        min = 4, max = 20, value = 10, step = 0.5))),
#     
#     
#     
#     column(width = 9, plotOutput("sc1trackoup"))     # End of column (6 space),
#     
#       
#       
#       
#   
#     # ), # End of column (6 space)
#   #    column(9,
#   #                    downloadButton("sc1trackoup.pdf", "Download PDF"),
#   #                    downloadButton("sc1trackoup.png", "Download PNG"), br(),
#   #                    div(style="display:inline-block",
#   #                        numericInput("sc1trackoup.h", "PDF / PNG height:", width = "138px",
#   #                                     min = 4, max = 20, value = 8, step = 0.5)),
#   #                    div(style="display:inline-block",
#   #                        numericInput("sc1trackoup.w", "PDF / PNG width:", width = "138px",
#   #                                     min = 4, max = 20, value = 10, step = 0.5)))
#   # 
#   )    # End of fluidRow (4 space)
#   
# ),     # End of tab (2 space)

##### END OF NEW CODE TO ADD #####


#### ### ### ### ### ### ### ### ### ### ### ### ### ("Peak/Motifs")### ### ### ### ### ###

### Tab1.a2: Metadata vs Metadata on dimRed 

### Tab1.a3: geneAccs vs geneAccs on dimRed 
tabPanel(
  HTML("MotifEnr vs MotifEnr"),
  h4("Motif enrichment vs motif enrichment on dimension reduction"),
  "In this tab, users can visualise two motif enrichment side-by-side ",
  "on low-dimensional representions.",
  br(),br(),
  
  
  fluidRow(
    column(
      6, 
   
      
      fluidPage(
        sidebarLayout(
          
          sidebarPanel(width = 6,
                       selectInput("sc2a3drX", "Dimension Reduction", 
                                   choices = sc2conf[-c(which(sc2conf$dimred==TRUE & grepl("_2",sc2conf$UI))),][dimred == TRUE]$UI,
                                   selected = sc2def$dimred[3]),
                       # selectInput("sc2a3drY", "Y-axis:", choices = sc2conf[dimred == TRUE]$UI,
                       #             selected = sc2def$dimred[4]),
                       selectInput("sc2a3inp1", "Motif name:", choices=NULL) %>%
                         helper(type = "inline", size = "m", fade = TRUE,
                                title = "Motif enrichment to colour cells by",
                                content = c("Select motif to colour cells by motif enrichment",
                                            paste0("- Motif enrichment are coloured in a ",
                                                   "White-Red colour scheme which can be ",
                                                   "changed in the plot controls"))),
                       actionButton("sc2a3tog1", "Toggle plot controls"),
                       conditionalPanel(
                         condition = "input.sc2a3tog1 % 2 == 1",
                         radioButtons("sc2a3col1", "Colour:",
                                      choices = c("White-Red", "Blue-Yellow-Red",
                                                  "Yellow-Green-Purple"),
                                      selected = "White-Red"),
                         radioButtons("sc2a3ord1", "Plot order:",
                                      choices = c("Max-1st", "Min-1st", "Original", "Random"),
                                      selected = "Max-1st", inline = TRUE)
                       ),br(),br(),
                       
                       actionButton("sc2a3togL", "Toggle to subset cells"),
                       conditionalPanel(
                         condition = "input.sc2a3togL % 2 == 1",
                         selectInput("sc2a3sub1", "Cell information to subset:",
                                     choices = sc2conf[grp == TRUE]$UI,
                                     selected = sc2def$grp1),
                         uiOutput("sc2a3sub1.ui"),
                         actionButton("sc2a3sub1all", "Select all groups", class = "btn btn-primary"),
                         actionButton("sc2a3sub1non", "Deselect all groups", class = "btn btn-primary")
                       ),br(),br(),
                       downloadButton("sc2a3oup1.pdf", "Download PDF", style = "width:123px;font-size: 13px;"),
                       downloadButton("sc2a3oup1.png", "Download PNG", style = "width:123px;font-size: 13px;"), br(),
                       div(style="display:inline-block",
                           numericInput("sc2a3oup1.h", "PDF / PNG height:", width = "123px",
                                        min = 4, max = 20, value = 6, step = 0.5)),
                       div(style="display:inline-block",
                           numericInput("sc2a3oup1.w", "PDF / PNG width:", width = "123px",
                                        min = 4, max = 20, value = 8, step = 0.5))            
                       
                       ),
          mainPanel = mainPanel(width = 6,uiOutput("sc2a3oup1.ui"))
        ))
      
      
      
  ), # End of column (6 space)
  column(
      6, 
      
      
      fluidPage(
        sidebarLayout(
    
          sidebarPanel(width = 6,
                       selectInput("sc2a3drX2", "Dimension Reduction", 
                                   choices = sc1conf[-c(which(sc1conf$dimred==TRUE & grepl("_2",sc1conf$UI))),][dimred == TRUE]$UI,
                                   selected = sc1def$dimred[5]),
                       # selectInput("sc2a3drY2", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI,
                       #             selected = sc1def$dimred[6]),
                       
                       selectInput("sc2a3inp2", "Motif name:", choices=NULL) %>%
                         helper(type = "inline", size = "m", fade = TRUE,
                                title = "Motif enrichment to colour cells by",
                                content = c("Select motif to colour cells by motif enrichment",
                                            paste0("- Motif enrichment are coloured in a ",
                                                   "White-Red colour scheme which can be ",
                                                   "changed in the plot controls"))),
                       actionButton("sc2a3tog2", "Toggle plot controls"),
                       conditionalPanel(
                         condition = "input.sc2a3tog2 % 2 == 1",
                         radioButtons("sc2a3col2", "Colour:",
                                      choices = c("White-Red", "Blue-Yellow-Red",
                                                  "Yellow-Green-Purple"),
                                      selected = "White-Red"),
                         radioButtons("sc2a3ord2", "Plot order:",
                                      choices = c("Max-1st", "Min-1st", "Original", "Random"),
                                      selected = "Max-1st", inline = TRUE)
                       ),br(),br(),
                       
                       actionButton("sc2a3tog0", "Toggle graphics controls"),
                       conditionalPanel(
                         condition = "input.sc2a3tog0 % 2 == 1",
                         fluidRow(
                           column(
                             6, sliderInput("sc2a3siz", "Point size:",
                                            min = 0, max = 4, value = 2.00, step = 0.25),
                             radioButtons("sc2a3psz", "Plot size:",
                                          choices = c("Small", "Medium", "Large"),
                                          selected = "Small", inline = TRUE),
                             radioButtons("sc2a3fsz", "Font size:",
                                          choices = c("Small", "Medium", "Large"),
                                          selected = "Small", inline = TRUE)
                           ),
                           column(
                             6, radioButtons("sc2a3asp", "Aspect ratio:",
                                             choices = c("Square", "Fixed", "Free"),
                                             selected = "Square", inline = TRUE),
                             checkboxInput("sc2a3txt", "Show axis text", value = FALSE)
                           )
                         )
                       ),br(),br(),
                       downloadButton("sc2a3oup2.pdf", "Download PDF", style = "width:123px;font-size: 13px;"),
                       downloadButton("sc2a3oup2.png", "Download PNG", style = "width:123px;font-size: 13px;"), br(),
                       div(style="display:inline-block",
                           numericInput("sc2a3oup2.h", "PDF / PNG height:", width = "123px",
                                        min = 4, max = 20, value = 6, step = 0.5)),
                       div(style="display:inline-block",
                           numericInput("sc2a3oup2.w", "PDF / PNG width:", width = "123px",
                                        min = 4, max = 20, value = 8, step = 0.5))  
                       
                       ),
          mainPanel = mainPanel(width = 6,
                                uiOutput("sc2a3oup2.ui")
                                )
        ))
      
    )  # End of column (6 space)
  )    # End of fluidRow (4 space)


),     # End of tab (2 space)

### Tab1.b2: Gene coAccsession plot 

### Tab1.c1: violinplot / boxplot 

### Tab1.c2: Proportion plot 

### Tab1.d1: Multiple gene Accs 
tabPanel(
  HTML("Motifs Heatmap"),
  h4("Motif enrichment heatmap"),
  "In this tab, users can visualise the motif enrichments of ",
  "multiple motifs grouped by categorical cell information (e.g. Sample / Cluster).", br(),
  # "The normalised enrichment are averaged, log-transformed and then plotted.",
  br(),br(),
  fluidPage(
    sidebarLayout(
      sidebarPanel(
        
        selectInput("sc2d1grp", "Group by:",
                    choices = 'Clusters',
                    selected = 'Clusters') %>%
          helper(type = "inline", size = "m", fade = TRUE,
                 title = "Cell information to group cells by",
                 content = c("Select categorical cell information to group cells by",
                             "- Single cells are grouped by this categorical covariate",
                             "- Plotted as the X-aqxis of the heatmap")),
        
        radioButtons("sc2d1plt", "Plot type:",
                     choices = c(
                       # "Bubbleplot",
                       "Heatmap"),
                     selected = "Heatmap", inline = TRUE),
        #       checkboxInput("sc2d1scl", "Scale gene accessibility", value = TRUE),
        checkboxInput("sc2d1row", "Cluster rows (genes)", value = FALSE),
        checkboxInput("sc2d1col", "Cluster columns (samples)", value = TRUE),
        br(),
        actionButton("sc2d1togL", "Toggle to subset cells"),
        conditionalPanel(
          condition = "input.sc2d1togL % 2 == 1",
          selectInput("sc2d1sub1", "Cell information to subset:",
                      choices = sc2conf[grp == TRUE]$UI[1:4],
                      selected = sc2def$grp1),
          uiOutput("sc2d1sub1.ui"),
          actionButton("sc2d1sub1all", "Select all groups", class = "btn btn-primary"),
          actionButton("sc2d1sub1non", "Deselect all groups", class = "btn btn-primary")
        ), br(), br(),
        actionButton("sc2d1tog", "Toggle graphics controls"),
        conditionalPanel(
          condition = "input.sc2d1tog % 2 == 1",
          radioButtons("sc2d1cols", "Colour scheme:",
                       choices = c("White-Red", "Blue-Yellow-Red",
                                   "Yellow-Green-Purple","comet","blueYellow"),
                       selected = "comet"),
          radioButtons("sc2d1psz", "Plot size:",
                       choices = c("Small", "Medium", "Large"),
                       selected = "Small", inline = TRUE),
          radioButtons("sc2d1fsz", "Font size:",
                       choices = c("Small", "Medium", "Large"),
                       selected = "Small", inline = TRUE)),br(),br(),
        downloadButton("sc2d1oup.pdf", "Download PDF"),
        downloadButton("sc2d1oup.png", "Download PNG"), br(),
        div(style="display:inline-block",
            numericInput("sc2d1oup.h", "PDF / PNG height:", width = "138px",
                         min = 4, max = 20, value = 10, step = 0.5)),
        div(style="display:inline-block",
            numericInput("sc2d1oup.w", "PDF / PNG width:", width = "138px",
                         min = 4, max = 20, value = 10, step = 0.5))
      ), mainPanel = mainPanel(
        h4(htmlOutput("sc2d1oupTxt")),
        plotOutput('sc2d1oup')
      )# End of column (6 space)
      
      
    )  # End of column (6 space)
  )  # End of fluidRow (4 space)
) ,     # End of tab (2 space)




### Tab1.x2: New tab for logo
tabPanel(
  HTML("Motifs logo Plot"),
  h4("Motifs logo Plot"),
  "In this tab, users can visualise the Motifs logo plot for marker Motif differential accessibility",
  br(),br(),
  fluidPage(
    sidebarLayout(
      sidebarPanel(
      selectInput("sc2logoinp", "Motif", choices=NULL),
      br(),
      # fluidRow( column(3,
                       downloadButton("sc2logooup.pdf", "Download PDF"),
                       downloadButton("sc2logooup.png", "Download PNG"), br(),
                       div(style="display:inline-block",
                           numericInput("sc2logooup.h", "PDF / PNG height:", width = "138px",
                                        min = 4, max = 20, value = 8, step = 0.5)),
                       div(style="display:inline-block",
                           numericInput("sc2logooup.w", "PDF / PNG width:", width = "138px",
                                        min = 4, max = 20, value = 10, step = 0.5)),
      
    ), mainPanel = mainPanel(
      plotOutput('sc2logooup')
    ) # End of column (6 space)
    # fluidRow(column(6, plotOutput("sc2logooup"))    
    )  # End of column (6 space),

  )    # End of fluidRow (4 space)
),     # End of tab (2 space)

##### END OF NEW CODE TO ADD #####
tabPanel(
  HTML("Genome Tracks"),
  h4("Genome Tracks"),
  "In this tab, users can visualise the genome browser for gene accessibility",
  br(),br(),
fluidPage(
  sidebarLayout(
    sidebarPanel(
      
      selectInput("sc1trackinp", "Gene Name", choices=NULL),
      
      br(),
      
      selectInput("sc1trackgrp", "Group Name", choices=NULL),
      
      br(),
      
      
      sliderInput("range_1", "Distance From Center (kb):", min = -250, max = 250, value = c(-50,50)),
      
      # actionButton("up","",icon("arrow-up")),
            splitLayout(cellWidths = c("50%","50%"),
                        numericInput("range_min_1", "Distance (-kb):", min = -250, max = 250, value = -50),
                        numericInput("range_max_1", "Distance (+kb):", min = -250, max = 250, value = 50)),
      downloadButton("sc1trackoup.pdf", "Download PDF"),
      downloadButton("sc1trackoup.png", "Download PNG"), br(),
      div(style="display:inline-block",
          numericInput("sc1trackoup.h", "PDF / PNG height:", width = "138px",
                       min = 4, max = 20, value = 8, step = 0.5)),
      div(style="display:inline-block",
          numericInput("sc1trackoup.w", "PDF / PNG width:", width = "138px",
                       min = 4, max = 20, value = 10, step = 0.5))
          
      
    ), mainPanel = mainPanel(plotOutput('sc1trackoup'))
  )
)
),



###########################################################
   
br(), 
# p(strong("Reference: "),"Pieper Lab Data Analysis ",style = "font-size: 125%;"), 
# p(em("This webpage was made using "), a("ShinyCell", 
#   href = "https://github.com/SGDDNB/ShinyCell",target="_blank")), 


br(),br(),br(),br(),br() 
))) 
 
 
 


 
