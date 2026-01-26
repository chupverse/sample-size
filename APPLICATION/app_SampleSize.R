library(shiny)
library(shinydashboard)
library(epiR)
library(shinythemes)
library(shinyjs)
library(drcarlate)
library(rpact)
library(shinycssloaders)
library(bsplus)


se <- function(width, alpha) # The standard error associated with the 1-alpha confidence interval
{
  fun <- function(x) { exp( qnorm(1-alpha/2, mean=0, sd=1) * x ) - exp(-1* qnorm(1-alpha/2, mean=0, sd=1) * x ) - width } 
  return(uniroot(fun, lower = 0.001, upper = 100)$root)
}

size.calib <- function(p0, width, alpha) # the minimum sample size to achieve this precision
{   
  (1-p0) / ((p0 * se(width=width, alpha=alpha)**2 ))
}




sample_survival <- function(p1, p0, power, r, alpha) {
  HR <- log(p1)/log(p0)
  phi <- 1/r
  m1 <- (1/phi)* (((1+phi*HR)/(1-HR))^2)* ((qnorm(p=1-(alpha/2),mean=0,sd=1)+qnorm(p=power,mean=0,sd=1))^2) / (1-p1 + phi*(1-p0))
  N.gp1 <- ceiling(m1)
  N.gp0 <- phi*N.gp1
  Ntotal <- N.gp1 + N.gp0
  res <- list(Ntotal, N.gp1, N.gp0)
  names(res) <- c("n.total","n.gp1","n.gp0")
  return(res)
}


SW <- function(ni, center, sequence, icc){
  aa <- -2*center*(sequence - 1/sequence)*icc*(1+sequence/2)
  bb <- 3*ni*(1-icc)*icc*(1+sequence) - 2*center*(sequence -1/sequence)*(1-icc)
  cc <- 3*ni*(1-icc)*(1-icc)
  m1 <- (-bb + sqrt(bb^2 - 4*aa*cc)) / (2*aa)
  m2 <- (-bb - sqrt(bb^2 - 4*aa*cc)) / (2*aa)
  m_sol <- max(m1,m2)
  Npat_center <- m_sol*(sequence+1)
  N_tot_SW <- Npat_center*center
  2*ceiling(N_tot_SW /2)
}

CRT <- function(ni,center,icc) {
  Npat_center <- (ni*(1-icc))/(center-(ni*icc))
  N_tot_CRT <- Npat_center*center
  2*ceiling(N_tot_CRT /2)
}


################################################################################################################################
################################################### UI #########################################################################
################################################################################################################################

ui <- dashboardPage(
  dashboardHeader(title = "SAMPLE SIZE CALCULATOR",
                  titleWidth = "100%"),
  dashboardSidebar(
    title = p(a(tags$img(src='LOGO.jpeg', width=230,align = "center"),
                target="_blank", href="https://www.chu-poitiers.fr", class="hidden-xs"),
              style="padding-left:0px !important"),
    sidebarMenu(
      id = "sidebarID",
      menuItem("DESCRIBE A POPULATION", tabName = "dp", id = "DPid", expandedName = "DP",
               radioButtons(inputId = "typeDP",
                            label = "Endpoint",
                            choices = c("Mean","Proportion"))),
      hidden(menuItem("hiddenChartsDP", tabName = "hiddenChartsDP")),
      menuItem("COMPARING MEANS",tabName = "comparingmean",id = "CMid",expandedName = "COMPARINGMEAN",
               radioButtons(inputId = "hypMean",
                            label = "Hypothesis",
                            choices = c("Superiority", "Non-inferiority")), 
               radioButtons(inputId = "meanDesign",
                            label = "Design",
                            choices = c("No intermediate analysis with individual randomization",
                                        "No intermediate analysis with parallel cluster randomization",
                                        "No intermediate analysis with stepped wedge randomization",
                                        "Intermediate analysis with individual randomization"))),
      hidden(menuItem("hiddenChartsMean", tabName = "hiddenChartsMean")),
      menuItem("COMPARING PROPORTION",tabName = "comparingprop",id = "CPid",expandedName = "COMPARINGPROP",
               radioButtons(inputId = "hypProp",
                            label = "Hypothesis",
                            choices = c("Superiority", "Non-inferiority")), 
               radioButtons(inputId = "propDesign",
                            label = "Design",
                            choices = c("No intermediate analysis with individual randomization", 
                                        "No intermediate analysis with parallel cluster randomization",
                                        "No intermediate analysis with stepped wedge randomization",
                                        "Intermediate analysis with individual randomization"))),
      hidden(menuItem("hiddenChartsProp", tabName = "hiddenChartsProp")),
      menuItem("COMPARING SURVIVAL",tabName = "comparingsurv",id = "CSid",expandedName = "COMPARINGSURV",
               radioButtons(inputId = "survDesign",
                            label = "Design",
                            choices = c("No intermediate analysis with individual randomization", 
                                        "No intermediate analysis with stepped wedge randomization"))),
      hidden(menuItem("hiddenChartsSurv", tabName = "hiddenChartsSurv")),
      menuItem("BINARY EVENT PREDICTION", tabName = "bep", id = "BEPid", expandedName = "BEP",
               radioButtons(inputId = "typePredBEP",
                            label = "Design",
                            choices = c("Construction of predictive model of binary event", "External validation of binary event"))),
      hidden(menuItem("hiddenChartsBEP", tabName = "hiddenChartsBEP"))
      # menuItem("COMPARING TWO ROC CURVES", tabName = "CTRC")
    )
  ),
  dashboardBody(
    useShinyjs(),
    tabItems(
      
# DESCRIBE A POPULATION ---------------------------------------------------------------------------------------------------
tabItem(tabName = "hiddenChartsDP",
        fluidRow(column(width = 12,
                        h3("DESCRIBE A POPULATION"))),
        br(),
        sidebarLayout(
          sidebarPanel(
            conditionalPanel(
              condition = "input.typeDP == 'Mean'",
              h4("Describe mean with related IC :"),
              br(),
              numericInput(inputId = "descExpMean",
                           label = "Expected standard deviation of the random variable",
                           value = NULL,
                           max = 100,
                           min = 0)%>%
                shinyInput_label_embed(
                  shiny_iconlink() %>%
                    bs_embed_tooltip("WARNING : Not the standard deviation of the mean",
                                     placement = "right")),
              numericInput(inputId = "widthDescMean",
                           label = "Desired width of the confidence interval",
                           value = NULL,
                           max = 100,
                           min = 0)
            ),
            conditionalPanel(
              condition = "input.typeDP == 'Proportion'",
              h4("Describe proportion with related IC :"),
              br(),
              numericInput(inputId = "descExpProp",
                           label = "Expected proportion (%)",
                           value = NULL,
                           max = 100,
                           min = 0),
              numericInput(inputId = "widthDescProp",
                           label = "Desired width of the confidence interval (%)",
                           value = NULL,
                           max = 100,
                           min = 0)
            ),
            numericInput(inputId = "descAlpha",
                         label = "Type I error rate (%)",
                         value = 5)
          ),
          mainPanel(
            box(
              tags$script("
              Shiny.addCustomMessageHandler('txt', function (txt) {
              navigator.clipboard.writeText(txt);
              });"),
              div(style = "position:absolute;right:1em;",
              actionButton("copy_link_description", "Copy",icon = icon("copy"))),
              br(),br(),
              htmlOutput("description")))
        ),
        fluidRow(column(width = 12,
                          tags$a("Chow, S.-C., Shao, J., Wang, H., & Lokhnygina, Y. (2017). Sample Size Calculations in Clinical Research (3rd ed.). Chapman and Hall/CRC.", 
                                 href = "https://ndl.ethernet.edu.et/bitstream/123456789/15304/1/9103.pdf")))),
      
# PAGE COMPARING MEANS --------------------------------------------------------------------------------------------------
      tabItem(tabName = "hiddenChartsMean",
              fluidRow(column(width = 12,
                              h3("COMPARING MEANS"))),
              br(),
              sidebarLayout(
                sidebarPanel(
                  textOutput("titleMean"),
                  tags$head(tags$style(HTML("
                            #titleMean {
                              font-size: 18px;
                            }
                            "))),
                  br(),
                    conditionalPanel(
                      condition = "input.hypMean == 'Superiority'",
                      radioButtons(inputId = "sidetestMean",
                                   label = "Two sided test",
                                   choices = c("Yes", "No")),
                      numericInput(inputId = "prevMean",
                                   label = "Expected mean in the experimental arm",
                                   value = NULL,
                                   step = 0.005)),
                    numericInput(inputId = "obsMean",
                                 label = "Expected mean in the control arm",
                                 value = NULL,
                                 step = 0.005),
                  numericInput(inputId = "powerMean",
                               label = "Power (%)",
                               value = 80)%>%
                    shinyInput_label_embed(
                      shiny_iconlink() %>%
                        bs_embed_tooltip("1 - type II error",
                                         placement = "right")),
                  numericInput(inputId = "alphaMean",
                               label = "Type I error rate (%)",
                               value = 5,
                               step = 0.5),
                  numericInput(inputId = "sigmaMean",
                               label = "Expected standard deviation in the two arms",
                               value = NULL,
                               step = 0.005),
                  conditionalPanel(
                    condition = "input.hypMean == 'Non-inferiority'",
                    numericInput(inputId = "deltaMean",
                                 label = "Absolute non-inferiority margin",
                                 value = NULL,
                                 step = 0,005)),
                  conditionalPanel(
                    condition = "input.meanDesign == 'No intermediate analysis with stepped wedge randomization'",
                    numericInput(inputId = "centersMean",
                                 label = "Number of centers",
                                 value = NULL,
                                 step = 1),
                    numericInput(inputId = "sequencesMean",
                                 label = "Number of sequences",
                                 value = NULL,
                                 step = 1)%>%
                      shinyInput_label_embed(
                        shiny_iconlink() %>%
                          bs_embed_tooltip("a sequence dictate the order (or timing) at which each cluster will switch to the intervention condition, number of period of trtA and trtB",
                                           placement = "right")),
                    numericInput(inputId = "iccMean",
                                 label = "Intraclass correlation coefficient",
                                 value = 0.05,
                                 step = 0.05)),
                  conditionalPanel(
                    condition = "input.meanDesign == 'Intermediate analysis with individual randomization'",
                    numericInput(inputId = "meanSlice",
                                 label = "Number of planned intermediate analysis",
                                 value = 2))
                ), # sidebar panel
                mainPanel(
                  box(
                    tags$script("
                    Shiny.addCustomMessageHandler('txt', function (txt) {
                    navigator.clipboard.writeText(txt);
                    });"),
                    div(style = "position:absolute;right:1em;",
                        actionButton("copy_link_mean", "Copy",icon = icon("copy"))),
                    br(),br(),
                    htmlOutput("mean"))
                  )
              ),
              conditionalPanel(
                condition = "input.meanDesign == 'No intermediate analysis with individual randomization'",
                fluidRow(column(width = 12,
                                tags$a("Chow, S.-C., Shao, J., Wang, H., & Lokhnygina, Y. (2017). 
                                       Sample Size Calculations in Clinical Research (3rd ed.). 
                                       Chapman and Hall/CRC. "),
                                       href = "https://ndl.ethernet.edu.et/bitstream/123456789/15304/1/9103.pdf"))
                ),
              conditionalPanel(
                condition = "input.meanDesign == 'No intermediate analysis with stepped wedge randomization'",
                fluidRow(column(width = 12,
                                tags$a("Hemming K and Taljaard M. Sample size calculations for 
                                       stepped wedge and cluster randomised trials: a unified 
                                       approach. J Clin Epidemiol. 2016 Jan;69:137-46",
                                       href = "https://www.sciencedirect.com/science/article/pii/S089543561500414X")))
              ),
              conditionalPanel(
                condition = "input.meanDesign == 'Intermediate analysis with individual randomization'",
                fluidRow(column(width = 12,
                                tags$a("Demets DL and Lan KG. Interim analysis: The alpha spending 
                                       function approach. Stat Med. 1994;13(13-14):1341–1352",
                                       href = "https://pubmed.ncbi.nlm.nih.gov/7973215/")))
              )
      ),

# PAGE COMPARING PROPORTIONS --------------------------------------------------------------------------------------------
      tabItem(tabName = "hiddenChartsProp",
               fluidRow(column(width = 12,
                               h3("COMPARING PROPORTIONS"))),
               br(),
               sidebarLayout(
                 sidebarPanel(
                   textOutput("titleProp"),
                   tags$head(tags$style(HTML("
                            #titleProp {
                              font-size: 18px;
                            }
                            "))),
                   br(),
                   conditionalPanel(
                     condition = "input.hypProp == 'Superiority'",
                     radioButtons(inputId = "sidetestProp",
                                  label = "Two sided test",
                                  choices = c("Yes", "No")),
                     numericInput(inputId = "prevProp",
                                  label = "Expected proportion of the experimental arm (%)",
                                  value = NULL,
                                  step = 0.005)),
                     numericInput(inputId = "obsProp",
                                  label = "Expected proportion of the control arm (%)",
                                  value = NULL,
                                  step = 0.005),
                   numericInput(inputId = "alphaProp",
                                label = "Type I error rate (%)",
                                value = 5,
                                max = 20,
                                min = 0,
                                step = 0.5),
                   numericInput(inputId = "powerProp",
                                label = "Power (%)",
                                value = 80,
                                max = 100,
                                min = 0)%>%
                     shinyInput_label_embed(
                       shiny_iconlink() %>%
                         bs_embed_tooltip("1 - type II error",
                                          placement = "right")),
                   conditionalPanel(
                     condition = "input.propDesign == 'No intermediate analysis with stepped wedge randomization'",
                     numericInput(inputId = "centersProp",
                                  label = "Number of centers",
                                  value = NULL,
                                  step = 1),
                     numericInput(inputId = "sequencesProp",
                                  label = "Number of sequences",
                                  value = NULL,
                                  step = 1),
                     numericInput(inputId = "iccProp",
                                  label = "Intraclass correlation coefficient",
                                  value = 0.05,
                                  step = 0.05)),
                   conditionalPanel(
                     condition = "input.hypProp == 'Non-inferiority'",
                     numericInput(inputId = "deltaProp",
                                  label = "Absolute non-inferiority margin (%)",
                                  value = NULL,
                                  step = 0.005)),
                   conditionalPanel(
                     condition = "input.propDesign == 'Intermediate analysis with individual randomization'",
                     numericInput(inputId = "propSlice",
                                  label = "Number of planned intermediate analysis",
                                  value = 2))
                 ),
                 mainPanel(
                   box(
                     tags$script("
                     Shiny.addCustomMessageHandler('txt', function (txt) {
                     navigator.clipboard.writeText(txt);
                     });"),
                     div(style = "position:absolute;right:1em;",
                         actionButton("copy_link_proportion", "Copy",icon = icon("copy"))),
                     br(),br(),
                     htmlOutput("proportion"))
                   )
               ), # sidebarLayout
              conditionalPanel(
                condition = "input.propDesign == 'No intermediate analysis with individual randomization'",
                fluidRow(column(width = 12,
                                tags$a("Chow, S.-C., Shao, J., Wang, H., & Lokhnygina, Y. (2017). 
                                       Sample Size Calculations in Clinical Research (3rd ed.). 
                                       Chapman and Hall/CRC. ",
                                       href = "https://ndl.ethernet.edu.et/bitstream/123456789/15304/1/9103.pdf")))
              ),
              conditionalPanel(
                condition = "input.propDesign == 'No intermediate analysis with stepped wedge randomization'",
                fluidRow(column(width = 12,
                                tags$a("Hemming K, Taljaard M. Sample size calculations for 
                                       stepped wedge and cluster randomised trials: a unified 
                                       approach. J Clin Epidemiol. 2016 Jan;69:137-46",
                                       href = "https://www.sciencedirect.com/science/article/pii/S089543561500414X")))
              ),
              conditionalPanel(
                condition = "input.propDesign == 'Intermediate analysis with individual randomization'",
                fluidRow(column(width = 12,
                                tags$a("Demets DL and Lan KG. Interim analysis: The alpha spending 
                                       function approach. Stat Med. 1994;13(13-14):1341–1352",
                                       href = "https://pubmed.ncbi.nlm.nih.gov/7973215/")))
              )
      ),# Comparing proportion


# PAGE COMPARING SURVIVAL --------------------------------------------------------------------------------------------
tabItem(tabName = "hiddenChartsSurv",
        fluidRow(column(width = 12,
                        h3("COMPARING SURVIVAL"))),
        br(),
        sidebarLayout(
          sidebarPanel(
            textOutput("titleSurv"),
            tags$head(tags$style(HTML("
                            #titleSurv {
                              font-size: 18px;
                            }
                            "))),
            br(),
            numericInput(inputId = "prevSurv",
                         label = "Expected survival proportion of the experimental arm (%)",
                         value = NULL,
                         step = 0.005),
            numericInput(inputId = "obsSurv",
                         label = "Expected survival proportion of the control arm (%)",
                         value = NULL,
                         step = 0.005),
            numericInput(inputId = "alphaSurv",
                         label = "Type I error rate (%)",
                         value = 5,
                         max = 20,
                         min = 0,
                         step = 0.5),
            numericInput(inputId = "powerSurv",
                         label = "Power (%)",
                         value = 80,
                         max = 100,
                         min = 0)%>%
              shinyInput_label_embed(
                shiny_iconlink() %>%
                  bs_embed_tooltip("1 - type II error",
                                   placement = "right")),
            conditionalPanel(
              condition = "input.survDesign == 'No intermediate analysis with stepped wedge randomization'",
              numericInput(inputId = "centersSurv",
                           label = "Number of centers",
                           value = NULL,
                           step = 1),
              numericInput(inputId = "sequencesSurv",
                           label = "Number of sequences",
                           value = NULL,
                           step = 1),
              numericInput(inputId = "iccSurv",
                           label = "Intraclass correlation coefficient",
                           value = 0.05,
                           step = 0.05))
          ),
          mainPanel(
            box(
              tags$script("
                     Shiny.addCustomMessageHandler('txt', function (txt) {
                     navigator.clipboard.writeText(txt);
                     });"),
              div(style = "position:absolute;right:1em;",
                  actionButton("copy_link_survival", "Copy",icon = icon("copy"))),
              br(),br(),
              htmlOutput("survival"))
          )
        ), # sidebarLayout
        
        fluidRow(column(width = 12,
                          tags$a("Machin, D., Campbell, M.J., Tan, S.B., Tan, S.H. (2009). 
                                      Sample Size Tables for Clinical Studies (3rd ed.). 
                                      Blackwell Publishing",
                                 href = "https://onlinelibrary.wiley.com/doi/book/10.1002/9781444300710"))),
        conditionalPanel(
          condition = "input.survDesign == 'No intermediate analysis with stepped wedge randomization'",
          fluidRow(column(width = 12,
                          tags$a("Hemming K, Taljaard M. Sample size calculations for 
                                       stepped wedge and cluster randomised trials: a unified 
                                       approach. J Clin Epidemiol. 2016 Jan;69:137-46",
                                 href = "https://www.sciencedirect.com/science/article/pii/S089543561500414X")))
        )
),# Comparing survival



# PAGE BINARY EVENT PREDICTION ------------------------------------------------------------------------------------------
tabItem(tabName = "hiddenChartsBEP",
         fluidRow(column(width = 12,
                         h3("BINARY EVENT PREDICTION"))),
         br(),
         sidebarLayout(
           sidebarPanel(
             conditionalPanel(
               condition = "input.typePredBEP == 'Construction of predictive model of binary event'",
               h4("Construction of predictive model :"),
               br(),
               numericInput(inputId = "predictorsToTest",
                            label = "Number of potential predictors",
                            value = NULL,
                            step = 1),
               numericInput(inputId = "shrinkageExpected",
                            label = "Expected shrinkage",
                            value = NULL,
                            step = 0.1)%>%
                 shinyInput_label_embed(
                   shiny_iconlink() %>%
                     bs_embed_tooltip("i.e. reduction in the effects of sampling variation",
                                      placement = "right")),
               numericInput(inputId = "R2",
                            label = "Expected predictive capacities in term of R2",
                            value = NULL,
                            step = 0.01)
             ),
             conditionalPanel(
               condition = "input.typePredBEP == 'External validation of binary event'",
               h4("External validation :"),
               br(),
               numericInput(inputId = "P0",
                            label = "Expected proportion of events (%)",
                            value = NULL,
                            step = 1),
               numericInput(inputId = "alphaExtVal",
                            label = "Type I error rate (%)",
                            value = 5,
                            step = 0.5),
               #The total length of the confidence of O/E
               numericInput(inputId = "widthExtVal",
                            label = "Expected width of the confidence interval (%)",
                            value = NULL,
                            step = 0.1)
             )
           ),
           mainPanel(
             box(
               tags$script("
               Shiny.addCustomMessageHandler('txt', function (txt) {
               navigator.clipboard.writeText(txt);
               });"),
               div(style = "position:absolute;right:1em;",
                   actionButton("copy_link_prediction", "Copy",icon = icon("copy"))),
               br(),br(),
               htmlOutput("prediction"))
             )
         ),
        conditionalPanel(
          condition = "input.typePredBEP == 'Construction of predictive model of binary event'",
          fluidRow(column(width = 12,
                          tags$a("Riley et al. Statistics in Medicine. equation 11,  
                                 2019;38:1276–1296",
                                 href = "https://onlinelibrary.wiley.com/doi/full/10.1002/sim.7992")))
        ),
        conditionalPanel(
          condition = "input.typePredBEP == 'External validation of binary event'",
          fluidRow(column(width = 12,
                          tags$a("Riley et al. Minimum sample size for external validation 
                                 of a clinical prediction model with a binary outcome. 
                                 Statistics in Medicine. 2021;19:4230-4251",
                                 href = "https://onlinelibrary.wiley.com/doi/full/10.1002/sim.9025")))
        )
), # BINARY EVENT PREDICTION



# PAGE COMPARING TWO ROC CURVES -----------------------------------------------------------------------------------------
tabItem(tabName = "CTRC",
         fluidRow(column(width = 12,
                         h3("COMPARING TWO ROC CURVES"))),
         br(),
         sidebarLayout(
           sidebarPanel(
             numericInput(inputId = "AUC1",
                          label = "AUC of the first diagnostic test",
                          value = NULL,
                          min = 0.5,
                          step = 0.001),
             numericInput(inputId = "AUC2",
                          label = "AUC of the second diagnostic test",
                          value = NULL,
                          min = 0.5,
                          step = 0.001),
             numericInput(inputId = "alphaAUC",
                          label = "Type I error rate (%)",
                          value = 5,
                          step = 0.5),
             numericInput(inputId = "powerAUC",
                          label = "Power (%)",
                          value = 80,
                          min = 0,
                          max = 100)
           ),
           mainPanel(
             box(tags$script("
                 Shiny.addCustomMessageHandler('txt', function (txt) {
                 navigator.clipboard.writeText(txt);
                 });"),
                 div(style = "position:absolute;right:1em;",
                     actionButton("copy_link_AUC", "Copy",icon = icon("copy"))),
                 br(),br(),
               htmlOutput("AUC")) 
             )
         )
) #ROC CURVES
    )
  )
)


################################################################################################################################
################################################### SERVER #####################################################################
################################################################################################################################


server <- function(input, output, session) {
  
  # for the sidbar --------------------------------------------------------
  observeEvent(input$sidebarItemExpanded, {
    if(input$sidebarItemExpanded == "DP"){
      updateTabItems(session, "sidebarID", selected = "hiddenChartsDP")
    }
    if(input$sidebarItemExpanded == "COMPARINGMEAN"){
      updateTabItems(session, "sidebarID", selected = "hiddenChartsMean")
    }
    if(input$sidebarItemExpanded == "COMPARINGPROP"){
      updateTabItems(session, "sidebarID", selected = "hiddenChartsProp")
    }
    if(input$sidebarItemExpanded == "COMPARINGSURV"){
      updateTabItems(session, "sidebarID", selected = "hiddenChartsSurv")
    }
    if(input$sidebarItemExpanded == "BEP"){
      updateTabItems(session, "sidebarID", selected = "hiddenChartsBEP")
    }
  })
  #------------------------------------------------------------------------
  
  # Title of the page
  output$titleMean <- renderText({paste0(input$hypMean, " trial. ",input$meanDesign, " :")})
  output$titleProp <- renderText({paste0(input$hypProp, " trial. ",input$propDesign, " :")})
  output$titleSurv <- renderText({paste0("Superiority trial. ",input$survDesign, " :")})
  
  
  
  sampleSizeDesc <- function(p,alpha,width){
    Z <- qnorm(1-alpha/2)
    (((2*Z)**2)*(p*(1-p)))/(width**2)
  }
  
  # convert to numeric value -----------------------------------------------------------------------------------------------------
  # "one-sided" <- 1
  # "tow-sided" <- 2
  stMean <- reactive( if (reactive(input$sidetestMean)()=="No") 1 else 2 )
  printStMean <- reactive( if (reactive(input$sidetestMean)()=="No") "one-sided" else "two-sided" )
  
  stProp <- reactive( if (reactive(input$sidetestProp)()=="No") 1 else 2 )
  printStProp <- reactive( if (reactive(input$sidetestProp)()=="No") "one-sided" else "two-sided" )
  
  
  varAUC <- function(auc){
    a <- norminv(auc)*1.414
    (0.0099*exp(-(a**2)/2))*(6*(a**2)+16)
  }
  
  compAUC <- function(auc1=input$AUC1, auc2=input$AUC2, alpha=input$alphaAUC/100, power=input$powerAUC/100){
    ((qnorm(1-alpha/2)*sqrt(2*varAUC((auc1+auc2)/2))+qnorm(power)*sqrt(varAUC(auc1)+varAUC(auc2)))**2)/((auc2-auc1)**2)
  }
  
  
  
  
  
  
  
  
  
  # sample size result of the computation ----------------------------------------------------------------------------------------
  
  #DESCRITPION -------------------------------------------------------------------------------------------------------------------
  resDesc <- reactive(
    if(reactive(input$typeDP)()=='Mean'){
      ceiling((2*qnorm(1-(input$descAlpha/100)/2)*input$descExpMean/input$widthDescMean)**2)
    }else{
      ceiling(sampleSizeDesc(p=input$descExpProp/100, alpha=input$descAlpha/100, width=input$widthDescProp/100))
    }
  )
  
  # MEAN -------------------------------------------------------------------------------------------------------------------------
  resMeanSeq <- reactive(
      if(reactive(input$meanDesign)()=='Intermediate analysis with individual randomization' & reactive(input$hypMean)()=='Superiority'){
        design <- getDesignGroupSequential(typeOfDesign = "OF", informationRates = seq(from=1/(input$meanSlice+1), to=1, by=1/(input$meanSlice+1)),
                                           alpha = input$alphaMean/100, beta = 1-input$powerMean/100, sided = stMean())
        designPlan <- getSampleSizeMeans(design, alternative = abs(input$prevMean-input$obsMean), stDev = input$sigmaMean,
                                         allocationRatioPlanned = 1)
        designPlan$numberOfSubjects
        
      }else if(reactive(input$meanDesign)()=='Intermediate analysis with individual randomization' & reactive(input$hypMean)()=='Non-inferiority'){
        design <- getDesignGroupSequential(typeOfDesign = "OF", informationRates = seq(from=1/(input$meanSlice+1), to=1, by=1/(input$meanSlice+1)),
                                           alpha = input$alphaMean/100, beta = 1-input$powerMean/100, sided = 1)
        designPlan <- getSampleSizeMeans(design, alternative = 0, stDev = input$sigmaMean,
                                         allocationRatioPlanned = 1, thetaH0 = -input$deltaMean)
        designPlan$numberOfSubjects
      }
  )
  
  

  
  resMean <- reactive(
    # MEAN SUPERIORITY NO INTERMEDIATE ANALYSIS AND INDIVIDUAL RANDOMIZATION  --------------------
    if(reactive(input$hypMean)()=='Superiority' & reactive(input$meanDesign)()=='No intermediate analysis with individual randomization'){ 
        epi.sscompc(N = NA, treat = input$prevMean, control = input$obsMean, 
                    sigma = input$sigmaMean, n = NA, power = input$powerMean/100, 
                    r = 1, design = 1, sided.test = stMean(), conf.level = 1-(input$alphaMean/100))$n.treat
    # MEAN SUPERIORITY NO INTERMEDIATE ANALYSIS AND PARALLEL CLUSTER RANDOMIZATION  -------------------
    }else if(reactive(input$hypMean)()=='Superiority' & reactive(input$meanDesign)()=='No intermediate analysis with parallel cluster randomization'){
      SampSize_I <- epi.sscompc(N = NA, treat = input$prevMean, control = input$obsMean, 
                                sigma = input$sigmaMean, n = NA, power = input$powerMean/100, 
                                r = 1, design = 1, sided.test = stMean(), conf.level = 1-(input$alphaMean/100))
      CRT(ni=SampSize_I$n.total, center=input$centersMean, icc=input$iccMean)/2
    # MEAN SUPERIORITY NO INTERMEDIATE ANALYSIS AND STEP WEDGE RANDOMIZATION -------------------------
    }else if(reactive(input$hypMean)()=='Superiority' & reactive(input$meanDesign)()=='No intermediate analysis with stepped wedge randomization'){
        SampSize_I <- epi.sscompc(N = NA, treat = input$prevMean, control = input$obsMean, 
                                  sigma = input$sigmaMean, n = NA, power = input$powerMean/100, 
                                  r = 1, design = 1, sided.test = stMean(), conf.level = 1-(input$alphaMean/100))
        SW(ni=SampSize_I$n.total, center=input$centersMean, sequence=input$sequencesMean, icc=input$iccMean)/2
      # MEAN SUPERIORITY INTERMEDIATE ANALYSIS AND INDIVIDUAL RANDOMIZATION ------------------------------------------------------
      }else if(reactive(input$hypMean)()=='Superiority' & reactive(input$meanDesign)()=='Intermediate analysis with individual randomization'){
        ceiling(resMeanSeq()[input$meanSlice+1]/2) 
      
      # MEAN NON-INFERIORITY NO INTERMEDIATE ANALYSIS AND INDIVIDUAL RANDOMIZATION ----------------
      }else if(reactive(input$hypMean)()=='Non-inferiority' & reactive(input$meanDesign)()=='No intermediate analysis with individual randomization'){
        epi.ssninfc(treat = input$obsMean, control = input$obsMean, sigma = input$sigmaMean,
                    delta = input$deltaMean, n = NA, power = input$powerMean/100, alpha = input$alphaMean/100, r = 1)$n.treat
      # MEAN NON-INFERIORITY NO INTERMEDIATE ANALYSIS AND PARALLEL CLUSTER RANDOMIZATION  -------------------
      }else if(reactive(input$hypMean)()=='Non-inferiority' & reactive(input$meanDesign)()=='No intermediate analysis with parallel cluster randomization'){
        SampSize_I <- epi.ssninfc(treat = input$obsMean, control = input$obsMean, sigma = input$sigmaMean,
                                  delta = input$deltaMean, n = NA, power = input$powerMean/100, alpha = input$alphaMean/100, r = 1)
        CRT(ni=SampSize_I$n.total, center=input$centersMean, icc=input$iccMean)/2
      # MEAN NON-INFERIORITY NO INTERMEDIATE ANALYSIS AND STEP WEDGE RANDOMIZATION ------------------------
      }else if(reactive(input$hypMean)()=='Non-inferiority' & reactive(input$meanDesign)()=='No intermediate analysis with stepped wedge randomization'){ 
        SampSize_I <- epi.ssninfc(treat = input$obsMean, control = input$obsMean, sigma = input$sigmaMean,
                                  delta = input$deltaMean, n = NA, power = input$powerMean/100, alpha = input$alphaMean/100, r = 1)
        SW(ni=SampSize_I$n.total, center=input$centersMean, sequence=input$sequencesMean, icc=input$iccMean)/2
      # MEAN NON-INFERIORITY Intermediate analysis with individual randomization -----------------------------------------------------------------------------
      }else if(reactive(input$hypMean)()=='Non-inferiority' & reactive(input$meanDesign)()=='Intermediate analysis with individual randomization'){
        ceiling(resMeanSeq()[input$meanSlice+1]/2)
      })
  
  
  # PROPORTION ------------------------------------------------------------------------------------------------------------------
  resPropSeq <- reactive(
    # PROPORTION INTERMEDIATE ANALYSIS WITH INDIVIDUAL RANDOMIZATION ------------------------
    if(reactive(input$propDesign)()=='Intermediate analysis with individual randomization' & reactive(input$hypProp)()=='Superiority'){
        design <- getDesignGroupSequential(typeOfDesign = "OF", informationRates = seq(from=1/(input$propSlice+1), to=1, by=1/(input$propSlice+1)),
                                           alpha = input$alphaProp/100, beta = 1-input$powerProp/100, sided = stProp())
        designPlan <- getSampleSizeRates(design,  pi1 = input$prevProp/100, pi2 = input$obsProp/100,
                                         allocationRatioPlanned = 1)
        designPlan$numberOfSubjects
      }else if(reactive(input$propDesign)()=='Intermediate analysis with individual randomization' & reactive(input$hypProp)()=='Non-inferiority'){
        design <- getDesignGroupSequential(typeOfDesign = "OF", informationRates = seq(from=1/(input$propSlice+1), to=1, by=1/(input$propSlice+1)),
                                           alpha = input$alphaProp/100, beta = 1-input$powerProp/100, sided = 1)
        designPlan <- getSampleSizeRates(design, pi1 = input$obsProp/100, pi2 = input$obsProp/100, thetaH0 = input$deltaProp/100)
        designPlan$numberOfSubjects
      })
  
  
  resProp <- reactive(
    # PROPORTION SUPERIORITY NO INTERMEDIATE ANALYSIS : INDIVIDUAL RANDOMIZATION -------------
    if(reactive(input$hypProp)()=='Superiority' & reactive(input$propDesign)()=='No intermediate analysis with individual randomization'){
        epi.sscohortc(N = NA, irexp1 = input$prevProp/100, irexp0 = input$obsProp/100, pexp = NA, n = NA, 
                      power = input$powerProp/100, r = 1, design = 1, sided.test = stProp(), 
                      finite.correction = FALSE, nfractional = FALSE, conf.level = 1-(input$alphaProp/100))$n.exp1
    # PROPORTION SUPERIORITY NO INTERMEDIATE ANALYSIS : PARALLEL CLUSTER RANDOMIZATION --------------
    }else if(reactive(input$hypProp)()=='Superiority' & reactive(input$propDesign)()=='No intermediate analysis with parallel cluster randomization'){
      SampSize_I <- epi.sscohortc(irexp1 = input$prevProp/100, irexp0 = input$obsProp/100, n = NA, r = 1,
                                  power = input$powerProp/100, sided.test = stProp(), 
                                  conf.level = 1-(input$alphaProp/100))
      CRT(ni=SampSize_I$n.total, center=input$centersProp, icc=input$iccProp)/2
    # PROPORTION SUPERIORITY NO INTERMEDIATE ANALYSIS : STEP WEDGE RANDOMIZATION ---------------------------
    }else if(reactive(input$hypProp)()=='Superiority' & reactive(input$propDesign)()=='No intermediate analysis with stepped wedge randomization'){
      SampSize_I <- epi.sscohortc(irexp1 = input$prevProp/100, irexp0 = input$obsProp/100, n = NA, r = 1,
                                    power = input$powerProp/100, sided.test = stProp(), 
                                    conf.level = 1-(input$alphaProp/100))
      SW(ni=SampSize_I$n.total, center=input$centersProp, sequence=input$sequencesProp, icc=input$iccProp)/2
    # PROPORTION SUPERIORITY INTERMEDIATE ANALYSIS AND INDIVIDUAL RANDOMIZATION -------------------------------------------------------------------------
    }else if(reactive(input$hypProp)()=='Superiority' & reactive(input$propDesign)()=='Intermediate analysis with individual randomization'){
        ceiling(resPropSeq()[input$propSlice+1]/2)
        
    # PROPORTION NON-INFERIORITY NO INTERMEDIATE ANALYSIS : INDIVIDUAL RANDOMIZATION ---------
    }else if(reactive(input$hypProp)()=='Non-inferiority' & reactive(input$propDesign)()=='No intermediate analysis with individual randomization'){
        epi.ssninfb(treat = input$obsProp/100, control = input$obsProp/100, delta = input$deltaProp/100, 
                    n = NA, r = 1, power = input$powerProp/100, alpha = input$alphaProp/100)$n.treat
    # PROPORTION NON-INFERIORITY NO INTERMEDIATE ANALYSIS : PARALLEL CLUSTER RANDOMIZATION -----------------
    }else if(reactive(input$hypProp)()=='Non-inferiority' & reactive(input$propDesign)()=='No intermediate analysis with parallel cluster randomization'){
      SampSize_I <- epi.ssninfb(treat = input$obsProp/100, control = input$obsProp/100, delta = input$deltaProp/100, 
                                n = NA, r = 1, power = input$powerProp/100, alpha = input$alphaProp/100)
      CRT(ni=SampSize_I$n.total, center=input$centersProp, icc=input$iccProp)/2
    # PROPORTION NON-INFERIORITY NO INTERMEDIATE ANALYSIS : STEP WEDGE RANDOMIZATION -----------------
    }else if(reactive(input$hypProp)()=='Non-inferiority' & reactive(input$propDesign)()=='No intermediate analysis with stepped wedge randomization'){
        SampSize_I <- epi.ssninfb(treat = input$obsProp/100, control = input$obsProp/100, delta = input$deltaProp/100, 
                                  n = NA, r = 1, power = input$powerProp/100, alpha = input$alphaProp/100)
        SW(ni=SampSize_I$n.total, center=input$centersProp, sequence=input$sequencesProp, icc=input$iccProp)/2
    # PROPORTION NON-INFERIORITY INTERMEDIATE ANALYSIS AND INDIVIDUAL RANDOMIZATION ---------------------------------------------------------------------  
    }else if(reactive(input$hypProp)()=='Non-inferiority' & reactive(input$propDesign)()=='Intermediate analysis with individual randomization'){
       ceiling(resPropSeq()[input$propSlice+1]/2)
    })
  
  
  
  resSurv <- reactive(
    # SURVIVAL SUPERIORITY NO INTERMEDIATE ANALYSIS : INDIVIDUAL RANDOMIZATION -------------
    if(reactive(input$survDesign)()=='No intermediate analysis with individual randomization'){
      sample_survival(p1=input$prevSurv/100, p0=input$obsSurv/100, 
                      power=input$powerSurv/100, r=1, alpha=input$alphaSurv/100)$n.gp1

    }
    # SURVIVAL SUPERIORITY NO INTERMEDIATE ANALYSIS : STEP WEDGE RANDOMIZATION ---------------------------
    else if(reactive(input$survDesign)()=='No intermediate analysis with stepped wedge randomization'){
      SampSize_I <- sample_survival(p1=input$prevSurv/100, p0=input$obsSurv/100, 
                                    power=input$powerSurv/100, r=1, alpha=input$alphaSurv/100)
      SW(ni=SampSize_I$n.total, center=input$centersSurv, sequence=input$sequencesSurv, icc=input$iccSurv)/2
    }
  )
  
  
  # PREDICTIVE ------------------------------------------------------------------------------------------------------------------
  resPred <- reactive(
    if(reactive(input$typePredBEP)()=='Construction of predictive model of binary event')
      ceiling(input$predictorsToTest/((input$shrinkageExpected-1)*log(1-input$R2/input$shrinkageExpected)))
    else if (reactive(input$typePredBEP)()=='External validation of binary event')
      ceiling(size.calib(p0=input$P0/100, width = input$widthExtVal/100, alpha = input$alphaExtVal/100))
  )
  
  # AUC -------------------------------------------------------------------------------------------------------------------------
  resAUC <- reactive(
    ceiling(compAUC())
  )
  
  # reactive example sentence ---------------------------------------------------------------------------------------------------
  z <- reactive(
    if(reactive(input$typeDP)()=='Mean'){
      paste0("This sample size is for a continuous endpoint descriptive study. In order to describe a mean for an outcome with an expected standard deviation of ",input$descExpMean," with a precision define by a ",input$widthDescMean," total width confidence interval and a ",input$descAlpha,"% two-sided type I error rate, the minimum sample size needed is ",resDesc()," patients")
    }else{
      paste0("This sample size is for a binary endpoint descriptive study. In order to demonstrate the expected proportion of event of ",input$descExpProp,"% with a precision define by a ",input$widthDescProp,"% width confidence interval and a ",input$descAlpha,"% two-sided type I error rate, the minimum sample size needed is ",resDesc()," patients")
    })
  
  
  a <- reactive(
      # MEAN SUPERIORITY NO INTERMEDIATE ANALYSIS AND INDIVIDUAL RANDOMIZATION  --------------------
      if(reactive(input$hypMean)()=='Superiority' & reactive(input$meanDesign)()=='No intermediate analysis with individual randomization'){
        if (resMean()<30) {
          warning("--> At least one group size is < 30, normality assumption is questionnable and sample 
                size calculation may not be valid.")
        }else
        paste0("This sample size is for a randomized superiority trial in two parallel groups (experimental treatment versus control treatment) with balanced randomization (ratio 1:1) for a continuous endpoint. The expected mean of the criteria is ",input$prevMean," units with experimental treatment compared to ",input$obsMean," units with control treatment. In order to highlight this absolute difference of ",abs(input$prevMean-input$obsMean),", with a standard deviation of ",input$sigmaMean,", with a ",printStMean()," alpha risk of ",input$alphaMean,"% and a power of ",input$powerMean,"%, the needed sample size is ",resMean()*2, " patients (",resMean()," patients in each group).")
      # MEAN SUPERIORITY NO INTERMEDIATE ANALYSIS AND PARALLEL CLUSTER RANDOMIZATION ---------------------
      }else if(reactive(input$hypMean)()=='Superiority' & reactive(input$meanDesign)()=='No intermediate analysis with parallel cluster randomization'){
        if (resMean()<30) {
          warning("--> At least one group size is < 30, normality assumption is questionnable and sample 
                size calculation may not be valid.")
        }else
          paste0("This sample size is for a parallel cluster randomized superiority trial in two parallel groups (experimental treatment versus control treatment) with balanced randomization (ratio 1:1) for a continuous endpoint. The expected mean of the criteria is ",input$prevMean," units with experimental treatment compared to ",input$obsMean," units with control treatment. In order to highlight this absolute difference of ",abs(input$prevMean-input$obsMean),", with a standard deviation of ",input$sigmaMean,", with a ",printStMean()," alpha risk of ",input$alphaMean,"%, a power of ",input$powerMean,"%, according to our parallel cluster design with ",input$centersMean," centers randomized and assuming an intraclass correlation coefficient of ",input$iccMean,", the needed sample size is ",resMean()*2, " patients (",resMean()," patients in each group)")
        # MEAN SUPERIORITY NO INTERMEDIATE ANALYSIS AND STEP WEDGE RANDOMIZATION ---------------------
      }else if(reactive(input$hypMean)()=='Superiority' & reactive(input$meanDesign)()=='No intermediate analysis with stepped wedge randomization'){
        if (resMean()<30) {
          warning("--> At least one group size is < 30, normality assumption is questionnable and sample 
                size calculation may not be valid.")
        }else
        paste0("This sample size is for a stepped wedge cluster randomized superiority trial in two groups (experimental treatment versus control treatment) with balanced randomization (ratio 1:1) for a continuous endpoint. The expected mean of the criteria is ",input$prevMean," units with experimental treatment compared to ",input$obsMean," units with control treatment. In order to highlight this absolute difference of ",abs(input$prevMean-input$obsMean),", with a standard deviation of ",input$sigmaMean,", with a ",printStMean()," alpha risk of ",input$alphaMean,"%, a power of ",input$powerMean,"%, according to our stepped wedge design with ",input$centersMean," centers randomized in ",input$sequencesMean," sequences and assuming an intraclass correlation coefficient of ",input$iccMean,", the needed sample size is ",resMean()*2, " patients (",resMean()," patients in each group)")
      # MEAN NON-INFERIORITY INTERMEDIATE ANALYSIS AND INDIVIDUAL RANDOMIZATION ------------------------------------------------------
      }else if(reactive(input$hypMean)()=='Superiority' & reactive(input$meanDesign)()=='Intermediate analysis with individual randomization'){
        paste0("This sample size is for a RCT with an expected mean of ",input$prevMean," units in patients in the experimental arm versus ",input$obsMean," units in the control arm. In order to demonstrate such a difference of ",abs(input$prevMean-input$obsMean)," units, with a standard deviation of ",input$sigmaMean,", a ",input$alphaMean,"% ",printStMean()," type I error rate and a power of ",input$powerMean,"%,  the final analysis should be carried out on ",resMean()*2," patients (",resMean()," patients per group). Intermediates analyses would be performed on ",paste(ceiling(head(resMeanSeq(),-2)), collapse = ", ")," and ",paste(tail(ceiling(head(resMeanSeq(),-1)),1), collapse = ", ")," patients respectively if their is nodecision of stopping the study.")
        
      # MEAN NON-INFERIORITY NO INTERMEDIATE ANALYSIS AND INDIVIDUAL RANDOMIZATION ----------------
      }else if(reactive(input$hypMean)()=='Non-inferiority' & reactive(input$meanDesign)()=='No intermediate analysis with individual randomization'){
        if (resMean()<30) {
          warning("--> At least one group size is < 30, normality assumption is questionnable and sample 
                size calculation may not be valid.")
        }else
        paste0("This sample size is for a randomized controlled non-inferiority trial in two parallel groups (experimental treatment versus control treatment) with balanced randomisation (ratio 1:1) for a continuous endpoint. The mean of the criteria is ",input$obsMean," with both treatment groups. Assuming an absolute non-inferiority margin of ",input$deltaMean,", with a standard deviation of ",input$sigmaMean,", with a one-sided type I error rate of ",input$alphaMean,"% and a power of ",input$powerMean,"%, the needed sample size is ",resMean()*2, "patients (",resMean()," patients in each group).")
      # MEAN NON-INFERIORITY NO INTERMEDIATE ANALYSIS AND PARALLEL CLUSTER RANDOMIZATION ---------------------
      }else if(reactive(input$hypMean)()=='Non-inferiority' & reactive(input$meanDesign)()=='No intermediate analysis with parallel cluster randomization'){
        if (input$prevMean<30 | input$obsMean<30) {
          warning("--> At least one group size is < 30, normality assumption is questionnable and sample 
                size calculation may not be valid.")
        }else
          paste0("This sample size for a parallel cluster randomised non-inferiority trial in two parallel groups (experimental treatment versus control treatment) with balanced randomisation (ratio 1:1) for a continuous endpoint.The mean of the criteria is ",input$obsMean," with both treatment groups. Assuming an absolute non-inferiority margin of ",input$deltaMean,", with a standard deviation of ",input$sigmaMean,", with a one-sided type I error rate of ",input$alphaMean,"%, a power of ",input$powerMean,"%, according to our parallel cluster design with ",input$centersMean,"centers randomized in two groups and assuming an intraclass correlation coefficient of ",input$iccMean,", the needed sample size is ",resMean()*2, " patients (",resMean()," patients in each group)")
      # MEAN NON-INFERIORITY NO INTERMEDIATE ANALYSIS AND STEP WEDGE RANDOMIZATION ---------------------
      }else if(reactive(input$hypMean)()=='Non-inferiority' & reactive(input$meanDesign)()=='No intermediate analysis with stepped wedge randomization'){
        if (input$prevMean<30 | input$obsMean<30) {
          warning("--> At least one group size is < 30, normality assumption is questionnable and sample 
                size calculation may not be valid.")
        }else
        paste0("This sample size for a stepped wedge cluster randomised non-inferiority trial in two groups (experimental treatment versus control treatment) with balanced randomisation (ratio 1:1) for a continuous endpoint.The mean of the criteria is ",input$obsMean," with both treamtent groups. Assuming an absolute non-inferiority margin of ",input$deltaMean,", with a standard deviation of ",input$sigmaMean,", with a one-sided type I error rate of ",input$alphaMean,"% and a power of ",input$powerMean,"%, according to our stepped wedge design with ",input$centersMean,"centers randomized in ",input$sequencesMean," sequences and assuming an intraclass correlation coefficient of ",input$iccMean,", the needed sample size is ",resMean()*2, " patients (",resMean()," patients in each group)")
      # MEAN NON-INFERIORITY INTERMEDIATE ANALYSIS AND INDIVIDUAL RANDOMIZATION -----------------------------------------------------------------------------
      }else if(reactive(input$hypMean)()=='Non-inferiority' & reactive(input$meanDesign)()=='Intermediate analysis with individual randomization'){
        paste0("This sample size is for a randomised controlled non-inferiority trial in two parallel groups (experimental treatment versus control treatment) with balanced randomization (ratio 1:1) for a continuous endpoint. Assuming an absolute non-inferiority margin of ",input$deltaMean,", with a standard deviation of ",input$sigmaMean,", with a one-sided alpha risk of ",input$alphaMean,"% and a power of ",input$powerMean,"%, the final analysis should be carried out on ",resMean()*2," patients (",resMean()," patients per group). The intermediate analyses would be performed on ",paste(ceiling(head(resMeanSeq(),-2)), collapse = ", ")," and ",paste(tail(ceiling(head(resMeanSeq(),-1)),1), collapse = ", ")," patients respectively, if their is no decision of stopping the study.")
      })
#----  
  
  b <- reactive(
    # PROPORTION SUPERIORITY NO INTERMEDIATE ANALYSIS AND INDIVIDUAL RANDOMIZATION  --------------------
    if(reactive(input$hypProp)()=='Superiority' & reactive(input$propDesign)()=='No intermediate analysis with individual randomization'){
      if (resProp()<30) {
        warning("--> At least one group size is < 30, normality assumption is questionnable and sample
                size calculation may not be valid.")
      }else
      paste0("This sample size is for a randomized controlled superiority trial in two parallel groups (experimental treatment versus control treatment) with balanced randomization (ratio 1:1) for a binary endpoint. The proportion of patient with the criteria is ",input$prevProp,"% with experimental treatment compared to ",input$obsProp,"% with control treatment. In order to highlight this absolute difference of ",abs(input$prevProp-input$obsProp),"%, with a ",printStProp()," type I error rate of ",input$alphaProp,"% and a power of ",input$powerProp,"%, the needed sample size is ",resProp()*2, " patients (",resProp()," patients in each group).")
    # PROPORTION SUPERIORITY NO INTERMEDIATE ANALYSIS AND PARALLEL CLUSTER RANDOMIZATION --------------------
    }else if(reactive(input$hypProp)()=='Superiority' & reactive(input$propDesign)()=='No intermediate analysis with parallel cluster randomization'){
      if (resProp()<30) {
        warning("--> At least one group size is < 30, normality assumption is questionnable and sample
                size calculation may not be valid.")
      }else
        paste0("This sample size is for a parallel cluster randomized superiority trial in two parallel groups (experimental treatment versus control treatment) with balanced randomization (ratio 1:1) for a binary endpoint. The expected propotion of event is ",input$prevProp,"% with experimental treatment compared to ",input$obsProp,"% with control treatment. In order to highlight this absolute difference of ",abs(input$prevProp-input$obsProp),"%, with a ",printStProp()," type I error rate of ",input$alphaProp,"%, a power of ",input$powerProp,"%, according to our parallel cluster design with ",input$centersProp," centers randomized in two groups and assuming an intraclass correlation coefficient of ",input$iccProp,", the needed sample size is ",resProp()*2, " patients (",resProp()," patients in each group)")
      # PROPORTION SUPERIORITY NO INTERMEDIATE ANALYSIS AND STEP WEDGE RANDOMIZATION --------------------
    }else if(reactive(input$hypProp)()=='Superiority' & reactive(input$propDesign)()=='No intermediate analysis with stepped wedge randomization'){
      if (resProp()<30) {
        warning("--> At least one group size is < 30, normality assumption is questionnable and sample
                size calculation may not be valid.")
      }else
      paste0("This sample size is for a stepped wedge cluster randomized superiority trial in two groups (experimental treatment versus control treatment) with balanced randomization (ratio 1:1) for a binary endpoint. The expected propotion of event is ",input$prevProp,"% with experimental treatment compared to ",input$obsProp,"% with control treatment. In order to highlight this absolute difference of ",abs(input$prevProp-input$obsProp),"%, with a ",printStProp()," type I error rate of ",input$alphaProp,"%, a power of ",input$powerProp,"%, according to our stepped wedge design with ",input$centersProp," centers randomized in ",input$sequencesProp," sequences and assuming an intraclass correlation coefficient of ",input$iccProp,", the needed sample size is ",resProp()*2, " patients (",resProp()," patients in each group)")
    # PROPORTION SUPERIORITY INTERMEDIATE ANALYSIS AND INDIVIDUAL RANDOMIZATION ------------------------------------------------------
    }else if(reactive(input$hypProp)()=='Superiority' & reactive(input$propDesign)()=='Intermediate analysis with individual randomization'){
      paste0("This sample size is for a RCT with an expected proportion of ",input$prevProp," units in patients in the experimental arm versus ",input$obsProp," units in the control arm. In order to demonstrate such adifference of ",abs(input$prevProp-input$obsProp)," units, with a standard deviation of ",input$sigmaProp,", a ",input$alphaProp,"% ",printStProp()," type I error rate and a power of ",input$powerProp,"%, the final analysis should be carried out on ",resProp()*2," patients (",resProp()," patients per group). The intermediate analyses would be performed on ",paste(ceiling(head(resPropSeq(),-2)), collapse = ", ")," and ",paste(tail(ceiling(head(resPropSeq(),-1)),1), collapse = ", ")," patients respectively if their is no decision of stopping the study.")
      
    # PROPORTION NON-INFERIORITY NO INTERMEDIATE ANALYSIS AND INDIVIDUAL RANDOMIZATION ----------------
    }else if(reactive(input$hypProp)()=='Non-inferiority' & reactive(input$propDesign)()=='No intermediate analysis with individual randomization'){
      if (resProp()<30) {
        warning("--> At least one group size is < 30, normality assumption is questionnable and sample 
                size calculation may not be valid.")
      }else
      paste0("This sample size is for a randomized controlled non-inferiority trial in two parallel groups (experimental treatment versus control treatment) with balanced randomization (ratio 1:1) for a binary endpoint. The expected percentage of events is ",input$obsProp,"% in patients in the control arm and no difference compared to the experimental arm. Assuming an absolute non-inferiority margin of ",input$deltaProp,"%, the needed sample size is ",resProp()*2, " patients (",resProp()," patients in each group) to achieve a ",input$alphaProp,"% one-sided type I error rate and a power of ",input$powerProp,"%.")
    # PROPORTION NON-INFERIORITY NO INTERMEDIATE ANALYSIS AND PARALLEL CLUSTER RANDOMIZATION ---------
    }else if(reactive(input$hypProp)()=='Non-inferiority' & reactive(input$propDesign)()=='No intermediate analysis with parallel cluster randomization'){
      if (resProp()<30) {
        warning("--> At least one group size is < 30, normality assumption is questionnable and sample 
                size calculation may not be valid.")
      }else
        paste0("This sample size is for a parallel cluster randomized non-inferiority trial in two parallel groups (experimental treatment versus control treatment) with balanced randomization (ratio 1:1) for a binary endpoint. The expected propotion of event is ",input$prevProp," % in patients in the control arm and no difference compared to the experimental arm. Assuming an absolute non-inferiority margin of ",input$deltaProp,"%, a one-sided type I error rate of ",input$alphaProp,"%, a power of ",input$powerProp,"%, according to our parallel cluster design with ",input$centersProp," centers randomized in two groups and assuming an intraclass correlation coefficient of ",input$iccProp,", the needed sample size is ",resProp()*2, " patients (",resProp()," patients in each group)")
    # PROPORTION NON-INFERIORITY NO INTERMEDIATE ANALYSIS AND STEP WEDGE RANDOMIZATION ---------
    }else if(reactive(input$hypProp)()=='Non-inferiority' & reactive(input$propDesign)()=='No intermediate analysis with stepped wedge randomization'){
      if (resProp()<30) {
        warning("--> At least one group size is < 30, normality assumption is questionnable and sample 
                size calculation may not be valid.")
      }else
      paste0("This sample size is for a stepped wedge cluster randomized non-inferiority trial in two groups (experimental treatment versus control treatment) with balanced randomization (ratio 1:1) for a binary endpoint. The expected propotion of event is ",input$prevProp," % in patients in the control arm and no difference compared to the experimental arm. Assuming an absolute non-inferiority margin of ",input$deltaProp,"%, a one-sided type I error rate of ", input$alphaProp,"%, a power of ",input$powerProp,"%, according to our stepped wedge design with ",input$centersProp," centers randomized in ",input$sequencesProp," sequences and assuming an intraclass correlation coefficient of ",input$iccProp,", the needed sample size is ",resProp()*2, " patients (",resProp()," patients in each group)")
    # PROPORTION NON-INFERIORITY INTERMEDIATE ANALYSIS AND INDIVIDUAL RANDOMIZATION -----------------------------------------------------------------------------  
    }else if(reactive(input$hypProp)()=='Non-inferiority' & reactive(input$propDesign)()=='Intermediate analysis with individual randomization'){
      paste0("This sample size is for a randomized controlled non-inferiority trial in two parallel groups (experimental treatment versus control treatment) with balanced randomization (ratio 1:1) for a binary endpoint. The expected percentage of events is ",input$obsProp,"% in patients in the control arm and no difference compared to the experimental arm. Assuming an absolute non-inferiority margin of ",input$deltaProp,"%, with a one-sided alpha risk of ",input$alphaProp,"% and a power of ",input$powerProp,"%, the final analysis should be carried out on ",resProp()*2," patients (",resProp()," patients per group). The intermediate analyses would be performed on ",paste(ceiling(head(resPropSeq(),-2)), collapse = ", ")," and ",paste(tail(ceiling(head(resPropSeq(),-1)),1), collapse = ", ")," patients respectively if their is no decision of stopping the study.")
    })
#---- 

  c <- reactive(
    if(reactive(input$typePredBEP)()=='Construction of predictive model of binary event'){
      paste0("This sample size is for developing a logistic regression model based on up to ",input$predictorsToTest," candidate predictors, with an anticipated R2 of at least ",input$R2,", and to target an expected shrinkage of ",input$shrinkageExpected," (equation 11 in Riley et al. Statistics in Medicine. 2019;38:1276–1296).")
    }else{
      paste0("This sample size is for external validation of a logistic regression model based with an expected outcome event proportions of ",input$P0,"%, with a alpha risk at ",input$alphaExtVal,"% and with a target confidence interval width of ",input$widthExtVal,"% (Riley et al.  Statistics in Medicine. 2021;19:4230-4251).")
    })

  d <- reactive(
    paste0("This sample size is for comparing accuracy of two diagnostic tests assuming a detection of ",(input$AUC2-input$AUC1)*100,"% difference in estimating two independent diagnostic systems (AUC1 = ",input$AUC1," and AUC2 = ",input$AUC2,") with ",(100-input$alphaAUC),"% confidence and ",input$powerAUC,"% power.")
  )  

#-----  
  e <- reactive(
    # SURVIVAL SUPERIORITY NO INTERMEDIATE ANALYSIS AND INDIVIDUAL RANDOMIZATION ----------------
    if (reactive(input$survDesign)()=='No intermediate analysis with individual randomization') {
      if (resSurv()<30) {
        warning("--> At least one group size is < 30, normality assumption is questionnable and sample
                size calculation may not be valid.")
      } else {
        paste0("This sample size is for a randomized controlled superiority trial in two parallel groups (experimental treatment versus control treatment) with balanced randomization (ratio 1:1) for a survival endpoint. The survival proportion of patient is ", input$prevSurv,"% with the experimental treatment compared to ", input$obsSurv,"% with the control treatment. In order to highlight this difference in survival curves, with a two-sided type I error rate of ", input$alphaSurv,"% and a power of ", input$powerSurv,"%, the needed sample size is ",resSurv()*2, " patients (",resSurv()," patients in each group).")
      }
    }
    # SURVIVAL SUPERIORITY NO INTERMEDIATE ANALYSIS AND STEPPED WEDGE RANDOMIZATION --------------
    else if (reactive(input$survDesign)()=='No intermediate analysis with stepped wedge randomization') {
      if (resSurv()<30) {
        warning("--> At least one group size is < 30, normality assumption is questionnable and sample
                size calculation may not be valid.")
      } else {
        paste0("This sample size is for a stepped wedge cluster randomized superiority trial in two groups (experimental treatment versus control treatment) with balanced randomization (ratio 1:1) for a survival endpoint. The survival proportion of patient is ", input$prevSurv,"% with the experimental treatment compared to ", input$obsSurv,"% with the control treatment. In order to highlight this difference in survival curves, with a two-sided type I error rate of ", input$alphaSurv,"%, a power of ", input$powerSurv,"%, according to our stepped wedge design with ", input$centersSurv," centers randomized in ", input$sequencesSurv," sequences and assuming an intraclass correlation coefficient of ", input$iccSurv,", the needed sample size is ",resSurv()*2, " patients (",resSurv()," patients in each group).")
      }
    }
  )
  
  
  # display the result ----------------------------------------------------------------------------------------------------------
  output$description <- renderText(paste("<p>The needed sample size is <b>", resDesc(), "</b> patients. </p>", z()))
  output$mean <- renderText(paste("<p>The needed sample size is <b>", resMean(), "</b> patients in each group. </p>", a()))
  output$proportion <- renderText(paste0("<p>The needed sample size is <b>", resProp(), "</b> patients in each group. </p>", b()))
  output$prediction <- renderText(paste0("<p>The needed sample size is <b>", resPred(), "</b> patients. </p>", c()))
  output$AUC <- renderText(paste0("<p>The needed sample size is <b>", resAUC(),"</b> patients.</p>", d()))
  output$survival <- renderText(paste0("<p>The needed sample size is <b>", resSurv(), "</b> patients in each group.</p>", e()))
  
  # copy the result -------------------------------------------------------------------------------------------------------------
  observeEvent(input$copy_link_description, {
    text <-  paste("The needed sample size is", resDesc(), "patients.", z())
    session$sendCustomMessage("txt", text)
  })
  observeEvent(input$copy_link_mean, {
    text <-  paste("The needed sample size is", resMean(), "patients in each group.", a())
    session$sendCustomMessage("txt", text)
  })
  observeEvent(input$copy_link_proportion, {
    text <-  paste("The needed sample size is", resProp(), "patients in each group.", b())
    session$sendCustomMessage("txt", text)
  })
  observeEvent(input$copy_link_prediction, {
    text <-  paste("The needed sample size is", resPred(), "patients.", c())
    session$sendCustomMessage("txt", text)
  })
  observeEvent(input$copy_link_AUC, {
    text <-  paste("The needed sample size is", resAUC(), "patients.", d())
    session$sendCustomMessage("txt", text)
  })
  observeEvent(input$copy_link_survival, {
    text <-  paste("The needed sample size is", resSurv(), "patients in each group.", e())
    session$sendCustomMessage("txt", text)
  })
}

# Run the application
shinyApp(ui = ui, server = server)
