#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(scales)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Subsampling TCR Data"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        selectInput("sex",
                     "Sex (M/F)",
                     choices = c('Male', 'Female')
                    ),
        
        numericInput("age",
                     "Age (Years):",
                     min = 0.0,
                     max = 100,
                     value = 2),
        
        
        numericInput("height",
                     "Height (M):",
                     min = 0.0,
                     max = 3.0,
                     value = 1.2),
         
         numericInput("weight",
                      "Weight (kg):",
                      min = 0.0,
                      max = 200,
                      value = 40),
         
         numericInput("CD3s",
                      "CD3 Count * 10^9/L:",
                      min = 0.0,
                      max = 100,
                      value = 1.56),
         
         numericInput("inputVolume",
                      "Input Volume For Protocol (mLs):",
                      min = 0.0,
                      max = 20,
                      value = 10),
         
         numericInput("usedPcnt",
                      "% of Input Volume Used For Sequencing",
                      min = 0.0,
                      max = 100,
                      value = 50),
         
         numericInput("sequencingPcnt",
                      "% of Input Reads Making It Through Sequencing",
                      min = 0.0,
                      max = 100,
                      value = 50),
         
         numericInput("decombinatorPcntIdentified",
                      "% of Sequences Identified By Decombinator",
                      min = 0.0,
                      max = 100,
                      value = 50),
         
         numericInput("nonproductivesPcnt",
                      "% of Sequences Identified As Nonproductive By Decombinator",
                      min = 0.0,
                      max = 100,
                      value = 7),
         
         numericInput("subsampledPcnt",
                      "Subsampling Depth (% of Sequenced Reads)",
                      min = 0.0,
                      max = 100,
                      value = 10),
         
         actionButton("drawPlot",
                      "Draw Plot")
      ),
   
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("distPlot")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  observeEvent( input$drawPlot, {
    output$distPlot <- renderPlot({
      # generate bins based on input$bins from ui.R
      #numSeqs in a child - Biology
      
      #colnames(dat) <- c('seqs', 'pcnt')
      
      CD3s <- input$CD3s * 10^6 #corrects for mLs
      
      sex <- input$sex
      
      age <- round(input$age)
      
      #Old method
      #TBV <- blood_vol_Nadler(input$height, input$weight, male = F)*1000
      
      weight <- input$weight
      
      height <- input$height
      
      if (age >= 0 & age <= 1){
        #M & F 0 - 1
        logTBV <- (0.7891*log10(weight)) + (0.004132*(height*100)) + 1.8117
        TBV <- 10 ^ logTBV
      }
      
      if (sex == 'Male' & age >= 2 & age <= 14 | sex == 'Female' & age >= 2 & age <= 6){
        #F 2 - 6, M 2- 14
        logTBV <- (0.6459*log10(weight)) + (0.002743*(height*100)) + 2.0324
        TBV <- 10 ^ logTBV
      }
      
      if (sex == 'Female' & age >= 7 & age <= 14){
        #F 7 - 14
        logTBV <- (0.6412*log10(weight)) + (0.001270*(height*100)) + 2.2169
        TBV <- 10 ^ logTBV
      }
      
      if (sex == 'Female' & age >= 15 & age <= 100){
        #F 15 - 100
        TBV <- ((35.5*(height*100)) + (2.27 * weight) - 3382)/0.6178
      }
      
      if (sex == 'Male' & age >= 15 & age <= 100){
        #M 14 - 100
        TBV <- ((13.1*(height*100)) + (18.05 * weight) - 480)/0.5723
      }
      
      
      #TBV <- 10 ^ logTBV
      
      
      numSeqs <- CD3s*TBV #Assuming 1 sequence per circulating CD3 cell
      
      inputFrac <- input$inputVolume/TBV
      
      usedFrac <- input$usedPcnt/100
      
      numSampledSeqs <- numSeqs*inputFrac*usedFrac
      numSampledSeqsPcnt <- numSampledSeqs/numSeqs*100
      sequencingFrac <- input$sequencingPcnt/100
      #print(c(numSampledSeqs, sequencingFrac))
      numSequendedSeqs <- numSampledSeqs * sequencingFrac
      numSequencedSeqsPcnt <- numSequendedSeqs/numSeqs*100
      
      decombinatorFrac <- input$decombinatorPcntIdentified/100
      
      nonproductivesFrac <- 1 - (input$nonproductivesPcnt/100)
      
      subsampledFrac <- input$subsampledPcnt/100
      
      numDecombinedSeqs <- numSequendedSeqs * decombinatorFrac * nonproductivesFrac * subsampledFrac 
      numDecombinedSeqsPcnt <- numDecombinedSeqs/numSeqs*100
      
      dat <- data.frame(c(numSeqs, numSampledSeqs, numSequendedSeqs, numDecombinedSeqs), c(100, numSampledSeqsPcnt, numSequencedSeqsPcnt, numDecombinedSeqsPcnt), c('Full TCR', 'Blood Sample', 'Sequenced', 'Pipeline'))
      colnames(dat) <- c('seqs', 'pcnt', 'lab')
      xbreaks <- dat$pcnt
      #print(xbreaks)
      xbreaks[xbreaks > 1] <- round(xbreaks[xbreaks > 1])
      breaks.majorx <- c(0.001,0.01,0.1,1,10,100)
      
      breaks.combx <- round(sort(c(breaks.majorx, xbreaks)),digits = 3)
      
      ybreaks <- dat$seqs
      beaks.majory <- c(1,10,100,1000,10000,100000,1000000,10000000,100000000,1000000000,10000000000)
      
      breaks.comby <- round(sort(c(beaks.majory, ybreaks)),digits = 3)
      
      # draw the histogram with the specified number of bins
      ggp <- ggplot(dat)+
        geom_line(aes(x = pcnt, y = seqs))+
        geom_text(aes(x = pcnt, y = seqs+(seqs*1.5), label = lab, colour = as.factor(lab)))+
        geom_point(aes(x = pcnt, y = seqs, colour = as.factor(lab)), size = 2)+
        geom_segment(aes(yend = seqs+(seqs*1), y=min(dat$seqs), x = pcnt,  xend = pcnt), linetype = 2)+
        geom_segment(aes(y=seqs, yend = seqs, x=min(dat$pcnt), xend = pcnt ), linetype = 2)+
        scale_x_log10(minor_breaks = xbreaks, breaks = breaks.combx)+
        scale_y_log10(minor_breaks = ybreaks, breaks = breaks.comby, labels = comma)+
        #expand_limits(x=500)+
        labs(y = 'Number of Sequences', x = '% of Estimated Total TCR')+
        guides(colour = F)+
        theme_classic()
      return(ggp)
    })
  })
}
# Run the application 
shinyApp(ui = ui, server = server)

