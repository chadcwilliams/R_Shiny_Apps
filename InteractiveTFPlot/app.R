##########################################################################
##Written by Chad C. Williams, PhD student in the Krigolson Lab, 2020   ##
##University of Victoria, British Columbia, Canada                      ##
##www.krigolsonlab.com                                                  ##
##www.chadcwilliams.com                                                 ##
##########################################################################

#### (1) Setup ####
library(ggplot2)
library(gganimate)
library(RColorBrewer)
library(reshape2)
library(jpeg)
library(grid)
library(gridExtra)
library(magick)
library(transformr)
library(ggplotify)
library(dplyr)
library(shiny)
#library(eegUtils) #This had some trouble, so I extracted the following functions from it:
source('ggplot2-extensions.R')
source('import_chans.R')
source('circ_rad_fun.R')

#### (2) Load and prepare data ####

#####  Load and Format Data
frequency_data_temp = read.csv('Frequency_Data.csv',header=FALSE) #Load data
frequency_data = rbind(rep(0,602),frequency_data_temp[1:30,],rep(0,602),frequency_data_temp[31:60,]) #Create data structure (filtered and unfiltered) with time variable 
colnames(frequency_data) = seq(-200,1000,by=2) #Change column names
frequency_data$Freq = as.factor(0:30) #Convert to factor
frequency_data$Filter = sort(rep(1:2,31)) #Create filter ID
data = melt(frequency_data,id=c('Freq','Filter')) #Convert to long format
colnames(data) = c('Freq','Filter','Time','Value') #Change column names
time_index = seq(-200,1000,by=2) #Determine time IDs
data$MeanValue = 0 #Setup column 
for (f in 1:30){ #Cycle through frequencies
  freq_temp = data[data$Freq == as.character(f),] #Determine rows of current frequency
  for (c in 1:601){ #Cycle through time
    if (c < 25){ #Determine if time is near beginning
      timing = seq(1,c) #Average up to this time point
    }else{ #Else the time is not near beginning
      timing = seq(c-24,c+24) #Average around this time point
    }
    data$MeanValue[data$Freq == as.character(f) & data$Time == as.character(time_index[c])] = mean(freq_temp$Value[timing],na.rm = TRUE) #Create grand averaged time-frequency data
  }
}

data$Freq = (as.numeric(as.character(data$Freq))) #Determine frequencies
data$Time = (as.numeric(as.character(data$Time))) #Determine times
colour_index = seq(-.6,.6,by=.1) #Determine colour scale

#Create Animated 2D Wavelet Plot
temp_WAV_data = melt(frequency_data[,1:603],id=c('Freq','Filter')) #Create data frame
colnames(temp_WAV_data) = c('Freq','Filter','Time','Value') #Rename data frame

#Correlational Data
corr_data = read.csv('BaseRate3_Corr_Data_Shiny.csv',header = FALSE) #Load data
colnames(corr_data) = c('Participant','Accuracy','Theta','Time','Frequency') #Name columns

#Topo Data
topo_data = read.csv('BaseRate3_Topo_Data_Shiny.csv',header = FALSE) #Load data
colnames(topo_data) = c('Electrode','Values','Time','Frequency') #Name columns
eeg_base = import_chans("BaseRate3_chanlocs2.txt") #Load channel locations
eeg_example = do.call("rbind", replicate(18030, eeg_base, simplify = FALSE)) #Create topographic dataframe
eeg_example$values = topo_data$Values #Re-assign column
eeg_example$Time = topo_data$Time #Re-assign column
eeg_example$Frequency = topo_data$Frequency #Re-assign column

#### Shiny ####
#### UI ####
ui = fluidPage(tags$head(tags$style(type="text/css", ".irs {max-width: 946px;}")), #Setup UI
 selectInput("filterSelect", label = " ", choices = list("Filtered" = 1, "Unfiltered" = 2), selected = 1), #Create dropdown interactive
 sliderInput(inputId = 'timeNum',label = 'Start and End Time',value = c(100,600),min = -200,max = 1000, step = 1, width ='95%'), #Create sliding bar interactive
 sliderInput(inputId = 'freqNum',label = 'Start and End Frequency',value = c(4,8),min = 1,max = 30, step = 1, width ='95%'), #Create sliding bar interactive
 plotOutput(outputId = 'allPlots', height = 650, width = 946), #Insert plots with specified dimensions
 mainPanel(p("\nSupplemental Figure 1. Interactive plot...")) #Add figure caption
 )

#### Server ####
server = function(input,output){ #Setup main code
  #Convert Variables
  startTime = reactive({((((input$timeNum[1])+200)/2)+1)}) #Determine start time
  endTime = reactive({((((input$timeNum[2])+200)/2)+1)}) #Determine end time
  startFreq = reactive({input$freqNum[1]+1}) #Determine start frequenct
  endFreq = reactive({input$freqNum[2]+1}) #Determine end frequency
  
  #Plot Wavelet Data
  output$allPlots = renderPlot({ #Setup plot
      #Create Wavelet Plot with Analysis Box
      WAV_data = temp_WAV_data[(temp_WAV_data$Filter==input$filterSelect),] #Extract data 
      wavPlot = ggplot(WAV_data, aes(as.numeric(Time), as.numeric(Freq), z=Value))+ #Setup plot
      geom_contour_filled(bins = 9)+ #Construct time-frequency plot
      theme_minimal()+theme(legend.position = "none",text = element_text(size=16))+ #Specify theme characteristics
      scale_fill_brewer(type = 'div',palette='RdBu',direction = -1)+ #Fill colours
      scale_x_continuous(breaks=c(1,seq(50,600,by=50)),labels=seq(-200,1000,by=100),expand=c(0,0))+xlab('Time (ms)')+ #Setup x axis
      scale_y_continuous(breaks=seq(3,31,2),labels=seq(2,30,2),expand = c(0,0))+ylab('Frequency (Hz)')+ #Setup y axis
      theme(plot.margin = margin(1,1,1,1, "cm"))+ #Determine margin size
      coord_cartesian(xlim =c(0,600),ylim = c(2,31))+ylab('Frequency (Hz)')+ #Setup limits and labels for axes
      geom_segment(aes(x=startTime(),xend=startTime(),y=startFreq(),yend=endFreq()),alpha=.5,size=.5)+ #Setup interactive box
      geom_segment(aes(x=endTime(),xend=endTime(),y=startFreq(),yend=endFreq()),alpha=.5,size=.5)+ #Setup interactive box
      geom_segment(aes(x=startTime(),xend=endTime(),y=startFreq(),yend=startFreq()),alpha=.5,size=.5)+ #Setup interactive box
      geom_segment(aes(x=startTime(),xend=endTime(),y=endFreq(),yend=endFreq()),alpha=.5,size=.5) #Setup interactive box

      #Determine Circle Bar Plot Data
      timeConstrainedData = data[(as.numeric(data$Time) >= input$timeNum[1]) & (as.numeric(data$Time) <= input$timeNum[2]) & (as.numeric(data$Filter) == input$filterSelect),] #Setup data structure
      circlePlotData = data.frame(Freq = 0,Time = mean(c(input$timeNum[1],input$timeNum[2])),Value = 0) #Determine time of interest
      for (currentFrequency in 1:30){ #Cycle through frequencies
        circlePlotData[currentFrequency+1,] = c(currentFrequency,mean(c(input$timeNum[1],input$timeNum[2])),mean(timeConstrainedData[timeConstrainedData$Freq == currentFrequency,4])) #Determined averaged data
      }
      scale = max(abs(c(min(circlePlotData$Value),max(circlePlotData$Value)))) #Determine max and min values of bars 
      #Adjust the scale of the plot dependent on current data
      if(scale <= .355){
        scale = .4
      } else if (scale <= .55) {
        scale = .6
      } else if (scale <= .75) {
        scale = .8
      } else {
        scale = 1.0
      }
      
      #Create Circle Bar Plot
      fsize = 4.5 #Determine font size
      circlePlot = ggplot(circlePlotData, aes(x=Freq, y=Value, fill=Value))+ #Setup plot
      geom_hline(yintercept=colour_index[seq(3,13,by=2)], linetype="dashed",alpha = .075)+ #Create y axis guides
      geom_bar(stat="identity")+ #Add bars to plot
      theme_minimal()+ #Change theme
      theme(plot.margin = margin(1,1,1,1, "cm"))+ #Setup margin size
      theme(axis.text = element_blank(),axis.title = element_blank(),panel.grid = element_blank(),legend.position = "none",plot.title = element_text(hjust = 0.5,vjust = -75),plot.margin = unit(rep(-1,4), "cm"))+ #Determine theme characteristics
      coord_polar(start = 0)+ #Determine orientation of plot
      scale_fill_gradient2(low='#313695',mid='#ffffff',high='#a50026')+ #Determine colour of bars
      ylim(-scale,scale-.01)+ #Determine y axis
      annotate("text", x = 31, y = c(colour_index[seq(3,13,by=2)],scale-.01), label = c(round(colour_index[seq(3,13,by=2)],1),'Power (dB)'),size = 5.5)+ #Add y axis label
      annotate("text", x = 1:31, y = scale-.1, label = c(1:30,NA),size = fsize)+ #Add frequency labels
      annotate("text", x = 15, y = scale-.01, label = 'Frequency (Hz)',size = 5.5) #Add x axis label
      
      #Determine Topographic Data
      constrainedTopoPlotData = eeg_example[(as.numeric(eeg_example$Time) >= input$timeNum[1]) & (as.numeric(eeg_example$Time) <= input$timeNum[2]) & (as.numeric(eeg_example$Frequency) >= input$freqNum[1]) &  (as.numeric(eeg_example$Frequency) <= input$freqNum[2]),] #Construct topographic data
      topoPlotData = eeg_example[1:29,] #Re-assign structure
      for (electrodeIndex in 1:29){ #Cycle through electrodes
        topoPlotData[electrodeIndex,10] = colMeans(as.data.frame(constrainedTopoPlotData[constrainedTopoPlotData$electrode==constrainedTopoPlotData$electrode[electrodeIndex],10])) #Average for each electrode
      }
      
      #Create Topographic Plot
      topoPlot = ggplot(topoPlotData,aes(x = x,y = y,fill = values,label = electrode))+ #Setup plot
        geom_topo(grid_res = 100,interp_limit = "head",chan_size = 2)+ #Construct plot
        scale_fill_distiller(palette = "RdBu")+ #Determine colours
        theme_void()+theme(legend.position = "none")+ #Change theme
        theme(plot.margin = margin(1.5,1.5,1.5,1.5, "cm"))+ #Change margin size
        coord_equal() #Determine coordinates
      
      #Create Correlational Data
      constrainedCorrData = corr_data[(as.numeric(corr_data$Time) >= input$timeNum[1]) & (as.numeric(corr_data$Time) <= input$timeNum[2]) & (as.numeric(corr_data$Frequency)>=input$freqNum[1]) &(as.numeric(corr_data$Frequency)<=input$freqNum[2]),] #Determine correlational structure
      corrData = constrainedCorrData[1:30,] #Re-assign structure
      for (participantIndex in 1:30){ #Cycle through participants
        corrData[participantIndex,3] = mean(constrainedCorrData[constrainedCorrData$Participant==participantIndex,3]) #Determine data
      }
      
      #Create Correlational Plot
      correlational_Plot = ggplot(aes(x=Accuracy,y=Theta),data=corrData)+ #Setup plot
        geom_point(colour = '#313695', alpha = .8)+ #Add scatterplot
        geom_smooth(method = 'lm',formula = 'y~x',se = FALSE, colour = '#a50026', alpha = .8)+ #Add regression line
        xlab('Percentage of Rational Responses')+ylab('Oscillatory Activity (dB)')+ylim(-4,4)+ #Add axis limits and labels
        theme_minimal()+theme(plot.margin = margin(1,1,1,1, "cm"))+ #Change theme
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size=16)) #Determine theme characteristics
  
      grid.arrange(wavPlot, correlational_Plot, circlePlot, topoPlot, ncol = 2) #Combine all plots 
  })
  
}

shinyApp(ui=ui,server=server) #Setup app to run
