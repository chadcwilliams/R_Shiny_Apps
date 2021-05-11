options(scipen = 999)
library(rsconnect)
library(BSDA)
library(rhandsontable) #Data Tables

#UI
ui = fluidPage(tags$head(tags$style(type = "text/css", ".irs {max-width: 946px;}")),
               sidebarLayout(
                   sidebarPanel(
                       selectInput(
                           "Test",
                           label = " ",
                           choices = list(
                               "Frequency Distribution" = 1,
                               "Descriptives" = 2,
                               "Single Sample Z-Test" = 3,
                               "Single Sample T-Test" = 4,
                               "Paired Sample T-Test" = 4
                           ),
                           selected = 1
                       ),
                       sliderInput(
                           inputId = 'num_of_participants',
                           label = 'Number of Participants',
                           value = 10,
                           min = 2,
                           max = 20,
                           step = 1,
                           width = '95%'
                       ),
                       actionButton('refresh', 'Refresh'),
                       actionButton('answers', 'Show Answers')
                   ),
                   mainPanel(
                       rHandsontableOutput("data_display"),
                       rHandsontableOutput("stats_display")
                   )
               ))

#Server
server = function(input, output) {
    stats = reactiveValues(data_table = NULL)
    observeEvent(input$refresh, {
        #Frequency Distribution
        if (input$Test == 1) {
            #Create Data
            data = sample(1:10, input$num_of_participants, replace = TRUE)
            #Create Stats
            freq_dist = table(data)
            rel_freq = prop.table(freq_dist)
            
            frequency_distribution = data.frame(
                Data = min(data):max(data),
                Frequency = NA,
                Relative_Frequency = NA,
                Cumulative_Frequency = NA,
                Cum_Rel_Freq = NA
            )
            
            frequency_dist = data.frame(
                Data = rownames(freq_dist),
                Frequency = as.integer(freq_dist),
                Relative_Frequency = as.numeric(rel_freq),
                Cumulative_Frequency = cumsum(freq_dist),
                Cum_Rel_Freq = cumsum(rel_freq)
            )
            for (counter in 1:dim(freq_dist)[1]){
                frequency_distribution[frequency_distribution$Data==frequency_dist$Data[counter],] = frequency_dist[counter,]
            }
            
            missing = which(is.na(frequency_distribution$Frequency))
            frequency_distribution$Frequency[missing]=0
            frequency_distribution$Relative_Frequency[missing] = 0.00
            for (counter in 1:dim(as.data.frame(missing))[1]){
                frequency_distribution$Cumulative_Frequency[missing[counter]] = frequency_distribution$Cumulative_Frequency[missing[counter]-1]
                frequency_distribution$Cum_Rel_Freq[missing[counter]] = frequency_distribution$Cum_Rel_Freq[missing[counter]-1]}
            frequency_distribution$Frequency = as.integer(frequency_distribution$Frequency)
            #Set outputs
            stats$data_table = frequency_distribution[order(nrow(frequency_distribution):1),]
            output$data_display = renderRHandsontable(rhandsontable(as.data.frame(data)))
            output$stats_display = renderRHandsontable({})
        }
        else if (input$Test == 2) {
            #Descriptives
            #Create Data
            data = sample(1:100, input$num_of_participants, replace = TRUE)
            #Create Stats
            descriptives = data_table = data.frame(
                Mean = mean(data),
                Variance = var(data),
                Standard_Deviation = sd(data)
            )
            #Set outputs
            stats$data_table = descriptives
            output$data_display = renderRHandsontable(rhandsontable(as.data.frame(data)))
            output$stats_display = renderRHandsontable({
            })
        }
        else if (input$Test == 3) {
            #Single Sample Z-Test
            #Create Data
            data = data.frame(
                Data = sample(1:100, input$num_of_participants, replace = TRUE),
                mu = c(
                    sample(1:100, 1),
                    rep(NA, input$num_of_participants - 1)
                ),
                sigma = c(
                    sample(1:10, 1),
                    rep(NA, input$num_of_participants - 1)
                )
            )
            data$mu[1] = as.integer(round(rnorm(1,mean(data$Data),sd(data$Data)/5)))
            t = z.test(data$Data, mu = data$mu[1], sigma.x = data$sigma[1])
            #Create Table
            descriptives = data_table = data.frame(
                Mean = mean(data$Data),
                Variance = var(data$Data),
                Standard_Deviation = sd(data$Data),
                Z_Value = round(as.numeric(t["statistic"]), digits = 2),
                Degrees_of_Freedom = as.numeric(input$num_of_participants - 1),
                P_Value = round(as.numeric(t["p.value"]), digits = 4)
            )
            #Set Outputs
            stats$data_table = descriptives
            output$data_display = renderRHandsontable(rhandsontable(as.data.frame(data)))
            output$stats_display = renderRHandsontable({
            })
        }
        else if (input$Test == 4) {
            #Single Sample T-Test
            #Create Data
            data = data.frame(
                Data = sample(1:100, input$num_of_participants, replace = TRUE),
                mu = c(
                    sample(1:100, 1),
                    rep(NA, input$num_of_participants - 1)
                )
            )
            data$mu[1] = as.integer(round(rnorm(1,mean(data$Data),sd(data$Data)/5)))
            #Create Stats
            t = t.test(data$Data, mu = data$mu[1])
            descriptives = data_table = data.frame(
                Mean = mean(data$Data),
                Variance = var(data$Data),
                Standard_Deviation = sd(data$Data),
                T_Value = round(as.numeric(t["statistic"]), digits = 2),
                Degrees_of_Freedom = as.numeric(input$num_of_participants - 1),
                P_Value = round(as.numeric(t["p.value"]), digits = 4)
            )
            #Set Outputs
            stats$data_table = descriptives
            output$data_display = renderRHandsontable(rhandsontable(as.data.frame(data)))
            output$stats_display = renderRHandsontable({
            })
        } else {
            #Repeated Measures T-Test
            run = 1
            while (run == 1) {
                data1 = sample(1:100, input$num_of_participants, replace = TRUE)
                data2 = sample(1:100, input$num_of_participants, replace = TRUE)
                M1 = mean(data1)
                V1 = var(data1)
                SD1 = sd(data1)
                M2 = mean(data2)
                V2 = var(data2)
                SD2 = sd(data2)
                t = t.test(data1, data2, paired = TRUE)
                if (input$decimal == 2) {
                    if (M1 %% 1 == 0 &
                        V1 %% 1 == 0 &
                        SD1 %% 1 == 0 &
                        M2 %% 1 == 0 &
                        V2 %% 1 == 0 &
                        SD2 %% 1 == 0) {
                        print(paste("Group 1 Data: ", toString(data1)))
                        print(paste("Group 2 Data: ", toString(data2)))
                        output$stats_display = renderRHandsontable({
                            
                        })
                        run = 0
                    }
                } else{
                    print(paste("Group 1 Data: ", toString(data1)))
                    print(paste("Group 2 Data: ", toString(data2)))
                    output$stats_display = renderRHandsontable({
                        
                    })
                    run = 0
                }
            }
        }
    })
    observeEvent(input$answers,
                 {
                     output$stats_display = renderRHandsontable(rhandsontable(stats$data_table))
                 })
}

shinyApp(ui = ui, server = server)


#Deploy App
##Published Version
#rsconnect::setAccountInfo(name='chadcwilliams', token='64ECB37313FD540CACEEC1410973047A', secret='dVExTcA9wF9YYBPLTe63lgVrNXj3EW/FYObYMAD6')

##Example Version
#rsconnect::setAccountInfo(name='chadcwilliams', token='DB339619E53122BF203F3772D0BE27B0', secret='Mls9Q8LjxNaxV+l+pro1Wc7f/uQ2LI8BiD3I+lwT')

#Deploy App
#rsconnect::deployApp(appName = "StatsProblemGenerator",' /Users/Chad/Documents/GitHub/Random R Script/StatCalcExamples/')
