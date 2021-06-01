####################################################################
####              Written by Chad C. Williams, 2021             ####
####                   www.chadcwilliams.com                    ####
####################################################################

options(scipen = 999) #Remove scientific notation
library(rsconnect) #Shiny
library(BSDA) #z-test function
library(rhandsontable) #Data tables
library(ggplot2) #Plotting

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
                               "Single Participant Z-Test" = 3,
                               "Correlation & Regression" = 4,
                               "Single Sample T-Test" = 5,
                               "Paired Sample T-Test" = 5
                           ),
                           selected = 1
                       ),
                       sliderInput(
                           inputId = 'num_of_participants',
                           label = 'Number of Participants',
                           value = 10,
                           min = 2,
                           max = 100,
                           step = 1,
                           width = '95%'
                       ),
                       sliderInput(
                           inputId = 'value_range',
                           label = 'Range of Values',
                           value = c(1, 10),
                           min = 1,
                           max = 100,
                           step = 1,
                           width = '95%'
                       ),
                       actionButton('refresh', 'Refresh'),
                       actionButton('distribution', 'Plot Data'),
                       actionButton('answers', 'Show Answers')
                   ),
                   mainPanel(fluidRow(
                       column(3, rHandsontableOutput("data_display")),
                       column(9, plotOutput('distribution_display'))
                   ),
                   fluidRow(column(
                       12, rHandsontableOutput("stats_display")
                   )))
               ))

#Server
server = function(input, output) {
    stats = reactiveValues(data_table = NULL)
    plotdata = reactiveValues(data = NULL)
    observeEvent(input$refresh, {
        #Frequency Distribution
        if (input$Test == 1) {
            #Create Data
            data = sample(
                input$value_range[1]:input$value_range[2],
                input$num_of_participants,
                replace = TRUE
            )
            plotdata$data = as.data.frame(data)
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
            for (counter in 1:dim(freq_dist)[1]) {
                frequency_distribution[frequency_distribution$Data == frequency_dist$Data[counter],] = frequency_dist[counter,]
            }
            
            missing = which(is.na(frequency_distribution$Frequency))
            frequency_distribution$Frequency[missing] = 0
            frequency_distribution$Relative_Frequency[missing] = 0.00
            for (counter in 1:dim(as.data.frame(missing))[1]) {
                frequency_distribution$Cumulative_Frequency[missing[counter]] = frequency_distribution$Cumulative_Frequency[missing[counter] -
                                                                                                                                1]
                frequency_distribution$Cum_Rel_Freq[missing[counter]] = frequency_distribution$Cum_Rel_Freq[missing[counter] -
                                                                                                                1]
            }
            frequency_distribution$Frequency = as.integer(frequency_distribution$Frequency)
            #Set outputs
            stats$data_table = frequency_distribution[order(nrow(frequency_distribution):1),]
            output$data_display = renderRHandsontable(rhandsontable(as.data.frame(data)))
            output$stats_display = renderRHandsontable({
                
            })
            output$distribution_display = renderPlot({
                
            })
        }
        else if (input$Test == 2) {
            #Descriptives
            #Create Data
            data = sample(
                input$value_range[1]:input$value_range[2],
                input$num_of_participants,
                replace = TRUE
            )
            plotdata$data = as.data.frame(data)
            #Force a Mode (By duplicating one of the numbers)
            if (length(unique(data)) == input$num_of_participants) {
                index = sample(1:input$num_of_participants, 2)
                data[index[1]] = data[index[2]]
            }
            #Setup Mode Function
            mod = function(data) {
                unique_x = unique(data)
                tabulate_x = tabulate(match(data, unique_x))
                unique_x[tabulate_x == max(tabulate_x)]
            }
            #Setup Semi-Interquartile Range Function
            siqr = function(data) {
                sorted_data = sort(data)
                if (length(sorted_data) %% 2 == 0) {
                    q1 = median(sorted_data[1:(length(sorted_data) / 2)])
                    q3 = median(sorted_data[((length(sorted_data) / 2) +
                                                 1):length(sorted_data)])
                } else{
                    q1 = median(sorted_data[1:((length(sorted_data) - 1) / 2)])
                    q3 = median(sorted_data[((length(sorted_data) - ((
                        length(sorted_data) - 1
                    ) / 2)) + 1):length(sorted_data)])
                }
                (q3 - q1) / 2
            }
            #Create Stats
            descriptives = data_table = data.frame(
                Mode = mod(data),
                Median = median(data),
                Mean = mean(data),
                Range = range(data)[2] - range(data)[1],
                SIQR = siqr(data),
                MAD = median(abs(data - median(data))),
                SS = sum((data - mean(data)) ^ 2),
                Var = sum((data - mean(data)) ^ 2) / input$num_of_participants,
                SD = sqrt(sum((
                    data - mean(data)
                ) ^ 2) / input$num_of_participants),
                SkewP = (3 * (mean(data) - median(data))) / sqrt(sum((
                    data - mean(data)
                ) ^ 2) / input$num_of_participants)
            )
            #Clear the Duplicate Values
            if (dim(descriptives)[1] > 1) {
                descriptives[2:dim(descriptives)[1], 2:dim(descriptives)[2]] = NA
            }
            #Set outputs
            stats$data_table = descriptives
            output$data_display = renderRHandsontable(rhandsontable(as.data.frame(data)))
            output$stats_display = renderRHandsontable({
                
            })
            output$distribution_display = renderPlot({
                
            })
        }
        else if (input$Test == 3) {
            #Single Participant Z-Test
            #Create Data
            data = data.frame(
                X = sample(
                    seq(input$value_range[1], input$value_range[2], by = .1),
                    1
                ),
                mu = sample(
                    seq(input$value_range[1], input$value_range[2], by = .1),
                    1
                ),
                sigma = rnorm(1, (
                    input$value_range[2] - input$value_range[1]
                ) / 5, .1)
            )
            data2 = data.frame(data = dnorm(
                seq((data$mu - (
                    4 * data$sigma
                )), (data$mu + (
                    4 * data$sigma
                )), length.out = 100),
                mean = data$mu,
                sd = data$sigma
            ))
            plotdata$data = data2
            
            #Create Table
            descriptives = data.frame(
                Z_Value = (data$X - data$mu) / data$sigma,
                P_Value_of_X_and_Below = round(pnorm(
                    round((data$X - data$mu) / data$sigma, digits = 2)
                ), digits = 4),
                P_Value_of_X_and_Above = round(pnorm(
                    round((data$X - data$mu) / data$sigma, digits = 2), lower.tail = F
                ), digits = 4)
            )
            #Set Outputs
            stats$data_table = descriptives
            output$data_display = renderRHandsontable(rhandsontable(as.data.frame(t(data))))
            output$stats_display = renderRHandsontable({
                
            })
            output$distribution_display = renderPlot({
                
            })
        }
        else if (input$Test == 4) {
            #Single Sample Z-Test
            #Create Data
            data = data.frame(
                X = sample(
                    seq(input$value_range[1], input$value_range[2], by = .1),
                    1
                ),
                mu = sample(
                    seq(input$value_range[1], input$value_range[2], by = .1),
                    1
                ),
                sigma = rnorm(1, (
                    input$value_range[2] - input$value_range[1]
                ) / 5, .1)
            )
            data2 = data.frame(data = dnorm(
                seq((data$mu - (
                    4 * data$sigma
                )), (data$mu + (
                    4 * data$sigma
                )), length.out = 100),
                mean = data$mu,
                sd = data$sigma
            ))
            plotdata$data = data2
            
            #Create Table
            descriptives = data.frame(
                Z_Value = (data$X - data$mu) / data$sigma,
                P_Value_of_X_and_Below = round(pnorm(
                    round((data$X - data$mu) / data$sigma, digits = 2)
                ), digits = 4),
                P_Value_of_X_and_Above = round(pnorm(
                    round((data$X - data$mu) / data$sigma, digits = 2), lower.tail = F
                ), digits = 4)
            )
            #Set Outputs
            stats$data_table = descriptives
            output$data_display = renderRHandsontable(rhandsontable(as.data.frame(t(data))))
            output$stats_display = renderRHandsontable({
                
            })
            output$distribution_display = renderPlot({
                
            })
        }
        else if (input$Test == 5) {
            #Single Sample T-Test
            #Create Data
            data = data.frame(
                Data = sample(
                    input$value_range[1]:input$value_range[2],
                    input$num_of_participants,
                    replace = TRUE
                ),
                mu = c(
                    sample(1:20, 1),
                    rep(NA, input$num_of_participants - 1)
                )
            )
            data$mu[1] = as.integer(round(rnorm(
                1, mean(data$Data), sd(data$Data) / 5
            )))
            plotdata$data = as.data.frame(data$Data)
            #Create Stats
            t = t.test(data$Data, mu = data$mu[1])
            descriptives = data_table = data.frame(
                Mean = mean(data$Data),
                Variance = var(data$Data),
                SD = sd(data$Data),
                T_Value = round(as.numeric(t["statistic"]), digits = 2),
                Degrees_of_Freedom = as.numeric(input$num_of_participants - 1),
                P_Value = round(as.numeric(t["p.value"]), digits = 4)
            )
            #Set Outputs
            stats$data_table = descriptives
            output$data_display = renderRHandsontable(rhandsontable(as.data.frame(data)))
            output$stats_display = renderRHandsontable({
                
            })
            output$distribution_display = renderPlot({
                
            })
        } else {
            #Repeated Measures T-Test
            run = 1
            while (run == 1) {
                data1 = sample(1:20, input$num_of_participants, replace = TRUE)
                data2 = sample(1:20, input$num_of_participants, replace = TRUE)
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
                     output$stats_display = renderRHandsontable(if (input$Test == 3) {
                         rhandsontable(stats$data_table) %>%
                             hot_col("P_Value_of_X_and_Below", format = "0.0000") %>%
                             hot_col("P_Value_of_X_and_Above", format = "0.0000")
                     }
                     else{
                         rhandsontable(stats$data_table)
                     })
                 })
    observeEvent(input$distribution,
                 {
                     output$distribution_display = renderPlot(if (input$Test == 3) {
                         ggplot(aes(x = 1:100, y = data), data = plotdata$data) +
                             geom_line() +
                             theme_void()
                     } else{
                         ggplot(aes(x = data), data = plotdata$data) +
                             geom_histogram(color = "#E27D60",
                                            fill = "#E8A87C",
                                            binwidth = 1) +
                             scale_x_continuous(
                                 breaks = 1:input$value_range[2],
                                 limits = c(input$value_range[1] - 1, input$value_range[2] +
                                                1),
                                 name = 'X Values'
                             ) +
                             ylab('Frequency Count') +
                             theme_classic()
                     })
                 })
}

shinyApp(ui = ui, server = server)
