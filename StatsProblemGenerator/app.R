####################################################################
####              Written by Chad C. Williams, 2021             ####
####                   www.chadcwilliams.com                    ####
####################################################################

options(scipen = 999) #Remove scientific notation
library(rsconnect) #Shiny
library(BSDA) #z-test function
library(rhandsontable) #Data tables
library(ggplot2) #Plotting
library(faux) #Creating correlated data (rnorm_multi)
library(rstatix) #Dependency of faux

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
                               "Correlation & Regression" = 4,
                               "Single Participant Z-Test" = 3,
                               "Single Sample Z-Test" = 5,
                               "Single Sample T-Test" = 6
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
                       column(6, rHandsontableOutput("data_display")),
                       column(6, plotOutput('distribution_display'))
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
            #Correlation & Regression
            #Create Data
            data = rnorm_multi(
                n = input$num_of_participants,
                mu = c(
                    sample(input$value_range[1]:input$value_range[2], 1),
                    sample(input$value_range[1]:input$value_range[2], 1)
                ),
                sd = c(rnorm(
                    1, (input$value_range[2] - input$value_range[1]) / 5, .1
                ), rnorm(
                    1, (input$value_range[2] - input$value_range[1]) / 5, .1
                )),
                r = sample(seq(-1, 1, .01), 1),
                varnames = c('X', 'Y')
            )
            data$X = as.integer(data$X)
            data$Y = as.integer(data$Y)
            data$X_SD = c(sqrt(sum((
                data$X - mean(data$X)
            ) ^ 2) / dim(data)[1]), rep(NA, dim(data)[1] - 1))
            data$Y_SD = c(sqrt(sum((
                data$Y - mean(data$Y)
            ) ^ 2) / dim(data)[1]), rep(NA, dim(data)[1] - 1))
            plotdata$data = data
            
            #Create Table
            descriptives = data.frame(
                X_Mean = mean(data$X),
                X_SD = sqrt(sum((
                    data$X - mean(data$X)
                ) ^ 2) / dim(data)[1]),
                Y_Mean = mean(data$Y),
                Y_SD = sqrt(sum((
                    data$Y - mean(data$Y)
                ) ^ 2) / dim(data)[1]),
                SP = sum((data$X - mean(data$X)) * (data$Y - mean(data$Y))),
                COV = sum((data$X - mean(data$X)) * (data$Y - mean(data$Y))) /
                    dim(data)[1],
                r = cor(data$X, data$Y)
            )
            descriptives$by = descriptives$r * (descriptives$Y_SD / descriptives$X_SD)
            descriptives$ay = descriptives$Y_Mean - (round(descriptives$by, 4) *
                                                         descriptives$X_Mean)
            descriptives$bx = descriptives$r * (descriptives$X_SD / descriptives$Y_SD)
            descriptives$ax = descriptives$X_Mean - (round(descriptives$bx, 4) *
                                                         descriptives$Y_Mean)
            descriptives$SD_XPrime = round((round(descriptives$X_SD, 2) *
                                                (round(
                                                    sqrt(1 - (abs(
                                                        round(cor(data$X, data$Y), 4)
                                                    ) ^ 2)), 4
                                                ))), 2)
            descriptives$SD_Yprime = round((round(descriptives$Y_SD, 2) *
                                                (round(
                                                    sqrt(1 - (abs(
                                                        round(cor(data$X, data$Y), 4)
                                                    ) ^ 2)), 4
                                                ))), 2)
            
            #Set Outputs
            stats$data_table = descriptives
            output$data_display = renderRHandsontable(rhandsontable(as.data.frame(data)))
            output$stats_display = renderRHandsontable({
                
            })
            output$distribution_display = renderPlot({
                
            })
        }
        else if (input$Test == 5) {
            #Single Sample Z-Test
            #Create Data
            data = data.frame(
                Data = sample(
                    input$value_range[1]:input$value_range[2],
                    input$num_of_participants,
                    replace = TRUE
                ),
                n = input$num_of_participants,
                mu = c(
                    sample(1:20, 1),
                    rep(NA, input$num_of_participants - 1)
                )
            )
            data$mu[1] = (round(rnorm(
                1, mean(data$Data), sd(data$Data) / 5
            )))
            data$sigma = c((round(runif(1, 0.5, 2.5),2)),
                           rep(NA, input$num_of_participants - 1))
            dir = runif(1)
            if (dir < .5) {
                data$direction = c('Two-Tail',
                                   rep(NA, input$num_of_participants - 1))
                direction = 1
            } else if (dir < .75) {
                data$direction = c('One-Tail (lower)',
                                   rep(NA, input$num_of_participants - 1))
                direction = 2
            }else{
                data$direction = c('One-Tail (higher)',
                                   rep(NA, input$num_of_participants - 1))
                direction = 3
            }
            
            data$p_alpha = c(.05, rep(NA, input$num_of_participants - 1))
            data$X_Mean = c(mean(data$Data),
                            rep(NA, input$num_of_participants - 1))
            plotdata$data = as.data.frame(data$Data)
            #Create Stats
            descriptives = data_table = data.frame(
                SE = data$sigma[1] / sqrt(input$num_of_participants),
                z_Obs = (data$X_Mean[1] - data$mu[1]) / (data$sigma[1] /
                                                             sqrt(input$num_of_participants)),
                z_Crit = if (direction == 1) {
                    '+-1.96'
                } else if (direction == 2) {
                    '-1.645'
                } else{
                    '+1.645'
                }
            )
            descriptives$p_obs = pnorm(round((data$X_Mean[1] - data$mu[1]) / (data$sigma[1] / sqrt(input$num_of_participants)),4))
            if (descriptives$z_Obs>=0){descriptives$p_obs=1-descriptives$p_obs}
            if (direction == 1) {descriptives$p_obs=descriptives$p_obs*2}
                
            descriptives$p_alpha = .05
            if (direction == 1){
                descriptives$H0 = if (descriptives$p_obs < .05){'Reject'}else{'Retain'}
                descriptives$H1 = if (descriptives$p_obs < .05){'Accept'}else{'Suspend'}}
            else if (direction == 2){
                if (descriptives$z_Obs<0){
                    descriptives$H0 = if (descriptives$p_obs < .05){'Reject'}else{'Retain'}
                    descriptives$H1 = if (descriptives$p_obs < .05){'Accept'}else{'Suspend'}}
                else{
                    descriptives$H0 = 'Retain'
                    descriptives$H1 = 'Suspend'}
            } else {
                if (descriptives$z_Obs>0){
                    descriptives$H0 = if (descriptives$p_obs < .05){'Reject'}else{'Retain'}
                    descriptives$H1 = if (descriptives$p_obs < .05){'Accept'}else{'Suspend'}}
                else{
                    descriptives$H0 = 'Retain'
                    descriptives$H1 = 'Suspend'}
            }
            
            
            #Set Outputs
            stats$data_table = descriptives
            output$data_display = renderRHandsontable(rhandsontable(as.data.frame(t(data[1, 2:dim(data)[2]]))))
            output$stats_display = renderRHandsontable({
                
            })
           output$distribution_display = renderPlot({
                
            })
        } else {
            #Single Sample Z-Test
            #Create Data
            data = data.frame(
                Data = sample(
                    input$value_range[1]:input$value_range[2],
                    input$num_of_participants,
                    replace = TRUE
                ),
                n = input$num_of_participants,
                mu = c(
                    sample(1:20, 1),
                    rep(NA, input$num_of_participants - 1)
                )
            )
            data$mu[1] = (round(rnorm(
                1, mean(data$Data), sd(data$Data) / 2.5
            )))
            data$SS = sum(((data$Data - mean(data$Data))^2))
            dir = runif(1)
            if (dir < .5) {
                data$direction = c('Two-Tail',
                                   rep(NA, input$num_of_participants - 1))
                direction = 1
            } else if (dir < .75) {
                data$direction = c('One-Tail (lower)',
                                   rep(NA, input$num_of_participants - 1))
                direction = 2
            }
            else{
                data$direction = c('One-Tail (higher)',
                                   rep(NA, input$num_of_participants - 1))
                direction = 3
            }
            
            data$p_alpha = c(.05, rep(NA, input$num_of_participants - 1))
            data$X_Mean = c(mean(data$Data),
                            rep(NA, input$num_of_participants - 1))
            plotdata$data = as.data.frame(data$Data)
            #Create Stats
            t=if (direction == 1){t.test(data$Data,mu=data$mu[1])
            }else if (direction == 2) {t.test(data$Data,mu=data$mu[1],alternative = 'less')
                    }else{t.test(data$Data,mu=data$mu[1],alternative = 'greater')}
            descriptives = data_table = data.frame(
                SS = sum(((data$Data - mean(data$Data))^2)),
                s = sd(data$Data),
                SE = sd(data$Data) / sqrt(input$num_of_participants),
                df = input$num_of_participants-1,
                t_Obs = as.numeric(t[['statistic']]),
                t_Crit = 
                if (direction == 1) {
                    paste('+-',toString(round(qt(p=.975, df=input$num_of_participants-1),2)))
                } else if (direction == 2) {
                    round(qt(p=.05, df=input$num_of_participants-1),2)
                } else{
                    round(qt(p=.95, df=input$num_of_participants-1),2)
                }
            )
            
            descriptives$p_obs = as.numeric(t[['p.value']])
            descriptives$p_alpha = .05
            if (direction == 1){
                descriptives$H0 = if (descriptives$p_obs < .05){'Reject'}else{'Retain'}
                descriptives$H1 = if (descriptives$p_obs < .05){'Accept'}else{'Suspend'}}
            else if (direction == 2){
                if (descriptives$t_Obs<0){
                    descriptives$H0 = if (descriptives$p_obs < .05){'Reject'}else{'Retain'}
                    descriptives$H1 = if (descriptives$p_obs < .05){'Accept'}else{'Suspend'}}
                else{
                    descriptives$H0 = 'Retain'
                    descriptives$H1 = 'Suspend'}
            } else {
                if (descriptives$t_Obs>0){
                    descriptives$H0 = if (descriptives$p_obs < .05){'Reject'}else{'Retain'}
                    descriptives$H1 = if (descriptives$p_obs < .05){'Accept'}else{'Suspend'}}
                else{
                    descriptives$H0 = 'Retain'
                    descriptives$H1 = 'Suspend'}
            }
            
            
            #Set Outputs
            stats$data_table = descriptives
            output$data_display = renderRHandsontable(rhandsontable(as.data.frame(t(data[1, 2:dim(data)[2]]))))
            output$stats_display = renderRHandsontable({
                
            })
            output$distribution_display = renderPlot({
                
            })
            
            
            
            
            
            
            
            
            
            
            
            
            
            }
    })
    observeEvent(input$answers,
                 {
                     output$stats_display = renderRHandsontable(if (input$Test == 3) {
                         rhandsontable(stats$data_table) %>%
                             hot_col("P_Value_of_X_and_Below", format = "0.0000") %>%
                             hot_col("P_Value_of_X_and_Above", format = "0.0000")
                     }else if (input$Test == 5 | input$Test == 6){
                         rhandsontable(stats$data_table) %>%
                             hot_col("p_obs", format = "0.0000")
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
                             geom_vline(xintercept = round((
                                 stats$data_table$P_Value_of_X_and_Below * 100
                             )) + .5, color = 'red') +
                             theme_void()
                     } else if (input$Test == 4) {
                         ggplot(aes(x = X, y = Y), data = plotdata$data) +
                             geom_point(size = 4, alpha = .5) +
                             geom_smooth(method = lm, se = F) +
                             geom_segment(
                                 y = min(plotdata$data$Y),
                                 x = (stats$data_table$ax + (
                                     stats$data_table$bx * min(plotdata$data$Y)
                                 )),
                                 yend = max(plotdata$data$Y),
                                 xend = (stats$data_table$ax + (
                                     stats$data_table$bx * max(plotdata$data$Y)
                                 )),
                                 color = 'red'
                             ) +
                             theme_classic() +
                             theme(text = element_text(size = 20))
                     }
                     else{
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
