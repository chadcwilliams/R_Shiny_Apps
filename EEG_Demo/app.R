
library(shiny)

ui = fluidPage(
  tags$head(tags$style(type = "text/css", ".irs {max-width: 946px;}")),
  sliderInput(
    inputId = 'delta_amp',
    label = 'Delta (2 Hz)',
    value = 1,
    min = 0,
    max = 10,
    step = 1,
    width = '95%'
  ),
  sliderInput(
    inputId = 'theta_amp',
    label = 'Theta (6 Hz)',
    value = 1,
    min = 0,
    max = 10,
    step = 1,
    width = '95%'
  ),
  sliderInput(
    inputId = 'alpha_amp',
    label = 'Alpha (11 Hz)',
    value = 1,
    min = 0,
    max = 10,
    step = 1,
    width = '95%'
  ),
  sliderInput(
    inputId = 'beta_amp',
    label = 'Beta (23 Hz)',
    value = 1,
    min = 0,
    max = 10,
    step = 1,
    width = '95%'
  ),
  plotOutput(
    outputId = 'allPlots',
    height = 650,
    width = 946
  ),
)

server = function(input, output) {
  output$allPlots = renderPlot({
    
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    library(grid)
    library(gridExtra)
    
    #Define Sine Waves
    delta = 2
    theta = 6
    alpha = 11
    beta = 23
    
    delta_x <- seq(0, delta * 2 * pi, length.out = 1000)
    theta_x <- seq(0, theta * 2 * pi, length.out = 1000)
    alpha_x <- seq(0, alpha * 2 * pi, length.out = 1000)
    beta_x <- seq(0, beta * 2 * pi, length.out = 1000)
    
    #Create Line plot data frames
    line_data = data.frame(
      Time = 1:1000,
      Delta = sin(delta_x) * input$delta_amp,
      Theta = sin(theta_x) * input$theta_amp,
      Alpha = sin(alpha_x) * input$alpha_amp,
      Beta = sin(beta_x) * input$beta_amp,
      EEG = sin(delta_x) * input$delta_amp + sin(theta_x) * input$theta_amp + sin(alpha_x) *
        input$alpha_amp + sin(beta_x) * input$beta_amp
    )
    
    #Create Line Plots
    delta_plot = ggplot(aes(x = Time, y = Delta), data = line_data) + scale_y_continuous(limits =
                                                                                           c(-10, 10)) + geom_line() + ggtitle('Delta (2 Hz)') + theme_void()
    theta_plot = ggplot(aes(x = Time, y = Theta), data = line_data) + scale_y_continuous(limits =
                                                                                           c(-10, 10)) + geom_line() + ggtitle('Theta (6 Hz)') + theme_void()
    alpha_plot = ggplot(aes(x = Time, y = Alpha), data = line_data) + scale_y_continuous(limits =
                                                                                           c(-10, 10)) + geom_line() + ggtitle('Alpha (11 Hz)') + theme_void()
    beta_plot = ggplot(aes(x = Time, y = Beta), data = line_data) + scale_y_continuous(limits =
                                                                                         c(-10, 10)) + geom_line() + ggtitle('Beta (23 Hz)') + theme_void()
    EEG_plot = ggplot(aes(x = Time, y = EEG), data = line_data) + geom_line() +
      ggtitle('Summed Sine Waves (EEG)') + theme_void()
    
    #Fast Fourier Transforms
    FFT = rep(0, 30)
    FFT[delta] = input$delta_amp / 2
    FFT[theta] = input$theta_amp / 2
    FFT[alpha] = input$alpha_amp / 2
    FFT[beta] = input$beta_amp / 2
    FFT_data = data.frame(Power = FFT,
                          Frequency = 1:30)
    
    FFT_plot = ggplot(aes(x = Frequency, y = Power), data = FFT_data) +
      geom_line() +
      scale_x_continuous(breaks = 1:30) +
      scale_y_continuous(limits = c(0, 10), expand = c(0, 0)) +
      ggtitle('Fast Fourier Transform Plot') +
      theme_classic()
    
    #Time-Frequency Wavelets
    WAV_data = as.data.frame(matrix(rep(0, len = 1000 * 30), nrow = 30))
    WAV_data[delta, 301:800] = input$delta_amp
    WAV_data[theta, 301:800] = input$theta_amp
    WAV_data[alpha, 301:800] = input$alpha_amp
    WAV_data[beta, 301:800] = input$beta_amp
    WAV_long = melt(WAV_data)
    WAV_long$Time = sort(rep(1:1000, 30))
    WAV_long$Frequency = rep(1:30, 1000)
    
    WAV_plot = ggplot(aes(
      x = as.numeric(Time),
      y = as.numeric(Frequency),
      z = as.numeric(value)
    ), data = WAV_long) +
      geom_contour_filled(bins = 9) +
      theme_minimal() + theme(text = element_text(size = 16)) +
      scale_fill_brewer(type = 'div',
                        palette = 'RdBu',
                        direction = -1) +
      scale_x_continuous(expand = c(0, 0)) + xlab('Time (ms)') +
      scale_y_continuous(
        breaks = seq(2, 31, 2),
        labels = seq(2, 30, 2),
        expand = c(0, 0)
      ) + ylab('Frequency (Hz)') +
      coord_cartesian(ylim = c(1, 30)) + ylab('Frequency (Hz)') +
      ggtitle('Time-Frequency Plot') +
      guides(fill = guide_legend(reverse = T))
    
    
    lay = rbind(c(1, 6, 6),
                c(2, 6, 6),
                c(3, 7, 7),
                c(4, 7, 7),
                c(5, 7, 7))
    grid.arrange(
      delta_plot,
      theta_plot,
      alpha_plot,
      beta_plot,
      EEG_plot,
      FFT_plot,
      WAV_plot,
      layout_matrix = lay
    )
  })
}

shinyApp(ui = ui, server = server)
