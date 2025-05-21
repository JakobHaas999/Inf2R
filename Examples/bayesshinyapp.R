library(shiny)

ui <- fluidPage(
  titlePanel("Bayesian Updating: Beta Prior and Binomial Likelihood"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Prior: Beta(a, b)"),
      numericInput("a", "a (Prior)", value = 1, min = 0.1, step = 0.1),
      numericInput("b", "b (Prior)", value = 1, min = 0.1, step = 0.1),
      hr(),
      h4("Data: Binomial(m+n, p)"),
      numericInput("m", "m (Observed successes)", value = 1, min = 0),
      numericInput("n", "n (Observed failures)", value = 10, min = 0),
      hr(),
      checkboxInput("proportional", "Show proportional posterior (unnormalized)", value = FALSE),
      hr(),
      h4("Interval for CDF"),
      sliderInput("interval", "Select interval (lower and upper bounds):",
                  min = 0, max = 1, value = c(0.083, 0.1), step = 0.01),
      verbatimTextOutput("intervalProb"),
      helpText("Posterior: Beta(a + m, b + n)")
    ),
    
    mainPanel(
      plotOutput("posteriorPlot"),
      br(),
      plotOutput("posteriorCDFPlot")
    )
  )
)

server <- function(input, output, session) {
  
  output$posteriorPlot <- renderPlot({
    
    p <- seq(0, 1, length.out = 500)
    prior <- dbeta(p, input$a, input$b)
    
    n_total <- input$m + input$n
    likelihood <- dbinom(input$m, size = n_total, prob = p)
    likelihood_scaled <- likelihood / max(likelihood) * max(prior)
    
    if (input$proportional) {
      posterior <- prior * likelihood
      posterior <- posterior / max(posterior) * max(prior)
    } else {
      posterior <- dbeta(p, input$a + input$m, input$b + input$n)
    }
    
    plot(p, prior, type = "l", lwd = 2, col = "blue",
         ylim = range(c(prior, likelihood_scaled, posterior)),
         xlab = "Probability p", ylab = "Density (up to proportionality)",
         main = ifelse(input$proportional, 
                       "Prior × Likelihood (Proportional Posterior)", 
                       "Prior, Likelihood, and Posterior"))
    lines(p, likelihood_scaled, lwd = 2, col = "orange")
    lines(p, posterior, lwd = 2, col = "darkgreen")
    
    legend("topright", legend = c("Prior (Beta)", "Likelihood (scaled Binomial)", 
                                  ifelse(input$proportional, "Proportional Posterior", "Posterior (Beta)")),
           col = c("blue", "orange", "darkgreen"), lwd = 2)
  })
  
  output$posteriorCDFPlot <- renderPlot({
    p <- seq(0, 1, length.out = 500)
    
    if (input$proportional) {
      plot(p, rep(NA, length(p)), type = "n", xlab = "Probability p", ylab = "CDF",
           main = "CDF only for normalized posterior")
      text(0.5, 0.5, "CDF not available for proportional posterior", cex = 1.2)
    } else {
      cdf <- pbeta(p, input$a + input$m, input$b + input$n)
      
      plot(p, cdf, type = "l", lwd = 2, col = "darkgreen",
           xlab = "Probability p", ylab = "CDF",
           main = "Posterior CDF (Beta Distribution)")
      
      # Add vertical lines
      abline(v = input$interval[1], col = "red", lwd = 2, lty = 2)
      abline(v = input$interval[2], col = "red", lwd = 2, lty = 2)
    }
  })
  
  output$intervalProb <- renderText({
    if (!input$proportional) {
      prob <- pbeta(input$interval[2], input$a + input$m, input$b + input$n) - 
        pbeta(input$interval[1], input$a + input$m, input$b + input$n)
      paste0("P(", round(input$interval[1], 3), " < p < ", 
             round(input$interval[2], 3), ") = ", round(prob, 4))
    } else {
      "Probability not available for proportional posterior."
    }
  })
}

shinyApp(ui, server)
