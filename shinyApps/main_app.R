
library(shiny)

server <- function(session, input, output) {
  
  print("main page")
  
}

shinyApp(ui=htmlTemplate("/srv/shiny-server/www/index.html"), server=server)
