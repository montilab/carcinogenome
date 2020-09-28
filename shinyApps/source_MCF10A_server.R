
library(shiny)

dat = readRDS("data/MCF10A/data.RDS")
runApp("app.R", host='0.0.0.0', port=4848)

