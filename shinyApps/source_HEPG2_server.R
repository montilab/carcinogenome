
library(shiny)

dat = readRDS("data/HEPG2/data.RDS")
runApp("app.R", host='0.0.0.0', port=8787)
