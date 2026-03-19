library(stringr)

# Source the package function for development/testing
source(file.path(dirname(dirname(dirname(dirname(
  rstudioapi::getSourceEditorContext()$path
)))), "R", "sg_translator.R"))
