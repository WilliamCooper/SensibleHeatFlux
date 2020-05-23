## This chunk loads some often used R packages and knitr options

library(knitr)
opts_chunk$set(echo = FALSE,
               include = FALSE,
               fig.lp = "fig:")
# note that fig.pos="center" gave errors, changed to fig.align
opts_chunk$set(
  fig.width = 6,
  fig.height = 5,
  fig.align = "center",
  digits = 4
)
thisFileName <- "ThisFileName"    ## change this
library(Ranadu, quietly = TRUE, warn.conflicts = FALSE)
library(scales)
library(reshape2)    ## used with ggplot facet plots
library(grid)
library(magrittr)    ## used for pipes (%>%)
library(dplyr)
options(stringsAsFactors = FALSE)