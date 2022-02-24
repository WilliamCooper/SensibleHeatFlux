#' Package Ranadu
#' 
#' This is a package of programs designed for use with the data archives
#' that contain measurements from research flights of the NCAR/EOL/RAF
#' aircraft, normally archived in netCDF format. These routines read those
#' archives, construct R data.frames containing selected variables, and
#' provide utility functions for working with the measurements.
"_PACKAGE"

if (!exists('VSpecEnv', envir=emptyenv())) {
  VSpecEnv <- new.env(parent = emptyenv())  # define if absent
}
# This is just to keep roxygen2 from complaining about no-visible-binding when these are OK:
utils::globalVariables(c(".x", "Time", "cf", "cospec", "fpf2", "fpf3", "ncospec", "ogive", "x", 
  "xc", "y", "ybar", "ymax", "ymin", "freq"))
