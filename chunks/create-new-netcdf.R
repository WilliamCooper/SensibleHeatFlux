## create-new-netcdf.R

fnew <- sub ('TEMP.nc', 'TC.nc', fnameW)
## beware: overwrites without warning!!
Z <- file.copy (fnameW, fnew, overwrite=TRUE)  ## BEWARE: overwrites without warning!!

# function to copy attributes from old variable (e.g., PITCH) to new one (e.g., PITCHKF)
copy_attributes <- function (atv, v, nfile) {
  for (i in 1:length(atv)) {
    aname <- names(atv[i])
    if (grepl ('name', aname)) {next}  # skips long and standard names
    if (grepl ('units', aname)) {next}
    if (grepl ('Dependencies', aname)) {next}
    if (grepl ('actual_range', aname)) {next}
    if (is.numeric (atv[[i]])) {
      ncatt_put (nfile, v, attname=aname, attval=as.numeric(atv[[i]]))
    } else {
      ncatt_put (nfile, v, attname=aname, attval=as.character (atv[[i]]))
    }
  }
}


