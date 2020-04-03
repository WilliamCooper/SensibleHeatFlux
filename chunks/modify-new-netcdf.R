## modify-new-netcdf.R

## variables needed for attributes 
VarList <- ATS
VarListRef <- VarList
FI <- DataFileInfo (fname)
VarList <- VarListRef

for (Var in VarList) {
  if (!(Var %in% FI$Variables)) {
    print (sprintf (' required variable %s not found in file %s; aborting...', Var, fname))
    exit()
  }
}

netCDFfile <- nc_open (fnew, write=TRUE) 
Rate <- 1
Dimensions <- attr (DXC, "Dimensions")
Dim <- Dimensions[["Time"]]
if ("sps25" %in% names (Dimensions)) {
  Rate <- 25
  Dim <- list(Dimensions[["sps25"]], Dimensions[["Time"]])
}
if ("sps50" %in% names (Dimensions)) {
  Rate <- 50
  Dim <- list(Dimensions[["sps50"]], Dimensions[["Time"]])
}
DATT <- DXC  ## save to ensure that attributes are preserved

## variables to add to the netCDF file:
VarOld <- ATS
VarNew <- paste0(ATS, 'C')
VarUnits <- vector('character', length(ATS))
VarStdName <- vector('character', length(ATS))
VarLongName <- vector('character', length(ATS))
for (i in 1 : length(ATS)) {
  VarUnits[i] <- attr(DXC[, ATS[i]], 'units')  
  VarStdName[i] <- attr(DXC[, ATS[i]], 'standard_name')
  VarLongName[i] <- paste0(attr(DXC[, ATS[i]], 'long_name'), ', corrected')
}

## create the new variables
varCDF <- list ()
for (i in 1:length(VarNew)) {
  print (sprintf ('new-netcdf %d%% done', as.integer(100*(i-1)/length(VarNew))))
  varCDF[[i]] <- ncvar_def (VarNew[i],  
                            units=VarUnits[i], 
                            dim=Dim, 
                            missval=as.single(-32767.), prec='float', 
                            longname=VarLongName[i])
  if (i == 1) {
    newfile <- ncvar_add (netCDFfile, varCDF[[i]])
  } else {
    newfile <- ncvar_add (newfile, varCDF[[i]])
  }
  ATV <- ncatt_get (netCDFfile, VarOld[i])
  copy_attributes (ATV, VarNew[i], newfile)
  ncatt_put (newfile, VarNew[i], attname="standard_name", 
             attval=VarStdName[i])
  if (Rate == 1) {
    ncvar_put (newfile, varCDF[[i]], DXC[, VarNew[i]])
  } else if (Rate == 25) {
    ncvar_put (newfile, varCDF[[i]], DXC[, VarNew[i]], count=c(25, nrow(DXC)/25))
  }
}
nc_close (newfile)

