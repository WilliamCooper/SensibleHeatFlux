## Using this script: from a unix shell:
## Rscript CorrectTemperature.R Project Flight
## (If run interactively from an R command window, e.g., in RStudio,
##  the script will request this input)
## The Flight should be in the form rf15h, not a simple number and
## without the trailing .nc. Two special "Flight" inputs are these:
##   NEXT -- processes the next unprocessed high-rate file, If a high-rate
##           file has an associated "..T.nc" file, it will be skipped.
##   ALL  -- process all unprocessed high-rate files in the project directory.

require(Ranadu, quietly = TRUE, warn.conflicts=FALSE)
# needed packages
library(zoo)
require(signal)
load(file='ARF.Rdata')    ## the filters
load(file='PAR.Rdata')  ## the response parameters

## Get file to process:

## Specify the Project and Flight:
Project <- 'SOCRATES'
Flight <- 'rf15h'  # should be high-rate
ProjectDir <- Project
getNext <- function(ProjectDir, Project) {
  Fl <- sort (list.files (sprintf ("%s%s/", DataDirectory (), ProjectDir),
    sprintf ("%srf..hT.nc", Project)), decreasing = TRUE)[1]
  if (is.na (Fl)) {
    Flight <- 'rf01h'
  } else {
    Flight <- sub (".*rf", '',  sub ("hT.nc", '', Fl))
    Flight <- as.numeric(Flight)+1
    Flight <- sprintf ('rf%02dh', Flight)
  }
  return (Flight)
}
if (!interactive()) {  ## can run interactively or via Rscript
  run.args <- commandArgs (TRUE)
  if (length (run.args) > 0) {
    if (nchar(run.args[1]) > 1) {
      Project <- run.args[1]
      ProjectDir <- Project
    }
  } else {
    print ("Usage: Rscript CorrectTemperature.R Project Flight")
    print ("Example: Rscript CorrectTemperature SOCRATES rf15h")
    stop("exiting...")
  }
  ## Flight
  if (length (run.args) > 1) {
    ALL <- FALSE
    if (run.args[2] == 'ALL') {
      ALL <- TRUE
    } else if (run.args[2] != 'NEXT') {
      Flight <- run.args[2]
    } else {
      ## Find max already-processed rf..h in data directory,
      ## Use as default if none supplied via command line:
      Flight <- getNext(ProjectDir, Project)
    }
  }
} else {
  x <- readline (sprintf ("Project is %s; CR to accept or enter new project name: ", Project))
  if (nchar(x) > 1) {
    Project <- x
    if (grepl ('HIPPO', Project)) {
      ProjectDir <- 'HIPPO'
    } else {
      ProjectDir <- Project
    }
  }
  x <- readline (sprintf ("Flight is %s; CR to accept, number 'ALL' or 'NEXT' for new flight name: ",
    Flight))
  ALL <- FALSE
  if (x == 'ALL') {
    ALL <- TRUE
  } else if (x == 'NEXT') {
    Flight <- getNext(ProjectDir, Project)
  } else if (nchar(x) > 0 && !is.numeric(x)) {
    Flight <- x
  } else if (nchar(x) > 0 && is.numeric(x)) {
    Flight <- sprintf('rf%02dh', as.numeric(x))
  } 
}
## A function to transfer attributes:
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

correctDynamicHeating <- function(D, AT) {
  platform <- attr(D, 'Platform')
  platform <- ifelse (grepl('677', platform), 'GV', 'C130')
  heated <- attr(D[, AT], 'long_name')
  heated <- ifelse(grepl('nheated', heated), FALSE, TRUE)
  if (grepl('RTF', RT)) {heated <- FALSE}
  # Select the right filter
  Rate <- attr(D, 'Rate')
  if (Rate == 25) {
    if (platform == 'GV') {
      if (heated) {
        AF <- ARH
        Lsh <- LshiftH
      } else {
        AF <- ARG
        Lsh <- Lshift
      }
    } else {
      if (heated) {
        AF <- ARH
        Lsh <- LshiftH
      } else {
        AF <- AR
        Lsh <- Lshift
      }
    }
  } else { ## assumed 1
    if (platform == 'GV') {
      if (heated) {
        AF <- ARH1
        Lsh <- LshiftH1
      } else {
        AF <- ARG1
        Lsh <- Lshift1
      }
    } else {
      if (heated) {
        AF <- ARH1
        Lsh <- LshiftH1
      } else {
        AF <- AR1
        Lsh <- Lshift1
      }
    }
  }
  # get the recovery factor from the attribute:
  rf.txt <- attr(D[, AT], 'RecoveryFactor')
  if (grepl('mach', rf.txt)) {
    rf <- gsub('mach', 'MACHX', rf.txt)
    rf <- gsub(' log', ' * log', rf)
    rf <- gsub(' \\(', ' * \\(', rf)
    rf <- with(D, eval(parse(text=rf)))
  } else {
    rf <- 0.985  ## placeholder -- if needed, add decoding here
  }
  # Is the associated recovery temperature present?
  dep <- attr(D[, AT], 'Dependencies')
  RT <- gsub(' .*', '', gsub('^. ', '', dep))
  # Get the dynamic-heating correction:
  if ('EWX' %in% names(D)) {
    ER <- SmoothInterp(D$EWX / D$PSXC, .Length = 0)
    CP <- SpecificHeats(ER)
  } else {
    CP <- SpecificHeats()
  }
  if (RT %in% names(D)) {
    X <- rf * D$MACHX^2 * CP[, 3] / (2 * CP[, 2])
    Q <- (273.15 + D[, RT]) * X / (1 + X)
  } else {
    Q <- rf * D$TASX^2 / 2010
  }
  # Interpolate to avoid missing values
  Q <- SmoothInterp(Q, .Length = 0)
  # Filter
  QF <- as.vector(signal::filter(AF, Q))
  QF <- ShiftInTime(QF, .shift=(-(Lsh + 1) * 40), .rate = Rate)
  # Change the measured air temperature by adding back Q and subtracting QF
  ATC <- D[, AT] + Q - QF
  return(ATC)
}

## The following applies either to recovery temperature or to air temperature after
## filtering the dynamic-heating term:
correctTemperature <- function(D, RT, responsePar = 
                                 list(a = 0.733, tau1 = 0.0308, tau2 = 0.447)){
  Rate <- attr(D, 'Rate')
  RTT <- D[, RT]
  # Protect against missing values by interpolating:
  RTT <- zoo::na.approx(as.vector(RTT), maxgap = 1000 * Rate, na.rm = FALSE, rule = 2)
  RTT[is.na(RTT)] <- 0
  heated <- attr(D[, RT], 'long_name')
  heated <- ifelse(grepl('nheated', heated), FALSE, TRUE)
  a <- responsePar$a
  tau1 <- responsePar$tau1
  tau2 <- responsePar$tau2
  if (heated) {
    DTMDT <- c(0, diff(RTT, 2), 0) * Rate / 2
    DTM2DT2 <- (c(diff(RTT), 0) - c(0, diff(RTT))) * Rate ^ 2
    RTC <- (tau1 + tau2) * DTMDT + RTT + tau1 * tau2 * DTM2DT2
    RTC <- zoo::na.approx (
      as.vector(RTC),
      maxgap = 1000 * Rate,
      na.rm = FALSE,
      rule = 2
    )
    CutoffPeriod <- 12.5
    RTC <- signal::filtfilt (signal::butter (3,
                          2 / CutoffPeriod), RTC)
  } else {
    ## RTT is the working solution; Ts is the support temperature
    Ts <- RTT
    Rate <- attr (D, 'Rate')
    # DTMDT <- c(0, diff(RT, 2), 0) * Rate / 2
    DTMDT <-  (c(0, 8 * diff(RTT, 2), 0) -
                    c(0, 0, diff(RTT, 4), 0, 0)) * Rate / 12
    fS <- function(y, i) {
      ((tau1 * DTMDT[i] + RTT[i] - (1 - a) * y) / a - y) / (Rate * tau2)
    }
    
    Ts <- rk4.integrate (fS, Ts[1], 1:nrow(D))
    RTC <- (1 / a) * (tau1 * DTMDT + RTT - (1 - a) * Ts)
  }
  return(RTC)
}
processFile <- function(ProjectDir, Project, Flight) {
  ## Find the available air_temperature variables:
  fname <- file.path(DataDirectory(), sprintf('%s/%s%s.nc',
    ProjectDir, Project, Flight))
  FI <- DataFileInfo(fname, LLrange = FALSE)
  TVARS <- FI$Measurands$air_temperature
  TVARS <- TVARS[-which ('ATX' == TVARS)]  # don't include ATX
  RVARS <- sub('^A', 'R', TVARS)
  ## get the old netCDF variables needed to calculate the modified variables
  VarList <- standardVariables(TVARS)
  ## Add the recovery temperatures if present; otherwise
  ## recalculate them:
  for (RV in RVARS) {
    if (RV %in% FI$Variables) {
      VarList <- c(VarList, RV)
    }
  }
  D <- getNetCDF (fname, VarList)
  Rate <- attr(D, 'Rate')

  ## Calculate the new variables:
  E <- SmoothInterp(D$EWX / D$PSXC, .Length = 0)
  D$Cp <- SpecificHeats (E)[, 1]
  D$DH <- D$TASX^2 / (2 * D$Cp)
  ## Recalculate recovery temperatures if missing:
  ## (Then, will have to modify how recovery factor is retrieved because not available...)
  retrievedRVARS <- rep(FALSE, length(RVARS))
  for (RV in RVARS) {
    if(!(RV %in% VarList)) {
      TV <- sub('^R', 'A', RV)
      prb <- 'HARCO'
      if (TV == 'ATH2') {prb <- 'HARCOB'}
      if (grepl('ATF', TV)) {prb <- 'UNHEATED'}
      retrievedRVARS[which(RV == RVARS)] <- TRUE
      D[, RV] <- D[, TV] + RecoveryFactor(D$MACHX, prb) * D$DH
    }
  }
  platform <- attr(D, 'Platform')
  platform <- ifelse (grepl('677', platform), 'GV', 'C130')
  # filter DH:
  CutoffPeriod <- rep(Rate * 1.0, length(TVARS)) # Standard is 1 s for DH filtering
  probe <- rep('HARCO', length(TVARS))  # used to determine the recovery factor
  PAR <- data.frame()
  for (i in 1:length(TVARS)) {
    PAR <- rbind(PAR, ParamSH)
  }
  # check for ATF
  ic <- which(grepl('ATF', TVARS))
  if (length(ic) > 0) {
    CutoffPeriod[ic] <- Rate * 0.5  # 0.5 s for ATFx
    probe[ic] <- 'UNHEATED'
    if (platform == 'GV') {
      PAR[ic, ] <- as.data.frame(ParamSF)
    } else {
      PAR[ic, ] <- as.data.frame(Param1)
    }
  }
  ic <- which(grepl('ATH2', TVARS))
  if (length(ic) > 0) {
    probe[ic] <- 'HARCOB'
    PAR[ic, ] <- as.data.frame(ParamSH)
  }
  if (Rate == 1) {  # Protection against script failure for a LRT file
    CutoffPeriod[CutoffPeriod == 1] <- 2.2
    CutoffPeriod[CutoffPeriod == 0.5] <- 2.0
  }

  DHM <- rep(D$DH, length(TVARS))
  dim(DHM) <- c(length(D$DH), length(TVARS))
  newTVARS <- paste0(TVARS, 'C')
  filteredTVARS <- paste0(TVARS, 'F')
  for (i in 1:length(TVARS)) {
    D[, filteredTVARS[i]] <- correctDynamicHeating(D, TVARS[i])
  }

  ## Calculate the corrected temperatures:
  for (i in 1:length(TVARS)) {
    ATC <- correctTemperature(D, filteredTVARS[i], responsePar = PAR[i,])
    D[, newTVARS[i]] <- ATC
  }
  ## Add the new variables:
  newTVARS <- c(newTVARS, filteredTVARS)

  ## create and open a copy of the old file for writing:
  fnew <- sub ('.nc', 'T.nc', fname)
  ## beware: overwrites without warning!!
  Z <- file.copy (fname, fnew, overwrite=TRUE)
  netCDFfile <- nc_open (fnew, write=TRUE)
  ## retrieve dimension info from the old file
  Dimensions <- attr (D, "Dimensions")
  Dim <- Dimensions[["Time"]]
  if ("sps25" %in% names (Dimensions)) {
    Rate <- 25
    Dim <- list(Dimensions[["sps25"]], Dimensions[["Time"]])
  }
  DATT <- D

  ## variables to add to the netCDF file:
  VarNew <- newTVARS
  VarOld <- c(TVARS, TVARS)
  VarUnits <- rep("deg_C", 2*length(TVARS))
  VarLongName <- c(rep("Ambient Temperature, corrected", length(TVARS)), 
                   rep("Ambient Temperature, filtered", length(TVARS)))
  VarStdName <- rep("air_temperature", 2*length(TVARS))
  Dependencies <- rep(sprintf('2 %s TASX', TVARS[1]), 2*length(TVARS))
  for (i in 1:(2*length(TVARS))) {
    Dependencies[i] <- sprintf('2 %s TASX', TVARS[i])
  }

  ## create the new variables
  varCDF <- list ()
  for (i in 1:length(VarNew)) {
    ## create the new variable and add it to the netCDF file
    varCDF[[i]] <- ncvar_def (VarNew[i],
      units=VarUnits[i],
      dim=Dim,
      missval=as.single(-32767.), prec="float",
      longname=VarLongName[i])
    if (i == 1) {
      newfile <- ncvar_add (netCDFfile, varCDF[[i]])
    } else {
      newfile <- ncvar_add (newfile, varCDF[[i]])
    }
    ## transfer attributes from the old variable and add new ones
    ATV <- ncatt_get (netCDFfile, VarOld[i])
    copy_attributes (ATV, VarNew[i], newfile)
    ncatt_put (newfile, VarNew[i], attname="standard_name",
      attval=VarStdName[i])
    ncatt_put (newfile, VarNew[i], attname="Dependencies",
      attval=Dependencies[i])
    ## add the measurements for the new variable
    if (Rate == 1) {
      ncvar_put (newfile, varCDF[[i]], D[, VarNew[i]])
    } else if (Rate == 25) {
      ncvar_put (newfile, varCDF[[i]], D[, VarNew[i]],
        count=c(25, nrow(D)/25))
    }
  }
  ## then close to write the new file
  nc_close (newfile)
}

if (ALL) {
  Fl <- sort (list.files (sprintf ("%s%s/", DataDirectory (), ProjectDir),
    sprintf ("%srf..h.nc", Project)), decreasing = TRUE)
  for (flt in Fl) {
    fcheck <- file.path(DataDirectory(), ProjectDir, '/', flt, fsep = '')
    fcheck <- sub('.nc', 'T.nc', fcheck)
    if (file.exists(fcheck)) {
      print (sprintf('processed file %s exists; skipping', flt))
    } else {
      print (sprintf('processing file %s', flt))
      processFile (ProjectDir, Project, sub('.*rf', 'rf', sub('.nc', '', flt)))
    }
  }
} else {
  processFile (ProjectDir, Project, Flight)
}

