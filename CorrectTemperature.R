## Using this script: from a unix shell:
## Rscript CorrectTemperature.R Project Flight
## (If run interactively from an R command window, e.g., in RStudio,
##  the script will request this input)
## The Flight should be in the form rf15h, not a simple number and
## without the trailing .nc. (Some untested code tries to adapt to
## single-number input or trailing .nc.)  Two special "Flight" inputs are these:
##   NEXT -- process the next unprocessed high-rate file, If a high-rate
##           file has an associated "..T.nc" file, it will be skipped.
##   NEXTLR - process the next unprocessed low-rate file.
##   ALL  -- process all unprocessed high-rate files in the project directory.
##   ALLLR - process all unprocessed low-rate files.
## Additional optional arguments are "RTN", "FFT" and "UH1", all FALSE by 
## default. There use is as follows:
##   RTN  -- For files that don't contain a recovery temperature, it must be
##           recalculated from the temperature and dynamic heating. This
##           argument causes the recalculated variable to be preserved in the
##           new netCDF file.
##   FFT  -- Processing via Fourier transformation is normally suppressed and
##           the associated "C2" variables below are not added to the new
##           netCDF file. If those variables are desired, they will be added
##           if the "FFT" argument is supplied.
##   UH1  -- Normally for 1-Hz files the change provided by filtering dynamic-
##           heating does not make a significant change for an unheated sensor,
##           so the filtered variables ("F" below) are omitted. This argument
##           will cause them to be included in the new netCDF file.
## Processed files will duplicate the original with the addition of these
## variables (where [name] is the original name; e.g., ATF1)
##  [name]C --  corrected temperature using method 1 (differential equations)
##  [name]C2 -- corrected temperature using method 2 (Fourier transforms)
##  [name]F  -- filtered dynamic heating applied to original, but otherwise
##              uncorrected
## The processed file with have "T" added to the name; e.g., SOCRATESrf15hT.nc .

require(Ranadu, quietly = TRUE, warn.conflicts=FALSE)
# needed packages
library(zoo)
require(signal)
load(file='AR.Rdata')    ## the filters
load(file='PAR.Rdata')  ## the response parameters

## Get file to process:

## Specify the Project and Flight:
Project <- 'SOCRATES'
Flight <- 'rf15h'  
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
getNextLR <- function(ProjectDir, Project) {
  Fl <- sort (list.files (sprintf ("%s%s/", DataDirectory (), ProjectDir),
                          sprintf ("%srf..T.nc", Project)), decreasing = TRUE)[1]
  if (is.na (Fl)) {
    Flight <- 'rf01'
  } else {
    Flight <- sub (".*rf", '',  sub ("T.nc", '', Fl))
    Flight <- as.numeric(Flight)+1
    Flight <- sprintf ('rf%02d', Flight)
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
    ALLLR <- FALSE
    if (run.args[2] == 'ALLLR') {
      ALLLR <- TRUE
      Flight <- 'ALLLR'
    } else if (run.args[2] == 'ALL') {
      ALL <- TRUE
      Flight <- 'ALL'
    } else if (run.args[2] != 'NEXT' && (run.args[2] != 'NEXTLR')) {
      Flight <- run.args[2]
      Flight <- sub('.nc$', '', Flight)
    } else {
      if (run.args[2] == 'NEXTLR') {
        Flight <- getNextLR(ProjectDir, Project)
      } else {
        ## Find max already-processed rf..h in data directory,
        ## Use as default if none supplied via command line:
        Flight <- getNext(ProjectDir, Project)
      }
    }
  }
  FFT <- FALSE
  RTN <- FALSE
  UH1 <- FALSE
  if (length (run.args) > 2) {
    if (any(run.args == 'FFT')) {FFT <- TRUE}
    if (any(run.args == 'RT')) {RTN <- TRUE}
    if (any(run.args == 'UH1')) {UH1 <- TRUE}
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
  ALLLR <- FALSE
  if (x == 'ALLLR') {
    ALLLR <- TRUE
    Flight <- 'ALLLR'
  } else if (x == 'ALL') {
    ALL <- TRUE
    Flight <- 'ALL'
  } else if (x == 'NEXTLR') {
    Flight <- getNextLR(ProjectDir, Project)
  } else if (x == 'NEXT') {
    Flight <- getNext(ProjectDir, Project)
  } else if (nchar(x) > 0 && !is.numeric(x)) {
    Flight <- x
    Flight <- sub('.nc$', '', Flight)
    if (!
        grepl('rf', Flight)) {Flight <- paste0('rf', Flight)}
    # if (!grepl('[0-9][0-9]h', Flight)) {Flight <- paste0(Flight, 'h')}
  } else if (nchar(x) > 0 && is.numeric(x)) {
    Flight <- sprintf('rf%02dh', as.numeric(x))
  }
  x <- readline ("Add reconstructed recovery temperatures if missing? (N/y):")
  RTN <- FALSE
  if ((x != '') && (x != 'N')) {RTN <- TRUE}
  x <- readline ("Include FFT-correctedtemperatures? (N/y):")
  FFT <- FALSE
  if ((x != '') && (x != 'N')) {FFT <- TRUE}
  x <- readline ("Include processing for unheated sensors if 1-Hz? (N/y):")
  UH1 <- FALSE
  if ((x != '') && (x != 'N')) {UH1 <- TRUE}
}

print (sprintf ('run parameters: Project = %s, Flight = %s, FFT = %s RTN = %s UH1 = %s',
                Project, Flight, FFT, RTN, UH1))
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

recoveryF <- function(D, TV) {  ## TV should be a temperature column in D
  # get the recovery factor from the attribute:
  rf.txt <- attr(TV, 'RecoveryFactor')
  if (is.null (rf.txt)) {  ## protection for omitted attribute
    rf <- 0.985
  } else {
    if (is.numeric(rf.txt)) {
      rf <- rf.txt
    } else {
      if (grepl('mach', rf.txt)) {
        rf <- gsub('mach', 'MACHX', rf.txt)
        rf <- gsub(' log', ' * log', rf)
        rf <- gsub(' \\(', ' * \\(', rf)
        rf <- with(D, eval(parse(text=rf)))
        rf <- SmoothInterp(rf, .Length = 0)
      } else {
        rf <- 0.985  ## placeholder -- if needed, add decoding here
      }
    }
  }
  return (rf)
}

correctDynamicHeating <- function(D, AT) {
  platform <- attr(D, 'Platform')
  platform <- ifelse (grepl('677', platform), 'GV', 'C130')
  heated <- attr(D[, AT], 'long_name')
  heated <- ifelse(grepl('nheated', heated), FALSE, TRUE)
  if (grepl('^ATR', AT)) {heated <- FALSE}  ## old name convention
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
    } else {  ## - C130 case
      if (heated) {
        AF <- ARH
        Lsh <- LshiftH
      } else {
        AF <- AR
        Lsh <- Lshift
      }
    }
  } else { ## assumed 1 Hz
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

  rf <- recoveryF(D, AT)
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
  QF <- ShiftInTime(QF, .shift=(-(Lsh) * 1000 / Rate), .rate = Rate)
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
  ## This indication is often missing from RTFx:
  if ((grepl('RTF', RT)) || (grepl('ATF', RT)) || (grepl('^ATR', RT)) || (grepl('^RTR', RT))) {
    heated <- FALSE
  }
  a <- responsePar$a
  tau1 <- responsePar$tau1
  tau2 <- responsePar$tau2
  if (heated) {
    # DTMDT <- c(0, diff(RTT, 2), 0) * Rate / 2
    # DTM2DT2 <- (c(diff(RTT), 0) - c(0, diff(RTT))) * Rate ^ 2
    DTMDT <-  (c(0, 8 * diff(RTT, 2), 0) -
                 c(0, 0, diff(RTT, 4), 0, 0)) * Rate / 12
    DTM2DT2 <- (-c(0,0, RTT)[1:nrow(D)] + 16*c(0,RTT)[1:nrow(D)] 
                - 30 * RTT + 16 * c(RTT[2:nrow(D)], 0) 
                - c(RTT[3:nrow(D)], 0, 0)) * (Rate^2) / 12
    RTC <- (tau1 + tau2) * DTMDT + RTT + tau1 * tau2 * DTM2DT2
    RTC <- zoo::na.approx (
      as.vector(RTC),
      maxgap = 1000 * Rate,
      na.rm = FALSE,
      rule = 2
    )
    CutoffPeriod <- ifelse(Rate == 25, Rate / 2, 5)
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


addFFTsoln <- function(D, RV, responsePar) {
  j1 <- which(D$TASX > 90)[1]
  j2 <- which(D$TASX > 90)
  j2 <- j2[length(j2)]
  DS <- D[j1:j2,]
  Rate <- attr(D, 'Rate')
  DS$TASX <- SmoothInterp(DS$TASX, .Length = 0)
  DS$QCF <- SmoothInterp(DS$QCF, .Length = 0)
  DS[, RV] <- SmoothInterp(DS[, RV], .Length = 0)
  RVC <- paste0(RV, 'C')
  DS[, RVC] <- DS[, RV]  ## column to hold the corrected values
  DD <- D
  DD[, RVC] <- DD[, RV]
  N <- ifelse(Rate == 25, 2^14, 2*9)  ## about 10.9 min at 25 Hz, 8.5 min at 1 Hz
  DT <- N / (Rate * 60)
  df <- Rate / N
  ## as needed for the Fourier transform (fft()):
  frq <- c(seq(0, Rate / 2, by = df), seq(-Rate / 2 + df, -df, by = df))
  N2 <- N %/% 2
  fmax <- 2
  nlim <- which(frq > fmax)[1]
  ## Proceed sequentially through the time series:
  irw <- nrow(DS)
  ir <- seq(1, irw - N, by = N2)
  i1 <- N2 %/% 2 + 1
  i2 <- i1 + N2 - 1
  rhozero <- 1013.25 * 100 / (287.05 * 288.15)
  rPar <- responsePar
  for (i in ir) {
    DSW <- DS[i:(i+N-1), ]
    f <- fft (DSW[, RV])
    ## get the mean values of air density and Mach number, for tau1 adjustment:
    MRHO <- MachNumber(DSW$PSXC, DSW$QCXC) * 
      DSW$PSXC * 100 / (287.05 * (273.15 + DSW$ATX)) / rhozero
    rPar$tau1 <- responsePar$tau1 * (0.3 / mean(MRHO, na.rm = TRUE)) ^ 0.68
    ## Modify the spectrum by the inverse of the response function:
    AFFT <- LTphase(frq, rPar)
    AFFT$frq <- frq
    AFFT$Phase <- AFFT$Phase * pi / 180
    fp <- fft(DSW[, RV])
    H <- complex (modulus = AFFT$Amp, argument = AFFT$Phase)
    xn <- Re(fft(fp / H, inverse = TRUE)) / N
    FFT <- xn
    if (i == 1) {
      DS[1:i2, RVC] <- FFT[1:i2]
    } else if (i == ir[length(ir)]) {
      DS[(i+i1-1):(i+N-1), RVC] <- FFT[i1:N]
    } else {
      ## save only the central part
      DS[(i + (i1:i2) - 1), RVC] <- FFT[i1:i2]
    }
  }
  rvc <- c(DD[1:(j1-1), RVC], DS[, RVC], DD[(j2+1):nrow(DD), RVC])
  return(rvc)
}

processFile <- function(ProjectDir, Project, Flight) {
  ## Find the available air_temperature variables:
  fname <- file.path(DataDirectory(), sprintf('%s/%s%s.nc',
                                              ProjectDir, Project, Flight))
  FI <- DataFileInfo(fname, LLrange = FALSE)
  Rate <- FI$Rate
  TVARS <- unlist(FI$Measurands$air_temperature)
  TVARS <- TVARS[-which ('ATX' == TVARS)]  # don't include ATX, AT_A, AT_A2
  if ('AT_A' %in% TVARS) {TVARS <- TVARS[-which('AT_A' == TVARS)]}
  if ('AT_A2' %in% TVARS) {TVARS <- TVARS[-which('AT_A2' == TVARS)]}
  if ('AT_VXL' %in% TVARS) {TVARS <- TVARS[-which('AT_VXL' == TVARS)]}
  if ('AT_VXL2' %in% TVARS) {TVARS <- TVARS[-which('AT_VXL2' == TVARS)]}
  if ('OAT' %in% TVARS) {TVARS <- TVARS[-which('OAT' == TVARS)]} #omit Ophir T
  if (!UH1 && (Rate != 25)) {  ## omit unheated-probe measurements
    for (i in length(TVARS):1) {  ## note reversed order
      if(grepl('nheated', FI$LongNames[which(TVARS[i] == FI$Variables)])) {
        TVARS <- TVARS[-i]
      }
      ## Some older projects (e.g., FRAPPE) have ATRx for the Rosemount
      if(grepl('^ATR', TVARS[i])) {
        TVARS <- TVARS[-i]
      }
    }
  }
  RVARS <- sub('^A', 'R', TVARS)
  ## get the old netCDF variables needed to calculate the modified variables
  VarList <- standardVariables(TVARS)
  ## Treat case, in some old projects, where EWX and MACHX are not present:
  if (!('EWX' %in% FI$Variables)) {
    VarList <- VarList[-which('EWX' == VarList)]
    VarList <- c(VarList, 'EDPC')
  }
  if (!('MACHX' %in% FI$Variables)) {
    VarList <- VarList[-which('MACHX' == VarList)]
    VarList <- c(VarList, 'XMACH2')
  }
  if (!('GGALT' %in% FI$Variables)) {
    VarList <- VarList[-which('GGALT' == VarList)]
    VarList <- c(VarList, 'GGALT_NTL')
  }
  ## Add the recovery temperatures if present; otherwise recalculate:
  RVS <- RVARS[RVARS %in% FI$Variables]
  ## The next lines deal with, for a 25-Hz file, the recovery temperature
  ## is present only at 1-Hz rate (e.g., WECANrf08h.nc)
  if (length(RVS) > 0) {
    if (Rate == 25) {
      St <- as.integer(gsub(':', '', sub('.* ', '', FI$Start)))
      Ed <- as.integer(gsub(':', '', sub('.* ', '', FI$Start + 100)))
      D <- getNetCDF(fname, RVS, St, Ed)
      for (RV in RVS) {  ## only add if present at 25 Hz in netCDF file
        if (attr(D[, RV], 'Dimensions')[[1]] == 'sps25') {
          VarList <- c(VarList, RV)
        }
      }
    } else {
      VarList <- c(VarList, RVS)
    }
  }
  D <- getNetCDF (fname, VarList)
  if ('EDPC' %in% VarList) {
    D$EWX <- D$EDPC
  }
  if ('XMACH2' %in% VarList) {
    D$MACHX <- sqrt(D$XMACH2)
  }
  if ('GGALT_NTL' %in% VarList) {
    D$GGALT <- D$GGALT_NTL
  }
  Rate <- attr(D, 'Rate')

  ## Calculate the new variables:
  E <- SmoothInterp(D$EWX / D$PSXC, .Length = 0)
  D$Cp <- SpecificHeats (E)[, 1]
  D$DH <- D$TASX^2 / (2 * D$Cp)
  ## Recalculate recovery temperatures if missing:
  retrievedRVARS <- rep(FALSE, length(RVARS))
  for (RV in RVARS) {
    if(!(RV %in% VarList)) {
      TV <- sub('^R', 'A', RV)
      prb <- 'HARCO'
      if (TV == 'ATH2') {prb <- 'HARCOB'}
      if (grepl('ATF', TV) || grepl('ATR', TV)) {prb <- 'UNHEATED'}
      retrievedRVARS[which(RV == RVARS)] <- TRUE
      D[, RV] <- SmoothInterp(D[, TV] + RecoveryFactor(D$MACHX, prb) * D$DH, .Length = 0)
      ## modify the attributes for saving in the new file:
      attr(D[, RV], 'long_name') <- paste0("Recovery Air Temperature, ", prb)
      attr(D[, RV], 'standard_name') <- NULL
      attr(D[, RV], 'actual_range') <- NULL
      attr(D[, RV], 'Dependencies') <- NULL
      attr(D[, RV], 'measurand') <- "recovery_temperature"
      attr(D[, RV], 'label') <- paste0("recovery temperature (", RV, ") [deg C]")
      attr(D[, RV], 'RecoveryFactor') <- NULL
      attr(D[, RV], 'DataQuality') <- "Reconstructed"
    }
  }
  platform <- attr(D, 'Platform')
  platform <- ifelse (grepl('677', platform), 'GV', 'C130')
  # filter DH:
  CutoffPeriod <- rep(Rate * 1.0, length(TVARS)) # Standard is 1 s for DH filtering
  probe <- rep('HARCO', length(TVARS))  # used to determine the recovery factor
  PAR <- list()
  if (length(TVARS) > 0) {
    for (i in 1:length(TVARS)) {
      PAR[[length(PAR) + 1]] <- ParamSH
    }
  }
  # check for ATF or ATR
  ic <- which(grepl('ATF', TVARS))
  if (length(ic) > 0) {
    CutoffPeriod[ic] <- Rate * 0.5  # 0.5 s for ATFx
    probe[ic] <- 'UNHEATED'
    if (platform == 'GV') {
      PAR[[ic]] <- ParamSF
    } else {
      PAR[[ic]] <- Param1
    }
  }
  ic <- which(grepl('ATR', TVARS))  ## old naming convention (e.g., FRAPPE)
  if (length(ic) > 0) {
    CutoffPeriod[ic] <- Rate * 0.5  # 0.5 s for ATFx
    probe[ic] <- 'UNHEATED'
    for (icc in ic) {
      if (platform == 'GV') {
        PAR[[icc]] <- ParamSF
      } else {
        PAR[[icc]] <- Param1
      }
    }
  }
  ic <- which(grepl('ATH.*2', TVARS))  ## .* allows for names like ATHL2
  if (length(ic) > 0) {
    probe[ic] <- 'HARCOB'
    PAR[[ic]] <- ParamSH
  }
  if (Rate == 1) {  # Protection against script failure for a LRT file
    CutoffPeriod[CutoffPeriod == 1] <- 2.2
    CutoffPeriod[CutoffPeriod == 0.5] <- 2.2
  }
  
  DHM <- rep(D$DH, length(TVARS))
  dim(DHM) <- c(length(D$DH), length(TVARS))
  newTVARS <- paste0(TVARS, 'C')
  newRVARS <- paste0(RVARS, 'C')
  newRVARS2 <- paste0(RVARS, 'C2')
  newTVARS2 <- paste0(TVARS, 'C2')
  filteredTVARS <- paste0(TVARS, 'F')
  for (i in 1:length(TVARS)) {
    D[, filteredTVARS[i]] <- correctDynamicHeating(D, as.character(TVARS[i]))
  }
  
  ## Calculate the corrected temperatures:
  for (i in 1:length(TVARS)) {
    D[, newTVARS[i]] <- correctTemperature(D, filteredTVARS[i], responsePar = PAR[[i]])
    D[, newRVARS[i]] <- correctTemperature(D, RVARS[i], responsePar = PAR[[i]])
    if (FFT) {
      D[, newRVARS2[i]] <- addFFTsoln(D, RVARS[i], responsePar = PAR[[i]])
    }
    # get the recovery factor from the attribute:
    rf <- recoveryF(D, TVARS[[i]])
    # Is the associated recovery temperature present?
    dep <- attr(D[, TVARS[i]], 'Dependencies')
    RT <- gsub(' .*', '', gsub('^. ', '', dep))
    # Get the dynamic-heating correction:
    if ('EWX' %in% names(D)) {
      ER <- SmoothInterp(D$EWX / D$PSXC, .Length = 0)
      CP <- SpecificHeats(ER)
    } else {
      CP <- SpecificHeats()
    }
    if (RT %in% VarList) {
      X <- rf * D$MACHX^2 * CP[, 3] / (2 * CP[, 2])
      Q <- (273.15 + D[, RT]) * X / (1 + X)
    } else {
      Q <- rf * D$TASX^2 / 2010
    }
    # Interpolate to avoid missing values
    Q <- SmoothInterp(Q, .Length = 0)
    ## This should be the same as newTVARS except for numerical effects
    D[, newTVARS2[i]] <- D[, newRVARS[i]] - Q  ## using diff.eq.
    if (FFT) {
      D[, newTVARS2[i]] <- D[, newRVARS2[i]] - Q  ## using FFT
    }
  }
  ## At this point, have all the new variables.
  ## If this is a 25-Hz file, the 1-Hz averages are better than those
  ## obtained by using this script at 1 Hz. The 1-Hz averages can be used 
  ## instead in the corrected 1-Hz file if they are available. Those averages
  ## can be obtained from the HRT file by averaging like this:
  if (Rate == 25) {
    ## Construct a new data.frame containing the 1-Hz averages.
    ## (By convention, 1-Hz values are the average of values beginning
    ## at the indicated time, so Time+0.5s or Var[ix+13] represents the
    ## value at Time[ix].)
    D1HZ <- data.frame(Time = D$Time[seq(1, nrow(D)-1, by = 25)])
    for (TV1 in c(newTVARS, newRVARS, filteredTVARS)) {
      vs <- zoo::na.approx (as.vector(D[, TV1]), na.rm=FALSE, rule = 2)
      rmean <- rollapply(vs, width=25, FUN=mean, fill=NA)
      rmean <- zoo::na.approx (as.vector(rmean), na.rm=FALSE)
      means <- mean(vs, na.rm=TRUE)
      vs[is.na(vs)] <- rmean[is.na(vs)]
      vs[is.na(vs)] <- means
      vs <- signal::filtfilt (signal::butter (3, 2/50), vs)
      D1HZ[, TV1] <- vs[seq(13, nrow(D), by = 25)]
    }
    save(D1HZ, file = sprintf('%s%s1Hztemp.Rdata', Project, Flight))
    ## also save as netCDF file:
    # D1HZ <- transferAttributes(D, D1HZ)
    # attr(D1HZ, 'Rate') <- 1
    # for (v in names(D1HZ)) {
    #   if (v == 'Time') {next}
    #   vv <- sub('.$', '', v)
    #   la <- attributes(D[, vv])
    #   la$Dimensions <- list('Time')
    #   attr(D1HZ[, v], 'Dimensions') <- list('Time')
    # }
    # makeNetCDF(D1HZ, sprintf('%s%sLRT.nc', Project, Flight))
  } else { ## This section uses the data.frame saved previously if available
    fchk <- sub('1Hz', 'h1Hz', sprintf('%s%s1Hztemp.Rdata', Project, Flight))
    if(file.exists(fchk)) {
      load(fchk)  ## loads D1HZ
      ## beware of time mis-match: get the indices where Time matches
      ix <- match(D$Time, D1HZ$Time)  ## will be NA for no match
      for (v in names(D1HZ)) {
        if (v == 'Time') {next}
        if (v %in% names(D)) {
          D[!is.na(ix), v] <- D1HZ[!is.na(ix), v] ## replace by saved value from 25-Hz run
        }
      }

    }
  }
  
  
  
  
  ## Add the new variables:
  if (FFT) {
    newTVARS <- c(newTVARS, filteredTVARS, newTVARS2)  ## omit newTVARS2?
  } else {
    newTVARS <- c(newTVARS, filteredTVARS) 
  }
  
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
  DATT <- D  ## Save in case this is needed
  
  ## variables to add to the netCDF file:
  VarNew <- newTVARS
  VarOld <- c(TVARS, TVARS, TVARS)
  VarUnits <- rep("deg_C", length(newTVARS))
  DataQ <- rep("corrected_for_sensor_response", length(newTVARS))
  VarLongName <- c(rep("Ambient Temperature, corrected", length(TVARS)), 
                   rep("Ambient Temperature, Q filtered", length(TVARS)),
                   rep("Ambient Temperature, FFT corrected", length(TVARS)))
  VarStdName <- rep("air_temperature", length(newTVARS))
  Dependencies <- rep('', length(newTVARS))
  for (i in 1:(length(newTVARS))) {
    Dependencies[i] <- sprintf('2 %s TASX', VarOld[i])
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
    ncatt_put (newfile, VarNew[i], attname="long_name",
               attval=VarLongName[i])
    ncatt_put (newfile, VarNew[i], attname="Dependencies",
               attval=Dependencies[i])
    ncatt_put (newfile, VarNew[i], attname="DataQuality",
               attval=DataQ[i])
    ## add the measurements for the new variable
    if (Rate == 1) {
      ncvar_put (newfile, varCDF[[i]], D[, VarNew[i]])
    } else if (Rate == 25) {
      ncvar_put (newfile, varCDF[[i]], D[, VarNew[i]],
                 count=c(25, nrow(D)/25))
    }
  }
  
  ## Add the corrected (newRVARS), FFT-corrected (newRVARS2) and 
  ## reconstructed (RVARS) recovery temperatures. Only add the 
  ## reconstructed recovery temperature if it was not already present.
  VarUnits <- "deg_C"
  DataQ <- "corrected_for_sensor_response"
  VarLongName <- c(rep("Recovery Temperature, corrected", length(newRVARS)),
                   rep("Recovery Temperature, FFT corrected", length(newRVARS2)))
  VarStdName <- "air_temperature"
  if (FFT) {
    newRVARS <- c(newRVARS, newRVARS2)
  }
  Dependencies <- c(rep('', length(newRVARS)))
  for (i in 1:length(newRVARS)) {
    Dependencies[i] <- sprintf('2 %s TASX', RVARS[i])
  }
  for (i in 1:length(newRVARS)) {
    ## create the new variable and add it to the netCDF file
    varCDF[[i + length(VarNew)]] <- ncvar_def (newRVARS[i],
                                               units=VarUnits,
                                               dim=Dim,
                                               missval=as.single(-32767.), prec="float",
                                               longname=VarLongName)
    newfile <- ncvar_add (newfile, varCDF[[i + length(VarNew)]])
    ## if the original variable is present, transfer attributes from the old 
    ## variable, then add new ones
    if (RVARS[i] %in% VarList) {
      ATV <- ncatt_get (netCDFfile, RVARS[i])
      copy_attributes (ATV, newRVARS[i], newfile)
    }
    ncatt_put (newfile, newRVARS[i], attname="standard_name",
               attval=VarStdName)
    ncatt_put (newfile, newRVARS[i], attname="long_name",
               attval=VarLongName[i])
    ncatt_put (newfile, newRVARS[i], attname="Dependencies",
               attval=Dependencies[i])
    ncatt_put (newfile, newRVARS[i], attname="DataQuality",
               attval=DataQ)
    ## add the measurements for the new variable
    if (Rate == 1) {
      ncvar_put (newfile, varCDF[[i+length(VarNew)]], D[, newRVARS[i]])
    } else if (Rate == 25) {
      ncvar_put (newfile, varCDF[[i+length(VarNew)]], D[, newRVARS[i]],
                 count=c(25, nrow(D)/25))
    }
  }
  ## Add reconstructed RVARS here if desired.
  if (RTN) {
    ## Special treatment of recovery temperatures that are present in the 25-Hz
    ## netCDF file only at 1 Hz, like RTH1 for WECAN
    ## -- the needed 25-Hz  measurement of RTH1 is recovered from ATH1 
    ##    (present at 25 Hz) and Q. Then, before writing the new
    ##    variable (if requested), the old variable is renamed to 'RTH1X' in
    ##    the netCDF file.  
    if (Rate == 25) {
      for (RV in RVARS) {
        if (attr(D[, RV], 'DataQuality') == 'Reconstructed') {
          if (RV %in% FI$Variables) {  ## avoid duplicate variable names
            newfile <- ncvar_rename(newfile, RV, paste0(RV, 'x'))
            print (sprintf ('renaming old variable %s to %s', RV, paste0(RV, 'x')))
          }
        }
      }
    }
    # if (attr(D[, RV], 'Dimensions')[[1]] == 'sps25') {
    if (Project == 'WECAN') {
      
    }
    VarUnits <- "deg_C"
    DataQ <- "reconstructed"
    VarLongName <- "Recovery Temperature, reconstructed"
    VarStdName <- "air_temperature"
    Dependencies <- c(rep('', length(RVARS)))
    for (i in 1:length(RVARS)) {
      Dependencies[i] <- sprintf('2 %s TASX', TVARS[i])
    }
    L <- length(VarNew) + length(newRVARS)
    for (i in 1:length(RVARS)) {
      ## create the new variable and add it to the netCDF file
      if ((RVARS[i] %in% VarList)) {next}
      varCDF[[i + L]] <- ncvar_def (RVARS[i],
                                    units=VarUnits,
                                    dim=Dim,
                                    missval=as.single(-32767.), prec="float",
                                    longname=VarLongName)
      newfile <- ncvar_add (newfile, varCDF[[i + L]])
      ncatt_put (newfile, RVARS[i], attname="standard_name",
                 attval=VarStdName)
      ncatt_put (newfile, RVARS[i], attname="long_name",
                 attval=VarLongName)
      ncatt_put (newfile, RVARS[i], attname="Dependencies",
                 attval=Dependencies[i])
      ncatt_put (newfile, RVARS[i], attname="DataQuality",
                 attval=DataQ)
      ## add the measurements for the new variable
      if (Rate == 1) {
        ncvar_put (newfile, varCDF[[i+L]], D[, RVARS[i]])
      } else if (Rate == 25) {
        ncvar_put (newfile, varCDF[[i+L]], D[, RVARS[i]],
                   count=c(25, nrow(D)/25))
      }
    }
  }
  
  ## then close to write the new file
  nc_close (newfile)
}

if (ALL) {
  Fl <- sort (list.files (sprintf ("%s%s/", DataDirectory (), ProjectDir),
                          sprintf ("%srf..h.nc", Project)), decreasing = FALSE)
  
} else if (ALLLR) {
  Fl <- sort (list.files (sprintf ("%s%s/", DataDirectory (), ProjectDir),
                          sprintf ("%srf...nc", Project)), decreasing = FALSE)
}
if (ALL || ALLLR) {
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

