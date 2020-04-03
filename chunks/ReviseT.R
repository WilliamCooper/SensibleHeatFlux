

print (sprintf ('ReviseT Processor is starting -- %s', Sys.time()))
require(Ranadu, quietly = TRUE, warn.conflicts = FALSE)
options(stringsAsFactors = FALSE)
setwd ('~/RStudio/SensibleHeatFlux')

##----------------------------------------------------------
## These are the run options; reset via command-line or UI:
GeneratePlots <- TRUE
ShowPlots <- FALSE  ## leave FALSE for batch runs or shiny-app runs
ALL <- FALSE
NEXT <- FALSE
StartTime <- 0
EndTime <- 400000
Project <- 'CSET'
ProjectDir <- Project
Flight <- 1
UpdateRT <- TRUE
FilterQ <- TRUE
FilterR <- FALSE
FilterN <- FALSE
thisFileName <- 'ReviseT'

Directory <- DataDirectory ()
getNext <- function(ProjectDir, Project) {
  Fl <-
    sort (list.files (
      sprintf ("%s%s/", DataDirectory (), ProjectDir),
      sprintf ("%srf..TC.nc", Project)
    ), decreasing = TRUE)[1]
  if (is.na (Fl)) {
    Flight <- 1
  } else {
    Flight <- sub (".*rf", '',  sub ("TC.nc", '', Fl))
    Flight <- as.numeric(Flight) + 1
  }
  return (Flight)
} 
if (!interactive()) {
  ## can run interactively or via Rscript; these are Rscript arguments
  run.args <- commandArgs (TRUE)
  if (length (run.args) > 0) {
    if (nchar(run.args[1]) > 1) {
      Project <- run.args[1]
      ProjectDir <- Project
    }
  } else {
    print ("Usage: Rscript ReviseT.R Project Flight StartTime EndTime UpdateRT[Y/n] FilterQ[Y/n] FilterR[N/y] FilterN[N/y]")
    print (" All must be provided in order, and defaults will be used beyond the last provided.")
    print (" UpdateRT integrates the diff. eqs. to estimate the corrected recovery temperature.")
    print (" FilterQ uses a numerical filter to update the dynamic-heating correction.")
    print ("    -- if 'n' numerical integration is used instead, usually with inferior result.")
    print ("FilterR tries to remove line-resonance effects from static and dynamic pressure.")
    print ("FilterN adds a noise filter to the pressure measurements to reduce noise.")
    print ("Example: Rscript ReviseT.R CSET 1")
    print ("    -- will use defaults for remaining arguments")
    stop("exiting... Run again with arguments.")
  }
  ## Flight
  if (length (run.args) > 1) {
    if (run.args[2] != 'NEXT') {
      Flight <- as.numeric (run.args[2])
    } else {
      ## Find max rf in data directory,
      ## Use as default if none supplied via command line:
      Flight <- getNext(ProjectDir, Project)
    }
  }
  if (length (run.args) > 2) {
    StartTime <- as.integer(run.args[3])
    if (StartTime < 0) {
      StartTime <- 0
    }
  }
  if (length (run.args) > 3) {
    EndTime <- as.integer(run.args[4])
    if (EndTime <= 0) {
      EndTime <- 400000
    }
  }
  if (length (run.args) > 4) {
    UpdateRT <- run.args[5]
    if (UpdateRT == 'y' || UpdateRT == 'Y' || UpdateRT == 'TRUE') {
      UpdateRT <- TRUE
    } else {
      UpdateRT <- FALSE
    }
  }
  if (length (run.args) > 5) {
    FilterQ <- run.args[6]
    if (FilterQ == 'y' || FilterQ == 'Y' || FilterQ == 'TRUE') {
      FilterQ <- TRUE
    } else {
      FilterQ <- FALSE
    }
  }
  if (length (run.args) > 6) {
    FilterR <- run.args[7]
    if (FilterR == 'y' || FilterR == 'Y' || FilterR == 'TRUE') {
      FilterR <- TRUE
    } else {
      FilterR <- FALSE
    }
  }
  if (length (run.args) > 7) {
    FilterN <- run.args[8]
    if (FilterN == 'y' || FilterN == 'Y' || FilterN == 'TRUE') {
      FilterN <- TRUE
    } else {
      FilterN <- FALSE
    }
  }
  } else {
  ## This is the interactive section
  x <-
    readline (sprintf ("Project is %s; CR to accept or enter new project name: ", Project))
  if (nchar(x) > 1) {
    Project <- x
    if (grepl ('HIPPO', Project)) {
      ProjectDir <- 'HIPPO'
    } else {
      ProjectDir <- Project
    }
  }
  x <-
    readline (
      sprintf (
        "Flight is %d; CR to accept, number or 'ALL' or 'NEXT' for new flight name: ",
        Flight
      )
    )
  if (x == 'ALL') {
    ALL <- TRUE
  } else if (x == 'NEXT') {
    Flight <- getNext(ProjectDir, Project)
  } else if (nchar(x) > 0 && !is.na(as.numeric(x))) {
    Flight <- as.numeric(x)
  }
  x <-
    readline (
      sprintf (
        "Start/End Times are %d %d; CR to accept, HHMMSS HHMMSS to change: ",
        StartTime,
        EndTime
      )
    )
  if (nchar(x) > 0) {
    StartTime <- as.integer(sub('[, ] *[0-9]* *', '', x))
    EndTime  <-  as.integer(sub('[0-9]*[, ]*', '', x))
  }
  x <- 
    readline (sprintf (
      "UpdateRT (recovery T) is %s; CR to accept, Y or T to enable, N or F to disable: ",
      UpdateRT
    ))
  if (nchar(x) > 0)
    UpdateRT <- ifelse ((grepl('^T', x) || grepl('^Y', x)), TRUE, FALSE)
  x <- 
    readline (sprintf (
      "FilterQ (dynamic heating) is %s; CR to accept, Y or T to enable, N or F to disable: ",
      FilterQ
    ))
  if (nchar(x) > 0)
    FilterQ <- ifelse ((grepl('^T', x) || grepl('^Y', x)), TRUE, FALSE)
  x <- 
    readline (sprintf (
      "FilterR (line resonance) is %s; CR to accept, Y or T to enable, N or F to disable: ",
      FilterR
    ))
  if (nchar(x) > 0)
    FilterR <- ifelse ((grepl('^T', x) || grepl('^Y', x)), TRUE, FALSE)
  x <- 
    readline (sprintf (
      "FilterN (ambient P noise) is %s; CR to accept, Y or T to enable, N or F to disable: ",
      FilterN
    ))
  if (nchar(x) > 0)
    FilterN <- ifelse ((grepl('^T', x) || grepl('^Y', x)), TRUE, FALSE)
}
print (
  sprintf (
    'run controls:  Project: %s;  Flight: %dTime Interval %d--%d', 
    Project,
    Flight,
    StartTime,
    EndTime))
print (
  sprintf (
    '  UpdateRT=%s, FilterQ=%s, FilterR=%s, FilterN=%s',
    UpdateRT,
    FilterQ,
    FilterR,
    FilterN 
  )
)

fname = sprintf("%s%s/%srf%02d.nc", Directory, ProjectDir, Project, Flight)

SaveRData2 <-
  sprintf("%s2.Rdata", thisFileName)  # In case these get saved...
SaveRData <- sprintf("%s.Rdata", thisFileName)
print (sprintf ('reading netCDF file %s', fname))
FI <- DataFileInfo(fname, LLrange = FALSE)
# Get the measurements of air_temperature:
ATS <- FI$Measurands$air_temperature
ATS <- ATS[-which('ATX' == ATS)]  # skip ATX; could update later
RTS <- sub('^AT', 'RT', ATS)
Probe <- rep('HARCO', length(ATS))
# Do the corresponding recovery temperatures exist? (if not, reconstruct them)
# This is useful because recovery temperature is sometimes omitted from archives
Vars <- c(ATS, RTS, 'TASX', 'MACHX', 'PSXC', 'QCXC', 'ATX', 'EWX')
Vars <-
  Vars[Vars %in% FI$Variables] # Only load those that are available
## Create a subset working file, overwriting if necessary:
fnameW <- sub('.nc', 'TEMP.nc', fname)
if (StartTime > 0 || EndTime < 400000) {
  suppressMessages(ncsubset(fname, fnameW, StartTime, EndTime))
} else {
  file.copy(fname, fnameW, overwrite = TRUE)
}
D <- getNetCDF(fnameW, Vars)
## Deal with NAs:
DS <- D
for (n in names(D)) {
  if (n == 'Time') {
    next
  }
  D[, n] <- SmoothInterp(D[, n], .maxGap = 5000, .Length = 0)
}
D <- transferAttributes(DS, D)
if (!('MACHX' %in% names(D))) {
  if ('EWX' %in% names(D)) {
    D$MACHX <- MachNumber(D$PSXC, D$QCXC, D$EWX)
  } else {
    D$MACHX <- MachNumber(D$PSXC, D$QCXC)
  }
}
for (AT in ATS) {
  if (grepl('Unheated', attr(D[, AT], 'long_name'))) {
    if (grepl('130', FI$Platform)) {
      Probe[which(AT == ATS)] <- 'FC'
    } else {
      Probe[which(AT == ATS)] <- 'F'
    }
  } else if (AT == 'ATH2') {
    Probe[which(AT == ATS)] <- 'H2'
  } else {
    Probe[which(AT == ATS)] <- 'H'
  }
}
print ('These are the probe types that will be processed:')
print (Probe)
## If some recovery temperatures are not present, reconstruct them:
for (i in 1:length(RTS)) {
  if (RTS[i] %in% names(D)) {
    print (sprintf ('%s is available', RTS[i]))
  } else {
    print (sprintf ('%s is not available; reconstructing it...', RTS[i]))
    D[, RTS[i]] <-
      D[, ATS[i]] + RecoveryFactor(D$MACHX, Probe[i]) * D$TASX ^ 2 / 2010
  }
}

## If the data file is 1 Hz, expand by interpolation
## (Tried this; replaced by adaptive-step procedure with rkck.integrate()
# Rate <- attr(D, 'Rate')
# if (Rate > 1) {
#   Xpand <- 1
#   DX <- D
# } else {
#   Xpand <- 11  ## Use odd number for later realignment
#   itime <- as.integer(D$Time-D$Time[1])
#   DX <- data.frame(Time=D$Time[1]+(0:(Xpand*length(itime)-1))/Xpand)
#   for (nm in names(D)) {
#     if (nm == 'Time') {next}
#     DX <- cbind(DX, data.frame(SmoothInterp(interp1(0:length(itime), D[, nm],
#           (0:(Xpand*length(itime)-1))/Xpand, method='linear'), .Length=Xpand)))
#   }
#   names(DX) <- names(D)
#   attr(DX, 'Rate') <- Xpand
# }
## Ready to calculate revised temperature. Define two functions:
# function to provide response parameters at Z=0.6 and adjust for input Z:
getParam <- function (S = 'H', Z = 0.3) {
  P <- switch(
    S,
    'F' = data.frame (a = 0.733, tau1 = 0.0308, tau2 = 0.45),
    'FC' = data.frame(a = 0.651, tau1 = 0.0295, tau2 = 1.04),
    'H' = data.frame (a = 0.051, tau1 = 0.23, tau2 = 1.04),
    'H2' = data.frame(a = 0.063, tau1 = 0.30, tau2 = 0.99),
    'S' = data.frame (a = 0.27,  tau1 = 0.74, tau2 = 0.74),
    data.frame (a = 0.5,   tau1 = 0.79, tau2 = 0.79) # default
  )
  tau1 <- P$tau1 * (0.3 / Z) ^ 0.6
  tau2 <- P$tau2 * (0.3 / Z) ^ 0.6
  a <- rep(P$a, length(tau1))
  return(data.frame(a = a, tau1 = tau1, tau2 = tau2))
}
## Arguments RTS, ATS, and SensorTypeS can all be vectors, to process all.
ReviseT <-
  function(.data,
           RTS = 'RTX',
           ## Recovery Ts
           ATS = 'ATH1',
           ## Air Ts
           SensorTypeS = 'H',
           UpdateRT = FALSE,
           FilterQ = TRUE) {
    Rate <- attr(.data, 'Rate')
    if (is.null(Rate)) {
      print ('Rate attribute missing so using 1-Hz processing...')
      Rate <- 1
    }
    DEX <-
      .data  ## to allow reuse of some code previously written with DEX
    load(file = 'AR.Rdata')  ## Retrieve the filter coefficients AR, ARH, ARG, AR1, ARH1, ARG1
    ## and Lshift, LshiftH, Lshift1, LshiftH1
    rhozero <-
      100 * 1013.25 / (287.05 * 288.15)  ## Reference conditions for air density
    if ('MACHX' %in% names(DEX)) {
      Mach <- DEX$MACHX
    } else {
      if ('EWX' %in% names(DEX)) {
        Mach <- MachNumber(DEX$PSXC, DEX$QCXC, DEX$EWX)
      } else {
        Mach <- MachNumber(DEX$PSXC, DEX$QCXC)
      }
    }
    Z <-
      Mach * 100 * DEX$PSXC / (287.05 * (273.15 + DEX$ATX)) / rhozero
    PRS <- rep('HARCO', length(SensorTypeS))
    rfs <- rep(0.95, length(SensorTypeS))
    for (SensorType in SensorTypeS) {
      PR <- switch(
        SensorType,
        'F' = 'UNHEATED',
        'FC' = 'UNHEATED',
        'S' = 'Rose',
        'H2' = 'HARCOB',
        'H' = 'HARCO'
      )
      itmp <- which(SensorType == SensorTypeS)[1]
      PRS[itmp] <- PR
      rf <- RecoveryFactor(Mach, PR)
      XXA <- rf * Mach ^ 2 / 5
      P <- getParam(SensorType, Z)
      a <- P$a
      tau1 <- P$tau1
      tau2 <- P$tau2
      ## Get the filter coefficients:
      ARC <- switch(
        SensorType,
        'F' = AR,
        'FC' = ARG,
        'H' = ARH,
        'H2' = ARHB,
        AR
      )
      Ls <- switch(
        SensorType,
        'F' = Lshift * 40,
        'FC' = Lshift * 40,
        'H' = LshiftH * 40,
        'H2' = LshiftH * 40,
        Lshift * 40
      )
      if (Rate == 1) {
        ARC <- switch(
          SensorType,
          'F' = AR1,
          'FC' = ARG1,
          'H' = ARH1,
          'H2' = ARHB1,
          AR1
        )
        Ls <- switch(
          SensorType,
          'F' = Lshift1 * 1000,
          'FC' = Lshift1 * 1000,
          'H' = LshiftH1 * 1000,
          'H2' = LshiftH1 * 1000,
          Lshift1 * 1000
        )
      }
      # DEX$RT <- DEX[, RT] # Now RT may be a vector
      NRW <- nrow(DEX)
      if (UpdateRT) {
        ## Update RT to correct for time response
        DEX$DTMDT <- c(0, diff(DEX[, RTS[itmp]], 2), 0) * Rate / 2
        ## Treat HARCOs differently:
        #         if (grepl('HARCO', PR)) {
        #           ## HARCO requires different treatment because a=0:
        #           DF1 <- diff(DEX[, RTS[itmp]])
        #           D2DTDT2 <- (c(DF1, 0) - c(0, DF1)) * Rate ^ 2
        #           ## The correction is very noisy before takeoff, so suppress
        #           ic <- DEX$TASX > 50
        #           RTC <- DEX[, RTS[itmp]]
        #           RTC[ic] <-
        #             (RTC + (tau1 + tau2) * DEX$DTMDT + tau1 * tau2 * D2DTDT2)[ic]
        #         } else {
        DEX$Ts <- DEX[, RTS[itmp]]  # just define a storage vector
        fS <- function(y, t) {
          # Eq. Ts3
          i <- min (as.integer(t), NRW)
          f <- t - i
          i1 <- min(i + 1, NRW)
          DTMDT <- (1 - f) * DEX$DTMDT[i] + f * DEX$DTMDT[i1]
          RT <- (1 - f) * DEX[i, RTS[itmp]] + f * DEX[i1, RTS[itmp]]
          ((1 / a[i]) * (tau1[i] * DTMDT + RT - (1 - a[i]) * y) - y) /
            (Rate * tau2[i])
        }
        DEX$Ts <- rk4.integrate (fS, DEX$Ts[1], 1:nrow(DEX))
        RTC <-
          (1 / a) * (tau1 * DEX$DTMDT + DEX[, RTS[itmp]] - (1 - a) * DEX$Ts)
        DEX$RTX <- RTC
      } else {
        DEX$RTX <- DEX[, RTS[itmp]]
      }
      DEX$Q <- (DEX$RTX + 273.15) * XXA / (1 + XXA)
      DEX$TsQ <- DEX$RTX
      DEX$Qp <- DEX$Q
      # fSQ <- function(y, i) {
      #   (DEX$Q[i] - y) / (Rate * tau2[i])
      # }
      fSQ <- function(y, t) {
        i <- min (as.integer(t), NRW)
        f <- t - i
        i1 <- min(i + 1, NRW)
        ((1 - f) * DEX$Q[i] + f * DEX$Q[i1] - y) / (Rate * tau2[i])
      }
      # fM <- function (y, i) {
      #   (a[i] * DEX$Q[i] + (1 - a[i]) * DEX$TsQ[i] - y) / (Rate * tau1[i])
      # }
      fM <- function (y, t) {
        i <- min (as.integer(t), NRW)
        f <- t - i
        i1 <- min(i + 1, NRW)
        R <- (a[i] * ((1 - f) * DEX$Q[i] + f * DEX$Q[i1]) +
                (1 - a[i]) * ((1 - f) * DEX$TsQ[i] + f * DEX$TsQ[i1]) - y) /
          (Rate * tau1[i])
        # print (sprintf ('fM evaluation at t=%.1f, i=%d, f=%.2f, y=%.2f, return=%.2f', t, i, f, y, R))
        return(R)
      }
      if (FilterQ) {
        ## Now add the filtered version:
        DEX$QF <- as.vector(signal::filter(ARC, DEX$Q))
        DEX$QF <- ShiftInTime(DEX$QF, .shift = -Ls, .rate = Rate)  
        DEX$ATXC <- DEX$RTX - DEX$QF
      } else {
        DEX$TsQ <- rk4.integrate (fSQ, DEX$Q[1], 1:nrow(DEX))
        DEX$Qp <- rk4.integrate (fM, DEX$Q[1], tv = 1:nrow(DEX))
        DEX$ATXC <- DEX$RTX - DEX$Qp
      }
      
      nm <- names(DEX)
      ATC <- paste0(ATS[itmp], 'C')
      nm[which('ATXC' == nm)] <- ATC
      names(DEX) <- nm
    }
    return (DEX)
  }

DXC <-
  ReviseT(D, RTS, ATS, SensorTypeS = Probe, UpdateRT = UpdateRT, FilterQ = FilterQ)

## If the original rate is 1, return to that rate.
## (Comment re averaging: The original measurements at 1 Hz represent
##  averages over the 1-s beginning at the listed times. The interpolation
##  process here begins at the listed time and adds "Xpand" values spaced
##  to the next second. Therefore they begin at the middle of the listed
##  time and should be displaced to represent that time before being
##  averaged again to 1 Hz. That is the reason for the shift below before
##  averaging.)
# if (Xpand > 1) {
#   ## Use averaging:
#   ixp <- Xpand %/% 2
#   for (AT in ATS) {
#     if (AT %in% names(DXC)) {
#       ATC <- paste0(AT, 'C')
#       A <- c(rep(DXC[1, ATC], ixp), DXC[1:(nrow(DXC)-ixp), ATC])
#       D[, ATC] <- zoo::rollapply (A, Xpand, FUN=mean, by=Xpand)
#     }
#   }
#   A <- c(rep(DXC$Qp[1], ixp), DXC$Qp[1:(nrow(DXC)-ixp)])
#   D$Qp <- zoo::rollapply (DXC$Qp, Xpand, FUN=mean, by=Xpand)  ## Also save Qprime
# }

## ----create-new-netcdf, cache=FALSE--------------------------------------

print (sprintf("making new netCDF file -- %s", Sys.time()))
source ('chunks/create-new-netcdf.R')

## ----modify-new-netcdf, include=TRUE-------------------------------------

source ('chunks/modify-new-netcdf.R')
## modify-new-netcdf.R
print (sprintf ("ReviseT Processor is finished -- %s", Sys.time()))