
## SPECIAL FOR SHOWING TRANSFER-FUNCTION MOD
## Also, doesn't set scale_x so it can be set later
VSpecC <-
  function (.data,
            .Variable = NA,
            VLabel = NA,
            col = NA,
            type = 'spectrum',
            corrected = FALSE,
            method = NA,
            xlim = NA,
            ylim = NA,
            Par = NA,
            # c(0.001, 15), ylim=c(0.0001,1),
            spans = 49,
            ae = 0.2,
            smoothBins = 0,
            segLength = 512,
            poles = 50,
            resolution = 0.0001,
            showErrors = 0,
            WavelengthScale = TRUE,
            ADD = NA,
            add = NA,
            EDR = FALSE,
            WACtheme = NA) {
    if (!is.data.frame(.data)) {
      # See if the first argument can be split into a data.frame and a variable:
      X <- substitute(.data)
      if (is.call(X)) {
        V <- try(eval(X), silent = TRUE)
        if (grepl('Error', V[[1]])) {
          V <- eval(plyr::as.quoted(X))  # eval(X) for names()
        }
        if (is.character(V[1])) {
          
        } else {
          V <- plyr::as.quoted(X)
          if (is.symbol(V[[1]])) {
            V <- vapply(V, deparse, 'character')
          }
        }
      } else {
        V <- plyr::as.quoted(X)
        if (is.symbol(V[[1]])) {
          V <- vapply(V, deparse, 'character')
        }
      }
      # print (c('first argument evaluates to', V))
      # Extract data.frame:
      .data <- get(V[[1]])
      .Variable <- V[[2]]
    }
    # print(str(.data))
    if (is.data.frame(.data)) {
      ## must be true, or exit. Needs to contain Time and TASX
      ## in addition to .Variable
      nm <- names(.data)
      V <- try(is.na(.Variable), silent = TRUE)
      if (grepl('Error', V[[1]])) {
        X <- substitute(.Variable)
        if (is.call(X)) {
          V <- try(eval(X), silent = TRUE)
          if (grepl('Error', V[[1]])) {
            V <- eval(plyr::as.quoted(X))  # eval(X) for names()
          }
          if (is.character(V[1])) {
            
          } else {
            V <- plyr::as.quoted(X)
            if (is.symbol(V[[1]])) {
              V <- vapply(V, deparse, 'character')
            }
          }
        } else {
          V <- plyr::as.quoted(X)
          if (is.symbol(V[[1]])) {
            V <- vapply(V, deparse, 'character')
          }
        }
        .Variable <- V
      }
      if (is.na(.Variable[1])) {
        nm <- nm[-which('Time' == nm)]
        nm <- nm[-which('TASX' == nm)]
        .Variable <- nm
        if (length(.Variable) > 3) {
          .Variable <- .Variable[1:3]
        }
      }
    } else {
      print('VSpec ERROR: first argument is not a data.frame.')
      return (NA)
    }
    if (is.null(attr(.data, 'Rate'))) {
      print ('VSpec warning: Rate attribute missing from data.frame, so using Rate=1')
      Rate <- 1
    } else {
      Rate <- attr(.data, 'Rate')
    }
    if (is.na(xlim[1])) {
      if (Rate == 1) {
        xlim <- c(0.001, 1)
      } else {
        xlim <- c(0.001, 15)
      }
    }
    if (is.na(ylim[1])) {
      ylim <- c(1.e-4, 1.)
    }
    if (!is.na(method)) {
      type <- method
    }  ## method over-rides if present
    for (.V in .Variable) {
      if (.V %in% names(.data)) {
        NV <- which(.V == .Variable)
        Z <- capture.output (v <- detrend (.data[, c('Time', .V)]))
        if (!is.na(VLabel[1]) &&
            length(VLabel) >= NV) {
          ## use this alternate name in legend
          V <- VLabel[NV]
        } else {
          V <- .V
        }
      } else {
        print(sprintf(
          'VSpec ERROR: Variable %s is not in the supplied data.frame',
          .V
        ))
        return (NA)
      }
      if (type != 'spectrum' &&
          type != 'Welch' && type != 'MEM' && type != 'mem') {
        print (sprintf ('type %s is unavailable; using type=spectrum', type))
        type <- 'spectrum'
      }
      siglim <- 1  ## 1-sigma error limits
      if (type == 'spectrum') {
        if (!is.null(spans[1])) {
          if (!(spans[1] %% 2)) {
            spans[1] <- spans[1] + 1
          }
          if (spans[1] <= 5) {
            spans <- NULL
          }
        }
        S <-
          spectrum (ts(
            SmoothInterp(v, .maxGap = 1000 * Rate, .Length = 0),
            frequency = Rate
          ),
          span = spans,
          plot = FALSE)
        freq <- S$freq
        if (corrected) {
          ARX <- LTphase(freq, Par)
          S$spec <- S$spec / ARX$Amp ^ 2
        }
        fpf <- 2 * S$spec * freq
      } else if (type == 'Welch') {
        ## bspec section
        ## force segLength to a power of 2
        segl <- segLength
        rsl <- log(segl) / log(2)
        ns <- round (rsl)
        if (2 ^ ns != segl) {
          if (2 ^ ns > segl) {
            segl <- 2 ^ (ns - 1)
          }
          else {
            segl <- 2 ^ (ns + 1)
          }
          segLength <- segl
          print (sprintf ('reset segLength to %d', segLength))
        }
        S2 <-
          bspec::welchPSD (
            ts(SmoothInterp(v, .Length = 0), frequency = Rate),
            seglength = segLength,
            windowfun = bspec::hammingwindow
          )
        # ci <- quantile.bspec(BSP, probs = c(0.025, 0.975),
        #   two.sided = FALSE)
        coverage <- 0.683
        tail <- 1 - coverage
        df <- 2 * 9 * S2$segments / 11 ##1.768849
        upper.quantile <-
          1 - tail * pchisq(df, df, lower.tail = FALSE)
        lower.quantile <- tail * pchisq(df, df)
        ci <- 1 / (qchisq(c(upper.quantile, lower.quantile), df) / df)
        df <- 1.46 * (S2$segments + 1)
        lower.limit <- qchisq (pnorm(-siglim), df) / df
        upper.limit <- qchisq (pnorm(siglim), df) / df
        # ci <- 0.5 + (ci-0.5) / sqrt(9 * S2$segments / 11)
        # print (sprintf ('ci2=%.3f -- %.3f segments %d', ci[1], ci[2], S2$segments))
        freq <- S2$frequency[-1]
        fpf <- S2$power[-1] * freq
      } else if (type == 'MEM' || type == 'mem') {
        ## MEM section
        A <- memCoef (v, poles)
        ld <- nrow(.data)
        fmin <- log (Rate / ld)
        fmax <- log (0.5 * Rate)
        bins <- as.integer (1 / resolution)
        df <- (fmax - fmin) / bins
        fdtl <- fmin + df * (0:bins)
        freq <- exp (fdtl)
        psComplex <- memEstimate (freq / Rate, A) / Rate
        ps <- 2 * Rate * Mod (psComplex) ^ 2
        fpf <- freq * ps
      }
      tasAverage <- mean(.data$TASX, na.rm = TRUE)
      if (EDR) {
        ps <- fpf / freq
        fpf <- (2 * pi / tasAverage) * (1.5 * ps) ^ 1.5 * freq ^ 2.5
      }
      if (smoothBins > 9) {
        bs1 <- binStats(data.frame(fpf, log(freq)), bins = smoothBins)
        bs1 <-
          rbind (bs1,
                 data.frame(
                   xc = bs1$xc[nrow(bs1)],
                   ybar = bs1$ybar[nrow(bs1)],
                   sigma = bs1$sigma[nrow(bs1)],
                   nb = 1
                 ))
        bs1 <- bs1[!is.na(bs1$ybar), ]
        freq <- exp(bs1$xc)
        fpf <- bs1$ybar
        bs1$sigma <- ifelse (bs1$nb > 2, bs1$sigma / sqrt(bs1$nb), NA)
        rna <- is.na(bs1$sigma)
        bs1$sigma[rna] <- bs1$ybar[rna] / 2
        # bs1 <<- bs1
      }
      if (.V == .Variable[1]) {
        DF <- data.frame(freq, fpf)
      } else if (.V == .Variable[2]) {
        DF$fpf2 <- fpf
      } else if (.V == .Variable[3]) {
        DF$fpf3 <- fpf
      }
    }
    
    VL <- .Variable
    if (!is.na(VLabel[1])) {
      VL <- VLabel
    }
    if (is.na(col[1])) {
      col = c("blue", "forestgreen", "black", "darkorange")
    }
    if (!is.na(add[1])) {
      ADD <- add
    }
    for (.V in .Variable) {
      NV <- which(.V == .Variable)
      if (is.na(ADD[1])) {
        if (NV == 1) {
          ## first call: redefine VSpecDF
          try(rm(list = names(VSpecEnv), envir = VSpecEnv), silent = TRUE)
          VSpecEnv$Variable <- .Variable
          assign('VSpecDF1', DF, envir = VSpecEnv)
          labx <- 'frequency [Hz]'
          if (EDR) {
            # laby <- sprintf('eddy dissipation rate for %s', .V)
            laby <-
              expression(paste("eddy dissipation rate [m" ^ "2", "s" ^ "-3", "]"))
          } else {
            laby <- sprintf('variance spectrum fP(f) for %s', .V)
          }
          g <- ggplot(data = DF)
          g <-
            g + geom_path (aes(
              x = freq,
              y = fpf,
              colour = VL[1]
            ),
            data = DF,
            na.rm = TRUE) +
            xlab(labx) + ylab (laby)
          .clinesVSpec <- col[1]
          names(.clinesVSpec) <- VL[1]
          VSpecEnv$clinesVSpec <- .clinesVSpec
        }
        if (NV == 2) {
          g <-
            g + geom_path (aes(
              x = freq,
              y = fpf2,
              colour = VL[2]
            ),
            data = DF,
            na.rm = TRUE)
          cl2 <- ifelse (length(col) >= 2, col[2], 'forestgreen')
          names(cl2) <- VL[2]
          .clinesVSpec <- c(VSpecEnv$clinesVSpec, cl2)
          VSpecEnv$clinesVSpec <- .clinesVSpec
        } else if (NV == 3) {
          g <-
            g + geom_path (aes(
              x = freq,
              y = fpf3,
              colour = VL[3]
            ),
            data = DF,
            na.rm = TRUE)
          cl3 <- ifelse (length(col) >= 3, col[3], 'black')
          names(cl3) <- VL[3]
          .clinesVSpec <- c(VSpecEnv$clinesVSpec, cl3)
          VSpecEnv$clinesVSpec <- .clinesVSpec
        }
      } else {
        ## assign name based on elements in clinesVSpec
        N <- length(VSpecEnv$clinesVSpec) + 1
        nc <- names(VSpecEnv$clinesVSpec)
        .clinesVSpec <- c(VSpecEnv$clinesVSpec, col[N])
        names(.clinesVSpec) <- c(nc, V)
        VSpecEnv$clinesVSpec <- .clinesVSpec
        VName <- sprintf('VSpecDF%d', N)
        assign(VName, DF, pos = VSpecEnv)
        if (N == 2) {
          VSpecEnv$VSpecVar2 <- V
          g <-
            ADD + geom_path (
              aes(
                x = freq,
                y = fpf,
                colour = VSpecEnv$VSpecVar2
              ),
              data = get(VName, envir = VSpecEnv),
              na.rm = TRUE
            )
        } else if (N == 3) {
          VSpecEnv$VSpecVar3 <- V
          g <-
            ADD + geom_path (
              aes(
                x = freq,
                y = fpf,
                colour = VSpecEnv$VSpecVar3
              ),
              data = get(VName, envir = VSpecEnv),
              na.rm = TRUE
            )
        } else if (N == 4) {
          VSpecEnv$VSpecVar4 <- V
          g <-
            ADD + geom_path (
              aes(
                x = freq,
                y = fpf,
                colour = VSpecEnv$VSpecVar4
              ),
              data = get(VName, envir = VSpecEnv),
              na.rm = TRUE
            )
        }
      }
    }
    
    g <-
      suppressMessages(g + scale_colour_manual (name = '', values = .clinesVSpec))
    # print (.clinesVSpec)
    if (showErrors > 0) {
      if (smoothBins > 9) {
        bse <-
          data.frame(
            x = exp(bs1$xc),
            ymin = bs1$ybar - showErrors * bs1$sigma,
            ymax = bs1$ybar + showErrors * bs1$sigma
          )
        bse$ymin[bse$ymin < ylim[1]] <- ylim[1]
      } else {
        coverage <-
          pnorm(showErrors) - pnorm(-showErrors)  ## 1-sigma, 0.68269
        tail <- 1 - coverage
        if (type == 'spectrum') {
          df <- S$df
        } else if (type == 'Welch') {
          df <- 1.46 * (S2$segments + 1)
        } else if (type == 'MEM') {
          df <- length(v) / poles
        }
        uq <- 1 - tail * pchisq(df, df, lower.tail = FALSE)
        lq <- tail * pchisq(df, df, lower.tail = TRUE)
        ci <- 1 / (qchisq(c(uq, lq), df) / df)
        lower.limit <- qchisq (pnorm(-showErrors), df) / df
        upper.limit <- qchisq (pnorm(showErrors), df) / df
        bse <-
          data.frame(x = freq,
                     ymin = lower.limit * fpf,
                     ymax = upper.limit * fpf)
      }
      # g <- g + geom_ribbon(data=bs1, aes(x=exp(xc), ymin=max(ylim[1], ybar-showErrors*sigma), ymax=ybar+showErrors*sigma),
      #   fill='cyan', alpha=0.25, show.legend=FALSE, inherit.aes=FALSE, na.rm=TRUE)
      g <- g + geom_ribbon(
        data = bse,
        aes(x = x, ymin = ymin, ymax = ymax),
        fill = 'gray50',
        alpha = 0.5,
        show.legend = FALSE,
        inherit.aes = FALSE,
        na.rm = TRUE
      )
      # bs1$xc <- exp(bs1$xc)
    }
    if (is.na(ADD[1])) {
      g <-
        g + # scale_x_log10(
        # breaks = trans_breaks("log10", function(x)
        # 10 ^ x, n = 4),
        #limits = xlim,
        # labels = trans_format("log10", math_format(10 ^ .x))
        # ) +
        scale_y_log10(
          breaks = trans_breaks("log10", function(x)
            10 ^ x, n = 4),
          #limits = ylim,
          labels = trans_format("log10", math_format(10 ^ .x))
        ) +
        annotation_logticks(sides = 'trbl') +
        coord_cartesian(xlim = xlim, ylim = ylim)
      # g <- g + theme(panel.grid.minor = element_line(colour = "black"))
      if (EDR) {
        ## add line showing highest-decade average EDR
        imx <- length(freq)
        imn <- which (freq > freq[imx] / 20)[1]
        aveEDR <- mean(fpf[imn:imx], na.rm = TRUE)
        ttl <- sprintf ('EDR=%.2e', aveEDR)
        DFL <- data.frame(x = c(freq[imn], freq[imx]), y = rep(aveEDR, 2))
        g <-
          g + geom_path(data = DFL,
                        aes(x = x, y = y),
                        lwd = 1.5,
                        colour = 'red')
        g <- g + ggtitle (bquote(.(ttl) ~ ' m' ^ 2 ~ 's' ^ -3))
        # g <- g + ggtitle(sprintf(' mean eddy dissipation rate %.2e m^2/s^3', aveEDR))
      } else {
        for (i in (-8:0)) {
          a = ae * 10. ^ (i * (2 / 3)) * tasAverage ^ (2 / 3)
          lw = ifelse(i == -4, 1.2, 0.5)
          DFL <-
            data.frame(x = xlim, y = c(a / xlim[1] ^ (2 / 3), a / xlim[2] ^ (2 / 3)))
          # print(DFL)
          g <-
            g + geom_path (
              data = DFL,
              aes(x = x, y = y),
              colour = 'darkorange',
              lwd = lw,
              lty = 3
            )
        }
      }
      if (WavelengthScale) {
        yl <- c(ylim[1] * 1.2, ylim[1] * 1.5)
        lclr <- 'slategrey'
        for (j1 in c(10, 100, 1000, 10000, 100000)) {
          DFL2 <- data.frame(x = rep(tasAverage / j1, 2), y = yl)
          g <-
            g + geom_path(
              data = DFL2,
              aes(x = x, y = y),
              colour = lclr,
              lwd = 1.0
            )
          if (j1 != 100000) {
            for (j2 in 2:9) {
              DFL2 <- data.frame(x = rep(tasAverage / (j1 * j2), 2), y = yl)
              g <-
                g + geom_path(
                  data = DFL2,
                  aes(x = x, y = y),
                  colour = lclr,
                  lwd = 0.6
                )
            }
          }
        }
        DFL2 <-
          data.frame (x = tasAverage * c(1 / 10, 1 / 100000), y = rep(yl[1], 2))
        g <-
          g + geom_path(data = DFL2,
                        aes(x = x, y = y),
                        colour = lclr,
                        lwd = 1.0)
        g <- g + annotate(
          "text",
          x = tasAverage * c(1 / 100000, 1 / 10000, 1 / 1000, 1 / 100, 1 /
                               10),
          y = rep(yl[2] * 1.5, 5),
          label = c("100 km", "10 km", "1 km", "0.1 km", " "),
          colour = lclr
        )
      }
      # g <- g + theme_WAC()
    }
    if (!is.na(WACtheme)) {
      g <- g + theme_WAC()
    }
    return(g)
  }
