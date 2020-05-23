## This is a special version of CohPhase with features not in the standard Ranadu version
CohP <-
  function (.data,
            .Var1,
            .Var2,
            col = 'blue',
            spans = 25,
            smoothBins = 50,
            plotType = 'ggplot',
            showErrors = 0,
            returnCospectrum = FALSE) {
    if (is.data.frame(.data)) {
      if (.Var1 %in% names(.data)) {
        Z <-
          capture.output (Vr <-
                            SmoothInterp(detrend (.data[, c('Time', .Var1)]), .Length = 0))
      } else {
        print(sprintf(
          'CohPhase ERROR: Variable %s is not in the supplied data.frame',
          .Var1
        ))
        return (NA)
      }
      if (.Var2 %in% names(.data)) {
        Z <-
          capture.output (VrC <-
                            SmoothInterp(detrend (.data[, c('Time', .Var2)]), .Length = 0))
      } else {
        print(sprintf(
          'CohPhase ERROR: Variable %s is not in the supplied data.frame',
          .Var2
        ))
        return (NA)
      }
    } else {
      print('CohPhase ERROR: first argument is not a data.frame.')
      return (NA)
    }
    if (is.null(attr(.data, 'Rate'))) {
      print ('CohPhase warning: Rate attribute missing from data.frame, so using Rate=1')
      Rate <- 1
    } else {
      Rate <- attr(.data, 'Rate')
    }
    vcv <-
      cbind(ts(Vr, frequency = Rate), ts(VrC, frequency = Rate))
    P <-
      spec.pgram(
        vcv,
        detrend = FALSE,
        fast = TRUE,
        plot = FALSE,
        spans = spans
      )
    df1 <- data.frame(P$coh, log(P$freq))
    df2 <- data.frame (P$phase, log(P$freq))
    df3 <- data.frame (P$spec[, 1], log(P$freq))
    df4 <- data.frame (P$spec[, 2], log(P$freq))
    pf1 <- binStats (df1, bins = smoothBins)
    pf2 <- binStats (df2, bins = smoothBins)
    pf3 <- binStats (df3, bins = smoothBins)
    pf4 <- binStats (df4, bins = smoothBins)
    pf1 <- pf1[!is.na (pf1$ybar), ]
    pf2 <- pf2[!is.na (pf2$ybar), ]
    pf3 <- pf3[!is.na (pf3$ybar), ]
    pf4 <- pf4[!is.na (pf4$ybar), ]
    # pf1$sigma[pf1$nb > 1] <- pf1$sigma[pf1$nb > 1] / sqrt(pf1$nb[pf1$nb > 2])
    pf1$sigma[pf1$nb <= 1] <- NA # pf1$ybar[pf1$nb <= 1] * 0.5
    # pf2$sigma[pf2$nb > 1] <- pf2$sigma[pf2$nb > 1] / sqrt(pf2$nb[pf2$nb > 2])
    pf2$sigma[pf2$nb <= 1] <- NA # pf2$ybar[pf2$nb <= 1] * 0.5is
    if (plotType != 'ggplot') {
      pf1 <- binStats (df1, bins = smoothBins, addBin = TRUE)
      pf2 <- binStats (df2, bins = smoothBins, addBin = TRUE)
      pf3 <- binStats (df3, bins = smoothBins, addBin = TRUE)
      pf4 <- binStats (df4, bins = smoothBins, addBin = TRUE)
      return(cbind(pf1, pf2, pf3, pf4)[, c(2, 1, 3, 4, 6, 7, 10)])
    } else {
      d2 <-
        data.frame(
          Time = exp(pf1$xc),
          coherence = pf1$ybar,
          phase = pf2$ybar * 180 / pi,
          clo = (pf1$ybar - showErrors * pf1$sigma),
          chi = pf1$ybar + showErrors * pf1$sigma,
          plo = (pf2$ybar - showErrors * pf2$sigma) * 180 / pi,
          phi = (pf2$ybar + showErrors * pf2$sigma) * 180 / pi
        )
      d2$clo[!is.na(d2$clo) & (d2$clo < 0)] <- 0
      labelP <- c('coherence', 'phase [degrees]')
      g <- ggplotWAC(
        d2[, c(1, 2, 3)],
        panels = 2,
        labelP = labelP,
        col = col,
        lwd = c(1.0),
        lty = c(1),
        xlab = 'freq'
      )
      g <-
        g + xlab('frequency [Hz]') + ylab (sprintf ('%s x %s', .Var1, .Var2))
      g <-
        g + scale_x_log10(
          breaks = trans_breaks("log10", function(x)
            10 ^ x, n = 4),
          labels = trans_format("log10", math_format(expr = 10 ^ .x))
        ) + xlab('frequency [Hz]')
      if (showErrors > 0 && smoothBins > 5) {
        da <- data.frame(d2[, c(1, 4, 5)])
        db <- data.frame(d2[, c(1, 6, 7)])
        names(da) <- c('Time', 'ymin', 'ymax')
        names(db) <- c('Time', 'ymin', 'ymax')
        da$PanelGroup <- labelP[1]
        db$PanelGroup <- labelP[2]
        d <- rbind(db, da)
        g <-
          g + geom_ribbon(
            data = d,
            aes(
              x = Time,
              ymin = ymin,
              ymax = ymax
            ),
            colour = 'grey',
            alpha = 0.15,
            inherit.aes = FALSE
          )
      }
      g <- g + theme_WAC(1) + theme(legend.position = 'none')
      if (returnCospectrum) {
        CS <-
          sqrt(P$coh[, 1] * P$spec[, 1] * P$spec[, 2] / (1 + tan(P$phase[, 1]) ^
                                                           2))
        v1 <- SmoothInterp(.data[, .Var1], .Length = 0)
        v2 <- SmoothInterp(.data[, .Var2], .Length = 0)
        v1 <- detrend(data.frame(Time = .data$Time, v1))
        v2 <- detrend(data.frame(Time = .data$Time, v2))
        ff1 <- fft(v1)
        ff2 <- fft(v2)
        G <- Re(ff1 * Conj(ff2)) / nrow(.data)
        GQ <- Im(ff1 * Conj(ff2)) / nrow(.data)
        N <- nrow(.data) %/% 2
        S1 <- Re(ff1 * Conj(ff1) / nrow(.data))
        S2 <- Re(ff2 * Conj(ff2) / nrow(.data))
        G <- G[2:(N + 1)]
        GQ <- GQ[2:(N + 1)]
        S1 <- S1[2:(N + 1)]
        S2 <- S2[2:(N + 1)]
        frq <- c(1:N) * Rate / nrow(.data)
        spec1 <- 2 * S1 / Rate
        spec2 <- 2 * S2 / Rate
        cospec <- 2 * G / Rate
        quad <- 2 * GQ / Rate
        # cospectrum - i * quadrature = (gain spectrum) * exp(i*(phase spectrum))
        # sqrt(cospectrum^2 + quadrature^2) is the amplitude or gain spectrum
        return(data.frame(
          freq = frq,
          cospec = cospec,
          quad = quad,
          spec1 = spec1,
          spec2 = spec2
        ))
      } else {
        return(g)
      }
    }
  }