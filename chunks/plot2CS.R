plot2CS <- function(CS,
                   Units,
                   spans = 149,
                   fL,
                   wavelengthLimit = 2000,
                   smoothBins = 100,
                   xlim = c(0.01, 15),
                   ylim = c(0.001, 80),
                   printTitle = TRUE,
                   CSprevious = NA,
                   plotFigure = TRUE,
                   plotRibbon = TRUE,
                   showNegative = TRUE,
                   ADD = FALSE) {
  if (ADD) {
    ## Add to an existing plot definition that is saved in a previous call
    if (exists('g.CS'))
      g <- g.CS  ## Saved from a previous call
  }
  ylab <- bquote(nu * " x flux cospectrum [" * .(Units) * "]")
  CSogive <- cumsum(CS$cospec) * CS$freq[1]
  CSogive <- CSogive[length(CSogive)] - CSogive
  CS$ogive <- CSogive
  CS$cospec <- SmoothInterp(CS$cospec, .Length = 0)  # treat NAs
  s25 <- spans %/% 25
  s10 <- spans %/% 10
  s3 <- spans %/% 3
  s25 <- s25 + (s25 + 1) %% 2
  s10 <- s10 + (s10 + 1) %% 2
  s3 <- s3 + (s3 + 1) %% 2
  CS$cospec <-
    zoo::rollapply(CS$cospec,
                   FUN = mean,
                   fill = 'extend',
                   width = s25)
  CS$cospec[CS$freq > 0.01] <-
    zoo::rollapply(CS$cospec,
                   FUN = mean,
                   fill = 'extend',
                   width = s10)[CS$freq > 0.01]
  CS$cospec[CS$freq > 0.1] <-
    zoo::rollapply(CS$cospec,
                   FUN = mean,
                   fill = 'extend',
                   width = s3)[CS$freq > 0.1]
  CS$cospec[CS$freq > 1] <-
    zoo::rollapply(CS$cospec,
                   FUN = mean,
                   fill = 'extend',
                   width = spans)[CS$freq > 1]
  FluxL <- CSogive[which(CS$freq > fL)[1]]
  Flux <- CSogive[which(CS$freq > 0.01)[1]]
  attr(CS, 'Flux') <- Flux
  attr(CS, 'FluxL') <- FluxL
  attr(CS, 'wavelengthLimit') <- wavelengthLimit
  ## Construct the plot:
  CS$ncospec <- -1 * CS$cospec
  ## Weight by frequency for log-abscissa plot:
  CS$cospec <- CS$cospec * CS$freq
  CS$ncospec <- CS$ncospec * CS$freq
  if (smoothBins > 5) {
    BS <-
      binStats(data.frame(CS$cospec, log(CS$freq)), bins = smoothBins)
    # lines(exp(BS$xc), BS$ybar, lwd=2, col='brown')
    BS$nybar <- -1 * BS$ybar
    BS$ybar[BS$ybar < 0] <- NA
    BS$nybar[BS$nybar < 0] <- NA
    BS$xc <- exp(BS$xc)
    # lines(exp(BS$xc), BS$nybar, lwd=2, col='magenta')
    attr(CS, 'smoothed data.frame') <- BS
    bse <-
      data.frame(
        x = BS$xc,
        ymin = BS$ybar - BS$sigma,
        ymax = BS$ybar + BS$sigma,
        yminN = BS$nybar - BS$sigma,
        ymaxN = BS$nybar + BS$sigma
      )
  }
  if (smoothBins > 5) {
    bse$ymin[bse$ymin < ylim[1]] <- ylim[1]
    bse$yminN[bse$yminN < ylim[1]] <- ylim[1]
  }
  if (ADD) {
    g <- g + geom_path(data = CS,
                       aes(
                         y = cospec,
                         colour = 'measured',
                         linetype = 'measured'
                       ))
    if (showNegative) {
      g <- g + geom_path(aes(
        y = ncospec,
        colour = '-cosp2',
        linetype = '-cosp2'
      ))
    }
    g <- g + geom_path(
      data = CS,
      aes(
        y = ogive,
        colour = 'exc-msrd',
        linetype = 'exc-msrd'
      ),
      lwd = 1.2
    )
  } else {
    g <- ggplot(data = CS, aes(x = freq))
    g <-
      g + geom_path(data = CS,
                    aes(
                      y = cospec,
                      colour = 'generated',
                      linetype = 'generated'
                    ))
    if (showNegative) {
      g <-
        g + geom_path(aes(
          y = ncospec,
          colour = '-cosp',
          linetype = '-cosp'
        ))
    }
    g <- g + geom_path(aes(
      y = ogive,
      colour = 'exceedance',
      linetype = 'exceedance'
    ),
    lwd = 1.2)
  }
  if (is.data.frame(CSprevious)) {
    g <- g + geom_path(
      data = CSprevious,
      aes(
        x = freq,
        y = ogive,
        colour = 'exc-corr',
        linetype = 'exc-corr'
      ),
      lwd = 1.2
    )
    if ('ogive2' %in% names(CSprevious)) {
      g <- g + geom_path(
        data = CSprevious,
        aes(
          x = freq,
          y = ogive2,
          colour = 'generated',
          linetype = 'generated'
        ),
        lwd = 1.3
      )
    }
  }
  if (smoothBins > 5) {
    pclr <- ifelse(ADD, 'measured', 'generated')
    pclr <- ifelse(ADD, 'forestgreen', 'blue')
    if (ADD) {
      g <-
        g + geom_point(data = BS,
                       aes(
                         x = xc,
                         y = ybar,
                         colour = 'measured'
                       ),
                       pch = 19)
    } else {
      g <-
        g + geom_point(data = BS,
                       aes(
                         x = xc,
                         y = ybar,
                         colour = 'generated'
                       ),
                       pch = 19)
    }
    if (showNegative) {
      g <-
        g + geom_point(
          data = BS,
          aes(x = xc, y = nybar),
          colour = 'darkred',
          pch = 19
        )
    }
    if (plotRibbon) {
      # GeomRibbon$handle_na <- function(data, params) {  data }
      g <- g + geom_ribbon(
        data = bse,
        aes(
          x = x,
          ymin = ymin,
          ymax = ymax
        ),
        fill = pclr,
        alpha = 0.2,
        show.legend = FALSE,
        inherit.aes = FALSE,
        na.rm = FALSE
      )
      if (showNegative) {
        g <- g + geom_ribbon(
          data = bse,
          aes(
            x = x,
            ymin = yminN,
            ymax = ymaxN
          ),
          fill = 'red',
          alpha = 0.2,
          show.legend = FALSE,
          inherit.aes = FALSE,
          na.rm = FALSE
        )
      }
    }
  }
  if (!ADD) {
    g <-
      g + geom_path(data = data.frame(x = rep(fL, 2), y = ylim),
                    aes(x = x, y = y),
                    linetype = 2)
    g <-
      g + scale_x_log10(
        breaks = trans_breaks("log10", function(x)
          10 ^ x, n = 2),
        labels = trans_format("log10", math_format(10 ^ .x))
      ) +
      scale_y_log10(
        breaks = trans_breaks("log10", function(x)
          10 ^ x, n = 4),
        labels = trans_format("log10", math_format(10 ^ .x))
      ) +
      annotation_logticks(sides = 'trbl') +
      coord_cartesian(xlim = xlim, ylim = ylim)
    g <- g + xlab('frequency [Hz]') + ylab(ylab)
    ttl <-
      bquote(
        'Total flux ' ~ .(format(Flux, digits = 3)) ~ .(Units) * '; partial <' *
          .(format((
            wavelengthLimit / 1000
          ), digits = 2)) ~ 'km:' ~ .(format(FluxL, digits = 3)) ~ .(Units)
      )
    if (printTitle) {
      g <- g + labs(title = ttl)
    }
    g <-
      g + theme_WAC(1) + theme(plot.title = element_text(size = 12)) +
      theme(legend.position = c(0.5, 0.91))
  } else {
    ## this is for ADD == TRUE
    if (showNegative) {
      g <- suppressWarnings(g + scale_colour_manual (
        name = '',
        values = c(
          'cosp' = 'blue',
          '-cosp' = 'red',
          'exceedance' = 'brown',
          'generated' = 'forestgreen',
          'measured' = 'forestgreen',
          '-cosp2' = 'darkorange',
          'exceedance2' = 'black'
        )
      ))
      g <- g + scale_linetype_manual (
        name = '',
        values = c(
          'cosp' = 1,
          '-cosp' = 1,
          'exceedance' = 1,
          'generated' = 4,
          'measured' = 1,
          '-cosp2' = 1,
          'exceedance2' = 2
        )
      )
    } else {
      g <- g + scale_linetype_manual (
        name = '',
        values = c(
          'generated' = 1,
          'exceedance' = 1,
          'exc-msrd' = 2,
          'measured' = 1,
          'exc-corr' = 2
        )
      )
      g <-
        g + scale_shape_manual(
          name = '',
          values = c(
            'generated' = 19,
            'exceedance' = NA,
            'exc-msrd' = NA,
            'measured' = 19,
            'exc-corr' = NA
          )
        )
      g <- g + scale_colour_manual (
        name = '',
        values = c(
          'generated' = 'blue',
          'exceedance' = 'brown',
          'exc-msrd' = 'brown',
          'measured' = 'forestgreen',
          'exc-corr' = 'black'
        ),
        guide = guide_legend(override.aes = list(
          shape = c('generated' = 19,
                    'measured' = 19)
        ))
      )
      g <- g + guides(
        colour = guide_legend(reverse = TRUE),
        linetype = guide_legend(reverse = TRUE),
        shape = guide_legend(reverse = TRUE)
      )
      g <-
        g + guides(shape = guide_legend(override.aes = list(shape = c(
          19, NA, NA, 19, NA
        ))))
    }
  }
  if (plotFigure)
    suppressWarnings(print(g))
  g.CS <<- g  ## Save for adding additional lines...
  return(CS)
}