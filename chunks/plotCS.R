plotCS <- function(CS,
                   Units,
                   spans = 149,
                   fL,
                   wavelengthLimit = 2000,
                   smoothBins = 100,
                   xlim = c(0.01, 15),
                   ylim = c(0.001, 80),
                   printTitle = TRUE,
                   CSprevious = NA,
                   plotRibbon = TRUE,
                   vwp = NA,
                   project = NA) {
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
    BS$nybar <- -1 * BS$ybar
    BS$ybar[BS$ybar < 0] <- NA
    BS$nybar[BS$nybar < 0] <- NA
    BS$xc <- exp(BS$xc)
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
  g <- ggplot(data = CS, aes(x = freq))
  g <-
    g + geom_path(aes(
      y = cospec,
      colour = 'cospectrum',
      linetype = 'cospectrum',
      size = 'cospectrum'
    ))
  g <-
    g + geom_path(aes(
      y = ncospec,
      colour = '-cospectrum',
      linetype = '-cospectrum',
      size = '-cospectrum'
    ))
  g <-
    g + geom_path(aes(
      y = ogive,
      colour = 'exceedance',
      linetype = 'exceedance',
      size = 'exceedance'
    )
  )
  if (is.data.frame(CSprevious)) {
    g <- g + geom_path(
      data = CSprevious,
      aes(x = freq, y = ogive,
          colour = 'exceedance'),
      lty = 2,
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
    g <-
      g + geom_point(
        data = BS,
        aes(x = xc, y = ybar, colour = 'cospectrum', shape = 'cospectrum')
      )
    g <-
      g + geom_point(
        data = BS,
        aes(x = xc, y = nybar, colour = '-cospectrum', shape = '-cospectrum')
      )
    if (plotRibbon) {
      # GeomRibbon$handle_na <- function(data, params) {  data }
      g <- g + geom_ribbon(
        data = bse,
        aes(
          x = x,
          ymin = ymin,
          ymax = ymax
        ),
        fill = 'blue',
        alpha = 0.2,
        show.legend = FALSE,
        inherit.aes = FALSE,
        na.rm = FALSE
      )
      g <-
        g + geom_ribbon(
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
  g <-
    g + geom_path(data = data.frame(x = rep(fL, 2), y = ylim),
                  aes(x = x, y = y),
                  linetype = 2)
  g <-
    g + scale_x_log10(
      breaks = trans_breaks("log10", function(x)
        10 ^ x, n = 4),
      #limits = xlim,
      labels = trans_format("log10", math_format(10 ^ .x))
    ) +
    scale_y_log10(
      breaks = trans_breaks("log10", function(x)
        10 ^ x, n = 4),
      #limits = ylim,
      labels = trans_format("log10", math_format(10 ^ .x))
    ) +
    annotation_logticks(sides = 'trbl') +
    coord_cartesian(xlim = xlim, ylim = ylim)
  g <- g + xlab('frequency [Hz]') + ylab(ylab)
  g <- g + scale_size_manual(
    name = '', 
    # guide = guide_legend(override.aes = list(linewidth = c(0.1, 0.1, 1.2) ) )
    values = c(
      'cospectrum' = 1,
      '-cospectrum' = 1,
      'exceedance' = 1.2
    )
  )
  g <- suppressWarnings(g + scale_colour_manual (
    name = '',
    values = c(
      'cospectrum' = 'blue',
      '-cospectrum' = 'red',
      'exceedance' = 'grey50'
    )
  ))
  g <- g + scale_linetype_manual (
    name = '',
    values = c(
      'cospectrum' = 'solid',
      '-cospectrum' = 'solid',
      'exceedance' = 'solid'
    )
  )
  g <- g + scale_shape_manual(
    name = '', 
    values = c(
      'cospectrum' = 19,
      '-cospectrum' = 19,
      'exceedance' = 32
    )
  )
  # g <- g + guides(linetype = guide_legend(override.aes = list(size = c(0.8, 0.8, 1.2))))
  # g <- g + guides(guide_legend(override.aes = list(linewidth = c(0.1, 0.1, 1.2) ) ))
  # g <-
  #   g + guides(col = guide_legend(reverse = TRUE),
  #              linetype = guide_legend(reverse = TRUE))
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
  if (!is.na(project)) {
    g <-
      g + annotate(
        'text',
        x = 4,
        y = 7,
        label = paste0('(', project, ')'),
        size = 4
      )
  }
  if (is.na(vwp[1])) {
    suppressWarnings(print(
      g + theme_WAC(1) + theme(plot.title = element_text(size = 12),
                               legend.text = element_text(size=14),
                               legend.position = c(0.5, 0.91))
    ))
  } else {
    suppressWarnings(print(
      g + theme_WAC(1) + 
        theme(
          # legend.key.size = unit(3,"line"),
          legend.position = c(0.67, 0.89),
          legend.text = element_text(size=11),
          plot.margin = unit(c(0.6, 0.6, 0.8, 0.6), "lines"),
          plot.title = element_text(size = 12),
          axis.title.y = 
            element_text (face='plain', size=12, color='blue', 
                          margin=margin(0,10,0,0), angle=90)
        ), vp = vwp)
    )
  }
  g <<- g  ## Save for adding uncorrected cospec for debugging...
  return(CS)
}