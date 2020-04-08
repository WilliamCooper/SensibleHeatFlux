
CompareTF <- function(f, a, tau2, tau1) {
  b <- sqrt(1 / (1 + (2 * pi * f * tau2) ^ 2))
  zeta <- -atan(2 * pi * f * tau2)
  C1 <- 1 / (1 + 4 * pi ^ 2 * f ^ 2 * tau1 ^ 2) *
    (-(a + (1 - a) * b * cos(zeta)) * 2 * pi * f * tau1 +
       (1 - a) * b * sin(zeta))
  C2 <- 1 / (1 + 4 * pi ^ 2 * f ^ 2 * tau1 ^ 2) *
    ((a + (1 - a) * b * cos(zeta)) +
       (1 - a) * b * sin(zeta) * 2 * pi * f * tau1)
  cTF <- sqrt(C1 ^ 2 + C2 ^ 2)
  phiTF <- atan2(C1, C2) * 180 / pi
  A2 <- (1-a)*tau2 / (tau2 - tau1)
  A1 <- 1-A2
  ## Reference value
  C1a <- -2*pi*f*((A1*tau1)/(1+(2*pi*f*tau1)^2)+(A2*tau2)/(1+(2*pi*f*tau2)^2))
  omega <- 2*pi*f
  ## checked alternate
  C1a <- -omega / ((1 + omega^2 * tau1^2) * (1 + omega^2 * tau2^2)) *
      (a*tau2*(-1+tau1*tau2*omega^2) + tau1 + tau2)
  # checked alternate
  # C1a <- -omega / ((1 + omega^2 * tau1^2) * (1 + omega^2 * tau2^2)) *
      # (a*(tau1*(1+omega^2 * tau2^2)-(tau2+tau1))+(tau2+tau1))
  # checked alternate
  C1a <- -omega / ((1 + omega^2 * tau1^2)) *
      (a*tau1 + (1-a) * (tau1+tau2) / (1 + omega^2 * tau2^2))
  C2a <- ((A1)/(1+(2*pi*f*tau1)^2)+(A2)/(1+(2*pi*f*tau2)^2))
  cTFa <- sqrt(C1a^2 + C2a^2)
  phiTFa <- atan2(C1a, C2a) * 180 / pi 
  AmpDiff <- cTF - cTFa
  PhaseDiff <- phiTF - phiTFa
  DF <- data.frame(Time=f, Amp=cTF, Amp2=cTFa, Phase=phiTF, Phase2= phiTFa, AmpD = AmpDiff, PhaseD = PhaseDiff)
  g <- DF %>% select(Time, Amp, Amp2, Phase, Phase2) %>% ggplotWAC(labx='Frequency', panels=2, labelP=c('Gain', 'Phase'), labelL=c('std', 'IF'))
  # g <- ggplotWAC(DF, labx='Frequency', panels=2, labelP=c('Gain', 'Phase'))
  g <- g + xlab('frequency [Hz]') + ylab(bquote('transfer function H(' * nu ~ ')'))
  g <-
    g + scale_x_log10(
      breaks = trans_breaks("log10", function(x)
        10 ^ x, n = 4),
      labels = trans_format("log10", math_format(expr = 10 ^ .x))
    ) + xlab('frequency [Hz]')
  g <- g + annotation_logticks(sides = 'tb')
  g <- g + theme_WAC(1) + theme(legend.position = c(0.8, 0.95))
#   print(g)
#   g <- DF %>% select(Time, AmpD, PhaseD) %>% 
#              ggplotWAC(labx='Frequency', panels=2, labelP=c('Gain', 'Phase'))
}
