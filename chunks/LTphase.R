# The Laplace-transform solution, given the parameters in P:
LTphase <- function(f, P) {
  ## f=frequency; P=Param
  tau1 <- P$tau1
  tau2 <- P$tau2
  a <- P$a
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
  return(list('Amp' = cTF, 'Phase' = phiTF))
}
