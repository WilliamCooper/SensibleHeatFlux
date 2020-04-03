## rk4.integrate
## [Do entire integration in one call instead, as for runge.kutta]
## [i.e., function (ya, dydt, tv, tol=0.005) and return entire vector]
rk4.integrate <- function (dydt, ystart, tv, tol = 0.005) {
  ## Try the full step and accept if error estimate is < tol;
  ## otherwise calculate the number of steps that will give the desired
  ## tolerance and divide the step into that number of smaller steps.
  L <- length(tv)
  yr <- rep(ystart, L)
  y <- ystart
  dt <- 1
  for (it in 2:L) {
    dt <- tv[it] - tv[it - 1]
    t <- tv[it - 1]
    RK4 <- rk4.step(y, t, dt, dydt)
    # print (sprintf ('rk4 return for y=%.2f and t=%.1f is %.2f with error estimate %.3f',
    #                 y1, t, RK4[[1]], RK4[[2]]))
    if (abs(RK4[[2]]) <= tol) {
      yr[it] <- y <- RK4[[1]]
    } else {
      N <- as.integer (1.1 * (abs(RK4[[2]]) / tol) ^ 0.3 + 1)
      # print (sprintf ('error estimate is %.3f so using %d steps', RK4[[2]], N))
      
      for (n in 1:N) {
        A <- rk4.step(y, t, dt / N, dydt)
        y <- A[[1]]
        if (abs(A[[2]]) > tol) {
          # print (sprintf ('Warning: error estimate still too large, t=%.2f, y=%.2f, err=%.2f', t, y, A[[2]]))
          # Try further reduced step size
          M <- as.integer (1.1 * (abs(A[[2]]) / tol) ^ 0.3 + 1)
          for (m in 1:M) {
            y <- rk4.step(y, t, dt / (N * M), dydt)[[1]]
            t <- t + dt / (N * M)
          }
        } else {
          t <- t + dt / N
        }
      }
      yr[it] <- y
    }
  }
  return(yr)
}
## These are the coefficients used by the Runge-Kutta Cash-Karp integration:
RKD <- data.frame(
  a2 = 0.2,
  a3 = 0.3,
  a4 = 0.6,
  a5 = 1,
  a6 = 0.875,
  b21 = 0.2,
  b31 = 3 / 40,
  b32 = 9 / 40,
  b41 = 0.3,
  b42 = -0.9,
  b43 = 1.2,
  b51 = -11 / 54,
  b52 = 2.5,
  b53 = -70 / 27,
  b54 = 35 / 27,
  b61 = 1631 / 55296,
  b62 = 175 / 512,
  b63 = 575 / 13824,
  b64 = 44275 / 110592,
  b65 = 253 / 4096,
  c1 = 37 / 378,
  c3 = 250 / 621,
  c4 = 125 / 594,
  c6 = 512 / 1771
)
RKD <-
  cbind(
    RKD,
    data.frame(
      dc1 = RKD$c1 - 2825 / 27648,
      dc3 = RKD$c3 - 18575 / 48384,
      dc4 = RKD$c4 - 13525 / 55296,
      dc5 = -277 / 14336,
      dc6 = RKD$c6 - 0.25
    )
  )
## Save this in a special environment that rk4.step() can access:
if (!exists('RKCKEnv', envir = emptyenv())) {
  # define if absent
  assign('RKCKEnv', new.env(parent = emptyenv()), envir = globalenv())
}
RKCKEnv$RKD <- RKD

rk4.step <- function (y1, t, dt, dydt) {
  RKD <- RKCKEnv$RKD
  dydt1 <- dydt(y1, t)
  trial2 <- dydt(y1 + dt * RKD$b21 * dydt1, t + RKD$a2 * dt)
  trial3 <- dydt(y1 + dt *
                   (RKD$b31 * dydt1 + RKD$b32 * trial2), t + RKD$a3 * dt)
  trial4 <- dydt(y1 + dt *
                   (RKD$b41 * dydt1 + RKD$b42 * trial2 + RKD$b43 * trial3),
                 t + RKD$a4 * dt)
  trial5 <- dydt(
    y1 + dt *
      (
        RKD$b51 * dydt1 + RKD$b52 * trial2 + RKD$b53 * trial3 +
          RKD$b54 * trial4
      ),
    t + RKD$a5 * dt
  )
  trial6 <- dydt(
    y1 + dt *
      (
        RKD$b61 * dydt1 + RKD$b62 * trial2 + RKD$b63 * trial3 +
          RKD$b64 * trial4 + RKD$b65 * trial5
      ),
    t + RKD$a6 * dt
  )
  yn <- y1 + dt *
    (RKD$c1 * dydt1 + RKD$c3 * trial3 + RKD$c4 * trial4 + RKD$c6 * trial6)
  err <-
    dt * (RKD$dc1 * dydt1 + RKD$dc3 * trial3 + RKD$dc4 * trial4 +
            RKD$dc5 * trial5 + RKD$dc6 * trial6)
  # print (sprintf ('rkck.integrate return is %.5f for input y=%.5f and t=%.2f', yn, y, t))
  return(list(yn, err))
}
