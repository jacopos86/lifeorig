#
#  This module implements
#  the Gillespie algorithm
#  for reaction kinetics
#  Stochastic simulation algorithm :
#  0. Initialize (t=t0, x=x0)
#  1. with system in state (x,t) compute aj(x), a0(x)
#  2. generate tau, j
#  3. next reaction t <- t + tau ; x <- x + vj
#  4. record (x,t) return to 1. 