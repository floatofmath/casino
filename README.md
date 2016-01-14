# The _casino_ r-package

This package collects a number of utility functions useful for
performing simulation studies. Especially we provide tools that
facilitate monte carlo optimization based on simulation results.

## Tools that facilitate general simulation studies

These consist mainly of wrappers around `(mc)lapply(2)` to run
simulations for a large number of scenarios. Simulations may be run in
parallel using fork-based parallellization.

## Tools that facilitate stochastic optimization

Currently most effort goes into implementing methods that facilitate
stochastic optimization where the goal is to optimiza a uitility
function that may only be evaluated by simulation.

