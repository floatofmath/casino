##' To use mclapply2 from bt88.03.704 you have to install this library (e.g. \code{devtools::install_github('floatofmath/bt88.03.704')}) unfortunately the statusbar implementation does not work with Rstudio.
##' 
##' @title Run a batch of simulations
##' @param scenarios Object with paramter settings for \code{delta}, \code{sigma} and \code{n1} in each row 
##' @param runs Number of simulation runs for each scenario
##' @param run_scenario Function that runs a single parameter scenario, takes at least the columns of \code{G} and \code{runs} as arguments
##' @param ... Additional named arguments for run_scenario
##' @param use_mclapply2 use mclapply2 from package bt88.03.704 which implements a status bar (see details)
##' @param multicore Should multicore parallelisation be used
##' @return Object with parameters and simulation results in each row
##' @author Florian Klinglmueller
##' @examples
##' simulate_ttest <- function(m,s,n,runs,alpha=.05){
##'   x <- matrix(rnorm(runs*n,m,s),ncol=runs)
##'   decisions <- apply(x,2,function(x) { t.test(x)$p.value <= alpha })
##'   c(n=n,power=mean(decisions))
##' }
##'
##' scenarios <- expand.grid(m=c(0,1),s=c(1,2))
##' simulate_batch(scenarios,100,simulate_ttest,n=10)
##' @export
simulate_batch <- function(scenarios,runs,run_scenario,...,use_mclapply2=FALSE,multicore=TRUE){
    ## Parallelization
    if(multicore) {
        require(parallel)
        if(use_mclapply2) {
            require(bt88.03.704)
            mcla <- mclapply2
        }
        mcla <- mclapply
    } 
    if(!multicore) {
        mcla  <-  lapply
    }
    params <- list(...)
    dplyr::bind_rows(mcla(1:nrow(scenarios),function(p) c(scenarios[p,],result=do.call('run_scenario',c(runs=runs,scenarios[p,],params)))))
}


