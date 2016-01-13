##' A strange function with error to maximize
##'
##' @title Disturbed function
##' @param m parameter \code{m}
##' @param s parameter \code{s}
##' @param n parameter \code{s} 
##' @param runs parameter \code{runs} 
##' @return numeric value of the disturbed function
##' @examples
##' par(mfrow=c(4,1))
##' plot(x <- seq(0,1,.001),simulate_poly(x,1,100,10^6))
##' plot(x <- seq(0,1,.001),simulate_poly(x,1.3,100,10^6))
##' plot(x <- seq(0,1,.001),simulate_poly(x,1.6,100,10^6))
##' plot(x <- seq(0,1,.001),simulate_poly(x,2,100,10^6))
##' @author Florian Klinglmueller
##' @export
simulate_poly <- function(m,s,n,runs){
   2*s*m-m^2-m^3-s + rnorm(length(m),sd=1000/sqrt(runs*n))
}


    
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
##'
##' scenarios <- expand.grid(m=c(0,1),s=c(1,2))
##' simulate_batch(scenarios,10000,simulate_poly,n=10)
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

##' This function selects rows from a dataframe of simulation results. For example using the default \code{functional} \code{which.max} the rows containing the maximum of column \code{what} within values of \code{group} are selected. The function adds column with name \code{select_what} that specifies based on which columns a row was selected.
##'
##' 
##' @title Select results
##' @param sim_results Simulation results from a previous simulation run
##' @param what Variable based on which the selection is performed
##' @param group Grouping variables within unique values of which results should be selected
##' @param functional Function that returns a unique vector position
##' @return
##'
##' @examples
##' scenarios <- expand.grid(m=c(0,1),s=c(1,2),n=10:20)
##' results <- simulate_batch(scenarios,10000,simulate_poly)
##' select_results(results,'n',"result")
##' @author Florian Klinglmueller
##' @import dplyr
##' @export
select_results <- function(sim_results,group,what,functional=which.max){
    sim_results %>%
        group_by_(.dots=group) %>% 
            slice_(lazyeval::interp(~functional(m),m=as.name(what))) %>% ungroup -> out
    out$select_what <- what
    out
}


##' Expand and refine parameter grid for simulation
##'
##' @title Add neighbouring parameter values
##' @param scenarios matrix with one parameter setting in each line
##' @param parameter list of parameters for which grid points should be added
##' @param radius width of interval to be added around each point in \code{scenarios}
##' @param points number of points to fill the interval with
##' @param drop should columns that start with \code{"result"} or \code{"select"} or  be dropped (which may have been added by \code{\link{simulate_batch}}, or \code{\link{select_results}})
##' @return object with the same columns as G
##' @examples
##' scenarios <- expand.grid(m=c(0,1),s=c(1,2),n=5:10)
##' neighbourhoods <- add_neighbours(scenarios,c('m','s'),c(.2,.4),4)
##' 
##' @author Florian Klinglmueller
##' @import data.table
##' @export
add_neighbours <- function(scenarios,parameter,radius,points,drop=T){
    by <- 2*radius/(points-1)
    scenarios <- as.data.table(scenarios)
    if(length(radius) != 1 | length(by) != 1){
        if(length(radius) != length(parameter) | {length(by)!= 1 & length(by) != length(parameter)})
            stop("Vector of interval widths has to have same length as number of columns to refine")
        if(length(by) != length(radius)){
            ## replicate the shorter one
            by <- by+radius*0
            radius <- radius+by*0
        }
    } else {
        radius <- rep(radius,length(parameter))
        by <- rep(by,length(parameter))
    }
    int <- lapply(1:length(by),function(i) seq(-radius[[i]],+radius[[i]],by[[i]]))
    for(j in 1:length(parameter)){
        .c <- substitute(`:=`(what,as.numeric(outer(int[[j]],scenarios[,what],'+'))),list(what=as.name(parameter[j])))
        scenarios <- scenarios[rep(1:.N,each=length(int[[j]]))][,eval(.c)]
    }
    out <- tbl_df(data.frame(scenarios))
    if(drop) out <- select(out,-matches('^result|^select'))
    return(out)
}

##' Fits a loess regression and adds candidates scenarios within the 95 percent confidence band around the fit. 
##'
##' @title smooth parameter grid
##' @param scenarios list of candadidate parameter settings
##' @param parameter parameter to smooth
##' @param groups explanatory variable
##' @param points number of points around each prediction
##' @param drop should columns that start with \code{"result"} or \code{"select"} or  be dropped (which may have been added by \code{\link{simulate_batch}}, or \code{\link{select_results}})
##' @param width width of the band around the smoothed fit in multiples of the standard error
##' @param ... additional options to \code{\link{stats::loess}} 
##' @return matrix of parameter settings with \code{points} times \code{nrow(scenarios)} rows
##' @author Florian Klinglmueller
##'
##' @examples
##' scenarios <- expand.grid(m=seq(0,1,.1),s=seq(1,2,.1),n=100)
##' results <- simulate_batch(scenarios,100000,simulate_poly)
##' candidates <- select_results(results,'s',"result")
##' smooth_criminals <- smooth_loess(candidates,'m','s',10)
##' 
##' @export
smooth_loess <- function(scenarios,parameter,groups,points,drop=T,width=6,...){
    model <- loess(substitute(y~x,list(y=as.name(parameter),x=as.name(groups))),data=scenarios,...)
    plx <- predict(model,se=T)
    upper  <- plx$fit + plx$se*width/2
    lower  <- plx$fit - plx$se*width/2
    N <- nrow(scenarios)
    newscenarios <- scenarios[rep(1:N,each=points),]
    newscenarios[[parameter]] <- unlist(lapply(1:N,function(i) seq(lower[[i]],upper[[i]],length.out=points)))
    select(newscenarios,-matches('^result|^select'))
}

##' Plot old, selected and new candidate scenarios.
##'
##' @title plot scenarios
##' @param old Original scenarios from which selection was made
##' @param selected Selected scenarios
##' @param new New scenarios based on selection
##' @param parameter Nuisance parameter
##' @param groups Parameter of interest
##' @return scatter plot of scenarios
##' @examples
##' scenarios <- expand.grid(m=seq(0,1,.01),s=seq(1,2,.01),n=1000)
##' results <- simulate_batch(scenarios,100000,simulate_poly)
##' candidates <- select_results(results,'s',"result")
##' smooth <- smooth_loess(candidates,'m','s',10)
##' grid <- add_neighbours(candidates,'m',.2,10)
##' plot_scenarios(results,candidates,smooth,'m','s')
##' plot_scenarios(results,candidates,grid,'m','s')
##' @author Florian Klinglmueller
##' @export
plot_scenarios <- function(old,selected,new,parameter,groups){
    ylim <- c(min(min(old[[parameter]]),min(new[[parameter]])),max(max(old[[parameter]]),max(new[[parameter]])))
    plot(new[[groups]],new[[parameter]],pch=3,ylim=ylim)
    points(selected[[groups]],selected[[parameter]],pch=8)
    points(old[[groups]],old[[parameter]],pch='.',col='black')
}

## smoothy <- smooth_loess(scens,'sigma','n1',10)
## smoochy <- smooth_loess(scens,'delta','n1',10)
## grizzly <- add_neighbours(scens,.5,.01,'sigma')
## plot_scenarios(maxsim,scens,smoothy,'sigma','n1')
## plot_scenarios(maxsim,scens,grizzly,'sigma','n1')
## plot_scenarios(maxsim,scens,smoochy,'delta','n1')

