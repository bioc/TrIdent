#' No pattern pattern-match
#'
#' A horizontal line at the mean or median coverage should be an optimal
#' pattern-match if the contig read coverage displays no sloping or block
#' patterns
#'
#' @param viralSubset A subset of the read coverage pileup that pertains only to
#'   the contig currently being assessed
#' @param searchMethod Search method to use. Either "grid" for the original grid
#'   search or "direct" for DIRECT global optimization.
#' @param DirectMaxEval Maximum number of DIRECT evaluations to make.
#' @param globalLocal Use global or local DIRECT search. Default is local. 
#' @return List
#' @keywords internal
noPattern <- function(viralSubset, searchMethod, DirectMaxEval,globalLocal) {
    if(searchMethod =="direct"){
     bestMatchInfo <- optimNoPattern(viralSubset, DirectMaxEval, globalLocal)
    }else if (searchMethod=="grid"){
  pattern1 <- rep(median(viralSubset[, 2]), nrow(viralSubset))
  pattern2 <- rep(mean(viralSubset[, 2]), nrow(viralSubset))
  diff1 <- mean(abs(viralSubset[, 2] - pattern1))
  diff2 <- mean(abs(viralSubset[, 2] - pattern2))
  diff <- ifelse(diff1 < diff2, diff1, diff2)
  value <- ifelse(diff1 < diff2, median(viralSubset[, 2]),
    mean(viralSubset[, 2])
  )
  bestMatchInfo <-
    list(
      diff,
      value,
      nrow(viralSubset),
      "NA",
      1,
      nrow(viralSubset),
      "NoPattern"
    )
    }
  return(bestMatchInfo)
}



#' Function for applying DIRECT global optimization algorithm
#'
#' Set upper and lower limits for DIRECT and apply it to the noPattern pattern
#'
#' @param viralSubset A subset of the read coverage pileup that pertains only to
#'   the contig currently being assessed
#' @param DirectMaxEval Maximum number of DIRECT evaluations to make.
#' @param globalLocal Use global or local DIRECT search. Default is local. 
#' @return List
#' @keywords internal
optimNoPattern <- function(viralSubset, DirectMaxEval, globalLocal){
    maxReadCov <- max(viralSubset[, 2])
    minReadCov <- min(viralSubset[, 2])
    wrapper <- function(dims){
        pattern <- rep(dims, nrow(viralSubset))
        diff <- mean(abs(viralSubset[, 2] - pattern))
        return(diff)
    }
    optim <- if(globalLocal=="local"){
            nloptr::directL(wrapper, minReadCov, maxReadCov, control = list(xtol_rel = 1e-4, maxeval = DirectMaxEval))}
            else{
            nloptr::direct(wrapper, minReadCov, maxReadCov, control = list(xtol_rel = 1e-4, maxeval = DirectMaxEval))}
    pattern <- rep(optim$par, nrow(viralSubset))
    diff <- mean(abs(viralSubset[, 2] - pattern))
    bestMatchInfo <-
        list(
            diff,
            optim$par,
            nrow(viralSubset),
            "NA",
            1,
            nrow(viralSubset),
            "NoPattern"
        )
    return(bestMatchInfo)
}

