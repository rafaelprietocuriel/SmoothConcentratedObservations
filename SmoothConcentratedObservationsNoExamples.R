#### Smooth functions evaluated in concentrated observations
#### Description
####
####     code created to run a smoothed function SmoothW on a data set.
####     The data contains at least two paired variables. The function runs
####     over a continuous variable, t, with a mark, x, and for the values
####     of x, returns the value of a function f evaluated over the marks
####     and for some intervals.
####
#### Usage
#### SmoothW(DB, fun, 
####         bw = 256,
####         part = 0.2,
####         iter = 100, 
####         extra = 0.1)
####
#### Arguments
#### DB - a data frame containing two arguments:
####            t = numeric, "time" or continuous variable which is the basis of the result  
####            x = marks (maybe more than one dimension)
#### fun - function to apply to the set x
#### bw - number of observations to return
#### part - a value between 0 and 1 which is the max length of the partition for refinement. 
#### larger values of part return a coarser partition. Smaller values 
#### take longer and give less smooth values
#### iter - iterations of the code
#### extra - an additional parameter which controls the extra width of the intervals.
#### The parameter extra is the proportion that the interval is extended.


SmoothW <- function(DB, fun,
              bw = 256,
              part = 0.2,
              iter = 100,
              extra = 0.1){
  u <- range(DB$t)
  l <- u[2] - u[1]
  
  #### the points in which the function is evaluated
  tCells <- seq(from = u[1] - extra*l, 
                to = u[2] + extra*l, 
                length.out = bw)
  
  #### the matrix has as basis only NA
  NT <- matrix(rep(NA, length(tCells) * iter), nrow = iter)
  #### for the first evaluation, it gives an overall value of the function
  NT[1, ] <- fun(DB)
  for (j in 2:iter){
    #### create a random partition
    Part <- seq(from = u[1]-runif(1)*.1 * l, 
                to = u[2]+runif(1)*.1 * l, 
                by = runif(1)* part *l)
    
    #### for the partition, identigy observations inside and evaluate the function
    for (k in 1:(length(Part)-1)){
      Qm <- Part[k]
      QM <- Part[k+1]
      
      #### observations inside the window
      D <- DB$t >= Qm & DB$t  < QM
      
      #### cells inside the window
      InCell <- tCells >= Qm & tCells  < QM
      
      #### only if cells are not empty and there are observations
      if (sum(D) > 0 & sum(InCell) > 0){
        NT[j,InCell] <- fun(DB[D,])
      }
    }
  }
  Swipe <- list(EvalTime = tCells,
                ObsMean = apply(NT, 2, mean, na.rm = TRUE),
                ObsMax = apply(NT, 2, max, na.rm = TRUE),
                ObsMin = apply(NT, 2, min, na.rm = TRUE),
                Obs05 = apply(NT, 2, quantile, probs = 0.05, na.rm = TRUE),
                Obs95 = apply(NT, 2, quantile, probs = 0.95, na.rm = TRUE))  
}

