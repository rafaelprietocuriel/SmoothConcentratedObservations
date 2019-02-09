#### Smooth functions evaluated in concentrated observations
#### Description
####
####     code created to run a smoothed function f over concentrated 
####     but continuous events, with a mark, and returns the value of f
####     evaluated over the marks and for some intervals
####
#### Usage
#### SmoothW(data, fun, 
####         bw = 256,
####         part = 0.2,
####         iter = 100, 
####         bw = 0.1)
####
#### Arguments
#### data - data frame containing two numeric arguments:
####            t = time  
####            x = marks
#### fun - function to apply to the set x
#### bw - number of observations to return
#### part - a value between 0 and 1 which is the max length of the partition for refinement. 
#### larger values of part return a coarser partition. Smaller values 
#### take longer and give less smooth values
#### iter - iterations of the code




f <- function(DB, fun,
              bw = 256,
              part = 0.2,
              iter = 100){
  u <- range(DB$time)
  l <- u[2] - u[1]
  #### the points in which the function is evaluated
  tCells <- seq(from = u[1] - .1*l, 
                to = u[2] + .1*l, 
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
      D <- DB$time >= Qm & DB$time  < QM
      
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




### test set
n <- 200
t <- runif(n)^8
x <- 16 * t + rnorm(n)
DB <- data.frame(x = x, time = t)

g <- function(DB){return(mean(DB$x))}
Res <- f(DB, fun = g, 
         part = 0.1,
         bw = 1000, iter = 200)
plot(DB$time, DB$x)
points(Res$EvalTime, Res$ObsMean, type = "l", col = 3, lwd = 2)
points(Res$EvalTime, Res$ObsMax, type = "l", col = 2)
points(Res$EvalTime, Res$ObsMin, type = "l", col = 2)
points(Res$EvalTime, Res$Obs05, type = "l", col = 2)
points(Res$EvalTime, Res$Obs95, type = "l", col = 2)




### test set 2
### the function returns a regression coeffiicient for observations within range
n <- 200
t <- runif(n)^8
x1 <- 16 * t + rnorm(n)
x2 <- -8 * x1 + 2* rnorm(n)
g <- function(DB){ 
  m <- lm(DB$x2 ~ DB$x1)
  return(m$coefficients[2])
  }
DB <- data.frame(x1 = x1, x2 = x2, time = t)


Res <- f(DB, fun = g)
plot(DB$time, DB$x1)
points(Res$EvalTime, Res$ObsMean, type = "l", col = 2, lwd = 2)
points(Res$EvalTime, Res$ObsMax, type = "l", col = 2)
points(Res$EvalTime, Res$ObsMin, type = "l", col = 2)
points(Res$EvalTime, Res$Obs05, type = "l", col = 2)
points(Res$EvalTime, Res$Obs95, type = "l", col = 2)

