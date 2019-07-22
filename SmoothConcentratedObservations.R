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


### Example 1
### It constructs a sample dataset, concentrated over t
n <- 200
t <- runif(n)^8
x <- 16 * t + rnorm(n)
DB <- data.frame(x = x, t = t)
g <- function(DB){return(mean(DB$x))}

Result <- SmoothW(DB, fun = g,
                  part = 0.1,
                  bw = 1000,
                  iter = 200,
                  extra = 0.01)
plot(DB$t, DB$x)
polygon(c(Result$EvalTime, rev(Result$EvalTime)),
        c(Result$ObsMin, rev(Result$ObsMax)), 
        col = rgb(.1, .2, .7, .2))
polygon(c(Result$EvalTime, rev(Result$EvalTime)),
        c(Result$Obs05, rev(Result$Obs95)), 
        col = rgb(.1, .2, .7, .4))
points(Result$EvalTime, Result$ObsMean, 
       type = "l", col = 2, lwd = 4)


### Example 2
### the function returns a regression coeffiicient for observations within range
n <- 1000
t <- runif(n)^3
x1 <- 16 * t + rnorm(n)
x2 <- -8 * x1 + 1* rnorm(n)
### the function g gives us, in this case, the coefficient of a linear model
g <- function(DB){ 
  m <- lm(DB$x2 ~ DB$x1)
  return(m$coefficients[2])
  }
DB <- data.frame(x1 = x1, x2 = x2, t = t)

Result <- SmoothW(DB, fun = g,
                  part = 0.4,
                  bw =  100,
                  iter = 200)
plot(DB$x1, DB$x2)

### notice that since the function g gives the coefficient of a linear model
plot(Result$EvalTime, Result$ObsMean, 
     ylim = c(-11, -5))
polygon(c(Result$EvalTime, rev(Result$EvalTime)),
        c(Result$ObsMin, rev(Result$ObsMax)), 
        col = rgb(.1, .2, .7, .2))
polygon(c(Result$EvalTime, rev(Result$EvalTime)),
        c(Result$Obs05, rev(Result$Obs95)), 
        col = rgb(.1, .2, .7, .4))
points(Result$EvalTime, Result$ObsMean, 
       type = "l", col = 2, lwd = 4)


### Example 3
### This example uses the package the swiss dataset, described here:
### https://stat.ethz.ch/R-manual/R-devel/library/datasets/html/swiss.html
DB <- data.frame(t = swiss$Catholic,
                 x = swiss$Education)
g <- function(DB){
  m <- lm(DB$x ~ DB$t)
  return(m$coefficients[2])
}

Result <- SmoothW(DB, fun = g,
                  part = 0.2,
                  bw =  100,
                  iter = 1000)

### notice that since the function g gives the coefficient of a linear model
plot(Result$EvalTime, Result$ObsMean, 
     ylim = c(-9, 7))
polygon(c(Result$EvalTime, rev(Result$EvalTime)),
        c(Result$ObsMin, rev(Result$ObsMax)), 
        col = rgb(.1, .2, .7, .2))
polygon(c(Result$EvalTime, rev(Result$EvalTime)),
        c(Result$Obs05, rev(Result$Obs95)), 
        col = rgb(.1, .2, .7, .4))
points(Result$EvalTime, Result$ObsMean, 
       type = "l", col = 2, lwd = 4)
