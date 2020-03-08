# 3-(a)
# translates trajectory to begin with time zero at the origin
trajectory_set <- function(x){
  # x is a n*3 marix representing the trajectory(x,y,t)
  n <- nrow(x)
  i <- 2
  while (i < n+1) {
    x[i,1] <- x[i,1]-x[1,1]
    x[i,2] <- x[i,2]-x[1,2]
    x[i,3] <- x[i,3]-x[1,3]
    i <- i+1
  }
  x[1,1] <- 0
  x[1,2] <- 0
  x[1,3] <- 0
  return(x)
}
# 3-(b)
# Computes the angle cita formed by the secant line
# The secant line connects the origin and the final position
angle_compute <- function(x){
  # x is a n*3 marix representing the trajectory(x,y,t)
  n <- nrow(x)
  a <- x[n,1]-x[1,1]
  b <- x[n,2]-x[1,2]
  a <- as.numeric(a)
  b <- as.numeric(b)
  longth <- sqrt(a^2+b^2)
  cita <- 0
  if(a==0 && b==0){
    cita <- 0
  } else if(a>=0 && b>=0){
    cita <- asin(b/longth)
  } else if(a<0 && b>=0){
    cita <- pi-asin(b/longth)
  } else if(a<0 && b<0){
    cita <- -pi-asin(b/longth)
  } else {
    cita <- asin(b/longth)
  }
  return(cita)
}
# 3-(c)
# Rotate the (x, y) coordinates of a trajectory 
# So that the final point lies along the positive x-axis
trajectory_rotate <- function(x){
  # x is a n*3 marix representing the trajectory(x,y,t)
  n <- nrow(x)
  t <- angle_compute(x)
  i <- 1
  a1 <- 0
  a2 <- 0
  if(t>=0){
    while (i<=n) {
      a1 <- x[i,1]*cos(t)+x[i,2]*sin(t)
      a2 <- x[i,2]*cos(t)-x[i,1]*sin(t)
      x[i,1] <- a1
      x[i,2] <- a2
      i <- i+1
    }
  }
  if(t<0){
    while (i<=n) {
      a1 <- x[i,1]*cos(abs(t))-x[i,2]*sin(abs(t))
      a2 <- x[i,2]*cos(abs(t))+x[i,1]*sin(abs(t))
      x[i,1] <- a1
      x[i,2] <- a2
      i <- i+1
    }
  }
  return(x)
}
# 3-(d)
# normalizes an n*3 trajectory matrix 
# to begin at the origin and end on the positive x-axis.
trajectory_normalize <- function(x){
  # x is a n*3 marix representing the trajectory(x,y,t)
  x = trajectory_set(x)
  x = trajectory_rotate(x)
  return(x)
}
# 3-(e)
# Compute some metrics describing the trajectory's curvature
metrics_calculate <- function(x){
  # x is a n*3 marix representing a normalized trajectory 
  # Compute the total distance
  x <- trajectory_normalize(x)
  distance <- 0
  n <- nrow(x)
  i<- 1
  while (i<n) {
    distance <- distance + sqrt((x[i+1,1]-x[i,1])^2+(x[i+1,2]-x[i,2])^2)
    i <- i+1
  }
  # compute the maximum absolute deviation
  max_deviation <- 0
  i <- 1
  while (i<n) {
    if(max_deviation<abs(x[i,2])){
      max_deviation <- abs(x[i,2])
    }
    i <- i+1
  }
  # compute the average absolute deviation
  i <- 2
  sum_deviation <- 0
  while (i<n) {
    sum_deviation <- sum_deviation + abs(x[i,2])
    i <- i+1
  }
  aver_deviation <- 0
  aver_deviation <- sum_deviation/n
  # Compute the absolute area under the trajectory and secant line 
  i <- 1
  area <- 0
  i <- 1
  for (i in 1:(n-1)){
    disi=x[(i+1),1]-x[i,1]
    area=area+(abs(x[i,2])+abs(x[i+1,2]))*disi/2
  }
  # return the result that we compute above
  metrics <- c(distance,max_deviation,aver_deviation,area)
  names(metrics) <- c('tot_dist','max_abs_dev','avg_abs_dev','AUC')
  return(metrics)
}