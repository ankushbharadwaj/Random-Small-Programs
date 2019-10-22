#bipolar disorder function
#this code represents tom chou's initial differential equations for bipolar disorder
#the values for all the parameters are lifted from tom chou's code

#set parameters equal to constant values
m_scale <- 10*(2/3)
t_scale <- (100/30)*(6.5/8)
parameters <- c(eta_m = m_scale/t_scale, eta_v = 1/t_scale , f = 2/m_scale, r = 2, k = 1/t_scale, k_3 = 3/(t_scale*m_scale^2))

#set the initial values of the state variables
state <- c(M = 0, V = 0)

#this function represents the ordinary differential equations from tom chou's paper
bpd<-function(t, state, parameters) {
  with(as.list(c(state,parameters)),{ 
    dM <- eta_m*(f*M + r - V) - k*M - k_3*M^3
    dV <- eta_v*(f*M + r - V)
    
    list(c(dM,dV))
  })
}

#arbitrarily decided to simulate the model from t=0 to t=72, where time is measured as hours
times <- seq(0,72,by=.01)

#call the package that includes the function that we will use to solve the set of differential equations
library(deSolve)

#solve the differential equations and set them equal to a variable
out1 <- ode(y = state, times = times, func = bpd, parms = parameters)

#############################################################################

#this file contains the code for the time dependent function for eta_m

#declare the function
eta_m_new<-function(t) {
  values <- c(rep.int(0,length(t))) 
  #first, create a vector values that repeats the value 0 such that values is a vector of equal length to the input t full of 0
  for (i in seq_along(t)) { 
    #loop through each element of t essentially 
    if (t[i] > 12) { 
      #if the value of t is greater than 12, we are using modulus to find by how much t is greater than 12 and using this value to calculate the value for eta_m at that time t
      t <- t %% 12
      values[i] <- c(((m_scale/t_scale))*(1+(1/18))^((t[i])))
      #ith element of values equals eta_m as a function of t, which is modeled here
    }
    else if (t[i] <= 12) { 
      #this conditional is for the case that the time value is less than or equal to 12, in which case we don't have to alter t before calculating the value of eta_m at that point of time 
      values[i] <- c(((m_scale/t_scale))*(1+(1/18))^((t[i])))
    }
    
  } 
  return(values) 
  #return the vector that includes the values of eta_m that corresponds to specific times
}

###########################################################################

#bipolar disorder function
#this is the one that includes the time dependent eta_m value rather than setting it equal to a value as a parameter
#most of this code is the same as their original function's code
#comments on the lines that have been altered 

m_scale <- 10*(2/3)
t_scale <- (100/30)*(6.5/8)
parameters <- c(eta_v = 1/t_scale , f = 2/m_scale, r = 2, k = 1/t_scale, k_3 = 3/(t_scale*m_scale^2))
#in the above line, we removed eta_m from the list of parameters, because it is not a parameter anymore

state <- c(M = 0, V = 0)

#in the function below, we change the mention of eta_m to eta_m(t), which calls the time dependent function for eta_m that we wrote whenever ODE is used to solve this set of differential equations
bpd<-function(t, state, parameters) {
  with(as.list(c(state,parameters)),{ 
    dM <- eta_m_new(t)*(f*M + r - V) - k*M - k_3*M^3
    dV <- eta_v*(f*M + r - V)
    
    list(c(dM,dV))
  })
}

times <- seq(0,72,by=.01)

library(deSolve)

#solve the differential equations and set them equal to a different variable
out <- ode(y = state, times = times, func = bpd, parms = parameters)

##############################################################################

#this file includes the code used in conjunction with each other to graph the results from solving the set of ordinary differential equations

#plots time vs mood
plot(times,out[,2], type = "l", main = "Time vs Mood", xlab = "Time - hours", ylab = "Mood", xlim = c(0,72), ylim = c(-10,10), col = "red")
par(new=T)
plot(times,out1[,2], type = "l", main = "Time vs Mood", xlab = "Time - hours", ylab = "Mood", xlim = c(0,72), ylim = c(-10,10), col = "blue") 
par(new=T)
legend(0,10,legend=c("Time Dependent eta_m", "Constant eta_m"), col = c("red","blue"),lty=1:2,cex=.8)

#plots time vs expectation
plot(times,out[,3], type = "l", main = "Time vs Expectation", xlab = "Time - hours", ylab = "Expectations", xlim = c(0,72), ylim = c(0,4), col = "red")
par(new=T)
plot(times,out1[,3], type = "l", main = "Time vs Expectation", xlab = "Time - hours", ylab = "Expectations", xlim = c(0,72), ylim = c(0,4), col = "blue")
par(new=T)
legend(0,4,legend=c("Time Dependent eta_m", "Constant eta_m"), col = c("red","blue"),lty=1:2,cex=.8)

#phase diagram
plot(out[,2],out[,3], pch = ".", main = "Mood vs Expectation", xlab = "Mood", ylab = "Expectation", xlim = c(-6,6), ylim = c(0,3), col = "red", type = "l")
par(new=T)
plot(out1[,2],out1[,3], pch = ".", main = "Mood vs Expectation", xlab = "Mood", ylab = "Expectation", xlim = c(-6,6), ylim = c(0,3), col = "blue", type = "l")
legend(-6,3,legend=c("Time Dependent eta_m", "Constant eta_m"), col = c("red","blue"),lty=1:2,cex=.8)

#this plots time vs eta_m
plot(t, eta_m_new(t), type = "l", main = "Time vs eta_m_new", xlab = "Time - hours", ylab = "eta_m", xlim = c(0,50), ylim =c(0,5), col = "red")
par(new=T)
plot(t, c(rep.int((m_scale/t_scale),length(t))) , col="blue", xlab = "Time - hours", ylab = "eta_m", xlim = c(0,50), ylim =c(0,5),type="l")
legend(0,1,legend=c("Time Dependent eta_m", "Constant eta_m"), col = c("red","blue"),lty=1:2,cex=.8)

