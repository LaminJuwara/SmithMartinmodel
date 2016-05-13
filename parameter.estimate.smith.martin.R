##########################################################################################################
#   parameter estimation for the Smith-Martin model
##########################################################################################################

#model definition
smith.martin.model <- function(t,Y,parms){
  with(as.list(c(parms,Y)), {
    if(t < deltat){
      lag1 <- init.vec[1:(ngen+1)]
    }
    else{
      lag1 <- lagvalue(t-deltat)[1:(ngen+1)]
    }
    temp <- exp(-d_0*deltat)
    dN <- rep(0,length(Y))
    dN[1] <- -(lambda_0 + d_0)*Y[1]
    dN[2] <- 2*lambda_0*temp*lag1[1] - (lambda+d)*Y[2]
    dN[3:(ngen+1)] <- 2*lambda*temp*lag1[2:ngen] - (lambda+d)*Y[3:(ngen+1)]
    
    dN[ngen+2] <- -lambda_0*temp*lag1[1] - d_0*Y[ngen+2]
    dN[ngen+3] <- 2*lambda_0*temp*lag1[1] - lambda*temp*lag1[2] - d*Y[ngen+3]
    dN[(ngen+4):(2*ngen+2)] <- 2*lambda*temp*lag1[2:ngen] - lambda*temp*lag1[3:(ngen+1)] - d*Y[(ngen+4):(2*ngen+2)]
    list(dN)
  })
}
#######################################################################################
#this function simulates the Smith-Martin model
#and returns the number of cells in each gen
#at each time point found in the vector t_pts.
simulate.smith.martin.model <- function(parms){
  out <- dede(init.vec,t_pts,smith.martin.model,parms)
  return (out[,(ngen+3):(2*ngen+3)])
}

#example run of the model
require(deSolve)
#cells in generations ranging from 0 to ngen will be considered

ngen <- 20
init.vec <- rep(0,2*ngen+2)
#number of cells in A phase in i'th generation is denoted NAi
#total number of cells in i'th generation is denoted Ni
names(init.vec) <- paste(rep(c('NA','N'),each=ngen+1),0:ngen,sep='')
#start with 1e6 cells in 0'th generation
init.vec[c('NA0','N0')] <- 50000
#define parameters of the model
parms <- c(lambda_0=1e-2, lambda=1e-2, d_0=1e-5, d=1e-5, deltat=3)
#set the time points at which data should be collected
t_pts <- seq(0,120,12)
#run the simulation
out <- simulate.smith.martin.model(parms)
round(out)
####################################################################################
#example run of the model

# define default parameters and time points
parms1<-c(lambda_0=1e-2, lambda=1e-2, d_0=1e-5, d=1e-5, deltat=3)
t_pts <- seq(0,120,12)
##################################################################
# Optimization of Smith-Martin model
##################################################################
##################################################################

# read the simulated data generated from agent based model
df <- read.table("./Desktop/datasets/agent_based_sim_output_1.txt",header = TRUE)

sum.cols <- function(idx, data) return (rowSums(data[,idx]))

# reshape data to have total number of cells in each division class in a single column
dataset <- matrix( rep(0, (length(t_pts)*(ngen+1))), nrow = length(t_pts))
dataset[1,1] <- init.vec[1]
dataset[-1,] <- sapply(split(2:43,rep(1:21,each=2)),sum.cols,df)

# generates a candidate dataset using the parameter to simulate model
f2.sm <- function(parms1){
  out <- dede(init.vec,t_pts,smith.martin.model,parms1)
  return (sapply(split(2:43,rep(1:21,each=2)),sum.cols,out))
}

# fits the simulated data from agent based model to model simulation
f1.sm <- function(parms1, dataset){
  sim.dataset <- f2.sm(parms1)
  return(sum((sim.dataset - dataset)^2))
}
# optim function determines the best parameter canditate
out.optim <- optim(parms1, f1.sm, gr=NULL, dataset)
print(out.optim)



#######################################################################
# plots Smith-Martin
#######################################################################
model.sim <- round(out)
agent.sim <- dataset
time <- t_pts
######################## plots #####################################
png(file="gen0.png",width=400,height=450,res = 120)
plot(c(0,120), c(0,50000), type = "n", xlab = "Hours", ylab = "Count",
     main = "Gen. 0",font.main=1)
# set up the axis
lines(time, model.sim[,1] , col="red", lwd=1.5)
lines(time, agent.sim[,1] , col="blue", lwd=1.5)
legend(40,49000, c('S-M', 'A-B'),
       lwd=c(1.5,1.5),col=c("red","blue"),lty=c(1,1))

dev.off()

png(file="gen1.png",width=400,height=450,res = 120)
plot(c(0,120), c(0,50000), type = "n", xlab = "Hours", ylab = "Count",
     main = "Gen. 1",font.main=1)
# set up the axis
lines(time, model.sim[,2] , col="red", lwd=1.5)
lines(time, agent.sim[,2] , col="blue", lwd=1.5)
legend(40,49000, c('S-M', 'A-B'),
       lwd=c(1.5,1.5),col=c("red","blue"),lty=c(1,1))

dev.off()
png(file="gen2.png",width=400,height=450,res = 120)
plot(c(0,120), c(0,50000), type = "n", xlab = "Hours", ylab = "Count",
     main = "Gen. 2",font.main=1)
# set up the axis
lines(time, model.sim[,3] , col="red", lwd=1.5)
lines(time, agent.sim[,3] , col="blue", lwd=1.5)
legend(40,49000, c('S-M', 'A-B'),
       lwd=c(1.5,1.5),col=c("red","blue"),lty=c(1,1))

dev.off()

png(file="gen3.png",width=400,height=450,res = 120)
plot(c(0,120), c(0,50000), type = "n", xlab = "Hours", ylab = "Count",
     main = "Gen. 3",font.main=1)
# set up the axis
lines(time, model.sim[,4] , col="red", lwd=1.5)
lines(time, agent.sim[,4] , col="blue", lwd=1.5)
legend(40,49000, c('S-M', 'A-B'),
       lwd=c(1.5,1.5),col=c("red","blue"),lty=c(1,1))

dev.off()

png(file="gen4.png",width=400,height=450,res = 120)
plot(c(0,120), c(0,50000), type = "n", xlab = "Hours", ylab = "Count",
     main = "Gen. 4",font.main=1)
# set up the axis
lines(time, model.sim[,5] , col="red", lwd=1.5)
lines(time, agent.sim[,5] , col="blue", lwd=1.5)
legend(40,49000, c('S-M', 'A-B'),
       lwd=c(1.5,1.5),col=c("red","blue"),lty=c(1,1))

dev.off()

png(file="gen5.png",width=400,height=450,res = 120)
plot(c(0,120), c(0,50000), type = "n", xlab = "Hours", ylab = "Count",
     main = "Gen. 5",font.main=1)
# set up the axis
lines(time, model.sim[,6] , col="red", lwd=1.5)
lines(time, agent.sim[,6] , col="blue", lwd=1.5)
legend(40,49000, c('S-M', 'A-B'),
       lwd=c(1.5,1.5),col=c("red","blue"),lty=c(1,1))

dev.off()
png(file="gen6.png",width=400,height=450,res = 120)
plot(c(0,120), c(0,50000), type = "n", xlab = "Hours", ylab = "Count",
     main = "Gen. 6",font.main=1)
# set up the axis
lines(time, model.sim[,7] , col="red", lwd=1.5)
lines(time, agent.sim[,7] , col="blue", lwd=1.5)
legend(40,49000, c('S-M', 'A-B'),
       lwd=c(1.5,1.5),col=c("red","blue"),lty=c(1,1))

dev.off()
png(file="gen7.png",width=400,height=450,res = 120)
plot(c(0,120), c(0,50000), type = "n", xlab = "Hours", ylab = "Count",
     main = "Gen. 7",font.main=1)
# set up the axis
lines(time, model.sim[,8] , col="red", lwd=1.5)
lines(time, agent.sim[,8] , col="blue", lwd=1.5)
legend(40,49000, c('S-M', 'A-B'),
       lwd=c(1.5,1.5),col=c("red","blue"),lty=c(1,1))

dev.off()
png(file="gen8.png",width=400,height=450,res = 120)
plot(c(0,120), c(0,50000), type = "n", xlab = "Hours", ylab = "Count",
     main = "Gen. 8",font.main=1)
# set up the axis
lines(time, model.sim[,9] , col="red", lwd=1.5)
lines(time, agent.sim[,9] , col="blue", lwd=1.5)
legend(40,49000, c('S-M', 'A-B'),
       lwd=c(1.5,1.5),col=c("red","blue"),lty=c(1,1))

dev.off()
#########################################################################################



