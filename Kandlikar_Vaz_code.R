
## Differential Equations in R
## An introduction with Lotka-Volterra models
## Gaurav Kandlikar (gkan@umd.edu)
## Marcel Vaz (mvaz@umd.edu)
## Prepared for ENTM 798v, Seminar in R
## Last Updated 4 December 2014

# install.packages("deSolve")
require("deSolve")

# exploring the contents of the deSolve package
ls("package:deSolve")
# this script will make use of the "ode" solver
str(ode)

# To use ode, we need to set up the function that we would like ode to solve for us.
# Define the function lv
lv <- function (time , init , params) {
# set everything up as a list to be accessed throughout the function
  with (as.list(c(time , init , params)), {

    # description of parameters to be included in params:
    # r1 = growth rate of Sp. 1; r2 = growth rate of Sp. 2
    # N = population size of Sp. 1; P = Population Size Sp. 2
    # a = competitive impact of Sp. 2 on Sp. 1; b = competitive impact of Sp 1 on Sp 2
    # K1 and K2 = carrying capacities of Sp. 1 and Sp. 2, respectively
    
    # Define the pair of ODEs:
    
    # Growth of species 1 is a function of species 2
    N1 <- ( r1 * N * (1 - (N+a*P) / K1) )
    # Growth of species 2 is a function of species 1
    N2 <- ( r2 * P * (1 - (P+b*N) / K2) )
    
    # Return a list of the current population size of Sp.1  and Sp. 2
    return (list (c (N1 , N2)))
  })
}



# Testing the function

# Randomly define a starting N and P (population sizes of the two competitors)
p <- runif(n = 2, min = 50, max = 150)
init<-c(N = floor(p[1]), P = floor (p[2])) 

# run the ODE for 100 time steps
time  <- seq (0, 100, by=1)
# define other parameters - species growth rates, carrying capacities and interspecific competition parameters
params <- c(r1 = 1, b = 1.5 , K1 = 50 , r2 = 1, a = .1 , K2 = 100)

# run the ODE 
lvout <- floor(as.data.frame(ode(func=lv,y=init,parms=params,times=time)))

str(lvout)


# We have the output of the ODE in a dataframe, but the output is best inspected visually.

# Set up the plotting area
par(mfrow=c(2,2),bg="white",oma=c(2,2,2,2))

# Panel A : write equation and parameters
plot(1:10,type="n",xlab="",ylab="",xaxt="n",yaxt="n",bty="n")
text(5,9,expression(frac(dN[1],dt) == r[1]*N[1]*bgroup("(",1-frac(N[1]+alpha*N[2],K[1]),")")),cex=1.5)
text(3,5.7,bquote(N[1]== .(floor(p[1]))),cex=1.5)
text(7,5.7,bquote(N[2]== .(floor(p[2]))),cex=1.5)
text(3,4.2,bquote(r[1]== .(params["r1"])),cex=1.5)
text(7,4.2,bquote(r[2]== .(params["r2"])),cex=1.5)
text(3,2.7,bquote(alpha== .(params["a"])),cex=1.5)
text(7,2.7,bquote(beta== .(params["b"])),cex=1.5)
text(3,1.2,bquote(K[1]== .(params["K1"])),cex=1.5)
text(7,1.2,bquote(K[2]== .(params["K2"])),cex=1.5)
  
# Panel B: plot the Zero Net Growth Isoclines based on parameters above.
plot(1,type="n",xlim=c(0,max(params["K1"],params["K2"]/params["b"])+20),
      ylim=c(0,max(params["K2"]  ,params["K1"]/params["a"])+20),
      xlab  ="Species 1",ylab="Species 2",main="ZNGIs for Sp.1 and Sp.2",
      xaxs="i",yaxs="i",cex.axis=1.5,cex.lab=1.5,cex.main=1.5)

lines(x = c(params["K1"],0),y = c(0,params["K1"]/params["a"]),lwd=2,col="blue")
lines(x = c(params["K2"]/params["b"], 0) ,y = c(0,params["K2"]),lwd=2,col="red")
# plot the starting population size onto the ZNGIs
points(x=p[1],y=p[2],cex=2,pch=20)

# Panel C: Population 1 vs Population 2
plot(lvout$P~lvout$N,type="o",col="blue",xlim=c(0,max(lvout$N)+20),
     ylim=c(0,max(lvout$P)+20), main = "Species 1 vs Species 2",xlab="Species 1", 
     ylab="Species 2",cex.axis=1.5,cex.lab=1.5,cex.main=1.5)

# Panel D: Population 1 and Population 2 vs time
plot(lvout$N~time,ylim=c(0,max(lvout$N)+20),ylab="Population size",type="l",lwd=2,col="blue",cex.axis=1.5,cex.lab=1.5,main="Population size vs time",cex.main=1.5)
points(lvout$P~time,col="red",type="l",lwd=2)

## ------------------------------------------------------------------------
 
## Replicating the process with Lotka-Volterra predator-prey interactions
## Set up a new function, lvpp
lvpp <- function(pp.time,pp.init,pp.params,prey_K=FALSE) {
  with (as.list(c(pp.time,pp.init,pp.params)), {
  
    # Parameters
    # N = prey population size; P = predator population size
    # r = intrinsic growth rate of prey
    # a = predation efficiency
    # b = conversion efficiency of prey into predator
    # d = intrinsic dseath rate of predator
    # prey_k = carrying capacity for prey; only used if user-defined
    
     if (exists("prey_k")) {
       dNdt <- ((r*N)*(1-(N/prey_k)))- (a*N*P) 
     }
      else {
        dNdt <- (r*N) - (a*N*P)
      }
    
      dPdt <- (b*a*N*P) - (d*P)
    
      return(list(c(dNdt,dPdt)))
  })
}


## ------------------------------------------------------------------------
# Set ODE parameters
pp.time <-seq(0,3000,by=.1)
pp.params <- c(r=0.1,a=0.005,d=0.01,b=.1,prey_k=500)

# Set pp.initial population sizes of Prey (N) and Predator (P)
pp.init = c(N=25,P=10)
# Run the ODE

lvppout<-floor(as.data.frame(ode(func=lvpp,y=pp.init,parms=pp.params,times=pp.time)))
str (lvppout)

# Set up plotting parameters
par(mfrow=c(1,2),bg="white")

# Panel A: Plot P vs N; draw in the starting N and P parameters, draw in the ZNGIs
plot(lvppout$P~lvppout$N,ylim=c(0,max(lvppout$P)+20),type="l",xlab="Prey population size",ylab="Predator population size")
points(x=pp.init["N"],y=pp.init["P"],col="red",pch=19)
abline(v=pp.params["d"]/(pp.params["b"]*pp.params["a"]))
# abline(h=pp.params["r"]/pp.params["a"])
abline (b=-(pp.params["r"]/(pp.params["prey_k"]*pp.params["a"])), a = pp.params["r"]/pp.params["a"])


# Panel B: Plot N & P vs pp.time
plot(lvppout$N~pp.time,type="l",xlab="pp.time",ylab="Population Size",ylim=c(0,max(max(lvout$N),max(lvout$P))+50))
points(lvppout$P~pp.time,col="red",type="l")
legend(x="topright",col=c("black","red"),lty=1,legend=c("Prey","Predator"),bty="n")
