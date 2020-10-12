# Code for reunning stochstic compartmental models





rand_vect_cont <- function(N, M, sd = 1) {
  vec <- rnorm(N, M/N, sd)
  vec / sum(vec) * M
}

#####################STOCHASTIC FIXED TIMESTEP SIR MODELS#################
#Runs a simple stochastic SIR model.
#Parameters
#  beta - the force of infection
#  gamma - the recovery rate
#  initial.state - the initial state of the system.
#               MUST BE WHOLE NUMBERS!
#  time.step -  what time step should we use
#  freq.dependent - should we use a frequency or densitiy
#             dependent model
#  final.only - speeds things up a little if we are only interested in
#           the final state of the epidemic
#  import - a vector showing number of importations into the sink
#  import.step.size - time step of the importations
#  import.step.size / step.size has to be an integer, same for step.size/import.step.size
#  max.time - how many days of data do we ask the model to generate
#Returns
#   a matrix with 4 columns: t, S, I, and R
runStochasticSEIR <- function(beta=.001,
                              sigma=.5,
                              gamma=.05,
                              initial.state= c(S=999, E=1, I1=0, I2=0, I3=0, R=0, I_trk=0, E_trk=0),
                              step.size = 1/4,
                              freq.dependent=FALSE,
                              final.only=FALSE,
                              max.time = NULL){

  #If the model is frequency dependent we modify beta
  #based on the total populations size
  beta.divisor <- ifelse(freq.dependent,
                           initial.state["S"] + 
                           initial.state["E"] +
                           initial.state["I1"] + 
                           initial.state["I2"] + 
                           initial.state["I3"] + 
                           initial.state["R"],
                           1)

  #create the parameter vector.
  param <- c(beta=beta/beta.divisor, sigma=sigma, gamma=gamma)

    
  #Since we are not using a fancy solver we will need to run this on our own.
  #note that the epidemic ends once there are no more susceptibles.
  t <- 0
  y <- initial.state

  sir.output <- matrix(ncol=9, nrow=1)
  colnames(sir.output) <- c("time", "S","E","I1","I2","I3","R","I_trk","E_trk")
  sir.output[1,] <- c(t,y)
  
  
  while ((y["E"]+y["I1"]+y["I2"]+y["I3"])>0 & t<=max.time ) {

    t <- t + step.size 

    delta <- stochastic.dx.dt(step.size, y, param)
    
    # how y changes for each step.size increment
    y <- y + delta
    
    #trick to speed up the code
    if (!final.only) {
      sir.output <- rbind(sir.output, c(t,y))
    }
  }

  if(final.only) {
    sir.output[1,] <- c(t,y)
  }

  return(sir.output)

}




#The dx.dt method for a stochastic SIR using simple
#rates (e.g., a constant probability of being moved
#from one comparment to the other)
#
#Parameters:
#  step.size -  the size of the step being taken compared
#         to that used to specify the parameters
#
#  y - the current state of the system.
#  param - the parameters of the system
#
#Returns:
#  a list representing the changes in each compartment of the system
stochastic.dx.dt <- function(step.size, y, param) {
  
  #calculate the probability of infection and recovery
  #in this time step
  p.expose  <- 1-exp(-step.size*param["beta"]*(y["I1"]+y["I2"]+y["I3"]))    
  p.infect  <- 1-exp(-step.size*param["sigma"])        
  p.recover <- 1-exp(-step.size*param["gamma"])      

  #Do our random stuff so we know how much to change.
  # The number of exposed cases follows a binomial
  # distribution where N=S and p=p.exposed
  expose.cases <- rbinom(1, y["S"], p.expose)
  
  # The number of incident cases follows an erlang/gamma
  # to get that, we use 3 compartments with binomial
  # distribution where N=E and p=p.infect, 
  incident1.cases <- rbinom(1, y["E"], p.infect)

  # The number of recovered cases follows a binomial
  # distribution where N=I and p=p.recover
  incident2.cases <- rbinom(1, y["I1"], p.recover)
  incident3.cases <- rbinom(1, y["I2"], p.recover)
  
  # The number of recovered cases follows a binomial
  # distribution where N=I and p=p.recover
  recovered.cases <- rbinom(1, y["I3"], p.recover)
  
  
  #Find the deltas for our compartments
  dS  <- -expose.cases
  dE  <- expose.cases - incident1.cases
  dI1 <- incident1.cases - incident2.cases
  dI2 <- incident2.cases - incident3.cases
  dI3 <- incident3.cases - recovered.cases
  dR  <- recovered.cases 
  dI_trk <- incident1.cases
  dE_trk <- expose.cases
  
  #Those susceptibles move to the recovered compartment
  return(c(dS,dE,dI1,dI2,dI3,dR,dI_trk,dE_trk))
}




