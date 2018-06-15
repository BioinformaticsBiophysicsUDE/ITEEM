if(!require("deSolve")) install.packages("deSolve", repos="http://ftp5.gwdg.de/pub/misc/cran/")
if(!require("gplots")) install.packages("gplots", repos="http://ftp5.gwdg.de/pub/misc/cran/")

GLVn <- function(t,pop,parameters) {
  with(as.list(c(pop, parameters)),{
    rp <- rep*pop
    dp <- (rp*(Ns - sum(pop)) + rp*(Interaction%*%pop) - (t(Interaction)%*%rp)*pop)/Ns - death*pop
    list(dp)
  })
}

GLV1 <- function(t,pop,parameters) {
  with(as.list(c(pop, parameters)),{
    dp <- rep*pop*(1-pop/Ns) - death*pop
    list(dp)
  })
}

# this function deletes the strains that their population is less than the threshold
DeleteRare <- function(pop,Interaction,rep,Nst,N.threshold) {
  if (Nst > 1) {
    ids <- which(pop>=N.threshold)
    Nst <- length(ids)
    return(list(pop[ids],Interaction[ids,ids],rep[ids],Nst))
  }
  else 
    return(list(pop,Interaction,rep,Nst))
}

# this function scales the population to Ns when it is higher than this value  
LimitPopulation <- function(pop,Ns) {
  if (sum(pop)>Ns) {
    return(pop*Ns/sum(pop))
  } else {
    return(pop)
  }
}

# this function produces the mutants
Producemutants <- function(pop,mutant,Interaction,rep,Nst,m,num.mutant,intra.sp.int) {
  
  Nst.temp <- Nst
  if (Nst > 1) {
    for (i in which(mutant>0)) {
      # here the noise is added to the parent strains to produce their mutants
      # in order to speed up the calculations each subsection of the interaction matrix is produced separately
      temp1 <- matrix(rep(Interaction[i,],mutant[i]),ncol = Nst.temp,byrow = T)
      temp1 <- temp1 + matrix(rnorm(n = Nst.temp*mutant[i], mean = 0, sd = m),ncol = Nst.temp)
      ids <- which(temp1>1)
      temp1[ids] <- 2-temp1[ids]
      ids <- which(temp1<0)
      temp1[ids] <- -temp1[ids]
      temp2 <- -t(temp1)+1
      temp3 <- matrix(rep(temp2[1,],mutant[i]),ncol = mutant[i],byrow = T) + matrix(rnorm(n = mutant[i]*mutant[i], mean = 0, sd = m),ncol = mutant[i])
      ids <- which(temp3>1)
      temp3[ids] <- 2-temp3[ids]
      ids <- which(temp3<0)
      temp3[ids] <- -temp3[ids]
      temp4 <- -t(temp3)+1
      temp3 <- temp3*upper.tri(temp3, diag = FALSE)
      temp4 <- temp4*lower.tri(temp4, diag = FALSE)
      temp5 <- temp3 + temp4
      diag(temp5) <- intra.sp.int
      Interaction <- rbind(cbind(Interaction,temp2),cbind(temp1,temp5))
      Nst.temp <- Nst.temp + mutant[i]
    }
  } else {
    temp1 <- matrix(rep(Interaction,mutant),ncol = Nst,byrow = T)
    temp1 <- temp1 + matrix(rnorm(n = Nst*mutant, mean = 0, sd = m),ncol = Nst)
    temp2 <- -t(temp1)+1
    temp3 <- matrix(rep(temp2[1,],mutant),ncol = mutant,byrow = T) + matrix(rnorm(n = mutant*mutant, mean = 0, sd = m),ncol = mutant)
    temp4 <- -t(temp3)+1
    temp3 <- temp3*upper.tri(temp3, diag = FALSE)
    temp4 <- temp4*lower.tri(temp4, diag = FALSE)
    temp5 <- temp3 + temp4
    diag(temp5) <- intra.sp.int
    Interaction <- rbind(cbind(Interaction,temp2),cbind(temp1,temp5))
  }
  return(Interaction)
}

#### parameters
delta <- 0.5            # trade-off parameter (0 < delta < 1)
mu <- 0.001             # mutation rate
m <- 0.02               # mutation width
Ns <- 10^5              # number of possible individuals (carrying capacity)
death <- 0              # death rate = 1/lambda
intra.sp.int <- 0.5     # intra specific interaction
N.threshold <- 1        # threshold for deleting rare strains 
evol.timestep <- 1      # evolutionary time step (mutations happen every evol.timestep)
ecol.timestep <- 0.1    # ecological time step (descrete time of solving population Lotka-Volterra dynamics) To speed up the simulation it is chosen to be high
T.write <- 10000        # every T.write generations the results are written in a file
first.step <- 0
last.step <- 1000000

#### initial conditions
Interaction <- intra.sp.int    # interaction matrix (initial condition with only 1 strain)
Nst <- 1              # number of strains (initial condition with only 1 strain)
pop <- Ns             # vector of population of strains (initial condition with a system full of individuals)

# number of evolutionary steps without writing 
num.internal.step <- T.write/evol.timestep
# trade-off exponent
s <- -log(1-delta,base=2)
# time serie for solving differential equation
times <- seq(ecol.timestep, evol.timestep, by = ecol.timestep)

time0 <- Sys.time()
print(delta)
for (main.step in as.integer((first.step/T.write)+1):as.integer(last.step/T.write)) {
  cat(main.step,",")
  for (step in 1:num.internal.step) {
    ########################
    # calculating the reproduction rate
    ########################
    if (Nst > 1) {
      if (s == 0) {
        rep=rep(1,Nst) # reproduction rate is equal to 1 when there is no trade-off
      } else {
        Comp.abil <- (apply(Interaction,MARGIN = 2,sum)-intra.sp.int)/Nst
        rep <- (1-Comp.abil^(1/s))^s # equation 3 from the appendix
      }
    } else {
      rep=1 # for only one strain reproduction rate is equal to 1 regardless of trade-off
    }
    
    ########################
    # Integration
    ########################
    parameters <- c(Ns,Nst,Interaction,rep,death)
    if (Nst>1) {
      out <- deSolve::ode(y = pop, times = times, func = GLVn, parms = parameters)
    } else {
      out <- deSolve::ode(y = pop, times = times, func = GLV1, parms = parameters)
    }
    
    pop <- out[evol.timestep/ecol.timestep,2:(1+Nst)]
    
    ########################
    # Solving problem of population boom
    ########################
    pop <- LimitPopulation(pop,Ns)
    
    ########################
    # deleting rare species
    ########################
    res <- DeleteRare(pop,Interaction,rep,Nst,N.threshold)
    pop <- res[[1]]
    Interaction <- res[[2]]
    rep <- res[[3]]
    Nst <- res[[4]]
    
    ########################
    # Adding mutants
    ########################
    if (Nst > 1) {
      mean.pop <- colMeans(out[,2:(1+Nst)])
    } else {
      mean.pop <- mean(out[,2:(1+Nst)])
    }
    # number of mutants for each strains is calculated approximately = (part of the population dynamics (Eq. 6 in the apendix) that increases the population of strains) * mutation rate * time interval 
    mutant <- round((mu*(rep*mean.pop*(max(min(Ns,Ns-sum(mean.pop)),0) + as.vector(Interaction%*%mean.pop))/Ns)))*evol.timestep
    num.mutant <- sum(mutant)
    pop <- pop - mutant
    Interaction <- Producemutants(pop,mutant,Interaction,rep,Nst,m,num.mutant,intra.sp.int)
    colnames(Interaction) <- NULL
    pop <- c(pop,rep(1,num.mutant))
    Nst <- Nst + num.mutant
  }
  time1 <- Sys.time()
  cat(time1-time0," | ")
  time0 <- time1
  current.step <- T.write*main.step
  save(list=c("current.step","last.step","pop","Nst","rep","Interaction","mu","m","delta","s","Ns","death","evol.timestep","ecol.timestep","N.threshold","T.write")
       ,file = paste("./",delta,"/Results_",current.step,"_",mu,"_",m,"_",delta,"_",Ns,"_",death,"_",evol.timestep,"_",ecol.timestep,"_",N.threshold,"_",T.write,"_",last.step,".Rdata",sep = ""))
}