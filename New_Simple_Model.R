# New model for phytotron experiment

# Parameters - not all used, params df copied from CORPSE script
Tref=293.15
Rugas=8.314472
initc <- c(2, 3, 1)
theta <- 0.33
thetasat <- 0.46
chem_types<-c('Fast','Slow','Necro')
params<-data.frame(
  "Vmaxref_Fast" = 9.0,
  "Vmaxref_Slow" = 0.25,
  "Vmaxref_Necro"= 4.5,
  "Ea_Fast" = 5e3,
  "Ea_Slow" = 30e3,
  "Ea_Necro"= 3e3,
  "kC_Fast" = 0.01,
  "kC_Slow" = 0.01,
  "kC_Necro"= 0.01,
  "gas_diffusion_exp" = 0.6,
  "minMicrobeC" = 1e-3,
  "Tmic"= 0.25,
  "et" = 0.6,
  "eup_Fast" = 0.6,
  "eup_Slow"= 0.1,
  "eup_Necro" = 0.6,
  "tProtected" = 100.0,
  "frac_N_turnover_min" = 0.2,
  "protection_rate_Fast" = 6, # was 0.7
  "protection_rate_Slow" = 3, # was 0.001
  "protection_rate_Necro" = 4.0,
  "nup_Fast" = 0.3,
  "nup_Slow" = 0.3,
  "nup_Necro" = 0.3,
  "CN_Microbe" = 7,
  "max_immobilization_rate" = 3.65,
  "substrate_diffusion_exp" = 1.5,
  "new_resp_units" = TRUE,
  "iN_loss_rate" = 5.0,
  "frac_turnover_slow" = 0.2
)

# Function to calculate Vmax - copied from CORPSE script

##Function to calculate Vmax of microbial decomposition 
##Vmax function, normalized to Tref=293.15 (T is in Kelvin)
Vmax<-function (T,params,Tref=293.15,Rugas=8.314472) {
  # Fast<-params$Vmaxref_Fast*exp(-params$Ea_Fast*(1.0/(Rugas*T)-1.0/(Rugas*Tref)))
  # Slow<-params$Vmaxref_Slow*exp(-params$Ea_Slow*(1.0/(Rugas*T)-1.0/(Rugas*Tref)))
  # Necro<-params$Vmaxref_Necro*exp(-params$Ea_Necro*(1.0/(Rugas*T)-1.0/(Rugas*Tref)))
  # Vmax<-data.frame(Fast,Slow,Necro)
  Vmax<-data.frame(params[paste('Vmaxref_',chem_types,sep='')]*exp(-params[paste('Ea_',chem_types,sep='')]*(1.0/(Rugas*T)-1.0/(Rugas*Tref))))
  names(Vmax)<-chem_types
  return(Vmax)
}


# Testing Vmax ~ temperature
temp <- matrix(0, ncol = 4, nrow = 21)
colnames(temp) <- c("Temperature", "Vmaxfast", "Vmaxslow", "Vmaxnecro")
temp[,1] <- c(290:310)

for (i in 1:length(temp[,1])) {
  temp[i,2] <- Vmax(temp[i,1], params, Tref=293.15, Rugas=8.314472)[[1]]
  temp[i,3] <- Vmax(temp[i,1], params, Tref=293.15, Rugas=8.314472)[[2]]
  temp[i,4] <- Vmax(temp[i,1], params, Tref=293.15, Rugas=8.314472)[[3]]
}
temp <- as.data.frame((temp))
plot(Vmaxfast ~ Temperature, data = temp) 
# Question 1: Should Vmax be linear with temperature?


# Testing decomposition rate (k) using CORPSE equation from tutorial, with 
# Vmax, kC, theta, and thetasat as constants
temp$decompfast <- temp$Vmaxfast * (
                    ((theta / thetasat)^3) *
                    ((1 - (theta / thetasat))^2.5) *
                    (initc[1] * (initc[1] / (initc[1] + params$kC_Fast)))
                   )

plot(decompfast ~ Temperature, data = temp) 
# showing decomposition rate as a linear transformation of Vmax....
# Question 2: Should decomp simply be a linear combination of Vmax and theta?

# Modeling daily decomposition as a fraction of annual k - annual k / 365
dayconvert <- 1/364
fastCpool <- as.vector(initc[1]) 
for (i in 1:364) {
  c <- fastCpool[i] - (temp$decompfast[1] * dayconvert)
  fastCpool <- c(fastCpool, c)
}
plot(fastCpool ~ c(1:365)) # plot daily pool size
result <- paste("Annual loss is",
                signif(temp$decompfast[1], 3),
                "kg C and cumulative daily loss is", 
                signif(initc[1] - tail(fastCpool, n = 1), 3),
                "kf C. With",
                signif(((tail(fastCpool, n = 1) / fastCpool[1]) * 100), 3), 
                "% of initial C remaining"
                )
result
mrfastCpool <- fastCpool / fastCpool[1] 
plot(mrfastCpool ~ c(1:365)) # plot % mass remaining
# C pool decays linearly, both as mass loss and % mass remaining
# Question 3: If decomposition rate is a function of pool size and the step is daily, 
# should we be calculating decomp rate daily, instead of dividing annual decomp by 365?


# Modeling daily decomposition rate, based on daily C pool size





