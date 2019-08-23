##############################################################
###### Cost-effectiveness of expanding the capacity of #######
########### opioid agonist treatment in Ukraine: #############
################# Dynamic modeling analysis ##################
###################### Addiction 2019 ########################
### Authors: Olga Morozova, Forrest W. Crawford, Ted Cohen,###
########### A David Paltiel, Frederick L. Altice #############
##############################################################
##################### Replication code #######################
##############################################################


# load packages
require(dplyr)
require(deSolve)


### rename city-specific variables to unify variable names across locations ###
## inputs:
# joint.same: sample from joint distribution of model parameters that are the same across locations
# joint.city: sample from joint distribution of model parameters that are different in different locations
## output:
# joint: sample from joint distribution of model parameters (same and city-specific) with unified variable names

smpl.joint.city <- function(joint.same, joint.city){
rename <- c('S0', 'E0', 'On0', 'Of0', 'A0', 'Q0', 'v0', 'w0',
            'lam', 'lam0', 'alp.f', 'alp0.f', 'irr.alp', 'rho.n', 'rho.f', 'dlt', 'gam', 
            'inj.On', 'inj.Q', 'inj.OAT')
colnames(joint.city) <- rename
joint <- as.data.frame(cbind(joint.same, joint.city))
joint$inj.Of <- with(joint, inj.On*inj.ratio.Of.On)
joint$inj.ratio.Of.On <- NULL
return(joint)
}








### compute summary statistics (marginal mean, 95% CI) for the model parameters ###
### inputs:
# smpl.joint: output of function 'smpl.joint.city': sample from joint distribution of model parameters with unified variable names
# city: 1=kyiv, 2=mykolaiv, 3=lviv
# rates: 
## 1 = base case peer effects assumption
## 2 = scenario when OAT patients contribute to opioid use rates
## 0 = no peer effects assumed besides capacity-dependence of OAT drop-out rate
### output:
# table reporting marginal means and 95%CI of model parameters
param.sumry <- function(smpl.joint, city, rates){
   if (city==1){Bs0 <- 829
               Bp0 <- 0}
   if (city==2){Bs0 <- 478
               Bp0 <- 30}
   if (city==3){Bs0 <- 142
               Bp0 <- 0}

smpl.joint$N0 <- with(smpl.joint, S0+E0+On0+Of0+A0+Q0+Bs0+Bp0)
smpl.joint$mu.s <- with(smpl.joint, mu0.s + mu1.s*(1-exp(-k.mu*((Bp0+Bs0)/(N0-S0-E0)))))

smpl.joint$cost.busy.s <- 306 * smpl.joint$cost.adjust.busy.s
smpl.joint$cost.busy.p <- 392 * smpl.joint$cost.adjust.busy.p
smpl.joint$cost.idle.s <- 66 * smpl.joint$cost.adjust.idle.s
smpl.joint$cost.idle.p <- 80 * smpl.joint$cost.adjust.idle.p

smpl.joint$cost.adjust.busy.s <- NULL
smpl.joint$cost.adjust.busy.p <- NULL
smpl.joint$cost.adjust.idle.s <- NULL
smpl.joint$cost.adjust.idle.p <- NULL

if (rates==1){
smpl.joint <- smpl.joint[c("lam", "lam0", "lam1", "alp.f", "alp0.f", "alp1.f", "irr.alp", "k.al",
                           "dlt", "gam", "rho.n", "rho0.n", "rho.f", "rho0.f", "rho1", 
                           "mu.s", "mu0.s", "mu1.s", "irr.drop", "k.mu", "u", "p", 
                           "m.On", "m.B", "m.A", 
                           "S0", "E0", "On0", "Of0", "Q0", "A0", "N0", "v0", "w0",
                           "hrqol.S", "hrqol.E", "hrqol.A", "hrqol.OnOfQ", "dif.hrqol.OAT.DU", 
                           "inj.On", "inj.Of", "inj.Q", "inj.OAT", 
                           "cost.busy.s", "cost.busy.p", "cost.idle.s", "cost.idle.p", "r")]

nms <- c("lam (t=0)", "lam0", "lam1", "alpha.f (t=0)", "alpha0.f", "alpha1.f", "irr.alpha", "k.al",
                           "delta", "gamma", "rho.n (t=0)", "rho0.n", "rho.f (t=0)", "rho0.f", "rho1", 
                           "mu.s (t=0)", "mu0.s", "mu1.s", "irr.mu", "k.mu", "u", "p", 
                           "m.DU", "m.OAT", "m.A", 
                           "S0", "E0", "On0", "Of0", "Q0", "A0", "N0", "v0", "w0",
                           "hrqol.S", "hrqol.E", "hrqol.A", "hrqol.DU", "dif.hrqol.OAT.DU", 
                           "inj.On", "inj.Of", "inj.Q", "inj.OAT", 
                           "cost.busy.s", "cost.busy.p", "cost.idle.s", "cost.idle.p", "r")
}

if (rates==2){
smpl.joint <- smpl.joint[c("lam", "lam0", "lam1", "z.B", "alp.f", "alp0.f", "alp1.f", "irr.alp", "k.al",
                           "dlt", "gam", "rho.n", "rho0.n", "rho.f", "rho0.f", "rho1", 
                           "mu.s", "mu0.s", "mu1.s", "irr.drop", "k.mu", "u", "p", 
                           "m.On", "m.B", "m.A", 
                           "S0", "E0", "On0", "Of0", "Q0", "A0", "N0", "v0", "w0",
                           "hrqol.S", "hrqol.E", "hrqol.A", "hrqol.OnOfQ", "dif.hrqol.OAT.DU", 
                           "inj.On", "inj.Of", "inj.Q", "inj.OAT", 
                           "cost.busy.s", "cost.busy.p", "cost.idle.s", "cost.idle.p", "r")]

nms <- c("lam (t=0)", "lam0", "lam1", "z.B", "alpha.f (t=0)", "alpha0.f", "alpha1.f", "irr.alpha", "k.al",
                           "delta", "gamma", "rho.n (t=0)", "rho0.n", "rho.f (t=0)", "rho0.f", "rho1", 
                           "mu.s (t=0)", "mu0.s", "mu1.s", "irr.mu", "k.mu", "u", "p", 
                           "m.DU", "m.OAT", "m.A", 
                           "S0", "E0", "On0", "Of0", "Q0", "A0", "N0", "v0", "w0",
                           "hrqol.S", "hrqol.E", "hrqol.A", "hrqol.DU", "dif.hrqol.OAT.DU", 
                           "inj.On", "inj.Of", "inj.Q", "inj.OAT", 
                           "cost.busy.s", "cost.busy.p", "cost.idle.s", "cost.idle.p", "r")
}

if (rates==0){
smpl.joint <- smpl.joint[c("lam", "alp.f", "irr.alp", 
                           "dlt", "gam", "rho.n", "rho.f",  
                           "mu.s", "mu0.s", "mu1.s", "irr.drop", "k.mu", "u", "p", 
                           "m.On", "m.B", "m.A", 
                           "S0", "E0", "On0", "Of0", "Q0", "A0", "N0", "v0", "w0",
                           "hrqol.S", "hrqol.E", "hrqol.A", "hrqol.OnOfQ", "dif.hrqol.OAT.DU", 
                           "inj.On", "inj.Of", "inj.Q", "inj.OAT", 
                           "cost.busy.s", "cost.busy.p", "cost.idle.s", "cost.idle.p", "r")]

nms <- c("lam", "alpha.f", "irr.alpha", 
                           "delta", "gamma", "rho.n", "rho.f",  
                           "mu.s (t=0)", "mu0.s", "mu1.s", "irr.mu", "k.mu", "u", "p", 
                           "m.DU", "m.OAT", "m.A", 
                           "S0", "E0", "On0", "Of0", "Q0", "A0", "N0", "v0", "w0",
                           "hrqol.S", "hrqol.E", "hrqol.A", "hrqol.DU", "dif.hrqol.OAT.DU", 
                           "inj.On", "inj.Of", "inj.Q", "inj.OAT", 
                           "cost.busy.s", "cost.busy.p", "cost.idle.s", "cost.idle.p", "r")
}

par.means <- par.CI.low <- par.CI.high <- c(ncol(smpl.joint))
for (i in 1:ncol(smpl.joint)){
   par.means[i] <- round(mean(smpl.joint[,i]), digits=3)
   par.CI.low[i] <- round(quantile(smpl.joint[,i], probs=c(0.025)), digits=3)
   par.CI.high[i] <- round(quantile(smpl.joint[,i], probs=c(0.975)), digits=3)
}
res <- as.data.frame(cbind(nms, par.means, par.CI.low, par.CI.high))
colnames(res) <- c('var.name', 'mean', '95%CI.low', '95%CI.high')
return(res)
}

















### prepare inputs for ODE system integration at the parameter values set to the marginal means of their joint distribution

### inputs:
# joint: output of function 'smpl.joint.city': sample from joint distribution of model parameters with unified variable names
# city: location code (1=kyiv, 2=mykolaiv, 3=lviv)
# Cs and Cp: OAT capacity at specialty / primary care
# ppy: periods per year
# rates: 
## 1 = base case peer effects assumption
## 2 = scenario when OAT patients contribute to opioid use rates
## 0 = no peer effects assumed besides capacity-dependence of OAT drop-out rate

### output: a list with the following elements:
# state0: vector on initial conditions for the model compartments
# params: vector of parameter values
# coef.costs: vector of costs for model compartments 
# coef.qaly: vector of health utility weights for model compartments 
# coef.inj: vector of injection frequency for model compartments  
# r: rate for time-discounting of costs and QALYs

ode.mod.inputs.means <- function(joint, city, Cs, Cp, ppy, rates){
## constants ##
if (city==1){Bs0 <- 829
         Bp0 <- 0
         nu <- 5353/ppy
         xi <- 0.0236/ppy
         m.S <- 0.0021/ppy
         m.E <- 0.0268/ppy}
if (city==2){Bs0 <- 478
         Bp0 <- 30
         nu <- 667/ppy
         xi <- 0.0281/ppy
         m.S <- 0.0029/ppy
         m.E <- 0.0379/ppy}
if (city==3){Bs0 <- 142
         Bp0 <- 0
         nu <- 1053/ppy
         xi <- 0.027/ppy
         m.S <- 0.0025/ppy
         m.E <- 0.0322/ppy}

## marginal means ##
# initial conditions
S0 <- mean(joint$S0) 
E0 <- mean(joint$E0) 
On0 <- mean(joint$On0)
Of0 <- mean(joint$Of0)
A0 <- mean(joint$A0) 
Q0 <- mean(joint$Q0) 
v0 <- mean(joint$v0) 
w0 <- mean(joint$w0) 
state0 <- c(S0,E0,On0,Of0,A0,Q0,Bs0,Bp0,v0,w0)

# model parameters
if (rates==1 | rates==2){
lam0 <- mean(joint$lam0)/ppy 
lam1 <- mean(joint$lam1)/ppy
al0.f <- mean(joint$alp0.f)/ppy 
al1.f <- mean(joint$alp1.f)/ppy 
al0.n <- mean(joint$alp0.f * joint$irr.alp)/ppy
al1.n <- mean(joint$alp1.f * joint$irr.alp)/ppy
rho0.n <- mean(joint$rho0.n)/ppy 
rho1.n <- rho1.f <- mean(joint$rho1)/ppy 
rho0.f <- mean(joint$rho0.f)/ppy }

if (rates==0){
lam0 <- mean(joint$lam)/ppy 
lam1 <- 0
al0.f <- mean(joint$alp.f)/ppy 
al1.f <- al1.n <- 0
al0.n <- mean(joint$alp0.f * joint$irr.alp)/ppy
rho0.n <- mean(joint$rho.n)/ppy 
rho1.n <- rho1.f <- 0 
rho0.f <- mean(joint$rho.f)/ppy}

if (rates==1 | rates==0){
z.B <- 0}

if (rates==2){
z.B <- mean(joint$z.B)}

mu0.s <- mean(joint$mu0.s)/ppy 
mu1.s <- mean(joint$mu1.s)/ppy 
mu0.p <- mean(joint$mu0.s * joint$irr.drop)/ppy 
mu1.p <- mean(joint$mu1.s * joint$irr.drop)/ppy 

dlt <- mean(joint$dlt)/ppy 
gam.n <- gam.f <- dlt.a <- mean(joint$gam)/ppy 

u <- mean(joint$u)
p <- mean(joint$p) 

m.On <- mean(joint$m.On)/ppy
m.Of <- m.On
m.A <- mean(joint$m.A)/ppy 
m.Q <- m.On
m.Bs <-  mean(joint$m.B)/ppy  
m.Bp <- m.Bs

k.al <- mean(joint$k.al)
k.mu <- mean(joint$k.mu)


# Coefficients for costs and outcomes for all compartments
# COSTS Order: S  E On  Of  A  Q  Bs  Bp  Cs-Bs  Cp-Bp #
cost.busy.s <- 306 * mean(joint$cost.adjust.busy.s)
cost.busy.p <- 392 * mean(joint$cost.adjust.busy.p)

cost.idle.s <- 66 * mean(joint$cost.adjust.idle.s)
cost.idle.p <- 80 * mean(joint$cost.adjust.idle.p)

coef.costs <- c(0, 0, 0, 0, 0, 0, cost.busy.s, cost.busy.p, cost.idle.s, cost.idle.p)/ppy
r <- mean(joint$r) 

# HRQoL utility weights #
qaly.S <- mean(joint$hrqol.S) 
qaly.E <- mean(joint$hrqol.E) 
qaly.On <- mean(joint$hrqol.OnOfQ) 
qaly.Of <- qaly.On
qaly.A <- mean(joint$hrqol.A) 
qaly.Q <- qaly.On
qaly.Bs <- qaly.On + mean(joint$dif.hrqol.OAT.DU) 
qaly.Bp <- qaly.Bs
# Order: S  E On  Of  A  Q  Bs  Bp
coef.qaly <- c(qaly.S, qaly.E, qaly.On, qaly.Of, qaly.A, qaly.Q, qaly.Bs, qaly.Bp)


# Injections #
inj.S <- inj.E <- inj.A <- 0
inj.On <- mean(joint$inj.On) 
inj.Of <- mean(joint$inj.Of) 
inj.Q <- mean(joint$inj.Q) 
inj.Bs <- mean(joint$inj.OAT) 
inj.Bp <- inj.Bs

coef.inj <- c(inj.S, inj.E, inj.On, inj.Of, inj.A, inj.Q, inj.Bs, inj.Bp)/ppy 

# collect model parameter values in a vector
params <- c(nu, xi, 
            lam0, lam1, 
            al0.n, al1.n, al0.f, al1.f, 
            mu0.s, mu1.s, mu0.p, mu1.p, 
            u, 
            dlt, 
            dlt.a, gam.n, gam.f, 
            rho0.n, rho1.n, rho0.f, rho1.f, 
            p, 
            m.S, m.E, m.On, m.Of, m.A, m.Q, m.Bs, m.Bp, 
            k.al, k.mu, z.B,
            Cs, Cp)

return(list(state0, params, coef.costs, coef.qaly, coef.inj, r))
}














### prepare inputs for ODE system integration at the parameter values set to a given draw from the sample of their joint distribution

### inputs:
# joint: output of function 'smpl.joint.city': sample from joint distribution of model parameters with unified variable names
# j: number of the line (draw) from joint distribution of model parameters
# city: location code (1=kyiv, 2=mykolaiv, 3=lviv)
# Cs and Cp: OAT capacity at specialty / primary care
# ppy: periods per year
# rates: 
## 1 = base case peer effects assumption
## 2 = scenario when OAT patients contribute to opioid use rates
## 0 = no peer effects assumed besides capacity-dependence of OAT drop-out rate

### output: a list with the following elements:
# state0: vector on initial conditions for the model compartments
# params: vector of parameter values
# coef.costs: vector of costs for model compartments 
# coef.qaly: vector of health utility weights for model compartments 
# coef.inj: vector of injection frequency for model compartments  
# r: rate for time-discounting of costs and QALYs

ode.mod.inputs.draw <- function(joint, j, city, Cs, Cp, ppy, rates){
## constants ##
if (city==1){Bs0 <- 829
         Bp0 <- 0
         nu <- 5353/ppy
         xi <- 0.0236/ppy
         m.S <- 0.0021/ppy
         m.E <- 0.0268/ppy}
if (city==2){Bs0 <- 478
         Bp0 <- 30
         nu <- 667/ppy
         xi <- 0.0281/ppy
         m.S <- 0.0029/ppy
         m.E <- 0.0379/ppy}
if (city==3){Bs0 <- 142
         Bp0 <- 0
         nu <- 1053/ppy
         xi <- 0.027/ppy
         m.S <- 0.0025/ppy
         m.E <- 0.0322/ppy}
 
# initial conditions
S0 <- joint$S0[j]
E0 <- joint$E0[j] 
On0 <- joint$On0[j]
Of0 <- joint$Of0[j]
A0 <- joint$A0[j] 
Q0 <- joint$Q0[j] 
v0 <- joint$v0[j] 
w0 <- joint$w0[j] 
state0 <- c(S0,E0,On0,Of0,A0,Q0,Bs0,Bp0,v0,w0)

# model parameters
if (rates==1 | rates==2){
lam0 <- joint$lam0[j]/ppy 
lam1 <- joint$lam1[j]/ppy
al0.f <- joint$alp0.f[j]/ppy 
al1.f <- joint$alp1.f[j]/ppy 
rho0.n <- joint$rho0.n[j]/ppy 
rho1.n <- rho1.f <- joint$rho1[j]/ppy 
rho0.f <- joint$rho0.f[j]/ppy }

if (rates==0){
lam0 <- joint$lam[j]/ppy 
lam1 <- 0
al0.f <- joint$alp.f[j]/ppy 
al1.f <- 0
rho0.n <- joint$rho.n[j]/ppy 
rho1.n <- rho1.f <- 0 
rho0.f <- joint$rho.f[j]/ppy}

if (rates==1 | rates==0){
z.B <- 0}

if (rates==2){
z.B <- joint$z.B[j]}

irr.alp <- joint$irr.alp[j] 
al0.n <- al0.f * irr.alp 
al1.n <- al1.f * irr.alp 

mu0.s <- joint$mu0.s[j]/ppy 
mu1.s <- joint$mu1.s[j]/ppy 
irr.mu <- joint$irr.drop[j] 
mu0.p <- mu0.s * irr.mu 
mu1.p <- mu1.s * irr.mu 

dlt <- joint$dlt[j]/ppy 
dlt.a <- gam.n <- gam.f <- joint$gam[j]/ppy 

u <- joint$u[j]
p <- joint$p[j] 

m.On <- joint$m.On[j]/ppy
m.Of <- m.On
m.A <- joint$m.A[j]/ppy 
m.Q <- m.On
m.Bs <-  joint$m.B[j]/ppy  
m.Bp <- m.Bs

k.al <- joint$k.al[j]
k.mu <- joint$k.mu[j]


# Coefficients for costs and outcomes for all compartments
# COSTS Order: S  E On  Of  A  Q  Bs  Bp  Cs-Bs  Cp-Bp #
cost.busy.s <- 306 * joint$cost.adjust.busy.s[j]
cost.busy.p <- 392 * joint$cost.adjust.busy.p[j]

cost.idle.s <- 66 * joint$cost.adjust.idle.s[j]
cost.idle.p <- 80 * joint$cost.adjust.idle.p[j]

coef.costs <- c(0, 0, 0, 0, 0, 0, cost.busy.s, cost.busy.p, cost.idle.s, cost.idle.p)/ppy
r <- joint$r[j] 

# HRQoL utility weights #
qaly.S <- joint$hrqol.S[j] 
qaly.E <- joint$hrqol.E[j] 
qaly.On <- joint$hrqol.OnOfQ[j] 
qaly.Of <- qaly.On
qaly.A <- joint$hrqol.A[j] 
qaly.Q <- qaly.On
qaly.Bs <- qaly.On + joint$dif.hrqol.OAT.DU[j] 
qaly.Bp <- qaly.Bs
# Order: S  E On  Of  A  Q  Bs  Bp
coef.qaly <- c(qaly.S, qaly.E, qaly.On, qaly.Of, qaly.A, qaly.Q, qaly.Bs, qaly.Bp)


# Injections #
inj.S <- inj.E <- inj.A <- 0
inj.On <- joint$inj.On[j] 
inj.Of <- joint$inj.Of[j] 
inj.Q <- joint$inj.Q[j] 
inj.Bs <- joint$inj.OAT[j] 
inj.Bp <- inj.Bs

coef.inj <- c(inj.S, inj.E, inj.On, inj.Of, inj.A, inj.Q, inj.Bs, inj.Bp)/ppy 

# collect model parameter values in a vector
params <- c(nu, xi, 
            lam0, lam1, 
            al0.n, al1.n, al0.f, al1.f, 
            mu0.s, mu1.s, mu0.p, mu1.p, 
            u, 
            dlt, 
            dlt.a, gam.n, gam.f, 
            rho0.n, rho1.n, rho0.f, rho1.f, 
            p, 
            m.S, m.E, m.On, m.Of, m.A, m.Q, m.Bs, m.Bp, 
            k.al, k.mu, z.B,
            Cs, Cp)

return(list(state0, params, coef.costs, coef.qaly, coef.inj, r))
}

























### rule for allocation of patients in the waiting list between types of facilities
## inputs:
# fs.s, fs.p: free slots in specialty and primary care
# tq: length of the waiting list
# p: probability of giving a preference to primary care
## output:
# number of patients allocated to specialty and primary care
allocate.new <- function(fs.s, fs.p, tq, p)
{
new.tx.p <- 0
new.tx.s <- 0
tfs <- fs.p + fs.s

if (tq >= tfs){
   new.tx.p <- fs.p
   new.tx.s <- fs.s}
if (tq < tfs){
            if (p*tq>=fs.p){
               new.tx.p <- fs.p
               new.tx.s <- tq-fs.p}
            if ((1-p)*tq>=fs.s){
               new.tx.s <- fs.s
               new.tx.p <- tq-fs.s}
            if (p*tq<fs.p & (1-p)*tq<fs.s){
               new.tx.p <- p*tq 
               new.tx.s <- (1-p)*tq}}
return (c(new.tx.s, new.tx.p))
}



### compute mlti-component rates 
f.lin.rate <- function(par0, par1, prop){
   return(par0 + par1*prop)
}

f.sat.rate <- function(par0, par1, prop, k){
   a <- 1-exp(-k*prop)
   return(par0 + par1*a)
}














### ODE system model ###  
## inputs:
# t: vector of time points for which output is wanted
# state: vector of initial conditions
# parameters: vector of parameter values
## output:
# list of derivatives in the ODE system
oat.model <- function(t, state, parameters){
# states
S <- state[1]
E <- state[2]
On <- state[3]
Of <- state[4]
A <- state[5]
Q <- state[6]
Bs <- state[7]
Bp <- state[8]
v <- state[9]
w <- state[10]

# parameters
nu <- parameters[1]
xi <- parameters[2]
lam0 <- parameters[3]
lam1 <- parameters[4]
al0.n <- parameters[5]
al1.n <- parameters[6]
al0.f <- parameters[7]
al1.f <- parameters[8]
mu0.s <- parameters[9]
mu1.s <- parameters[10]
mu0.p <- parameters[11]
mu1.p <- parameters[12]
u <- parameters[13]
dlt <- parameters[14]
dlt.a <- parameters[15]
gam.n <- parameters[16]
gam.f <- parameters[17]
rho0.n <- parameters[18]
rho1.n <- parameters[19]
rho0.f <- parameters[20]
rho1.f <- parameters[21]
p <- parameters[22]
m.S <- parameters[23]
m.E <- parameters[24]
m.On <- parameters[25]
m.Of <- parameters[26]
m.A <- parameters[27]
m.Q <- parameters[28]
m.Bs <- parameters[29] 
m.Bp <- parameters[30]
k.al <- parameters[31]
k.mu <- parameters[32]
z.B <- parameters[33]
Cs <- parameters[34]
Cp <- parameters[35]

# proportions and multi-component transition rates
N <- S+E+On+Of+A+Q+Bs+Bp
du.prop <- (On+Of+Q+(Bs+Bp)*z.B)/(N-E)
slot.prop <- (Cs+Cp)/(N-(S+E))
lam <- f.lin.rate(lam0,lam1,du.prop)
al.n <- f.sat.rate(al0.n, al1.n, slot.prop, k.al)
al.f <- f.sat.rate(al0.f, al1.f, slot.prop, k.al)
mu.s <- f.sat.rate(mu0.s, mu1.s, slot.prop, k.mu)
mu.p <- f.sat.rate(mu0.p, mu1.p, slot.prop, k.mu)
rho.n <- f.lin.rate(rho0.n, rho1.n, du.prop)
rho.f <- f.lin.rate(rho0.f, rho1.f, du.prop)
fs.s <- Cs-Bs
fs.p <- Cp-Bp
new.tx <- allocate.new(fs.s, fs.p, tq=Q, p)
new.tx.s <- new.tx[1]
new.tx.p <- new.tx[2]
new.tx <- new.tx.s+new.tx.p

# transition equations
dS <- -S*(lam+xi+m.S) + nu
dE <- -E*m.E + S*xi
dOn <- -On*(al.n+gam.n+m.On) + lam*S + v*dlt*Q + w*rho.n*A
dOf <- -Of*(al.f+gam.f+m.Of) + (1-u)*(mu.s*Bs+mu.p*Bp) + (1-v)*dlt*Q + (1-w)*rho.f*A
dA <- -A*(w*rho.n+(1-w)*rho.f+m.A) + gam.n*On + gam.f*Of + dlt.a*Q + u*(mu.s*Bs+mu.p*Bp)
dQ <- -Q*(dlt+dlt.a+m.Q) - new.tx + al.n*On + al.f*Of
dBs <- -Bs*(mu.s+m.Bs) + new.tx.s
dBp <- -Bp*(mu.p+m.Bp) + new.tx.p

dv <- (-v*(al.n*On + al.f*Of) + al.n*On)/(Q*(1-dlt-dlt.a-m.Q) - new.tx + al.n*On + al.f*Of)
dw <- (-w*(gam.n*On+gam.f*Of+dlt.a*Q+u*(mu.s*Bs+mu.p*Bp))+gam.n*On+v*dlt.a*Q)/(A*(1-w*rho.n-(1-w)*rho.f-m.A)+gam.n*On+gam.f*Of+dlt.a*Q+u*(mu.s*Bs+mu.p*Bp))
if ((v+dv)<0 | (v+dv)>1) {dv <- 0}
if ((w+dw)<0 | (w+dw)>1) {dw <- 0}

return(list(c(dS,dE,dOn,dOf,dA,dQ,dBs,dBp,dv,dw)))
} 















### detailed results of the ODE system integration: temporal dynamics

## inputs:
# state0: vector on initial conditions for the model compartments
# times: vector of time points for which output is wanted
# model: ODE system model
# params: vector of parameter values
# coef.costs: vector of costs for model compartments 
# coef.qaly: vector of health utility weights for model compartments 
# coef.inj: vector of injection frequency for model compartments  
# r: rate for time-discounting of costs and QALYs
# ppy: periods per year

## output: a data frame with rows corresponding to 'times' and colums reporting:
# values of all model compartments at given times
# values of modeled population, POUD, active POUD sizes at given times
# proportion in OAT among POUD and capacity / POUD at given times
# total cost, accumulated QALYs, number of injections and number of drug use initiations per period 

ode.int.detail <- function(state0, times, model, params, coef.costs, coef.qaly, coef.inj, r, ppy){
Cs <- params[34]
Cp <- params[35]
lam0 <- params[3]
lam1 <- params[4]
z.B <- params[33]

out <- rk4(state0, times, oat.model, params, verbose = F)
out <- as.data.frame(out)
colnames(out) <- c('times', 'S', 'E', 'On', 'Of', 'A', 'Q', 'Bs', 'Bp', 'v', 'w')
out$N <- with(out, S+E+On+Of+A+Q+Bs+Bp)
out$POUD <- with(out, N-(S+E))
out$actPOUD <- with(out, POUD-A-Bs-Bp)
out$prop.OAT <- with (out, (Bs+Bp)/POUD)
out$prop.capacity <- with(out, (Cs+Cp)/POUD)

# add 4 columns: cost, qaly, injections and drug use initiations per period #
out$dsct <- (((1+r)^(1/ppy))^out$times)^(-1)
cc <- coef.costs
qc <- coef.qaly
ic <- coef.inj
out$cost <- with(out, dsct*(S*cc[1]+E*cc[2]+On*cc[3]+Of*cc[4]+A*cc[5]+Q*cc[6]+Bs*cc[7]+Bp*cc[8]+(Cs-Bs)*cc[9]+(Cp-Bp)*cc[10]))
out$qaly <- with(out, dsct*(S*qc[1]+E*qc[2]+On*qc[3]+Of*qc[4]+A*qc[5]+Q*qc[6]+Bs*qc[7]+Bp*qc[8])/ppy)
out$inj <- with(out, S*ic[1]+E*ic[2]+On*ic[3]+Of*ic[4]+A*ic[5]+Q*ic[6]+Bs*ic[7]+Bp*ic[8])
out$du.init <- with(out, S*(lam0+lam1*((On+Of+Q+(Bs+Bp)*z.B)/(N-E))))
return(out)
}














### summary results of the ODE system integration

## inputs:
# state0: vector on initial conditions for the model compartments
# times: vector of time points for which output is wanted
# model: ODE system model
# params: vector of parameter values
# coef.costs: vector of costs for model compartments 
# coef.qaly: vector of health utility weights for model compartments 
# coef.inj: vector of injection frequency for model compartments  
# r: rate for time-discounting of costs and QALYs
# ppy: periods per year

## output: a vector of:
# specialty and primary care capacity values that define a strategy
# modeling horizon (years)
# capacity utilization at the end of modeling horizon
# in OAT among POUD at the end of modeling horizon
# total time-discounted cost of a given strategy over a given modeling horizon
# total time-discounted number of accumulated QALYs for a given strategy over a given modeling horizon
# total number of injections for a given strategy over a given modeling horizon
# total number of drug use initiations for a given strategy over a given modeling horizon

ode.int.sum <- function(state0, times, model, params, coef.costs, coef.qaly, coef.inj, r, ppy){
out <- ode.int.detail(state0, times, model, params, coef.costs, coef.qaly, coef.inj, r, ppy)
Cs <- params[34]
Cp <- params[35]
iend <- nrow(out)
years <- nrow(out)/365
util.end <- (out$Bs[iend]+out$Bp[iend])/(Cs+Cp)
inOAT.end <- (out$Bs[iend]+out$Bp[iend])/out$POUD[iend]

smry <- c(format(Cs, scientific = F),
          format(Cp, scientific = F),
          format(years, scientific = F),
          round(util.end, digits=2),
          round(inOAT.end, digits=3),
          round(sum(out$cost), digits=0), 
          round(sum(out$qaly), digits=0), 
          round(sum(out$inj), digits=0), 
          round(sum(out$du.init), digits=0))
names(smry) <- c('CS','CP','modeling.horizon', 'capacity.utilization.end', 'OAT.coverage.end', 
                'total.cost','total.qaly','total.injections','total.drug.use.initiations') 
return(smry)
}














