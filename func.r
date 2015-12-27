
EAKF <- function(distr.tseries, N1, N2, distr.dist, ntrn, nens, nfor, trn.init = 1){
	# Ensemble adjustment Kalman filter
	# Parts of code:
	#	1. -> initialize starting ensemble members
	#	2. -> train Kalman filter on first ntrn incidence data
	#	3. -> forecast using trained Kalman filter
	#
	# ODE subroutine (method for SEIR prediction step):
	#	-> Runge-Kutta 4th order
	# 	-> something else?
	
	# Inputs :
	#	ntrn -> number of training steps
	#	nens -> number of ensembles
	#	nfor -> number of forecast steps (after training period)
	#	tseries -> weekly incidence data (cols are different districts)
	#	N1, N2, distr.dist -> matrices for connectivity matrix
	
	# check viability of inputs:
	if (!all( dim(distr.dist) == dim(N1) & dim(N1) == dim(N2)))
		warning('Wrong N1, N2, distr.dist dimensions.')
	if (dim(distr.dist)[1] != dim(distr.tseries)[1])
		warning('Wrong tseries dimension.')
	if (ntrn > dim(distr.tseries)[2])
		warning('More training data required then available.')
	
	# init variables
	nobs <- dim(distr.tseries)[2]		# number of observations
	ndist <- dim(distr.dist)[2]	# number of districs
	nvar <- 4*ndist + 5		# number of variables
					# for district: S,E,I,R (not recovered)
					#	global: Z, D, tau1, tau2, ro
	obs.mean <- 3
	obs.sd <- 2
	
	
	# Part 1. - INIT
	# In the original paper this was done by sampling from latin hypercube
	# Here it will be done less sophisticated first
	
	# variables are orderes as Ss, Es, Is, Rs, global variables are last
	xprior <- array(0,c(nvar,nens))
	xpost  <- array(0,c(nvar,nens))
	
	# this is to be done for forecasts
	# fcast=array(0,c(3,nens,nfor));	# S, E, I
	
	# generate random intial draw
	# S0:
	xprior[1:ndist,] <- matrix(runif(nens*ndist, 
					 	min = 0.4*N1[,1], max = N1[,1]), ncol = nens)
	# E0
	xprior[1:ndist+ndist,] <- matrix(rtruncnorm(nens*ndist, 
						a=0, b=Inf, mean = obs.mean, sd = 2*obs.sd), ncol = nens)
	# I0
	xprior[1:ndist+2*ndist,] <- matrix(rtruncnorm(nens*ndist, 
						a=0, b=Inf, mean = obs.mean, sd = 2*obs.sd), ncol = nens)
	# R0 (these are not recovered !)
	xprior[1:ndist+3*ndist,] <- matrix(runif(nens*ndist, 
						min = 0.5, max = 3.5), ncol = nens)
	# Z
	xprior[1:ndist+3*ndist+1,] <- matrix(runif(nens, 
						min = 2, max = 14), ncol = nens)
	# D
	xprior[1:ndist+3*ndist+2,] <- matrix(runif(nens,
						min = 5, max = 14), ncol = nens)
	# tau1, tau2
	xprior[1:ndist+3*ndist+3,] <- matrix(runif(nens*2, 
						min = 0.3, max = 0.7), ncol = nens)
	# ro
	xprior[1:ndist+3*ndist+4,] <- matrix(runif(nens, 
						   min = 2, max = 8), ncol = nens)
	
	View(xprior)
	# Part 2. - Train filter
	for (i in 1:ntrn+trn.init-1){
		# propagate SEIR - NSDE method here
		# update filter
		var.obs <- apply(xpost,1,var)	# here goes inflation
		var.prior <- apply(xprior,1,var)
		xprior.mean <- apply(xprior,1,mean)
		z <- distr.tseries[i,]	# new incidence data
		
		den <- var.obs + var.prior
		# eq. (4.4):
		xpost.mean <- (var.obs*xprior.mean + var.prior*z)/den
		# eq. (4.5):
		xpost <- xpost.mean + sqrt(var.obs/den)*(xprior-xprior.mean)
	}
}


rk4 <- function(t1, t2, dt, S, E, I, param, N1, N2, distr.dist){
	# propagate for one ensemble
	
	t.vec <- seq(t1, t2, by=dt)
	for (t in t.vec){
		
	}
}

camelCase <- function(x) {
	# "I like pizza" to "I Like Pizza"
	s <- tolower(x)
	s <- strsplit(s, " ")[[1]]
	paste(toupper(substring(s, 1,1)), substring(s, 2),
	      sep="", collapse="")
}