
EAKF <- function(distr.tseries, N, N1, N2, distr.dist, ntrn, nens, nfor, trn.init = 24){
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
	#	trn.init -> initial week (24th obs. is where all districts had ebola)
	
	# check viability of inputs:
	if (!all( dim(distr.dist) == dim(N1) & dim(N1) == dim(N2)))
		warning('Wrong N1, N2, distr.dist dimensions.')
	if (dim(distr.dist)[1] != dim(distr.tseries)[1])
		warning('Wrong tseries dimension.')
	if (ntrn > dim(distr.tseries)[2])
		warning('More training data required then observations available.')
	
	# init variables
	nobs <- dim(distr.tseries)[2]	# number of observations
	ndist <- dim(distr.tseries)[1]	# number of districs
	nvar <- 4*ndist + 5		# number of variables
					# for district: S,E,I,R (not recovered)
					#	global: Z, D, tau1, tau2, ro
	
	trn.range <- 1:ntrn+trn.init-1		# training range for tseries
	for.range <- 1:nfor+ntrn+trn.init-1	# forecast range 
	obs.mean <- apply(distr.tseries[,trn.range],1,mean)	# mean of whole trn.range
	obs.sd   <- apply(distr.tseries[,trn.range],1,sd)	# std. dev. of whole trn.range (for init. draw)
	
	#
	# Part 1. - INIT
	#
	
	# In the original paper this was done by sampling from latin hypercube
	# Here it will be done less sophisticated first
	
	# variables are orderes as Ss, Es, Is, Rs, global variables are last
	xprior <- array(0,c(nvar,nens))
	xpost  <- array(0,c(nvar,nens))
	xobs  <- array(0,c(nvar,nens))
	
	# this is to be done for forecasts
	# fcast=array(0,c(3,nens,nfor));	# S, E, I
	
	# generate random intial draw
	# S0:
	xprior[1:ndist,] <- matrix(runif(nens*ndist, 
					 	min = 0.4*N, max = N), ncol = nens)
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
	xprior[4*ndist+1,] <- matrix(runif(nens, 
						min = 2, max = 14), ncol = nens)
	# D
	xprior[4*ndist+2,] <- matrix(runif(nens,
						min = 5, max = 14), ncol = nens)
	# tau1, tau2
	xprior[4*ndist+3:4,] <- matrix(runif(nens*2, 
						min = 0.3, max = 0.7), ncol = nens)
	# ro
	xprior[4*ndist+5,] <- matrix(runif(nens, 
						   min = 2, max = 8), ncol = nens)

	# Create initial connectivity matrices
	Cns <- array(0,c(ndist,ndist,nens))
	for (j in 1:nens){
		Cns[,,j] <- create.Cm(N1 = N1, N2 = N2, distr.dist = distr.dist, 
					   xprior[4*ndist+3, j],
					   xprior[4*ndist+4, j],
					   xprior[4*ndist+5, j])	# last 3 are tau1, tau2, ro
	}
	
	#
	# Part 2. - Train filter
	#
	
	print('Training:')
	for (i in trn.range){
		print(i)
		# propagate SEIR using NSDE method (for each ensemble)
		for (j in 1:nens){
			xobs[,j] <- rk4(t1 = i,t2 = i+1, dt = 1/7, N = N, xprior.n = xprior[,j], Cn = Cns[,,j])
		}
		# update filter
		var.obs <- apply(xobs,1,var)*1.03	# here will go inflation OEV
		var.prior <- apply(xprior,1,var)
		xprior.mean <- apply(xprior,1,mean)
		z <- c(numeric(2*ndist),distr.tseries[,i],numeric(nvar-3*ndist))	# new incidence data
											# Note: we observe only new I
		
		den <- var.obs + var.prior
		# eq. (4.4):
		xpost.mean <- (var.obs*xprior.mean + var.prior*z)/den		# vector (nvar,1)
		# eq. (4.5):
		xpost <- (xpost.mean + sqrt(var.obs/den)) %*% t(rep(1,nens))*
				(xprior-xprior.mean %*% t(rep(1,nens)))		# matrix (nvar, nens)

		xprior <- xpost

		# update Cns array
		for (j in 1:nens){
			Cns[,,j] <- create.Cm(N1 = N1, N2 = N2, distr.dist = distr.dist, 
						xprior[4*ndist+3, j],
						xprior[4*ndist+4, j],
						xprior[4*ndist+5, j])	# last 3 are tau1, tau2, ro
		}
	}
	#
	# Part 3. - Forecast
	#
# 	
# 	xfor <- array(0,c(nvar,nens,nfor+1))
# 	xfor[,,1] <- xprior
# 	n = 2
# 	print('Forecasting:')
# 	for (i in for.range){
# 		print(i)
# 		for (j in 1:nens){
# 			xfor[,j,n] <- rk4(t1 = i, t2 = i+1, dt = 1/7, N = N, 
# 					  xprior.n = xfor[,j,n-1], Cn = Cns[,,j])
# 		}
# 		n = n + 1
# 	}
# 	return(xfor)
}

rk4 <- function(t1, t2, dt, N, xprior.n, Cn){
	nvar <- length(xprior.n)
	ndist <- dim(Cn)[1]
	SEI <- xprior.n[1:(3*ndist)]
	Z <- xprior.n[nvar-4]
	D <- xprior.n[nvar-3]
	beta <- xprior.n[1:ndist+3*ndist]/D
	N.hat <- Cn %*% N
	print('###### SEI:')
	print(SEI)
	# propagate for one ensemble
	xpost.n <- xprior.n
	t.vec <- seq(t1, t2, by=dt)
	for (t in t.vec){
		k1s <- SEIRfunc(SEI, N.hat, beta, alpha, Z, D, Cn)
		k2s <- SEIRfunc(SEI + 0.5*dt*k1s, N.hat, beta, alpha, Z, D, Cn)
		k3s <- SEIRfunc(SEI + 0.5*dt*k2s, N.hat, beta, alpha, Z, D, Cn)
		k4s <- SEIRfunc(SEI +   1*dt*k3s, N.hat, beta, alpha, Z, D, Cn)
		SEI <- SEI + dt*(k1s + 2*k2s + 2*k3s + k4s)/6
	}
	
	xpost.n[1:(3*ndist)] <- SEI
	return(xpost.n)
}

SEIRfunc <- function(SEI, N.hat, beta, alpha,  Z, D, Cn){
	ndist <- length(SEI)/3
	S <- SEI[1:ndist] 
	E <- SEI[1:ndist+ndist] 
	I <- SEI[1:ndist+2*ndist]

	h1 <- S*sum(t(Cn) %*% beta*I/N.hat)		# I.hat := I
	h2 <- E/Z
	SEI[1:ndist] <- -h1 - alpha			# S
	SEI[1:ndist+ndist] <- h1 - h2 + alpha	# E
	SEI[1:ndist+2*ndist] <- h2 - I/D	# I
	SEI <- check.bounds(SEI.new = SEI, ndist = ndist, N.hat =  N.hat)
	return(SEI)
}

camelCase <- function(x) {
	# "I like pizza" to "ILikePizza"
	s <- tolower(x)
	s <- strsplit(s, " ")[[1]]
	paste(toupper(substring(s, 1,1)), substring(s, 2),
	      sep="", collapse="")
}

check.bounds <- function(SEI.new, ndist, N.hat){
	# posterior corrections for xnew (xnew is matrix)
	N.hat3 <-
	# SEI => 0, otherwise put 0 instead (non-negative populations)
	SEI.new[SEI.new < 0] <- 0
	
	# SEI < N.hat, otherwise put N.hat instead (can not have in district more infected people than present)
	SEI.new[SEI.new > rep(N.hat,3)] <- rep(N.hat,3)
	
	return(SEI.new)
}