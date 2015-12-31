
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
	
	first.nonzeros <- apply(distr.tseries,1, function(x) head(x[x!=0],1))
	obs.mean <- apply(distr.tseries[,trn.range],1,mean)	# mean of whole trn.range
	obs.sd   <- apply(distr.tseries[,trn.range],1,sd)	# std. dev. of whole trn.range (for init. draw)
	
# 	tmp=rep(0,nobs)
# 	for (i in 4:nobs){
# 		tmp[i]=mean(distr.tseries[,(i-3):(i-1)]);
# 	}
# 	obs.vars = sum((obs.mean-tmp)^2)/(nobs-1);
	
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
						a=0, b=Inf, mean = first.nonzeros, sd = 2*obs.sd), ncol = nens)
	# I0
	xprior[1:ndist+2*ndist,] <- matrix(rtruncnorm(nens*ndist, 
						a=0, b=Inf, mean = first.nonzeros, sd = 2*obs.sd), ncol = nens)
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
		if (any(is.na(xprior[4*ndist+3:4, j]))){
			stop('tau1/2 is NA !')
		}
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
			if (any(is.na(xprior[,j]))){
				stop('xprior[,j] is NA')
			}
			xobs[,j] <- rk4(t1 = i,t2 = i+1, dt = 1, N = N, xprior.n = xprior[,j], Cn = Cns[,,j])
		}
		#xobs <- check.bounds.xnew(xobs, ndist, N)

		# update filter
		#var.obs <- apply(xobs,1,var)	# here will go inflation OEV
		var.obs <- obs.vars
		var.prior <- apply(xprior,1,var)
		xprior.mean <- apply(xprior,1,mean)
		z <- c(numeric(2*ndist), distr.tseries[,i], numeric(nvar-3*ndist))	# new incidence data
											# Note: we observe only new I
		den <- var.obs + var.prior +0.3
		

		# eq. (4.4):
		xpost.mean <- var.obs/den*xprior.mean + var.prior/den*z		# vector (nvar,1)

		# eq. (4.5):
		ones <- rep(1,nens)
		xpost <- xpost.mean %*% t(ones) + (sqrt(var.obs/den) %*% t(ones))*
				(xprior-xprior.mean %*% t(rep(1,nens)))		# matrix (nvar, nens)
		
		if (any(is.na(xpost))){
			stop('xpost is NA')
		}
		xprior <- xpost

		# update Cns array
		for (j in 1:nens){
			if (any(is.na(xprior[4*ndist+3:4, j]))){
				stop('tau1/2 is NA !')
			}
			Cns[,,j] <- create.Cm(N1 = N1, N2 = N2, distr.dist = distr.dist, 
						xpost[4*ndist+3, j],
						xpost[4*ndist+4, j],
						xpost[4*ndist+5, j])	# last 3 are tau1, tau2, ro
		}
	}
	return(xpost)
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
	if (any(is.na(N.hat))){
		stop('N.hat is NA')
	}
	if (any(is.na(SEI))){
		stop('SEI is NA')
	}
	# propagate for one ensemble
	xpost.n <- xprior.n
	t.vec <- seq(t1, t2, by=dt)
	for (t in t.vec){
		k1s <- SEIRfunc(SEI = SEI       , N.hat = N.hat, N = N, beta = beta, alpha = alpha, Z = Z, D = D, Cn = Cn)
		k2s <- SEIRfunc(SEI + 0.5*dt*k1s, N.hat = N.hat, N = N, beta = beta, alpha = alpha, Z = Z, D = D, Cn = Cn)
		k3s <- SEIRfunc(SEI + 0.5*dt*k2s, N.hat = N.hat, N = N, beta = beta, alpha = alpha, Z = Z, D = D, Cn = Cn)
		k4s <- SEIRfunc(SEI +   1*dt*k3s, N.hat = N.hat, N = N, beta = beta, alpha = alpha, Z = Z, D = D, Cn = Cn)
		SEI <- SEI + dt*(k1s + 2*k2s + 2*k3s + k4s)/6
	}
	
	xpost.n[1:(3*ndist)] <- SEI
	return(xpost.n)
}

SEIRfunc <- function(SEI, N.hat, N, beta, alpha,  Z, D, Cn){
	ndist <- length(SEI)/3
	S <- SEI[1:ndist] 
	E <- SEI[1:ndist+ndist] 
	I <- SEI[1:ndist+2*ndist]

	h1 <- S*sum(t(Cn) %*% (beta*I/N.hat))		# I.hat := I
	h2 <- E/Z 				# Z != 0 !
	SEI[1:ndist] <- -h1 - alpha		# S
	SEI[1:ndist+ndist] <- h1 - h2 + alpha	# E
	SEI[1:ndist+2*ndist] <- h2 - I/D	# I	D != 0 !
	if (any(is.na(SEI))){
		View(Cn)
		print('----- h1:')
		print(h1)
		print('----- beta:')
		print(beta)
		print('----- N.hat:')
		print(N.hat)
		print('----- S:')
		print(S)
		print('----- E:')
		print(E)
		print('----- I:')
		print(I)
		print('-----')
		stop('Na in SEIR func')
	}
	#SEI <- check.bounds.SEI(SEI, N)
	return(SEI)
}

camelCase <- function(x) {
	# "I like pizza" to "ILikePizza"
	s <- tolower(x)
	s <- strsplit(s, " ")[[1]]
	paste(toupper(substring(s, 1,1)), substring(s, 2),
	      sep="", collapse="")
}

check.bounds.xnew <- function(xnew, ndist, N){
	# posterior corrections for xnew (xnew is matrix)
	
	nvar <- dim(xnew)[1]
	nens <- dim(xnew)[2]

	Ss <- xnew[1:ndist,]
	Es <- xnew[1:ndist+ndist,]
	Is <- xnew[1:ndist+2*ndist,]
	Rs <- xnew[1:ndist+3*ndist,]
	Zs <- xnew[4*ndist+1,]
	Ds <- xnew[4*ndist+2,]
	taus <- xnew[4*ndist + 3:4,]
	ros <- xnew[4*ndist+5,]

	# Parameter bounds
# 	Rs[Rs < 0.5] = 0.5
# 	Rs[Rs > 3.5] = 3.5
# 	Zs[Zs < 2] = 2
# 	Zs[Zs > 14] = 14
# 	Ds[Ds < 5] = 5
# 	Ds[Ds > 14] = 14
# 	taus[taus < 0.3] = 0.3
# 	taus[taus < 0.7] = 0.3
# 	ros[ros < 2] = 2
# 	ros[ros > 8] = 8
	
	# S
	Ns <- N %*%t(rep(1,nens))
	Ss[Ss > Ns] <- Ns[Ss > Ns]-10
	
	# E
	E.median <- apply(Es,1,median) %*% t(rep(1,nens))
	Es[Es > Ns] <- E.median[Es > Ns]
	
	# I
	I.median <- apply(Is,1,median) %*% t(rep(1,nens))
	Is[Is > Ns] <- I.median[Is > Ns]
	
	xnew <- rbind(Ss,Es,Is,Rs,Zs, Ds, taus, ros)
	
	return(xnew)
}

check.bounds.SEI <- function(SEI, N){
	# SEI is vector
	ndist <- length(SEI)/3
	
	S <- SEI[1:ndist]
	E <- SEI[1:ndist+ndist]
	I <- SEI[1:ndist+2*ndist]
	
	# S
	S[S > N] <- N[S > N] - 10
	S[S <= 0] <- 1
	
	# E
	E[E > N] <- N[E > N] - 10
	E[E <= 0] <- 1
	
	# I
	I[I > N] <- N[I > N] - 10
	I[I <= 0] <- 1
	
	SEI <- rbind(S,E,I)
	return(SEI)
}