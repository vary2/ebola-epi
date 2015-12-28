clearall = T
if (clearall) rm(list = ls(all.names = TRUE))
echo = T
#
# 	Prepare incidence data
#

# Incidence district-level data for Western African countries (GU, LI, SL)
wa <- read.csv('../data/WA-WHO-2015-12-16.csv')	
	
# Omit cumulative rows and rows without location
# Note: These are summary rows in the original csv WHO file
wa <- wa[wa$Indicator.type!='Cumulative' & wa$Location != '',]

# Take confirmed records only from Patient database source
# Note: Data from WHO include various sources
wa <- wa[wa$Ebola.data.source == 'Patient database' & wa$Case.definition == 'Confirmed',]

# Change NA values to 0s
wa[is.na(wa$Cases),]$Cases <- numeric(dim(wa[is.na(wa$Cases),])[1])

# Order file by time
wa <- wa[order(wa$Time),]

# Districts naming convenvtion
wa$Location <- unlist(lapply(as.character(wa$Location), camelCase))

# Prepare smaller datasets specifically for GU, LI, SL
wa <- wa[,c('Country', 'Location', 'Cases', 'Time')]
nweeks <- sum(wa$Location == 'Bo')

gu <- wa[wa$Country == 'Guinea',]
li <- wa[wa$Country == 'Liberia',]
sl <- wa[wa$Country == 'Sierra Leone',]

#
# 	Prepare districts/regions data
#
sl.name <- 'Sierra Leone'
# short name is needed for googleAPI (should be fixed by proper URLencode-this works for now)
sl.name.short <- 'SierraLeone'

# Sierra Leone districts information (Area, Population)
sl.distr.inf <- read.csv('../data/sl.distr.inf.csv')
sl.distr.inf <- sl.distr.inf[sl.distr.inf$Country == sl.name,]
sl.distr.n <- dim(sl.distr.inf)[1]

# Create distric tseries data
sl.distr.tseries <- array(0,ncol = nweeks, nrow = sl.distr.n)
i = 1
for (distr in sl.distr.inf$District){
	sl.distr.tseries[i,] <- sl[sl$Location == distr,]$Cases
	i = i + 1
}

#
#	Create connectivity matrix
#

# Sierra Leone districts adjacency matrix
sl.distr.adj <- read.csv('../data/sl.distr.adj.csv',header = T)
row.names(sl.distr.adj) <- colnames(sl.distr.adj)

# create edgelist from adjacency matrix
sl.distr.el <- which(sl.distr.adj == 1, arr.ind = T,useNames = F)

# create matrix of interdistrict distances
# distances of two districts is approx. by proxies (distr. capitols)
sl.distr.dist <- as.matrix(sl.distr.adj)
sl.distr.time <- as.matrix(sl.distr.adj)

# load function for Google API
source('fetchGoogle.r',echo = F)

# iteration through edges - to get travel distances and times between distr. capitols
for (i in 1:nrow(sl.distr.el)){
	# district index
	distr1 <- sl.distr.el[i,1]
	distr2 <- sl.distr.el[i,2]
	# district capitols
	capital1 <- sl.distr.inf[sl.distr.inf$Country == sl.name,]$Capital[distr1]
	capital2 <- sl.distr.inf[sl.distr.inf$Country == sl.name,]$Capital[distr2]
	
	#
	if (echo) print(paste('Getting: ',paste(capital1,sl.name.short,sep = ','),' <-> ',
			    		paste(capital2,sl.name.short,sep = ',')))
	
	dist <- getDistance(paste(capital1,sl.name.short,sep = ','),
			    paste(capital2,sl.name.short,sep = ','), 
			    googleKey)
	if (echo) print(paste('Got: distance=',dist[1],'; time=',dist[2],sep = ''))
	# road distance in kilometers
	sl.distr.dist[distr1,distr2] <- dist[1]
	sl.distr.dist[distr2,distr1] <- dist[1]
	# time travel in hours
	sl.distr.time[distr1,distr2] <- dist[2]
	sl.distr.time[distr2,distr1] <- dist[2]
}

# Based on the publication:
r.ave <- sqrt(sum(sl.distr.inf[sl.distr.inf$Country == sl.name,'Area'])/pi)
sl.distr.dist <- sl.distr.dist/r.ave
diag(sl.distr.dist) <- sqrt(sl.distr.inf[sl.distr.inf$Country == sl.name,'Area']/pi)

# Create N1, N2 matrices - to be used in construction of connectivity matrix
n = length(sl.distr.inf[sl.distr.inf$Country == sl.name,]$Population)
sl.pop.tot <- sum(sl.distr.inf[sl.distr.inf$Country == sl.name,]$Population)
r.ave <- sqrt(sum(sl.distr.inf[sl.distr.inf$Country == sl.name,]$Area)/pi)
sl.N1 <- matrix(rep(sl.distr.inf[sl.distr.inf$Country == sl.name,]$Population,n), ncol = n)
sl.N2 <- matrix(rep(sl.distr.inf[sl.distr.inf$Country == sl.name,]$Population,n), ncol = n, byrow = T)
sl.N <- sl.N1[,1]
sl.N1 <- sl.N1 / sl.pop.tot
sl.N2 <- sl.N2 / sl.pop.tot


# function for constructing connectivity matrix for certain tau1, tau2, ro
create.Cm <- function(N1, N2, distr.dist, tau1, tau2, ro){
	Cm <- N1^tau1 * N2^tau2 / distr.dist^ro
	Cm[is.infinite(Cm)] <- 0
	Cm <- Cm %*% diag(1 / rowSums(Cm))
	return(Cm)
}

