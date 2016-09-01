###################################################################################################
## original code by Robert Dorazio
## modified by Vira Koshkina
## Integrated Species Distribution Models: Combining presence-only data and presence-absence data with imperfect detection
## utility functions for the file POPA-functions
## 19/08/2016
###################################################################################################

#likelihood functions
# Utility functions
logit = function(pp) { log(pp) - log(1-pp) }
expit = function(eta) {1/(1+exp(-eta))}

# function that calculates a probability of occupancy for each location in the dataset X
predict=function(mymodel, X){
	beta=mymodel$coefs[1:ncol(X),'Value']

	lambda=exp(X %*% beta)
	# psi =1- exp(-lambda*area)
	return(lambda)
}


#Function that fits IPP model
po.ipp=function(X.po, W.po,X.back, W.back){

	beta.names=colnames(X.back)
	beta.names[1]='beta0'

	alpha.names=colnames(W.back)
	alpha.names[1]='alpha0'

	par.names=c(beta.names,	alpha.names)



	minrecipCondNum = 1e-6
	paramGuess = c(rep(.1, ncol(X.po)), rep(-.1, ncol(W.po)))


	fit.po = optim(par=paramGuess, fn=negLL.po, method='BFGS', hessian=FALSE
								 , X.po=X.po, W.po=W.po,X.back=X.back,W.back=W.back ) # params for likelyhood function


	# calculating se with Hessian matrix
	recipCondNum.po = NA
	se.po = rep(NA, length(fit.po$par))
	if (fit.po$convergence==0) {
		hess = ObsInfo.po(fit.po$par, X.po, W.po,X.back, W.back)
		ev = eigen(hess)$values
		recipCondNum.po = ev[length(ev)]/ev[1]
		if (recipCondNum.po>minrecipCondNum) {
			vcv = chol2inv(chol(hess))
			se.po = sqrt(diag(vcv))
		}
	}

	#printing PO results
	tmp= data.frame(par.names,fit.po$par,se.po)
	names(tmp)=c('Parameter name', 'Value', 'Standard error')
	p=NULL
	p$coefs=tmp
	p$convergence=fit.po$convergence
	p$optim_message=fit.po$message
	p$value=fit.po$value
	# print("Estimated parameters beta and alpha", quote=FALSE)
	# print(p)
	return(p)
}

#Function that fits Mackenzie model
pa.mackenzie=function(X.pa,W.pa,y.pa){

	beta.names=colnames(X.pa)
	beta.names[1]='beta0'
	# find sites with at least one detection
	y.pa.pres = y.pa[rowSums(y.pa)>=1,] #detection/non detection matrix for sites with detection in at least one of the surveys


	alpha.names.pa=NULL
	for (i in 1:(dim(W.pa)[3])){
		alpha.names.pa[i]=paste("alpha",as.character(i-1), ".pa", sep="")}

	par.names.pa=c(beta.names,alpha.names.pa)


	#Analyzing Presence-Absence data ------------------------------------------------

	minrecipCondNum = 1e-6

	paramGuess = c(rep(.2, dim(X.pa)[2]), rep(.1, dim(W.pa)[3]))
	fit.pa = NA
	fit.pa = optim(par=paramGuess, fn=negLL.pa, method='BFGS', hessian=TRUE,y.pa.pres=y.pa.pres,y.pa=y.pa, X.pa=X.pa, W.pa=W.pa)

	# calculating se with Hessian matrix
	recipCondNum.pa = NA
	se.pa = rep(NA, length(fit.pa$par))
	if (fit.pa$convergence==0) {
		hess = fit.pa$hessian
		ev = eigen(hess)$values
		recipCondNum.pa = ev[length(ev)]/ev[1]
		if (recipCondNum.pa>minrecipCondNum) {
			vcv = chol2inv(chol(hess))
			se.pa = sqrt(diag(vcv))
		}
	}

	#print PA results
	tmp=data.frame(par.names.pa,fit.pa$par,se.pa)
	names(tmp)=c('Parameter name', 'Value', 'Standard error')
	p=NULL
	p$coefs=tmp
	p$convergence=fit.pa$convergence
	p$optim_message=fit.pa$message
	p$value=fit.pa$value
	return(p)

}

#Function that fits Combined data model
poANDpa.combined=function(X.po, W.po,X.back, W.back,X.pa,W.pa,y.pa){

	beta.names=colnames(X.back)
	beta.names[1]='beta0'

	alpha.names=colnames(W.back)
	alpha.names[1]='alpha0'

	alpha.names.pa=NULL
	for (i in 1:(dim(W.pa)[3])){
		alpha.names.pa[i]=paste("alpha",as.character(i-1), ".pa", sep="")}

	par.names=c(beta.names,	alpha.names, alpha.names.pa)

	y.pa.pres = y.pa[rowSums(y.pa)>=1,] #detection/non detection matrix for sites with detection in at least one of the surveys
	minrecipCondNum = 1e-6
	paramGuess = c(rep(0, dim(X.po)[2]),rep(0, dim(W.po)[2]), rep(0, dim(W.pa)[3]))
	fit.poANDpa = optim(par=paramGuess, fn=negLL.poANDpa, method='BFGS', hessian=TRUE
											,y.pa.pres=y.pa.pres,y.pa=y.pa, X.po=X.po, W.po=W.po, X.back=X.back, W.back=W.back, X.pa=X.pa, W.pa=W.pa )

	# calculating se with Hessian matrix
	recipCondNum.poANDpa = NA
	se.poANDpa = rep(NA, length(fit.poANDpa$par))
	if (fit.poANDpa$convergence==0) {
		hess = fit.poANDpa$hessian
		ev = eigen(hess)$values
		recipCondNum.poANDpa = ev[length(ev)]/ev[1]
		if (recipCondNum.poANDpa>minrecipCondNum) {
			vcv = chol2inv(chol(hess))
			se.poANDpa = sqrt(diag(vcv))
		}
	}

	#print Combined Data results
	tmp=data.frame(par.names,fit.poANDpa$par,se.poANDpa)
	names(tmp)=c('Parameter name', 'Value', 'Standard error')

	p=NULL
	p$coefs=tmp
	p$convergence=fit.poANDpa$convergence
	p$optim_message=fit.poANDpa$message
	p$value=fit.poANDpa$value
	return(p)
}

# negative loglikelihood function for Poisson point process
negLL.pp = function(param) {

	beta = param[1:dim(X.pp)[2]]
	lambda = exp(X.back %*% beta)
	mu = lambda * area.back

	logL.pp = sum(X.pp %*% beta) - sum(mu)

	(-1)*sum(logL.pp)
}

# negative loglikelihood function for thinned Poisson point process
negLL.po = function(param, X.po, W.po,X.back, W.back) {

	beta = param[1:dim(X.po)[2]]
	alpha = param[(dim(X.po)[2]+1):(dim(X.po)[2]+dim(W.po)[2])]
 	# dim(X.back)
 	# length(beta)
 	# length(area.back)
	lambda = exp(X.back %*% beta)
	mu = lambda * area.back
	p = expit(W.back %*% alpha)

	logL.po = sum(X.po %*% beta) + sum(W.po %*% alpha) - sum(log(1 + exp(W.po %*% alpha))) - sum(mu*p)

	(-1)*sum(logL.po)
}

# Observed hessian matrix of negative loglikelihood function for thinned Poisson point process
ObsInfo.po = function(param, X.po,W.po,X.back, W.back) {

	beta = param[1:dim(X.back)[2]]
	alpha = param[(dim(X.back)[2]+1):(dim(X.back)[2]+dim(W.back)[2])]

	lambda = exp(X.back %*% beta)
	mu = lambda * area.back
	p = expit(W.back %*% alpha)

	p.po = expit(W.po %*% alpha)

	nxcovs = length(beta)
	nwcovs = length(alpha)

	nparams = nxcovs + nwcovs
	Hmat = matrix(nrow=nparams, ncol=nparams)

	#  beta partials
	for (i in 1:nxcovs) {
		for (j in 1:i) {
			Hmat[i,j] = sum(X.back[,i] * X.back[,j] * mu * p)
			Hmat[j,i] = Hmat[i,j]
		}
	}

	# alpha partials
	for (i in 1:nwcovs) {
		for (j in 1:i) {
			Hmat[nxcovs+i, nxcovs+j] = sum(W.back[,i] * W.back[,j] * mu * p * ((1-p)^3) * (1 - exp(2 * W.back %*% alpha)) ) + sum(W.po[,i] * W.po[,j] * p.po * (1-p.po))
			Hmat[nxcovs+j, nxcovs+i] = Hmat[nxcovs+i, nxcovs+j]
		}
	}

	# alpha-beta partials
	for (i in 1:nwcovs) {
		for (j in 1:nxcovs) {
			Hmat[nxcovs+i, j] = sum(X.back[,j] * W.back[,i] * mu * p * (1-p))
			Hmat[j, nxcovs+i] = Hmat[nxcovs+i, j]
		}
	}

	Hmat
}

# Expected hessian matrix of negative loglikelihood function for thinned Poisson point process
FisherInfo.po = function(param) {

	beta = param[1:dim(X.back)[2]]
	alpha = param[(dim(X.back)[2]+1):(dim(X.back)[2]+dim(W.back)[2])]

	lambda = exp(X.back %*% beta)
	mu = lambda * area.back
	p = expit(W.back %*% alpha)


	nxcovs = length(beta)
	nwcovs = length(alpha)

	nparams = nxcovs + nwcovs
	Hmat = matrix(nrow=nparams, ncol=nparams)

	#  beta partials
	for (i in 1:nxcovs) {
		for (j in 1:i) {
			Hmat[i,j] = sum(X.back[,i] * X.back[,j] * mu * p)
			Hmat[j,i] = Hmat[i,j]
		}
	}

	# alpha partials
	for (i in 1:nwcovs) {
		for (j in 1:i) {
			Hmat[nxcovs+i, nxcovs+j] = sum(W.back[,i] * W.back[,j] * mu * p * ((1-p)^3) * (1 - exp(2 * W.back %*% alpha)) ) + sum(W.back[,i] * W.back[,j] * p * (1-p) * mu * p)
			Hmat[nxcovs+j, nxcovs+i] = Hmat[nxcovs+i, nxcovs+j]
		}
	}

	# alpha-beta partials
	for (i in 1:nwcovs) {
		for (j in 1:nxcovs) {
			Hmat[nxcovs+i, j] = sum(X.back[,j] * W.back[,i] * mu * p * (1-p))
			Hmat[j, nxcovs+i] = Hmat[nxcovs+i, j]
		}
	}

	Hmat
}

# negative loglikelihood function for Mackenzie model
negLL.pa = function(param, y.pa.pres, y.pa,X.pa,W.pa) {

	beta = param[1:dim(X.pa)[2]]
	alpha = param[(dim(X.pa)[2]+1):(dim(X.pa)[2]+dim(W.pa)[3])]

	#temp --------------------------
	area.pa=1

	lambda.pa = exp(X.pa %*% beta)
	psi =1- exp(-lambda.pa*area.pa)


	mean(lambda.pa)
	mean(psi)

	p.pa = matrix(nrow=dim(W.pa)[1], ncol=J.pa)

	for (j in 1:J.pa) {
		p.pa[, j] = expit(as.matrix(W.pa[,j,], nrow=dim(W.pa)[1]) %*% alpha)
	}


	# prob of detection for sites with presence at least in one of the surveys, and with no presence detected
	p.pa.pres=p.pa[rowSums(y.pa)>=1,]
	p.pa.non.pres=p.pa[rowSums(y.pa)==0,]

	# prob of occupancy for sites with presence at least in one of the surveys, and with no presence detected
	psi.pres=psi[rowSums(y.pa)>=1,]
	psi.non.pres=psi[rowSums(y.pa)==0,]


	#If there is only one site with no observed animals R automatically turns p.pa.non.pres in a vector while we need it in a form of a matrix with 1 row and J.pa (number of surveys rows) for the function rowProds to work
	if (length(p.pa.non.pres)==J.pa) {dim(p.pa.non.pres)=c(1,J.pa)}

	pa.pres=sum(log(psi.pres)+rowSums(y.pa.pres*log(p.pa.pres)+(1-y.pa.pres)*log(1-p.pa.pres)))


	pa.non.pres=0
	if (!is.null(psi.non.pres))	{
		pa.non.pres=sum(log(psi.non.pres*rowProds(1-p.pa.non.pres)+1-psi.non.pres))
	}

	-(pa.pres+pa.non.pres)

}



negLL.poANDpa = function(param,y.pa.pres,y.pa, X.po, W.po, X.back, W.back, X.pa, W.pa )  {

	param.po = param[1:(dim(X.po)[2]+dim(W.po)[2])]
	param.pa = param[c(1:dim(X.po)[2], (dim(X.po)[2]+dim(W.po)[2]+1):(dim(X.po)[2]+dim(W.po)[2]+dim(W.pa)[3]))]
	negLL.po(param.po, X.po, W.po,X.back, W.back ) + negLL.pa(param.pa,y.pa.pres,y.pa,X.pa,W.pa)
}

