library(parallel)
suppressMessages(library(hash))

options(arn=-1)

args <- commandArgs()

inputFile <-args[6];
no_cores  <-args[7];

score <- function(score) {
	return(-10*log(score)/log(10))
	}

unscore <- function(phred) {
	return(10^(-phred/10))
	}

dataProb <- function(target, nref, nalt, error) {
	if (is.na(error)) {return(1)}
	if (target == 0) {return(dbinom(nalt, nref+nalt, 1-error))}
	if (target == 1) {return(dbinom(nalt, nref+nalt, 0.5))}
	if (target == 2) {return(dbinom(nref, nref+nalt, 1-error))}
	}

genProb <- function(target, nref, nalt, error, p) {
	prior <- c(p^2, 2*p*(1-p), (1-p)^2)
	priorD <- 0
	for(ti in c(0,1,2)) {
		priorD <- priorD + dataProb(ti, nref, nalt, error)*prior[ti+1]
		}
	return(dataProb(target, nref, nalt, error)*prior[target+1]/priorD)
	}

call <- function(nref, nalt, error, p) {
	if (is.na(error)) {return(NA)}
	probs <- c(genProb(0,nref,nalt,error,p), genProb(1,nref,nalt,error,p), genProb(2,nref,nalt,error,p))
	if (probs[1] == "NaN") {return(NA)}
	if (probs[2] == "NaN") {return(NA)}
	if (probs[3] == "NaN") {return(NA)}
	probMax <- max(probs)
	for (ci in c(0,1,2)) {
		if (probs[ci + 1] == probMax) {
			return(ci)
			}
		}
	}

RescaleProb <- function(target, prob, p) {
	if (target == 0) {
		return( (prob - p*p)/(1 - p*p) )
		}
	if (target == 1) {
		return( (prob - 2*p*(1-p))/(1 - 2*p*(1-p)) )
		}
	if (target == 2) {
		return( (prob - 2*(1-p)*(1-p))/(1 - (1-p)*(1-p)) )
		}
	}

indicatorG <- function(target, nref, nalt, error, p) {
	if (is.na(error)) {return(0)}
	if (is.na(call(nref, nalt, error, p))) {return(0)}
	if (call(nref, nalt, error, p) == target) {
		prob <- RescaleProb(target, genProb(target, nref, nalt, error, p), p)
		if (prob < 0) {
			return(0)
			} else {
			phred <- -10*log(1 - prob)/log(10)
			if (phred > 20) {
				return(1)
				} else {return(0)}
			}
		} else {
		return(0)
		}
	}

approxSens <- function(cover, p, error) {
	return(exp(-log((p**(1/250))*(error**0.125)*0.8823648)*cover))
	}

coverSet <- function(cover, conf_pp) {
	k_min	<- cover
	k_max	<- cover
	k_m	<- cover
	if (k_m == 0) {} else {
		while (TRUE) {
			if (k_min > 0) {
				k_min <- k_min - 1
				pp_min <- ppois(k_min, k_m)
				if (k_min == 0) {pp_min <- 0}
				if (ppois(k_max, k_m) - pp_min > conf_pp) {
					break
					}
				}
			k_max <- k_max + 1
			pp_min <- ppois(k_min, k_m)
			if (k_min == 0) {pp_min <- 0}
			if (ppois(k_max, k_m) - pp_min > conf_pp) {
				break
				}
			}
		}
	k_normal <- 0
	if (k_m == 0) {k_normal <- 1} else {
		for (k in (k_min:k_max)) {
			k_normal <- k_normal + dpois(k, k_m)
			}
		}
	cSet <- matrix(,0,2)
	for (k in k_min:k_max) {
		cSet <- rbind(cSet, c(k, dpois(k, k_m)/ k_normal))
		}
	return(cSet)
	}

singleRun <- function(element) {
	p	<- as.numeric(element[1])
	error	<- as.numeric(element[3])
	cover	<- as.numeric(element[2])
	name	<- as.character(element[4])
	sens_i	<- 0
	if (cover > 0) {
		if (approxSens(cover, p, error) < 1e-6) {
			return(list(name, p, p))
			}
		} else {
		return(list(name, 0, p))
		}
	cSet <- coverSet(cover, 0.99)
	for (ki in (1:nrow(cSet))) {
		k	<- cSet[ki,1]
		k_m	<- cover
		p_k	<- cSet[ki,2]
		for (j in (0:k)) {
			nref	<- k - j
			nalt	<- j
			errorL	<- unscore(error)
			sens_i	<- sens_i +
				p_k*dataProb(0, k - j, j, errorL) *
					indicatorG(0, k-j, j, errorL, p) *
					2 * p^2 +
				p_k*dataProb(1, k - j, j, errorL) *
					indicatorG(1, k-j, j, errorL, p) *
					2 * p * (1 - p)
			}
		}
	return(list(name, sens_i/2, p))
	}

#--------------------------------------------
#--- HEAD
#--------------------------------------------

input <- scan(inputFile, what="", sep="\n", quiet=TRUE)
givenData <- matrix(,0,3)       # allele frequency, km - coverage, error rate
givenName <- c()
dataList  <- list()
n <- 1
cover <- hash()
for (i in (1:(length(input)/3))) {
	name <- as.character(unlist(strsplit(input[i*3-2], "\t"))[1])
	p <- as.numeric(unlist(strsplit(input[i*3-2], "\t"))[3])
	x <- c()
	x1 <- strsplit(input[i*3-1], "\t")
	if (length(unlist(x1)) > 2) {
		x1 <- as.numeric(unlist(x1)[3:length(unlist(x1))])
		} else {
		x1 <- NULL
		}
	x2 <- strsplit(input[i*3], "\t")
	if (length(unlist(x2)) > 2) {
		x2 <- as.numeric(unlist(x2)[3:length(unlist(x2))])
		if (length(c(x1,x2)) > 0) {
			if (indicatorG(0, length(x1), length(x2), unscore(mean(c(x1,x2))), p) > 0) {
				x2 <- c()
				}
			if (indicatorG(1, length(x1), length(x2), unscore(mean(c(x1,x2))), p) > 0) {
				x2 <- c()
				}
			}
		} else {
		x2 <- NULL
		}
	element <- c()
	if (length(c(x1,x2)) == 0) {
		element <- c(p, 0, NA)
		givenName <- c(givenName, name)
		} else {
		element <- c(p, length(c(x1,x2)), mean(c(x1,x2)))
		givenName <- c(givenName, name)
		}
	givenData <- rbind(givenData, c(element))
	dataList[[n]] <- list(element[1], element[2], element[3], name)
	cover[[name]] <- element[2]
	n <- n + 1
	}

c1 <- makeCluster(no_cores, type="FORK")
result <- parLapply(c1, dataList, singleRun)
stopCluster(c1)

sens <- 0
for(i in (1:length(result))) {
	cat(as.character(result[[i]][1])," ",
		as.numeric(result[[i]][2])," ",
		as.numeric(result[[i]][3])," ",
		cover[[as.character(result[[i]][1])]],"\n"
		)
	sens <- sens + as.numeric(result[[i]][2])
	}
cat(sens,"\n")























