options(arn=-1)

args <- commandArgs()

inputFile <-args[6];

input <- scan(inputFile, what="", sep="\n")
givenData <- matrix(,0,3)	# allele frequency, km - coverage, error rate
givenName <- c()
for (i in (1:(length(input)/3))) {
	name <- unlist(strsplit(input[i*3-2], "\t"))[1]
	p <- as.numeric(unlist(strsplit(input[i*3-2], "\t"))[3])
	x <- c()
	x1 <- strsplit(input[i*3-1], "\t")
	if (length(unlist(x1)) > 2) {
		x1 <- as.numeric(unlist(x1)[3:length(unlist(x1))])
		x1 <- x1[!x1 %in% boxplot.stats(x1, coef=0.5)$out]
		} else {
		x1 <- NULL
		}
	x2 <- strsplit(input[i*3], "\t")
	if (length(unlist(x2)) > 2) {
		x2 <- as.numeric(unlist(x2)[3:length(unlist(x2))])
		x2 <- x2[!x2 %in% boxplot.stats(x2, coef=0.5)$out]
		} else {
		x2 <- NULL
		}
	if (length(c(x1,x2)) == 0) {
		givenData <- rbind(givenData, c(p, 0, NA))
		givenName <- c(givenName, name)
		} else {
		givenData <- rbind(givenData, c(p, length(c(x1,x2)), mean(c(x1,x2))))
		givenName <- c(givenName, name)
		}
	}

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

sens <- 0
#print(givenData)
for (i in (1:nrow(givenData))) {
	sens_i	<- 0
	p	<- givenData[i,1]
	if (givenData[i,2] > 0) {
		if (approxSens(givenData[i,2], givenData[i,1], givenData[i,3]) < 1e-6) {
			sens_i <- 2*givenData[i,1]
			sens <- sens + sens_i
			cat(givenName[i],sens_i/2,p^2 + p * (1-p),"\n")
			next
			}
		}
	cSet <- coverSet(givenData[i,2], 0.99)
	for (ki in (1:nrow(cSet))) {
		k	<- cSet[ki,1]
		k_m	<- givenData[i,2]
		p_k	<- cSet[ki,2]
#		cat(k, "(",k_m,") - ", p_k,"\n")
		for (j in (0:k)) {
			nref	<- k - j
			nalt	<- j
			error	<- unscore(givenData[i,3])
#			cat(i," - ",p,"\t",k, "\t", j, "\t", indicatorG(0, nref, nalt, error, p),"-",dataProb(0, nref, nalt, error),"\t",indicatorG(1, nref, nalt, error, p),"-",dataProb(1, nref, nalt, error),"\n")
#			cat(nref,"\t",nalt,"\t",error,"\t",p,"\t",genProb(0, nref, nalt, error, p),"\t",genProb(1, nref, nalt, error, p),"\n")
			sens_i <- sens_i + 
				p_k*dataProb(0, k-j, j, error) * 
					indicatorG(0, k-j, j, error, givenData[i,1]) *
					2 * p^2 + 
				p_k*dataProb(1, k-j, j, error) *
					indicatorG(1, k-j, j, error, givenData[i,1]) *
					2 * p * (1-p)
			}
		}
	cat(givenName[i],sens_i/2,p^2 + p * (1-p),"\n")	
	sens <- sens + sens_i
	}
sens <- sens/2
cat(sens,"\n")
























