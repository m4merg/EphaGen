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

indicatorG <- function(nref, nalt, error) {
	if (is.na(error)) {return(0)}
	prob = poisson.test(nalt, (nref + nalt)*error)$p.value;
	phred <- -10*log(prob)/log(10)
	if (phred > 100) {return(1)} else {return(0)}
	}

approxSens <- function(error) {
	return(((error)/550)**(-2))
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
		if (cover > approxSens(error)) {
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
		A_range <- ceiling(k/2) - round(k*0.05) + 1
		for (j in (round(k*0.05):ceiling(k/2))) {
			nref	<- k - j
			nalt	<- j
			errorL	<- unscore(error)
			sens_i	<- sens_i +
				p_k*(1/A_range)*
					indicatorG(k-j, j, errorL) * p
			}
		}
	return(list(name, sens_i, p))
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
			if (indicatorG(length(x1), length(x2), unscore(mean(c(x1,x2)))) > 0) {
				x2 <- c()
				}
			if (indicatorG(length(x1), length(x2), unscore(mean(c(x1,x2)))) > 0) {
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























