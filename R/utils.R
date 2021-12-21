nearest_simulated_norms <- function(n, deg) {
	
	# This extracts the most suitable variable mres.norm.* to be used
	
	if (n > 2150) stop("The largest permitted n is 2150.")
	
	len <- as.character(round(round(min(max(1, n/100), 21))*100))

	if (deg == 0) dg <- "c" else if (deg == 1) dg <- "l"
	
	var.name <- paste("mres.norm.", dg, ".", len, sep="")
	
	get(var.name)
	
}




order_chron <- function(nsp.obj) {
	
	# Order intervals of significance chronologically.
	# nsp.obj - quantity returned by one of the nsp_* functions.

	d <- dim(nsp.obj)
	if (d[2]) {
		nsp.obj.ord <- order(nsp.obj[1,])
		nsp.obj <- nsp.obj[,nsp.obj.ord]
	}	
	
	nsp.obj
	
}



all_dyadic_scans_array <- function(x) {
	
	d <- dim(x)
	n <- d[1]
	
	if (n) {
		
		add.scales <- floor(logb(n, 2))
		shifts <- rep(0, add.scales+1)
		res <- array(x, c(d[1], d[2], add.scales+1))
		if (add.scales) for (j in 1:add.scales) {
			res[1:(n-2^j+1),,(j+1)] <- 2^(-1/2) * (res[1:(n-2^j+1),,j] + res[(2^(j-1)+1):(n-2^j+1+2^(j-1)),,j])
			res[(n-2^j+2):(n),,(j+1)] <- 0
			shifts[j+1] <- 2^j-1
		}		
		
	}
	else {
		res <- array(0, c(d[1], d[2], 0))	
		shifts <- integer(0)
	}
	
	list(res=res, shifts=shifts)
	
}


iter_random_checks_scan_array <- function(ind, ads.array, M, thresh, overlap = FALSE, buffer = 0) {
	
	s <- ind[1]
	e <- ind[2]
	n <- e - s + 1

	indices <- ((ads.array$shifts+1) <= (n/2))
	
	ads.array$res <- ads.array$res[s:e,,indices,drop=F]
	
	ads.array$shifts <- ads.array$shifts[indices]
		
	if (n > 1) {
		
		next.int <- random_checks_scan_2stage_array(c(1,n), ads.array, M, thresh)
		
		if (!is.na(next.int$selected.ind))  {
		
			if (!overlap) {
			
			if (next.int$selected.val[1,1]-buffer >= 2) left <- iter_random_checks_scan_array(c(1, next.int$selected.val[1,1]-buffer), ads.array, M, thresh, overlap, buffer) else left <- matrix(NA, 3, 0)
			if (n - next.int$selected.val[1,2]-buffer >= 1) {
				
				right <- iter_random_checks_scan_array(c(next.int$selected.val[1,2]+buffer, n), ads.array, M, thresh, overlap, buffer)
				if (dim(right)[2]) right <- right + c(rep(next.int$selected.val[1,2]-1+buffer, 2), 0)
				
			}
			else right <- matrix(NA, 3, 0)
			}
			
			else {
				

			if (floor(mean(next.int$selected.val[1,1:2]))-buffer >= 2) left <- iter_random_checks_scan_array(c(1, floor(mean(next.int$selected.val[1,1:2]))-buffer), ads.array, M, thresh, overlap, buffer) else left <- matrix(NA, 3, 0)
			if (n - floor(mean(next.int$selected.val[1,1:2]))-buffer >= 2) {
				right <- iter_random_checks_scan_array(c(floor(mean(next.int$selected.val[1,1:2]))+1+buffer, n), ads.array, M, thresh, overlap, buffer)
				if (dim(right)[2]) right <- right + c(rep(floor(mean(next.int$selected.val[1,1:2]))+buffer, 2), 0)
				
			}
			else right <- matrix(NA, 3, 0)
				
				
			}
			
			
			return(cbind(t(next.int$selected.val), left, right))
			
			
		}
		
		else(return(matrix(NA, 3, 0)))
		
		
	}

	else(return(matrix(NA, 3, 0)))
	
}



random_checks_scan_2stage_array <- function(ind, ads.array, M, thresh) {
		
	s1 <- random_checks_scan_array_1by1(ind, ads.array, M, thresh)
		
	if (!is.na(s1$selected.ind)) {
		
		s <- s1$selected.val[1,1] + ind[1] - 1
		e <- s1$selected.val[1,2] + ind[1] - 1
		
		s2 <- random_checks_scan_array_1by1(c(s,e), ads.array, M, thresh)

		if (!is.na(s2$selected.ind)) {
			
			replacement <- s2$selected.val
			replacement[1,1:2] <- replacement[1,1:2] + s - ind[1]
			s1$selected.val <- replacement
			
		}

	}
	
	s1	
	
}


random_checks_scan_array_1by1 <- function(ind, ads.array, M, thresh) {

	s <- ind[1]
	e <- ind[2]
	n <- e - s + 1
	
	if (n > 1) {

		indices <- ((ads.array$shifts+1) <= (n/2))
	
		ads.array$res <- ads.array$res[s:e,,indices,drop=F]
	
		ads.array$shifts <- ads.array$shifts[indices]
	
		M <- min(M, (n-1)*n/2)

		ind <- grid_intervals_sorted(n, M)

		M <- dim(ind)[2]

		res <- matrix(0, M, 3)

		res[,1:2] <- t(ind)

		zero.check <- TRUE
		j <- 1
		
		while (zero.check && (j <= M)) {
			
			res[j,3] <- check_interval_array(res[j,1:2], ads.array, thresh)
			zero.check <- (res[j,3] == 0)
			j <- j + 1
			
		}

		if (zero.check) {
			
			selected.ind <- NA
			selected.val <- matrix(0, 0, 3)
			
		}

		else {
			
			selected.ind <- j-1
			selected.val <- res[selected.ind,,drop=FALSE]
			
		}

	
	}

	else {
		
		selected.val <- matrix(0, 0, 3)
		selected.ind <- NA
		M <- 0
		
	}
	

	list(selected.ind = selected.ind, selected.val = selected.val, M.eff=max(1, M))

}


all_intervals_flat <- function(n) {
	
	if (n == 2) ind <- matrix(1:2, 2, 1) else {
		M <- (n-1)*n/2	
		ind <- matrix(0, 2, M)
		ind[1,] <- rep(1:(n-1), (n-1):1)
		ind[2,] <- 2:(M+1) - rep(cumsum(c(0, (n-2):1)), (n-1):1)
	}
	ind

}


all_intervals_sorted <- function(n) {
	
	d <- all_intervals_flat(n)
	d.ord <- order(d[2,] - d[1,])
	d[,d.ord, drop=FALSE]
	
}


grid_intervals_sorted <- function(n, M) {
	
	if ((n==2) || (M == 0)) ind <- matrix(c(1, n), 2, 1)
	
	else if (M >= (n-1)*n/2) ind <- all_intervals_sorted(n)
	
	else {
		k <- 1
		while (k*(k-1)/2 < M) k <- k+1
		ind2 <- all_intervals_sorted(k)
		ind2.mx <- max(ind2)
		ind <- round((ind2 - 1) * ((n-1) / (ind2.mx-1)) + 1)
	}	
	
	ind	
}


check_interval_array <- function(ind, ads.array, thresh) {
	
	s <- ind[1]
	e <- ind[2]
	n <- e - s + 1

	indices <- ((ads.array$shifts+1) <= (n/2))
	
	ads.array$res <- ads.array$res[s:e,,indices,drop=F]
	
	ads.array$shifts <- ads.array$shifts[indices]

	dm <- dim(ads.array$res)

	f.con.rhs.core <- matrix(0, 0, 1+2*(dm[2]-1))

	scales <- length(ads.array$shifts)	
	
	for (i in 1:scales) {
		
		f.con.rhs.current <- ads.array$res[1:(n-ads.array$shifts[i]),,i]
		f.con.rhs.current <- cbind(f.con.rhs.current, -f.con.rhs.current[,-1])
		
		f.con.rhs.core <- rbind(f.con.rhs.core, f.con.rhs.current)		
		
	}

	f.con.rhs.core <- rbind(f.con.rhs.core, -f.con.rhs.core)

	f.obj <- c(1, rep(0, 2*(dm[2]-1)))
	f.rhs <- f.con.rhs.core[,1]
	f.con <- f.con.rhs.core
	f.con[,1] <- 1	
	d <- dim(f.con.rhs.core)

	f.dir <- rep(">=", d[1])
	linf <- lpSolve::lp("min", f.obj, f.con, f.dir, f.rhs)$solution[1]
	linf.t <- linf * (linf > thresh)
	linf.t

}


lp_selfnorm <- function(ind, ads.array, selfnorm.array) {
	
	s <- ind[1]
	e <- ind[2]
	n <- e - s + 1

	indices <- ((ads.array$shifts+1) <= (n/2))
	
	ads.array$res <- ads.array$res[s:e,,indices,drop=F]
	
	ads.array$shifts <- ads.array$shifts[indices]

	dm <- dim(ads.array$res)

	f.con.rhs.core <- matrix(0, 0, 1+2*(dm[2]-1))

	scales <- length(ads.array$shifts)	
	
	for (i in 1:scales) {
				
		f.con.rhs.current <- ads.array$res[1:(n-ads.array$shifts[i]),,i]  / selfnorm.array[1:(n-ads.array$shifts[i]),i]
		f.con.rhs.current[!is.finite(f.con.rhs.current)] <- 0
		f.con.rhs.current <- cbind(f.con.rhs.current, -f.con.rhs.current[,-1])
		
		f.con.rhs.core <- rbind(f.con.rhs.core, f.con.rhs.current)		
		
	}

	f.con.rhs.core <- rbind(f.con.rhs.core, -f.con.rhs.core)

	f.obj <- c(1, rep(0, 2*(dm[2]-1)))
	f.rhs <- f.con.rhs.core[,1]
	f.con <- f.con.rhs.core
	f.con[,1] <- 1	
	d <- dim(f.con.rhs.core)

	f.dir <- rep(">=", d[1])
	lp.sol <- lpSolve::lp("min", f.obj, f.con, f.dir, f.rhs)$solution
	lp.sol
	
}


create_selfnorm_array_Vn2est <- function(resid, Vn2est, eps, c = exp(1 + 2 * eps)) {
	
	m <- length(resid)
	
	zz <- all_dyadic_scans_array(matrix(resid^2, m, 1))
	zz$res <- zz$res[,,]
	
	zz.norm <- all_dyadic_scans_array(matrix(1, m, 1))
	zz.norm$res <- zz.norm$res[,,]
		
	(1 + eps) * sqrt(zz$res / zz.norm$res) * log(c * pmax(1, Vn2est / (zz$res * zz.norm$res)))^(1/2 + eps)
	
}


linreg_resid <- function(ind, ads.array) {
	
	s <- ind[1]
	e <- ind[2]

	lmmat <- ads.array$res[s:e,,1]
	
	res <- as.numeric(stats::lm(lmmat[,1] ~ lmmat[,-1]-1)$resid)
	
	if (sum(res^2) == 0) res <- (lmmat[,1] - mean(lmmat[,1]))
	
	if (sum(res^2) == 0) res <- lmmat[,1]
		
	res		
	
}


check_interval_array_selfnorm <- function(ind, ads.array, thresh, Vn2est, eps = 0.03, c = exp(1 + 2 * eps)) {
	
	resid <- linreg_resid(ind, ads.array)
	resid.sna <- create_selfnorm_array_Vn2est(resid, Vn2est, eps, c)
	a <- lp_selfnorm(ind, ads.array, resid.sna)
	sol <- a[1] * (a[1] > thresh)
	sol	
	
}



random_checks_scan_array_selfnorm_1by1 <- function(ind, ads.array, M, thresh, Vn2est, eps = 0.03, c = exp(1 + 2 * eps)) {

	s <- ind[1]
	e <- ind[2]
	n <- e - s + 1
	
	if (n > 1) {

		indices <- ((ads.array$shifts+1) <= (n/2))
	
		ads.array$res <- ads.array$res[s:e,,indices,drop=F]
	
		ads.array$shifts <- ads.array$shifts[indices]
	
		M <- min(M, (n-1)*n/2)

		ind <- grid_intervals_sorted(n, M)

		M <- dim(ind)[2]

		res <- matrix(0, M, 3)

		res[,1:2] <- t(ind)

		zero.check <- TRUE
		j <- 1
		
		while (zero.check && (j <= M)) {
			
			res[j,3] <- check_interval_array_selfnorm(res[j,1:2], ads.array, thresh, Vn2est, eps, c)
			zero.check <- (res[j,3] == 0)
			j <- j + 1
			
		}

		if (zero.check) {
			
			selected.ind <- NA
			selected.val <- matrix(0, 0, 3)
			
		}

		else {
			
			selected.ind <- j-1
			selected.val <- res[selected.ind,,drop=FALSE]
			
		}


	}

	else {
		
		filtered.res <- selected.val <- matrix(0, 0, 3)
		selected.ind <- NA
		M <- 0
		
	}

	list(selected.ind = selected.ind, selected.val = selected.val, M.eff=max(1, M))

}



rcs2sas <- function(ind, ads.array, M, thresh, Vn2est, eps = 0.03, c = exp(1 + 2 * eps)) {
		
	s1 <- random_checks_scan_array_selfnorm_1by1(ind, ads.array, M, thresh, Vn2est, eps, c)
		
	if (!is.na(s1$selected.ind)) {
		
		s <- s1$selected.val[1,1] + ind[1] - 1
		e <- s1$selected.val[1,2] + ind[1] - 1
		
		s2 <- random_checks_scan_array_selfnorm_1by1(c(s,e), ads.array, M, thresh, Vn2est, eps, c)

		if (!is.na(s2$selected.ind)) {
			
			replacement <- s2$selected.val
			replacement[1,1:2] <- replacement[1,1:2] + s - ind[1]
			s1$selected.val <- replacement
			
		}

	}
	
	s1	
	
}


ircs2sas <- function(ind, ads.array, M, thresh, Vn2est, eps = 0.03, c = exp(1 + 2 * eps), overlap = FALSE) {
	
	s <- ind[1]
	e <- ind[2]
	n <- e - s + 1

	indices <- ((ads.array$shifts+1) <= (n/2))
	
	ads.array$res <- ads.array$res[s:e,,indices,drop=F]
	
	ads.array$shifts <- ads.array$shifts[indices]
		
	if (n > 1) {
		
		next.int <- rcs2sas(c(1,n), ads.array, M, thresh, Vn2est, eps, c)
		
		if (!is.na(next.int$selected.ind))  {
		
			if (!overlap) {
			
			if (next.int$selected.val[1,1] >= 2) left <- ircs2sas(c(1, next.int$selected.val[1,1]), ads.array, M, thresh, Vn2est, eps, c, overlap) else left <- matrix(NA, 3, 0)
			if (n - next.int$selected.val[1,2] >= 1) {
				
				right <- ircs2sas(c(next.int$selected.val[1,2], n), ads.array, M, thresh, Vn2est, eps, c, overlap)
				if (dim(right)[2]) right <- right + c(rep(next.int$selected.val[1,2]-1, 2), 0)
				
			}
			else right <- matrix(NA, 3, 0)
			}
			
			else {
				

			if (floor(mean(next.int$selected.val[1,1:2])) >= 2) left <- ircs2sas(c(1, floor(mean(next.int$selected.val[1,1:2]))), ads.array, M, thresh, Vn2est, eps, c, overlap) else left <- matrix(NA, 3, 0)
			if (n - floor(mean(next.int$selected.val[1,1:2])) >= 2) {
				right <- ircs2sas(c(floor(mean(next.int$selected.val[1,1:2]))+1, n), ads.array, M, thresh, Vn2est, eps, c, overlap)
				if (dim(right)[2]) right <- right + c(rep(floor(mean(next.int$selected.val[1,1:2])), 2), 0)
				
			}
			else right <- matrix(NA, 3, 0)
				
				
			}
			
			
			return(cbind(t(next.int$selected.val), left, right))
			
		}
		
		else(return(matrix(NA, 3, 0)))
		
	}

	else(return(matrix(NA, 3, 0)))
	
}


max_holder <- function(e, eps, c = exp(1 + 2 * eps)) {
 	
 	n <- length(e)
	
	eps.cum <- c(0, cumsum(e))
		
	max.stat <- 0
	
	for (i in 0:(n-1)) for (j in (i+1):n) {
		
		scan.stat <- abs(eps.cum[j+1] - eps.cum[i+1]) / sqrt(j-i) / log(c * n / (j-i))^(1/2+eps)

		if (scan.stat > max.stat) max.stat <- scan.stat
				
	}

	max.stat
	
}


est_var <- function(y, x, power = 1/2, min.size = 20, estVn2 = FALSE) {
	
	n <- length(y)
	w.size <- min(n, max(round(n^power), min.size))

	how.many <- n - w.size + 1
	
	res <- rep(0, how.many)
	
	for (i in 1:how.many) {
		
		resp <- y[i:(i+w.size-1)]
		covs <- x[i:(i+w.size-1),]
		
		res[i] <- summary(stats::lm(resp ~ covs))$sigma
		
	}	

	if (estVn2) est <- n / (n - w.size + 1) * sum(res^2)
	else est <- stats::median(res)
	
	est
	
}
