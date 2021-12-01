# Time to human readable format 
proc_time <- function(e) {
    if (e < 60) {
	sprintf("%.1fs", e)
    } else if (e <= 3600) {
	m <- e / 60
	s <- (m - floor(m)) * 60
	sprintf("%.0fm %.0fs", m, s)
    } else {
	h <- e / 3600
	m <- (h - floor(h)) * 60
	sprintf("%.0fh %.0fm", h, m)
    }
}


# ----------------------------------------
# Comparison to MCMC
# ----------------------------------------
metrics <- c("l2", "l1", "tpr", "fdr", "auc", "coverage.n0", "length.n0", 
	     "coverage.0", "length.0", "elapsed")
for (d in 1:5) {
    for (s in 1:2) {
	for (m in 1:2) {
	    load(file=sprintf("../../RData/comparison/c_%d_%d_%d.RData", d, s, m))
	    x <- get(sprintf("c_%d_%d_%d", d, s, m))
	    x.m <- apply(x[ , metrics], 2, function(x) mean(unlist(x), na.rm=T))
	    x.sd <- apply(x[ , metrics], 2, function(x) sd(unlist(x), na.rm=T))

	    for (m in seq_along(metrics)) {
		if (metrics[m] == "elapsed") {
		    cat(" ", proc_time(x.m[m]),
			paste0(" (", proc_time(x.sd[m]), ")"))
		} else {
		    cat(sprintf(" %.3f (%.3f)", x.m[m], x.sd[m]))
		}
		cat(ifelse(m == length(metrics), "\\\\", "&"))
	    }
	    cat("\n")
	}
    }
}


# ----------------------------------------
# Comparison to other methods
# ----------------------------------------
metrics <- c("l2", "l1", "tpr", "fdr", "auc", "elapsed")
for (d in 1:5) { # dgp
    for (s in 1:6) { # size
	for (m in c(3, 1, 2)) { # method

	    if (d == 5 & any(s %in% c(3,4))) next
	    
	    load(file=sprintf("../../RData/simulations/r_%d_%d_%d.RData", d, s, m))
	    x <- get(sprintf("r_%d_%d_%d", d, s, m))
	    x.m <- apply(x[ , metrics], 2, function(x) mean(x, na.rm=T))
	    x.sd <- apply(x[ , metrics], 2, function(x) sd(x, na.rm=T))

	    for (i in seq_along(metrics)) {
		if (metrics[i] == "elapsed") {
		    x.m <- apply(x[ , metrics], 2, function(x) median(x, na.rm=T))
		    cat(" ", proc_time(x.m[i]),
			paste0(" (", proc_time(x.sd[i]), ")"))
		} else {
		    cat(sprintf(" %.3f (%.3f)", x.m[i], x.sd[i]))
		}
		cat(ifelse(i == length(metrics), "\\\\", "&"))
	    }
	    cat("\n")
	}
    }
    cat("\n")
}

# coverage (only offered by SVB)
metrics <- c("coverage.n0", "length.n0", "coverage.0", "length.0")
for (d in 1:5) { # setting
    for (s in 1:6) { # size

	if (d == 5 & any(s %in% c(3,4))) next

	load(file=sprintf("../../RData/simulations/r_%d_%d_%d.RData", d, s, 3))
	x <- get(sprintf("r_%d_%d_%d", d, s, 3))

	x.m <- apply(x[ , metrics], 2, mean)
	x.sd <- apply(x[ , metrics], 2, sd)

	for (m in seq_along(metrics))
	    cat(sprintf(" %.3f (%.3f)", x.m[m], x.sd[m]),
		ifelse(m == length(metrics), "\\\\", "&"))
	cat("\n")
    }
    cat("\n")
}

