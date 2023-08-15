# ------------------------------
# EM algorithm
# ------------------------------
log_likelihood <- function(delta, psi, W, classes) {
    n <- dim(W[[1]])[1]
    R <- length(delta)
    M <- length(classes)
    res <- 0
    for (i in 1:n) {
        tmp1 <- 0
        for (k in 1:R) {
            tmp2 <- delta[k]
            for (l in 1:M) {
                for (j in 1:classes[l]) {
                    tmp2 <- tmp2 * psi[[k]][[l]][j] ^ W[[l]][i, j]
                }
            }
            tmp1 <- tmp1 + tmp2
        }
        res <- res + log(tmp1)
    }
    return(res)
}


Marg.prob <- function(Y, r.rank, classes) {
    n <- dim(Y)[1]
    M <- dim(Y)[2]
    
    W <- list()
    for (l in 1:M) {
        W[[l]] <- matrix(0, nrow = n, ncol = classes[l])
    }
    for (i in 1:n) {
        for (l in 1:M) {
            W[[l]][i, Y[i, l] + 1] <- 1
        }
    }

    # -----------------
    # Initialization
    # -----------------
    delta <- rep(1, r.rank)
    delta <- delta / length(delta)
    psi <- list()
    for (k in 1 : r.rank) {
        psi[[k]] <- list()
        for (l in 1 : M) {
            psi[[k]][[l]] <- runif(classes[l])
            psi[[k]][[l]] <- psi[[k]][[l]] / sum(psi[[k]][[l]])
        }
    }
    q <- matrix(1 / r.rank, nrow = n, ncol = r.rank)

    objPre <- Inf
    obj <- log_likelihood(delta, psi, W, classes)
    cat("--------------------\n")
    cat("Initial: ", obj, "\n")

    while (abs(objPre - obj) > 1e-6) {
        objPre <- obj
        # ----------------
        # E-step
        # ----------------
        for (i in 1:n) {
            for (k in 1:r.rank) {
                num <- delta[k]
                for (l in 1:M) {
                    num <- num * psi[[k]][[l]][Y[i, l] + 1]
                }
                
                den <- 0
                for (s in 1:r.rank) {
                    tmp <- delta[s]
                    for (l in 1:M) {
                        tmp <- tmp * psi[[s]][[l]][Y[i, l] + 1]
                    }
                    den <- den + tmp
                }
                
                q[i, k] <- num / den
            }
        }
        
        # ----------------
        # Q-step
        # ----------------
        for (k in 1:r.rank) {
            delta[k] <- 1 / n * sum(q[, k])
            
            for (l in 1:M) {
                for (j in 1:classes[l]) {
                    num <- 0
                    for (i in 1:n) {
                        num <- num + q[i, k] * W[[l]][i, j]
                    }
                    den <- sum(q[, k])
                    psi[[k]][[l]][j] <- num / den
                }
            }
        }
        obj <- log_likelihood(delta, psi, W, classes)
    }
    cat("Final Eval: ", obj, "\n")


    y_all <- as.matrix(do.call(expand.grid, lapply(classes, function(x) 0:(x-1))))
    est.prob <- rep(0, dim(y_all)[1])
    for (j in 1:(dim(y_all)[1]) ){
        for (k in 1:r.rank) {
            tmp <- delta[k]
            for (l in 1:M) {
                tmp <- tmp * psi[[k]][[l]][y_all[j, l] + 1]
            }
            est.prob[j] <- est.prob[j] + tmp
        }
    }
    return(est.prob)
}