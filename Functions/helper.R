library(mvtnorm)
library(Matrix)
library(mda)
library(matrixcalc)
library(psych)
library(MASS)
library(pracma)
library(msda)
library(sparseLDA)
library(Rcpp)

source("Marginal.R")
sourceCpp("alpha.cpp")

# -----------------------------------
# Two covariance matrix structures
# -----------------------------------
covCS <- function(p, rho) {
    a <- matrix(rho, ncol = p, nrow = p)
    b <- diag(rep(1-rho, p))
    return (a + b)
}

covAR <- function(p, rho) {
    times <- 1:p
    sigma <- 1
    H <- abs(outer(times, times, "-"))
    V <- sigma * rho^H
    p <- nrow(V)
    V[cbind(1:p, 1:p)] <- V[cbind(1:p, 1:p)] * sigma
    return (V)
}

# ------- Candidate functions ---------
fcn <- function(yy, beta, all.linear = TRUE) {
    res <- c(beta %*% yy)
    if (!all.linear){
        res <- sin(pi / 6 * res)
    }
    return(res) 
}

# -----------------------------------
# Discrete Kernel Functions
# -----------------------------------
ker_1 <- function(J_1, J_2, classes) {
    if (sum(J_1 == J_2) == length(classes)) {
        res <- sum((J_1 == J_2) * sqrt(classes)) + 1e-2
    } else {
        res <- sum((J_1 == J_2) * sqrt(classes))
    }
    return(res)
}

ker_2 <- function(J_1, J_2, classes) {
    # --------------------
    # Eskin Measure
    # --------------------
    res <- 0
    for (j in 1:length(J_1)) {
        res <- res + classes[j] ^ 2 / (classes[j] ^ 2 + 2) +
            (J_1[j] == J_2[j]) * 2 / (classes[j] ^ 2 + 2)
    }
    return(res)
}

ker_3 <- function(J_1, J_2, classes, xtrain) {
    # --------------------
    # IOF
    # --------------------
    freq <- list()
    for (k in 1:length(classes)) {
        freq[[k]] <- rep(0, classes[k])
        for (j in 1:classes[k]) {
            freq[[k]][j] <- sum(xtrain[, k] == j - 1)
        }
    }
    res <- 0
    for (j in 1:length(J_1)) {
        res <- res + 1 / (1 + log(freq[[j]][J_1[j] + 1]) * log(freq[[j]][J_2[j] + 1])) +(J_1[j] == J_2[j]) * (log(freq[[j]][J_1[j] + 1]) * log(freq[[j]][J_2[j] + 1])) / (1 + log(freq[[j]][J_1[j] + 1]) * log(freq[[j]][J_2[j] + 1])) 
    }
    if (sum(J_1 == J_2) == length(classes)) {
        res <- res + 1e-2
    }
    return(res)
}

Kx.mat <- function(xtrain, classes, kernel) {
    # ------------------------------------------------
    # Calculate K matrix, with a factor of 1/n.
    # ------------------------------------------------
    ntrain <- dim(xtrain)[1]
    K <- matrix(0, nrow = ntrain, ncol = ntrain)
    if (kernel == 1) {
        for (j in 1 : ntrain) {
            for (l in 1 : ntrain) {
                K[j, l] <- ker_1(xtrain[j, ], xtrain[l, ], classes)
            }
        }
    } else if (kernel == 2) {
        for (j in 1 : ntrain) {
            for (l in 1 : ntrain) {
                K[j, l] <- ker_2(xtrain[j, ], xtrain[l, ])
            }
        }
    } else if (kernel == 3) {
        for (j in 1 : ntrain) {
            for (l in 1 : ntrain) {
                K[j, l] <- ker_3(xtrain[j, ], xtrain[l, ], classes)
            }
        }
    } 
    return(K / ntrain)
}

Kx.vec <- function(x, xtrain, classes, kernel) {
    # ------------------------------------------------------------------
    # For a particular J, calculate vector (K(J_1,J),...,K(J_n,J)).
    # ------------------------------------------------------------------
    res <- rep(0, dim(xtrain)[1])
    if (kernel == 1) {
        for (i in 1 : dim(xtrain)[1]) {
            res[i] <- ker_1(x, xtrain[i, ], classes)
        }
    } else if (kernel == 2) {
        for (i in 1 : dim(xtrain)[1]) {
            res[i] <- ker_2(x, xtrain[i, ])
        }
    } else if (kernel == 3) {
        for (i in 1 : dim(xtrain)[1]) {
            res[i] <- ker_3(x, xtrain[i, ], classes)
        }
    }
    return(res)
}


# -------------------------
# Auxiliary Functions
# -------------------------
findIndex <- function(x, classes) {
    # ------------------------------------------------
    # Given a J_i, find its index in expand.grid(...)
    # ------------------------------------------------
    dec <- cumprod(classes)
    res <- x[1]
    for (j in 2 : length(x)) {
        res <- res + x[j] * dec[j - 1]
    }
    return(res + 1)
}

getQ <- function(x, classes, kernel = 1) {
    n <- dim(x)[1]
    ntmp <- prod(classes)
    Q <- matrix(0, nrow = n, ncol = ntmp)
    Xtmp <- matrix(NA, nrow = prod(classes), ncol = length(classes))
    
    for(i in 1:n) {
        ind <- findIndex(x[i, ], classes)
        Q[i, ind] <- 1
        Xtmp[ind, ] <- x[i, ]
    }
    
    Xtmp <- na.omit(Xtmp)
    K <- Kx.mat(Xtmp, classes, kernel)
    
    if (any(colSums(Q) == 0)) {
        Q <- Q[,-which(colSums(Q) == 0)]
    }
    return(list("Q" = Q, "K" = K, "Xtilde" = Xtmp))
}


# ------------------------------------------
# Helper functions for Nonconvex Estimator
# ------------------------------------------
update.Omega.nonconvex <- function(y, alpha, K, Q, eta, diag.cov = FALSE) {
    n <- dim(y)[1]
    p <- dim(y)[2]
    ntilde <- dim(alpha)[1]
    if (!diag.cov){
        if (eta == 0) {
            return (chol2inv(chol(1 / n * crossprod(y - sqrt(ntilde) * Q %*% K %*% alpha))))
        } else {
            S <- 1 / n * crossprod(y - sqrt(ntilde) * Q %*% K %*% alpha)
            decomp <- eigen(S)
            D <- diag(decomp$values)
            V <- decomp$vectors
            Q <- 0.5 / eta * V %*% tcrossprod(-D + (D ^ 2 + 4 * eta * diag(p)) ^ 0.5, V)
            return(Q)
        }
    } else {
        sampleCov <- 1 / n * crossprod(y - sqrt(ntilde) * Q %*% K %*% alpha)
        return(diag(diag(sampleCov) ^ -1))
    }
}


# ------------------------------------------
# Helper functions for Convex Estimator
# ------------------------------------------
eval.obj.convex <- function(y, Beta, omega, K, Q, gamma, eta) { # FIXME: 
    # ------------------------------------------------
    # Evaluation of objective function
    # ------------------------------------------------
    n <- dim(y)[1]
    p <- dim(y)[2]
    ntilde <- dim(K)[1]
    
    res <- 0
    tt <- y %*% omega - sqrt(ntilde) * Q %*% K %*% Beta
    res <- res + 1 / n * tr(tt %*% solve(omega) %*% t(tt))
    # res <- res - log(det(omega))
    res <- res - determinant(omega, logarithm = TRUE)$modulus[1]
    res <- res + gamma * sum(sqrt(colSums(Beta^2)))
    res <- res + eta * sum(omega ^ 2) / 2
    return(as.numeric(res))
}

eval.f <- function(y, Beta, omega, K, Q) { 
    n <- dim(y)[1]
    ntilde <- dim(K)[1]
    res <- -2 * sqrt(ntilde) / n * tr(Q %*% K %*% Beta %*% t(y)) + 
        ntilde / n * tr(Q %*% K %*% Beta %*% tcrossprod(solve(omega), Beta) %*% tcrossprod(K, Q))
    return(res)
}

eval.f.hat <- function(y, Beta, Beta0, omega, K, Q, stepsize) { 
    res <- eval.f(y, Beta0, omega, K, Q) + 
        tr(crossprod(get_grad_Beta(y, Beta0, K, Q, omega), (Beta - Beta0))) +
        1 / (2 * stepsize) * sum((Beta - Beta0) ^ 2)
    return(res)
}

get_grad_Beta <- function(y, Beta, K, Q, omega) { 
    n <- dim(y)[1]
    ntilde <- dim(K)[1]
    res <- - 2 * sqrt(ntilde) / n * tcrossprod(K, Q) %*% y + 
        2 * ntilde / n * K %*% crossprod(Q) %*% K %*% Beta %*% solve(omega)
    return(res)
}

prox <- function(v, lambda) {  
    v_norm <- sqrt(sum(v ^ 2))
    if (v_norm >= lambda) {
        return((1 - lambda / v_norm) * v)
    } else {
        return(0)
    }
}

update.Beta <- function(y, K, Q, omega, gamma, eta, stepsize = 2, BetaPre = NULL) {
    p <- dim(y)[2]
    ntilde <- dim(K)[1]
    if (!is.null(BetaPre)) {
        Beta <- BetaPre
    } else {
        Beta <- matrix(0, nrow = ntilde, ncol = p)
        BetaPre <- matrix(0, nrow = ntilde, ncol = p)
    }

    scaling <- 0.5
    objPre <- Inf
    obj <- eval.obj.convex(y, Beta, omega, K, Q, gamma, eta)
    num.iter <- 0
    
    # ------------------------------------------
    # Accelerated proximal gradient method
    # ------------------------------------------
    while (abs(objPre - obj) > 1e-4 & num.iter <= 1000) {
        Beta0 <- Beta + num.iter / (num.iter + 3) * (Beta - BetaPre)
        objPre <- obj
        grad.f <- get_grad_Beta(y, Beta0, K, Q, omega)
        while (TRUE) {
            z <- matrix(0, nrow = ntilde, ncol = p)
            for (j in 1:p) {
                z[, j] <- prox(Beta0[, j] - stepsize * grad.f[, j], stepsize * gamma)
            }
            if (eval.f(y, z, omega, K, Q) <= eval.f.hat(y, z, Beta0, omega, K, Q, stepsize)) {
                break
            }
            stepsize <- stepsize * scaling
        }
        
        BetaPre <- Beta
        Beta <- z
        obj <- eval.obj.convex(y, Beta, omega, K, Q, gamma, eta)
        num.iter <- num.iter + 1
    }
    return(list('Beta' = Beta, 'stepsize' = stepsize))
}

get_grad_Omega <- function(y, Beta, omega, K, Q, eta) {
    n <- dim(y)[1]
    p <- dim(y)[2]
    ntilde <- dim(Beta)[1]
    H <- sqrt(ntilde) * Q %*% K %*% Beta

    res <- matrix(0, nrow = p, ncol = p)
    omega_inv <- solve(omega)
    res <- crossprod(y) / n - omega_inv %*% crossprod(H) %*% omega_inv / n - omega_inv + eta * omega
    return(res)
}

conv.proj <- function(X){
    eig <- eigen(X)
    dv <- eig$values
    dv[dv < 1e-2] <- 1e-2
    dv[dv > 1e2] <- 1e2
    return (eig$vectors %*% diag(dv) %*% t(eig$vectors))
}

update.Omega.convex <- function(y, Beta, omega, K, Q, gamma, eta, vareps = 1e-2) {
    n <- dim(y)[1]
    p <- dim(y)[2]
    ntilde <- dim(K)[1]

    obj.orig <- eval.obj.convex(y, Beta, omega, K, Q, gamma, eta)
    ss <- 1
    num.iter <- 0

    while (num.iter <= 1000) {
        tmp_grad <- get_grad_Omega(y, Beta, omega, K, Q, eta)
        linesearch <- TRUE

        while (linesearch) {
            search.point <- omega - ss * tmp_grad
            omega.new <- conv.proj(search.point)
            obj.new <- eval.obj.convex(y, Beta, omega.new, K, Q, gamma, eta)
            if (obj.new <= (obj.orig + sum(tmp_grad * (omega.new - omega)) + (1/(2 * ss)) * sum((omega - omega.new) ^ 2))) {
                linesearch <- FALSE
                obj.orig <- obj.new
            } else {
                ss <- ss / 2
            }
        }
        if (max(abs(omega - omega.new)) > 1e-6) {
            omega <- omega.new
        } else {
            break
        }

        num.iter <- num.iter + 1
    }
    return(omega)
}



