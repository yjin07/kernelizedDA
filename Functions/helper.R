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

source("/blue/amolstad/y.jin/KernelLDA/Functions/Marginal.R")
sourceCpp("/blue/amolstad/y.jin/KernelLDA/Functions/alpha.cpp")

# -----Two covariance Structure---------
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

# --------Kernel Functions--------------------------
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
# --------------------------------------------------

# --------Auxiliary Functions-----------------------
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
# --------------------------------------------------


# ----------------Cross Validation------------------
cv.kernel.nonconvex <- function(x, y, xval, yval, x_all, classes, eta.vec, max.rank = 4, 
                ngamma = 35, delta = 0.001, kernel = 1, max.iter = 100, diag.cov = FALSE) 
{
    n <- dim(x)[1]
    p <- dim(y)[2]
    M <- length(classes)

    # ------------------
    #   Centralization
    # ------------------
    col_means <- colMeans(y)
    y <- t(t(y) - col_means) # FIXME:
    # yval <-  t(t(yval) - col_means)

    # ------------------
    # Find seq of gamma
    # ------------------
    getElement <- getQ(x, classes, kernel = kernel)
    K <- getElement$K
    Q0 <- getElement$Q
    xtilde <- getElement$Xtilde
    ntilde <- dim(xtilde)[1]
    Q <- (diag(n) - matrix(1 / n, nrow = n, ncol = n)) %*% Q0 # FIXME:

    mat <- K %*% crossprod(Q, y) %*% ginv(crossprod(y))
    vec.val <- sqrt(colSums(mat ^ 2))

    gamma.max <- max(vec.val) * 2 * sqrt(ntilde) * 2
    gamma.min <- delta * gamma.max
    gamma.vec <- 10^seq(log10(gamma.max), log10(gamma.min), length=ngamma)
    neta <- length(eta.vec)
    # --ENDS HERE------
    
    # ------------------------
    # Initialization of Omega
    # ------------------------
    SampleMean <- matrix(0, nrow = dim(x_all)[1], ncol = p)
    counts <- rep(0, dim(x_all)[1])
    for (ii in 1:n) {
        index <- findIndex(x[ii,], classes)
        counts[index] <- counts[index] + 1
        SampleMean[index, ] <- SampleMean[index, ] + y[ii, ]
    }

    for (jj in 1 : dim(x_all)[1]) {
        SampleMean[jj, ] <- SampleMean[jj, ] / counts[jj]
    }

    freq0 <- counts / n

    # --- Calculate Sample Covariance 1 ---
    Tmp <- na.omit(cbind(x_all, SampleMean))
    SampleCov <- matrix(0, ncol = p, nrow = p)
    ind <- rep(0, dim(Tmp)[1])
    for (j in 1:dim(Tmp)[1]) {
     ind[j] <- findIndex(Tmp[j, 1:M], classes)
    }

    for (i in 1 : n) {
        index <- which(ind == findIndex(x[i,], classes))
        SampleCov <- SampleCov + tcrossprod(y[i, ] - Tmp[index, (M+1):(M+p)])
    }

    SampleCov <- SampleCov / n
    # ----------------ENDS HERE------------------

    # ----------Marginal Probability-------------
    prob.list  <- list()
    for (rr in 1 : max.rank) {
        prob.list[[rr]] <- Marg.prob(x, rr, classes)
        cat("Rank", rr, "decomp finished.\n")
    }
    # ----------------ENDS HERE------------------
    
    total <- prod(classes)
    y_all <- x_all
    cvErr1 <- array(0, dim = c(neta, ngamma, max.rank))
    cvErr2 <- array(0, dim = c(neta, ngamma, max.rank))
    cvErr3 <- array(0, dim = c(neta, ngamma, max.rank))

    # ---------------------------
    # Train-validation split
    # ---------------------------
    xtest <- xval
    ytest <- yval
    xtrain <- x
    ytrain <- y
    ntrain <- dim(xtrain)[1]
    ntest <- dim(xtest)[1]
    alpha.list <- list()
    Omega.list <- list()

    for (ieta in 1 : neta) {
        eta <- eta.vec[ieta]
        cat("\n\n\n\n=======================================\n")
        cat(ieta,"-th eta:, ", eta, "\n", sep = "")
        cat("=======================================\n")

        alpha.list[[ieta]] <- list()
        Omega.list[[ieta]] <- list()

        if (!diag.cov) {
            if (eta == 0){
                Omega.initial <- chol2inv(chol(SampleCov))
            } else {
                decomp <- eigen(SampleCov)
                D <- diag(decomp$values)
                V <- decomp$vectors
                Omega.initial <- 0.5 / eta * V %*% tcrossprod(-D + (D ^ 2 + 4 * eta * diag(p)) ^ 0.5, V)
            }
        } else {
            Omega.initial <- diag(diag(SampleCov) ^ -1)
        }

        for (j in 1 : ngamma) {
            cat(j,"-th gamma:, ", gamma.vec[j], "\n", sep = "")
            gamma <- gamma.vec[j]
    
            Omega <- Omega.initial
            llpre <- Inf
            if (j == 1) {
                starting <- update_alpha(ytrain, K, Q, Omega, gamma, eta)
            } else {
                starting <- update_alpha(ytrain, K, Q, Omega, gamma, eta, alphaPre = alpha)
            }
            
            alpha <- starting$alpha
            stepsize <- starting$stepsize
            llcur <- eval_obj(ytrain, alpha, Omega, K, Q, gamma, eta)
            iter <- 1
            
            # ------------------------
            # Method 1
            # ------------------------
            while (abs(llcur - llpre) > 1e-6 & iter <= max.iter) {
                # Record
                llpre <- llcur

                # Update
                Omega <- update.Omega.nonconvex(ytrain, alpha, K, Q, eta, diag.cov = diag.cov)
                tmp <- update_alpha(ytrain, K, Q, Omega, gamma, eta, stepsize, alphaPre = alpha)
                alpha <- tmp$alpha
                stepsize <- tmp$stepsize
                llcur <-  eval_obj(ytrain, alpha, Omega, K, Q, gamma, eta)
                iter <- iter + 1
                # cat("+++Iter: ", iter, "Eval: ", llcur, "\n")
            }
            cat("Number of iterations: ", iter, "\n")
            # cat("========================================\n")
            
            Omega.list[[ieta]][[j]] <- Omega
            alpha.list[[ieta]][[j]] <- alpha
            
            est.cov <- chol2inv(chol((Omega + t(Omega)) / 2))
            est.cov <- (est.cov + t(est.cov)) / 2
            est.mean <- matrix(0, nrow = dim(y_all)[1], ncol = p)

            for (ii in 1 : total) {
                est.mean[ii, ] <-  Kx.vec(y_all[ii,], xtilde, classes, kernel) %*% alpha / sqrt(ntilde) + col_means - colMeans(Q0 %*% K %*% alpha) * sqrt(ntilde) 
            }


            Pmat <- matrix(0, nrow = ntest, ncol = prod(classes))
            prm0 <- Sys.time()
            for (kk in 1 : max.rank) {
                log_lld <- 0
                for (jj in 1 : total) {
                    ycent <- t(t(ytest) - est.mean[jj, ])
                    Pmat[, jj] <-  - diag(ycent %*% tcrossprod(Omega, ycent)) + 2 * log(prob.list[[kk]][jj])
                }

                for (ii in 1:ntest) {
                    ind.true <- findIndex(xtest[ii,], classes)
                    log_lld <- log_lld + log(dmvnorm(ytest[ii, ], mean = est.mean[ind.true, ], sigma = est.cov)) + log(prob.list[[kk]][ind.true])
                }

                pred.classes <- apply(Pmat, 1, which.max)
                x.pred <- matrix(0, ncol = M, nrow = ntest)
                x.pred <- y_all[pred.classes, ]
                cvErr1[ieta, j, kk] <- sum(rowSums(x.pred == xtest) == M) / ntest
                cvErr2[ieta, j, kk] <- mean(colSums(x.pred == xtest) / ntest)
                cvErr3[ieta, j, kk] <- log_lld
                cat("----->Rank of ", kk, " : ", round(cvErr1[ieta, j, kk], 3), round(cvErr2[ieta, j, kk], 3),round(cvErr3[ieta, j, kk], 3),"------\n", sep = " ")
            }
        }

    }

    # -----------
    # Accuracy
    # -----------
    inds <- which(cvErr1 == max(cvErr1), arr.ind = TRUE)
    sparsity <- rep(0, dim(inds)[1])
    for(lll in 1:dim(inds)[1]){
        sparsity[lll] <- sum(colSums(alpha.list[[inds[lll,1]]][[inds[lll,2]]]!=0))
    }
    best.inds <- inds[which(sparsity == min(sparsity))[1],]
    t1_acc <- best.inds[1]
    t2_acc <- best.inds[2]
    t3_acc <- best.inds[3]
    eta.tuned.acc <- eta.vec[t1_acc]
    gamma.tuned.acc <- gamma.vec[t2_acc]
    rank.tuned.acc <- t3_acc

    # ------------------
    # Hamming Distance
    # ------------------
    inds <- which(cvErr2 == max(cvErr2), arr.ind = TRUE)
    sparsity <- rep(0, dim(inds)[1])
    for(lll in 1:dim(inds)[1]){
        sparsity[lll] <- sum(colSums(alpha.list[[inds[lll,1]]][[inds[lll,2]]]!=0))
    }
    best.inds <- inds[which(sparsity == min(sparsity))[1],]
    t1_ham <- best.inds[1]
    t2_ham <- best.inds[2]
    t3_ham <- best.inds[3]
    eta.tuned.ham <- eta.vec[t1_ham]
    gamma.tuned.ham <- gamma.vec[t2_ham]
    rank.tuned.ham <- t3_ham
    
    # Validation
    inds <- which(cvErr3 == max(cvErr3), arr.ind = TRUE)
    sparsity <- rep(0, dim(inds)[1])
    for(lll in 1:dim(inds)[1]){
        sparsity[lll] <- sum(colSums(alpha.list[[inds[lll,1]]][[inds[lll,2]]]!=0))
    }
    best.inds <- inds[which(sparsity == min(sparsity))[1],]
    t1_lld <- best.inds[1]
    t2_lld <- best.inds[2]
    t3_lld <- best.inds[3]
    eta.tuned.lld <- eta.vec[t1_lld]
    gamma.tuned.lld <- gamma.vec[t2_lld]
    rank.tuned.lld <- t3_lld

    cat("\n\nBest Acc on Val: ", max(cvErr1), "\n")
    cat("Eta (Acc): ", eta.tuned.acc, "\n")
    cat("Gamma (Acc): ", gamma.tuned.acc, "\n")
    cat("Rank (Acc): ", rank.tuned.acc, "\n")
    cat("Eta (hamming): ", eta.tuned.ham, "\n")
    cat("Gamma (hamming): ", gamma.tuned.ham, "\n")
    cat("Rank (hamming): ", rank.tuned.ham, "\n")
    cat("Eta (lld): ", eta.tuned.lld, "\n")
    cat("Gamma (lld): ", gamma.tuned.lld, "\n")
    cat("Rank (lld): ", rank.tuned.lld, "\n")

    tuned <- rbind(c(eta.tuned.acc, gamma.tuned.acc,rank.tuned.acc),
                  c(eta.tuned.ham, gamma.tuned.ham,rank.tuned.ham),
                  c(eta.tuned.lld, gamma.tuned.lld,rank.tuned.lld))
    
    return (list("tune" = tuned, "alpha" = alpha.list[[t1_acc]][[t2_acc]], "Omega" = Omega.list[[t1_acc]][[t2_acc]], "xtilde" = xtilde, "freq" = prob.list[[t3_acc]], "sample.freq" = freq0, "sample.cov" = SampleCov, "Q0" = Q0, "K" = K,
    "alpha.lld" = alpha.list[[t1_lld]][[t2_lld]], "Omega.lld" = Omega.list[[t1_lld]][[t2_lld]], "freq.lld" = prob.list[[t3_lld]],
    "alpha.ham" = alpha.list[[t1_ham]][[t2_ham]], "Omega.ham" = Omega.list[[t1_ham]][[t2_ham]], "freq.ham" = prob.list[[t3_ham]]
    ))
}



# -----------------------------------
#          Convex Estimator
# -----------------------------------
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
    # BetaPre <- matrix(0, nrow = ntilde, ncol = p)
    obj <- eval.obj.convex(y, Beta, omega, K, Q, gamma, eta)
    num.iter <- 0
    
    # Accelerated proximal gradient method
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
        # cat("Iter: ", num.iter, "eval: ", obj, "\n")
    }
    # cat("\nNumber of iterations when updating Beta: ", num.iter, "\n")
    # cat("Final eval: ", obj, "\n")
    # cat("Time costs: ", (proc.time() - prm)[1], "\n")
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


update_Omega <- function(y, Beta, omega, K, Q, gamma, eta, vareps = 1e-2) {
    n <- dim(y)[1]
    p <- dim(y)[2]
    ntilde <- dim(K)[1]

    obj.orig <- eval.obj.convex(y, Beta, omega, K, Q, gamma, eta)
    # cat("\nOriginal obj: ", obj.orig, "\n")
    ss <- 1
    num.iter <- 0
    # omega.Pre <- diag(p)

    while (num.iter <= 1000) {
        tmp_grad <- get_grad_Omega(y, Beta, omega, K, Q, eta)
        linesearch <- TRUE

        while (linesearch) {
            search.point <- omega - ss * tmp_grad
            omega.new <- conv.proj(search.point)
            obj.new <- eval.obj.convex(y, Beta, omega.new, K, Q, gamma, eta)
            # cat("\nstepsize:",ss, ":", obj.new, "\n")
            if (obj.new <= (obj.orig + sum(tmp_grad * (omega.new - omega)) + (1/(2 * ss)) * sum((omega - omega.new) ^ 2))) {
                linesearch <- FALSE
                obj.orig <- obj.new
            } else {
                ss <- ss / 2
            }
        }

        # cat("Iteration at", num.iter, "eval: ", obj.orig, "\n")
        if (max(abs(omega - omega.new)) > 1e-6) {
            omega <- omega.new
        } else {
            break
        }

        num.iter <- num.iter + 1
    }

    # cat("\nNumber of iterations when updating Omega: ", num.iter, "\n")
    return(omega)
}
# --------------------------------------------------


# ----------------Cross Validation------------------
cv.kernel.convex <- function(x, y, xval, yval, x_all, classes, eta.vec, max.rank = 4, ngamma = 35, 
                eta = 0, delta = 0.001, kernel = 1, max.iter = 100, diag.cov = FALSE) {
    n <- dim(x)[1]
    p <- dim(y)[2]
    M <- length(classes)

    # ------------------
    #   Centralization
    # ------------------
    col_means <- colMeans(y)
    y <- t(t(y) - col_means)

    # ------------------
    # Find seq of gamma
    # ------------------
    getElement <- getQ(x, classes, kernel = kernel)
    K <- getElement$K
    Q0 <- getElement$Q
    xtilde <- getElement$Xtilde
    ntilde <- dim(xtilde)[1]
    Q <- (diag(n) - matrix(1 / n, nrow = n, ncol = n)) %*% Q0 # FIXME:

    mat <- K %*% crossprod(Q, y)
    vec.val <- sqrt(colSums(mat ^ 2))

    gamma.max <- max(vec.val) * 2 * sqrt(ntilde) * 2
    gamma.min <- delta * gamma.max
    gamma.vec <- 10^seq(log10(gamma.max), log10(gamma.min), length=ngamma)
    neta <- length(eta.vec)
    # --ENDS HERE------
    
    # ------------------------
    # Initialization of Omega
    # ------------------------
    SampleMean <- matrix(0, nrow = dim(x_all)[1], ncol = p)
    counts <- rep(0, dim(x_all)[1])
    for (ii in 1:n) {
        index <- findIndex(x[ii,], classes)
        counts[index] <- counts[index] + 1
        SampleMean[index, ] <- SampleMean[index, ] + y[ii, ]
    }

    for (jj in 1 : dim(x_all)[1]) {
        SampleMean[jj, ] <- SampleMean[jj, ] / counts[jj]
    }

    freq0 <- counts / n

    Tmp <- na.omit(cbind(x_all, SampleMean))
    SampleCov <- matrix(0, ncol = p, nrow = p)
    ind <- rep(0, dim(Tmp)[1])
    for (j in 1:dim(Tmp)[1]) {
     ind[j] <- findIndex(Tmp[j, 1:M], classes)
    }

    for (i in 1 : n) {
        index <- which(ind == findIndex(x[i,], classes))
        SampleCov <- SampleCov + tcrossprod(y[i, ] - Tmp[index, (M+1):(M+p)])
    }

    SampleCov <- SampleCov / n

    # if (!diag.cov) {
    #     if (eta == 0){
    #         Omega.initial <- chol2inv(chol(SampleCov))
    #     } else {
    #         decomp <- eigen(SampleCov)
    #         D <- diag(decomp$values)
    #         V <- decomp$vectors
    #         Omega.initial <- 0.5 / eta * V %*% tcrossprod(-D + (D ^ 2 + 4 * eta * diag(p)) ^ 0.5, V)
            
    #     }
    # } else {
    #     Omega.initial <- diag(diag(SampleCov) ^ -1)
    # }
    # --ENDS HERE------

    # ----------Marginal Probability-------------
    prm1 <- proc.time()
    prob.list  <- list()
    for (rr in 1 : max.rank) {
        prob.list[[rr]] <- Marg.prob(x, rr, classes)
        cat("Rank", rr, "decomp finished.\n")
    }
    prob.list[[max.rank + 1]] <- freq0
    cat("------------------------------------------\n")
    cat("Computation time: ", proc.time() - prm1, "\n")
    cat("------------------------------------------\n")
    # ----------------ENDS HERE------------------
    
    total <- prod(classes)
    y_all <- x_all
    cvErr1 <- array(0, dim = c(neta, ngamma, max.rank))
    cvErr2 <- array(0, dim = c(neta, ngamma, max.rank))
    cvErr3 <- array(0, dim = c(neta, ngamma, max.rank))

    # ---------------------------
    # Train-validation split
    # ---------------------------
    xtest <- xval
    ytest <- yval
    xtrain <- x
    ytrain <- y
    ntrain <- dim(xtrain)[1]
    ntest <- dim(xtest)[1]
    # class.err <- array(Inf, dim = c(neta, ngamma, max.rank, ntest))

    Beta.list <- list()
    Omega.list <- list()

    for (ieta in 1 : neta) {
        eta <- eta.vec[ieta]
        cat("\n\n\n\n=======================================\n")
        cat(ieta,"-th eta:, ", eta, "\n", sep = "")
        cat("=======================================\n")

        Beta.list[[ieta]] <- list()
        Omega.list[[ieta]] <- list()
    
        if (!diag.cov) {
            if (eta == 0){
                Omega.initial <- chol2inv(chol(SampleCov))
            } else {
                decomp <- eigen(SampleCov)
                D <- diag(decomp$values)
                V <- decomp$vectors
                Omega.initial <- 0.5 / eta * V %*% tcrossprod(-D + (D ^ 2 + 4 * eta * diag(p)) ^ 0.5, V)
            }
        } else {
            Omega.initial <- diag(diag(SampleCov) ^ -1)
        }
    
        for (j in 1 : ngamma) {
            cat(j,"-th gamma:, ", gamma.vec[j], "\n", sep = "")
            gamma <- gamma.vec[j]
            # -------------------------------------------------------
            # Use sample covariance for initialization of Omega
            # -------------------------------------------------------
            Omega <- Omega.initial
            llpre <- Inf
            if (j == 1) {
                starting <- update.Beta(ytrain, K, Q, Omega, gamma, eta)
            } else {
                starting <- update.Beta(ytrain, K, Q, Omega, gamma, eta, BetaPre = Beta)
            }
            
            Beta <- starting$Beta
            stepsize <- starting$stepsize
            llcur <- eval.obj.convex(ytrain, Beta, Omega, K, Q, gamma, eta)
            iter <- 1
            
            # ------------------------
            # Method 1
            # ------------------------
            while (abs(llcur - llpre) > 1e-6 & iter <= max.iter) {
                # Record
                llpre <- llcur

                # Update
                Omega <- update_Omega(ytrain, Beta, Omega, K, Q, gamma, eta)
                tmp <- update.Beta(ytrain, K, Q, Omega, gamma, eta, stepsize, BetaPre = Beta)
                Beta <- tmp$Beta
                stepsize <- tmp$stepsize
                llcur <-  eval.obj.convex(ytrain, Beta, Omega, K, Q, gamma, eta)
                iter <- iter + 1
                # cat("Iter: ", iter, "Eval: ", llcur, "\n")
            }

            Omega.list[[ieta]][[j]] <- Omega
            Beta.list[[ieta]][[j]] <- Beta
            
            est.cov <- chol2inv(chol((Omega + t(Omega)) / 2))
            est.cov <- (est.cov + t(est.cov)) / 2
            est.mean <- matrix(0, nrow = dim(y_all)[1], ncol = p)

            for (ii in 1 : total) {
                est.mean[ii, ] <-  Kx.vec(y_all[ii,], xtilde, classes, kernel)%*% Beta %*% est.cov / sqrt(ntilde) + col_means - colMeans(Q0 %*% K %*% Beta %*% est.cov) * sqrt(ntilde)
            }

            # counts <- counts / ntrain
            Pmat <- matrix(0, nrow = ntest, ncol = prod(classes))
            # prm0 <- Sys.time()
            for (kk in 1 : max.rank) {
                log_lld <- 0
                for (jj in 1 : total) {
                    ycent <- t(t(ytest) - est.mean[jj, ])
                    Pmat[, jj] <-  - diag(ycent %*% tcrossprod(Omega, ycent)) + 2 * log(prob.list[[kk]][jj])
                }

                for (ii in 1:ntest) {
                    ind.true <- findIndex(xtest[ii,], classes)
                    log_lld <- log_lld + log(dmvnorm(ytest[ii, ], mean = est.mean[ind.true, ], sigma = est.cov)) + log(prob.list[[kk]][ind.true])
                }

                pred.classes <- apply(Pmat, 1, which.max)
                x.pred <- matrix(0, ncol = M, nrow = ntest)
                x.pred <- y_all[pred.classes, ]
                # Accuracy
                cvErr1[ieta, j, kk] <- sum(rowSums(x.pred == xtest) == M) / ntest
                # # Accuracy with 1-SE rule
                # class.err[ieta, j, kk, ] <- rowSums(x.pred == xtest) == M
                # Hamming Distance
                cvErr2[ieta, j, kk] <- mean(colSums(x.pred == xtest) / ntest)
                # Validation
                cvErr3[ieta, j, kk] <- log_lld
                cat("----->Rank of ", kk, " : ", round(cvErr1[ieta, j, kk], 3), round(cvErr2[ieta, j, kk], 3),round(cvErr3[ieta, j, kk], 3),"------\n", sep = " ")
            }
            # cat("\nTime cost: ", Sys.time() - prm0, "\n")
        }
    }

    # -----------
    # Accuracy
    # -----------
    inds <- which(cvErr1 == max(cvErr1), arr.ind = TRUE)
    sparsity <- rep(0, dim(inds)[1])
    for(lll in 1:dim(inds)[1]){
        sparsity[lll] <- sum(colSums(Beta.list[[inds[lll,1]]][[inds[lll,2]]]!=0))
    }
    best.inds <- inds[which(sparsity == min(sparsity))[1],]
    t1_acc <- best.inds[1]
    t2_acc <- best.inds[2]
    t3_acc <- best.inds[3]
    eta.tuned.acc <- eta.vec[t1_acc]
    gamma.tuned.acc <- gamma.vec[t2_acc]
    rank.tuned.acc <- t3_acc

    # ------------------
    # Hamming Distance
    # ------------------
    inds <- which(cvErr2 == max(cvErr2), arr.ind = TRUE)
    sparsity <- rep(0, dim(inds)[1])
    for(lll in 1:dim(inds)[1]){
        sparsity[lll] <- sum(colSums(Beta.list[[inds[lll,1]]][[inds[lll,2]]]!=0))
    }
    best.inds <- inds[which(sparsity == min(sparsity))[1],]
    t1_ham <- best.inds[1]
    t2_ham <- best.inds[2]
    t3_ham <- best.inds[3]
    eta.tuned.ham <- eta.vec[t1_ham]
    gamma.tuned.ham <- gamma.vec[t2_ham]
    rank.tuned.ham <- t3_ham
    
    # Validation
    inds <- which(cvErr3 == max(cvErr3), arr.ind = TRUE)
    sparsity <- rep(0, dim(inds)[1])
    for(lll in 1:dim(inds)[1]){
        sparsity[lll] <- sum(colSums(Beta.list[[inds[lll,1]]][[inds[lll,2]]]!=0))
    }
    best.inds <- inds[which(sparsity == min(sparsity))[1],]
    t1_lld <- best.inds[1]
    t2_lld <- best.inds[2]
    t3_lld <- best.inds[3]
    eta.tuned.lld <- eta.vec[t1_lld]
    gamma.tuned.lld <- gamma.vec[t2_lld]
    rank.tuned.lld <- t3_lld

    cat("\n\nBest Acc on Val: ", max(cvErr1), "\n")
    cat("Eta (Acc): ", eta.tuned.acc, "\n")
    cat("Gamma (Acc): ", gamma.tuned.acc, "\n")
    cat("Rank (Acc): ", rank.tuned.acc, "\n")
    cat("Eta (hamming): ", eta.tuned.ham, "\n")
    cat("Gamma (hamming): ", gamma.tuned.ham, "\n")
    cat("Rank (hamming): ", rank.tuned.ham, "\n")
    cat("Eta (lld): ", eta.tuned.lld, "\n")
    cat("Gamma (lld): ", gamma.tuned.lld, "\n")
    cat("Rank (lld): ", rank.tuned.lld, "\n")
    
    tuned <- rbind(c(eta.tuned.acc, gamma.tuned.acc,rank.tuned.acc),
                  c(eta.tuned.ham, gamma.tuned.ham,rank.tuned.ham),
                  c(eta.tuned.lld, gamma.tuned.lld,rank.tuned.lld))
    
    return (list("tune" = tuned, "Beta" = Beta.list[[t1_acc]][[t2_acc]], "Omega" = Omega.list[[t1_acc]][[t2_acc]], "xtilde" = xtilde, "freq" = prob.list[[t3_acc]], "sample.freq" = freq0, "sample.cov" = SampleCov, "Q0" = Q0, "K" = K,
    "Beta.lld" = Beta.list[[t1_lld]][[t2_lld]], "Omega.lld" = Omega.list[[t1_lld]][[t2_lld]], "freq.lld" = prob.list[[t3_lld]],
    "Beta.ham" = Beta.list[[t1_ham]][[t2_ham]], "Omega.ham" = Omega.list[[t1_ham]][[t2_ham]], "freq.ham" = prob.list[[t3_ham]]
    ))
}


