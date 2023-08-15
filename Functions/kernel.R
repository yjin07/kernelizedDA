source("helper.R")

# ----------------------------------------------
#  Nonconvex Estimator
# ----------------------------------------------
# x: n x M matrix of responses for training
# y: n x p matrix of predictors for training
# xval: nval x M matrix of responses for validation
# yval: nval x p matrix of predictors for validation
# classes: a vector of length M, with each component being the number of categories for each response
# ----------------------------------------------

cv.kernel.nonconvex <- function(x, y, xval, yval, 
                                classes, eta.vec, max.rank = 4, 
                                ngamma = 35, delta = 0.001, 
                                kernel = 1, max.iter = 100, 
                                diag.cov = FALSE) 
{
    # -----------------------
    # Preliminaries
    # -----------------------
    n <- dim(x)[1]
    p <- dim(y)[2]
    M <- length(classes)
    neta <- length(eta.vec)
    total <- prod(classes)
    x_all <- as.matrix(do.call(expand.grid, lapply(classes, function(x) 0:(x-1))))

    # -----------------------
    # Centralization
    # -----------------------
    col_means <- colMeans(y)
    y <- t(t(y) - col_means) 

    # -----------------------
    # Set gamma.vec 
    # -----------------------
    getElement <- getQ(x, classes, kernel = kernel)
    K <- getElement$K
    Q0 <- getElement$Q
    xtilde <- getElement$Xtilde
    ntilde <- dim(xtilde)[1]
    Q <- (diag(n) - matrix(1 / n, nrow = n, ncol = n)) %*% Q0
    mat <- K %*% crossprod(Q, y) %*% ginv(crossprod(y))
    vec.val <- sqrt(colSums(mat ^ 2))
    gamma.max <- max(vec.val) * 2 * sqrt(ntilde) * 2
    gamma.min <- delta * gamma.max
    gamma.vec <- 10^seq(log10(gamma.max), log10(gamma.min), length=ngamma)
    
    # ---------------------------------
    # Calculation of sample covariance
    # ---------------------------------
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

    # -------------------------------------------
    # Calcalate the marginal probability of x
    # -------------------------------------------
    prob.list  <- list()
    for (rr in 1 : max.rank) {
        prob.list[[rr]] <- Marg.prob(x, rr, classes)
        cat("Rank", rr, "decomp finished.\n")
    }

    
    # ------------------------------------
    # Arrays and lists to store results 
    # ------------------------------------
    cvErr1 <- array(0, dim = c(neta, ngamma, max.rank))
    cvErr2 <- array(0, dim = c(neta, ngamma, max.rank))
    cvErr3 <- array(0, dim = c(neta, ngamma, max.rank))
    alpha.list <- list()
    Omega.list <- list()

    # ---------------------------
    # Train-validation split
    # ---------------------------
    xtest <- xval
    ytest <- yval
    xtrain <- x
    ytrain <- y
    ntrain <- dim(xtrain)[1]
    ntest <- dim(xtest)[1]


    for (ieta in 1 : neta) {
        eta <- eta.vec[ieta]
        cat("\n\n\n\n=======================================\n")
        cat(ieta,"-th eta:, ", eta, "\n", sep = "")
        cat("=======================================\n")
        alpha.list[[ieta]] <- list()
        Omega.list[[ieta]] <- list()

        # -------------------------
        # Initialization of Omega
        # -------------------------
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
            
            while (abs(llcur - llpre) > 1e-6 & iter <= max.iter) {
                llpre <- llcur
                Omega <- update.Omega.nonconvex(ytrain, alpha, K, Q, eta, diag.cov = diag.cov)
                tmp <- update_alpha(ytrain, K, Q, Omega, gamma, eta, stepsize, alphaPre = alpha)
                alpha <- tmp$alpha
                stepsize <- tmp$stepsize
                llcur <-  eval_obj(ytrain, alpha, Omega, K, Q, gamma, eta)
                iter <- iter + 1
            }
            
            Omega.list[[ieta]][[j]] <- Omega
            alpha.list[[ieta]][[j]] <- alpha

            # ------------------------
            # Mean vector estimation
            # ------------------------
            est.cov <- chol2inv(chol((Omega + t(Omega)) / 2))
            est.cov <- (est.cov + t(est.cov)) / 2
            est.mean <- matrix(0, nrow = dim(x_all)[1], ncol = p)
            for (ii in 1 : total) {
                est.mean[ii, ] <-  Kx.vec(x_all[ii,], xtilde, classes, kernel) %*% alpha / sqrt(ntilde) + col_means - colMeans(Q0 %*% K %*% alpha) * sqrt(ntilde) 
            }

            # -----------------------------
            # Prediction on validation set
            # -----------------------------
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
                x.pred <- x_all[pred.classes, ]
                # --------------------
                # Prediction Accuracy
                # --------------------
                cvErr1[ieta, j, kk] <- sum(rowSums(x.pred == xtest) == M) / ntest
                # --------------------
                # Hamming distance
                # --------------------
                cvErr2[ieta, j, kk] <- mean(colSums(x.pred == xtest) / ntest)
                # --------------------
                # Likelihood
                # --------------------
                cvErr3[ieta, j, kk] <- log_lld
                cat("----->Rank of ", kk, " : ", round(cvErr1[ieta, j, kk], 3), round(cvErr2[ieta, j, kk], 3),round(cvErr3[ieta, j, kk], 3),"------\n", sep = " ")
            }
        }

    }

    # ----------------------
    # Tuning via Accuracy
    # ----------------------
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

    # -----------------------------
    # Tuning via Hamming Distance
    # -----------------------------
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
    
    # -------------------------
    # Tuning via Likelihood
    # -------------------------
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


# ----------------------------------------------
#     Convex Estimator
# ----------------------------------------------
# x: n x M matrix of responses for training
# y: n x p matrix of predictors for training
# xval: nval x M matrix of responses for validation
# yval: nval x p matrix of predictors for validation
# classes: a vector of length M, with each component being the number of categories for each response
# ----------------------------------------------

cv.kernel.convex <- function(x, y, xval, yval,
                            classes, eta.vec, max.rank = 4, 
                            ngamma = 35, delta = 0.001, 
                            kernel = 1, max.iter = 100, 
                            diag.cov = FALSE) 
{
    # -----------------------
    # Preliminaries
    # -----------------------
    n <- dim(x)[1]
    p <- dim(y)[2]
    M <- length(classes)
    neta <- length(eta.vec)
    total <- prod(classes)
    x_all <- as.matrix(do.call(expand.grid, lapply(classes, function(x) 0:(x-1))))

    # ------------------
    #   Centralization
    # ------------------
    col_means <- colMeans(y)
    y <- t(t(y) - col_means)

    # ------------------
    # Set gamma.vec 
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
    
    # ---------------------------------
    # Calculation of sample covariance
    # ---------------------------------
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

    # -------------------------------------------
    # Calcalate the marginal probability of x
    # -------------------------------------------
    prob.list  <- list()
    for (rr in 1 : max.rank) {
        prob.list[[rr]] <- Marg.prob(x, rr, classes)
        cat("Rank", rr, "decomp finished.\n")
    }
    

    # ------------------------------------
    # Arrays and lists to store results 
    # ------------------------------------
    cvErr1 <- array(0, dim = c(neta, ngamma, max.rank))
    cvErr2 <- array(0, dim = c(neta, ngamma, max.rank))
    cvErr3 <- array(0, dim = c(neta, ngamma, max.rank))
    Beta.list <- list()
    Omega.list <- list()

    # ---------------------------
    # Train-validation split
    # ---------------------------
    xtest <- xval
    ytest <- yval
    xtrain <- x
    ytrain <- y
    ntrain <- dim(xtrain)[1]
    ntest <- dim(xtest)[1]


    for (ieta in 1 : neta) {
        eta <- eta.vec[ieta]
        cat("\n\n\n\n=======================================\n")
        cat(ieta,"-th eta:, ", eta, "\n", sep = "")
        cat("=======================================\n")
        Beta.list[[ieta]] <- list()
        Omega.list[[ieta]] <- list()
    
        # -------------------------
        # Initialization of Omega
        # -------------------------
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
                starting <- update.Beta(ytrain, K, Q, Omega, gamma, eta)
            } else {
                starting <- update.Beta(ytrain, K, Q, Omega, gamma, eta, BetaPre = Beta)
            }
            
            Beta <- starting$Beta
            stepsize <- starting$stepsize
            llcur <- eval.obj.convex(ytrain, Beta, Omega, K, Q, gamma, eta)
            iter <- 1
            
            while (abs(llcur - llpre) > 1e-6 & iter <= max.iter) {
                llpre <- llcur
                Omega <- update.Omega.convex(ytrain, Beta, Omega, K, Q, gamma, eta)
                tmp <- update.Beta(ytrain, K, Q, Omega, gamma, eta, stepsize, BetaPre = Beta)
                Beta <- tmp$Beta
                stepsize <- tmp$stepsize
                llcur <-  eval.obj.convex(ytrain, Beta, Omega, K, Q, gamma, eta)
                iter <- iter + 1
            }

            Omega.list[[ieta]][[j]] <- Omega
            Beta.list[[ieta]][[j]] <- Beta

            # ------------------------
            # Mean vector estimation
            # ------------------------
            est.cov <- chol2inv(chol((Omega + t(Omega)) / 2))
            est.cov <- (est.cov + t(est.cov)) / 2
            est.mean <- matrix(0, nrow = dim(x_all)[1], ncol = p)
            for (ii in 1 : total) {
                est.mean[ii, ] <-  Kx.vec(x_all[ii,], xtilde, classes, kernel)%*% Beta %*% est.cov / sqrt(ntilde) + col_means - colMeans(Q0 %*% K %*% Beta %*% est.cov) * sqrt(ntilde)
            }

            # -----------------------------
            # Prediction on validation set
            # -----------------------------
            Pmat <- matrix(0, nrow = ntest, ncol = prod(classes))
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
                x.pred <- x_all[pred.classes, ]
                # --------------------
                # Prediction Accuracy
                # --------------------
                cvErr1[ieta, j, kk] <- sum(rowSums(x.pred == xtest) == M) / ntest
                # --------------------
                # Hamming distance
                # --------------------
                cvErr2[ieta, j, kk] <- mean(colSums(x.pred == xtest) / ntest)
                # --------------------
                # Likelihood
                # --------------------
                cvErr3[ieta, j, kk] <- log_lld

                cat("----->Rank of ", kk, " : ", round(cvErr1[ieta, j, kk], 3), round(cvErr2[ieta, j, kk], 3),round(cvErr3[ieta, j, kk], 3),"------\n", sep = " ")
            }
        }
    }

    # ----------------------
    # Tuning via Accuracy
    # ----------------------
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

    # -----------------------------
    # Tuning via Hamming Distance
    # -----------------------------
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
    
    # -------------------------
    # Tuning via Likelihood
    # -------------------------
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
