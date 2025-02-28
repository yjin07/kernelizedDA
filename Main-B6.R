Sys.time()
source("../Functions/helper-all-revised.R")
library(glmnet)
library(TULIP)
uu <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

set.seed(uu + 1000)
combs <- expand.grid(p = rep(c(50, 100, 150), each = 100), rho = c(0.8, 1.0, 1.2, 1.4))
final <- list()

# -------------------------
# Simulation setup 
# -------------------------
M  <-  6
p  <-  combs$p[uu]
c  <-  2
ntest  <-  1000
ntrain  <-  200
nval <- ntrain
n <- ntest + ntrain + nval
classes <- rep(c, M)
rho  <-  combs$rho[uu]
cMat <- covAR(p, 0.7)
Omega.true <- solve(cMat)
beta <- matrix(0, nrow = 10, ncol = M)
inds <- sample(1 : (10 * M), 5 * M, replace = FALSE)
beta[inds] <- 2
y_all <- as.matrix(do.call(expand.grid, lapply(classes, function(x) 0:(x-1))))
kernel <- 1

# -----------------
# Store Result
# -----------------
result <- matrix(0, nrow = 9, ncol = M + 5)
meandiff <- matrix(0, nrow = 2, ncol = 2)

# -------------------------------
# Random Generation of Y 
# -------------------------------
Y <- matrix(0, ncol = M, nrow = n)
probs <- runif(prod(classes))
probs <- probs / sum(probs)
draw <- sample(1 : prod(classes), n, prob = probs, replace = T)
for (i in 1 : n) {
    Y[i,] <- y_all[draw[i],]
}

# -------------------------------
# Random Generation of the Mean 
# -------------------------------
nn <- dim(y_all)[1]
nonMat <- matrix(0, nrow = nn, ncol =  10)
for (i in 1 : nn) {
    nonMat[i,] <- fcn(y_all[i, ], beta, all.linear = TRUE)
}
nonMat <- nonMat / rho
mean.ind <- sample(1:p, 10, replace = F)
meanMat <- matrix(0, nrow = n, ncol = p)
meanMat[, mean.ind] <- nonMat[draw, ]; meanMat <- meanMat %*% cMat # TODO:

# ---------------------
# True MEAN
# ---------------------
truemean <- matrix(0, nrow = nn, ncol = p)
truemean[, mean.ind] <- nonMat; truemean <- truemean %*% cMat # TODO:

R <- matrix(0, nrow = n, ncol = p)
for (i in 1 : n) {
    R[i, ] <- rmvnorm(1, mean = meanMat[i, ], sigma = cMat)
}

# -------------------------------------------
# Split of training, validation and testing
# -------------------------------------------
index <- sample(1:n, ntrain + nval, replace = F)
Xtest <- R[-index, ]
Ytest <- Y[-index, ]
X <- R[index, ]
Y <- Y[index, ]
index <- sample(1 : (ntrain + nval), ntrain, replace = F)
Xval <- X[-index, ]
Yval <- Y[-index, ]
X <- X[index, ]
Y <- Y[index, ]

# --------------------------------------
# Check number of occured categories
# --------------------------------------
tmp <- Y
cc <- rep(0, dim(tmp)[1])
for(i in 1:length(cc)) {
    cc[i] <- findIndex(tmp[i, ], classes)
}
count <- length(levels(factor(cc)))
cat("----------------\ncount: ", count,"\n----------------\n")

samplemean <- matrix(0, nrow = dim(y_all)[1], ncol = p)
for (kk in 1:prod(classes)) {
    if (is.null(dim(X[which(cc==kk), ]))) {
        samplemean[kk, ] <- X[which(cc==kk), ]
    } else {
        samplemean[kk, ] <- colMeans(X[which(cc == kk), ])
    }
}
na.rows <- na.action(na.omit(samplemean))

# ---------------------------------
# Extract important variables 
# ---------------------------------
theta <- solve(cMat, t(truemean))
enorms <- apply(theta, 1, function(x) sqrt(sum(x ^ 2)))
imps <- sort(union(mean.ind, which(enorms > 1e-9)))
final[['VOI']] <- imps
if (2 * length(imps) > p) {
    combined <- 1 : p
} else {
    noise <- sample((1:p)[-imps], length(imps), replace = F)
    combined <- union(imps, noise)
}

# ----------------------
# Method 1: Oracle
# ----------------------
Pmat <- matrix(0, nrow = ntest, ncol = prod(classes))
Omega <- solve(cMat)
for (jj in 1:prod(classes)) {
    Xcent <- t(t(Xtest) - truemean[jj, ])
    Pmat[, jj] <- - diag(Xcent %*% tcrossprod(Omega, Xcent)) + 2 * log(probs[jj])
} 
pred.classes <- apply(Pmat, 1, which.max)
Y.pred <- matrix(0, ncol = M, nrow = ntest)
Y.pred <- y_all[pred.classes, ] 
result[1, 1 : M] <- colSums(Y.pred == Ytest)/ntest
result[1, M + 1] <- sum(rowSums(Y.pred == Ytest) == M) / ntest


# -----------------------------------
# Method 2: KLDA-M (Nonconvex)
# -----------------------------------
eta.vec <- 2 ^ seq(2, -4, length = 5)
res <- cv.kernel.nonconvex(Y, X, Yval, Xval, classes, eta.vec, delta = 1e-6)
final[['nonconvex']] <- res
est.mean <- matrix(0, nrow = dim(y_all)[1], ncol = p)
for (i in 1:prod(classes)) {
    est.mean[i, ] <-  Kx.vec(y_all[i,], res$xtilde, classes, kernel)%*% res$alpha / sqrt(dim(res$xtilde)[1]) + colMeans(X) - colMeans(res$Q0 %*% res$K %*% res$alpha) * sqrt(dim(res$xtilde)[1])
}
# --------------------------------------------------------------------------
tt <- 1
if (is.null(na.rows)) {
    meandiff[tt, 1] <- sum((samplemean - truemean) ^ 2) / nn
    meandiff[tt, 2] <- sum((est.mean - truemean) ^ 2) / nn
} else {
    meandiff[tt, 1] <- sum((samplemean[-na.rows, ] - truemean[-na.rows, ]) ^ 2) / dim(samplemean[-na.rows, ])[1]
    meandiff[tt, 2] <- sum((est.mean[-na.rows, ] - truemean[-na.rows, ]) ^ 2) / dim(samplemean[-na.rows, ])[1]
}
# --------------------------------------------------------------------------
freq <- res$freq
Pmat <- matrix(0, nrow = ntest, ncol = prod(classes))
for (jj in 1:prod(classes)) {
    Xcent <- t(t(Xtest) - est.mean[jj, ])
    Pmat[, jj] <- - diag(Xcent %*% tcrossprod(res$Omega, Xcent)) + 2 * log(freq[jj])
} 
pred.classes <- apply(Pmat, 1, which.max)
Y.pred <- matrix(0, ncol = M, nrow = ntest)
Y.pred <- y_all[pred.classes, ]
result[2, 1 : M] <- colSums(Y.pred == Ytest)/ntest
result[2, M + 1] <- sum(rowSums(Y.pred == Ytest) == M) / ntest
enorms <- apply(res$alpha, 2, function(x) sqrt(sum(x ^ 2)))
imp.mod <- which(enorms > 0)
result[2, M + 2] <- length(intersect(imps, imp.mod)) / length(imps)
result[2, M + 3] <- length(intersect((1:p)[-imps], (1:p)[-imp.mod])) / (p - length(imps))
result[2, M + 4] <- length(intersect(mean.ind, imp.mod)) / length(mean.ind)
result[2, M + 5] <- length(intersect((1:p)[-mean.ind], (1:p)[-imp.mod])) / (p - length(mean.ind))
final[['pred.nonconvex']] <- Y.pred

# -----------------------------------
# Method 3: KLDA-D (Convex)
# -----------------------------------
eta.vec <- 2 ^ seq(2, -4, length = 5)
res <- cv.kernel.convex(Y, X, Yval, Xval, classes, eta.vec, delta = 1e-6)
final[['convex']] <- res
est.mean <- matrix(0, nrow = dim(y_all)[1], ncol = p)
est.cov <- solve((res$Omega + t(res$Omega)) / 2)
for (i in 1:prod(classes)) {
    est.mean[i, ] <-  Kx.vec(y_all[i,], res$xtilde, classes, kernel) %*% res$Beta %*% est.cov / sqrt(dim(res$xtilde)[1]) + colMeans(X) - colMeans(res$Q0 %*% res$K %*% res$Beta %*% est.cov) * sqrt(dim(res$xtilde)[1])
}
tt <- 2
# --------------------------------------------------------------------------
if (is.null(na.rows)) {
    meandiff[tt, 1] <- sum((samplemean - truemean) ^ 2) / nn
    meandiff[tt, 2] <- sum((est.mean - truemean) ^ 2) / nn
} else {
    meandiff[tt, 1] <- sum((samplemean[-na.rows, ] - truemean[-na.rows, ]) ^ 2) / dim(samplemean[-na.rows, ])[1]
    meandiff[tt, 2] <- sum((est.mean[-na.rows, ] - truemean[-na.rows, ]) ^ 2) / dim(samplemean[-na.rows, ])[1]
}
freq <- res$freq
Pmat <- matrix(0, nrow = ntest, ncol = prod(classes))
for (jj in 1:prod(classes)) {
    Xcent <- t(t(Xtest) - est.mean[jj, ])
    Pmat[, jj] <- - diag(Xcent %*% tcrossprod(res$Omega, Xcent)) + 2 * log(freq[jj])
} 
pred.classes <- apply(Pmat, 1, which.max)
Y.pred <- matrix(0, ncol = M, nrow = ntest)
Y.pred <- y_all[pred.classes, ]
result[3, 1 : M] <- colSums(Y.pred == Ytest)/ntest
result[3, M + 1] <- sum(rowSums(Y.pred == Ytest) == M) / ntest
enorms <- apply(res$Beta, 2, function(x) sqrt(sum(x ^ 2)))
imp.mod <- which(enorms > 0)
result[3, M + 2] <- length(intersect(imps, imp.mod)) / length(imps)
result[3, M + 3] <- length(intersect((1:p)[-imps], (1:p)[-imp.mod])) / (p - length(imps))
result[3, M + 4] <- length(intersect(mean.ind, imp.mod)) / length(mean.ind)
result[3, M + 5] <- length(intersect((1:p)[-mean.ind], (1:p)[-imp.mod])) / (p - length(mean.ind))
final[['pred.convex']] <- Y.pred


# --------------------------------
# Method 4: Sep-Logistic
# --------------------------------
xtrain <- X
ytrain <- Y
xval <- Xval
yval <- Yval
lambda_vals <- 10 ^ seq(log10(1e3), log10(1e-3), length = 35)
Y.pred0 <- matrix(0, nrow = ntest, ncol = M)
imp.Logistic <- NULL
for (k in 1:M) {
    accuracy <- rep(0, length(lambda_vals))
    for (j in 1:length(lambda_vals)) {
        fit <- glmnet(x = xtrain, y = ytrain[, k], family = 'multinomial', 
                    alpha = 1, lambda = lambda_vals[j])
        pred <- as.numeric(predict(fit, newx = xval, type = "class"))
        accuracy[j] <- sum(pred == yval[,k]) / length(yval[,k])
    }

    cat("For", k, "-th response: ", which.max(accuracy), "\n")
    lambda_opt <- lambda_vals[which.max(accuracy)]

    fit_R <- glmnet(xtrain, ytrain[,k], family = 'multinomial', 
                    alpha = 1, lambda = lambda_opt)
    for (jj in 1:length(fit_R$beta)) {
        imp.Logistic <- union(imp.Logistic, which(fit_R$beta[[jj]] != 0))
    }            
    Y.pred0[, k] <- as.numeric(predict(fit_R, Xtest, type = 'class'))
}
result[4, 1 : M] <- colSums(Y.pred0 == Ytest)/ntest
result[4, M + 1] <- sum(rowSums(Y.pred0 == Ytest) == M) / ntest
result[4, M + 2] <- length(intersect(imps, imp.Logistic)) / length(imps)
result[4, M + 3] <- length(intersect((1:p)[-imps], (1:p)[-imp.Logistic])) / (p - length(imps))
final[['pred.sep.logistic']] <- Y.pred0

# ----------------------
# Method 5: Sep-MSDA
# ----------------------
xtest <- Xval
ytest <- Yval
xtrain <- X
ytrain <- Y
Y.pred.msda <- matrix(0, nrow = ntest, ncol = M)
imp.MSDA <- NULL
for (jj in 1: M) {
    obj <- dsda(x = xtrain, y = ytrain[,jj] + 1)
    res <- predict(obj, xtest) - 1
    nlambda <- dim(res)[2]
    acc <- rep(0, nlambda)
    for (i in 1 : nlambda) {
        acc[i] <- sum(res[, i] == ytest[,jj]) / dim(res)[1]
    }
    obj <- dsda(x = X, y = Y[,jj] + 1, lambda = obj$lambda[which.max(acc)])
    Y.pred.msda[,jj] <- (predict(obj, Xtest) -  1)
    cat("Response: ", jj, ":", which(obj$beta[-1] != 0), "\n")
    imp.MSDA <- union(imp.MSDA, which(obj$beta[-1] != 0))
}
result[5, 1 : M] <- colSums(Y.pred.msda == Ytest)/ntest
result[5, M + 1] <- sum(rowSums(Y.pred.msda == Ytest) == M) / ntest
result[5, M + 2] <- length(intersect(imps, imp.MSDA)) / length(imps)
result[5, M + 3] <- length(intersect((1:p)[-imps], (1:p)[-imp.MSDA])) / (p - length(imps))
final[['pred.sep.msda']] <- Y.pred.msda



# ---------------------------------
# Method 6: Combined Logistic
# ---------------------------------
tmp <- Y
cc <- rep(0, dim(tmp)[1])
for(i in 1:length(cc)) {
    cc[i] <- findIndex(tmp[i, ], classes)
}
count <- length(levels(factor(cc)))
cat("\nNum of appeared categories: ", count,"\n")
Yf <- cc

tmp <- Yval
cc <- rep(0, dim(tmp)[1])
for(i in 1:length(cc)) {
    cc[i] <- findIndex(tmp[i, ], classes)
}
count <- length(levels(factor(cc)))
cat("\nNum of appeared categories: ", count,"\n")
Yvalf <- cc

tmp <- Ytest
cc <- rep(0, dim(tmp)[1])
for(i in 1:length(cc)) {
    cc[i] <- findIndex(tmp[i, ], classes)
}
count <- length(levels(factor(cc)))
cat("\nNum of appeared categories: ", count,"\n")
Ytestf <- cc

# ------------------------------------------------------------
# Drop categories whose number of appearance is less than 2
# ------------------------------------------------------------
counts <- table(Yf)
selected <- names(counts[counts >= 2])
Yf_s <- Yf[Yf %in% selected]
cat("\nThe responses left after we drop those with observations less than 2: \n")
X_s <- X[Yf %in% selected, ]
X_train <- X_s
Y_train <- Yf_s
X_val <- Xval
Y_val <- Yvalf

mod <- glmnet(x = X_train, y = Y_train, family = "multinomial",
                alpha = 1)
pred <- predict(mod, newx = X_val, type = "class")
accuracy <- rep(0, dim(pred)[2])
for (i in 1:(dim(pred)[2])) {
    accuracy[i] <- sum(as.numeric(pred[,i]) == Y_val) / length(Y_val)
}

cat("======================================\n")
cat(which.max(accuracy), "out of", length(accuracy), "\n")
cat("======================================\n")
best.ind <- which.max(accuracy)
pred <- as.numeric(predict(mod, newx = Xtest, type = "class")[,best.ind])
imp.comb.Logistic <- NULL
for (jj in 1:length(mod$beta)) {
    imp.comb.Logistic <- union(imp.comb.Logistic, which(mod$beta[[jj]][,best.ind] != 0))
}
accuracy <- sum(pred == Ytestf) / ntest
Ypred <- matrix(0, nrow = dim(Ytest)[1], ncol = dim(Ytest)[2])
for (ii in 1 : (dim(Ypred)[1])) {
    Ypred[ii, ] <- y_all[pred[ii], ]
}
result[6, 1 : M] <- colSums(Ypred == Ytest)/ntest
result[6, M + 1] <- sum(rowSums(Ypred == Ytest) == M) / ntest
result[6, M + 2] <- length(intersect(imps, imp.comb.Logistic)) / length(imps)
result[6, M + 3] <- length(intersect((1:p)[-imps], (1:p)[-imp.comb.Logistic])) / (p - length(imps))
final[['pred.comb.logistic']] <- Ypred

# ---------------------------------
# Method 7: Combined-MSDA
# ---------------------------------
xval <- X_val
yval <- Y_val
xtrain <- X_train
ytrain <- Y_train
xtest <- Xtest
ytest <- Ytestf
imp.MSDA <- NULL

nntrain <- length(ytrain)
nnval <- length(yval)
nntest <- length(ytest)
y_comb0 <- c(ytrain, yval, ytest)
y_comb <- as.numeric(factor(y_comb0, levels = unique(y_comb0)))
ytrain <- y_comb[1:nntrain]
yval <- y_comb[(nntrain + 1):(nntrain + nnval)]
ytest <- y_comb[-(1:(nntrain + nnval))]

obj <- msda(x = xtrain, y = ytrain, nlambda = 35)
res <- predict(obj, xval)
nlambda <- dim(res)[2]
acc <- rep(0, nlambda)
for (i in 1 : nlambda) {
    acc[i] <- sum(res[, i] == yval) / dim(res)[1]
}
imp.comb.MSDA <- NULL
obj <- msda(x = xtrain, y = ytrain, lambda = obj$lambda[which.max(acc)])
imp.comb.MSDA <- sort(union(imp.comb.MSDA, as.numeric(which(obj$beta$'1'[,1] != 0, arr.ind = T))))
Y.pred.msda <- predict(obj, xtest)
Y.original <- rep(NA, nntest)
for (jj in 1 : nntest) {
    Y.original[jj] <- y_comb0[which(y_comb == Y.pred.msda[jj])][1]
}
Ypred <- matrix(0, nrow = dim(Ytest)[1], ncol = dim(Ytest)[2])
for (ii in 1 : (dim(Ypred)[1])) {
    Ypred[ii, ] <- y_all[Y.original[ii], ]
}
result[7, 1 : M] <- colSums(Ypred == Ytest)/ntest
result[7, M + 1] <- sum(rowSums(Ypred == Ytest) == M) / ntest
result[7, M + 2] <- length(intersect(imps, imp.comb.MSDA)) / length(imps)
result[7, M + 3] <- length(intersect((1:p)[-imps], (1:p)[-imp.comb.MSDA])) / (p - length(imps))
final[['pred.comb.msda']] <- Ypred
final[['meandiff']] <- meandiff

# ----------------------------------------
# Method 8: Sep-MDA with all predictors
# ----------------------------------------
log_info("Start fitting S-MDA regression...")
xtest <- Xval
ytest <- Yval
xtrain <- X
ytrain <- Y
Y.pred.mda <- matrix(0, nrow = ntest, ncol = M)
imp.MDA <- NULL
imps

for (jj in 1:M) {
    response <- ytrain[, jj] + 1
    obj <- mda(response ~ xtrain, subclasses = 27)
    Y.pred.mda[, jj] <- as.numeric(predict(obj, Xtest)) - 1
}

colSums(Y.pred.mda == Ytest)/ntest
sum(rowSums(Y.pred.mda == Ytest) == M) / ntest
result[8, 1 : M] <- colSums(Y.pred.mda == Ytest)/ntest
result[8, M + 1] <- sum(rowSums(Y.pred.mda == Ytest) == M) / ntest


# ----------------------------------------------
# Method 9: Sep-MDA with important predictors
# ----------------------------------------------
Y.pred.mda <- matrix(0, nrow = ntest, ncol = M)
imp.MDA <- NULL
imps

for (jj in 1:M) {
    response <- ytrain[, jj] + 1
    obj <- mda(response ~ xtrain[, imps], subclasses = 27)
    Y.pred.mda[, jj] <- as.numeric(predict(obj, Xtest[, imps])) - 1
}

colSums(Y.pred.mda == Ytest)/ntest
sum(rowSums(Y.pred.mda == Ytest) == M) / ntest
result[9, 1 : M] <- colSums(Y.pred.mda == Ytest)/ntest
result[9, M + 1] <- sum(rowSums(Y.pred.mda == Ytest) == M) / ntest





result
final[['acc']] <- result
saveRDS(final, file = paste("Results-B6/Res_",uu,".RDS", sep = ""))
