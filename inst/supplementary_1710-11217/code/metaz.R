## NOT A POLISHED PIECE OF PUBLIC-USE SOFTWARE!  Provided "as is",
## licensed under GPL2 or higher.  NO WARRANTY OF FITNESS FOR ANY
## PURPOSE!
##
## Claudia Di Caterina [aut] dicaterina@stat.unipd.it
##
## 26 October 2017
##
## This code still has lots of rough edges, some of which are
## indicated by the embedded comments.


############################# ADJUSTED Z-TESTS for META-REGRESSION ############################

## object: metaLik object

## profile log-lik for psi
profpsi <- function(psi, ystar, Xstar, X, sigma2) {
  tau2 <- psi
  w <- 1/(sigma2 + tau2)
  beta <- coef(lm.wfit(x = Xstar, y = ystar, w = w))
  res <- ystar - Xstar %*% beta
  W.X <- sqrt(w)*X
  XWX <- crossprod(W.X)
  - 0.5*sum(log(w)) + 0.5*sum(w*res^2)
}

## bias of psi
bias <- function(pars, X, sigma2, y, fit = NULL) {
  if (is.null(fit)) {
    fit <- fitfun(pars, X, sigma2, y, info = TRUE)
  }
  with(fit, {
    - sum(w*h)/sum(w^2)
  })
}

fitFun <- function(object) {
  fitM <- fitfun(pars=as.vector(object$mle), X=object$X, sigma2=object$sigma2, y=object$y, info = TRUE)
  parML <- as.vector(object$mle)
  X <- model.matrix(object)
  p <- ncol(X)
  betahat <- object$beta.mle
  psihat <- object$tau2.mle
  sigma2 <- object$sigma2
  W <- diag(1/(sigma2+psihat))
  W2 <- W %*% W
  F <- cbind(rbind(t(X)%*%W%*%X, rep(0,p)), c(rep(0,p), sum(W2)/2))
  vcovMat <- chol2inv(chol(F))
  rownames(vcovMat) <- colnames(vcovMat) <- attr(object$mle, "names")
  biases <- -c(numeric(p), bias(fit=fitM))
  list(y = object$y,
       X = X,
       p = p,
       K = nrow(X),
       sigma2 = object$sigma2,
       W = W,
       W2 = W2,
       weights = object$weights,
       coefs = parML,
       betahat = betahat,
       psihat = psihat,
       vcovMat = vcovMat,
       se = sqrt(diag(vcovMat)),
       bias = biases)
}



fitFunR <- function(object, beta0=0, arg=1) {
  #  fit <- fitFun(object)
  y <- object$y
  X <- model.matrix(object)
  sigma2 <- object$sigma2
  p <- ncol(X)
  XR<-X[,-arg]
  parR <- rep(NA, length(object$mle))
  parR[arg] <- beta0
  betahatR <- rep(NA, length(object$beta.mle))
  betahatR[arg] <- beta0
  if(ncol(XR)==0) {
    psihatR <- optimize(profpsi, interval=c(0, 1000), ystar=y-as.matrix(X[,arg])%*%beta0, Xstar=as.matrix(XR), X=X, sigma2=sigma2)$minimum
    parR[-arg] <- psihatR
  }
  else {
    mR<-metaLik(y ~ offset(beta0*X[,arg]) + XR - 1, sigma2=object$sigma2, data = data.frame(y = y, XR = XR, X = X))
    parR[-arg] <- mR$mle
    betahatR[-arg] <- mR$beta.mle
    psihatR <- mR$tau2.mle
  }
  WR <- diag(1/(sigma2+psihatR))
  W2R <- WR %*% WR
  F <- cbind(rbind(t(X)%*%WR%*%X, rep(0,p)), c(rep(0,p), sum(W2R)/2))
  vcovMatR <- chol2inv(chol(F))
  rownames(vcovMatR) <- colnames(vcovMatR) <- attr(object$mle, "names")
  bR <- bias(pars=parR, X=X, sigma2=sigma2, y=y)
  biasR <- -c(numeric(p), bR)
  list(y = y,
       X = X,
       p = p,
       K = nrow(X),
       sigma2 = sigma2,
       W = WR,
       W2 = W2R,
       weights = object$weights,
       coefs = parR,
       betahat = betahatR,
       psihat = psihatR,
       vcovMat = vcovMatR,
       se = sqrt(diag(vcovMatR)),
       bias = biasR)
}

biases <- function(object, beta0 = 0, arg = 1, at=c("ML", "Restr")) {
  if (at=="ML") {
    fit <- fitFun(object)
    bias <- fit$bias
  }
  if (at=="Restr") {
    fit <- fitFunR(object, beta0=beta0, arg=arg)
    bias <- fit$bias
  }
  return(as.vector(bias))
}

## The derivatives of the standard errors wrt to the parameter arg
dse <- function(object, beta0 = beta0, arg = 1, at=c("ML", "Restr")) {
  if (at=="ML") {
    fit <- fitFun(object)
  }
  if (at=="Restr") {
    fit <- fitFunR(object, beta0=beta0, arg=arg)
  }
  with(fit, {
    Wstar <- diag(-(sigma2+psihat)^(-2))
    vstar <- -sum((sigma2+psihat)^(-3))

    dinfoB <- lapply(seq.int(p), function(var) matrix(0, nrow=p+1, ncol=p+1))
    dinfoP <- lapply(seq.int(1), function(var) cbind(rbind(t(X)%*%Wstar%*%X, rep(0,p)), c(rep(0,p), vstar)))

    dinfo <- c(dinfoB, dinfoP)

    sapply(seq.int(p+1), function(var) {
      -0.5*diag(vcovMat%*%dinfo[[var]]%*%vcovMat)[arg]/se[arg]
    })
  })
}

## The second derivatives of the standard errors wrt to the parameter arg
ddse <- function(object, beta0 = beta0, arg = 1, at=c("ML", "Restr")) {
  if (at=="ML") {
    fit <- fitFun(object)
  }
  if (at=="Restr") {
    fit <- fitFunR(object, beta0=beta0, arg=arg)
  }
  with(fit, {
    Wstar <- diag(-(sigma2+psihat)^(-2))
    Wstarstar <- diag(2*(sigma2+psihat)^(-3))
    vstar <- -sum((sigma2+psihat)^(-3))
    vstarstar <- 3*sum((sigma2+psihat)^(-4))

    dinfoB <- lapply(seq.int(p), function(var) matrix(0, nrow=p+1, ncol=p+1))
    dinfoP <- lapply(seq.int(1), function(var) cbind(rbind(t(X)%*%Wstar%*%X, rep(0,p)), c(rep(0,p), vstar)))

    dinfo <- function(s) {
      c(dinfoB, dinfoP)[[s]]
    }

    ddinfoPP <- cbind(rbind(t(X)%*%Wstarstar%*%X, rep(0,p)), c(rep(0,p), vstarstar))

    ddinfo <- function(s,t) {
      if((is.element(s, 1:p)) & (is.element(t, 1:p))) ddinfo <- matrix(0, p+1, p+1)
      if((is.element(s, 1:p)) & (t==p+1)) ddinfo <- matrix(0, p+1, p+1)
      if((s==p+1) & (is.element(t, 1:p))) ddinfo <- matrix(0, p+1, p+1)
      if((s==p+1) & (t==p+1)) ddinfo <- ddinfoPP
      return(ddinfo)
    }

    der <- function(s, t) {
      - 0.25 * diag(vcovMat%*%dinfo(s)%*%vcovMat)[arg] *
        diag(vcovMat%*%dinfo(t)%*%vcovMat)[arg]/se[arg]^3 +
        0.5 * diag(vcovMat%*%dinfo(s)%*%vcovMat%*%dinfo(t)%*%vcovMat)[arg]/se[arg] +
        0.5 * diag(vcovMat%*%dinfo(t)%*%vcovMat%*%dinfo(s)%*%vcovMat)[arg]/se[arg] -
        0.5 * diag(vcovMat%*%ddinfo(s, t)%*%vcovMat)[arg]/se[arg]
    }
    outer(seq.int(p+1), seq.int(p+1), Vectorize(der))
  })
}

tau <- function(object, beta0 = 0, arg = 1) {
  fit <- fitFun(object)
  with(fit, {
    (coefs[arg] - beta0)*se[arg]^(-1)
  })
}

dtau <- function(object, beta0 = 0, arg = 1, at=c("ML", "Restr")) {
  if (at=="ML") {
    fit <- fitFun(object)
  }
  if (at=="Restr") {
    fit <- fitFunR(object, beta0=beta0, arg=arg)
  }
  tau.val <- tau(object, beta0=beta0, arg=arg)
  dse.val <- dse(object, beta0=beta0, arg=arg, at=at)
  with(fit, {
    evec <- rep(0, p+1)
    evec[arg] <- 1
    (evec - tau.val*dse.val)/se[arg]
  })
}

ddtau <- function(object, beta0 = 0, arg = 1, at=c("ML", "Restr")) {
  if (at=="ML") {
    fit <- fitFun(object)
  }
  if (at=="Restr") {
    fit <- fitFunR(object, beta0=beta0, arg=arg)
  }
  tau.val <- tau(object, beta0=beta0, arg=arg)
  dtau.val <- dtau(object, beta0=beta0, arg=arg, at=at)
  dse.val <- dse(object, beta0=beta0, arg=arg, at=at)
  ddse.val <- ddse(object, beta0=beta0, arg=arg, at=at)
  with(fit, {
    -(dtau.val %*% t(dse.val) +
        dse.val %*% t(dtau.val) +
        tau.val * ddse.val)/se[arg]
  })
}

zstatR <- function(object, beta0 = 0, arg = 1, correction = TRUE, at=c("ML", "Restr")) {
  tau <- tau(object, beta0=beta0, arg=arg)
  nam <- names(tau)
  if (!correction) {
    tau
  }
  else {
    if (at=="ML") {
      fit <- fitFun(object)
    }
    if (at=="Restr") {
      fit <- fitFunR(object, beta0=beta0, arg=arg)
    }
    biasvec <- c(biases(object, beta0=beta0, arg=arg, at=at))

    dtau <-  dtau(object, beta0=beta0, arg=arg, at=at)
    traces <- sum(diag(fit$vcovMat%*%ddtau(object, beta0=beta0, arg=arg, at=at)))
    tauc <- tau - (dtau %*% biasvec + 0.5 * traces)
    tauc <- as.double(tauc)
    names(tauc) <- nam
    tauc
  }
}

score_test<-function(object, beta0 = 0, arg = 1) {
  fit <- fitFun(object)
  fitR <- fitFunR(object, beta0=beta0, arg=arg)
  with(fitR, {
    R <- y - X%*%betahat
    scoreB <- as.vector(t(X)%*%W%*%R)
    scoreP <- 0.5*(t(R)%*%W2%*%R)-sum(W)/2
    score <- c(scoreB, scoreP)
    stest <- t(score) %*% vcovMat %*% score
    as.vector(sign(fit$coefs[arg] - beta0) * sqrt(stest))
  })
}
