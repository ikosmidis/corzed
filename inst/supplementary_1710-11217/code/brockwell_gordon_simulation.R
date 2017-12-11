## NOT A POLISHED PIECE OF PUBLIC-USE SOFTWARE!  Provided "as is",
## licensed under GPL2 or higher.  NO WARRANTY OF FITNESS FOR ANY
## PURPOSE!
##
## Claudia Di Caterina [aut] dicaterina@stat.unipd.it,
## Ioannis Kosmidis [ctb], i.kosmidis@ucl.ac.uk
##
## 26 October 2017
##
## This code still has lots of rough edges, some of which are
## indicated by the embedded comments.

path <- "."

## functionsMPL.R is from the supplementary material of
## https://doi.org/10.1093/biomet/asx001
source(paste(path,"functionsMPL.R", sep = "/"))
## metaz.R implements the location-adjusted Wald test
source(paste(path, "metaz.R", sep = "/"))

library("plyr")
library("metaLik")
library("metatest")
library("parallel")
library("plyr")
library("doMC")

registerDoMC(3)

## as in Zhen and Lin...
generate.sigma2s <- function(K){
  sigma2s <- 0.25 * rchisq(K, df = 1)
  ok <- (sigma2s > 0.009 & sigma2s < 0.6)
  while(sum(ok) < K){
    tmp <- 0.25 * rchisq(K-sum(ok), df = 1)
    sigma2s[!ok] <- tmp
    ok <- (sigma2s > 0.009 & sigma2s < 0.6)
  }
  sigma2s
}

simulate.BG <- function(beta = 0.5,
                        var.re = 0.01,
                        sigma2s = runif(10, 0, 1)) {
  N <- length(sigma2s)
  list(K = N, truebeta = beta, truepsi = var.re, betahat = beta + rnorm(N)*sqrt(var.re + sigma2s),
       betahat.se = sqrt(sigma2s), X = matrix(1, nrow = N))
}

## sample sizes to consider
Ks <- c(seq(5, 50, by = 5), 100, 200)

## variance components to consider
truepsis <- seq(0, 0.1, length = 11)
## other simulation constants
truebeta <- 0.5

## compute X matrices and generate sigma2
sigma2 <- Xmat <- as.list(rep(0, length(Ks)))
## generate sigmas and set Xs
set.seed(123)
for (j in seq_along(Ks)) {
  sigma2[[j]] <- generate.sigma2s(Ks[j])
  Xmat[[j]] <- matrix(1, nrow = Ks[j])
}

## simulation size per sample size
nsimu <- 10000

## nominal significance levels
alpha <- c(0.01, 0.05, 0.1, 0.2, 0.4, 0.8)

sizeML <- sizeBC <- sizeMLl <- sizeBCl <- sizeMLr <- sizeBCr <- sizeplr <- sizeplrl <- sizeplrr <- sizeDL <-
  sizeDLr <- sizeDLl <- sizeBart <- sizeBartl <- sizeBartr <- sizeKH <- sizeKHl <- sizeKHr <-
  sizes <- sizer <- sizerstar <- sizesl <- sizesr <- sizerl <- sizerr <- sizerstarr <- sizerstarl <-
  array(NA,dim=c(length(alpha), length(truepsis), length(Ks)), dimnames=list("alpha"=alpha, "psi"=truepsis, "K"=Ks))

qML <- qBC <- qBCest <- qr <- qs <- qrstar <- qBCR <- qplr <- qDL <- qBart <- qKH <- array(NA,dim=c(length(alpha), length(truepsis),
                                                                                                           length(Ks)), dimnames=list("alpha"=alpha,
                                                                                                                                      "psi"=truepsis, "K"=Ks))

zMLlist <- zBClist <- zBCestlist <- slist <- rlist <- rstarlist <- zBCRlist <- errlist <- plrlist <- DLlist <- Bartlist <- ZLlist <- KHlist <-
  as.list(numeric(length(Ks)))
names(zMLlist) <- names(zBClist) <- names(zBCestlist) <- names(slist) <- names(rlist) <- names(zBCRlist) <- names(rstarlist) <-
  names(plrlist) <- names(DLlist) <- names(Bartlist) <- names(KHlist) <- Ks

## statistics computed on each dataset
simuMeta1 <- function(dati, X, trace=TRUE) {
  if (dati$id %% trace ==0) cat(dati$id, "\n")
  p <- ncol(X)
  ## Refit the full model
  mod <- metaLik(betahat ~ 1, data = dati, sigma2 = betahat.se^2)
  zML <- zstatR(mod, beta0 = truebeta, arg = 1, correction = FALSE, at="ML") # t
  zBC <- zstatR(mod, beta0 = truebeta, arg = 1, correction = TRUE, at="ML") # t*
  s <- score_test(mod, beta0 = truebeta, arg = 1) # score
  stats <- perform_tests(mod$y, mod$X, mod$sigma2, null = truebeta, what = 1, B = 1,
                         silent = FALSE, maxiter = 10000)[, c("statistics")]
  r <- stats["LR"]
  plr <- stats["PLR"]
  rstar <- stats["Skovgaard"]
  Bart <- stats["Bartlett"]
  DL <- stats["DL"]
  KH <- stats["KH"]
  data.frame(zML=zML, zBC=zBC,
             ## zBCR=zBCR, zBCest=zBCest,
             s=s, r=r, rstar=rstar, plr=plr,
             Bart=Bart, DL=DL, KH=KH)
}


simuMetaAll<-function(data_all, trace=TRUE) {

  X <- attr(data_all, "X")
  beta <- attr(data_all, "beta")
  sigma2 <- attr(data_all, "sigma2")
  psi <- attr(data_all, "psi")

  Nsim <- length(data_all)

  out = ldply(.data=data_all, .fun=simuMeta1, X=X, trace=trace, .parallel=TRUE, .inform=TRUE)

  data.frame(out)
}

for(K in Ks) {
  cat(paste("K =", K), "\n")

  out<-array(NA, dim=c(length(truepsis), 9, nsimu), dimnames=list("psi"=truepsis, "test"=c("zML", "zBC",
                                                                                            ## "zBCR", "zBCest",
                                                                                            "s", "r", "rstar",
                                                                                            "plr", "Bart",
                                                                                            ## "ZL",
                                                                                            "DL", "KH"),
                                                                   "sim"=seq.int(nsimu)))

  for(psi in truepsis) {
    cat(paste("psi =", psi), "\n")

    set.seed(123)
    data.sim <- as.list(numeric(nsimu))

    for (i in seq.int(nsimu)) {
      dat <- simulate.BG(beta = truebeta, var.re = psi, sigma2s = sigma2[[which(Ks==K)]])
      data.sim[[i]]<-list(id=i, beta=dat$truebeta, K=dat$K, psi=dat$truepsi, betahat=dat$betahat, X=dat$X, betahat.se=dat$betahat.se)
    }

    attr(data.sim, "beta")=truebeta
    attr(data.sim, "psi")=psi
    attr(data.sim, "sigma2")=sigma2[[which(Ks==K)]]
    attr(data.sim, "X")=Xmat[[which(Ks==K)]]
    attr(data.sim, "Nsim")=nsimu

    resMeta<-simuMetaAll(data.sim, trace=100)

    out[which(truepsis==psi),,] <- t(as.matrix(resMeta))

  }
  zMLlist[[as.character(K)]] <- out[,"zML",]
  zBClist[[as.character(K)]] <- out[,"zBC",]
  slist[[as.character(K)]] <- out[,"s",]
  rlist[[as.character(K)]] <- out[,"r",]
  rstarlist[[as.character(K)]] <- out[,"rstar",]
  plrlist[[as.character(K)]] <- out[,"plr",]
  Bartlist[[as.character(K)]] <- out[,"Bart",]
  DLlist[[as.character(K)]] <- out[,"DL",]
  KHlist[[as.character(K)]] <- out[,"KH",]

  ### empirical size
  ## t
  # two-sided
  sizeML[,,which(Ks==K)]<- t(sapply(alpha, function(a) {
    quant <- qnorm(1 - a/2)
    apply(out[,"zML",], 1, function(x) mean(((x > quant) | (x < -quant)), na.rm = TRUE))
  }))
  # right-sided
  sizeMLr[,,which(Ks==K)]<- t(sapply(alpha, function(a) {
    quant <- qnorm(1 - a)
    apply(out[,"zML",], 1, function(x) mean(x > quant, na.rm = TRUE))
  }))
  # left-sided
  sizeMLl[,,which(Ks==K)]<- t(sapply(alpha, function(a) {
    quant <- qnorm(a)
    apply(out[,"zML",], 1, function(x) mean(x < quant, na.rm = TRUE))
  }))

  ## t*
  sizeBC[,,which(Ks==K)]<- t(sapply(alpha, function(a) {
    quant <- qnorm(1 - a/2)
    apply(out[,"zBC",], 1, function(x) mean(((x > quant) | (x < -quant)), na.rm = TRUE))
  }))
  sizeBCr[,,which(Ks==K)]<- t(sapply(alpha, function(a) {
    quant <- qnorm(1 - a)
    apply(out[,"zBC",], 1, function(x) mean(x > quant, na.rm = TRUE))
  }))
  sizeBCl[,,which(Ks==K)]<- t(sapply(alpha, function(a) {
    quant <- qnorm(a)
    apply(out[,"zBC",], 1, function(x) mean(x < quant, na.rm = TRUE))
  }))

  # score stat
  sizes[,,which(Ks==K)]<- t(sapply(alpha, function(a) {
    quant <- qnorm(1 - a/2)
    apply(out[,"s",], 1, function(x) mean(((x > quant) | (x < -quant)), na.rm = TRUE))
  }))
  sizesr[,,which(Ks==K)]<- t(sapply(alpha, function(a) {
    quant <- qnorm(1 - a)
    apply(out[,"s",], 1, function(x) mean(x > quant, na.rm = TRUE))
  }))
  sizesl[,,which(Ks==K)]<- t(sapply(alpha, function(a) {
    quant <- qnorm(a)
    apply(out[,"s",], 1, function(x) mean(x < quant, na.rm = TRUE))
  }))

  # likelihood ratio stat
  sizer[,,which(Ks==K)]<- t(sapply(alpha, function(a) {
    quant <- qnorm(1 - a/2)
    apply(out[,"r",], 1, function(x) mean(((x > quant) | (x < -quant)), na.rm=TRUE))
  }))
  sizerr[,,which(Ks==K)]<- t(sapply(alpha, function(a) {
    quant <- qnorm(1 - a)
    apply(out[,"r",], 1, function(x) mean(x > quant, na.rm=TRUE))
  }))
  sizerl[,,which(Ks==K)]<- t(sapply(alpha, function(a) {
    quant <- qnorm(a)
    apply(out[,"r",], 1, function(x) mean(x < quant, na.rm = TRUE))
  }))

  # modified likelihood ratio stat
  sizerstar[,,which(Ks==K)]<- t(sapply(alpha, function(a) {
    quant <- qnorm(1 - a/2)
    apply(out[,"rstar",], 1, function(x) mean(((x > quant) | (x < -quant)), na.rm=TRUE))
  }))
  sizerstarr[,,which(Ks==K)]<- t(sapply(alpha, function(a) {
    quant <- qnorm(1 - a)
    apply(out[,"rstar",], 1, function(x) mean(x > quant, na.rm=TRUE))
  }))
  sizerstarl[,,which(Ks==K)]<- t(sapply(alpha, function(a) {
    quant <- qnorm(a)
    apply(out[,"rstar",], 1, function(x) mean(x < quant, na.rm=TRUE))
  }))

  # penalized likelihood ratio stat
  sizeplr[,,which(Ks==K)]<- t(sapply(alpha, function(a) {
    quant <- qnorm(1 - a/2)
    apply(out[,"plr",], 1, function(x) mean(((x > quant) | (x < -quant)), na.rm=TRUE))
  }))
  sizeplrr[,,which(Ks==K)]<- t(sapply(alpha, function(a) {
    quant <- qnorm(1 - a)
    apply(out[,"plr",], 1, function(x) mean(x > quant, na.rm=TRUE))
  }))
  sizeplrl[,,which(Ks==K)]<- t(sapply(alpha, function(a) {
    quant <- qnorm(a)
    apply(out[,"plr",], 1, function(x) mean(x < quant, na.rm=TRUE))
  }))

  # Bartlett-corrected likelihood ratio stat
  sizeBart[,,which(Ks==K)]<- t(sapply(alpha, function(a) {
    quant <- qnorm(1 - a/2)
    apply(out[,"Bart",], 1, function(x) mean(((x > quant) | (x < -quant)), na.rm=TRUE))
  }))
  sizeBartr[,,which(Ks==K)]<- t(sapply(alpha, function(a) {
    quant <- qnorm(1 - a)
    apply(out[,"Bart",], 1, function(x) mean(x > quant, na.rm=TRUE))
  }))
  sizeBartl[,,which(Ks==K)]<- t(sapply(alpha, function(a) {
    quant <- qnorm(a)
    apply(out[,"Bart",], 1, function(x) mean(x < quant, na.rm=TRUE))
  }))

  # DerSimonian & Laird method
  sizeDL[,,which(Ks==K)]<- t(sapply(alpha, function(a) {
    quant <- qnorm(1 - a/2)
    apply(out[,"DL",], 1, function(x) mean(((x > quant) | (x < -quant)), na.rm=TRUE))
  }))
  sizeDLr[,,which(Ks==K)]<- t(sapply(alpha, function(a) {
    quant <- qnorm(1 - a)
    apply(out[,"DL",], 1, function(x) mean(x > quant, na.rm=TRUE))
  }))
  sizeDLl[,,which(Ks==K)]<- t(sapply(alpha, function(a) {
    quant <- qnorm(a)
    apply(out[,"DL",], 1, function(x) mean(x < quant, na.rm=TRUE))
  }))

  # modified DL
  sizeKH[,,which(Ks==K)]<- t(sapply(alpha, function(a) {
    quant <- qnorm(1 - a/2)
    apply(out[,"KH",], 1, function(x) mean(((x > quant) | (x < -quant)), na.rm=TRUE))
  }))
  sizeKHr[,,which(Ks==K)]<- t(sapply(alpha, function(a) {
    quant <- qnorm(1 - a)
    apply(out[,"KH",], 1, function(x) mean(x > quant, na.rm=TRUE))
  }))
  sizeKHl[,,which(Ks==K)]<- t(sapply(alpha, function(a) {
    quant <- qnorm(a)
    apply(out[,"KH",], 1, function(x) mean(x < quant, na.rm=TRUE))
  }))

  ### empirical quantiles
  ## t
  qML[,,which(Ks==K)]<-t(sapply(alpha, function(a) {
    apply(out[,"zML",], 1, function(x) quantile(x, probs=a, na.rm = TRUE))
  }))

  ## t*
  qBC[,,which(Ks==K)]<-t(sapply(alpha, function(a) {
    apply(out[,"zBC",], 1, function(x) quantile(x, probs=a, na.rm = TRUE))
  }))

  # score stat
  qs[,,which(Ks==K)]<-t(sapply(alpha, function(a) {
    apply(out[,"s",], 1, function(x) quantile(x, probs=a, na.rm = TRUE))
  }))

  # likelihood ratio stat
  qr[,,which(Ks==K)]<-t(sapply(alpha, function(a) {
    apply(out[,"r",], 1, function(x) quantile(x, probs=a, na.rm = TRUE))
  }))

  # modified likelihood ratio stat
  qrstar[,,which(Ks==K)]<-t(sapply(alpha, function(a) {
    apply(out[,"rstar",], 1, function(x) quantile(x, probs=a, na.rm = TRUE))
  }))

  # penalized likelihood ratio stat
  qplr[,,which(Ks==K)]<-t(sapply(alpha, function(a) {
    apply(out[,"plr",], 1, function(x) quantile(x, probs=a, na.rm = TRUE))
  }))

  # Bartlett-corrected likelihood ratio stat
  qBart[,,which(Ks==K)]<-t(sapply(alpha, function(a) {
    apply(out[,"Bart",], 1, function(x) quantile(x, probs=a, na.rm = TRUE))
  }))

  # DerSimonian & Laird method
  qDL[,,which(Ks==K)]<-t(sapply(alpha, function(a) {
    apply(out[,"DL",], 1, function(x) quantile(x, probs=a, na.rm = TRUE))
  }))

  # modified DL
  qKH[,,which(Ks==K)]<-t(sapply(alpha, function(a) {
    apply(out[,"KH",], 1, function(x) quantile(x, probs=a, na.rm = TRUE))
  }))

}

save(truebeta, Ks, truepsis, sizeML, sizeMLl, sizeMLr, sizeBC, sizeBCl, sizeBCr, sizeKH, sizeKHl, sizeKHr, qKH,
     zMLlist, zBClist, sizer, sizerr, sizerl, sizes, sizesl, sizesr, sizeBart, sizeBartl, sizeBartr, qBart, Bartlist,
     sizerstar, sizerstarl, sizerstarr, slist, rlist, rstarlist, qML, qBC, qr, nsimu, sizeDL, sizeDLl, sizeDLr, qDL, DLlist, KHlist,
     qs, qrstar, alpha, sigma2, Xmat, sizeplr, sizeplrl, sizeplrr, qplr, plrlist,
     file = "results/brockwell_gordon_simulation.rda")
