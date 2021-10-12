# This is from https://github.com/kolesarm/Robust-Small-Sample-Standard-Errors, October 2021, with minor adjustment for formatting when there
# is only one regressor (gives errors otherwise, but we need precisely this for comparison)
dfadjustSE<-function (model, clustervar = NULL, ell = NULL, IK = TRUE, tol = 1e-09,
    rho0 = FALSE)
{
    Q <- as.matrix(qr.Q(model$qr))
    R <- qr.R(model$qr)
    n <- NROW(Q)
    u <- stats::residuals(model)
    K <- model$rank
    rho <- sig <- NA
    sandwich <- function(meat) backsolve(R, t(backsolve(R, meat)))


    if (is.null(clustervar)) {
        diaghat <- try(stats::hatvalues(model), silent = TRUE)
        AQ <- (1 - diaghat >= tol) * (1/sqrt(pmax(1 - diaghat,
            tol))) * Q
        HC2 <- crossprod(u * AQ)
        HC1 <- n/(n - K) * crossprod(u * Q)
        df0 <- function(ell) {
            a <- drop(AQ %*% backsolve(R, ell, transpose = TRUE))
            B <- a * Q
            (sum(a^2) - sum(B^2))^2/(sum(a^4) - 2 * sum((a *
                B)^2) + sum(crossprod(B)^2))
        }
    } else {
        if (!is.factor(clustervar))
            stop("'clustervar' must be a factor")
        ord <- order(clustervar)
        u <- u[ord]
        Q <- as.matrix(Q[ord, ])
        clustervar <- clustervar[ord]
        S <- nlevels(clustervar)
        uj <- apply(u * Q, 2, function(x) tapply(x, clustervar,
            sum))
        HC1 <- S/(S - 1) * (n - 1)/(n - K) * crossprod(uj)
        AQf <- function(s) {
            Qs <- Q[clustervar == s, , drop = FALSE]
            e <- eigen(crossprod(Qs))
            Ds <- e$vectors %*% ((1 - e$values >= tol) * (1/sqrt(pmax(1 -
                e$values, tol))) * t(e$vectors))
            Qs %*% Ds
        }
        AQ <- lapply(levels(clustervar), AQf)
        AQ <- do.call(rbind, AQ)
        uj <- apply(u * AQ, 2, function(x) tapply(x, clustervar,
            sum))
        HC2 <- crossprod(uj)
        if (IK) {
            ssr <- sum(u^2)
            den <- sum(tapply(u, clustervar, length)^2) - n
            rho <- 0
            if (den > 0)
                rho <- (sum(tapply(u, clustervar, sum)^2) - ssr)/den
            if (rho0)
                rho <- max(rho, 0)
            sig <- max(ssr/n - rho, 0)
        }
        df0 <- function(ell) {
            a <- drop(AQ %*% backsolve(R, ell, transpose = TRUE))
            as <- as.vector(tapply(a^2, clustervar, sum))
            B <- apply(a * Q, 2, function(x) tapply(x, clustervar,
                sum))
            if (!IK) {
                (sum(as) - sum(B^2))^2/(sum(as^2) - 2 * sum(as *
                  B^2) + sum(crossprod(B)^2))
            }
            else {
                D <- as.vector(tapply(a, clustervar, sum))
                Fm <- apply(Q, 2, function(x) tapply(x, clustervar,
                  sum))
                GG <- sig * (diag(as) - tcrossprod(B)) + rho *
                  tcrossprod(diag(D) - tcrossprod(B, Fm))
                sum(diag(GG))^2/sum(GG^2)
            }
        }
    }
    Vhat <- sandwich(HC2)
    VhatStata <- sandwich(HC1)
    if (!is.null(ell)) {
        se <- drop(sqrt(crossprod(ell, Vhat) %*% ell))
        dof <- df0(ell)
        seStata <- drop(sqrt(crossprod(ell, VhatStata) %*% ell))
        beta <- sum(ell * model$coefficients)
    }
    else {
        se <- sqrt(diag(Vhat))
        dof <- vapply(seq(K), function(k) df0(diag(K)[, k]),
            numeric(1))
        seStata <- sqrt(diag(VhatStata))
        beta <- model$coefficients
    }
    r <- cbind(Estimate = beta, `HC1 se` = seStata, `HC2 se` = se,
        `Adj. se` = se * stats::qt(0.975, df = dof)/stats::qnorm(0.975),
        df = dof)
    rownames(r) <- names(beta)
    colnames(Vhat) <- rownames(Vhat) <- names(model$coefficients)
    structure(list(vcov = Vhat, coefficients = r, rho = rho,
        sig = sig), class = "dfadjustSE")
}
