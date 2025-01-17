## These are the R functions needed to run ADELLE.
library(mpoly)

pregrid.points <- function(minP, np) {
	  ## Create a grid of p-values to be used when computing the l-value grid.
    ## The grid of p-values are spaced uniformly on the log scale between
    ## minP and 1.0, with the total number of points being np+1. It's useful
    ## when picking the number of points to think about the number of grid
    ## points you want for each power of 10 in the p-values. So if the minP
    ## for the grid is 10^-20, if you might want 500 points per power of 10,
    ## you would get np = 10000.
    ##
    ## Input:
    ## minP: The smallest p-value in the grid. This should be the smallest
    ##       p-value you are likely to see.
    ## np: The number of grid points. In fact the number of grid points is np+1
    ##
    ## Output: A vector of p-values to use for the l-value grid

    if (minP > 1 || minP <= 0) {
        stop("minP must be between 0 and 1")
    }
    H <- np + 1
    p.gr <- minP^((H - 1:H) / (H - 1))
    ## The last element of p.gr is exactly 1, but this can result in NaN
    ## in some later computations, so let's make that last point the closest
    ## double to 1 that is still less than 1.
    p.gr[np + 1] <- 1 - .Machine$double.neg.eps
    p.gr
}

precomp.lvalues <- function(pval, N.ord, cor.mat) {
    ## Compute the l-values for each of the p-values in pval, for each
    ## of the p-values being the 1...N.ord order statistic. This will create
    ## a matrix of l-values with each column being a p-value and the rows
    ## being the index of the order statistic.
    ## Input:
    ## pval: vector of p-values, sorted from small to large
    ## N.ord: sorted p-values are the order statistic with index 1...N.ord
    ## cor.mat: correlation matrix of the Z scores that gave pval
    ##
    ## Output:
    ## lv.grid: Matrix of precomputed l-values. Each column is for one of the
    ##          pvalues in the input pval vector. The first row gives the l-value
    ##          for that p-value if it were the first order statistic. The
    ##          second row gives the l-value if the p-value were the second
    ##          order statistic, etc. There are N.ord rows.
    ## pvals: The input pval vector. These are the p-values for each column of
    ##        the grid for lv.grid.
    ## gammas: The gamma parameter for the beta-binomial for each p-value
    ## prep: The output of prep.lvalue

    lval.prep <- prep.lvalue(cor.mat)
    N.tests <- lval.prep$N.tests
    if (N.ord > N.tests) {
        stop("N.ord cannot be larger than the number of rows in the correlation matrix")
    }
    gam.p <- precomp.gamma(pval, lval.prep)

    ## M.mat will hold the pmf for each of the probability parameter values
    ## (each column) for each order statistic (each row). I'll get the cdf
    ## by summing the columns (cumsum to do the summation once)
    M.mat <- lv.pmf.mat(pval, gam.p, N.tests)
    lv.grid <- apply(M.mat, 2, function(x) rev(cumsum(rev(x)))[1:N.ord])
    list(lv.grid = lv.grid, pvals = pval, gammas = gam.p, prep = lval.prep)
}

get.lvalue <- function(pval, ord.idx, precomp) {
    ## Given a p-value and which order statistic it is, return its l-value
    ## by looking it up in the grid. If it is not exactly a p-value grid
    ## point, linearly interpolate the l-value from the two l-values that
    ## are immediately larger and smaller. If the p-value is smaller or larger
    ## than the bounds of the grid, compute the l-value. This only computes
    ## the l-value for p-values whose order statistic is less than or equal
    ## to the maximum order statistic given in the precomputation.
    ##
    ## Input:
    ## pval: The p-value for which to compute an l-value. May be a vector
    ## ord.idx: Which order statistic the p-value (or each p-value) is
    ## precomp: The output from precomp.lvalues
    ##
    ## Output:
    ## Vector of l-values, one for each input p-value

    if (length(pval) != length(ord.idx)) {
        stop("pval and ord.idx must be of the same length")
    }

    lv.grid <- precomp$lv.grid
    n.ord <- nrow(lv.grid)
    n.pgrid <- ncol(lv.grid)
    n.pval <- length(pval)
    if (any(ord.idx < 1) || any(ord.idx > n.ord)) {
        stop("Illegal order statistic index")
    }

    ## I assume the grid goes to a max pvalue of 1, so no p-value
    ## can be too big for the grid, in principle. If a p-value is
    ## too big, just replace it with the largest grid p-value minus
    ## a small amount to keep it off the boundary
    too.big <- which(pval > max(precomp$pvals) - .Machine$double.neg.eps)
    if (length(too.big) > 0) {
        pval[too.big] <- max(precomp$pvals) - .Machine$double.neg.eps
    }

    ## Extract out p-values that are too small for the grid. We will compute
    ## their l-values directly. Remaining p-values will get their l-values
    ## from the grid.
    ## A p-value of zero can be problematic. Replace these with a super
    ## small p-value.
    pval <- ifelse(pval < .Machine$double.xmin, .Machine$double.xmin, pval)
    too.small <- which(pval < precomp$pvals[1])
    if (length(too.small) > 0) {
        small.p <- pval[too.small]
        small.ord <- ord.idx[too.small]
        not.too.small <- (1:n.pval)[-too.small]
        p.grid <- pval[not.too.small]
        ord.grid <- ord.idx[not.too.small]

        lval.small <- lvalues.single(small.p, small.ord, precomp$prep)$l.values
    } else {
        not.too.small <- 1:n.pval
        p.grid <- pval
        ord.grid <- ord.idx
    }

    if (length(p.grid) > 0) {
        ## Get the p-value grid points that bound pval
        ## idx <- sapply(p.grid,
        ##               function(p) which.min(p > precomp$pvals) - 1)
        idx <- findInterval(p.grid, precomp$pvals)
        idxp1 <- idx + 1
        hi <- precomp$pvals[idx]
        hip1 <- precomp$pvals[idxp1]
        ## Lineaerly interpolate the l-value
        ## Taking the diagonal selects the elements (ord.idx[1], idx[1]),
        ## (ord.idx[2], idx[2]), etc.
        if (length(ord.grid) > 1) {
            ## lvi <- diag(lv.grid[ord.grid, idx])
            ## lvip1 <- diag(lv.grid[ord.grid, idxp1])
            vi <- (idx - 1) * n.ord + ord.grid
            vip1 <- (idxp1 - 1) * n.ord + ord.grid
            lvi <- lv.grid[vi]
            lvip1 <- lv.grid[vip1]
        } else {
            lvi <- lv.grid[ord.grid, idx]
            lvip1 <- lv.grid[ord.grid, idxp1]
        }
        lval.grid <- (p.grid - hi) * (lvip1 - lvi) / (hip1 - hi) + lvi
    }

    lval <- rep(NA, n.pval)
    if (length(not.too.small > 0)) {
        lval[not.too.small] <- lval.grid
    }
    if (length(too.small) > 0) {
        lval[too.small] <- lval.small
    }

    lval
}

shrink.C <- function(ycor, ny, cor.eig = NULL, tune = 0, ntest = 4) {
	  ## Given a correlation matrix made from a sample of Y's, apply shrinkage
    ## to get an estimate of the (true) population correlation matrix. The
    ## measure for accuracy is how closely the estimated population correlation
    ## gives a sample correlation matrix with the same mean r^2 (i.e. mean of
    ## the off diagonal elements) as the input correlation matrix. To do this
    ## the sample correlation matrix is shrunk towards the identity and
    ## ntest data sets, each with ny samples, are generated from which
    ## sample correlation matrices are computed. We get the mean r^2 for each
    ## of the ntest sample correlation matrices, and then the mean of these
    ## means. If this overall mean is larger than the mean r^2 of the original
    ## sample correlation matrix, set the tune parameter lower (more negative)
    ## and try again. If the overall mean is too small, set the tune parameter
    ## larger. Note that the method always picks a shrinkage parameter based
    ## on the largest eigenvalue of the original correlation matrix, the
    ## dimension of ycor, and the number of samples ny, so tune is an
    ## adjustment from this value.
    ##
    ## Input:
    ## ycor, sample correlation matrix
    ## ny, number of samples used to create ycor
    ## tune, parameter to tune the amount of shrinkage, negative values shrink
    ##       more towards the identity, positive values toward the sample
    ##       correlation matrix
    ## cor.eig, eigen decomposition of ycor
    ## ntest, number of times to sample a new correlation matrix to see if
    ##        the generated sample correlation matrices are similar to ycor
    ##
    ## Output:
    ## est.cor: The new shrunken correlation matrix
    ## est.r2: Vector of mean r^2 across ntest samples
    ## est.mean.r2: The mean of the values in est.r2
    ## y.mean.r2: The mean r^2 in the original sample correlation matrix ycor
    ## weight: The final tuning parameter. A value of one would give ycor while
    ##         a value of zero would give the identity matrix

    if (is.null(cor.eig)) {
        cor.eig <- eigen(ycor, symmetric = TRUE)
    }
    eigv1 <- cor.eig$values[1]
    p <- nrow(ycor)
    gam <- p / ny
    ## tune.simple comes from what the shrinkage parameter would have to be
    ## such that the largest shrunken eigenvalue is equal to the true largest
    ## eigenvalue when the true largest eigenvalue is much larger than one.
    ## When the true largest eigenvalue is big the sample eigenvalue is
    ## approximately the true eigenvalue plus gamma.
    tune.simple <- (eigv1 - 1 - gam) / (eigv1 - 1)
    tune.adj <- tune.simple + tune
    if (tune.adj < 0 || tune.adj > 1) {
        stop(paste0("tune must be larger than -", tune.simple, " and smaller than ",
                    1 - tune.simple))
    }

    shrunk.cor <- tune.adj * ycor
    diag(shrunk.cor) <- 1
    shrunk.eig <- eigen(shrunk.cor, symmetric = TRUE)

    mean.r2 <- numeric(ntest)
    for (i in 1:ntest) {
        y.samp <- rmvnorm.maa(ny, shrunk.eig)
        y.samp.cor <- crossprod(scale(y.samp)) / (ny - 1)
        s.offd <- trimat2vec(y.samp.cor)
        mean.r2[i] <- mean(s.offd^2)
    }
    y.offd <- trimat2vec(ycor)
    m.r2 <- mean(y.offd^2)
    mean.mean <- mean(mean.r2)
    ## sem <- sd(mean.r2) / sqrt(ntest)
    ## z.m <- (mean.mean - m.r2) / sem

    list(est.cor = shrunk.cor,
         est.eig = shrunk.eig,
         est.r2 = mean.r2,
         est.mean.r2 = mean.mean,
         y.mean.r2 = m.r2,
         weight = tune.adj)
}

make.emp.null <- function(M, Sig.z, precomp, block.size = 1000, verbose = TRUE,
                          n.cores = parallel::detectCores()/2) {
    ## Create an empiric distribution of ELL statistics. This is done by
    ## generating M multivariate normal vectors with covariance matrix
    ## given by Sig.z. Sig.z should be a correlation matrix for the
    ## generated vector to be a vector of Z scores, though this is not
    ## checked. For each generated vector, compute the l-values for the
    ## elements and retain the smallest l-value and which order statistic
    ## it came from. Set the seed before calling this function.
    ##
    ## Input:
    ## M: Number of replicates
    ## Sig.z: Covariance matrix for the Z scores, size D x D
    ## precomp: Output from doing the precomputation of the l-value grid
    ## block.size: How many replicates to do at one time. A larger number
    ##             wil be somewhat faster, but take up more memory
    ## verbose: Logical. If TRUE print a message each time a block is done.
    ## n.cores: Number of cores to use when computing the statistic across
    ##          replicates in each block. This does not affect the number of
    ##          cores used when generating the Z scores, which is determined
    ##          by the BLAS library used.
    ##
    ## Output is a list with these objects:
    ## adelle.dist: Vector of length M with sorted values of the ELL statistic
    ## ord.dist: Which order statistic each ELL statistic was. These are
    ##           sorted to be in the same order as in adelle.dist

    if (is(Sig.z, "eigen")) {
        if (nrow(Sig.z$vectors) != precomp$prep$N.tests) {
            stop("Number of rows of Sigma must match the number of tests.")
        }
        sig.eig <- Sig.z
    } else {
        if (nrow(Sig.z) != precomp$prep$N.tests) {
            stop("Number of rows of Sigma must match the number of tests.")
        }
        sig.eig <- eigen(Sig.z, symmetric = TRUE)
    }
    N.ord <- nrow(precomp$lv.grid)

    ## Do replicates in blocks of size block.size. This is to avoid making
    ## such a huge matrix of Z scores that it fills up memory
    if (block.size > M) {
        block.size <- M
    }
    n.blocks <- M %/% block.size
    n.rem <- M %% block.size

    adelle.dist <- rep(NA_real_, M)
    ord.dist <- rep(NA_real_, M)
    i.block <- 1
    first.idx <- 1
    while (i.block <= n.blocks) {
        min.dist <- comp.dist(block.size, sig.eig, N.ord, precomp, n.cores)
        adelle.dist[first.idx:(i.block * block.size)] <- min.dist[1,]
        ord.dist[first.idx:(i.block * block.size)] <- min.dist[2,]
        first.idx <- i.block * block.size + 1
        if (verbose) {
            cat("Done with block ", i.block, "\n")
        }
        i.block <- i.block + 1
    }

    if (n.rem > 0) {
        min.dist <- comp.dist(n.rem, sig.eig, N.ord, precomp, n.cores)
        adelle.dist[first.idx:(first.idx + n.rem - 1)] <- min.dist[1,]
        ord.dist[first.idx:(first.idx + n.rem - 1)] <- min.dist[2,]
    }

    sort.idx <- order(adelle.dist)
    list(adelle.dist = adelle.dist[sort.idx], ord.dist = ord.dist[sort.idx])
}

adelle.emp.pval <- function(lval, l.dist) {
    ## Get the p-value for the input l-value based on the empirical distribution
    ## given. Smaller l-values are considered more significant and the input
    ## empiric distribution is ordered from smallest to largest.
    ##
    ## Input:
    ## lval: The l-value to get the significance of. May be a vector.
    ## l.dist: Vector with ordered l-values from the null distribution
    ##
    ## Output:
    ## Vector of length equal to lval with the empiric p-value for each l-value

    L <- length(l.dist)
    M <- length(lval)

    l.order <- order(lval)
    l.sort <- lval[l.order]
    r <- which(order(c(l.sort, l.dist)) <= M) - 1:M
    pv.sort <- (r + 1) / (L + 1)
    pv.unsort <- pv.sort[order(l.order)]
    return(pv.unsort)
}

adelle.lval.thresh <- function(thr, l.dist) {
    ## Get the l-value that corresponds to a given significance threshold.
    ## That is, return the largest l-value that corresponds to a p-value
    ## that is less than or equal to the threshold specified. Use the given
    ## empiric distribution as the distribution of the l-values. Thresholds
    ## may not be larger than 1. If a threshold is <= 1/(L+1) return NA.
    ##
    ## Input:
    ## thr: A desired significance threshold. May be a vector of thresholds
    ## l.dist: Vector of ordered l-values with higher significance
    ##         corresponding to smaller l-values
    ##
    ## Output:
    ## Vector of length equal to thr with the l-value that corresponds to
    ## a p-value no larger than thr. If thr <= 1/(L+1) return NA

    if (any(thr > 1)) {
        stop("Thresholds may not be larger than 1")
    }

    L <- length(l.dist)

    r <- ceiling(thr * (L + 1) - 1)
    r <- ifelse(r <= 0, NA, r)
    l.dist[r]
}

########################################################################
##
## Various helper functions called by the main Adelle functions
##
##

precomp.gamma <- function(pval, lval.prep) {
    ## Given a vector of p-values for correlated tests, compute the parameter
    ## gamma for each of the p-values.
    ##
    ## Input:
    ## pval: Vector of p-values
    ## lval.prep: The output of prep.lvalue
    ##
    ## Output:
    ## gamma: The values of the gamma parameters for the N.ord smallest p-values
    ## smallp.idxs: Indices of pval that have the N.ord smallest p-values

    zt <- qnorm(0.5 * pval, lower.tail = FALSE)
    h.mat <- lval.prep$Herm.funcs(zt)^2
    c.val <- drop(4 * dnorm(zt)^2 * (h.mat %*% (lval.prep$r.means / lval.prep$fac.vals)) /
                  (pval * (1 - pval)))

    c.val / (1 - c.val) # = gamma
}

lv.pmf.mat <- function(pval, gam.p, N.tests) {
	  ## Make the matrix that has the pmf needed to compute the l-values.
    ## The columns are for each parameter value (pval and gam.p) and
    ## the rows for each possible value of the order statistic (from
    ## 1 to N.tests). By summing the correct entries in a column, you
    ## get the l-value, which is 1 - cdf. This algorithm takes care to
    ## check whether the binomial or beta-binomial should be used.
    ## Return the matrix.

    M.mat <- matrix(NA_real_, nrow = N.tests, ncol = length(pval))
    ## When gamma is zero, or close enough to it, we use a binomial rather
    ## than the beta-binomial to compute the l-value
    binom.idx <- which(gam.p < 10 * .Machine$double.xmin)
    betab.idx <- which(gam.p >= 10 * .Machine$double.xmin)

    ## binom.M.mat has each column the binomial pmf for that column's probability
    ## with each row the possible values of the order statistic index (1 to N.tests)
    if (length(binom.idx > 0)) {
        binom.p <- pval[binom.idx]
        M.mat[, binom.idx] <- sapply(1:length(binom.p),
                                     function(j) dbinom(1:N.tests, N.tests, binom.p[j]))
    }

    ## Now get the matrix with the beta-binomial pmfs
    if (length(betab.idx) > 0) {
        betab.p <- pval[betab.idx]
        bb.gam <- gam.p[betab.idx]
        log.denom <- lbeta(betab.p / bb.gam, (1 - betab.p) / bb.gam)
        log.choose <- lchoose(N.tests, 1:N.tests)
        ## log.M.mat <- matrix(NA_real_, nrow = N.tests, ncol = N.ord)
        bb.log.M.mat <- sapply(1:length(betab.p),
                               function(j)
                                   lbeta(betab.p[j] / bb.gam[j] + 1:N.tests,
                                   (1 - betab.p[j]) / bb.gam[j] + N.tests - 1:N.tests)
                               )
        ## log.M.mat has one column for each parameter value pair (pval and gam.p)
        ## With the rows corresponding to the index of the order statistic

        M.mat[, betab.idx] <- exp(t(t(bb.log.M.mat) - log.denom) + log.choose)
    }
    M.mat
}

lvalues.single <- function(pval, ord.idx, lval.prep) {
    ## Compute the l-values for each of the p-values, order index pair in
    ## pval, ord.idx. That is, l-values are computed for the pairs
    ## (pval_1, ord.idx_1), (pval_2, ord.idx_2), ....
    ## The output is a vector of the l-values for each of the input pairs.
    ##
    ## Input:
    ## pval: vector of p-values
    ## ord.idx: vector of order statistic indices, one for each pval
    ## lval.prep: Output of prep.lvalue using the correlation matrix of Z's
    ##
    ## Output:
    ## l.values: Vector of l-values for each of the p-values
    ## gammas: The gamma parameter for the beta-binomial for each p-value

    N.tests <- lval.prep$N.tests
    if (length(pval) != length(ord.idx)) {
        stop("Length of pval must equal the length of ord.idx.")
    }

    gam.p <- precomp.gamma(pval, lval.prep)

    M.mat <- lv.pmf.mat(pval, gam.p, N.tests)
    lv.mat <- as.matrix(apply(M.mat, 2, function(x) rev(cumsum(rev(x)))[1:max(ord.idx)]))

    if (length(ord.idx) > 1) {
        l.values <- diag(lv.mat[ord.idx, ])
    } else {
        l.values <- lv.mat[ord.idx, ]
    }

    list(l.values = l.values, gammas = gam.p)
}

trimat2vec <- function(m, diagonal = FALSE) {
    ## m is a lower triangular (or symmetric) matrix
    ## diag = keep diagonal (TRUE) or not (FALSE)
    ## This gives the same output as m[lower.tri(m)] but faster
    ## for large matrices (nrow = a few thousand)

    n <- dim(m)[1]
    if (diagonal) {
        v <- vector(mode="numeric", 0.5*n*(n+1))

        begin <- 1
        for ( i in 1:n) {
            end <- begin + n - i
            v[ begin : end ] <- m[ i:n , i ]
            begin <- end + 1
        }
    } else {
        v <- vector(mode="numeric", 0.5 * n * (n - 1))

        begin <- 1
        for ( i in 1:(n-1)) {
            end <- begin + n - i - 1
            v[ begin : end ] <- m[ (i + 1):n , i ]
            begin <- end + 1
        }
    }
    v
}

prep.lvalue <- function(cor.mat,
                        max.terms = 20,
                        N.tests = nrow(cor.mat),
                        n.cores = parallel::detectCores()/2) {
    ## Prepare quantities need to compute the l-value for a SNP.
    ## cor.mat: Correlation matrix of the tests, or a single value
    ##          equal to the correlation of all non-diagonal elements
    ##          in the matrix.
    ## max.terms: The number of terms used in the sum over Hermite polynomials
    ## N.tests: Number of tests (i.e. traits). If cor.mat is a single number then
    ##    an integer equal to the number of rows in the correlation matrix
    ## n.cores: Number of cores to use
    ##
    ## Output:
    ## A list with functions for the Hermite polynomials, means of the off
    ## diagonals of the correlation matrix raised to an even power, and
    ## other quantities needed to compute the Hermite polynomial representation
    ## of a multivariate Gaussian.

    if (is.null(N.tests)) {
        stop("Enter number of tests for N.tests.")
    }

    sum.idxs <- 2 * (1:max.terms)
    herm.idxs <- sum.idxs - 1
    fac.vals <- factorial(sum.idxs)

    Herm.funcs <- as.function(hermite(herm.idxs))
    ## Herm.funcs is a function that takes a single numeric input
    ## and returns a vector with the value of all the Hermite
    ## polynomials (probabilist's version) from 1 to (2 * max.terms - 1)
    ## computed at that value.

    if (length(cor.mat) == 1) {
        cor.vec <- cor.mat
        r.means <- sapply(sum.idxs, function(x) cor.vec^x)
    } else {
        ## This extracts the upper triangle (no diagonal) of the matrix
        ## into a vector
        ## cor.vec <- cor.mat[col(cor.mat) > row(cor.mat)]
        cor.vec <- trimat2vec(cor.mat) ## A bit faster than the above line
        ## r.means <- sapply(sum.idxs, function(x) mean(cor.vec^x))
        r.means <- unlist(parallel::mclapply(sum.idxs, function(x) mean(cor.vec^x),
                                             mc.cores = n.cores))
    }

    list(#cor.mat = cor.mat,
        sum.idxs = sum.idxs,
        herm.idxs = herm.idxs,
        fac.vals = fac.vals,
        r.means = r.means,
        Herm.funcs = Herm.funcs,
        N.tests = N.tests)
}

rmvnorm.maa <- function(M, Sig.z, eps = 1e-12) {
    ## Generate a matrix of random multivariate normal vectors. The mean
    ## is zero and covariance matrix is given by Sig.z. Set the seed before
    ## calling this function.
    ##
    ## Input:
    ## M - Number of replicates
    ## Sig.z - Covariance matrix of dim D x D, D is # of traits. Can
    ##         input the eigen decomposition of Sig.z instead.
    ## eps - If there are negative eigenvalues the absolute value of the
    ##       smallest must be greater than eps * (largest eigenvalue)
    ##
    ## Output:
    ## An M x D matrix with each row a replicate from the MVN distribution.


    if (is(Sig.z, "eigen")) {
        sig.eig <- Sig.z
    } else {
        sig.eig <- eigen(Sig.z, symmetric = TRUE)
    }

    D <- length(sig.eig$values)
    if (sig.eig$values[D] < 0) {
        if (abs(sig.eig$values[D] / sig.eig$values[1]) > eps) {
            stop("Sig.z matrix is not sufficiently close to positive definite.")
        }
        sig.eig$values <- sig.eig$values + 2 * abs(sig.eig$values[D])
    }

    ## L.mat <- sig.eig$vectors %*% (sqrt(sig.eig$values) * t(sig.eig$vectors))
    L.mat <- sqrt(sig.eig$values) * t(sig.eig$vectors)
    ZZ <- matrix(rnorm(D * M), nrow = M, ncol = D, byrow = TRUE)
    Z <- ZZ %*% L.mat
    Z
}

comp.dist <- function(n.repl, sig.eig, N.ord, precomp,
                      n.cores = parallel::detectCores()/2) {
	  ## Compute the ELL null distribution for n.repl replicates. Output
    ## a matrix with two rows, the minimum l-value for each replicate in
    ## row one, and which order statistic gave the minimum l-value in
    ## row two.
    Z <- rmvnorm.maa(n.repl, sig.eig)
    Z <- abs(Z)
    pmat <- t(2 * pnorm(Z, lower.tail = FALSE))
    ## pmat now has each column as a replicate rather than rows
    ## So, we want the minimum l-value from each column and which
    ## order statistic it came from.
    ## min.dist <- sapply(1:ncol(pmat), function(i) { # non-parallel version
    min.dist <- parallel::mclapply(1:ncol(pmat), function(i) {
        p.sort <- sort(pmat[, i])[1:N.ord]
        l.rep <- get.lvalue(p.sort, 1:N.ord, precomp)
        which.ord <- which.min(l.rep)
        min.l <- l.rep[which.ord]
        c(min.l, which.ord)
    }, mc.cores = n.cores)
    ## Next line necessary for mclapply
    min.dist <- matrix(unlist(min.dist), nrow = 2, ncol = length(min.dist))
    ## First row is min l-value, second row is which order stat
    min.dist
}
