
processAll <- function() {
    A <- matrix(c(-128, 80, -272, 128), ncol=2)
    dt <- 1e-3
    x0 <- c(1,1)
    nSteps <- 800

    AforI <- dt*A + diag(2)
    x <- matrix(NA, nrow=nrow(A), ncol=nSteps)
    x[,1] <- x0
    for(n in 2:nSteps) {
        x[,n] <-AforI%*%x[,n-1]
    }
    plot(x[1,], x[2,])
    browser()
}

processAll()
