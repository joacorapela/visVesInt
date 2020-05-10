fit_FA_MARSS_model <- function(firingRatesDataBlock, stateInputs, obsInputs, dimLat, stateInputMemorySecs, obsInputMemorySecs, sRate, maxIter, silentLevel=2, controlFA=list(trace=TRUE, nstart=5)) {


    dimObs <- nrow(firingRatesDataBlock)
    nObs <- ncol(firingRatesDataBlock)

    # begin estimate initial conditions
    firingRatesDataBlockForFA <- t(as.matrix(firingRatesDataBlock))
    initialConds <- estimateKFInitialCondFA(z=firingRatesDataBlockForFA, nFactors=dimLat, control=controlFA)
    # end estimate initial conditions

    # begin create model
    B1List <- c()
    for(j in 1:dimLat) {
        for(i in 1:dimLat) {
            B1List <- c(B1List, list(sprintf("b%d%d", i, j)))
        }
    }
    B1  <- matrix(B1List, nrow=dimLat)
    U1  <- "unequal"
    Q1  <- "diagonal and equal"
    Z1List <- c()
    for(j in 1:dimLat) {
        for(i in 1:dimObs) {
            Z1List <- c(Z1List, list(sprintf("z%d%d", i, j)))
        }
    }
    Z1  <- matrix(Z1List, nrow=dimObs)
    if(!is.nan(stateInputs)) {
        C1  <- "unconstrained"
        c1 <- stateInputs
    } else {
        C1  <- "zero"
        c1 <- NA
    }
    A1  <- "unequal"
    R1  <- "diagonal and unequal"
    pi1 <- "unequal"
    V01 <- "diagonal and equal"
    if(!is.nan(obsInputs)) {
        D1  <- "unconstrained"
        d1 <- obsInputs
    } else {
        D1  <- "zero"
        d1 <- NA
    }

    model.list <- list(B=B1, U=U1, C=C1, c=c1, Q=Q1, Z=Z1, A=A1, D=D1, d=d1, R=R1, x0=pi1, V0=V01)
    # end create model

    # begin set initial conditions
    B0 <- matrix(as.vector(initialConds$A), ncol=1)
    Z0 <- matrix(as.vector(initialConds$C), ncol=1)
    R0 <- matrix(initialConds$sigmaDiag, ncol=1)
    control <- list(maxit=maxIter)

    inits <- list(B=B0, Z=Z0, R=R0)
    # end set initial conditions

    kem <- MARSS(as.matrix(firingRatesDataBlock), model=model.list, inits=inits, control=control, silent=silentLevel)
    return(kem)
}
