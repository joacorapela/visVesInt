fit_FA_MARSSmodel <- function(dataBlock) {
    if(!is.nan(stateInputMemorySecs)) {
        stateVelInputs <- buildVelocityInputs(velocities=as.numeric(dataBlock[1,]), inputMemorySecs=stateInputMemorySecs, sRate=sRate)
    } else {
        stateVelInputs <- "zero"
    }
    if(!is.nan(obsInputMemorySecs)) {
        obsVelInputs <- buildVelocityInputs(velocities=as.numeric(dataBlock[1,]), inputMemorySecs=obsInputMemorySecs, sRate=sRate)
    } else {
        obsVelInputs <- "zero"
    }

    dimLat <- nFactors
    dimObs <- nrow(zs)
    nObs <- ncol(zs)

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
    if(!is.nan(stateInputMemorySecs)) {
        C1  <- "unconstrained"
    } else {
        C1  <- "zero"
    }
    A1  <- "unequal"
    R1  <- "diagonal and unequal"
    pi1 <- "unequal"
    V01 <- "diagonal and equal"
    if(!is.nan(obsInputMemorySecs)) {
        D1  <- "unconstrained"
    } else {
        D1  <- "zero"
    }

    model.list <- list(B=B1, U=U1, C=C1, c=stateVelInputs, Q=Q1, Z=Z1, A=A1, D=D1, d=obsVelInputs, R=R1, x0=pi1, V0=V01)
    # end create model

    # begin set initial conditions
    zsForFA <- t(as.matrix(zs))
    initialConds <- estimateKFInitialCondFA(z=zsForFA, nFactors=nFactors)

    B0 <- matrix(as.vector(initialConds$A), ncol=1)
    Z0 <- matrix(as.vector(initialConds$C), ncol=1)
    R0 <- matrix(initialConds$sigmaDiag, ncol=1)
    control <- list(maxit=maxIter)

    inits <- list(B=B0, Z=Z0, R=R0)
    # end set initial conditions

    kem <- MARSS(as.matrix(zs), model=model.list, inits=inits, control=control, silent=2)
    kfRes <- MARSSkf(kem)

    kem0 <- kem
    kem0$par <- kem0$start
    kfRes0 <- MARSSkf(kem0)

    resultsFilenamePattern <- "results/expt%02d_MARSS_condition%s_nLat%02d_stateInputMemory%.02f_obsInputMemory%.02f.RData"
    resultsFilename <- sprintf(resultsFilenamePattern, exptN, conditionNoBlanks, nFactors, stateInputMemorySecs, obsInputMemorySecs)
    results <- list(kem=kem, kfRes=kfRes, kfRes0=kfRes0, zs=zs, trialStartTimes=trialStartTimes)
}
