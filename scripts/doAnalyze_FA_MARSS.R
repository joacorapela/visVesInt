
require(MASS)
require(MARSS)
require(ggplot2)
require(plotly)
require(reshape2)
source("~/dev/learning/kalmanFilter/code/src/plotTrueInitialAndEstimatedMatrices.R")
source("~/dev/learning/kalmanFilter/code/src/plotTrueInitialAndEstimatedVectors.R")
source("~/dev/learning/kalmanFilter/code/src/estimateKFInitialCondFA.R")

source("getDataSubset.R")
source("buildExpCellByTimeConcatenatedTrialsDataBlock.R")
source("circleFun.R")

buildVelocityInputs <- function(velocities, inputMemorySecs, sRate) {
    inputMemorySamples <- round(inputMemorySecs*sRate)
    velInputs <- c()
    N <- length(velocities)
    for(i in 0:inputMemorySamples) {
        velInputs <- rbind(velInputs, c(rep(0, time=i), velocities[1:(N-i)]))
    }
    return(velInputs)
}

processAll <- function() {
    nFactors <- 6
    sRate <- 10
    stateInputMemorySecs <- NaN
    obsInputMemorySecs <- 0.0
    maxIter <- 100
    # exptN <- NA
    exptN <- 6
    cellN <- NA
    trialN <- NA
    # condition <- NA
    # condition <- "Visual + Vestibular"
    condition <- "Visual"
    # condition <- "Vestibular"
    region <- NA
    # region <- c("RSPg")
    layer <- NA
    # layer <- "Layer 2/3"
    cellType <- NA
    # cellType <- "Wide"
    sdGamma <- 1000.0
    sdXHat00 <- 1000.0
    sdSigmaX00 <- 1000.0

    dataFilename <- "../../data/011620/Spike_Freq_Dataframe_Hugo.csv"

    conditionNoBlanks <- gsub(" ", "_", condition)
    conditionNoBlanks <- gsub("_+_", "_", conditionNoBlanks)

    data <- read.csv(dataFilename)
    dataSubset <- getDataSubset(data=data, exptN=exptN, cellN=cellN, trialN=trialN, condition=condition, region=region, layer=layer, cellType=cellType)
    # dataSubsetAvgTrials <- aggregate(x=data.frame("AvgFiringRate"=dataSubset$Firing.Rate), by=dataSubset[, c("Expt..", "Cell..", "Bin..", "Position..deg.", "Speed..deg.s.", "Trial.Condition", "Region", "Layer", "Cell.Type")], FUN=mean)
    # zs <- buildExpTrialBySpeedFiringRateDataBlock(data=dataSubsetAvgTrials, firingRateColLabel="AvgFiringRate")
    res <- buildExpCellByTimeConcatenatedTrialsDataBlock(data=dataSubset, dataColLabel="Firing.Rate")
    trialStartTimes <- cumsum(res$trialNSamples)/sRate
    zs <- res$dataBlock
    zs <- sqrt(zs)
    res <- buildExpCellByTimeConcatenatedTrialsDataBlock(data=dataSubset, dataColLabel="Speed..deg.s.")
    if(!is.nan(stateInputMemorySecs)) {
        stateVelInputs <- buildVelocityInputs(velocities=as.numeric(res$dataBlock[1,]), inputMemorySecs=stateInputMemorySecs, sRate=sRate)
    } else {
        stateVelInputs <- "zero"
    }
    if(!is.nan(obsInputMemorySecs)) {
        obsVelInputs <- buildVelocityInputs(velocities=as.numeric(res$dataBlock[1,]), inputMemorySecs=obsInputMemorySecs, sRate=sRate)
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
    save(results, file=resultsFilename)

    figFilenamePattern <- "figures//expt%02d//exp%02d_condition%s_nLat%02d_stateInputMemory%.02f_obsInputMemory%.02f_MARSS_%s.html"

    AFigFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, nFactors, stateInputMemorySecs, obsInputMemorySecs, "A")
    plotTrueInitialAndEstimatedMatrices(initial=initialConds$A, estimated=matrix(coef(kem)$B, nrow=dimLat), title="A", figFilename=AFigFilename)
    # begin plot params
    uFigFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, nFactors, stateInputMemorySecs, obsInputMemorySecs, "u")
    plotTrueInitialAndEstimatedVectors(estimated=matrix(coef(kem)$U, nrow=dimLat), title="u", figFilename=uFigFilename)
    #
    CFigFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, nFactors, stateInputMemorySecs, obsInputMemorySecs, "C")
    plotTrueInitialAndEstimatedMatrices(initial=initialConds$C, estimated=matrix(coef(kem)$Z, nrow=dimObs), title="C", figFilename=CFigFilename)
    #
    aFigFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, nFactors, stateInputMemorySecs, obsInputMemorySecs, "a")
    plotTrueInitialAndEstimatedVectors(estimated=matrix(coef(kem)$A, nrow=dimObs), title="a", figFilename=aFigFilename)
    #
    GammaFigFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, nFactors, stateInputMemorySecs, obsInputMemorySecs, "Gamma")
    plotTrueInitialAndEstimatedVectors(estimated=coef(kem)$Q, title="Gamma Variance", figFilename=GammaFigFilename)
    #
    SigmaFigFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, nFactors, stateInputMemorySecs, obsInputMemorySecs, "Sigma")
    plotTrueInitialAndEstimatedVectors(initial=initialConds$sigmaDiag, estimated=coef(kem)$R, title="Sigma Diagonal", figFilename=SigmaFigFilename)
    #
    x0FigFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, nFactors, stateInputMemorySecs, obsInputMemorySecs, "x0")
    plotTrueInitialAndEstimatedVectors(estimated=coef(kem)$x0, title="x0", figFilename=x0FigFilename)
    #
    V0FigFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, nFactors, stateInputMemorySecs, obsInputMemorySecs, "V00")
    plotTrueInitialAndEstimatedVectors(estimated=coef(kem)$x0, title="V00", figFilename=x0FigFilename)
    #
    if(!is.nan(stateInputMemorySecs)) {
        CInputFigFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, nFactors, stateInputMemorySecs, obsInputMemorySecs, "CInput")
        plotTrueInitialAndEstimatedMatrices(estimated=t(matrix(coef(kem)$C), nrow=dimObs), title="C Input", figFilename=CInputFigFilename)
    }
    if(!is.nan(obsInputMemorySecs)) {
        DInputFigFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, nFactors, stateInputMemorySecs, obsInputMemorySecs, "DInput")
        plotTrueInitialAndEstimatedMatrices(estimated=t(matrix(coef(kem)$D, nrow=dimObs)), title="D Input", figFilename=DInputFigFilename)
    }
    # end plot params

    # begin plot state dynamics spectrum
    A <- matrix(kem$start$B, ncol=nFactors)
    dfCircle <- circleFun(c(0,0), radius=1, npoints = 200)
    #geom_path will do open circles, geom_polygon will do filled circles
    eigvals <- eigen(x=A, symmetric=FALSE, only.values=TRUE)$values
    dfEigvals <- data.frame(x=Re(eigvals), y=Im(eigvals))
    p <- ggplot()
    p <- p + geom_path(data=dfCircle, mapping=aes(x=x, y=y), color="gray")
    p <- p + geom_point(data=dfEigvals, mapping=aes(x=x, y=y, color="eigenvalue"))
    p <- p + geom_vline(xintercept=0, linetype="solid")
    p <- p + geom_hline(yintercept=0, linetype="solid")
    p <- p + xlab("Real")
    p <- p + ylab("Imaginary")
    p <- p + theme(legend.title = element_blank())
    p <- ggplotly(p)
    AeigenFigFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, nFactors, stateInputMemorySecs, obsInputMemorySecs, "Aeigen")
    htmlwidgets::saveWidget(as_widget(p), file.path(normalizePath(dirname(AeigenFigFilename)), basename(AeigenFigFilename)))
    # end plot state dynamics spectrum
    #
    # begin plot predicted observations
    #
    # begin compute one-lag ahead observation predictions stats
    Z <- matrix(coef(kem)$Z, nrow=dimObs)
    D <- matrix(coef(kem)$D, nrow=dimObs)
    if(!is.nan(obsInputMemorySecs)) {
        ytt1 <- Z%*%kfRes$xtt1+as.numeric(coef(kem)$A)+D%*%obsVelInputs
    } else {
        ytt1 <- Z%*%kfRes$xtt1+as.numeric(coef(kem)$A)
    }
    Wtt1 <- array(NA, dim=c(dimObs, dimObs, nObs))
    if(kem$call$model$R=="diagonal and unequal") {
        R <- diag(coef(kem)$R)
    } else {
        error("Functionality not yet implemented for R!=<diagonal and unequal>")
    }
    for(n in 1:nObs) {
        Wtt1[,,n] <- Z%*%kfRes$Vtt1[,,n]%*%t(Z)+R
    }
    # end compute one-lag ahead observation predictions stats
    #
    for(i in 1:dimObs) {
        df <- data.frame(time=(1:ncol(ytt1))/sRate, predMean=ytt1[i,], predMeanLower=ytt1[i,]-1.96*sqrt(Wtt1[i,i,]), predMeanUpper=ytt1[i,]+1.96*sqrt(Wtt1[i,i,]), observation=as.numeric(zs[i,]))
        p <- ggplot(data=df, mapping=aes(x=time))
        p <- p + geom_line(aes(y=predMean, col="prediction"))
        p <- p + geom_point(aes(y=observation, col="observation"))
        p <- p + geom_ribbon(mapping=aes(ymin=predMeanLower, ymax=predMeanUpper), alpha=0.3, fill="red")
        p <- p + scale_colour_manual(name="", values=c("prediction"="red", "observation"="black"))
        p <- p + geom_vline(xintercept=trialStartTimes, linetype="dotted")
        p <- p + xlab("Time (sec)")
        p <- p + ylab("Square-Root Transformed Firing Rate")
        p <- p + ggtitle(rownames(zs)[i])
        p <- p + theme(legend.title = element_blank())
        p <- ggplotly(p)
        suffix <- sprintf("pred%s", rownames(zs)[i])
        predFigFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, nFactors, stateInputMemorySecs, obsInputMemorySecs, suffix)
        htmlwidgets::saveWidget(as_widget(p), file.path(normalizePath(dirname(predFigFilename)), basename(predFigFilename)))
        # print(p)
    }
    # end plot predicted observations
    #
    # begin compute nmses
    nmse <- array(rep(NA, dimObs))
    rownames(nmse) <- rownames(zs)
    for(i in 1:dimObs) {
        zsi <- as.numeric(zs[i,])
        nmse[i] <- mean((ytt1[i,]-zsi)^2)/var(zsi)
    }
    # end compute nmses
    #
    # begin plot nmses
    sortRes <- sort(nmse, index.return=TRUE)
    sortedNMSE <- sortRes$x
    sortedCellNames <- factor(rownames(nmse)[sortRes$ix], levels=rownames(nmse)[sortRes$ix])
    df <- data.frame(x=sortedCellNames, y=sortedNMSE)
    p <- ggplot(data=df, aes(x=x, y=y))
    p <- p + geom_bar(stat="identity")
    p <- p + ylab("NMSE")
    p <- p + xlab("")
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    p <- ggplotly(p)
    nmseFigFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, nFactors, stateInputMemorySecs, obsInputMemorySecs, "NMSE")
    htmlwidgets::saveWidget(as_widget(p), file.path(normalizePath(dirname(nmseFigFilename)), basename(nmseFigFilename)))
    # end plot nmses
    #
    # begin plot latents
    alphaRibbon <- 0.3
    df <- data.frame()
    for(i in 1:nrow(kfRes$xtT)) {
        dataBlock <- data.frame(sample=1:length(kfRes$xtT[i,]),
                                mean=kfRes$xtT[i,],
                                sd=sqrt(kfRes$VtT[i,i,]),
                                latentID=rep(i, length(kfRes$xtT[i,])),
                                latentType=rep("estimated",
                                               length(kfRes$xtT[i,])))
        df <- rbind(df, dataBlock)
    }
    for(i in 1:nrow(kfRes0$xtT)) {
        dataBlock <- data.frame(sample=1:length(kfRes0$xtT[i,]),
                                mean=kfRes0$xtT[i,],
                                sd=sqrt(kfRes0$VtT[i,i,]),
                                latentID=rep(i, length(kfRes0$xtT[i,])),
                                latentType=rep("initial",
                                               length(kfRes0$xtT[i,])))
        df <- rbind(df, dataBlock)
    }
    p <- ggplot(df, aes(x=sample/sRate, y=mean,
                        ymin=mean-1.96*sd,
                        ymax=mean+1.96*sd,
                        color=factor(latentID),
                        fill=factor(latentID),
                        linetype=factor(latentType)))
    p <- p + geom_line()
    p <- p + geom_ribbon(alpha=0.3)
    p <- p + geom_hline(yintercept=0)
    p <- p + geom_vline(xintercept=0)
    p <- p + geom_vline(xintercept=trialStartTimes, linetype="dotted")
    p <- p + ylab("Kalman Smoother State Estimate")
    p <- p + xlab("Time")
    p <- p + theme(legend.title = element_blank())
    p <- ggplotly(p)
    latentsFigFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, nFactors, stateInputMemorySecs, obsInputMemorySecs, "Latents")
    htmlwidgets::saveWidget(as_widget(p), file.path(normalizePath(dirname(latentsFigFilename)), basename(latentsFigFilename)))
    # print(p)
    # end plot latents

    browser()
}

processAll()
