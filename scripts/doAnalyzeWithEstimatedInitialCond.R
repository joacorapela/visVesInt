
require(MASS)
require(ggplot2)
require(mvtnorm)
require(ramcmc)                                                       
require(gridExtra)

source("~/dev/learning/kalmanFilter/code/src/lsolve.R")
source("~/dev/learning/kalmanFilter/code/src/chol_downdate_higherOrder.R")
source("~/dev/learning/kalmanFilter/code/src/squareRootKF.R")
source("~/dev/learning/kalmanFilter/code/src/smoothLDS.R")
source("~/dev/learning/kalmanFilter/code/src/emEstimationSquareRootKF.R")
source("~/dev/learning/kalmanFilter/code/src/estimateKFInitialCond.R")
source("getDataSubset.R")
source("buildExpTrialBySpeedFiringRateDataBlock.R")

processAll <- function() {
    nFactors <- 25
    exptN <- NA
    cellN <- NA
    trialN <- NA
    condition <- "Visual + Vestibular"
    # condition <- "Visual"
    # condition <- "Vestibular"
    coordinatesToPlot <- c(1, 5)
    # region <- "RSPg"
    region <- c("RSPg")
    # layer <- NA
    layer <- "Layer 2/3"
    cellType <- "Wide"
    sdXHat00 <- 1000.0
    sdSigmaX00 <- 1000.0

    sRate <- 10
    nIter <- 50
    dataFilename <- "../../data/011620/Spike_Freq_Dataframe_Hugo.csv"
    figFilename <- "figures/kfEstimates.png"

    data <- read.csv(dataFilename)
    dataSubset <- getDataSubset(data=data, exptN=exptN, cellN=cellN, trialN=trialN, condition=condition, region=region, layer=layer, cellType=cellType)
    dataSubsetAvgTrials <- aggregate(x=data.frame("AvgFiringRate"=dataSubset$Firing.Rate), by=dataSubset[, c("Expt..", "Cell..", "Bin..", "Position..deg.", "Speed..deg.s.", "Trial.Condition", "Region", "Layer", "Cell.Type")], FUN=mean)
    zs <- buildExpTrialBySpeedFiringRateDataBlock(data=dataSubsetAvgTrials, firingRateColLabel="AvgFiringRate")
    zsForFA <- as.matrix(zs)
    # initialCond <- estimateKFInitialCond(z=zsForFA, nFactors=nFactors)
    browser()
    initialCond <- estimateKFInitialCond(z=zsForFA, nFactors=nFactors)
    nObs <- ncol(zs)
    A0 <- initialCond$A
    Gamma0 <- initialCond$Gamma
    SRSigmaW0 <- chol(Gamma0)
    C0 <- initialCond$C
    Sigma0 <- initialCond$Sigma
    D0 <- matrix(0, nrow=nrow(C0), ncol=1)
    Sigma0 <- initialCond$Sigma
    SRSigmaV0 <- chol(Sigma0)
    B0 <- matrix(0, nrow=nrow(A0), ncol=1)
    xHat00 <- rep(0, nrow(A0))
    SRSigmaX00 <- sdSigmaX00*diag(rep(1, length(xHat00)))
    nTimeSamples <- ncol(zs)
    us <- matrix(0, nrow=1, ncol=nTimeSamples)

    t <- (1:ncol(zs))/sRate
    eRes <- emEstimationSquareRootKF(zs=zs, A0=A0, SRSigmaW0=SRSigmaW0, C0=C0, SRSigmaV0=SRSigmaV0, B0=B0, D0=D0, xHat00=xHat00, SRSigmaX0=SRSigmaX00, nIter=nIter, varsToEstimate=list(initialStateMean=TRUE, initialStateCovariance=TRUE, transitionMatrix=TRUE, transitionCovariance=TRUE, observationMatrix=TRUE, observationCovariance=TRUE))

    df <- data.frame(x=1:length(eRes$logLikelihoods), 
                     y=eRes$logLikelihoods)
    p1 <- ggplot(df, aes(x=x, y=y))
    p1 <- p1 + geom_line()
    p1 <- p1 + geom_point()
    p1 <- p1 + xlab("Time")
    p1 <- p1 + ylab("Log Likelihood")

    fRes0 <- squareRootKF(A=A0, B=B0, C=C0, D=D0, xHat0=xHat00, SRSigmaX0=SRSigmaX00, SRSigmaW=SRSigmaW0, SRSigmaV=SRSigmaV0, us=us, zs=zs)
    sRes0 <- smoothLDS(A=A0, mu=fRes0$xHat, V=fRes0$SigmaXHat, P=fRes0$SigmaX[2:length(fRes0$SigmaX)])

    fRes <- squareRootKF(A=eRes$A, B=B0, C=eRes$C, D=D0, xHat0=eRes$xHat0, SRSigmaX0=chol(x=eRes$V0), SRSigmaW=chol(x=eRes$Gamma), SRSigmaV=chol(x=eRes$Sigma), us=us, zs=zs)
    sRes <- smoothLDS(A=eRes$A, mu=fRes$xHat, V=fRes$SigmaXHat, P=fRes$SigmaX[2:length(fRes$SigmaX)])

    df <- data.frame(x=c(t, t, t, t), y=c(sRes0$muHat[coordinatesToPlot[1],], sRes0$muHat[coordinatesToPlot[2],], sRes$muHat[coordinatesToPlot[1],], sRes$muHat[coordinatesToPlot[2],]), type=factor(c(rep("initial cond. latent1", nTimeSamples), rep("initial cond. latent2", nTimeSamples), rep("estimated latent1", nTimeSamples), rep("estimated latent2", nTimeSamples))))
    p2 <- ggplot(data=df, aes(x=x, y=y, colour=type))
    p2 <- p2 + geom_point() + geom_line()
    p2 <- p2 + scale_colour_manual(values=c("initial cond. latent1"="orchid", "initial cond. latent2"="orchid4", "estimated latent1"="orangered", "estimated latent2"="darkred"))
    p2 <- p2 + xlab("Time")
    p2 <- p2 + ylab("")
    p2 <- p2 + theme(legend.title = element_blank()) 
    #
    grid.arrange(p1, p2, ncol = 1)
    p <- arrangeGrob(p1, p2, ncol=1)
    ggsave(filename=figFilename, plot=p)

    browser()
}

processAll()
