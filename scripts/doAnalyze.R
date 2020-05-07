
require(ggplot2)
require(mvtnorm)
require(ramcmc)                                                       
require(gridExtra)

source("~/dev/learning/kalmanFilter/code/src/lsolve.R")
source("~/dev/learning/kalmanFilter/code/src/chol_downdate_higherOrder.R")
source("~/dev/learning/kalmanFilter/code/src/squareRootKF.R")
source("~/dev/learning/kalmanFilter/code/src/smoothLDS.R")
source("~/dev/learning/kalmanFilter/code/src/emEstimationSquareRootKF.R")
source("getDataSubset.R")
source("buildExpTrialBySpeedFiringRateDataBlock.R")

processAll <- function() {
    exptN <- 11
    cellN <- NA
    trialN <- NA
    condition <- "Visual + Vestibular"
    region <- "RSPg"
    layer <- "Layer 2/3"
    cellType <- "Wide"
    sdW0 <- .01
    sdV0 <- .01
    sdXHat00 <- .01
    sdSigmaX00 <- .01
    # A0 <- rbind(c(1.0, -5e-4), c(.05, 1.0))
    A0 <- rbind(c(1.0, 0), c(0, 1.0))
    B0 <- matrix(0, nrow=nrow(A0), ncol=1)
    SRSigmaW0 <- sdW0*diag(rep(1, nrow(A0)))
    xHat00 <- rep(0, nrow(A0))
    SRSigmaX00 <- sdSigmaX00*diag(rep(1, length(xHat00)))
    sRate <- 10
    nIter <- 50
    dataFilename <- "../../data/011620/Spike_Freq_Dataframe_Hugo.csv"
    figFilename <- "figures/kfEstimates.png"

	data <- read.csv(dataFilename)
    dataSubset <- getDataSubset(data=data, exptN=exptN, cellN=cellN, trialN=trialN, condition=condition, region=region, layer=layer, cellType=cellType)
    dataSubsetAvgTrials <- aggregate(x=data.frame("AvgFiringRate"=dataSubset$Firing.Rate), by=dataSubset[, c("Expt..", "Cell..", "Bin..", "Position..deg.", "Speed..deg.s.", "Trial.Condition", "Region", "Layer", "Cell.Type")], FUN=mean)
    zs <- buildExpTrialBySpeedFiringRateDataBlock(data=dataSubsetAvgTrials, firingRateColLabel="AvgFiringRate")
    nObs <- nrow(zs)
    # C0 <- cbind(rep(0, nObs), rep(1, nObs))
    C0 <- matrix(rnorm(2*nObs), nrow=nObs)
    D0 <- matrix(0, nrow=nObs, ncol=1)
    SRSigmaV0 <- sdV0*diag(rep(1, nObs))
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

    df <- data.frame(x=c(t, t, t, t), y=c(sRes0$muHat[1,], sRes0$muHat[2,], sRes$muHat[1,], sRes$muHat[2,]), type=factor(c(rep("initial cond. latent1", nTimeSamples), rep("initial cond. latent2", nTimeSamples), rep("estimated latent1", nTimeSamples), rep("estimated latent2", nTimeSamples))))
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
