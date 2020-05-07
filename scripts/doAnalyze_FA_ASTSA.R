
require(MASS)
require(astsa)
require(pcaMethods)
require(ggplot2)
require(plotly)
require(reshape2)
require(mvtnorm)
require(ramcmc)
require(gridExtra)
source("~/dev/learning/kalmanFilter/code/src/squareRootKF.R")
source("~/dev/learning/kalmanFilter/code/src/smoothLDS.R")
source("~/dev/learning/kalmanFilter/code/src/plotTrueInitialAndEstimatedMatrices.R")
source("~/dev/learning/kalmanFilter/code/src/plotTrueInitialAndEstimatedVectors.R")
source("~/dev/learning/kalmanFilter/code/src/estimateKFInitialCondFA.R")

source("getDataSubset.R")
source("buildExpTrialBySpeedFiringRateDataBlock.R")

processAll <- function() {
    nFactors <- 2
    maxIter <- 1000
    tol <- 1e-6
    # exptN <- NA
    exptN <- 6
    cellN <- NA
    trialN <- NA
    # condition <- NA
    condition <- "Visual + Vestibular"
    # condition <- "Visual"
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
    coordinatesToPlot <- c(1, 2)

    sRate <- 10
    nIter <- 50
    dataFilename <- "../../data/011620/Spike_Freq_Dataframe_Hugo.csv"
    figFilename <- "figures/kfEstimates.png"

    data <- read.csv(dataFilename)
    dataSubset <- getDataSubset(data=data, exptN=exptN, cellN=cellN, trialN=trialN, condition=condition, region=region, layer=layer, cellType=cellType)
    # dataSubsetAvgTrials <- aggregate(x=data.frame("AvgFiringRate"=dataSubset$Firing.Rate), by=dataSubset[, c("Expt..", "Cell..", "Bin..", "Position..deg.", "Speed..deg.s.", "Trial.Condition", "Region", "Layer", "Cell.Type")], FUN=mean)
    # zs <- buildExpTrialBySpeedFiringRateDataBlock(data=dataSubsetAvgTrials, firingRateColLabel="AvgFiringRate")
    zs <- buildExpTrialBySpeedFiringRateDataBlock(data=dataSubset, firingRateColLabel="Firing.Rate")
    zs <- sqrt(zs)
    zsForFA <- t(as.matrix(zs))
    initialCond <- estimateKFInitialCondFA(z=zsForFA, nFactors=nFactors)
    # initialCond <- estimateKFInitialCondPPCA(z=zsForFA, nFactors=nFactors)
    nObs <- ncol(zs)
    A0 <- initialCond$A
    Gamma0 <- sdGamma*diag(rep(1, ncol(A0)))
    SRSigmaW0 <- chol(Gamma0)
    C0 <- initialCond$C
    Sigma0 <- diag(initialCond$sigmaDiag)
    D0 <- matrix(0, nrow=nrow(C0), ncol=1)
    SRSigmaV0 <- chol(Sigma0)
    B0 <- matrix(0, nrow=nrow(A0), ncol=1)
    xHat00 <- rep(0, ncol(A0))
    V00 <- sdSigmaX00*diag(rep(1, length(xHat00)))
    SRSigmaX00 <- chol(V00)
    nTimeSamples <- ncol(zs)
    us <- matrix(0, nrow=1, ncol=nTimeSamples)

    emRes <- EM0(num=ncol(zs), y=t(zs), A=C0, mu0=xHat00, Sigma0=V00, Phi=A0, cQ=SRSigmaW0, cR=SRSigmaV0, max.iter=maxIter, tol=tol)

    figFilenamePattern <- "figures//exp%02d_ASTSA_%s.html"
    df <- data.frame(x=1:length(emRes$like), 
                     y=emRes$like)
    p <- ggplot(df, aes(x=x, y=y))
    p <- p + geom_line()
    p <- p + geom_point()
    p <- p + xlab("Time")
    p <- p + ylab("Log Likelihood")
    p <- ggplotly(p)
    llFigFilename <- sprintf(figFilenamePattern, exptN, "LogLik")
    htmlwidgets::saveWidget(as_widget(p), file.path(normalizePath(dirname(llFigFilename)), basename(llFigFilename)))
    print(p)

    AFigFilename <- sprintf(figFilenamePattern, exptN, "A")
    plotTrueInitialAndEstimatedMatrices(initial=A0, estimated=emRes$Phi, title="A", figFilename=AFigFilename)
    #
    CFigFilename <- sprintf(figFilenamePattern, exptN, "C")
    plotTrueInitialAndEstimatedMatrices(initial=C0, title="C", figFilename=CFigFilename)
    #
    GammaFigFilename <- sprintf(figFilenamePattern, exptN, "Gamma")
    plotTrueInitialAndEstimatedMatrices(initial=Gamma0, estimated=emRes$Q, title="Gamma", figFilename=GammaFigFilename)
    #
    SigmaFigFilename <- sprintf(figFilenamePattern, exptN, "Sigma")
    plotTrueInitialAndEstimatedMatrices(initial=Sigma0, estimated=emRes$R, title="Sigma", figFilename=SigmaFigFilename)
    #
    V0FigFilename <- sprintf(figFilenamePattern, exptN, "V00")
    plotTrueInitialAndEstimatedMatrices(initial=V00, estimated=emRes$Sigma0, title="V00", figFilename=V0FigFilename)
    #
    xHat0FigFilename <- sprintf(figFilenamePattern, exptN, "xHat0")
    plotTrueInitialAndEstimatedVectors(initial=xHat00, estimated=emRes$mu0, title="xHat0", figFilename=xHat0FigFilename)
    #
    fRes0 <- squareRootKF(A=A0, B=B0, C=C0, D=D0, xHat0=xHat00, SRSigmaX0=SRSigmaX00, SRSigmaW=SRSigmaW0, SRSigmaV=SRSigmaV0, us=us, zs=zs)
    sRes0 <- smoothLDS(A=A0, mu=fRes0$xHat, V=fRes0$SigmaXHat, P=fRes0$SigmaX[2:length(fRes0$SigmaX)])
    #
    fRes <- squareRootKF(A=emRes$Phi, B=B0, C=C0, D=D0, xHat0=emRes$mu0, SRSigmaX0=chol(x=emRes$Sigma0), SRSigmaW=chol(x=emRes$Q), SRSigmaV=chol(x=emRes$R), us=us, zs=zs)
    sRes <- smoothLDS(A=emRes$Phi, mu=fRes$xHat, V=fRes$SigmaXHat, P=fRes$SigmaX[2:length(fRes$SigmaX)])
    #
    t <- (1:ncol(zs))/sRate
    df <- data.frame(x=c(t, t, t, t), y=c(sRes0$muHat[coordinatesToPlot[1],], sRes0$muHat[coordinatesToPlot[2],], sRes$muHat[coordinatesToPlot[1],], sRes$muHat[coordinatesToPlot[2],]), type=factor(c(rep("initial cond. latent1", nTimeSamples), rep("initial cond. latent2", nTimeSamples), rep("estimated latent1", nTimeSamples), rep("estimated latent2", nTimeSamples))))
    p <- ggplot(data=df, aes(x=x, y=y, colour=type))
    p <- p + geom_point() + geom_line()
    p <- p + scale_colour_manual(values=c("initial cond. latent1"="orchid", "initial cond. latent2"="orchid4", "estimated latent1"="orangered", "estimated latent2"="darkred"))
    p <- p + xlab("Time")
    p <- p + ylab("")
    p <- p + theme(legend.title = element_blank())
    p <- ggplotly(p)
    latentsFigFilename <- sprintf(figFilenamePattern, exptN, "Latents")
    htmlwidgets::saveWidget(as_widget(p), file.path(normalizePath(dirname(latentsFigFilename)), basename(latentsFigFilename)))
    print(p)

    browser()
}

processAll()
