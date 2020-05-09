
require(MARSS)
require(ggplot2)
require(plotly)
require(reshape2)
source("../src/utils/signalProcessing/computeNMSE.R")
source("../src/utils/kalmanFilter/buildVelocityInputs.R")
source("../src/plot/kalmanFilter/getPlotTrueInitialAndEstimatedMatrices.R")
source("../src/plot/kalmanFilter/getPlotTrueInitialAndEstimatedVectors.R")
source("../src/utils/plot/circleFun.R")

processAll <- function() {
    dimLat <- 8
    sRate <- 10
    stateInputMemorySecs <- NaN
    obsInputMemorySecs <- 0.1
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

    conditionNoBlanks <- gsub(" ", "_", condition)
    conditionNoBlanks <- gsub("_\\+_", "_", conditionNoBlanks)

    resultsFilenamePattern <- "results/expt%02d/expt%02d_MARSS_condition%s_nLat%02d_stateInputMemory%.02f_obsInputMemory%.02f.RData"
    resultsFilename <- sprintf(resultsFilenamePattern, exptN, exptN, conditionNoBlanks, dimLat, stateInputMemorySecs, obsInputMemorySecs)

    results <- get(load(resultsFilename))
    kem <- results$kem
    velocitiesDataBlock <- results$velocitiesDataBlock
    initialConds <- list(A=matrix(kem$start$B, ncol=dimLat), C=matrix(kem$start$Z, ncol=dimLat), sigmaDiag=as.vector(kem$start$R))
    kfRes <- MARSSkf(kem)
    kem0 <- kem
    kem0$par <- kem0$start
    kfRes0 <- MARSSkf(kem0)

    firingRatesDataBlock <- results$firingRatesDataBlock
    trialStartTimes <- results$trialStartTimes

    dimObs <- nrow(firingRatesDataBlock)
    nObs <- ncol(firingRatesDataBlock)

    figFilenamePattern <- "figures//expt%02d//expt%02d_condition%s_nLat%02d_stateInputMemory%.02f_obsInputMemory%.02f_MARSS_%s.%s"

    # begin plot true params

    APNGFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, dimLat, stateInputMemorySecs, obsInputMemorySecs, "A", "png")
    AHTMLFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, dimLat, stateInputMemorySecs, obsInputMemorySecs, "A", "html")
    p <- getPlotTrueInitialAndEstimatedMatrices(initial=initialConds$A, estimated=matrix(coef(kem)$B, nrow=dimLat), title="A")
    ggsave(plot=p, filename=APNGFilename)
    pPlotly <- ggplotly(p)
    htmlwidgets::saveWidget(as_widget(pPlotly), file.path(normalizePath(dirname(AHTMLFilename)), basename(AHTMLFilename)))
    #
    uPNGFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, dimLat, stateInputMemorySecs, obsInputMemorySecs, "u", "png")
    uHTMLFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, dimLat, stateInputMemorySecs, obsInputMemorySecs, "u", "html")
    p <- getPlotTrueInitialAndEstimatedVectors(estimated=matrix(coef(kem)$U, nrow=dimLat), title="u", xlab="Latent Index")
    ggsave(plot=p, filename=uPNGFilename)
    pPlotly <- ggplotly(p)
    htmlwidgets::saveWidget(as_widget(pPlotly), file.path(normalizePath(dirname(uHTMLFilename)), basename(uHTMLFilename)))
    #
    CPNGFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, dimLat, stateInputMemorySecs, obsInputMemorySecs, "C", "png")
    CHTMLFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, dimLat, stateInputMemorySecs, obsInputMemorySecs, "C", "html")
    p <- getPlotTrueInitialAndEstimatedMatrices(initial=t(initialConds$C), estimated=t(matrix(coef(kem)$Z, nrow=dimObs)), title="C", xlab="Latent Index")
    ggsave(plot=p, filename=CPNGFilename)
    pPlotly <- ggplotly(p)
    htmlwidgets::saveWidget(as_widget(pPlotly), file.path(normalizePath(dirname(CHTMLFilename)), basename(CHTMLFilename)))
    #
    aPNGFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, dimLat, stateInputMemorySecs, obsInputMemorySecs, "a", "png")
    aHTMLFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, dimLat, stateInputMemorySecs, obsInputMemorySecs, "a", "html")
    p <- getPlotTrueInitialAndEstimatedVectors(estimated=matrix(coef(kem)$A, nrow=dimObs), title="a", xlab="Neuron Index")
    ggsave(plot=p, filename=aPNGFilename)
    pPlotly <- ggplotly(p)
    htmlwidgets::saveWidget(as_widget(pPlotly), file.path(normalizePath(dirname(aHTMLFilename)), basename(aHTMLFilename)))
    #
    GammaPNGFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, dimLat, stateInputMemorySecs, obsInputMemorySecs, "Gamma", "png")
    GammaHTMLFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, dimLat, stateInputMemorySecs, obsInputMemorySecs, "Gamma", "html")
    p <- getPlotTrueInitialAndEstimatedVectors(estimated=coef(kem)$Q, title="Gamma Variance")
    ggsave(plot=p, filename=GammaPNGFilename)
    pPlotly <- ggplotly(p)
    htmlwidgets::saveWidget(as_widget(pPlotly), file.path(normalizePath(dirname(GammaHTMLFilename)), basename(GammaHTMLFilename)))
    #
    SigmaPNGFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, dimLat, stateInputMemorySecs, obsInputMemorySecs, "Sigma", "png")
    SigmaHTMLFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, dimLat, stateInputMemorySecs, obsInputMemorySecs, "Sigma", "html")
    p <- getPlotTrueInitialAndEstimatedVectors(initial=initialConds$sigmaDiag, estimated=coef(kem)$R, title="Sigma Diagonal")
    ggsave(plot=p, filename=SigmaPNGFilename)
    pPlotly <- ggplotly(p)
    htmlwidgets::saveWidget(as_widget(pPlotly), file.path(normalizePath(dirname(SigmaHTMLFilename)), basename(SigmaHTMLFilename)))
    #
    x0PNGFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, dimLat, stateInputMemorySecs, obsInputMemorySecs, "x0", "png")
    x0HTMLFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, dimLat, stateInputMemorySecs, obsInputMemorySecs, "x0", "html")
    p <- getPlotTrueInitialAndEstimatedVectors(estimated=coef(kem)$x0, title="x0", xlab="Latent Index")
    ggsave(plot=p, filename=x0PNGFilename)
    pPlotly <- ggplotly(p)
    htmlwidgets::saveWidget(as_widget(pPlotly), file.path(normalizePath(dirname(x0HTMLFilename)), basename(x0HTMLFilename)))
    #
    V0PNGFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, dimLat, stateInputMemorySecs, obsInputMemorySecs, "V0", "png")
    V0HTMLFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, dimLat, stateInputMemorySecs, obsInputMemorySecs, "V0", "html")
    p <- getPlotTrueInitialAndEstimatedVectors(estimated=coef(kem)$V0, title="V0")
    ggsave(plot=p, filename=V0PNGFilename)
    pPlotly <- ggplotly(p)
    htmlwidgets::saveWidget(as_widget(pPlotly), file.path(normalizePath(dirname(V0HTMLFilename)), basename(V0HTMLFilename)))
    #
    if(!is.nan(stateInputMemorySecs)) {
        CInputPNGFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, dimLat, stateInputMemorySecs, obsInputMemorySecs, "CInput", "png")
        CInputHTMLFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, dimLat, stateInputMemorySecs, obsInputMemorySecs, "CInput", "html")
        p <- getPlotTrueInitialAndEstimatedMatrices(estimated=t(matrix(coef(kem)$C), nrow=dimObs), title="C Input")
        ggsave(plot=p, filename=CInputPNGFilename)
        pPlotly <- ggplotly(p)
        htmlwidgets::saveWidget(as_widget(pPlotly), file.path(normalizePath(dirname(CInputHTMLFilename)), basename(CInputHTMLFilename)))
    }
    if(!is.nan(obsInputMemorySecs)) {
        DInputPNGFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, dimLat, stateInputMemorySecs, obsInputMemorySecs, "DInput", "png")
        DInputHTMLFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, dimLat, stateInputMemorySecs, obsInputMemorySecs, "DInput", "html")
        p <- getPlotTrueInitialAndEstimatedMatrices(estimated=t(matrix(coef(kem)$D, nrow=dimObs)), title="D Input")
        ggsave(plot=p, filename=DInputPNGFilename)
        pPlotly <- ggplotly(p)
        htmlwidgets::saveWidget(as_widget(pPlotly), file.path(normalizePath(dirname(DInputHTMLFilename)), basename(DInputHTMLFilename)))
    }
    # end plot true params

    # begin plot state dynamics spectrum
    A <- matrix(kem$start$B, ncol=dimLat)
    dfCircle <- circleFun(c(0,0), radius=1, npoints = 200)
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
    AeigenPNGFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, dimLat, stateInputMemorySecs, obsInputMemorySecs, "Aeigen", "png")
    AeigenHTMLFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, dimLat, stateInputMemorySecs, obsInputMemorySecs, "Aeigen", "html")
    ggsave(plot=p, filename=AeigenPNGFilename)
    pPlotly <- ggplotly(p)
    htmlwidgets::saveWidget(as_widget(pPlotly), file.path(normalizePath(dirname(AeigenHTMLFilename)), basename(AeigenHTMLFilename)))
    # end plot state dynamics spectrum

    # begin plot predicted observations

    # begin compute one-lag ahead observation predictions stats
    if(!is.nan(obsInputMemorySecs)) {
        obsInputs <- buildVelocityInputs(velocities=as.numeric(velocitiesDataBlock[1,]), inputMemorySecs=obsInputMemorySecs, sRate=sRate)
    } else {
        obsInputs <- NA
    }
    Z <- matrix(coef(kem)$Z, nrow=dimObs)
    D <- matrix(coef(kem)$D, nrow=dimObs)
    if(!is.nan(obsInputMemorySecs)) {
        ytt1 <- Z%*%kfRes$xtt1+as.numeric(coef(kem)$A)+D%*%obsInputs
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

    for(i in 1:dimObs) {
        df <- data.frame(time=(1:ncol(ytt1))/sRate, predMean=ytt1[i,], predMeanLower=ytt1[i,]-1.96*sqrt(Wtt1[i,i,]), predMeanUpper=ytt1[i,]+1.96*sqrt(Wtt1[i,i,]), observation=as.numeric(firingRatesDataBlock[i,]))
        p <- ggplot(data=df, mapping=aes(x=time))
        p <- p + geom_line(aes(y=predMean, col="prediction"))
        p <- p + geom_point(aes(y=observation, col="observation"))
        p <- p + geom_ribbon(mapping=aes(ymin=predMeanLower, ymax=predMeanUpper), alpha=0.3, fill="red")
        p <- p + scale_colour_manual(name="", values=c("prediction"="red", "observation"="black"))
        p <- p + geom_vline(xintercept=trialStartTimes, linetype="dotted")
        p <- p + xlab("Time (sec)")
        p <- p + ylab("Square-Root Transformed Firing Rate")
        p <- p + ggtitle(rownames(firingRatesDataBlock)[i])
        p <- p + theme(legend.title = element_blank())
        suffix <- sprintf("pred%s", rownames(firingRatesDataBlock)[i])
        predPNGFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, dimLat, stateInputMemorySecs, obsInputMemorySecs, suffix, "png")
        predHTMLFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, dimLat, stateInputMemorySecs, obsInputMemorySecs, suffix, "html")
        ggsave(plot=p, filename=predPNGFilename)
        pPlotly <- ggplotly(p)
        htmlwidgets::saveWidget(as_widget(pPlotly), file.path(normalizePath(dirname(predHTMLFilename)), basename(predHTMLFilename)))
        # print(p)
    }
    # end plot predicted observations

    # begin compute nmses
    nmse <- array(rep(NA, dimObs))
    rownames(nmse) <- rownames(firingRatesDataBlock)
    for(i in 1:dimObs) {
        firingRatesDataBlocki <- as.numeric(firingRatesDataBlock[i,])
        nmse[i] <- mean((ytt1[i,]-firingRatesDataBlocki)^2)/var(firingRatesDataBlocki)
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
    nmsePNGFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, dimLat, stateInputMemorySecs, obsInputMemorySecs, "NMSE", "png")
    nmseHTMLFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, dimLat, stateInputMemorySecs, obsInputMemorySecs, "NMSE", "html")
    ggsave(plot=p, filename=nmsePNGFilename)
    pPlotly <- ggplotly(p)
    htmlwidgets::saveWidget(as_widget(pPlotly), file.path(normalizePath(dirname(nmseHTMLFilename)), basename(nmseHTMLFilename)))
    # end plot nmses

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
    latentsPNGFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, dimLat, stateInputMemorySecs, obsInputMemorySecs, "Latents", "png")
    latentsHTMLFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, dimLat, stateInputMemorySecs, obsInputMemorySecs, "Latents", "html")
    ggsave(plot=p, filename=latentsPNGFilename)
    pPlotly <- ggplotly(p)
    htmlwidgets::saveWidget(as_widget(pPlotly), file.path(normalizePath(dirname(latentsHTMLFilename)), basename(latentsHTMLFilename)))
    # end plot latents

    browser()
}

processAll()
