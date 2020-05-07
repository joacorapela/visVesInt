
require(ggplot2)
require(plotly)
require(reshape2)
source("~/dev/learning/kalmanFilter/code/src/plotTrueInitialAndEstimatedMatrices.R")
source("~/dev/learning/kalmanFilter/code/src/plotTrueInitialAndEstimatedVectors.R")
source("circleFun.R")

processAll <- function() {
    nFactors <- 6
    sRate <- 10
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

    resultsFilenamePattern <- "results/expt%02d_MARSS_condition%s_nLat%02d.RData"
    resultsFilename <- sprintf(resultsFilenamePattern, exptN, conditionNoBlanks, nFactors)

    results <- get(load(resultsFilename))
    kem <- results$kem
    kfRes <- results$kfRes
    kfRes0 <- results$kfRes0
    zs <- results$zs
    trialStartTimes <- results$trialStartTimes

    conditionNoBlanks <- gsub(" ", "_", condition)

    dimLat <- nFactors
    dimObs <- nrow(zs)
    nObs <- ncol(zs)

    figFilenamePattern <- "figures//expt%02d//exp%02d_condition%s_nLat%02d_MARSS_%s.html"

    # begin plot true params
    AFigFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, nFactors, "A")
    plotTrueInitialAndEstimatedMatrices(initial=A, estimated=matrix(coef(kem)$B, nrow=dimLat), title="A", figFilename=AFigFilename)
    #
    uFigFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, nFactors, "u")
    plotTrueInitialAndEstimatedVectors(estimated=matrix(coef(kem)$U, nrow=dimLat), title="u", figFilename=uFigFilename)
    #
    CFigFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, nFactors, "C")
    plotTrueInitialAndEstimatedMatrices(initial=matrix(kem$start$Z, ncol=nFactors), estimated=matrix(coef(kem)$Z, nrow=dimObs), title="C", figFilename=CFigFilename)
    #
    aFigFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, nFactors, "a")
    plotTrueInitialAndEstimatedVectors(estimated=matrix(coef(kem)$A, nrow=dimObs), title="a", figFilename=aFigFilename)
    #
    GammaFigFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, nFactors, "Gamma")
    plotTrueInitialAndEstimatedVectors(estimated=coef(kem)$Q, title="Gamma Variance", figFilename=GammaFigFilename)
    #
    SigmaFigFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, nFactors, "Sigma")
    plotTrueInitialAndEstimatedVectors(initial=kem$start$R, estimated=coef(kem)$R, title="Sigma Diagonal", figFilename=SigmaFigFilename)
    #
    x0FigFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, nFactors, "x0")
    plotTrueInitialAndEstimatedVectors(estimated=coef(kem)$x0, title="x0", figFilename=x0FigFilename)
    #
    V0FigFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, nFactors, "V00")
    plotTrueInitialAndEstimatedVectors(estimated=coef(kem)$x0, title="V00", figFilename=x0FigFilename)
    # end plot true params

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
    AeigenFigFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, nFactors, "Aeigen")
    htmlwidgets::saveWidget(as_widget(p), file.path(normalizePath(dirname(AeigenFigFilename)), basename(AeigenFigFilename)))
    # end plot state dynamics spectrum

    # begin plot predicted observations

    # begin compute one-lag ahead observation predictions stats
    Z <- matrix(coef(kem)$Z, ncol=dimLat)
    ytt1 <- Z%*%kfRes$xtt1+as.numeric(coef(kem)$A)
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
        predFigFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, nFactors, suffix)
        htmlwidgets::saveWidget(as_widget(p), file.path(normalizePath(dirname(predFigFilename)), basename(predFigFilename)))
        # print(p)
    }
    # end plot predicted observations

    # begin compute nmses
    nmse <- array(rep(NA, dimObs))
    rownames(nmse) <- rownames(zs)
    for(i in 1:dimObs) {
        zsi <- as.numeric(zs[i,])
        nmse[i] <- mean((ytt1[i,]-zsi)^2)/var(zsi)
    }
    # end compute nmses

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
    predFigFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, nFactors, "NMSE")
    htmlwidgets::saveWidget(as_widget(p), file.path(normalizePath(dirname(predFigFilename)), basename(predFigFilename)))
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
    p <- ggplotly(p)
    latentsFigFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, nFactors, "Latents")
    htmlwidgets::saveWidget(as_widget(p), file.path(normalizePath(dirname(latentsFigFilename)), basename(latentsFigFilename)))
    # print(p)
    # end plot latents

    browser()
}

processAll()
