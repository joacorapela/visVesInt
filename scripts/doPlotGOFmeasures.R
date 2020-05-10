
require(ggplot2)
require(plotly)
require(reshape2)

processAll <- function() {
    dimLat <- 2:10
    obsInputMemorySecs <- rep(NaN, times=length(dimLat))
    xVal <- dimLat
    xlab <- "Number of Latents"
    # obsInputMemorySecs <- (0:5)*0.1
    # dimLat <- rep(8, times=length(obsInputMemorySecs))
    # xVal <- obsInputMemorySecs
    # xlab <- "Observation Speed Memory (sec)"
    ylab <- "Goodness of Fit Measure"
    stateInputMemorySecs <- rep(NaN, times=length(dimLat))
    exptN <- 6
    condition <- "Visual + Vestibular"
    # condition <- "Visual"
    # condition <- "Vestibular"
    resultsFilenamePattern <- "results/expt%02d/expt%02d_MARSS_condition%s_nLat%02d_stateInputMemory%.02f_obsInputMemory%.02f.RData"
    figFilenamePattern <- "figures//expt%02d//expt%02d_condition%s_nLat%s_stateInputMemory%s_obsInputMemory%s_MARSS_%s.%s"

    stopifnot(length(dimLat)==length(stateInputMemorySecs) & length(stateInputMemorySecs)==length(obsInputMemorySecs))

    conditionNoBlanks <- gsub(" ", "_", condition)
    conditionNoBlanks <- gsub("_\\+_", "_", conditionNoBlanks)
    # logLik <- rep(NA, length(dimLat))
    AIC <- rep(NA, length(dimLat))
    AICc <- rep(NA, length(dimLat))
    for(i in 1:length(dimLat)) {
        resultsFilename <- sprintf(resultsFilenamePattern, exptN, exptN, conditionNoBlanks, dimLat[i], stateInputMemorySecs[i], obsInputMemorySecs[i])
        res <- get(load(resultsFilename))
        # logLik[i] <- res$kem$logLik
        AIC[i] <- res$kem$AIC
        AICc[i] <- res$kem$AICc
    }

    # df <- data.frame(x=xVal, logLik=logLik, AIC=AIC, AICc=AICc)
    df <- data.frame(x=xVal, AIC=AIC, AICc=AICc)
    mdf <- melt(df, id="x")
    p <- ggplot(data=mdf, mapping=aes(x=x, y=value, color=variable))
    p <- p + geom_line()
    p <- p + geom_point()
    p <- p + xlab(xlab)
    p <- p + ylab(ylab)
    p <- p + theme(legend.title = element_blank())
    pngFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, paste(dimLat, collapse="_"), paste(stateInputMemorySecs, collapse="_"), paste(obsInputMemorySecs, collapse="_"), "GOF", "png")
    htmlFilename <- sprintf(figFilenamePattern, exptN, exptN, conditionNoBlanks, paste(dimLat, collapse="_"), paste(stateInputMemorySecs, collapse="_"), paste(obsInputMemorySecs, collapse="_"), "GOF", "html")
    ggsave(plot=p, filename=pngFilename)
    p <- ggplotly(p)
    htmlwidgets::saveWidget(as_widget(p), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))

    browser()
}

processAll()
