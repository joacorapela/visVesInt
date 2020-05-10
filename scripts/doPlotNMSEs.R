
require(MARSS)
require(ggplot2)
require(plotly)
require(reshape2)
source("../src/utils/kalmanFilter/buildVelocityInputs.R")
source("../src/utils/signalProcessing/computeNMSE.R")

processAll <- function() {
    exptN <- 6
    sRate <- 10
    bestModelsInfoFilenamePattern <- "results/expt%02d/bestModelsInfo.csv"
    figFilenamePattern <- "figures//expt%02d//exp%02d_MARSS_%s.%s"

    bestModelsInfoFilename <- sprintf(bestModelsInfoFilenamePattern, exptN)
    bestModelsInfo <- read.csv(file=bestModelsInfoFilename)
    resultsFilenamePattern <- "results/expt%02d/expt%02d_MARSS_condition%s_nLat%02d_stateInputMemory%.02f_obsInputMemory%.02f.RData"

    nModels <- nrow(bestModelsInfo)

    nmse <- c()
    for(i in 1:nModels) {
        condition <- as.character(bestModelsInfo[i,"condition"])
        conditionNoBlanks <- gsub(" ", "_", condition)
        conditionNoBlanks <- gsub("_\\+_", "_", conditionNoBlanks)
        dimLat <- as.numeric(bestModelsInfo[i,"dimLat"])
        obsInputMemorySecs <- as.numeric(bestModelsInfo[i,"obsInputMemorySecs"])
        stateInputMemorySecs <- as.numeric(bestModelsInfo[i,"stateInputMemorySecs"])
        resultsFilename <- sprintf(resultsFilenamePattern, exptN, exptN, conditionNoBlanks, dimLat, stateInputMemorySecs, obsInputMemorySecs)
        res <- get(load(resultsFilename))
        kem <- res$kem
        firingRatesDataBlock <- as.matrix(res$firingRatesDataBlock)
        velocitiesDataBlock <- res$velocitiesDataBlock
        kfRes <- MARSSkf(kem)

        if(!is.nan(obsInputMemorySecs)) {
            obsVelInputs <- buildVelocityInputs(velocities=as.numeric(velocitiesDataBlock[1,]), inputMemorySecs=obsInputMemorySecs, sRate=sRate)
        } else {
            obsVelInputs <- "zero"
        }

        dimObs <- nrow(firingRatesDataBlock)
        Z <- matrix(coef(kem)$Z, nrow=dimObs)
        D <- matrix(coef(kem)$D, nrow=dimObs)
        if(!is.nan(obsInputMemorySecs)) {
            predictions <- Z%*%kfRes$xtt1+as.numeric(coef(kem)$A)+D%*%obsVelInputs
        } else {
            predictions <- Z%*%kfRes$xtt1+as.numeric(coef(kem)$A)
        }
        nmsei <- computeNMSEs(observations=firingRatesDataBlock,
                              predictions=predictions)
        nmsei <- matrix(nmsei, ncol=1)
        colnames(nmsei) <- condition
        rownames(nmsei) <- rownames(firingRatesDataBlock)
        nmse <- cbind(nmse, nmsei)
    }

    nmseM <- melt(nmse)
    p <- ggplot(data=nmseM, aes(x=Var1, y=value, fill=Var2))
    p <- p + geom_bar(stat="identity", position="dodge")
    p <- p + ylab("NMSE")
    p <- p + xlab("")
    p <- p + theme(legend.title = element_blank())
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))

    pngFigFilename <- sprintf(figFilenamePattern, exptN, exptN, "conditionsNMSE", "png")
    ggsave(plot=p, filename=pngFigFilename)
    #
    p <- ggplotly(p)
    htmlFigFilename <- sprintf(figFilenamePattern, exptN, exptN, "conditionsNMSE", "html")
    htmlwidgets::saveWidget(as_widget(p), file.path(normalizePath(dirname(htmlFigFilename)), basename(htmlFigFilename)))

    browser()
}

processAll()
