
require(ggplot2)
source("../src/utils/kalmanFilter/getDataSubset.R")

processAll <- function() {
    nSamplesPerTrial <- 70
    xlab <- "Trial Sample"
    ylab <- "Average across Trials of Square-Root Transformed Firing Rates"
    vlinestyle <- "dotted"
    errorBarCol <- "darkgray"
    titlePattern <- "Expt %02d, %s, Cell %02d, Mean Firing Rate %.02f (spikes/sec)"
	dataFilename <- "../../../data/011620/Spike_Freq_Dataframe_Hugo.csv"
    dirnamePattern <- "figures/expt%02d"
    figFilenamePattern <- "%s/tunningCurveExpt%02dCondition%sCell%02d.%s"

    sem <- function(x) {sd(x)/sqrt(length(x))}
	data <- read.csv(dataFilename)
    data[,"Firing.Rate"] <- sqrt(data[,"Firing.Rate"])
    sData <- data[order(data$Expt.., data$Trial.Condition, data$Cell.., data$Trial..),]
    trialSamples <- 1:nSamplesPerTrial
    sData <- cbind(trialSample=trialSamples, sData)
    meanData <- aggregate(Firing.Rate ~ Expt.. + Trial.Condition + Cell.. + trialSample + Position..deg.  + Speed..deg.s., sData, mean)
    names(meanData)[names(meanData)=="Firing.Rate"] <- "meanFiringRate"
    semData <- aggregate(Firing.Rate ~ Expt.. + Trial.Condition + Cell.. + trialSample + Position..deg.  + Speed..deg.s., sData, sem)
    names(semData)[names(semData)=="Firing.Rate"] <- "semFiringRate"
    df <- cbind(meanData, semFiringRate=semData["semFiringRate"])
    exptNs <- unique(df$Expt..)
    conditions <- unique(df$Trial.Condition)
    for(exptN in exptNs) {
        dirName <- sprintf(dirnamePattern, exptN)
        if(!dir.exists(dirName)) {
            dir.create(dirName)
        }
        for(condition in conditions) {
            conditionNoBlanks <- gsub(" ", "_", condition)
            conditionNoBlanks <- gsub("_\\+_", "_", conditionNoBlanks)
            dfForExpAndCond <- subset(x=df, Expt..==exptN & Trial.Condition==condition)
            maxMeanFiringRate <- max(dfForExpAndCond$meanFiringRate+1.96*dfForExpAndCond$semFiringRate)
            minMeanFiringRate <- min(dfForExpAndCond$meanFiringRate-1.96*dfForExpAndCond$semFiringRate)
            ylim <- c(minMeanFiringRate, maxMeanFiringRate)
            cellNs <- unique(dfForExpAndCond$Cell..)
            for(cellN in cellNs) {
                show(sprintf("Processing experiment %02d, condition %s, cell %02d", exptN, condition, cellN))
                dfForExpCondAndCell <- subset(x=dfForExpAndCond, Cell..==cellN)
                # medianBin <- median(dfForExpAndCell$Bin..)
                p <- ggplot(data=dfForExpCondAndCell, mapping=aes(x=trialSample, y=meanFiringRate))
                p <- p + geom_line()
                p <- p + geom_point()
                p <- p + geom_errorbar(mapping=aes(ymin=meanFiringRate-1.96*semFiringRate, ymax=meanFiringRate+1.96*semFiringRate), color=errorBarCol)
                p <- p + geom_vline(xintercept=round(nSamplesPerTrial/2), linetype=vlinestyle)
                p <- p + ylim(ylim)
                p <- p + xlab(xlab)
                p <- p + ylab(ylab)
                p <- p + ggtitle(sprintf(titlePattern, exptN, condition, cellN, mean(dfForExpCondAndCell$meanFiringRate)))
                pngFilename <- sprintf(figFilenamePattern, dirName, exptN, conditionNoBlanks, cellN, "png")
                htmlFilename <- sprintf(figFilenamePattern, dirName, exptN, conditionNoBlanks, cellN, "html")
                ggsave(plot=p, filename=pngFilename)
                pPlotly <- ggplotly(p)
                htmlwidgets::saveWidget(as_widget(pPlotly), file.path(normalizePath(dirname(htmlFilename)), basename(htmlFilename)))
                print(p)
                # browser()
            }
        }
    }

    browser()
}

processAll()

