
require(ggplot2)
source("getDataSubset.R")

processAll <- function() {
    xlab <- "Bin"
    ylab <- "Avg Firing Rate Across Trials"
    vlinestyle <- "dotted"
    errorBarCol <- "darkgray"
    titlePattern <- "Experiment %02d, Cell %02d, Mean Firing Rate %.02f (spikes/sec)"
	dataFilename <- "../../data/011620/Spike_Freq_Dataframe_Hugo.csv"
    dirnamePattern <- "figures/exp%02d"
    figFilenamePattern <- "%s/tunningCurveExp%02dCell%02d.png"

    sem <- function(x) {sd(x)/sqrt(length(x))}
	data <- read.csv(dataFilename)
    meanData <- aggregate(Firing.Rate ~ Expt.. + Cell.. + Bin.., data, mean)
    names(meanData)[names(meanData)=="Firing.Rate"] <- "meanFiringRate"
    semData <- aggregate(Firing.Rate ~ Expt.. + Cell.. + Bin.., data, sem)
    names(semData)[names(semData)=="Firing.Rate"] <- "semFiringRate"
    df <- cbind(meanData, semFiringRate=semData["semFiringRate"])
    exptNs <- unique(df$Expt..)
    for(exptN in exptNs) {
        dirName <- sprintf(dirnamePattern, exptN)
        if(!dir.exists(dirName)) {
            dir.create(dirName)
        }
        dfForExp <- subset(x=df, Expt..==exptN)
        cellNs <- unique(dfForExp$Cell..)
        for(cellN in cellNs) {
            show(sprintf("Processing experiment %02d, cell %02d", exptN, cellN))
            dfForExpAndCell <- subset(x=dfForExp, Cell..==cellN)
            maxMeanFiringRate <- max(dfForExpAndCell$meanFiringRate)
            medianBin <- median(dfForExpAndCell$Bin..)
            ylim <- c(0, maxMeanFiringRate)
            p <- ggplot(data=dfForExpAndCell, mapping=aes(x=Bin.., y=meanFiringRate))
            p <- p + geom_line()
            p <- p + geom_point()
            p <- p + geom_errorbar(mapping=aes(ymin=meanFiringRate-semFiringRate, ymax=meanFiringRate+semFiringRate), color=errorBarCol)
            p <- p + geom_vline(xintercept=medianBin, linetype=vlinestyle)
            # p <- p + ylim(ylim)
            p <- p + xlab(xlab)
            p <- p + ylab(ylab)
            p <- p + ggtitle(sprintf(titlePattern, exptN, cellN, mean(dfForExpAndCell$meanFiringRate)))
            ggsave(filename=sprintf(figFilenamePattern, dirName, exptN, cellN))
            print(p)
        }
    }

    browser()
}

processAll()

