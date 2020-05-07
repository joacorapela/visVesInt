
require(ggplot2)
source("getDataSubset.R")

processAll <- function() {
    xlab <- "Trial Number"
    ylab <- "Avg Firing Rate Per Trial"
    titlePattern <- "Experiment %02d, #Cells %02d"
	dataFilename <- "../../data/011620/Spike_Freq_Dataframe_Hugo.csv"
    figFilenamePattern <- "figures/avgFiringRateVsTrialExp%02d.png"

	data <- read.csv(dataFilename)
    avgData <- aggregate(Firing.Rate ~ Expt.. + Cell.. + Trial.., data, mean)
    exptNs <- unique(avgData$Expt..)
    maxFiringRate <- max(avgData$Firing.Rate)
    ylim <- c(0, maxFiringRate)
    for(exptN in exptNs) {
        avgDataForExp <- subset(x=avgData, Expt..==exptN)
        nCells <- length(unique(avgDataForExp$Cell..))
        p <- ggplot(data=avgDataForExp, mapping=aes(x=Trial.., y=Firing.Rate, color=factor(Cell..)))
        p <- p + geom_line()
        p <- p + ylim(ylim)
        p <- p + xlab(xlab)
        p <- p + ylab(ylab)
        p <- p + ggtitle(sprintf(titlePattern, exptN, nCells))
        p <- p + theme(legend.position="none")
        ggsave(filename=sprintf(figFilenamePattern, exptN))
        print(p)
    }

    browser()
}

processAll()

