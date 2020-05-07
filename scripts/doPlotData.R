
require(ggplot2)
source("getDataSubset.R")

getPlotTrial <- function(firingRates, speeds, sRate=10) {
    N <- length(firingRates)
    times <- (1:N)/sRate
    df <- data.frame(x=c(times, times), y=c(firingRates, speeds), type=factor(c(rep("firing rate", N), rep("speed", N))))
    p <- ggplot(data=df, map=aes(x=x, y=y, col=type))
    p <- p + geom_line()
}

processAll <- function() {
    exptN <- 13
    cellN <- 18
    trialN <- 1:5
    # condition <- "Visual + Vestibular"
    condition <- "Visual"
    region <- "V1"
    layer <- "Layer 6"
    cellType <- "Wide"
	data <- read.csv("../../data/011620/Spike_Freq_Dataframe_Hugo.csv")
    dataSubset <- getDataSubset(data=data, exptN=exptN, cellN=cellN, trialN=trialN, condition=condition, region=region, layer=layer, cellType=cellType)
    #
    p <- ggplot(data=dataSubset, map=aes(x=Speed..deg.s., y=Firing.Rate))
    p <- p + geom_point()
    print(p)
    browser()

    exptN <- 10
    cellN <- NA
    trialN <- NA
    condition <- "Visual + Vestibular"
    region <- "RSPg"
    layer <- "Layer 6"
    cellType <- "Wide"
	data <- read.csv("../../data/011620/Spike_Freq_Dataframe_Hugo.csv")
    dataSubset <- getDataSubset(data=data, exptN=exptN, cellN=cellN, trialN=trialN, condition=condition, region=region, layer=layer, cellType=cellType)
    firingRatesDataBlock <- buildExpTrialBySpeedFiringRateDataBlock(data=dataSubset)
    # i=3; p <- getPlotTrial(firingRates=as.numeric(firingRatesDataBlock[i,]), speeds=as.numeric(colnames(firingRatesDataBlock))); print(p)
}

processAll()

