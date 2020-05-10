
require(MASS)
require(MARSS)
require(reshape2)
source("../src/stats/kalmanFilter/estimateKFInitialCondFA.R")
source("../src/utils/kalmanFilter/getDataSubset.R")
source("../src/utils/kalmanFilter/buildVelocityInputs.R")
source("../src/utils/kalmanFilter/buildExpCellByTimeConcatenatedTrialsDataBlock.R")
source("../src/stats/kalmanFilter/fit_FA_MARSS_model.R")

processAll <- function() {
    obsInputMemorySecs <- 0.1*(0:5)
    dimLat <- rep(5, times=length(obsInputMemorySecs))
    stateInputMemorySecs <- rep(NaN, times=length(dimLat))
    sRate <- 10
    maxIter <- 100
    # exptN <- NA
    exptN <- 6
    cellN <- NA
    trialN <- NA
    # condition <- NA
    # condition <- "Visual + Vestibular"
    # condition <- "Visual"
    condition <- "Vestibular"
    region <- NA
    # region <- c("RSPg")
    layer <- NA
    # layer <- "Layer 2/3"
    cellType <- NA
    # cellType <- "Wide"
    sdGamma <- 1000.0
    sdXHat00 <- 1000.0
    sdSigmaX00 <- 1000.0
    dataFilename <- "../../../data/011620/Spike_Freq_Dataframe_Hugo.csv"
    #
    stopifnot(length(dimLat)==length(stateInputMemorySecs) & length(stateInputMemorySecs)==length(obsInputMemorySecs))
    #
    conditionNoBlanks <- gsub(" ", "_", condition)
    conditionNoBlanks <- gsub("_+_", "_", conditionNoBlanks)
    #
    data <- read.csv(dataFilename)
    dataSubset <- getDataSubset(data=data, exptN=exptN, cellN=cellN, trialN=trialN, condition=condition, region=region, layer=layer, cellType=cellType)
    res <- buildExpCellByTimeConcatenatedTrialsDataBlock(data=dataSubset, dataColLabel="Firing.Rate")
    trialStartTimes <- cumsum(res$trialNSamples)/sRate
    firingRatesDataBlock <- res$dataBlock
    firingRatesDataBlock <- sqrt(firingRatesDataBlock)
    velocitiesDataBlock <- buildExpCellByTimeConcatenatedTrialsDataBlock(data=dataSubset, dataColLabel="Speed..deg.s.")$dataBlock
    #
    dimObs <- nrow(firingRatesDataBlock)
    nObs <- ncol(firingRatesDataBlock)
    #
    for(i in 1:length(dimLat)) {
        show(sprintf("Processing number of latents=%d, state input memory=%f sec, observation input memory=%f sec", dimLat[i], stateInputMemorySecs[i], obsInputMemorySecs[i]))
        if(!is.nan(stateInputMemorySecs[i])) {
            stateInputs <- buildVelocityInputs(velocities=as.numeric(velocitiesDataBlock[1,]), inputMemorySecs=stateInputMemorySecs[i], sRate=sRate)
        } else {
            stateInputs <- NA
        }
        if(!is.nan(obsInputMemorySecs[i])) {
            obsInputs <- buildVelocityInputs(velocities=as.numeric(velocitiesDataBlock[1,]), inputMemorySecs=obsInputMemorySecs[i], sRate=sRate)
        } else {
            obsInputs <- NA
        }
        kem <- fit_FA_MARSS_model(firingRatesDataBlock=firingRatesDataBlock, stateInputs=stateInputs, obsInputs=obsInputs, dimLat=dimLat[i], stateInputMemorySecs=stateInputMemorySecs[i], obsInputMemorySecs=obsInputMemorySecs[i], sRate=sRate, maxIter=maxIter)
        show(sprintf("number of latents=%d, state input memory=%f sec, observation input memory=%f sec: logLik=%f, AIC=%f, AICc=%f", dimLat[i], stateInputMemorySecs[i], obsInputMemorySecs[i], kem$logLik, kem$AIC, kem$AICc))
        #
        resultsFilenamePattern <- "results/expt%02d/expt%02d_MARSS_condition%s_nLat%02d_stateInputMemory%.02f_obsInputMemory%.02f.RData"
        resultsFilename <- sprintf(resultsFilenamePattern, exptN, conditionNoBlanks, dimLat[i], stateInputMemorySecs[i], obsInputMemorySecs[i])
        results <- list(kem=kem, firingRatesDataBlock=firingRatesDataBlock, velocitiesDataBlock=velocitiesDataBlock, trialStartTimes=trialStartTimes)
        save(results, file=resultsFilename)
    }
}

processAll()
