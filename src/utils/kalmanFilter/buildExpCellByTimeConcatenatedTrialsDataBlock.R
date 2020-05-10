
buildExpCellByTimeConcatenatedTrialsDataBlock <- function(data, exptColLabel="Expt..", cellColLabel="Cell..", trialColLabel="Trial..", binColLabel="Bin..", dataColLabel="Firing.Rate") {
    dataBlock <- c()
    colNames <- c()
    rowNames <- c()
    exptCellAndTrialIDs <- unique(data[,c(exptColLabel, cellColLabel, trialColLabel)])
    exptAndCell <- unique(exptCellAndTrialIDs[, c(exptColLabel, cellColLabel)])
    trials <- unique(exptCellAndTrialIDs[, trialColLabel])
    for(i in 1:nrow(exptAndCell)) {
        exptID <- exptAndCell[i, exptColLabel]
        cellID <- exptAndCell[i, cellColLabel]
        rowName <- sprintf("Expt%02dCell%02d", exptID, cellID)
        trialNSamples <- rep(NA, length(trials))
        exptCellDataBlock <- c()
        for(j in 1:length(trials)) {
            trialID <- trials[j]
            indices <- which(data[,exptColLabel]==exptID &
                             data[,cellColLabel]==cellID &
                             data[,trialColLabel]==trialID)
            exptCellTrialData <- data[indices, dataColLabel]
            trialNSamples[j] <- length(exptCellTrialData)
            # make sure I have selected only one trial (i.e., bin numbers are not repeated)
            stopifnot(trialNSamples[j]==2*length(unique(data[indices, binColLabel])))
            #
            exptCellDataBlock <- c(exptCellDataBlock, exptCellTrialData)
        }
        rowNames <- c(rowNames, rowName)
        dataBlock <- rbind(dataBlock, exptCellDataBlock)
    }
    dataBlock <- data.frame(dataBlock)
    colnames(dataBlock) <- colNames
    rownames(dataBlock) <- rowNames
    answer <- list(dataBlock=dataBlock, trialNSamples=trialNSamples)
    return(answer)
}
