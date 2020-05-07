
source("getDataSubset.R")

processAll <- function() {
    # [1] "X"               "Expt.."          "Cell.."          "Bin.."          
    # [5] "Trial.."         "Position..deg."  "Speed..deg.s."   "Firing.Rate"    
    # [9] "Trial.Condition" "Region"          "Layer"           "Cell.Type"      

	data <- read.csv("../../data/011620/Spike_Freq_Dataframe_Hugo.csv")
    dataSummaryFilename <- "results/dataSummary.RData"
    dataSummary <- data.frame()
    exptNs <- unique(data$Expt..)
    conditions <- unique(data$Trial.Condition)
    regions <- unique(data$Region)
    layers <- unique(data$Layer)
    cellTypes <- unique(data$Cell.Type)
    #
    for(exptN in exptNs) {
        for(condition in conditions) {
            for(region in regions) {
                for(layer in layers) {
                    for(cellType in cellTypes) {
                        dataSubset <- getDataSubset(data=data, exptN=exptN, trialN=NA, condition=condition, region=region, layer=layer, cellType=cellType)
                        nCells = length(unique(dataSubset$Cell..))
                        meanFiringRate = mean(dataSubset$Firing.Rate)
                        dataSummary <- rbind(dataSummary, data.frame("Expt.."=exptN, 
                                                             "Trial.Condition"=condition,
                                                             "Region"=region,
                                                             "Layer"=layer,
                                                             "Cell.Type"=cellType,
                                                             "nCells"=nCells,
                                                             "meanFiringRate"=meanFiringRate))
                    }
                }
            }
        }
    }
    save(dataSummary, file=dataSummaryFilename)
    browser()
}

processAll()

