
require(MASS)
require(ggplot2)
require(mvtnorm)
require(ramcmc)
require(gridExtra)
require(corrplot)

source("~/dev/learning/kalmanFilter/code/src/lsolve.R")
source("~/dev/learning/kalmanFilter/code/src/chol_downdate_higherOrder.R")
source("~/dev/learning/kalmanFilter/code/src/squareRootKF.R")
source("~/dev/learning/kalmanFilter/code/src/smoothLDS.R")
source("~/dev/learning/kalmanFilter/code/src/emEstimationSquareRootKF.R")
# source("~/dev/learning/kalmanFilter/code/src/estimateKFInitialCond.R")
source("getDataSubset.R")
source("buildExpTrialBySpeedFiringRateDataBlock.R")

processAll <- function() {
    exptN <- 6
    cellN <- NA
    trialN <- NA
    condition <- "Visual + Vestibular"
    # condition <- "Visual"
    # condition <- "Vestibular"
    # region <- NA
    region <- "RSPg"
    # region <- c("RSPg")
    layer <- NA
    # layer <- "Layer 2/3"
    cellType <- NA
    # cellType <- "Wide"
    nFactors <- 4
    scores <- "regression"
    rotation <- "varimax"

    dataFilename <- "../../data/011620/Spike_Freq_Dataframe_Hugo.csv"

    data <- read.csv(dataFilename)
    dataSubset <- getDataSubset(data=data, exptN=exptN, cellN=cellN, trialN=trialN, condition=condition, region=region, layer=layer, cellType=cellType)
    zs <- buildExpTrialBySpeedFiringRateDataBlock(data=dataSubset, firingRateColLabel="Firing.Rate")
    zs <- t(as.matrix(zs))
    corrplot(cor(zs), method="ellipse")
    fa <- factanal(~zs, factors=nFactors, rotation=rotation, scores=scores)
    theLoadings <- fa$loadings
    theUniqueness <- fa$uniqueness
    print(fa, digits=2, cutoff=.6, sort=TRUE)
    print("Loadings")
    print(theLoadings)
    print("Variances")
    print(theUniqueness)

    browser()
}

processAll()
