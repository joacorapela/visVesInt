
getDataSubset <- function(data, exptN=7, cellN=35, trialN=NA, condition="Visual + Vestibular", region="RSPg", layer="Layer 5", cellType="Wide") {
    if(!is.na(exptN)) {
        indices <- which(data$Expt..%in%exptN)
        data <- data[indices,]
    }
    if(!is.na(cellN[1])) {
        indices <- which(data$Cell..%in%cellN)
        data <- data[indices,]
    }
    if(!is.na(trialN[1])) {
        indices <- which(data$Trial..%in%trialN)
        data <- data[indices,]
    }
    if(!is.na(condition[1])) {
        indices <- which(data$Trial.Condition%in%condition)
        data <- data[indices,]
    }
    if(!is.na(region[1])) {
        indices <- which(data$Region%in%region)
        data <- data[indices,]
    }
    if(!is.na(layer[1])) {
        indices <- which(data$Layer%in%layer)
        data <- data[indices,]
    }
    if(!is.na(cellType)) {
        indices <- which(data$Cell.Type%in%cellType)
        data <- data[indices,]
    }
    show(sprintf("Number of cells %d", length(unique(data$Cell..))))
    return(data)
}

