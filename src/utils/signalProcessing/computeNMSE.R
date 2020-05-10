
computeNMSEs <- function(observations, predictions) {
    dimObs <- nrow(observations)
    nmse <- array(rep(NA, dimObs))
    rownames(nmse) <- rownames(observations)
    for(i in 1:dimObs) {
        nmse[i] <- mean((observations[i,]-predictions[i,])^2)/var(observations[i,])
    }
    return(nmse)
}

