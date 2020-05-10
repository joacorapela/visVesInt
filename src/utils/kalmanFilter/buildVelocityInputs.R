
buildVelocityInputs <- function(velocities, inputMemorySecs, sRate) {
    inputMemorySamples <- round(inputMemorySecs*sRate)
    velInputs <- c()
    N <- length(velocities)
    for(i in 0:inputMemorySamples) {
        velInputs <- rbind(velInputs, c(rep(0, time=i), velocities[1:(N-i)]))
    }
    return(velInputs)
}

