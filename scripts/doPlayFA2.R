
require(corrplot)

processAll <- function() {
    bioData <- read.csv("http://msekce.karlin.mff.cuni.cz/~maciak/NMST539/bioData.csv", header = T)
    corrplot(cor(bioData[,2:18]), method="ellipse")

    browser()
}

processAll()
