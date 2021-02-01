expdat <- read.csv("exp1_data.csv", header=TRUE)
f <- aggregate(soilrespyrg ~ suillus + fert, data = expdat, FUN = mean)
f
