library(ggplot2)

wdir="/Users/liuya/PycharmProjects/nano-compare/src/rplot_func"
infn='performance-results.csv'

setwd(wdir)

perf_data <- read.csv(file = infn)

ggplot(perf_data, aes(x=Location, y=Accuracy))+
    geom_bar(stat='identity', fill="forest green")+
    ylab("Accuracy")


