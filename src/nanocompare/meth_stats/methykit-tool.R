#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#
#BiocManager::install("methylKit")

library(methylKit)
file.list=list(system.file("extdata", "test1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "test2.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
                system.file("extdata", "control2.myCpG.txt", package = "methylKit") )


# read the files to a methylRawListDB object: myobjDB
# and save in databases in folder methylDB
myobjDB=methRead(file.list,
                 sample.id=list("test1","test2","ctrl1","ctrl2"),
                 assembly="hg18",
                 treatment=c(1,1,0,0),
                 context="CpG",
                 dbtype = "tabix",
                 dbdir = "methylDB"
)

print(myobjDB[[1]]@dbpath)

getMethylationStats(myobj[[2]],plot=TRUE,both.strands=FALSE)
