require(gdata)
require(ProjectTemplate)
require(data.table)
load.project()
sheets<-c("Tcell","Treg","Bcell","DCMonoNK","Thelper")
sheetind<-1:5
manual<-sapply(sheetind,function(sheetindex)xlsx::read.xlsx("data/To FlowCAP/CA c3 v2 ALL SITES.xlsx",sheetIndex=sheetindex,check.names=FALSE))
names(manual)<-sheets


 #Import and munge the manual gates.
 Tcell <- prepareManualTcellGates(manual[["Tcell"]])
 Bcell <- prepareManualBcellGates(manual[["Bcell"]])
 Mono <- prepareManualMonocyteGates(manual[["DCMonoNK"]])
 Treg <- prepareManualTregGates(manual[["Treg"]])
 Thelper <- prepareManualThelperGates(manual[["Thelper"]])


## Import automated gates from OpenCyto
OpenCytoGateFiles <- list.files(path="data/AutomatedGating/",pattern="OpenCyto",full=TRUE)
OCGates <- lapply(OpenCytoGateFiles,function(x)read.csv(x,check.names=FALSE))
names(OCGates) <- do.call(c,lapply(strsplit(basename(OpenCytoGateFiles),"-"),function(x)gsub("\\.csv","",x[[3]])))


## Import automated gates from flowDensity
flowDensityGateFiles <- list.files(path="data/AutomatedGating/",pattern="Brinkman",full=TRUE)
FDGates <- lapply(flowDensityGateFiles[c(1,3,4)],function(x)read.csv(x,check.names=FALSE))
names(FDGates) <- do.call(c,lapply(strsplit(basename(flowDensityGateFiles[c(1,3,4)]),"-"),function(x)x[[2]]))
FDGates[[4]] <- read.table(flowDensityGateFiles[[2]],sep="\t",header=TRUE)
names(FDGates)[4] <- do.call(c,lapply(strsplit(basename(flowDensityGateFiles[c(2)]),"-"),function(x)x[[2]]))



                  ##Rename the levels of a factor using mapvalue or revalue from plyr
require(plyr)
