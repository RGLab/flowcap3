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
Tcell <- rename(melt(Tcell,id=c("Sample","Center","File")),c(variable="population",value="proportion"))
Bcell <- prepareManualBcellGates(manual[["Bcell"]])
Bcell <- rename(melt(Bcell,id=c("Sample","Center","File")),c(variable="population",value="proportion"))


