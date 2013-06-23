require(gdata)
require(ProjectTemplate)
require(data.table)
load.project()
sheets<-c("Tcell","Treg","Bcell","DCMonoNK","Thelper")
sheetind<-1:5
manual<-sapply(sheetind,function(sheetindex)xlsx::read.xlsx("data/To FlowCAP/CA c3 v2 ALL SITES.xlsx",sheetIndex=sheetindex,check.names=FALSE))
names(manual)<-sheets

##' Clean up the T-cell centralized manual gating table from Marc
##'
##' This function will read and clean up the manual gating table from Marc for the T-cells.
##' Unfortunately all the panels are formatted differently and need to be handled separately.
##' 
##' @title prepareManualTcellGates
##' @param manual A data.table of the manual gates for T-cells are passed in, as returned by read.xlsx with check.names=FALSE
##' @return A data frame of just the T-cell statistics
##' @author Greg Finak
prepareManualTcellGates <- function(manual){
                                        #Clean up column names for the T-cells
  Tcell <- manual
  colnames(Tcell)<-(gsub("\\n"," ",colnames(Tcell))) #remove newlines from column headers
  Tcell<-subset(Tcell,!File%in%c("mean","SD","CV")) # remove the estimates of the mean, CV and SD.

                                        #What's in the remaining NA columns? The first column is "Center", the rest we can toss out.
  Tcell<-(Tcell[,-grep("NA",colnames(Tcell))[-1L]]) #drop irrelevant columns
  colnames(Tcell)[grep("NA",colnames(Tcell))] <- "Center" #Second column is center
  colnames(Tcell) <- gsub("\\.([0-9])$","\\.Alternate",colnames(Tcell)) # rename the populations that seem to have been calculated an alternate way (some other parent)
  Tcell <- Tcell[-1L,] # drop the first row, it is empty

                                        #Fill in blank rows for sample and center
  Tcell[is.na(Tcell[,1]),1]<-""
  Tcell[,1]<-factor(fill(Tcell[,1]))
  Tcell[,2] <- as.character(Tcell[,2])
  Tcell[is.na(Tcell[,2]),2] <- ""
  Tcell[,2] <- factor(fill(Tcell[,2]))

                                        #Remove missing data from centers that were excluded
  Tcell<-na.omit(Tcell)
  invisible(Tcell)
}

##' Clean up the B-cell gates
##'
##' Reads and cleans up the manual gates from Marc.
##' @title prepareManualBcellGates
##' @param manual a data.frame read in by read.xlsx with check.names=FALSE
##' @return data.frame of the B-cell statistics.
##' @author Greg Finak
prepareManualBcellGates <- function(manual){
  Bcell <- manual
                                        #Parse rows one and two as a new header
  h1 <- colnames(Bcell)
  h2 <- Bcell[1,]

  Bcell <- Bcell[-1L,]#Drop the first row, which  we have stored.
  colnames(Bcell) <- gsub("\\n"," ",as.matrix(h2)) #Remove newlines from column names
  colnames(Bcell)[2] <- c("Center") #name the column which contains the center
  Bcell <- subset(Bcell,!File%in%c("mean","CV","SD")) #Drop the mean, SD and CV
  tmp <- as.character(Bcell[,1]) #Fill in missing entries
  tmp[is.na(tmp)] <- ""
  Bcell[,1] <- factor(fill(tmp))
  tmp <- as.character(Bcell[,2])
  tmp[is.na(tmp)] <- ""
  Bcell[,2] <- factor(fill(tmp))

                                        #drop NA.3
  Bcell <- Bcell[,-grep("NA.3",colnames(Bcell))]
                                        #NA.1 and NA.2 separate the populatoins that are parented by h1
                                        #prefix the populations with the appropriate parent name
  cur<-eval(as.call(c(`seq`,as.list(grep("NA",colnames(Bcell))+c(1,-1)))))
  colnames(Bcell)[cur] <- paste(h1[11],colnames(Bcell)[cur],sep="/")
  colnames(Bcell)[16:18] <- paste(h1[16],colnames(Bcell)[16:18],sep="/")
                                        #remove NA columns we don't need
  Bcell <- Bcell[,-grep("NA",colnames(Bcell))];
  Bcell <- Bcell[,-4]
  Bcell <- na.omit(Bcell)
  Bcell
}

                                        #Import and munge the manual gates.
Tcell <- prepareManualTcellGates(manual[["Tcell"]])
Tcell <- rename(melt(Tcell,id=c("Sample","Center","File")),c(variable="population",value="proportion"))
Bcell <- prepareManualBcellGates(manual[["Bcell"]])
Bcell <- rename(melt(Bcell,id=c("Sample","Center","File")),c(variable="population",value="proportion"))
