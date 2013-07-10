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

CGates <- list(Tcell=Tcell,Treg=Treg,Bcell=Bcell,DCMonoNK=Mono,Thelper=Thelper)

## Import automated gates from OpenCyto
OpenCytoGateFiles <- list.files(path="data/AutomatedGating/",pattern="OpenCyto",full=TRUE)
OCGates <- lapply(OpenCytoGateFiles,function(x)read.csv(x,check.names=FALSE))
names(OCGates) <- do.call(c,lapply(strsplit(basename(OpenCytoGateFiles),"-"),function(x)gsub("\\.csv","",x[[3]])))


## Import automated gates from flowDensity
flowDensityGateFiles <- list.files(path="data/AutomatedGating/",pattern="Brinkman",full=TRUE)
FDGates <- lapply(flowDensityGateFiles[c(1,3,4,5)],function(x)read.csv(x,check.names=FALSE))
names(FDGates) <- do.call(c,lapply(strsplit(basename(flowDensityGateFiles[c(1,3,4,5)]),"-"),function(x)x[[2]]))
FDGates[[5]] <- read.table(flowDensityGateFiles[[2]],sep="\t",header=TRUE)
names(FDGates)[5] <- do.call(c,lapply(strsplit(basename(flowDensityGateFiles[c(2)]),"-"),function(x)x[[2]]))



## First, order the Panels for each Method
FDGates <- FDGates[names(FDGates)[c(2,4,1,5,3)]]
OCGates <- OCGates[names(OCGates)[c(3,5,1,2,4)]]
CGates <- CGates[sheets]
## Standardize the panel names
names(FDGates) <- sheets
names(OCGates) <- sheets
names(CGates) <- sheets


## What are the column names in each table?
## Drop the first column from all FDGates tables. It is a row number, we don't need it.
for(i in 1:length(FDGates)){
  if(!is.null(FDGates[[i]]))
    FDGates[[i]] <- FDGates[[i]][,-1L]
}

oldnames.FD <- lapply(FDGates,colnames)
oldnames.OC <- lapply(OCGates,colnames)
oldnames.C <- lapply(CGates,colnames)

## For flowDensity gates, compute the proportions where they'remissing
for(i in 1:length(FDGates)){
  if(!any(grepl("Proportion",colnames(FDGates[[i]])))&!is.null(colnames(FDGates[[i]]))){
    FDGates[[i]]$Proportion <- as.numeric(as.character(FDGates[[i]]$`Cell count`))/as.numeric(as.character(FDGates[[i]]$`Parent Cell count`))

  }
}
## For OCgates, add an empty Phenotyp column
for(i in 1:length(OCGates)){
  if(!is.null(OCGates[[i]])){
    OCGates[[i]]$Phenotype <- NA
  }
}
standard.names <- c("Sample","Center","File","Population","Proportion","Phenotype","Count","Parent")

FDmap.1 <- list(Sample="Sample",Center="Center",File="Filename",Population="Population",Proportion="Proportion",Phenotype="Phenotype",Count="Cell count",Parent="Parent population")
FDmap.2 <- list(Sample="Sample",Center="Center",File="Filename",Population="Population",Proportion="Proportion.of.Parent",Phenotype="Phenotype",Count="Cell.Count",Parent="Parent.Population")
OCmap <- list(Sample="Sample",Center="Center",File="Filename",Population="Population",Proportion="Proportion",Phenotype="Phenotype",Count="Count",Parent="Parent")
Cmap <- list(Sample="Sample",Center="Center",File="File",Population="population",Proportion="proportion")


## Organize the central gating data frame
for(i in 1:length(CGates)){
  CGates[[i]] <- CGates[[i]][,do.call(c,Cmap)]
  colnames(CGates[[i]]) <- names(Cmap)
}

## Organize the OpenCyto gating data frame
for(i in seq_along(OCGates)){
  if(!is.null(OCGates[[i]])){
    OCGates[[i]] <- OCGates[[i]][,do.call(c,OCmap)]
    colnames(OCGates[[i]]) <- names(OCmap)
  }
}
## Organize the flowDensity gating data frame
for(i in seq_along(FDGates)){
  if(all(do.call(c,FDmap.1)%in%colnames(FDGates[[i]]))){
    FDmap <- FDmap.1
  }else{
    FDmap <- FDmap.2
  }
  if(!is.null(colnames(FDGates[[i]]))){
    FDGates[[i]] <- FDGates[[i]][,do.call(c,FDmap)]
    colnames(FDGates[[i]]) <- names(FDmap)
  }
}


## Sample Columns should be factors
for(i in 1:length(OCGates)){
  if(!is.null(colnames(OCGates[[i]]))){
    ## Standardize the Sample column
    if(class(OCGates[[i]]$Sample)!="factor"){
      OCGates[[i]]$Sample <- factor(OCGates[[i]]$Sample)
      cat(levels(OCGates[[i]]$Sample),"\n")
    }
    ## Standardize the File column
  }
}



#Fix up the FDGates sample column
for(i in seq_along(FDGates)){
  if(!is.null(colnames(FDGates[[i]]))){
    FDGates[[i]]$Sample <- factor(gsub("1228","12828",as.numeric(gsub(".*(1349|1369|12828|1228).*","\\1",FDGates[[i]]$Sample))))
   }
}


## Drop Centers with missing data for manual gates
for(i in seq_along(CGates)){
  if(!is.null(colnames(CGates[[i]]))){
    CGates[[i]]$Center <- revalue(CGates[[i]]$Center,replace=c(BIIR="Baylor"))
    CGates[[i]] <- subset(CGates[[i]],Center%in%getLevels(FDGates,"Center")[[1]])
    CGates[[i]]$Center <- factor(CGates[[i]]$Center)
  }
}

pops <- cbind(getLevels(FDGates,"Population"),
getLevels(OCGates,"Population"),
getLevels(CGates,"Population"))
colnames(pops) <- c("FD","OC","C")



exchangeNV <- function(x){
  vals <- as.character(names(x))
  nms <- as.vector(do.call(c,x))
  l <- as.list(vals)
  names(l) <- nms
  l
}

## FD population names.. need to uniquely identify the population..
## Tcells
selected <- "Tcell"
pops[selected,]
popname <- as.character(FDGates[[selected]]$Phenotype) 
popname[(popname=="")] <- as.character(FDGates[[selected]]$Population[popname==""])
FDGates[[selected]]$Population<-factor(popname)
OCGates[[selected]]$Population <- factor(paste(OCGates[[selected]]$Parent,OCGates[[selected]]$Population))

pops <- cbind(getLevels(FDGates,"Population"),
getLevels(OCGates,"Population"),
getLevels(CGates,"Population"))
colnames(pops) <- c("FD","OC","C")

Tcell.central.map <- as.list(pops[[selected,"C"]])
Tcell.FD.map <- as.list(pops[[selected,"FD"]])
Tcell.OC.map <- as.list(pops[[selected,"OC"]])
names(Tcell.central.map) <- c("Lymphocytes","CD3","CD4","CD4 Activated","CD4 Naive","CD4 Central Memory","CD4 Effector Memory","CD4 Effector","CD8","CD8 Activated","CD8 Naive","CD8 Central Memory","CD8 Effector Memory", "CD8 Effector")
names(Tcell.FD.map) <- c("remove","CD3","CD4 Effector Memory","CD4 Effector","CD4 Central Memory","CD4 Naive","remove","remove","remove","CD4 Activated","CD8","CD4","CD8 Effector Memory","CD8 Effector","CD8 Central Memory","CD8 Naive","remove","remove","remove","CD8 Activated","remove","Lymphocytes","remove")
names(Tcell.OC.map) <- c("CD8","CD4","CD8 Effector Memory","CD8 Effector","CD8 Central Memory","CD8 Naive","CD4 Effector Memory","CD4 Effector","CD4 Central Memory","CD4 Naive","CD8 Activated","CD4 Activated","CD3","Lymphocytes")

#Remap Tcells
CGates[[selected]]$Population <-mapvalues(CGates[[selected]]$Population,from=do.call(c,Tcell.central.map),to=names(Tcell.central.map))
FDGates[[selected]]$Population <-
  mapvalues(FDGates[[selected]]$Population,from=do.call(c,Tcell.FD.map),to=names(Tcell.FD.map))
OCGates[[selected]]$Population <-
  mapvalues(OCGates[[selected]]$Population,from=do.call(c,Tcell.OC.map),to=names(Tcell.OC.map))


## Merge the Tcell data
#allCacheFiles <-
#  list.files(".", pattern="^tmp_R_cache_TCELLS\\.RData$", full.name=TRUE)
#         file.remove(allCacheFiles)
TCELLS <- subset(rbind(cbind(CGates[[selected]][,1:5],Method="Manual"),
                       cbind(FDGates[[selected]][,1:5],Method="flowDensity"),
                 cbind(OCGates[[selected]][,1:5],Method="OpenCyto")),!Population%in%"remove")

## Compute the CV for each sample and method and choose the best automated method
TCELL.CV <- ddply(TCELLS,.(Sample,Method,Population),summarize,CV=sd(Proportion,na.rm=TRUE)/mean(Proportion,na.rm=TRUE))
TCELL.CV <- melt(ddply(cast(TCELL.CV,value="CV",id=c("Sample","Population","Method"),Sample+Population~Method),.(Sample,Population),transform,Automated=min(flowDensity,OpenCyto)))
setnames(TCELL.CV,c("variable","value"),c("Method","CV"))

pdf("Comparison/TCELLS.pdf",width=15,height=5)
ggplot(subset(TCELL.CV,!(Population%in%"Lymphocytes")&Method%in%c("Manual","Automated")))+geom_bar(aes(y=CV,x=Method,fill=Method),stat="identity")+facet_grid(Sample~Population)+theme(axis.text.x=element_text(angle=90,hjust=1),strip.text.x=element_text(size=8))
dev.off()



## Tregs
selected <- "Treg"
pops[selected,]
OCGates[[selected]][,c("Population","Parent")]
Treg.central.map <- as.list(pops[[selected,"C"]])
Treg.FD.map <- as.list(pops[[selected,"FD"]])
Treg.OC.map <- as.list(pops[[selected,"OC"]])
names(Treg.central.map) <- c("Lymphocytes","CD3","CD4","Lo127Hi25","Naive","Memory","Total Treg","Activated","remove","remove","remove","remove","remove","remove","remove")
names(Treg.OC.map) <- c("Activated","remove","remove","remove","remove","remove","remove","remove","remove","Lo127Hi25","CD3","CD4","Lymphocytes","Memory","Naive","Total Treg")


## Remap Tregs
CGates[[selected]]$Population <- mapvalues(CGates[[selected]]$Population,from=do.call(c,Treg.central.map),to=names(Treg.central.map))
OCGates[[selected]]$Population <- mapvalues(OCGates[[selected]]$Population,from=do.call(c,Treg.OC.map),to=names(Treg.OC.map))
## Merge Tregs
TREG <- subset(rbind(cbind(CGates[[selected]][,1:5],Method="Manual"),
cbind(OCGates[[selected]][,1:5],Method="OpenCyto")),!Population%in%"remove")


## Compute CV for TREGS
TREG.CV <- ddply(TREG,.(Sample,Method,Population),summarize,CV=sd(Proportion,na.rm=TRUE)/mean(Proportion,na.rm=TRUE))

pdf("Comparison/TREG.pdf",width=15,height=5)
ggplot(subset(TREG.CV,!(Population%in%"Lymphocytes")&Method%in%c("Manual","Automated")))+geom_bar(aes(y=CV,x=Method,fill=Method),stat="identity")+facet_grid(Sample~Population)+theme(axis.text.x=element_text(angle=90,hjust=1),strip.text.x=element_text(size=8))
dev.off()


## Bcells
selected <- "Bcell"
pops[selected,]
Bcell.central.map <- as.list(pops[[selected,"C"]])
Bcell.FD.map <- as.list(pops[[selected,"FD"]])
Bcell.OC.map <- as.list(pops[[selected,"OC"]])
names(Bcell.central.map) <- c("Lymphocytes","CD3","CD19","remove","CD20","Naive","Memory IgD+","Memory IgD-","Transitional","Plasmablasts","remove","remove")
names(Bcell.FD.map) <- c("remove","CD19","CD20","CD3","remove","Memory IgD-","Naive","Memory IgD+","remove","Lymphocytes","Plasmablasts","remove","Transitional")
names(Bcell.OC.map) <- c("CD19","CD20","CD3","remove","Memory IgD-","Naive","Memory IgD+","Lymphocytes","Plasmablasts","Transitional")
## Remap
remap <- function(G,m){
  G$Population <- mapvalues(G$Population,from=do.call(c,m),to=names(m))
  G
}
## Remap Bcells
CGates[[selected]]$Population <- remap(CGates[[selected]],Bcell.central.map)$Population
OCGates[[selected]]$Population <- remap(OCGates[[selected]],Bcell.OC.map)$Population
FDGates[[selected]]$Population <- remap(FDGates[[selected]],Bcell.FD.map)$Population
## Merge Bcells
BCELL <- subset(rbind(cbind(CGates[[selected]][,1:5],Method="Manual"),
cbind(FDGates[[selected]][,1:5],Method="flowDensity"),
      cbind(OCGates[[selected]][,1:5],Method="OpenCyto")),!Population%in%c("remove","CD3")|is.na(Population))
BCELL$Population <- factor(BCELL$Population)
## Compute CV for Bcells
BCELL.CV <- ddply(BCELL,.(Sample,Method,Population),summarize,CV=sd(Proportion,na.rm=TRUE)/mean(Proportion,na.rm=TRUE))
BCELL.CV <- melt(ddply(cast(BCELL.CV,value="CV",id=c("Sample","Population","Method"),Sample+Population~Method),.(Sample,Population),transform,Automated=min(flowDensity,OpenCyto)))
setnames(BCELL.CV,c("variable","value"),c("Method","CV"))

pdf("Comparison/BCELLS.pdf",width=15,height=5)
ggplot(subset(BCELL.CV,!(Population%in%"Lymphocytes")&!(Method%in%"Automated")))+geom_bar(aes(y=CV,x=Method,fill=Method),stat="identity")+facet_grid(Sample~Population)+theme(axis.text.x=element_text(angle=90,hjust=1),strip.text.x=element_text(size=8))
dev.off()



## Monocytes
selected <- "DCMonoNK"
pops[selected,]
Mono.central.map <- as.list(pops[[selected,"C"]])
Mono.FD.map <- as.list(pops[[selected,"FD"]])
Mono.OC.map <- as.list(pops[[selected,"OC"]])


## Thelper
selected <- "Thelper"
## Fix the Thelper population names
FDGates[[selected]]$Population <- factor(paste(FDGates[[selected]]$Phenotype,FDGates[[selected]]$Population))
OCGates[[selected]]$Population <- factor(paste(OCGates[[selected]]$Parent,OCGates[[selected]]$Population))
pops <- cbind(getLevels(FDGates,"Population"),
getLevels(OCGates,"Population"),
getLevels(CGates,"Population"))
colnames(pops) <- c("FD","OC","C")

pops[selected,]
Thelper.central.map <- as.list(pops[[selected,"C"]])
Thelper.FD.map <- as.list(pops[[selected,"FD"]])
Thelper.OC.map <- as.list(pops[[selected,"OC"]])
names(Thelper.central.map) <- c("Lymphocytes","CD3","CD4","CD4 Activated","CD4 Th1","CD4 Th2","CD4 Th17","CD8","CD8 Activated","CD8 Th1","CD8 Th2","CD8 Th17")
names(Thelper.FD.map) <- c("remove","CD3","remove","remove","remove","CD4 Activated","CD4","CD8","CD4 Th2","CD4 Th17","CD4 Th1","remove","remove","remove","remove","CD8 Activated","CD8 Th2","CD8 Th17","CD8 Th1","remove","remove","Lymphocytes","remove")
names(Thelper.OC.map) <- c("CD8","CD4","remove","remove","remove","CD8 Activated","remove","remove","remove","CD4 Activated","CD8 Th2","CD8 Th17","CD8 Th1","remove","CD4 Th2","CD4 Th17", "CD4 Th1","remove","CD3","remove")

## Map
CGates[[selected]]$Population <- remap(CGates[[selected]],Thelper.central.map)$Population
OCGates[[selected]]$Population <- remap(OCGates[[selected]],Thelper.OC.map)$Population
FDGates[[selected]]$Population <- remap(FDGates[[selected]],Thelper.FD.map)$Population

## Merge
THELPER <- subset(rbind(cbind(CGates[[selected]][,1:5],Method="Manual"),
cbind(FDGates[[selected]][,1:5],Method="flowDensity"),
      cbind(OCGates[[selected]][,1:5],Method="OpenCyto")),!Population%in%c("remove")|is.na(Population))
THELPER$Population <- factor(THELPER$Population)


## Compute CV for Thelper
THELPER.CV <- ddply(THELPER,.(Sample,Method,Population),summarize,CV=sd(Proportion,na.rm=TRUE)/mean(Proportion,na.rm=TRUE))
THELPER.CV <- melt(ddply(cast(THELPER.CV,value="CV",id=c("Sample","Population","Method"),Sample+Population~Method),.(Sample,Population),transform,Automated=min(flowDensity,OpenCyto)))
setnames(THELPER.CV,c("variable","value"),c("Method","CV"))

pdf("Comparison/THELPER.pdf",width=15,height=5)
ggplot(subset(THELPER.CV,!(Population%in%"Lymphocytes")&!(Method%in%"Automated")))+geom_bar(aes(y=CV,x=Method,fill=Method),stat="identity")+facet_grid(Sample~Population)+theme(axis.text.x=element_text(angle=90,hjust=1),strip.text.x=element_text(size=8))
dev.off()

save("THELPER","TCELLS","BCELL","TREG",file="data/MergedTables.rda")
