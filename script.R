computerFlag="cluster"
computerFlag=""

if (computerFlag=="cluster") {
} else {
    dirSrc="/Users/royr/UCSF/"
    dirSrc2=dirSrc
    setwd(paste(dirSrc2,"JoeWiemels/leukMeth",sep=""))
}

####################################################################
####################################################################
## Section 1


cohort="_aml"
cohort="_leuk"
cohort="_allGuthSet1Set2"
#cohort="_tcgaGbm"
cohort="_birthDefect"
cohort="_allGuthSet1"
cohort="_allGuthSet2"

subsetFlag="_ctrl"
subsetFlag="_case"

subsetFlag=""
#subsetFlag="noHisp"
subsetFlag="hisp"
subsetFlag="noHispWt"

subsetName=ifelse(subsetFlag=="",subsetFlag,paste(subsetFlag,"Subset",sep=""))

### R code from vignette source 'vignettes/minfi/inst/doc/minfi.Rnw'

options(width=70)

if (F) {
    require(minfi)
    #require(minfiData)
    library(IlluminaHumanMethylation450kmanifest)
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
}

if (computerFlag=="cluster") {
    if (cohort=="_allGuthSet1") {
        baseDir="data/set1/idat/"
        targets=read.table("data/set1/clin_allGuthSet1_20160928.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        tbl=read.table("data/meth_yrbirth.csv",sep=",",h=T,quote="",comment.char="",as.is=T,fill=T)
        targets$yob=NA
        j=match(targets$subjectID,tbl$subjectid); j1=which(!is.na(j)); j2=j[j1]
        targets$yob[j1]=tbl$ch_birth_year[j2]
    } else if (cohort=="_allGuthSet2") {
        baseDir="data/set2/idat/"
        #targets=read.table("data/set2/clinGuthrieReplJune2012.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        targets=read.table("data/set2/clin_guthrieSet2_20140619.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        tbl=read.table("data/meth_yrbirth.csv",sep=",",h=T,quote="",comment.char="",as.is=T,fill=T)
        targets$yob=NA
        j=match(targets$subjectID,tbl$subjectid); j1=which(!is.na(j)); j2=j[j1]
        targets$yob[j1]=tbl$ch_birth_year[j2]
    } else if (cohort=="_allGuthSet1Set2") {
        baseDir="data/idat/"
        targets=read.table("data/set1set2/clin_guthrieSet1Set2_20140619.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        targets$id=paste("X",targets$guthrieId,sep="")
    } else if (cohort=="_leuk") {
        baseDir="data/set1/idat/"
        
        clin=read.delim("data/set1/i.LEU.v2.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        names(clin)[match(c("Plate","Sentrix_ID","Sentrix_Position"),names(clin))]=c("Batch","Beadchip","Position")
        clin$id=sapply(clin$Sample,function(x) {if (is.na(as.integer(substr(x,1,1)))) x else paste("X",x,sep="")},USE.NAMES=F)
        clin$group="Leukemia"
        clin$group2="Leuk"
        
        tbl1=read.table("data/0708011 Sample_Sheet (Fetal blood).csv",sep=",",h=T,quote="",comment.char="",as.is=T,fill=T,skip=14)
        names(tbl1)[match(c("Sample_Name","Sample_Well","Sample_Plate","Sample_Group","Pool_ID","Sentrix_ID","Sentrix_Position"),names(tbl1))]=c("id1","Sample_Well","Sample_Plate","group","Pool_ID","Beadchip","Position")
        tbl1$id=tbl1$id1
        j=which(!is.na(as.integer(tbl1$id1)))
        tbl1$group[j]=paste("S",tbl1$group[j],sep="")
        tbl1$group2=tbl1$group
        tbl1$group2[which(tbl1$group=="Allcell-B")]="B"
        tbl1$group2[which(tbl1$group=="Allcell-non-B")]="nB"
        tbl1$id=paste(tbl1$group2,"_",tbl1$id1,sep="")
        tbl1$sex=NA
        tbl1[tbl1$id%in%tbl1$id[duplicated(tbl1$id)],]
        
        names(clin)[names(clin)%in%names(tbl1)]
        #names(clin)[!names(clin)%in%names(tbl1)]
        #names(tbl1)[!names(tbl1)%in%names(clin)]
        
        k=match(names(clin),names(tbl1)); k1=which(!is.na(k)); k2=k[k1]
        targets=rbind(clin[,k1],tbl1[!duplicated(tbl1$id),k2])
        x=gsub("_+","_",gsub("(","_",gsub(" |-|)","_",targets$id),fixed=T))
        j=which(substr(x,nchar(x),nchar(x))=="_")
        if (length(j)!=0) x[j]=substr(x[j],1,nchar(x[j])-1)
        x[!(x%in%colnames(yMVals))]
        targets$id=x
        targets$Basename=paste("/",targets$Beadchip,"/",targets$Beadchip,"_",targets$Position,sep="")
        
        if (F) {
            targets$tmp=paste(targets$Beadchip,"_",targets$Position,sep="")
            targets$tmp=paste(targets$Beadchip,"/",targets$Beadchip,"_",targets$Position,"_Grn.idat",sep="")
            x=dir(path="data/set1/idat",pattern="_Grn.idat",recursive=T)
            table(targets$tmp%in%x)
            
            tbl1=read.table("data/0708011 Sample_Sheet (Fetal blood).csv",sep=",",h=T,quote="",comment.char="",as.is=T,fill=T,skip=14)
            names(tbl1)[match(c("Sample_Name","Sample_Well","Sample_Plate","Sample_Group","Pool_ID","Sentrix_ID","Sentrix_Position"),names(tbl1))]=c("id1","Sample_Well","Sample_Plate","group","Pool_ID","Beadchip","Position")
            tbl1$id=tbl1$id1
            j=which(!is.na(as.integer(tbl1$id1)))
            tbl1$id[j]=paste("X",tbl1$id1[j],"_S",tbl1$group[j],sep="")
            
            targets=tbl1
            targets$Basename=paste("/",targets$Beadchip,"/",targets$Beadchip,"_",targets$Position,sep="")
            
            targets$tmp=paste(targets$Beadchip,"/",targets$Beadchip,"_",targets$Position,"_Grn.idat",sep="")
            x=dir(path="data/set1/idat",pattern="_Grn.idat",recursive=T)
            table(targets$tmp%in%x)
        }
    } else if (cohort=="_birthDefect") {
        baseDir="data/birthDefect/idat/"
        targets=read.table("data/birthDefect/clin_birthDefect_20160607.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    } else if (cohort=="_tcgaGbm") {
        baseDir="data/tcgaGbm/idat/"
        x=dir(baseDir)
        x=x[grep("_Grn",x)]
        targets=as.data.frame(t(sapply(x,function(x) {
            y=sub("_Grn.idat","",x,fixed=T)
            y=strsplit(y,"_")[[1]]
            c(paste("X",y[1],"_",y[2],sep=""),y)
        },USE.NAMES=F)),stringsAsFactors=F)
        names(targets)=c("id","Beadchip","Position")
        rm(x)
    } else {
        baseDir="data/aml/I169_DataFiles_part4_110311/"
        datadir="data/"
        #tbl1=read.table("data/List of 110 pilot ANB specimens to Joe 12-15-2011.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        #tbl2=read.table(paste(datadir,"AML-Sample-Layout.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,skip=1)
        #tbl3=read.table(paste(datadir,"aml/i.AML.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        tbl4=read.table(paste(datadir,"aml/11072011 Sample_Sheet (AML).csv",sep=""),sep=",",h=T,quote="",comment.char="",as.is=T,fill=T,skip=14)
        names(tbl4)[match(c("Sample_Name","Sample_Well","Sample_Plate","Sample_Group","Pool_ID","Sentrix_ID","Sentrix_Position"),names(tbl4))]=c("guthrieID","Sample_Well","Sample_Plate","sampleGroup","Pool_ID","Beadchip","Position")
        targets=read.table(paste(datadir,"aml/i.AML.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        names(targets)=c("guthrieId","caco","sex","subjectID")
        j=match(targets$guthrieID,tbl4$guthrieID)
        targets=cbind(targets,tbl4[,c("sampleGroup","Beadchip","Position")])
    }
} else {
    if (cohort=="_allGuthSet1") {
        baseDir="docs/all/set1/idat/"
        #baseDir="docs/all/set1/idat/6057833042/"
        targets=read.table("docs/all/set1/clin_allGuthSet1_20160928.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        targets=targets[,which(!names(targets)%in%c("Position"))]
        names(targets)[match(c("Position1"),names(targets))]=c("Position")
    } else if (cohort=="_allGuthSet2") {
        baseDir="docs/all/set2/idat/"
        targets=read.table("docs/all/set2/clinGuthrieReplJune2012.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        tbl=read.table("docs/all/meth_yrbirth.csv",sep=",",h=T,quote="",comment.char="",as.is=T,fill=T)
        targets$yob=NA
        j=match(targets$subjectID,tbl$subjectid); j1=which(!is.na(j)); j2=j[j1]
        targets$yob[j1]=tbl$ch_birth_year[j2]
    } else {
    }
}
if (cohort%in%c("_allGuthSet2","_allGuthSet1Set2")) {
    # "0855G" - T-ALL
    targets=targets[which(!is.na(targets$guthrieId) & !targets$guthrieId%in%c("1458G","1762G","0635G","1588G","0855G")),]
}
targets$Basename=paste("/",targets$Beadchip,"/",targets$Beadchip,"_",targets$Position,sep="")
if (subsetFlag=="noHispWt") {targets=targets[which(targets$int_ch_ethnicity==2),]}

if (F) {
    #list.files(baseDir)
    #list.files(file.path(baseDir, "6042308046"))
    #targets=read.450k.sheet(baseDir)
    if (!cohort%in%c("_tcgaGbm")) {
        x=paste("/",sub("_Grn.idat","",list.files(baseDir,recursive=T,pattern="_Grn.idat"),fixed=T),sep="")
        table(targets$Basename%in%x)
        table(x%in%targets$Basename)
        targets=targets[targets$Beadchip%in%list.files(baseDir),]
    }
}

####################################################################
####################################################################

########Please keep the names of the graphs and tables as indicated in the script below because this ########will be used as a code to link tables/graphs to relevant scripts
########
#rm(list=ls())

### the libraries needed for the analysis
### please check that you have the latest versions, update if necessary
#library(BiocGenerics)
library(minfi)
#library(limma)
#library(quadprog)

################ Uncomment and run if necessary

#source("https://bioconductor.org/biocLite.R")
#biocLite("FlowSorted.Blood.450k")
#library(FlowSorted.Blood.450k)

#source("https://bioconductor.org/biocLite.R")
#biocLite("FlowSorted.CordBlood.450k")
library(FlowSorted.CordBlood.450k)

#source("https://bioconductor.org/biocLite.R")
#biocLite("IlluminaHumanMethylation450kmanifest")
library(IlluminaHumanMethylation450kmanifest)

#source("https://bioconductor.org/biocLite.R")
#biocLite("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)


#### Read RGSet file ####
#baseDir <- c("idats")              ##### provide the path to the idat folder
output_dir="output/"               ##### provide the path for results/outputs folder
#pdata <-read.metharray.sheet(baseDir)    ##### the csv file should be in the idat folder
pdata=targets
rm(targets)
##### factor character columns if it's not the case in pdata (e.g. pdata$Gender<-factor(pdata$Gender) )
pdata$Gender=as.character(pdata$sex)
pdata$Gender[which(pdata$sex==1)]="M"
pdata$Gender[which(pdata$sex==2)]="F"

#RGSetraw <- read.metharray.exp(base = baseDir, targets = pdata)
load("rgsetRaw_allGuthSet2.RData")
RGSetraw <- rgsetRaw
rm(rgsetRaw)
#RGSetraw <- read.metharray.exp(base = baseDir, targets = pdata, recursive=T)
####!!!##### remove samples from RGSetraw after “QC Steps”
####### RGSetraw<-RGsetcounts[,-c(“low quality samples; use their indexes”)]   ###We call this Line 1 ####### for later reference.

j=match(sampleNames(RGSetraw),paste(pdata$Beadchip,"_",pdata$Position,sep="")); j1=which(!is.na(j)); j2=j[j1]
RGSetraw=RGSetraw[,j1]
pdata=pdata[j2,]
#save(RGSetraw,file=paste("RGSetraw",cohort,".RData",sep=""))

colId1=c("caco","Gender","yob","ch_ageref","nRBC","CD8T","CD4T", "NK","Bcell","Mono","Gran",paste("epistr",1:5,sep=""),"Plate","Position")
colId1=c("caco","Gender","yob","ch_ageref","Plate","Position")
j=1:nrow(pdata)
for (k in colId1) {
    j=j[which(!is.na(pdata[j,k]))]
}
RGSetraw=RGSetraw[,j]
pdata=pdata[j,]

########### estimate cell counts
counts <- estimateCellCounts(RGSetraw,compositeCellType = "CordBlood",cellTypes=c("nRBC","CD8T","CD4T", "NK","Bcell","Mono","Gran"),meanPlot=FALSE)
save(counts, file=paste0(output_dir,"Test.CordBlood.CellType.nRBC.Rdata"))

#####
rawbetas <- getBeta(RGSetraw)

probNAcont <- which(apply(rawbetas,1,function(i) sum(!is.finite(i)))>0)   ##### check for NAs and save for later

write.table(rawbetas, file=paste(output_dir, "betas.csv",sep="/"), row.names=T, sep="\t")

save.image("tmp_1.RData")

##################QC Steps
######change aesthetics (e.g. sampGroups, sampNames) according to the specific need

MSet <- preprocessRaw(RGSetraw)

save.image("tmp_2.RData")

estimate_mode <- function(x) {
    d <- density(x,na.rm=T)
    d$x[which.max(d$y)]
}
x=apply(rawbetas,2,median,na.rm=T)
x2=apply(rawbetas,2,estimate_mode)
print('pdata[,c("guthrieId","Beadchip","Position")][(x2>.2 & x2<.8),]')
print(pdata[,c("guthrieId","Beadchip","Position")][(x2>.2 & x2<.8),])
"
guthrieId   Beadchip Position
417     1762G 9761749095   R05C01
"
x=rep("Good",length(x2)); x[(x2>.2 & x2<.8)]="Bad"

qc <- getQC(MSet)
jpeg(paste(output_dir, "plotQC.jpg",sep="/"), width=800, height=800)
plotQC(qc,badSampleCutoff=11)
invisible(dev.off())

print('pdata[,c("guthrieId","Beadchip","Position")][c(182,212:214),]')
print(pdata[,c("guthrieId","Beadchip","Position")][c(182,212:214),])
"
guthrieId   Beadchip Position
182     1458G 9702496074   R01C02
212     0911G 9702496086   R04C02
213     0380G 9702496086   R05C01
214     1338G 9702496086   R05C02
"

jpeg(paste(output_dir, "methylated.jpg",sep="/"), width=800, height=800)
boxplot(log(RGSetraw@assayData$Red+1)[,order(colnames(rawbetas))], main="methylated (Red channel)", las=2, cex.axi=0.8)
invisible(dev.off())

jpeg(paste(output_dir, "unmethylated.jpg",sep="/"), width=800, height=800)
boxplot(log(RGSetraw@assayData$Green+1)[,order(colnames(rawbetas))], main="unmethylated (Green channel)", las=2, cex.axi=0.8)
invisible(dev.off())

jpeg(paste(output_dir, "densityPlot.jpg",sep="/"), width=800, height=800)
#densityPlot(rawbetas, sampGroups = pdata$Cohort)
densityPlot(rawbetas, sampGroups = x)
invisible(dev.off())

jpeg(paste(output_dir, "densityBeanPlot.jpg",sep="/"), width=800, height=800)
#densityBeanPlot(rawbetas, sampGroups = pdata$Cohort)
densityBeanPlot(rawbetas, sampGroups = x)
invisible(dev.off())


jpeg(paste(output_dir, "mdsPlot.jpg",sep="/"), width=800, height=800)
if (dim(pdata)[1] < 50 ) {
    mdsPlot(rawbetas,numPositions=1000,sampGroups=as.character(pdata$Gender),
    legendPos = "bottomleft",sampNames=1:215)
} else {
    mdsPlot(rawbetas,numPositions=1000,sampGroups=as.character(pdata$Gender))
}
invisible(dev.off())

####if you have bad quality samples, go back to beginning and re-run all the script from Line 1 onwards

########## check Gender
GMsetEx <- mapToGenome(RGSetraw)
estSex <- getSex(GMsetEx)

#pdata$Gender <- ifelse(pdata$Gender=="Male","M","F")     #####  code gender by “M” and “F” and not by “1” and “2” or anything else.
#print("sum(estSex$predictedSex!=pdata$Gender)")
#print(sum(estSex$predictedSex!=pdata$Gender))  ###  this will indicate the number of mismatches
print("table(predSex=estSex$predictedSex,ObsSex=pdata$Gender,exclude=NULL)")
print(table(predSex=estSex$predictedSex,ObsSex=pdata$Gender,exclude=NULL))  ###  this will indicate the number of mismatches
## Exclude 0635G, 1588G

#######mdsPlots are done elsewhere in the code, before and after normalization, to further assess for #######gender mismatches.  In case mismatches exists, compare mdsPlots generated using the genders #######reported in the questionnaire data versus those predicted by getsex.

####################### Normalize with funnorm
Gender <- which(colnames(pdata)=="Gender")

suppressWarnings(RGset <- preprocessFunnorm(RGSetraw, sex=pdata[,Gender]))

########################remove unwanted probes (cross reactive probes)
crp_probes <- "data/Crossreactive_probes.csv"       ###add the path to the csv file
cross_reactive <- read.table(crp_probes, header=F, sep="\t")$V1
RGset <- RGset[ ! featureNames(RGset)%in%cross_reactive,]

####################### remove probes with too many NAs; otherwise, replace by the mean value of the other samples for the respective probes
CpGNAexcl = function(betas,CpGlimit=0.05) {
    NArow <- apply(betas,1,function(i) sum(is.na(i)))
    exclprob <- rownames(betas)[NArow>CpGlimit*ncol(betas)]
    return(exclprob)
}

betas <- getBeta(RGset)
probNAcont <- probNAcont[names(probNAcont)%in%setdiff(names(probNAcont), as.character(cross_reactive))]
for (i in names(probNAcont)) {
    betas[i,is.na(rawbetas[i,])] <- mean(rawbetas[i,],na.rm=T)
}
remunwprob <- CpGNAexcl(rawbetas)    ##### threshold is 5% missing values
betas <- betas[-which(rownames(betas)%in%remunwprob),]

x2=apply(betas,2,estimate_mode)
print('pdata[,c("guthrieId","Beadchip","Position")][(x2>.2 & x2<.8),]')
print(pdata[,c("guthrieId","Beadchip","Position")][(x2>.2 & x2<.8),])
x=rep("Good",length(x2)); x[(x2>.2 & x2<.8)]="Bad"
if (F) {
    x[which(pdata$guthrieId%in%c("1762G","0635G","1588G"))]="Ugly"
    x[which(pdata$guthrieId%in%c("1458G","0911G","0380G","1338G"))]="Ugly2"
    j=which(pdata$guthrieId%in%c("1458G","0911G","0380G","1338G"))
    x[j]=pdata$guthrieId[j]
    x[which(is.na(pdata$guthrieId))]="NA"
    x=rep("Good",length(x2))
    j=which(pdata$guthrieId%in%c("1458G","0911G","0380G","1338G"))
    x[j]=pdata$guthrieId[j]
}

#######################save
write.table(betas, file=paste(output_dir, "betas.filtered.csv",sep="/"), row.names=T, sep="\t")
write.table(pdata, file=paste(output_dir, "pdata.txt",sep="/"), row.names=F, sep="\t")

######some QC graphs after normalization
jpeg(paste(output_dir, "densityPlot2.jpg",sep="/"), width=800, height=800)
#densityPlot(betas, sampGroups = pdata$Sample_Plate)
densityPlot(betas, sampGroups = x)
invisible(dev.off())

jpeg(paste(output_dir, "densityBeanPlot2.jpg",sep="/"), width=800, height=800)
densityBeanPlot(betas, sampGroups = pdata$Cov1)
invisible(dev.off())

jpeg(paste(output_dir, "mdsPlot2.jpg",sep="/"), width=800, height=800)
if (dim(pdata)[1] < 50 ) {
    mdsPlot(betas,numPositions=1000,sampGroups=pdata$Gender,
    legendPos = "bottomleft",sampNames=1:215)
} else {
    mdsPlot(betas,numPositions=1000,sampGroups=pdata$Gender)
}
invisible(dev.off())


save.image("tmp_3.RData")

############################## Examples on how to factor variables, modify variables, modify variable names..etc  if necessary
#colnames(pdata)[14]<-"Sentrix_Position"
#colnames(pdata)[15]<-"Sentrix_ID"
#pdata$Gender<-ifelse(pdata$Gender=="Male","M","F")
#pdata<-merge(pdata,yearofBdata,by="RetrievalDetail_ID",sort=F)
#for (i in c(3,5,10,14,15,25))
#{
#  pdata[,i]<-factor(pdata[,i])
#}
#rm(i)

if (!subsetFlag%in%c("hisp","noHispWt")) {
    fName1="_moba"
    dirClin="data/"
    clin2=read.table(paste(dirClin,"epistructure",cohort,fName1,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    id=pdata$id[!pdata$id%in%clin2$id]
    if (length(id)!=0) {
        tmp=clin2[1:length(id),]
        for (k in 1:ncol(tmp)) tmp[,k]=NA
        tmp$id=id
        clin2=rbind(clin2,tmp)
    }
    pdata=cbind(pdata,clin2[match(paste(pdata$Beadchip,pdata$Position,sep="_"),clin2$id),grep("epistr",names(clin2))])
    rm(clin2,tmp,id)
}

for (k in 1:ncol(pdata)) {
    if (is.factor(pdata[,k])) pdata[,k]=as.character(pdata[,k])
}
pdata$caco=as.character(pdata$caco)
pdata$Plate=as.character(pdata$Plate)
for (k in which(names(pdata)%in%c("caco","Gender","yob","ch_ageref","nRBC","CD8T","CD4T", "NK","Bcell","Mono","Gran",paste("epistr",1:5,sep=""),"Plate","Position"))) {
    #for (k in 1:ncol(pdata)) {
    if (is.factor(pdata[,k])) pdata[,k]=as.character(pdata[,k])
    #if (is.numeric(pdata[,k])) {
    cat(colnames(pdata)[k],": ",class(pdata[,k]),"\n")
    #}
}


## ----------------------------------------------
if (T) {
    if (computerFlag=="") {
        load(file="ann.RData")
    } else {
        if (computerFlag=="cluster") {
            ann=read.delim(paste("data/","HumanMethylation450_15017482_v.1.2.csv",sep=""),header=TRUE, sep=",",quote="",comment.char="",as.is=T,fill=T, skip=7)
            #snpVec=read.table(paste("data/CpGs to exclude_FINAL.txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
            snpVec=read.table(paste("data/list_to_exclude_Sept_24.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        } else {
            ann=read.delim(paste("docs/yuanyuan/HumanMethylation450_15017482_v.1.2.csv",sep=""),header=TRUE, sep=",",quote="",comment.char="",as.is=T,fill=T, skip=7)
            #snpVec=read.table(paste("docs/ShwetaChoudhry/CpGs to exclude_FINAL.txt",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
            snpVec=read.table(paste("docs/SemiraGonsethNussle/list_to_exclude_Sept_24.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        }
        ann[which(ann[,"CHR"]=="X"),"CHR"]="23"
        ann[which(ann[,"CHR"]=="Y"),"CHR"]="24"
        ann[,"CHR"]=as.integer(ann[,"CHR"])
        ann=ann[,-match(c("AddressA_ID","AlleleA_ProbeSeq","AddressB_ID","AlleleB_ProbeSeq", "Next_Base",  "Color_Channel","Forward_Sequence","SourceSeq"),colnames(ann))]
        for (k in 1:ncol(ann)) if (class(ann[,k])=="factor") ann[,k]=as.character(ann[,k])
        i=match(rownames(betas),ann[,"IlmnID"])
        table(is.na(i))
        ann=ann[i,]
        
        snpVec=snpVec[,1]
        ann$snp=0; ann$snp[which(ann$IlmnID%in%snpVec)]=1
    }
    
    
    ann$geneSym=sapply(toupper(ann$UCSC_RefGene_Name),function(x) {
        strsplit(x,";")[[1]][1]
    },USE.NAMES=F)
    ann$geneSym[is.na(ann$geneSym)]=""
    
    #keep=!ann$CHR%in%c(23,24)
    #keep=!ann$CHR%in%c(23,24) & apply(meth,1,function(x) {any(!is.na(x))})
    #keep=ann$snp==0 & !ann$CHR%in%c(23,24) & apply(meth,1,function(x) {any(!is.na(x))})
    #keep=ann$snp==0 & ann$CHR%in%1:22 & apply(meth,1,function(x) {any(!is.na(x))})
}
betas <- betas[which(ann$CHR%in%1:22),]

#############################
#############################
########################### PCs before SVA
library(lumi)
library(sva)

lgBetas <- t(betas)
nPCs <- 10

invlgBetas <- 1/lgBetas
invlgBetas[!is.finite(invlgBetas)] <- 0
pca <- prcomp(invlgBetas,scale=T)

# plot proportion of var explained
propVarExpl <- summary(pca)$importance[2,]

jpeg(paste(output_dir, "componentsfilteredbeforeSVA.jpg",sep="/"))    ###PC plot
barplot(propVarExpl[1:10],ylab='propn var expl',xlab='PC')
invisible(dev.off())
#pdataBatch <- pdata[,c(3,15,14,5,10,25,17:23)]    #### put indexes for technical and biological variables
pdataBatch <- cbind(pdata,counts)
pdataBatch <- pdataBatch[,which(names(pdataBatch)%in%c("caco","Gender","yob","ch_ageref","nRBC","CD8T","CD4T", "NK","Bcell","Mono","Gran",paste("epistr",1:5,sep=""),"Plate","Position"))]    #### put indexes for technical and biological variables

variables <- colnames(pdataBatch)
# assemble potential batch covariates plus other covs of interest
nVar = length(variables)
batch_data <- data.frame(row.names = rownames(pdataBatch))
variables <- colnames(pdataBatch)
for (i in variables) {
    if (is.numeric(pdataBatch[,i]) && length(levels(as.factor(pdataBatch[,i]))) < 10 ) {
        pdataBatch[,i] <- as.factor(pdataBatch[,i])
    }
    batch_data <- cbind(batch_data,pdataBatch[i])
}
batch_data <- cbind(batch_data,pca$x[,1:nPCs])
covs <- names(batch_data)[1:nVar] # covs to consider
PCnames <- names(batch_data)[(nVar+1):length(batch_data)]

# matrix for results (pVals)
batch_corr_pVals <- matrix(nrow=length(covs),ncol=nPCs,dimnames=list(covs,PCnames))
#batch_data$Gender <- factor(batch_data$Gender)
for (name in PCnames) {
    # two factor vars - wilcox signed rank sum tesst
    for (variable in variables) {
        if (is.character(batch_data[,variable])) {
            if (length(levels(as.factor(sort(batch_data[,variable])))) > 1) {
                if (length(levels(as.factor(sort(batch_data[,variable])))) == 2) {
                    batch_corr_pVals[variable,name] = wilcox.test(get(name)~batch_data[, variable],data=batch_data)$p.value
                } else {
                    batch_corr_pVals[variable,name] = kruskal.test(get(name)~as.factor(batch_data[, variable]),data=batch_data)$p.value
                }
            }
        } else if (is.numeric(batch_data[,variable])) {
            #print(variable)
            fit=lm(get(name)~batch_data[, variable],data=batch_data)
            batch_corr_pVals[variable,name] = summary.lm(fit)$coefficients[2,4]
        } else {
            cat(variable,": ",class(batch_data[,variable]),"\n")
        }
    }
}

batch <- list(pVals=batch_corr_pVals, data=batch_data)
################## Tables of PCs*covariates
write.table(cbind(variable=rownames(batch$pVals),batch$pVals), file=paste(output_dir,"batchfilteredbeforeSVA.txt",sep="/"), row.names=T, sep="\t")

########################## SVA
meth <- beta2m(betas)    #### transformation from beta-values to m-values
meth[!is.finite(meth)] <- min(meth[is.finite(meth)])
#variables <- colnames(pdata)[c(10,3,14)]   ######variable of interest (Case/Control) + technical covariates
variables <- c("caco","Plate","Position")   ######variable of interest (Case/Control) + technical covariates

formula1 <- as.formula(paste("~ " ,paste(variables,collapse="+")))
design <- model.matrix(formula1,data=pdata)

#formula0 <- as.formula(paste("~ 0+", paste(variables[-1],collapse="+")))
formula0 <- as.formula(paste("~ ", paste(variables[-1],collapse="+")))
model0 <- model.matrix(formula0, data=pdata)
sva <- sva(meth,design,model0)

meth_adjSVA <- t(residuals(lm(t(meth)~sva$sv)))
betas_adjSVA <- m2beta(meth_adjSVA)   #### back transformation from m-values to beta-values

########################### PCs after SVA
lgBetas <- t(betas_adjSVA)
nPCs <- 10

invlgBetas <- 1/lgBetas
invlgBetas[!is.finite(invlgBetas)] <- 0
pca <- prcomp(invlgBetas,scale=T)

# plot proportion of var explained
propVarExpl <- summary(pca)$importance[2,]

jpeg(paste(output_dir, "componentsfilteredpostSVA.jpg",sep="/"))
barplot(propVarExpl[1:10],ylab='propn var expl',xlab='PC')
invisible(dev.off())
#pdataBatch <- pdata[,c(3,15,14,5,10,25,17:23)]    ##### use same variable indexes in PCs as before SVA
pdataBatch <- cbind(pdata,counts)
pdataBatch <- pdataBatch[,which(names(pdataBatch)%in%c("caco","Gender","yob","ch_ageref","nRBC","CD8T","CD4T", "NK","Bcell","Mono","Gran",paste("epistr",1:5,sep=""),"Plate","Position"))]    #### put indexes for technical and biological variables

variables <- colnames(pdataBatch)
# assemble potential batch covariates plus other covs of interest
nVar <- length(variables)
batch_data <- data.frame(row.names = rownames(pdataBatch))
variables <- colnames(pdataBatch)
for (i in variables) {
    if (is.numeric(pdataBatch[,i]) && length(levels(as.factor(pdataBatch[,i]))) < 10 ){
        pdataBatch[,i]<-as.factor(pdataBatch[,i])
    }
    batch_data <- cbind(batch_data,pdataBatch[i])
}
batch_data <- cbind(batch_data,pca$x[,1:nPCs])
covs <- names(batch_data)[1:nVar] # covs to consider
PCnames <- names(batch_data)[(nVar+1):length(batch_data)]

# matrix for results (pVals)
batch_corr_pVals <- matrix(nrow=length(covs),ncol=nPCs,dimnames=list(covs,PCnames))
#batch_data$Gender <- factor(batch_data$Gender)
for (name in PCnames) {
    # two factor vars - wilcox signed rank sum tesst
    for (variable in variables) {
        if (is.character(batch_data[,variable])) {
            if (length(levels(as.factor(sort(batch_data[,variable])))) > 1) {
                if (length(levels(as.factor(sort(batch_data[,variable])))) == 2) {
                    batch_corr_pVals[variable,name] <- wilcox.test(get(name)~batch_data[, variable],data=batch_data)$p.value
                } else {
                    batch_corr_pVals[variable,name] <- kruskal.test(get(name)~as.factor(batch_data[, variable]),data=batch_data)$p.value
                }
            }
        } else if (is.numeric(batch_data[,variable])) {
            print(variable)
            fit <- lm(get(name)~batch_data[, variable],data=batch_data)
            batch_corr_pVals[variable,name] <- summary.lm(fit)$coefficients[2,4]
        } else {
            cat(variable,": ",class(batch_data[,variable]))
        }
    }
}

batch <- list(pVals=batch_corr_pVals, data=batch_data)
write.table(batch$pVals, file=paste(output_dir,"batchfilteredpostSVA.txt",sep="/"), row.names=T, sep="\t")

##########################################
############## e.g. of choice of indexes for exclusion
# toremove<- which(phenotypes$sex!=phenotypes$Gender)



#################### Regression Models


library(data.table)# to process results
library(MASS) # rlm function for robust linear regression
library(sandwich) #Huberís estimation of the standard error
library(lmtest) # to use coeftest
#library(parallel) # to use multicore approach - part of base R

#Function
removeOutliers<-function(probes) {
    require(matrixStats)
    if(nrow(probes) < ncol(probes)) warning("expecting probes are rows (long dataset)")
    rowIQR <- rowIQRs(probes, na.rm = T)
    row2575 <- rowQuantiles(probes, probs = c(0.25, 0.75), na.rm = T)
    maskL <- probes < row2575[,1] - 3 * rowIQR
    maskU <- probes > row2575[,2] + 3 * rowIQR
    initial_NAs <- rowSums(is.na(probes))
    probes[maskL] <- NA
    removed_lower <- rowSums(is.na(probes))-initial_NAs
    probes[maskU] <- NA
    removed_upper <- rowSums(is.na(probes))-removed_lower-initial_NAs
    N_for_probe <- rowSums(!is.na(probes))
    Log <- data.frame(initial_NAs,removed_lower,removed_upper,N_for_probe)
    return(list(probes, Log))
}


meth <- betas_adjSVA   #use betas after SVA
#Remove outliers from METH (methylation data where probes are rows and samples are columns)
system.time(OutlierResults <- removeOutliers(meth))
meth.2 <- OutlierResults[[1]]
Log <- OutlierResults[[2]]
save(Log,file="output/Outlier_logBW.Rdata") #save log

tdat <- t(meth.2)     ##### this matrix will be used for the regression models

save(tdat,pdata,counts,file="data.RData")
write.table(pdata,file=paste("pdata.txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
write.table(cbind(id=rownames(counts),counts),file=paste("counts.txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
write.table(cbind(cpgId=rownames(meth.2),as.data.frame(meth.2)),file=paste("meth.2.txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
#write.table(cbind(cpgId=rownames(meth.2)[1:101],as.data.frame(meth.2[1:101,])),file=paste("meth.2.txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)

meth.2 <- t(tdat)     ##### this matrix will be used for the regression models
write.table(cbind(cpgId=rownames(meth.2),as.data.frame(meth.2)),file=paste("meth.2.txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)


######################################################### SUMMARIZE PROBES

descriptives <- function(x) {
    tmp <- c(min(x,na.rm=T),quantile(x,probs=c(.1,.25,.5),na.rm=T),mean(x,na.rm=T), median(x,na.rm=T),sd(x,na.rm=T),quantile(x,probs=c(.75,.90),na.rm=T),max(x,na.rm=T),sum(is.na(x)))
    names(tmp)[c(1,5:7,10:11)]<-c("Min.","Mean","Median","SD","Max.","NA")
    return(tmp)
}

desc <- t(apply(tdat,2,descriptives))
write.table(desc, file = paste(output_dir,"/MOBA3_Descriptives.txt", sep = ""),
sep = "\t", col.names = T, row.names = T, append = F, quote=FALSE)

####################################################################
####################################################################
## Section 2

## Run linReg_moba.R with linReg_moba.txt script

#######################
#######################
## Section 3
## Annotation

if (F) {
    load("ann.RData")
    datadir=""; fileList="tmp.txt"
    load("annAll.RData")
    datadir="results/moba/output/stat/"
    fileList=dir(datadir,pattern="ind.res_model")
    #fileList=fileList[1]
    f1=read.table(paste(datadir,fileList[1],sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    if (length(fileList)>1) {
        k=c()
        for (fId in 2:length(fileList)) {
            f2=read.table(paste(datadir,fileList[fId],sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
            if (any(is.na(f1$probeID)!=is.na(f2$probeID)) | any(f1$probeID!=f2$probeID,na.rm=T)) {
                k=c(k,fId)
            }
        }
        cat("Files with diff lists of loci: ",fileList[k],"\n\n",sep="")
    }
    
    table(is.na(match(f1$probeID,ann$IlmnID)))
    annThis=ann[match(f1$probeID,ann$IlmnID),]
    
    
    library(FDb.InfiniumMethylation.hg19)
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    
    ###### tes is a data.frame for annotating the CpGs to the nearest gene. It should have at least 2 columns, one with CpG ID and another with nearest gene name.
    genome="hg19"
    getTxdb <- function(genome="hg19") {
        if      (grepl("hg",genome)) { species="Hsapiens"
        } else if (grepl("mm",genome)) { species="Mmusculus"
        } else { return ("no reference") }
        
        suppressMessages(pkg<-paste("TxDb", species, "UCSC", genome ,"knownGene", sep="."))
        print(pkg)
        library(pkg, character.only = TRUE)
        txdb <- eval(parse(text = pkg))
        return(txdb)
    }
    txdb <- getTxdb(genome)
    columns(txdb)
    #genes <- genes(txdb)
    #(ranges(genes)@start+ranges(genes)@width-1)
    genes <- genes(txdb,columns=c("TXNAME","TXSTART","TXEND"))
    
    datadir=""
    datadir="results/moba/misc/"
    ann2 <- read.table(paste(datadir,"annotation_ucsc_hg19_knownGene.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,skip=1)
    names(ann2) <- c("ucscId","geneSym")
    #table(unique(unlist(res$TXNAME))%in%ann2$ucscId)
    geneSym=rep("",length(genes))
    length(geneSym)
    for (i in 1:length(genes$TXNAME)) {
        if (i%%5000==0) print(i)
        k=match(genes$TXNAME[[i]],ann2$ucscId); k=k[!is.na(k)]
        if (length(k)!=0) geneSym[i]=paste(unique(ann2$geneSym[k]),collapse=" ")
    }
    #########################
    
    
    hm <- get450k()
    overlap <- nearest(hm,genes)
    
    #tes <- data.frame(probeID=hm@ranges@NAMES,nearestgene=genes[overlap]$gene_name)
    #tesm <- merge(all.results1f,tes,by="probeID",sort=F)
    tes <- data.frame(IlmnID=hm@ranges@NAMES,nearestgene=geneSym[overlap])
    tesm <- merge(annThis,tes,by="IlmnID",sort=F)
    N <- 50
    N <- 5
    ndx <- order(tesm$P_VAL)[1:N]
    tesm[ndx,]
    i=which(tesm$geneSym!="")
    mean(tesm$geneSym[i]==tesm$nearestgene[i])
    save(tesm,file="tesm.RData")
}

#######################
#######################

if (T) {
    if (F) {
    cat("\n\n============================ ",modelId,"\n\n",sep="")


    ## ------------------------
    nProbe=10001
    nProbe=-1

    load("data.RData")
        tdat=read.table(paste("meth.2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=nProbe)
        rownames(tdat)=tdat$cpgId
        tdat=tdat[,-1]
        tdat=as.matrix(tdat)
        tdat=t(tdat)
        
        pdata=read.table(paste("pdata.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        
        counts=read.table(paste("counts.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        rownames(counts)=counts$id
        counts=as.matrix(counts[,-1])
        
        save(tdat,pdata,counts,file="data.RData")
    }

    betas=t(tdat)
    rm(tdat,counts)

    ## ------------------------
    if (F) {
        library(FDb.InfiniumMethylation.hg19)
        library(TxDb.Hsapiens.UCSC.hg19.knownGene)

        datadir=""
        datadir="results/moba/output/stat/"
        #all.results1=read.table(paste("ind.res_model",gsub("-","_",modelId),".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        all.results1=read.table(paste(datadir,"ind.res_model",gsub("-","_",modelId),".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)


        output_dir="tmp"

        ## Calculate lambda
        lambda <- qchisq(median(all.results1$P_VAL,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)

        jpeg(paste(output_dir, "/qqplotmodel",gsub("-","_",modelId),".jpg",sep=""))    #### change qqplot name according to regression #### model code
        #suppressWarnings(qq(all.results1$P_VAL,main=paste("QQ plot: ","lambda=",lambda,sep="")))
        suppressWarnings(qqnorm(all.results1$P_VAL,main=paste("QQ plot: ","lambda=",lambda,sep="")))
        invisible(dev.off())

        qvalues <- p.adjust(all.results1$P_VAL,method="BH")
        pvBonf <- p.adjust(all.results1$P_VAL,method="bonferroni")
        qvaluesInter <- p.adjust(all.results1$Inter_P_VAL,method="BH")
        group1 <- which(pdata$caco=="1")   ### define first group of
        group2 <- which(pdata$caco=="0")
        # selection <- sort(c(group1,group2))
        deltabeta <- rowMeans(betas[,group1],na.rm = T)-rowMeans(betas[,group2],na.rm = T) #(Case - Control)

        group1M <- which(pdata$caco=="1" & pdata$Gender=="M")   ### define first group of
        group2M <- which(pdata$caco=="0" & pdata$Gender=="M")
        deltabetaM <- rowMeans(betas[,group1M],na.rm = T)-rowMeans(betas[,group2M],na.rm = T) #(Case - Control)

        group1F <- which(pdata$caco=="1" & pdata$Gender=="F")   ### define first group of
        group2F <- which(pdata$caco=="0" & pdata$Gender=="F")
        deltabetaF <- rowMeans(betas[,group1F],na.rm = T)-rowMeans(betas[,group2F],na.rm = T) #(Case - Control)
        all.results1f <- cbind(all.results1,qvalues,pvBonf,qvaluesInter,deltabeta,deltabetaM,deltabetaF)

        ###### tes is a data.frame for annotating the CpGs to the nearest gene. It should have at least 2 columns, one with CpG ID and another with nearest gene name.
        genome="hg19"
        getTxdb <- function(genome="hg19") {
            if      (grepl("hg",genome)) { species="Hsapiens"
            } else if (grepl("mm",genome)) { species="Mmusculus"
            } else { return ("no reference") }
            
            suppressMessages(pkg<-paste("TxDb", species, "UCSC", genome ,"knownGene", sep="."))
            print(pkg)
            library(pkg, character.only = TRUE)
            txdb <- eval(parse(text = pkg))
            return(txdb)
        }
        columns(txdb)
        txdb <- getTxdb(genome)
        #genes <- genes(txdb)
        #(ranges(genes)@start+ranges(genes)@width-1)
        genes <- genes(txdb,columns=c("TXNAME","TXSTART","TXEND"))

        datadir="results/moba/misc/"
        datadir=""
        ann2 <- read.table(paste(datadir,"annotation_uscs_hg19_knownGene.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,skip=1)
        names(ann2) <- c("ucscId","geneSym")
        table(unique(unlist(genes$TXNAME))%in%ann2$ucscId)
        geneSym=rep("",length(genes))
        length(geneSym)
        for (i in 1:length(genes$TXNAME)) {
            if (i%%5000==0) print(i)
            k=match(genes$TXNAME[[i]],ann2$ucscId); k=k[!is.na(k)]
            if (length(k)!=0) geneSym[i]=paste(unique(ann2$geneSym[k]),collapse=" ")
        }
        #########################


        hm <- get450k()
        overlap <- nearest(hm,genes)

        #tes <- data.frame(probeID=hm@ranges@NAMES,nearestgene=genes[overlap]$gene_name)
        tes <- data.frame(probeID=hm@ranges@NAMES,nearestgene=geneSym[overlap])
    }
    tes=readRDS("probeannotation.rds")
    tesm <- merge(all.results1f,tes,by="probeID",sort=F)
    save(tesm,file="tesm.RData")
    N <- 50
    N <- 5
    ndx <- order(tesm$P_VAL)[1:N]
    tesm[ndx,]
}

#########################
## SNP
## Not yet working

library(BSgenome.Hsapiens.UCSC.hg19.masked)
#library(SNPlocs.Hsapiens.dbSNP142.GRCh37)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

lib="SNPlocs.Hsapiens.dbSNP142.GRCh37"
lib="SNPlocs.Hsapiens.dbSNP144.GRCh37"

installed.SNPs()
available.SNPs()
genome <- BSgenome.Hsapiens.UCSC.hg19.masked
SNPlocs_pkgname(genome)
genome2 <- injectSNPs(genome, lib)
genome2 # note the extra "with SNPs injected from ..." line
SNPlocs_pkgname(genome2)
snpcount(genome2)
#head(snplocs(genome2, "chr1"))
res=snplocs(genome, "chr22")
res2=snplocs(genome2, "chr22")


i=match(paste("rs",res$RefSNP_id,sep=""),ann$Probe_SNPs); i2=which(!is.na(i)); i1=i[i2]
summary(ann$MAPINFO[i1]-res$loc[i2])
i=match(paste("rs",res$RefSNP_id,sep=""),ann$Probe_SNPs_10); i2=which(!is.na(i)); i1=i[i2]
summary(ann$MAPINFO[i1]-res$loc[i2])
i=1:3
cbind(ann[i1[i],c("IlmnID","CHR","MAPINFO","Probe_SNPs")],res[i2[i],])


chrInfo=data.frame(chr=1:24,chrName=paste("chr",1:24,sep=""),stringsAsFactors=F)
chrInfo$chrName[which(chrInfo$chrName=="chr23")]="chrX"
chrInfo$chrName[which(chrInfo$chrName=="chr24")]="chrY"

i2=1:2

annThis$snp50=annThis$snp="No"
timeStamp=Sys.time()
print(format(timeStamp, "%x %X"))
for (chr in 1) {
    #for (chr in 1:nrow(chrInfo)) {
    res=snplocs(genome2, chrInfo$chrName[chr])[i2,]
    i1=which(annThis$CHR==chr)
    for (i in i1) {
        if (any(res$loc>(annThis$MAPINFO[i]-50) & res$loc<(annThis$MAPINFO[i]+50),na.rm=T)) annThis$snp50="Yes"
    }
}
timeStamp=c(timeStamp,Sys.time())
print(format(timeStamp[2], "%x %X"))
print(diff(timeStamp))


annThis$snp50=annThis$snp="No"
timeStamp=Sys.time()
print(format(timeStamp, "%x %X"))
for (chr in 1) {
    res=snplocs(genome2, chrInfo$chrName[chr])[i2,]
    i1=which(annThis$CHR==chr)
    annThis$snp[i1]=sapply(annThis$MAPINFO[i1],function(x) {
        ifelse(any(res$loc>(x-50) & res$loc<(x+50),na.rm=T),"Yes","No")
    },USE.NAMES=F)
}
timeStamp=c(timeStamp,Sys.time())
print(format(timeStamp[2], "%x %X"))
print(diff(timeStamp))




############# Comparison of total methylation levels, globally and regionally
## Not used
## Not working

library(rtracklayer)
library(ChIPseeker)


#dir.create(path)
path="tmp"

cgi_distribution <- function(GR, cgi.all) {
    GR.hyper <- GR[GR$b>0]
    GR.hypo <- GR[GR$b<0]
    print("mergeByOverlaps")
    DMRs.cgi <- mergeByOverlaps(cgi.all,GR)
    print("duplicated")
    DMRs.cgi <- DMRs.cgi[!duplicated(DMRs.cgi$dmr), ]
    DMRs.cgi.hyper <- DMRs.cgi[DMRs.cgi$b>0,]
    DMRs.cgi.hypo <- DMRs.cgi[DMRs.cgi$b<0,]
    print("000")
    df <- data.frame(matrix(0,5,3))
    rownames(df) <- c("N_Shore", "N_Shelf", "CGI", "S_Shelf", "S_Shore")
    for (x in rownames(df)){
        df[x,1] <- sum(DMRs.cgi$feature == x)
        df[x,2] <- sum(DMRs.cgi.hyper$feature==x)
        df[x,3] <- sum(DMRs.cgi.hypo$feature==x)
    }
    print("001")
    Open_sea <- c ( length(GR) - dim(DMRs.cgi)[1],
    length(GR.hyper) - dim(DMRs.cgi.hyper)[1] ,
    length(GR.hypo) - dim(DMRs.cgi.hypo)[1] )
    df <- rbind(df,Open_sea)
    rownames(df)[6] <- "Open_sea"
    colnames(df) <- c("All","hyper","hypo")
    print("002")
    df.perc <- df
    df.perc$All <- df$All/sum(df$All)*100
    df.perc$hyper <- df$hyper/sum(df$hyper)*100
    df.perc$hypo <- df$hypo/sum(df$hypo)*100
    df.perc <- as.matrix(df.perc)
    print("003")
    return(df.perc)
}

feature_distribution <- function(GR, txdb) {
    if (genome=="hg19" | genome=="hg38") {annotation="org.Hs.eg.db"}
    if (genome=="mm9" | genome=="mm10") {annotation="org.Mm.eg.db"}
    
    GR.hyper <- GR[GR$b>0]
    GR.hypo <- GR[GR$b<0]
    
    suppressMessages(peakAnno <- annotatePeak(GR, tssRegion=c(-2000, 2000), TxDb=txdb, annoDb=annotation))
    suppressMessages(peakAnno_hyper <- annotatePeak(GR.hyper, tssRegion=c(-2000, 2000), TxDb=txdb, annoDb=annotation))
    suppressMessages(peakAnno_hypo <- annotatePeak(GR.hypo, tssRegion=c(-2000, 2000), TxDb=txdb, annoDb=annotation))
    df.perc <- cbind(peakAnno@annoStat[,2], peakAnno_hyper@annoStat[,2],peakAnno_hypo@annoStat[,2])
    colnames(df.perc) <- c("All","hyper","hypo")
    rownames(df.perc) <- peakAnno@annoStat[,1]
    return(df.perc)
}


cgi.all <- get_cgi_annotation(genome)
txdb <- getTxdb(genome)
GBM$dmr <- table$name
GBM$b <- table$meanbetafc
save(GBM,file="/tmp/gbm.r")

print("GGI distribution")
df1 <- cgi_distribution(GBM, cgi.all)
print("feature distribution")
df2 <- feature_distribution(GBM, txdb)
print(df1)
print(df2)

colors1 <- c("Cyan","SandyBrown","ForestGreen","SandyBrown","Cyan", "DarkBlue")
colors2 <- c("DarkBlue","Blue","DarkGreen","DarkGreen","Green","LightGreen","Gold","Wheat","Grey","Black")

jpeg(paste(path, "/distribution",gsub("-","_",modelId),".jpg",sep=""), width=800, height=800)
par(mfrow = c(2,1), mar=c(5, 5, 4, 14))
barplot(df1, main="DMR Distribution according to CpG islands (CGI)", xlab="", col=colors1, legend = rownames(df1), horiz=T, args.legend=list(x="topright",bty = "n", inset=c(-0.30, 0)))
barplot(df2, main="DMR Distribution according to feature annotation", xlab="", col=colors2, legend = rownames(df2), horiz=T, args.legend=list(x="topright",bty = "n", inset=c(-0.40, 0)))
dev.off()

####################################################################
####################################################################
## Section 4
## Create table 1, table 2

computerFlag="cluster"
computerFlag=""

if (computerFlag=="cluster") {
} else {
    dirSrc="/Users/royr/UCSF/"
    dirSrc2=dirSrc
    setwd(paste(dirSrc2,"JoeWiemels/leukMeth",sep=""))
}
## ------------------------

output_dir="tmp"
if (computerFlag=="cluster") {
    dirAnn=""
    dirRes="./"
} else {
    dirAnn="results/moba/misc/"
    dirRes="results/moba/output/stat/"
}


load(paste(dirAnn,"tesm.RData",sep=""))


## ------------------------
load("tmp_3.RData")
load("data.RData")
load("tesm.RData")
tmp=rep(NA,nrow(tesm))
ann2=data.frame(IlmnID=tesm$IlmnID,snp=tmp,snp50=tmp,stringsAsFactors=F)

if (F) {
    nProbe=-1
    tdat=read.table(paste("meth.2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=nProbe)
    rownames(tdat)=tdat$cpgId
    tdat=tdat[,-1]
    tdat=as.matrix(tdat)
    tdat=t(tdat)
    
    pdata=read.table(paste("pdata.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    
    counts=read.table(paste("counts.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    rownames(counts)=counts$id
    counts=as.matrix(counts[,-1])
}

## ------------------------


fileList=dir(dirRes,pattern="ind.res_model")
fileList=fileList[1:2]
modelList=c("1-i","1-ii","2a-i","2a-ii","2b-i","2b-ii","3b-i","3b-ii","3c-i","3c-ii","4b-i","4b-ii","4c-i","4c-ii")
modelList=modelList[c(1,14)]
modelList=modelList[c(1,10)]
tmp=rep(NA,length(modelList)); tmpC=rep("",length(modelList))
tbl1 <- data.frame(model=tmpC,numSam=tmp,numProbe=tmp,numSig=tmp,lambda=tmp,stringsAsFactors=F)
tbl1$numSam=nrow(tdat)
tbl1$numProbe=ncol(tdat)
for (fId in 1:length(modelList)) {
    #modelId=sub("ind.res_model","",sub(".txt","",fileList[fId],fixed=T),fixed=T)
    modelId=modelList[fId]
    
    cat("\n\n============================ ",modelId,"\n\n",sep="")
    
    
    all.results1=read.table(paste(dirRes,paste("ind.res_model",gsub("-","_",modelId),".txt",sep=""),sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    cat("Statictic type:")
    print(table(all.results1$statisticType),exclude=NULL)
    pvaluesInter=qvaluesInter=deltabeta=deltabetaM=deltabetaF=rep(NA,nrow(all.results1))
    
    ## Calculate lambda
    lambda <- qchisq(median(all.results1$P_VAL,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
    
    jpeg(paste(output_dir, "/qqplot_model",gsub("-","_",modelId),".jpg",sep=""))    #### change qqplot name according to regression #### model code
    #suppressWarnings(qq(all.results1$P_VAL,main=paste("QQ plot: ","lambda=",lambda,sep="")))
    suppressWarnings(qqnorm(all.results1$P_VAL,main=paste("QQ plot: ","lambda=",lambda,sep="")))
    invisible(dev.off())
    
    qvalues <- p.adjust(all.results1$P_VAL,method="BH")
    pvBonf <- p.adjust(all.results1$P_VAL,method="bonferroni")
    if ("Inter_P_VAL"%in%names(all.results1)) {
        pvaluesInter <- all.results1$Inter_P_VAL
        qvaluesInter <- p.adjust(all.results1$Inter_P_VAL,method="BH")
    }
    group1 <- which(pdata$caco=="1")   ### define first group of
    group2 <- which(pdata$caco=="0")
    # selection <- sort(c(group1,group2))
    deltabeta <- rowMeans(betas[,group1],na.rm = T)-rowMeans(betas[,group2],na.rm = T) #(Case - Control)
    
    if (F) {
        group1M <- which(pdata$caco=="1" & pdata$Gender=="M")   ### define first group of
        group2M <- which(pdata$caco=="0" & pdata$Gender=="M")
        deltabetaM <- rowMeans(betas[,group1M],na.rm = T)-rowMeans(betas[,group2M],na.rm = T) #(Case - Control)
        
        group1F <- which(pdata$caco=="1" & pdata$Gender=="F")   ### define first group of
        group2F <- which(pdata$caco=="0" & pdata$Gender=="F")
        deltabetaF <- rowMeans(betas[,group1F],na.rm = T)-rowMeans(betas[,group2F],na.rm = T) #(Case - Control)
        all.results1f <- cbind(all.results1,qvalues,pvBonf,qvaluesInter,deltabeta,deltabetaM,deltabetaF)
    }
    
    tbl2=cbind(all.results1[,c("probeID","BETA","SE","z.value","P_VAL")],qvalues,pvBonf,pvaluesInter,qvaluesInter,deltabeta,tesm$nearestgene,ann2[,c("snp","snp50")])
    
    tbl1$model[fId]=modelId
    tbl1$numSig[fId]=sum(all.results1$qvalues<0.05,na.rm=T)
    tbl1$lambda[fId]=lambda
    
    write.table(paste(c("probeID","Beta","SE","t value","P_VAL","qvalues","Bonferonni","Inter_P_VAL","QvaluesInter","Deltabeta","nearestgene","SNP","Near SNP"),collapse=","),file=paste("table2_model",modelId,"_cbs.csv",sep=""), sep=",", col.names=F, row.names=F, quote=F, append=F)
    write.table(tbl2,file=paste("table2_model",modelId,"_cbs.csv",sep=""), sep=",", col.names=F, row.names=F, quote=F, append=T)
}
write.table(paste(c("Model name","Number of samples","Number of probes","Number of significant CpGs (FDR<0.05)","lambda"),collapse=","),file=paste("table1_cbs.csv",sep=""), sep=",", col.names=F, row.names=F, quote=F, append=F)
write.table(tbl1,file=paste("table1_cbs.csv",sep=""), sep=",", col.names=F, row.names=F, quote=F, append=T)


####################################################################
####################################################################
