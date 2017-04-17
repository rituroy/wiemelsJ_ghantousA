## Run script.R section 1 then refactor_forMoba.R first

####################################################################

computerFlag=""
computerFlag="cluster"

if (computerFlag=="cluster") {
} else {
    dirSrc="/Users/royr/UCSF/"
    dirSrc2=dirSrc
    setwd(paste(dirSrc2,"JoeWiemels/leukMeth",sep=""))
}

####################################################################
####################################################################

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
        #targets$guthrieId=paste("X",targets$guthrieId,sep="")
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
    targets=targets[which(!is.na(targets$guthrieId) & !targets$guthrieId%in%c("1458G","1762G","0635G","1588G")),]
}
targets$Basename=paste("/",targets$Beadchip,"/",targets$Beadchip,"_",targets$Position,sep="")

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


load("tmp_3.RData")

####################################################################
####################################################################
## ---------------------------------

computerFlag=""
computerFlag="cluster"

## ---------------------------------

nProbe=101
nProbe=10001
nProbe=-1

## ---------------------------------
subsetName2=""

normFlag=""
normFlag="_funNorm"
normFlag="_bmiq"

datType="_leuk"; subsetName2=""
datType="_aml"; subsetName2=""
datType="_allGuthSet1"; subsetName2=""
datType="_allGuthSet2"; subsetName2=""

setFlag=ifelse(subsetName2=="",tolower(sub("allGuth","",datType)),subsetName2)

subsetFlag="case"
subsetFlag="ctrl"

subsetFlag=""

#subsetFlag="noHisp"
subsetFlag="hisp"
subsetFlag="noHispWt"

covFlag=""
mediationFlag=F
varFlag=""
covPCFlag="_covPrinComp1234"
winsorTail=0.05

library(limma)
library(FactoMineR)

#for (subsetName2 in c("_set1","_set2")) {
#for (datType in c("_allGuthSet1","_allGuthSet2")) {
		
## ---------------------------------

## All smoke variables from set1 & set2
## All smoke variables in set2 are in set1

## ---------------------------------

## ---------------------------------

#for (subsetFlag in c("noHisp","hisp")) {

#	for (varFlag in c("_hyperdipCtrl","_telamlCtrl","_noHypTelamlCtrl")) {
#for (varFlag in c("_dfeFoodCat","_dfeSupCat","_dfeNat","_dfeTot","_dfeFood","_dfeFort","_dfeNatCat","_dfeTotCat")) {
#for (varFlag in paste("_smoke_",varSmoke$varOut,sep="")) {

##############################################

if (computerFlag=="cluster") {
	setwd("/home/royr/project/JoeWiemels")
} else {
	dirSrc="/Users/royr/UCSF/"
	dirSrc2=dirSrc
	setwd(paste(dirSrc2,"JoeWiemels/leukMeth",sep=""))
}

subsetName=subsetFlag
if (subsetFlag!="") {
	#varFlag=paste(varFlag,"_",subsetFlag,"Subset",sep="")
	subsetName=paste("_",subsetFlag,"Subset",sep="")
	switch(subsetFlag,
		   "case"={varName=sub(")"," cases)",varName)},
		   "ctrl"={varName=sub(")"," controls)",varName)},
		   "hisp"={varName=sub(")"," hispanics)",varName)},
		   "noHisp"={varName=sub(")"," non-hispanics)",varName)},
		   "hyperdip"={varName=sub(")"," hyperdiploids)",varName)},
		   "telaml"={varName=sub(")"," tel/aml1s)",varName)},
		   "noHypTelaml"={varName=sub(")"," non-hyperdiploids/non-tel/aml1s)",varName)}
		   )
}

heading=paste(c(varFlag,", ",subsetFlag,", ",covFlag,", ",covPCFlag,", ",datType,subsetName2,", ",normFlag),collapse="")
cat("\n\n============================ Linear Regression, Epistructure ===========================\n\n")
cat("\n\n============================",varFlag,", ",subsetFlag,", ",covFlag,", ",covPCFlag,", ",datType,subsetName2,", ",normFlag,"===========================\n\n")

##############################################

timeStamp=Sys.time()
print(format(timeStamp, "%x %X"))

if (computerFlag=="cluster") {
	dirMeth="data/"
	dirClin=dirMeth
	dirBW="data/"
	dirCom=dirMeth
	dirRefactor=dirEpistructure=dirMeth
	switch(datType,
		"_allGuthSet2"={
		   dirMeth=dirClin="data/set2/"
		   fNameMeth=paste("beta",normFlag,ifelse(normFlag=="_funNorm","_set2",datType),sep="")
		   fNameClin="clin_guthrieSet2_20140619"
           fNameClin="clin_allGuthSet2_20160928"
		},
		"_allGuthSet1"={
		   dirMeth=dirClin="data/set1/"
		   fNameMeth=paste("beta",normFlag,ifelse(normFlag=="_funNorm","_set1",datType),sep="")
		   fNameClin="final"
           fNameClin="clin_allGuthSet1_20160928"
		   dirMethLeuk=dirMeth
		   fNameMethLeuk="i.LEU.v2"
        },
        "_allGuthSet1Set2"={
            dirMeth=dirClin="data/set1set2/"
            fNameMeth=paste("beta",normFlag,ifelse(normFlag=="_funNorm","_set1set2",datType),sep="")
            fNameClin="clin_guthrieSet1Set2_20140619"
            fNameClin="clin_guthrieSet1Set2_20151022"
        },
		"_leuk"={
		   dirMeth=dirClin=dirClin2="data/set1/"
		   fNameMeth=paste("beta",normFlag,datFlag,sep="")
		   fNameClin="i.LEU.v2"
		   fNameClin2="0708011 Sample_Sheet (Fetal blood)"
		},
		"_allGuthSet1Set2Combat"={
		   dirMeth=dirClin="data/set1set2/"
		   fNameMeth="combatAdjustedBeta_allGuthSet1Set2_set"
		   fNameClin="clin_guthrieSet1Set2_20140619"
		},
		"_aml"={
		   dirClin=dirMeth=dirRefactor="data/aml/"
		   fNameMeth=paste("beta",normFlag,"_aml",sep="")
		   fNameClin=paste("clin_aml_20150114",sep="")
		},
		"_ivorra"={
		   dirMeth=dirClin="data/ivorra2014/"
		   fNameMeth="beta_ivorra"
		   fNameClin="clin_ivorra"
		}
	)
} else {
	dirBW="docs/birthWeight/"
	dirCom="docs/all/"
    #dirRefactor="docs/SemiraGonsethNussle/refactor/"
    dirEpistructure="docs/SemiraGonsethNussle/epistructure/"
    switch(datType,
		"_allGuthSet2"={
		   dirMeth=dirClin="docs/all/set2/"
           dirRefactor=dirClin
           fNameMeth=paste("beta",normFlag,ifelse(normFlag=="_funNorm","_set2",datType),sep="")
		   fNameClin="clin_guthrieSet2_20140619"
           fNameClin="clin_allGuthSet2_20160928"
	   },
	   "_allGuthSet1"={
			dirMeth=dirClin="docs/all/set1/"
            dirRefactor=dirClin
            fNameMeth=paste("beta",normFlag,ifelse(normFlag=="_funNorm","_set1",datType),sep="")
			fNameClin="final"
            fNameClin="clin_allGuthSet1_20160928"
			dirMethLeuk="docs/all/set1/LEU.data/"
			fNameMethLeuk="i.LEU.v2"
		},
		"_leuk"={
		   dirMeth="docs/all/set1/"
           dirRefactor=dirClin
           fNameMeth=paste("beta",normFlag,datFlag,sep="")
		   dirClin="docs/all/set1/LEU.data/"
		   fNameClin="i.LEU.v2"
		   dirClin2="docs/all/set1/preBcell/"
		   fNameClin2="0708011 Sample_Sheet (Fetal blood)"
		},
		"_aml"={
		   dirClin=dirMeth="docs/aml/"
		   dirRefactor="docs/aml/refactor/"
		   fNameMeth=paste("beta",normFlag,"_aml",sep="")
		   fNameClin=paste("clin_aml_20150114",sep="")
	   },
	   "_ivorra"={
		   dirMeth=dirClin="docs/misc/ivorra2014/"
           dirRefactor=dirClin
           fNameMeth="beta_ivorra"
		   fNameClin="clin_ivorra"
		}
	)
}

## ----------------------------------------------
fName1="_moba"

load(paste(dirRefactor,"Refactor_dat",ifelse(subsetFlag=="","","_"),subsetFlag,datType,fName1,".RData",sep=""))
colnames(Refactor_dat)=paste("prinComp",1:6,sep="")
j=match(colnames(betas),rownames(Refactor_dat)); j1=which(!is.na(j)); j2=j[j1]
betas=betas[,j1]
pdata=cbind(pdata[j1,],Refactor_dat[j2,])

if (T) {
if (computerFlag=="") {
	load(file="ann.RData")
    crp_probes <- "docs/Crossreactive_probes.csv"       ###add the path to the csv file
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
    crp_probes <- "data/Crossreactive_probes.csv"       ###add the path to the csv file
}


ann$geneSym=sapply(toupper(ann$UCSC_RefGene_Name),function(x) {
    strsplit(x,";")[[1]][1]
},USE.NAMES=F)
ann$geneSym[is.na(ann$geneSym)]=""


########################remove unwanted probes (cross reactive probes)
cross_reactive <- read.table(crp_probes, header=T, sep="\t")$TargetID
#RGset <- RGset[ ! featureNames(RGset)%in%cross_reactive,]


#keep=!ann$CHR%in%c(23,24)
keep=!ann$CHR%in%c(23,24) & apply(betas,1,function(x) {any(!is.na(x))})
keep=ann$snp==0 & !ann$CHR%in%c(23,24) & apply(betas,1,function(x) {any(!is.na(x))})
keep=ann$snp==0 & ann$CHR%in%1:22 & apply(betas,1,function(x) {any(!is.na(x))})
keep=ann$snp==0 & ann$CHR%in%1:22 & !ann$IlmnID%in%cross_reactive & apply(betas,1,function(x) {any(!is.na(x))})
}

################################################
# Regression
################################################

save.image("tmp.RData")
load("tmp.RData")
rownames(pdata)=pdata$id
ann=data.frame(IlmnID=rownames(betas),id=rownames(betas),stringsAsFactors=T)

prId=read.table(paste(dirEpistructure,"epistructure_reference_sites.txt",sep=""), sep="\t", h=F, quote="", comment.char="",as.is=T,fill=T)
prId=prId[,1]

keep=which(ann$IlmnID%in%prId)

model1=""
model1=paste(model1,"+",paste("prinComp",strsplit(sub("_covPrinComp","",covPCFlag),"")[[1]],collapse="+",sep=""),sep="")
model1=sub("+","~",model1,fixed=T)
model1=as.formula(model1)
Xunadj=model.matrix(model1,data=pdata)
samId=match(rownames(Xunadj),pdata$id)

lm.unadj=eBayes(lmFit(betas[keep,samId], Xunadj),robust=ifelse(is.null(winsorTail),F,T), winsor.tail.p=winsorTail)
res=residuals(lm.unadj,y=betas[keep,samId])

head(lm.unadj$coef)
head(lm.unadj$p.value)


table(is.na(res))
res_nona=na.omit(res)
res_PCA=PCA(res_nona, graph = F)

epistr1=res_PCA$var$coord[, c(1)]
epistr2=res_PCA$var$coord[, c(2)]
epistr3=res_PCA$var$coord[, c(3)]
epistr4=res_PCA$var$coord[, c(4)]
epistr5=res_PCA$var$coord[, c(5)]
epistr=cbind(epistr1, epistr2, epistr3, epistr4, epistr5)

tbl=cbind(id=rownames(epistr),as.data.frame(epistr))
write.table(tbl,file=paste("epistructure",ifelse(subsetFlag=="","","_"),subsetFlag,datType,fName1,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
