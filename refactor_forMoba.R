## Run script.R section 1 first

####################################################################

cohort="_allGuthSet2"
cohort="_noHispWt_allGuthSet2"


baseDir="data/set2/idat/"
#pdata=read.table("data/set2/clinGuthrieReplJune2012.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
pdata=read.table("data/set2/clin_guthrieSet2_20140619.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
tbl=read.table("data/meth_yrbirth.csv",sep=",",h=T,quote="",comment.char="",as.is=T,fill=T)
pdata$yob=NA
j=match(pdata$subjectID,tbl$subjectid); j1=which(!is.na(j)); j2=j[j1]
pdata$yob[j1]=tbl$ch_birth_year[j2]


if (length(grep("_allGuthSet2",cohort))==1) {
    # "0855G" - T-ALL
    pdata=pdata[which(!is.na(pdata$guthrieId) & !pdata$guthrieId%in%c("1458G","1762G","0635G","1588G","0855G")),]
}
if (length(grep("_noHispWt",cohort))==1) {
    pdata=pdata[which(pdata$int_ch_ethnicity==2),]
}
pdata$Basename=paste("/",pdata$Beadchip,"/",pdata$Beadchip,"_",pdata$Position,sep="")

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

# Functions
#------------------
# return the variable position in dataset

number <- function(data) {
    data_matrix <- matrix(c(1:dim(data)[2],(names(data))),,2)
    data_matrix
}

# Remove NAs if more than 3% of the row
remove_NAs <- function(x){
    if (table(is.na(x))/ length(x) < 0.97) {TRUE} else {FALSE}
}

# Impute mean values on missing data

impute_mean <- function(x) {
    ifelse (is.na(x), mean(x, na.rm = TRUE), x)
}

## ---------------------------------

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
cat("\n\n============================ Linear Regression, Refactor ===========================\n\n")
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
    cross_reactive <- read.table(crp_probes, header=T, sep="\t")$TargetID
    save(ann,cross_reactive,file="ann.RData")
}

################################################

fName1="_moba"
load("ann.RData")

betaThis=betas
clinThis=pdata
rm(betas,pdata)

# Load outliers found by Elior
# ------------------------


if (T) {
    ## Find sample outliers
    
    #load("tmp.RData")
    i=match(rownames(betaThis),ann$IlmnID); i1=which(!is.na(i)); i2=i[i1]
    betaThis=betaThis[i1,]
    ann=ann[i2,]
    
    IDs_to_keep=1:ncol(betaThis)
	x=apply(betaThis,1,function(x) mean(!is.na(x)))
	CpG_to_keep <- which(ann$snp==0 & ann$CHR%in%1:22 & !ann$IlmnID%in%cross_reactive & x>=0.03)

	betaThis <- betaThis[CpG_to_keep,IDs_to_keep]
	betaThis=t(betaThis)
	R_est <- PCA(betaThis, graph = F, ncp =6)
	save(R_est,file=paste("R_est_init",ifelse(subsetFlag=="","","_"),subsetFlag,datType,fName1,".RData",sep=""))
    
    limThis=c(-200,200)
    limThis=c(-250,250)
    limThis=c(-300,300)
	Refactor_dat <- R_est$ind$coord[,c(1:6)]
	png("pca_1.png")
	plot(R_est,xlim=c(-2000,1000),ylim=c(-2000,1000))
	abline(c(0,1),lty="dotted")
	dev.off()
    ## Outliers are those outside the dotted lines
	png("pca_1_1.png")
	plot(R_est$ind$coord[,1],R_est$ind$coord[,2],xlim=c(-2000,1000),ylim=c(-2000,1000))
	abline(c(0,1),lty="dotted"); abline(h=0,lty="dotted"); abline(v=0,lty="dotted")
	abline(h=limThis,lty="dotted")
	abline(v=limThis,lty="dotted")
	dev.off()
	png("pca_1_2.png")
	plot(R_est,xlim=c(-2000,1000),ylim=c(-2000,1000),label="none")
	abline(c(0,1),lty="dotted")
	dev.off()

    #load(file=paste("R_est_init",ifelse(subsetFlag=="","","_"),subsetFlag,datType,fName1,".RData",sep=""))
    table((abs(R_est$ind$coord[,1])>limThis[2] | abs(R_est$ind$coord[,2])>limThis[2]))
    ## Outliers
	paste(rownames(R_est$ind$coord)[which(abs(R_est$ind$coord[,1])>limThis[2] | abs(R_est$ind$coord[,2])>limThis[2])],collapse=",")
}



#load("tmp.RData")
i=match(rownames(betaThis),ann$IlmnID); i1=which(!is.na(i)); i2=i[i1]
betaThis=betaThis[i1,]
ann=ann[i2,]

## Exclude outliers
if (F) {
    if (datType=="_allGuthSet1") {
        IDs_to_keep=which(!colnames(betaThis)%in%c("X1218G","X0580G","X0404G","X1336G","X1225G","X1202G","X1162G","X1224G","X0527G","X0665G","X0451G","X0557G","X1664G","X0498G","X0419G","X1446G","X2006G","X1577G","X0486G","X1419G","X1739G","X1540G","X1736G","X1988G","X1718G","X0525G","X2221G","X1749G","X0420G","X1769G","X0511G","X0563G","X1650G","X1237G","X1328G","X1297G","X0993G","X1114G","X1042G","X1116G","X1870G","X1240G","X1967G","X1011G","X1439G","X1906G","X1070G","X1933G","X0935G","X1244G","X1639G","X1384G","X2045G","X1246G","X1579G","X1205G","X0489G","X0280G","X0183G","X0924G","X0470G","X0162G","X0481G","X1288G","X1578G","X1037G","X1157G","X0384G","X0283G","X0716G","X1323G","X1965G","X0804G","X1392G","X1221G","X0284G","X0194G","X0264G","X0285G","X0459G","X1141G","X1098G","X0594G","X0586G","X1280G","X0401G","X0136G","X0504G","X0332G","X0381G","X0584G","X0541G"))
    } else if (datType=="_allGuthSet1Set2_ctrlSubset") {
        IDs_to_keep=which(!colnames(betaThis)%in%c("X1218G","X0580G","X1336G","X1204G","X1225G","X1202G","X1252G","X1162G","X1224G","X0527G","X0665G","X0465G","X0471G","X1477G","X1562G","X0486G","X1718G","X0525G","X1769G","X1650G","X0311G","X0469G","X1248G","X1639G","X1384G","X1579G","X0885G","X1341G","X1008G","X1205G","X1117G","X1405G","X1288G","X1578G","X1260G","X0874G","X1700G","X1326G","X1164G","X1359G","X1392G","X1147G","X1221G","X1141G","X1098G","X0594G","X0586G","X1280G","X1404G","X0552G","X2046G","X1770G","X1824G","X2029G","X0247G","X1658G","X0080G","X1360G","X1024G","X1173G","X0229G","X1021G","X1994G","X0132G","X0091G","X0137G","X0197G","X0531G","X1881G","X1781G","X0784G","X0184G","X1285G","X0053G","X1230G","X1920G","X0380G","X0231G","X0876G","X0125G","X0815G","X0036G","X0643G","X1628G","X1092G","X1901G","X0827G","X1169G","X1807G","X1884G","X1777G","X0513G","X0270G","X0437G","X1834G","X1946G","X1277G","X1843G","X1859G","X0173G","X0107G","X1123G","X1996G","X2188G","X0065G","X0225G","X0149G","X1184G"))
    } else {
        IDs_to_keep=which(!colnames(betaThis)%in%c("X1288G","X1466G","X0153G","X0077G","X0201G","X1191G","X1264G","X1219G","X0873G","X0927G","X1075G","X1766G","X0691G"))
    }
}
load(file=paste("R_est_init",ifelse(subsetFlag=="","","_"),subsetFlag,datType,fName1,".RData",sep=""))
limThis=c(-300,300)
IDs_to_keep=which(abs(R_est$ind$coord[,1])<=limThis[2] & abs(R_est$ind$coord[,2])<=limThis[2])

x=apply(betaThis,1,function(x) mean(!is.na(x)))
CpG_to_keep <- which(ann$snp==0 & ann$CHR%in%1:22 & !ann$IlmnID%in%cross_reactive & x>=0.03)
betaThis <- betaThis[CpG_to_keep,IDs_to_keep]
ann=ann[CpG_to_keep,]
clinThis=clinThis[IDs_to_keep,]


# # Removing CpG sites with >3% of NAs
# # --------------------------------
# # set 1
# remove_NAs <- function(x){
	# if (table(is.na(x))/ length(x) < 0.97) {TRUE} else {FALSE}
	# }

# res_set1 <- apply(set1, 1, remove_NAs)

# set1_nona <- set1[(res_set1) == F,]; dim(set1_nona)

# # set 2

# res_set2 <- apply(set2, 1, remove_NAs)

# set2_nona <- set2[(res_set2) == F,]; dim(set2_nona)

# Imputing mean methylation levels on the other NAs
# --------------------------------------------

# set 2
betaThis <- t(apply(betaThis, 1, impute_mean))


# reFACTor calculations
# -------------------------
# Reduction of the number of CpG sites (acc. to Elior Rahmani)

# x1 <- set1
x2 <- betaThis

# mean_x1 <- apply(x1, 1, mean)
# names(mean_x1) <- rownames(x1)
mean_x2 <- apply(x2, 1, mean)
names(mean_x2) <-rownames(x2)

x2_proc <- mean_x2[mean_x2 < 0.9 & mean_x2 > 0.1]

betaThis <- x2[which(rownames(x2) %in% names(x2_proc)),]

#save(betaThis, file = paste("beta_matrix",ifelse(subsetFlag=="","","_"),subsetFlag,datType,fName1,".RData",sep=""))
#load(file = paste("beta_matrix",ifelse(subsetFlag=="","","_"),subsetFlag,datType,fName1,".RData",sep=""))

# refactor from Elior''s matlab code

# Running a standard PCA...
# ---------------------------------

O_prime <- (betaThis)

# z score
zscore <- function(x){
	(x - mean(x)) / sd(x)
	}

z_norm_O_prime <- t(apply(O_prime, 1, zscore))

K = 6

t0 <- Sys.time()
res_PCA <- princomp(z_norm_O_prime)
Sys.time() - t0


save(res_PCA, file = paste("res_PCA_for_reFACTOR",ifelse(subsetFlag=="","","_"),subsetFlag,datType,fName1,".RData",sep=""))

# Compute a low rank approximation of input data and rank sites...
# ---------------------------------

x <- res_PCA$scores[,c(1:K)] %*% t((res_PCA$loadings[,c(1:K)]))

# An = bsxfun(@minus,O',mean(O',1));
An <- z_norm_O_prime

# Bn=bsxfun(@minus,x,mean(x,1));
Bn <- t(apply(x, 1, zscore))

# An=bsxfun(@times,An,1./sqrt(sum(An.^2,1)));

low_rank <- function(data){
	data2 <- (data^2)
	data3 <- sum(data2)
	data4 <- sqrt(data3)
	data5 <- 1/data4
	data*data5
}
An2 <- t(apply(An, 1, low_rank))

# Bn=bsxfun(@times,Bn,1./sqrt(sum(Bn.^2,1)));

Bn2 <- t(apply(Bn, 1, low_rank))


# Find the distance of each site from its low rank approximation
# ---------------------------------

res_AnBn <- An2-Bn2
# save(res_AnBn, file = paste("For_cluser_res_AnBn.RData",fName1,sep=""))

#load(paste("For_cluser_res_AnBn",fName1,".RData",sep=""))

# distances = sum((An-Bn).^2,1).^0.5 ;
# [~,ranked_list] = sort(distances);

res_interm <- apply((res_AnBn)^2, 1, FUN = sum)
distances <- sqrt(res_interm)
ranked_list <- sort(distances)

# Compute ReFACTor components...
# ---------------------------------
t = 500
# sites = ranked_list(1:t);
sites <- ranked_list[1:t]

# [~,R_est] = pca(zscore(O(sites,:)''));

t_O_sites <- t(O_prime[which(rownames(z_norm_O_prime) %in% names(sites)),])

t_O_sites_zs <- (apply(t_O_sites, 2, zscore))

R_est <- PCA(t_O_sites_zs, graph = F, ncp =6)

Refactor_dat <- R_est$ind$coord[,c(1:6)]

save(Refactor_dat, file = paste("Refactor_dat",ifelse(subsetFlag=="","","_"),subsetFlag,datType,fName1,".RData",sep=""))

#load(file = paste("Refactor_dat",ifelse(subsetFlag=="","","_"),subsetFlag,datType,fName1,".RData",sep=""))

