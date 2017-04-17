

###################################################  Models

argv=commandArgs(TRUE)
modelId=argv[1]

nProbe=101
nProbe=-1

if (F) {
    modelId="2b-ii"
    nProbe=101
}

cat("modelId: ",modelId,"\n\n",sep="")

switch(modelId,
"1-i"={
    compName="Methylation = Case/control"
    testFlag="RLMtestmodif"
    colId1="caco"
},
"1-ii"={
    compName="Methylation = Case/control + WBC"
    testFlag="RLMtestmodif"
    colId1=c("caco","nRBC","CD8T","CD4T", "NK","Bcell","Mono")
},
"2a-i"={
    compName="Methylation = Case/control + Gender"
    testFlag="RLMtestmodif"
    colId1=c("caco","Gender",paste("epistr",1:5,sep=""))
},
"2a-ii"={
    compName="Methylation = Case/control + Gender + WBC"
    testFlag="RLMtestmodif"
    colId1=c("caco","Gender","nRBC","CD8T","CD4T", "NK","Bcell","Mono",paste("epistr",1:5,sep=""))
},
"2b-i"={
    compName="Methylation = Case/control + Gender + Case/control*Gender"
    testFlag="RLMtestinteraction"
    colId1=c("caco","Gender",paste("epistr",1:5,sep=""))
},
"2b-ii"={
    compName="Methylation = Case/control + Gender + Case/control*Gender + WBC"
    testFlag="RLMtestinteraction"
    colId1=c("caco","Gender","nRBC","CD8T","CD4T", "NK","Bcell","Mono",paste("epistr",1:5,sep=""))
},
"3b-i"={
    compName="Methylation = Case/control + Gender + YofB"
    testFlag="RLMtestmodif"
    colId1=c("caco","Gender","yob",paste("epistr",1:5,sep=""))
},
"3b-ii"={
    compName="Methylation = Case/control + Gender + YofB + WBC"
    testFlag="RLMtestmodif"
    colId1=c("caco","Gender","yob","nRBC","CD8T","CD4T", "NK","Bcell","Mono",paste("epistr",1:5,sep=""))
},
"3c-i"={
    compName="Methylation = Case/control + Gender + Case/control*Gender + YofB"
    testFlag="RLMtestinteraction"
    colId1=c("caco","Gender","yob",paste("epistr",1:5,sep=""))
},
"3c-ii"={
    compName="Methylation = Case/control + Gender + Case/control*Gender + YofB + WBC"
    testFlag="RLMtestinteraction"
    colId1=c("caco","Gender","yob","nRBC","CD8T","CD4T", "NK","Bcell","Mono",paste("epistr",1:5,sep=""))
},
"4b-i"={
    compName="Methylation = Case/control + Gender + AgeofD"
    testFlag="RLMtestmodif"
    colId1=c("caco","Gender","ch_ageref",paste("epistr",1:5,sep=""))
},
"4b-ii"={
    compName="Methylation = Case/control + Gender + AgeofD + WBC"
    testFlag="RLMtestmodif"
    colId1=c("caco","Gender","ch_ageref","nRBC","CD8T","CD4T", "NK","Bcell","Mono",paste("epistr",1:5,sep=""))
},
"4c-i"={
    compName="Methylation = Case/control + Gender + Case/control*Gender + AgeofD"
    testFlag="RLMtestinteraction"
    colId1=c("caco","Gender","ch_ageref",paste("epistr",1:5,sep=""))
},
"4c-ii"={
    compName="Methylation = Case/control + Gender + Case/control*Gender + AgeofD + WBC"
    testFlag="RLMtestinteraction"
    colId1=c("caco","Gender","ch_ageref","nRBC","CD8T","CD4T", "NK","Bcell","Mono",paste("epistr",1:5,sep=""))
}
)
#colId1=colId1[-grep("epistr",colId1)]

#c("caco","Gender","yob","ch_ageref","nRBC","CD8T","CD4T", "NK","Bcell","Mono","Gran",paste("epistr",1:5,sep=""))


library(MASS) # rlm function for robust linear regression
library(sandwich) #Huberís estimation of the standard error
library(lmtest) # to use coeftest


#load("data.RData")
tdat=read.table(paste("meth.2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=nProbe)
rownames(tdat)=tdat$cpgId
tdat=tdat[,-1]
tdat=as.matrix(tdat)
tdat=t(tdat)

pdata=read.table(paste("pdata.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

counts=read.table(paste("counts.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
rownames(counts)=counts$id
counts=as.matrix(counts[,-1])

probeVec=1:ncol(tdat)

cat("\n\n============================ ",modelId,": ",testFlag,"\n",compName,"\n\n",sep="")

#xMat <- cbind(pdata,counts)[,c("caco","Gender","yob","ch_ageref","nRBC","CD8T","CD4T", "NK","Bcell","Mono","Gran",paste("epistr",1:5,sep=""),"Plate","Position")]    ##### use same variable indexes in PCs as before SVA
if (length(colId1)==1) {
    xMat <- cbind(pdata,counts)[,c(colId1,"Position")]
} else {
    xMat <- cbind(pdata,counts)[,colId1]
}

#### for RLMtestmodif use this
if (testFlag=="RLMtestmodif") {
    modelThis=paste("tdat[,probeVec[p]]~",paste(colId1,collapse="+"),sep="")
    cat(modelThis,"\n\n")
    modelThis=as.formula(modelThis)
    ind.res=matrix(nrow=length(probeVec),ncol=4,dimnames=list(colnames(tdat)[probeVec],c("Estimate", "Std. Error", "z value", "Pr(>|z|)")))
    res2=data.frame(isT=rep("",length(probeVec)),stringsAsFactors=F)
    
    timeStamp=Sys.time()
    print(format(timeStamp, "%x %X"))
    for (p in 1:length(probeVec)) {
        mod <- tryCatch(lm(modelThis, data=xMat),error = function(e) e)
        if (class(mod) == "try-error") {
            #print(paste("error thrown by column", methcol))
            #invisible(rep(NA, 3))
        } else {
            cf <- coeftest(mod, vcov=vcovHC(mod, type="HC0"))
            if ("t value"%in%colnames(cf)) {
                ind.res[p,] <- cf[2, c("Estimate", "Std. Error", "t value", "Pr(>|t|)")]
                res2$isT="t"
            } else {
                ind.res[p,] <- cf[2, c("Estimate", "Std. Error", "z value", "Pr(>|z|)")]
                res2$isT="z"
            }
        }
        
    }
    timeStamp=c(timeStamp,Sys.time())
    print(format(timeStamp[2], "%x %X"))
    print(diff(timeStamp))
    colId2 <- c("probeID","BETA","SE","z value", "P_VAL","statisticType")
} else {
    modelThis=paste("tdat[,probeVec[p]]~",paste(colId1[1:2],collapse="*"),sep="")
    if (ncol(xMat)>2) {
        modelThis=paste(modelThis,"+",paste(colId1[3:length(colId1)],collapse="+"),sep="")
    }
    cat(modelThis,"\n\n")
    modelThis=as.formula(modelThis)
    ind.res=matrix(nrow=length(probeVec),ncol=5,dimnames=list(colnames(tdat)[probeVec],c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "Pr(>|z|)_int")))
    res2=data.frame(isT=rep("",length(probeVec)),stringsAsFactors=F)
    
    timeStamp=Sys.time()
    print(format(timeStamp, "%x %X"))
    for (p in 1:length(probeVec)) {
        mod <- tryCatch(lm(modelThis, data=xMat),error = function(e) e)
        if (class(mod) == "try-error") {
            #print(paste("error thrown by column", methcol))
            #invisible(rep(NA, 3))
        } else {
            cf <- coeftest(mod, vcov=vcovHC(mod, type="HC0"))
            if ("t value"%in%colnames(cf)) {
                ind.res[p,] <- c(cf[2, c("Estimate", "Std. Error", "t value", "Pr(>|t|)")],cf[nrow(cf), "Pr(>|t|)"])
                res2$isT="t"
            } else {
                ind.res[p,] <- c(cf[2, c("Estimate", "Std. Error", "z value", "Pr(>|z|)")],cf[nrow(cf), "Pr(>|z|)"])
                res2$isT="z"
            }
        }
        
    }
    timeStamp=c(timeStamp,Sys.time())
    print(format(timeStamp[2], "%x %X"))
    print(diff(timeStamp))
    colId2 <- c("probeID","BETA","SE","z value", "P_VAL","Inter_P_VAL","statisticType")
}

ind.res=as.data.frame(ind.res,stringsAsFactors=F)
ind.res=cbind(probeID=rownames(ind.res),ind.res,isT=res2$isT,stringsAsFactors=F)
names(ind.res)=colId2
names(ind.res)=colId2

#save(ind.res,file=paste("ind.res_model",gsub("-","_",modelId),".RData",sep=""))
write.table(ind.res,file=paste("ind.res_model",gsub("-","_",modelId),".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
