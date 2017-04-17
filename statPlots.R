
####################################################################
####################################################################
## Volcano plot, etc.

computerFlag="cluster"
computerFlag=""

if (computerFlag=="cluster") {
} else {
    dirSrc="/Users/royr/UCSF/"
    dirSrc2=dirSrc
    setwd(paste(dirSrc2,"JoeWiemels/leukMeth",sep=""))
}
## ------------------------

output_dir="tmp2"
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
#modelList=modelList[c(1,14)]
#modelList=modelList[c(14)]



####################################################################
####################################################################
####################################################################
####################################################################
## ----------------------------------------------

plotFlag=""
plotFlag="_onePlot"

outlierFlag=T
outlierFlag=F

colIdEst="BETA"; colIdPV=c("P_VAL","qvalues"); 	pThres=0.05
colIdEst="BETA"; colIdPV=c("P_VAL","P_VAL"); 	pThres=0.05
colList=c("BETA","P_VAL","qvalues")
colName=c("coef","pv","fdrBH")


colIdEst="BETA"; colIdPV=c("pvaluesInter","Inter_P_VAL"); 	pThres=0.05
colIdEst="BETA"; colIdPV=c("pvaluesInter","qvaluesInter"); 	pThres=0.05
colList=c("BETA","P_VAL","qvalues","pvaluesInter","qvaluesInter")
colName=c("coef","pv","fdrBH","pvInter","fdrBHInter")

ann=tesm
ann$keep=1

## -------------------
for (compId in modelList) {
    
    modelId=compId
    
    fName1=paste("_model",gsub("-","_",modelId),sep="")
    
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
    
    
    

    all.results1=read.table(paste(dirRes,paste("ind.res_model",gsub("-","_",modelId),".txt",sep=""),sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    if ("Inter_P_VAL"%in%colList & !"Inter_P_VAL"%in%names(all.results1)) next
    #cat("Statictic type:")
    #print(table(all.results1$statisticType),exclude=NULL)
    pvaluesInter=qvaluesInter=deltabeta=deltabetaM=deltabetaF=rep(NA,nrow(all.results1))

    cat("\n\n============================ ",modelId,"\n",compName,"\n\n",sep="")

    if ("Inter_P_VAL"%in%names(all.results1)) {
        
        ## Calculate lambda
        lambda <- qchisq(median(all.results1$P_VAL,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
        jpeg(paste(output_dir, "/qqplot_model",gsub("-","_",modelId),".jpg",sep=""))    #### change qqplot name according to regression #### model code
        #suppressWarnings(qq(all.results1$P_VAL,main=paste("QQ plot: ","lambda=",lambda,sep="")))
        suppressWarnings(qqnorm(all.results1$P_VAL,main=paste("QQ plot: ","lambda=",lambda,sep="")))
        invisible(dev.off())
    }
    
    qvalues <- p.adjust(all.results1$P_VAL,method="BH")
    pvBonf <- p.adjust(all.results1$P_VAL,method="bonferroni")
    if ("Inter_P_VAL"%in%names(all.results1)) {
        pvaluesInter <- all.results1$Inter_P_VAL
        qvaluesInter <- p.adjust(all.results1$Inter_P_VAL,method="BH")
    }
    group1 <- which(pdata$caco=="1")   ### define first group of
    group2 <- which(pdata$caco=="0")
    # selection <- sort(c(group1,group2))
    #deltabeta <- rowMeans(betas[,group1],na.rm = T)-rowMeans(betas[,group2],na.rm = T) #(Case - Control)

    stat2=cbind(all.results1[,c("probeID","BETA","SE","z.value","P_VAL")],qvalues,pvBonf,pvaluesInter,qvaluesInter,deltabeta,tesm$nearestgene)
    names(stat2)[match(c("probeID"),names(stat2))]="cpgId"

    i=match(stat2$cpgId,ann$IlmnID)
    stat2=stat2[which(ann$keep[i]==1),]
    iA2=match(stat2$cpgId,ann$IlmnID)
    
    ann2=ann[iA2,]
    
    ####################################################################
    
    i=which(ann$keep[iA2]==1)
    stat=stat2
    fName=fName1
    
    x1=matrix(0,nrow=2,ncol=2,dimnames=list(c("dn","up"),paste(colIdPV[2],c(">=","<"),pThres,sep="")))
    ii=i[which(stat[i,colIdEst]!=0)]
    x2=table(stat[ii,colIdEst]>0,stat[ii,colIdPV[2]]<pThres)
    x1[match(rownames(x2),c("FALSE","TRUE")),match(colnames(x2),c("FALSE","TRUE"))]=x2
    print(x1)
    
    if (plotFlag=="_onePlot") {
        png(paste("plots",fName,".png",sep=""),width=3*240, height=1*240)
        par(mfcol=c(1,3))
    }
    
    if (plotFlag=="") {
        png(paste("qqplot",fName,".png",sep=""))
        header=compName
    } else {
        header=""
    }
    pvs <- sort(na.exclude(stat[i,colIdPV[1]]))
    qqplot(-log10(runif(length(pvs),0,1)),-log10(pvs),xlab="Expected -log10(p-values) by random",ylab="Observed -log10(p-values)",pch=19,cex.axis=1.5,cex.lab=1.5,main=header)
    abline(0,1)
    if (plotFlag=="") {
        dev.off()
    }
    
    if (plotFlag=="") {
        png(paste("histogram",fName,".png",sep=""))
        header=compName
    } else {
        header=""
    }
    hist(stat[i,colIdPV[1]],xlab="P-value",pch=19,cex.axis=1.5,cex.lab=1.5,main=header)
    if (plotFlag=="") {
        dev.off()
    }
    #	title(main=sub(" ReFACTor","\nReFACTor",compName))
    title(main=compName)
    
    if (plotFlag=="") {
        png(paste("volcanoPlot",fName,"_",colIdPV[1],pThres,".png",sep=""))
        header=compName
    } else {
        header=""
    }
    iThis=i
    if (!outlierFlag) {
        x=quantile(abs(stat[iThis,colIdEst]),probs=.95,na.rm=T)
        iThis=iThis[which(abs(stat[iThis,colIdEst])<=x)]
    }
    plot(stat[iThis,colIdEst],-log10(stat[iThis,colIdPV[1]]),xlab="Estimate",ylab="-log10(p-value)",pch=19,cex.axis=1.5,cex.lab=1.5,main=header)
    ii=iThis[stat[iThis,colIdPV[2]]<pThres]
    points(stat[ii,colIdEst],-log10(stat[ii,colIdPV[2]]),col="red")
    if (plotFlag=="") {
        dev.off()
    }
    
    
    if (plotFlag=="_onePlot") {
        #		par(mfrow=c(1,1))
        #		title(main=compName)
        dev.off()
    }
    
    if (F) {
        
        png(paste("plots",fName,".png",sep=""),width=3*240, height=1*240)
        par(mfcol=c(1,3))
        plot(1:6)
        plot(1:6)
        title(main=sub(" ReFACTor","\nReFACTor",compName))
        plot(1:6)
        dev.off()
    }
    
    ii=i[order(stat[i,colIdPV[2]])]
    ii=ii[stat[ii,colIdPV[2]]<pThres]
    tbl=stat[ii,colList]
    names(tbl)=colName
    tbl=cbind(ann[iA2,c("IlmnID","CHR","MAPINFO","Relation_to_UCSC_CpG_Island","nearestgene")][ii,],tbl)
    write.table(tbl, file=paste("stat",fName,"_",colIdPV[2],pThres,".txt",sep=""), append=F,col.names=T,row.names=F, sep="\t",quote=F)
    
    ####################################################################
    ## Gene level summarization of p-values
    ####################################################################
    
    pvFlag=""
    
    k=which(colnames(stat2)==colIdPV[2])
    kk=grep("Mean",colnames(stat2)[k])
    if (length(kk)!=0) k=k[-kk]
    
    fName2=fName
    # pThres=.001; prId=which(ann2$CHR%in%1:22 & !ann2$IlmnID%in%snpVec & !is.na(stat2[,k]) & stat2[,k]<pThres)
    # prId=which(ann2$CHR%in%1:22 & !ann2$IlmnID%in%snpVec & !is.na(stat2[,k]) & stat2[,k]<pThres)
    prId=which(ann2$keep==1 & !is.na(stat2[,k]) & stat2[,k]<pThres)
    
    if (length(prId)!=0) {
        annotSel=ann2[prId,]
        pval=stat2[,k][prId]
        coef=stat2[,colIdEst[1]][prId]
        signGiven=stat2$sLEU[prId]
        
        
        if (F) {
            # Assign polycomb status to each CpG
            tmp1 = strsplit(annotSel$UCSC_RefGene_Accession,";")
            names(tmp1)=paste(1:length(tmp1),":",sep="")
            tmp2 = unlist(tmp1)
            tmp3 = as.numeric(gsub("[:][[:digit:]]*$","",names(tmp2)))
            tmp4 = data.frame(UCSC_REFGENE_ACCESSION=tmp2,
            rowid=tmp3, stringsAsFactors=FALSE)
            load("PolycombComplete-120109.RData")
            tmp5 = merge(polycombTab,tmp4)
            tmp6 = unique(tmp5$rowid)
            annotSel$PcG = rep(0, dim(annotSel)[1])
            annotSel$PcG[tmp6] = 1
        }
        
        # Get an index of CpGs by gene region
        tmp1 = strsplit(annotSel$UCSC_RefGene_Name,";")
        names(tmp1)=paste(1:length(tmp1),":",sep="")
        tmp2 = unlist(tmp1)
        tmp3 = as.numeric(gsub("[:][[:digit:]]*$","",names(tmp2)))
        tmp4 = strsplit(annotSel$UCSC_RefGene_Group,";")
        names(tmp4)=paste(1:length(tmp4),":",sep="")
        tmp5 = unlist(tmp4)
        tmp6 = as.numeric(gsub("[:][[:digit:]]*$","",names(tmp5)))
        all(tmp3==tmp6)
        
        if (F) {
            library(qvalue)
            qval = qvalue(pval)
            qval$pi0
            qThresh = max(pval[qval$qv<=0.05])
            qThresh
        }
        
        GeneAnnotation = unique(data.frame(UCSC_REFGENE_NAME=tmp2,UCSC_REFGENE_GROUP=tmp5,
        rowid=tmp3, stringsAsFactors=FALSE))
        GeneIndex = split(GeneAnnotation$rowid,with(GeneAnnotation,paste(UCSC_REFGENE_NAME,UCSC_REFGENE_GROUP,sep=":")))
        GeneIndexN = sapply(GeneIndex, length)
        
        if (length(GeneIndexN)!=0) {
            medPval = sapply(GeneIndex, function(u) median(pval[u],na.rm=T))
            propHit = sapply(GeneIndex, function(u) mean(pval[u]<pThres,na.rm=T))
            
            #		isPcG = sapply(GeneIndex, function(u) min(annotSel$PcG[u]))
            #		isNearPcG = sapply(GeneIndex, function(u) max(annotSel$PcG[u]))
            
            tmp1 = sapply(GeneIndex, function(u) u[which.min(annotSel$MAPINFO[u])])
            tmp2 = sapply(GeneIndex, function(u) u[which.max(annotSel$MAPINFO[u])])
            GeneSym1 = annotSel$UCSC_RefGene_Name[tmp1]
            GeneSym2 = annotSel$UCSC_RefGene_Name[tmp2]
            
            annotSelMap <- as.numeric(annotSel$MAPINFO)
            tmpF <- function(u) {
                sgn = sign(coef[u])
                pv = pval[u]
                sgnChar = ifelse(pv>pThres, ".", ifelse(sgn<0,"-","+"))
                paste(sgnChar[order(annotSel$CHR[u], annotSelMap[u])],collapse="")
            }
            # check
            #			if (length(GeneIndex)>2) print(tmpF(GeneIndex[[3]]))
            
            hypohyper =sapply(GeneIndex, tmpF)
            
            combineSigns <- function(u) {
                sgn = signGiven[u]
                sgnChar = ifelse(sgn==0, ".", ifelse(sgn<0,"-","+"))
                paste(sgnChar[order(annotSel$CHR[u], annotSelMap[u])],collapse="")
            }
            hypohyper2 =sapply(GeneIndex, combineSigns)
            
            tmp = strsplit(names(GeneIndexN),":")
            #GeneResults = data.frame(Gene=sapply(tmp,function(u)u[1]),Region=sapply(tmp,function(u)u[2]), nCpG=GeneIndexN, Sign=hypohyper, SignLeuk=hypohyper2, medPval, propHit, stringsAsFactors=FALSE, GeneSymsFirst=GeneSym1, GeneSymsLast=GeneSym2)
            GeneResults = data.frame(Gene=sapply(tmp,function(u)u[1]),Region=sapply(tmp,function(u)u[2]), nCpG=GeneIndexN, Sign=hypohyper, medPval, propHit, stringsAsFactors=FALSE, GeneSymsFirst=GeneSym1, GeneSymsLast=GeneSym2)
            
            rownames(GeneResults)=NULL
            ord = order(medPval)
            
            if (pThres>1) {
                fName=paste("geneSummary_top500_",ifelse(pvFlag=="",colIdPV[2],"pvPerm"),fName2,".txt",sep="")
                write.table(GeneResults[ord[1:500],which(!names(GeneResults)%in%c("isPcG","isNearPcG"))], file=fName, append=FALSE,col.names=T,row.names=FALSE, sep="\t",quote=FALSE)
                
                fName=paste("geneSummary_med",ifelse(pvFlag=="",colIdPV[2],"PvPerm"),".001",fName2,".txt",sep="")
                write.table(GeneResults[ord[medPval[ord]<.001],which(!names(GeneResults)%in%c("isPcG","isNearPcG"))], file=fName, append=FALSE,col.names=T,row.names=FALSE, sep="\t",quote=FALSE)
                
                fName=paste("geneSummary_med",ifelse(pvFlag=="",colIdPV[2],"PvPerm"),".05",fName2,".txt",sep="")
                write.table(GeneResults[ord[medPval[ord]<.05],which(!names(GeneResults)%in%c("isPcG","isNearPcG"))], file=fName, append=FALSE,col.names=T,row.names=FALSE, sep="\t",quote=FALSE)
            } else {
                fName=paste("geneSummary_",ifelse(pvFlag=="",colIdPV[2],"pvPerm"),pThres,fName2,".txt",sep="")
                write.table(GeneResults[ord,which(!names(GeneResults)%in%c("isPcG","isNearPcG"))], file=fName, append=FALSE,col.names=T,row.names=FALSE, sep="\t",quote=FALSE)
            }
        }
        
        tbl=GeneResults[ord,which(!names(GeneResults)%in%c("isPcG","isNearPcG"))]
        
        i=prId[grep("FAM5C",ann2$UCSC_RefGene_Name[prId])]
        i=i[grep("5'UTR",ann2$UCSC_RefGene_Group[i])]
        tbl[which(tbl$Gene=="FAM5C"),]
        stat2[i,]
    }
    
    ####################################################################
    
}


####################################################################
####################################################################
####################################################################
####################################################################
