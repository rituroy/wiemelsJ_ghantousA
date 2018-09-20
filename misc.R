
####################################################################
####################################################################
## Statistics for manuscript

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
    dirRes="results/moba/allSamples/output/stat/"
    dirRes="results/moba/noHispWt/output/stat/"
}

## ----------------------------------------------
if (computerFlag=="cluster") {
    ann=read.delim(paste("data/HumanMethylation450_15017482_v.1.2.csv",sep=""),header=TRUE, sep=",",quote="",comment.char="",as.is=T,fill=T, skip=7)
    snpVec=read.table(paste("data/list_to_exclude_Sept_24.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
} else {
    ann=read.delim(paste("docs/yuanyuan/HumanMethylation450_15017482_v.1.2.csv",sep=""),header=TRUE, sep=",",quote="",comment.char="",as.is=T,fill=T, skip=7)
    snpVec=read.table(paste("docs/SemiraGonsethNussle/list_to_exclude_Sept_24.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    #dmr=read.table(paste("docs/SeungTae/leukemia.DMRs/leukemia.DMRs.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
}
ann[which(ann[,"CHR"]=="X"),"CHR"]="23"
ann[which(ann[,"CHR"]=="Y"),"CHR"]="24"
ann[,"CHR"]=as.integer(ann[,"CHR"])
ann=ann[,-match(c("AddressA_ID","AlleleA_ProbeSeq","AddressB_ID","AlleleB_ProbeSeq", "Next_Base",  "Color_Channel","Forward_Sequence","SourceSeq"),colnames(ann))]
for (k in 1:ncol(ann)) if (class(ann[,k])=="factor") ann[,k]=as.character(ann[,k])

snpVec=snpVec[,1]
ann$snp=0; ann$snp[which(ann$IlmnID%in%snpVec)]=1
ann$keep=0
ann$keep[which(ann$snp==0 & ann$CHR%in%1:22)]=1

## ------------------------
# 03/02/18
if (computerFlag=="cluster") {
    dirClin=""
} else {
    dirClin="docs/moba/forManuscript/"
}

candGene=read.table(paste(dirClin,"CpGlist_DMR.TXT",sep=""),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
candGene=gsub("\"| |`","",candGene[,1])
#candGene=gsub(" ","",gsub("\"","",candGene[,1]))
candGene=candGene[which(substr(candGene,1,1)!="#")]
k1=grep("<-",candGene)
k2=c(k1[2:length(k1)]-1,length(candGene))
out=rep("",length(k1))
for (k in 1:length(k1)) {
    out[k]=paste(candGene[k1[k]:k2[k]],collapse="")
}
x=c()
for (k in 1:length(out)) {
    kk=grep("<-",out[k])
    if (length(kk)!=1) x=c(x,k)
}
candGene=NULL
for (k in 1:length(out)) {
    x=as.list(strsplit(sub(")","",out[k]),"<-c(",fixed=T)[[1]])
    #x[[1]]=gsub("-","_",x[[1]])
    x[[2]]=strsplit(x[[2]],",")[[1]]
    tbl2=data.frame(geneSym=rep(x[[1]],length(x[[2]])),cpgId=x[[2]],stringsAsFactors=F)
    candGene=rbind(candGene,tbl2)
}

write.table(candGene,file="CpGlist_DMR_format2.txt", sep="\t", col.names=T, row.names=F, quote=F)

candGeneDmr=candGene

## ---------------------------------
# 08/13/18
datadir="docs/moba/forManuscript/"
candGeneDmrRlm=read.table(paste(datadir,"CpGlist_DMR_RLM.txt",sep=""),sep="\n",h=F,quote="",comment.char="",as.is=T,fill=T)
candGeneReg=read.table(paste(datadir,"Annot_Regul_reg.txt",sep=""),sep=" ",h=T,quote="",comment.char="",as.is=T,fill=T)
candGeneDen=read.table(paste(datadir,"Annot_Density.txt",sep=""),sep=" ",h=T,quote="",comment.char="",as.is=T,fill=T)
dim(candGeneDmrRlm)
dim(candGeneReg)
dim(candGeneDen)

names(candGeneReg)=names(candGeneDen)=c("cpgId","loc")

tbl=candGeneDmrRlm[,1]
tbl=gsub("\"","",tbl)
for (k in 1:ncol(candGeneReg)) candGeneReg[,k]=gsub("\"","",candGeneReg[,k])
for (k in 1:ncol(candGeneDen)) candGeneDen[,k]=gsub("\"","",candGeneDen[,k])

k=grep("# ",tbl)
k1=k+1; k2=c(k[2:length(k)]-1,length(tbl))
candGeneDmrRlm=NULL
for (k in 1:length(k1)) {
    modelThis=sub("# ","",tbl[k1[k]-1])
    tbl2=NULL
    for (i in k1[k]:k2[k]) {
        x=tbl[i]
        y=sub(")","",x)
        y=strsplit(y,"<-c(",fixed=T)[[1]]
        geneSym=y[1]
        cpgId=strsplit(y[2],", ")[[1]]
        model=rep(modelThis,length(cpgId))
        geneSym=rep(geneSym,length(cpgId))
        tbl2=rbind(tbl2,data.frame(model,geneSym,cpgId,stringsAsFactors=F))
    }
    candGeneDmrRlm=rbind(candGeneDmrRlm,tbl2)
}

####################################################################
####################################################################
## Get data before SVA

if (F) {
    computerFlag=""
    computerFlag="cluster"

    ## ------------------------
    if (computerFlag=="cluster") {
        dirClin=""
    } else {
        dirClin="docs/moba/forManuscript/"
    }

    subsetFlag="noHispWt"
    subsetFlag="allSamples"

    nProbe=101
    nProbe=-1

    subsetList=c("allSamples","noHispWt")
    for (subsetFlag in subsetList) {
        tdat=read.table(paste(dirClin,subsetFlag,"/betas.filtered.csv",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=nProbe)
        pdata=read.table(paste(dirClin,subsetFlag,"/pdata.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

        colnames(tdat)=sub("X","",gsub(".","",colnames(tdat),fixed=T))
        rownames(tdat)=gsub("\"","",rownames(tdat))
        colnames(pdata)=sub("X","",gsub(".","",colnames(pdata),fixed=T))
        tdat=as.matrix(tdat)
        for (k in 1:ncol(pdata)) {
            if (is.character(pdata[,k])) pdata[,k]=gsub("\"","",pdata[,k])
        }
        table(colnames(tdat)==paste(pdata$Beadchip,"_",pdata$Position,sep=""))

        save(tdat,pdata,file=paste(dirClin,subsetFlag,"/data.RData",sep=""))
    }
}

####################################################################
####################################################################
## Candidate gene summary
datFlag=""
datFlag="_beforeSVA"

candGeneList="_dmr"
candGeneList=c("_dmrRlm","_reg")
candGeneList=c("_dmrRlm","_reg","_den")

varList=c("caco","sex")
varList=c("caco","sex","cacoSex")
varList=c("caco","cacoSex","cacoEthn","cacoSexEthn")
varList=c("cacoSex")
varList=c("caco","cacoSex","cacoEthn","cacoSexEthn")

subsetList=c("noHispWt")
subsetList=c("allSamples","noHispWt")

if (datFlag=="_beforeSVA") datadir=paste("allSamples/",sep="") else datadir=paste("results/moba/allSamples/misc/",sep="")
load(paste(datadir,"data.RData",sep=""))
varInfo1=list(caco=sort(unique(pdata$caco)),sex=sort(unique(pdata$sex)),"ethn"=sort(unique(pdata$int_ch_ethnicity)),"cacoSex"=sort(unique(paste(pdata$caco,pdata$sex))),"cacoEthn"=sort(unique(paste(pdata$caco,pdata$int_ch_ethnicity))),"cacoSexEthn"=sort(unique(paste(pdata$caco,pdata$sex,pdata$int_ch_ethnicity))))
varInfo2=list(caco=c("ctrl","case"),sex=c("male","female"),"ethn"=c("hispanic","nonHispanicWhite","nonHispanicOther"),"cacoSex"=paste(c("male","female"),"_",rep(c("ctrl","case"),each=2),sep=""),"cacoEthn"=paste(c("hispanic","nonHispanicWhite","nonHispanicOther"),"_",rep(c("ctrl","case"),each=3),sep=""),"cacoSexEthn"=paste(c("hispanic","nonHispanicWhite","nonHispanicOther"),"_",rep(paste(c("male","female"),"_",rep(c("ctrl","case"),each=2),sep=""),each=3),sep=""))

for (candGeneFlag in candGeneList) {
    switch(candGeneFlag,
        "_dmr"={
            header="DMR"
            candGene=candGeneDmr
            colInfo=data.frame(id=c("id","gene"),value=c("cpgId","geneSym"),stringsAsFactors=F)
            colId2="gene"; colId1=colInfo$value[colInfo$id=="gene"]
            colId2="cpgId"; colId1=colInfo$value[colInfo$id=="id"]
        },
        "_dmrRlm"={
            header="DMR-RLM"
            candGene=candGeneDmrRlm
            candGene$model=paste(unique(candGene$model),collapse=" or ")
            candGene=rbind(candGeneDmrRlm,candGene)
            colInfo=data.frame(id=c("group","id","gene"),value=c("model","cpgId","geneSym"),stringsAsFactors=F)
            colId2="cpgId"; colId1=colInfo$value[colInfo$id=="id"]
        },
        "_reg"={
            header="Regulatory region"
            candGene=candGeneReg
            colInfo=data.frame(id=c("id","gene"),value=c("cpgId","loc"),stringsAsFactors=F)
            colId2="loc"; colId1=colInfo$value[colInfo$id=="gene"]
        },
        "_den"={
            header="Density"
            candGene=candGeneDen
            colInfo=data.frame(id=c("id","gene"),value=c("cpgId","loc"),stringsAsFactors=F)
            colId2="loc"; colId1=colInfo$value[colInfo$id=="gene"]
        }
    )
    for (vId in 1:length(varList)) {
        tbl1=tbl2=NULL
        for (subsetFlag in subsetList) {
            cat("\n\n",subsetFlag,"\n")
            if (datFlag=="_beforeSVA") datadir=paste(subsetFlag,"/",sep="") else datadir=paste("results/moba/",subsetFlag,"/misc/",sep="")
            load(paste(datadir,"data.RData",sep=""))
            pdata$cacoSex=paste(pdata$caco,pdata$sex)
            pdata$cacoEthn=paste(pdata$caco,pdata$int_ch_ethnicity)
            pdata$cacoSexEthn=paste(pdata$caco,pdata$sex,pdata$int_ch_ethnicity)
            if (datFlag=="") {
                tdat=t(tdat)
            }
            #print("table(pdata$beadchip_position==colnames(tdat))")
            #print(table(paste(pdata$Beadchip,"_",pdata$Position,sep="")==colnames(tdat)))
            grp=pdata[,varList[vId]]
            grpUniq=sort(unique(grp))
            switch(varList[vId],
                "caco"={ttl=c("ctrl","case")
                },
                "sex"={ttl=c("male","female")
                },
                "ethn"={ttl=c("hispanic","nonHispanicWhite","nonHispanicOther")
                },
                "cacoSex"={ttl=paste(c("male","female"),"_",rep(c("ctrl","case"),each=2),sep="")
                },
                "cacoEthn"={ttl=paste(c("hispanic","nonHispanicWhite","nonHispanicOther"),"_",rep(c("ctrl","case"),each=3),sep="")
                },
                "cacoSexEthn"={ttl=paste(c("hispanic","nonHispanicWhite","nonHispanicOther"),"_",rep(paste(c("male","female"),"_",rep(c("ctrl","case"),each=2),sep=""),each=3),sep="")
                }
            )
            ttl=ttl[match(grpUniq,varInfo1[[varList[vId]]])]
            cat("\n\nNo. of samples:\n")
            for (gId in 1:length(grpUniq)) {
                cat(ttl[gId],": ",length(grp[which(grp==grpUniq[gId])]),"\n",sep="")
            }
            grp2=candGene[,colId1]
            grp2Uniq=unique(grp2)
            #if (candGeneFlag%in%c("_dmr","_dmrRlm")) {
            if (candGeneFlag%in%c("_dmr","_dmrRlm","_reg","_den")) {
                n=length(grp2Uniq)+1
                tmp=rep(NA,n); tmpC=rep("",n)
                tbl=data.frame(cpgId=c("global",grp2Uniq),stringsAsFactors=F)
                if (colId2=="cpgId") {
                    nm="gene"
                    tbl=cbind(gene=tmpC,tbl,stringsAsFactors=F)
                } else {
                    nm=c()
                }
                nm=c(nm,colId2,paste(c("mean","sd"),"_",rep(ttl,each=2),sep=""))
            } else {
                n=1
                tmp=rep(NA,n); tmpC=rep("",n)
                nm=c("gene","cpgId",paste(c("mean","sd"),"_",rep(ttl,each=2),sep=""))
                tbl=data.frame(gene=tmpC,cpgId=c("global"),stringsAsFactors=F)
            }
            for (gId in 1:length(grpUniq)) tbl=cbind(tbl,tmp,tmp)
            if (F) {
                if (varList[vId]%in%c("cacoSex","cacoEthn")) {
                    tbl=data.frame(gene=tmpC,cpgId=c("global",grp2Uniq),tmp,tmp,tmp,tmp,tmp,tmp,tmp,tmp,stringsAsFactors=F)
                } else {
                    tbl=data.frame(gene=tmpC,cpgId=c("global",grp2Uniq),tmp,tmp,tmp,tmp,stringsAsFactors=F)
                }
            }
            names(tbl)=nm
            #names(tbl)[1]=colId2
            if ("gene"%in%names(tbl)) {
                i=match(candGene$cpgId,tbl$cpgId); i1=which(!is.na(i)); i2=i[i1]
                tbl$gene[i2]=candGene[i1,colInfo$value[colInfo$id=="gene"]]
            }
            k=which(tbl[,colId2]=="global")
            for (gId in 1:length(grpUniq)) {
                p=which(names(tbl)==paste("mean_",ttl[gId],sep=""))
                x=c(tdat[,which(grp==grpUniq[gId])])
                tbl[k,p]=mean(x,na.rm=T)
                tbl[k,p+1]=sd(x,na.rm=T)
            }
            #if (candGeneFlag%in%c("_dmr","_dmrRlm")) {
            if (candGeneFlag%in%c("_dmr","_dmrRlm","_reg","_den")) {
                i=match(candGene$cpgId,rownames(tdat)); i1=which(!is.na(i)); i2=i[i1]
                for (k in which(tbl[,colId2]%in%candGene[i1,colId1])) {
                    for (gId in 1:length(grpUniq)) {
                        p=which(names(tbl)==paste("mean_",ttl[gId],sep=""))
                        x=c(tdat[i2,][which(candGene[i1,colId1]==tbl[k,colId2]),which(grp==grpUniq[gId])])
                        tbl[k,p]=mean(x,na.rm=T)
                        tbl[k,p+1]=sd(x,na.rm=T)
                    }
                }
            }
            switch(subsetFlag,
                "allSamples"={
                    tbl1=tbl
                },
                "noHispWt"={
                    tbl2=tbl
                }
            )
        }
        fName=paste("summaryTbl",datFlag,candGeneFlag,"_",varList[vId],"_ccls.txt",sep="")
        colId=2:ncol(tbl1)
        tbl=tbl1
        nm=c("",subsetList[1],rep("",length(colId)-1))
        if (length(subsetList)>1) {
            colId=2:ncol(tbl2)
            tbl=cbind(tbl,tbl2[,colId])
            nm=c(nm,subsetList[2],rep("",length(colId)-1))
        }
        write.table(paste(nm,collapse="\t"),file=fName, sep="\t", col.names=F, row.names=F, quote=F,append=F)
        #nm=c("gene",rep(names(tbl1)[colId],length(subsetList)))
        nm=names(tbl)
        write.table(paste(nm,collapse="\t"),file=fName, sep="\t", col.names=F, row.names=F, quote=F,append=T)
        write.table(tbl,file=fName, sep="\t", col.names=F, row.names=F, quote=F,append=T)
    }
}

####################################################################
####################################################################
## NOT USED
## Clinical variable summary

if (F) {
    fName1=""

    varList="sex"

    subsetList=c("allSamples","noHispWt")
    for (subsetFlag in subsetList) {
        switch(subsetFlag,
            "allSamples"={
                subsetName="All samples"
            },
            "noHispWt"={
                subsetName="Non-hispanic white"
            }
        )
        cat("\n\n",subsetName,"\n")
        if (fName1=="_beforeSVA") {
            datadir=paste(subsetFlag,"/",sep="")
        } else {
            datadir=paste("results/moba/",subsetFlag,"/misc/",sep="")
        }
        load(paste(datadir,"data.RData",sep=""))
        for (vId in 1:length(varList)) {
            grp=pdata[,varList[vId]]
            switch(varList[vId],
                "sex"={
                    grp=as.character(grp)
                    grp[which(grp=="1")]="Male"
                    grp[which(grp=="2")]="Female"
                }
            )
            grpUniq=unique(grp)
            for (gId in 1:length(grpUniq)) {
                x=grp[which(grp==grpUniq[gId])]
                cat(grpUniq,": ","Mean ",mean(x),", SD ",sd(x),"\n")
            }
        }
    }
}

####################################################################
####################################################################
