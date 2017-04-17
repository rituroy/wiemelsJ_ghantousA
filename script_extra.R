## Run script.R
## This is just extra stuff

####################################################################
####################################################################
## Section 2
## This hs been replaced by linReg_moba.R

###################################################  Models

################################## For models without interaction, use RLMtestmodif function.
########### For models with interaction, use RMLtestinteraction
########### Label the qqplots and significant CpG sets identified below using the coding given for the ########### regression models above (Section VIII)

RLMtestmodif = function(meth_matrix,methcol,exposure, covmatrix=NULL) {
    expr<-c()
    if (!is.null(covmatrix)) {
        for (i in 1:ncol(covmatrix)) {
            assign(paste0("X",i),covmatrix[,i])
            expr<-paste0(expr,"+X",i)
        }
    }
    eval(parse(text=paste0("mod <- try(rlm(meth_matrix[, methcol]~exposure",expr,",maxit=200))")))
    if(class(mod) == "try-error") {
        print(paste("error thrown by column", methcol))
        invisible(rep(NA, 3))
    } else {
        cf = coeftest(mod, vcov=vcovHC(mod, type="HC0"))
    }
    cf[2, c("Estimate", "Std. Error", "z value", "Pr(>|z|)")]
}


RLMtestinteraction = function(meth_matrix,methcol,exposure, covmatrix=NULL) {
    expr <- c()
    if (!is.null(covmatrix)) {
        for (i in 1:ncol(covmatrix)) {
            assign(paste0("X",i),covmatrix[,i])
            expr<-paste0(expr,"+X",i)
        }
    }
    expr <- paste0(expr,"+exposure*X1")
    eval(parse(text=paste0("mod <- try(rlm(meth_matrix[, methcol]~exposure",expr,",maxit=200))")))
    if(class(mod) == "try-error"){
        print(paste("error thrown by column", methcol))
        invisible(rep(NA, 3))
    } else {
        cf <- coeftest(mod, vcov=vcovHC(mod, type="HC0"))
    }
    c(cf[2, c("Estimate", "Std. Error","z value","Pr(>|z|)")],cf[4, "Pr(>|z|)"])
}

##covmatrix <- data.frame(pdata[,5]) indicate indexes of the covariates; if none, then delete the term #########”covmatrix” from mclapply below
system.time(ind.res <- mclapply(setNames(seq_len(ncol(tdat)), dimnames(tdat)[[2]]), RLMtestinteraction, meth_matrix=tdat,
exposure <- pdata[,10],covmatrix,mc.cores=12))


system.time(ind.res <- mclapply(setNames(seq_len(ncol(tdat)), dimnames(tdat)[[2]]), RLMtestinteraction, meth_matrix=tdat,
exposure <- pdata[,10],covmatrix,mc.cores=12))


#### for RLMtestmodif use this
setattr(ind.res, 'class', 'data.frame')
setattr(ind.res, "row.names", c(NA_integer_,5))
setattr(ind.res, "names", make.names(names(ind.res), unique=TRUE))
probelistnames <- names(ind.res)

all.results1 <- t(data.table(ind.res))
all.results1 <- data.table(all.results1)
all.results1[, probeID := probelistnames]
setnames(all.results1, c("BETA","SE","z value", "P_VAL", "probeID")) # rename columns
setcolorder(all.results1, c("probeID","BETA","SE","z value", "P_VAL"))
rm(probelistnames, ind.res)
#######################

#### for RLMtestinteraction use this
setattr(ind.res, 'class', 'data.frame')
setattr(ind.res, "row.names", c(NA_integer_,6))
setattr(ind.res, "names", make.names(names(ind.res), unique=TRUE))
probelistnames <- names(ind.res)

all.results1 <- t(data.table(ind.res))
all.results1 <- data.table(all.results1)
all.results1[, probeID := probelistnames]
setnames(all.results1, c("BETA","SE","z value", "P_VAL", "Inter_P_VAL","probeID")) # rename columns
setcolorder(all.results1, c("probeID","BETA","SE","z value", "P_VAL", "Inter_P_VAL"))
rm(probelistnames, ind.res)
#######################


####################################################################
####################################################################
