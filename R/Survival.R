


classify <- function(x,q= 2/3, ...){
    factor(ifelse(x > quantile(x, q, ...),"High","Low"),levels=c("Low","High"))
}

simpleCoxPH <- function(frml, data=clinical){
    coxresult <- summary(coxph(frml, data=data))
    simpleresult <- data.frame(variable=row.names(coxresult$coefficients),coxresult$coefficients[,c(2,4,5),drop=F],coxresult$conf.int[,3:4,drop=F])
    colnames(simpleresult) <- c("variable","HR","z","p.val","lower95","upper95")
    return(simpleresult)
}

plotSimpleCoxPH<-function(overView, maxHR=3, plotTitle="",wrapformula=as.formula("~Outcome"), wrapcols=1, forcelevels=NULL, scl="free",
                          colours = c("grey40","coral4","coral2","coral")){
    require(ggplot2)

    overView$lower95<-sapply(overView$lower95,function(x){ifelse(x<1,(-1/x)+1,x-1)})
    overView$upper95<-sapply(overView$upper95,function(x){ifelse(x<1,(-1/x)+1,x-1)})
    overView$HR<-sapply(overView$HR,function(x){ifelse(x<1,(-1/x)+1,x-1)})

    lbl=c(paste0("1/",seq(maxHR,1,-1)+1),seq(0,maxHR,1)+1)
    overView$significant<-factor(ifelse(overView$p.val<0.001,"***",
                                        ifelse(overView$p.val<0.01, "**",
                                               ifelse(overView$p.val<0.05, "*","NS"))),levels=c("NS","*","**","***"))

    if(!is.null(forcelevels)){overView$variable = factor(overView$variable, levels = forcelevels)}

    overView<-overView[order(overView$HR,decreasing = F),]

    p<-ggplot(overView, aes(x=variable, y=HR,col=significant))+
        geom_point(position=position_dodge(0.9),stat="identity")+
        geom_errorbar(aes(ymin = lower95, ymax = upper95))+
        coord_flip(ylim=c(-maxHR,maxHR))+
        scale_y_continuous(breaks=seq(-maxHR,maxHR,1),labels=lbl)+
        scale_colour_manual(values=colours)+
        geom_hline(yintercept = 0, col="darkgrey")+
        theme(panel.margin = unit(2, "lines"))+
        ggtitle(plotTitle)+xlab("")+ylab("Hazard Ratio (95% CI)")+
        facet_wrap(wrapformula, ncol=wrapcols, strip.position = ifelse(wrapcols==1,"right","top"), scales = scl )
    return(p)
}

HR_plot <- function(in.df=clinical.train, probCol="normProb", survCol="OverallSurvival", minimum_group_size=2, add=F, extraFormula = "", varofinterest=1, symmetric=F, showStats=F, ...){
    in.df <- in.df[order(in.df[,probCol]),]

    # top and bottom n patients are always proficient / deficient
    probs=in.df[-c(1:minimum_group_size,(nrow(in.df)-minimum_group_size):nrow(in.df)),probCol]

    result=lapply(probs, function(pr){
        in.df$Class = factor(ifelse(in.df[,probCol] <= pr,"Proficient","Deficient"), levels=c("Proficient","Deficient"))
        res = simpleCoxPH(as.formula(paste0(survCol,"~ Class",extraFormula)), data = in.df)
        res$cutoff = pr
        return(res[varofinterest,,drop=F])
    })
    result = plyr::ldply(result,rbind)

    # display hazard ratios <1 on the same scale as those >1
    if(symmetric){
        result$HR<-sapply(result$HR,function(x){ifelse(x<1,(-1/x)+1,x-1)})
    }

    pchs = ifelse(result$p.val<0.05,19,1)
    sizes= ifelse(result$p.val<0.05,0.8,0)

    if(!add){
        plot(y=result$HR, x=result$cutoff, pch=pchs, type="l", cex=sizes, ylab="Hazard Ratio", xlab="cut-off used",
             yaxt = ifelse(symmetric,'n','s'), xaxt="n", ...)
        axis(1, at = seq(0,1,0.1), labels=seq(0,1,0.1))
        if(symmetric){
            # for symmmetric plots, suppress standard y axis labels and replace with my own:
            axis(2, at = seq(-200, 200, by = 1), las=2, labels= c(paste0("1/",seq(201,2,-1)),1:201) )
        }
        lines(y=result$HR, x=result$cutoff, pch=pchs, type="p", cex=sizes, ...)
        abline(v=quantile(in.df[,probCol],probs = seq(0,1,0.25)), col="grey30", lty=2)
        abline(h=seq(-50,50), col="grey90", lty=2)
        abline(h=ifelse(symmetric,0,1), col="grey30", lty=1)
    }else{
        lines(y=result$HR, x=result$cutoff, pch=pchs, type="l", cex=sizes, ...)
        lines(y=result$HR, x=result$cutoff, pch=pchs, type="p", cex=sizes, ...)
    }

    # return nr of significant results and median HR
    if(mean(result$p.val<0.05)==0 & showStats){
        return(data.frame(fractionSignificant=0,medianHR=1))
    }else if(showStats){
        return(data.frame(fractionSignificant=mean(result$p.val<0.05), medianHR=median(result$HR[which(result$p.val<0.05)])))
    }
}

rankclassify = function(x, na.rm=T, q){
    # q = "proportion high"; If you do this with quantiles then everything gets shifted if you have a lot of zeros or something? Only matters for Sig3 here
    normProb = norm_01(rank(x))
    factor(ifelse(normProb <= q,"Low","High"),levels=c("Low","High"))
}

plotSurvCutoff = function(in.df,survCol,sigCol,cutoff=0.5,...){
    in.df$mysplit = rankclassify(in.df[,sigCol], q=cutoff)
    plotSurvCategory(in.df=in.df, survCol = survCol, sigCol="mysplit", ...)
}



# for using reference levels according to the radplat cohort (so that both cohorts have the same references in the end)
table_and_cox <- function(var=reportcols[1], clinical=clinical.all, clinref = NULL, outcome="OverallSurvival"){
    if(is.null(clinref)){clinref=clinical}

    df1 = data.frame(Variable = var, table(clinical[,var]))

    # add "Unknowns" to the table
    if(sum(is.na(clinical[,var]))>0){
        df1=rbind(df1,c(var,NA,sum(is.na(clinical[,var]))))
    }

    # add reference level of the reference cohort.
    df_ref = data.frame(Variable = var, table(clinref[,var]))
    clinical[,var]=relevel(factor(clinical[,var]), ref= as.character(df_ref$Var1[which.max(df_ref$Freq)]))

    df2 = simpleCoxPH(as.formula(paste0(outcome, " ~ ",var)), clinical)
    df2$Var1 = gsub(var,"",df2$variable)
    df3 = merge(df1,df2[c("Var1","HR","p.val")], by="Var1", all.x=T)[,c(2,1,3:5)]
    colnames(df3)[2]="Value"
    return(df3)
}


survival2category<-function(response, time=2.5*365.25){
    # assumes expr and response are in the same order!
    response<-as.matrix(response)
    message(paste0("Response Dimensions: ", paste(dim(response),collapse=", ")))
    class<-ifelse(response[,2]==0 & response[,1]>time,"Negative",
                  ifelse(response[,2]==0 & response[,1]<time,"Censored",
                         ifelse(response[,2]==1 & response[,1]>time,"Negative","Positive")))
    class<-factor(class)

    return(class)
}


plotInteractionKM = function(df, first, second, q1=0.5, q2=0.5, outcomes.surv = outcomes, showCoxPH=F, shortnames=NULL, ...){

    require(SurvivalPlots)

    if(is.null(shortnames)){
        shortnames=c(first,second)
    }

    df = subset(df, !is.na(df[,first]) & !is.na(df[,second]))
    df$interaction = interaction(paste0(shortnames[1] ,as.character(classify(df[,first ],q=q1, na.rm=T))),
                                 paste0(shortnames[2] ,as.character(classify(df[,second],q=q2, na.rm=T))))

    df[,first ] = rankclassify(df[,first ],q=q1, na.rm=T)
    df[,second] = rankclassify(df[,second],q=q2, na.rm=T)
    interactions=list()

    par(mfrow=c(1,4))
    for(outcome in outcomes.surv){
        plotSurvCategory(df, sigCol = "interaction", survCol=outcome, title = outcome, ...)

        interact = simpleCoxPH(as.formula(paste0(outcome, "~",first,"*",second)))
        interactions[[outcome]] = interact

        text(y=0.2,x=0.5*365.25,labels=paste0("interaction p=",round.pval(interact[3,4])), adj=0)
    }

    if(showCoxPH){
        interactions = plyr::ldply(interactions)
        colnames(interactions)[1]="Outcome"
        return(interactions)
    }
}

