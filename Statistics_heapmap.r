library(openxlsx)
library(reshape2)
library(pheatmap)
library(RColorBrewer)


data<-list(PLacylreport=read.xlsx("~/Documents/Excel_rawdata/Disease_acylfinal.xlsx",1),
             PLreport=read.xlsx("~/Documents/Excel_rawdata/DiseaseMolfinal.xlsx",1))


#PLacylreport = AfterAcyl


colnames(data[["PLacylreport"]])[1]<-"ID"
colnames(data[["PLreport"]])[1]<-"ID"

row.names(data[["PLacylreport"]])<-data[["PLacylreport"]]$ID
row.names(data[["PLreport"]])<-data[["PLreport"]]$ID


expdesign<-list(skinny=paste0("TMT",c(126:128)),
                obese=paste0("TMT",c(129:131)))


for(m in c("PLacylreport","PLreport")){
 # dir.create(paste0("~Document/results/",m))
  for(d in c("raw")){ ## correctbytotal can be added to enable correction
    if(d=="raw"){
      workingdata<-data[[m]]
    }else{
      workingdata<-data[[m]]
      workingdata[unlist(expdesign)]<-workingdata[unlist(expdesign)]*workingdata$PreInt
      workingdata[unlist(expdesign)]<-t(t(workingdata[unlist(expdesign)])/
                                          apply(log10(workingdata[unlist(expdesign)]),2,
                                                function(x)sum(x[!is.infinite(x)])))
      ##the two ts is because R process through rows not columns
    }
    for(g in unlist(expdesign)){
      if(any(unlist(workingdata[,g])==0)){
        workingdata[unlist(workingdata[,g])==0,g]<-NA
      }
    }
    workingdata[unlist(expdesign)]<-t(apply(log10(workingdata[unlist(expdesign)]),1,function(x)x-mean(x,na.rm = T))) ##normalize with log mean
    testresult<-data.frame(log10diff=NULL,t.test.p=NULL) ## Not enough treatment for mann-whitney
    for(x in 1:nrow(workingdata)){
      log10diff<-mean(unlist(workingdata[x,expdesign[["obese"]]]),na.rm = T)-
        mean(unlist(workingdata[x,expdesign[["skinny"]]]),na.rm = T)
      t.test.p<-try(t.test(unlist(workingdata[x,expdesign[["skinny"]]]),
                       unlist(workingdata[x,expdesign[["obese"]]]),paired=F,na.rm = T,)$p.value,silent = T)
      if(!is.numeric(t.test.p)){
        t.test.p<-NA
      }
      testresult<-rbind(testresult,c(log10diff,t.test.p))
    }
    
    colnames(testresult)<-c("log10diff","t.test.p")
    row.names(testresult)<-row.names(workingdata)
    
    
    testresult$bonferroni <- p.adjust(testresult$t.test.p,method = "bonferroni")
    testresult$fdr <- p.adjust(testresult$t.test.p,method = "fdr")
    
    writexlsx<-cbind(row.names(testresult),testresult)
    colnames(writexlsx)[1]<-m
    write.xlsx(writexlsx,
               paste0("~/Documents/results/PLacyltest_",m,"_",d,".xlsx"))
    
    
    group<-pheatmap(as.matrix(workingdata[unlist(expdesign)]),silent = T)
    
    cluster<-hclust(dist((as.matrix(workingdata[unlist(expdesign)]))))
    
    k=ifelse(m=="PL",4,
             ifelse(m=="Isomer",3,
                    ifelse(m=="PLacylreport",5,4)))
    
    cut_avg <- cutree(cluster, k = k)
    clusterdf<-data.frame(ID=row.names(workingdata)[group$tree_row$order])
    clusterdf$clusters<-cut_avg[clusterdf$ID]
    write.xlsx(clusterdf,
               paste0("~/Documents/results/PLacyltestcluster_",m,"_",d,".xlsx"))
    
    annotation_col = data.frame(
      groups = factor(c(rep("skinny",3),rep("obese",3))))
    rownames(annotation_col) = paste0("TMT",c(126:131))
    
    annotation_row = data.frame(
      LipidClusters = factor(paste0("C",cut_avg))
    )
    rownames(annotation_row) = names(cut_avg)
    
    ann_colors = list(
      groups = rainbow(2),
      LipidClusters = rainbow(k))
    names(ann_colors[["groups"]])<-c("skinny", "obese")
    names(ann_colors[["LipidClusters"]])<-unique(paste0("C",cut_avg))
    
    pdf(paste0("~/Documents/results/PLacyltestclustermap_",m,"_",d,".pdf"),height = 10,width=7)
    pheatmap(as.matrix(workingdata[unlist(expdesign)]), 
             annotation_col = annotation_col, 
             annotation_row = annotation_row, 
             annotation_colors = ann_colors, 
             row_split = annotation_row$GeneClass,
             column_split = annotation_col$CellType,
             col=brewer.pal(9,"GnBu"),na_col = "grey70")
    dev.off()
  }

    
}
