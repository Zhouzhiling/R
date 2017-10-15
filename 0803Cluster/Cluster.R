library(cluster);
library(norm);
library(SC3);
library(stats);

rm(list=ls(all=TRUE));

Input <- read.table(file.choose(),head=T,sep="\t",fill=TRUE);
Input = as.matrix(Input);


fun.Calculate <- function(TextNum,KeyWord,CellThresholdUpper,CellThresholdLower,GeneThreshold,ClusterNum)
{
  #TextNum = 3;
  #KeyWord = 2;
  #CellThresholdUpper = 1000000;
  #CellThresholdLower = 100000;
  #GeneThreshold = 5;
  #ClusterNum = 500;
  lncRNA = matrix(as.integer(Input[,-c(1:TextNum)]),nrow=nrow(Input),ncol=ncol(Input)-TextNum);
  rownames(lncRNA)= Input[c(1:nrow(Input)),KeyWord]
  colnames(lncRNA)= colnames(Input)[-c(1:TextNum)]
  #rownames(lncRNA)
  #colnames(lncRNA)
  if(!is.numeric(lncRNA))
  {
    lncRNA[is.na(lncRNA)] <- 0
  }
  #is.numeric(lncRNA)
  
  #去除低质量细胞
  lncRNA <- lncRNA[,colSums(lncRNA)>CellThresholdLower];
  lncRNA <- lncRNA[,colSums(lncRNA)<CellThresholdUpper];
  
  #去除低质量基因
  lncRNA <- lncRNA[rowSums(lncRNA)>GeneThreshold,];
  
  #去除低表达基因：基因至少在5个细胞中表达
  lncRNA.filter <- lncRNA[ sapply(1:dim(lncRNA)[1],function(X) sum(lncRNA[X,]>0))>GeneThreshold,] 
  
  
  #标准化： RC/size_factor
  #for(i in 1:nrow(lncRNA.filter))
  #{
  #  lncRNA.filter.scale[i,] <- lncRNA.filter[i,]/rowSums(lncRNA.filter)[i];			#按行合并
  #}
  
  #标准化
  lncRNA.filter.scale <- scale(lncRNA.filter.scale, center=rep(0,ncol(lncRNA.filter)),scale=F)
  
  #根据方差排序
  varance <- c();
  lncRNA.filter.rank <- matrix(rep(-1,nrow(lncRNA.filter.scale)*ncol(lncRNA.filter.scale)),nrow = nrow(lncRNA.filter.scale),ncol = ncol(lncRNA.filter.scale));
  colnames(lncRNA.filter.rank) <- colnames(lncRNA.filter.scale);
  rownames(lncRNA.filter.rank) <- rownames(lncRNA.filter.scale);
  #for(i in 1:nrow(lncRNA.filter.scale))
  for(i in 1:nrow(lncRNA.filter.scale))
  varance[i] <- var(lncRNA.filter.scale[i,]);
  
  
  orderResult <- order(varance)
  for(i in 1:nrow(lncRNA.filter.scale))
    {
    rownames(lncRNA.filter.rank)[orderResult[i]] <- rownames(lncRNA.filter.scale)[i];
      lncRNA.filter.rank[orderResult[i],] <- lncRNA.filter.scale[i,];
  }
  lncRNA.filter.tocluster <- (lncRNA.filter.rank[c(1:ClusterNum),]);

  
  ##计算最佳k值
  ##Calculate the best K, written into outputClustNum
  #begin = 3; 
  #length = 15;
  #count = 5;
  #end = begin + length - 1;
  #SCNow=-1;
  #result = c();
  # result[begin:end] = -1;
  #for(i in begin:end) {
  #  # Silhouette coefficient
  #  tmp = c();
  #  tmp[1:count] = 0;
  #  for(j in 1:count) {
  #    kcluster = pam(lncRNA.filter.tocluster, i);
  #    tmp[j] = kcluster$silinfo$avg.width;
  #  }
  # result[i] = mean(tmp);
  # if(result[i] > SCNow)
  #  {
  #   SCNow <- result[i];
  #   outputClustNum <- i;
  # }
  #}
  
  
  #SC3
  library(SC3);
  library(scater);
  install.packages("scater");
  lncRNA.clu<-read.table("path here",head=T)
  lncRNA.clu_ann <- read.table("path here",sep="\t")
  rownames(lncRNA.clu_ann) <- lncRNA.clu_ann[,1];
  lncRNA.clu_ann <- lncRNA.clu_ann[,-c(1)]
  head(lncRNA.clu_ann)
  lncRNA.clu<-log2(lncRNA.clu+1)
  lncRNA.clu_pd <- new("AnnotatedDataFrame", data = lncRNA.clu_ann)
  lncRNA.clu_sceset<-newSCESet(fpkmData=as.matrix(lncRNA.clu),phenoData=lncRNA.clu_pd,logExprsOffset = 0)
  is_exprs(lncRNA.clu_sceset) <- exprs(lncRNA.clu_sceset) > 0
  lncRNA.clu_sceset <- calculateQCMetrics(lncRNA.clu_sceset)
  lncRNA.clu_sceset <- sc3(lncRNA.clu_sceset, ks = 2:10, biology = TRUE, n_cores = 1)
  write.table(lncRNA.clu_sceset,file="path here",sep=",");
  
  
  #setwd(ClusterResultPath);
  #jpeg(file="NumberOfCluster.jpeg");
  #plot(result, type="o",main=outputClustNum, xlab="Number of Cluster", ylab="Silhouette Coefficient");
  #dev.off();
  
  #Calculate the cluster
  #cl = kmeans(lncRNA.filter.tocluster, outputClustNum)
  #plot(lncRNA.filter.tocluster, col = cl$cluster)
  
  #ResultOutput <- matrix(rep("",nrow(lncRNA.filter.tocluster)*2),nrow = nrow(lncRNA.filter.tocluster) ,ncol = 2)
  #rownames(ResultOutput)= rownames(lncRNA.filter.tocluster);
  #ResultOutput[,2] <- cl$cluster;
  #ResultOutput <- ResultOutput[,c(2)];
  #write.table(ResultOutput,file="C:/Users/Milkway/Desktop/R/0729/ClusterResult.csv",sep=",")
  
  #DESeq
  #install.packages("E:/source/DESeq_1.28.0.zip", 
  #                 destdir = 'E:/source',
  #                 lib = 'E:/libs', 
  #                 repos = NULL)
  #library('DESeq', lib.loc = 'E:/libs')
  #讲数据放到表中
  #描述实验设置重复
  #ncol(tmpdelete)
  #floor(ncol(tmpdelete)/2)
  #ceiling(ncol(tmpdelete)/2)
  
  library('DESeq');
  flag <- ceiling(ncol(lncRNA.filter)/2) - floor(ncol(lncRNA.filter)/2);
  if(flag==1) {
    design<- rep (c("P","Mo"),each=ceiling(ncol(lncRNA.filter)/2));
    design<-design[-(ncol(lncRNA.filter)+1)];
  } else {
    design<- rep (c("P","Mo"),each=ncol(lncRNA.filter)/2); }
  
  de  <- newCountDataSet(lncRNA.filter, design)
  de  <- estimateSizeFactors(de)
  de  <- estimateDispersions(de)
  res  <- nbinomTest(de,"P","Mo")
  #计算结果
  #sum(na.omit(res$padj<0.05))
  #输出差异显著的值。
  #最后再加一部把数据输出出来
  write.table(res,file="path here",sep=",");
  

fun.Calculate(3,2,1000000,100000,5,500);
#参数一n：数据的前n列是文本信息，处理时需要删除
#参数二m：可以区分每个基因的主键（需要是唯一的）
#参数三k1：去除低质量细胞 total reads < k
#参数四k2：去除低质量细胞 total reads > k
#参数五p：去除低表达基因 total count < p
#参数六q：选择方差为前q的基因进行聚类








