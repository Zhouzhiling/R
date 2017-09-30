#读入txt文件（以制表符为分割）
library(ggplot2)
rm(list=ls())
#sample <- read.table(file.choose(),head=T,sep="\t")
#fpkm <- read.table("/BIGDATA/zju_hwouyang_1/ZZL/fpkm_matrix.txt",head=T,sep="\t")
fpkm <- read.table(file.choose(),head=T,sep="\t")

row.names(fpkm)<-fpkm[,1]
fpkmNum <- fpkm[,-1]
fpkmNum <- t(as.matrix(fpkmNum))
is.numeric(fpkmNum)


#0930new added
#选择MAD最大的前1/4基因
addTage <- fpkmNum[1,];
i=1;
while(i<length(addTage))
{
  addTage[i] <- mad(fpkmNum[,i]);
  i = i+1;
}
fpkmWithTag <- t(data.frame(addTage,t(fpkmNum)))
fpkmWithTagg <-fpkmWithTag[,order(fpkmWithTag[1,],decreasing = TRUE)]

chosenGene = ceiling(ncol(fpkmWithTag) / 4);
fpkmNumS <- fpkmWithTagg[-1,]
fpkmNumS <- fpkmNumS[,1:chosenGene]

#自己电脑跑先少一点

DDT <- t(fpkmNumS) %*% fpkmNumS
CovM <- cov(DDT)

#计算特征向量
Aeigen=eigen(CovM)
tezhenzhi <- Aeigen$values
#Aeigen$vectors
max(tezhenzhi)
maxplace = -1;
i = 1;
while(maxplace == -1)
{
  if(tezhenzhi[i] == max(tezhenzhi))
  {
    maxplace = i;
  }
  i <- i + 1;
}
maxplace
wvector <- Aeigen$vectors[,maxplace]
pResult <-  array(0,nrow(fpkmNumS))

for(j in 1:nrow(fpkmNumS))
{
  tmp <- 0;
  for(i in 1:ncol(fpkmNumS))
  {
    tmp = tmp + wvector[i] * fpkmNumS[j,i];
  }
  pResult[j] = tmp;
}
plot(pResult)

#准备名字
cellName <- rownames(fpkmNumS)
cellName <- substr(cellName,8,10)
cellName

#然后需要新建一个列表来存需要画出来的结果
result <- matrix(pResult,3,ncol = length(pResult)*2,byrow=TRUE)

result[1,1:length(pResult)] <- cellName
result[1,(length(pResult)+1):ncol(result)] <- cellName
result[2,1:length(pResult)] <- pResult
result[2,(length(pResult)+1):ncol(result)] <- pResult
result[c(2,3),] <- as.numeric(result[c(2,3),])

average0 = 0;
average16 = 0;
average17 = 0;
count0 = count16 = count17 = 0;
for(i in 1:length(pResult))
{
  if(result[1,i]=="E16")
  {
    result[3,i] = 1;
    count16 = count16 + 1;
    average16 = average16 + as.numeric(result[2,i]);
  }
  else
  {
    if(result[1,i]=="E17")
    {
      result[3,i] = 0;
      count17 = count17 + 1;
      average17 = average17 +  as.numeric(result[2,i]);
    }
    else
    {
      result[1,i]= "P0";
      result[3,i] = 2;
      count0 = count0 + 1;
      average0 = average0 +  as.numeric(result[2,i]);
    }
  }
}

average0 = average0/count0;
average16 = average16/count16;
average17 = average17/count17;

for(i in (length(pResult)+1):ncol(result))
{
  if(result[1,i]=="P0_")
    result[1,i]= "P0";
  result[3,i] = -0.5;
}

avg <- matrix(c("P0",average0,-1,"E16",average16,-1,"E17",average17,-1),nrow = 3,ncol = 3);

toPlot <- as.matrix(data.frame(avg,result))

setwd("C:/Users/37908/Desktop/")
#jpeg(file="/BIGDATA/zju_hwouyang_1/ZZL/myplot.jpeg",width=1400,height=500)
jpeg(file="myplot.jpeg",width=800,height=500)
qplot(toPlot[2,],toPlot[3,],color=toPlot[1,],size=I(4),alpha=I(1/2),xlab = "Pseudotime",ylab = "")


dev.off()

#qplot(result[2,],result[3,],color=result[1,])
#qplot(result[2,],0,color=result[1,])

#Aeigen$vectors
