library(ggplot2)
rm(list=ls())
#sample <- read.table(file.choose(),head=T,sep="\t")
#fpkm <- read.table("/BIGDATA/zju_hwouyang_1/ZZL/fpkm_matrix.txt",head=T,sep="\t")
fpkm <- read.table(file.choose(),head=T,sep="\t")

sam_mat <- as.matrix(fpkm)
colnames(sam_mat[,1])

data <- rep(0,nrow(sam_mat)*3)

tot <- matrix(data, nrow = nrow(sam_mat) , ncol = 3, byrow = F)
row.names(tot) <- row.names(sam_mat)
colnames(tot) <- c(0,16,17)

for(i in 1 : nrow(sam_mat))
{
  for(j in 1:ncol(sam_mat))
  {
    if(j<=96)
    {
      tot[i,1] = tot[i,1] + sam_mat[i,j]
    }
    else if(j<=192)
    {
      tot[i,2] = tot[i,2] + sam_mat[i,j]
    }
    else
    {
      tot[i,3] = tot[i,3] + sam_mat[i,j]
    }
  }
}

setwd("C:/Users/37908/Desktop/")
jpeg(file="output.jpg",width=1400,height=500)
par(mfrow=c(3,5))
tot_mat <- data.frame(tot)
for(i in 1:15)
plot(x=c(0,16,19),y=tot[i,],  type="l",xlab="time",ylab="amount",main=row.names(tot_mat)[i],cex.main=2)
dev.off()

data <- rep(0,nrow(sam_mat) * ncol(sam_mat)*3)
more <- matrix(data, nrow = nrow(sam_mat) * ncol(sam_mat) , ncol = 3, byrow = F)

k = 1;
for(i in 1 : nrow(sam_mat))
{
  for(j in 1:ncol(sam_mat))
  {
    more[k,1] <- rownames(sam_mat)[i];
    more[k,2] <- sam_mat[i,j];
    more[k,3] <- colnames(sam_mat)[j]
    #if(grep('E16',colnames(sam_mat)[j]) > 0)
    #{
    #  more[k,3] <- "16";
    #}
    #else if(grep('E17',colnames(sam_mat)[j]) > 0)
    #{
    #  more[k,3] <- "19";
    #}
    #else
    #{
    #  more[k,3] <- "0";
    #}
      k = k+1;
  }
}
more <- data.frame(more)

ggplot(data= NULL, aes(x = more$X2, y = more$X3)) +  #开始绘图
  geom_point(aes(color = more$X1));
  #添加点
#  annotate("text",x =13 , y = 20,parse = T,
#           label = "x[1] == x[2]") #添加注释

  
  
  
  
  
  
  
  


#data <- rep(" ",278*3)
#more <- matrix(data, nrow = 278 , ncol = 3, byrow = F)

#par(mfrow=c(3,5))

#for(i in 1:nrow(sam_mat))

  i=15;
  k = 1;
  data <- rep(" ",278*3)
  more <- matrix(data, nrow = 278 , ncol = 3, byrow = F)
  
  for(j in 1:ncol(sam_mat))
  {
    #more[k,1] <- rownames(sam_mat)[i];
    more[k,2] <- sam_mat[i,j];
    #more[k,3] <- colnames(sam_mat)[j]
    if(!identical(grep('E16',colnames(sam_mat)[j]), integer(0)))
    {
      more[k,3] <- 16;
    }
    if(!identical(grep('E17',colnames(sam_mat)[j]), integer(0)))
    {
      more[k,3] <- 19;
    }
    if(!identical(grep('P0',colnames(sam_mat)[j]), integer(0)))
    {
      more[k,3] <- 0;
    }
    k = k+1;
  }
  more <- data.frame(more)
  setwd("C:/Users/37908/Desktop/")
  jpeg(file=paste(i,".jpg"),width=1400,height=500)
  ggplot(data= NULL, aes(x = more$X3, y = more$X2)) +  #开始绘图
    geom_point(aes(color = more$X1),stat = "sum") + theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank())
  dev.off();
  

  
