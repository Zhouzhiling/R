# 目的
在Pseudotime下，根据基因量随时间变化的趋势给其分类。

# 源文件
fpkm_matrix.txt

#处理过程

> * 按照pseudotime升序给基因排序
> * 删除pseudotime异常高和低的若干个时间点
> * 以上下三分位数把pseudotime分成三段，把每一个基因在三段时间内的平均量计算出来
> * 删除三段时间的平均量均为0的基因,得到剩余的19157中基因记录
> * 预设变化趋势为九种，分别在三个时间点的变化为：降降，降升，降稳，升降，升升，升稳，升降，稳升，稳降，稳稳
>  * 其中对稳的定义为，两个时间段的平均值的变化量不足10%。
>  * 对升的定义为，后一个时间段的平均值的相比前一个升高了超过10%。
>  *对降的定义为，后一个时间段的平均值的相比前一个降低了超过10%。
> * 计算并把基因名放入对应的变化趋势下，导出得到分类结果。

# R代码
```R
rm(list=ls())
src <- read.table(file.choose(),head=T,sep="\t")
src
src <- src[,order(src[1,])]

#setwd("C:/Users/37908/Desktop/")
#jpeg(file="1027.jpg",width=20000,height=10000)
#plot(src[1,])
#dev.off()

src_trim <-as.matrix(src[,-rep(1:26)])
low_quant <- ceiling(ncol(src_trim) / 3)
high_quant <- ceiling(ncol(src_trim) / 3) * 2
low_quant_num <- src_trim[1,low_quant]
high_quant_num <- src_trim[1,high_quant]

persutime <- src_trim[1,]
src_trim <- src_trim[-1,]
gene_name <- row.names(src_trim)
flag <- matrix(0,nrow=nrow(src_trim),ncol=3)

count1 <- 0
count2 <- 0
count3 <- 0

#归类求各个时期个个基因表达量的和
for(j in 1:length(persutime))
{
  if(persutime[j] < low_quant_num)
    for(i in 1:nrow(src_trim))
    {
      if(persutime[j] < low_quant_num)
      {
        flag[i,1] = flag[i,1]+src_trim[i,j]
        count1 <- count1 + 1
      }
    }
  
  else if(persutime[j] < high_quant_num)
  {
    for(i in 1:nrow(src_trim))
    {
      flag[i,2] = flag[i,2]+src_trim[i,j]
      count2 <- count2 + 1
    }
  }
  
  else
    for(i in 1:nrow(src_trim))
    {
      {
        flag[i,3] = flag[i,3]+src_trim[i,j]
        count3 <- count3 + 1
      }
    }
}

row.names(flag) <- gene_name


#把全为0的删掉先
tmp_flag <- matrix(0,nrow=nrow(flag),ncol=3)
tmp_gene_frame <- matrix("",nrow=nrow(flag),ncol=1)
count_no_zero <- 0
#is.numeric(flag)
for(k in 1:nrow(flag))
{
  if((flag[k,1]!=0)||(flag[k,2]!=0)||(flag[k,3]!=0))
  {
    tmp_flag[count_no_zero,] <- flag[k,]
    tmp_gene_frame[count_no_zero,1] <- gene_name[k]
    count_no_zero <- count_no_zero + 1
  }
}


#基因三阶段数量
flag_no_zero <- tmp_flag[rep(1:19157),]
#基因名称
tmp_gene_frame <- tmp_gene_frame[rep(1:19157),1]
row.names(flag_no_zero) <- tmp_gene_frame


#得到和之后求平均值
flag_avg <- matrix(0,nrow=nrow(flag_no_zero),ncol=3)
for(i in 1:nrow(flag_no_zero))
{
  {
    flag_avg[i,1] = flag_no_zero[i,1] * 1000000 / count1
    flag_avg[i,2] = flag_no_zero[i,2] * 1000000 / count2
    flag_avg[i,3] = flag_no_zero[i,3] * 1000000 / count3
  }
}
row.names(flag_avg) <- tmp_gene_frame

#取样看个例子w
setwd("C:/Users/37908/Desktop/")
jpeg(file="output.jpg",width=2500,height=500)
par(mfrow=c(3,5))
flag_avg_mat <- data.frame(flag_avg)
for(i in 1:15)
  plot(x=c(0,4,19),y=flag_avg_mat[i,],  type="l",xlab="time",ylab="amount",main=row.names(flag_avg_mat)[i],cex.main=2)
dev.off()

#导出成表格，备用查询具体细节
write.csv(flag_avg,file="C:/Users/37908/Desktop/classify.csv",row.names=TRUE)

#?write.csv

class <- matrix("",nrow=nrow(flag_avg) , ncol=9 )
# d means decrease
# k means keep stable
# i means increase
dd_num <- 0
di_num <- 0
dk_num <- 0
id_num <- 0
ii_num <- 0
ik_num <- 0
kd_num <- 0
ki_num <- 0
kk_num <- 0


#确认是否是数字类型
#is.numeric(flag_avg)

error = 0.1

for(i in 1:nrow(flag_avg))
# for(i in 1:8) 
{
  # 从第一到第二阶段降
  if(flag_avg[i,2] < flag_avg[i,1]*(1-error))
  {
    
    if(flag_avg[i,2] > flag_avg[i,3]*(1+error) )
    {
      dd_num <- dd_num + 1 
      #第一类,降降
      class[dd_num,1] <- tmp_gene_frame[i]
    }
    
    if(flag_avg[i,3] > flag_avg[i,2]*(1+error) )
    {
      di_num <- di_num + 1 
      #第二类，降升
      class[di_num,2] <- tmp_gene_frame[i]
    }
    
    if(flag_avg[i,3]<= flag_avg[i,2]*(1+error) && flag_avg[i,3]>=flag_avg[i,2]*(1-error))
    {
      dk_num <- dk_num + 1
      #第三类，降平
      class[dk_num,3] <- tmp_gene_frame[i]
    }
  }
  
  
  # 从第一到第二阶段升
  if(flag_avg[i,2] > flag_avg[i,1]*(1+error))
  {
    if(flag_avg[i,3] < flag_avg[i,2]*(1-error) )
    {
      id_num <- id_num + 1 
      #第四类,升降
      class[id_num,4] <- tmp_gene_frame[i]
    }
    
    if(flag_avg[i,3] > flag_avg[i,2]*(1+error) )
    {
      ii_num <- ii_num + 1 
      #第五类，升升
      class[ii_num,5] <- tmp_gene_frame[i]
    }
    
    if(flag_avg[i,3]<= flag_avg[i,2]*(1+error) && flag_avg[i,3]>=flag_avg[i,2]*(1-error))
    {
      ik_num <- ik_num + 1
      #第六类，升平
      class[ik_num,6] <- tmp_gene_frame[i]
    }
  }
  
  #从第一到第二阶段平
  if(flag_avg[i,2] <= flag_avg[i,1]*(1+error) && flag_avg[i,2] >= flag_avg[i,1]*(1-error))
  {
    if(flag_avg[i,3] < flag_avg[i,2]*(1-error) )
    {
      kd_num <- kd_num + 1 
      #第七类,平降
      class[kd_num,7] <- tmp_gene_frame[i]
    }
    if(flag_avg[i,3] > flag_avg[i,2]*(1+error) )
    {
      ki_num <- ki_num + 1 
      #第八类，平升
      class[ki_num,8] <- tmp_gene_frame[i]
    }
    if(flag_avg[i,3]<= flag_avg[i,2]*(1+error) && flag_avg[i,3]>=flag_avg[i,2]*(1-error))
    {
      kk_num <- kk_num + 1
      #第九类，平平
      class[kk_num,9] <- tmp_gene_frame[i]
    }
  }
  }

colnames(class) <- c("↘↘","↘↗","↘→","↗↘","↗↗","↗→","→↘","→↗","→→ ")
class <- class[rep(1:5432),]
write.csv(class,file="C:/Users/37908/Desktop/result.csv",row.names=TRUE)
```




# 输出结果
### *detail_of_each_gene.csv*
每行纪录了每一个基因在三个时间段中的表达量，第一列为基因名，后三列是时间增序排列下的表达量
### *classification_results.csv*
归类结果，每一列为一种变化趋势，同一列下是相同变化趋势的所有基因名称
