用于处理的代码是Cluster.r，其中主要函数是fun.Calculate <- function(TextNum,KeyWord,CellThresholdUpper,CellThresholdLower,GeneThreshold,ClusterNum)

#参数一n：数据的前n列是文本信息，处理时需要删除
#参数二m：可以区分每个基因的主键（需要是唯一的）
#参数三k1：去除低质量细胞 total reads < k
#参数四k2：去除低质量细胞 total reads > k
#参数五p：去除低表达基因 total count < p
#参数六q：选择方差为前q的基因进行聚类


100，101，111、129、157行的路径需要在运行前修改

聚类的结果输出在文件ClusterResult.csv里

DESeq的处理结果在文件DESeqResult.csv里

还存在的问题是：
不知道SC3聚类的结果是直接输出吗，如果有格式要求的话需要在111行那里改
