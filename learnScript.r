library(ggplot2);
ggplot(data=NULL,aes(x=x,y=y)) + geom_point(color="darkred") + annotate("text",label="Pearson = 0.9645")



set.seed(1410)
dsmall <- diamonds[sample(nrow(diamonds),100),]
qplot(carat,price,data=dsmall,color=color,xlab="Carat(?)",ylab="Price($)",main="Price and Carat with color")

#散点图
qplot(carat,price,data=dsmall,alpha = I(1/4),geom=c("point","smooth"),span=0.8)
#alpha表示每个点的透明度，这里四个重合的点构成一个不透明的点
#span决定拟合曲线的平滑度，1为完全平滑，0为完全曲折

qplot(carat,price,data=diamonds,geom=c("point","smooth"))

#根据颜色分类图
qplot(color,price/carat,data=diamonds,geom="jitter",alpha=I(1/5))
qplot(color,price/carat,data=diamonds,geom="jitter",alpha=I(1/50))
qplot(color,price/carat,data=diamonds,geom="jitter",alpha=I(1/50),color=color,ylab=expression(frac(price,carat)))

#直方图,binwidth控制跨度
qplot(carat,data=diamonds,geom="histogram",binwidth=0.1,fill=color)
#密度曲线图,adjust控制平滑度
qplot(carat,data=diamonds,geom="density",adjust=1,color=color)

#条形图
#按照颜色的数目进行分类
qplot(color,data = diamonds,geom="bar")
#用carat作为权重相加得到高度
qplot(color,data = diamonds,geom="bar",weight=carat)

#时间序列中的线条图
qplot(date,unemploy / pop,data=economics,geom="line")

#时间序列中的路径图
year <- function(x) as.POSIXlt(x)$year + 1900
qplot(unemploy / pop, uempmed, data = economics,geom=c("point","path"))
qplot(unemploy / pop, uempmed, data = economics,geom="path",colour=year(date))


#facets=color ~ . 表示生成color行，1列的窗格
qplot(carat,data=diamonds,facets=color ~ .,geom="histogram",binwidth=0.1,xlim=c(0,3))
#将密度而不是频数映射到y轴，使得比较不同组的分布时不受样本量的影响
qplot(carat,..density..,data=diamonds,facets=color ~ .,geom="histogram",binwidth=0.1,xlim=c(0,3))


qplot(displ,hwy,data=mpg,colour=factor(cyl))
