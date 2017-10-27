rm(list=ls())
library(reshape2)
library(ggplot2)

df <- read.table("C:/Users/37908/Desktop/qt.txt",header = TRUE)
head(df)

df_melt <- melt(df, measure.vars = c("Healthy","Diseases"))

#setwd("C:/Users/37908/Desktop/")
#jpeg(file="output.jpg",width=500,height=500)
ggplot(data=df_melt, aes(x=M, y=value, fill=variable))+
  geom_bar(position="identity",stat = "identity", alpha = 0.7)+
  scale_fill_manual(values=c("#18B3B7", "#F35E5A")) +
  scale_y_continuous(name = "Abundance",
                     breaks = seq(0, 30, 5),
                     limits=c(0, 30)) +
  scale_x_discrete(name = "") +
  ggtitle("fafsdafsad")+
  theme_bw()+
  theme(plot.title=element_text(size = 20, hjust = 0.5),
        legend.title=element_blank())+
  #legend.text  = element_text(color="134", size=16, face="bold"))+
  coord_flip()
#dev.off()
