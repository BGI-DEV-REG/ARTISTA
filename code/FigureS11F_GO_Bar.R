library(ggplot2)
library(viridis)

setwd('../data/FigureS11_data/FigureS11F')
s1 = read.csv('Selected_GO_brown_0517.txt',sep='\t',header=T)
s1 = s1[which(s1$Count >= 5),]
s1$LogP_1 = -log10(s1$p.adjust)
colnames(s1)
s1 = s1[order(s1$Count),] 
s1$Description = factor(s1$Description,levels = s1$Description)

p = ggplot(data=s1, aes(x=Description,y=Count, fill=LogP_1)) + geom_bar(stat="identity", width=0.9) + 
  coord_flip()  + ylab("")+theme(axis.text.y=element_text(color="black", size=12)) +
  scale_fill_viridis(option = "D", direction = -1,limits = c(3,7))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+ theme_bw() + theme(panel.grid=element_blank())+
  scale_y_continuous(limits=c(0,10), breaks=seq(0,10,2))

ggsave(plot=p,file='GO_ brown_selected.pdf')

####
s1 = read.csv('Selected_KEGG_brown_0517.txt',sep='\t',header=T)
s1 = s1[which(s1$Count >= 2),]
s1$LogP_1 = -log10(s1$p.adjust)
colnames(s1)
s1 = s1[order(s1$Count),] 
s1$Description = factor(s1$Description,levels = s1$Description)

ggplot(data=s1, aes(x=Description,y=Count, fill=LogP_1)) + geom_bar(stat="identity", width=0.9) + 
  coord_flip()  + ylab("")+theme(axis.text.y=element_text(color="black", size=12)) +scale_fill_viridis(option = "D", direction = -1,limits = c(0, 10))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+ theme_bw() + theme(panel.grid=element_blank())

ggsave('KEGG_brown_slected.pdf')
