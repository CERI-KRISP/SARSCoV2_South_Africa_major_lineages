

source("functions.plotting.R")
require(tidyverse)
library("readxl")
library("lubridate")
require("ggalt")
library("wesanderson")


#################################################
###### PANEL 2B

data2<-read_excel('data_original/SA_final_dataset_n=1365_v2_cluster_annotated_final.xlsx')
data2$cluster<-factor(data2$cluster,levels = c("remaining","cluster1","cluster3","cluster4","cluster5"))
data2$days<-as.Date(cut(data2$date,breaks = "day",start.on.monday = FALSE))
data2$date<-as.Date(cut(data2$date,breaks = "week",start.on.monday = FALSE))
data2$first_case<-as.Date(cut(data2$first_case,breaks = "day",start.on.monday = FALSE))
data2$first_cluster<-as.Date(cut(data2$first_cluster,breaks = "day",start.on.monday = FALSE))
data2$cluster_first_case<-as.Date(cut(data2$cluster_first_case, breaks = "day", start.on.monday = FALSE))
data2$cluster_last_case<-as.Date(cut(data2$cluster_last_case,breaks = "day",start.on.monday = FALSE))
data2<- data2 %>% filter(lineage_cluster!="others")
data2<- subset(data2, !is.na(lineage_cluster))
##this will force the lineages to come first on the y axis =)
data2$lineage_cluster[which(data2$lineage_cluster==unique(data2$lineage_cluster)[1])]<- paste(" ",unique(data2$lineage_cluster)[1]," ")
data2$lineage_cluster[which(data2$lineage_cluster==unique(data2$lineage_cluster)[2])]<- paste(" ",unique(data2$lineage_cluster)[2]," ")
data2$lineage_cluster[which(data2$lineage_cluster==unique(data2$lineage_cluster)[3])]<- paste(" ",unique(data2$lineage_cluster)[3]," ")


lineages_color<- wes_palette("Zissou1")[c(1,3,5)]
panelB1<-ggplot(data=data2) + ggstyleTILE +
  geom_segment(aes(x=min(data2$date), y = division, xend=max(data2$date), yend=division, group=division), colour="grey88", size=10) +
  geom_point(aes(x=days, y=division, fill=lineage_cluster), position = position_jitter(width=0.3, height=0.3), shape=21, col='grey33', size=2.1, alpha=1)+
  geom_segment(aes(x = cluster_first_case, y = lineage_cluster, xend = cluster_last_case, yend = lineage_cluster, colour = lineage_cluster, group=lineage_cluster), size=1.2) +
  ylab('')+ xlab('month')+ ggtitle('South Africa lineages by division') +
  scale_color_manual("Lineage",values=lineages_color) +
  scale_fill_manual("Lineage",values=lineages_color) +
  theme(legend.position="none") +
  scale_y_discrete(breaks=sort(unique(data2$division)), labels=sort(unique(data2$division))) + ##hides the lineages
  geom_label(aes(x=cluster_first_case, y=lineage_cluster, label=lineage_cluster,  fill = lineage_cluster, group=lineage_cluster), color='white', fontface = "bold", nudge_y=0.3, nudge_x=4, size=3)


data3<- subset(subset(data2, division=='KwaZulu-Natal'), !is.na(lineage_cluster))

lineages_color<- wes_palette("Zissou1")[c(1,3,5)]
panelB2<-ggplot(data=data3) + ggstyleTILE +
  geom_segment(aes(x=min(data2$date), y = location, xend=max(data2$date), yend=location, group=location), colour="grey88", size=5) +
  geom_point(aes(x=days, y=location, fill=lineage_cluster), position = position_jitter(width=0.3, height=0.3), shape=21, col='grey33', size=2.1, alpha=1)+
  geom_segment(aes(x = cluster_first_case, y = lineage_cluster, xend = cluster_last_case, yend = lineage_cluster, colour = lineage_cluster, group=lineage_cluster), size=1.2) +
  ylab('')+ xlab('month')+ ggtitle('KwaZulu-Natal (SA) lineages by location') +
  scale_color_manual("Lineage",values=lineages_color) +
  scale_fill_manual("Lineage",values=lineages_color) +
  theme(legend.position="none") +
  scale_y_discrete(breaks=sort(unique(data2$location)), labels=sort(unique(data2$location))) + ##hides the lineages
  geom_label(aes(x=cluster_first_case, y=lineage_cluster, label=lineage_cluster,  fill = lineage_cluster, group=lineage_cluster), color='white', fontface = "bold", nudge_y=0.5, nudge_x=4, size=3)

  top<- cowplot::plot_grid(gempty, panelB1, panelB2, labels=c("A","B","C"), ncol=3, rel_widths=c(1,1,1))

##############################################################
##############################################################
### PANEL C

data2<-read_excel('data_original/SA_final_dataset_n=1365_v2_cluster_annotated_final.xlsx')
data2$cluster<-factor(data2$cluster,levels = c("remaining","cluster1","cluster3","cluster4","cluster5"))
data2$days<-as.Date(cut(data2$date,breaks = "day",start.on.monday = FALSE))
data2$date<-as.Date(cut(data2$date,breaks = "week",start.on.monday = FALSE))
data2$first_case<-as.Date(cut(data2$first_case,breaks = "day",start.on.monday = FALSE))
data2$first_cluster<-as.Date(cut(data2$first_cluster,breaks = "day",start.on.monday = FALSE))
data2$cluster_first_case<-as.Date(cut(data2$cluster_first_case, breaks = "day", start.on.monday = FALSE))
data2$cluster_last_case<-as.Date(cut(data2$cluster_last_case,breaks = "day",start.on.monday = FALSE))
data4<- data2 %>% mutate(Collection_date=as.POSIXct(date))
data4<- subset(data4, !is.na(lineage_cluster))
data4<- subset(data4, !is.na(strain))

panelC<- ggplot(data=data4, mapping = aes(x = date, fill=lineage_cluster))+ ggstyleTILE +
        geom_bar(position='fill',width=5,color='grey33')+
        scale_x_date(date_labels = "%b",date_breaks = "month")+
        scale_fill_manual(name="Lineage",values=c(lineages_color,"grey"))+
        xlab('month')+ ylab('proportion') + ggtitle('Proportion of genomes by lineage') +
        theme(legend.position="none")

############################################################
############################################################
### PANEL D

data2<-read_excel('data_original/SA_final_dataset_n=1365_v2_cluster_annotated_final.xlsx')
data2$cluster<-factor(data2$cluster,levels = c("remaining","cluster1","cluster3","cluster4","cluster5"))
data2$days<-as.Date(cut(data2$date,breaks = "day",start.on.monday = FALSE))
data2$date<-as.Date(cut(data2$date,breaks = "week",start.on.monday = FALSE))
data2$first_case<-as.Date(cut(data2$first_case,breaks = "day",start.on.monday = FALSE))
data2$first_cluster<-as.Date(cut(data2$first_cluster,breaks = "day",start.on.monday = FALSE))
data2$cluster_first_case<-as.Date(cut(data2$cluster_first_case, breaks = "day", start.on.monday = FALSE))
data2$cluster_last_case<-as.Date(cut(data2$cluster_last_case,breaks = "day",start.on.monday = FALSE))
data5<- subset(subset(data2, !is.na(lineage_cluster)))
data5<- subset(data5, !is.na(strain))
data5<- subset(data5, division!="South Africa")


panelD<-ggplot(data=data5)+ ggstyleTILE +
  geom_bar(aes(x = division,fill=lineage_cluster), width=0.8,color='grey33')+
  scale_fill_manual("Lineage", values = c(lineages_color,"grey"))+
  xlab('')+ ylab('number of genomes')+ ggtitle("South Africa genomes by division") +
  coord_flip() + theme(legend.position="none")


############################################################
############################################################
### PANEL E

data6<- subset(data5, division=='KwaZulu-Natal')
data6<- subset(data6, location!='ND')

panelE<-ggplot(data=data6)+ ggstyleTILE +
  geom_bar(aes(x = location,fill=lineage_cluster), width=0.8,color='grey33')+
  scale_fill_manual("Lineage", values = c(lineages_color,"grey"))+
  xlab('')+ ylab('number of genomes')+ ggtitle("KwaZulu-Natal (SA) genomes by location") +
  coord_flip() + theme(legend.position="none")

bottom<- cowplot::plot_grid(panelC, panelD, panelE, labels=c("D","E","F"), ncol=3, rel_widths=c(1.7,2,1.9))


######################################################
######################################################
######### ALL PANELS

gempty<- ggplot()


all<- cowplot::plot_grid(top, bottom, nrow=2, rel_heights=c(1,1))

pdf("output/fig2.pdf",w=14,h=8,bg='white')
print(all)
a<-dev.off()
