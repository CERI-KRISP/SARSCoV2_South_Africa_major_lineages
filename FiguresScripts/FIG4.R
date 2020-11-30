library(ggplot2)
library(ape)
library(repr)
library("readxl")
library('gridExtra')
#library(tidyverse)
library(dplyr)
library(hrbrthemes)
library(ggpubr)
library(cowplot)
library(ggthemes)
library(viridis)
library(ggrepel)
library("ggsci")
library(ggalt)
library("Hmisc")
library(ggtree)
library(tidyverse)
library(treeio)


### Fig 4A

tree <- read.beast("figures_original/BEAST_results/trees/B.1.106_corrected.tree")
p <- ggtree(tree, mrsd="2020-07-21", as.Date=TRUE,color='gray30',size=1) + theme_tree2()


p <- p+ ggstyleTILE +
  #geom_tiplab(aes(
  #subset=(grepl('SouthAfrica',label,fixed=TRUE)==TRUE)),linesize=2,color='black',align=F)+
  #scale_fill_manual(name='Location',values=c("cadetblue1", "cadetblue2",'cadetblue3','cadetblue4','steelblue1','steelblue2','steelblue3','steelblue4','lightcyan1','lightcyan2','lightcyan3'))+
  scale_fill_simpsons(name='Location')+
  #geom_tippoint(aes(
  #  subset=(grepl('KZN',label,fixed=TRUE)==TRUE)),size=4, align=F, color='black', fill='cadetblue',shape=21)+
  geom_tippoint(aes(
    subset=(grepl('ETH',label,fixed=TRUE)==TRUE),fill='ETH'),size=3, align=F, color='black' ,shape=21)+
  geom_tippoint(aes(
    subset=(grepl('DC25',label,fixed=TRUE)==TRUE), fill='DC25'),size=3, align=F, color='black',shape=21)+
  geom_tippoint(aes(
    subset=(grepl('DC43',label,fixed=TRUE)==TRUE), fill='DC43'),size=3, align=F, color='black',shape=21)+
  geom_tippoint(aes(
    subset=(grepl('DC29',label,fixed=TRUE)==TRUE),fill='DC29'),size=3, align=F, color='black',shape=21)+
  geom_tippoint(aes(
    subset=(grepl('DC28',label,fixed=TRUE)==TRUE),fill='DC28'),size=3, align=F, color='black',shape=21)+
  geom_tippoint(aes(
    subset=(grepl('DC21',label,fixed=TRUE)==TRUE),fill='DC21'),size=3, align=F, color='black',shape=21)+
  geom_tippoint(aes(
    subset=(grepl('DC22',label,fixed=TRUE)==TRUE), fill='DC22'),size=3, align=F, color='black',shape=21)+
  geom_tippoint(aes(
    subset=(grepl('DC27',label,fixed=TRUE)==TRUE), fill='DC27'),size=3, align=F, color='black',shape=21)+
  geom_tippoint(aes(
    subset=(grepl('DC24',label,fixed=TRUE)==TRUE), fill='DC24'),size=3, align=F, color='black',shape=21)+
  geom_tippoint(aes(
    subset=(grepl('DC23',label,fixed=TRUE)==TRUE), fill='DC23'),size=3, align=F, color='black',shape=21)+
  geom_tippoint(aes(
    subset=(grepl('DC26',label,fixed=TRUE)==TRUE), fill='DC26'),size=3, align=F, color='black',shape=21)+
  geom_tippoint(aes(
    subset=(grepl('CH3',label,fixed=TRUE)==TRUE), shape='CH3'),size=3, align=F, fill='olivedrab4',color='black')+
  geom_tippoint(aes(
    subset=(grepl('CH1',label,fixed=TRUE)==TRUE),shape='CH1'),size=3, align=F, fill='olivedrab4',color='black')+
  scale_x_date(date_labels = "%b",date_breaks = "1 month")+
  scale_shape_manual(name = "Outbreaks",
                     values = c(25, 22)) +
  theme(axis.text=element_text(size=10))+
  expand_limits(y = 70)+
  ggtitle("B.1.106") +xlab('month')
p


#### Fig 4B

data2<-read_excel('data tables/SA_final_dataset_n=1365_v2_cluster_annotated_final.xlsx')

data2$cluster<-factor(data2$cluster,levels = c("remaining","cluster1","cluster3","cluster4","cluster5"))

data2$days<-as.Date(cut(data2$date,
                        breaks = "day",
                        start.on.monday = FALSE))

data2$date<-as.Date(cut(data2$date,
                        breaks = "week",
                        start.on.monday = FALSE))

data2$first_case<-as.Date(cut(data2$first_case,
                              breaks = "day",
                              start.on.monday = FALSE))

data2$first_cluster<-as.Date(cut(data2$first_cluster,
                                 breaks = "day",
                                 start.on.monday = FALSE))


data2$cluster_first_case<-as.Date(cut(data2$cluster_first_case,
                                      breaks = "day",
                                      start.on.monday = FALSE))

data2$cluster_last_case<-as.Date(cut(data2$cluster_last_case,
                                     breaks = "day",
                                     start.on.monday = FALSE))



P_106_KZN_dates<-data2 %>%
  mutate(Collection_date=as.POSIXct(date)) %>% #convert date to date
  ggplot(data=subset(data2, division=='KwaZulu-Natal'), mapping = aes(x = date, fill=(Pangolin=='B.1.106')))+
  ggstyleTILE +
  geom_bar(position='fill',width=5,color='black')+
  xlab("Sampling Date")+
  scale_x_date(date_labels = "%b",date_breaks = "1 month")+
  scale_fill_manual(name='Lineages',values=c('white','black'),labels = c("others", "B.1.106"))+
  xlab('month')+
  ylab('proportion of B.1.106 in KZN')


pdf("output/fig4.1.pdf", w=4.5, h=7, bg='white')
print(cowplot::plot_grid(p,P_106_KZN_dates,ncol=1,labels=c("A","B"),align='vh'))
a<-dev.off()
