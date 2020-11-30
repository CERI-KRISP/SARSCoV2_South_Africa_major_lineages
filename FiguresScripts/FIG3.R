
source("functions.plotting.R")
require(tidyverse)
library("readxl")
library("lubridate")
require("ggalt")
library("wesanderson")

#####################################################
#####################################################
### PANEL 3A


data2<-read_excel('data_original/SA_final_dataset_n=1365_v2_cluster_annotated_final.xlsx')
data2$cluster<-factor(data2$cluster,levels = c("remaining","cluster1","cluster3","cluster4","cluster5"))
data2$days<-as.Date(cut(data2$date, breaks = "day", start.on.monday = FALSE))
data2$date<-as.Date(cut(data2$date, breaks = "week", start.on.monday = FALSE))
data2$first_case<-as.Date(cut(data2$first_case,   breaks = "day",start.on.monday = FALSE))
data2$first_cluster<-as.Date(cut(data2$first_cluster, breaks = "day", start.on.monday = FALSE))
data2$cluster_first_case<-as.Date(cut(data2$cluster_first_case,breaks = "day",start.on.monday = FALSE))
data2$cluster_last_case<-as.Date(cut(data2$cluster_last_case,breaks = "day",start.on.monday = FALSE))
data2<- subset(data2, !is.na(lineage_cluster))

lineages_color<- c("darkolivegreen1",wes_palette("Zissou1")[c(1,3)],wes_palette("Zissou1")[c(5)])
lineages_color<- c(lineages_color,"grey")

data2$lineage_cluster<- factor(data2$lineage_cluster, levels=c("B.1.106", "B.1.1.54", "B.1.1.56", "C.1", "others"))

panel3A <- ggplot(data=data2, aes(lineage_cluster, totalMutations, fill=lineage_cluster))+ ggstyleTILE +
  geom_violin()+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),geom = "pointrange", color = "black" )+
  scale_fill_manual("Lineage",values=lineages_color) +
  xlab('')+ ylab('number of mutations')+ggtitle("Lineages and mutations")+
  theme(legend.position="none", axis.text.x = element_text(angle=45, hjust=1))



#######################################################
########################################################
### PANEL B

codeMutationTags<- function(count, thre){
    over<- names(count[which(as.numeric(count)>=thre)])
    over<- strsplit(over, "(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])", perl=TRUE)
    over<- matrix(unlist(over),ncol=3,byrow=T)
    positions<- over[,2]
    over<- cbind(over[,2], ": ", over[,1], " to ", over[,3])
    over<- apply(over, MARG=1, FUN=function(X) { paste(X, collapse='') })
    over<- data.frame(x_pos=as.numeric(positions), y_pos=thre/2, tag=as.character(over), stringsAsFactors=FALSE)
    return(over)
}

B1154 <- read.csv(file='data_original/mutation_table_B.1.1.54_V2.csv', sep=",", header=T)
B1156 <- read.csv(file='data_original/mutation_table_B.1.1.56_V2.csv', sep=",", header=T)
C1 <- read.csv(file='data_original/mutation_table_C.1_V2.csv', sep=",", header=T)
B1106 <- read.csv(file='data_original/mutation_table_B.1.106.csv', sep=",", header=T)

B1154_tags<- codeMutationTags(table(B1154$mutation), thre=288)
B1156_tags<- codeMutationTags(table(B1156$mutation), thre=94)
C1_tags<- codeMutationTags(table(C1$mutation), thre=135)
B1106_tags<- codeMutationTags(table(B1106$mutation), thre=61)

B1154table<-table(B1154$location)
B1154table<- data.frame(count=as.numeric(B1154table),location=as.numeric(names(B1154table)))
nx<- 400
ss<- 2.4
panelB1 <- ggplot(B1154table, aes(x = location, y=count))+ ggstyleTILE +
  geom_bar(fill="black", col="black",stat="identity") +
  xlab("genome site")+ ylab("number of genomes")+ ggtitle("B.1.1.54 mutational map") +
  scale_x_continuous(limits = c(0,29903),breaks = seq(0, 29903, by = 5000))+
  geom_hline(yintercept=288, linetype="dashed", color = "grey33") +
  ##sadly, this is the easier way
  annotate("text", fontface =2,x=B1154_tags$x_pos[1]-nx, y=B1154_tags$y_pos[1], label=B1154_tags$tag[1], angle=90, size=ss)+
  annotate("text", fontface =2,x=B1154_tags$x_pos[2]-nx, y=B1154_tags$y_pos[2], label=B1154_tags$tag[2], angle=90, size=ss)+
  annotate("text", fontface =2,x=B1154_tags$x_pos[3]+nx, y=B1154_tags$y_pos[3], label=B1154_tags$tag[3], angle=90, size=ss)+
  annotate("text", fontface =2,x=B1154_tags$x_pos[4]-nx, y=B1154_tags$y_pos[4], label=B1154_tags$tag[4], angle=90, size=ss,col='tomato3') +
  annotate("text", fontface =2,x=B1154_tags$x_pos[5]-nx, y=B1154_tags$y_pos[5], label=B1154_tags$tag[5], angle=90, size=ss)+
  annotate("text", fontface =2,x=B1154_tags$x_pos[6]-nx, y=B1154_tags$y_pos[6], label=B1154_tags$tag[6], angle=90, size=ss)+
  annotate("text", fontface =2,x=B1154_tags$x_pos[7]-nx*3, y=B1154_tags$y_pos[7], label=B1154_tags$tag[7], angle=90, size=ss) +
  annotate("text", fontface =2,x=B1154_tags$x_pos[8]-nx*5, y=B1154_tags$y_pos[8], label=B1154_tags$tag[8], angle=90, size=ss) +
  annotate("text", fontface =2,x=B1154_tags$x_pos[9]-nx, y=B1154_tags$y_pos[9], label=B1154_tags$tag[9], angle=90, size=ss,col='tomato3')

B1156table<-table(B1156$location)
B1156table<- data.frame(count=as.numeric(B1156table),location=as.numeric(names(B1156table)))
nx<- 400
ss<- 2.4
panelB2 <- ggplot(B1156table, aes(x = location, y=count))+ ggstyleTILE +
  geom_bar(fill="black", col="black",stat="identity") +
  xlab("genome site")+ ylab("number of genomes")+ ggtitle("B.1.1.56 mutational map") +
  scale_x_continuous(limits = c(0,29903),breaks = seq(0, 29903, by = 5000))+
  geom_hline(yintercept=94, linetype="dashed", color = "grey33") +
  ##sadly, this is the easier way
  annotate("text", fontface =2, x=B1156_tags$x_pos[1]+nx, y=B1156_tags$y_pos[1], label=B1156_tags$tag[1], angle=90, size=ss)+
  annotate("text", fontface =2, x=B1156_tags$x_pos[2]-nx, y=B1156_tags$y_pos[2], label=B1156_tags$tag[2], angle=90, size=ss) +
  annotate("text", fontface =2, x=B1156_tags$x_pos[3]-nx, y=B1156_tags$y_pos[3], label=B1156_tags$tag[3], angle=90, size=ss, col='tomato3') +
  annotate("text", fontface =2, x=B1156_tags$x_pos[4]+nx, y=B1156_tags$y_pos[4], label=B1156_tags$tag[4], angle=90, size=ss)+
  annotate("text", fontface =2, x=B1156_tags$x_pos[5]-nx, y=B1156_tags$y_pos[5], label=B1156_tags$tag[5], angle=90, size=ss) +
  annotate("text", fontface =2, x=B1156_tags$x_pos[6]-nx, y=B1156_tags$y_pos[6], label=B1156_tags$tag[6], angle=90, size=ss) +
  annotate("text", fontface =2, x=B1156_tags$x_pos[7]+nx, y=B1156_tags$y_pos[7], label=B1156_tags$tag[7], angle=90, size=ss)+
  annotate("text", fontface =2, x=B1156_tags$x_pos[8]+nx*2.5, y=B1156_tags$y_pos[8], label=B1156_tags$tag[8], angle=90, size=ss)

C1table<-table(C1$location)
C1table<- data.frame(count=as.numeric(C1table),location=as.numeric(names(C1table)))
nx<- 400
ss<- 2.4
panelB3 <- ggplot(C1table, aes(x = location,y=count))+ ggstyleTILE +
  geom_bar(fill="black", col="black",stat="identity") +
  xlab("genome site")+ ylab("number of genomes")+ ggtitle("C.1 mutational map") +
  scale_x_continuous(limits = c(0,29903),breaks = seq(0, 29903, by = 5000))+
  geom_hline(yintercept=135, linetype="dashed", color = "grey33") +
  ##sadly, this is the easier way
  annotate("text", fontface =2, x=C1_tags$x_pos[1]-nx, y=C1_tags$y_pos[1], label=C1_tags$tag[1], angle=90, size=ss) +
  annotate("text", fontface =2, x=C1_tags$x_pos[2]-nx, y=C1_tags$y_pos[2], label=C1_tags$tag[2], angle=90, size=ss, col='tomato3') +
  annotate("text", fontface =2, x=C1_tags$x_pos[3]-nx, y=C1_tags$y_pos[3], label=C1_tags$tag[3], angle=90, size=ss) +
  annotate("text", fontface =2, x=C1_tags$x_pos[4]-nx, y=C1_tags$y_pos[4], label=C1_tags$tag[4], angle=90, size=ss, col='tomato3') +
  annotate("text", fontface =2, x=C1_tags$x_pos[5]+nx, y=C1_tags$y_pos[5], label=C1_tags$tag[5], angle=90, size=ss, col='tomato3') +
  annotate("text", fontface =2, x=C1_tags$x_pos[6]+nx, y=C1_tags$y_pos[6], label=C1_tags$tag[6], angle=90, size=ss) +
  annotate("text", fontface =2, x=C1_tags$x_pos[7]-nx, y=C1_tags$y_pos[7], label=C1_tags$tag[7], angle=90, size=ss) +
  annotate("text", fontface =2, x=C1_tags$x_pos[8]-nx, y=C1_tags$y_pos[8], label=C1_tags$tag[8], angle=90, size=ss, col='tomato3') +
  annotate("text", fontface =2, x=C1_tags$x_pos[9]-nx, y=C1_tags$y_pos[9], label=C1_tags$tag[9], angle=90, size=ss, col='tomato3') +
  annotate("text", fontface =2, x=C1_tags$x_pos[10]-nx, y=C1_tags$y_pos[10], label=C1_tags$tag[10], angle=90, size=ss) +
  annotate("text", fontface =2, x=C1_tags$x_pos[11]-nx*3, y=C1_tags$y_pos[11], label=C1_tags$tag[11], angle=90, size=ss) +
  annotate("text", fontface =2, x=C1_tags$x_pos[12]-nx*5, y=C1_tags$y_pos[12], label=C1_tags$tag[12], angle=90, size=ss)

B1106table<-table(B1106$location)
B1106table<- data.frame(count=as.numeric(B1106table),location=as.numeric(names(B1106table)))
nx<- 400
ss<- 2.4
panelB4 <- ggplot(B1106table, aes(x = location,y=count))+ ggstyleTILE +
  geom_bar(fill="black", col="black",stat="identity") +
  xlab("genome site")+ ylab("number of genomes")+ ggtitle("B.1.106 mutational map") +
  scale_x_continuous(limits = c(0,29903),breaks = seq(0, 29903, by = 5000))+
  geom_hline(yintercept=61, linetype="dashed", color = "grey33") +
  ##sadly, this is the easier way
  annotate("text", fontface =2, x=B1106_tags$x_pos[1]-nx, y=B1106_tags$y_pos[1], label=B1106_tags$tag[1], angle=90, size=ss) +
  annotate("text", fontface =2, x=B1106_tags$x_pos[2]-nx, y=B1106_tags$y_pos[2], label=B1106_tags$tag[2], angle=90, size=ss) +
  annotate("text", fontface =2, x=B1106_tags$x_pos[3]-nx, y=B1106_tags$y_pos[3], label=B1106_tags$tag[3], angle=90, size=ss, col='tomato3') +
  annotate("text", fontface =2, x=B1106_tags$x_pos[4]+nx, y=B1106_tags$y_pos[4], label=B1106_tags$tag[4], angle=90, size=ss) +
  annotate("text", fontface =2, x=B1106_tags$x_pos[5]+nx, y=B1106_tags$y_pos[5], label=B1106_tags$tag[5], angle=90, size=ss)


panelBe<- ggplot()


#######################################################
#######################################################
#### PANEL 3C

mutation_metadata<-read_excel('data_original/world_mutations_dataset_16Sept.xlsx')
mutation_data<-read_excel('data_original/mutations_world.xlsx')
all_mutations_world<-merge(mutation_metadata,mutation_data)
all_mutations_world$date<-as.Date(cut(all_mutations_world$date, breaks = "week", start.on.monday = FALSE))
all_mutations_world<- all_mutations_world %>% mutate(sampling_date=as.POSIXct(date))
all_mutations_world12503<- subset(all_mutations_world, T12503C=='C')
all_mutations_world29721<- subset(all_mutations_world, C29721T=='T')
all_mutations_world22675<- subset(all_mutations_world, C22675T=='T')

all_time_range<- range(all_mutations_world12503$date, all_mutations_world29721$date, all_mutations_world22675$date)

panelC1<- ggplot(data=all_mutations_world12503) + ggstyleTILE +
          geom_bar(aes(x = date, fill=(country=='South Africa')), width=5, color='black', size=0.3)+
          xlab("month") + ylab('num. sequences') + ggtitle("12503: T to C") +
          scale_x_date(date_labels = "%b",date_breaks = "month", limits=all_time_range)+
          scale_fill_manual("Country",values=c('cyan','grey33'),labels = c("rest of the world", "South Africa")) +
          theme(legend.position=c(0.35,0.7))

panelC2<- ggplot(data=all_mutations_world29721) + ggstyleTILE +
          geom_bar(aes(x = date, fill=(country=='South Africa')), width=5, color='black', size=0.3)+
          xlab("month") + ylab('num. sequences') + ggtitle("29721: C to T") +
          scale_x_date(date_labels = "%b",date_breaks = "month", limits=all_time_range)+
          scale_fill_manual("Country",values=c('cyan','grey33'),labels = c("rest of the world", "South Africa")) +
          theme(legend.position=c(0.35,0.7))

panelC3<- ggplot(data=all_mutations_world22675) + ggstyleTILE +
          geom_bar(aes(x = date, fill=(country=='South Africa')), width=5, color='black', size=0.3)+
          xlab("month") + ylab('num. sequences') + ggtitle("22675: C to T") +
          scale_x_date(date_labels = "%b",date_breaks = "month", limits=all_time_range)+
          scale_fill_manual("Country",values=c('cyan','grey33'),labels = c("rest of the world", "South Africa")) +
          theme(legend.position=c(0.35,0.7))


panelC<- cowplot::plot_grid(panelC1,panelC2,panelC3,nrow=3)

#####################################################
#####################################################
#####################################################

left<- cowplot::plot_grid(panel3A,panelC1,panelC2,panelC3,nrow=4,labels=c("A","C","",""))
right<- cowplot::plot_grid(panelBe,panelB4,panelB1,panelB2,panelB3,nrow=5,rel_heights=c(3,3,3,3,3),labels=c("B","","",""))

all<- cowplot::plot_grid(left,right,ncol=2,rel_widths=c(1,2))

pdf("output/fig3.pdf",w=9,h=10,bg='white')
print(all)
a<-dev.off()
