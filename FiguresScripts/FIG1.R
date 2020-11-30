

source("functions.plotting.R")
require(tidyverse)
library("readxl")
library("lubridate")

####################################################
####################################################
## PANEL 1A

    ##tulio: genome sampling
    data2<-read_excel('data_original/SA_final_dataset_n=1365_v2_cluster_annotated_final.xlsx')
    data2$cluster<-factor(data2$cluster,levels = c("remaining","cluster1","cluster3","cluster4","cluster5"))
    data2$days<-as.Date(cut(data2$date, breaks = "day", start.on.monday = FALSE))
    data2$date<-as.Date(cut(data2$date, breaks = "week", start.on.monday = FALSE))
    data2$first_case<-as.Date(cut(data2$first_case, breaks = "day", start.on.monday = FALSE))
    data2$first_cluster<-as.Date(cut(data2$first_cluster,breaks = "day",start.on.monday = FALSE))
    data2$cluster_first_case<-as.Date(cut(data2$cluster_first_case, breaks = "day", start.on.monday = FALSE))
    data2$cluster_last_case<-as.Date(cut(data2$cluster_last_case,breaks = "day",start.on.monday = FALSE))
    data2<- data2 %>% mutate(Collection_date=as.POSIXct(days))
    data2<- subset(data2, !is.na(strain))


    ggplot()+
    geom_rug(data= data2, aes(x=days),color='dodgerblue3',alpha=0.2,outside =FALSE,length = unit(0.03, "npc"), show.legend=TRUE)


    ##Re estimates from the github repo
    re<- read.csv("data/ZAF-estimates.csv", sep=',')
    reSA<- re %>% filter(region=="ZAF") %>% filter(data_type=="Confirmed cases") %>% filter(estimate_type=="Cori_slidingWindow")
    reSA$date<- as.Date(reSA$date, format="%Y-%m-%d")

    ##tulio: covid incidence
    covid<- read.csv("data/covid19za_provincial_cumulative_timeline_confirmed_16Sept.csv", sep=',')
    covid<- covid %>% select(date, SA_daily)
    covid$date<- as.Date(covid$date, format="%d/%m/%Y")

    ##match dates between Re and cases data
    reSA$cases<- NA
    match_dates<- match(covid$date, reSA$date)
    match_dates<- which(!is.na(match_dates))
    reSA$cases[match_dates]<- covid$SA_daily[which(!is.na(match_dates))]

    ##define level dates TODO:: make these the official ones
    level5<- as.Date(c("2020/03/09","2020/05/01"), format="%Y/%m/%d")
    level4<- as.Date(c("2020/05/01","2020/06/01"), format="%Y/%m/%d")
    level3<- as.Date(c("2020/06/01","2020/08/26"), format="%Y/%m/%d")
    level2<- as.Date(c("2020/08/26","2020/09/09"), format="%Y/%m/%d")

    # collevels<- wes_palette("Royal2")[c(2,3,1,5,4)]
    collevels<- c("grey80","white","grey55","grey30","black")

    ##make the plot
    f<- 3500
    pp<- 0.4
    panelA<- ggplot(data=reSA) +  ggstyleTILE +
        geom_hline(yintercept=0, color='black', size=0.3) +
        geom_hline(yintercept=1, color='red3', linetype=2) +
        geom_segment(aes(x = level5[1], y = 0-pp, xend = level5[2], yend = 0-pp), color=collevels[5], size=2.5) +
          geom_segment(aes(x = level4[1], y = 0-pp, xend = level4[2], yend = 0-pp), color=collevels[4], size=2.5) +
          geom_segment(aes(x = level3[1], y = 0-pp, xend = level3[2], yend = 0-pp), color=collevels[3], size=2.5) +
          geom_segment(aes(x = level2[1], y = 0-pp, xend = level2[2], yend = 0-pp), color=collevels[1], size=2.5) +
        geom_text(aes(x=mean(level5), y=0.2-pp), label="LEVEL 5", color=collevels[5], label.size = 0, size=3) +
          geom_text(aes(x=mean(level4), y=0.2-pp), label="LEVEL 4", color=collevels[4], label.size = 0, size=3) +
          geom_text(aes(x=mean(level3), y=0.2-pp), label="LEVEL 3", color=collevels[3], label.size = 0, size=3) +
          geom_text(aes(x=mean(level2), y=0.2-pp), label="LEVEL 2", color=collevels[1], label.size = 0, size=3) +
        geom_ribbon(aes(x=date, ymin=median_R_lowHPD, ymax=median_R_highHPD), fill='tomato', alpha=0.2) +
        geom_line(aes(x=date, y=median_R_mean, col="re"), size=1) +
        geom_line(aes(x=date, y=median_R_lowHPD, col="re"), size=0.2) +
        geom_line(aes(x=date, y=median_R_highHPD, col="re"), size=0.2) +
        geom_line(aes(x=date, y=cases/f, col="cases"), size=1) +
        geom_rug(data= data2, aes(x=days, color="rug"),alpha=0.2,outside =FALSE,length = unit(0.03, "npc"), show.legend=TRUE) +
        scale_y_continuous(sec.axis = sec_axis(~.*f, name = "cases"), limits=c(0-pp,4.5)) +
        ggtitle("Epidemic and genomic data in South Africa") + ylab("Re") + xlab("month") +
        scale_color_manual("",values=c("black","red3","blue"), labels=c("Cases","Re","Genomes"))  +
        theme(legend.position=c(0.3,0.8),legend.text= element_text(size=12, family="Helvetica")) +
        scale_shape_manual("", values='_', labels="Level 5") +
        scale_x_date(labels = date_format("%b"), date_breaks = "1 month")



####################################################
####################################################
##### Panel 1B

        ##tulio: introduction data
        data5<-read_excel('data_original/SA_introductions_V2_dated.xlsx')
        data5$days<-as.Date(cut(data5$date, breaks = "day", start.on.monday = FALSE))
        data5$date<-as.Date(cut(data5$date, breaks = "week", start.on.monday = FALSE))
        data5$days<-as.Date(cut(data5$date, breaks = "day", start.on.monday = FALSE))
        data5$date<-as.Date(cut(data5$date, breaks = "week", start.on.monday = FALSE))
        # data5<- data5 %>% subset(data5, !is.na(Origin)) always true
        data5<- data5 %>% mutate(date=as.Date(date, format="%Y-%m-%d"))

        all<-data.frame(dataset=c(rep('South Africa',length(data2$strain)),rep('Outside',length(data5$date))),date=c(data2$date,data5$date),origin=c(data2$virus,data5$Origin_region))

        # region_colors<- c("#FF1026","#009076","#E9C341","#25201D","grey88","#597DC5")
        region_colors<- c("#FF1026","#009076","#BBA18B","#25201D","white","#597DC5")

        panelB<- ggplot(data=data5)+  ggstyleTILE +
          geom_bar(aes(x = date,fill=Origin_region), color='grey33',width=5)+
          scale_fill_manual("Origin", values=region_colors)+
          scale_x_date(date_labels = "%b",date_breaks = "month")+
          xlab('dates')+ ylab('count')+ ggtitle("Introduction into South Africa") +
          theme(legend.position=c(0.2,0.7),legend.text= element_text(size=12, family="Helvetica"))


####################################################
####################################################
##### Panel 1C

        panelC<-ggplot(all,aes(x=date,fill=dataset))+  ggstyleTILE +
          geom_bar(width=5,size=1)+
          scale_fill_manual("Origin",values=c('black','grey'))+
          scale_x_date(date_labels = "%b",date_breaks = "month")+
          xlab('month')+ ylab('count')+
          ggtitle("Sampling in South Africa") +
          theme(legend.position=c(0.2,0.8),legend.text= element_text(size=12, family="Helvetica"))


        top<- cowplot::plot_grid(panelA,panelB,panelC,labels=c("A","B","C"),ncol=3)
        bottom<- cowplot::plot_grid(gempty,labels=c("D"),ncol=1)
        all<- cowplot::plot_grid(top,bottom,nrow=2)

        pdf("output/fig1.pdf",w=14,h=9,bg='white')
        print(all)
        a<-dev.off()
