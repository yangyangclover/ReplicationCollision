# library
library(dplyr)
library(ggplot2)
library(ggh4x)
library(reshape)
library(ggpubr)

# function
plot_largetd_across_tumors_function <- function(simple_sv,
                                                cohort,
                                                output_path,
                                                allsample=NA,
                                                histology=NA,
                                                largetd_group){
  if (!dir.exists(output_path)){
    dir.create(output_path,showWarnings = T,recursive = T)
  }
  
  names(simple_sv)[grep("primaryTumorLocation",names(simple_sv))] <- "histology_abbreviation"
  
  #  large dup:  td10-16
  simple_sv$large_td <- ifelse(simple_sv$sig_group %in% paste0("td_",seq(10,16,1)),1,0)
  sample_largetd <- simple_sv %>% group_by(sample,large_td) %>% count()  %>% filter(large_td==1)
  
  # td3  and td4
  if (largetd_group=="pcawg"){
    simple_sv$td3 <- ifelse(simple_sv$sig_group %in% paste0("td_",seq(10,11,1)),1,0)
    simple_sv$td4 <- ifelse(simple_sv$sig_group %in% paste0("td_",seq(12,16,1)),1,0)
  }
  if (largetd_group=="hartwig"){
    simple_sv$td3 <- ifelse(simple_sv$sig_group %in% paste0("td_",seq(10,12,1)),1,0)
    simple_sv$td4 <- ifelse(simple_sv$sig_group %in% paste0("td_",seq(13,16,1)),1,0)
  }
  sample_td3 <- simple_sv %>% group_by(sample,td3) %>% count()  %>% filter(td3==1)
  sample_td4 <- simple_sv %>% group_by(sample,td4) %>% count()  %>% filter(td4==1)
  
  # add large td num to sample_largetd_all
  if (is.na(allsample[1])){
    sample_largetd_all  <- data.frame("sample"=unique(simple_sv$sample))
  } else {
    sample_largetd_all <- data.frame("sample"=allsample)
  }
  sample_largetd_all <- merge(sample_largetd_all,
                              sample_largetd[,c("sample","n")],
                              by.x = "sample",
                              by.y = "sample",
                              all.x = T)
  names(sample_largetd_all)[ncol(sample_largetd_all)] <- "largetd.num"
  sample_largetd_all$largetd.num <- ifelse(is.na(sample_largetd_all$largetd.num),0,sample_largetd_all$largetd.num)
  
  # add td3/td4 to sample_largetd_all
  sample_largetd_all <- merge(sample_largetd_all,
                              sample_td3[,c("sample","n")],
                              by.x = "sample",
                              by.y = "sample",
                              all.x = T)
  names(sample_largetd_all)[ncol(sample_largetd_all)] <- "td3"
  sample_largetd_all$td3 <- ifelse(is.na(sample_largetd_all$td3),0,sample_largetd_all$td3)
  sample_largetd_all <- merge(sample_largetd_all,
                              sample_td4[,c("sample","n")],
                              by.x = "sample",
                              by.y = "sample",
                              all.x = T)
  names(sample_largetd_all)[ncol(sample_largetd_all)] <- "td4"
  sample_largetd_all$td4 <- ifelse(is.na(sample_largetd_all$td4),0,sample_largetd_all$td4)
  
  # add histology
  histology_num <- histology %>% group_by(WGS,histology_abbreviation)%>% count()%>% group_by(histology_abbreviation) %>% count()
  histology_num$fullname <- paste0(histology_num$histology_abbreviation," (",histology_num$n,")")

  
  sample_histology <- histology %>% group_by(WGS,histology_abbreviation) %>% count()%>% select(c("WGS","histology_abbreviation"))
  sample_simplesv <-  simple_sv %>% group_by(sample) %>% count() 
  sample_histology <- merge(sample_histology,sample_simplesv,
                            by.x = "WGS",by.y = "sample",all.x = T)
  names(sample_histology)[ncol(sample_histology)] <- "simple.sv.num"
  sample_histology$simple.sv.num <- ifelse(is.na(sample_histology$simple.sv.num),0,sample_histology$simple.sv.num)
  
  sample_histology <- merge(sample_histology,
                            histology_num[,c("histology_abbreviation","fullname")])
  sample_largetd_all <- merge(sample_largetd_all,
                              sample_histology,
                              by.x = "sample",
                              by.y = "WGS",
                              all.x = T)
  
  
  # tumortype largetd frequency
  tumortype_largetd_num <- sample_largetd_all %>%
    group_by(fullname) %>%
    dplyr::summarize(large_td.num = sum(largetd.num>0, na.rm = TRUE))
  tumortype_largetd_freq10 <- sample_largetd_all %>% filter(largetd.num >= 10) %>% 
    group_by(fullname) %>%
    dplyr::summarize(freq10 = sum(largetd.num>0, na.rm = TRUE))
  tumortype_largetd_num <- merge(tumortype_largetd_num,tumortype_largetd_freq10,all.x = T)
  tumortype_largetd_num <- merge(tumortype_largetd_num,histology_num[,c("fullname","n")])
  tumortype_largetd_num$freq10 <- ifelse(is.na(tumortype_largetd_num$freq10),0,tumortype_largetd_num$freq10)
  tumortype_largetd_num$freq10.percentage <- tumortype_largetd_num$freq10/tumortype_largetd_num$n
  tumortype_largetd_num$freq10.percentage2 <- paste(sprintf("%.0f",tumortype_largetd_num$freq10.percentage*100),"%",sep="")
  tumortype_largetd_num$fullname2 <- paste0(
    substr(tumortype_largetd_num$fullname,1,nchar(tumortype_largetd_num$fullname)-1),
    ", ", tumortype_largetd_num$freq10.percentage2,")")
  
  # add fullnames
  sample_largetd_all <- merge(sample_largetd_all,
                              tumortype_largetd_num[,c("fullname","fullname2")],
                              all.x=T)
  
  # filter
  keep_tumortype <- histology_num %>% filter(n>=20 & histology_abbreviation !="Unknown") %>% pull(histology_abbreviation)
  sample_largetd_all <- sample_largetd_all[sample_largetd_all$histology_abbreviation  %in% keep_tumortype,]
  
  sample_largetd_all_long <- melt(sample_largetd_all,measure.vars= c("td3","td4"))
  
  # sort
  tumortype_order <- tumortype_largetd_num[order(tumortype_largetd_num$freq10.percentage,decreasing = T),] %>% pull(fullname2)
  sample_largetd_all <- sample_largetd_all[order(sample_largetd_all$largetd.num,decreasing = T),]
  sample_largetd_all$sample <- factor(sample_largetd_all$sample,levels = unique(sample_largetd_all$sample))
  sample_largetd_all$fullname2 <- factor(sample_largetd_all$fullname2, levels = tumortype_order)
  sample_largetd_all_long$sample <- factor(sample_largetd_all_long$sample,levels = unique(sample_largetd_all$sample))
  sample_largetd_all_long$fullname2 <- factor(sample_largetd_all_long$fullname2, levels = tumortype_order)
  
  # plot
  p1  <-  ggplot(data=sample_largetd_all,aes(x=sample,y=largetd.num)) + 
    geom_point(size=0.05)  +
    scale_y_log10(
      # limits = c(1,1000)
                  ) + 
    facet_grid2(. ~ fullname2,
                scales = "free_x",
                space = "free_x",
                switch = "both",
                strip = strip_vanilla(clip = "off")) + 
    ylab("Number of large TDs") + 
    theme_bw()+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=6,colour = "black"),
          axis.title.y = element_text(size=6,colour = "black"),
          axis.title.x = element_blank(),
          axis.line.y = element_line(size=0.1,colour = "black"), 
          axis.line.x.top = element_line(size=0.1,colour = "black"), 
          axis.line.x.bottom = element_line(size=0.1,colour = "black"), 
          axis.ticks.x =  element_blank(),
          axis.ticks.length.y =unit(0.8, "mm"),
          panel.background = element_blank(),
          # panel.margin.x=unit(0.25, "cm"),
          panel.spacing.x = unit(0.2, "cm"),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          strip.text.x = element_text(size=6, angle = 90, colour = "black",hjust = 1,vjust = 0.5),
          strip.background =  element_blank(),
          strip.placement = "top",
          plot.margin=grid::unit(c(0,2,0,2), "mm")
    )
  
  p2  <- ggplot(data=sample_largetd_all_long,aes(x=sample,y=value,fill=variable)) + 
    geom_bar(stat='identity',position='fill')  +
    scale_y_continuous(labels = scales::percent)+
    scale_fill_manual(values = c("#F1C44B","#E89080"),
                      breaks = c("td3","td4"),
                      labels  = c("TD3","TD4")) +
    facet_grid2(. ~ fullname2,
                scales = "free_x",
                space = "free_x",
                switch = "both",
                strip = strip_vanilla(clip = "off")) +
    ylab("")+
    xlab("")+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=6,colour = "black",vjust = 1),
          axis.title.y = element_text(size=4,colour = "black",vjust = 1),
          axis.title.x = element_text(size=4,colour = "black",vjust = 1),
          axis.line.y = element_blank(), 
          # axis.line.x.top = element_line(size=0.1,colour = "black"),
          axis.ticks.x =  element_blank(),
          axis.ticks.y =  element_blank(),
          legend.key = element_rect(fill = "transparent"),
          legend.key.size = unit(0.15, "cm"),
          legend.text = element_text(size=6,colour = "black"),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.margin=margin(t=-20),
          legend.box.margin=margin(0,0,0,0),
          panel.background = element_blank(),
          # panel.margin.x=unit(0.25, "cm"),
          # panel.margin.y=unit(0, "cm"),
          panel.spacing.x = unit(0.2, "cm"),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          strip.text.x = element_blank(),
          strip.background = element_blank(),
          strip.placement = "bottom",
          plot.margin=grid::unit(c(-3,2,1,2), "mm")
    )
  
  p <-  ggarrange(p1,p2,
            ncol = 1, nrow = 2,
            heights = c(1, 0.35))
  
  ggsave(plot=p,file.path(output_path,paste0(cohort,"_large_td_across_tumor_types.pdf")),width = 165,height = 90,units = c("mm"))
}




##### RUN
pcawg_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature"
output_path <- file.path(pcawg_path,"scratch/large_dup/pcawg/008-largetd/freq")

# pcawg
simple_sv_path <- file.path(pcawg_path,"scratch/large_dup/pcawg/004_signature_matrix","pcawg_simple_sv_with_49svsig_simple.txt")
cohort <- "pcawg"
pcawg2583_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/pcawg2583_sampleid_match.csv"
pcawg2583 <- read.csv(pcawg2583_path, fileEncoding="UTF-8-BOM")
simple_sv <- read.delim(simple_sv_path)
# simple_sv <- simple_sv[simple_sv$sample %in% unique(pcawg2583$sample),]
plot_largetd_across_tumors_function(simple_sv=simple_sv,
                                    cohort=cohort,
                                    output_path=output_path,
                                    allsample=unique(pcawg2583$WGS),
                                    histology=pcawg2583,
                                    largetd_group="pcawg")

# hartwig
simple_sv_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/pub.simple_sv.siggroup.tumortype.bedpe"
pub.sample.name_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/hartwig2367_sampleid_match.csv"
cohort <- "hartwig"
simple_sv <- read.delim(simple_sv_path)
pub.sample.name <- read.csv(pub.sample.name_path,fileEncoding = "UTF-8-BOM")
plot_largetd_across_tumors_function(simple_sv=simple_sv,
                                    cohort=cohort,
                                    output_path=output_path,
                                    allsample=unique(pub.sample.name$WGS),
                                    histology=pub.sample.name,
                                    largetd_group="hartwig")

