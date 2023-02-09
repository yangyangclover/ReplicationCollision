### plot hotspots of large TD from GISTIC ###
### library
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(metR)
library(ggrepel)
library(ggh4x)
library(ggpubr)
library(patchwork)
library(grid)
library(egg)


### Function
hotspot_function <- function(gisticpath,scorepath,top.num,output_file,titlename,chrname){
  centromeres <- read.delim("D:/Data/REF/hg19/centromeres/hg19.centromeres.bed",header = F)
  centromeres$V1 <- gsub("chr","",centromeres$V1)
  names(centromeres) <- c("Chromosome","Pos1","Pos2","Cen")
  
  chrlength <- read.delim("D:/Data/REF/hg19/genome_length/genome_length_hg19.tsv")
  names(chrlength) <- c("Chromosome","Length")
  # read  file
  gistic <- read.csv(gisticpath)
  score <- read.delim(scorepath)
  # remove na rows
  gistic <- gistic %>% filter(!is.na(q.value)) %>% as.data.frame()
  # add chr, start, end
  gistic$chr <- apply(gistic,1,function(x) unlist(strsplit(x["region"],":"))[1])
  gistic$start <- apply(gistic,1,function(x) as.numeric(unlist(strsplit(unlist(strsplit(x["region"],":"))[2],"-"))[1]))
  gistic$end <- apply(gistic,1,function(x) as.numeric(unlist(strsplit(unlist(strsplit(x["region"],":"))[2],"-"))[2]))
  # change to Grange
  gistic_grange <- GRanges(seqnames = gistic[,"chr"], ranges = IRanges(gistic[,"start"], gistic[,"end"]))
  gistic_grange$q.value <- gistic$q.value
  gistic_grange$region <- gistic$region
  gistic_grange <- unique(gistic_grange)
  gistic_grange$q.order <- seq(1,gistic_grange@elementMetadata@nrows,1) # add q.value order
  # keep oncogene/fusion gene
  gistic_oncogene <- gistic[grep("[oncogene|fusion]",gistic$Role.in.Cancer),]
  gistic_oncogene$gene <- ifelse(gistic_oncogene$Role.in.Cancer=="fusion",
                                 paste0(gistic_oncogene$Translocation.Partner,"-",gistic_oncogene$gene," fusion"),
                                 gistic_oncogene$gene)
  # add oncogene to gistic_grange
  gistic_grange$oncogene <- apply(data.frame(gistic_grange), 1, function(x) gistic_oncogene %>% filter(region==x["region"]) %>% pull(gene) %>% toString)
  # plot
  gistic_plot <- as.data.frame(gistic_grange)
  gistic_plot$seqnames <- gsub("chr","",gistic_plot$seqnames)
  gistic_plot$seqnames <- factor(gistic_plot$seqnames, levels = c(seq(1,22,1),"X","Y"))
  gistic_plot$log10.q.value <- -log10(gistic_plot$q.value)
  # gistic_plot$oncogene <- apply(gistic_plot,1,function(x) paste0(unlist(strsplit(x["oncogene"],", ")),collapse = "\n",sep=""))
  gistic_plot$oncogene <- apply(gistic_plot,1,function(x){
    all.gene <- unlist(strsplit(x["oncogene"],", "))
    first.gene <- all.gene[1]
    first.gene <- ifelse(length(all.gene)>1,paste0(first.gene," etc."),first.gene)
  })
  names(gistic_plot)[1] <- "Chromosome"
  score_plot <- score[score$Type=="Amp",]
  score_plot$Chromosome <- factor(score_plot$Chromosome, levels = c(seq(1,22,1),"X","Y"))
  centromeres$chrtitle1 <- ifelse(as.numeric(centromeres$Chromosome)%%2 %in% 1,centromeres$Chromosome,NA)
  centromeres$chrtitle2 <- ifelse(as.numeric(centromeres$Chromosome)%%2 %in% 0,centromeres$Chromosome,NA)
  centromeres$Chromosome <- factor(centromeres$Chromosome, levels = c(seq(1,22,1),"X","Y"))
  centromeres <- centromeres[centromeres$Chromosome %in% c(seq(1,22,1)),]
  chrlength$Chromosome <- factor(chrlength$Chromosome, levels = c(seq(1,22,1),"X","Y"))
  chrlength <- chrlength[chrlength$Chromosome %in% c(seq(1,22,1)),]
  if (chrname==TRUE){
    p <- ggplot()  + 
      scale_y_continuous(trans = "log2",
                         breaks = c(0.25,10,50),
                         labels = c(0.25,10,50),
                         limits = c(0.006,200),
                         sec.axis = sec_axis(~., name=titlename)
      ) +
      scale_x_continuous(expand = c(0,0))+
      geom_rect(data=score_plot,aes(fill = Chromosome),xmin = -Inf,xmax = Inf,
                ymin = -Inf,ymax = log2(50),alpha = 0.3) +
      scale_fill_manual(breaks = c(seq(1,22,1)),
                        values = rep(c("white","grey90"),11)) +
      geom_linerange(data=centromeres,aes(x=Pos1,ymin=0.006,ymax=52), alpha=0.5,linetype="dotted",color="black",size=0.1)+
      geom_linerange(data=chrlength,aes(x=Length,ymin=0.006,ymax=52), alpha=0.5,linetype="solid",color="white",size=0.1)+
      geom_linerange(data=chrlength,aes(x=0,ymin=0.006,ymax=52), alpha=0.5,linetype="solid",color="white",size=0.1)+
      geom_line(data=score_plot,aes(x=(Start+End)/2,y = X.log10.q.value.),size=0.01)+
      geom_point(data= gistic_plot %>% filter(q.order <= top.num & oncogene != ""),
                 aes(x=(start+end)/2,y=log10.q.value),
                 size=0.5,
                 color="red")+
      facet_grid2(.~Chromosome,
                  scales = "free_x",
                  space = "free",
                  strip = strip_vanilla(clip = "off")) + 
      geom_text(data= gistic_plot %>% filter(q.order <= top.num),
                      aes(x=(start+end)/2,y=log10.q.value,label=oncogene),
                      size = 2,
                      color = 'red',vjust=-1
                      )+
      geom_text(data=centromeres,aes(x=Pos1,y=80,label=chrtitle1),size=2.45)+
      geom_text(data=centromeres,aes(x=Pos1,y=200,label=chrtitle2),size=2.45)+
      coord_cartesian(clip = "off") +
      ylab(ifelse(chrname==TRUE,"-log10(q value)",""))+
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.y.left = element_text(color="black",size=7),
            axis.title.y.right = element_text(color="black",size=7,angle = 0,hjust = 0,vjust = 0.5),
            axis.text.y.left = element_text(color="black",size=7),
            axis.text.y.right = element_blank(),
            axis.ticks.length.x = unit(0,"mm"),
            axis.ticks.length.y.right = unit(0,"mm"),
            plot.title = element_text(color="black",size=7),
            # strip.text.x = element_text(size = 7, colour = "black", angle = 0),
            strip.text.x = element_blank(),
            legend.position = "none",
            panel.spacing.x = unit(0, "mm"),
            strip.background = element_blank(),
            strip.placement = "outside",
            plot.margin = unit(c(5, 5, 5, 5), "mm"),
            panel.background = element_blank())
    
  }
  if (chrname==FALSE){
    p <- ggplot()  + 
      scale_y_continuous(trans = "log2",
                         breaks = c(0.25,10,50),
                         labels = c(0.25,10,50),
                         limits = c(0.006,50),
                         sec.axis = sec_axis(~., name=titlename)
      ) +
      scale_x_continuous(expand = c(0,0))+
      geom_rect(data=score_plot,aes(fill = Chromosome),xmin = -Inf,xmax = Inf,
                ymin = -Inf,ymax = log2(50),alpha = 0.3) +
      scale_fill_manual(breaks = c(seq(1,22,1)),
                        values = rep(c("white","grey90"),11)) +
      geom_linerange(data=centromeres,aes(x=Pos1,ymin=0.006,ymax=50), alpha=0.5,linetype="dotted",color="black",size=0.1)+
      geom_linerange(data=chrlength,aes(x=Length,ymin=0.006,ymax=52), alpha=0.5,linetype="solid",color="white",size=0.1)+
      geom_linerange(data=chrlength,aes(x=0,ymin=0.006,ymax=52), alpha=0.5,linetype="solid",color="white",size=0.1)+
      geom_line(data=score_plot,aes(x=(Start+End)/2,y = X.log10.q.value.),size=0.01)+
      geom_point(data= gistic_plot %>% filter(q.order <= top.num & oncogene != ""),
                 aes(x=(start+end)/2,y=log10.q.value),
                 size=0.5,
                 color="red")+
      facet_grid2(.~Chromosome,
                  scales = "free_x",
                  space = "free",
                  strip = strip_vanilla(clip = "off")) + 
      geom_text(data= gistic_plot %>% filter(q.order <= top.num),
                aes(x=(start+end)/2,y=log10.q.value,label=oncogene),
                size = 2,
                color = 'red',vjust=-1
      )+
      coord_cartesian(clip = "off") +
      ylab(ifelse(chrname==TRUE,"-log10(q value)",""))+
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.y.left = element_blank(),
            axis.title.y.right = element_text(color="black",size=7,angle = 0,hjust = 0,vjust = 0.5),
            axis.text.y.left = element_text(color="black",size=7),
            axis.text.y.right = element_blank(),
            axis.ticks.length.x = unit(0,"mm"),
            axis.ticks.length.y.right = unit(0,"mm"),
            plot.title = element_text(color="black",size=7),
            # strip.text.x = element_text(size = 7, colour = "black", angle = 0),
            strip.text.x = element_blank(),
            legend.position = "none",
            panel.spacing.x = unit(0, "mm"),
            strip.background = element_blank(),
            strip.placement = "outside",
            plot.margin = unit(c(0, 5, 5, 5), "mm"),
            panel.background = element_blank())
    
  }
  
  ggsave(output_file,width = 23,height = 18/15)
  return(p)
  
}


#########  RUN 
### PCAWG PANCAN ###
# parameter1
gisticpath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/pcawg/pancan/conf_90/amp_genes.conf_90_long.csv"
# parameter 2
scorepath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/pcawg/pancan/conf_90/scores.gistic"
# parameter3
top.num <- 10
# parameter4
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/pcawg/pancan/conf_90/pcawg_pancan_largetd_gistic_hotspot10.pdf"
hotspot_function(gisticpath,scorepath,top.num,output_file,"PCAWG\nPancancer",TRUE)
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/largetd.gistic.hotspot.merge/pcawg_pancan_largetd_gistic_hotspot10.pdf"
pcawg.pan.10 <- hotspot_function(gisticpath,scorepath,top.num,output_file,"PCAWG\nPancancer",TRUE)

### PCAWG Breast-AdenoCA ###
# parameter1
gisticpath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/pcawg/Breast-AdenoCA/conf_90/amp_genes.conf_90_long.csv"
# parameter 2
scorepath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/pcawg/Breast-AdenoCA/conf_90/scores.gistic"
# parameter3
top.num <- 10
# parameter4
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/pcawg/Breast-AdenoCA/conf_90/pcawg_Breast-AdenoCA_largetd_gistic_hotspot10.pdf"
hotspot_function(gisticpath,scorepath,top.num,output_file,"PCAWG\nBreast-AdenoCA",FALSE)
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/largetd.gistic.hotspot.merge/pcawg_Breast-AdenoCA_largetd_gistic_hotspot10.pdf"
pcawg.breast.10 <- hotspot_function(gisticpath,scorepath,top.num,output_file,"PCAWG\nBreast-AdenoCA",FALSE)

### PCAWG Eso-AdenoCA ###
# parameter1
gisticpath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/pcawg/Eso-AdenoCA/conf_90/amp_genes.conf_90_long.csv"
# parameter 2
scorepath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/pcawg/Eso-AdenoCA/conf_90/scores.gistic"
# parameter3
top.num <- 10
# parameter4
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/pcawg/Eso-AdenoCA/conf_90/pcawg_Eso-AdenoCA_largetd_gistic_hotspot10.pdf"
hotspot_function(gisticpath,scorepath,top.num,output_file,"PCAWG Eso-AdenoCA",FALSE)
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/largetd.gistic.hotspot.merge/pcawg_Eso-AdenoCA_largetd_gistic_hotspot10.pdf"
hotspot_function(gisticpath,scorepath,top.num,output_file,"PCAWG Eso-AdenoCA",FALSE)


### PCAWG Liver-HCC ###
# parameter1
gisticpath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/pcawg/Liver-HCC/conf_90/amp_genes.conf_90_long.csv"
# parameter 2
scorepath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/pcawg/Liver-HCC/conf_90/scores.gistic"
# parameter3
top.num <- 10
# parameter4
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/pcawg/Liver-HCC/conf_90/pcawg_Liver-HCC_largetd_gistic_hotspot10.pdf"
hotspot_function(gisticpath,scorepath,top.num,output_file,"PCAWG Liver-HCC",FALSE)
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/largetd.gistic.hotspot.merge/pcawg_Liver-HCC_largetd_gistic_hotspot10.pdf"
hotspot_function(gisticpath,scorepath,top.num,output_file,"PCAWG Liver-HCC",FALSE)

### PCAWG Ovary-AdenoCA ###
# parameter1
gisticpath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/pcawg/Ovary-AdenoCA/conf_90/amp_genes.conf_90_long.csv"
# parameter 2
scorepath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/pcawg/Ovary-AdenoCA/conf_90/scores.gistic"
# parameter3
top.num <- 10
# parameter4
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/pcawg/Ovary-AdenoCA/conf_90/pcawg_Ovary-AdenoCA_largetd_gistic_hotspot10.pdf"
hotspot_function(gisticpath,scorepath,top.num,output_file,"PCAWG\nOvary-AdenoCA",FALSE)
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/largetd.gistic.hotspot.merge/pcawg_Ovary-AdenoCA_largetd_gistic_hotspot10.pdf"
pcawg.ov.10 <- hotspot_function(gisticpath,scorepath,top.num,output_file,"PCAWG\nOvary-AdenoCA",FALSE)

### PCAWG Prost-AdenoCA ###
# parameter1
gisticpath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/pcawg/Prost-AdenoCA/conf_90/amp_genes.conf_90_long.csv"
# parameter 2
scorepath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/pcawg/Prost-AdenoCA/conf_90/scores.gistic"
# parameter3
top.num <- 10
# parameter4
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/pcawg/Prost-AdenoCA/conf_90/pcawg_Prost-AdenoCA_largetd_gistic_hotspot10.pdf"
hotspot_function(gisticpath,scorepath,top.num,output_file,"PCAWG\nProst-AdenoCA",FALSE)
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/largetd.gistic.hotspot.merge/pcawg_Prost-AdenoCA_largetd_gistic_hotspot10.pdf"
pcawg.prost.10 <- hotspot_function(gisticpath,scorepath,top.num,output_file,"PCAWG\nProst-AdenoCA",FALSE)

### PCAWG Stomach-AdenoCA ###
# parameter1
gisticpath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/pcawg/Stomach-AdenoCA/conf_90/amp_genes.conf_90_long.csv"
# parameter 2
scorepath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/pcawg/Stomach-AdenoCA/conf_90/scores.gistic"
# parameter3
top.num <- 10
# parameter4
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/pcawg/Stomach-AdenoCA/conf_90/pcawg_Stomach-AdenoCA_largetd_gistic_hotspot10.pdf"
hotspot_function(gisticpath,scorepath,top.num,output_file,"PCAWG Stomach-AdenoCA",FALSE)
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/largetd.gistic.hotspot.merge/pcawg_Stomach-AdenoCA_largetd_gistic_hotspot10.pdf"
hotspot_function(gisticpath,scorepath,top.num,output_file,"PCAWG Stomach-AdenoCA",FALSE)

### PCAWG Uterus-AdenoCA ###
# parameter1
gisticpath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/pcawg/Uterus-AdenoCA/conf_90/amp_genes.conf_90_long.csv"
# parameter 2
scorepath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/pcawg/Uterus-AdenoCA/conf_90/scores.gistic"
# parameter3
top.num <- 10
# parameter4
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/pcawg/Uterus-AdenoCA/conf_90/pcawg_Uterus-AdenoCA_largetd_gistic_hotspot10.pdf"
hotspot_function(gisticpath,scorepath,top.num,output_file,"PCAWG\nUterus-AdenoCA",FALSE)
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/largetd.gistic.hotspot.merge/pcawg_Uterus-AdenoCA_largetd_gistic_hotspot10.pdf"
pcawg.uterus.10 <- hotspot_function(gisticpath,scorepath,top.num,output_file,"PCAWG\nUterus-AdenoCA",FALSE)

### Hartwig PANCAN ###
# parameter1
gisticpath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/hartwig/pancan/conf_90/amp_genes.conf_90_long.csv"
# parameter 2
scorepath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/hartwig/pancan/conf_90/scores.gistic"
# parameter3
top.num <- 10
# parameter4
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/hartwig/pancan/conf_90/hartwig_pancan_largetd_gistic_hotspot10.pdf"
hotspot_function(gisticpath,scorepath,top.num,output_file,"Hartwig\nPancancer",FALSE)
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/largetd.gistic.hotspot.merge/hartwig_pancan_largetd_gistic_hotspot10.pdf"
hartwig.pan.10 <- hotspot_function(gisticpath,scorepath,top.num,output_file,"Hartwig\nPancancer",FALSE)

### Hartwig Breast ###
# parameter1
gisticpath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/hartwig/Breast/conf_90/amp_genes.conf_90_long.csv"
# parameter 2
scorepath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/hartwig/Breast/conf_90/scores.gistic"
# parameter3
top.num <- 10
# parameter4
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/hartwig/Breast/conf_90/hartwig_Breast_largetd_gistic_hotspot10.pdf"
hotspot_function(gisticpath,scorepath,top.num,output_file,"Hartwig\nBreast",FALSE)
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/largetd.gistic.hotspot.merge/hartwig_Breast_largetd_gistic_hotspot10.pdf"
hartwig.breast.10 <- hotspot_function(gisticpath,scorepath,top.num,output_file,"Hartwig\nBreast",FALSE)

### Hartwig Esophagus ###
# parameter1
gisticpath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/hartwig/Esophagus/conf_90/amp_genes.conf_90_long.csv"
# parameter 2
scorepath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/hartwig/Esophagus/conf_90/scores.gistic"
# parameter3
top.num <- 10
# parameter4
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/hartwig/Esophagus/conf_90/hartwig_Esophagus_largetd_gistic_hotspot10.pdf"
hotspot_function(gisticpath,scorepath,top.num,output_file,"Hartwig Esophagus",FALSE)
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/largetd.gistic.hotspot.merge/hartwig_Esophagus_largetd_gistic_hotspot10.pdf"
hotspot_function(gisticpath,scorepath,top.num,output_file,"Hartwig Esophagus",FALSE)

### Hartwig Ovary ###
# parameter1
gisticpath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/hartwig/Ovary/conf_90/amp_genes.conf_90_long.csv"
# parameter 2
scorepath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/hartwig/Ovary/conf_90/scores.gistic"
# parameter3
top.num <- 10
# parameter4
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/hartwig/Ovary/conf_90/hartwig_Ovary_largetd_gistic_hotspot10.pdf"
hotspot_function(gisticpath,scorepath,top.num,output_file,"Hartwig\nOvary",FALSE)
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/largetd.gistic.hotspot.merge/hartwig_Ovary_largetd_gistic_hotspot10.pdf"
hartwig.ov.10 <- hotspot_function(gisticpath,scorepath,top.num,output_file,"Hartwig\nOvary",FALSE)

### Hartwig Prostate ###
# parameter1
gisticpath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/hartwig/Prostate/conf_90/amp_genes.conf_90_long.csv"
# parameter 2
scorepath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/hartwig/Prostate/conf_90/scores.gistic"
# parameter3
top.num <- 10
# parameter4
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/hartwig/Prostate/conf_90/hartwig_Prostate_largetd_gistic_hotspot10.pdf"
hotspot_function(gisticpath,scorepath,top.num,output_file,"Hartwig\nProstate",FALSE)
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/largetd.gistic.hotspot.merge/hartwig_Prostate_largetd_gistic_hotspot10.pdf"
hartwig.prost.10 <- hotspot_function(gisticpath,scorepath,top.num,output_file,"Hartwig\nProstate",FALSE)

### Hartwig Urothelial-tract ###
# parameter1
gisticpath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/hartwig/Urothelial-tract/conf_90/amp_genes.conf_90_long.csv"
# parameter 2
scorepath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/hartwig/Urothelial-tract/conf_90/scores.gistic"
# parameter3
top.num <- 10
# parameter4
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/hartwig/Urothelial-tract/conf_90/hartwig_Urothelial-tract_largetd_gistic_hotspot10.pdf"
hotspot_function(gisticpath,scorepath,top.num,output_file,"Hartwig\nUrothelial-tract",FALSE)
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/largetd.gistic.hotspot.merge/hartwig_Urothelial-tract_largetd_gistic_hotspot10.pdf"
hartwig.urothelial.10 <- hotspot_function(gisticpath,scorepath,top.num,output_file,"Hartwig\nUrothelial-tract",FALSE)

### Hartwig Uterus ###
# parameter1
gisticpath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/hartwig/Uterus/conf_90/amp_genes.conf_90_long.csv"
# parameter 2
scorepath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/hartwig/Uterus/conf_90/scores.gistic"
# parameter3
top.num <- 10
# parameter4
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/hartwig/Uterus/conf_90/hartwig_Uterus_largetd_gistic_hotspot10.pdf"
hotspot_function(gisticpath,scorepath,top.num,output_file,"Hartwig Uterus",FALSE)
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/largetd.gistic.hotspot.merge/hartwig_Uterus_largetd_gistic_hotspot10.pdf"
hotspot_function(gisticpath,scorepath,top.num,output_file,"Hartwig Uterus",FALSE)

### POG570 PANCAN ###
# parameter1
gisticpath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/pog570/pancan/conf_90/amp_genes.conf_90_long.csv"
# parameter 2
scorepath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/pog570/pancan/conf_90/scores.gistic"
# parameter3
top.num <- 10
# parameter4
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/pog570/pancan/conf_90/pog570_pancan_largetd_gistic_hotspot10.pdf"
hotspot_function(gisticpath,scorepath,top.num,output_file,"POG570\nPancancer",FALSE)
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/largetd.gistic.hotspot.merge/pog570_pancan_largetd_gistic_hotspot10.pdf"
pog570.pan.10 <-hotspot_function(gisticpath,scorepath,top.num,output_file,"POG570\nPancancer",FALSE)

### POG570 BRCA ###
# parameter1
gisticpath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/pog570/BRCA/conf_90/amp_genes.conf_90_long.csv"
# parameter 2
scorepath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/pog570/BRCA/conf_90/scores.gistic"
# parameter3
top.num <- 10
# parameter4
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/pog570/BRCA/conf_90/pog570_BRCA_largetd_gistic_hotspot10.pdf"
hotspot_function(gisticpath,scorepath,top.num,output_file,"POG570\nBRCA",FALSE)
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/largetd.gistic.hotspot.merge/pog570_BRCA_largetd_gistic_hotspot10.pdf"
pog570.breast.10 <- hotspot_function(gisticpath,scorepath,top.num,output_file,"POG570\nBRCA",FALSE)

###  BRCA-EU ###
# parameter1
gisticpath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/BRCA-EU/BRCA-EU/conf_90/amp_genes.conf_90_long.csv"
# parameter 2
scorepath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/BRCA-EU/BRCA-EU/conf_90/scores.gistic"
# parameter3
top.num <- 10
# parameter4
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/BRCA-EU/BRCA-EU/conf_90/BRCA-EU_largetd_gistic_hotspot10.pdf"
hotspot_function(gisticpath,scorepath,top.num,output_file,"ICGC\nBRCA-EU",FALSE)
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/largetd.gistic.hotspot.merge/BRCA-EU_largetd_gistic_hotspot10.pdf"
brcaeu.10 <- hotspot_function(gisticpath,scorepath,top.num,output_file,"ICGC\nBRCA-EU",FALSE)

###  OV-AU ###
# parameter1
gisticpath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/OV-AU/OV-AU/conf_90/amp_genes.conf_90_long.csv"
# parameter 2
scorepath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/OV-AU/OV-AU/conf_90/scores.gistic"
# parameter3
top.num <- 10
# parameter4
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/OV-AU/OV-AU/conf_90/OV-AU_largetd_gistic_hotspot10.pdf"
hotspot_function(gisticpath,scorepath,top.num,output_file,"ICGC\nOV-AU",FALSE)
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/largetd.gistic.hotspot.merge/OV-AU_largetd_gistic_hotspot10.pdf"
ovau.10 <- hotspot_function(gisticpath,scorepath,top.num,output_file,"ICGC\nOV-AU",FALSE)

###  PACA-AU ###
# parameter1
gisticpath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/PACA-AU/PACA-AU/conf_90/amp_genes.conf_90_long.csv"
# parameter 2
scorepath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/PACA-AU/PACA-AU/conf_90/scores.gistic"
# parameter3
top.num <- 10
# parameter4
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/PACA-AU/PACA-AU/conf_90/PACA-AU_largetd_gistic_hotspot10.pdf"
hotspot_function(gisticpath,scorepath,top.num,output_file,"ICGC PACA-AU",FALSE)
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/largetd.gistic.hotspot.merge/PACA-AU_largetd_gistic_hotspot10.pdf"
hotspot_function(gisticpath,scorepath,top.num,output_file,"ICGC PACA-AU",FALSE)

###  LIRI-JP ###
# parameter1
gisticpath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/LIRI-JP/LIRI-JP/conf_90/amp_genes.conf_90_long.csv"
# parameter 2
scorepath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/LIRI-JP/LIRI-JP/conf_90/scores.gistic"
# parameter3
top.num <- 10
# parameter4
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/LIRI-JP/LIRI-JP/conf_90/LIRI-JP_largetd_gistic_hotspot10.pdf"
hotspot_function(gisticpath,scorepath,top.num,output_file,"ICGC\nLIRI-JP",FALSE)
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/largetd.gistic.hotspot.merge/LIRI-JP_largetd_gistic_hotspot10.pdf"
hotspot_function(gisticpath,scorepath,top.num,output_file,"ICGC\nLIRI-JP",FALSE)

###  PRAD-UK ###
# parameter1
gisticpath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/PRAD-UK/PRAD-UK/conf_90/amp_genes.conf_90_long.csv"
# parameter 2
scorepath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/PRAD-UK/PRAD-UK/conf_90/scores.gistic"
# parameter3
top.num <- 10
# parameter4
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/PRAD-UK/PRAD-UK/conf_90/PRAD-UK_largetd_gistic_hotspot10.pdf"
hotspot_function(gisticpath,scorepath,top.num,output_file,"ICGC\nPRAD-UK",FALSE)
output_file <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/largetd.gistic.hotspot.merge/PRAD-UK_largetd_gistic_hotspot10.pdf"
praduk.10 <- hotspot_function(gisticpath,scorepath,top.num,output_file,"ICGC\nPRAD-UK",FALSE)

### merge
s <- plot_spacer() + theme(panel.background = element_blank())


p2 <- ggarrange(pcawg.pan.10,pcawg.breast.10,pcawg.ov.10,
                pcawg.prost.10,pcawg.uterus.10,
                hartwig.pan.10,hartwig.breast.10,hartwig.ov.10,
                hartwig.prost.10,hartwig.urothelial.10,
                pog570.pan.10,pog570.breast.10,
                brcaeu.10,ovau.10,praduk.10,
                heights = c(1.2,rep(1,14)),
                nrow = 15)
ggsave("D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.gistic/largetd.gistic.hotspot.merge/top10.mainplots.1col.pdf",
       p2,height = 230,width=180,units = "mm",dpi=500)
