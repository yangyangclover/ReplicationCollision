### library
library(dplyr)
library(ggplot2)
library(grid)

### function
svenrichment_plot_function <- function(svenrichment,mysig,myexpression,plottitle){
  signamenew  <- c("Del0","Del1","Del2","Del3","Del4",
                   "TD0","TD1","TD2","TD3","TD4",
                   "Fragile sites",
                   "Fb inv",
                   "Unbal inv",
                   "Large mixed",
                   "Unbal tra")
  names(signamenew) <- c("del0","del1","del2","del3","del4",
                         "dup0","dup1","dup2","dup3","dup4",
                         "fragile",
                         "fbinv",
                         "unbalinv",
                         "largemixed",
                         "unbaltra")
  
  svenrichment$signaturenew <- signamenew[as.character(svenrichment$signature)]
  
  svenrichment_mysig <- svenrichment %>% filter(expression==myexpression & signature != "all")
  svenrichment_mysig$q <- p.adjust(svenrichment_mysig$raw.p,method = "bonferroni",n=nrow(svenrichment_mysig))
  svenrichment_mysig$logq <- -log10(svenrichment_mysig$q)
  svenrichment_mysig$logq <- ifelse(is.infinite(svenrichment_mysig$logq),309,svenrichment_mysig$logq)
  svenrichment_mysig$color <- "grey"
  svenrichment_mysig$color <-  ifelse(svenrichment_mysig$or >1 & svenrichment_mysig$q < 0.05,"red",svenrichment_mysig$color)
  svenrichment_mysig$color <-  ifelse(svenrichment_mysig$or <1  & svenrichment_mysig$q < 0.05,"blue",svenrichment_mysig$color)
  svenrichment_mysig$color <- factor(svenrichment_mysig$color, levels = c("red", "blue", "grey"))
  svenrichment_mysig$n <- svenrichment_mysig$a + svenrichment_mysig$b
  
  maxlogq <- max(svenrichment_mysig$logq)
  # sort
  svenrichment_mysig$signaturenew <- factor(svenrichment_mysig$signaturenew,levels = rev(signamenew))
  # x.lable
  x.breaks <- if(maxlogq<100){seq(0,300,100)}else{seq(0,300,100)}
  
  num.font <- 6
  hjust.text <- 1
  yposition1 <- -25
  ggplot(data = svenrichment_mysig,
         aes(y=signaturenew,x=logq,fill=color))+
    geom_bar(stat="identity") + 
    # geom_vline(xintercept=-log10(0.05), linetype="dashed", color = "black")+
    scale_fill_manual(values = c("red", "blue", "grey"),
                      breaks = c("red", "blue", "grey"),
                      labels = c(paste0("Enriched in ","observed SVs"),
                                 paste0("Depleted in ","observed SVs"),
                                 "Not significant"),
                      drop=FALSE)+
    
    annotation_custom(textGrob("n", gp=gpar(fontsize=num.font, fontface="bold"), rot = 0,vjust = 0,hjust = hjust.text),ymin=-.25,ymax=-.25,xmin=-3,xmax=-3) +
    annotation_custom(textGrob(paste(svenrichment_mysig[1,c("n")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(svenrichment_mysig),ymax=nrow(svenrichment_mysig),xmin=yposition1,xmax=yposition1) +
    annotation_custom(textGrob(paste(svenrichment_mysig[2,c("n")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(svenrichment_mysig)-1,ymax=nrow(svenrichment_mysig)-1,xmin=yposition1,xmax=yposition1) +
    annotation_custom(textGrob(paste(svenrichment_mysig[3,c("n")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(svenrichment_mysig)-2,ymax=nrow(svenrichment_mysig)-2,xmin=yposition1,xmax=yposition1) +
    annotation_custom(textGrob(paste(svenrichment_mysig[4,c("n")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(svenrichment_mysig)-3,ymax=nrow(svenrichment_mysig)-3,xmin=yposition1,xmax=yposition1) +
    annotation_custom(textGrob(paste(svenrichment_mysig[5,c("n")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(svenrichment_mysig)-4,ymax=nrow(svenrichment_mysig)-4,xmin=yposition1,xmax=yposition1) +
    annotation_custom(textGrob(paste(svenrichment_mysig[6,c("n")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(svenrichment_mysig)-5,ymax=nrow(svenrichment_mysig)-5,xmin=yposition1,xmax=yposition1) +
    annotation_custom(textGrob(paste(svenrichment_mysig[7,c("n")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(svenrichment_mysig)-6,ymax=nrow(svenrichment_mysig)-6,xmin=yposition1,xmax=yposition1) +
    annotation_custom(textGrob(paste(svenrichment_mysig[8,c("n")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(svenrichment_mysig)-7,ymax=nrow(svenrichment_mysig)-7,xmin=yposition1,xmax=yposition1) +
    annotation_custom(textGrob(paste(svenrichment_mysig[9,c("n")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(svenrichment_mysig)-8,ymax=nrow(svenrichment_mysig)-8,xmin=yposition1,xmax=yposition1) +
    annotation_custom(textGrob(paste(svenrichment_mysig[10,c("n")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(svenrichment_mysig)-9,ymax=nrow(svenrichment_mysig)-9,xmin=yposition1,xmax=yposition1) +
    annotation_custom(textGrob(paste(svenrichment_mysig[11,c("n")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(svenrichment_mysig)-10,ymax=nrow(svenrichment_mysig)-10,xmin=yposition1,xmax=yposition1) +
    annotation_custom(textGrob(paste(svenrichment_mysig[12,c("n")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(svenrichment_mysig)-11,ymax=nrow(svenrichment_mysig)-11,xmin=yposition1,xmax=yposition1) +
    annotation_custom(textGrob(paste(svenrichment_mysig[13,c("n")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(svenrichment_mysig)-12,ymax=nrow(svenrichment_mysig)-12,xmin=yposition1,xmax=yposition1) +
    annotation_custom(textGrob(paste(svenrichment_mysig[14,c("n")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(svenrichment_mysig)-13,ymax=nrow(svenrichment_mysig)-13,xmin=yposition1,xmax=yposition1) +
    
    
    ggtitle(plottitle)+
    scale_x_continuous(breaks = c(x.breaks),
                       labels = c(x.breaks))+
    coord_fixed(xlim = c(yposition1*7, 300*1.1),
                ratio=60)+
    xlab("-log10(adj. p value)") + 
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 6,colour = "black"),
      axis.text.y = element_text(size = 6,colour = "black", angle=0,hjust=0,vjust = 0.5),
      axis.text.x = element_text(size = 6,colour = "black",hjust=0.5,vjust = 0.5),
      # axis.ticks.length.y = unit(1.8, "cm"),
      axis.ticks.length.y = unit(0, "mm"),
      axis.ticks.length.x = unit(0.5, "mm"),
      axis.ticks.x = element_line(colour = "black"),
      axis.ticks.y = element_blank(),
      
      legend.position = "bottom",
      legend.justification = "right",
      legend.direction = "vertical",
      legend.margin=margin(-10,0,0,0),
      legend.box.margin=margin(0,0,0,0),
      legend.title = element_blank(),
      legend.text = element_text(size = 6,colour = "black",margin = margin(0,0,0,0)),
      legend.key = element_rect(fill = "transparent"),
      legend.key.size = unit(0.2, 'cm'),
      legend.spacing.y = unit(0, 'cm'),
      legend.spacing.x = unit(0, 'cm'),
      strip.text = element_text(size = 6),
      plot.title = element_text(size = 6,colour = "black",hjust=1,margin = margin(0,0,0,0)),
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = unit(c(0,0,0,0), "mm")
    )+
    guides(fill = guide_legend(nrow = 3, byrow = TRUE,border = "black",box.lwd = 0))
  ggsave(file.path(plot_path,paste0(plottitle,".pdf")),width =1.25,height = 1.75)
  
}




############################ RUN ####################
enrichment_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv_enrichment"
plot_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv_enrichment/plot"
if(!dir.exists(plot_path)){
  dir.create(plot_path,showWarnings = F,recursive = T)
}


hartwig_path <- file.path(enrichment_path,"hartwig_sv_enrichment_across_expression.csv")
hartwig <- read.csv(hartwig_path)
svenrichment_plot_function(svenrichment=hartwig,myexpression = "match1.vs.notgene",plottitle="hartwig.match1.vs.notgene")
svenrichment_plot_function(svenrichment=hartwig,myexpression = "match2.vs.notgene",plottitle="hartwig.match2.vs.notgene")
svenrichment_plot_function(svenrichment=hartwig,myexpression = "match3.vs.notgene",plottitle="hartwig.match3.vs.notgene")
svenrichment_plot_function(svenrichment=hartwig,myexpression = "match4.vs.notgene",plottitle="hartwig.match4.vs.notgene")
svenrichment_plot_function(svenrichment=hartwig,myexpression = "match12.vs.notgene",plottitle="hartwig.match12.vs.notgene")
svenrichment_plot_function(svenrichment=hartwig,myexpression = "match34.vs.notgene",plottitle="hartwig.match34.vs.notgene")
svenrichment_plot_function(svenrichment=hartwig,myexpression = "match234.vs.notgene",plottitle="hartwig.match234.vs.notgene")
svenrichment_plot_function(svenrichment=hartwig,myexpression = "match1234.vs.notgene",plottitle="hartwig.match1234.vs.notgene")
svenrichment_plot_function(svenrichment=hartwig,myexpression = "gene.vs.notgene",plottitle="hartwig.gene.vs.notgene")
svenrichment_plot_function(svenrichment=hartwig,myexpression = "normaltcga1.vs.notgene",plottitle="hartwig.normaltcga1.vs.notgene")
svenrichment_plot_function(svenrichment=hartwig,myexpression = "normaltcga2.vs.notgene",plottitle="hartwig.normaltcga2.vs.notgene")
svenrichment_plot_function(svenrichment=hartwig,myexpression = "normaltcga3.vs.notgene",plottitle="hartwig.normaltcga3.vs.notgene")
svenrichment_plot_function(svenrichment=hartwig,myexpression = "normaltcga4.vs.notgene",plottitle="hartwig.normaltcga4.vs.notgene")
svenrichment_plot_function(svenrichment=hartwig,myexpression = "normaltcga12.vs.notgene",plottitle="hartwig.normaltcga12.vs.notgene")
svenrichment_plot_function(svenrichment=hartwig,myexpression = "normaltcga34.vs.notgene",plottitle="hartwig.normaltcga34.vs.notgene")
svenrichment_plot_function(svenrichment=hartwig,myexpression = "normaltcga234.vs.notgene",plottitle="hartwig.normaltcga234.vs.notgene")
svenrichment_plot_function(svenrichment=hartwig,myexpression = "normaltcga1234.vs.notgene",plottitle="hartwig.normaltcga1234.vs.notgene")

pcawg_path <- file.path(enrichment_path,"pcawg_sv_enrichment_across_expression.csv")
pcawg <- read.csv(pcawg_path)
svenrichment_plot_function(svenrichment=pcawg,myexpression = "match1.vs.notgene",plottitle="pcawg.match1.vs.notgene")
svenrichment_plot_function(svenrichment=pcawg,myexpression = "match2.vs.notgene",plottitle="pcawg.match2.vs.notgene")
svenrichment_plot_function(svenrichment=pcawg,myexpression = "match3.vs.notgene",plottitle="pcawg.match3.vs.notgene")
svenrichment_plot_function(svenrichment=pcawg,myexpression = "match4.vs.notgene",plottitle="pcawg.match4.vs.notgene")
svenrichment_plot_function(svenrichment=pcawg,myexpression = "match12.vs.notgene",plottitle="pcawg.match12.vs.notgene")
svenrichment_plot_function(svenrichment=pcawg,myexpression = "match34.vs.notgene",plottitle="pcawg.match34.vs.notgene")
svenrichment_plot_function(svenrichment=pcawg,myexpression = "match234.vs.notgene",plottitle="pcawg.match234.vs.notgene")
svenrichment_plot_function(svenrichment=pcawg,myexpression = "match1234.vs.notgene",plottitle="pcawg.match1234.vs.notgene")
svenrichment_plot_function(svenrichment=pcawg,myexpression = "gene.vs.notgene",plottitle="pcawg.gene.vs.notgene")
svenrichment_plot_function(svenrichment=pcawg,myexpression = "normaltcga1.vs.notgene",plottitle="pcawg.normaltcga1.vs.notgene")
svenrichment_plot_function(svenrichment=pcawg,myexpression = "normaltcga2.vs.notgene",plottitle="pcawg.normaltcga2.vs.notgene")
svenrichment_plot_function(svenrichment=pcawg,myexpression = "normaltcga3.vs.notgene",plottitle="pcawg.normaltcga3.vs.notgene")
svenrichment_plot_function(svenrichment=pcawg,myexpression = "normaltcga4.vs.notgene",plottitle="pcawg.normaltcga4.vs.notgene")
svenrichment_plot_function(svenrichment=pcawg,myexpression = "normaltcga12.vs.notgene",plottitle="pcawg.normaltcga12.vs.notgene")
svenrichment_plot_function(svenrichment=pcawg,myexpression = "normaltcga34.vs.notgene",plottitle="pcawg.normaltcga34.vs.notgene")
svenrichment_plot_function(svenrichment=pcawg,myexpression = "normaltcga234.vs.notgene",plottitle="pcawg.normaltcga234.vs.notgene")
svenrichment_plot_function(svenrichment=pcawg,myexpression = "normaltcga1234.vs.notgene",plottitle="pcawg.normaltcga1234.vs.notgene")

