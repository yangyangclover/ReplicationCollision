### library
library(dplyr)
library(ggplot2)
library(ggpubr)
library(grid)

### functions
plot_function <- function(plot_df,
                          titlename,
                          corrected=TRUE,
                          savepath){
  if (nrow(plot_df)==0){
    out <- NA
  } else {
    
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
    
    plot_df$signaturenew <- signamenew[as.character(plot_df$signature)]
    
    # adjust p
    if(corrected==FALSE){
      plot_df$q <-  plot_df$raw.p
    } else  {
      plot_df$q <-  p.adjust(plot_df$raw.p,method = "bonferroni",n=nrow(plot_df))
    }
    
    print(nrow(plot_df))
    print(titlename)
    # log10(p)
    plot_df$logq <- ifelse(plot_df$or>1,-log10(plot_df$q),log10(plot_df$q))
    # color
    plot_df$color <- "grey"
    plot_df$color <-  ifelse(plot_df$logq > -log10(0.05),"red",plot_df$color)
    plot_df$color <-  ifelse(plot_df$logq < log10(0.05),"blue",plot_df$color)
    plot_df$color <- factor(plot_df$color, levels = c("red", "blue", "grey"))
    # abs p
    plot_df$logq <- abs(plot_df$logq)
    #  round p
    plot_df$logq <-  round(plot_df$logq,2)
    plot_df$p.text <- ifelse(plot_df$logq>=-log10(0.05),plot_df$logq,"")
    plot_df$p.text <- ifelse(is.na(plot_df$p.text),"",plot_df$p.text)
    plot_df$p.max <- ifelse(plot_df$logq>=8,8,plot_df$logq)
    plot_df$p.max <- ifelse(is.na(plot_df$p.max),0,plot_df$p.max)
    #  round or
    plot_df$or <-  round(plot_df$or,2)
    # order chrss
    plot_df$cluster <- plot_df$signaturenew
    plot_df$cluster <- factor(plot_df$cluster,levels=rev(signamenew))
    # sum a+b
    plot_df$ab <- plot_df$a + plot_df$b
    # nunmber position
    yposition1 <- -1.5
    yposition2 <- -0.8
    #  maxy
    maxy  <-  21
    maxlogq <- ceiling(max(plot_df$logq))
    maxlogq <- ifelse(maxlogq>8,maxlogq,8)
    # num  font size
    num.font <- 6
    num.font2 <- 6
    # other value
    hjust.text <- 1
    # label
    descip <- as.character(unique(plot_df$test))
    if (descip %in% c("rep.vs.unzip")){
      red.label <- "replicated strand"
      blue.label <- "unreplicated strand"
    }
    if (descip %in% c("rep.vs.unzip.oppo")){
      red.label <- "replicated strand"
      blue.label <- "unreplicated strand"
    }
    if (descip %in% c("rep.vs.unzip.same")){
      red.label <- "replicated strand"
      blue.label <- "unreplicated strand"
    }
    
    
    if (nrow(plot_df)==13){
      out=ggplot(data = plot_df, aes(y = cluster, x = ifelse(logq>maxy,maxy,logq), fill = color)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = c("red", "blue", "grey"),
                          breaks = c("red", "blue", "grey"),
                          labels = c(paste0("Enriched in ",red.label),
                                     paste0("Enriched in ",blue.label),
                                     "Not significant"),
                          drop=FALSE) +
        scale_x_continuous(breaks = seq(0,maxlogq,2),
                           labels = c(seq(0,maxlogq,2))
        )+
        xlab("-log10(adj. p value)") +
        ylab("")+
        ggtitle(titlename) +
        coord_fixed(xlim = c(-4, maxlogq),clip = "off",ratio = 2)+
        geom_vline(xintercept=-log10(0.05), linetype="dashed", color = "black")+
        # # geom_vline(xintercept=8, linetype="solid", color = "black")+
        annotation_custom(textGrob("n", gp=gpar(fontsize=num.font, fontface="bold"), rot = 0,vjust = 0,hjust = hjust.text),ymin=-.25,ymax=-.25,xmin=-3,xmax=-3) +
        annotation_custom(textGrob(paste(plot_df[1,c("ab")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(plot_df),ymax=nrow(plot_df),xmin=yposition1,xmax=yposition1) +
        annotation_custom(textGrob(paste(plot_df[2,c("ab")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(plot_df)-1,ymax=nrow(plot_df)-1,xmin=yposition1,xmax=yposition1) +
        annotation_custom(textGrob(paste(plot_df[3,c("ab")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(plot_df)-2,ymax=nrow(plot_df)-2,xmin=yposition1,xmax=yposition1) +
        annotation_custom(textGrob(paste(plot_df[4,c("ab")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(plot_df)-3,ymax=nrow(plot_df)-3,xmin=yposition1,xmax=yposition1) +
        annotation_custom(textGrob(paste(plot_df[5,c("ab")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(plot_df)-4,ymax=nrow(plot_df)-4,xmin=yposition1,xmax=yposition1) +
        annotation_custom(textGrob(paste(plot_df[6,c("ab")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(plot_df)-5,ymax=nrow(plot_df)-5,xmin=yposition1,xmax=yposition1) +
        annotation_custom(textGrob(paste(plot_df[7,c("ab")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(plot_df)-6,ymax=nrow(plot_df)-6,xmin=yposition1,xmax=yposition1) +
        annotation_custom(textGrob(paste(plot_df[8,c("ab")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(plot_df)-7,ymax=nrow(plot_df)-7,xmin=yposition1,xmax=yposition1) +
        annotation_custom(textGrob(paste(plot_df[9,c("ab")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(plot_df)-8,ymax=nrow(plot_df)-8,xmin=yposition1,xmax=yposition1) +
        annotation_custom(textGrob(paste(plot_df[10,c("ab")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(plot_df)-9,ymax=nrow(plot_df)-9,xmin=yposition1,xmax=yposition1) +
        annotation_custom(textGrob(paste(plot_df[11,c("ab")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(plot_df)-10,ymax=nrow(plot_df)-10,xmin=yposition1,xmax=yposition1) +
        annotation_custom(textGrob(paste(plot_df[13,c("ab")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(plot_df)-11,ymax=nrow(plot_df)-11,xmin=yposition1,xmax=yposition1) +
        annotation_custom(textGrob(paste(plot_df[12,c("ab")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(plot_df)-12,ymax=nrow(plot_df)-12,xmin=yposition1,xmax=yposition1) +
        
        theme(
          axis.title.x = element_text(size = num.font,colour = "black"),
          axis.title.y = element_text(size = num.font,colour = "black"),
          axis.text.x = element_text(size = num.font,colour = "black"),
          axis.text.y = element_text(size = num.font,colour = "black",hjust=0,vjust = 0.5),
          # axis.ticks.length.y = unit(1.8, "cm"),
          axis.ticks.length.x = unit(0.5, "mm"),
          axis.ticks.length.y = unit(0, "mm"),
          axis.ticks.y = element_blank(),
          
          legend.position = "bottom",
          legend.justification = "right",
          legend.direction = "vertical",
          legend.margin=margin(-10,0,0,0),
          legend.box.margin=margin(0,0,0,0),
          legend.title = element_blank(),
          legend.text = element_text(size = num.font,colour = "black",margin = margin(0,0,0,0)),
          legend.key = element_rect(fill = "transparent"),
          legend.key.size = unit(0.2, 'cm'),
          legend.spacing.y = unit(0, 'cm'),
          legend.spacing.x = unit(0, 'cm'),
          strip.text = element_text(size = num.font),
          plot.title = element_text(size = num.font,colour = "black",hjust=0.5,margin = margin(0,0,0,0)),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(0,-10,0,-10), "mm")
          # axis.line = element_line(colour = "black")
        )+
        guides(fill = guide_legend(nrow = 3, byrow = TRUE,border = "black",box.lwd = 0))
    }
    
    if (nrow(plot_df)==14){
      out=ggplot(data = plot_df, aes(y = cluster, x = ifelse(logq>maxy,maxy,logq), fill = color)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = c("red", "blue", "grey"),
                          breaks = c("red", "blue", "grey"),
                          labels = c(paste0("Enriched in ",red.label),
                                     paste0("Enriched in ",blue.label),
                                     "Not significant"),
                          drop=FALSE) +
        scale_x_continuous(breaks = seq(0,maxlogq,2),
                           labels = c(seq(0,maxlogq,2))
        )+
        xlab("-log10(adj. p value)") +
        ylab("")+
        ggtitle(titlename) +
        coord_fixed(xlim = c(-4, maxlogq),clip = "off",ratio = 2)+
        geom_vline(xintercept=-log10(0.05), linetype="dashed", color = "black")+
        # # geom_vline(xintercept=8, linetype="solid", color = "black")+
        annotation_custom(textGrob("n", gp=gpar(fontsize=num.font, fontface="bold"), rot = 0,vjust = 0,hjust = hjust.text),ymin=-.25,ymax=-.25,xmin=-3,xmax=-3) +
        annotation_custom(textGrob(paste(plot_df[1,c("ab")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(plot_df),ymax=nrow(plot_df),xmin=yposition1,xmax=yposition1) +
        annotation_custom(textGrob(paste(plot_df[2,c("ab")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(plot_df)-1,ymax=nrow(plot_df)-1,xmin=yposition1,xmax=yposition1) +
        annotation_custom(textGrob(paste(plot_df[3,c("ab")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(plot_df)-2,ymax=nrow(plot_df)-2,xmin=yposition1,xmax=yposition1) +
        annotation_custom(textGrob(paste(plot_df[4,c("ab")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(plot_df)-3,ymax=nrow(plot_df)-3,xmin=yposition1,xmax=yposition1) +
        annotation_custom(textGrob(paste(plot_df[5,c("ab")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(plot_df)-4,ymax=nrow(plot_df)-4,xmin=yposition1,xmax=yposition1) +
        annotation_custom(textGrob(paste(plot_df[6,c("ab")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(plot_df)-5,ymax=nrow(plot_df)-5,xmin=yposition1,xmax=yposition1) +
        annotation_custom(textGrob(paste(plot_df[7,c("ab")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(plot_df)-6,ymax=nrow(plot_df)-6,xmin=yposition1,xmax=yposition1) +
        annotation_custom(textGrob(paste(plot_df[8,c("ab")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(plot_df)-7,ymax=nrow(plot_df)-7,xmin=yposition1,xmax=yposition1) +
        annotation_custom(textGrob(paste(plot_df[9,c("ab")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(plot_df)-8,ymax=nrow(plot_df)-8,xmin=yposition1,xmax=yposition1) +
        annotation_custom(textGrob(paste(plot_df[10,c("ab")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(plot_df)-9,ymax=nrow(plot_df)-9,xmin=yposition1,xmax=yposition1) +
        annotation_custom(textGrob(paste(plot_df[11,c("ab")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(plot_df)-10,ymax=nrow(plot_df)-10,xmin=yposition1,xmax=yposition1) +
        annotation_custom(textGrob(paste(plot_df[12,c("ab")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(plot_df)-11,ymax=nrow(plot_df)-11,xmin=yposition1,xmax=yposition1) +
        annotation_custom(textGrob(paste(plot_df[14,c("ab")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(plot_df)-12,ymax=nrow(plot_df)-12,xmin=yposition1,xmax=yposition1) +
        annotation_custom(textGrob(paste(plot_df[13,c("ab")], sep="\n", collapse=" "), gp=gpar(fontsize=num.font, fontface="plain"), rot = 0,hjust = hjust.text),ymin=nrow(plot_df)-13,ymax=nrow(plot_df)-13,xmin=yposition1,xmax=yposition1) +
       
        theme(
          axis.title.x = element_text(size = num.font,colour = "black"),
          axis.title.y = element_text(size = num.font,colour = "black"),
          axis.text.x = element_text(size = num.font,colour = "black"),
          axis.text.y = element_text(size = num.font,colour = "black",hjust=0,vjust = 0.5),
          # axis.ticks.length.y = unit(1.8, "cm"),
          axis.ticks.length.x = unit(0.5, "mm"),
          axis.ticks.length.y = unit(0, "mm"),
          axis.ticks.y = element_blank(),
          
          legend.position = "bottom",
          legend.justification = "right",
          legend.direction = "vertical",
          legend.margin=margin(-10,0,0,0),
          legend.box.margin=margin(0,0,0,0),
          legend.title = element_blank(),
          legend.text = element_text(size = num.font,colour = "black",margin = margin(0,0,0,0)),
          legend.key = element_rect(fill = "transparent"),
          legend.key.size = unit(0.2, 'cm'),
          legend.spacing.y = unit(0, 'cm'),
          legend.spacing.x = unit(0, 'cm'),
          strip.text = element_text(size = num.font),
          plot.title = element_text(size = num.font,colour = "black",hjust=0.5,margin = margin(0,0,0,0)),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(0,-10,0,-20), "mm")
          # axis.line = element_line(colour = "black")
        )+
        guides(fill = guide_legend(nrow = 3, byrow = TRUE,border = "black",box.lwd = 0))
    }
    
  }
  ggsave(savepath,width = ifelse(maxlogq*0.12<1.8,1.8,maxlogq*0.12),height = 2)
  return(out)
}
#############################################################
shrink <- 25000
large_dup_path <- file.path("D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup")
plot_path <- file.path(large_dup_path,"pcawg/006-features/plot/single",paste0("shrink",shrink))
if (!dir.exists(plot_path)){
  dir.create(plot_path,showWarnings = F,recursive = T)
}

pcawg_path <- file.path(large_dup_path,"pcawg/006-features/result_test")
validation_path <- file.path(large_dup_path,"validation")

input_path=file.path(pcawg_path,paste0("sv.rt.NA.consistent8",".shrink",shrink,".simple.fishertest.csv"))
input <- read.csv(input_path)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "all"),
  titlename = paste0("pcawg.consistent8.rep.vs.unrep\n(all)"),
  savepath = file.path(plot_path, paste0("pcawg.consistent8",".shrink",shrink,".allsv", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip.oppo", expression == "all"),
  titlename = paste0("pcawg.consistent8.rep.vs.unrep\n(head-on)"),
  savepath = file.path(plot_path, paste0("pcawg.consistent8",".shrink",shrink,".head-on", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip.same", expression == "all"),
  titlename = paste0("pcawg.consistent8.rep.vs.unrep\n(co-direction)"),
  savepath = file.path(plot_path, paste0("pcawg.consistent8",".shrink",shrink,".co-direction", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "match1"),
  titlename = paste0("pcawg.consistent8.rep.vs.unrep\n(match1)"),
  savepath = file.path(plot_path, paste0("pcawg.consistent8",".shrink",shrink,".match1", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "match234"),
  titlename = paste0("pcawg.consistent8.rep.vs.unrep\n(match234)"),
  savepath = file.path(plot_path, paste0("pcawg.consistent8",".shrink",shrink,".match234", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "normaltcga1"),
  titlename = paste0("pcawg.consistent8.rep.vs.unrep\n(normaltcga1)"),
  savepath = file.path(plot_path, paste0("pcawg.consistent8",".shrink",shrink,".normaltcga1", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "normaltcga234"),
  titlename = paste0("pcawg.consistent8.rep.vs.unrep\n(normaltcga234)"),
  savepath = file.path(plot_path, paste0("pcawg.consistent8",".shrink",shrink,".normaltcga234", ".pdf"))
)


input_path=file.path(pcawg_path,paste0("sv.rt.Liver-HCC.HepG2",".shrink",shrink,".simple.fishertest.csv"))
input <- read.csv(input_path)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "all"),
  titlename = paste0("PCAWG.Liver-HCC.HepG2.rep.vs.unrep\n(all)"),
  savepath = file.path(plot_path, paste0("pcawg.Liver-HCC.HepG2",".shrink",shrink,".allsv", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "match1"),
  titlename = paste0("PCAWG.Liver-HCC.HepG2.rep.vs.unrep\n(match1)"),
  savepath = file.path(plot_path, paste0("pcawg.Liver-HCC.HepG2",".shrink",shrink,".match1", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "match234"),
  titlename = paste0("PCAWG.Liver-HCC.HepG2.rep.vs.unrep\n(match234)"),
  savepath = file.path(plot_path, paste0("pcawg.Liver-HCC.HepG2",".shrink",shrink,".match234", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "normaltcga1"),
  titlename = paste0("PCAWG.Liver-HCC.HepG2.rep.vs.unrep\n(normaltcga1)"),
  savepath = file.path(plot_path, paste0("pcawg.Liver-HCC.HepG2",".shrink",shrink,".normaltcga1", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "normaltcga234"),
  titlename = paste0("PCAWG.Liver-HCC.HepG2.rep.vs.unrep\n(normaltcga234)"),
  savepath = file.path(plot_path, paste0("pcawg.Liver-HCC.HepG2",".shrink",shrink,".normaltcga234", ".pdf"))
)

input_path=file.path(pcawg_path,paste0("sv.rt.Ovary-AdenoCA.consistent8",".shrink",shrink,".simple.fishertest.csv"))
input <- read.csv(input_path)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "all"),
  titlename = paste0("PCAWG.Ovary-AdenoCA.consistent8.rep.vs.unrep\n(all)"),
  savepath = file.path(plot_path, paste0("pcawg.Ovary-AdenoCA.consistent8",".shrink",shrink,".allsv", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "match1"),
  titlename = paste0("PCAWG.Ovary-AdenoCA.consistent8.rep.vs.unrep\n(match1)"),
  savepath = file.path(plot_path, paste0("pcawg.Ovary-AdenoCA.consistent8",".shrink",shrink,".match1", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "match234"),
  titlename = paste0("PCAWG.Ovary-AdenoCA.consistent8.rep.vs.unrep\n(match234)"),
  savepath = file.path(plot_path, paste0("pcawg.Ovary-AdenoCA.consistent8",".shrink",shrink,".match234", ".pdf"))
)


input_path=file.path(pcawg_path,paste0("sv.rt.Breast-AdenoCA.MCF7",".shrink",shrink,".simple.fishertest.csv"))
input <- read.csv(input_path)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "all"),
  titlename = paste0("PCAWG.Breast-AdenoCA.MCF7.rep.vs.unrep\n(all)"),
  savepath = file.path(plot_path, paste0("pcawg.Breast-AdenoCA.MCF7",".shrink",shrink,".allsv", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "match1"),
  titlename = paste0("PCAWG.Breast-AdenoCA.MCF7.rep.vs.unrep\n(match1)"),
  savepath = file.path(plot_path, paste0("pcawg.Breast-AdenoCA.MCF7",".shrink",shrink,".match1", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "match234"),
  titlename = paste0("PCAWG.Breast-AdenoCA.MCF7.rep.vs.unrep\n(match234)"),
  savepath = file.path(plot_path, paste0("pcawg.Breast-AdenoCA.MCF7",".shrink",shrink,".match234", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "normaltcga1"),
  titlename = paste0("PCAWG.Breast-AdenoCA.MCF7.rep.vs.unrep\n(normaltcga1)"),
  savepath = file.path(plot_path, paste0("pcawg.Breast-AdenoCA.MCF7",".shrink",shrink,".normaltcga1", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "normaltcga234"),
  titlename = paste0("PCAWG.Breast-AdenoCA.MCF7.rep.vs.unrep\n(normaltcga234)"),
  savepath = file.path(plot_path, paste0("pcawg.Breast-AdenoCA.MCF7",".shrink",shrink,".normaltcga234", ".pdf"))
)
input_path=file.path(pcawg_path,paste0("sv.rt.Ovary-AdenoCA.consistent8",".shrink",shrink,".simple.fishertest.csv"))
input <- read.csv(input_path)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "all"),
  titlename = paste0("PCAWG.Ovary-AdenoCA.consistent8.rep.vs.unrep\n(all)"),
  savepath = file.path(plot_path, paste0("pcawg.Ovary-AdenoCA.consistent8",".shrink",shrink,".allsv", ".pdf"))
)


input_path=file.path(pcawg_path,paste0("sv.rt.Stomach-AdenoCA.consistent8",".shrink",shrink,".simple.fishertest.csv"))
input <- read.csv(input_path)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "all"),
  titlename = paste0("PCAWG.Stomach-AdenoCA.consistent8.rep.vs.unrep\n(all)"),
  savepath = file.path(plot_path, paste0("pcawg.Stomach-AdenoCA.consistent8",".shrink",shrink,".allsv", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "normaltcga1"),
  titlename = paste0("PCAWG.Stomach-AdenoCA.consistent8.rep.vs.unrep\n(normaltcga1)"),
  savepath = file.path(plot_path, paste0("pcawg.Stomach-AdenoCA.consistent8",".shrink",shrink,".normaltcga1", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "normaltcga234"),
  titlename = paste0("PCAWG.Stomach-AdenoCA.consistent8.rep.vs.unrep\n(normaltcga234)"),
  savepath = file.path(plot_path, paste0("pcawg.Stomach-AdenoCA.consistent8",".shrink",shrink,".normaltcga234", ".pdf"))
)

input_path=file.path(pcawg_path,paste0("sv.rt.Prost-AdenoCA.consistent8",".shrink",shrink,".simple.fishertest.csv"))
input <- read.csv(input_path)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "all"),
  titlename = paste0("PCAWG.Prost-AdenoCA.consistent8.rep.vs.unrep\n(all)"),
  savepath = file.path(plot_path, paste0("pcawg.Prost-AdenoCA.consistent8",".shrink",shrink,".allsv", ".pdf"))
)


input_path=file.path(validation_path,"PG570","genomic_features/result_test",paste0("sv.rt.NA.consistent8",".shrink",shrink,".fishertest.csv"))
input <- read.csv(input_path)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "all"),
  titlename = paste0("PG570.consistent8.rep.vs.unrep\n(all)"),
  savepath = file.path(plot_path, paste0("POG570.consistent8",".shrink",shrink,".allsv", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "match1"),
  titlename = paste0("PG570.consistent8.rep.vs.unrep\n(match1)"),
  savepath = file.path(plot_path, paste0("POG570.consistent8",".shrink",shrink,".match1", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "match234"),
  titlename = paste0("PG570.consistent8.rep.vs.unrep\n(match234)"),
  savepath = file.path(plot_path, paste0("POG570.consistent8",".shrink",shrink,".match234", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "normaltcga1"),
  titlename = paste0("PG570.consistent8.rep.vs.unrep\n(normaltcga1)"),
  savepath = file.path(plot_path, paste0("POG570.consistent8",".shrink",shrink,".normaltcga1", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "normaltcga234"),
  titlename = paste0("PG570.consistent8.rep.vs.unrep\n(normaltcga234)"),
  savepath = file.path(plot_path, paste0("POG570.consistent8.normaltcga234", ".pdf"))
)

input_path=file.path(validation_path,"PG570","genomic_features/result_test",paste0("sv.rt.BRCA.MCF7",".shrink",shrink,".fishertest.csv"))
input <- read.csv(input_path)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "all"),
  titlename = paste0("PG570.BRCA.MCF7.rep.vs.unrep\n(all)"),
  savepath = file.path(plot_path, paste0("POG570.BRCA.MCF7",".shrink",shrink,".allsv", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "match1"),
  titlename = paste0("PG570.BRCA.MCF7.rep.vs.unrep\n(match1)"),
  savepath = file.path(plot_path, paste0("POG570.BRCA.MCF7",".shrink",shrink,".match1", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "match234"),
  titlename = paste0("PG570.BRCA.MCF7.rep.vs.unrep\n(match234)"),
  savepath = file.path(plot_path, paste0("POG570.BRCA.MCF7",".shrink",shrink,".match234", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "normaltcga1"),
  titlename = paste0("PG570.BRCA.MCF7.rep.vs.unrep\n(normaltcga1)"),
  savepath = file.path(plot_path, paste0("POG570.BRCA.MCF7",".shrink",shrink,".normaltcga1", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "normaltcga234"),
  titlename = paste0("PG570.BRCA.MCF7.rep.vs.unrep\n(normaltcga234)"),
  savepath = file.path(plot_path, paste0("POG570.BRCA.MCF7",".shrink",shrink,".normaltcga234", ".pdf"))
)

input_path=file.path(validation_path,"Hartwig","genomic_features/result_test",paste0("sv.rt.NA.consistent8",".shrink",shrink,".fishertest.csv"))
input <- read.csv(input_path)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "all"),
  titlename = paste0("Hartwig.consistent8.rep.vs.unrep\n(all)"),
  savepath = file.path(plot_path, paste0("Hartwig.consistent8",".shrink",shrink,".allsv", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip.oppo", expression == "all"),
  titlename = paste0("Hartwig.consistent8.rep.vs.unrep\n(head-on)"),
  savepath = file.path(plot_path, paste0("Hartwig.consistent8",".shrink",shrink,".head-on", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip.same", expression == "all"),
  titlename = paste0("Hartwig.consistent8.rep.vs.unrep\n(co-direction)"),
  savepath = file.path(plot_path, paste0("Hartwig.consistent8",".shrink",shrink,".co-direction", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "match1"),
  titlename = paste0("Hartwig.consistent8.rep.vs.unrep\n(match1)"),
  savepath = file.path(plot_path, paste0("Hartwig.consistent8",".shrink",shrink,".match1", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "match234"),
  titlename = paste0("Hartwig.consistent8.rep.vs.unrep\n(match234)"),
  savepath = file.path(plot_path, paste0("Hartwig.consistent8",".shrink",shrink,".match234", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "normaltcga1"),
  titlename = paste0("Hartwig.consistent8.rep.vs.unrep\n(normaltcga1)"),
  savepath = file.path(plot_path, paste0("Hartwig.consistent8",".shrink",shrink,".normaltcga1", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "normaltcga234"),
  titlename = paste0("Hartwig.consistent8.rep.vs.unrep\n(normaltcga234)"),
  savepath = file.path(plot_path, paste0("Hartwig.consistent8",".shrink",shrink,".normaltcga234", ".pdf"))
)


input_path=file.path(validation_path,"Hartwig","genomic_features/result_test",paste0("sv.rt.Breast.MCF7",".shrink",shrink,".fishertest.csv"))
input <- read.csv(input_path)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "all"),
  titlename = paste0("Hartwig.Breast.MCF7.rep.vs.unrep\n(all)"),
  savepath = file.path(plot_path, paste0("Hartwig.Breast.MCF7",".shrink",shrink,".allsv", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip.oppo", expression == "all"),
  titlename = paste0("Hartwig.Breast.MCF7.rep.vs.unrep\n(head-on)"),
  savepath = file.path(plot_path, paste0("Hartwig.Breast.MCF7",".shrink",shrink,".head-on", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip.same", expression == "all"),
  titlename = paste0("Hartwig.Breast.MCF7.rep.vs.unrep\n(co-direction)"),
  savepath = file.path(plot_path, paste0("Hartwig.Breast.MCF7",".shrink",shrink,".co-direction", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "match1"),
  titlename = paste0("Hartwig.Breast.MCF7.rep.vs.unrep\n(match1)"),
  savepath = file.path(plot_path, paste0("Hartwig.Breast.MCF7",".shrink",shrink,".match1", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "match234"),
  titlename = paste0("Hartwig.Breast.MCF7.rep.vs.unrep\n(match234)"),
  savepath = file.path(plot_path, paste0("Hartwig.Breast.MCF7",".shrink",shrink,".match234", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "normaltcga1"),
  titlename = paste0("Hartwig.Breast.MCF7.rep.vs.unrep\n(normaltcga1)"),
  savepath = file.path(plot_path, paste0("Hartwig.Breast.MCF7",".shrink",shrink,".normaltcga1", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "normaltcga234"),
  titlename = paste0("Hartwig.Breast.MCF7.rep.vs.unrep\n(normaltcga234)"),
  savepath = file.path(plot_path, paste0("Hartwig.Breast.MCF7",".shrink",shrink,".normaltcga234", ".pdf"))
)

input_path=file.path(validation_path,"Hartwig","genomic_features/result_test",paste0("sv.rt.Ovary.consistent8",".shrink",shrink,".fishertest.csv"))
input <- read.csv(input_path)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "all"),
  titlename = paste0("Hartwig.Ovary.consistent8.rep.vs.unrep\n(all)"),
  savepath = file.path(plot_path, paste0("Hartwig.Ovary.consistent8",".shrink",shrink,".allsv", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "match1"),
  titlename = paste0("Hartwig.Ovary.consistent8.rep.vs.unrep\n(match1)"),
  savepath = file.path(plot_path, paste0("Hartwig.Ovary.consistent8",".shrink",shrink,".match1", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "match234"),
  titlename = paste0("Hartwig.Ovary.consistent8.rep.vs.unrep\n(match234)"),
  savepath = file.path(plot_path, paste0("Hartwig.Ovary.consistent8",".shrink",shrink,".match234", ".pdf"))
)

input_path=file.path(validation_path,"Hartwig","genomic_features/result_test",paste0("sv.rt.Esophagus.consistent8",".shrink",shrink,".fishertest.csv"))
input <- read.csv(input_path)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "all"),
  titlename = paste0("Hartwig.Esophagus.consistent8.rep.vs.unrep\n(all)"),
  savepath = file.path(plot_path, paste0("Hartwig.Esophagus.consistent8",".shrink",shrink,".allsv", ".pdf"))
)

input_path=file.path(validation_path,"Hartwig","genomic_features/result_test",paste0("sv.rt.Prostate.consistent8",".shrink",shrink,".fishertest.csv"))
input <- read.csv(input_path)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "all"),
  titlename = paste0("Hartwig.Prostate.consistent8.rep.vs.unrep\n(all)"),
  savepath = file.path(plot_path, paste0("Hartwig.Prostate.consistent8",".shrink",shrink,".allsv", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip.oppo", expression == "all"),
  titlename = paste0("Hartwig.Prostate.consistent8.rep.vs.unrep\n(head-on)"),
  savepath = file.path(plot_path, paste0("Hartwig.Prostate.consistent8",".shrink",shrink,".head-on", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip.same", expression == "all"),
  titlename = paste0("Hartwig.Prostate.consistent8.rep.vs.unrep\n(co-direction)"),
  savepath = file.path(plot_path, paste0("Hartwig.Prostate.consistent8",".shrink",shrink,".co-direction", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "match1"),
  titlename = paste0("Hartwig.Prostate.consistent8.rep.vs.unrep\n(match1)"),
  savepath = file.path(plot_path, paste0("Hartwig.Prostate.consistent8",".shrink",shrink,".match1", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "match234"),
  titlename = paste0("Hartwig.Prostate.consistent8.rep.vs.unrep\n(match234)"),
  savepath = file.path(plot_path, paste0("Hartwig.Prostate.consistent8",".shrink",shrink,".match234", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "normaltcga1"),
  titlename = paste0("Hartwig.Prostate.consistent8.rep.vs.unrep\n(normaltcga1)"),
  savepath = file.path(plot_path, paste0("Hartwig.Prostate.consistent8",".shrink",shrink,".normaltcga1", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "normaltcga234"),
  titlename = paste0("Hartwig.Prostate.consistent8.rep.vs.unrep\n(normaltcga234)"),
  savepath = file.path(plot_path, paste0("Hartwig.Prostate.consistent8",".shrink",shrink,".normaltcga234", ".pdf"))
)

# input_path=file.path(validation_path,"Hartwig","genomic_features/result_test","sv.rt.Liver.HepG2.fishertest.csv")
# input <- read.csv(input_path)
# plot_function(
#   plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "all"),
#   titlename = paste0("Hartwig.Liver.HepG2.rep.vs.unrep\n(all)"),
#   savepath = file.path(plot_path, paste0("Hartwig.Liver.HepG2.allsv", ".pdf"))
# )

input_path=file.path(validation_path,"BRCA-EU","genomic_features/result_test",paste0("sv.rt.NA.MCF7",".shrink",shrink,".fishertest.csv"))
input <- read.csv(input_path)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "all"),
  titlename = paste0("BRCA-EU.MCF7.rep.vs.unrep\n(all)"),
  savepath = file.path(plot_path, paste0("BRCA-EU.MCF7",".shrink",shrink,".allsv", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip.oppo", expression == "all"),
  titlename = paste0("BRCA-EU.MCF7.rep.vs.unrep\n(head-on)"),
  savepath = file.path(plot_path, paste0("BRCA-EU.MCF7",".shrink",shrink,".head-on", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip.same", expression == "all"),
  titlename = paste0("BRCA-EU.MCF7.rep.vs.unrep\n(co-direction)"),
  savepath = file.path(plot_path, paste0("BRCA-EU.MCF7",".shrink",shrink,".co-direction", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "match1"),
  titlename = paste0("BRCA-EU.MCF7.rep.vs.unrep\n(match1)"),
  savepath = file.path(plot_path, paste0("BRCA-EU.MCF7",".shrink",shrink,".match1", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "match234"),
  titlename = paste0("BRCA-EU.MCF7.rep.vs.unrep\n(match234)"),
  savepath = file.path(plot_path, paste0("BRCA-EU.MCF7",".shrink",shrink,".match234", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "normaltcga1"),
  titlename = paste0("BRCA-EU.MCF7.rep.vs.unrep\n(normaltcga1)"),
  savepath = file.path(plot_path, paste0("BRCA-EU.MCF7",".shrink",shrink,".normaltcga1", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "normaltcga234"),
  titlename = paste0("BRCA-EU.MCF7.rep.vs.unrep\n(normaltcga234)"),
  savepath = file.path(plot_path, paste0("BRCA-EU.MCF7",".shrink",shrink,".normaltcga234", ".pdf"))
)

input_path=file.path(validation_path,"OV-AU","genomic_features/result_test",paste0("sv.rt.NA.consistent8",".shrink",shrink,".fishertest.csv"))
input <- read.csv(input_path)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "all"),
  titlename = paste0("OV-AU.consistent8.rep.vs.unrep\n(all)"),
  savepath = file.path(plot_path, paste0("OV-AU.consistent8",".shrink",shrink,".allsv", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "match1"),
  titlename = paste0("OV-AU.consistent8.rep.vs.unrep\n(match1)"),
  savepath = file.path(plot_path, paste0("OV-AU.consistent8",".shrink",shrink,".match1", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "match234"),
  titlename = paste0("OV-AU.consistent8.rep.vs.unrep\n(match234)"),
  savepath = file.path(plot_path, paste0("OV-AU.consistent8",".shrink",shrink,".match234", ".pdf"))
)

input_path=file.path(validation_path,"PRAD-CA","genomic_features/result_test",paste0("sv.rt.NA.consistent8",".shrink",shrink,".fishertest.csv"))
input <- read.csv(input_path)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "all"),
  titlename = paste0("PRAD-CA.consistent8.rep.vs.unrep\n(all)"),
  savepath = file.path(plot_path, paste0("PRAD-CA.consistent8",".shrink",shrink,".allsv", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "normaltcga1"),
  titlename = paste0("PRAD-CA.consistent8.rep.vs.unrep\n(normaltcga1)"),
  savepath = file.path(plot_path, paste0("PRAD-CA.consistent8",".shrink",shrink,".normaltcga1", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "normaltcga234"),
  titlename = paste0("PRAD-CA.consistent8.rep.vs.unrep\n(normaltcga234)"),
  savepath = file.path(plot_path, paste0("PRAD-CA.consistent8",".shrink",shrink,".normaltcga234", ".pdf"))
)

input_path=file.path(validation_path,"PRAD-UK","genomic_features/result_test",paste0("sv.rt.NA.consistent8",".shrink",shrink,".fishertest.csv"))
input <- read.csv(input_path)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "all"),
  titlename = paste0("PRAD-UK.consistent8.rep.vs.unrep\n(all)"),
  savepath = file.path(plot_path, paste0("PRAD-UK.consistent8",".shrink",shrink,".allsv", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "normaltcga1"),
  titlename = paste0("PRAD-UK.consistent8.rep.vs.unrep\n(normaltcga1)"),
  savepath = file.path(plot_path, paste0("PRAD-UK.consistent8",".shrink",shrink,".normaltcga1", ".pdf"))
)
plot_function(
  plot_df = input %>% dplyr::filter(test == "rep.vs.unzip", expression == "normaltcga234"),
  titlename = paste0("PRAD-UK.consistent8.rep.vs.unrep\n(normaltcga234)"),
  savepath = file.path(plot_path, paste0("PRAD-UK.consistent8",".shrink",shrink,".normaltcga234", ".pdf"))
)
