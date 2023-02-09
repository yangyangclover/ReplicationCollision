# library
library(preprocessCore)
library(gdata)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(stringr)

rt_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/genomic_features/rt_cellline"
direction_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/genomic_features/rt_cellline/direction"
summary_path <- file.path(rt_path,"rt_cell_line_summary.csv")
pub_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/rptiming"

if (!exists(direction_path)){
  dir.create(direction_path,showWarnings = F,recursive = T)
}

### read summary
mychr <- paste0("chr",c(seq(1,22,1),"X","Y"))
summary <- read.csv(summary_path)


########### functions ################
### quantile normalization
add_raw_rt_function <- function(name){
  mychr <- paste0("chr",c(seq(1,22,1),"X","Y"))
  file_name <- as.character(summary[(summary$name == name)  %in% "TRUE", "file"])
  file <- read.table(file.path(rt_path, file_name))
  file <- file[file$V1 %in% mychr,]
  normalize.df <- if(exists("normalize.df"))  {
    cbindX(normalize.df,data.frame(file$V4))
  } else  {
    data.frame(file$V4)
  }
  names(normalize.df)[ncol(normalize.df)] <- name
  return(normalize.df)
}



# add normalized  value to orginal  file
rt_direction_function <- function(name = "HEK293_final",
                                  range_size = 50000,
                                  outname = "HEK293_final_direction_quantile_norm_18cellline_range50000"
                                  ){
  
  file_name <- as.character(summary[(summary$name == name)  %in% "TRUE", "file"])
  mychr <- paste0("chr",c(seq(1,22,1),"X","Y"))
  file <- read.table(file.path(rt_path, file_name))
  file <- file[file$V1 %in% mychr,]
  file[,"normalized_rt"] <- normalize.mat_normalized[,name][!is.na(normalize.mat_normalized[,name])]
  file$mid <-  floor((file$V2 + file$V3)/2)
  file <- file[order(file$V1,file$mid),]
  
  file_simple <- file
  names(file_simple) <- c("chr","start","end","rawrt","rt1","mid")
  file_simple$group <- floor(file_simple$mid/range_size)
  write.table(file_simple[,c(1,2,3,5)],file.path(direction_path,paste0(gsub("direction","rt",outname),".txt")),quote = F,sep = "\t",row.names = F,col.names = T)
  
  file$group <- floor(file$mid/range_size)
  

  peak_valley_all <- data.frame()
  file_chr_new_all <- data.frame()
  for (chri in unique(file$V1)){ #unique(file$V1)
    print(chri)
    file_chr <- file[file$V1==chri,] #do chr by chr
    
    ###  smooth part
    # row that min mid
    min.mid.row <- which.min(file_chr$mid)
    max.mid.row <- which.max(file_chr$mid)
    boundary1.row <- data.frame(
      "V1" = chri,
      "V2" = file_chr[min.mid.row,"V2"],
      "V3" = file_chr[min.mid.row,"V3"],
      "normalized_rt" = file_chr[min.mid.row,"normalized_rt"],
      "mid" = file_chr[min.mid.row,"mid"]
    )
    boundary2.row <- data.frame(
      "V1" = chri,
      "V2" = file_chr[max.mid.row,"V2"],
      "V3" = file_chr[max.mid.row,"V3"],
      "normalized_rt" = file_chr[max.mid.row,"normalized_rt"],
      "mid" = file_chr[max.mid.row,"mid"]
    )
    v2.vector <- c()
    v3.vector <- c()
    rt.vector <- c()
    mid.vector <- c()
    for (i in unique(file_chr$group)) {
      rows <- which(file_chr$group == i)
      line1 <- head(rows, 1)
      line2 <- tail(rows, 1)
      v2.vector <- c(v2.vector, file_chr[line1, "V2"])
      v3.vector <- c(v3.vector, file_chr[line1, "V3"])
      rt.vector <-
        c(rt.vector, mean(file_chr[line1:line2, "normalized_rt"]))
      #mid.vector <- c(mid.vector,floor((file_chr[line1,"V2"]+file_chr[line2,"V3"])/2))
      mid.vector <- c(mid.vector, (i + 1) * (range_size) - range_size /
                        2)
      file_chr_new <- data.frame(
        "V1" = chri,
        "V2" = v2.vector,
        "V3" = v3.vector,
        "normalized_rt" = rt.vector,
        "mid" = mid.vector
      )
    }
    ### reorder 
    file_chr_new <- file_chr_new[order(file_chr_new$mid),]
    ### add boundary
    boundary1.row.mid <- boundary1.row[1,"mid"]
    boundary2.row.mid <- boundary2.row[1,"mid"]
    if (boundary1.row.mid < file_chr_new[1,"mid"]){
      file_chr_new <- rbind(boundary1.row,file_chr_new)
    }
    if (boundary2.row.mid > file_chr_new[nrow(file_chr_new),"mid"]){
      file_chr_new <- rbind(file_chr_new,boundary2.row)
    }
    

    file_chr_new$normalized_rt.next <- c(file_chr_new$normalized_rt[2:nrow(file_chr_new)],NA)
    file_chr_new$mid.next <- c(file_chr_new$mid[2:nrow(file_chr_new)],NA)
    file_chr_new$size <- file_chr_new$mid.next-file_chr_new$mid
    file_chr_new$slope <- (file_chr_new$normalized_rt.next-file_chr_new$normalized_rt)/file_chr_new$size
    # merge all smoothed chromosome
    file_chr_new_all <- rbind(file_chr_new_all,file_chr_new)
  }
  
  names(file_chr_new_all) <- c("chr","V2","V3","normalized_rt","start","normalized_rt.next","end","size","slope")
  #write.table(peak_valley_all,file.path(direction_path,paste0(outname,".txt")),sep = "\t",row.names = F,quote = F)
  write.table(file_chr_new_all,file.path(direction_path,paste0(paste0(gsub("direction","bin",outname)),".txt")),sep = "\t",row.names = F,quote = F)
  return(file_chr_new_all)
}

bin2window_function <- function(rtpath = file.path(direction_path,"HCT116_wb_bin_quantile_norm_14cellline_range50000"),slope_cutoff=5e-7){
  file <- read.table(paste0(rtpath,".txt"),header = T)
  mychr=paste0("chr",c(seq(1,22,1),"X","Y"))
  file <- file[file$chr %in% mychr, ]
  
  file_chr_new_all <- data.frame()
  for (chri in unique(file$chr)){
    file_chr <- file[file$chr==chri,]
    file_chr$direction <- ifelse(abs(file_chr$slope)<slope_cutoff,"flat",
                                 ifelse(file_chr$slope<0,"right","left"))
    
    right.row <- grep("right",file_chr$direction)
    right.row.pair <- unname(tapply(right.row, cumsum(c(1, diff(right.row)) != 1), range))
    right.row.pair1 <- sapply(right.row.pair,"[[",1)
    right.row.pair2 <- sapply(right.row.pair,"[[",2)
    
    left.row <- grep("left",file_chr$direction)
    left.row.pair <- unname(tapply(left.row, cumsum(c(1, diff(left.row)) != 1), range))
    left.row.pair1 <- sapply(left.row.pair,"[[",1)
    left.row.pair2 <- sapply(left.row.pair,"[[",2)
    
    flat.row <- grep("flat",file_chr$direction)
    flat.row.pair <- unname(tapply(flat.row, cumsum(c(1, diff(flat.row)) != 1), range))
    flat.row.pair1 <- sapply(flat.row.pair,"[[",1)
    flat.row.pair2 <- sapply(flat.row.pair,"[[",2)
    
    file_chr_new <- rbind(data.frame(
      "chr"=chri,
      "start"=file_chr[right.row.pair1,"start"],
      "end"=file_chr[right.row.pair2,"end"],
      "rt1"=file_chr[right.row.pair1,"normalized_rt"],
      "rt2"=file_chr[right.row.pair2,"normalized_rt.next"],
      "direction"="right"
    ),
    data.frame(
      "chr"=chri,
      "start"=file_chr[left.row.pair1,"start"],
      "end"=file_chr[left.row.pair2,"end"],
      "rt1"=file_chr[left.row.pair1,"normalized_rt"],
      "rt2"=file_chr[left.row.pair2,"normalized_rt.next"],
      "direction"="left"
    ),
    data.frame(
      "chr"=chri,
      "start"=file_chr[flat.row.pair1,"start"],
      "end"=file_chr[flat.row.pair2,"end"],
      "rt1"=file_chr[flat.row.pair1,"normalized_rt"],
      "rt2"=file_chr[flat.row.pair2,"normalized_rt.next"],
      "direction"="flat"
    )
    )
    
    file_chr_new <- file_chr_new[order(file_chr_new$start),]
    file_chr_new_all <- rbind(file_chr_new_all,file_chr_new)
  }
  write.table(file_chr_new_all,file.path(paste0(paste0(gsub("bin","window",rtpath)),"_slope",slope_cutoff,".txt")),sep = "\t",row.names = F,quote = F)
}



plot_direction_bin_funciton  <- function(name=HCT116_wb_bin,
                                     color=HCT116_wb,
                                     mychr="chr1",
                                     range=c(50000000,75000000),
                                     titlename="test",shrink=25000,shrinkplot=TRUE){
  # read orginal file  add  normalized  rt value
  file = name
  
  # file_plot is the region need to be plotted
  file_plot <- file[file$chr==mychr & file$start >= range[1] & file$start <= range[2],]
  # color is the color file
  color$slope_simple <- ifelse(color$direction=="left" ,1,NA)
  color$slope_simple <- ifelse(color$direction=="right" ,-1,color$slope_simple)
  # add slope_simple to file_plot
  file_plot$slope_simple <- unlist(apply(file_plot,1,function(x) 
    color[grep(TRUE,color$chr==x[1] & color$start<=as.numeric(x["start"]) & color$end>as.numeric(x["start"])),"slope_simple"][1]
  ))
  file_plot$slope_simple <- ifelse(is.na(file_plot$slope_simple),0,file_plot$slope_simple)
  # add boundary
  color$boundary1 <- ifelse(color$direction %in% c("left","right"),color$start+shrink,NA)
  color$boundary2 <- ifelse(color$direction %in% c("left","right"),color$end-shrink,NA)
  color.boundary <- color[color$chr==mychr,]
  
  out <-
    ggplot()  +
    geom_line(data = file_plot,
              aes(x = start,
                  y = normalized_rt,
                  color = slope_simple),size = 0.018)  +
    labs(x = mychr, y = "Normalized\nReplication Timing",title = titlename) +
    scale_x_continuous(breaks = c(range[1],range[2]))+
    scale_color_gradient2(
      low = "blue",
      midpoint = 0,
      mid = "black",
      high = "green"
    ) +
    theme_bw()+
    theme(
      plot.title = element_text(hjust = 0.5,size=6),
      legend.position = "none",
      axis.line = element_line(size=0.018),
      axis.text = element_text(size=6,colour = "black"),
      axis.title = element_text(size=6,colour = "black"),
      panel.background = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  if (shrinkplot==TRUE){
    out <- out + geom_rect(aes(xmin=color.boundary$boundary1,
                               xmax=color.boundary$boundary2,
                               ymin=-Inf,
                               ymax=Inf),alpha=0.1)
  } else {
    out <- out
  }
  return(out)
}



in_domain_funciton <- function(input,file){
  chri  <- input[1]
  starti <- as.numeric(input[4])
  endi  <- as.numeric(input[5])
  tell <- file$chr  == chri & file$start<= starti & file$end >= endi
  if (sum(tell)>0){
    return(file[grep(TRUE,tell),"direction"])
  } else  {
    return("Splited")
  }
}
plot_gene_function <- function(mychr=mychr,range=myrange,window=cl.50000.window,shrink=25000){
  coding.gene.plot <- hg19_coding[hg19_coding$V1==mychr & hg19_coding$V4>=range[1] & hg19_coding$V5<= range[2],]
  
  if (nrow(coding.gene.plot)==0){
    return(NA)
  }else {
    coding.gene.plot$y <- rep(c(1,9,5,3,7,2,8,4,6),1000)[1:nrow(coding.gene.plot)]
    coding.gene.plot$direction <- apply(coding.gene.plot, 1, in_domain_funciton,window)
    # add boundary
    window$boundary1 <- ifelse(window$direction %in% c("left","right"),window$start+shrink,NA)
    window$boundary2 <- ifelse(window$direction %in% c("left","right"),window$end-shrink,NA)
    window <- window[window$chr %in% mychr,]
    
    if (nrow(coding.gene.plot)>20){
      ggplot()+
        geom_linerange(data=coding.gene.plot,aes(xmin=V4,xmax=V5,y=y,color=direction),size=2)+
        geom_rect(aes(xmin=window$boundary1, 
                      xmax=window$boundary2, 
                      ymin=-Inf,
                      ymax=Inf),alpha=0.1)+
        xlim(range[1],range[2])+
        ylab("Gene")+
        scale_color_manual(values = c("green","blue","pink","black"),
                           breaks = c("left","right","flat","Splited"))+
        theme_bw()+
        theme(axis.title.x = element_blank(),
              axis.title.y = element_text(size=6,colour = "black"),
              axis.text.y = element_blank(),
              axis.text.x = element_text(size=6,colour = "black"),
              axis.ticks.length.y = unit(0,"mm"),
              panel.background = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank(),
              legend.position = "none"
              # panel.grid.major.x = element_blank(),
              # panel.grid.minor.x = element_blank()
        )
    }else{
      ggplot()+
        geom_linerange(data=coding.gene.plot,aes(xmin=V4,xmax=V5,y=y,color=direction),size=2)+
        geom_rect(aes(xmin=window$boundary1, 
                      xmax=window$boundary2, 
                      ymin=-Inf,
                      ymax=Inf),alpha=0.1)+
        geom_text(data=coding.gene.plot,aes(x=V4,y=y,label=gene_name),hjust=1,size=2)+
        xlim(range[1],range[2])+
        ylab("Gene")+
        scale_color_manual(values = c("green","blue","pink","black"),
                           breaks = c("left","right","flat","Splited"))+
        theme_bw()+
        theme(axis.title.x = element_blank(),
              axis.title.y = element_text(size=6,colour = "black"),
              axis.text.y = element_blank(),
              axis.text.x = element_text(size=6,colour = "black"),
              axis.ticks.length.y = unit(0,"mm"),
              panel.background = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank(),
              legend.position = "none"
              # panel.grid.major.x = element_blank(),
              # panel.grid.minor.x = element_blank()
        )
    }
  }

}


##################### run ##################
# add all rt raw value together
rm(normalize.df)
normalize.df <- add_raw_rt_function(name="HCT116_final")
normalize.df <- add_raw_rt_function(name="HEK293_final")
normalize.df <- add_raw_rt_function(name="RPE1_hTERT_final")
normalize.df <- add_raw_rt_function(name="U2OS_final")
normalize.df <- add_raw_rt_function(name="A13_final")
normalize.df <- add_raw_rt_function(name="BG01_NPC_final")
normalize.df <- add_raw_rt_function(name="CyT49_Neural_Crest_final")
normalize.df <- add_raw_rt_function(name="CyT49_Liver_final")
normalize.df <- add_raw_rt_function(name="CyT49_Panc_final")
normalize.df <- add_raw_rt_function(name="CyT49_Definitive_Endoderm_final")
normalize.df <- add_raw_rt_function(name="HeLaS3_1")
normalize.df <- add_raw_rt_function(name="BG02_ESC_1") #BG02_ESC_1 Bg02es_final
normalize.df <- add_raw_rt_function(name="MCF7_1") #MCF7_1 MCF7_final
normalize.df <- add_raw_rt_function(name="HepG2_1") #HepG2_1 HepG2_final

# quantile normalization
normalize.mat <- as.matrix(normalize.df)
normalize.mat_normalized <- data.frame(normalize.quantiles(as.matrix(normalize.mat)))
names(normalize.mat_normalized) <- names(normalize.df)

### generate bin file at sepecific bin size
HCT116_wb <- rt_direction_function(name="HCT116_final",outname = "HCT116_wb_direction_quantile_norm_14cellline_range50000",range_size = 50000)
HEK293_wb <- rt_direction_function(name="HEK293_final",outname = "HEK293_wb_direction_quantile_norm_14cellline_range50000",range_size = 50000)
HeLaS3_wb <- rt_direction_function(name="HeLaS3_1",outname = "HeLaS3_wb_direction_quantile_norm_14cellline_range50000",range_size = 50000)
RPE1_hTERT_wb <- rt_direction_function(name="RPE1_hTERT_final",outname = "RPE1_hTERT_wb_direction_quantile_norm_14cellline_range50000",range_size = 50000)
U2OS_wb <- rt_direction_function(name="U2OS_final",outname = "U2OS_wb_direction_quantile_norm_14cellline_range50000",range_size = 50000)
A13_wb <- rt_direction_function(name="A13_final",outname = "A13_wb_direction_quantile_norm_14cellline_range50000",range_size = 50000)
BG01_NPC_wb <- rt_direction_function(name="BG01_NPC_final",outname = "BG01_NPC_wb_direction_quantile_norm_14cellline_range50000",range_size = 50000)
CyT49_Neural_Crest_wb <- rt_direction_function(name="CyT49_Neural_Crest_final",outname = "CyT49_Neural_Crest_wb_direction_quantile_norm_14cellline_range50000",range_size = 50000)
CyT49_Liver_wb <- rt_direction_function(name="CyT49_Liver_final",outname = "CyT49_Liver_wb_direction_quantile_norm_14cellline_range50000",range_size = 50000)
CyT49_Panc_wb <- rt_direction_function(name="CyT49_Panc_final",outname = "CyT49_Panc_wb_direction_quantile_norm_14cellline_range50000",range_size = 50000)
CyT49_Definitive_Endoderm_wb <- rt_direction_function(name="CyT49_Definitive_Endoderm_final",outname = "CyT49_Definitive_Endoderm_wb_direction_quantile_norm_14cellline_range50000",range_size = 50000)
BG02_ESC_wb <- rt_direction_function(name="BG02_ESC_1",outname = "BG02_ESC_wb_direction_quantile_norm_14cellline_range50000",range_size = 50000) #BG02_ESC_1 Bg02es_final
MCF7_wb <- rt_direction_function(name="MCF7_1",outname = "MCF7_wb_direction_quantile_norm_14cellline_range50000",range_size = 50000) #MCF7_1 MCF7_final
HepG2_wb <- rt_direction_function(name="HepG2_1",outname = "HepG2_wb_direction_quantile_norm_14cellline_range50000",range_size = 50000) #HepG2_1 HepG2_final

### generate window file with sepecific slope cutoff
slope_cutoff=5e-7
bin2window_function(rtpath = file.path(direction_path,"HCT116_wb_bin_quantile_norm_14cellline_range50000"),slope_cutoff=slope_cutoff)
bin2window_function(rtpath = file.path(direction_path,"HEK293_wb_bin_quantile_norm_14cellline_range50000"),slope_cutoff=slope_cutoff)
bin2window_function(rtpath = file.path(direction_path,"HeLaS3_wb_bin_quantile_norm_14cellline_range50000"),slope_cutoff=slope_cutoff)
bin2window_function(rtpath = file.path(direction_path,"RPE1_hTERT_wb_bin_quantile_norm_14cellline_range50000"),slope_cutoff=slope_cutoff)
bin2window_function(rtpath = file.path(direction_path,"U2OS_wb_bin_quantile_norm_14cellline_range50000"),slope_cutoff=slope_cutoff)
bin2window_function(rtpath = file.path(direction_path,"A13_wb_bin_quantile_norm_14cellline_range50000"),slope_cutoff=slope_cutoff)
bin2window_function(rtpath = file.path(direction_path,"BG01_NPC_wb_bin_quantile_norm_14cellline_range50000"),slope_cutoff=slope_cutoff)
bin2window_function(rtpath = file.path(direction_path,"CyT49_Neural_Crest_wb_bin_quantile_norm_14cellline_range50000"),slope_cutoff=slope_cutoff)
bin2window_function(rtpath = file.path(direction_path,"CyT49_Liver_wb_bin_quantile_norm_14cellline_range50000"),slope_cutoff=slope_cutoff)
bin2window_function(rtpath = file.path(direction_path,"CyT49_Panc_wb_bin_quantile_norm_14cellline_range50000"),slope_cutoff=slope_cutoff)
bin2window_function(rtpath = file.path(direction_path,"CyT49_Definitive_Endoderm_wb_bin_quantile_norm_14cellline_range50000"),slope_cutoff=slope_cutoff)
bin2window_function(rtpath = file.path(direction_path,"BG02_ESC_wb_bin_quantile_norm_14cellline_range50000"),slope_cutoff=slope_cutoff)
bin2window_function(rtpath = file.path(direction_path,"MCF7_wb_bin_quantile_norm_14cellline_range50000"),slope_cutoff=slope_cutoff)
bin2window_function(rtpath = file.path(direction_path,"HepG2_wb_bin_quantile_norm_14cellline_range50000"),slope_cutoff=slope_cutoff)

### plot replication timing
mycellline <- c("HCT116","HEK293","RPE1_hTERT","U2OS","A13",
                             "BG01_NPC","CyT49_Neural_Crest","CyT49_Liver","CyT49_Panc","CyT49_Definitive_Endoderm",
                             "HeLaS3","BG02_ESC","MCF7","HepG2")
for (cellline in mycellline){
  myrange <- c(50000000,75000000)
  mychr <- "chr1"
  slope_cutoff=5e-7
  options(scipen = 3)
  ### publication
  binsize=50000
  cl.50000.window <- read.delim(file.path(direction_path,paste0(cellline,"_wb_window_quantile_norm_14cellline_range",binsize,"_slope",slope_cutoff,".txt")))
  cl.50000.bin <- read.delim(file.path(direction_path,paste0(cellline,"_wb_bin_quantile_norm_14cellline_range",binsize,".txt")))
  cl.50000 <- plot_direction_bin_funciton(name=cl.50000.bin,color=cl.50000.window,mychr=mychr,range= myrange,titlename=paste0(cellline,"_bin",binsize,"_slope",slope_cutoff),shrink=25000,shrinkplot=FALSE)
  options(scipen = 10)
  ggsave(file.path(plot_path,paste0(cellline,".",mychr,"-",myrange[1],"-",myrange[2],".pub.pdf")),height = 50,width = 55,units = "mm",dpi = 300)
  
}
