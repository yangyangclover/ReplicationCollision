# library
library(dplyr)
library(stringr)

#########################################################################################################
############# This part is for bins #########
#########################################################################################################
### function
in_domain_funciton <- function(input,file,halfboundary){
  chri  <- input[1]
  starti <- as.numeric(input[2])
  endi  <- as.numeric(input[3])
  tell <- file$chr  == chri & ((file$start+halfboundary)<= starti) & ((file$end-halfboundary) >= endi)
  if (sum(tell)>0){
    return(file[grep(TRUE,tell),"direction"])
  } else  {
    return("Splited")
  }
}

bin_slope_function <- function(input=input,
                               rtpath = file.path(direction_path,"HCT116_wb_bin_quantile_norm_14cellline_range50000.txt"),
                               celllinename="HCT116",
                               halfboundary=25000){
  file <- read.table(rtpath,header = T)
  mychr=paste0("chr",c(seq(1,22,1),"X","Y"))
  file <- file[file$chr %in% mychr, ]
  
  file_chr_new_all <- data.frame()
  # add direction
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
  input[,celllinename] <- apply(input, 1, in_domain_funciton,file_chr_new_all,halfboundary)
  return(input)
}



### path
hg19_coding_path <- file.path("D:/Data/REF/hg19/coding_gene/coding_gene.gtf")
rt_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/genomic_features/rt_cellline"
direction_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/genomic_features/rt_cellline/direction"
outpath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/genomic_features/rt_cellline/direction/codinggene"

###parameters
mycellline  <- c("HCT116","HEK293","RPE1_hTERT","U2OS","A13",
                 "BG01_NPC","CyT49_Neural_Crest","CyT49_Liver","CyT49_Panc","CyT49_Definitive_Endoderm",
                 "HeLaS3","BG02_ESC","MCF7","HepG2") 
mychr <- paste0("chr",c(seq(1,22,1),"X","Y"))
slope_cutoff  <- 5e-7

### read files and annotate
hg19_coding <- read.delim(hg19_coding_path,sep = "\t",stringsAsFactors = F,header = F) 
hg19_coding <- hg19_coding %>% filter(V1 %in% mychr)
hg19_coding$gene_name <- str_extract(hg19_coding$V9,"gene_name.*; transcript_type") %>% gsub("gene_name ","",.) %>% gsub("; transcript_type","",.)
hg19_coding$ensembl_gene_id <- str_extract(hg19_coding$V9,"gene_id.*; transcript_id") %>% gsub("gene_id ","",.) %>% gsub("; transcript_id","",.)
hg19_coding$ensembl_gene_id <- substr(hg19_coding$ensembl_gene_id,1,15)

### direction
halfboundary <- 25000
input <- hg19_coding[,c("V1","V4","V5","V7","gene_name","ensembl_gene_id")]
input <- bin_slope_function(input=input,rtpath = file.path(direction_path,"HCT116_wb_bin_quantile_norm_14cellline_range50000.txt"),celllinename="HCT116",halfboundary=halfboundary)
input <- bin_slope_function(input=input,rtpath=file.path(direction_path,"HEK293_wb_bin_quantile_norm_14cellline_range50000.txt"),celllinename="HEK293",halfboundary=halfboundary)
input <- bin_slope_function(input=input,rtpath=file.path(direction_path,"HeLaS3_wb_bin_quantile_norm_14cellline_range50000.txt"),celllinename="HeLaS3",halfboundary=halfboundary)
input <- bin_slope_function(input=input,rtpath=file.path(direction_path,"RPE1_hTERT_wb_bin_quantile_norm_14cellline_range50000.txt"),celllinename="RPE1_hTERT",halfboundary=halfboundary)
input <- bin_slope_function(input=input,rtpath=file.path(direction_path,"U2OS_wb_bin_quantile_norm_14cellline_range50000.txt"),celllinename="U2OS",halfboundary=halfboundary)
input <- bin_slope_function(input=input,rtpath=file.path(direction_path,"A13_wb_bin_quantile_norm_14cellline_range50000.txt"),celllinename="A13",halfboundary=halfboundary)
input <- bin_slope_function(input=input,rtpath=file.path(direction_path,"BG01_NPC_wb_bin_quantile_norm_14cellline_range50000.txt"),celllinename="BG01_NPC",halfboundary=halfboundary)
input <- bin_slope_function(input=input,rtpath=file.path(direction_path,"CyT49_Neural_Crest_wb_bin_quantile_norm_14cellline_range50000.txt"),celllinename="CyT49_Neural_Crest",halfboundary=halfboundary)
input <- bin_slope_function(input=input,rtpath=file.path(direction_path,"CyT49_Liver_wb_bin_quantile_norm_14cellline_range50000.txt"),celllinename="CyT49_Liver",halfboundary=halfboundary)
input <- bin_slope_function(input=input,rtpath=file.path(direction_path,"CyT49_Panc_wb_bin_quantile_norm_14cellline_range50000.txt"),celllinename="CyT49_Panc",halfboundary=halfboundary)
input <- bin_slope_function(input=input,rtpath=file.path(direction_path,"CyT49_Definitive_Endoderm_wb_bin_quantile_norm_14cellline_range50000.txt"),celllinename="CyT49_Definitive_Endoderm",halfboundary=halfboundary)
input <- bin_slope_function(input=input,rtpath=file.path(direction_path,"BG02_ESC_wb_bin_quantile_norm_14cellline_range50000.txt"),celllinename="BG02_ESC",halfboundary=halfboundary)
input <- bin_slope_function(input=input,rtpath=file.path(direction_path,"MCF7_wb_bin_quantile_norm_14cellline_range50000.txt"),celllinename="MCF7",halfboundary=halfboundary)
input <- bin_slope_function(input=input,rtpath=file.path(direction_path,"HepG2_wb_bin_quantile_norm_14cellline_range50000.txt"),celllinename="HepG2",halfboundary=halfboundary)

### add consistent and write
cellline_num <- c(8,10)[1]
input$consistent <-
  apply(input[, mycellline], 1, function(x)
    if (sum(x %in% "left") >= cellline_num) {
      return("left")
    } else if (sum(x %in% "right") >= cellline_num) {
      return("right")
    }  else {
      return(NA)
    })
sum(input$consistent %in% c("left","right"))
filename <- paste0("coding_geng_all_rtcellline_direction_","slope",slope_cutoff,"_cellline",cellline_num,"_shrink",halfboundary,"_flatonbin",".tsv")
write.table(input,file.path(outpath,filename),sep = "\t",row.names = F,quote = F)
