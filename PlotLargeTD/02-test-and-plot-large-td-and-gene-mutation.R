# library
library(dplyr)
library(ggplot2)
library(reshape)
library(tidyr)
library(stringr)
require(openxlsx)
library(cowplot)


# path
coding_gene <- read.delim("D:/Data/REF/hg19/coding_gene/coding_gene_position_hg19.bed",header = F)
names(coding_gene) <- c("chr","start","end","ensg","gene")
coding_gene <- coding_gene[coding_gene$gene  !=  "",] #remove non coding
coding_gene <- coding_gene[duplicated(coding_gene$gene)==FALSE,] #remove duplicated
ensg_gene_match  <-  coding_gene$gene
names(ensg_gene_match) <- coding_gene$ensg
plotpath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.snv/plot"

# function
plot_function <- function(sv,allsv,snv_t,sample,mygene=c("TP53","CDK12"),savepath,ptext,name,plotwide){
  all_sample <- data.frame(unique(sample[,c("WGS","histology_abbreviation")]))
  names(all_sample)[1] <- "sample"
  ### large td
  sample_largetd <- sv %>% filter(sig_group %in% paste0("td_",seq(10,16,1)))%>% group_by(sample) %>% dplyr::count()
  names(sample_largetd)[2] <- "largetd"
  all_sample <- merge(all_sample,sample_largetd,
                      by.x = "sample",by.y = "sample",all.x = T)
  all_sample$largetd <- ifelse(is.na(all_sample$largetd),0,all_sample$largetd)
  
  ### simple sv
  sample_simplesv <- sv %>% group_by(sample) %>% dplyr::count()
  names(sample_simplesv)[2] <- "simplesv"
  all_sample <- merge(all_sample,sample_simplesv,
                      by.x = "sample",by.y = "sample",all.x = T)
  all_sample$simplesv <- ifelse(is.na(all_sample$simplesv),0,all_sample$simplesv)
  all_sample$simplesv  <- log10(all_sample$simplesv)
  
  ### all sv
  sample_allsv <- allsv %>% group_by(sample) %>% dplyr::count()
  names(sample_allsv)[2] <- "allsv"
  all_sample <- merge(all_sample,sample_allsv,
                      by.x = "sample",by.y = "sample",all.x = T)
  all_sample$allsv <- ifelse(is.na(all_sample$allsv),0,all_sample$allsv)
  all_sample$allsv  <- log10(all_sample$allsv)
  
  ### merge
  all_sample  <- merge(all_sample,snv_t,by.x = "sample",by.y = "sample",all.x = T)
  
  ### -to.
  all_sample$histology_abbreviation <- gsub("-","",all_sample$histology_abbreviation)
  all_sample$histology_abbreviation <- gsub(" ","",all_sample$histology_abbreviation)
  all_sample$histology_abbreviation <- gsub("/","",all_sample$histology_abbreviation)
  
  ### plot
  sample_order <- all_sample[order(all_sample$largetd,decreasing = T),"sample"]
  all_sample_plot <- all_sample[order(all_sample$largetd,decreasing = T),c("sample","largetd",mygene)]
  all_sample_plot$sample <- factor(all_sample_plot$sample,levels = sample_order)
  all_sample_plot[,mygene] <- apply(as.data.frame(all_sample_plot[,mygene]), 2, as.factor)
  all_sample_plot_long <- gather(all_sample_plot,gene,type,mygene,factor_key = T)
  all_sample_plot_long$sample <- factor(all_sample_plot_long$sample,levels = sample_order)
  
  p1 <- ggplot(data=all_sample_plot)+
    geom_bar(aes(x=sample,y=largetd),stat = "identity",color="black",fill="black")+
    annotate(geom="text", x=nrow(all_sample_plot)/2, y=max(all_sample_plot$largetd)/2, label=ptext,color="red",size=6*0.35)+
    ylab("Number of large TDs")+
    scale_y_continuous(expand = c(0,0))+
    theme_bw() +
    theme(axis.title.y = element_text(size=6,color="black"),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size=6,color="black"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line = element_line(colour = "black",size=0.1189),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.margin=grid::unit(c(2,2,0,0), "mm"))
  p2 <- ggplot(data=all_sample_plot_long)+
    geom_bar(aes(x=sample,y=1,fill=type),
             stat = "identity", 
             width = 1)+
    facet_grid(gene~.,switch = "y",scales='free')+
    scale_fill_manual(breaks = c(0,1,2),
                      values = c("white","white","red"))+
    coord_cartesian(clip = "off")+
    ylab(mygene)+
    xlab(name)+
    theme_bw()+
    theme(axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.title.x = element_text(angle = 0,size=6,color="black"),
          axis.ticks = element_blank(),
          # axis.line.x.bottom  = element_line(colour = "black",size=0.1189),
          # axis.line.x.top = element_line(colour = "black",size=0.1189),
          # axis.line.y.right = element_line(colour = "black",size=0.1189),
          # axis.line.y.left = element_line(colour = "black",size=0.1189),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # panel.border = element_rect(color = "blue",fill = NA,size=0.1189),
          legend.position = "none",
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0,size=6,color="black",hjust=0.5,vjust = 0.5),
          strip.background =element_blank(),
          panel.spacing = unit(1, "mm"),
          plot.margin=grid::unit(c(0,2,0,0), "mm"))
  plot_grid(p1, p2, align = "v", ncol = 1, axis = "l", rel_heights = c(2.5, 0.5+0.75*length(mygene)))
  ggsave(savepath,width = plotwide,height = sum(2.5, 0.5+0.75*length(mygene))*4.5,units = "mm")
}

test_function <- function(sv,allsv,snv_t,mutload,sample,out_path,mutcutoff,prefix){
  if  (!dir.exists(out_path)){
    dir.create(out_path,showWarnings = F,recursive = T)
  }
  
  
  all_sample <- data.frame(unique(sample[,c("WGS","histology_abbreviation")]))
  names(all_sample)[1] <- "sample"
  ### large td
  sample_largetd <- sv %>% filter(sig_group %in% paste0("td_",seq(10,16,1)))%>% group_by(sample) %>% dplyr::count()
  names(sample_largetd)[2] <- "largetd"
  all_sample <- merge(all_sample,sample_largetd,
                      by.x = "sample",by.y = "sample",all.x = T)
  all_sample$largetd <- ifelse(is.na(all_sample$largetd),0,all_sample$largetd)
  
  ### simple sv
  sample_simplesv <- sv %>% group_by(sample) %>% dplyr::count()
  names(sample_simplesv)[2] <- "simplesv"
  all_sample <- merge(all_sample,sample_simplesv,
                      by.x = "sample",by.y = "sample",all.x = T)
  all_sample$simplesv <- ifelse(is.na(all_sample$simplesv),0,all_sample$simplesv)
  all_sample$simplesv  <- log10(all_sample$simplesv)
  
  ### all sv
  sample_allsv <- allsv %>% group_by(sample) %>% dplyr::count()
  names(sample_allsv)[2] <- "allsv"
  all_sample <- merge(all_sample,sample_allsv,
                      by.x = "sample",by.y = "sample",all.x = T)
  all_sample$allsv <- ifelse(is.na(all_sample$allsv),0,all_sample$allsv)
  all_sample$allsv  <- log10(all_sample$allsv)
  
  ### mutation load
  mutload$mut_load  <- log10(mutload$mut_load)
  all_sample  <- merge(all_sample,mutload,by.x = "sample",by.y = "sample",all.x = T)
  
  ### merge
  all_sample  <- merge(all_sample,snv_t,by.x = "sample",by.y = "sample",all.x = T)
  
  ### -to.
  all_sample$histology_abbreviation <- gsub("-","",all_sample$histology_abbreviation)
  all_sample$histology_abbreviation <- gsub(" ","",all_sample$histology_abbreviation)
  all_sample$histology_abbreviation <- gsub("/","",all_sample$histology_abbreviation)
  
  ### test
  num.gene.test <- ncol(all_sample) - 7
  wilcox_df <- data.frame()
  fisher_df <- data.frame()
  glm_df <- data.frame()
  lm_df  <- data.frame()
  lm2_df <- data.frame()
  
  ### pick gene mut at least cutoff number
  allgenes  <- names(all_sample)[8:ncol(all_sample)]
  driver.gene.candidates <-  allgenes[(apply(all_sample[,allgenes], 2, function(x) sum(x==2,na.rm = T)) >= mutcutoff) == TRUE]
  # passenger.gene.candidates <-  allgenes[(apply(all_sample[,allgenes], 2, function(x) sum(x>=1,na.rm = T)) >= mutcutoff) == TRUE]
  driver.gene.candidates.num <- length(driver.gene.candidates)
  # passenger.gene.candidates.num <- length(passenger.gene.candidates)
  
  
  
  for  (i  in  1:num.gene.test){
    genei  <- names(all_sample)[i+7]
    driver_tell <- ifelse(all_sample[,genei]==2,1,0)
    driver_tell_factor <- factor(driver_tell,levels = c(1,0,NA))
    driver.num <- sum(driver_tell,na.rm = T)
    nodriver.num <- sum(driver_tell %in% 0,na.rm = T)
    passenger_tell <- ifelse(all_sample[,genei] >=1,1,0)
    passenger_tell_factor <- factor(passenger_tell,levels = c(1,0,NA))
    passenger.num <- sum(passenger_tell,na.rm = T)
    nopassenger.num <- sum(passenger_tell %in% 0,na.rm = T)
    
    

    #### linear regression  model
    linear_df  <- cbind(data.frame("largetd"=log10(all_sample[,"largetd"]+1),
                                     "driver"=ifelse(all_sample[,genei]==2,1,0),
                                     "passenger"=ifelse(all_sample[,genei]>=1,1,0)
    ),
    all_sample[,c("simplesv","allsv","mut_load")])


    if (genei %in% driver.gene.candidates & driver.num>0 & nodriver.num>0){
      driver_mutload <- lm(largetd ~ driver + mut_load, data = linear_df)
      driver_mutload.p <- signif(summary(driver_mutload)[["coefficients"]][2,4],2)
      # driver_mutload.p.bfrn <-  signif(p.adjust(summary(driver_mutload)[["coefficients"]][2,4],method = "bonferroni",n=driver.gene.candidates.num),2)
      driver_mutload.estimate <- signif(summary(driver_mutload)[["coefficients"]][2,1],2)
    }  else {
      driver_mutload <- NA
      driver_mutload.p <- NA
      # driver_mutload.p.bfrn <- NA
      driver_mutload.estimate <- NA
    }


    lm_df <- rbind(lm_df,data.frame(
      "gene"=genei,
      "test"="linear",
      "driver.count"=driver.num,
      "driver.raw.p.with.mutload"=driver_mutload.p,
      "driver_mutload.estimate"=driver_mutload.estimate
    ))

    
    #### linear regression model with tumor type
    if (length(unique(all_sample$histology_abbreviation))>1){
      linear_df  <- cbind(data.frame("sample"=all_sample[,"sample"],
                                     "largetd"=log10(all_sample[,"largetd"]+1),
                                     "histology_abbreviation"=all_sample[,"histology_abbreviation"],
                                     "driver"=ifelse(all_sample[,genei]==2,1,0),
                                     "passenger"=ifelse(all_sample[,genei]>=1,1,0)
      ),
      all_sample[,c("simplesv","allsv","mut_load")])
      allhis <- unique(all_sample$histology_abbreviation)
      linear_df$hist_mark <- 1
      linear_df_wide  <- spread(linear_df,key = histology_abbreviation,value = hist_mark)
      linear_df_wide[is.na(linear_df_wide)] <- 0
      
      if (genei %in% driver.gene.candidates & driver.num>0 & nodriver.num>0){
        myformula <- as.formula(paste0("largetd ~ driver + mut_load +",paste0(allhis,collapse = "+")))
        driver_mutload <- lm(myformula, data = linear_df_wide)
        driver_mutload.p <- signif(summary(driver_mutload)[["coefficients"]][2,4],2)
        # driver_mutload.p.bfrn <-  signif(p.adjust(summary(driver_mutload)[["coefficients"]][2,4],method = "bonferroni",n=driver.gene.candidates.num),2)
        driver_mutload.estimate <- signif(summary(driver_mutload)[["coefficients"]][2,1],2)
      }  else {
        driver_mutload <- NA
        driver_mutload.p <- NA
        # driver_mutload.p.bfrn <- NA
        driver_mutload.estimate <- NA
      }
      

      lm2_df <- rbind(lm2_df,data.frame(
        "gene"=genei,
        "test"="linear",
        "driver.count"=driver.num,
        "driver.raw.p.with.mutload"=driver_mutload.p,
        "driver_mutload.estimate"=driver_mutload.estimate
      ))
      
      
    }
  }

  
  all.sample.number <- nrow(sample)

  
  ten.of.100 <- all.sample.number*10/100
  five.of.100 <- all.sample.number*5/100
  four.of.100 <- all.sample.number*4/100
  three.of.100 <- all.sample.number*3/100
  two.of.100 <- all.sample.number*2/100
  
  ten.of.100.sample <- lm_df$driver.count >= ten.of.100
  five.of.100.sample <- lm_df$driver.count >= five.of.100
  four.of.100.sample <- lm_df$driver.count >= four.of.100
  three.of.100.sample <- lm_df$driver.count >= three.of.100
  two.of.100.sample <- lm_df$driver.count >= two.of.100
  lm_df[,"driver.bfrn.p.with.mutload"] <- p.adjust(lm_df$driver.raw.p.with.mutload,method = "bonferroni",n=driver.gene.candidates.num)
  lm_df[ten.of.100.sample,"driver.bfrn.p.with.mutload.10of100"] <- p.adjust(lm_df$driver.raw.p.with.mutload[ten.of.100.sample],
                                                                            method = "bonferroni",
                                                                            n=sum(ten.of.100.sample))
  lm_df[five.of.100.sample,"driver.bfrn.p.with.mutload.5of100"] <- p.adjust(lm_df$driver.raw.p.with.mutload[five.of.100.sample],
                                                                   method = "bonferroni",
                                                                   n=sum(five.of.100.sample))
  lm_df[four.of.100.sample,"driver.bfrn.p.with.mutload.4of100"] <- p.adjust(lm_df$driver.raw.p.with.mutload[four.of.100.sample],
                                                                            method = "bonferroni",
                                                                            n=sum(four.of.100.sample))
  lm_df[three.of.100.sample,"driver.bfrn.p.with.mutload.3of100"] <- p.adjust(lm_df$driver.raw.p.with.mutload[three.of.100.sample],
                                                                            method = "bonferroni",
                                                                            n=sum(three.of.100.sample))
  lm_df[two.of.100.sample,"driver.bfrn.p.with.mutload.2of100"] <- p.adjust(lm_df$driver.raw.p.with.mutload[two.of.100.sample],
                                                                            method = "bonferroni",
                                                                            n=sum(two.of.100.sample))
  lm_df$driver.fdr.p.with.mutload <- p.adjust(lm_df$driver.raw.p.with.mutload,method = "fdr",n=driver.gene.candidates.num)
  lm_df[ten.of.100.sample,"driver.fdr.p.with.mutload.10of100"] <- p.adjust(lm_df$driver.raw.p.with.mutload[ten.of.100.sample],
                                                                           method = "fdr",
                                                                           n=sum(ten.of.100.sample))
  lm_df[five.of.100.sample,"driver.fdr.p.with.mutload.5of100"] <- p.adjust(lm_df$driver.raw.p.with.mutload[five.of.100.sample],
                                                                            method = "fdr",
                                                                            n=sum(five.of.100.sample))
  lm_df[four.of.100.sample,"driver.fdr.p.with.mutload.4of100"] <- p.adjust(lm_df$driver.raw.p.with.mutload[four.of.100.sample],
                                                                            method = "fdr",
                                                                            n=sum(four.of.100.sample))
  lm_df[three.of.100.sample,"driver.fdr.p.with.mutload.3of100"] <- p.adjust(lm_df$driver.raw.p.with.mutload[three.of.100.sample],
                                                                             method = "fdr",
                                                                             n=sum(three.of.100.sample))
  lm_df[two.of.100.sample,"driver.fdr.p.with.mutload.2of100"] <- p.adjust(lm_df$driver.raw.p.with.mutload[two.of.100.sample],
                                                                           method = "fdr",
                                                                           n=sum(two.of.100.sample))
  
  lm_df <- lm_df[order(lm_df$driver.raw.p.with.mutload,decreasing = F),c("gene","test","driver.count","driver.raw.p.with.mutload",
                    "driver.bfrn.p.with.mutload",
                    "driver.bfrn.p.with.mutload.10of100",
                    "driver.bfrn.p.with.mutload.5of100","driver.bfrn.p.with.mutload.4of100",
                    "driver.bfrn.p.with.mutload.3of100","driver.bfrn.p.with.mutload.2of100",
                    "driver.fdr.p.with.mutload",
                    "driver.fdr.p.with.mutload.10of100",
                    "driver.fdr.p.with.mutload.5of100","driver.fdr.p.with.mutload.4of100",
                    "driver.fdr.p.with.mutload.3of100","driver.fdr.p.with.mutload.2of100",
                    "driver_mutload.estimate"
                    )]
  

  write.csv(lm_df,file.path(out_path,paste0(prefix,"_large.td_and_snv_linear_test.csv")),quote = F,row.names = F)
  
  if (nrow(lm2_df)>0){
    lm2_df[,"driver.bfrn.p.with.mutload"] <- p.adjust(lm2_df$driver.raw.p.with.mutload,method = "bonferroni",n=driver.gene.candidates.num)
    lm2_df[ten.of.100.sample,"driver.bfrn.p.with.mutload.10of100"] <- p.adjust(lm2_df$driver.raw.p.with.mutload[ten.of.100.sample],
                                                                               method = "bonferroni",
                                                                               n=sum(ten.of.100.sample))
    lm2_df[five.of.100.sample,"driver.bfrn.p.with.mutload.5of100"] <- p.adjust(lm2_df$driver.raw.p.with.mutload[five.of.100.sample],
                                                                              method = "bonferroni",
                                                                              n=sum(five.of.100.sample))
    lm2_df[four.of.100.sample,"driver.bfrn.p.with.mutload.4of100"] <- p.adjust(lm2_df$driver.raw.p.with.mutload[four.of.100.sample],
                                                                              method = "bonferroni",
                                                                              n=sum(four.of.100.sample))
    lm2_df[three.of.100.sample,"driver.bfrn.p.with.mutload.3of100"] <- p.adjust(lm2_df$driver.raw.p.with.mutload[three.of.100.sample],
                                                                               method = "bonferroni",
                                                                               n=sum(three.of.100.sample))
    lm2_df[two.of.100.sample,"driver.bfrn.p.with.mutload.2of100"] <- p.adjust(lm2_df$driver.raw.p.with.mutload[two.of.100.sample],
                                                                             method = "bonferroni",
                                                                             n=sum(two.of.100.sample))
    lm2_df$driver.fdr.p.with.mutload <- p.adjust(lm2_df$driver.raw.p.with.mutload,method = "fdr",n=driver.gene.candidates.num)
    lm2_df[ten.of.100.sample,"driver.fdr.p.with.mutload.10of100"] <- p.adjust(lm2_df$driver.raw.p.with.mutload[ten.of.100.sample],
                                                                              method = "fdr",
                                                                              n=sum(ten.of.100.sample))
    lm2_df[five.of.100.sample,"driver.fdr.p.with.mutload.5of100"] <- p.adjust(lm2_df$driver.raw.p.with.mutload[five.of.100.sample],
                                                                             method = "fdr",
                                                                             n=sum(five.of.100.sample))
    lm2_df[four.of.100.sample,"driver.fdr.p.with.mutload.4of100"] <- p.adjust(lm2_df$driver.raw.p.with.mutload[four.of.100.sample],
                                                                             method = "fdr",
                                                                             n=sum(four.of.100.sample))
    lm2_df[three.of.100.sample,"driver.fdr.p.with.mutload.3of100"] <- p.adjust(lm2_df$driver.raw.p.with.mutload[three.of.100.sample],
                                                                              method = "fdr",
                                                                              n=sum(three.of.100.sample))
    lm2_df[two.of.100.sample,"driver.fdr.p.with.mutload.2of100"] <- p.adjust(lm2_df$driver.raw.p.with.mutload[two.of.100.sample],
                                                                            method = "fdr",
                                                                            n=sum(two.of.100.sample))
    lm2_df <- lm2_df[order(lm2_df$driver.raw.p.with.mutload,decreasing = F),c("gene","test","driver.count","driver.raw.p.with.mutload",
                      "driver.bfrn.p.with.mutload",
                      "driver.bfrn.p.with.mutload.10of100",
                      "driver.bfrn.p.with.mutload.5of100","driver.bfrn.p.with.mutload.4of100",
                      "driver.bfrn.p.with.mutload.3of100","driver.bfrn.p.with.mutload.2of100",
                      "driver.fdr.p.with.mutload",
                      "driver.fdr.p.with.mutload.10of100",
                      "driver.fdr.p.with.mutload.5of100","driver.fdr.p.with.mutload.4of100",
                      "driver.fdr.p.with.mutload.3of100","driver.fdr.p.with.mutload.2of100",
                      "driver_mutload.estimate"
                      # "passenger.raw.p.with.mutload","passenger.bfrn.p.with.mutload","passenger.fdr.p.with.mutload","passenger_mutload.estimate"
    )]
    write.csv(lm2_df,file.path(out_path,paste0(prefix,"_large.td_and_snv_linear_tumortypevariable_test.csv")),quote = F,row.names = F)
  }
}


##################################### RUN #############################
############# BRCA-EU #####################
sample_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/BRCA-EU/brcaeu532_sampleid_match.csv"
sv_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/BRCA-EU/real_and_random_withsig.txt"
out_path  <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/BRCA-EU/plot/largetd/snv"
snv_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/BRCA-EU/snv/sample_gene.txt"
allsv_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/BRCA-EU/shatterseek/sv_with_complex_annotation.txt"
pub_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.snv/BRCA-EU"
mutload_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/BRCA-EU/snv/mut_load.tsv"
### read
sv <- read.delim(sv_path)
allsv <- read.delim(allsv_path)
sample <- read.csv(sample_path,fileEncoding = "UTF-8-BOM")
sv$sample <- gsub(".sv.nocplx.type.bedpe","",sv$sample)
snv <- read.delim(snv_path,header = T,row.names = 1)
mutload <- read.delim(mutload_path)
# filter gene
snv <- snv[rownames(snv) %in% coding_gene$ensg,]
rownames(snv) <- as.character(ensg_gene_match[rownames(snv)])
snv_t <- data.frame(t(snv))
snv_t$sample <- row.names(snv_t)
# run
test_function(sv,allsv,snv_t,mutload,sample,out_path,mutcutoff=2,"BRCA-EU")
plot_function(sv,allsv,snv_t,sample,"PIK3CA",file.path(plotpath,"BRCA-EU.large.td.and.PIK3CA.pdf"),"FDR = 4.3e-13","BRCA-EU (n=532)",plotwide=40)

############# POG570 #####################
sample_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/pog570_sampleid_match.csv"
sv_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/real_and_random_withsig.txt"
out_path  <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/plot/largetd/snv"
snv_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/snv/sample_gene.txt"
allsv_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/shatterseek/sv_with_complex_annotation.txt"
pub_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.snv/POG570"
mutload_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/snv/mut_load.tsv"
### read
sv <- read.delim(sv_path)
allsv <- read.delim(allsv_path)
sample <- read.csv(sample_path,fileEncoding = "UTF-8-BOM")
sv$sample <- gsub(".sv.nocplx.type.bedpe","",sv$sample)
snv <- read.delim(snv_path,header = T,row.names = 1)
mutload <- read.delim(mutload_path)
# filter gene
snv <- snv[rownames(snv) %in% coding_gene$gene,]
snv_t <- data.frame(t(snv))
snv_t$sample <- gsub("X","",row.names(snv_t))
# run
test_function(sv,allsv,snv_t,mutload,sample,out_path,mutcutoff=2,"pancan")
test_function(sv,allsv,snv_t,mutload,sample[sample$histology_abbreviation=="BRCA",],
              out_path,mutcutoff=2,"BRCA")
plot_function(sv,allsv,snv_t,sample,"TP53",file.path(plotpath,"POG570.large.td.and.TP53.pdf"),"FDR = 0.015","POG570 Pan-cancer (n=570)",plotwide=40)
plot_function(sv,allsv,snv_t,sample[sample$histology_abbreviation=="BRCA",],"PIK3CA",file.path(plotpath,"POG570.BRCA.and.PIK3CA.pdf"),"FDR = 0.005","POG570 Breast (n=144)",plotwide=40)

############# PCAWG #####################
sample_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/pcawg2583_sampleid_match.csv"
sv_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/004_signature_matrix/pcawg_simple_sv_with_49svsig_simple.txt"
out_path  <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/008-largetd/snv"
snv_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/009-snv/sample_gene.txt"
allsv_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/002_starfish/pcawg_all_sv_starfish_annotation.txt"
pub_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.snv/PCAWG"
mutload_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/009-snv/mut_load.tsv"
### read
sv <- read.delim(sv_path)
allsv <- read.delim(allsv_path)
sample <- read.csv(sample_path,fileEncoding = "UTF-8-BOM")
sv$sample <- gsub(".sv.nocplx.type.bedpe","",sv$sample)
snv <- read.delim(snv_path,header = T,row.names = 1)
mutload <- read.delim(mutload_path)
# filter gene
snv <- snv[rownames(snv) %in% coding_gene$gene,]
snv_t <- data.frame(t(snv))
snv_t$sample <- gsub("X","",row.names(snv_t))
snv_t$sample <- gsub("\\.","-",snv_t$sample )
# run
test_function(sv,allsv,snv_t,mutload,sample,out_path,mutcutoff=2,"pancan")
test_function(sv,allsv,snv_t,mutload,sample[sample$histology_abbreviation=="Liver-HCC",],out_path,mutcutoff=2,"Liver-HCC")
test_function(sv,allsv,snv_t,mutload,sample[sample$histology_abbreviation=="Breast-AdenoCA",],out_path,mutcutoff=2,"Breast-AdenoCA")
test_function(sv,allsv,snv_t,mutload,sample[sample$histology_abbreviation=="Ovary-AdenoCA",],out_path,mutcutoff=2,"Ovary-AdenoCA")
test_function(sv,allsv,snv_t,mutload,sample[sample$histology_abbreviation=="Stomach-AdenoCA",],out_path,mutcutoff=2,"Stomach-AdenoCA")
test_function(sv,allsv,snv_t,mutload,sample[sample$histology_abbreviation=="Eso-AdenoCA",],out_path,mutcutoff=2,"Eso-AdenoCA")
test_function(sv,allsv,snv_t,mutload,sample[sample$histology_abbreviation=="Skin-Melanoma",],out_path,mutcutoff=2,"Skin-Melanoma")
test_function(sv,allsv,snv_t,mutload,sample[sample$histology_abbreviation=="Uterus-AdenoCA",],out_path,mutcutoff=2,"Uterus-AdenoCA")
test_function(sv,allsv,snv_t,mutload,sample[sample$histology_abbreviation=="Prost-AdenoCA",],out_path,mutcutoff=2,"Prost-AdenoCA")
# plot
plot_function(sv,allsv,snv_t,sample,"TP53",file.path(plotpath,"PCAWG.pancaner.large.td.and.TP53.pdf"),"FDR = 3.0e-17","PCAWG Pan-cancer (n=2,583)",plotwide=40)
plot_function(sv,allsv,snv_t,sample[sample$histology_abbreviation=="Breast-AdenoCA",],"PIK3CA",file.path(plotpath,"PCAWG.Breast-AdenoCA.large.td.and.PIK3CA.pdf"),"FDR = 0.053","PCAWG Breast (n=195)",plotwide=40)
plot_function(sv,allsv,snv_t,sample[sample$histology_abbreviation=="Uterus-AdenoCA",],c("ARID1A","PTEN"),file.path(plotpath,"PCAWG.Uterus-AdenoCA.large.td.and.ARID1A.PTEN.pdf"),"FDR\nARID1A = 0.025\nPTEN = 0.034","PCAWG Uterus (n=44)",plotwide=40)
plot_function(sv,allsv,snv_t,sample[sample$histology_abbreviation=="Ovary-AdenoCA",],"CDK12",file.path(plotpath,"PCAWG.Ovary-AdenoCA.large.td.and.CDK12.pdf"),"FDR = 7.5e-4","PCAWG Ovary (n=110)",plotwide=40)
plot_function(sv,allsv,snv_t,sample[sample$histology_abbreviation=="Prost-AdenoCA",],"SPOP",file.path(plotpath,"PCAWG.Prost-AdenoCA.large.td.and.SPOP.pdf"),"FDR = 0.0025","PCAWG Prostate (n=199)",plotwide=40)

############# Hartwig #####################
sample_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/hartwig2367_sampleid_match.csv"
sv_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/real_and_random_withsig.txt"
out_path  <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/plot/largetd/snv"
snv_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/snv/sample_gene.txt"
allsv_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/shatterseek/sv_with_complex_annotation.txt"
pub_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.snv/Hartwig"
mutload_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/snv/mut_load.tsv"
### read
sv <- read.delim(sv_path)
allsv <- read.delim(allsv_path)
sample <- read.csv(sample_path,fileEncoding = "UTF-8-BOM")
sv$sample <- gsub(".sv.nocplx.type.bedpe","",sv$sample)
snv <- read.delim(snv_path,header = T)
mutload <- read.delim(mutload_path)
# filter gene
snv <- snv[snv$gene != "" & snv$gene !="UGT2A1",]
rownames(snv) <-  snv$gene
snv <- snv[,names(snv)!="gene"]
snv <- snv[rownames(snv) %in% coding_gene$gene,]
snv_t <- data.frame(t(snv))
snv_t$sample <- row.names(snv_t)
# run
test_function(sv,allsv,snv_t,mutload,sample,out_path,mutcutoff=2,"pancan")
test_function(sv,allsv,snv_t,mutload,sample[sample$histology_abbreviation ==  "Breast",],out_path,mutcutoff=2,"breast")
test_function(sv,allsv,snv_t,mutload,sample[sample$histology_abbreviation ==  "Ovary",],out_path,mutcutoff=2,"ovary")
test_function(sv,allsv,snv_t,mutload,sample[sample$histology_abbreviation ==  "Prostate",],out_path,mutcutoff=2,"prostate")
test_function(sv,allsv,snv_t,mutload,sample[sample$histology_abbreviation ==  "Urothelial tract",],out_path,mutcutoff=2,"urothelial-tract")
test_function(sv,allsv,snv_t,mutload,sample[sample$histology_abbreviation ==  "Esophagus",],out_path,mutcutoff=2,"esophagus")
test_function(sv,allsv,snv_t,mutload,sample[sample$histology_abbreviation ==  "Uterus",],out_path,mutcutoff=2,"uterus")
# plot
plot_function(sv,allsv,snv_t,sample,"TP53",file.path(plotpath,"Hartiwg.pancaner.large.td.and.TP53.pdf"),"FDR = 1.5e-16","Hartiwg Pan-cancer (n=2,367)",plotwide=40)
plot_function(sv,allsv,snv_t,sample[sample$histology_abbreviation ==  "Breast",],"PIK3CA",file.path(plotpath,"Hartiwg.Breast.large.td.and.PIK3CA.pdf"),"FDR = 9.1e-5","Hartiwg Breast (n=471)",plotwide=40)
plot_function(sv,allsv,snv_t,sample[sample$histology_abbreviation ==  "Prostate",],"CDK12",file.path(plotpath,"Hartiwg.Prostate.large.td.and.CDK12.pdf"),"FDR = 1.9e-7","Hartiwg Prostate (n=206)",plotwide=40)
plot_function(sv,allsv,snv_t,sample[sample$histology_abbreviation ==  "Uterus",],c("ARID1A","PTEN"),file.path(plotpath,"Hartiwg.Uterus.large.td.and.ARID1A.PTEN.pdf"),"FDR\nARID1A = 0.029\nPTEN = 0.036","Hartiwg Uterus (n=45)",plotwide=40)
