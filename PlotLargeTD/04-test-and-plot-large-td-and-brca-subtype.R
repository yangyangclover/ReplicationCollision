### library
library(dplyr)
library(ggplot2)
library(ggpubr)

### path
sample_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/pcawg2583_sampleid_match.csv"
sv_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/004_signature_matrix/pcawg_simple_sv_with_49svsig_simple.txt"
idmatch_path <- "D:/Data/PCAWG/sample_info/pcawg_summary.tsv"
tcga.uuid2barcode_path <- "D:/Data/PCAWG/sample_info/pc_annotation-tcga_uuid2barcode.tsv"
brca.subtype_path <- "D:/Data/PCAWG/published/brca.subtype/41467_2021_27079_MOESM5_ESM.csv"
brca.subtype_path2 <- "D:/Data/PCAWG/published/brca.subtype/1-s2.0-S0092867415011952-mmc2.csv"
plot_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/large.td.and.brcasubtype"

### run
sv <- read.delim(sv_path)
sample <- read.csv(sample_path,fileEncoding = "UTF-8-BOM")
tcga.uuid2barcode <- read.delim(tcga.uuid2barcode_path)
idmatch <- read.delim(idmatch_path,fileEncoding = "UTF-8-BOM")
idmatch <- idmatch[,c("submitter_donor_id","tumor_wgs_aliquot_id","dcc_project_code")]
newidmatch <- data.frame()
for (i in 1:nrow(idmatch)){
  newidmatch <- rbind(newidmatch,data.frame(
    "submitter_donor_id"=idmatch[i,1],
    "tumor_wgs_aliquot_id"=unlist(strsplit(idmatch[i,2],",")),
    "dcc_project_code"=idmatch[i,3]
  ))
}
brca.subtype <- read.csv(brca.subtype_path,fileEncoding = "UTF-8-BOM")
brca.subtype2 <- read.csv(brca.subtype_path2,fileEncoding = "UTF-8-BOM")

sample <- merge(sample,idmatch[,c("tumor_wgs_aliquot_id","dcc_project_code")],by.x = "WGS",by.y = "tumor_wgs_aliquot_id",all.x = T)
all_sample <- data.frame(unique(sample[,c("WGS","histology_abbreviation")]))
names(all_sample)[1] <- "sample"
### large td
sample_largetd <- sv %>% filter(sig_group %in% paste0("td_",seq(10,16,1)))%>% group_by(sample) %>% dplyr::count()
names(sample_largetd)[2] <- "largetd"
all_sample <- merge(all_sample,sample_largetd,
                    by.x = "sample",by.y = "sample",all.x = T)
all_sample$largetd <- ifelse(is.na(all_sample$largetd),0,all_sample$largetd)

### add tcgaid
all_sample <- merge(all_sample,newidmatch,by.x = "sample",by.y = "tumor_wgs_aliquot_id",all.x = T)
### add tcga barcode
all_sample <- merge(all_sample,tcga.uuid2barcode[,c("barcode","uuid")],by.x="submitter_donor_id",by.y = "uuid",all.x = T)
### add brca subtype
all_sample <- merge(all_sample,brca.subtype[,c("Sample","Subtype")],by.x = "barcode",by.y = "Sample",all.x = T)
### add brca subtype2
all_sample <- merge(all_sample,brca.subtype2[,c("Case.ID","PAM50","Cell.cycle.score")],by.x = "barcode",by.y = "Case.ID",all.x = T)

# Subtype   PAM50     n
# <chr>     <chr> <int>
#   1 HER2+     Her2      5
# 2 HER2+     LumA      3
# 3 HER2+     LumB      6
# 4 HER2+     NA        1
# 5 HR-/HER2- Basal    30
# 6 HR-/HER2- Her2      1
# 7 HR-/HER2- NA        5
# 8 HR+/HER2- Basal     2
# 9 HR+/HER2- Her2      1
# 10 HR+/HER2- LumA      7
# 11 HR+/HER2- LumB      8
# 12 HR+/HER2- NA        4
# 13 NA        Basal     4
# 14 NA        Her2      5
# 15 NA        LumA      1
# 16 NA        LumB      2

### plot and stat
all_sample_subtype <- all_sample[is.na(all_sample$PAM50)==F,]
all_sample_subtype$PAM50 <- factor(all_sample_subtype$PAM50,levels = c("Basal","Her2","LumA","LumB"))


###plot
all_sample_subtype$largetd <- ifelse(all_sample_subtype$largetd>0,all_sample_subtype$largetd,all_sample_subtype$largetd + 0.1)
# add sample number
subtype.count <- all_sample_subtype %>% group_by(PAM50) %>% count() %>% as.data.frame()
all_sample_subtype <- merge(all_sample_subtype,subtype.count)
all_sample_subtype$subtype.num <- paste0(all_sample_subtype$PAM50,"\n","(",all_sample_subtype$n,")")
# my comparison
my_comparisons <- list( c(unique(all_sample_subtype$subtype.num)[1],unique(all_sample_subtype$subtype.num)[2]),
                        c(unique(all_sample_subtype$subtype.num)[1],unique(all_sample_subtype$subtype.num)[3]),
                        c(unique(all_sample_subtype$subtype.num)[1],unique(all_sample_subtype$subtype.num)[4]),
                        c(unique(all_sample_subtype$subtype.num)[2],unique(all_sample_subtype$subtype.num)[3]),
                        c(unique(all_sample_subtype$subtype.num)[2],unique(all_sample_subtype$subtype.num)[4]),
                        c(unique(all_sample_subtype$subtype.num)[3],unique(all_sample_subtype$subtype.num)[4])
)
# plot
ggboxplot(all_sample_subtype, x = "subtype.num", y = "largetd",
          color = "subtype.num", palette = "jco",size=0.5,width = 0.5)+ 
  stat_compare_means(comparisons = my_comparisons,size=2)+ # Add pairwise comparisons p-value
  stat_compare_means(
   label.y = log10(30000),
                     size=2)+    # Add global p-value
  ylab("Number of Large TDs")+
  scale_y_log10()+
  theme(axis.text = element_text(size=6,color="black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6,color="black"),
        legend.position = "none")
ggsave(file.path(plot_path,"PCAWG.BRCA.subtype.log10largtd.twosides.pdf"),width = 3,height = 2.5)


ggboxplot(all_sample_subtype, x = "PAM50", y = "largetd",
          color = "PAM50", palette = "jco",size=0.5,width = 0.5)+ 
  stat_compare_means(comparisons = my_comparisons,size=2)+ # Add pairwise comparisons p-value
  stat_compare_means(
    label.y = 500,
    size=2)+    # Add global p-value
  ylab("Number of Large TD")+
  theme(axis.text = element_text(size=6,color="black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6,color="black"),
        legend.position = "none")
ggsave(file.path(plot_path,"PCAWG.BRCA.subtype.largtd.twosides.pdf"),width = 3,height = 2.5)

ggboxplot(all_sample_subtype, x = "PAM50", y = "largetd",
          color = "PAM50", palette = "jco",size=0.5,width = 0.5)+ 
  stat_compare_means(comparisons = my_comparisons,size=2,method.args = list(alternative = "greater"))+ # Add pairwise comparisons p-value
  stat_compare_means(
    label.y = log10(30000),
    size=2)+    # Add global p-value
  ylab("Number of Large TD")+
  scale_y_log10()+
  theme(axis.text = element_text(size=6,color="black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6,color="black"),
        legend.position = "none")
ggsave(file.path(plot_path,"PCAWG.BRCA.subtype.log10largtd.greater.pdf"),width = 3,height = 2.5)

ggboxplot(all_sample_subtype, x = "PAM50", y = "largetd",
          color = "PAM50", palette = "jco",size=0.5,width = 0.5)+ 
  stat_compare_means(comparisons = my_comparisons,size=2,method.args = list(alternative = "greater"))+ # Add pairwise comparisons p-value
  stat_compare_means(
    label.y = 500,
    size=2)+    # Add global p-value
  ylab("Number of Large TD")+
  theme(axis.text = element_text(size=6,color="black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=6,color="black"),
        legend.position = "none")
ggsave(file.path(plot_path,"PCAWG.BRCA.subtype.largtd.greater.pdf"),width = 3,height = 2.5)

