### library
library(dplyr)
library(ggplot2)
require("ggrepel")
library(data.table)
library("doParallel") 
library("foreach") 
library(ggpubr)
library(stringr)

### function
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

tau.function <- function(drugi,myalternative){
  test = "Kendall"
  tau.test <- cor.test(largetd.and.drug.auc$largetd,largetd.and.drug.auc[,drugi], method="kendall", alternative = myalternative)
  tau.p <- tau.test$p.value
  tau <- tau.test$estimate
  tao.test.df <- data.frame(
    "Drug"=drugi,
    "Test"=test,
    "Tau.p"=tau.p,
    "Tau"=tau
  )
  return(tao.test.df)
}
rho.function <- function(drugi,myalternative){
  # Spearman rank correction test
  test = "Spearman"
  rho.test <- cor.test(largetd.and.drug.auc$largetd,largetd.and.drug.auc[,drugi], method="spearman", alternative = myalternative)
  rho.p <- rho.test$p.value
  rho <- rho.test$estimate
  rho.test.df <- data.frame(
    "Drug"=drugi,
    "Test"=test,
    "Rho.p"=rho.p,
    "Rho"=rho
  )
  return(rho.test.df)
}
t.function <- function(drugi,myalternative){
  # T-test
  test = "Student's t-Test"
  t.test <- t.test(largetd.and.drug.auc[largetd.and.drug.auc$withlargetd==2,drugi],
                   largetd.and.drug.auc[largetd.and.drug.auc$withlargetd==0,drugi],
                   alternative = myalternative)
  t.p <- t.test$p.value
  x.estimate <- t.test$estimate[1]
  y.estimate <- t.test$estimate[2]
  t.test.df <- data.frame(
    "Drug"=drugi,
    "Test"=test,
    "T.p"=t.p,
    "T"=x.estimate-y.estimate,
    "x.estimate"=x.estimate,
    "y.estimate"=y.estimate
  )
  return(t.test.df)
}

volcano_plot_function <- function(test.df_plot=t.test.df,
                                  p.cutoff=0.05,
                                  estimate.cutoff=0.1){
  test.df_plot$color <- ifelse(test.df_plot[,3] < p.cutoff,
                               ifelse(test.df_plot[,4] < -estimate.cutoff, "neg",
                                      ifelse(test.df_plot[,4] > estimate.cutoff, "pos","notsig")),
                               "notsig")
  test.df_plot$annotation <- ifelse(test.df_plot[,3] < 0.025,
                                    ifelse(test.df_plot[,4] < -estimate.cutoff,
                                           paste0(test.df_plot$Drug,"(",test.df_plot$moa,")"),
                                           ifelse(test.df_plot[,4] > estimate.cutoff, paste0(test.df_plot$Drug,"(",test.df_plot$moa,")"),"")
                                    ),
                                    ""
  )
  test.df_plot$annotation.simple <- ifelse(test.df_plot[,3] < 0.025,
                                    ifelse(test.df_plot[,4] < -estimate.cutoff,
                                           firstup(paste0(test.df_plot$Drug)),
                                           ifelse(test.df_plot[,4] > estimate.cutoff, firstup(paste0(test.df_plot$Drug)),"")
                                    ),
                                    ""
  )
  estimate.max <- ceiling(max(abs(test.df_plot[,4]))/0.05)*0.05
  
  p <- ggplot(
    test.df_plot, aes(x = test.df_plot[,4], y = -log10(test.df_plot[,3]),color=color)) +
    geom_point(alpha=0.8, size=0.3) +
    scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
    scale_x_continuous(breaks = c(-estimate.max,0,estimate.max),
                       labels = c(-estimate.max,0,estimate.max),
                       limits = c(-estimate.max,estimate.max))+
    geom_vline(xintercept=c(-estimate.cutoff,estimate.cutoff),lty=4,col="black",lwd=0.1189) +
    geom_hline(yintercept = -log10(p.cutoff),lty=4,col="black",lwd=0.1189) +
    geom_text_repel(
      aes(label = annotation.simple),
      box.padding = 0.1,
      point.padding = 1e-26,
      max.overlaps = 20,
      size = 2,
      segment.size=0.05
    )+
    
    labs(x=names(test.df_plot)[4],
         y="-log10 (p-value)")+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5,size=6), 
          axis.title = element_text(size=6,color="black"),
          axis.text = element_text(size=6,colour = "black"),
          axis.ticks.length = unit(0.5,"mm"),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          legend.position="none", 
          legend.title = element_blank())
  return(p)
}
p.adjust_function <- function(df,col=3){
  df$p.bfrn <- p.adjust(df[,col],method = "bonferroni",n=nrow(df))
  df$p.fdr <- p.adjust(df[,col],method = "BH")
  return(df)
}

### path
ccle.drug.auc.path <- "D:/Data/CCLE/raw.data.of.10.1038/PRISM.drug.repurposing.screen/secondary-screen-dose-response-curve-parameters.csv"
largetd.path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/CCLE/plot/sv/ccle_sample_largetd_annotation.csv"
plot.path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/ccle/drug"
moa.interested.path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/ccle/moa.interested.txt"

### read
ccle.drug.auc <- read.csv(ccle.drug.auc.path) %>% filter(is.na(ccle_name)==FALSE)  %>% as.data.frame()
ccle.drug.auc.HTS002 <- ccle.drug.auc %>% filter(screen_id == "HTS002") %>% as.data.frame()
ccle.drug.auc.MTS006 <- ccle.drug.auc %>% filter(screen_id == "MTS006") %>% as.data.frame()
ccle.drug.auc.MTS010 <- ccle.drug.auc %>% filter(screen_id == "MTS010") %>% as.data.frame()
ccle.drug.auc.test <- ccle.drug.auc.HTS002

largetd <- read.csv(largetd.path,row.names = 1)
moa.interested <- read.delim(moa.interested.path,header=F)
drug.of.moa.interested <- unique(ccle.drug.auc.test[ccle.drug.auc.test$moa %in% unique(moa.interested$V1),c("name","moa")])
# write.table(drug.of.moa.interested,
#        file="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/ccle/drug.of.moa.interested.txt",
#        sep = "\t",row.names = F,quote = F)

### drug file, from long to wide
cellline <- unique(ccle.drug.auc.test$ccle_name)
drug <- unique(ccle.drug.auc.test$name)
cellline.drug.df <- data.frame(matrix(data=NA,nrow = length(cellline),ncol = length(drug)))
rownames(cellline.drug.df) <- cellline
colnames(cellline.drug.df) <- drug
for (i in 1:nrow(ccle.drug.auc.test)){
  celli <- ccle.drug.auc.test$ccle_name[i]
  drugi <- ccle.drug.auc.test$name[i]
  auci <- ccle.drug.auc.test$auc[i]
  cellline.drug.df[celli,drugi] <- auci
}

### merge auc and largetd
#Higher AUC implies less sensitivity to the compound
# ccle.drug.auc %>% filter(moa=="PARP inhibitor") %>% select(name) %>% unique() # PARP inhibitor
largetd.and.drug.auc <- merge(largetd,cellline.drug.df,by = 'row.names', all.x = TRUE)

### Test
system.time({
  cl<- makeCluster(8)      
  registerDoParallel(cl)
  tau.test.df <- foreach(drugi = drug,
                                .combine = rbind) %dopar% tau.function(drugi,"two.sided")
  
  stopCluster(cl)
})
system.time({
  cl<- makeCluster(8)      
  registerDoParallel(cl)
  rho.test.df <- foreach(drugi = drug,
                         .combine = rbind) %dopar% rho.function(drugi,"two.sided")
  
  stopCluster(cl)
})
system.time({
  cl<- makeCluster(8)      
  registerDoParallel(cl)
  t.test.df <- foreach(drugi = drug,
                         .combine = rbind) %dopar% t.function(drugi,"two.sided")
  
  stopCluster(cl)
})

### add moa
tau.test.df <- merge(tau.test.df,unique(ccle.drug.auc[,c("name","moa")]),by.x = "Drug",by.y = "name",all.x = T)
write.table(tau.test.df,file.path(plot.path,"Kendall.rank.test.txt"),row.names = F,quote = F,sep = "\t")
rho.test.df <- merge(rho.test.df,unique(ccle.drug.auc[,c("name","moa")]),by.x = "Drug",by.y = "name",all.x = T)
write.table(rho.test.df,file.path(plot.path,"Spearman.rank.test.txt"),row.names = F,quote = F,sep = "\t")
t.test.df <- merge(t.test.df,unique(ccle.drug.auc[,c("name","moa")]),by.x = "Drug",by.y = "name",all.x = T)
write.table(t.test.df,file.path(plot.path,"Ttest.test.txt"),row.names = F,quote = F,sep = "\t")

### SUB FDR
sub <- c("twosides","greater","less")[1]
drug.of.moa.interested.path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/publication/plot/ccle/drug.of.moa.interested2.txt"

drug.of.moa.interested <- read.delim(drug.of.moa.interested.path)
tau.test.df <- read.delim(file.path(plot.path,sub,"Kendall.rank.test.txt"))
rho.test.df <- read.delim(file.path(plot.path,sub,"Spearman.rank.test.txt"))
t.test.df <- read.delim(file.path(plot.path,sub,"Ttest.test.txt"))

tau.test.df.sub <- tau.test.df[tau.test.df$Drug %in% drug.of.moa.interested$name,]
rho.test.df.sub <- rho.test.df[rho.test.df$Drug %in% drug.of.moa.interested$name,]
t.test.df.sub <- t.test.df[t.test.df$Drug %in% drug.of.moa.interested$name,]

tau.test.df.sub <- p.adjust_function(tau.test.df.sub,col=3)
rho.test.df.sub <- p.adjust_function(rho.test.df.sub,col=3)
t.test.df.sub <- p.adjust_function(t.test.df.sub,col=3)

write.table(tau.test.df.sub,file.path(plot.path,sub,"Kendall.rank.test.interested.txt"),row.names = F,quote = F,sep = "\t")
write.table(rho.test.df.sub,file.path(plot.path,sub,"Spearman.rank.test.interested.txt"),row.names = F,quote = F,sep = "\t")
write.table(t.test.df.sub,file.path(plot.path,sub,"Ttest.test.interested.txt"),row.names = F,quote = F,sep = "\t")

### Number calculator
sum(unique(largetd.and.drug.auc$Row.names) %in% ccle.drug.auc.test$ccle_name)
# [1] 203 cell lines with large TD are in the screen
length(unique(ccle.drug.auc.test$ccle_name))
# [1] 480 cell lines in the screen
length(unique(ccle.drug.auc.test$name))
# [1] 1394 drugs in the screen

### plot volcano
sub <- c("twosides","greater","less")[1]
tau.test.df <- read.delim(file.path(plot.path,sub,"Kendall.rank.test.interested.txt"))
rho.test.df <- read.delim(file.path(plot.path,sub,"Spearman.rank.test.interested.txt"))
t.test.df <- read.delim(file.path(plot.path,sub,"Ttest.test.interested.txt"))

p <- volcano_plot_function (test.df_plot=tau.test.df,p.cutoff=0.05,estimate.cutoff=0.1)
ggsave(plot=p,file.path(plot.path,"Kendall.rank.test.volcano.plot.pdf"),width = 34,height = 34,units = "mm")
p <- volcano_plot_function (test.df_plot=rho.test.df,p.cutoff=0.05,estimate.cutoff=0.1)
ggsave(plot=p,file.path(plot.path,"Spearman.rank.test.volcano.plot.pdf"),width = 34,height = 34,units = "mm")
p <- volcano_plot_function (test.df_plot=t.test.df,p.cutoff=0.05,estimate.cutoff=0.1)
ggsave(plot=p,file.path(plot.path,"Ttest.test.volcano.plot.pdf"),width = 34,height = 34,units = "mm")

### plot large td and drug
largetd.and.drug.auc$withlargetd <- factor(largetd.and.drug.auc$withlargetd,levels = c(0,1,2))
mydrug <- c("bendamustine",
            "famciclovir",
            "MK-1775",
            "nelarabine",
            "PJ-34"
)

for (mydrugi in mydrug){
  plot.name <- paste0("largetd.and.",mydrugi,".group.pdf")
  ggplot(data=largetd.and.drug.auc,aes(x=withlargetd,y=largetd.and.drug.auc[,mydrugi]))+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(width = 0.2,size=0.1)+
    scale_x_discrete(labels = c('<=5','6-15','>=16'))+
    ggtitle(firstup(mydrugi))+
    ylab("AUC")+
    xlab("Large TD")+
    theme_bw()+
    theme(axis.title.x = element_text(size=6,color = "black"),
          axis.title.y = element_text(size=6,color = "black"),
          axis.text.x = element_text(size=6,color = "black"),
          axis.text.y = element_text(size=6,color = "black"),
          plot.title = element_text(hjust = 0.5,size=6,color = "black"),
          panel.background = element_blank(),
          panel.grid = element_blank())
  ggsave(file.path(plot.path,plot.name),width = 50,height = 50,units = "mm")
}

for (mydrugi in mydrug){
  plot.name <- paste0("largetd.and.",mydrugi,".scatter.pdf")
  # ggplot(data=largetd.and.drug.auc,aes(x=largetd,y=largetd.and.drug.auc[,mydrugi]))+
  #   geom_point()+
  mydrugi.max <- max(largetd.and.drug.auc[,mydrugi],na.rm = T)
  
  ggscatter(largetd.and.drug.auc, x = "largetd", y = mydrugi,
              add = "reg.line",  # Add regressin line
              add.params = list(color = "blue", fill = "lightgray",size=0.1189), # Customize reg. line
              conf.int = TRUE, # Add confidence interval
            size=0.001
    ) +
    stat_cor(method = "kendall",
             alternative = "two.sided",
             cor.coef.name = "tau",
             label.x = 0, label.y = mydrugi.max*4.8/5,
             size = 2)+

    ggtitle(firstup(mydrugi))+
    ylab("AUC")+
    xlab("Number of large TDs")+
    theme_bw()+
    theme(axis.title.x = element_text(size=6,color = "black"),
          axis.title.y = element_text(size=6,color = "black"),
          axis.text.x = element_text(size=6,color = "black"),
          axis.text.y = element_text(size=6,color = "black"),
          axis.ticks.length = unit(0.5,"mm"),
          plot.title = element_text(hjust = 0.5,size=6,color = "black"),
          panel.background = element_blank(),
          panel.grid = element_blank())
  ggsave(file.path(plot.path,plot.name),width = 35,height = 35,units = "mm")
}
