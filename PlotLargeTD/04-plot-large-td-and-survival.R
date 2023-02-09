# 2022-09-27
# must have time and death column
# draw_data was used to plot survival plot:
# sample                                  tumor           death   time      have_td
# 0009b464-b376-4fbc-8a56-da538269a02f    Ovary-AdenoCA   1       1972      Have_TD


library(dplyr)
library(survival)
library(survminer)
library(cowplot)
library(egg)



single_psurvival_risk_cumevents_table_theme <- theme(
  axis.line.x=element_line(colour="black", size=0.1189),
  axis.line.y=element_line(colour="black", size=0.1189),
  axis.ticks.x=element_line(colour="black", size=0.1189),
  axis.ticks.y=element_line(colour="black", size=0.1189),
  axis.ticks.length.y.left=unit(0.706, "mm"),
  axis.ticks.length.x.bottom=unit(0.706, "mm"),
  axis.text.x=element_text(colour="black", size=unit(6, "picas")),
  # axis.text.y=element_text(colour="black", size=unit(6, "picas")),
  axis.title.x=element_text(colour="black", size=unit(6, "picas")),
  axis.title.y=element_text(colour="black", size=unit(6, "picas")),
  plot.title=element_text(colour="black", size=unit(6, "picas"))
)


draw_survial_plot <- function(input_data, input_survfit, input_category_column_name, input_category_order, input_color_vector, p_value, output_file_name){
  # input_color_vector is in the order of input_category_order
  category_order=input_category_order
  legendLabel=paste(category_order, "\n(N=", table(input_data[[input_category_column_name]])[category_order], ")", sep="")
  psurvival=ggsurvplot(input_survfit, data=input_data, font.main=c(2),
                       palette=input_color_vector,
                       font.x = c(6),
                       font.y = c(6),
                       font.tickslab = c(6),font.legend=c(6),legend.title = "",
                       legend = c("bottom"),
                       censor.size=1,
                       size=0.2,
                       pval.size=2,
                       fontsize=2,
                       legend.labs=legendLabel,
                       pval=p_value,
                       xlab="Days",
                       surv.plot.height=0.9,
                       risk.table=TRUE,
                       cumevents=TRUE)
  psurvival$plot=psurvival$plot+theme(axis.line.x=element_line(colour="black", size=0.1189),
                                      axis.line.y=element_line(colour="black", size=0.1189),
                                      axis.ticks.x=element_line(colour="black", size=0.1189),
                                      axis.ticks.y=element_line(colour="black", size=0.1189),
                                      axis.ticks.length.y.left=unit(0.706, "mm"),
                                      axis.ticks.length.x.bottom=unit(0.706, "mm"))+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))
  psurvival$table=psurvival$table+single_psurvival_risk_cumevents_table_theme
  psurvival$cumevents=psurvival$cumevents+single_psurvival_risk_cumevents_table_theme
  pdf(output_file_name, width=5, height=4, onefile=FALSE)
  print(psurvival)
  dev.off()
}


my_ggsurvplot_facet<-function(pval.size = 5, ...)
{
  newcall <- bquote(
    p <- p + geom_text(data = pvals.df, aes(x = pval.x, y = pval.y, 
                                            label = pval.txt), size = .(pval.size), hjust = 0)
  )
  
  body(ggsurvplot_facet)[[20]][[3]][[8]] <- newcall
  ggsurvplot_facet(...)
}

draw_survial_plot_with_facet <- function(input_data, input_survfit, input_category_column_name, input_fact_column_name, input_category_order, input_color_vector, p_value, output_file_name){
  # input_color_vector is in the order of input_category_order
  category_order=input_category_order
  legendLabel=paste(category_order, "\n(N=", table(input_data[[input_category_column_name]])[category_order], ")", sep="")
  # fail to change axis.line size as draw_survial_plot
  # legendLable is for all category, not for each facet
  psurvival=my_ggsurvplot_facet(input_survfit, data=input_data, facet.by=input_fact_column_name,
                                font.main=c(2),
                                palette=input_color_vector,
                                font.x = c(6),
                                font.y = c(6),
                                font.tickslab = c(6),
                                font.legend=c(6),
                                legend.title="",
                                legend="bottom",
                                censor.size=3.5,
                                size=0.5,
                                pval=TRUE,
                                pval.size=2,
                                legend.labs=legendLabel,
                                xlab="Days",
                                surv.plot.height=0.9,
                                panel.labs.font=list(size=6),
                                panel.labs.background=list(color=NA, fill=NA),
                                short.panel.labs=TRUE)
  psurvival$theme=psurvival$theme+theme(axis.line.x=element_line(colour="black", size=0.1189),
                                        axis.line.y=element_line(colour="black", size=0.1189),
                                        axis.ticks.x=element_line(colour="black", size=0.1189),
                                        axis.ticks.y=element_line(colour="black", size=0.1189),
                                        axis.ticks.length.y.left=unit(0.706, "mm"),
                                        axis.ticks.length.x.bottom=unit(0.706, "mm"))
  
  temp_surv_formula_string=as.formula(paste("Surv(time=time, event=death)~ ",input_category_column_name, " + ", input_fact_column_name, sep=""))
  fit_for_risk_and_cumevents <- surv_fit(temp_surv_formula_string, data=input_data)
  # bug: there is bug when there are more than 24 variable (https://github.com/kassambara/survminer/issues/526 and https://github.com/tidyverse/ggplot2/issues/4409), 
  # use palette to solve the bug
  facet_number=length(unique(input_data[[input_fact_column_name]]))
  risk_and_cumevents_table<-ggsurvplot(fit_for_risk_and_cumevents, data=input_data, risk.table=T, cumevents=T, palette=rep("black", facet_number*2))
  risk_plot<-ggplot(risk_and_cumevents_table$table$data, aes_string("time", input_category_column_name))+
    geom_text(aes(label = n.risk), colour="black", size=2)+
    facet_wrap(as.formula(paste("~", input_fact_column_name, sep="")))+ylab("Number at risk")+
    facet_psurvival_risk_cumevents_table_theme
  cumevent_plot<-ggplot(risk_and_cumevents_table$cumevents$data, aes_string("time", input_category_column_name))+
    geom_text(aes(label = n.event), colour="black", size=2)+
    facet_wrap(as.formula(paste("~", input_fact_column_name, sep="")))+
    ylab("Cumulative number of events")+facet_psurvival_risk_cumevents_table_theme
  
  merge_plot<-ggarrange(psurvival, risk_plot, cumevent_plot, heights = c(3,1,1), ncol=1)
  
  pdf(output_file_name, width=10, height=18, onefile=FALSE) # widht=8, height=7 when no tables
  print(merge_plot)
  dev.off()
}

draw_survial_plot_with_facet_and_bar_plot <- function(input_data, input_survfit, input_category_column_name, input_fact_column_name, input_category_order, input_color_vector, p_value, output_file_name){
  # input_color_vector is in the order of input_category_order
  category_order=input_category_order
  legendLabel=paste(category_order, "\n(N=", table(input_data[[input_category_column_name]])[category_order], ")", sep="")
  # fail to change axis.line size as draw_survial_plot
  # legendLable is for all category, not for each facet
  psurvival=my_ggsurvplot_facet(input_survfit, data=input_data, facet.by=input_fact_column_name,ncol=4,
                                font.main=c(2),
                                palette=input_color_vector,
                                font.x = c(6),
                                font.y = c(6),
                                font.tickslab = c(6),
                                font.legend=c(6),
                                legend.title="",
                                legend="bottom",
                                censor.size=3.5,
                                size=0.5,
                                pval=TRUE,
                                pval.size=2,
                                legend.labs=legendLabel,
                                xlab="Days",
                                surv.plot.height=0.9,
                                panel.labs.font=list(size=6),
                                panel.labs.background=list(color=NA, fill=NA),
                                short.panel.labs=TRUE)
  psurvival$theme=psurvival$theme+theme(axis.line.x=element_line(colour="black", size=0.1189),
                                        axis.line.y=element_line(colour="black", size=0.1189),
                                        axis.ticks.x=element_line(colour="black", size=0.1189),
                                        axis.ticks.y=element_line(colour="black", size=0.1189),
                                        axis.ticks.length.y.left=unit(0.706, "mm"),
                                        axis.ticks.length.x.bottom=unit(0.706, "mm"))
  
  sample_count<-as.data.frame(table(input_data[[input_category_column_name]],input_data[[input_fact_column_name]]))
  colnames(sample_count)=c("category", "facet", "count")
  #sample_count$label=paste(sample_count$category, "\n(N=", sample_count$count, ")", sep="")
  sample_count_plot<-ggplot(sample_count, aes(x=category, y=count))+geom_bar(stat="identity")+
    facet_wrap(~facet, ncol=4, scale="free_y")+
    theme(axis.text.x.bottom=element_text(angle=70, vjust=0.6))+
    ylab("Sample Count")+facet_psurvival_risk_cumevents_table_theme
  
  merge_plot<-ggarrange(psurvival, sample_count_plot, widths=c(1.3,1), ncol=2)
  
  pdf(output_file_name, width=8, height=7, onefile=FALSE) # height=2.5 for Non_PCAWG samples
  print(merge_plot)
  dev.off()
}


pairwise_comparison <- function(input_data, input_category_column_name, input_category_order, input_title){
  draw_data_reference <- input_data
  draw_data_reference[[input_category_column_name]] = factor(draw_data_reference[[input_category_column_name]], 
                                                             levels=input_category_order)
  fit.coxph <- coxph(as.formula(paste("Surv(time=time, event=death)~",input_category_column_name,sep="")), data=draw_data_reference)
  hazardPlot <- ggforest(fit.coxph, data=draw_data_reference,fontsize=0.7,
                         main=input_title)
  print(hazardPlot)
  pairwise_survdiff(as.formula(paste("Surv(time=time, event=death)~",input_category_column_name,sep="")), data=input_data, 
                    p.adjust.method="none")
  
}

# 1. PCAWG 
setwd("D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/xiaoming")
# 1. pcawg survival data
sur_pcawg=read.csv("PCAWG/pcawg_survival.csv",header=T,stringsAsFactors=F)
sur_pcawg$time=ifelse(!is.na(sur_pcawg$donor_survival_time),sur_pcawg$donor_survival_time,sur_pcawg$donor_interval_of_last_followup)
sur_pcawg$death=ifelse(sur_pcawg$donor_vital_status=="alive",0,ifelse(sur_pcawg$donor_vital_status=="deceased",1,NA))
sur_pcawg$sample=sur_pcawg$aliquot_id
sur_pcawg=sur_pcawg[c("sample","death","time")]

# 2. all analysis samples (2583 samples)
all_sample=read.csv("PCAWG/pcawg2583_sampleid_match.csv",header=T,stringsAsFactors=F)
colnames(all_sample)=c("sample", "RNA", "tumor")
all_sample=all_sample[,c("sample", "tumor")]

# 3. Large TD samples
large_td_data=read.delim("PCAWG/largetd_for_bedtools.bed", header=F, stringsAsFactors=F)
colnames(large_td_data)=c("chr", "start", "end", "sample")


# 4. TP53 mutation
tp53_maf <- read.delim("PCAWG/tp53.maf", header=T, sep="\t", stringsAsFactors=F)
# Variant_Classification: 3'UTR, Frame_Shift_Del, Frame_Shift_Ins, IGR, In_Frame_Del, 
#   In_Frame_Ins, Intron, Missense_Mutation, Nonsense_Mutation, Silent, 
# select Frame_Shift_Del, Frame_Shift_Ins, 
#   Missense_Mutation, Nonsense_Mutation, Splice_Site
tp53_maf <- tp53_maf[tp53_maf$Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins", 
                                                            "Missense_Mutation", "Nonsense_Mutation", "
                                                            "),]
tp53_mutated_samples=unique(tp53_maf$Tumor_Sample_Barcode)

# 5. merge data
merge_data=merge(all_sample, sur_pcawg, by.x="sample", by.y="sample", all.x=T)
merge_data$have_td="No_TD"
merge_data[merge_data$sample %in% large_td_data$sample,]$have_td="Have_TD"
merge_data$tp53="WT"
merge_data[merge_data$sample %in% tp53_mutated_samples,]$tp53="Mut"
# filter unsure data
merge_data=merge_data[!is.na(merge_data$time) & ! is.na(merge_data$death),]
merge_data=merge_data[merge_data$time>=0,] # time of 19233fd1-5229-466e-acf3-5882165758e0 in THCA and 93b51c61-6eea-4228-a102-840a2e118522 in Kidney-ChRCC less than 0
# Cervix-AdenoCA (0 No_TD, 2 Have_TD), Myeloid-MDS (2 No_TD, 0 Have_TD),  Myeloid-MPN (23 No_TD, 0 Have_TD)


category_order=c("No_TD", "Have_TD")
merge_data$have_td=factor(merge_data$have_td, levels=category_order)
color_vector=c("#7b94d1", "#ea6293")

# 6. thrshold
td_count=large_td_data %>% count(sample)
td_threshold=10
filter_have_td_sample=td_count[td_count$n<td_threshold,]$sample
filter_merge_data=merge_data[! merge_data$sample %in% filter_have_td_sample, ]
# 6.1. all samples
fit1 <- survfit(Surv(time=time, event=death) ~ have_td, data=filter_merge_data)
my_pvalue<-surv_pvalue(fit1)
draw_survial_plot(input_data=filter_merge_data, input_survfit=fit1, input_category_column_name="have_td", 
                  input_category_order=category_order, input_color_vector=color_vector, p_value=my_pvalue$pval, 
                  output_file_name=paste("pcawg_pancan_at_least_", td_threshold, "_large_td_both_TP53_mut_and_WT.pdf"))

# 6.2. PCAWG pancan TP53 Mut samples
temp_data=filter_merge_data[filter_merge_data$tp53=="Mut",]
fit1 <- survfit(Surv(time=time, event=death) ~ have_td, data=temp_data)
my_pvalue<-surv_pvalue(fit1)
draw_survial_plot(temp_data, fit1, "have_td", category_order, color_vector, my_pvalue$pval, paste("pcawg_pancan_at_least_", td_threshold, "_large_td_TP53_mut.pdf"))

# 6.3. TP53 WT samples
temp_data=filter_merge_data[filter_merge_data$tp53=="WT",]
fit1 <- survfit(Surv(time=time, event=death) ~ have_td, data=temp_data)
my_pvalue<-surv_pvalue(fit1)
draw_survial_plot(temp_data, fit1, "have_td", category_order, color_vector, my_pvalue$pval, paste("pcawg_pancan_at_least_", td_threshold, "_large_td_TP53_WT.pdf"))

# 6.4. PCAWG PRAD TP53 WT
temp_data=filter_merge_data[filter_merge_data$tp53=="WT" & filter_merge_data$tumor=="Prost-AdenoCA",]
fit1 <- survfit(Surv(time=time, event=death) ~ have_td, data=temp_data)
my_pvalue<-surv_pvalue(fit1)
draw_survial_plot(temp_data, fit1, "have_td", category_order, color_vector, my_pvalue$pval, paste("pcawg_Prost-AdenoCA_at_least_", td_threshold, "_large_td_TP53_WT.pdf"))

# 6.5. PCAWG Skin-Melanoma TP53 WT
temp_data=filter_merge_data[filter_merge_data$tp53=="WT" & filter_merge_data$tumor=="Skin-Melanoma",]
fit1 <- survfit(Surv(time=time, event=death) ~ have_td, data=temp_data)
my_pvalue<-surv_pvalue(fit1)
draw_survial_plot(temp_data, fit1, "have_td", category_order, color_vector, my_pvalue$pval, paste("pcawg_Skin-Melanoma_at_least_", td_threshold, "_large_td_TP53_WT.pdf"))

# 6.6. PCAWG Eso-AdenoCA TP53 mut
temp_data=filter_merge_data[filter_merge_data$tp53=="Mut" & filter_merge_data$tumor=="Eso-AdenoCA",]
fit1 <- survfit(Surv(time=time, event=death) ~ have_td, data=temp_data)
my_pvalue<-surv_pvalue(fit1)
draw_survial_plot(temp_data, fit1, "have_td", category_order, color_vector, my_pvalue$pval, paste("pcawg_Eso-AdenoCA_at_least_", td_threshold, "_large_td_TP53_mut.pdf"))

