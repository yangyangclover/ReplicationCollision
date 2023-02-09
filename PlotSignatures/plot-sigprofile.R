# library
library(reshape)
library(magrittr)
library(ggplot2)
library(dplyr)
library(ggh4x)

# functions
plot_sigprofiler_function <- function(phat_path,
                                      output_signature_dir,
                                      class_num,
                                      color_path,
                                      figure_name,
                                      signaturenum,
                                      signame_list,
                                      order){
  if (!dir.exists(output_signature_dir)) {
    dir.create(output_signature_dir, recursive = TRUE)
  }
  
  ####################### annotate row ##############################
  if  (class_num=="49"){
    yinter <- c(2,11,30)+0.5
  }
  ############################ READ FILES ##########################
  phat <- read.delim(phat_path,row.names = 1)
  colnames(phat) <- gsub(paste0("SBS",class_num),"",colnames(phat))
  
  
  ############################# plot signatures ########################
  phat$group  <- row.names(phat)
  # phat$group_fullname <- ifelse(phat$group %in% names(del_dup_rename),
  #                               del_dup_rename[phat$group],
  #                               phat$group)
  phat <- melt(phat, value.name="percent", variable.name="sig", na.rm=TRUE)
  phat_matrix_color<-read.csv(color_path,header=TRUE, fileEncoding="UTF-8-BOM")
  phat_color <- merge(phat,phat_matrix_color,by.x="group",by.y="sig")
  signew=dplyr::group_by(phat_color, group) %>% dplyr::mutate(percentage = value/sum(value))
  # signew=dplyr::group_by(phat_color, variable) %>% dplyr::mutate(percentage = value/sum(value))
  signew_data=as.data.frame(signew)
  
  signew_data_order <- signew_data[order(signew_data$order,signew_data$group), ]
  signew_data_order$svtype = factor(signew_data_order$svtype, levels = c("DEL","DUP","INV","TRA","C_A","C_G","C_T","T_A","T_C","T_G"))
  
  if  (class_num=="49"){
    signew_data_order$group_fullname = factor(
      signew_data_order$group_fullname,
      levels = rev(unique(phat_matrix_color$group_fullname))
    )
  }
  ################### plot signatures ###################
  if  (!is.na(signame_list)){
    names(signame_list) <- order
    signew_data_order$signame <- signame_list[as.character(signew_data_order$variable)]
    signew_data_order$signame = factor(signew_data_order$signame,levels = signame_list)
  } else{
    signew_data_order$signame <- signew_data_order$variable
  }
  
  
  
  signature_plot <- ggplot(data = signew_data_order, aes(y=group_fullname,x=percentage,fill=group,
                                                         # color=group,
                                                         alpha=1)) + 
    geom_col(width = 0.9) +
    geom_hline(yintercept = yinter,
              color="grey")+
    scale_fill_manual(values=as.vector(signew_data_order$sigcolor),
                      breaks = as.vector(signew_data_order$group)) +
    # scale_color_manual(values=as.vector(signew_data_order$sigcolor),
    #                    breaks = as.vector(signew_data_order$group)) +
    scale_x_continuous(breaks = c(0,0.5),
                       labels = c(0,0.5)) +
    facet_wrap2(~signame, 
                nrow = 1,
                scales = "free_x",
                strip = strip_vanilla(clip = "off")
                # labeller=labeller(variable = ss_labels)
    ) +
    ggtitle(figure_name)+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 3.91,color = "black",vjust = 0.5,hjust = 0),
          axis.title = element_blank(),
          axis.ticks.x =  element_blank(),
          panel.background = element_blank(),
          legend.position = "none",
          strip.text.x = element_text(size=6, angle = 0, colour = "black",hjust = 0,vjust = 0),
          strip.background = element_blank(),
          strip.placement = "outside",
          panel.spacing.x=unit(0, "lines"),
          plot.title = element_text(size=6,color = "black"),
          plot.margin=grid::unit(c(2,5,2,2), "mm")
    )
  
  signature_plot
  ggsave(file.path(output_signature_dir,paste0(figure_name,".pdf")), height = 73, width =(30+6*as.numeric(signaturenum))*1.4,units = c("mm"))
}


###################### RUN ###################
# PCAWG
plot_sigprofiler_function(
  phat_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/005_sigprofile/PCAWG_SV49_simple_clusterasone/CH49/All_Solutions/SBS49_13_Signatures/Signatures/SBS49_S13_Signatures.txt",
  output_signature_dir="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/005-sigprofiler-plot/PCAWG_SV49_simple_clusterasone",
  class_num="49",
  color_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/49sv_simple_matrix_color_svtype.csv",
  figure_name="PCAWG_SV49_simple_13_Signatures",
  signaturenum="13",
  signame_list=c("Del1","Del2","Del3","Del4","TD1","TD2","TD3","TD4","Fragile sites","Fb inv","Unbal inv","Large mixed","Unbal tra"),
  order=c("B","C","E","I","A","G","F","K","M","D","L","J","H")
)

plot_sigprofiler_function(
  phat_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/005_sigprofile/PCAWG_SV49_simple_clusterasone/CH49/All_Solutions/SBS49_12_Signatures/Signatures/SBS49_S12_Signatures.txt",
  output_signature_dir="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/005-sigprofiler-plot/PCAWG_SV49_simple_clusterasone",
  class_num="49",
  color_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/49sv_simple_matrix_color_svtype.csv",
  figure_name="PCAWG_SV49_simple_12_Signatures",
  signaturenum="12",
  signame_list=c("Del1","Del2","Del3","TD1","TD2","TD3","TD4","Fragile sites","Fb inv","Unbal inv","Large mixed","Unbal tra"),
  order=c("A","B","D","C","G","E","K","L","F","J","I","H")
)

plot_sigprofiler_function(
  phat_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/005_sigprofile/PCAWG_SV49_simple_clusterasone/CH49/All_Solutions/SBS49_14_Signatures/Signatures/SBS49_S14_Signatures.txt",
  output_signature_dir="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/005-sigprofiler-plot/PCAWG_SV49_simple_clusterasone",
  class_num="49",
  color_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/49sv_simple_matrix_color_svtype.csv",
  figure_name="PCAWG_SV49_simple_14_Signatures",
  signaturenum="14",
  signame_list=c("Del1","Del2","Del3","Del4","TD1","TD2","TD3","TD4","Fragile sites","Fb inv","Unbal inv","Large mixed","Unbal tra","Mid mixed"),
  order=c("A","C","F","J","B","H","E","K","N","D","L","M","I","G")
)
# Hartwig
plot_sigprofiler_function(
  phat_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/sigprofiler/Pub_Hartwig_SV49_simple_clusterasone/CH49/All_Solutions/SBS49_14_Signatures/Signatures/SBS49_S14_Signatures.txt",
  output_signature_dir="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/sigprofiler_plot",
  class_num="49",
  color_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/49sv_simple_matrix_color_svtype.csv",
  figure_name="Hartwig_SV49_simple_14_Signatures",
  signaturenum="14",
  signame_list=c("Del0","Del1","Del2","Del3","Del4","TD0","TD1","TD2","TD3","TD4","Fragile sites","Fb inv","Large mixed","Unbal tra"),
  order=c("A","B","E","C","L","N","G","H","F","J","M","D","I","K")
)
plot_sigprofiler_function(
  phat_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/sigprofiler/Pub_Hartwig_SV49_simple_clusterasone/CH49/All_Solutions/SBS49_13_Signatures/Signatures/SBS49_S13_Signatures.txt",
  output_signature_dir="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/sigprofiler_plot",
  class_num="49",
  color_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/49sv_simple_matrix_color_svtype.csv",
  figure_name="Hartwig_SV49_simple_13_Signatures",
  signaturenum="13",
  signame_list=c("Del0","Del1","Del2","Del3","Del4","TD0","TD1","TD2","TD3","TD4","Fragile sites","Fb inv","Tra"),
  order=c("A","C","E","B","J","M","G","H","F","I","L","D","K")
)
plot_sigprofiler_function(
  phat_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/sigprofiler/Pub_Hartwig_SV49_simple_clusterasone/CH49/All_Solutions/SBS49_15_Signatures/Signatures/SBS49_S15_Signatures.txt",
  output_signature_dir="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/sigprofiler_plot",
  class_num="49",
  color_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/49sv_simple_matrix_color_svtype.csv",
  figure_name="Hartwig_SV49_simple_15_Signatures",
  signaturenum="15",
  signame_list=c("Del0","Del1","Del2","Del3","Del4","Del5","TD0","TD1","TD2","TD3","TD4","Fragile sites","Fb inv","Large mixed","Unbal tra"),
  order=c("A","D","I","F","B","K", "N","G","H","E","L", "M","C","O","J")
)

# POG570
plot_sigprofiler_function(
  phat_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/sigprofiler/PCAWG_VALIDATION_PG570_SV49/CH49/All_Solutions/SBS49_13_Signatures/Signatures/SBS49_S13_Signatures.txt",
  output_signature_dir="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/sigprofiler_plot",
  class_num="49",
  color_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/49sv_simple_matrix_color_svtype.csv",
  figure_name="POG570_SV49_simple_13_Signatures",
  signaturenum="13",
  signame_list=c("Del0","Del1","Del2","Del3","TD0","TD1","TD2","TD3", "Fragile Sites","Fb inv","Unbal inv","Large mixed","Unbal tra"),
  order=c("A","E","D","G","M","C","J","I","K","F","H","L","B")
)
plot_sigprofiler_function(
  phat_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/sigprofiler/PCAWG_VALIDATION_PG570_SV49/CH49/All_Solutions/SBS49_14_Signatures/Signatures/SBS49_S14_Signatures.txt",
  output_signature_dir="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/sigprofiler_plot",
  class_num="49",
  color_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/49sv_simple_matrix_color_svtype.csv",
  figure_name="POG570_SV49_simple_14_Signatures",
  signaturenum="14",
  signame_list=c("Del0","Del1","Del2","Del3","Del4","TD0","TD1","TD2","TD3", "Fragile Sites","Fb inv","Unbal inv","Large mixed","Unbal tra"),
  order=c("A","E","D","F","L","N","C","G","J","K","H","I","M","B")
)
plot_sigprofiler_function(
  phat_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/sigprofiler/PCAWG_VALIDATION_PG570_SV49/CH49/All_Solutions/SBS49_15_Signatures/Signatures/SBS49_S15_Signatures.txt",
  output_signature_dir="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/sigprofiler_plot",
  class_num="49",
  color_path="D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/49sv_simple_matrix_color_svtype.csv",
  figure_name="POG570_SV49_simple_15_Signatures",
  signaturenum="15",
  signame_list=c("Del0","Del1","Del2","Del3","Del4","Del5","TD0","TD1","TD2","TD3", "Fragile Sites","Fb inv","Unbal inv","Large mixed","Unbal tra"),
  order=c("A","E","J","D","F","M","O","C","G","I","L","H","K","N","B")
)
