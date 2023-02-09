# library
library(jtools)
library(ggplot2)
library(ggstance)
library(ggpubr)
require(scales)
library(dplyr)

library(facetscales)

### functions
dot2na <- function(x){
  x=as.numeric(as.character(ifelse(x==".",NA,x)))
}

rep2bi  <- function(x){
  x=ifelse(x=="replicated",1,x)
  x=ifelse(x=="unzip",0,x)
  x=as.numeric(as.character(x))
}

na2zero <- function(x){
  x=as.numeric(as.character(ifelse(is.na(x),0,x)))
}

# preprocessing_pcawg
preprocess_pcawg_function <- function(features_path,feature.detail_path,logistic_path,real_sig_path,sample_path,input_col_need){
  if (!dir.exists(logistic_path)){
    dir.create(logistic_path,recursive = F,showWarnings = F)
  }
  
  ### read  files
  features <- read.delim(features_path)
  names(features)
  
  feature.detail <- read.csv(feature.detail_path,header = T,row.names = 1)
  real_sig <- read.delim(real_sig_path,header = T)
  
  sample <- read.csv(sample_path,fileEncoding = "UTF-8-BOM")
  
  ### define  features  and  clusters
  myfeatures_ctn <-  as.character(feature.detail[feature.detail$DataType=="continue","valuename"]) #continue value
  myfeatures_bi <-  as.character(feature.detail[feature.detail$DataType=="Binary","valuename"]) #binary
  myfeatures_rawvalue <- as.character(feature.detail[feature.detail$RawValue=="Yes","valuename"]) #raw value
  myfeature <- c(myfeatures_ctn,myfeatures_bi)
  myfeatures_ctn_full <-  as.character(feature.detail[feature.detail$DataType=="continue","fullname"])
  myfeatures_bi_full <-  as.character(feature.detail[feature.detail$DataType=="Binary","fullname"])
  names(myfeature) <- c(myfeatures_ctn_full,myfeatures_bi_full)
  
  
  #  delete breakpoints in >=2 genes
  duplicated_id <- features[duplicated(features[,"order"]),"order"]
  features_new <- features[!features[,"order"] %in% duplicated_id,]
  #  delete breakpoints not in any genes
  features_new <- features_new[!(features_new[,input_col_need[1]]=="."|features_new[,input_col_need[2]]=="."),]
  
  #  delete breakpoints in >=2 genes
  duplicated_id <- real_sig[duplicated(real_sig$order),"order"]
  real_sig_new <- real_sig[!real_sig$order %in% duplicated_id,]
  #  delete breakpoints not in any genes
  real_sig_new <- real_sig_new[!(real_sig_new[,input_col_need[1]] %in% "."|real_sig_new[,input_col_need[2]] %in% "."),]
  
  # remove SVs two ends in one gene
  total.id <- max(real_sig$order)
  brk1.id <- seq(1,total.id/2,1)
  brk2.id <- seq((total.id/2 + 1),total.id,1)
  
  brk1.row <- match(x = brk1.id, table = as.vector(real_sig_new[,"order"]), nomatch = NA)
  brk2.row <- match(x = brk2.id, table = as.vector(real_sig_new[,"order"]), nomatch = NA)
  
  brk1.gene <- real_sig_new[brk1.row,input_col_need[1]]
  brk2.gene <- real_sig_new[brk2.row,input_col_need[1]]
  
  sv.removed <- grep(TRUE,brk1.gene==brk2.gene)
  row.removed <- c(brk1.row[sv.removed],brk2.row[sv.removed])
  
  real_sig_new[row.removed,"endsinonegene"] <- 1
  real_sig_new$endsinonegene <- ifelse(grepl("td_|del_",real_sig_new$sig_group),real_sig_new$endsinonegene,NA)
  
  
  ### add signature
  if ("order" %in% names(features)){
    features_new <- merge(features_new,
                          unique(real_sig_new[is.na(real_sig_new$endsinonegene),c("sample","order","sig_group")]),
                          by.x = c("sample","order"),
                          by.y = c("sample","order"),
                          all.x = F,
                          all.y = T,
                          sort = F)
  }
  
  
  
  ### filter sample
  features_new <- features_new[features_new$sample %in% unique(sample$WGS),]
  
  #mysigs
  features_new$signature <- NA
  features_new$signature <-  ifelse(features_new$sig_group  %in% paste0("del_",seq(1,3,1)),"Del1",features_new$signature)
  features_new$signature <-  ifelse(features_new$sig_group  %in% paste0("del_",seq(4,6,1)),"Del2",features_new$signature)
  features_new$signature <-  ifelse(features_new$sig_group  %in% paste0("del_",seq(7,11,1)),"Del3",features_new$signature)
  features_new$signature <-  ifelse(features_new$sig_group  %in% paste0("del_",seq(12,15,1)),"Del4",features_new$signature)
  features_new$signature <-  ifelse(features_new$sig_group  %in% paste0("td_",seq(1,5,1)),"TD1",features_new$signature)
  features_new$signature <-  ifelse(features_new$sig_group  %in% paste0("td_",seq(6,9,1)),"TD2",features_new$signature)
  features_new$signature <-  ifelse(features_new$sig_group  %in% paste0("td_",seq(10,11,1)),"TD3",features_new$signature)
  features_new$signature <-  ifelse(features_new$sig_group  %in% paste0("td_",seq(12,16,1)),"TD4",features_new$signature)
  features_new$signature <-  ifelse(features_new$sig_group  %in% c("td_fragile","del_fragile"),"Fragile site",features_new$signature)
  features_new$signature <-  ifelse(features_new$sig_group  %in% c("fb_inv"),"Fb inv",features_new$signature)
  features_new$signature <-  ifelse(features_new$sig_group  %in% c("unbal_inv_1","unbal_inv_2","unbal_inv_3"),"Unbal inv",features_new$signature)
  features_new$signature <-  ifelse(features_new$sig_group  %in% c("unbal_tra"),"Unbal tra",features_new$signature)
  features_new$signature <-  ifelse(features_new$sig_group  %in% c("recip_tra",
                                                                   "recip_inv_1","recip_inv_2","recip_inv_3","recip_inv_4","recip_inv_5",
                                                                   paste0("del_",seq(16,18,1)),
                                                                   paste0("td_",seq(17,18,1))),"Large mixed",features_new$signature)
  
  ### feature type
  if (sum(names(features_new) %in% myfeatures_bi)>1){
    col <- names(features_new)[names(features_new) %in% myfeatures_bi]
    features_new[,col] <-  apply(features_new[,col] , 2, rep2bi)
    features_new[,col] <- lapply(features_new[,col], factor)
    features_new[,myfeatures_rawvalue] <-  apply(features_new[,myfeatures_rawvalue], 2, na2zero)
  } else {
    col <- names(features_new)[names(features_new) %in% myfeatures_bi]
    features_new[,col] <-  ifelse(features_new[,col]=="replicated",1, features_new[,col])
    features_new[,col] <-  ifelse(features_new[,col]=="unzip",0, features_new[,col])
    features_new[,col] <- as.factor(as.numeric(as.character(features_new[,col])))
    features_new[,myfeatures_rawvalue] <-  apply(features_new[,myfeatures_rawvalue] , 2, na2zero)
  }
  return(features_new)
}

# preprocessing_hartwig
preprocess_hartwig_function <- function(features_path,feature.detail_path,logistic_path,real_sig_path,sample_path,input_col_need){
  if (!dir.exists(logistic_path)){
    dir.create(logistic_path,recursive = F,showWarnings = F)
  }
  
  ### read  files
  features <- read.delim(features_path)
  names(features)
  
  feature.detail <- read.csv(feature.detail_path,header = T,row.names = 1)
  real_sig <- read.delim(real_sig_path,header = T)
  
  sample <- read.csv(sample_path,fileEncoding = "UTF-8-BOM")
  
  ### define  features  and  clusters
  myfeatures_ctn <-  as.character(feature.detail[feature.detail$DataType=="continue","valuename"]) #continue value
  myfeatures_bi <-  as.character(feature.detail[feature.detail$DataType=="Binary","valuename"]) #binary
  myfeatures_rawvalue <- as.character(feature.detail[feature.detail$RawValue=="Yes","valuename"]) #raw value
  myfeature <- c(myfeatures_ctn,myfeatures_bi)
  myfeatures_ctn_full <-  as.character(feature.detail[feature.detail$DataType=="continue","fullname"])
  myfeatures_bi_full <-  as.character(feature.detail[feature.detail$DataType=="Binary","fullname"])
  names(myfeature) <- c(myfeatures_ctn_full,myfeatures_bi_full)
  
  #  delete breakpoints in >=2 genes
  duplicated_id <- features[duplicated(features[,"order"]),"order"]
  features_new <- features[!features[,"order"] %in% duplicated_id,]
  #  delete breakpoints not in any genes
  features_new <- features_new[!(features_new[,input_col_need[1]]=="."|features_new[,input_col_need[2]]=="."),]
  
  #  delete breakpoints in >=2 genes
  duplicated_id <- real_sig[duplicated(real_sig$order),"order"]
  real_sig_new <- real_sig[!real_sig$order %in% duplicated_id,]
  #  delete breakpoints not in any genes
  real_sig_new <- real_sig_new[!(real_sig_new[,input_col_need[1]] %in% "."|real_sig_new[,input_col_need[2]] %in% "."),]
  
  # remove SVs two ends in one gene
  total.id <- max(real_sig$order)
  brk1.id <- seq(1,total.id/2,1)
  brk2.id <- seq((total.id/2 + 1),total.id,1)
  
  brk1.row <- match(x = brk1.id, table = as.vector(real_sig_new[,"order"]), nomatch = NA)
  brk2.row <- match(x = brk2.id, table = as.vector(real_sig_new[,"order"]), nomatch = NA)
  
  brk1.gene <- real_sig_new[brk1.row,input_col_need[1]]
  brk2.gene <- real_sig_new[brk2.row,input_col_need[1]]
  
  sv.removed <- grep(TRUE,brk1.gene==brk2.gene)
  row.removed <- c(brk1.row[sv.removed],brk2.row[sv.removed])
  
  real_sig_new[row.removed,"endsinonegene"] <- 1
  real_sig_new$endsinonegene <- ifelse(grepl("td_|del_",real_sig_new$sig_group),real_sig_new$endsinonegene,NA)
  
  
  ### add signature
  if ("order" %in% names(features)){
    features_new <- merge(features_new,
                          unique(real_sig_new[is.na(real_sig_new$endsinonegene),c("sample","order","sig_group")]),
                          by.x = c("sample","order"),
                          by.y = c("sample","order"),
                          all.x = F,
                          all.y = T,
                          sort = F)
  }
  
  
  
  ### filter sample
  features_new <- features_new[features_new$sample %in% unique(sample$WGS),]
  
  #mysigs
  features_new$signature <- NA
  features_new$signature <-  ifelse(features_new$sig_group  %in% paste0("del_",1),"Del0",features_new$signature)
  features_new$signature <-  ifelse(features_new$sig_group  %in% paste0("del_",seq(2,3,1)),"Del1",features_new$signature)
  features_new$signature <-  ifelse(features_new$sig_group  %in% paste0("del_",seq(4,6,1)),"Del2",features_new$signature)
  features_new$signature <-  ifelse(features_new$sig_group  %in% paste0("del_",seq(7,11,1)),"Del3",features_new$signature)
  features_new$signature <-  ifelse(features_new$sig_group  %in% paste0("del_",seq(12,15,1)),"Del4",features_new$signature)
  features_new$signature <-  ifelse(features_new$sig_group  %in% paste0("td_",1),"TD0",features_new$signature)
  features_new$signature <-  ifelse(features_new$sig_group  %in% paste0("td_",seq(2,5,1)),"TD1",features_new$signature)
  features_new$signature <-  ifelse(features_new$sig_group  %in% paste0("td_",seq(6,9,1)),"TD2",features_new$signature)
  features_new$signature <-  ifelse(features_new$sig_group  %in% paste0("td_",seq(10,12,1)),"TD3",features_new$signature)
  features_new$signature <-  ifelse(features_new$sig_group  %in% paste0("td_",seq(13,16,1)),"TD4",features_new$signature)
  features_new$signature <-  ifelse(features_new$sig_group  %in% c("td_fragile","del_fragile"),"Fragile site",features_new$signature)
  features_new$signature <-  ifelse(features_new$sig_group  %in% c("fb_inv","unbal_inv_1","unbal_inv_2"),"Fb inv",features_new$signature)
  features_new$signature <-  ifelse(features_new$sig_group  %in% c("unbal_tra"),"Unbal tra",features_new$signature)
  features_new$signature <-  ifelse(features_new$sig_group  %in% c("unbal_inv_3",
                                                                   "recip_tra",
                                                                   "recip_inv_1","recip_inv_2","recip_inv_3","recip_inv_4","recip_inv_5",
                                                                   paste0("del_",seq(16,18,1)),
                                                                   paste0("td_",seq(17,18,1))),"Large mixed",features_new$signature)
  
  
  
  ### feature type
  if (sum(names(features_new) %in% myfeatures_bi)>1){
    col <- names(features_new)[names(features_new) %in% myfeatures_bi]
    features_new[,col] <-  apply(features_new[,col] , 2, rep2bi)
    features_new[,col] <- lapply(features_new[,col], factor)
    features_new[,myfeatures_rawvalue] <-  apply(features_new[,myfeatures_rawvalue], 2, na2zero)
  } else {
    col <- names(features_new)[names(features_new) %in% myfeatures_bi]
    features_new[,col] <-  ifelse(features_new[,col]=="replicated",1, features_new[,col])
    features_new[,col] <-  ifelse(features_new[,col]=="unzip",0, features_new[,col])
    features_new[,col] <- as.factor(as.numeric(as.character(features_new[,col])))
    features_new[,myfeatures_rawvalue] <-  apply(features_new[,myfeatures_rawvalue] , 2, na2zero)
  }
  return(features_new)
}

# interaction_function
interaction_function(sv=features[features$signature=="Del1",],
                     expr.method=transi,
                     rt.method=rti,
                     name="Del1",
                     rtvalue="RTValue_BG02_ESC",
                     rpstrand="RepStrandconsistent8_shrink25000",
                     inttell=inttelli,
                     randomtell=randomtelli,
                     rttell=rttelli)
interaction_function <- function(sv,expr.method,rt.method,name,rtvalue,rpstrand,inttell,randomtell,rttell){
  print(name)
  if (randomtell==TRUE){
    sv <- sv
  } else {
    sv <- sv[sv$treatment=="real",]
  }
  #################################
  ########## expr.method ##########
  ################################
  ########### binary  ########### 
  if (expr.method=="50th"){
    sv$match.expression <- ifelse(is.na(sv[,"match_exp"]),NA,
                                  ifelse(sv[,"match_exp"]<50,0,1
                                  )
    )
    sv$match.expression <- factor(sv$match.expression)
  }
  ########### quartile ########### 
  if (expr.method=="25th"){
    sv$match.expression <- ifelse(is.na(sv[,"match_exp"]),NA,
                                  ifelse(sv[,"match_exp"]<25,1,
                                         ifelse(sv[,"match_exp"]<50,2,
                                                ifelse(sv[,"match_exp"]<75,3,
                                                       4)
                                         )
                                  )
    )
    sv$match.expression <- factor(sv$match.expression)
  }
  ########### 1vs4 ########### 
  if (expr.method=="1vs4"){
    sv$match.expression <- ifelse(is.na(sv[,"match_exp"]),NA,
                                  ifelse(sv[,"match_exp"]<25,1,
                                         ifelse(sv[,"match_exp"]<50,2,
                                                ifelse(sv[,"match_exp"]<75,3,
                                                       4)
                                         )
                                  )
    )
    sv <- sv[sv$match.expression %in% c(1,4),]
    sv$match.expression <- factor(sv$match.expression)
  }
  ########### 1vs234 ########### 
  if (expr.method=="1vs234"){
    sv$match.expression <- ifelse(is.na(sv[,"match_exp"]),NA,
                                  ifelse(sv[,"match_exp"]<25,0,1
                                  )
    )
    sv$match.expression <- factor(sv$match.expression)
  }
  ########### continue ###########
  if (expr.method=="continuous"){
    sv$match.expression <- as.numeric(floor(sv[,"match_exp"]))
  }
  
  #################################
  ########## rt.method ##########
  ################################
  ########### binary  ########### 
  median <- median(sv[,rtvalue],na.rm = T)
  rtquantile <- quantile(sv[,rtvalue],na.rm = T)
  if (rt.method=="50th"){
    sv$RTvalue <- ifelse(sv[,rtvalue]>=median,1,0)
    sv$RTvalue <- factor(sv$RTvalue,levels = c(0,1))
  }
  ########### quartile ########### 
  if (rt.method=="25th"){
    sv$RTvalue <- ifelse(sv[,rtvalue] < rtquantile[2],1,
                         ifelse(sv[,rtvalue] < rtquantile[3],2,
                                ifelse(sv[,rtvalue] < rtquantile[4],3,4
                                )
                         )
    )
    sv$RTvalue <- factor(sv$RTvalue,levels = c(1,2,3,4))
  }
  ########### 1vs4 ########### 
  if (rt.method=="1vs4"){
    sv$RTvalue <- ifelse(sv[,rtvalue] < rtquantile[2],1,
                         ifelse(sv[,rtvalue] < rtquantile[3],2,
                                ifelse(sv[,rtvalue] < rtquantile[4],3,4
                                )
                         )
    )
    sv <- sv[sv$RTvalue %in% c(1,4),]
    sv$RTvalue <- factor(sv$RTvalue,levels = c(1,4))
  }
  ########### 1vs234 ########### 
  if (rt.method=="1vs234"){
    sv$RTvalue <- ifelse(sv[,rtvalue] < rtquantile[2],1,
                         ifelse(sv[,rtvalue] < rtquantile[3],2,
                                ifelse(sv[,rtvalue] < rtquantile[4],3,4
                                )
                         )
    )
    sv$RTvalue <- as.factor(ifelse(is.na(sv$RTvalue),NA,
                                   ifelse(sv$RTvalue %in% c(2,3,4),1,0)))   
  }
  ########### continue ###########
  if (rt.method=="continuous"){
    sv$RTvalue <- as.numeric(floor(sv[,rtvalue]))
  }
  
  # glm
  if (randomtell==TRUE){
    fmla.cov <- paste(c(testfeatures[testfeatures %in%  c("match_exp") == FALSE],"treatment"), collapse= "+")
  } else {
    fmla.cov <- paste(c(testfeatures[testfeatures %in%  c("match_exp") == FALSE]), collapse= "+")
  }
  
  if (inttell==TRUE){
    fmla.int <- paste0("match.expression * " , "RTvalue")
  } else {
    fmla.int <- paste0("match.expression + " , "RTvalue")
  }
  
  if (rttell == TRUE){
    fmla.int <- fmla.int
  } else {
    fmla.int <- paste0("match.expression")
  }
  
  fmla <- as.formula(paste(rpstrand, "~",
                           ifelse(fmla.cov=="",
                                  fmla.int,
                                  paste0(c(fmla.cov,fmla.int),collapse="+")
                           )))
  glm.test <- tryCatch(glm(formula = fmla, family = binomial(link="logit"), data = sv),error = function(e) e)
  
  if (class(glm.test)[2]=="error"){
    # out_plot_names <- c("term","estimate","std.error","statistic","p.value", "pstar","Coeff","OR","2.5 %","97.5 %","fullname","signature")
    # out_plot <- data.frame(matrix(data=1,ncol = length(out_plot_names),nrow = length(testfeatures)))
    # names(out_plot) <- out_plot_names
    out_plot <- NA
  } else {
    myconfint <- tryCatch(confint(glm.test),error = function(e) e)
    myor <- tryCatch(coef(glm.test),error = function(e) e)
    out <- summary(glm.test)
    out_plot  <- as.data.frame(tidy(glm.test))
    # add pstar
    out_plot$pstar <- ""
    out_plot$pstar <- ifelse(out_plot$p.value< 0.05,"*",out_plot$pstar)
    out_plot$pstar <- ifelse(out_plot$p.value< 0.01,"**",out_plot$pstar)
    out_plot$pstar <- ifelse(out_plot$p.value< 0.001,"***",out_plot$pstar)
    out_plot$Coeff <- out_plot$term
    if (class(myconfint)[2] %in% "array" & class(myor) %in% "numeric"){
      # add OR 2.5% 97.5%
      out_or <- cbind(OR = myor, myconfint)
      out_or[is.na(out_or)] <- 0
      out_or <- round(exp(out_or), 4)
      out_plot <- cbind(out_plot,out_or)
    } else {
      out_plot$OR <- NA
      out_plot$`2.5 %` <- NA
      out_plot$`97.5 %` <- NA
    }
    out_plot <- out_plot[out_plot$Coeff!="(Intercept)",]
    # add fullname
    feature.name.df <- read.csv(feature.detail_path,encoding = "UTF-8") %>% select("valuename","fullname")
    rownames(feature.name.df) <- feature.name.df$valuename
    order2 <- feature.name.df[testfeatures,"fullname"]
    valuename.new <- c("treatmentreal","match.expression1","RTvalue1","match.expression1:RTvalue1")
    fullname.new <- c("SV (random-observed)","Expression (unexpressed-expressed)","Replication Timing (late-early)","Expression:Replication Timing")
    order1 <- fullname.new
    feature.name.df <- rbind(feature.name.df,
                               data.frame(
                                 "valuename"=valuename.new,
                                 "fullname"=fullname.new
                               ))
    out_plot <- merge(out_plot,feature.name.df[,c("valuename","fullname")],
                        by.x = "term",
                        by.y = "valuename",
                        all.x = T)
    # order full name
    fullname.order <- rev(unique(c(order1,order2)[c(order1,order2) %in% out_plot$fullname]))
    out_plot$fullname <- factor(out_plot$fullname,levels = fullname.order)
    out_plot.order <- unlist(lapply(fullname.order, function(x) {
        grep(gsub("\\(", ":", x),
             gsub("\\(", ":", out_plot$fullname))
      }))
    out_plot <- out_plot[rev(out_plot.order),]
    out_plot$signature <- name
    # add name
    out_plot$signature <- name
    
  }
  
  # return
  return(out_plot)
}


interaction_plot_function <- function(table,savepath){
  sig.order <- c(paste0("Del",seq(0,4,1)),
                 paste0("TD",seq(0,4,1)),
                 "Fragile site","Fb inv","Unbal inv",
                 "Large mixed","Unbal tra")
  table$signature <- factor(table$signature,levels = sig.order)
  table$OR <- ifelse(is.na(table$OR),1,table$OR)
  table$`2.5 %` <- ifelse(is.na(table$`2.5 %`),1,table$`2.5 %`)
  table$`97.5 %` <- ifelse(is.na(table$`97.5 %`),1,table$`97.5 %`)
  # change 0 to 0.001
  options(digits=1)
  scales_x <- list(
    "Del0"=scale_x_log10(),
    "Del1"=scale_x_log10(),
    "Del2"=scale_x_log10(),
    "Del3"=scale_x_log10(),
    "Del4"=scale_x_log10(),
    "TD0"=scale_x_log10(),
    "TD1"=scale_x_log10(),
    "TD2"=scale_x_log10(),
    "TD3"=scale_x_log10(),
    "TD4"=scale_x_log10(),
    "Fragile site"=scale_x_log10(),
    "Fb inv"=scale_x_log10(),
    "Unbal inv"=scale_x_log10(),
    "Large mixed"=scale_x_log10(),
    "Unbal tra"=scale_x_log10()
  )
  plot_out <- ggplot(data = table)+
    facet_grid_sc(.~signature,
               # scales = "free_x",
               scales=list(y="fixed",
                           x=scales_x
                           )
                 )+
    # scale_x_log10()+
    geom_pointrange(
      aes(y=fullname,xmin=`2.5 %`,xmax=`97.5 %`,x=OR),
      size=0.25,
      shape=18,
      lwd=0.5 /(.pt*72.27/96),
      fatten = 0.05,
      color="black")+
    geom_vline(aes(xintercept=1),linetype="dotted",color="grey30",size=0.25 /(.pt*72.27/96))+
    coord_cartesian(clip = "off")+
    geom_text(inherit.aes = FALSE,
              data = table,
              aes(x = 1,y = fullname,label = pstar),
              vjust = 0.5,
              hjust=ifelse(table$OR<1,1,0),
              color=ifelse(table$estimate < 0,"blue","red"),
              size=3)+
    theme(axis.text.x.bottom  = element_text(size=5,color="black"),
          axis.title.x.bottom = element_blank(),
          axis.text.x.top  = element_blank(),
          axis.title.x.top = element_text(size=6,color="black"),
          axis.text.y  = element_text(size=6,color="black",hjust = 0),
          axis.ticks.x.top = element_blank(),
          axis.ticks.x.bottom = element_line(size=0.25 /(.pt*72.27/96),color="black"),
          axis.ticks.length.x.bottom = unit(0.5,"mm"),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.x.bottom = element_line(size=0.25 /(.pt*72.27/96),color="black"),
          panel.background = element_blank(),
          strip.text = element_text(size=6,color="black"),
          strip.background = element_blank(),
          panel.spacing.y = unit(0,"mm"),
          plot.margin = unit(c(1, 0, 1, 0), "mm"),
          legend.position = "none")
  ggsave(savepath,width = 160,height = 80,units="mm",dpi = 300)
  savepath <- gsub(".pdf",".2.pdf",savepath)
  ggsave(savepath,width = 350,height = 80,units="mm",dpi = 300)
  return(plot_out)
}

### logistic for features separately 
testfeatures <- c(
  "GCContent",
  "CpGIsland_dist",
  "centromere_dist",
  "telomere_dist",
  
  "alu_dist",
  "l1_dist",
  "l2_dist",
  "ltr_dist",
  "mir_dist",
  "transposon_dist",
  "simple_repeat_dist",
  'low_complexity_dist',
  
  "short_tandem_repeats_dist",
  "Z_DNA_dist",
  "a_phased_repeats_dist",
  "direct_repeats_dist",
  "inverted_repeats_dist",
  "mirror_repeats_dist",
  "G_quadruplex_dist",
  
  "Gm12878H3k27acStdSig",
  "Gm12878H3k27me3StdSig",
  "Gm12878H3k36me3StdSig",
  "Gm12878H3k4me1StdSig",
  "Gm12878H3k4me2StdSig",
  "Gm12878H3k4me3StdSig",
  "Gm12878H3k79me2StdSig",
  "Gm12878H3k9acStdSig",
  "Gm12878H3k9me3StdSig",
  "Gm12878H4k20me1StdSig",
  
  "Nucleosome_Occupancy",
  "LAD_density",
  "TadGm12878_dist"
)

###############################################################################################
################################### Run PCAWG #################################################
###############################################################################################
features_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/real_random_brkpt_withorder_withexpression_withfeatures.txt"
feature.detail_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/feature_detail.csv"
logistic_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/logistic"
real_sig_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-random/real_random_brkpt_withorder_withexpression_withsig.txt"
sample_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/pcawg2583_sampleid_match.csv"

features <- preprocess_pcawg_function(features_path,feature.detail_path,logistic_path,real_sig_path,sample_path,input_col_need = c("V12","V13"))
### PCAWG PANCANCER
count <- 0
for (transi in c("continuous","50th","25th","1vs4","1vs234")[c(5)]){
  for (rti in c("continuous","50th","25th","1vs4","1vs234")[2]){
    for (randomtelli in c(FALSE,TRUE)){
      for (inttelli in c(FALSE,TRUE)[2]){
        for (rttelli in c(FALSE,TRUE)[2]){
          if (inttelli==TRUE & rttelli==FALSE){
          }else {
            count <- count + 1
            intdel1 <- interaction_function(features[features$signature=="Del1",],transi,rti,"Del1","RTValue_BG02_ESC","RepStrandconsistent8_shrink25000",inttelli,randomtelli,rttelli)
            intdel2 <- interaction_function(features[features$signature=="Del2",],transi,rti,"Del2","RTValue_BG02_ESC","RepStrandconsistent8_shrink25000",inttelli,randomtelli,rttelli)
            intdel3 <- interaction_function(features[features$signature=="Del3",],transi,rti,"Del3","RTValue_BG02_ESC","RepStrandconsistent8_shrink25000",inttelli,randomtelli,rttelli)
            intdel4 <- interaction_function(features[features$signature=="Del4",],transi,rti,"Del4","RTValue_BG02_ESC","RepStrandconsistent8_shrink25000",inttelli,randomtelli,rttelli)
            inttd1 <- interaction_function(features[features$signature=="TD1",],transi,rti,"TD1","RTValue_BG02_ESC","RepStrandconsistent8_shrink25000",inttelli,randomtelli,rttelli)
            inttd2 <- interaction_function(features[features$signature=="TD2",],transi,rti,"TD2","RTValue_BG02_ESC","RepStrandconsistent8_shrink25000",inttelli,randomtelli,rttelli)
            inttd3 <- interaction_function(features[features$signature=="TD3",],transi,rti,"TD3","RTValue_BG02_ESC","RepStrandconsistent8_shrink25000",inttelli,randomtelli,rttelli)
            inttd4 <- interaction_function(features[features$signature=="TD4",],transi,rti,"TD4","RTValue_BG02_ESC","RepStrandconsistent8_shrink25000",inttelli,randomtelli,rttelli)
            intfragile <- interaction_function(features[features$signature=="Fragile site",],transi,rti,"Fragile site","RTValue_BG02_ESC","RepStrandconsistent8_shrink25000",inttelli,randomtelli,rttelli)
            intfbinv <- interaction_function(features[features$signature=="Fb inv",],transi,rti,"Fb inv","RTValue_BG02_ESC","RepStrandconsistent8_shrink25000",inttelli,randomtelli,rttelli)
            intunbalinv <- interaction_function(features[features$signature=="Unbal inv",],transi,rti,"Unbal inv","RTValue_BG02_ESC","RepStrandconsistent8_shrink25000",inttelli,randomtelli,rttelli)
            intunbaltra <- interaction_function(features[features$signature=="Unbal tra",],transi,rti,"Unbal tra","RTValue_BG02_ESC","RepStrandconsistent8_shrink25000",inttelli,randomtelli,rttelli)
            intlargemixed <- interaction_function(features[features$signature=="Large mixed",],transi,rti,"Large mixed","RTValue_BG02_ESC","RepStrandconsistent8_shrink25000",inttelli,randomtelli,rttelli)
            # save
            print("save")
            filename <- paste0("expr-",transi,".rt-",rti,".int",sum(inttelli),".random",sum(randomtelli),".rt",sum(rttelli))
            filename <- paste0("pcawg.pancancer.RTValue_BG02_ESC.RepStrandconsistent8.",filename,".csv")
            table <- rbind(intdel1,intdel2,
                           intdel3,intdel4,
                           inttd1,inttd2,inttd3,inttd4,
                           intfragile,
                           intfbinv,intunbalinv,intunbaltra,intlargemixed
                           )
            table <- table[is.na(table$term)==FALSE,]
            write.csv(table,file.path(logistic_path,filename),quote = F,row.names = F)
            # plot
            filename <- gsub(".csv",".pdf",filename)
            interaction_plot_function(table,file.path(logistic_path,filename))
          }
        }
      }
    }
  }
}

###############################################################################################
################################### Run Hartwig ###############################################
###############################################################################################
features_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/real_random_brkpt_withorder_withexpression_withfeatures.txt"
feature.detail_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/feature_detail.csv"
logistic_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/logistic"
real_sig_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/real_random_brkpt_withorder_withexpression_withsig.txt"
sample_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/hartwig2367_sampleid_match.csv"

features <- preprocess_hartwig_function(features_path,feature.detail_path,logistic_path,real_sig_path,sample_path,input_col_need = c("V20","V21"))
### Hartwig PANCANCER
count <- 0
for (transi in c("continuous","50th","25th","1vs4","1vs234")[c(5)]){
  for (rti in c("continuous","50th","25th","1vs4","1vs234")[2]){
    for (randomtelli in c(FALSE,TRUE)){
      for (inttelli in c(FALSE,TRUE)[2]){
        for (rttelli in c(FALSE,TRUE)[2]){
          if (inttelli==TRUE & rttelli==FALSE){
          }else {
            count <- count + 1
            intdel0 <- interaction_function(features[features$signature=="Del0",],transi,rti,"Del0","RTValue_BG02_ESC","RepStrandconsistent8_shrink25000",inttelli,randomtelli,rttelli)
            intdel1 <- interaction_function(features[features$signature=="Del1",],transi,rti,"Del1","RTValue_BG02_ESC","RepStrandconsistent8_shrink25000",inttelli,randomtelli,rttelli)
            intdel2 <- interaction_function(features[features$signature=="Del2",],transi,rti,"Del2","RTValue_BG02_ESC","RepStrandconsistent8_shrink25000",inttelli,randomtelli,rttelli)
            intdel3 <- interaction_function(features[features$signature=="Del3",],transi,rti,"Del3","RTValue_BG02_ESC","RepStrandconsistent8_shrink25000",inttelli,randomtelli,rttelli)
            intdel4 <- interaction_function(features[features$signature=="Del4",],transi,rti,"Del4","RTValue_BG02_ESC","RepStrandconsistent8_shrink25000",inttelli,randomtelli,rttelli)
            inttd0 <- interaction_function(features[features$signature=="TD0",],transi,rti,"TD0","RTValue_BG02_ESC","RepStrandconsistent8_shrink25000",inttelli,randomtelli,rttelli)
            inttd1 <- interaction_function(features[features$signature=="TD1",],transi,rti,"TD1","RTValue_BG02_ESC","RepStrandconsistent8_shrink25000",inttelli,randomtelli,rttelli)
            inttd2 <- interaction_function(features[features$signature=="TD2",],transi,rti,"TD2","RTValue_BG02_ESC","RepStrandconsistent8_shrink25000",inttelli,randomtelli,rttelli)
            inttd3 <- interaction_function(features[features$signature=="TD3",],transi,rti,"TD3","RTValue_BG02_ESC","RepStrandconsistent8_shrink25000",inttelli,randomtelli,rttelli)
            inttd4 <- interaction_function(features[features$signature=="TD4",],transi,rti,"TD4","RTValue_BG02_ESC","RepStrandconsistent8_shrink25000",inttelli,randomtelli,rttelli)
            intfragile <- interaction_function(features[features$signature=="Fragile site",],transi,rti,"Fragile site","RTValue_BG02_ESC","RepStrandconsistent8_shrink25000",inttelli,randomtelli,rttelli)
            intfbinv <- interaction_function(features[features$signature=="Fb inv",],transi,rti,"Fb inv","RTValue_BG02_ESC","RepStrandconsistent8_shrink25000",inttelli,randomtelli,rttelli)
            intunbaltra <- interaction_function(features[features$signature=="Unbal tra",],transi,rti,"Unbal tra","RTValue_BG02_ESC","RepStrandconsistent8_shrink25000",inttelli,randomtelli,rttelli)
            intlargemixed <- interaction_function(features[features$signature=="Large mixed",],transi,rti,"Large mixed","RTValue_BG02_ESC","RepStrandconsistent8_shrink25000",inttelli,randomtelli,rttelli)
            # save
            print("save")
            filename <- paste0("expr-",transi,".rt-",rti,".int",sum(inttelli),".random",sum(randomtelli),".rt",sum(rttelli))
            filename <- paste0("Hartwig.pancancer.RTValue_BG02_ESC.RepStrandconsistent8.interaction.",filename,".csv")
            table <- rbind(intdel0,
                           intdel1,intdel2,
                           intdel3,intdel4,
                           inttd0,
                           inttd1,inttd2,
                           inttd3,inttd4,
                           intfragile,
                           intfbinv,
                           intunbaltra,intlargemixed
            )
            table <- table[is.na(table$term)==FALSE,]
            write.csv(table,file.path(logistic_path,filename),quote = F,row.names = F)
            # plot
            filename <- gsub(".csv",".pdf",filename)
            interaction_plot_function(table,file.path(logistic_path,filename))
          }
        }
      }
    }
  }
}



###############################################################################################
################################### Run all BRCA ###############################################
###############################################################################################
# hartwig brca
features_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/real_random_brkpt_withorder_withexpression_withfeatures.txt"
feature.detail_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/feature_detail.csv"
logistic_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/logistic"
real_sig_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/real_random_brkpt_withorder_withexpression_withsig.txt"
sample_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/hartwig2367_sampleid_match.csv"
features <- preprocess_hartwig_function(features_path,feature.detail_path,logistic_path,real_sig_path,sample_path,input_col_need = c("V20","V21"))
features_hartwig <- features[features$histology_abbreviation=="Breast" ,c("sample",testfeatures,"treatment","match_exp","normal_tcga_exp" ,"RTValue_MCF7","RepStrandMCF7_shrink25000","signature") ]
# pcawg brca
features_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/real_random_brkpt_withorder_withexpression_withfeatures.txt"
feature.detail_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/feature_detail.csv"
logistic_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/logistic"
real_sig_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-random/real_random_brkpt_withorder_withexpression_withsig.txt"
sample_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/pcawg2583_sampleid_match.csv"
features <- preprocess_pcawg_function(features_path,feature.detail_path,logistic_path,real_sig_path,sample_path,input_col_need = c("V12","V13"))
features_pcawg <- features[features$histology_abbreviation=="Breast-AdenoCA" ,c("sample",testfeatures,"treatment","match_exp","normal_tcga_exp" ,"RTValue_MCF7","RepStrandMCF7_shrink25000","signature") ]
# brca-eu
features_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/BRCA-EU/real_random_brkpt_withorder_withexpression_withfeatures.txt"
feature.detail_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/feature_detail.csv"
logistic_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/BRCA-EU/genomic_features/logistic"
real_sig_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/BRCA-EU/real_random_brkpt_withorder_withexpression_withsig.txt"
sample_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/BRCA-EU/brcaeu532_sampleid_match.csv"
features <- preprocess_pcawg_function(features_path,feature.detail_path,logistic_path,real_sig_path,sample_path,input_col_need = c("V20","V21"))
features_brcaeu <- features[,c("sample",testfeatures,"treatment","match_exp","normal_tcga_exp" ,"RTValue_MCF7","RepStrandMCF7_shrink25000","signature") ]
# pog570 brca
features_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/real_random_brkpt_withorder_withexpression_withfeatures.txt"
feature.detail_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/feature_detail.csv"
logistic_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/logistic"
real_sig_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/real_random_brkpt_withorder_withexpression_withsig.txt"
sample_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/pog570_sampleid_match.csv"
features <- preprocess_hartwig_function(features_path,feature.detail_path,logistic_path,real_sig_path,sample_path,input_col_need = c("V20","V21"))
features_pog570 <- features[features$histology_abbreviation=="BRCA",c("sample",testfeatures,"treatment","match_exp","normal_tcga_exp","RTValue_MCF7","RepStrandMCF7_shrink25000","signature") ]

features_brca <- rbind(
  features_hartwig,
  features_pcawg,
  features_brcaeu,
  features_pog570)

### Merged BRCA
count <- 0
for (transi in c("continuous","50th","25th","1vs4","1vs234")[c(5)]){
  for (rti in c("continuous","50th","25th","1vs4","1vs234")[2]){
    for (randomtelli in c(FALSE,TRUE)){
      for (inttelli in c(FALSE,TRUE)[2]){
        for (rttelli in c(FALSE,TRUE)[2]){
          if (inttelli==TRUE & rttelli==FALSE){
          }else {
            count <- count + 1
            intdel0 <- interaction_function(features_brca[features_brca$signature=="Del0",],transi,rti,"Del0","RTValue_MCF7","RepStrandMCF7_shrink25000",inttelli,randomtelli,rttelli)
            intdel1 <- interaction_function(features_brca[features_brca$signature=="Del1",],transi,rti,"Del1","RTValue_MCF7","RepStrandMCF7_shrink25000",inttelli,randomtelli,rttelli)
            intdel2 <- interaction_function(features_brca[features_brca$signature=="Del2",],transi,rti,"Del2","RTValue_MCF7","RepStrandMCF7_shrink25000",inttelli,randomtelli,rttelli)
            intdel3 <- interaction_function(features_brca[features_brca$signature=="Del3",],transi,rti,"Del3","RTValue_MCF7","RepStrandMCF7_shrink25000",inttelli,randomtelli,rttelli)
            intdel4 <- interaction_function(features_brca[features_brca$signature=="Del4",],transi,rti,"Del4","RTValue_MCF7","RepStrandMCF7_shrink25000",inttelli,randomtelli,rttelli)
            inttd0 <- interaction_function(features_brca[features_brca$signature=="TD0",],transi,rti,"TD0","RTValue_MCF7","RepStrandMCF7_shrink25000",inttelli,randomtelli,rttelli)
            inttd1 <- interaction_function(features_brca[features_brca$signature=="TD1",],transi,rti,"TD1","RTValue_MCF7","RepStrandMCF7_shrink25000",inttelli,randomtelli,rttelli)
            inttd2 <- interaction_function(features_brca[features_brca$signature=="TD2",],transi,rti,"TD2","RTValue_MCF7","RepStrandMCF7_shrink25000",inttelli,randomtelli,rttelli)
            inttd3 <- interaction_function(features_brca[features_brca$signature=="TD3",],transi,rti,"TD3","RTValue_MCF7","RepStrandMCF7_shrink25000",inttelli,randomtelli,rttelli)
            inttd4 <- interaction_function(features_brca[features_brca$signature=="TD4",],transi,rti,"TD4","RTValue_MCF7","RepStrandMCF7_shrink25000",inttelli,randomtelli,rttelli)
            intfragile <- interaction_function(features_brca[features_brca$signature=="Fragile site",],transi,rti,"Fragile site","RTValue_MCF7","RepStrandMCF7_shrink25000",inttelli,randomtelli,rttelli)
            intfbinv <- interaction_function(features_brca[features_brca$signature=="Fb inv",],transi,rti,"Fb inv","RTValue_MCF7","RepStrandMCF7_shrink25000",inttelli,randomtelli,rttelli)
            intunbalinv <- interaction_function(features_brca[features_brca$signature=="Unbal inv",],transi,rti,"Unbal inv","RTValue_MCF7","RepStrandMCF7_shrink25000",inttelli,randomtelli,rttelli)
            intunbaltra <- interaction_function(features_brca[features_brca$signature=="Unbal tra",],transi,rti,"Unbal tra","RTValue_MCF7","RepStrandMCF7_shrink25000",inttelli,randomtelli,rttelli)
            intlargemixed <- interaction_function(features_brca[features_brca$signature=="Large mixed",],transi,rti,"Large mixed","RTValue_MCF7","RepStrandMCF7_shrink25000",inttelli,randomtelli,rttelli)
            # save
            print("save")
            filename <- paste0("expr-",transi,".rt-",rti,".int",sum(inttelli),".random",sum(randomtelli),".rt",sum(rttelli))
            filename <- paste0("AllBRCA.RTValue_MCF7.RepStrandMCF7.interaction.",filename,".csv")
            table <- rbind(intdel0,
                           intdel1,intdel2,
                           intdel3,intdel4,
                           inttd0,
                           inttd1,inttd2,
                           inttd3,inttd4,
                           intfragile,
                           intfbinv,intunbalinv,
                           intunbaltra,intlargemixed
            )
            table <- table[is.na(table$term)==FALSE,]
            write.csv(table,file.path(logistic_path,filename),quote = F,row.names = F)
            # plot
            filename <- gsub(".csv",".pdf",filename)
            interaction_plot_function(table,file.path("D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/logistic",
                                                      filename))
          }
        }
      }
    }
  }
}
