# function
test_bias_function <- function(raw,
                               input_path,
                               tumortype,
                               input_col_need,
                               outpath,
                               assignment,
                               keep){
  if (!is.na(tumortype)){
    raw2 <- raw[raw$histology_abbreviation==tumortype,]
  } else {
    raw2 <- raw
  }
  
  raw2 <- raw2[raw2$sample %in% unique(keep$WGS),]
  
  
  input <- read.delim(input_path,header = T)
  
  ################################################# add features #################################################
  #  delete breakpoints in >=2 genes
  duplicated_id <- input[duplicated(input[,input_col_need[1]]),input_col_need[1]]
  input_new <- input[!input[,input_col_need[1]] %in% duplicated_id,]
  #  delete breakpoints not in any genes
  input_new <- input_new[!(input_new[,input_col_need[2]]=="."|input_new[,input_col_need[3]]=="."),]
  
  #  delete breakpoints in >=2 genes
  duplicated_id <- raw2[duplicated(raw2$order),"order"]
  raw_new <- raw2[!raw2$order %in% duplicated_id,]
  #  delete breakpoints not in any genes
  raw_new <- raw_new[!(raw_new[,input_col_need[2]] %in% "."|raw_new[,input_col_need[3]] %in% "."),]
  
  # uniform  colnames
  names(raw_new)[names(raw_new)=="chrom1"] <- "chrom"
  names(raw_new)[names(raw_new)=="start1"] <- "start"
  names(raw_new)[names(raw_new)=="strand1"] <- "strand"
  

  # remove SVs two ends in one gene
  total.id <- max(raw$order)
  brk1.id <- seq(1,total.id/2,1)
  brk2.id <- seq((total.id/2 + 1),total.id,1)
  
  brk1.row <- match(x = brk1.id, table = as.vector(raw_new[,"order"]), nomatch = NA)
  brk2.row <- match(x = brk2.id, table = as.vector(raw_new[,"order"]), nomatch = NA)
  
  brk1.gene <- raw_new[brk1.row,input_col_need[2]]
  brk2.gene <- raw_new[brk2.row,input_col_need[2]]
  
  
  sv.removed <- grep(TRUE,brk1.gene==brk2.gene) #only for deletion duplication
  row.removed <- c(brk1.row[sv.removed],brk2.row[sv.removed])

  raw_new[row.removed,"endsinonegene"] <- 1
  raw_new$endsinonegene <- ifelse(grepl("td_|del_",raw_new$sig_group),raw_new$endsinonegene,NA)
  
  raw_feature <- merge(raw_new[is.na(raw_new$endsinonegene),],
                       unique(input_new[,c("V1","V2","V4",input_col_need[1],tail(names(input_new),4))]),
                       by.x = c("chrom","start","strand","order"),
                       by.y = c("V1","V2","V4",input_col_need[1]),
                       all.x = T,
                       sort = F)
  
  
  ################################################# test #################################################
  t.real <- raw_feature$treatment %in% "real"
  t.random <- raw_feature$treatment %in% "random"
  t.replicated <- raw_feature$replicated %in% "replicated"
  t.unrep <- raw_feature$replicated %in% "unreplicated"
  t.oppo  <-  raw_feature$direction %in% "oppo"
  t.same  <-  raw_feature$direction %in% "same"
  t.5 <- raw_feature$prime %in% "5"
  t.3 <- raw_feature$prime %in% "3"
  
  t.match1  <- raw_feature$match_exp %in% seq(0,24,1)
  t.match2  <- raw_feature$match_exp %in% seq(25,49,1)
  t.match3  <- raw_feature$match_exp %in% seq(50,74,1)
  t.match4  <- raw_feature$match_exp %in% seq(75,100,1)
  t.match234 <- raw_feature$match_exp %in% seq(25,100,1)
  t.match12 <- raw_feature$match_exp %in% seq(0,49,1)
  t.match34 <- raw_feature$match_exp %in% seq(50,100,1)
  
  t.normalpcawg1  <- raw_feature$normal_pcawg_exp %in% seq(0,24,1)
  t.normalpcawg2  <- raw_feature$normal_pcawg_exp %in% seq(25,49,1)
  t.normalpcawg3  <- raw_feature$normal_pcawg_exp %in% seq(50,74,1)
  t.normalpcawg4  <- raw_feature$normal_pcawg_exp %in% seq(75,100,1)
  t.normalpcawg234  <- raw_feature$normal_pcawg_exp %in% seq(25,100,1)
  t.normalpcawg12  <- raw_feature$normal_pcawg_exp %in% seq(0,49,1)
  t.normalpcawg34  <- raw_feature$normal_pcawg_exp %in% seq(50,100,1)
  
  t.tumorpcawg1  <- raw_feature$tumor_pcawg_exp %in% seq(0,24,1)
  t.tumorpcawg2  <- raw_feature$tumor_pcawg_exp %in% seq(25,49,1)
  t.tumorpcawg3  <- raw_feature$tumor_pcawg_exp %in% seq(50,74,1)
  t.tumorpcawg4  <- raw_feature$tumor_pcawg_exp %in% seq(75,100,1)
  t.tumorpcawg234  <- raw_feature$tumor_pcawg_exp %in% seq(25,100,1)
  t.tumorpcawg12  <- raw_feature$tumor_pcawg_exp %in% seq(0,49,1)
  t.tumorpcawg34  <- raw_feature$tumor_pcawg_exp %in% seq(50,100,1)
  
  t.normaltcga1  <- raw_feature$normal_tcga_exp %in% seq(0,24,1)
  t.normaltcga2  <- raw_feature$normal_tcga_exp %in% seq(25,49,1)
  t.normaltcga3  <- raw_feature$normal_tcga_exp %in% seq(50,74,1)
  t.normaltcga4  <- raw_feature$normal_tcga_exp %in% seq(75,100,1)
  t.normaltcga234  <- raw_feature$normal_tcga_exp %in% seq(25,100,1)
  t.normaltcga12  <- raw_feature$normal_tcga_exp %in% seq(0,49,1)
  t.normaltcga34  <- raw_feature$normal_tcga_exp %in% seq(50,100,1)
  
  t.tumortcga1  <- raw_feature$tumor_tcga_exp %in% seq(0,24,1)
  t.tumortcga2  <- raw_feature$tumor_tcga_exp %in% seq(25,49,1)
  t.tumortcga3  <- raw_feature$tumor_tcga_exp %in% seq(50,74,1)
  t.tumortcga4  <- raw_feature$tumor_tcga_exp %in% seq(75,100,1)
  t.tumortcga234  <- raw_feature$tumor_tcga_exp %in% seq(25,100,1)
  t.tumortcga12  <- raw_feature$tumor_tcga_exp %in% seq(0,49,1)
  t.tumortcga34  <- raw_feature$tumor_tcga_exp %in% seq(50,100,1)
  
  expression_list <- list(t.match1,t.match2,t.match3,t.match4,t.match12,t.match34,t.match234,
                          t.normalpcawg1,t.normalpcawg2,t.normalpcawg3,t.normalpcawg4,t.normalpcawg12,t.normalpcawg34,t.normalpcawg234,
                          t.tumorpcawg1,t.tumorpcawg2,t.tumorpcawg3,t.tumorpcawg4,t.tumorpcawg12,t.tumorpcawg34,t.tumorpcawg234,
                          t.normaltcga1,t.normaltcga2,t.normaltcga3,t.normaltcga4,t.normaltcga12,t.normaltcga34,t.normaltcga234,
                          t.tumortcga1,t.tumortcga2,t.tumortcga3,t.tumortcga4,t.tumortcga12,t.tumortcga34,t.tumortcga234,
                          TRUE)
  names(expression_list) <- c("match1","match2","match3","match4","match12","match34","match234",
                              "normalpcawg1","normalpcawg2","normalpcawg3","normalpcawg4","normalpcawg12","normalpcawg34","normalpcawg234",
                              "tumorpcawg1","tumorpcawg2","tumorpcawg3","tumorpcawg4","tumorpcawg12","tumorpcawg34","tumorpcawg234",
                              "normaltcga1","normaltcga2","normaltcga3","normaltcga4","normaltcga12","normaltcga34","normaltcga234",
                              "tumortcga1","tumortcga2","tumortcga3","tumortcga4","tumortcga12","tumortcga34","tumortcga234",
                              "all")
  
  
  if  (assignment=="pcawg"){
    t.del1 <-  raw_feature$sig_group  %in% paste0("del_",seq(1,3,1))
    t.del2 <-  raw_feature$sig_group  %in% paste0("del_",seq(4,6,1))
    t.del3 <-  raw_feature$sig_group  %in% paste0("del_",seq(7,11,1))
    t.del4 <-  raw_feature$sig_group  %in% paste0("del_",seq(12,15,1))

    t.dup1  <-  raw_feature$sig_group  %in% paste0("td_",seq(1,5,1))
    t.dup2  <-  raw_feature$sig_group  %in% paste0("td_",seq(6,9,1))
    t.dup3  <-  raw_feature$sig_group  %in% paste0("td_",seq(10,11,1))
    t.dup4 <- raw_feature$sig_group  %in% paste0("td_",seq(12,16,1))

    t.fragile <- raw_feature$sig_group  %in% c("td_fragile","del_fragile")
    t.fb <- raw_feature$sig_group  %in% c("fb_inv")
    t.unbalinv <- raw_feature$sig_group  %in% c("unbal_inv_1","unbal_inv_2","unbal_inv_3")
    t.unbaltra <- raw_feature$sig_group  %in% c("unbal_tra")
    t.largemixed <- raw_feature$sig_group  %in% c("recip_tra",
                                             "recip_inv_1","recip_inv_2","recip_inv_3","recip_inv_4","recip_inv_5",
                                             paste0("del_",seq(16,18,1)),
                                             paste0("td_",seq(17,18,1))
    )
    
    mysig_list <- list(t.del1,t.del2,t.del3,t.del4, 
                       t.dup1,t.dup2,t.dup3,t.dup4,
                       t.fragile,
                       t.fb,
                       t.unbalinv,
                       t.unbaltra,
                       t.largemixed
    )
    names(mysig_list) <- c("del1","del2","del3","del4",
                           "dup1","dup2","dup3","dup4",
                           "fragile",
                           "fbinv",
                           "unbalinv",
                           "unbaltra",
                           "largemixed")
  }
  
  if  (assignment=="hartwig"){
    t.del0 <-  raw_feature$sig_group  %in% paste0("del_",1)
    t.del1 <-  raw_feature$sig_group  %in% paste0("del_",seq(2,3,1))
    t.del2 <-  raw_feature$sig_group  %in% paste0("del_",seq(4,6,1))
    t.del3 <-  raw_feature$sig_group  %in% paste0("del_",seq(7,11,1))
    t.del4 <-  raw_feature$sig_group  %in% paste0("del_",seq(12,15,1))
    
    t.dup0  <-  raw_feature$sig_group  %in% paste0("td_",1)
    t.dup1  <-  raw_feature$sig_group  %in% paste0("td_",seq(2,5,1))
    t.dup2  <-  raw_feature$sig_group  %in% paste0("td_",seq(6,9,1))
    t.dup3  <-  raw_feature$sig_group  %in% paste0("td_",seq(10,12,1))
    t.dup4 <- raw_feature$sig_group  %in% paste0("td_",seq(13,16,1))
    
    t.fragile <- raw_feature$sig_group  %in% c("td_fragile","del_fragile")
    t.fb <- raw_feature$sig_group  %in% c("fb_inv","unbal_inv_1","unbal_inv_2")
    t.unbaltra <- raw_feature$sig_group  %in% c("unbal_tra")
    t.largemixed <- raw_feature$sig_group  %in% c("unbal_inv_3",
                                                  "recip_tra",
                                                  "recip_inv_1","recip_inv_2","recip_inv_3","recip_inv_4","recip_inv_5",
                                                  paste0("del_",seq(16,18,1)),
                                                  paste0("td_",seq(17,18,1))
    )
    
    mysig_list <- list(t.del0,t.del1,t.del2,t.del3,t.del4, 
                       t.dup0,t.dup1,t.dup2,t.dup3,t.dup4,
                       t.fragile,
                       t.fb,
                       t.unbaltra,
                       t.largemixed
    )
    names(mysig_list) <- c("del0","del1","del2","del3","del4",
                           "dup0","dup1","dup2","dup3","dup4",
                           "fragile",
                           "fbinv",
                           "unbaltra",
                           "largemixed")
  }
  
  # test
  test_df <- data.frame()
  for (j in 1:length(expression_list)){
    exp_name <- names(expression_list)[j]
    exp <- expression_list[[j]]
    for (i in 1:length(mysig_list)){
      signame <- names(mysig_list)[i]
      t.sig <- mysig_list[[i]]
      
      #relicated vs. unrep
      a <- sum((t.real&t.replicated)&t.sig&exp)
      b <- sum((t.real&t.unrep)&t.sig&exp)
      c <- sum((t.random&t.replicated)&t.sig&exp)
      d <- sum((t.random&t.unrep)&t.sig&exp)
      test.matrix <- matrix(data = c(a,b,c,d),nrow = 2)
      test <- fisher.test(test.matrix)
      test_df <- rbind(test_df,
                       data.frame("signature"=signame,
                                  "test"="rep.vs.unrep",
                                  "expression"=exp_name,
                                  "a"=a,
                                  "b"=b,
                                  "c"=c,
                                  "d"=d,
                                  "raw.p"=test$p.value,
                                  "or"=test$estimate,
                                  "adjust.p"=p.adjust(test$p.value,method = "bonferroni",n=length(mysig_list))
                       )
      )
      
      #relicated vs. unrep at oppo
      a <- sum((t.real&t.replicated&t.oppo)&t.sig&exp)
      b <- sum((t.real&t.unrep&t.oppo)&t.sig&exp)
      c <- sum((t.random&t.replicated&t.oppo)&t.sig&exp)
      d <- sum((t.random&t.unrep&t.oppo)&t.sig&exp)
      test.matrix <- matrix(data = c(a,b,c,d),nrow = 2)
      test <- fisher.test(test.matrix)
      test_df <- rbind(test_df,
                       data.frame("signature"=signame,
                                  "test"="rep.vs.unrep.oppo",
                                  "expression"=exp_name,
                                  "a"=a,
                                  "b"=b,
                                  "c"=c,
                                  "d"=d,
                                  "raw.p"=test$p.value,
                                  "or"=test$estimate,
                                  "adjust.p"=p.adjust(test$p.value,method = "bonferroni",n=length(mysig_list))
                       )
      )
      
      #relicated vs. unrep at same
      a <- sum((t.real&t.replicated&t.same)&t.sig&exp)
      b <- sum((t.real&t.unrep&t.same)&t.sig&exp)
      c <- sum((t.random&t.replicated&t.same)&t.sig&exp)
      d <- sum((t.random&t.unrep&t.same)&t.sig&exp)
      test.matrix <- matrix(data = c(a,b,c,d),nrow = 2)
      test <- fisher.test(test.matrix)
      test_df <- rbind(test_df,
                       data.frame("signature"=signame,
                                  "test"="rep.vs.unrep.same",
                                  "expression"=exp_name,
                                  "a"=a,
                                  "b"=b,
                                  "c"=c,
                                  "d"=d,
                                  "raw.p"=test$p.value,
                                  "or"=test$estimate,
                                  "adjust.p"=p.adjust(test$p.value,method = "bonferroni",n=length(mysig_list))
                       )
      )
    }
  }
  
  for (i in 1:length(mysig_list))  {
    signame <- names(mysig_list)[i]
    t.sig <- mysig_list[[i]]
    
    # Q34/Q12 in match
    a=sum(t.real&t.replicated&t.sig&t.match34)
    b=sum(t.real&t.replicated&t.sig&t.match12)
    c=sum(t.real&t.unrep&t.sig&t.match34)
    d=sum(t.real&t.unrep&t.sig&t.match12)
    test.matrix <- matrix(data = c(a,b,c,d),nrow = 2)
    test <- fisher.test(test.matrix)
    test_df <- rbind(test_df,
                     data.frame("signature"=signame,
                                "test"="rep.vs.unrep",
                                "expression"="match34.vs.match12",
                                "a"=a,
                                "b"=b,
                                "c"=c,
                                "d"=d,
                                "raw.p"=test$p.value,
                                "or"=test$estimate,
                                "adjust.p"=p.adjust(test$p.value,method = "bonferroni",n=length(mysig_list))
                     )
    )
    # Q234/Q1 in match
    a=sum(t.real&t.replicated&t.sig&t.match234)
    b=sum(t.real&t.replicated&t.sig&t.match1)
    c=sum(t.real&t.unrep&t.sig&t.match234)
    d=sum(t.real&t.unrep&t.sig&t.match1)
    test.matrix <- matrix(data = c(a,b,c,d),nrow = 2)
    test <- fisher.test(test.matrix)
    test_df <- rbind(test_df,
                     data.frame("signature"=signame,
                                "test"="rep.vs.unrep",
                                "expression"="match234.vs.match1",
                                "a"=a,
                                "b"=b,
                                "c"=c,
                                "d"=d,
                                "raw.p"=test$p.value,
                                "or"=test$estimate,
                                "adjust.p"=p.adjust(test$p.value,method = "bonferroni",n=length(mysig_list))
                     )
    )
    # Q34/Q12 in pcawg normal
    a=sum(t.real&t.replicated&t.sig&t.normalpcawg34)
    b=sum(t.real&t.replicated&t.sig&t.normalpcawg12)
    c=sum(t.real&t.unrep&t.sig&t.normalpcawg34)
    d=sum(t.real&t.unrep&t.sig&t.normalpcawg12)
    test.matrix <- matrix(data = c(a,b,c,d),nrow = 2)
    test <- fisher.test(test.matrix)
    test_df <- rbind(test_df,
                     data.frame("signature"=signame,
                                "test"="rep.vs.unrep",
                                "expression"="normalpcawg34.vs.normalpcawg12",
                                "a"=a,
                                "b"=b,
                                "c"=c,
                                "d"=d,
                                "raw.p"=test$p.value,
                                "or"=test$estimate,
                                "adjust.p"=p.adjust(test$p.value,method = "bonferroni",n=length(mysig_list))
                     )
    )
    # Q234/Q1 in pcawg normal
    a=sum(t.real&t.replicated&t.sig&t.normalpcawg234)
    b=sum(t.real&t.replicated&t.sig&t.normalpcawg1)
    c=sum(t.real&t.unrep&t.sig&t.normalpcawg234)
    d=sum(t.real&t.unrep&t.sig&t.normalpcawg1)
    test.matrix <- matrix(data = c(a,b,c,d),nrow = 2)
    test <- fisher.test(test.matrix)
    test_df <- rbind(test_df,
                     data.frame("signature"=signame,
                                "test"="rep.vs.unrep",
                                "expression"="normalpcawg234.vs.normalpcawg1",
                                "a"=a,
                                "b"=b,
                                "c"=c,
                                "d"=d,
                                "raw.p"=test$p.value,
                                "or"=test$estimate,
                                "adjust.p"=p.adjust(test$p.value,method = "bonferroni",n=length(mysig_list))
                     )
    )
    # Q34/Q12 in tcga normal
    a=sum(t.real&t.replicated&t.sig&t.normaltcga34)
    b=sum(t.real&t.replicated&t.sig&t.normaltcga12)
    c=sum(t.real&t.unrep&t.sig&t.normaltcga34)
    d=sum(t.real&t.unrep&t.sig&t.normaltcga12)
    test.matrix <- matrix(data = c(a,b,c,d),nrow = 2)
    test <- fisher.test(test.matrix)
    test_df <- rbind(test_df,
                     data.frame("signature"=signame,
                                "test"="rep.vs.unrep",
                                "expression"="normaltcga34.vs.normaltcga12",
                                "a"=a,
                                "b"=b,
                                "c"=c,
                                "d"=d,
                                "raw.p"=test$p.value,
                                "or"=test$estimate,
                                "adjust.p"=p.adjust(test$p.value,method = "bonferroni",n=length(mysig_list))
                     )
    )
    # Q234/Q1 in tcga normal
    a=sum(t.real&t.replicated&t.sig&t.normaltcga234)
    b=sum(t.real&t.replicated&t.sig&t.normaltcga1)
    c=sum(t.real&t.unrep&t.sig&t.normaltcga234)
    d=sum(t.real&t.unrep&t.sig&t.normaltcga1)
    test.matrix <- matrix(data = c(a,b,c,d),nrow = 2)
    test <- fisher.test(test.matrix)
    test_df <- rbind(test_df,
                     data.frame("signature"=signame,
                                "test"="rep.vs.unrep",
                                "expression"="normaltcga234.vs.normaltcga1",
                                "a"=a,
                                "b"=b,
                                "c"=c,
                                "d"=d,
                                "raw.p"=test$p.value,
                                "or"=test$estimate,
                                "adjust.p"=p.adjust(test$p.value,method = "bonferroni",n=length(mysig_list))
                     )
    )
  }
  
  write.csv(test_df,file.path(outpath),row.names = F,quote = F)
  
}
 

###################################### RUN #######################################
# PCAWG
raw_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-random/real_random_brkpt_withorder_withexpression_withsig.txt"
raw <- read.delim(raw_path,header = T)
pcawg2583_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/pcawg2583_sampleid_match.csv"
pcawg2583 <- read.csv(pcawg2583_path, fileEncoding="UTF-8-BOM")
raw <-  raw[raw$sample %in% unique(pcawg2583$WGS),]

test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/sv.rt.consistent8.shrink0.tsv",
                   tumortype = NA,
                   input_col_need = c("V7","V12","V13"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv.rt.NA.consistent8.shrink0.simple.fishertest.csv",
                   assignment = "pcawg",
                   keep = pcawg2583)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/sv.rt.consistent8.shrink25000.tsv",
                   tumortype = NA,
                   input_col_need = c("V7","V12","V13"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv.rt.NA.consistent8.shrink25000.simple.fishertest.csv",
                   assignment = "pcawg",
                   keep = pcawg2583)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/sv.rt.MCF7.shrink0.tsv",
                   tumortype = "Breast-AdenoCA",
                   input_col_need = c("V7","V12","V13"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv.rt.Breast-AdenoCA.MCF7.shrink0.simple.fishertest.csv",
                   assignment = "pcawg",
                   keep = pcawg2583)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/sv.rt.MCF7.shrink25000.tsv",
                   tumortype = "Breast-AdenoCA",
                   input_col_need = c("V7","V12","V13"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv.rt.Breast-AdenoCA.MCF7.shrink25000.simple.fishertest.csv",
                   assignment = "pcawg",
                   keep = pcawg2583)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/sv.rt.HepG2.shrink0.tsv",
                   tumortype = "Liver-HCC",
                   input_col_need = c("V7","V12","V13"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv.rt.Liver-HCC.HepG2.shrink0.simple.fishertest.csv",
                   assignment = "pcawg",
                   keep = pcawg2583)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/sv.rt.HepG2.shrink25000.tsv",
                   tumortype = "Liver-HCC",
                   input_col_need = c("V7","V12","V13"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv.rt.Liver-HCC.HepG2.shrink25000.simple.fishertest.csv",
                   assignment = "pcawg",
                   keep = pcawg2583)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/sv.rt.consistent8.shrink0.tsv",
                   tumortype = "Eso-AdenoCA",
                   input_col_need = c("V7","V12","V13"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv.rt.Eso-AdenoCA.consistent8.shrink0.simple.fishertest.csv",
                   assignment = "pcawg",
                   keep = pcawg2583)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/sv.rt.consistent8.shrink25000.tsv",
                   tumortype = "Eso-AdenoCA",
                   input_col_need = c("V7","V12","V13"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv.rt.Eso-AdenoCA.consistent8.shrink25000.simple.fishertest.csv",
                   assignment = "pcawg",
                   keep = pcawg2583)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/sv.rt.consistent8.shrink0.tsv",
                   tumortype = "Ovary-AdenoCA",
                   input_col_need = c("V7","V12","V13"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv.rt.Ovary-AdenoCA.consistent8.shrink0.simple.fishertest.csv",
                   assignment = "pcawg",
                   keep = pcawg2583)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/sv.rt.consistent8.shrink25000.tsv",
                   tumortype = "Ovary-AdenoCA",
                   input_col_need = c("V7","V12","V13"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv.rt.Ovary-AdenoCA.consistent8.shrink25000.simple.fishertest.csv",
                   assignment = "pcawg",
                   keep = pcawg2583)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/sv.rt.consistent8.shrink0.tsv",
                   tumortype = "Panc-AdenoCA",
                   input_col_need = c("V7","V12","V13"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv.rt.Panc-AdenoCA.consistent8.shrink0.simple.fishertest.csv",
                   assignment = "pcawg",
                   keep = pcawg2583)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/sv.rt.consistent8.shrink25000.tsv",
                   tumortype = "Panc-AdenoCA",
                   input_col_need = c("V7","V12","V13"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv.rt.Panc-AdenoCA.consistent8.shrink25000.simple.fishertest.csv",
                   assignment = "pcawg",
                   keep = pcawg2583)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/sv.rt.consistent8.shrink0.tsv",
                   tumortype = "Stomach-AdenoCA",
                   input_col_need = c("V7","V12","V13"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv.rt.Stomach-AdenoCA.consistent8.shrink0.simple.fishertest.csv",
                   assignment = "pcawg",
                   keep = pcawg2583)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/sv.rt.consistent8.shrink25000.tsv",
                   tumortype = "Stomach-AdenoCA",
                   input_col_need = c("V7","V12","V13"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv.rt.Stomach-AdenoCA.consistent8.shrink25000.simple.fishertest.csv",
                   assignment = "pcawg",
                   keep = pcawg2583)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/sv.rt.consistent8.shrink0.tsv",
                   tumortype = "Uterus-AdenoCA",
                   input_col_need = c("V7","V12","V13"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv.rt.Uterus-AdenoCA.consistent8.shrink0.simple.fishertest.csv",
                   assignment = "pcawg",
                   keep = pcawg2583)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/sv.rt.consistent8.shrink25000.tsv",
                   tumortype = "Uterus-AdenoCA",
                   input_col_need = c("V7","V12","V13"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv.rt.Uterus-AdenoCA.consistent8.shrink25000.simple.fishertest.csv",
                   assignment = "pcawg",
                   keep = pcawg2583)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/sv.rt.consistent8.shrink0.tsv",
                   tumortype = "Prost-AdenoCA",
                   input_col_need = c("V7","V12","V13"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv.rt.Prost-AdenoCA.consistent8.shrink0.simple.fishertest.csv",
                   assignment = "pcawg",
                   keep = pcawg2583)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/sv.rt.consistent8.shrink25000.tsv",
                   tumortype = "Prost-AdenoCA",
                   input_col_need = c("V7","V12","V13"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv.rt.Prost-AdenoCA.consistent8.shrink25000.simple.fishertest.csv",
                   assignment = "pcawg",
                   keep = pcawg2583)

test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/sv.rt.consistent10.shrink0.tsv",
                   tumortype = NA,
                   input_col_need = c("V7","V12","V13"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv.rt.NA.consistent10.shrink0.simple.fishertest.csv",
                   assignment = "pcawg",
                   keep = pcawg2583)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/sv.rt.consistent10.shrink25000.tsv",
                   tumortype = NA,
                   input_col_need = c("V7","V12","V13"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv.rt.NA.consistent10.shrink25000.simple.fishertest.csv",
                   assignment = "pcawg",
                   keep = pcawg2583)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/sv.rt.consistent10.shrink0.tsv",
                   tumortype = "Eso-AdenoCA",
                   input_col_need = c("V7","V12","V13"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv.rt.Eso-AdenoCA.consistent10.shrink0.simple.fishertest.csv",
                   assignment = "pcawg",
                   keep = pcawg2583)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/sv.rt.consistent10.shrink25000.tsv",
                   tumortype = "Eso-AdenoCA",
                   input_col_need = c("V7","V12","V13"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv.rt.Eso-AdenoCA.consistent10.shrink25000.simple.fishertest.csv",
                   assignment = "pcawg",
                   keep = pcawg2583)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/sv.rt.consistent10.shrink0.tsv",
                   tumortype = "Ovary-AdenoCA",
                   input_col_need = c("V7","V12","V13"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv.rt.Ovary-AdenoCA.consistent10.shrink0.simple.fishertest.csv",
                   assignment = "pcawg",
                   keep = pcawg2583)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/sv.rt.consistent10.shrink25000.tsv",
                   tumortype = "Ovary-AdenoCA",
                   input_col_need = c("V7","V12","V13"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv.rt.Ovary-AdenoCA.consistent10.shrink25000.simple.fishertest.csv",
                   assignment = "pcawg",
                   keep = pcawg2583)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/sv.rt.consistent10.shrink0.tsv",
                   tumortype = "Panc-AdenoCA",
                   input_col_need = c("V7","V12","V13"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv.rt.Panc-AdenoCA.consistent10.shrink0.simple.fishertest.csv",
                   assignment = "pcawg",
                   keep = pcawg2583)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/sv.rt.consistent10.shrink25000.tsv",
                   tumortype = "Panc-AdenoCA",
                   input_col_need = c("V7","V12","V13"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv.rt.Panc-AdenoCA.consistent10.shrink25000.simple.fishertest.csv",
                   assignment = "pcawg",
                   keep = pcawg2583)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/sv.rt.consistent10.shrink0.tsv",
                   tumortype = "Stomach-AdenoCA",
                   input_col_need = c("V7","V12","V13"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv.rt.Stomach-AdenoCA.consistent10.shrink0.simple.fishertest.csv",
                   assignment = "pcawg",
                   keep = pcawg2583)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/sv.rt.consistent10.shrink25000.tsv",
                   tumortype = "Stomach-AdenoCA",
                   input_col_need = c("V7","V12","V13"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv.rt.Stomach-AdenoCA.consistent10.shrink25000.simple.fishertest.csv",
                   assignment = "pcawg",
                   keep = pcawg2583)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/sv.rt.consistent10.shrink0.tsv",
                   tumortype = "Uterus-AdenoCA",
                   input_col_need = c("V7","V12","V13"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv.rt.Uterus-AdenoCA.consistent10.shrink0.simple.fishertest.csv",
                   assignment = "pcawg",
                   keep = pcawg2583)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/sv.rt.consistent10.shrink25000.tsv",
                   tumortype = "Uterus-AdenoCA",
                   input_col_need = c("V7","V12","V13"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv.rt.Uterus-AdenoCA.consistent10.shrink25000.simple.fishertest.csv",
                   assignment = "pcawg",
                   keep = pcawg2583)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/sv.rt.consistent10.shrink0.tsv",
                   tumortype = "Prost-AdenoCA",
                   input_col_need = c("V7","V12","V13"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv.rt.Prost-AdenoCA.consistent10.shrink0.simple.fishertest.csv",
                   assignment = "pcawg",
                   keep = pcawg2583)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result/sv.rt.consistent10.shrink25000.tsv",
                   tumortype = "Prost-AdenoCA",
                   input_col_need = c("V7","V12","V13"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv.rt.Prost-AdenoCA.consistent10.shrink25000.simple.fishertest.csv",
                   assignment = "pcawg",
                   keep = pcawg2583)

# Hartwig
raw_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/real_random_brkpt_withorder_withexpression_withsig.txt"
raw <- read.delim(raw_path,header = T)
hartwig2367_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/hartwig2367_sampleid_match.csv"
hartwig2367 <- read.csv(hartwig2367_path, fileEncoding="UTF-8-BOM")

test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.consistent8.shrink0.tsv",
                   tumortype = NA,
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.NA.consistent8.shrink0.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.consistent8.shrink25000.tsv",
                   tumortype = NA,
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.NA.consistent8.shrink25000.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.MCF7.shrink0.tsv",
                   tumortype = "Breast",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.Breast.MCF7.shrink0.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.MCF7.shrink25000.tsv",
                   tumortype = "Breast",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.Breast.MCF7.shrink25000.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.HepG2.shrink0.tsv",
                   tumortype = "Liver",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.Liver.HepG2.shrink0.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.HepG2.shrink25000.tsv",
                   tumortype = "Liver",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.Liver.HepG2.shrink25000.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.consistent8.shrink0.tsv",
                   tumortype = "Ovary",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.Ovary.consistent8.shrink0.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367)  
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.consistent8.shrink25000.tsv",
                   tumortype = "Ovary",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.Ovary.consistent8.shrink25000.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367) 
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.consistent8.shrink0.tsv",
                   tumortype = "Uterus",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.Uterus.consistent8.shrink0.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367) 
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.consistent8.shrink25000.tsv",
                   tumortype = "Uterus",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.Uterus.consistent8.shrink25000.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367) 
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.consistent8.shrink0.tsv",
                   tumortype = "Pancreas",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.Pancreas.consistent8.shrink0.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367)  
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.consistent8.shrink25000.tsv",
                   tumortype = "Pancreas",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.Pancreas.consistent8.shrink25000.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367)  
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.consistent8.shrink0.tsv",
                   tumortype = "Esophagus",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.Esophagus.consistent8.shrink0.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367)  
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.consistent8.shrink25000.tsv",
                   tumortype = "Esophagus",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.Esophagus.consistent8.shrink25000.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367) 
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.consistent8.shrink0.tsv",
                   tumortype = "Stomach",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.Stomach.consistent8.shrink0.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367) 
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.consistent8.shrink25000.tsv",
                   tumortype = "Stomach",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.Stomach.consistent8.shrink25000.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367) 
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.consistent8.shrink0.tsv",
                   tumortype = "Gastroesophageal",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.Gastroesophageal.consistent8.shrink0.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367) 
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.consistent8.shrink25000.tsv",
                   tumortype = "Gastroesophageal",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.Gastroesophageal.consistent8.shrink25000.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367) 
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.consistent8.shrink0.tsv",
                   tumortype = "Prostate",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.Prostate.consistent8.shrink0.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.consistent8.shrink25000.tsv",
                   tumortype = "Prostate",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.Prostate.consistent8.shrink25000.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.consistent10.shrink0.tsv",
                   tumortype = NA,
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.NA.consistent10.shrink0.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.consistent10.shrink25000.tsv",
                   tumortype = NA,
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.NA.consistent10.shrink25000.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.consistent10.shrink0.tsv",
                   tumortype = "Ovary",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.Ovary.consistent10.shrink0.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.consistent10.shrink25000.tsv",
                   tumortype = "Ovary",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.Ovary.consistent10.shrink25000.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.consistent10.shrink0.tsv",
                   tumortype = "Uterus",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.Uterus.consistent10.shrink0.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367)  
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.consistent10.shrink25000.tsv",
                   tumortype = "Uterus",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.Uterus.consistent10.shrink25000.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367)  
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.consistent10.shrink0.tsv",
                   tumortype = "Pancreas",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.Pancreas.consistent10.shrink0.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367)  
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.consistent10.shrink25000.tsv",
                   tumortype = "Pancreas",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.Pancreas.consistent10.shrink25000.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367)  
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.consistent10.shrink0.tsv",
                   tumortype = "Esophagus",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.Esophagus.consistent10.shrink0.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.consistent10.shrink25000.tsv",
                   tumortype = "Esophagus",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.Esophagus.consistent10.shrink25000.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.consistent10.shrink0.tsv",
                   tumortype = "Stomach",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.Stomach.consistent10.shrink0.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367) 
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.consistent10.shrink25000.tsv",
                   tumortype = "Stomach",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.Stomach.consistent10.shrink25000.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367) 
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.consistent10.shrink0.tsv",
                   tumortype = "Gastroesophageal",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.Gastroesophageal.consistent10.shrink0.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.consistent10.shrink25000.tsv",
                   tumortype = "Gastroesophageal",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.Gastroesophageal.consistent10.shrink25000.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.consistent10.shrink0.tsv",
                   tumortype = "Prostate",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.Prostate.consistent10.shrink0.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367)  
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/sv.rt.consistent10.shrink25000.tsv",
                   tumortype = "Prostate",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/genomic_features/result_test/sv.rt.Prostate.consistent10.shrink25000.fishertest.csv",
                   assignment = "hartwig",
                   keep = hartwig2367) 

# POG570
raw_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/real_random_brkpt_withorder_withexpression_withsig.txt"
raw <- read.delim(raw_path,header = T)
pog570_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/pog570_sampleid_match.csv"
pog570 <- read.csv(pog570_path, fileEncoding="UTF-8-BOM")

test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/sv.rt.consistent8.shrink0.tsv",
                   tumortype = NA,
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/result_test/sv.rt.NA.consistent8.shrink0.fishertest.csv",
                   assignment = "hartwig",
                   keep = pog570)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/sv.rt.consistent8.shrink25000.tsv",
                   tumortype = NA,
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/result_test/sv.rt.NA.consistent8.shrink25000.fishertest.csv",
                   assignment = "hartwig",
                   keep = pog570)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/sv.rt.MCF7.shrink0.tsv",
                   tumortype = "BRCA",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/result_test/sv.rt.BRCA.MCF7.shrink0.fishertest.csv",
                   assignment = "hartwig",
                   keep = pog570)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/sv.rt.MCF7.shrink25000.tsv",
                   tumortype = "BRCA",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/result_test/sv.rt.BRCA.MCF7.shrink25000.fishertest.csv",
                   assignment = "hartwig",
                   keep = pog570)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/sv.rt.consistent8.shrink0.tsv",
                   tumortype = "OV",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/result_test/sv.rt.OV.consistent8.shrink0.fishertest.csv",
                   assignment = "hartwig",
                   keep = pog570)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/sv.rt.consistent8.shrink25000.tsv",
                   tumortype = "OV",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/result_test/sv.rt.OV.consistent8.shrink25000.fishertest.csv",
                   assignment = "hartwig",
                   keep = pog570)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/sv.rt.consistent8.shrink0.tsv",
                   tumortype = "UCEC",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/result_test/sv.rt.UCEC.consistent8.shrink0.fishertest.csv",
                   assignment = "hartwig",
                   keep = pog570)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/sv.rt.consistent8.shrink25000.tsv",
                   tumortype = "UCEC",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/result_test/sv.rt.UCEC.consistent8.shrink25000.fishertest.csv",
                   assignment = "hartwig",
                   keep = pog570)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/sv.rt.consistent8.shrink0.tsv",
                   tumortype = "PANC",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/result_test/sv.rt.PANC.consistent8.shrink0.fishertest.csv",
                   assignment = "hartwig",
                   keep = pog570)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/sv.rt.consistent8.shrink25000.tsv",
                   tumortype = "PANC",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/result_test/sv.rt.PANC.consistent8.shrink25000.fishertest.csv",
                   assignment = "hartwig",
                   keep = pog570)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/sv.rt.consistent8.shrink0.tsv",
                   tumortype = "ESCA",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/result_test/sv.rt.ESCA.consistent8.shrink0.fishertest.csv",
                   assignment = "hartwig",
                   keep = pog570)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/sv.rt.consistent8.shrink25000.tsv",
                   tumortype = "ESCA",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/result_test/sv.rt.ESCA.consistent8.shrink25000.fishertest.csv",
                   assignment = "hartwig",
                   keep = pog570)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/sv.rt.consistent8.shrink0.tsv",
                   tumortype = "STAD",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/result_test/sv.rt.STAD.consistent8.shrink0.fishertest.csv",
                   assignment = "hartwig",
                   keep = pog570)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/sv.rt.consistent8.shrink25000.tsv",
                   tumortype = "STAD",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/result_test/sv.rt.STAD.consistent8.shrink25000.fishertest.csv",
                   assignment = "hartwig",
                   keep = pog570)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/sv.rt.consistent10.shrink0.tsv",
                   tumortype = NA,
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/result_test/sv.rt.NA.consistent10.shrink0.fishertest.csv",
                   assignment = "hartwig",
                   keep = pog570)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/sv.rt.consistent10.shrink25000.tsv",
                   tumortype = NA,
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/result_test/sv.rt.NA.consistent10.shrink25000.fishertest.csv",
                   assignment = "hartwig",
                   keep = pog570)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/sv.rt.consistent10.shrink0.tsv",
                   tumortype = "OV",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/result_test/sv.rt.OV.consistent10.shrink0.fishertest.csv",
                   assignment = "hartwig",
                   keep = pog570)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/sv.rt.consistent10.shrink25000.tsv",
                   tumortype = "OV",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/result_test/sv.rt.OV.consistent10.shrink25000.fishertest.csv",
                   assignment = "hartwig",
                   keep = pog570)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/sv.rt.consistent10.shrink0.tsv",
                   tumortype = "UCEC",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/result_test/sv.rt.UCEC.consistent10.shrink0.fishertest.csv",
                   assignment = "hartwig",
                   keep = pog570)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/sv.rt.consistent10.shrink25000.tsv",
                   tumortype = "UCEC",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/result_test/sv.rt.UCEC.consistent10.shrink25000.fishertest.csv",
                   assignment = "hartwig",
                   keep = pog570)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/sv.rt.consistent10.shrink0.tsv",
                   tumortype = "PANC",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/result_test/sv.rt.PANC.consistent10.shrink0.fishertest.csv",
                   assignment = "hartwig",
                   keep = pog570)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/sv.rt.consistent10.shrink25000.tsv",
                   tumortype = "PANC",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/result_test/sv.rt.PANC.consistent10.shrink25000.fishertest.csv",
                   assignment = "hartwig",
                   keep = pog570)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/sv.rt.consistent10.shrink0.tsv",
                   tumortype = "ESCA",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/result_test/sv.rt.ESCA.consistent10.shrink0.fishertest.csv",
                   assignment = "hartwig",
                   keep = pog570)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/sv.rt.consistent10.shrink25000.tsv",
                   tumortype = "ESCA",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/result_test/sv.rt.ESCA.consistent10.shrink25000.fishertest.csv",
                   assignment = "hartwig",
                   keep = pog570)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/sv.rt.consistent10.shrink0.tsv",
                   tumortype = "STAD",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/result_test/sv.rt.STAD.consistent10.shrink0.fishertest.csv",
                   assignment = "hartwig",
                   keep = pog570)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/sv.rt.consistent10.shrink25000.tsv",
                   tumortype = "STAD",
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PG570/genomic_features/result_test/sv.rt.STAD.consistent10.shrink25000.fishertest.csv",
                   assignment = "hartwig",
                   keep = pog570)

# BRCA-EU
raw_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/BRCA-EU/real_random_brkpt_withorder_withexpression_withsig.txt"
raw <- read.delim(raw_path,header = T)
brcaeu532_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/BRCA-EU/brcaeu532_sampleid_match.csv"
brcaeu532 <- read.csv(brcaeu532_path, fileEncoding="UTF-8-BOM")

test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/BRCA-EU/genomic_features/sv.rt.MCF7.shrink0.tsv",
                   tumortype = NA,
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/BRCA-EU/genomic_features/result_test/sv.rt.NA.MCF7.shrink0.fishertest.csv",
                   assignment = "pcawg",
                   keep = brcaeu532)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/BRCA-EU/genomic_features/sv.rt.MCF7.shrink25000.tsv",
                   tumortype = NA,
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/BRCA-EU/genomic_features/result_test/sv.rt.NA.MCF7.shrink25000.fishertest.csv",
                   assignment = "pcawg",
                   keep = brcaeu532)

# BRCA-FR
raw_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/BRCA-FR/real_random_brkpt_withorder_withexpression_withsig.txt"
raw <- read.delim(raw_path,header = T)
brcafr72_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/BRCA-FR/brcafr72_sampleid_match.csv"
brcafr72 <- read.csv(brcafr72_path, fileEncoding="UTF-8-BOM")

test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/BRCA-FR/genomic_features/sv.rt.MCF7.shrink0.tsv",
                   tumortype = NA,
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/BRCA-FR/genomic_features/result_test/sv.rt.NA.MCF7.shrink0.fishertest.csv",
                   assignment = "pcawg",
                   keep = brcafr72)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/BRCA-FR/genomic_features/sv.rt.MCF7.shrink25000.tsv",
                   tumortype = NA,
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/BRCA-FR/genomic_features/result_test/sv.rt.NA.MCF7.shrink25000.fishertest.csv",
                   assignment = "pcawg",
                   keep = brcafr72)

# OV-AU
raw_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/OV-AU/real_random_brkpt_withorder_withexpression_withsig.txt"
raw <- read.delim(raw_path,header = T)
ovau115_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/OV-AU/ovau115_sampleid_match.csv"
ovau115 <- read.csv(ovau115_path, fileEncoding="UTF-8-BOM")

test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/OV-AU/genomic_features/sv.rt.consistent8.shrink0.tsv",
                   tumortype = NA,
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/OV-AU/genomic_features/result_test/sv.rt.NA.consistent8.shrink0.fishertest.csv",
                   assignment = "pcawg",
                   keep = ovau115)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/OV-AU/genomic_features/sv.rt.consistent8.shrink25000.tsv",
                   tumortype = NA,
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/OV-AU/genomic_features/result_test/sv.rt.NA.consistent8.shrink25000.fishertest.csv",
                   assignment = "pcawg",
                   keep = ovau115)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/OV-AU/genomic_features/sv.rt.consistent10.shrink0.tsv",
                   tumortype = NA,
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/OV-AU/genomic_features/result_test/sv.rt.NA.consistent10.shrink0.fishertest.csv",
                   assignment = "pcawg",
                   keep = ovau115)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/OV-AU/genomic_features/sv.rt.consistent10.shrink25000.tsv",
                   tumortype = NA,
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/OV-AU/genomic_features/result_test/sv.rt.NA.consistent10.shrink25000.fishertest.csv",
                   assignment = "pcawg",
                   keep = ovau115)

# LIRI-JP
raw_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/LIRI-JP/real_random_brkpt_withorder_withexpression_withsig.txt"
raw <- read.delim(raw_path,header = T)
lirijp255_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/LIRI-JP/lirijp255_sampleid_match.csv"
lirijp255 <- read.csv(lirijp255_path, fileEncoding="UTF-8-BOM")

test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/LIRI-JP/genomic_features/sv.rt.HepG2.shrink0.tsv",
                   tumortype = NA,
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/LIRI-JP/genomic_features/result_test/sv.rt.NA.HepG2.shrink0.fishertest.csv",
                   assignment = "pcawg",
                   keep = lirijp255)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/LIRI-JP/genomic_features/sv.rt.HepG2.shrink25000.tsv",
                   tumortype = NA,
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/LIRI-JP/genomic_features/result_test/sv.rt.NA.HepG2.shrink25000.fishertest.csv",
                   assignment = "pcawg",
                   keep = lirijp255)

# PACA-AU
raw_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PACA-AU/real_random_brkpt_withorder_withexpression_withsig.txt"
raw <- read.delim(raw_path,header = T)
pacaau178_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PACA-AU/pacaau178_sampleid_match.csv"
pacaau178 <- read.csv(pacaau178_path, fileEncoding="UTF-8-BOM")

test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PACA-AU/genomic_features/sv.rt.consistent8.shrink0.tsv",
                   tumortype = NA,
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PACA-AU/genomic_features/result_test/sv.rt.NA.consistent8.shrink0.fishertest.csv",
                   assignment = "pcawg",
                   keep = pacaau178)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PACA-AU/genomic_features/sv.rt.consistent8.shrink25000.tsv",
                   tumortype = NA,
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PACA-AU/genomic_features/result_test/sv.rt.NA.consistent8.shrink25000.fishertest.csv",
                   assignment = "pcawg",
                   keep = pacaau178)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PACA-AU/genomic_features/sv.rt.consistent10.shrink0.tsv",
                   tumortype = NA,
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PACA-AU/genomic_features/result_test/sv.rt.NA.consistent10.shrink0.fishertest.csv",
                   assignment = "pcawg",
                   keep = pacaau178)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PACA-AU/genomic_features/sv.rt.consistent10.shrink25000.tsv",
                   tumortype = NA,
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PACA-AU/genomic_features/result_test/sv.rt.NA.consistent10.shrink25000.fishertest.csv",
                   assignment = "pcawg",
                   keep = pacaau178)

# PRAD-CA
raw_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PRAD-CA/real_random_brkpt_withorder_withexpression_withsig.txt"
raw <- read.delim(raw_path,header = T)
pradca327_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PRAD-CA/pradca327_sampleid_match.csv"
pradca327 <- read.csv(pradca327_path, fileEncoding="UTF-8-BOM")

test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PRAD-CA/genomic_features/sv.rt.consistent8.shrink0.tsv",
                   tumortype = NA,
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PRAD-CA/genomic_features/result_test/sv.rt.NA.consistent8.shrink0.fishertest.csv",
                   assignment = "pcawg",
                   keep = pradca327)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PRAD-CA/genomic_features/sv.rt.consistent8.shrink25000.tsv",
                   tumortype = NA,
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PRAD-CA/genomic_features/result_test/sv.rt.NA.consistent8.shrink25000.fishertest.csv",
                   assignment = "pcawg",
                   keep = pradca327)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PRAD-CA/genomic_features/sv.rt.consistent10.shrink0.tsv",
                   tumortype = NA,
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PRAD-CA/genomic_features/result_test/sv.rt.NA.consistent10.shrink0.fishertest.csv",
                   assignment = "pcawg",
                   keep = pradca327)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PRAD-CA/genomic_features/sv.rt.consistent10.shrink25000.tsv",
                   tumortype = NA,
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PRAD-CA/genomic_features/result_test/sv.rt.NA.consistent10.shrink25000.fishertest.csv",
                   assignment = "pcawg",
                   keep = pradca327)

# PRAD-UK
raw_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PRAD-UK/real_random_brkpt_withorder_withexpression_withsig.txt"
raw <- read.delim(raw_path,header = T)
praduk269_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PRAD-UK/praduk269_sampleid_match_unique.csv"
praduk269 <- read.csv(praduk269_path, fileEncoding="UTF-8-BOM")
praduk269 <- praduk269[praduk269$keep=="yes",1:3]

test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PRAD-UK/genomic_features/sv.rt.consistent8.shrink0.tsv",
                   tumortype = NA,
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PRAD-UK/genomic_features/result_test/sv.rt.NA.consistent8.shrink0.fishertest.csv",
                   assignment = "pcawg",
                   keep = praduk269)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PRAD-UK/genomic_features/sv.rt.consistent8.shrink25000.tsv",
                   tumortype = NA,
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PRAD-UK/genomic_features/result_test/sv.rt.NA.consistent8.shrink25000.fishertest.csv",
                   assignment = "pcawg",
                   keep = praduk269)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PRAD-UK/genomic_features/sv.rt.consistent10.shrink0.tsv",
                   tumortype = NA,
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PRAD-UK/genomic_features/result_test/sv.rt.NA.consistent10.shrink0.fishertest.csv",
                   assignment = "pcawg",
                   keep = praduk269)
test_bias_function(raw,
                   input_path = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PRAD-UK/genomic_features/sv.rt.consistent10.shrink25000.tsv",
                   tumortype = NA,
                   input_col_need = c("V15","V20","V21"),
                   outpath = "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/PRAD-UK/genomic_features/result_test/sv.rt.NA.consistent10.shrink25000.fishertest.csv",
                   assignment = "pcawg",
                   keep = praduk269)
