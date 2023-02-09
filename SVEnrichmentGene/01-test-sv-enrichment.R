### function
fisher_test_df_function <- function(x1,x2,y1,y2,sig,signame,testname,expressionname,orginaldf){
  a=sum(x1&y1&sig)
  b=sum(x1&y2&sig)
  c=sum(x2&y1&sig)
  d=sum(x2&y2&sig)
  fisher_test <- fisher.test(matrix(c(a,b,c,d),nrow=2,ncol=2))
  fisher_test_df <-  data.frame("signature"=signame,
                                "test"=testname,
                                "expression"=expressionname,
                                "a"=a,
                                "b"=b,
                                "c"=c,
                                "d"=d,
                                "raw.p"=fisher_test$p.value,
                                "or"=fisher_test$estimate)
  orginaldf <- rbind(orginaldf,fisher_test_df)
  return(orginaldf)
}

svenrichment_function <-  function(raw,genecolname,testnormal,assignment){
  raw_input <- raw
  
  # delete duplicated, if one breakpoint in 2 genes, only keep one
  raw_input <- raw_input[!duplicated(raw_input$order),]
  
  if (assignment=="pcawg"){
    t.del1 <-  raw_input$sig_group  %in% paste0("del_",seq(1,3,1))
    t.del2 <-  raw_input$sig_group  %in% paste0("del_",seq(4,6,1))
    t.del3 <-  raw_input$sig_group  %in% paste0("del_",seq(7,11,1))
    t.del4 <-  raw_input$sig_group  %in% paste0("del_",seq(12,15,1))
    t.dup1  <-  raw_input$sig_group  %in% paste0("td_",seq(1,5,1))
    t.dup2  <-  raw_input$sig_group  %in% paste0("td_",seq(6,9,1))
    t.dup3  <-  raw_input$sig_group  %in% paste0("td_",seq(10,11,1))
    t.dup4 <- raw_input$sig_group  %in% paste0("td_",seq(12,16,1))
    t.fragile <- raw_input$sig_group  %in% c("td_fragile","del_fragile")
    t.fb <- raw_input$sig_group  %in% c("fb_inv")
    t.unbalinv <- raw_input$sig_group  %in% c("unbal_inv_1","unbal_inv_2","unbal_inv_3")
    t.unbaltra <- raw_input$sig_group  %in% c("unbal_tra")
    t.largemixed <- raw_input$sig_group  %in% c("recip_tra",
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
                       t.largemixed,
                       TRUE
    )
    names(mysig_list) <- c("del1","del2","del3","del4",
                           "dup1","dup2","dup3","dup4",
                           "fragile",
                           "fbinv",
                           "unbalinv",
                           "unbaltra",
                           "largemixed",
                           "all")
  }
  
  if (assignment=="hartwig"){
    t.del0 <-  raw_input$sig_group  %in% paste0("del_",1)
    t.del1 <-  raw_input$sig_group  %in% paste0("del_",seq(2,3,1))
    t.del2 <-  raw_input$sig_group  %in% paste0("del_",seq(4,6,1))
    t.del3 <-  raw_input$sig_group  %in% paste0("del_",seq(7,11,1))
    t.del4 <-  raw_input$sig_group  %in% paste0("del_",seq(12,15,1))
    
    t.dup0  <-  raw_input$sig_group  %in% paste0("td_",1)
    t.dup1  <-  raw_input$sig_group  %in% paste0("td_",seq(2,5,1))
    t.dup2  <-  raw_input$sig_group  %in% paste0("td_",seq(6,9,1))
    t.dup3  <-  raw_input$sig_group  %in% paste0("td_",seq(10,12,1))
    t.dup4 <- raw_input$sig_group  %in% paste0("td_",seq(13,16,1))
    
    t.fragile <- raw_input$sig_group  %in% c("td_fragile","del_fragile")
    t.fb <- raw_input$sig_group  %in% c("fb_inv","unbal_inv_1","unbal_inv_2")
    t.unbaltra <- raw_input$sig_group  %in% c("unbal_tra")
    t.largemixed <- raw_input$sig_group  %in% c("unbal_inv_3",
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
                       t.largemixed,
                       TRUE
    )
    names(mysig_list) <- c("del0","del1","del2","del3","del4",
                           "dup0","dup1","dup2","dup3","dup4",
                           "fragile",
                           "fbinv",
                           "unbaltra",
                           "largemixed",
                           "all")
  }
  
  
  t.real <- raw_input$treatment %in% "real"
  t.random <- raw_input$treatment %in% "random"
  t.match_q1234gene <-  raw_input$match_exp %in% seq(0,100,1)
  t.match_q1gene <-  raw_input$match_exp %in% seq(0,24,1)
  t.match_q2gene <-  raw_input$match_exp %in% seq(25,49,1)
  t.match_q3gene <-  raw_input$match_exp %in% seq(50,74,1)
  t.match_q4gene <-  raw_input$match_exp %in% seq(75,100,1)
  t.match_q12gene <-  raw_input$match_exp %in% seq(0,49,1)
  t.match_q34gene <-  raw_input$match_exp %in% seq(50,100,1)
  t.match_q234gene <-  raw_input$match_exp %in% seq(25,100,1)
  t.normaltcga_q1234gene <-  raw_input$normal_tcga_exp %in% seq(0,100,1)
  t.normaltcga_q1gene <-  raw_input$normal_tcga_exp %in% seq(0,24,1)
  t.normaltcga_q2gene <-  raw_input$normal_tcga_exp %in% seq(25,49,1)
  t.normaltcga_q3gene <-  raw_input$normal_tcga_exp %in% seq(50,74,1)
  t.normaltcga_q4gene <-  raw_input$normal_tcga_exp %in% seq(75,100,1)
  t.normaltcga_q12gene <-  raw_input$normal_tcga_exp %in% seq(0,49,1)
  t.normaltcga_q34gene <-  raw_input$normal_tcga_exp %in% seq(50,100,1)
  t.normaltcga_q234gene <-  raw_input$normal_tcga_exp %in% seq(25,100,1)
  t.nogene <- raw_input[,genecolname]  %in%  "."
  t.gene <- !t.nogene
  
  test_df <- data.frame()
  
  
  for (i in 1:length(mysig_list))  {
    signame <- names(mysig_list)[i]
    t.sig <- mysig_list[[i]]
  
  test_df <- fisher_test_df_function(t.real,t.random,t.match_q1234gene,t.nogene,t.sig,signame,"sv.enrichment.observed.vs.random","match1234.vs.notgene",test_df)
  test_df <- fisher_test_df_function(t.real,t.random,t.match_q1gene,t.nogene,t.sig,signame,"sv.enrichment.observed.vs.random","match1.vs.notgene",test_df)
  test_df <- fisher_test_df_function(t.real,t.random,t.match_q2gene,t.nogene,t.sig,signame,"sv.enrichment.observed.vs.random","match2.vs.notgene",test_df)
  test_df <- fisher_test_df_function(t.real,t.random,t.match_q3gene,t.nogene,t.sig,signame,"sv.enrichment.observed.vs.random","match3.vs.notgene",test_df)
  test_df <- fisher_test_df_function(t.real,t.random,t.match_q4gene,t.nogene,t.sig,signame,"sv.enrichment.observed.vs.random","match4.vs.notgene",test_df)
  test_df <- fisher_test_df_function(t.real,t.random,t.match_q12gene,t.nogene,t.sig,signame,"sv.enrichment.observed.vs.random","match12.vs.notgene",test_df)
  test_df <- fisher_test_df_function(t.real,t.random,t.match_q34gene,t.nogene,t.sig,signame,"sv.enrichment.observed.vs.random","match34.vs.notgene",test_df)
  test_df <- fisher_test_df_function(t.real,t.random,t.match_q234gene,t.nogene,t.sig,signame,"sv.enrichment.observed.vs.random","match234.vs.notgene",test_df)
  test_df <- fisher_test_df_function(t.real,t.random,t.gene,t.nogene,t.sig,signame,"sv.enrichment.observed.vs.random","gene.vs.notgene",test_df)
  
  if (testnormal==TRUE){
    test_df <- fisher_test_df_function(t.real,t.random,t.normaltcga_q1234gene,t.nogene,t.sig,signame,"sv.enrichment.observed.vs.random","normaltcga1234.vs.notgene",test_df)
    test_df <- fisher_test_df_function(t.real,t.random,t.normaltcga_q1gene,t.nogene,t.sig,signame,"sv.enrichment.observed.vs.random","normaltcga1.vs.notgene",test_df)
    test_df <- fisher_test_df_function(t.real,t.random,t.normaltcga_q2gene,t.nogene,t.sig,signame,"sv.enrichment.observed.vs.random","normaltcga2.vs.notgene",test_df)
    test_df <- fisher_test_df_function(t.real,t.random,t.normaltcga_q34gene,t.nogene,t.sig,signame,"sv.enrichment.observed.vs.random","normaltcga3.vs.notgene",test_df)
    test_df <- fisher_test_df_function(t.real,t.random,t.normaltcga_q4gene,t.nogene,t.sig,signame,"sv.enrichment.observed.vs.random","normaltcga4.vs.notgene",test_df)
    test_df <- fisher_test_df_function(t.real,t.random,t.normaltcga_q12gene,t.nogene,t.sig,signame,"sv.enrichment.observed.vs.random","normaltcga12.vs.notgene",test_df)
    test_df <- fisher_test_df_function(t.real,t.random,t.normaltcga_q34gene,t.nogene,t.sig,signame,"sv.enrichment.observed.vs.random","normaltcga34.vs.notgene",test_df)
    test_df <- fisher_test_df_function(t.real,t.random,t.normaltcga_q234gene,t.nogene,t.sig,signame,"sv.enrichment.observed.vs.random","normaltcga234.vs.notgene",test_df)
  }
  }
  

  return(test_df)
}




################################# RUN #####################
output_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-features/result_test/sv_enrichment"

# PCAWG
raw_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/006-random/real_random_brkpt_withorder_withexpression_withsig.txt"
raw <- read.delim(raw_path,header = T)
pcawg2583_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/pcawg/pcawg2583_sampleid_match.csv"
pcawg2583 <- read.csv(pcawg2583_path, fileEncoding="UTF-8-BOM")
raw <-  raw[raw$sample %in% unique(pcawg2583$WGS),]
pcawg_svenrichment <- svenrichment_function(raw,genecolname="V13",testnormal=TRUE,assignment="pcawg")
write.csv(pcawg_svenrichment,file.path(output_path,"pcawg_sv_enrichment_across_expression.csv"),row.names = F,quote = F)

#Hartwig
raw_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/real_random_brkpt_withorder_withexpression_withsig.txt"
raw <- read.delim(raw_path,header = T)
hartwig2367_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/scratch/large_dup/validation/Hartwig/hartwig2367_sampleid_match.csv"
hartwig2367 <- read.csv(hartwig2367_path, fileEncoding="UTF-8-BOM")
raw <-  raw[raw$sample %in% unique(hartwig2367$WGS),]
hartwig_svenrichment <- svenrichment_function(raw,genecolname="V20",testnormal=TRUE,assignment="hartwig")
write.csv(hartwig_svenrichment,file.path(output_path,"hartwig_sv_enrichment_across_expression.csv"),row.names = F,quote = F)
