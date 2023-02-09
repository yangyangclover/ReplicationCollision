# path
raw_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/pubcode/01-test-strand-bias/02-annotate.strand/example.input/real_random_brkpt_withorder.txt"
real_rand_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/pubcode/01-test-strand-bias/02-annotate.strand/example.input/real_and_random.txt"
fragile_site_path <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/pubcode/01-test-strand-bias/02-annotate.strand/example.input/fragile_site.csv"
outpath <- "D:/OneDrive/OneDrive - The University of Chicago/office/GitHub/pcawg-sv-signature/pubcode/01-test-strand-bias/02-annotate.strand/example.output/real_random_brkpt_withorder_withsig.txt"

# read
fragile_site <- read.csv(fragile_site_path)
real_random <- read.delim(real_rand_path)

del_groups_left <- c(0,
                     1000,2500,5000,7500,
                     10000,25000,50000,75000,
                     100000,250000,500000,750000,
                     1000000,2500000,5000000,7500000,
                     10000000)
del_groups_right <- c(1000,2500,5000,7500,
                      10000,25000,50000,75000,
                      100000,250000,500000,750000,
                      1000000,2500000,5000000,7500000,
                      10000000,1000000000)
dup_groups_left <- c(0,
                     1000,2500,5000,7500,
                     10000,25000,50000,75000,
                     100000,250000,500000,750000,
                     1000000,2500000,5000000,7500000,
                     10000000)
dup_groups_right <- c(1000,2500,5000,7500,
                      10000,25000,50000,75000,
                      100000,250000,500000,750000,
                      1000000,2500000,5000000,7500000,
                      10000000,1000000000)

recip_inv_groups_left <- c(0,5000,50000,500000,5000000)
recip_inv_groups_right <- c(5000,50000,500000,5000000,1000000000)

unbal_inv_groups_left <- c(50000,500000,5000000)
unbal_inv_groups_right <- c(500000,5000000,1000000000)

real_random <- real_random[real_random$treatment=="real",]
for (i in 1:nrow(real_random)){
  chr1 <-  as.character(real_random$chrom1[i])
  chr2 <-  as.character(real_random$chrom2[i])
  pos1 <- real_random$start1[i]
  pos2 <- real_random$start2[i]
  group <- as.character(real_random$FULL_SV_TYPE[i])
  if (group %in% c("del","del_templated_ins")){
    fragile_tell <- (fragile_site$chrom==chr1 & fragile_site$start<= pos1 & fragile_site$end>= pos1)|(fragile_site$chrom==chr2 & fragile_site$start<= pos2 & fragile_site$end>= pos2)
    if (sum(fragile_tell)>0){
      group <- "del_fragile"
    } else {
      size <- abs(pos2-pos1)
      tell <- size >= del_groups_left & size < del_groups_right
      subgroup <- grep("TRUE",tell)
      group <- paste0("del","_",subgroup) 
    }
  }
  if (group %in% c("td")){
    fragile_tell <- (fragile_site$chrom==chr1 & fragile_site$start<= pos1 & fragile_site$end>= pos1)|(fragile_site$chrom==chr2 & fragile_site$start<= pos2 & fragile_site$end>= pos2)
    if (sum(fragile_tell)>0){
      group <- "td_fragile"
    } else {
      size <- abs(pos2-pos1)
      tell <- size >= dup_groups_left & size < dup_groups_right
      subgroup <- grep("TRUE",tell)
      group <- paste0(group,"_",subgroup) 
    }
  }
  if (group %in% c("unbal_inv")){
    size <- abs(pos2-pos1)
    tell <- size >= unbal_inv_groups_left & size < unbal_inv_groups_right
    subgroup <- grep("TRUE",tell)
    group <- paste0(group,"_",subgroup) 
  }
  if (group %in% c("recip_inv_del",
                   "recip_inv_bal",
                   "recip_inv_templated_ins")){
    size <- abs(pos2-pos1)
    tell <- size >= recip_inv_groups_left & size < recip_inv_groups_right
    subgroup <- grep("TRUE",tell)
    group <- paste0("recip_inv","_",subgroup) 
  }
  if (group %in% c("recip_tra_del",
                   "recip_tra_bal",
                   "recip_tra_templated_ins")){
    group <- "recip_tra"
  }
  real_random[i,"sig_group"] <- group
}

# read files
raw <- read.delim(raw_path,header = T)
real_sig <- real_random


# add signature
raw <- merge(raw,
             real_sig[,c("sample","sv_id","sig_group")],
             by.x = c("sample","sv_id"),
             by.y = c("sample","sv_id"),
             all.x = T,
             sort = F)
# write
write.table(raw,
            outpath,
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)
