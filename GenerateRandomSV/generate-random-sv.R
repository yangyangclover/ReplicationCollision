# This is a script used to generate random SVs for validation
# Yang Yang
# 10/25/2021

# library
library(dplyr,lib.loc="/gpfs/data/lyang-lab/users/yangyang/R/lib")
library(BSgenome.Hsapiens.UCSC.hg19,lib.loc="/gpfs/data/lyang-lab/users/yangyang/R/lib")
library(BSgenome.Hsapiens.UCSC.hg19.masked,lib.loc="/gpfs/data/lyang-lab/users/yangyang/R/lib")

args <-commandArgs(TRUE)

sampleid <- args[1]
annotation_path <- args[2]
randomdir <- args[3]
mappablepath <- args[4]
num.regions <- as.numeric(args[5])


if (dir.exists(randomdir)==FALSE){
  dir.create(randomdir,recursive = T,showWarnings = FALSE)
}
############################ read files #######################
all_annotation <- read.delim(annotation_path,header=T)

## ===================== Load unmappable region =====================
genome <- getBSgenome(BSgenome.Hsapiens.UCSC.hg19.masked, masked=TRUE)
# mappable file was downloaded from ucsc:https://bismap.hoffmanlab.org/
# hg19 k100 single read
mappable <- read.delim(file.path(mappablepath),header = F)
names(mappable) <- c("chr","start","end","k100","n1","n2")
# generate unmappable region
unmappable <- data.frame(matrix(NA,nrow = 0,ncol = 3))
# generate mappable file chr by chr
for (i in unique(mappable$chr)){
  chrregion <- mappable[mappable$chr == i,]
  unmappable_add <- data.frame("start"=c(0,chrregion$end),"end"=c(chrregion$start,seqlengths(genome)[i]))
  unmappable_add$chr <- i
  unmappable <- rbind(unmappable,unmappable_add)
}
unmappable<-unmappable[,c("chr","start","end")]


## ==================== FUNCTION generate random sv for intra sv ===================
# function: HowManyPairsFail
HowManyPairsFail <- function(site,unmappable_chr) {
  tell1 <- sum((site[1] > unmappable_chr$start)+(site[1] < unmappable_chr$end)  == 2) # site in the region will be 2, if not in any region, sum will be 0
  tell2 <- sum((site[2] > unmappable_chr$start)+(site[2] < unmappable_chr$end)  == 2) # site in the region will be 2, if not in any region, sum will be 0
  tell <- tell1 + tell2
  # if tell is 0, means all SV pairs are passed (not in unmappable region)
  if (tell > 0 ) {
    out <- "FAILE"
  } else {
    out <- "PASS"
  }
  return(out)
}
HowManyTRAFail <- function(site,unmappable_chr1,unmappable_chr2){
  tell1 <- sum((site[1] > unmappable_chr1$start)+(site[1] < unmappable_chr1$end)  == 2) # site in the region will be 2, if not in any region, sum will be 0
  tell2 <- sum((site[2] > unmappable_chr2$start)+(site[2] < unmappable_chr2$end)  == 2) # site in the region will be 2, if not in any region, sum will be 0
  tell <- tell1 + tell2
  # if tell is 0, means all SV pairs are passed (not in unmappable region)
  if (tell > 0 ) {
    out <- "FAILE"
  } else {
    out <- "PASS"
  }
  return(out)
}



# function: generate_randomsv
# input file
mysv <- all_annotation[all_annotation$sample ==  sampleid,]
#set.seed(123456)

# generate intra
mysv_intra <- mysv[as.character(mysv$chrom1) == as.character(mysv$chrom2),]

new.regions.intra <- data.frame(matrix(NA, nrow = 0, ncol = ncol(mysv_intra)))
names(new.regions.intra) <- names(mysv_intra)

mysv_intra_new <- data.frame(matrix(NA,nrow = 0,ncol=ncol(mysv_intra)))
names(mysv_intra_new) <- names(mysv_intra)

# generate random sv line by line
if (nrow(mysv_intra)>=1){
  for (i in 1:nrow(mysv_intra)) {
    size <- abs(mysv_intra$start2[i]-mysv_intra$start1[i])
    chrom <- paste0("chr",mysv_intra$chrom1[i])
    unmappable_chr <- unmappable[unmappable$chr == chrom, ]
    judge_sv <- apply(cbind(mysv_intra[i,"start1"],mysv_intra[i,"start2"]), 1,  HowManyPairsFail,unmappable_chr)
    if (judge_sv == "PASS") {
      mysv_intra_new <- rbind(mysv_intra_new,mysv_intra[i,])
      judge <- FALSE
      while (judge  == FALSE){
        # generate site1 randomly
        site1 <- sample(seqlengths(genome)[chrom]-size,num.regions,replace = F)
        # site2 is site1 + size
        site2 <- site1 + size
        site <- cbind(site1,site2)
        judge <- sum(apply(site, 1, HowManyPairsFail,unmappable_chr)  == "PASS") == num.regions
      }
      new.regions.add <- mysv_intra[rep(i,each=num.regions),]
      new.regions.add$start1 <- site1
      new.regions.add$start2 <- site2
      new.regions.add$end1 <- site1 + 1
      new.regions.add$end2 <- site2 + 1
      new.regions.intra <- rbind(new.regions.intra,new.regions.add)
    }
  }
}
# generate inter
mysv_inter <- mysv[as.character(mysv$chrom1) != as.character(mysv$chrom2),]
new.regions.inter <- data.frame(matrix(NA,nrow = 0,ncol = ncol(mysv_inter)))
names(new.regions.inter) <- names(mysv_inter)
mysv_inter_new <- data.frame(matrix(NA,nrow = 0,ncol=ncol(mysv_inter)))
names(mysv_inter_new) <- names(mysv_inter)
# generate random sv line by line
if (nrow(mysv_inter)>=1){
  for (i in 1:nrow(mysv_inter)) {
    chrom1 <- paste0("chr",mysv_inter$chrom1[i])
    chrom2 <- paste0("chr",mysv_inter$chrom2[i])
    unmappable_chr1 <- unmappable[unmappable$chr == chrom1, ]
    unmappable_chr2 <- unmappable[unmappable$chr == chrom2, ]
    
    judge_sv <- apply(cbind(mysv_inter[i,"start1"],mysv_inter[i,"start2"]), 1,  HowManyTRAFail,unmappable_chr1,unmappable_chr2)
    if (judge_sv == "PASS"){
      mysv_inter_new <- rbind(mysv_inter_new,mysv_inter[i,])
      judge <- FALSE
      while (judge  == FALSE){
        # generate site1 randomly
        site1 <- sample(seqlengths(genome)[chrom1],num.regions,replace = F)
        # generate site2 randomly
        site2 <- sample(seqlengths(genome)[chrom2],num.regions,replace = F)
        site <- cbind(site1,site2)
        judge <- sum(apply(site, 1, HowManyTRAFail,unmappable_chr1,unmappable_chr2)  == "PASS") == num.regions
      }
      
      new.regions.add <- mysv_inter[rep(i,each=num.regions),]
      new.regions.add$start1 <- site1
      new.regions.add$start2 <- site2
      new.regions.add$end1 <- site1 + 1
      new.regions.add$end2 <- site2 + 1
      new.regions.inter <- rbind(new.regions.inter,new.regions.add)
    }
  }
}

random_sv <- rbind(new.regions.intra,new.regions.inter)
random_sv$treatment <- "random"


random_sv$size <-  ifelse(as.character(random_sv$chrom1)==as.character(random_sv$chrom2),abs(random_sv$start1-random_sv$start2),0)

unrandom_sv <- rbind(mysv_intra_new,mysv_inter_new)
unrandom_sv$treatment <- "real"
unrandom_sv$size <-  ifelse(as.character(unrandom_sv$chrom1)==as.character(unrandom_sv$chrom2),abs(unrandom_sv$start1-unrandom_sv$start2),0)

write.table(random_sv,file.path(randomdir,paste0(sampleid,"_rand.tsv")),sep = "\t",row.names = F,quote = F)
write.table(unrandom_sv,file.path(randomdir,paste0(sampleid,"_real.tsv")),sep = "\t",row.names = F,quote = F)
