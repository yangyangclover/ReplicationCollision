### LIBRARY
library(parallel)

### FUNCTIONS
strand_function <- function(input,rtfile,cellline="HCT116"){
  names(rtfile)[1:4] <- c("chr","start","end","strand")
  print("raw file")
  print(head(input))
  #add  direction
  input <- merge(input,
                       rtfile[,c("gene_name","ensembl_gene_id",cellline)],
                       by.x="V21",
                       by.y="ensembl_gene_id",
                       sort=F,
                       all.x=T)

  print("merged")
  print(head(input))
  ######################### add direction ##############
  input$direction  <- 
    ifelse(
      is.na(input[cellline]),
      NA,
      ifelse((
        as.character(input$V19) %in% "+" &
          as.character(input[,cellline])  %in% "right"
      ) |
        (as.character(input$V19) %in% "-" &
           as.character(input[,cellline])  %in% "left"),
      "same",
      "oppo"
      )
    )
  print("direction added")
  print(head(input))
  print(unique(input$direction))
  ######################## add prime ###########################
  input$prime <-
    ifelse(input$V19 %in% c("-","+"),
           ifelse((
             input$V19 %in% "+" &
               input$V4 %in% "+"
           ) |
             (
               input$V19 %in% "-" &
                 input$V4 %in% "-"
             ),
           5,
           3
           ),
           NA
           )
  print("prime added")
  print(head(input))
  print(unique(input$prime))
  ######################## add replicated / unreplicated annotation ###########################
  #  filters
  oppo_filter  <-  input$direction %in% c("oppo")
  same_filter  <-  input$direction %in% c("same")
  prime3_filter <- input$prime  %in% c(3)
  prime5_filter <- input$prime  %in% c(5)

  left_filter <- as.character(input[,cellline])  %in% "left"
  right_filter <- as.character(input[,cellline])  %in% "right"
  bp_or1_filter <- input$V4 %in% "+"
  bp_or0_filter <- input$V4 %in% "-"

  # replicated
  replicated_tell <- ifelse((left_filter & bp_or0_filter) | (right_filter & bp_or1_filter),"replicated","")
  # unreplicated
  unreplicated_tell <- ifelse((left_filter & bp_or1_filter) | (right_filter & bp_or0_filter), "unreplicated","")

  # paste replicated annotation
  input$replicated=paste0(replicated_tell,unreplicated_tell)
  # order
  input <- input[order(input$V15,decreasing = F),]
  return(input)
}

############################################ RUN
args <- commandArgs(T)

input_path <- args[1] #/path/to/real_random_brkpt_withorder.gene.bed
rtpath  <- args[2] #/path/to/coding_geng_all_rtcellline_direction_slope5e-07_cellline8_shrink25000_flatonbin.tsv
cellline <- args[3] #consistent/other cell line name
ths <- as.numeric(as.character(args[4])) #thread to speed up
output <- args[5] #output path


# read files
input=read.delim(input_path ,header=F )
rtfile=read.delim(rtpath)
rtfile_names=names(rtfile)

# annotate
if (!cellline %in% names(rtfile)){
    print("check your cellline name")
  } else {

clus <- makeCluster(ths)
envir3 <- environment(strand_function )
clusterExport(clus, varlist = ls(envir3), envir = envir3) 

input <- strand_function(input, rtfile,cellline)
write.table(input,output,quote=FALSE,sep="\t",col.names = T,row.names = F)
}
