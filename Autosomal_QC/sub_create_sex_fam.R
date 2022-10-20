
#################### update family file with sex information ###################
library(data.table)
library(tidyverse)
library(optparse)
#########################################################################################################
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Path to plink files, it assumes a bed, bim and fam file with the same file name", metavar="character"),
  make_option(c("-f", "--famfile"), type="character", default=NULL, 
              help="complete path of the pedigree file"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="complete path of the pedigree file" )
); 
opt_parser  <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
out <- file.path(opt$output)
################################################################################
ped<-fread(opt$famfile,data.table=F,colClasses = c(rep("character", 6)))
fam<-fread(opt$input,data.table=F,colClasses = c(rep("character", 6))
fam$sex<-ped$V5[match(fam$V2,ped$V2)]
newfam<-fam[,c(1,2,3,4,7,6)]
write.table(newfam,paste0(out,"/sex_info.fam"),quote = F,row.names = F,col.names = F)
#### done
