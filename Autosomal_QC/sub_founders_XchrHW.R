###################################
### Stats on founders and Chromosome X Hardy-weinberg fiter
### date: 12-10-2022
### version: 0.01
### authors: EL - RAG
###################################
### notes
## make sure the .fam file in all the plink files contains information ofthe actual pedigrees so plink can identify the founders
library(data.table)
library(tidyverse)
library(optparse)
library(gridExtra)
library(grid)
library(scales)

##cluster test
# opt <- list()
# opt$plink <- "/groups/umcg-llnext/tmp01/umcg-elopera/ugli_next/merged/merged_commonsnps2/"
# opt$xplink <- "/groups/umcg-llnext/tmp01/umcg-elopera/ugli_next/merged/merged_commonsnps2/chr_X"
# opt$out <- "/groups/umcg-llnext/tmp01/umcg-elopera/ugli_next/merged/"
# opt$plinkexe <- "PLINK/1.9-beta6-20190617"
# opt$pedfamfile < "/groups/umcg-llnext/tmp01/umcg-elopera/ugli_next/merged/merged_commonsnps2/chr_2.fam"

## functions
HW.MAF.dist.plot  <- function(maf.dat, hw.dat, out, name=""){
  chr.labels<-unique(maf.dat$CHR)
  maf.dist.plot.chr <- ggplot(maf.dat, aes(x=MAF))+
    stat_density(aes(color= CHR), position= "identity", geom= "line")+
    ggtitle("MAF distribution per chromosome")+
    scale_color_manual(labels=(chr.labels), values=rainbow(23))+
    theme_bw()+
    theme(text=element_text(size=10, family="Helvetica"))
  
  maf.dist.plot.all <- ggplot(maf.dat, aes(x=MAF))+
    stat_density(color="black",  position= "identity", geom= "line")+
    ggtitle("MAF distribution")+
    theme_bw()+
    theme(text=element_text(size=10, family="Helvetica"))
  
  hw.dist.plot.chr <- ggplot(hw.dat, aes(x=-log10(P)))+
    stat_density(aes(color= CHR), position= "identity", geom= "line")+
    scale_color_manual(labels=names(chr.labels), values=rainbow(23))+
    ggtitle("HW pVal distribution  \n per chromosome")+
    xlab("-log10(HW-P)")+
    xlim(c(0, 20))+
    theme_bw()+
    theme(text=element_text(size=10, family="Helvetica"))
  
  hw.maf.dist.plot.chr <- ggplot(hw.dat[which(hw.dat$MAF > 0.01),], aes(x=-log10(P)))+
    stat_density(aes(color= CHR), position= "identity", geom= "line")+
    geom_vline(xintercept = 6)+
    ggtitle("HW pVal distribution of \n SNPs with a MAF > 0.01")+
    scale_color_manual(labels=names(chr.labels), values=rainbow(23))+
    xlim(c(0, 20))+
    xlab("-log10(HW-P)")+
    theme_bw()+
    theme(text=element_text(size=10, family="Helvetica"))
  
  hw.maf.dist.plot.all <- ggplot(hw.dat[which(hw.dat$MAF > 0.01),], aes(x=-log10(P)))+
    stat_density(position= "identity", geom= "line")+
    geom_vline(xintercept = 6)+
    xlim(c(0, 20))+
    ggtitle("HW pVal distribution of \n SNPs with a MAF > 0.01")+
    xlab("-log10(HW-P)")+
    theme_bw()+
    theme(text=element_text(size=10, family="Helvetica"))
  
  
  # layout matrix for grid.arrange
  lay <- rbind(c(1,1,1,1,2,2),
               c(3,3,4,4,5,5)
  )
  
  ### printing out the maf.hw distributions
  maf.hw.title <- paste0("MAF and HW distributions from Founders","\n", name," ", date())
  maf.hw.conclude<-paste0(sum(hw.dat$P<0.0001)," markers excluded ", (nrow(hw.dat)/3), " remaining.")
  maf.hw.tiff.file <- file.path(out, paste0(name, "MAFHW.density.tiff"))
  tiff(maf.hw.tiff.file, width = 3000, height = 2100, units = "px", res = 300, compression = "lzw")
  grid.arrange(maf.dist.plot.chr, maf.dist.plot.all, 
               hw.dist.plot.chr, hw.maf.dist.plot.chr, 
               hw.maf.dist.plot.all, 
               layout_matrix = lay, top=maf.hw.title, 
               bottom=textGrob(maf.hw.conclude, gp=gpar(fontsize=9,font=8)))
  dev.off()
}


#########################################################################################################
option_list = list(
  make_option(c("-p", "--plink"), type="character", default=NULL, 
              help="Path to plink files index, it assumes a bed, bim and fam file with the same file name", metavar="character"),
  make_option(c("-x", "--xplink"), type="character", default=NULL, 
              help="Path to plink files index, it assumes a bed, bim and fam file with the same file name", metavar="character"),
  make_option(c("-e", "--plinkexe"), metavar="character", default="PLINK/1.9-beta6-20190617", 
              help="complete name of the plink module to load" ),
  make_option(c("-f", "--pedfamfile"), metavar="character", default=NULL, 
              help=" .fam file with complete pedigree information with parents and sex in. [ORDER MUST BE THE SAME OF THE CURRENT .fam FILES] " ),
  make_option(c("-o", "--out"), metavar="character", default="./famCheck_genotypeQC", 
              help="Output path to save report")
); 

opt_parser  <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$pedfamfile)){
  print_help(opt_parser)
  stop("a pedigree file to with parents and sex info must be supplied to replace the .fam file", call.=FALSE)
}

out.founder <- file.path(opt$out, "7_FounderStats")
out.plot <- file.path(opt$out, "plots")
dir.create(out.founder)

##### Run plink command to merge all chromosomes

all.chr.list <- list.files(opt$plink, pattern = ".bim", full.names = T)
all.chr.list <- gsub(all.chr.list, pattern = ".bim$", replacement = "")

all.chr.Xchr.list <- c(all.chr.list, opt$xplink)
all.chr.Xchr.list <- gsub(all.chr.Xchr.list, pattern = ".bim$", replacement = "")
all.chr.Xchr.list.file <- file.path(out.founder, "plink.all.chr.Xchr.list")
write.table(as.data.frame(all.chr.Xchr.list), file= all.chr.Xchr.list.file, quote = FALSE, row.names = FALSE, col.names = FALSE)

###merge all chromosomes

merged.Xchr.plink <- file.path(out.founder, "merged.Xchr")
plink.merge.chr.Xchr.call <- paste0(
  "ml ", opt$plinkexe,"\n",
  " plink --merge-list ", all.chr.Xchr.list.file, " --out ", merged.Xchr.plink, " \n"
) 
system(plink.merge.chr.Xchr.call)

##### Extract only females for X chr founder analysis. 
##### With phenotypes from the pairing file generate a plink file for the X chromosome that would include on females. 
females.in.fam.file <- file.path(out.founder, "females.in.fam.txt")
new.fam<-fread(paste0(out.founder,"/merged.Xchr.fam"),data.table=F,header=F, colClasses = c(rep("character", 6)) )
names(new.fam)<-c("FID","IID","FATHER","MOTHER","Sex","PHENO")

females.in.fam <- new.fam[which(new.fam$Sex == 2), c("FID","IID")]
write.table(females.in.fam, file=females.in.fam.file, row.names = FALSE, quote = FALSE, col.names = FALSE, sep="\t")
xchr.founders.plink.file <- file.path(out.founder, "xchr.founders")
xchr.founder.plink.call <- paste0(
  "ml ", opt$plinkexe,"\n",
  "plink --bfile ", merged.Xchr.plink, " \\",
  "--fam ",opt$pedfamfile, " \\",
  "--chr X \\",
  "--make-bed \\",
  "--keep ", females.in.fam.file, " \\",
  "--no-pheno"," \\",
  "--out ", xchr.founders.plink.file, " \n"
)
system(xchr.founder.plink.call)

############################################################
############################################################
##### Run plink command to get stats within founders using 

## autosomes chromosomes only
plink.founderStats.call <- paste0(
  "ml ", opt$plinkexe,"\n",
  "plink ", 
  "--bfile ", merged.Xchr.plink,  " ", 
  "--fam ",opt$pedfamfile, " \\",
  "--freq ", 
  "--hardy ",
  "--filter-founders ",
  "--out ", out.founder, "/FounderStats \n"
)

system(plink.founderStats.call)


plink.founderStats.Xchr.call <- paste0(
  "ml ", opt$plinkexe,"\n",
  "plink ", 
  "--bfile ", xchr.founders.plink.file,  " ", 
  "--fam ",opt$pedfamfile, " \\",
  "--freq ", 
  "--hardy ",
  "--filter-founders ",
  "--out ", out.founder, "/Xchr.FounderStats \n"
)

system(plink.founderStats.Xchr.call)

## x chromosome only, including only females. 

maf.file <- file.path(out.founder, "FounderStats.frq")
maf.dat <- fread(maf.file, data.table = FALSE)
maf.dat$CHR <-as.factor(maf.dat$CHR)

hw.file <- file.path(out.founder, "FounderStats.hwe")
hw.dat <- fread(hw.file, data.table = FALSE)
hw.dat$CHR <-as.factor(hw.dat$CHR)
hw.dat$MAF <- maf.dat$MAF[match(hw.dat$SNP,maf.dat$SNP)]
chr.labels<-unique(hw.dat$CHR)

xchr.maf.file <- file.path(out.founder, "Xchr.FounderStats.frq")
xchr.maf.dat <- fread(xchr.maf.file, data.table = FALSE)
xchr.maf.dat$CHR <-as.factor(xchr.maf.dat$CHR)

xchr.hw.file <- file.path(out.founder, "Xchr.FounderStats.hwe")
xchr.hw.dat <- fread(xchr.hw.file, data.table = FALSE)
xchr.hw.dat$CHR <-as.factor(xchr.hw.dat$CHR)
xchr.hw.dat$MAF <- xchr.maf.dat$MAF[match(xchr.hw.dat$SNP,xchr.maf.dat$SNP)]


HW.MAF.dist.plot(maf.dat= maf.dat, hw.dat=hw.dat, out=out.plot, name="Autosome_founders")

xchr.maf.dist.plot.all <- ggplot(xchr.maf.dat, aes(x=MAF))+
  stat_density(color="black",  position= "identity", geom= "line")+
  ggtitle("MAF distribution")+
  theme_bw()+
  theme(text=element_text(size=10, family="Helvetica"))

xchr.hw.maf.dist.plot.all <- ggplot(xchr.hw.dat[which(xchr.hw.dat$MAF > 0.01),], aes(x=-log10(P)))+
  stat_density(position= "identity", geom= "line")+
  geom_vline(xintercept = 6)+
  xlim(c(0, 20))+
  ggtitle("HW pVal distribution of \n SNPs with a MAF > 0.01")+
  xlab("-log10(HW-P)")+
  theme_bw()+
  theme(text=element_text(size=10, family="Helvetica"))

### printing out the maf.hw distributions
xchr.maf.hw.title <- paste0("MAF and HW distributions from Founders","\n", "X_chr"," ", date())
xchr.maf.hw.conclude <-paste0(sum(xchr.hw.dat$P<0.0001)," markers excluded ", (nrow(xchr.hw.dat)), " remaining.")
xchr.maf.hw.tiff.file <- file.path(out.plot, paste0("X_chr", "MAFHW.density.tiff"))
tiff(xchr.maf.hw.tiff.file, width = 2000, height = 1500, units = "px", res = 300, compression = "lzw")
grid.arrange(xchr.maf.dist.plot.all, xchr.hw.maf.dist.plot.all, 
             top=xchr.maf.hw.title, nrow=1,
             bottom=textGrob(xchr.maf.hw.conclude, gp=gpar(fontsize=9,font=8)))
dev.off()

final.list.excluded.snps <- unique(c(xchr.hw.dat$SNP[which(xchr.hw.dat$P<0.0001)],
                                     hw.dat$SNP[which(hw.dat$P<0.0001)],
                                     xchr.hw.dat$SNP[which(xchr.hw.dat$P<0.0001)],
                                     xchr.maf.dat$SNP[which(xchr.maf.dat$MAF==0.000)]))

final.list.excluded.snps.file <- file.path(out.founder, "founder.final.list.excluded.snps")
write.table(final.list.excluded.snps, file=final.list.excluded.snps.file,col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

###### done ########


