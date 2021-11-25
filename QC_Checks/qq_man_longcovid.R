#!/usr/bin/env Rscript

### required libaries
#packs <- c("optparse","ggplot2","tidyverse","data.table","R.utils")
packs <- c("optparse","qqman","data.table","R.utils")

for(p in packs){
   if(!require(p,character.only=T,quietly=T)){
      install.packages(p,repos=c(CRAN="http://cran.r-project.org"))
   }
   library(p,character.only=T,quietly=T)
}


### set input arguments
option_list = list(
   make_option(c("-f", "--file"), type="character", default=NULL,
      help="dataset file name", metavar="character"),
   make_option(c("-o", "--out"), type="character",
      help="output file name [default= %default]", metavar="character")
);

### perform input check and set input file name
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, positional_arguments=0);

infile <- opt$options$file

# read first row of file to get header
infile_colnames <- colnames(fread(infile, header=T, nrows=0))

# check column names to find format
regenie_colnames <- c("CHROM","GENPOS","A1FREQ","INFO","LOG10P")
saige_colnames   <- c("CHR","POS","AF_Allele2","imputationInfo","p.value")

if(all(regenie_colnames %in% infile_colnames)) {
   filetype="REGENIE"
} else if (all(saige_colnames %in% infile_colnames)) {
   filetype="SAIGE"
} else {
   cat(paste0("ERROR :: Minimal set of columns not present in input file: ", infile,"\n"))
   cat("         Please ensure that at least the following columns are present:\n")
   cat("         -> REGENIE: CHROM GENPOS A1FREQ INFO and LOG10P\n")
   cat("         -> SAIGE: CHR POS AF_Allele2 imputationInfo p.value\n\n")
   quit(status=1)
}
# report filetype
cat(paste0("INFO :: Input file format is: ",filetype,"\n"))

# set output file prefix using option (or using filename if not set)
if(!is.null(opt$options$out)){
      output_prefix=opt$options$out
} else {
   output_prefix=gsub('\\.txt','',gsub('\\.gz','',tail(unlist(strsplit(infile,"/",fixed=T)),n=1)))
   cat(paste0("INFO :: Output prefix not set, defaulting to: ",output_prefix,"\n"))
}


### read in data, extracting only required fields and convert fields to standard form
newcolnames <- c("CHR","POS","EAF","INFO")

cat(paste0("INFO :: Reading input file: ",infile,"\n"))
if(filetype=="REGENIE") {
   data <- fread(infile,header=T,select=regenie_colnames)
   colnames(data)[which(colnames(data) %in% regenie_colnames[!grepl('LOG10P',regenie_colnames)])] <- newcolnames
   data$P <- 10^(-data$LOG10P)
   #data <- subset(data,select=-c(LOG10P))
} else {
   data <- fread(infile,header=T,select=saige_colnames)
   colnames(data)[which(colnames(data) %in% saige_colnames[!grepl('p.value',saige_colnames)])] <- newcolnames
   colnames(data)[colnames(data)=="p.value"] <- "P"
   data$LOG10P <- -log10(data$P)
}


### correct chromosome code to number and filter out any that can't be converted after substitutions
cat("INFO :: Converting chromosome codes to numeric format\n")
data$CHR <- gsub("chr","",data$CHR)
data$CHR <- gsub("X|chrX","23",data$CHR)
data$CHR <- gsub("Y|chrY","24",data$CHR)
data$CHR <- gsub("MT|chrMT|M|chrM","25",data$CHR)

data$CHR <- as.numeric(data$CHR)
data <- data[!is.na(data$CHR)]


### add dummy SNP column because qqman is broken and needs it, despite claiming not to...
### (bug may only exist in some versions of qqman)
data$SNP <- paste(data$CHR,data$POS,sep=":")


### filter input results based on EAF and INFO
data <- subset(data,!is.na(EAF)&EAF>0.001&EAF<0.999&!is.na(INFO)&INFO>0.8)

### plot QQ
quants <- c(0.7,0.5,0.1,0.01,0.001)
subdata <- data[!is.na(data$P)&is.numeric(data$P)]
lambda  <- round(quantile((qchisq(1-subdata$P,1)),probs=quants)/qchisq(quants,1),3)
cat("INFO :: Generating QQ plot\n")
png(paste0(output_prefix,"_P_qqplot.png"))
qq(subdata$P, main=paste("\nlambda ", quants, ": ", lambda, sep="" ))
dev.off()

sink(paste0(output_prefix,"_P_qquantiles.txt"))
cat(paste0(quants,": ",lambda,"\n"))
sink()


### plot standard manhattan
cat("INFO :: Subsetting P-values < 0.01 for manhattan plot\n")
subdata <- subdata[subdata$P<0.01&subdata$P>0]
cat(paste("INFO :: Plotting manhattan with",nrow(subdata),"variants\n"))
#print(summary(subdata$P))
png(paste0(output_prefix,"_P_manhattan.png"),width=1000,height=400)
manhattan(data.table(subdata[,c("SNP","POS","P","CHR"),with=F]),bp="POS",ylim=c(2,max(subdata$LOG10P)+1))
dev.off()



