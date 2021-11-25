# read input arguments: 1 = VCF file, 2 = fam file
args <- commandArgs(TRUE)
vcf_fil <- args[1]
fam_fil <- args[2]
gen_bld <- args[3]

# get list of installed packages
inst.pkgs <- installed.packages()[,"Package"]

# install BiocManager if not installed
if(!requireNamespace("BiocManager", quietly = TRUE)) {
   install.packages("BiocManager")
}

# install package through biocLite, if not installed
if(!"VariantAnnotation" %in% inst.pkgs) {
   BiocManager::install("VariantAnnotation")
}

if(!"ggplot2" %in% inst.pkgs) {
   install.packages("ggplot2")
}
if(!"gridExtra" %in% inst.pkgs) {
   install.packages("gridExtra")
}


# load VariantAnnotation library
suppressMessages(library("VariantAnnotation"))
suppressMessages(library("ggplot2"))
suppressMessages(library("gridExtra"))


### READ IN AND CHECK EXTRACTED VARIANTS
# read the vcf file
vcf_dat <- readVcf(vcf_fil,gen_bld)

# read the sample data
fam_dat <- read.table(fam_fil,sep=" ")[c("V1","V5")]
fam_dat$IID <- as.character(fam_dat$V1)
fam_dat$SEX_M1_F2 <- fam_dat$V5
fam_dat <- fam_dat[c("IID","SEX_M1_F2")]

# now check if VCF contains genotype dosages
if ("DS" %in% names(geno(vcf_dat))) {
   use="DS"
   vcf_ds <- geno(vcf_dat)$DS
   vcf_vars <- rownames(vcf_ds)
} else if("GP" %in% names(geno(vcf_dat))) {
   cat("\nWARNING :: Genotype format \"DS\" (dosages) not found in VCF file.\n")
   cat("           Genotype format \"GP\" (genotype probabilities) will be used.\n")
   use="GP"
   vcf_gp <- geno(vcf_dat)$GP
   vcf_vars <- rownames(vcf_gp)
} else {
   cat("\nERROR :: Genotype formats \"DS\" and \"GP\" not found in VFC file.\n")
   cat("         Unable to check X chromosome ploidy... exiting.\n")
   quit(status=1)
}

### loop over variants and plot probabilities
for(var in vcf_vars) {

   # get REF and ALT alleles
   var_ref <- as.character(rowRanges(vcf_dat)[match(var,names(rowRanges(vcf_dat)))]$REF)
   var_alt <- as.character(unlist(rowRanges(vcf_dat)[match(var,names(rowRanges(vcf_dat)))]$ALT))

   if (use=="DS") {

      # create data frame with this variant's dosages
      vcf_ds_var <- data.frame(IID=names(vcf_ds[dimnames(vcf_ds)[[1]]==var,]),DOSAGE=as.numeric(vcf_ds[dimnames(vcf_ds)[[1]]==var,]))
      # match in sex and separate into female/male
      vcf_ds_var$SEX_M1_F2 <- fam_dat$SEX_M1_F2[match(vcf_ds_var$IID,fam_dat$IID)]
      vcf_ds_var_fem <- subset(vcf_ds_var,SEX_M1_F2==2)
      vcf_ds_var_mal <- subset(vcf_ds_var,SEX_M1_F2==1)

      # create plot objects for each
      p_female <- ggplot(vcf_ds_var_fem,aes(x=DOSAGE)) + geom_histogram(binwidth=0.1) + labs(title="Female",x=paste(var,var_alt,"Dosage")) + coord_cartesian(xlim=c(0,2)) + theme_bw()
      p_male <- ggplot(vcf_ds_var_mal,aes(x=DOSAGE)) + geom_histogram(binwidth=0.1) + labs(title="Male",x=paste(var,var_alt,"Dosage")) + coord_cartesian(xlim=c(0,2)) + theme_bw()
      plot_ <- grid.arrange(p_female,p_male, ncol=2)

      # save plot
      ggsave(file=paste0("dosage_distribution_by_sex___",gsub(":","_",var),".png"),plot=plot_,device="png",width=8,height=4,units="in")

   } else if(use=="GP") {

      # create data frame with this variant's genotype probabilities
      vcf_gp_var <- data.frame(IID=names(vcf_gp[dimnames(vcf_gp)[[1]]==var,,1]),PROBRR=as.numeric(vcf_gp[dimnames(vcf_gp)[[1]]==var,,1]),PROBRA=as.numeric(vcf_gp[dimnames(vcf_gp)[[1]]==var,,2]),PROBAA=as.numeric(vcf_gp[dimnames(vcf_gp)[[1]]==var,,3]))
      # match in sex
      vcf_gp_var$SEX_M1_F2 <- fam_dat$SEX_M1_F2[match(vcf_gp_var$IID,fam_dat$IID)]
      # flatten the data
      vcf_gp_var_flat <- data.frame(IID=rep(vcf_gp_var$IID,3),SEX_M1_F2=rep(vcf_gp_var$SEX_M1_F2,3),PROB=c(vcf_gp_var$PROBRR,vcf_gp_var$PROBRA,vcf_gp_var$PROBAA),GENO=c(rep(c(paste0(var_ref,var_ref),paste0(var_ref,var_alt),paste0(var_alt,var_alt)),each=length(vcf_gp_var$IID))))
      # now separate into male and female
      vcf_gp_var_flat_fem <- subset(vcf_gp_var_flat,SEX_M1_F2==2)
      vcf_gp_var_flat_mal <- subset(vcf_gp_var_flat,SEX_M1_F2==1)

      # create plot objects
      p_female <- ggplot(vcf_gp_var_flat_fem,aes(x=PROB)) + geom_histogram(aes(fill=GENO),col="black",binwidth=0.05) + labs(title="Female",x="Genotype probability") + coord_cartesian(xlim=c(0,1)) + theme_bw()
      p_male <- ggplot(vcf_gp_var_flat_mal,aes(x=PROB)) + geom_histogram(aes(fill=GENO),col="black",binwidth=0.05) + labs(title="Male",x="Genotype probability") + coord_cartesian(xlim=c(0,1)) + theme_bw()
      plot_ <- grid.arrange(p_female,p_male, ncol=2)

      # save the plot to file
      ggsave(file=paste0("genotype_probabilities_by_sex___",gsub(":","_",var),".png"),plot=plot_,device="png",width=8,height=4,units="in")

   }

}



