#!/bin/bash
#

# take different actions depending on input arguments
case $# in
   0) echo -e "\nCorrect usage:\n\t $0 [chr X .vcf.gz filepath] [chr X .fam file]\n"
      exit 0
      ;;
   2) vcfin=$1
      famin=$2
      ;;
   *) echo -e "\nERROR :: Two arguments expected but $# given."
      echo -e "\nCorrect usage:\n\t $0 [chr X .vcf.gz filepath] [chr X .fam file]\n"
      exit 1
      ;;
esac

# check input VCF exists and isn't empty
if [ ! -s $vcfin ]
then
   echo -e "\nERROR :: Chromosome X VCF file \"$vcfin\" does not exist or is empty\n"
   exit 1
fi

# check input VCF exists and isn't empty
if [ ! -s $famin ]
then
   echo -e "\nERROR :: Chromosome X .fam file \"$famin\" does not exist or is empty\n"
   exit 1
fi

# check whether input file is bgzipped or not and use the appropriate cat command
if ( file `readlink -f $vcfin` | grep -q compressed )
then
   catcmd="zcat"
else
   catcmd="cat"
fi

# check whether fam file is of correct format
nnf=`awk '{print NF}' $famin | sort | uniq | wc -l`
if (( nnf == 1 ))
then
   nf=`awk '{print NF}' $famin | sort | uniq`
   if (( nf != 6 ))
   then
      echo -e "\nERROR :: Chromosome X .fam file \"$famin\" does not contain six columns."
      exit 1
   fi
else
   echo -e "\nERROR :: Chromosome X .fam file \"$famin\" line lengths are not equal."
   echo -e "           Please check your fam file and rerun.\n"
   exit 1
fi

# check whether fam file contains sex information (column 5)
nmal=`awk '$5==1' $famin | wc -l`
nfem=`awk '$5==2' $famin | wc -l`

if (( nmal == 0 )) && (( nfem == 0 ))
then
   echo -e "\nERROR :: Chromosome X .fam file \"$famin\" does not contain sex information."
   echo -e "           Please check code males as 1 and females as 2 in column 5.\n"
   exit 1
fi
if (( nmal < 10 ))
then
   echo -e "\nWARNING :: Chromosome X .fam file \"$famin\" contains <10 male samples."
   echo -e "             If this is correct, please ignore this warning."
fi
if (( nfem < 10 ))
then
   echo -e "\nWARNING :: Chromosome X .fam file \"$famin\" contains <10 female samples."
   echo -e "             If this is correct, please ignore this warning."
fi

# check whether IDs in .fam file match those in vcf
nmat=`awk 'NR==FNR{a[$1]++;next}($1 in a){print}' $famin <($catcmd $vcfin | awk '$0~/^#CHROM/{for(i=10;i<=NF;i++){print $i};exit}') | wc -l`
nall=`$catcmd $vcfin | awk '$0~/^#CHROM/{for(i=10;i<=NF;i++){print $i};exit}' | wc -l`

if (( nmat == 0 ))
then
   echo -e "\nERROR :: No match between sample IDs in input .fam file and in input"
   echo -e "           VCF file. Please check both files and, if necessary, edit"
   echo -e "           your .fam file to match the VCF sample IDs."
   exit 1
elif (( `echo "$nmat/$nall<0.9" | bc` == 1 ))
then
   echo -e "\nWARNING :: Less than 90% match between .fam sample IDs and VCF sample"
   echo -e "             IDs. Does the .fam file correspond to the same set of samples"
   echo -e "             as the VCF file?"
fi

### EXTRACT VARIANTS
# extract variants of interest - output file is "chrX_3vars.vcf"
echo "INFO :: Starting variant extraction from X chromosome VCF file"
awk 'NR==FNR{if(NR==1){next};split($2,y,":");split($3,z,":");a["X:"y[2]":"y[3]":"y[4]]++;a["X:"y[2]":"y[4]":"y[3]]++;a["chrX:"y[2]":"y[3]":"y[4]]++;a["chrX:"y[2]":"y[4]":"y[3]]++;a["23:"y[2]":"y[3]":"y[4]]++;a["23:"y[2]":"y[4]":"y[3]]++;a["chr23:"y[2]":"y[3]":"y[4]]++;a["chr23:"y[2]":"y[4]":"y[3]]++;a["X:"z[2]":"z[3]":"z[4]]++;a["X:"z[2]":"z[4]":"z[3]]++;a["chrX:"z[2]":"z[3]":"z[4]]++;a["chrX:"z[2]":"z[4]":"z[3]]++;a["23:"z[2]":"z[3]":"z[4]]++;a["23:"z[2]":"z[4]":"z[3]]++;a["chr23:"z[2]":"z[3]":"z[4]]++;a["chr23:"z[2]":"z[4]":"z[3]]++;next}$0~/^#/{print > "chrX_3vars.vcf";next}($3 in a){print > "chrX_3vars.vcf"; print "Variant "$3" extracted"}' check_variants.txt <($catcmd $vcfin)

# check whether any variants were found
nvars=`awk '$0!~/^#/' chrX_3vars.vcf | wc -l`
case $nvars in
   0) echo -e "\nERROR :: No variants extracted using either b37 or b38 positions."
      echo -e "           Please check the input VCF file and ensure that the input"
      echo -e "           variants listed in \"check_variants.txt\" are present in"
      echo -e "           the VCF file. If using your own variants to the file"
      echo -e "           \"check_variants.txt\", please ensure that their MAF is"
      echo -e "           at least 10% in your population and that the b37 and b38"
      echo -e "           positions and alleles are correct."
      exit 1
      ;;
   1) echo -e "\nWARNING :: Only one variant listed in \"check_variants.txt\" was"
      echo -e "             found in your X chromosome VCF. We recommend testing"
      echo -e "             the pliody of at least three variants. You can edit the"
      echo -e "             file \"check_variants.txt\" to include your own variants"
      echo -e "             if needed."
      ;;
   2) echo -e "\nWARNING :: Only two variant listed in \"check_variants.txt\" was"
      echo -e "             found in your X chromosome VCF. We recommend testing"
      echo -e "             the pliody of at least three variants. You can edit the"
      echo -e "             file \"check_variants.txt\" to include your own variants"
      echo -e "             if needed."
      ;;
   *) echo -e "\nINFO :: At least three variants listed in \"check_variants.txt\""
      echo -e "          were found in your X chromosome VCF file."
      ;;
esac

### Check whether variants are b37 or b38
genbld=`awk 'BEGIN{b37c=0;b38c=0}NR==FNR{if(NR==1){next};split($2,y,":");split($3,z,":");a["X:"y[2]":"y[3]":"y[4]]++;a["X:"y[2]":"y[4]":"y[3]]++;a["chrX:"y[2]":"y[3]":"y[4]]++;a["chrX:"y[2]":"y[4]":"y[3]]++;a["23:"y[2]":"y[3]":"y[4]]++;a["23:"y[2]":"y[4]":"y[3]]++;a["chr23:"y[2]":"y[3]":"y[4]]++;a["chr23:"y[2]":"y[4]":"y[3]]++;b["X:"z[2]":"z[3]":"z[4]]++;b["X:"z[2]":"z[4]":"z[3]]++;b["chrX:"z[2]":"z[3]":"z[4]]++;b["chrX:"z[2]":"z[4]":"z[3]]++;b["23:"z[2]":"z[3]":"z[4]]++;b["23:"z[2]":"z[4]":"z[3]]++;b["chr23:"z[2]":"z[3]":"z[4]]++;b["chr23:"z[2]":"z[4]":"z[3]]++;next}$0~/^#/{next}($3 in a){b37c++}($3 in b){b38c++}END{if(b37c>b38c){print "hg19"}else if(b38c>b37c){print "hg38"}}' check_variants.txt chrX_3vars.vcf`

### Call R script to run checks
Rscript plot_vcf_genotypes_by_sex.R chrX_3vars.vcf $famin $genbld
