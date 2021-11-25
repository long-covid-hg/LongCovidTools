#!/bin/bash
#

# take different actions depending on input arguments
case $# in
   0) echo -e "\nCorrect usage:\n\t $0 [chr X .vcf.gz filepath] [chr X .fam file] [genome build: \"hg19\" or \"hg38\"]\n"
      exit 0
      ;;
   3) vcfin=$1
      famin=$2
      genbld=$3
      ;;
   *) echo -e "\nERROR :: Three arguments expected but $# given."
      echo -e "\nCorrect usage:\n\t $0 [chr X .vcf.gz filepath] [chr X .fam file] [genome build: \"hg19\" or \"hg38\"]\n"
      exit 1
      ;;
esac

# check input VCF exists and isn't empty
if [ ! -s $vcfin ]
then
   echo -e "\nERROR :: Chromosome X VCF file \"$vcfin\" does not exist or is empty\n"
   exit 1
fi

# check input .fam exists and isn't empty
if [ ! -s $famin ]
then
   echo -e "\nERROR :: Chromosome X .fam file \"$famin\" does not exist or is empty\n"
   exit 1
fi

# check genome build is valid - only accept hg19 or hg38
if [ "$genbld" == "hg19" ]
then
   nonparstartpos=2699521
   nonparendpos=154931043
elif [ "$genbld" == "hg38" ]
then
   nonparstartpos=2781480
   nonparendpos=155701382
else
   echo -e "\nERROR :: Invalid genome build given. Please specify either"
   echo -e "         \"hg19\" or \"hg38\".\n"
fi

# check whether input file is bgzipped or not and use the appropriate cat command
if ( file `readlink -f $vcfin` | grep -q compressed )
then
   catcmd="zcat"
else
   catcmd="cat"
fi

# finally check that this is chromosome 23
# look for chr codes X, 23, chrX and chr23
chrcheck=`$catcmd $vcfin | awk 'BEGIN{a["X"]++;a[23]++;a["chrX"]++;a["chr23"]++}$0!~/^#/&&(!($1 in a)){print $1}' | sort | uniq`
if [[ ! -z "${chrcheck// }" ]]
then
   echo -e "\nERROR :: Non X-chromosome code(s) identified in input VCF. These are:"
   for chrc in $chrcheck
   do
      echo -e "         --> $chrc"
   done
   echo -e "\n         Please check that your VCF contains only chromosome X variants and"
   echo -e "         that the chromosome codes are consistent. Accepted codes for this"
   echo -e "         chromosome are \"X\", \"23\", \"chrX\" and \"chr23\".\n"
   exit 1
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
elif (( `echo "$nmat/$nall<0.999999" | bc` == 1 ))
then
   echo -e "\nERROR :: Less than 100% match between .fam sample IDs and VCF sample"
   echo -e "           IDs. Does the .fam file correspond to the same set of samples"
   echo -e "           as the VCF file?"
fi

# build new VCF name
basename=`echo $vcfin | awk '{split($1,a,"/");print a[length(a)]}' | sed 's/\.gz//g;s/\.vcf//g'`
newname="${basename}_het.males.to.alt.hom.vcf.gz"

# get fields containing male genotypes (and necessary columns)
awk -v npsp=$nonparstartpos -v npep=$nonparendpos 'BEGIN{OFS="\t"}NR==FNR{if($5==1){a[$1]++};next}$0~/^##/{print;next}$0~/^#CHR/{print;for(i=10;i<=NF;i++){if($i in a){c[i]++}};next}($2<npsp||$2>npep){print;next}{for(i=10;i<=NF;i++){if((i in c)&&($i=="0/1"||$i=="1/0")){$i="1/1"}else if((i in c)&&($i=="0|1"||$i=="1|0")){$i="1|1"}};print}' $famin <($catcmd $vcfin) | gzip --best > $newname

# print message to show completion
echo -e "\nINFO :: Conversion completed. New genotypes written to:"
echo -e "        \"$newname\".\n"

### EOF ###
