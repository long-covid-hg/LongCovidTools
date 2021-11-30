#!/bin/bash
#

# get script name
scriptname="$(basename "$(test -L "$0" && readlink "$0" || echo "$0")")"

# store number of input args
nargs=$#

# create usage function
usage() {
   echo -e "\nINFO :: Usage:\n\n$scriptname -b [genome build] -i [plink binary fileset prefix] -p [plink executable path] -r [R or Rscript executable path]\nFlags \"-b\" and \"-i\" and their arguments are required. Flags \"-p\" and \"-r\" are recommended but optional.\n\nAccepted formats for the \"-b\" argument are \"GRCh37\", \"hg19\", \"b37\", \"GRCh38\", \"hg38\" or \"b38\".\n\nFor the \"-i\" flag, you must point to the prefix of your pre-imputed genotype files that are in plink .bed/.bim/.fam format. All three files must have the same prefix and if not stored in folder where this script is run, the prefix must include the full path of the fileset.\n\nIf the script is unable to locate your default plink installation, please specific the full path of the plink executable using the \"-p\" flag.\n\nThe script will try to detect your default installation of R (or Rscript). If this R installation does not have the required packages installed (\"ggplot2\"), or if no default installation is located, please specify the path of either the R or Rscript executable using the \"-r\" flag and ensure that this version has the required packages installed.\n" 1>&2; exit 1
}

#######################
###                 ###
### CHECK ARGUMENTS ###
###                 ###
#######################

# get arguments
while getopts “:b::i::p::r:” opt; do
  case $opt in
    b) if [ -z "$OPTARG" -o "${OPTARG:0:1}" = "-" ]
       then
	  echo -e "\nERROR :: Flag \"-b\" requires an argument.\n" 1>&2
	  exit 1
       else
          BUILD=$OPTARG
       fi ;;
    i) if [ -z "$OPTARG" -o "${OPTARG:0:1}" = "-" ]
       then
          echo -e "\nERROR :: Flag \"-i\" requires an argument.\n" 1>&2
          exit 1
       else
          INPUT=$OPTARG
       fi ;;
    p) if [ -z "$OPTARG" -o "${OPTARG:0:1}" = "-" ]
       then
          echo -e "\nERROR :: Flag \"-p\" requires an argument.\n" 1>&2
          exit 1
       else
          PLINK=$OPTARG
       fi ;;
    r) if [ -z "$OPTARG" -o "${OPTARG:0:1}" = "-" ]
       then
          echo -e "\nERROR :: Flag \"-r\" requires an argument.\n" 1>&2
          exit 1
       else
          RSCRIPT=$OPTARG
       fi ;;
    *) echo -e "\nERROR :: Unrecognised flag. Valid flags are \"-b\", \"-i\"," 1>&2
       echo -e "         \"-p\" and \"-r\". All flags must be given an argument.\n"
       exit 1 ;;
    :) echo -e "\nERROR :: Unrecognised flag used or argument missing.\n" 1>&2
       exit 1 ;;
  esac
done

# if no arguments included, print usage
(( nargs == 0 )) && usage

####################
###              ###
### CHECK INPUTS ###
###              ###
####################

# set start time to report duration at the end
START=$(date +%s)

### check whether input files exist and have size > 0
for ext in bed bim fam
do
   if [ ! -s "$INPUT.$ext" ]
   then
      echo -e "\nERROR :: Input file \"$INPUT.$ext\" not found or is empty.\n"
      exit 1
   fi
done

### check if plink is callable if defined (try to search for it, if not)
if [ -z "$PLINK" ]
then
   if ! command -v plink &> /dev/null
   then
       echo -e "\nERROR :: PLINK executable could not be detected automatically."
       echo -e "         Please specify PLINK's executable location using the \"-p\" flag.\n"
       exit 1
   else
      pcmd=`command -v plink`
      if [[ "$pcmd" == alias* ]] # if alias, extract command and remove first/last single and double quotes
      then
         plink=`echo $pcmd | awk '{split($0,a,"plink=");print a[2]}'`
         plink="${plink%\"}"; plink="${plink%\'}"; plink="${plink#\"}"; plink="${plink#\'}"
      else # otherwise, detected plink is an executable
         plink=$pcmd
      fi
   fi
else # test plink executable provided to -p flag
   $PLINK --help > /dev/null 2>&1 || { echo -e "\nERROR :: Specified PLINK executable not found or failed to run.\n         Please check the location given as an argument to the \"-p\" flag.\n"; exit 1; }
   plink=$PLINK
fi

# check version of PLINK
pversion=`$plink | head -n 1 | awk '{str=$2;gsub("v","",str);split(str,a,"[a-z]");print a[1]}'`
if ! echo $pversion | awk '$1<2{exit 1}'
then
   echo -e "\nERROR :: PLINK executable needs to be v2.00 or later - current version identified"
   echo -e "         as \"$pversion\". Please pass the full path of the PLINK2 executable to the \"-p\""
   echo -e "         flag. You can obtain a copy at https://www.cog-genomics.org/plink/2.0/\n"
   exit 1
fi

### check if either R or Rscript are callable
if [ -z "$RSCRIPT" ]
then
   if ! command -v Rscript &> /dev/null
   then
       if ! command -v R &> /dev/null
       then
          echo -e "\nERROR :: An R or Rscript executable could not be detected automatically.\n"
          echo -e "         Please specify the R or Rscript executable location using the \"-r\" flag."
          exit 1
       else
          pcmd=`command -v R`
          if [[ "$pcmd" == alias* ]] # if alias, extract command and remove first/last single and double quotes
          then
             Rscript=`echo $pcmd | awk '{split($0,a,"R=");print a[2]}'`
             Rscript="${Rscript%\"}"; Rscript="${Rscript%\'}"; Rscript="${Rscript#\"}"; Rscript="${Rscript#\'}"
          else # otherwise, detected plink is an executable
             Rscript=$pcmd
          fi
       fi
   else
      pcmd=`command -v Rscript`
      if [[ "$pcmd" == alias* ]] # if alias, extract command and remove first/last single and double quotes
      then
         Rscript=`echo $pcmd | awk '{split($0,a,"Rscript=");print a[2]}'`
         Rscript="${Rscript%\"}"; Rscript="${Rscript%\'}"; Rscript="${Rscript#\"}"; Rscript="${Rscript#\'}"
      else # otherwise, detected plink is an executable
         Rscript=$pcmd
      fi
   fi
else # test plink executable provided to -p flag
   $RSCRIPT --help > /dev/null 2>&1 || { echo -e "\nERROR :: Specified R or Rscript executable not found or failed to run.\n         Please check the location given as an argument to the \"-r\" flag.\n"; exit 1; }
   Rscript=$RSCRIPT
fi

### Test whether identified executable is R or Rscript and set the Rcmd variable according
if echo -e "1" | $Rscript -s -e 'x<-scan(file("stdin")); x' &> /dev/null
then
   Rcmd="$Rscript CMD BATCH" # command is just "R"
elif $Rscript -e "print(1)" &> /dev/null
then
   Rcmd="$Rscript" # command is "Rscript"
else
   echo -e "\nERROR :: Unable to determine whether the R command is \"R\" or \"Rscript\"."
   echo -e "         Please pass the full location of your R or Rscript executable as an"
   echo -e "         argument to the \"-r\" flag. If you have already done this, please"
   echo -e "         ensure that your R installation is functional and is v3.5.0 or newer.\n"
   exit 1
fi

#echo -e "\nDEBUG :: R command is: $Rcmd\n"

### check if R has required packages
if ! $Rscript -e "library(ggplot2)" &> /dev/null
then
   echo -e "\nERROR :: Detected or specified R or Rscript installation does not have the required"
   echo -e "         package \"ggplot2\" installed. Please either install \"ggplot2\" to your"
   echo -e "         default R or Rscript installation, or specify the location of an alternative"
   echo -e "         R or Rscript executable, using the \"-r\" flag, that has \"ggplot2\" installed.\n"
   exit 1
fi

### Check genome build flag
if [[ "$BUILD" =~ ^(GRCh37|hg19|b37)$ ]]
then
   genomebuild="b37"
elif [[ "$BUILD" =~ ^(GRCh38|hg38|b38)$ ]]
then
   genomebuild="b38"
else
   echo -e "\nERROR :: Unrecognised genome build specified. Please pass one of the following to"
   echo -e "         the \"-b\" flag: \"GRCh37\", \"hg19\", \"b37\", \"GRCh38\", \"hg38\" or \"b38\".\n"
   exit 1
fi

### Check plink input files by silently generating frequency information
### (if successful, keep counts data for PC projection
$plink --bfile $INPUT --freq counts --silent --out $INPUT > /dev/null 2>&1 || { perror=$( ( $plink --bfile $INPUT --freq --silent) 2>&1); echo -e "\nERROR :: PLINK encountered an error when checking the input files."; echo -e "         PLINK's error message is:\n"; echo "         =========================================="; echo "$perror" | awk '{print "         "$0}'; echo "         =========================================="; echo -e "\n         Please check your input files and the argument provided to the \"-i\" flag.\n"; exit 1; }

### Check other required files are present
for file in 1KG.b37.afreq.gz 1KG.b38.afreq.gz 1KG.eigenvec.gz 1KG.eigenvec.var.gz
do
   if [ -z "$file" ]
   then
      echo "\nERROR :: File \"$file\" is missing from this directory. Please ensure that the following"
      echo "         files are all present in this directory:\"1KG.eigenvec.gz\", \"1KG.eigenvec.var.gz\","
      echo "         \"1KG.b37.afreq.gz\" and \"1KG.b38.afreq.gz\".\n"
      exit 1
   fi
done


#####################
###               ###
### CALCULATE PCS ###
###               ###
#####################

# extract only unique biallelic SNPs from input dataset, by dropping sites with multiple variants
$plink --bfile $INPUT --extract <(awk 'BEGIN{a["A"]++;a["C"]++;a["G"]++;a["T"]++}NR==FNR{d[$1]++;next}(!($1":"$4 in d))&&($5 in a)&&($6 in a){print $2}' <(awk 'BEGIN{print "DUMMY"}a[$1":"$4]++{print $1":"$4}' $INPUT.bim) $INPUT.bim) --silent --make-bed --out $INPUT.filtered

# create new .bim file to match 1KG variant weights file
awk 'function nn(x){split(x,xx,"");str=n[xx[1]];for(i=2;i<=length(xx);i++){str=str""n[xx[i]]};return str}BEGIN{OFS="\t";n["A"]=1;n["C"]=2;n["G"]=3;n["T"]=4}{id=$1"_"$4"_"$5"_"$6}nn($6)<nn($5){id=$1"_"$4"_"$6"_"$5}{$2=id;print}' $INPUT.filtered.bim > $INPUT.filtered.CPRA.bim

# set correct ID column depending on genome build
if [[ "$genomebuild" == "b37" ]]
then
   idcol=1
else
   idcol=2
fi

# set max number of PCs to calculate
maxpcs=20

# run PC projection for input data (PLINK 2)
echo -e "\nINFO :: Projecting 1000G PCs into input data."
gunzip -c 1KG.eigenvec.var.gz > 1KG.eigenvec.var
gunzip -c 1KG.$genomebuild.afreq.gz > 1KG.afreq
maxpccol=$((maxpcs+4))
$plink --bed $INPUT.filtered.bed --bim $INPUT.filtered.CPRA.bim --fam $INPUT.filtered.fam --read-freq 1KG.afreq --silent --score 1KG.eigenvec.var $idcol 3 header-read no-mean-imputation variance-normalize --score-col-nums 5-$maxpccol --out $INPUT.filtered > /dev/null 2>&1
#done
rm 1KG.eigenvec.var 1KG.afreq

# merge with 1KG PCs and assign population name "STUDY"
maxpccol=$((maxpcs+2))
{
zcat 1KG.eigenvec.gz | cut -f1-$maxpccol
cut -f1,2,5- $INPUT.filtered.sscore | awk 'BEGIN{OFS="\t"}NR>1{$1="STUDY";print}'
} | gzip --best > $INPUT.1KG.eigenvec.gz

# plot using R
echo -e "INFO :: Creating combined cohort/1000G PC plots."
if [[ "$Rcmd" == *BATCH ]]
then
   $Rcmd --no-save --no-restore "--args $INPUT.1KG.eigenvec.gz $maxpcs" plot_PCs_by_pop.R Rplot.out && rm Rplot.out
else
   $Rcmd plot_PCs_by_pop.R $INPUT.1KG.eigenvec.gz $maxpcs > /dev/null 2>&1
fi

# clean up
#rm $INPUT.filtered.* $INPUT.filtered.CPRA.bim

# print exit message
END=$(date +%s)
DIFF=$(( $END - $START ))
hours=$(($DIFF/3600))
minutes=$((($DIFF/60)%60))
seconds=$(($DIFF%60))
(( hours < 10 )) || hh=$hours && hh=`printf "%0*d\n" 2 $hours`
mm=`printf "%0*d\n" 2 $minutes`
ss=`printf "%0*d\n" 2 $seconds`

echo -e "INFO :: Genetic PC projection and plotting complete in $hh:$mm:$ss.\n"
