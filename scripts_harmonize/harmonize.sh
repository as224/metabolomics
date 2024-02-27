
# commands to transform raw gwas file to harmonized format (after preprocessing.sh)
# run: bash harmonize.sh -i gwas -t trait --pvalue --bed --eaf 
# run: bash harmonize.sh -i pmid29875488_IL1_eur.merge.final.gz -t alanine --pvalue --bed --eaf 

# 1. format columns 
# 2. filter p>0.05 -> if pvalue (& remove unfiltered file)
# 3. if needed to filter unvalid EAFs
# 4. create .bed format file of filtered and unfiltered data for liftover tool

#!/bin/bash

# set default values
filter_pvalue=false # only filters for pvalue, if wanted
apply_bed_filter=false # only create bed, if wanted
filter_eaf=false # exclude columns where the EAF is NA 
gwas=""
trait=""

# Check the used parameters
while [[ "$#" -gt 0 ]]; do
  case $1 in
    --pvalue)
      filter_pvalue=true
      shift
      ;;
    --bed)
      filter_bed=true
      shift
      ;;
    --eaf)
      filter_eaf=true
      shift
      ;;
    -i|--input)
      gwas="$2"
      shift 2
      ;;
    -t|--trait)
      trait="$2"
      shift 2
      ;;
    *)
      echo "Unknown Parameter: $1"
      exit 1
      ;;
  esac
done

# Check if all required parameters are provided
if [ -z "$gwas" ] || [ -z "$trait" ]; then
  echo "Usage: $0 [--pvalue] [--bed] [--eaf] -i <gwas> -t <trait>"
  exit 1
fi

# format columns -> adapt that to the specific gwas file
awk -F"\t" -v t="$trait" '{print $1"\t"$17"\t"$18"\t"toupper($2)"\t"toupper($3)"\t"$4"\t"$8"\t"exp($8)"\t"$9"\t"$10"\t"t"\t"$16}' <(zcat < $gwas) > ${trait}_unfiltered.tsv # are Alleles in upper case

# py script if needed to filter unvalid EAFs
if [ "$filter_eaf" = true ]; then
  python3 filter_EAF.py "${trait}_unfiltered.tsv" "${trait}_EAF_unfiltered.tsv"
fi


# filter pvalues if needed
if [ "$filter_pvalue" = true ]; then
    awk -F"\t" 'NR==1; NR>1 && $10 !~ /NA/ && $10<=0.05{print $0}' ${trait}_unfiltered.tsv > ${trait}_pval.tsv # adjust file names if conducted EAF filtering
    if [ "$filter_bed" = true ]; then
        awk -F"\t" 'NR>1{print "chr"$2"\t"$3"\t"$3"\t"$1}' ${trait}_pval.tsv > ${trait}_pval.bed # create .bed format file for liftover tool
    fi
    rm ${trait}_unfiltered.tsv # clean-up
else
    if [ "$filter_bed" = true ]; then
        awk -F"\t" 'NR>1{print "chr"$2"\t"$3"\t"$3"\t"$1}' ${trait}_unfiltered.tsv > ${trait}_unfiltered.bed # create .bed format file for liftover tool
    fi
fi

# upload .bed file to "https://genome.ucsc.edu/cgi-bin/hgLiftOver"
# select hg18 as original and hg38 as new assembly
# download results

# (edit and run manually on console) 
# merge result from liftover website with gwas-file to get new positions (manually edit the name of the website result file)
#awk -F"\t" 'FNR==NR{a[$4]=$2; next}{print $1"\t"$2"\t"a[$1]"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' hglft_glycerol.bed glycerol_pval.tsv > glycerol_pval.final 
#awk -F"\t" 'FNR==NR{a[$4]=$2; next}{print $1"\t"$2"\t"a[$1]"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' hglft_glycerol_unfiltered.bed glycerol_unfiltered.tsv > glycerol_unfiltered.final