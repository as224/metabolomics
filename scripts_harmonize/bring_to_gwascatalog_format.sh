
# commands to transform raw gwas file to harmonized format (after preprocessing.sh)
# run: bash bring.sh M32198.gz acetylcarnitine
# bash bring_to_gwascatalog_format.sh Shin_selected_2/M15122.metal.pos.txt.gz glycerol
# M15122	glycerol

gwas=$1 	# eg. M12345.gz
trait=$2	# eg. alanine

# format columns
awk -F"\t" -v t="$trait" '{print $1"\t"$17"\t"$18"\t"toupper($2)"\t"toupper($3)"\t"$4"\t"$8"\t"exp($8)"\t"$9"\t"$10"\t"t"\t"$16}' <(zcat < $gwas) > ${trait}_unfiltered.tsv # Alleles in upper case additionally

# py script if needed to filter unvalid EAFs
# python3 filter_EAF.py ${trait}_unfiltered.tsv ${trait}_EAF_unfiltered.tsv

# filter p>0.05 
awk -F"\t" 'NR==1; NR>1 && $10 !~ /NA/ && $10<=0.05{print $0}' ${trait}_unfiltered.tsv > ${trait}_pval.tsv # adjust file names if conducted EAF filtering

# create .bed format file of filtered and unfiltered data for liftover tool
awk -F"\t" 'NR>1{print "chr"$2"\t"$3"\t"$3"\t"$1}' ${trait}_pval.tsv > ${trait}_pval.bed 
awk -F"\t" 'NR>1{print "chr"$2"\t"$3"\t"$3"\t"$1}' ${trait}_unfiltered.tsv > ${trait}_unfiltered.bed

# upload .bed file to "https://genome.ucsc.edu/cgi-bin/hgLiftOver"
# select hg18 as original and hg38 as new assembly
# download results

# (edit and run manually on console) 
# merge result from liftover website with gwas-file to get new positions (manually edit the name of the website result file)
#awk -F"\t" 'FNR==NR{a[$4]=$2; next}{print $1"\t"$2"\t"a[$1]"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' hglft_glycerol.bed glycerol_pval.tsv > glycerol_pval.final 
#awk -F"\t" 'FNR==NR{a[$4]=$2; next}{print $1"\t"$2"\t"a[$1]"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' hglft_glycerol_unfiltered.bed glycerol_unfiltered.tsv > glycerol_unfiltered.final


# clean-up
# rm ${trait}_unfiltered.tsv
