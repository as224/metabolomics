# script to preprocess harmonized file.gz from gwas catalog to use with merger.pl script later
# e.g. http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007515/harmonised/29632382-GCST007515-EFO_0001360.h.tsv.gz
# call with "bash preprocess.sh pmid29875488_IL1_eur.h.tsv.gz"
# M00054.metal.pos.txt.gz
# adjust name, trait and sample size


file=${1?parameter missing - file name.}

name="pmid24816252_trp_eur"  # use pubmed id, trait and population as name
trait="trp" # "M00054"
# samplesize="7824"	# todo sample size von denen nehmen?

# remove unnecessary columns (adjust columns)
zcat < $file | cut -f1,2,3,4,8,9,10,16,17,18 > $name.tsv
#zcat < $file > $name.tsv

# reorder columns
awk -F"\t" '{print $1"\t"$9"\t"$10"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' $name.tsv > ${name}_fix.tsv && mv ${name}_fix.tsv $name.tsv


# (optional) reorder columns if pvalue & std.error are switched in previous step (correct order: first std.err then pvalue)
#awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$10"\t"$9}' $name > ${name}_fix && mv ${name}_fix $name

echo -e "\nColumns must match (if not -> adjust cut fields in previous step (and reorder columns))"
echo -e "hm_rsid\thm_chrom\thm_pos\thm_other_allele\thm_effect_allele\thm_beta\thm_odds_ratio\thm_effect_allele_frequency\tstandard_error\tp_value"
head -n1 $name.tsv
echo ""


# Alternative: filter pvalue with python script
# python3 filter_pval.py $name.tsv ${name}_pval_python

# pval cutoff 0.05
# awk -F"\t" 'NR==1; NR>1 && $10 !~ /NA/ && $10<=0.05 { print $0 }' $name.tsv > ${name}_pval # for the new files, we don't have a given odds_ratio -> one column less
awk -F"\t" 'NR==1; NR>1 && $9 !~ /NA/ && $9<=0.05 { print $0 }' $name.tsv > ${name}_pval # 173347

# create "header" file, reorder columns & attach to header
echo -e "SNP\tchr\tbp_38\teffect_allele\tnon_effect_allele\tEAF\teffect_direction\tOR\tSE\tpvalue\ttrait\tsample_size" > ${name}_pval_final
awk -F"\t" -v t=$trait -v s=$samplesize 'NR>1 {print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$8"\t"$6"\t"exp($6)"\t"$9"\t"$10"\t"t"\t"s}' ${name}_pval >> ${name}_pval_final

# create same file but without pvalue filtering if needed
# echo -e "SNP\tchr\tbp_38\teffect_allele\tnon_effect_allele\tEAF\teffect_direction\tOR\tSE\tpvalue\ttrait\tsample_size" > ${name}_no_pval
# awk -F"\t" -v t=$trait 'NR>1 {print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$8"\t"$6"\t"exp($6)"\t"$9"\t"$10"\t"t"\t"s}' ${name}.tsv >> ${name}_no_pval


# clean-up (optional)
# rm ${name}_pval
# mv ${name}_pval_final ${name}_afterSH.tsv
