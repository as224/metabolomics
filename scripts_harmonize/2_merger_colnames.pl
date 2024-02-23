#! /usr/bin/perl -w

use strict;
use diagnostics;
use File::Basename;

# Script adds allele frequency from dbSNP to studyfiles and/or fixes MAF incase it was provided
# Usage on single file:
# perl ~/2_merger_colnames.pl /home/icb/studyfilefolder/pmid123.tsv /lustre/groups/sysmet01/projects/snipa/data/grch38/refset/1kgpp3v5/all/vep_new/external_files/dbsnp/ /lustre/groups/sysmet01/projects/studyfile_pipeline/condition_files_2/result_folder/ &
# perl 2_merger_colnames.pl  Shin_afterSH/tryptophan_unfiltered.final /home/ output
# Useage on whole folder (with nohup):
#nohup perl ~/2_merger_colnames.pl /home/icb/studyfilefolder/ /lustre/groups/sysmet01/projects/snipa/data/grch38/refset/1kgpp3v5/all/vep_new/external_files/dbsnp/ &

my $harmonized_format = 1;	# gwas catalog harmonized format (after preprocessing)?
my $add_allele_freq = 0;		# add allele frequency?
my $debug = 0;	# debug messages?


# check number of command line arguments passed
die("Did not provide enough (2 minimum) arguments to run\n") if $#ARGV<1;

# study file as input argument 1
my $sfile = $ARGV[0];
print "\nPath to folder where study file(s) are located: $sfile\n" if -d $sfile;
print "\nLocation of file to process: $sfile\n" if -f $sfile;

# dbSNP file as input argument 2
my $dbsnp = $ARGV[1]; # dbsnp file / chr
print "\nPath to dbSNP files split by chromosome: $dbsnp\n" if $add_allele_freq == 1;

# output file path as input argument 3
my $outfile_path = $ARGV[2];
$outfile_path = "${sfile}out/" if !defined($outfile_path) && -d $sfile;
#$outfile_path = join("/",dirname($sfile),"out/") if !defined($outfile_path) && -f $sfile;
die("You must provide an outfile path if you want to process a single file!\n") if !defined($outfile_path) && -f $sfile;
print "\nPath to output: $outfile_path\n";

# variable declarations
my ($ch, $p, $rs, $rsid, $ea, $oa, $eaf, $effect, $ss, $se, $pval, $gb, $trait, $effect_sign, $or, $def, $ndef, $whole_string, $whole_string_reverse_alleles);
my @info = ();
my %seen;

my @chromosomes_to_iterate = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23);#23=="X"

foreach my $chr_i (@chromosomes_to_iterate){

$chr_i="X" if $chr_i==23;

# LOAD dbSNP file into hash (memory)
print "\nLoading dbSNP file for chromosome $chr_i into memory...\n" if $add_allele_freq == 1;
my %dbsnp_hash = ();
#open(IN, "zcat /lustre/groups/sysmet01/projects/snipa/data/grch38/refset/1kgpp3v5/all/vep_new/external_files/dbsnp/chr${dbsnp}_snv_mnv_indel_freq_VCfix_dbSNP155.vcf.gz |") or die $!;
open(IN, "zcat ${dbsnp}chr$chr_i.dbSNP156.38.ensembl.gz |") or die $! if $add_allele_freq == 1;
my $line_limit=0; 

if($add_allele_freq == 1){

while(<IN>){
	chomp($_);
	my $line = $_;
	#	last if ++$line_limit==10000;	#limit entries for testing purpose
	
	@info = split("\t", $_);
	
	my $dbsnp_chr = $info[0];
	my $pos = $info[1];
	my $ra = $info[3];
	my $alt = $info[4];
	my @aa = split(",",$alt);
	my $rs = $info[2];
	next if !defined($rs);
	my @af_split; my $af_source;
	my $af = $info[5];

	my $allele_counter = 1;# keeps track of which how many alleles are being saved per entry
	#TODO change to use rsid instead of position
	#TODO maybe install a switch to use position if rsid is not available?
	
	# for every alternative allele save different hash entry
	foreach (@aa){
		$whole_string="$dbsnp_chr:$info[1]:$_:$ra";

		if(defined($rs) && defined($af)){#check if rsid and allele-freq are present in dbSNP file; if not -> put in "undef"
		@af_split = split(/\|/,$af);
		my @topmed_A = grep(/^TOPMED:/,@af_split);   # find TOPMED entry
		my @gnomAD_A = grep(/^GnomAD:/,@af_split);   # find GnomAD entry
		if(scalar @topmed_A > 0){
			$af_source = $topmed_A[0];# First prio -> TOPMED
		} elsif (scalar @gnomAD_A > 0){
			$af_source = $gnomAD_A[0];# Second prio -> gnomAD
		} else {
			$af_source = $af_split[0];# Otherwise -> take first
		}
		$dbsnp_hash{$whole_string} = "$rs;$af_source;$allele_counter";

		} elsif(!defined($rs) && defined($af)) {
			@af_split = split(/\|/,$af);
			my @topmed_A = grep(/^TOPMED:/,@af_split);   # is TOPMED entry there?
			my @gnomAD_A = grep(/^GnomAD:/,@af_split);   # is GnomAD entry there?
			if(scalar @topmed_A > 0){
				$af_source = $topmed_A[0];# First prio -> TOPMED
			} elsif (scalar @gnomAD_A > 0){
				$af_source = $gnomAD_A[0];# Second prio -> gnomAD
			} else {
				$af_source = $af_split[0];# Otherwise -> take first
			}
			$dbsnp_hash{$whole_string} = "undef;$af_source;$allele_counter";

			} elsif(defined($rs) && !defined($af)) {
				$dbsnp_hash{$whole_string} = "$rs;undef;$allele_counter";
			} elsif(!defined($rs) && !defined($af)) {
				$dbsnp_hash{$whole_string} = "undef;undef;$allele_counter";
			}
			$allele_counter++;
	}
}
close(IN);

#my $count = 1;
#foreach my $key (keys %dbsnp_hash) {
#print $key, " = \t", $dbsnp_hash{$key}, "\n";
#$count++;
#last if $count==10;
#}

}

my @files;
if( -f $sfile ){
	print "You provided a single file for processing.\n" if $debug==1;
	push(@files, $sfile);

} elsif( -d $sfile ){
	print "You provided a folder of files for processing.\n" if $debug==1;
	
	# read files in provided folder
	opendir(DH, $sfile);
	@files = readdir(DH);
	closedir(DH);
}


# Iterate over every provided file
foreach my $file (@files){

# skip . and ..
next if($file =~ /^\.$/);
next if($file =~ /^\.\.$/);
#next if($file =~ /^lof/);

if($debug==1){
	print "File in process: $file\n";
	print "Outfile here: $outfile_path/$file.merge\n" if -d $sfile;
	print "Logfile here: $outfile_path/$file.merge.log\n" if -d $sfile;
	print "Outfile here: $file.merge\n" if -f $sfile;
	print "Logfile here: $file.merge.log\n" if -f $sfile;
}

if(-f $sfile){
	# result file (careful: it will always append to the end of a file if a previous result version exists -> move/rm previous result files)
	open(OUT, ">>", "$file.merge") or die $!;
	# log file (same appending happening here)
	open(LOG, ">>", "$file.merge.log") or die $!;
} else {
	# result file (careful: it will always append to the end of a file if a previous result version exists -> move/rm previous result files)
	open(OUT, ">>", "$outfile_path/$file.merge") or die $!;
	# log file (same appending happening here)
	open(LOG, ">>", "$outfile_path/$file.merge.log") or die $!;
}

if($harmonized_format==1){
	if(-f $sfile){
		open(IN, "<", "$file") or die $!;
	} else {
		open(IN, "<", "${sfile}${file}") or die $!;
	}
} else {
	if(-f $sfile){
		open(IN, "zcat ${file} |") or die $!;
	} else {
		open(IN, "zcat ${sfile}${file} |") or die $!;
	}
}

# some counters for the log file
$def = 0;
$ndef = 0;
my $case1 = 0; my $case2 = 0; my $case3 = 0; my $case4 = 0; my $case5 = 0; my $reverse_match = 0;
my $case1_r = 0; my $case2_r = 0; my $case3_r = 0; my $case4_r = 0; my $case5_r = 0; my $case1_no = 0; my $case2_no = 0; my $case3_no = 0;

# header for result file
#print OUT "chr\tpos\trsid\tEAF\ta1\ta2\teffect\tsampleSize\tstdErr\tp-val\tgenomeBuild\ttrait\teffect_sign\tOR\tsource\n" if !defined($seen{$file});
print OUT "SNP\tchr\tbp_38\teffect_allele\tnon_effect_allele\tEAF\tbeta_effect\tOR\tSE\tpvalue\ttrait\tsample_size\tEAF_source\n" if !defined($seen{$file});
$seen{$file}=1;
print LOG "Merging $file chr$chr_i with dbSNP...\n" if $add_allele_freq == 1;
print LOG "Checking alleles in $file chr$chr_i ...\n" if $add_allele_freq == 0;


# Processing Study Files:**
# The script processes either a single study file or all files in a folder.
# It iterates over each chromosome, loads dbSNP information, and processes each study file.
# For each study file, it opens an output file and a log file.

while(<IN>){
	if ($. == 1) {
        next;  # Skip the first line (header)
    }
	chomp($_);
	my $line = $_;
	@info = split("\t", $_);

	my $eaf_source;
	if($harmonized_format==0){ # custom column format of input files
		$ch = $info[0];
		$p = $info[1];
		$ea = $info[3];
		$oa = $info[4];
		$effect = $info[5];
		$ss = $info[6];
		$se = $info[7];
		$pval= $info[8];
		$gb = $info[9];
		$trait = $info[10];
		$effect_sign = $info[11];
		$or = $info[12];
		$rsid = $info[13];
	} else {  # file format of harmonized gwas catalog files
		$rsid = $info[0];
		$ch = $info[1];
		$p = $info[2];
		$ea = $info[3];
		$oa = $info[4];
		$eaf = $info[5];
		$effect = $info[6];
		$or = $info[7];
		$se = $info[8];
		$pval= $info[9];
		$trait = $info[10]; 	# edit here or change later
		$ss = $info[11]; 		# edit here or change later
		$eaf_source = "Study"; 	# edit here or change later
		#$ss = "NA"; 		# edit here or change later
	}

	next if($ch eq "NA" || $ch eq "chr" || $ch eq "hm_chrom");
	# next if($pval>0.05);
	if ($pval !~ /^[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?$/) {
    print "Non-numeric P-value found: $pval\n"; #ADDED hier code hinzugefügt
	} #elsif ($pval > 0.05) {
    	#next; # Skip the current iteration if $pval is greater than 0.05
	#}

	$ch = 23 if $ch eq "X"; # dbsnp gives chrX as 23
	$chr_i = 23 if $chr_i eq "X"; #LN: changes also the chr_i back to numeric
	next if($ch != $chr_i);

	$whole_string = "$ch:$p:$ea:$oa"; # string to search dbSNP hash
	$whole_string_reverse_alleles = "$ch:$p:$oa:$ea"; # string to search dbSNP hash with switched alleles

# Processing Individual Study Entries:
# reads each line from the study file, extracts relevant information, and checks for matches in the dbSNP hash.
# handles different cases based on whether dbSNP information is found.
# calculates new effect sizes and allele frequencies if needed.


 if($add_allele_freq == 1){ # -> add allele freq? if yes -> search hash for match
	if(defined($dbsnp_hash{$whole_string})){# -> study entry found in dbSNP

		# get contents of dbSNP hash saved earlier to extract source/MAF
		# dbsnp hash content: 1:12345:A:T => rs123;1000Genomes:0.9208,0.07917|GnomAD:0.9939,0.006087;1
		my $hash_content = $dbsnp_hash{$whole_string};

		my @hash_content_split = split(";", $hash_content);
		my $dbsnp_rs = $hash_content_split[0];
		my $dbsnp_af = $hash_content_split[1];
		my $dbsnp_ac = $hash_content_split[2];

		my $source; my $freq; my $dbsnp_ra_freq; my $dbsnp_aa_freq;
		if($dbsnp_af ne "undef"){
			my @af_split = split(":", $dbsnp_af);
			$source = $af_split[0];
			$freq = $af_split[1];
			my @freq_split = split(",", $freq);
			$dbsnp_ra_freq = $freq_split[0];
			$dbsnp_aa_freq = $freq_split[$dbsnp_ac];#get AF for correct alternative allele (in case there are more than one)
		} else {
			$source = "undef";
			$freq = "NA";
			$dbsnp_ra_freq = "NA";
			$dbsnp_aa_freq = "NA";#get AF for correct alternative allele (in case there are more than one)
		}

		if($dbsnp_ra_freq ne "NA" && $dbsnp_aa_freq ne "NA" && $dbsnp_aa_freq ne "."){
			# check if dbsnp reference allele is major allele; if not switch alleles and switch effect sign
			if($dbsnp_ra_freq ne "." && $dbsnp_ra_freq>=0.5){# ref > 0.5 -> alleles stay
				$or = sprintf("%.5f", $or);
				print OUT "$rsid\t$ch\t$p\t$ea\t$oa\t$dbsnp_aa_freq\t$effect\t$or\t$se\t$pval\t$trait\t$ss\t$source\n";
				$case1++;
			} elsif ($dbsnp_ra_freq<0.5 && $dbsnp_aa_freq>=0.5) { # ref<0.5 and aa>=0.5 -> switch alleles and sign
				$effect = -$effect;
				$or = exp($effect); #calculate OR with effect
				$or = sprintf("%.5f", $or);
				print OUT "$rsid\t$ch\t$p\t$oa\t$ea\t$dbsnp_ra_freq\t$effect\t$or\t$se\t$pval\t$trait\t$ss\t$source\n";
				$case2++;
			} elsif ($dbsnp_ra_freq<0.5 && $dbsnp_aa_freq<0.5 && $dbsnp_ra_freq>=$dbsnp_aa_freq) { #ref and aa <0.5 but ref>=aa -> alleles stay
				$or = sprintf("%.5f", $or);
				print OUT "$rsid\t$ch\t$p\t$ea\t$oa\t$dbsnp_aa_freq\t$effect\t$or\t$se\t$pval\t$trait\t$ss\t$source\n";
				$case3++;
			} elsif ($dbsnp_ra_freq<0.5 && $dbsnp_aa_freq<0.5 && $dbsnp_ra_freq<$dbsnp_aa_freq) { #ref and aa <0.5 and ref<aa -> alleles switch
				$effect = -$effect;
				$or = exp($effect); #calculate OR with effect
				$or = sprintf("%.5f", $or);
				print OUT "$rsid\t$ch\t$p\t$oa\t$ea\t$dbsnp_ra_freq\t$effect\t$or\t$se\t$pval\t$trait\t$ss\t$source\n";
				$case4++;
			}
		} else {
			$or = sprintf("%.5f", $or);
			print OUT "$dbsnp_rs\t$ch\t$p\t$ea\t$oa\t$dbsnp_aa_freq\t$effect\t$effect_sign\t$or\t$se\t$pval\t$trait\t$ss\t$gb\t$source\n";
			$case5++;
		}
		$def++;

	} elsif($add_allele_freq == 1 && defined($dbsnp_hash{$whole_string_reverse_alleles})){# -> study entry found in dbSNP but with reversed alleles

		$reverse_match++;

		# switch $ea and $oa
		my $ea_copy=$ea;
		$ea=$oa;
		$oa=$ea_copy;

		# get contents of dbSNP hash saved earlier to extract source/MAF
		# dbsnp hash content: 1:12345:T:A => rs123;1000Genomes:0.9208,0.07917|GnomAD:0.9939,0.006087;1
		my $hash_content = $dbsnp_hash{$whole_string_reverse_alleles};

		my @hash_content_split = split(";", $hash_content);
		my $dbsnp_rs = $hash_content_split[0];
		my $dbsnp_af = $hash_content_split[1];
		my $dbsnp_ac = $hash_content_split[2];

		my $source; my $freq; my $dbsnp_ra_freq; my $dbsnp_aa_freq;
		if($dbsnp_af ne "undef"){
			my @af_split = split(":", $dbsnp_af);
			$source = $af_split[0];
			$freq = $af_split[1];
			my @freq_split = split(",", $freq);
			$dbsnp_ra_freq = $freq_split[0];
			$dbsnp_aa_freq = $freq_split[$dbsnp_ac];#get AF for correct alternative allele (in case there are more than one)
		} else {
			$source = "undef";
			$freq = "NA";
			$dbsnp_ra_freq = "NA";
			$dbsnp_aa_freq = "NA";
		}

		if($dbsnp_ra_freq ne "NA" && $dbsnp_aa_freq ne "NA" && $dbsnp_aa_freq ne "."){
			# check if dbsnp reference allele is major allele; if not switch alleles and switch effect sign
				if($dbsnp_ra_freq ne "." && $dbsnp_ra_freq>=0.5){# ref > 0.5 -> alleles stay
					$or = sprintf("%.5f", $or);
					print OUT "$rsid\t$ch\t$p\t$ea\t$oa\t$dbsnp_aa_freq\t$effect\t$or\t$se\t$pval\t$trait\t$ss\t$source\n";
					$case1_r++;
				} elsif ($dbsnp_ra_freq<0.5 && $dbsnp_aa_freq>=0.5) { # ref<0.5 and aa>=0.5 -> switch alleles and sign
					$effect = -$effect;
					$or = exp($effect); #calculate OR with effect
					$or = sprintf("%.5f", $or);
					print OUT "$rsid\t$ch\t$p\t$oa\t$ea\t$dbsnp_ra_freq\t$effect\t$or\t$se\t$pval\t$trait\t$ss\t$source\n";
					$case2_r++;
				} elsif ($dbsnp_ra_freq<0.5 && $dbsnp_aa_freq<0.5 && $dbsnp_ra_freq>=$dbsnp_aa_freq) { #ref and aa <0.5 but ref>=aa -> alleles stay
					$or = sprintf("%.5f", $or);
					print OUT "$rsid\t$ch\t$p\t$ea\t$oa\t$dbsnp_aa_freq\t$effect\t$or\t$se\t$pval\t$trait\t$ss\t$source\n";
					$case3_r++;
				} elsif ($dbsnp_ra_freq<0.5 && $dbsnp_aa_freq<0.5 && $dbsnp_ra_freq<$dbsnp_aa_freq) { #ref and aa <0.5 and ref<aa -> alleles switch
					$effect = -$effect;
					$or = exp($effect); #calculate OR with effect
					$or = sprintf("%.5f", $or);
					print OUT "$rsid\t$ch\t$p\t$oa\t$ea\t$dbsnp_ra_freq\t$effect\t$or\t$se\t$pval\t$trait\t$ss\t$source\n";
					$case4_r++;
				}
		} else {
			$or = sprintf("%.5f", $or);
			print OUT "$rsid\t$ch\t$p\t$ea\t$oa\t$dbsnp_aa_freq\t$effect\t$or\t$se\t$pval\t$trait\t$ss\t$source\n";
			$case5_r++;
		}
		$def++;

	} elsif ($add_allele_freq == 1){ # no match in dbSNP found -> print line with NoMatch flag
		$or = sprintf("%.5f", $or);
		print OUT "NA\t$ch\t$p\t$ea\t$oa\tNA\t$effect\t$or\t$se\t$pval\t$trait\t$ss\tNoMatch\n";
		$ndef++;
	}
 } else {	# when no allele freq needs to be added -> just switch alleles
		# check if dbsnp reference allele is major allele; if not switch alleles and switch effect sign
		if($eaf ne "." && $eaf ne "" && $eaf ne "NA" && $eaf<0.5){ # eaf<0.5 -> alleles stay
			#		print OUT "$rsid\t$ch\t$p\t$ea\t$oa\t$eaf\t$effect\t$or\t$se\t$pval\t$trait\t$ss\tStudy\n";
			print OUT "$rsid\t$ch\t$p\t$ea\t$oa\t$eaf\t$effect\t$or\t$se\t$pval\t$trait\t$ss\t$eaf_source\n";
			$case1_no++;
		} elsif ($eaf ne "." && $eaf ne "" && $eaf ne "NA" && $eaf>=0.5) { # eaf>=0.5 -> switch alleles, effect sign and recalc OR
			$effect = -$effect;
			$or = exp($effect); #calculate OR with effect
			$or = sprintf("%.5f", $or);
			$eaf = 1-$eaf;
			#			print OUT "$rsid\t$ch\t$p\t$oa\t$ea\t$eaf\t$effect\t$or\t$se\t$pval\t$trait\t$ss\tStudy\n";
			print OUT "$rsid\t$ch\t$p\t$oa\t$ea\t$eaf\t$effect\t$or\t$se\t$pval\t$trait\t$ss\t$eaf_source\n";
			$case2_no++;
		} else { # you should only land in here if $eaf has some non-numeric value
			$or = sprintf("%.5f", $or);
			#	print LOG "$line\tProvided eaf is not numeric\n";
			#			print OUT "$rsid\t$ch\t$p\t$ea\t$oa\tNA\t$effect\t$or\t$se\t$pval\t$trait\t$ss\tNoEAF\n";
			print OUT "$rsid\t$ch\t$p\t$ea\t$oa\tNA\t$effect\t$or\t$se\t$pval\t$trait\t$ss\t$eaf_source\n";
			$case3_no++;
		}
 }
}
close(IN);
close(OUT);


# Debugging and information about the processing are logged 

if($add_allele_freq == 1){

	my $perc = ($def + $ndef) != 0 ? ($def * 100) / ($def + $ndef) : 0;
	my $perc_round = sprintf("%.2f", $perc);

	print LOG "\nStudy entries for $file chr$chr_i found in dbSNP: $def\tNot found: $ndef\tPercentage: $perc_round%\n" if $add_allele_freq == 1;
	print LOG "Regular allele matches:\n";
	print LOG "Ref>0.5 (no switch):\t$case1\nRef<0.5 & AA>=0.5 (switch):\t$case2\nRef<0.5 & AA<0.5 & Ref>=AA (no switch):\t$case3\nRef<0.5 & AA<0.5 & Ref<AA (switch):\t$case4\ndbSNP file did not provide freq:\t$case5\n";
	print LOG "Switched allele matches:\t$reverse_match entries total\n";
	print LOG "Ref>0.5 (no switch):\t$case1_r\nRef<0.5 & AA>=0.5 (switch):\t$case2_r\nRef<0.5 & AA<0.5 & Ref>=AA (no switch):\t$case3_r\nRef<0.5 & AA<0.5 & Ref<AA (switch):\t$case4_r\ndbSNP file did not provide freq:\t$case5_r\n";
	print LOG "\n===========\n";
} else {
	print LOG "EAF<0.5 (no switch):\t$case1_no\nEAF>=0.5 (switch):\t$case2_no\nFile did not provide allele freq:\t$case3_no\n";
}

}

}
