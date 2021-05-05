#!/usr/local/bin/perl
use warnings;
use strict 'vars';
use Getopt::Long;
use Pod::Usage;
use Cwd;


###############################################################################
# PARSE COMMAND LINE OPTIONS
###############################################################################


=head1 OPTIONS

        -help   This message
        -g	GWAS summary stat file
        -l	Lead SNP file (simple file)


        Example:  perl ~/Projects/Scripts/Blueprint/GWAS_overlap/GWAS_WP10_overlap_summary_stat_phase2_GWAS_lead.pl  -g /nfs/team151_data03/PublicData/GWAS_summary_stats/GWAS_unified_generic_formatting/Multiple_sclerosis_IMSGC_2013_NatGen_Immunochip.txt -l /lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/QTL_summary_stat/mono_gene_nor_combat_peer_10_all_summary.txt.20052016.simple.txt

=cut

#### COMMAND LINE OPTIONS

my ($i_help, $GWAS_Input, $Lead_Input, $chr, $pos, $CellType, $tag, $lead_val);   # in- and output files


GetOptions('help'               =>  \$i_help,
           'g|gwas=s'			=>	\$GWAS_Input,
           'l|lead=s'			=>	\$Lead_Input);

check input
pod2usage(-exitstatus => 1, -verbose => 1) if $i_help;


&usage unless $Lead_Input;
die "ERROR: $Lead_Input do/does not exist.\n" if (!(-e $Lead_Input));

&usage unless $GWAS_Input;
die "ERROR: $GWAS_Input do/does not exist.\n" if (!(-e $GWAS_Input));

######################################## PATHS ##########################################
#																						#
#  			Paths and Env. variable for all required softwares and files				#
#########################################################################################

my $LD = "/lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/Re_Analysis/VCF_QC/Processed_VCF/Final_VCFs/ld/tags/rsid/r08";
# my $uk10_1kgp_LD="/lustre/scratch113/projects/uk10k/users/vi1/uk10k+1kg-phase3/LD/combined";
my $uk10_1kgp_LD="/lustre/scratch119/humgen/projects/uk10k/users/vi1/LD/uk10k+1kg-phase3";

#########################################################################################


my $Lead_file = `basename $Lead_Input`;
chomp($Lead_file);

if($Lead_file =~/(\w{4})\_(\w{3,5})\_.*/) {
$CellType = $1;
$tag = $2;
}
print "$CellType\t$tag\n";

my $GWAS_file = `basename $GWAS_Input`;
chomp($GWAS_file);

my $tmp = "tmp\_$CellType\_$tag\_$GWAS_file";

`mkdir -p $tmp`;

# `zcat $GWAS_Input| sed '1d' | awk '{if(\$7<=5e-8) print \$2":"\$3"\t"\$1}'|sort -rk2|  awk '!seen[\$1]++' > $tmp/GWAS.sig`;
## Please note that here the $GWAS_Input is not the summary stat, rather it is published GWAS loci /lustre/scratch114/projects/hematopoiesis/Blueprint/Analysis/kk8/GWAS_WP10_QTLs_overlap/GWAS_loci/*
`cat $GWAS_Input| sed '1d' | awk '{print \$2":"\$3"\t"\$1}'|sort -rk2|  awk '!seen[\$1]++' > $tmp/GWAS.sig`;

# `cat $QTL_Input | awk '{if(\$7<=0.05) print \$1}' | awk '!seen[\$0]++' > $tmp/WP10_QTL.sig` ;

`cat $Lead_Input | awk '{if(\$7<=0.05) print \$1}' | awk '!seen[\$0]++' > $tmp/WP10_QTL.lead` ; ##  CHANGE the output file as "$tmp/lead_snps.txt" if you want to active next line, although you don't have to as next line is just for double check.
# `perl -ne 'print if (\$seen{\$_} .= \@ARGV) =~ /10\$/' $tmp/lead_snps.txt $tmp/WP10_QTL.sig > $tmp/WP10_QTL.lead`; ## $tmp/lead_snps.txt and $tmp/WP10_QTL.lead would not be different, I have just checked twice.
# `awk 'FNR==NR {a[\$0]++; next} !a[\$0]' $tmp/lead_snps.txt $tmp/WP10_QTL.sig > $tmp/WP10_QTL.nolead`;



for ($chr=1; $chr<=22; $chr++) {

# 	`cat $LD/chr$chr | sort -uk 1b,1 > $tmp/chr$chr\_LD.txt` ;

	`cat $tmp/WP10_QTL.lead | grep ^$chr\: > $tmp/chr$chr\_WP10_QTL.lead`;
# 	`cat $tmp/WP10_QTL.nolead | grep ^$chr\: > $tmp/chr$chr\_WP10_QTL.nolead`;

	`awk 'NR==FNR {h[\$1] = \$2; next} {print \$1,h[\$1]}' $LD/chr$chr $tmp/chr$chr\_WP10_QTL.lead| awk '{if(\$2 != "") print \$0}'| awk '{print \$1"\t"\$1"\\n"\$2"\t"\$1}'| awk '{gsub(/,/, "\t"\$2"\\n"); print }'| awk '{print \$1"\t1\t"\$2}' > $tmp/chr$chr\_WP10_QTL_LD.lead`;
# 	`awk 'NR==FNR {h[\$1] = \$2; next} {print \$1,h[\$1]}' $LD/chr$chr $tmp/chr$chr\_WP10_QTL.nolead| awk '{if(\$2 != "") print \$0}'| awk '{print \$1"\t"\$1"\\n"\$2"\t"\$1}'| awk '{gsub(/,/, "\t"\$2"\\n"); print }'| awk '{print \$1"\t0\t"\$2}' > $tmp/chr$chr\_WP10_QTL_LD.nolead`;


# 	`join $tmp/chr$chr\_WP10_QTL.lead $tmp/chr$chr\_LD.txt | tr "," "\n" | tr " " "\n" |  awk '{print \$1"\t1"}' > $tmp/chr$chr\_WP10_QTL_LD.lead`;
# 	`join $tmp/chr$chr\_WP10_QTL.nolead $tmp/chr$chr\_LD.txt | tr "," "\n" | tr " " "\n"  | awk '{print \$1"\t0"}' > $tmp/chr$chr\_WP10_QTL_LD.nolead`;


	`cat $tmp/chr$chr\_WP10_QTL_LD.lead| awk '{print \$1}' > $tmp/chr$chr\_WP10_QTL_LD.lead.tmp`;
# 	`cat $tmp/chr$chr\_WP10_QTL_LD.nolead| awk '{print \$1}' > $tmp/chr$chr\_WP10_QTL_LD.nolead.tmp`;

	`awk 'FNR==NR {a[\$0]++; next} !a[\$0]' $tmp/chr$chr\_WP10_QTL_LD.lead.tmp $tmp/chr$chr\_WP10_QTL.lead |  awk '{print \$1"\t1\t"\$1}'> $tmp/chr$chr\_extra.lead` ;
# 	`awk 'FNR==NR {a[\$0]++; next} !a[\$0]' $tmp/chr$chr\_WP10_QTL_LD.nolead.tmp $tmp/chr$chr\_WP10_QTL.nolead |  awk '{print \$1"\t0\t"\$1}'> $tmp/chr$chr\_extra.nolead` ;



	## `cat $tmp/chr$chr\_WP10_QTL_LD.lead $tmp/chr$chr\_extra.lead $tmp/chr$chr\_WP10_QTL_LD.nolead $tmp/chr$chr\_extra.nolead > $tmp/chr$chr\_WP10_QTL_LD_final.txt`;

	## This ordering is very IMPORTANT: Keep the lead SNPs with 2nd column assigned 1 on the top of the file -- "$tmp/chr$chr\_WP10_QTL_LD_final.txt";
	## REASON: Many lead SNP or r^2 0.8 LD (with lead) SNPs are also found in non_lead SNP r^0.8 LD, because a SNP might be LD with a lead SNP but also
	## the same SNP might be non-LD with other non-lead SNP, so when we look for non-lead SNP this SNP appear both in lead and non-lead LD
	## files - "$tmp/chr$chr\_WP10_QTL_LD.lead" and "$tmp/chr$chr\_WP10_QTL_LD.nolead". Thus we want to delete the non-lead entry, since the SNP has the evidence to be a LD SNP with
	## at least one lead SNP. We are deleting the duplicate entry below by using "awk '!seen[\$1]++'"

	`cat $tmp/chr$chr\_WP10_QTL_LD.lead $tmp/chr$chr\_extra.lead |tr "_" ":"| awk -F":" '{print \$1":"\$2,\$4":"\$5"_"\$6"_"\$7}' | awk '{print \$1"\t"\$3"\t"\$4}' | awk '!seen[\$0]++'  > $tmp/chr$chr\_WP10_QTL_LD_final_lead_LD.sig`;


	`cat $tmp/chr$chr\_WP10_QTL_LD_final_lead_LD.sig | cut -f3| sort -u| awk -F"_" '{print \$1"\t"\$1"_"\$2"_"\$3}' > $tmp/chr$chr\_WP10_QTL_LD_final_lead_LD_copy.sig`;
	`zcat $uk10_1kgp_LD/rsid.EUR.chr$chr.r08.tags.gz| awk -v OFS="\t" '\$1=\$1'| sed '1d'| awk '{print \$8"\t"\$2":"\$3}'| awk '{gsub(/\\|/, "\t"\$2"\\n"); print }'| awk '{if(\$1=="NONE") print \$2"\t"\$2; else print \$0}'|tr "\t" ":"| awk -F":" '{if(\$5!="" && \$6!="") print \$5":"\$6"\t"\$1":"\$2; else print \$3":"\$4"\t"\$1":"\$2}' > $tmp/chr$chr\_uk10k_1kgp.txt`;
	`awk 'NR==FNR {h[\$1] = \$0; next} {print \$0,h[\$1]}' $tmp/chr$chr\_WP10_QTL_LD_final_lead_LD_copy.sig $tmp/chr$chr\_uk10k_1kgp.txt | awk '{if(\$3!="") print \$2"\t1\t"\$4}' > $tmp/chr$chr\_WP10_QTL_LD_final_lead_uk10k_1kgp_LD.sig`;
	`cat $tmp/chr$chr\_WP10_QTL_LD_final_lead_LD.sig $tmp/chr$chr\_WP10_QTL_LD_final_lead_uk10k_1kgp_LD.sig | awk '!seen[\$0]++'| sort -k3 > $tmp/chr$chr\_WP10_QTL_LD_final_lead.sig`;


	`awk 'NR==FNR {h[\$1] = \$2; next} {OFS="\t"; print \$1,h[\$1],\$2,\$3}' $tmp/GWAS.sig $tmp/chr$chr\_WP10_QTL_LD_final_lead.sig | awk -F"\t" '{if(\$2 != "") print \$0}' > $tmp/chr$chr\_$CellType\_$tag\_WP10_GWAS_overlap_SNPs_lead.txt`;

	# `cat $tmp/chr$chr\_WP10_QTL_LD.nolead $tmp/chr$chr\_extra.nolead| tr "_" ":"| awk -F":" '{print \$1":"\$2,\$4":"\$5"_"\$6"_"\$7}' | awk '{print \$1"\t"\$3"\t"\$4}' | awk '!seen[\$0]++'  > $tmp/chr$chr\_WP10_QTL_LD_final_nolead.sig`;
	# `awk 'NR==FNR {h[\$1] = \$2; next} {OFS="\t"; print \$1,h[\$1],\$2,\$3}' $tmp/GWAS.sig $tmp/chr$chr\_WP10_QTL_LD_final_nolead.sig | awk -F"\t" '{if(\$2 != "") print \$0}' > $tmp/chr$chr\_$CellType\_$tag\_WP10_GWAS_overlap_SNPs_nolead.txt`;


}

`cat $tmp/*_$CellType\_$tag\_WP10_GWAS_overlap_SNPs_lead.txt | awk '!seen[\$0]++' > $CellType\_$tag\_WP10_GWAS_overlap_SNPs.txt`;
# `cat $tmp/*_$CellType\_$tag\_WP10_GWAS_overlap_SNPs_lead.txt $tmp/*_$CellType\_$tag\_WP10_GWAS_overlap_SNPs_nolead.txt | awk '!seen[\$0]++' > $CellType\_$tag\_WP10_GWAS_overlap_SNPs.txt`;



########## REMOVING TEMP DIRECTORY ###############
`rm -rf $tmp`;


