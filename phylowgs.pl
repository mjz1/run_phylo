#!/hpf/tools/centos6/perl/5.20.1/bin/perl
#Parse battenberg CNV output and Shlien lab mutect RDATA to create phylowgs input and run for an individual sample

# To test subsampling:
# perl /home/mjz1/bin/run_phylowgs/phylowgs.pl -n kics_32_273811_274026 -m /hpf/largeprojects/adam/projects/kics/data/wgs_ssms/0032/N_-_274026+T_-_273811/273811_annotated_filtered_clipped.rda -c /hpf/largeprojects/adam/projects/icgc_tcga_datasets/RNAmp/kics/data/battenberg_3.3.2cgp_2.2.8bberg/kics_32_273811_274026_battenberg/kics_32_273811_274026_subclones.txt -o /hpf/largeprojects/adam/matthew/test_subsamp -p 0.79553 -s male

use Getopt::Long;
use Data::Dumper;
use List::MoreUtils qw/ uniq/;
use File::Basename;
use File::Path qw/make_path/;
use TorquePBS;
use Moose;
use strict;
use warnings;

my ($sample_name, $mut, $cnv, $outdir, $cellularity);

# Set default options
my $gender = "auto";

my $subsamp = 5000;



GetOptions(
	'n=s' => \$sample_name,
	'm=s' => \$mut,
	'c=s' => \$cnv,
	'o=s' => \$outdir,
	'p=s' => \$cellularity,
	's:s' => \$gender,
	'b:s' => \$subsamp
	);

if (!defined($sample_name)) {
	print "Please provide sample_name\n";
	help();
	exit;
}

if (!defined($mut)) {
	print "Please provide mutect rdata frame\n";
	help();
	exit;
}

if (!defined($outdir)) {
	print "Please provide output directory\n";
	help();
	exit;
}

if (!defined($cnv)) {
	print "Please provide cnv data \n";
	help();
	exit;
}


if (! -e ($outdir)) {
	print "Creating output directory: $outdir\n";
	make_path($outdir);
}

my @jobs;

my $tmpfile = "$outdir/$sample_name".".tmp";
my $vcffile = "$outdir/$sample_name".".vcf";
my $cnv_output = "$outdir/$sample_name".".cnvs";
my $cnvs_final = "$outdir/$sample_name"."_cnv_data.txt";
my $variants_final = "$outdir/$sample_name"."_ssm_data.txt";
my $params_json = "$outdir/$sample_name"."_params.json";
my $top_k_trees = "$outdir/$sample_name".".top_k_trees";
my $clonal_freq = "$outdir/$sample_name".".clonalFrequencies";
my $tmpdir = "$outdir/tmp";
my $nonsubsamp_variants = "$outdir/$sample_name"."_ssm_data_nonsubsamp.txt";
my $trees_zip = "$outdir/trees.zip";
my $posthoc_out = "$outdir/posthoc_assign.json";
my $python = "/hpf/tools/centos6/python/2.7.12/bin/python2";

if (!-e $tmpdir) {
	make_path($tmpdir);
}


# Convert mutect rdata file into tmp file
print "Converting mutect rda to tmpfile...\n";
system("Rscript /hpf/largeprojects/adam/local/bin/mutect_to_vcftmp.R --mutect $mut --outfile $tmpfile\n") == 0 or die "Failed to convert mutect rda file\n";

print "Converting tmpfile to mutectvcf4.1...\n";

# Convert tmp file into mutect vcf
open (my $vcfout, '>', $vcffile);

#print VCF4.1 header

print $vcfout "##fileformat=VCFv4.1\n##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP Membership\">\n##FORMAT=<ID=TD,Number=.,Type=Integer,Description=\"Tumor allelic depths for the ref and alt alleles in the order listed\">\n##FORMAT=<ID=ND,Number=.,Type=Integer,Description=\"Normal allelic depths for the ref and alt alleles in the order listed\">\n##INFO=<ID=TR,Number=1,Type=Integer,Description=\"Approximate tumor read depth; some reads may have been filtered\">\n##INFO=<ID=NR,Number=1,Type=Integer,Description=\"Approximate normal read depth; some reads may have been filtered\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sample_name\n";

#read in mutation data and print out to vcf
# Mutation counter
my $cnt = 0;
open (my $mut_f, '<', $tmpfile);
while (my $row = <$mut_f>) {
	chomp $row;
	next if ($row =~ /annovar/);
	$cnt++;
	my @data = split(/\t/, $row);
	my $snv_id = $data[0];
	my $chr = $data[1];
	my $pos = $data[2];
	my $id = ".";
	my $ref = $data[3];
	my $alt = $data[4];
	my $t_ref_counts = $data[5];
	my $t_alt_counts = $data[6];
	my $n_ref_counts = $data[7];
	my $n_alt_counts = $data[8];
	my $t_depth = ($t_ref_counts + $t_alt_counts);
	my $n_depth = ($n_ref_counts + $n_alt_counts);
	my $qual = ".";
	my $filter = ".";
	my $info = "SOMATIC";
	my $format = "TD:ND:TR:NR";
	my $sample_data = "$t_ref_counts".","."$t_alt_counts".":"."$n_ref_counts".","."$n_alt_counts".":"."$t_depth".":"."$n_depth";

	# if ($chr =~ "X") {
	# 	$chr = 23;
	# }

	# if ($chr =~ "Y") {
	# 	$chr = 24;
	# }

	print $vcfout "$chr\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\t$format\t$sample_data\n";
}

close $mut_f;
close $vcfout;


# Run phylowgs cnv preparser parse_cnvs.py
print "Parsing Battenberg CNVs...\n";
system("module load phylowgs/bc4e098; $python /hpf/tools/centos6/phylowgs/bc4e098/parser/parse_cnvs.py -f battenberg-smchet -c $cellularity --cnv-output $cnv_output $cnv\n") == 0 or die "Failed to parse battenberg CNVs: $!\n";

# Run phylowgs input creation script
print "Preparing PhyloWGS input...\n";
if ($cnt > $subsamp) {
	print "Subsampling $subsamp mutations from $cnt for $sample_name.\n";
	system("module load phylowgs/bc4e098;$python /hpf/tools/centos6/phylowgs/bc4e098/parser/create_phylowgs_inputs.py -s $subsamp --cnvs $sample_name=$cnv_output --vcf-type $sample_name=mutect_tcga $sample_name=$vcffile --output-cnvs $cnvs_final --output-variants $variants_final --nonsubsampled-variants $nonsubsamp_variants --output-params $params_json --sex $gender --verbose") == 0 or die "Failed to run create_phylowgs_inputs.py\n";
} else {
	system("module load phylowgs/bc4e098;$python /hpf/tools/centos6/phylowgs/bc4e098/parser/create_phylowgs_inputs.py --cnvs $sample_name=$cnv_output --vcf-type $sample_name=mutect_tcga $sample_name=$vcffile --output-cnvs $cnvs_final --output-variants $variants_final --output-params $params_json --sex $gender --verbose") == 0 or die "Failed to run create_phylowgs_inputs.py\n";
}

# Run phylowgs
print "Running PhyloWGS...\n";
system("module load phylowgs/bc4e098;cd $outdir/;$python /hpf/tools/centos6/phylowgs/bc4e098/evolve.py -k $top_k_trees -f $clonal_freq -t $tmpdir $variants_final $cnvs_final\n") == 0 or die "Failed to run evolve.py\n";

# Write results
my $results_dir = "$outdir/results";
my $summ_json = "$results_dir/$sample_name".".summ.json.gz";
my $muts_json = "$results_dir/$sample_name".".muts.json.gz";
my $mutass = "$results_dir/$sample_name".".mutass.zip";

mkdir $results_dir;
print "Writing PhyloWGS results...\n";
system("module load phylowgs/bc4e098;$python /hpf/tools/centos6/phylowgs/bc4e098/write_results.py --include-ssm-names $sample_name $trees_zip $summ_json $muts_json $mutass") == 0 or die "Failed to write phylowgs results\n";


# Write human readable reports using morris lab smchet challenge code
my $output_dir = "$outdir/outputs";
mkdir $output_dir;

print "Writing PhyloWGS report...\n";
system("PYTHONPATH=/hpf/tools/centos6/phylowgs/bc4e098/;$python /home/mjz1/bin/smchet-challenge/create-smchet-report/write_report.py $summ_json $muts_json $mutass $output_dir") == 0 or die "Failed to write phylowgs report\n";

# if ($cnt > $subsamp) {
# 	print "Post hoc assigning non subsampled mutations $sample_name.\n";

# 	my $ssm_ids = `cut -f1 $nonsubsamp_variants | tail -n +2 | head -n 200 | tr '\n' ' '`;

# 	system("module load phylowgs/bc4e098;PYTHONPATH=/hpf/tools/centos6/phylowgs/bc4e098/;python2 /hpf/tools/centos6/phylowgs/bc4e098/misc/post_assign_ssm.py --cnvs $cnv_output $nonsubsamp_variants $trees_zip $ssm_ids > $posthoc_out") == 0 or die "Failed to run post_assign_ssm.py\n";
# 	exit;
# }


print "All done!\n";

sub help {
	print "Usage: $0 -n sample_name -m mut_rda -c battenberg_subclones -o outdir -p purity [-s sex] [-b subsamp_n]\n";
}
