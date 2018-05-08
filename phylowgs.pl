#!/hpf/tools/centos6/perl/5.20.1/bin/perl
#Parse battenberg CNV output and Shlien lab mutect RDATA to create phylowgs input and run for an individual sample

# MULTISAMPLE BRANCH -- Adapt to allow for multi-sample runs


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

my ($sample_name, $mut, $cnv, $outdir, $cellularity, $run_name);

my $python = "/hpf/tools/centos6/python/2.7.11/bin/python2";

# Set default options
my $gender = "auto";

my $subsamp = 5000;
my $priority_bed = "/home/mjz1/bin/run_phylowgs/comsic.v85.all.grch37.bed";

# evolve.py options - currently set low for testing purposes
my $burn = 100; # 1000 default
my $mcmc = 250; # 2500 default
my $metrhast = 500; # 5000 default

GetOptions(
	'n=s' => \$sample_name,
	'r=s' => \$run_name,
	'm=s' => \$mut,
	'c=s' => \$cnv,
	'o=s' => \$outdir,
	'p=s' => \$cellularity,
	's:s' => \$gender,
	'b:s' => \$subsamp,
	'i:s' => \$priority_bed
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

if (!defined($run_name)) {
	print "Please provide a run name\n";
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

# Split each array value
my @samples = split(/,/, $sample_name);
my @muts = split(/,/, $mut);
my @cnvs = split(/,/, $cnv);
my @cellularitys = split(/,/, $cellularity);


# Detect how many samples
my $multi_num = @samples;
my $multi_index = $multi_num - 1; # Indexing is one less for loop later
my $multi_flag = 0;
if ($multi_num eq 1) {
	print "$multi_num sample submitted...\n";
} else {
	print "$multi_num samples submitted: ";
	print join("; ", @samples);
	print "\nRunning PhyloWGS in multi-sample mode...\n";
	$multi_flag = 1;
}

my $tmpdir = "$outdir/tmp";

if (!-e $tmpdir) {
	make_path($tmpdir);
}

# Step 1: Convert mutect RDAs into vcf files, and parse CNVS for each sample individually to prepare inputs for step 2

# Convert mutect rdata file into tmp file
my @vcfs_parse;
my @cnvs_parse;
foreach my $i (0..$multi_index) {
	my $vcffile = "$tmpdir/$samples[$i]".".vcf";
	push (@vcfs_parse, $vcffile);

	# Parse sample CNVs using phylo's provided parser
	my $cnv_output = "$outdir/$samples[$i]".".cnvs";
	push (@cnvs_parse, $cnv_output);

	if (-e $vcffile) {
		print "VCF file already found for $samples[$i]...\n";
	} else {
		print "Converting mutect rdas to vcfs: $samples[$i]\n";

		my $tmpfile = "$tmpdir/$samples[$i]".".tmp";
		my $rda = $muts[$i];

		system("module unload R; module load R/3.4.4; /hpf/tools/centos6/R/3.4.4/bin/Rscript /hpf/largeprojects/adam/local/bin/mutect_to_vcftmp.R --mutect $rda --outfile $tmpfile") == 0 or die "Failed to convert mutect rda file\n";

		# Convert tmp file into mutect vcf
		open (my $vcfout, '>', $vcffile);

		# print VCF4.1 header
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
	}

	### FIX ME ### --- parse cnvs has integer error
	if (-e $cnv_output) {
		print "Parsed CNVs already detected: $samples[$i]...\n";
	} else {
		print "Parsing battenberg CNVs: $samples[$i];\n";
		system("module unload python;module load phylowgs/bc4e098; $python /hpf/tools/centos6/phylowgs/bc4e098/parser/parse_cnvs.py -f battenberg-smchet -c $cellularitys[$i] --cnv-output $cnv_output $cnvs[$i]\n") == 0 or die "Failed to parse battenberg CNVs: $!\n";
	}

}

# Step 2: Run create_phylowgs_inputs.py to create ssm_data and cnv_data files
my $cnvs_final = "$outdir/$run_name"."_cnv_data.txt";
my $variants_final = "$outdir/$run_name"."_ssm_data.txt";
my $nonsubsamp_variants = "$outdir/$run_name"."_ssm_data_nonsubsamp.txt";
my $params_json = "$outdir/$run_name"."_params.json";


# Run phylowgs input creation script
if (-e $variants_final) {
	print "SSM data located: $variants_final...Skipping...\n";
} else {
	print "Preparing PhyloWGS input...\n";

	# Create priority SSMs file (newline seperated <chr>_<locus>)
	

	# Construct the command
	my @cnv_arg;
	my @vcftype_arg;
	my @vcf_arg;
	for my $i (0..$multi_index) {
		push(@cnv_arg, join(" ", "--cnv", join("=", $samples[$i], $cnvs_parse[$i])));
		push(@vcftype_arg, join(" ", "--vcf-type", join("=", $samples[$i], "mutect_tcga")));
		push(@vcf_arg, join("=", $samples[$i], $vcfs_parse[$i]));
	}

	my $create_phylo_inputs_cmd = join(" ", "module unload python; module load phylowgs/bc4e098;$python /hpf/tools/centos6/phylowgs/bc4e098/parser/create_phylowgs_inputs.py -s $subsamp", join(" ", @cnv_arg, @vcftype_arg, @vcf_arg), "--output-cnvs $cnvs_final --output-variants $variants_final --nonsubsampled-variants $nonsubsamp_variants --output-params $params_json --sex $gender --verbose");

	system("$create_phylo_inputs_cmd") == 0 or die "Failed to run create_phylowgs_inputs.py\n";
}


# Step 3: Run PhyloWGS (evolve.py)
my $top_k_trees = "$outdir/$run_name".".top_k_trees";
my $clonal_freq = "$tmpdir/$run_name".".clonalFrequencies";

if (-e $top_k_trees) {
	print "PhyloWGS already run...Skipping...\n";
} else {
	print "Running PhyloWGS...\n";
	system("module unload python; module load phylowgs/bc4e098;cd $outdir/;$python /hpf/tools/centos6/phylowgs/bc4e098/evolve.py -k $top_k_trees -f $clonal_freq --params $params_json -B $burn -s $mcmc -i $metrhast -t $tmpdir $variants_final $cnvs_final\n") == 0 or die "Failed to run evolve.py\n";
}


# Write results
my $results_dir = "$outdir/results";
my $summ_json = "$results_dir/$run_name".".summ.json.gz";
my $muts_json = "$results_dir/$run_name".".muts.json.gz";
my $mutass = "$results_dir/$run_name".".mutass.zip";
my $trees_zip = "$outdir/trees.zip";

mkdir $results_dir;
if (-e $mutass) {
	print "PhyloWGS results already found...Skipping...\n";
} else {
	print "Writing PhyloWGS results...\n";
	system("module unload python; module load phylowgs/bc4e098;$python /hpf/tools/centos6/phylowgs/bc4e098/write_results.py --include-ssm-names $sample_name $trees_zip $summ_json $muts_json $mutass") == 0 or die "Failed to write phylowgs results\n";
}

# # Write human readable reports using morris lab smchet challenge code
if ($multi_flag == 0) {
	my $output_dir = "$outdir/outputs";
	if (!-e $output_dir) {
		mkdir $output_dir;
	}

	if (-e "$output_dir/3B.txt.gz") {
		print "Output already written...ALL DONE...\n";
	} else {
		print "Writing PhyloWGS report...\n";
		system("module purge; module load phylowgs/bc4e098; PYTHONPATH=/hpf/tools/centos6/phylowgs/bc4e098/:/hpf/tools/centos6/python/2.7.11/lib/python2.7; $python /home/mjz1/bin/smchet-challenge/create-smchet-report/write_report.py $summ_json $muts_json $mutass $output_dir") == 0 or die "Failed to write phylowgs report\n";
	}
}

# Step 4: Run post-hoc assignment
# First check for non-subsampled variants
my $count = `wc -l < $nonsubsamp_variants`;
die "wc failed: $?" if $?;
chomp($count);

if ($count <= 1) {
	print "No subsampling performed. No post-hoc assignment necessary.\n";
} elsif ($count > 1) {
	print "$count additional mutations to be used for post-hoc assignment...\n"
}

# Get the SSM ids for post hoc assignment
my $ssm_ids = `cut -f1 $nonsubsamp_variants | tail -n +2 | head -n 1 | tr '\n' ' '`;

# print "$ssm_ids\n";

foreach my $i (0..$multi_index) {
	my $posthoc_out = "$results_dir/$samples[$i].posthoc.json";
	system("module load phylowgs/bc4e098;PYTHONPATH=/hpf/tools/centos6/phylowgs/bc4e098/;python2 /hpf/tools/centos6/phylowgs/bc4e098/misc/post_assign_ssm.py --cnvs $cnvs_parse[$i] $nonsubsamp_variants $trees_zip $ssm_ids > $posthoc_out") == 0 or die "Failed to run post_assign_ssm.py\n";
}

# Step 5: Parse JSON files and produce final outputs
###

print "All done!\n";

sub help {
	print "Usage: $0 -n sample_name -m mut_rda -c battenberg_subclones -o outdir -p purity [-s sex] [-b subsamp_n]\n";
}