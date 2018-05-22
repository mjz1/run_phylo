#!/hpf/tools/centos6/perl/5.20.1/bin/perl
# Wrapper script to launch the phylowgs pipeline.

# MULTISAMPLE BRANCH -- Adapt to allow for multi-sample runs

# Nick's multisample input test file: 
# /hpf/largeprojects/adam/projects/lfs/analysis/cnv/phylowgs/phylowgs_subsample_2018_04-21/multisample.test.input.tab
# Battenberg dir: 
# /hpf/largeprojects/adam/projects/lfs/analysis/cnv/battenberg/battenberg-2018-04-09

# ~/bin/run_phylowgs/run.phylowgs.pl -i /hpf/largeprojects/adam/matthew/phylo_mutlisample_test_2/multisample.test.input.tab -b /hpf/largeprojects/adam/projects/lfs/analysis/cnv/battenberg/battenberg-2018-04-09 -o /hpf/largeprojects/adam/matthew/phylo_mutlisample_test_3/ -n 100


use Getopt::Long;
use Data::Dumper;
use List::MoreUtils qw/ uniq/;
use File::Basename;
use Time::HiRes qw( time );
use File::Path qw/make_path/;
use File::Find::Rule;
use TorquePBS;
use Moose;
use strict;
use warnings;


my ($sample_info, $outdir, $bberg_dir);
my $subsamp = 5000;
my $priority_bed = "/home/mjz1/bin/run_phylowgs/comsic.v85.all.grch37.bed";


GetOptions(
	'i=s' => \$sample_info,
	'o=s' => \$outdir,
	'b=s' => \$bberg_dir,
	'n:s' => \$subsamp,
	'p:s' => \$priority_bed
	);

if (!defined($sample_info)) {
	print "Please provide sample info tab file:\n\tFORMAT:\nSAMPLE_NAME\tTUMOUR_BAM\tNORMAL_BAM\tGENDER (XX or XY)\tMUTECT_RDA\tMULTI_SAMPLE_FLAG\n";
	help();
	exit;
}

if (!defined($outdir)) {
	print "Please provide output directory\n";
	help();
	exit;
}

if (!defined($bberg_dir)) {
	print "Please provide battenberg directory\n";
	help();
	exit;
}


if (! -e ("$outdir")) {
	print "Creating output directory: $outdir\n";
	mkdir "$outdir";
}

my @jobs;

my %sample_hash;

print "Reading through sample manifest...\n";

open (my $info, '<', $sample_info);
while (my $row = <$info>) {
	chomp $row;
	my @fields = split(/\t/, $row);
	next if ($row !~ /\.bam/); #ensure we read the first line with alignment files

	my $sample_name = $fields[0];
	$sample_name =~ s/-/_/g;

	my $gender = $fields[3];
	my $mut = $fields[4];
	my $multi_flag = $fields[5];

	if (!-e $mut) {
		print "Mutect rda not found for $sample_name...skipping...\n";
		next;
	}
	my $sex;

	if ($gender eq "XX") {
		$sex = "female"
	} elsif ($gender eq "XY") {
		$sex = "male"
	} else {
		$sex = "auto"
	}

	my $bberg_sample_dir = "$bberg_dir/$sample_name"."_battenberg";
	my $subclones_out = (File::Find::Rule->file()->name("*_subclones.txt")->maxdepth("1")->in($bberg_sample_dir))[0];
	my $cellularity_out = (File::Find::Rule->file()->name("*_cellularity_ploidy.txt")->maxdepth("1")->in($bberg_sample_dir))[0];

	# Add check for final non empty cn out file:
	if (!-e $subclones_out or -z $subclones_out or !-e $cellularity_out) {
		print "Battenberg data not detected for $sample_name...Skipping...\n";
		next;
	}

	# Get the cellularity value
	my $cellularity = `cut -f1 $cellularity_out | tail -n 1`;
	chomp $cellularity;

	# Define sample directory
	my $sample_dir = "$outdir/$sample_name"."_phylowgs";


	my @sample_args = ($mut, $subclones_out, $cellularity, $sex);

	# # Push all sample data into a hash -- multi sample flag use for key when it is present, otherwise use sample name
	if (!defined($multi_flag)) {
		@{$sample_hash{'single'}{$sample_name}} = @sample_args;
	} else {
		push(@{$sample_hash{'multi'}{$multi_flag}{$sample_name}}, @sample_args);
	}
}
close $info;

# print Dumper \%sample_hash;

print "Submitting samples to phylowgs...\n";

# Need to retrieve parameters and construct submission command

for my $mode (sort keys %sample_hash) {
	if ($mode =~ /single/) {
		foreach my $sample (sort keys %{$sample_hash{'single'}}) {
			# Retrieve single sample parameters into array
			my @params = @{$sample_hash{'single'}{$sample}};

			# Define individual parameters
			my $mut = $params[0];
			my $subclones_out = $params[1];
			my $sample_dir = "$outdir/$sample"."_phylowgs.single";

			if (!-e $sample_dir) {
				make_path($sample_dir);
			}

			my $cellularity = $params[2];
			my $sex = $params[3];

			my $cmd = "phylowgs.pl -r $sample -n $sample -m $mut -c $subclones_out -o $sample_dir -p $cellularity -s $sex -b $subsamp -i $priority_bed";

			# Below is removed -- not fully tested but previously run samples should resume
			# # If run not completed remove any files present in sample dir
			# if (system("touch $sample_dir/touch; rm -rf $sample_dir/*") != 0) {
			# 	print "Unable to delete previous run data $sample\n";
			# 	next;
			# }


			# Submit to cluster
			print("Submitting $sample to PhyloWGS in single-sample mode\n");
			my $phylo = TorquePBS->new(
				jobname => "$sample.phylowgs",
				command => "$cmd",
				log_dir => "$sample_dir/log",
				root_dir => "$sample_dir",
				script_dir => "$sample_dir/scripts",
				memory => 32,
				parallel => 1,
				template_dir => '/hpf/largeprojects/adam/local/lib/perl5/auto/share/dist/TorquePBS/templates/',
				template => 'submit_to_pbs.template',
				queue => 'long'
				);
			$phylo->create_shell;
			my $phylo_job = $phylo->submit_shell;
		}
	} elsif ($mode =~ /multi/) {
		foreach my $multi_sample (sort keys %{$sample_hash{'multi'}}) {
			# Define output directory for the multisample run
			my $sample_dir = "$outdir/$multi_sample"."_phylowgs.multi";

			if (!-e $sample_dir) {
				make_path($sample_dir);
			}

			# Create arrays for each command
			my (@muts_a, @samples_a, @subclones_out_a, @cellularity_a);

			my $sex;

			# Loop over each of the samples
			foreach my $sample (sort keys %{$sample_hash{'multi'}{$multi_sample}}) {
				my @samp_array = @{$sample_hash{'multi'}{$multi_sample}{$sample}};
				push (@muts_a, $samp_array[0]);
				push (@subclones_out_a, $samp_array[1]);
				push (@cellularity_a, $samp_array[2]);
				push (@samples_a, $sample);
				$sex = $samp_array[3];
			}

			# Construct submission command
			my $cmd = join(" ", "phylowgs.pl -r $multi_sample -n", join(",", @samples_a), "-m", join(",", @muts_a), "-c", join(",", @subclones_out_a), "-o $sample_dir", "-p", join(",", @cellularity_a), "-s", $sex, "-b $subsamp -i $priority_bed");

			# # If run not completed remove any files present in sample dir
			# if (system("touch $sample_dir/touch; rm -rf $sample_dir/*") != 0) {
			# 	print "Unable to delete previous run data $multi_sample\n";
			# 	next;
			# }


			# Submit to cluster
			print("Submitting $multi_sample to PhyloWGS in multi-sample mode\n");
			my $phylo = TorquePBS->new(
				jobname => "$multi_sample.phylowgs",
				command => "$cmd",
				log_dir => "$sample_dir/log",
				root_dir => "$sample_dir",
				script_dir => "$sample_dir/scripts",
				memory => 32,
				parallel => 1,
				template_dir => '/hpf/largeprojects/adam/local/lib/perl5/auto/share/dist/TorquePBS/templates/',
				template => 'submit_to_pbs.template',
				queue => 'long'
				);
			$phylo->create_shell;
			my $phylo_job = $phylo->submit_shell;	
		}
	}
}



	# my $file_check = "$sample_dir/results/3B.txt.gz";

	# # Check for output files:...
	# if (-e $file_check) {
	# 	print "PhyloWGS output already detected for $sample_name...Skipping...\n";
	# 	next;
	# }




print "ALL DONE!\n";


sub help {
	print "Usage: $0 -i SAMPLE_INFO  -b BBERG_DIR -o OUTDIR [-n subsamp_n]\n";
}