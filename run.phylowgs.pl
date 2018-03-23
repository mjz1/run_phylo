#!/hpf/tools/centos6/perl/5.20.1/bin/perl
# Wrapper script to launch the phylowgs pipeline.

use Getopt::Long;
use Data::Dumper;
use List::MoreUtils qw/ uniq/;
use File::Basename;
use Time::HiRes qw( time );
use File::Find::Rule;
use TorquePBS;
use Moose;
use strict;
use warnings;


my ($sample_info, $outdir, $bberg_dir);


GetOptions(
	'i=s' => \$sample_info,
	'o=s' => \$outdir,
	'b=s' => \$bberg_dir
	);

if (!defined($sample_info)) {
	print "Please provide sample info tab file:\n\tFORMAT:\nSAMPLE_NAME\tTUMOUR_BAM\tNORMAL_BAM\tGENDER (XX or XY)\tMUTECT_RDA\n";
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

open (my $info, '<', $sample_info);
while (my $row = <$info>) {
	chomp $row;
	my @fields = split(/\t/, $row);
	next if ($row !~ /\.bam/); #ensure we read the first line with alignment files

	my $sample_name = $fields[0];
	$sample_name =~ s/-/_/g;

	my $gender = $fields[3];
	my $mut = $fields[4];

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

	my $sample_dir = "$outdir/$sample_name"."_phylowgs";
	my $file_check = "$sample_dir/outputs/3B.txt.gz";

	# Check for output files:...
	if (-e $file_check) {
		print "PhyloWGS output already detected for $sample_name...Skipping...\n";
		next;
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

	# Submit to cluster
	print("Submitting $sample_name to PhyloWGS\n");
	my $phylo = TorquePBS->new(
		jobname => "$sample_name.phylowgs",
		command => "phylowgs.pl -n $sample_name -m $mut -c $subclones_out -o $sample_dir -p $cellularity -s $sex",
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

close $info;
print "ALL DONE!\n";


sub help {
	print "Usage: $0 -i SAMPLE_INFO  -b BBERG_DIR -o OUTDIR\n";
}