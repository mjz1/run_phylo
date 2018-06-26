#!/hpf/tools/centos6/perl/5.20.1/bin/perl
# Wrapper script to launch the phylowgs pipeline.

# MULTISAMPLE BRANCH -- Adapt to allow for multi-sample runs

# Nick's multisample input test file: 
# /hpf/largeprojects/adam/projects/lfs/analysis/cnv/phylowgs/phylowgs_subsample_2018_04-21/multisample.test.input.tab
# Battenberg dir: 
# /hpf/largeprojects/adam/projects/lfs/analysis/cnv/battenberg/battenberg-2018-04-09

# perl ~/bin/run_phylowgs/run.phylowgs.pl -i /hpf/largeprojects/adam/matthew/phylo_mutlisample_test_2/multisample.test.input.tab -b /hpf/largeprojects/adam/projects/lfs/analysis/cnv/battenberg/battenberg-2018-04-09 -o /hpf/largeprojects/adam/matthew/phylo_mutlisample_test_3/


use Getopt::Long;
use Pod::Usage;
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

# Opt help and manual
my ($opt_help, $opt_man);

# multievolve.py options
my $burn = 1000; # 1000 default
my $mcmc = 2500; # 2500 default
my $metrhast = 5000; # 5000 default
my $num_chains = 10;

# June 1 master phylowgs commit: phylowgs/262325b
my $phylo_ver = "phylowgs/262325b/";

my $priority_bed = "/home/mjz1/bin/run_phylowgs/comsic.v85.all.grch37.bed";
my $regions = "normal_and_abnormal_cn";


GetOptions(
	'i=s' => \$sample_info,
	'o=s' => \$outdir,
	'b=s' => \$bberg_dir,
	'subsamp:s' => \$subsamp,
	'p:s' => \$priority_bed,
	'num_chains:s' => \$num_chains,
	'regions:s' => \$regions,
	'burn:s' => \$burn,
	'mcmc:s' => \$mcmc,
	'metrhast:s' => \$metrhast,
	'help!' => \$opt_help,
	'man!' => \$opt_man
	) or pod2usage(-verbose => 1) && exit;
pod2usage(-verbose => 1) && exit if defined $opt_help;
pod2usage(-verbose => 2) && exit if defined $opt_man;


# Assign memory per thread
my $memuse = 6 * $num_chains;

if (!defined($sample_info)) {
	print "Please provide sample info tab file:\n\tFORMAT:\nSAMPLE_NAME\tTUMOUR_BAM\tNORMAL_BAM\tGENDER (XX or XY)\tMUTECT_RDA\tMULTI_SAMPLE_FLAG\n";
	pod2usage(-verbose => 1) && exit;
	exit;
}

if (!defined($outdir)) {
	print "Please provide output directory\n";
	pod2usage(-verbose => 1) && exit;
	exit;
}

if (!defined($bberg_dir)) {
	print "Please provide battenberg directory\n";
	pod2usage(-verbose => 1) && exit;
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

			my $cmd = "phylowgs.pl -r $sample -n $sample -m $mut -c $subclones_out -o $sample_dir -p $cellularity --gender $sex --subsamp $subsamp -i $priority_bed --num_chains $num_chains --regions $regions --burn $burn --mcmc $mcmc --metrhast $metrhast";

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
				memory => $memuse,
				parallel => $num_chains,
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
			my $cmd = join(" ", "phylowgs.pl -r $multi_sample -n", join(",", @samples_a), "-m", join(",", @muts_a), "-c", join(",", @subclones_out_a), "-o $sample_dir", "-p", join(",", @cellularity_a), "--gender", $sex, "--subsamp $subsamp -i $priority_bed --regions $regions --burn $burn --mcmc $mcmc --metrhast $metrhast --num_chains $num_chains");

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
				memory => $memuse,
				parallel => $num_chains,
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

__END__

=head1 NAME

run.phylowgs.pl - Perform phylogenetic reconstruction from tumour sequencing data. Convenience wrapper for submitting multiple samples

=head1 SYNOPSIS

run.phylowgs.pl [options]

  Required parameters:
    -i                          Sample info tab file
    -o                          Output directory
    -b                          Battenberg output directory

   Optional parameters:
    -subsamp                    Number of subsampled mutations [5000]
    -p                          Priority SSMs bed file [/home/mjz1/bin/run_phylowgs/comsic.v85.all.grch37.bed]
    -regions                    Which regions to use variants from. [normal_and_abnormal_cn]
    -burn                       Number of burnin samples [1000]
    -mcmc                       Number of mcmc iterations [2500]
    -metrhast                   Number of MH iterations [5000]
    -num_chains                 Number of chains for mulevolve.py [10]

   Other:
    -help     -h  Brief help message.
    -man      -m  Full documentation.

=head1 OPTIONS

=over 8

=item B<-i>

Sample info tab file with the following columns: [SAMPLE_NAME,TUMOUR_BAM,NORMAL_BAM,GENDER (XX or XY),MUTECT_RDA,MULTI_SAMPLE_ID]

=item B<-o>

Directory to write output to.

=item B<-b>

Folder containing battenberg results for listed samples

=item B<-regions>

Which regions to use variants from. Allowable values are {normal_cn,normal_and_abnormal_cn,all}. (default: normal_and_abnormal_cn) 
                        
=item B<-burn>

Number of burnin samples (default 1000)

=item B<-mcmc>

Number of mcmc samples (default 2500)

=item B<-metrhast>

Number of MH iterations (default 5000)

=item B<-num_chains>

Number of chains for B<multievolve.py> (default 10)

=back

=head1 DESCRIPTION

B<run.phylowgs.pl> will attempt to run all steps of the phylowgs pipeline, starting with SSM and CNV calls to producing JSON results files.

Automatic wrapper for B<phylowgs.pl> for use with a sample tab file.

=cut