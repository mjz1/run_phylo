# Run PhyloWGS

https://github.com/morrislab/phylowgs

## Prerequisites

To run the PhyloWGS pipeline you require CNV calls from Battenberg (the subclones.txt output) and mutation calls from the lab ssm pipeline.

## Running
### Single sample mode
For a single sample run: 

```{bash}
phylowgs.pl -n sample_name -m mut_rda -c battenberg_subclones -o outdir -p purity [-s sex] [-b subsamp_n]

sex in 'male' or 'female' format. 
```

### Multi-sample mode
PhyloWGS has a multi-sample mode where linked samples (e.g metastatic or multi-region or both) can be analyzed concurrently to produce a single phylogeny to explain all tumours.
To run multiple samples, provide comma seperated sample names/files where shown below. Additionally use the `-r` option to specify the merged sample name:

```{bash}
phylowgs.pl -n {sample_name.1,sample_name.2...sample_name.x} -m {mut_rda.1,mut_rda.2...mut_rda.x} -c {battenberg_subclones.1,battenberg_subclones.2...battenberg_subclones.x} -o outdir -p {purity.1,purity.2...purity.x} [-s sex] [-b subsamp_n]

sex in 'male' or 'female' format. 
```

### Convenience wrapper
A wrapper to run multiple samples will work with the same manifest given to run.battenberg.pl, with at least one additional column appended indicating the location of the sample mutect rda calls. It requires this `sample_info.txt`, the `battenberg results directory` and an `output directory`. A final column can be used to indicate linked samples, which will automatically be run in multi-sample mode.

```{bash}
run.phylowgs.pl -i SAMPLE_INFO -b BBERG_DIR -o OUTDIR [-n subsamp_n]

Please provide sample info tab file:
        FORMAT:
SAMPLE_NAME     TUMOUR_BAM      NORMAL_BAM      GENDER (XX or XY)       MUTECT_RDA        MULTI_SAMPLE
```

The `-b` for single sample or `-n` for wrapper runs can be specified on the command line to indicate the number of mutations you want to subsample in cases where there are any mutations. This is currently set to a default value 5,000, which allows the program to complete within 100 hours. 

## Output
Several outputs from this pipeline are produced for single sample mode (multi-sample in progress).

The PhyloWGS preparsing for this pipeline creates a convenient tab file of CNVs and cellular prevalance (better than battenberg output which is very convoluted). This file is located in the sample directory and is `[sample_name].cnvs`

The outputs directory within the sample directory contains the rest of the most important results:

* `1A.txt`: Inferred cellularity of the sample
* `1B.txt`: Best guess for number of cancerous populations
* `1C.txt`: For each cluster, 3 columns: the cluster number (starting from 1),  the typical number of ssms assigned to it and the phi value of the cluster.  
* `2A.txt`: Assignment of SSMs to clusters. Row i will list the cluster index for SSM number i (i.e., the SSM with the identifier `s<i - 1>`).
* `2B.txt`: Full NxN co-clustering matrix
* `3A.txt`: For each of the best guess clusters, 2 columns: cluster ID and the cluster ID of its parent (0 is root node)
* `3B.txt`: NxN Ancestor-decedent matrix. Entry i,j = The probability that i is in node that is an ancestor of node containing j. 

There are other methods for viewing results. Refer to PhyloWGS page for this.