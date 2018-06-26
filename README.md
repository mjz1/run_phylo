# Run PhyloWGS

https://github.com/morrislab/phylowgs

## Prerequisites

To run the PhyloWGS pipeline you require CNV calls from Battenberg (the subclones.txt output) and mutation calls from the lab ssm pipeline.

## Running
```{bash}
phylowgs.pl [options]

      Required parameters:
        -n                          Sample name
        -r                          Run name
        -m                          Mutect RDA file
        -c                          Battenberg cnv file
        -o                          Output directory
        -p                          Sample purity

       Optional parameters:
        -gender                     Patient gender
        -subsamp                    Number of subsampled mutations [5000]
        -i                          Priority SSMs bed file [/home/mjz1/bin/run_phylowgs/comsic.v85.all.grch37.bed]
        -regions                    Which regions to use variants from. [normal_and_abnormal_cn]
        -burn                       Number of burnin samples [1000]
        -mcmc                       Number of mcmc iterations [2500]
        -metrhast                   Number of MH iterations [5000]
        -num_chains                 Number of chains for mulevolve.py [10]
```

### Multi-sample mode
PhyloWGS has a multi-sample mode where linked samples (e.g metastatic or multi-region or both) can be analyzed concurrently to produce a single phylogeny to explain all tumours.
To run multiple samples, provide comma seperated sample names/files where shown below. Additionally use the `-r` option to specify the merged sample name:

```{bash}
phylowgs.pl -n {s1,s2,s3} -m {m1,m2,m3} -c {c1,c2,c3} -p {p1,p2,p3} -r run_name [OPTIONS]
```

For full options list: 
```{bash}
phylowgs.pl -man
```

### Convenience wrapper
A wrapper to run multiple samples will work with the same manifest given to run.battenberg.pl, with at least one additional column appended indicating the location of the sample mutect rda calls. It requires this `sample_info.txt`, the `battenberg results directory` and an `output directory`. A final column can be used to indicate linked samples, which will automatically be run in multi-sample mode.

```{bash}
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
```

For full options list: 
```{bash}
run.phylowgs.pl -man
```


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