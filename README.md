# Run PhyloWGS

https://github.com/morrislab/phylowgs

## Prerequisites

To run the PhyloWGS pipeline you require CNV calls from Battenberg (the subclones.txt output) and mutation calls from the lab ssm pipeline.

## Running
For a single sample run
```{bash}
phylowgs.pl -n sample_name -m mut_rda -c battenberg_subclones -o outdir -p purity [-s sex]

sex in 'male' or 'female' format. 
```
A wrapper to run multiple samples will work with the same manifest given to run.battenberg.pl, with a final column appended indicating the location of the sample mutect rda calls. It requires this `sample_info.txt`, the `battenberg results directory` and an `output directory`.

```{bash}
run.phylowgs.pl -i SAMPLE_INFO -b BBERG_DIR -o OUTDIR

Please provide sample info tab file:
        FORMAT:
SAMPLE_NAME     TUMOUR_BAM      NORMAL_BAM      GENDER (XX or XY)       MUTECT_RDA
```

## Output
Several useful outputs from this pipeline are produced.

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