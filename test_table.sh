#!/bin/bash

module load python/3.5.2 R/3.3.2shlib

python phylowgs_table.py -s /hpf/largeprojects/adam/projects/lfs/analysis/cnv/phylowgs/phylowgs_multisample_plus_multichain_2018-06-26/4856_phylowgs.multi -m /hpf/largeprojects/adam/projects/lfs/analysis/ssm/4856_1.wgs/4856_1_annotated_filtered_clipped.rda /hpf/largeprojects/adam/projects/lfs/analysis/ssm/4856_3.wgs/4856_3_annotated_filtered_clipped.rda
