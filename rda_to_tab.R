#!/usr/bin/Rscript

# Rscript to read in mutect data and convert into temporary tab file that can be converted into a vcf for parse_cnvs.py for phylowgs

library('getopt')
library('plyr')

usage <- function() {
  usage.text <-
    '\nUsage: rda_to_tab.R --mutect mutect_rdata --outfile output_file\n'
  return(usage.text)
}

params = matrix(
  c(
    'mutect', 'm', 1, 'character',
    'outfile', 'o', 1, 'character'
  ),
  ncol = 4,
  byrow = TRUE
)


opt = getopt(params)

if (is.null(opt$mutect)) { stop(usage()) }
if (is.null(opt$outfile)) { stop(usage()) }


data <- get(load(opt$mutect))


write.table(x = data, file = opt$outfile, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
