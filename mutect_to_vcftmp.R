# Rscript to read in mutect data and convert into temporary tab file that can be converted into a vcf for parse_cnvs.py for phylowgs

library('getopt')
library('plyr')

usage <- function() {
  usage.text <-
    '\nUsage: mutect_to_vcftmp.R --mutect mutect_rdata --outfile output_file\n'
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

tmp <- data[, c("snvid", "annovar_chr", "annovar_start", "annovar_ref", "annovar_alt", "gatk_tumour_ref_count", "gatk_tumour_alt_count", "gatk_normal_ref_count", "gatk_normal_alt_count", "tumour_name")]

sample <- unique(tmp$tumour_name)[1]


write.table(x = tmp, file = opt$outfile, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)