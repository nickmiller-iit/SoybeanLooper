#
# Subsample lines from a BED FILE
#
# Set seed to ensure the same subsample
#
#
# Order of command line args: input_BED_file num_lines_to_sample seed

library(dplyr)

args <- commandArgs(trailingOnly = T)

if (length(args) < 2){
  stop("At least 2 arguments must be supplied: input_BED_file num_lines_to_sample <optional: seed>")
}

inFile <- args[1]

sampleSize <- as.integer(args[2])

if (length(args) > 2){
  sd <- as.integer(args[3])
  set.seed(sd)
}

message <- paste("input file =", inFile, "number of lines to sample =", sampleSize)

if (exists("sd")){
  message <- paste(message, "random seed =", sd)
}

write(message, file = stderr())

inbed <- read.table(inFile)

outbed <- sample_n(tbl = inbed,
                   size = sampleSize)

write.table(outbed,
            file = stdout(),
            quote = F,
            sep = '\t',
            row.names = F,
            col.names = F
            )

