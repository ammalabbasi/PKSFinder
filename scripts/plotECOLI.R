#!/usr/bin/env Rscript
library(circlize)
library(Cairo)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
coverage_file <- args[1]
cytoband_file <- args[2]
output_pdf <- args[3]


# Read coverage data - adjust if necessary
coverage <- read.csv(coverage_file, sep="\t", stringsAsFactors = FALSE, header=0)
colnames(coverage) <- c('chr', 'start', 'end', 'value')
coverage <- coverage[, c('chr', 'start', 'end', 'value')]

# Read cytoband data
cytoband.df <- read.csv(cytoband_file, sep = "\t", stringsAsFactors = FALSE)
cytoband.df$start <- as.numeric(cytoband.df$start)
cytoband.df$end <- as.numeric(cytoband.df$end)
cytoband.df <- cytoband.df[, c('chrom', 'start', 'end', 'Name', 'gieStain')]



coverage <- coverage[((coverage$chr == 'NC_017628.1') & (coverage$start >= min(cytoband.df$start)) & (coverage$end <= max(cytoband.df$end))), ]
   
if(nrow(coverage)<1){
coverage <- read.csv(coverage_file, sep="\t", stringsAsFactors = FALSE, header=0)
colnames(coverage) <- c('chr', 'start', 'end', 'value')
coverage <- coverage[, c('chr', 'start', 'end', 'value')]
coverage$value <- rep(0,nrow(coverage))
}


# Extract sample name from coverage filename (optional)
sample_name <- tools::file_path_sans_ext(basename(coverage_file))

# Open PDF device
CairoPDF(output_pdf, width = 5, height = 5)

circos.clear()
circos.par(gap.after=3, gap.before=3, ADD = TRUE)
circos.initializeWithIdeogram(cytoband.df)

# You can uncomment this if you want labels inside the ideogram
# circos.labels(cytoband.df$chrom, x = as.numeric(cytoband.df$start), 
#               labels = cytoband.df$Name, side = "inside", facing="clockwise", niceFacing=TRUE)

if(length(unique(coverage$value)) == 1){
  circos.genomicTrack(coverage, ylim = c(0, 1), numeric.column = 4,
    panel.fun = function(region, value, ...) {
      circos.genomicLines(region, value, type = "h", numeric.column = 1, col = "#008000", track.height = 0.3)
    })  
} else {
  circos.genomicTrack(coverage, numeric.column = 4,
    panel.fun = function(region, value, ...) {
      circos.genomicLines(region, value, type = "h", numeric.column = 1, col = "#008000", track.height = 0.3)
    })
}

max_y_value <- max(coverage$value) + 100
y_axis_breaks <- seq(0, max_y_value, by = 1e4)

circos.yaxis(
  side = "left",
  at = y_axis_breaks,
  sector.index = unique(coverage$chr),
  labels.cex = 0.1,
  labels.niceFacing = TRUE,
  labels.col = "black"
)

title(sample_name,  cex = 1.5, side = "top", adj = 0.5)

circos.clear()
dev.off()
