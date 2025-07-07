args <- commandArgs(trailingOnly = TRUE)
coverage_file <- args[1]      # Full path to the coverage.bedgraph file
cytoband_file <- args[2]      # Cytoband file
output_pdf <- args[3]         # Output PDF file

library(circlize)

# Read cytoband file once
cytoband.df <- read.table(cytoband_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cytoband.df$start <- as.numeric(cytoband.df$start)
cytoband.df$end <- as.numeric(cytoband.df$end)
cytoband.df <- cytoband.df[, c('chrom', 'start', 'end', 'Name', 'gieStain')]

if (nrow(cytoband.df) == 0) {
    stop("Cytoband file is empty or invalid.")
}



# Extract sample name for title
base_name <- basename(coverage_file)
sample <- sub(".coverage\\.bedgraph$", "", base_name)

# Read coverage data
coverage <- read.table(coverage_file, sep = "\t", header=0, stringsAsFactors = FALSE)
colnames(coverage) <- c( 'chr', 'start','end', 'value')
coverage = coverage[,c('chr', 'start', 'end','value' )]

# Filter coverage
coverage <- coverage[coverage$chr == 'NC_017628.1' & coverage$end <= max(cytoband.df$end), ]
if (nrow(coverage) == 0) {
    stop("coverage file is empty or invalid.")
}

# Plot
pdf(output_pdf, width = 5, height = 5)
circos.clear()
circos.par(gap.after = 3)

circos.initializeWithIdeogram(cytoband.df)

circos.labels(cytoband.df$chrom, x = cytoband.df$start, labels = cytoband.df$Name,
              side = "inside", facing = "clockwise", niceFacing = TRUE)

circos.genomicTrack(coverage, numeric.column = 4,
    panel.fun = function(region, value, ...) {
        circos.genomicLines(region, value, type = "h", numeric.column = 1,
                            col = c("#008000"), track.height = 0.3)
    }
)

# Adjust y-axis parameters
max_y_value <- max(coverage$value, na.rm = TRUE)
y_axis_breaks <- seq(0, max_y_value, by = max_y_value / 5)



circos.yaxis(side = "left", at = y_axis_breaks, sector.index = unique(coverage$chr),
             labels.cex = 0.1, labels.niceFacing = TRUE, labels.col = "black")

title(sample, cex = 1.5)
circos.clear()
dev.off()

