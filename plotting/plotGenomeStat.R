#!/usr/bin/env Rscript


library("optparse")

option_list = list(
  make_option(c("-s", "--file_stat"), type="character", default=NULL, 
              help="file name with the statistic (3-col tab file with header: chromosome, position, stat)", metavar="character"),
  make_option(c("-f", "--file_fai"), type="character", default=NULL, 
              help="file name with chromosome lengths (fai)", metavar="character"),
  make_option(c("-t", "--plot_style"), type="character", default="single-line", 
              help="plot style: single-line, multiple-line or circos [default= %default]", metavar="character"),
  make_option(c("-m", "--max_chromosomes"), type="numeric", default=30, 
              help="maximum number of chromosomes to plot [default= %default]", metavar="numeric"),
  make_option(c("-n", "--name"), type="character", default="Statistic", 
              help="name of the statistic [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file_stat)) {
  print_help(opt_parser)
  stop("Input file with stat missing", call.=FALSE)
} else if (is.null(opt$file_fai)) {
  print_help(opt_parser)
  stop("Input fai file missing", call.=FALSE)
}

options(error = quote({dump.frames(to.file=TRUE); q()}))

library(dplyr)
library(tidyr)
library(ggplot2)
library(circlize)

#fname_stats <- 'stats.out'
#fname_fai <- 'S288C_reference_sequence_R64_mod.fasta.fai'
#max_chroms <- 30
#stat_name <- 'Distance'
#style <- 'single-line'

# Arguments
fname_stats <- opt$file_stat
fname_fai <- opt$file_fai
max_chroms <- opt$max_chromosomes
stat_name <- opt$name
style <- opt$plot_style 



# Reading fai file
df <- read.csv(fname_fai, sep='\t', header = FALSE)
#head(df)
#dim(df)

# Filtering upt maximum number of chromosomes
df2 <- df[1:max_chroms,1:2] %>% drop_na()
lengths <- c(0,as.vector(df2$V2)[1:length(df2$V2)-1])

additions <- c()
for (i in 1:length(lengths)) {
  subset_sum <- sum(lengths[1:i])
  additions[i] <- subset_sum
}
#additions
df2$addition <- additions
df2$chrom_mid <- df2$addition + (df2$V2/2)
#head(df2)

# Reading stats file
ds <- read.csv(fname_stats, sep='\t', header = TRUE)
#head(ds)
#dim(ds)
colnames(ds) <- c("chrom","mid","stat")
#head(ds)

# Keeping only filtered chromosomes in stat dataset
chroms <- as.vector(df2$V1)
ds2 <- ds %>% filter(chrom %in% chroms)
#tail(ds2)

# Keeping only filtered chromosomes (same as in stat dataset) in fai dataset
df3 <- df2 %>% filter(V1 %in% as.vector(unique(ds2$chrom)))
#tail(df3)
chroms <- as.vector(df3$V1)

# Converting postions
ds3 <- merge(ds2, df3, by.x = "chrom", by.y = "V1", sort = FALSE)
ds3$mid2 <- ds3$mid + ds3$addition
#ds3$chrom_mid <- ds3$addition + (ds3$V2/2)
ds4 <- ds3 %>% drop_na()
#head(ds4)
#dim(ds3)
#dim(ds4)
ds4$chrom <- factor(ds4$chrom, levels = chroms)

# Chromosome label postions on the plot
chromosome_mids <- unique(ds4$chrom_mid)
chromosomes <- unique(ds4$chrom)


# Plotting statistic according to the style
if (style == 'single-line') {
### Single line plot
png('plotGenomeStat_singleline.png',w=1200,h=400,res=150)
plot1 <- ggplot(ds4) + aes(x = mid2, y = stat, colour = chrom) +
  geom_point(size=1.5) +
  scale_colour_manual(values = rep(c("steelblue4","indianred2"),(length(chroms)/2) + 1)) +
  scale_x_continuous(breaks = chromosome_mids, labels = chromosomes) +
  labs(y = stat_name) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        axis.text.x = element_text(size=10,angle=45,hjust=1),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        legend.position = "none")
print(plot1)
dev.off()
} else if (style == 'multiple-line') {
### Multiple lines plot
png('plotGenomeStat_multipleline.png',w=900,h=1000,res=150)
plot2 <- ggplot(ds4) + aes(x = mid, y = stat, colour = chrom) +
  geom_point(size=1) +
  #scale_colour_manual(values = rep(c("steelblue4","indianred2"),(length(chroms)/2) + 1)) +
  scale_colour_manual(values = rep(c("steelblue4"),length(chroms))) +
  #scale_x_continuous(breaks = chromosome_mids, labels = chromosomes) +
  facet_wrap(~chrom,ncol=1,strip.position = "right") +
  labs(y = stat_name, x = 'Chromosome position (bp)') +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.title = element_text(size=12),
        legend.position = "none",
        strip.background = element_rect(fill=NA),
        strip.text.y.right = element_text(angle = 0))
print(plot2)
dev.off()
} else if (style == 'circos') {
### Circos plot
png('plotGenomeStat_circos.png',w=1000,h=1000,res=150)
circos.clear()
circos.par("track.height" = 0.2,"start.degree" = 90)
circos.initialize(ds4$chrom, x = ds4$mid)
circos.track(ds4$chrom, y = ds4$stat,
              panel.fun = function(x, y) {
                sector.index = get.cell.meta.data("sector.index")
                xcenter = get.cell.meta.data("xcenter")
                ycenter = get.cell.meta.data("ycenter")
                circos.text(xcenter, max(ds4$stat) + mm_y(5), sector.index, cex=0.75)
                })
col = rep(c("steelblue4","indianred2"),length(chroms)/2)
circos.trackPoints(ds4$chrom, ds4$mid, ds4$stat, col = col, pch = 16, cex = 0.5)
dev.off()
}



