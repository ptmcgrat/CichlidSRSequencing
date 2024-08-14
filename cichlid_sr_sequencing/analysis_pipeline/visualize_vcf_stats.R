# I think the easiest way to invoke this annoying script will be to copy it into the filtering_stats directory and pipe all outputs into a new directory called stats_visualizations or something 
# NOTE: This script will not run on Utaka due to issues installing tidyverse. Get the files on to Mzebra and run it there. 
library('tidyverse')
setwd(getwd())

# Create a stats_visualization dir if it doesn't already exist in the filtering_stats dir
if (!(dir.exists('stats_figures'))){
    dir.create(file.path('stats_figures'))
}
out_dir = paste(getwd(), '/stats_figures', sep='')

# for now I will hard_code in that all filenames will have the prefix "612_CVAnalysis" since that's the cohort I will run this on the first time I use the script
# NOTE: need to update script to allow for different prefix names

# stats for mean depth per sample (column)
ind_depth <- read_delim("CVAnalysis.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)
a <- ggplot(ind_depth, aes(depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
ggsave(path = out_dir, filename = 'mean_depth_per_sample.png')
sink('mean_depth_per_column_stats.txt')
print(summary(ind_depth$depth))
print(quantile(ind_depth$depth, c(.05, .10, .25, .50, .75, .90, .95)))
sink()

# stats for depth for each variant (depth in each row)
var_depth <- read_delim("CVAnalysis.ldepth.mean", delim = "\t",
           col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
b <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
b + theme_light() + xlim(0,50)
ggsave(path = out_dir, filename = 'mean_depth_per_row.png')
sink('mean_depth_per_row_stats.txt')
print(summary(var_depth$mean_depth))
print(quantile(var_depth$mean_depth, c(.05, .10, .25, .50, .75, .90, .95)))
sink()

# stats for depth per each depth, not averaged (by row)
var_site_depth <- read_delim("CVAnalysis.ldepth", delim = "\t",
           col_names = c("chr", "pos", "sum_depth", "sumsq_depth"), skip = 1)
c <- ggplot(var_site_depth, aes(sum_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
c + theme_light() + xlim(0, 25000) # the x_lim may need to increase as more samples are included
ggsave(path = out_dir, filename = 'total_depth_per_row.png')
sink('sum_depth_per_row_stats.txt')
print(summary(var_site_depth$sum_depth))
print(quantile(var_site_depth$sum_depth, c(.05, .10, .25, .50, .75, .90, .95)))
sink()

# stats for how much data is missing per sample (column)
ind_miss  <- read_delim("CVAnalysis.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
d <- ggplot(ind_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
d + theme_light()
ggsave(path = out_dir, filename = 'missing_data_per_sample.png')

# stats for missing data by variant (row)
var_miss <- read_delim("CVAnalysis.lmiss", delim = "\t",
            col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
e <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
e + theme_light()
ggsave(path = out_dir, filename = 'missing_data_per_variant.png')

# stats for allele frequency of variants
var_freq <- read_delim("CVAnalysis.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
f <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
f + theme_light() + xlim(0, 0.025)
ggsave(path = out_dir, filename = 'mean_allele_frequency_per_variant.png')
sink('allele_freq_stats.txt')
print(summary(var_freq$maf))
print(quantile(var_freq$maf, c(.05, .10, .25, .50, .75, .90, .95)))
sink()

# stats for heterozygosity & inbreeding coefficient
ind_het <- read_delim("CVAnalysis.het", delim = "\t",
           col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
g <- ggplot(ind_het, aes(f)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
g + theme_light()
ggsave(path = out_dir, filename = 'heterozygosity_per_individual.png')
sink('het_inbreeding_stats.txt')
print(summary(ind_het$f))
print(quantile(ind_het$f, c(.05, .10, .25, .50, .75, .90, .95)))
sink()

# stats for variant quality (row)
var_qual <- read_delim("CVAnalysis.lqual", delim = "\t",
           col_names = c("chr", "pos", "qual"), skip = 1)
h <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
h + theme_light() + xlim(0,1000)
ggsave(path = out_dir, filename = 'quality_per_variant.png')
sink('variant_quality_stats.txt')
print(summary(var_qual$qual))
print(quantile(var_qual$qual, c(.05, .10, .25, .50, .75, .90, .95)))
sink()

print('Script run successful')
