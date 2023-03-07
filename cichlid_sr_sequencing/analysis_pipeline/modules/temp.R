# commandArgs will take in values like sys.argv in python and I can work with these values in the code below:
args <- commandArgs(trailingOnly = TRUE)
# load libraries
# library(gdsfmt)
library(tidyverse)

# load in eigenvector and eigenvalue files
pca <- read_table("/Users/kmnike/Desktop/EigenFiles/LG1.eigenvec", col_names = FALSE)
eigenval <- scan("/Users/kmnike/Desktop/EigenFiles/LG1.eigenval")
# pca table has extra column of sample names. Remove this column:
pca <- pca[,-1]

# set names for metadata information
# Rename column 1 to "sample" and remaining columns as PC1-20:
names(pca)[1] <- "sample"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# Generate new names from metadata file
new_names <- read.csv('/Users/kmnike/Desktop/EigenFiles/csv_files/metadata_for_R.csv')

# Reorder the structure to match the order of samples in pca eignevector file
sample_ids <- pca[,1] # get the sample_ids from pca (i.e. the vcf file)
# Reorder the metadata_subset data frame to match the order of the sample IDs in the eigenvec data
metadata_ids_ordered <- new_names[match(sample_ids$sample, new_names$SampleID),] # reorganizes the metadata file to reflect order of samples in vcf file

# Set names of samples equal to metadata_id_ordered names
pca$sample=metadata_ids_ordered$metadata_id

# assign ego groups to the samples based on the grep value found in the metadata_id
eco_group <- rep(NA, length(pca$sample))
eco_group[grep("Mbuna", pca$sample)] <- "Mbuna"
eco_group[grep("Utaka", pca$sample)] <- "Utaka"
eco_group[grep("Deep_Benthic", pca$sample)] <- "Deep_Benthic"
eco_group[grep("Shallow_Benthic", pca$sample)] <- "Shallow_Benthic"
eco_group[grep("Diplotaxodon", pca$sample)] <- "Diplotaxodon"
eco_group[grep("Rhamphochromis", pca$sample)] <- "Rhamphochromis"

# define a pca tibble
pca <- as_tibble(data.frame(pca, eco_group))
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# Plot the PCA with ggplot
b <- ggplot(pca, aes(PC1, PC2, col = eco_group), ) + geom_point(size = 3)
b <- b + scale_colour_manual(values = c("red", "blue", "green", "brown", "purple", "pink"))
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
b <- b + labs(title = args[3])
# Save the output as png to the PCA_output directory
ggsave('test.png', device='png')

