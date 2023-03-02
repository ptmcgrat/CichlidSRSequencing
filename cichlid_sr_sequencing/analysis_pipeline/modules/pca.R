# commandArgs will take in values like sys.argv in python and I can work with these values in the code below:
args <- commandArgs(trailingOnly = TRUE)
# load libraries
library(gdsfmt)
library(tidyverse)

# load in eigenvector and eigenvalue files
pca <- read_table("test.eigenvec", col_names = FALSE)
eigenval <- scan("test.eigenval")
pca <- pca[,-1]

# set names for metadata information
names(pca)[1] <- "sample"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
new_names <- read.csv(args[1])
pca$sample=new_names$metadata_id
eco_group <- rep(NA, length(pca$sample))
eco_group[grep("Mbuna", pca$sample)] <- "Mbuna"
eco_group[grep("Utaka", pca$sample)] <- "Utaka"
eco_group[grep("Deep_Benthic", pca$sample)] <- "Deep_Benthic"
eco_group[grep("Shallow_Benthic", pca$sample)] <- "Shallow_Benthic"
eco_group[grep("Diplotaxodon", pca$sample)] <- "Diplotaxodon"
eco_group[grep("Rhampochromis", pca$sample)] <- "Rhampochromis"

# define a pca tibble
pca <- as_tibble(data.frame(pca, eco_group))
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# Plot the PCA with ggplot
b <- ggplot(pca, aes(PC1, PC2, col = eco_group)) + geom_point(size = 3)
b <- b + scale_colour_manual(values = c("red", "blue", "green", "brown", "purple", "pink"))
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

# Save the output as png to the PCA_output directory
ggsave(sprintf("%s%s_pca.png", args[2], args[3]), plot=b, device='png')