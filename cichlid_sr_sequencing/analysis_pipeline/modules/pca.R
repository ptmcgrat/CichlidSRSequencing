print('hi from R')
# args <- commandArgs(trailingOnly = TRUE)
metadata_file <- readLines(args[1])
# print(metadata_file)
library(gdsfmt)
library(tidyverse)
pca <- read_table("./test.eigenvec", col_names = FALSE)
eigenval <- scan("./test.eigenval")
pca <- pca[,-1]
# set names
names(pca)[1] <- "sample"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
new_names <- read.csv('filtered_samples.csv')
pca$sample=new_names$metadata_id
eco_group <- rep(NA, length(pca$sample))
eco_group[grep("Mbuna", pca$sample)] <- "Mbuna"
eco_group[grep("Utaka", pca$sample)] <- "Utaka"
eco_group[grep("Deep_Benthic", pca$sample)] <- "Deep_Benthic"
eco_group[grep("Shallow_Benthic", pca$sample)] <- "Shallow_Benthic"
eco_group[grep("Diplotaxodon", pca$sample)] <- "Diplotaxodon"
eco_group[grep("Rhampochromis", pca$sample)] <- "Rhampochromis"
pca <- as_tibble(data.frame(pca, eco_group))
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
b <- ggplot(pca, aes(PC1, PC2, col = eco_group)) + geom_point(size = 3)
b <- b + scale_colour_manual(values = c("red", "blue", "green", "brown", "purple", "pink"))
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
ggsave('test_out', plot=b, device='png')
