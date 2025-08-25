#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
BiocManager::install("readr")
library(readr)
library(dplyr)

setwd("C:/Users/shann/OneDrive - The University of Manchester/MScRP2/R/Third_trimester/yang_mitoxall")

#Load expression data
X_data <- read.csv("C:/Users/shann/OneDrive - The University of Manchester/MScRP2/R/Third_trimester/yang_X.csv",row.names=1)
cell_obs<-read.csv("C:/Users/shann/OneDrive - The University of Manchester/MScRP2/R//Third_trimester/yang_obs.csv",row.names=1)
gene_var<-read.csv("C:/Users/shann/OneDrive - The University of Manchester/MScRP2/R/Third_trimester/yang_var.csv",row.names=1)

#### CTB_1 ####  
set.seed(3)

# Filter cells by cluster name
CTB_1 <- cell_obs[cell_obs$unique_cell_type == "CTB_1", ]

# subset X_data by cell IDs only in CTB_1
CTB_1_X <- X_data[rownames(X_data) %in% rownames(CTB_1), ]

# calculate standard deviation for each column (gene) in CTB_1 expression data
sd.scores0<-apply(CTB_1_X,2,sd) #2 = selects columns not rows
CTB_1_X_sd<-CTB_1_X[,which(sd.scores0>0)] # only keep genes with SD>0; removes genes with no variation

# idenfity MT genes that are also velocity genes
CTB_1_mito <- rownames(gene_var)[gene_var$mito_gene == "True"]

# Run the correlation
CTB_1_cor_data <- cor(CTB_1_X_sd[, intersect(CTB_1_mito, colnames(CTB_1_X_sd))], # mito genes   
                      CTB_1_X_sd[, !colnames(CTB_1_X_sd) %in% CTB_1_mito]) # non-mito genes

#Binarise correlation values using sd of correlation matrix
bin_CTB_1 <-abs(CTB_1_cor_data) # takes absolute value of correlation matrix
bin_CTB_1[which(bin_CTB_1>sd(CTB_1_cor_data))]<-1 # binarises values, where values greater than SD of correlation matrix are set to 1
bin_CTB_1[which(bin_CTB_1!=1)]<-0 # all other values are set to 0; this creates the hypergraph incidence matrix

#Matrix multiplication to generate hypergraph adjacency matrix
hyp_CTB_1<-bin_CTB_1 %*% t(bin_CTB_1) # adjacency matrix; multiplies the matrix by its transpose

#ranking
CTB_1_rowsum_result<- rowSums(hyp_CTB_1) # sum each row in adjacency matrix
hist(CTB_1_rowsum_result) # histogram of row sums to visualise distribution
CTB_1_row_sums <- as.data.frame(CTB_1_rowsum_result) # converts row sums vector to a data frame

colnames(CTB_1_row_sums)[1] <- "rowsum" # rename column to rowsum
CTB_1_row_sums$unique_cell_type <- "CTB_1"  # adds a column with cell type label
CTB_1_row_sums$rank <- rank(CTB_1_row_sums$rowsum, ties.method = "average") # ranks genes based on row sums, using average ranking for tied values 
CTB_1_row_sums$norm_rank <- CTB_1_row_sums$rank / nrow(CTB_1_row_sums) # normalises ranks by dividing by total number of rows, creating a value between 0 and 1

CTB_1_row_sums <- CTB_1_row_sums[order(CTB_1_row_sums$rank, decreasing = TRUE), ] # sorts data frame in descending order of rank

write.csv(CTB_1_row_sums, "C:/Users/Team Knowhow/OneDrive - The University of Manchester/MScRP2/R/Third_trimester/yang_mitoxall/CTB_1_row_sums_results_ranked.csv", row.names=TRUE)

## plotting the heatmaps ##

#initially plot without row_split = 2,column_split = 2, so you can see your data
hm_CTB_1_1<-Heatmap(hyp_CTB_1, heatmap_width = unit(15, "cm"),heatmap_height = unit(15, "cm"),show_row_names = FALSE, show_column_names = FALSE, heatmap_legend_param = list(title = "Adjacency Score"), use_raster = T) #plot the heatmap, specifying the number of column and row clusters to split it into, raster improves performance at cost of resolution
draw(hm_CTB_1_1) #above comand will save the output to an object. use this command to plot the heatmap. plot(hm_FGR_1) also works

#then plot again when you can select the clusters in the data (represented by higher value colours based on the heatmap colourscale)
hm_CTB_1_2<-Heatmap(hyp_CTB_1, heatmap_width = unit(15, "cm"),heatmap_height = unit(15, "cm"),show_row_names = FALSE, show_column_names = FALSE, row_split = 2, column_split = 2, heatmap_legend_param = list(title = "Adjacency Score"),use_raster = T) #plot the heatmap, specifying the number of column and row clusters to split it into, raster improves performance at cost of resolution
ht <- draw(hm_CTB_1_2) #above comand will save the output to an object. 
row_order_CTB_1_2 <- row_order(ht)

# Save the heatmap as a PNG
png("C:/Users/Team Knowhow/OneDrive - The University of Manchester/MScRP2/R/Third_trimester/yang_mitoxall/CTB_1_1_heatmap.png", width = 2300, height = 2000, res = 300)  
draw(hm_CTB_1_1)
dev.off()

# Save the heatmap as a PNG
png("C:/Users/Team Knowhow/OneDrive - The University of Manchester/MScRP2/R/Third_trimester/yang_mitoxall/CTB_1_2_heatmap.png", width = 2300, height = 2000, res = 300)  
draw(hm_CTB_1_2)
dev.off()

## calculating the galois for the cluster ##

# extract row orders
row_order_CTB_1_2<-row_order(ht) #extract the row order from the heatmap into a list; each list element represents a cluster from the image and the number represents the row from the original data, arranged in the same order as the heatmap.
CTB_1_clust<-rownames(hyp_CTB_1)[row_order_CTB_1_2[[1]]] #extract row names from row numbers extracted in step above.

## calculate the galois correspondacnce
CTB_1_gal<-bin_CTB_1[match(CTB_1_clust,rownames(bin_CTB_1)),] #subset the original binary matrix to the rows from the cluster defined above

gal_98<-CTB_1_gal[,colSums(CTB_1_gal)>quantile(colSums(CTB_1_gal), 0.98)]
write.csv(gal_98, file = "C:/Users/Team Knowhow/OneDrive - The University of Manchester/MScRP2/R/Third_trimester/yang_mitoxall/galois/CTB_1_gal_98.csv", row.names = TRUE)
