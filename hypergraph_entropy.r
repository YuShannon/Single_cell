# install.packages("BiocManager")
# BiocManager::install("dplyr")
library(dplyr)     # Used for data manipulation
# BiocManager::install("entropy")
library(entropy)   # Used for entropy calculation
# BiocManager::install("ggplot2")
library(ggplot2)   # Used for creating the violin plot

setwd("C:/Users/shann/OneDrive - The University of Manchester/MScRP2/R/Third_trimester/yang_allxallentropy")

#Load expression data
X_data <- read.csv("C:/Users/shann/OneDrive - The University of Manchester/MScRP2/R/Third_trimester/yang_X.csv",row.names=1)
cell_obs<-read.csv("C:/Users/shann/OneDrive - The University of Manchester/MScRP2/R//Third_trimester/yang_obs.csv",row.names=1)
gene_var<-read.csv("C:/Users/shann/OneDrive - The University of Manchester/MScRP2/R/Third_trimester/yang_var.csv",row.names=1)

#### CTB_1 random ####

set.seed(3)

CTB_1 <- cell_obs[cell_obs$unique_cell_type == "CTB_1", ] # filters cells by cluster name
CTB_1_X <- X_data[rownames(X_data) %in% rownames(CTB_1), ] # subsets X_data by cell IDs only in CTB_1

# calculate SD for each column (gene)
sd.scores0<-apply(CTB_1_X,2,sd) # 2 = selects columns not rows
CTB_1_X_sd<-CTB_1_X[,which(sd.scores0>0)] # only keeps genes with SD>0; removes genes with no variation

CTB_1_result_df <- data.frame() # initialises empty df to store entropy results
cor_data_CTB_1 <- cor(CTB_1_X_sd, CTB_1_X_sd) # all genes x all genes correlation matrix 

num_columns <- ncol(CTB_1_X_sd) # stores the total n.genes (columns) in expression matrix for downstream sampling

for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  
  CTB_1_genes<- colnames(CTB_1_X_sd)[sample(1:num_columns,size=100)]
  # creates a vector of indices from 1 to total n.genes; randomly samples 100 of these indices; selects corresponding column names from CTB_1_X_sd
  
  # create submatrix with 100 sampled genes (rows) against all other genes (columns)
  cor_data_CTB_1_2 <- cor_data_CTB_1[na.omit(match(CTB_1_genes,rownames(cor_data_CTB_1))), # (row) finds position of 100 genes are in corr matrix
                                     -na.omit(match(CTB_1_genes,colnames(cor_data_CTB_1)))] # (column) excludes the 100 genes to avoid self-correlations
  
  # create vector of 1000 gene names
  CTB_1_sub_sample <- colnames(cor_data_CTB_1_2)[sample(1:(num_columns-100),size=1000)] 
  # creates seq of indices for remaining genes to sample from, and select corresponding gene names from cor_data_CTB_1_2
  
  # create final submatrix of corr matrix (keep rows, subset columns)
  cor_data_CTB_1_2 <- cor_data_CTB_1_2[,na.omit(match(CTB_1_sub_sample,colnames(cor_data_CTB_1_2)))] 
  # finds positions of the 1000 genes in cor_data_CTB_1_2 that match gene names in CTB_1_sub_sample
  
  bin_CTB_1<-abs(cor_data_CTB_1_2) # takes absolute value 
  bin_CTB_1[which(bin_CTB_1>sd(cor_data_CTB_1_2))]<-1 # applies threshold to the absolute corr matrix; values greater than SD are set to 1
  bin_CTB_1[which(bin_CTB_1!=1)]<-0 # rest are set to 0; this is the hypergraph incidence matrix
  hyp_CTB_1<-bin_CTB_1 %*% t(bin_CTB_1) # matrix multiplication of incidence matrix and its transpose; creates adjacency matrix 
  
  CTB_1_entropy_result<- entropy(hyp_CTB_1) # calculates entropy
  CTB_1_result_df <- rbind(CTB_1_result_df, CTB_1_entropy_result) # appends result to df
}

write.csv(CTB_1_result_df, "C:/Users/shann/OneDrive - The University of Manchester/MScRP2/R/Third_trimester/yang_allxallentropy/CTB_1_entropy_results_10000.csv", row.names=FALSE)

#### CTB_1 NuMP ####

set.seed(3) # Sets random seed for reproducibility; ensures sampling is the same every time the code is run

# Filter CTB_1 cells
CTB_1 <- cell_obs[cell_obs$all == "CTB_1", ] # Filters cell_obs for CTB_1
CTB_1_X <- X_data[rownames(X_data) %in% rownames(CTB_1), ] # Subsets X_data to include only CTB_1 cells, based on cell IDs (rownames)

# Keep genes with non-zero standard deviation
sd.scores0 <- apply(CTB_1_X, 2, sd) # Calculates standard deviation for each gene (column)
CTB_1_X_sd <- CTB_1_X[, which(sd.scores0 > 0)] # Filters out genes that have zero variability across CTB_1 cells, keeping only informative genes

# Identify mitochondrial genes present in data
MT_genes <- rownames(gene_var)[gene_var$mito_gene == "True"] # Selects mito genes from gene_var
MT_genes_in_data <- intersect(MT_genes, colnames(CTB_1_X_sd)) # Defines mito genes present in X data
rest_genes <- setdiff(colnames(CTB_1_X_sd), MT_genes_in_data) # Defines rest of transcriptome, exclusing mito genes; returns genes in x that are not in y for setdiff(x,y)

# Calculate MT genes vs ALL genes correlation matrix
cor_data_CTB_1 <- cor(CTB_1_X_sd[, MT_genes_in_data], CTB_1_X_sd) # Calculates Pearson correlation of MT genes vs all genes

# Loop for hypergraph entropy
CTB_1_result_df <- data.frame() # Initialises an empty data frame to store entropy results per iteration


for (i in 1:10000) { # Starts a loop that runs 10,000 times; each iteration simulates a "random experiment" to estimate entropy
  print(i)
  
  # Sample 100 MT genes
  sampled_MT <- sample(MT_genes_in_data, size = 100)
  
  # Sample 1000 rest-of-transcriptome genes
  sampled_rest <- sample(rest_genes, size = 1000)
  
  # Subset the correlation matrix: sampled MT genes (rows) × sampled rest genes (columns)
  cor_data_CTB_1_2 <- cor_data_CTB_1[na.omit(match(sampled_MT, rownames(cor_data_CTB_1))), 	 # Extract part of correlation matrix that contains the 100 MT genes
                                     na.omit(match(sampled_rest, colnames(cor_data_CTB_1)))] # and their correlations with the 1000 random non-MT genes sampled in this iteration
  
  
  # Binarize: apply absolute value and threshold; result is a binary incidence matrix (MT genes × rest genes)
  bin_CTB_1 <- abs(cor_data_CTB_1_2) # Takes absolute values of correlations
  bin_CTB_1[bin_CTB_1 > sd(cor_data_CTB_1_2)] <- 1 # Sets values above global SD to 1 (strong association)
  bin_CTB_1[bin_CTB_1 != 1] <- 0 # Everything else becomes 0
  
  # Create hypergraph projection (MT genes × MT genes adjacency)
  hyp_CTB_1 <- bin_CTB_1 %*% t(bin_CTB_1) # Projects the binary MT × rest matrix into a MT × MT adjacency matrix
  # This tells you how often two MT genes are jointly connected to the same rest genes (i.e., co-correlated)
  
  # Calculate entropy
  CTB_1_entropy_result <- entropy(hyp_CTB_1)
  
  # Store result
  CTB_1_result_df <- rbind(CTB_1_result_df, CTB_1_entropy_result) # Appends the entropy value to the results table
}

# Save results
write.csv(CTB_1_result_df, "C:/Users/shann/OneDrive - The University of Manchester/MScRP2/R/CTB_1_MT_entropy_results_10000.csv", row.names = FALSE)
