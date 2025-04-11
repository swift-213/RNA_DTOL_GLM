library(readr)
library(dplyr)
library(DESeq2)
library(edgeR)

Calculating_TPM <- function(feature_counts_raw) {
    #convert gene length to kilobases
    feature_counts_raw <- feature_counts_raw %>%
        mutate(Length_kb = Length/1000)
    #determine reads per kilobase
    feature_counts_raw <- feature_counts_raw %>%    
        mutate(RPK = gene_counts/Length_kb)    
    #get that scaling factor
    million_scaling_factor <- (sum(feature_counts_raw$RPK)/1000000)
    #determine transcripts per million
    feature_counts_raw <- feature_counts_raw %>%    
        mutate(TPM = RPK/million_scaling_factor)

    return(feature_counts_raw)
}


#reference_exon_only_feature_counts
reference_feature_counts = read.table("/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/featureCounts/iyBomPrat1.1_ERR7113577_GCA_930367275.1_pcrdups_aligned_exon_annotated_featurecounts.txt", header=TRUE)
#renaming the really long gene counts name to be just 'Gene_counts'
colnames(reference_feature_counts)[colnames(reference_feature_counts) == colnames(reference_feature_counts)[7]] <- "Gene_counts"
#checking it worked
colnames(reference_feature_counts)

#alternate exon only feature counts
alternate_feature_counts = read.table("/media/s1929681/Seagate/Frankie_DTOL_lep_project/RNA_seq_practise/featureCounts/iyBomPrat1.1_ERR7113577_GCA_930367225.1_pcrdups_aligned_exon_annotated_featurecounts.txt", header=TRUE)
#renaming the really long gene counts columns
colnames(alternate_feature_counts)[colnames(alternate_feature_counts) == colnames(alternate_feature_counts)[7]] <- "Gene_counts"

#the alinged exons
filtered_aligned_exons <- read_csv("~/Desktop/filtered_aligned_exons.csv")
#splitting the gene id columns to be just the name
filtered_aligned_exons$ref_gene_id_parsed <- sapply(strsplit(filtered_aligned_exons$ref_gene_id, split=":"), `[`, 2)
filtered_aligned_exons$alt_gene_id_parsed <- sapply(strsplit(filtered_aligned_exons$alt_gene_id, split=":"), `[`, 2)

colnames(filtered_aligned_exons)


#we now need to subset the alt and ref feature counts to keep only the ref and alt aligned exons!
filtered_aligned_exons$ref_gene_id_parsed
colnames(reference_feature_counts)
reference_feature_counts$Geneid == filtered_aligned_exons$ref_gene_id_parsed



reference_feature_counts_aligned_only <- reference_feature_counts[reference_feature_counts$Geneid %in% filtered_aligned_exons$ref_gene_id_parsed, ]
alternate_feature_counts_aligned_only <- alternate_feature_counts[alternate_feature_counts$Geneid %in% filtered_aligned_exons$alt_gene_id_parsed, ]

columns = c('reference_gene_ID', 'reference_count', 'alternate_count')

counts = data.frame(matrix(nrow = nrow(reference_feature_counts_aligned_only), ncol = length(columns)))
colnames(counts) = columns

for (i in 1:nrow(reference_feature_counts_aligned_only)){
    reference_gene_id = reference_feature_counts_aligned_only$Geneid[i]
    alt_gene_id = filtered_aligned_exons[filtered_aligned_exons$ref_gene_id_parsed == reference_gene_id,]$alt_gene_id_parsed[1]
    alternate_feature_counts_aligned_only[alternate_feature_counts_aligned_only$Geneid == alt_gene_id,]
    counts$reference_gene_ID[i] = reference_gene_id
    counts$reference_count[i] = reference_feature_counts_aligned_only$Gene_counts[i]
    counts$alternate_count[i] = alternate_feature_counts_aligned_only[alternate_feature_counts_aligned_only$Geneid == alt_gene_id,]$Gene_counts
}

#setting the gene id as the row name
rownames(counts) <- counts[[1]]
counts <- counts[, -1] 

#filter counts
#filter for any genes that have no aligned reads
counts_filtered = counts[rowSums(counts)> 0,]
#filter for reads that have more the x reads associated
counts_filtered = counts_filtered[rowSums(counts_filtered) >10, ]

nrow(counts_filtered)
#need to workoutif the counts are stat significant from each other.
#negative binomial test


#binomial test
# Run binomial test for each gene
pvals <- apply(counts_filtered, 1, function(x) {
  if (sum(x) == 0) return(NA)  # Skip genes with zero counts
  binom.test(x[1], sum(x), p = 0.5)$p.value  # Test reference bias
})


binomial_model <- glm(cbind(reference_count, sum(reference_count, alternate_count) - reference_count) ~ 1, 
                      family = binomial(link = "logit"), 
                      data = counts_filtered)

# Calculate residual deviance and degrees of freedom
deviance <- binomial_model$deviance
df <- binomial_model$df.residual

# Dispersion statistic (ratio of deviance to degrees of freedom)
dispersion_statistic <- deviance / df

# Print the result
dispersion_statistic

library(MASS)

AIC(binomial_model, nb_model)

summary(binomial_model)
# Adjust p-values for multiple testing (Benjamini-Hochberg correction)
pvals_adj <- p.adjust(pvals, method = "BH")

pvals_adj_filtered <- pvals_adj[pvals_adj != 'NA']

# Create results table
results <- data.frame(Gene = rownames(counts_filtered), Reference = counts_filtered[,1], Alternate = counts_filtered[,2], P_value = pvals, Adjusted_P = pvals_adj)
colnames(results)
# Filter for significant ASE genes (adjusted p < 0.05)
results$significant_genes <- ifelse(results$Adjusted_P < 0.05, "Significant", "Not_singificant")
nrow(significant_genes)
# Save results
write.csv(results, "ASE_binomial_results.csv")


col_vec <- ifelse(results$significant_genes == "Significant", "red", "black")

# Base R scatter plot
plot(counts_filtered$reference_count, counts_filtered$alternate_count,
     col = col_vec, pch = 16,  # pch=16 gives filled circles
     xlab = "Reference Allele Count", ylab = "Alternate Allele Count",
     main = "Allele-Specific Expression (ASE)",
     xlim = c(0, 15000),
     ylim = c(0, 10000))
     xlim = c(0, max(counts_filtered$reference_count) + 10),
     ylim = c(0, max(counts_filtered$alternate_count) + 10))


#negative binomial
library(MASS)
library(tidyr)

colnames(counts_filtered)

#takes the row names and makes it a column again
counts_filtered <- counts_filtered %>%
  tibble::rownames_to_column(var = "Gene_name")

#check for overdispersion
model_test <- glm.nb(Counts ~ Allele, family = poisson, data = long_data)
dispersion <- sum(residuals(model_test, type = "pearson")^2) / model_test$df.residual
print(dispersion)  # Should be ≈1 for Poisson to be valid


#creating the melted table
long_data <- pivot_longer(counts_filtered, cols = c("reference_count", "alternate_count"), 
                          names_to = "Allele", values_to = "Counts")

nb_model <- glm.nb(Counts ~ Allele, data = long_data)

# Calculate Pearson residuals
pearson_residuals <- residuals(nb_model, type = "pearson")

# Calculate sum of squared Pearson residuals
sum_squared_residuals <- sum(pearson_residuals^2)

# Calculate degrees of freedom
df_residual <- nb_model$df.residual

# Dispersion statistic (ratio of sum of squared residuals to degrees of freedom)
dispersion_statistic <- sum_squared_residuals / df_residual

# Print the result
dispersion_statistic


# Fit Poisson model (no overdispersion)
poisson_model <- glm(Counts ~ Allele, family = poisson, data = long_data)

# Compare AIC of both models
AIC(nb_model, poisson_model)


summary(nb_model)

# Fit null model (without Allele effect)
null_model <- glm.nb(Counts ~ 1, data = long_data)

# Compare full and null models using Likelihood Ratio Test
anova(null_model, nb_model, test = "Chisq")


#i think the negtative binomial doesn't work here as we don't have enough power?
counts_filtered$p_value <- sapply(counts_filtered$Gene_name, function(gene) {
  gene_data <- subset(long_data, Gene_name == gene)
  
  # Ensure the gene has at least two rows
  if (nrow(gene_data) < 2) return(NA)
  
  # Try fitting the full model
  model <- tryCatch(glm.nb(Counts ~ Allele, data = gene_data), error = function(e) NULL)
  if (!inherits(model, "glm")) return('model_issue')  # If model fitting failed, return NA
  
  # Try fitting the null model
  null_model <- tryCatch(glm.nb(Counts ~ 1, data = gene_data), error = function(e) NULL)
  if (!inherits(null_model, "glm")) return('null_model_issue')
  
  # Perform likelihood ratio test
  p_value <- tryCatch(anova(null_model, model, test = "Chisq")$`Pr(>Chi)`[2], error = function(e) NA)
  
  return(p_value)
})


sum(is.na(counts_filtered$p_value))
which(is.na(counts_filtered$p_value))

#comparing binomial to negtive binomial
binomial_model <- glm(cbind(counts_filtered$reference_count, counts_filtered$alternate_count) ~ 0.5, family = binomial(link = "logit"), data = counts_filtered)
library(MASS)
neg_binomial_model <- glm.nb(successes ~ predictor, data = your_data)


AIC(binomial_model, neg_binomial_model)

lrtest(binomial_model, neg_binomial_model)
