##Gamma model
library(readr)
library(glmmTMB)

all_indels_df <- read_delim("/mnt/loki/martin/frankie/RNA_seq_glm/2_R_files/2_output_datafiles/all_indel_nearest_gene_imbalance_data_order.tsv", 
                                           delim = "\t", escape_double = FALSE, 
                                           trim_ws = TRUE)

##removing indels with no genes within 100 kb

all_indels_lep_subset = subset(all_indels_df, all_indels_df$Order == 'Lepidoptera')
all_indel_no_na = subset(all_indels_df, complete.cases(all_indels_df))

##adding in our odds of imbalance ratio. We add 1 not 0.5 to be conservative on the log scale for the small values e.g 1/0, log(1.5/0.5) = 1.09 vs log(2/1) = 0.69 . We don't want to overestimate imbalance
all_indel_no_na$pos_odds <- (pmax(all_indel_no_na$number_of_reads_aligning_to_reference_haplotype, all_indel_no_na$number_of_reads_aligning_to_alternate_haplotype)+1)/(pmin(all_indel_no_na$number_of_reads_aligning_to_reference_haplotype, all_indel_no_na$number_of_reads_aligning_to_alternate_haplotype)+1)

##adding in normal odds ratio - adding 1 so when we log we don't get infinity values
all_indel_no_na$odds <- (all_indel_no_na$number_of_reads_aligning_to_reference_haplotype+1)/(all_indel_no_na$number_of_reads_aligning_to_alternate_haplotype+1)
all_indel_no_na$log_odds <- log(all_indel_no_na$odds)

#defining_weight_column
all_indel_no_na$depth_weight <- log(
  all_indel_no_na$number_of_reads_aligning_to_reference_haplotype + 
  all_indel_no_na$number_of_reads_aligning_to_alternate_haplotype + 1
)

#rescaling the dostance to TSS column -> too small values make the model unstable
all_indel_no_na$distance_to_TSS_kb <- all_indel_no_na$distance_to_TSS / 1000

## my indel length is super RH skewed so we log and then pray on the interpretation
all_indel_no_na$indel_length_log <- log(all_indel_no_na$indel_length)

##remove zeros - if there are any but i don't think there should be !
summary(all_indel_no_na$pos_odds)

##removing_crazy_outlier
all_indel_no_na = subset(all_indel_no_na, all_indel_no_na$gene_name != 'gene:ENSRTBG00005004283')

##scaling the pos_odds
all_indel_no_na$pos_odds_scaled <- all_indel_no_na$pos_odds / mean(all_indel_no_na$pos_odds)


gamma_model_all_individuals <- glmmTMB(
  pos_odds_scaled ~ indel_length_log + distance_to_TSS_kb + gene_location + (1 | ID / indel_chrom),
  family = Gamma(link = "log"),
  weights = depth_weight,
  data = all_indel_no_na,
  control = glmmTMBControl(
    optCtrl = list(iter.max = 2000, eval.max = 2000),
    parallel = 4
  )
)

#save the model
saveRDS(gamma_model_all_individuals, "/mnt/loki/martin/frankie/RNA_seq_glm/2_R_files/2_output_datafiles/gamma_model_all_individuals_lep_subset.rds")
#gamma_model_all_individuals <- readRDS("/mnt/loki/martin/frankie/RNA_seq_glm/2_R_files/2_output_datafiles/gamma_model_all_individuals.rds")

#analyse the model
#summary(gamma_model_all_individuals)
#gamma_model_all_individuals$fit$convergence
#diagnose(gamma_model_all_individuals)


##for one individual
#all_indel_no_na_subset = subset(all_indel_no_na, all_indel_no_na$ID =='ilAglIoxx1.1')

#gamma_model_one_individual <- glmmTMB(
#  pos_odds ~ indel_length_log + distance_to_TSS_kb + gene_location + (1 | indel_chrom) + (1| Order),
#  family = Gamma(link = "log"),
#  weights = depth_weight,
#  data = all_indel_no_na_subset
#)

#save the model
##saveRDS(gamma_model_one_individual, "/mnt/loki/martin/frankie/RNA_seq_glm/2_R_files/2_output_datafiles/gamma_model_one_individual.rds")
##gamma_model_one_individual <- readRDS("/mnt/loki/martin/frankie/RNA_seq_glm/2_R_files/2_output_datafiles/gamma_model_one_individual.rds")

#analyse the model
##summary(gamma_model_one_individual)
##gamma_model_one_individual$fit$convergence
##diagnose(gamma_model_one_individual)max(pair_counts_imbalance$V1,pair_counts_imbalance$V2)





