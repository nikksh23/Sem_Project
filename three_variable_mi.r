
library(tidyverse)


# func discretize data
discretize_data <- function(seq, num_bins) {
  breaks <- quantile(seq, probs = seq(0, 1, length.out = num_bins + 1), na.rm = TRUE)
  discretized_seq <- cut(seq, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  return(discretized_seq)
}

# Func calculate entropy
calculate_entropy <- function(seq) {
  probs <- table(seq) / length(seq)
  entropy <- -sum(probs * log2(probs + .Machine$double.eps), na.rm = TRUE)
  return(entropy)
}

# Func joint entropy of two
calculate_joint_entropy <- function(seq1, seq2) {
  joint_table <- table(seq1, seq2) / length(seq1) + .Machine$double.eps
  joint_entropy <- -sum(joint_table * log2(joint_table), na.rm = TRUE)
  return(joint_entropy)
}

# Func joint entropy of three variables
calculate_triple_entropy <- function(seq1, seq2, seq3) {
  joint_table <- table(seq1, seq2, seq3) / length(seq1) + .Machine$double.eps
  joint_entropy <- -sum(joint_table * log2(joint_table), na.rm = TRUE)
  return(joint_entropy)
}

# Func three-variable Mutual Information (MMI)
calculate_mmi <- function(data, A, B, C) {
  # Extract sequences for A, B, and C
  seqA <- data[[A]]
  seqB <- data[[B]]
  seqC <- data[[C]]
  
  # individual entropies
  HA <- calculate_entropy(seqA)
  HB <- calculate_entropy(seqB)
  HC <- calculate_entropy(seqC)
  
  # pairwise joint entropies
  HAB <- calculate_joint_entropy(seqA, seqB)
  HBC <- calculate_joint_entropy(seqB, seqC)
  HCA <- calculate_joint_entropy(seqC, seqA)
  
  # joint entropy for all three variables
  HABC <- calculate_triple_entropy(seqA, seqB, seqC)
  
  # Pairwise MI and MMI
  MI12 <- HA + HB - HAB  # MI between seqA and seqB
  MI23 <- HB + HC - HBC  # MI between seqB and seqC
  MI13 <- HA + HC - HCA  # MI between seqA and seqC
  MMI <- HA + HB + HC - (HAB + HBC + HCA) + HABC
  
  return(list(MI12 = MI12, MI23 = MI23, MI13 = MI13, MMI = MMI))
}

data <- read_csv("filtered_genes.csv")

# Func average of quadruplets
average_quadruplets <- function(df, column_prefix) {
  df %>%
    select(starts_with(column_prefix)) %>%
    rowMeans() %>%
    as.numeric()
}

#  average for each sample
data <- data %>%
  mutate(
    CPZ_avg = average_quadruplets(., "CPZ"),
    CFIX_avg = average_quadruplets(., "CFIX"),
    AMK_avg = average_quadruplets(., "AMK"),
    NM_avg = average_quadruplets(., "NM"),
    DOXY_avg = average_quadruplets(., "DOXY"),
    CP_avg = average_quadruplets(., "CP"),
    AZM_avg = average_quadruplets(., "AZM"),
    TP_avg = average_quadruplets(., "TP"),
    ENX_avg = average_quadruplets(., "ENX"),
    CPFX_avg = average_quadruplets(., "CPFX")
  )



# relevant columns for MI
selected_data <- data %>%
  select(CPZ_avg, CFIX_avg, AMK_avg, NM_avg, DOXY_avg, CP_avg, AZM_avg, TP_avg, ENX_avg, CPFX_avg)


# Discretize the data
num_bins <- 10  # Set number of bins for discretization
discretized_data <- selected_data %>%
  mutate(across(everything(), ~ discretize_data(.x, num_bins)))

#  MMI and Pairwise MI for unique combinations of three columns
unique_combinations <- combn(colnames(discretized_data), 3, simplify = FALSE)  # Generate unique combinations
results <- tibble(A = character(), B = character(), C = character(), MI12 = numeric(), MI23 = numeric(), MI13 = numeric(), MMI = numeric())

for (combo in unique_combinations) {
  mi_values <- calculate_mmi(discretized_data, combo[1], combo[2], combo[3])
  results <- results %>% add_row(A = combo[1], B = combo[2], C = combo[3], MI12 = mi_values$MI12, MI23 = mi_values$MI23, MI13 = mi_values$MI13, MMI = mi_values$MMI)
}

# View all results
#print(results)

#results to a CSV file
write.csv(results, "three_variable_mi_filtered_genes_old_cutoff.csv", row.names = FALSE)

library(ggplot2)

# combination identifier for plotting
results <- results %>%
    mutate(Combination = paste(A, B, C, sep = "_"))

# Bar plot for MMI values
ggplot(results, aes(x = Combination, y = MMI)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels
    labs(title = "Three-Variable Mutual Information (MMI)", 
         x = "Variable Combinations", 
         y = "MMI Value")


