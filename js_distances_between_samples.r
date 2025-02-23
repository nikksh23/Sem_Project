# Load necessary libraries
library(tidyverse)
library(ggplot2)
library(reshape2)
library(scales)  # for scientific notation

# Read the data from CSV
data <- read_csv("Japan_Lab_Transcriptomic_Data1.csv")

# Function to calculate average of quadruplets
average_quadruplets <- function(df, column_prefix) {
    df %>%
        select(starts_with(column_prefix)) %>%
        rowMeans(na.rm = TRUE) %>%
        as.numeric()
}

# Calculate average for each sample
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

# Calculate average for parents
data <- data %>%
    mutate(
        Parent_avg = rowMeans(select(., Parent, `Parent (2nd)`), na.rm = TRUE)
    )

# Remove original quadruplet columns and 'Parent (2nd)' column
data <- data %>%
    select(Gene, ends_with("_avg"), Parent_avg)

# Check the structure of the dataframe
print(str(data))

# Remove zero values from all columns except the first one (gene names)
remove_zeros <- function(df) {
    for (col in 2:ncol(df)) {
        df <- df[df[, col] != 0, ]
    }
    return(df)
}

# Apply the function to remove zero values
processed_df <- remove_zeros(data)

# Check the processed dataframe
print(head(processed_df))

# Function to bin data
bin_data <- function(data) {
    bins <- seq(0, max(data, na.rm = TRUE) + 1, by = 1)
    binned_data <- cut(data, breaks = bins, right = FALSE, include.lowest = TRUE)
    return(binned_data)
}

# Function to calculate JS distance
calc_js <- function(data1, data2) {
    binned_data1 <- bin_data(data1)
    binned_data2 <- bin_data(data2)
    freq_data1 <- table(binned_data1)
    freq_data2 <- table(binned_data2)
    prob_data1 <- freq_data1 / sum(freq_data1)
    prob_data2 <- freq_data2 / sum(freq_data2)
    
    common_bins <- intersect(names(prob_data1), names(prob_data2))
    prob_data1 <- prob_data1[common_bins]
    prob_data2 <- prob_data2[common_bins]
    
    epsilon <- 1e-10
    prob_data1 <- prob_data1 + epsilon
    prob_data2 <- prob_data2 + epsilon
    
    relative_entropy1 <- sum(prob_data1 * log2(prob_data1 / prob_data2), na.rm = TRUE)
    relative_entropy2 <- sum(prob_data2 * log2(prob_data2 / prob_data1), na.rm = TRUE)
    js_dist <- (relative_entropy1 + relative_entropy2) / 2
    return(js_dist)
}

# Calculate JS distances between all pairs of samples
samples <- names(processed_df)[-1]
js_matrix <- matrix(0, nrow = length(samples), ncol = length(samples),
                    dimnames = list(samples, samples))

for (i in 1:length(samples)) {
    for (j in i:length(samples)) {
        js_dist <- calc_js(processed_df[[samples[i]]], processed_df[[samples[j]]])
        js_matrix[i, j] <- js_dist
        js_matrix[j, i] <- js_dist  # JS distance is symmetric
    }
}

# Convert the JS distances matrix to a data frame
js_distances_df <- melt(js_matrix, varnames = c("Sample1", "Sample2"), value.name = "JS_Distance")

# Apply -log10 transformation to JS Distance values
js_distances_df <- js_distances_df %>%
    mutate(JS_Distance = -log10(JS_Distance))

# Print the JS distances dataframe
print(js_distances_df)

# Save the JS distances to a CSV file
write.csv(js_distances_df, "js_distances_between_samples_log10.csv", row.names = FALSE)

# Create a heatmap of JS distances with values in scientific notation
ggplot(js_distances_df, aes(x = Sample1, y = Sample2, fill = JS_Distance)) +
    geom_tile() +
    geom_text(aes(label = JS_Distance, digits = 2), color = "white", size = 3) +
    scale_fill_gradient2(low = "blue", high = "red", midpoint = median(js_distances_df$JS_Distance)) +
    labs(title = "Heatmap of JS Distances Between Samples",
         x = "Sample",
         y = "Sample") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the heatmap to a file
ggsave("js_distances_heatmap_log10.png", width = 10, height = 8)

cat("JS distances between all samples have been calculated and saved to 'js_distances_between_samples.csv' and the heatmap has been saved to 'js_distances_heatmap.png'.\n")