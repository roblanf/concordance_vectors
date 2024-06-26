#!/usr/bin/env Rscript

# Load necessary libraries
library(tidyverse)
library(ggplot2)
library(boot)

bootstrap_ci_counts <- function(counts, R = 1000, alpha = 0.05) {
    # Convert counts to a data frame for bootstrapping
    data <- tibble(category = rep(1:length(counts), round(counts)))
    
    # Bootstrap function to compute counts for each category
    bootstrap_counts <- function(data, indices) {
        resample <- data[indices, ]
        resampled_counts <- resample %>%
            group_by(category) %>%
            summarise(count = n(), .groups = 'drop')
        
        # Ensure all categories are present in the resampled counts
        resampled_counts <- resampled_counts %>%
            complete(category = 1:length(counts), fill = list(count = 0))
        
        return(resampled_counts$count)
    }
    
    # Perform bootstrap resampling
    results <- boot(data = data, statistic = bootstrap_counts, R = R)
    
    # Calculate 95% confidence intervals for each category
    ci_results <- tibble(
        lower_ci_count = apply(results$t, 2, quantile, probs = alpha / 2),
        upper_ci_count = apply(results$t, 2, quantile, probs = 1 - alpha / 2)
    )
    
    return(ci_results)
}

# Define the create_heatmap function with bootstrap CIs
create_heatmap <- function(branch_id, conc_vectors) {
    # Filter and select data based on branch ID
    clade_data <- conc_vectors %>%
        filter(ID == branch_id) %>%
        select(gene_psi1, gene_psi2, gene_psi3, gene_psi4, site_psi1, site_psi2, site_psi3, site_psi4, quartet_psi1, quartet_psi2, quartet_psi3, quartet_psi4,
               gene_psi1_N, gene_psi2_N, gene_psi3_N, gene_psi4_N, site_psi1_N, site_psi2_N, site_psi3_N, site_psi4_N, quartet_psi1_N, quartet_psi2_N, quartet_psi3_N, quartet_psi4_N)
    
    # Reshape data from wide to long format for proportions
    long_data <- clade_data %>%
        select(-ends_with("_N")) %>%
        pivot_longer(cols = everything(), names_to = "type_psi", values_to = "value") %>%
        separate(type_psi, into = c("type", "psi"), sep = "_psi") %>%
        mutate(type = factor(type, levels = c("quartet", "site", "gene")))
    
    # Reshape data from wide to long format for counts
    long_data_N <- clade_data %>%
        select(ends_with("_N")) %>%
        pivot_longer(cols = everything(), names_to = "type_psi_N", values_to = "count") %>%
        separate(type_psi_N, into = c("type", "psi_N"), sep = "_psi") %>%
        mutate(type = factor(type, levels = c("quartet", "site", "gene")))
    
    # Calculate bootstrap 95% CIs on the concordance vectors
    long_data_N <- long_data_N %>%
        group_by(type) %>%
        group_modify(~ {
            counts <- .x$count
            total_counts <- sum(counts)
            ci <- bootstrap_ci_counts(counts)
            ci <- ci %>%
                mutate(lower_CI = (lower_ci_count / total_counts) * 100,
                       upper_CI = (upper_ci_count / total_counts) * 100)
            .x <- .x %>% bind_cols(ci)
            return(.x)
        }) %>%
        ungroup()
    
    # Add bootstrap CIs to the long_data table in percentage terms
    long_data <- long_data %>%
        mutate(psi = paste0(psi, "_N")) %>%
        left_join(long_data_N %>% select(type, psi_N, lower_CI, upper_CI), by = c("type", "psi" = "psi_N")) %>%
        mutate(psi = str_replace(psi, "_N", ""),
               main_label = scales::number(value, accuracy = 0.001),
               ci_label = paste0("(", scales::number(lower_CI, accuracy = 0.001), ", ", 
                                 scales::number(upper_CI, accuracy = 0.001), ")"))
    
    # Generate the heatmap
    plot <- ggplot(long_data, aes(x = psi, y = type, fill = value)) +
        geom_tile(color = "white") + # Add tiles with white borders
        geom_text(aes(label = main_label), color = "black", size = 6, vjust = -0.5) + # Add main number text labels
        geom_text(aes(label = ci_label), color = "black", size = 3, vjust = 1.5) + # Add CI text labels
        scale_fill_gradient(low = "white", high = "red", limits = c(0, 100), name = "proportion") + # Color gradient
        scale_x_discrete(position = 'top', labels = c(bquote(Psi[1]), bquote(Psi[2]), bquote(Psi[3]), bquote(Psi[4]))) +
        theme_minimal() + # Minimal theme
        ggtitle(paste("Concordance factors for branch ID", branch_id), subtitle = "Values are percentages, numbers in brackets are bootstrap 95% CIs") +
        theme(plot.title = element_text(size = 24),
              axis.title = element_blank(), # Remove axis titles
              axis.text = element_text(size = 16), # Increase axis text size
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              line = element_blank(),
              legend.position = "bottom") + # Place legend at the bottom
        guides(fill = guide_colourbar(
            title = "Concordance or Discordance factor (%)",
            title.position = "top",
            title.hjust = 0.5,
            barwidth = unit(8, "cm"),
            barheight = unit(0.5, "cm"),
            label.position = "bottom",
            label.theme = element_text(size = 14),
            ticks.colour = "black"
        ))
    
    output_data <- long_data %>%
        select(type, psi, value, lower_CI, upper_CI)
    
    return(list(plot = plot, data = output_data))
}

# Main script
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Please provide both the concordance vectors CSV file and a branch ID as command line arguments.")
}

csv_file <- args[1]
branch_id <- as.numeric(args[2])

# Load the input data
concordance_vectors <- read_csv(csv_file)

output <- create_heatmap(branch_id, concordance_vectors)

# Save the plot to a PDF file
pdf_filename <- paste0("concordance_table_", branch_id, ".pdf")
ggsave(pdf_filename, plot = output$plot, width = 8, height = 5)

# Save the long_data table to a CSV file
csv_filename <- paste0("concordance_table_", branch_id, ".csv")
write_csv(output$data, csv_filename)
