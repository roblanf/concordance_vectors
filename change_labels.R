# Load necessary libraries
library(ape)
library(tidyverse)

# Read the tree file
treeo <- read.tree("gcf.cf.branch")
tree <- treeo

# Read the CSV file
concordance_data <- read_csv("concordance_vectors.csv")

# Prepare the label information
concordance_data <- concordance_data %>%
    mutate(new_label = paste0(ID, ":g", gene_psi1, ",s", site_psi1, ",q", quartet_psi1))

# Create a lookup table for the new labels
label_lookup <- setNames(concordance_data$new_label, as.character(concordance_data$ID))

# Assign new labels to node labels
tree$node.label <- sapply(tree$node.label, function(x) {
    new_label <- label_lookup[as.character(x)]
    if(is.na(new_label)) return(x) else return(new_label)
})

# fix the names
names(tree$node.label) = treeo$node.label

write.nexus(tree, file = "id_gcf_scf_qcf.nex")

