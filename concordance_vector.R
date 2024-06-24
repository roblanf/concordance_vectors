library(tidyverse)

# Read the data files
scf <- read_table("scfl.cf.stat", comment = "#")
gcf <- read_table("gcf.cf.stat", comment = "#")
cbl <- read_table("coalescent_bl.cf.stat", comment = "#")

# Process the scf file
scv <- scf %>%
  rename(site_psi1 = sCF, 
         site_psi1_N = sCF_N, 
         site_N = sN,
         length_subs_per_site = Length) %>%
  mutate(site_psi4 = 0, 
         site_psi4_N = 0,
         quartet_psi4 = 0, 
         quartet_psi4_N = 0) %>%
  rowwise() %>%
  mutate(site_psi2 = max(sDF1, sDF2),
         site_psi3 = min(sDF1, sDF2),
         site_psi2_N = max(sDF1_N, sDF2_N),
         site_psi3_N = min(sDF1_N, sDF2_N)) %>%
  ungroup() %>%
  mutate(
    quartet_psi1 = round(as.numeric(str_extract(Label, "(?<=q1=)[^;]+")), 2),
    quartet_psi2 = round(pmax(as.numeric(str_extract(Label, "(?<=q2=)[^;]+")),
                              as.numeric(str_extract(Label, "(?<=q3=)[^;]+"))), 2),
    quartet_psi3 = round(pmin(as.numeric(str_extract(Label, "(?<=q2=)[^;]+")),
                              as.numeric(str_extract(Label, "(?<=q3=)[^;]+"))), 2),
    quartet_psi1_N = round(as.numeric(str_extract(Label, "(?<=f1=)[^;]+")), 2),
    quartet_psi2_N = round(pmax(as.numeric(str_extract(Label, "(?<=f2=)[^;]+")),
                                as.numeric(str_extract(Label, "(?<=f3=)[^;]+"))), 2),
    quartet_psi3_N = round(pmin(as.numeric(str_extract(Label, "(?<=f2=)[^;]+")),
                                as.numeric(str_extract(Label, "(?<=f3=)[^;]+"))), 2),
    quartet_psi1_pp = round(as.numeric(str_extract(Label, "(?<=pp1=)[^;]+")), 1),
    quartet_psi2_pp = round(pmax(as.numeric(str_extract(Label, "(?<=pp2=)[^;]+")),
                                 as.numeric(str_extract(Label, "(?<=pp3=)[^;]+"))), 1),
    quartet_psi3_pp = round(pmin(as.numeric(str_extract(Label, "(?<=pp2=)[^;]+")),
                                 as.numeric(str_extract(Label, "(?<=pp3=)[^;]+"))), 1),
    quartet_N = round(as.numeric(str_extract(Label, "(?<=EN=)[^\\]]+")), 1)
  ) %>%
  select(ID, site_psi1, site_psi2, site_psi3, site_psi4, site_psi1_N, site_psi2_N, site_psi3_N, site_psi4_N, site_N, 
        quartet_psi1, quartet_psi2, quartet_psi3, quartet_psi4, quartet_psi1_N, quartet_psi2_N, quartet_psi3_N, quartet_psi4_N, quartet_psi1_pp, quartet_psi2_pp, quartet_psi3_pp, quartet_N, length_subs_per_site)

# Process the gcf file
gcv <- gcf %>%
  rename(gene_psi1 = gCF, 
         gene_psi1_N = gCF_N, 
         gene_N = gN,
         gene_psi4 = gDFP,
         gene_psi4_N = gDFP_N) %>%
  rowwise() %>%
  mutate(gene_psi2 = max(gDF1, gDF2),
         gene_psi3 = min(gDF1, gDF2),
         gene_psi2_N = max(gDF1_N, gDF2_N),
         gene_psi3_N = min(gDF1_N, gDF2_N)) %>%
  ungroup() %>%
  select(ID, gene_psi1, gene_psi2, gene_psi3, gene_psi4, gene_psi1_N, gene_psi2_N, gene_psi3_N, gene_psi4_N, gene_N)

# Process the cbl file
cbl <- cbl %>%
  rename(length_coalescent = Length) %>%
  select(ID, length_coalescent)

# Combine the processed data
concordance_vectors <- scv %>%
  left_join(gcv, by = "ID") %>%
  left_join(cbl, by = "ID") %>%
  select(ID, gene_psi1, gene_psi2, gene_psi3, gene_psi4, gene_psi1_N, gene_psi2_N, gene_psi3_N, gene_psi4_N, gene_N,
         site_psi1, site_psi2, site_psi3, site_psi4, site_psi1_N, site_psi2_N, site_psi3_N, site_psi4_N, site_N,
         quartet_psi1, quartet_psi2, quartet_psi3, quartet_psi4, quartet_psi1_N, quartet_psi2_N, quartet_psi3_N, quartet_psi4_N, quartet_N, quartet_psi1_pp, quartet_psi2_pp, quartet_psi3_pp,
         length_subs_per_site, length_coalescent)

# Write to CSV
write_csv(concordance_vectors, "concordance_vectors.csv")
