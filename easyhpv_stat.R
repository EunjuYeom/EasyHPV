
library(dplyr)
library(writexl)
library(stringr)
library(tidyr)
library(argparse)

parser <- ArgumentParser(description="HPV integration analysis sample")

parser$add_argument("-n", "--name", help="Sample name", required=TRUE)

args <- parser$parse_args()


Name <- args$name


data <- read.table(sprintf("%s_result_true_filtered_annot_50M.txt", Name, Name), sep="\t", header=T)

integrated_read_count_df <- data %>%
  distinct(read,HPV,gene, .keep_all = TRUE) %>%
  group_by(HPV, gene) %>%
  summarise(integrated_read = n(), .groups = "drop") %>%
  arrange(desc(integrated_read))

write.table(
  integrated_read_count_df, 
  file = sprintf("%s_filtered_HPV_gene_counts_50M.txt", Name, Name),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
  )
  
integrated_read_count_df_2 <- data %>%
  distinct(read,HPV,gene,exon_num,type, .keep_all = TRUE) %>%
  group_by(HPV, gene, exon_num, type) %>% 
  summarise(integrated_read_2 = n(), .groups = "drop") %>%
  arrange(desc(integrated_read_2))
  
write.table(
  integrated_read_count_df_2,
  file = sprintf("%s_filtered_HPV_pos_counts_50M.txt", Name, Name),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
  )
  