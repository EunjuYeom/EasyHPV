
library(dplyr)
library(writexl)
library(stringr)
library(tidyr)
library(argparse)

parser <- ArgumentParser(description="HPV integration analysis sample")

parser$add_argument("-n", "--name", help="Sample name", required=TRUE)

args <- parser$parse_args()


Name <- args$name

data <- read.table(sprintf("%s_results_read_count.txt", Name, Name), sep="\t", header=F)
colnames(data) <- c("read", "chr", "pos", "CIGAR", "MAPQ","HPV", "F_supporting_kmer", "R_supporting_kmer", "num_of_reads")
data <- data %>% mutate(`T/F` = "T")

data <- data %>%
  filter(!str_detect(HPV, "^Alpha|^Beta|^Gamma|^First")) %>%
  filter(!str_detect(HPV, "(taxid 10566)|(taxid 173087)"))

data <- data %>%
  mutate(`T/F` = ifelse(str_detect(CIGAR, "H"), "F", `T/F`)) %>%
  mutate(`T/F` = ifelse(`F_supporting_kmer` <= 20 & `R_supporting_kmer` < 30, "F", `T/F`)) %>%
  mutate(`T/F` = ifelse(`F_supporting_kmer` < 30 & `R_supporting_kmer` <= 20, "F", `T/F`)) %>%
  mutate(`T/F` = ifelse(`num_of_reads` == 1 & `CIGAR` == '*', "F", `T/F`)) %>%
  group_by(read) %>%
  mutate(`T/F` = ifelse(any(`T/F` == "F"), "F", `T/F`)) %>%
  ungroup() %>%
  #mutate(`T/F` = ifelse((`T/F` != "F"), "T", `T/F`)) %>%
  mutate(`T/F` = ifelse(str_detect(CIGAR, "\\d+S\\d+M\\d+S"), "F", `T/F`)) %>% 
  rowwise() %>%
  mutate(`T/F` = ifelse(sum(as.numeric(unlist(str_extract_all(CIGAR, "[0-9]+(?=M)")))) < 40, "F", `T/F`))%>%
  mutate(`T/F` = ifelse(`MAPQ` <= 50,"F",`T/F`)) %>%
  mutate(`T/F` = ifelse(`num_of_reads` == 1 & sum(as.numeric(unlist(str_extract_all(CIGAR, "[0-9]+(?=M)")))) > 100, "F", `T/F`)) %>%
  ungroup()
  
  

data_true <- data %>%
  filter(`T/F` == "T")
  
data_true_final <- data_true %>%
  filter(!str_detect(CIGAR, "^\\d+S\\d+M\\d+S$")) %>%
  filter(CIGAR != "*")
  
write.table(data,
            file = sprintf("%s_TF_result_50M.txt", Name, Name),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
  

write.table(data_true_final,
            file = sprintf("%s_result_true_filtered_50M.txt", Name, Name),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

