
library(dplyr)
library(writexl)
library(stringr)
library(tidyr)
library(argparse)

parser <- ArgumentParser(description="HPV integration analysis sample")

parser$add_argument("-n", "--name", help="Sample name", required=TRUE)

args <- parser$parse_args()


Name <- args$name

data <- read.table(sprintf("%s_results_read_count.txt", Name), sep="\t", header=F)
colnames(data) <- c("read", "chr", "pos", "CIGAR", "MAPQ","HPV", "F_supporting_kmer", "R_supporting_kmer", "num_of_reads")
data <- data %>% mutate(`T/F` = "T")

data <- data %>%
  filter(!str_detect(HPV, "^Alpha|^Beta|^Gamma|^First")) %>%
  filter(!str_detect(HPV, "(taxid 10566)|(taxid 173087)"))

data <- data %>%
  mutate(`T/F` = ifelse(str_detect(CIGAR, "H"), "F", `T/F`)) %>%
  mutate(`T/F` = ifelse(F_supporting_kmer <= 20 & R_supporting_kmer < 30, "F", `T/F`)) %>%
  mutate(`T/F` = ifelse(F_supporting_kmer < 30 & R_supporting_kmer <= 20, "F", `T/F`)) %>%
  group_by(read) %>%
  mutate(`T/F` = ifelse(any(`T/F` == "F"), "F", `T/F`)) %>%
  ungroup() %>%
  mutate(`T/F` = ifelse(str_detect(CIGAR, "\\d+S\\d+M\\d+S"), "F", `T/F`)) %>% 
  mutate(`T/F` = ifelse(num_of_reads == 1 & CIGAR == '*', "F", `T/F`)) %>%
  rowwise() %>%
  mutate(`T/F` = ifelse(sum(as.numeric(unlist(str_extract_all(CIGAR, "[0-9]+(?=M)")))) < 30, "F", `T/F`))%>%
  mutate(`T/F` = ifelse(MAPQ <= 50,"F",`T/F`)) %>%
  mutate(`T/F` = ifelse(num_of_reads == 1 & sum(as.numeric(unlist(str_extract_all(CIGAR, "[0-9]+(?=M)")))) > 100, "F", `T/F`)) %>%
  ungroup()
  
  

data_true <- data %>%
  filter(`T/F` == "T")
  
data_true_final <- data_true %>%
  filter(!str_detect(CIGAR, "^\\d+S\\d+M\\d+S$")) %>%
  filter(CIGAR != "*")

breakpoint_cal <- function(pos,CIGAR) {
  sm = str_extract_all(CIGAR, "\\d+[SM]")[[1]]
  for (element in sm) {
    n_num = as.numeric(str_extract(element, "\\d+"))
    n_type = str_extract(element, "[SM]")
    if (n_num >= 5) {
      if (n_type == "S") return (pos-1)
      if (n_type == "M") return (pos+n_num)
    }
  }
  return(NA)
}

breakpoint_cal_2 <- function(br_group) {
  if(length(unique(br_group$breakpoint)) > 1) {
    S_included <- br_group %>% filter(str_detect(CIGAR, "S"))
    if (nrow(S_included) > 0) {
      final_br <- S_included$breakpoint[1]
      br_group$breakpoint <- final_br
    }
  }
  return(br_group)
}

data_true_final_br <- data_true_final %>%
  rowwise() %>%
  mutate(breakpoint=breakpoint_cal(pos,CIGAR)) %>%
  ungroup() %>%
  group_by(read) %>%
  group_modify(~ breakpoint_cal_2(.x)) %>%
  ungroup()
    
    


write.table(data,
            file = sprintf("%s_TF_result_50M.txt", Name),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
  

write.table(data_true_final_br,
            file = sprintf("%s_result_true_filtered_50M.txt", Name),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

