import os
import pandas as pd
import pyranges as pr
import argparse
import sys
import re
import warnings
import pysam
from datetime import datetime

    

def mapping(file1, file2, name, ref, threads, path):
  
    current_dir = path
    
    #os.path.basename() : extract file name
    ref_fa = os.path.basename(ref)
    ref_name = os.path.basename(ref).split(".fa")[0]
    ref_dir = os.path.dirname(ref)
    
    index_file = os.path.join(ref_dir, f"{ref_fa}.bwt")
    
    if os.path.exists(index_file) :
      
      pass
    
    else :
      
      cmd1 = f"bwa index {ref}"
      
      os.system(cmd1)
    
    sample_dir = os.path.join(current_dir, name)
    os.makedirs(sample_dir, exist_ok=True)
    os.chdir(sample_dir)    
    
    cmd4 = f"bwa mem -t {threads} {ref} {file1} {file2} > {name}.sam"
    cmd5 = f"samtools view -Sb {name}.sam > {name}.bam"
    cmd6 = f"samtools sort {name}.bam -o {name}_sorted.bam"
    cmd7 = f"picard MarkDuplicates I={name}_sorted.bam O={name}_sorted_dedup.bam M={name}_dup_metrics.txt REMOVE_DUPLICATES=true"
    cmd8 = f"samtools index {name}_sorted_dedup.bam"
    cmd9 = f"samtools flagstat {name}_sorted_dedup.bam > {name}_flagstat.txt" 
    
    os.system(cmd4)
    os.system(cmd5)
    os.system(cmd6)
    os.system(cmd7)
    os.system(cmd8)
    os.system(cmd9)
    
#----------------------------------------------------------------------------------------#
def genotyping(file1, file2, name, HPV_ref, kraken2_dir, threads, path): 

    current_dir = path
    
    db_dir = os.path.dirname(os.path.normpath(HPV_ref))
    final_file = os.path.join(db_dir, "HPV.txt")
       
    os.chdir(db_dir) 
    
    if not os.path.exists(final_file) :
      
      #taxonomy file 
      cmd10 = f"{kraken2_dir}kraken2-build --download-taxonomy --db {db_dir} --threads {threads}"
      os.system(cmd10)
    
      #library file 
      for file in os.listdir(HPV_ref) :
        if file.endswith(".fa"):
          HPV_file = os.path.join(HPV_ref,file)
          HPV_name = os.path.basename(file).split(".fa")[0]
          cmd11 = f"{kraken2_dir}kraken2-build --add-to-library {HPV_file} -db {db_dir} --threads {threads}"
          os.system(cmd11)
        
      
      #library build & inspect .k2d 
      cmd12 = f"{kraken2_dir}kraken2-build --build --db {db_dir} --threads {threads}"
      cmd13 = f"{kraken2_dir}kraken2-inspect --threads {threads} --db {db_dir} > {db_dir}/HPV.txt"
        
      os.system(cmd12)
      os.system(cmd13)
    
    os.chdir(current_dir)
    sample_dir = os.path.join(current_dir, name)
    os.makedirs(sample_dir, exist_ok=True)
    os.chdir(sample_dir)
    
    #run kraken2
    cmd14=f"{kraken2_dir}kraken2 --use-names --db {db_dir} --threads {threads} \
            --paired {file1} {file2} > {name}_kraken2.kraken \
            --report {name}_kraken2.report --gzip-compressed"
    os.system(cmd14)

    #HPV genotype file 
    df_1=pd.read_csv(f"{name}_kraken2.report", sep="\t", header = None, names=["perc", "read_1", "read_2", "classification","taxid", "species"])
    
    df_1["species"]=df_1["species"].str.strip()
    
    HPV_rows = []
    
    for i, row in df_1.iterrows() :
      names = row["species"]
      if ("human" in names or "Human" in names) and (row["read_2"] !=0) and (row["classification"] != "S") and (row["classification"].startswith("S")) :
        HPV_rows.append({"HPV" : names, "Reads" : row["read_2"]})
      
    df_2 = pd.DataFrame(HPV_rows)
    
    with open(f"{name}_flagstat.txt", "r") as flagstat:
      line = flagstat.readlines()
      total_mapped_reads = int(line[6].split("+")[0])
    
    df_2["Total_mapped_reads"] = total_mapped_reads
    df_2["Normalized_reads"] = df_2["Reads"]/df_2["Total_mapped_reads"]*(10**8)
    
    df_2["HPV_positive"] = "F"
    df_2.loc[df_2["Normalized_reads"] >= 100.9, "HPV_positive"] = "T"
                  
    df_2.to_csv(f"{name}_kraken2_genotype_results.csv", sep = "\t", index = False)


def integration(file1, file2, name, HPV_ref, annot_bed, path): 

    current_dir = path

    db_dir = os.path.dirname(HPV_ref)
    easyhpv_path = os.path.join(db_dir,"../../")
    easyhpv_path = os.path.abspath(easyhpv_path)
    Int_R_path = os.path.join(easyhpv_path,"easyhpv_int.R")
    
    
    sample_dir = os.path.join(current_dir, name)
    os.makedirs(sample_dir, exist_ok=True)
    os.chdir(sample_dir)

    #extract kraken classified reads
    with open(f"{name}_kraken2.kraken", "r") as kraken :
      with open(f"{name}_classified_read_extract.csv","w") as result:  
        for line in kraken :
          category=line.strip().split("\t")
          classified = category[0]
          read = category[1]
          classification = category[2] 
          length = category[3] 
          info = category[4]
          
          joined="\t".join([classified, read, classification, length, info])
          
          if classified == "C" :
            result.write(joined + "\n")
    
    #find and filter integration 
    included_chromo=["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"] 

    read_list=set()
  
    with open(f"{name}_classified_read_extract.csv", "r") as extracted_reads:
      for line in extracted_reads :
        read_info=line.split()[1]
        read_list.add(read_info)
    
    
    
    with pysam.AlignmentFile(f"{name}_sorted_dedup.bam", "rb") as bam_file, open(f"{name}_int_reads.txt", "w") as out_file:
      for read in bam_file :
        if read.query_name in read_list : 
          chromo = read.reference_name
          if chromo in included_chromo :
            out_file.write(read.to_string()+"\n")
                 
      
    
    #add kmer 
    df_3 = pd.read_csv(f"{name}_int_reads.txt", sep="\t", header=None, usecols=range(6), dtype={2:str})
  
    df_3.columns = ["read","flag","chr","pos","MAPQ","CIGAR"]
  
    kmer_df = pd.read_csv(f"{name}_classified_read_extract.csv", sep = "\t", header = None, names = ["classification", "read", "HPV", "length", "kmer"])
    
    merged_df = df_3.merge(kmer_df[["read", "HPV", "kmer"]], on = "read")
    
    def process_kmer(row):
      taxid = row["HPV"].split("taxid ")[1].strip(")")
      kmer_counts = {"F_total" : 0, "R_total" : 0}
      
      if "|:|" in row["kmer"]:
        F_kmers, R_kmers = row["kmer"].split(" |:| ")
        kmer_counts["F_total"] = sum(int(x.split(":")[1]) for x in F_kmers.split() if x.startswith(taxid))
        kmer_counts["R_total"] = sum(int(x.split(":")[1]) for x in R_kmers.split() if x.startswith(taxid))
        
      else : 
        F_kmers = row["kmer"]
        kmer_counts["F_total"] = sum(int(x.split(":")[1]) for x in F_kmers.split() if x.startswith(taxid))
        
      return pd.Series(kmer_counts)
      
    merged_df[["F_supporting_kmer", "R_supporting_kmer"]] = merged_df.apply(process_kmer, axis=1)
    
    merged_df = merged_df[["read", "chr", "pos", "CIGAR","MAPQ","HPV","F_supporting_kmer", "R_supporting_kmer"]]
    
    merged_df.to_csv(f"{name}_results.txt", sep="\t", index=False, header=None)
    
    
    #add_num_reads
    df_4=pd.read_csv(f"{name}_results.txt", sep="\t", header = None, names=["read", "chr", "pos", "CIGAR","MAPQ", "HPV", "F_supporting_kmer", "R_supporting_kmer"])
    
    read_counts = df_4["read"].value_counts().rename("num_of_reads")
    
    df_5 = df_4.merge(read_counts, left_on = "read", right_index = True)
    
    df_5.to_csv(f"{name}_results_read_count.txt", sep = "\t", index = False, header = False)
    
    
    #true_false
    cmd15 = f"Rscript {easyhpv_path}/easyhpv_int.R -n {name}"
    os.system(cmd15)
    
    #annotation 
    
    df_6 = pd.read_csv(f"{name}_result_true_filtered_50M.txt", sep="\t",dtype= {"read" : str, "chr":"category","pos": int, "CIGAR": str, "MAPQ": int, "HPV": "category", "F_supporting_kmer":int,"R_supporting_kmer" : int, "num_of_reads":int, "T/F" : "category", "breakpoint":int})
    
    annot_ref = pd.read_csv(f"{annot_bed}", sep="\t", dtype ={"chr":"category", "start":int, "end":int, "NM":str,"gene":"category", "type":"category", "strand":"category", "exon_num":"category"})
    
    df_6["Start"] = df_6["breakpoint"]
    df_6["End"] = df_6["breakpoint"]+1
    pr_df = pr.PyRanges(df_6.rename(columns={"chr":"Chromosome"}))
    
    annot_ref=annot_ref.rename(columns={"chr" : "Chromosome", "start" : "Start", "end" : "End"})
    pr_annot = pr.PyRanges(annot_ref)
    
    joined=pr_df.join(pr_annot)
    
    
    joined_df = joined.df
    joined_df["Position"] = joined_df["Chromosome"].astype(str)+":"+joined_df["Start"].astype(str)
    
    joined_final = joined_df[["read", "Position", "HPV", "CIGAR","MAPQ", "gene", "exon_num", "F_supporting_kmer", "R_supporting_kmer", "type"]]
    
    final_df = joined_final.groupby(["read", "Position","HPV","gene","exon_num"])
    
    final_result=[]
    
    for key, group in final_df :
      combined_type = ",".join(sorted(set(group["type"])))
      
      row = group.iloc[0][["read", "Position", "HPV", "CIGAR", "MAPQ", "gene", "exon_num", "F_supporting_kmer", "R_supporting_kmer"]].to_dict()
      
      row["type"]=combined_type
      final_result.append(row)
      
      
    final_annot_df = pd.DataFrame(final_result)
    final_annot_df = final_annot_df.sort_values(by = "Position", ascending = True).reset_index(drop=True) 
    final_annot_df.to_csv(f"{name}_result_true_filtered_annot_50M.txt", sep = "\t", index = False, header = True)
    
    #stat
    cmd16 = f"Rscript {easyhpv_path}/easyhpv_stat.R -n {name}"
    os.system(cmd16) 

def run(file1, file2, name, ref, HPV_ref, kraken2_dir, threads, annot_bed, path):
    
    mapping(file1, file2, name, ref, threads, path)
    genotyping(file1, file2, name, HPV_ref, kraken2_dir, threads, path)
    integration(file1, file2, name, HPV_ref, annot_bed, path) 
#----------------------------------------------------------------------------------------#
if __name__ == '__main__' : 
    path = os.getcwd()
    warnings.simplefilter(action='ignore', category=FutureWarning)
    parser = argparse.ArgumentParser(description='EasyHPV Usage')
    parser.add_argument("-r1", "--R1", dest = "R1", action = "store") 
    parser.add_argument("-r2", "--R2", dest = "R2", action = "store") 
    parser.add_argument("-o","--output", dest = "output", action = "store", help = "Output name")
    parser.add_argument("-r","--reference", dest = "ref", action = "store", help = "Human reference fasta file directory")
    parser.add_argument("-hr", "--hpv-reference", dest = "hpv_reference", action = "store", help = "HPV reference fasta file directory")
    parser.add_argument("-k", "--kraken2-directory", dest = "kraken2_directory", action = "store", help = "Kraken2 tool directory")
    parser.add_argument("-m", "--mapping-only", dest="mapping", action="store", choices=['y', 'n'], default='n')
    parser.add_argument("-g", "--genotyping-only", dest="genotyping", action="store", choices=['y', 'n'], default='n')
    parser.add_argument("-i", "--integration-only", dest="integration_only", action="store", choices=['y', 'n'], default='n')
    parser.add_argument("-b", "--bed_file", dest="annotation_bed", action="store", help = "Bed file for annotation")
    parser.add_argument("-t", "--threads", dest="threads", action="store", default='2')
    args = parser.parse_args()

    if args.mapping == 'y' :
    
        mapping(args.R1, args.R2, args.output, args.ref, args.threads, path)
        
    elif args.genotyping == 'y' : 
        
        genotyping(args.R1, args.R2, args.output, args.hpv_reference, args.kraken2_directory, args.threads, path)

    elif args.integration_only == 'y' :
        
        integration(args.R1,args.R2, args.output, args.hpv_reference, args.annotation_bed, path)
        
    else :
        run(args.R1, args.R2, args.output, args.ref, args.hpv_reference, args.kraken2_directory, args.threads, args.annotation_bed, path)
        



