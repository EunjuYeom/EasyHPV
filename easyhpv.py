import os
import pandas as pd
import argparse
import sys
import re
import warnings
from datetime import datetime

    

def mapping(file1, file2, name, ref):

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
    
    os.makedirs(name, exist_ok=True)
    os.chdir(name)    
    
    cmd4 = f"bwa mem {ref} {file1} {file2} > {name}.sam"
    cmd5 = f"samtools view -Sb {name}.sam > {name}.bam"
    cmd6 = f"samtools sort {name}.bam -o {name}_sorted.bam"
    cmd7 = f"picard MarkDuplicates I={name}_sorted.bam O={name}_sorted_dedup.bam M={name}_dup_metrics.txt REMOVE_DUPLICATES=true"
    cmd8 = f"samtools index {name}_sorted_dedup.bam"
    
    os.system(cmd4)
    os.system(cmd5)
    os.system(cmd6)
    os.system(cmd7)
    os.system(cmd8)
    
    

#----------------------------------------------------------------------------------------#
if __name__ == '__main__' : 
    path = os.getcwd()
    warnings.simplefilter(action='ignore', category=FutureWarning)
    parser = argparse.ArgumentParser(description='EasyHPV Usage')
    parser.add_argument("-r1", "--R1", dest = "R1", action = "store") 
    parser.add_argument("-r2", "--R2", dest = "R2", action = "store") 
    parser.add_argument("-o","--output", dest = "output", action = "store", help = "Output name")
    parser.add_argument("-r","--reference", dest = "ref", action = "store", help = "Reference fasta file directory")
    parser.add_argument("-m", "--mapping-only", dest="mapping", action="store", choices=['y', 'n'], default='n')
    parser.add_argument("-k", "--kraken-only", dest="kraken", action="store", choices=['y', 'n'], default='n')
    parser.add_argument("-f", "--filter-only", dest="filter_only", action="store", choices=['y', 'n'], default='n')
    args = parser.parse_args()

    if args.mapping == 'y' :
    
        mapping(args.R1, args.R2, args.output, args.ref)
        
    elif args.kraken == 'y' : 
        
        kraken(args.R1, args.output)

    elif args.filter_only == 'y' :
        
        filter_only(args.R1, args.output)
        
    else :
        run(args.R1, args.output)
        
#python3 /labmed/99.YEJ/EasyHPV/easyhpv.py -r1 

#samtools, bwa, picard, pandas, python

