# EasyHPV
EasyHPV is a pipeline detecting HPV genotypes and human genome integration breakpoints. 

### Requirements
Kraken2 

Samtools

bwa

picard

pandas

python

pyranges

pysam

### Steps
1. Map fastq file to human reference
2. Map fastq file to HPVs using Kraken
3. Generate HPV genotype result file 
4. Extract reads that are both mapped to human and HPV references
5. Filter true integrated HPV reads
6. Generate HPV integration result file

## Options
| Options | Description |
|---------|-------------|
|-r1, --R1|Read 1 fastq file|
|-r2, --R2|Read 2 fastq file|
|-o, --output|Output name|
|-r, --reference|Human reference fasta file|
|-hr, --hpv-reference|HPV reference fasta file directory|
|-k, --kraken2-directory|Kraken2 tool directory|
|-b, --bed_file|Bed file for annotation|
|-t, --threads|Threads|
|-m, --mapping-only|If only mapping to human reference step `-m y`|
|-g, --genotyping-only|If only genotyping by Kraken2 step `-g y`|
|-i, --integration-only|If only integration step `-i y`|


### Example Run Code
`python3 easyhpv.py -r1 {fastq1} -r2 {fastq2} -o {output name} -r {human refence file} -hr {HPV reference fasta file directory} -k {Kraken2 tool directory} -t {threads} -b {bed file}`

### Example Bed file 
Require following columns
|chr|start|end|NM|gene|type|strand|exon_num|
|--|---|---|--|----|----|------|-------|
|chr1|65419|65433|NM_001005484|OR4F5|exon|+|exon_1|


