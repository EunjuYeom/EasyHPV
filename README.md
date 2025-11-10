# EasyHPV
EasyHPV is a pipeline detecting HPV genotypes and human genome integration breakpoints. 

## Installation
Install EasyHPV codes and files with following code

`git clone https://github.com/EunjuYeom/EasyHPV.git`

then, install required tools using environment.yaml file

`conda env update --file environment.yaml`

If you are not using conda, install the following tools manually. 
*Kraken2 is required to be manually installed even when using conda environment.
### Requirements
[kraken2](https://github.com/DerrickWood/kraken2)

[samtools](https://github.com/samtools/samtools)

[bwa](https://github.com/lh3/bwa)

[picard](https://github.com/broadinstitute/picard)

pandas

python

pyranges

pysam

## Pipeline workflow
<img width="988" height="220" alt="image" src="https://github.com/user-attachments/assets/a4095d00-2cea-4071-bd06-49107f620232" />


1. Map fastq file to human reference
2. Map fastq file to HPVs using Kraken
3. Generate HPV genotype result file 
4. Extract reads that are both mapped to human and HPV references
5. Filter true integrated HPV reads
6. Generate HPV integration result file

## Options
| Options | Description |
|---------|-------------|
|-r1, --R1|Read 1 fastq file (must be gzipped)|
|-r2, --R2|Read 2 fastq file (must be gzipped)|
|-o, --output|Output name|
|-r, --reference|Human reference fasta file|
|-hr, --hpv-reference|HPV reference fasta file directory|
|-k, --kraken2-directory|Kraken2 tool directory|
|-b, --bed_file|Bed file for annotation|
|-t, --threads|Threads|
|-m, --mapping-only|If only mapping to human reference step `-m y` (default: n)|
|-g, --genotyping-only|If only genotyping by Kraken2 step `-g y` (default: n)|
|-i, --integration-only|If only integration step `-i y` (default: n)|


### Example Run Code
`python3 easyhpv.py -r1 {fastq1} -r2 {fastq2} -o {output name} -r {human refence file} -hr {HPV reference fasta file directory} -k {Kraken2 tool directory} -t {threads} -b {bed file}`

### Example Bed file 
Require following columns
|chr|start|end|NM|gene|type|strand|exon_num|
|--|---|---|--|----|----|------|-------|
|chr1|65419|65433|NM_001005484|OR4F5|exon|+|exon_1|


