# Code for Inbreeding of Exmoor ponies

---
* Title: Code for Inbreeding of Exmoor ponies
* Author: Alexandra Brumwell
---

Code notebook for honours project on Exmoor pony populations genetics

--------------------------------------------------------------------------------------------
## Introduction
This aims to look at inbreeding and mutational load of Exmoor ponies

### Downloading a reference genome
I will be using a complete horse genome the EquCab3 reference (Genbank ID: NC_009144.3).


## FASTQ to BAM pipeline

### 1. Mapping FASTQs
```bash
cat /shared5/Alex/Inbreeding_project/List_Sep_2022_FASTQ.txt | awk '{print $2}' | while read file ; do

  # Defining variables:
  location=$(echo ${file} | awk -F "/" '{$NF=""}1' OFS="/") ; # variable of file path
  name=$(echo ${file} | awk -F "/" '{print $NF}' | sed 's/_1_val_1.fq.gz//'); # variable of sample name

  # Mapping
  bwa mem /shared5/Alex/Inbreeding_project/horse_WG_reference/GCF_002863925.1_EquCab3.0_genomic.fna \ # reference genome in fasta format
  ${location}${name}_1_val_1.fq.gz \ # read 1
  ${location}${name}_2_val_2.fq.gz \ # read 2
  2> /shared5/Alex/Inbreeding_project/BAMs/log_files/${name}.bwamem.log \ # creates a log file containing errors
  > /shared5/Alex/Inbreeding_project/BAMs/${name}_horseref.sam &  # output file location and name
done
```

### 1.1 Validating SAM Files
```bash
cat /shared5/Alex/Inbreeding_project/List_Sep_2022_FASTQ.txt | awk '{print $1}' | while read name ; do
  java -jar /shared5/Alex/picard-2.27.5/picard.jar ValidateSamFile -I /shared5/Alex/Inbreeding_project/BAMs/${name}*sam -O  /shared5/Alex/Inbreeding_project/BAMs/${name}_sam_horseref_validate.txt --MODE SUMMARY -R /shared5/Alex/Inbreeding_project/horse_WG_reference/GCF_002863925.1_EquCab3.0_genomic.fna &
done
```



### 2. Converting to BAM and sorting
Now that the SAM files have been created, I convert them to BAM format so they take less space and sort them by reading name using `samtools view` and then `samtools sort`:

```bash
cat /shared5/Alex/Inbreeding_project/List_Sep_2022_FASTQ.txt | awk '{print $1}' | while read name ; do
  samtools view -h -b -u /shared5/Alex/Inbreeding_project/BAMs/${name}_horseref.sam |  samtools sort -o /shared5/Alex/Inbreeding_project/BAMs/${name}_sorted_horseref.bam &
done
```



### 3. Marking duplicates
REDO WITH SAMTOOLS
I use Picardtools MarkDuplicates as it removes interchromsomal duplicates and overall described as better:
```bash
cat /shared5/Alex/Inbreeding_project/List_Sep_2022_FASTQ.txt | awk '{print $1}' | while read name ; do
  java -Xmx4g -jar /shared5/Alex/picard-2.27.5/picard.jar MarkDuplicates  \
    -I /shared5/Alex/Inbreeding_project/BAMs/${name}_sorted_horseref.bam \ # input file
    -O /shared5/Alex/Inbreeding_project/BAMs/${name}_rmdp_sorted_horseref.bam \ #output file
    -M /shared5/Alex/Inbreeding_project/BAMs/log_files/${name}_markdup_metrics.txt \ # text file with marked_dup_metrics
    --REMOVE_DUPLICATES true --ASSUME_SORT_ORDER coordinate --VALIDATION_STRINGENCY SILENT & \
done
```
