# Code for 5 largest Exmoor sequences

---
* Title: Code for WG Exmoor analysis
* Author: Alexandra Brumwell
---

Code notebook for honours project on Exmoor pony populations genetics

--------------------------------------------------------------------------------------------
## Introduction
This little sideproject consisted of mapping the whole genome of the 5 largest Exmoor sequences I had from the neww batch of samples.

1.	First I create a list with the 5 largest genomes.
```bash
ll -h /shared5/Alex/Exmoor_sequencing_data/Sep_2022/01.RawData/*/*1.fq.gz | cut -d " " -f 10 > Exmoor_top5_list.txt

# Exmoor_list.txt content:
/shared5/Alex/Exmoor_sequencing_data/Sep_2022/01.RawData/E_102004/E_102004_EKDN220034353-1A_H5GL5DSX5_L2_1.fq.gz
/shared5/Alex/Exmoor_sequencing_data/Sep_2022/01.RawData/E_107013/E_107013_EKDN220034352-1A_H5GL5DSX5_L2_1.fq.gz
```

2.	Trimming FASTQs
```bash
cat /shared5/Alex/Exmoor_sequencing_data/Sep_2022/Exmoor_top5_list.txt) | while read file; do
  #defining variables
  name=$(echo $file| cut -d "/" -f 7 | cut -d "_" -f 1-5); location=$(echo $file | cut -d "/" -f 1-6);

  #trimming
  trim_galore –-paired ${location}/${name}_1.fq.gz \ #read 1
  ${location}/${name}_2.fq.gz \ #read 2
  -o /shared5/Alex/new_Exmoor_genomes/FASTQs_trimmed &
done
```

3.	Summary statistics of FASTQs
A quick way to see the size of the FASTQs by order is:
```bash
cat /shared5/Alex/Exmoor_sequencing_data/Sep_2022/Exmoor_top5_list.txt) | while read file; do
  name=$(echo $file| cut -d "/" -f 7 | cut -d "_" -f 1-5);
  ll -Sh ${name}*_1.fq.gz | cut -d " " -f 5 ;
done
```
Although the proper way is to count the number of reads per FASTQ by using `zcat` and then dividing by 4:
```bash
echo "$(zcat <file>) | wc -l) / 4 | bc"

#same but looping over files
cat /shared5/Alex/Exmoor_sequencing_data/Sep_2022/Exmoor_list.txt | while read file; do
  paste <(echo $file) <( echo $(zcat $file | wc -l)/4|bc) ;
done >> fastq_reads.txt
```

4. Mapping
```bash
cat /shared5/Alex/Exmoor_sequencing_data/Sep_2022/Exmoor_top5_list.txt | while read name; do
  bwa mem /shared5/Alex/Donkey_ref/GCF_016077325.2_ASM1607732v2_genomic.fna.gz /shared5/Alex/Exmoor_sequencing_data/Sep_2022/FASTQs_trimmed/${name}_1_val_1.fq.gz /shared5/Alex/Exmoor_sequencing_data/Sep_2022/FASTQs_trimmed/${name}_2_val_2.fq.gz > /shared5/Alex/Exmoor_sequencing_data/Sep_2022/BAM_mappeddonkey/${name}.sam &
done
```

5. Converting to BAM and sorting by read name
```bash
cat /shared5/Alex/Exmoor_sequencing_data/Sep_2022/Exmoor_top5_list.txt | while read name; do
  samtools view -q 30 -u /shared5/Alex/Exmoor_sequencing_data/Sep_2022/BAM_mappeddonkey/${name}.sam | samtools sort -n -o /shared5/Alex/Exmoor_sequencing_data/Sep_2022/BAM_mappeddonkey/${name}_sorted.bam &
done
```

6. Markduplicates
```bash
cat /shared5/Alex/Exmoor_sequencing_data/Sep_2022/Exmoor_top5_list.txt | while read name; do
  samtools fixmate -m /shared5/Alex/Exmoor_sequencing_data/Sep_2022/BAM_mappeddonkey /${name}_sorted.bam - | samtools sort - | samtools markdup - /shared5/Alex/Exmoor_sequencing_data/Sep_2022/BAM_mappeddonkey/${name}_sorted_markdup.bam &
done
```

7. Indexing
```bash
for file in *markdup*; do
  samtools index $file
done
```
8. Calculate BAM depth
```bash
cat /shared5/Alex/Exmoor_sequencing_data/Sep_2022/Exmoor_top5_list.txt | while read name; do
  qualimap bamqc -bam /shared5/Alex/Exmoor_sequencing_data/Sep_2022/BAM_mappeddonkey/${name}_sorted_markdup.bam --java-mem-size=200G &
done
```
