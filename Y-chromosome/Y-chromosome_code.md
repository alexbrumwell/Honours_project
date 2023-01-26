# Code for Y-chromosomal analysis of Exmoor ponies

---
* Title: Code for Y-chromosomal analyses
* Author: Alexandra Brumwell
---

Code notebook for honours project on Exmoor pony populations genetics

--------------------------------------------------------------------------------------------
## Introduction
I will be studying the Y-chromosome in male Exmoor ponies and horses. I will be mapping these sequnces to a reference horse Y-chromosome. The sequences I  have are: 22 Exmoor sequences, 34 horse sequences and 4 Przewalski sequences.

Prior to starting any data processing or analyses I create a spreadsheet containing the metadata of the sequences like genome ID, date of sequencing, breed, file name, sex, etc. I save all this information in a [spreadsheet](ADD LINK TO SPREADSHEET IN GOOGLEDRIVE) in Excel or Googledrive.

##### Example of metadata table:

| Individual | File_name | Breed | Pony number | Pony name | Sequencing_date | Sex | Notes | Original_FASTQ_file_path |
|------------|-----------|-------|-------------|-----------|-----------------|-----|-------|--------------------------|
| ERR111 | ERR111_EYB1000 | Exmoor | 111/111 | Cheddar | Mar2020 | F | From unknown maternal line | /Path/to/FASTQ/file |



## Downloading/obtaining a reference genome
We need to have a reference genome to map the FASTQs on to.
There are several different Y-chromosome assemblies.

The first I will be using is a non-repetitive Male specific region from Wallner et al., 2017. It is 1.6 Mbp in total.
To download it, I use NCBI and search in Genbank for its ID (PRJNA353919). Then I open the associated FTP directory and copy the link for the ".fna.gz" file. In the cluster, I use `wget -P /chosen_directory/ <copied_link>` to download the assembly directly into my directory. This reference file is called "GCA_002166905.2_LipY764_genomic.fna".

I also want to map to a number of specific Y-chromosome genes, in order to avoid repetitive areas. So I read the list of genes they used in Felhel et al., 2019: AMELY, SRY, TSPY, USP9Y, UTY, ZFY. Check table (file:///C:/Users/asus/Downloads/03011table1.html) from Paria et al., 2011 to find Genbank IDs to download the genes. --> In the end, Barbara said that this was unecessary as there is little to no variation in horse Ychromosome genes.



## Creation of sample list
First I check the pedigree of the indivudals that each sequence is from to see which are male. Then I make a list of all the male Exmoor, horse and Przewalski sequences and save this in a text file called *List_ychromosome.txt*. I add the file path for read 1 of each FASTQ as a correspoding second column

```bash
# Example of list:
S21029_EDSW200015548-1a_HJFWVDSXY_L4    /shared5/Alex/Exmoor_sequencing_data/Sep_2020/S21029_EDSW200015548-1a_HJFWVDSXY_L4_1_val_1.fq.gz
S21092_EDSW200015546-1a_HJFWVDSXY_L4    /shared5/Alex/Exmoor_sequencing_data/Sep_2020/S21092_EDSW200015546-1a_HJFWVDSXY_L4_1_val_1.fq.gz
S21109_EDSW200015550-1a_HJFWVDSXY_L4    /shared5/Alex/Exmoor_sequencing_data/Sep_2020/S21109_EDSW200015550-1a_HJFWVDSXY_L4_1_val_1.fq.gz
```

## FASTQtoBAM pipeline
### 1. Mapping FASTQs
**Important to note**: The FASTQs I am working with are already trimmed.

Sidenote: Some papers (Wallner et al., 2017; Felkhel et al., 2019) on horse Y-chromosome use bwa aln but i think it's because they have many reptetitive regions, so I will continue using bwa mem.

I map using `bwa mem` with the `-M` option that "Marks shorter split hits as secondary (for Picard compatibility)". It esentially marks uncertain reads as secondary.

```bash
cat /shared5/Alex/Y-chromosome_project/List_ychromosome.txt | awk '{print $2}' | while read file ; do

  # Defining variables:
  location=$(echo ${file} | awk -F "/" '{$NF=""}1' OFS="/") ; # variable of file path
  name=$(echo ${file} | awk -F "/" '{print $NF}' | sed 's/_1_val_1.fq.gz//'); # variable of sample name

  # Mapping
  bwa mem -M /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/GCA_002166905.2_LipY764_genomic.fna \ # reference genome in fasta format
  ${location}${name}_1_val_1.fq.gz \ # read 1
  ${location}${name}_2_val_2.fq.gz \ # read 2
  2> /shared5/Alex/Y-chromosome_project/BAMs_LipY764/log_files/${name}.bwamem.log
  > /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764.sam &  # output file location and name
done
```
***Note on code formatting:*** *\"\\" at the end of a line means that the command continues in the following line of code. In other words, the five lines of code in the `bwa mem` command are usually written in one line but, in this case, it is more visual to see each argument on a different line*



## 2. Converting to BAM and sorting
Now that the SAM files have been created, I convert them to BAM format so they take less space and sort them by reading name using `samtools view` and then `samtools sort`:

I am using options based on Felkel et al. (2019).

```bash
cat /shared5/Alex/Y-chromosome_project/List_ychromosome.txt | awk '{print $1}' | while read name ; do
  samtools view -h -b -F 4 -u /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764.sam | \ # "-q" Skip alignments with MAPQ
  samtools sort -o /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted.bam & # this creates a log file containing standard errors
done
```
> `samtools view` options:
* `-h`: *Include the header in the output*
* `-b`: *Output in BAM format*
* `-F`: *Exclusion flag*. In this case the flag is set at "4" so it will exclude unmapped reads.
* `-u`: *Output uncompressed data*. Preferred when piping to another samtools command.

> `samtools sort`options:
* `-o`: *Write the final sorted output to specified file name*


## 3. Marking duplicates
I use Picardtools MarkDuplicates as it removes interchromsomal duplicates and overall desribed as better:
```bash
cat /shared5/Alex/Y-chromosome_project/List_ychromosome.txt | awk '{print $1}' | while read name ; do
  java -Xmx4g jar /shared5/Alex/picard-2.27.5/picard.jar MarkDuplicates  \
    -I /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted.bam \ # input file
    -O /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted_rmdp.bam \ #output file
    -M /shared5/Alex/Y-chromosome_project/BAMs_LipY764/log_files/${name}_markdup_metrics.txt \ # text file with marked_dup_metrics
    --REMOVE_DUPLICATES true --ASSUME_SORT_ORDER coordinate --VALIDATION_STRINGENCY SILENT & \
done
```

## 4. Merging BAM files
Since some of the sequences I have belong to the same individual, I have to merge the corresponding BAM files together.

`samtools merge` can be used for this:
```bash
samtools merge E_23279_merged_LipY764_sorted_rmdp.bam E_23279_EKDN220034363-1A_H5J5MDSX5_L2_LipY764_sorted_rmdp.bam E_23279_EKDN220034363-1A_H7VKMDSX5_L3_LipY764_sorted_rmdp.bam;

samtools merge E_23416_merged_LipY764_sorted_rmdp.bam E_23416_EKDN220034368-1A_H5J5MDSX5_L2_LipY764_sorted_rmdp.bam E_23416_EKDN220034368-1A_H7VKMDSX5_L1_LipY764_sorted_rmdp.bam E_23416_EKDN220034368-1A_H7VL2DSX5_L2_LipY764_sorted_rmdp.bam;

samtools merge E_320005_merged_LipY764_sorted_rmdp.bam E_320005_EKDN220034373-1A_H5J5MDSX5_L2_LipY764_sorted_rmdp.bam E_320005_EKDN220034373-1A_H7VKMDSX5_L1_LipY764_sorted_rmdp.bam E_320005_EKDN220034373-1A_H7VL2DSX5_L4_LipY764_sorted_rmdp.bam;

samtools merge E_32023_merged_LipY764_sorted_rmdp.bam E_32023_EKDN220034372-1A_H5J5MDSX5_L2_LipY764_sorted_rmdp.bam E_32023_EKDN220034372-1A_H7VKMDSX5_L1_LipY764_sorted_rmdp.bam
```
And I create a new list (List_Y-chromosome_merged.txt) with the updated IDs of the sequences. E.g. The merged file IDs have now become "s479021_merged" instead of "s2479021_EDSW210003783-1a_H3WNKDSX2_L4"

## 5. Quality filtering
I use quality filters form Felkhel et al. (2019):

```bash
cat /shared5/Alex/Y-chromosome_project/List_Y-chromosome_merged.txt | awk '{print $1}' | while read name ; do
  samtools view -h -b -F 4 -q 20 /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted_rmdp.bam -o /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted_rmdp_MQ20.bam &
done
```

## 6. Indexing
Now that the BAMs has been created, sorted and duplicates have been marked, they need to be indexed for future analyses.

```bash
cat /shared5/Alex/Y-chromosome_project/List_Y-chromosome_merged.txt | awk '{print $1}' | while read name ; do
  samtools index /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted_rmdp_MQ20.bam &
done
```

## 7. Validating BAM files
It is good practice to make sure that the BAMs have been created correctly.

There are two tools to check the validity of a BAM file: `samtools quickcheck` and `ValidateSamFile` from Picardtools. Samtools is very quick but not that thorough whereas Picard will report all errors encountered when checking a BAM file. As described by GATK, "*ValidateSamFile is most useful for troubleshooting errors encountered with other tools that may be caused by improper formatting, faulty alignments, incorrect flag values, etc.*

Generally I use `ValidateSamFile` on the final BAMs I have created and `samtools quickcheck` on the files that I have created in previous steps (i.e. SAMs or the sorted BAMs).

```bash
#Samtools quickecheck:
samtools quickcheck <input.bam>

#ValidateSamFile:
java -jar /shared5/Alex/picard-2.27.5/picard.jar ValidateSamFile \
  -I <input.bam> \ # input file
  -O <output.textfile> \ # output textfile containing any errors
  -R <reference_sequence> \ # optional
  --MODE SUMMARY # outputs a summary table listing the numbers of all 'errors' and 'warnings'.
```

## 8. Calculating BAM coverage
I use `qualimap bamqc` to calculate the coverage of the BAM files [(qualimap bamqc manual)](http://qualimap.conesalab.org/doc_html/analysis.html#bamqc).

*Important to note:* There is no need to specify an output as qualimap creates an output folder containing the stats results.

```bash
cat /shared5/Alex/Y-chromosome_project/List_Y-chromosome_merged.txt | awk '{print $1}' | while read name ; do
  qualimap bamqc -bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted_rmdp_MQ20.bam --java-mem-size=200G 2>/shared5/Alex/Y-chromosome_project/BAMs_LipY764/log_files/${name}_qualimap.log & # Sets desired memory size
done
```

## 9.1 Prelimnary consensus FASTA creation

```bash
#creating depth list
cat /shared5/Alex/Y-chromosome_project/List_Y-chromosome_merged.txt | awk '{print $1}' | while read name; do
  coverage=$(cat /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}*stats/genome_results.txt | grep "mean coverageData"  | awk '{print $4}' | sed 's/X//' | sed 's/,//' | cut -d "." -f 1);
  paste <(echo ${name}) <(echo ${coverage}) <(echo "${coverage} * 2" | bc);
done > /shared5/Alex/Y-chromosome_project/Y-chromosome_ID_depth_list_new.txt

#consensus FASTAs
cat /shared5/Alex/Y-chromosome_project/Y-chromosome_ID_depth_list_new.txt | while read line; do
  #Defining variables:
  name=$(echo ${line} | awk '{print $1}');
  number=$(echo ${line} | awk '{print $(NF-1)}');
  maxnumber=$(echo ${line} | awk '{print $(NF)}');
  #consensus FASTA
  angsd -i /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted_rmdp_MQ20.bam -minMapQ 30 -minQ 20 -setMinDepth ${number} -setMaxDepth ${maxnumber} -remove_bads 1 -doFasta 2 -doCounts 1 -ref /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/GCA_002166905.2_LipY764_genomic.fna -out /shared5/Alex/Y-chromosome_project/FASTAs_LipY764/${name}_LipY764 &
done

#first we unzip the compressed FASTAs
gunzip *fa.gz

#change headers
cat /shared5/Alex/Mitochondrial_project/Mitochondrial_ID_depth_list_horseref.txt | awk '{print $1}' | while read name; do
  sed -i "s/>NC_001640.1/>${name}_>NC_001640.1/g" /shared5/Alex/Mitochondrial_project/FASTAs_horseref/${name}_mtdna_horseref.fa
done
```



## 9. Variant calling
My original pipeline used ANGSD to call consensus FASTAs but most papers use GATK Haplotype caller.
Felkhel et al., (2019) does:
```bash
java -jar GenomeAnalysisTK.jar -T HaplotypeCaller –I sample_mapped_sort_rmdup_MQ20.bam
-R Lip_Y_Assembly_reapr.fa > sample.g.vcf
java -jar GenomeAnalysisTK.jar -T CombineGVCFs -V sample1.g.vcf -V sample2.g.vcf
–V sampleX.g.vcf -o cohort.g.vcf
```
In Wallner et al. (2017) they use vcftools for variant calling:
"    Variants for were called with Samtools.

    samtools-0.1.18/bcf-tools

    samtools-0.1.18/samtools view -h filename_trimmed_mapped_sorted_nodupSAMTminusS_MQ20.bam | samtools view -bS - | samtools mpileup -A -u -f reference.fa - | bcftools view -vcg - | vcfutils.pl varFilter -d 2 -D 200 > filename_trimmed_only_mapped_sorted_nodupSAMTminusS_MQ20.vcf"


Whilst the Canid Y-chromosome paper did: "We next recalled Y-chromosome genotypes in the male samples using the GATK haplotype caller using the EMIT_ALL_SITES flag. Following methods previously used for identifying repetitive sequence unsuitable for read mapping in primate Y-chromosomes, we used mapping and depth statistics from the VCF info field to identify callable regions on the canine Y-chromosome (Fig. 1) [26, 27]. Once a callable region was identified, we applied site level filtering to the remaining sequence: dropping maximum likelihood heterozygotes, missing sites, and positions with an MQ0/raw depth ratio > 0.10. A second depth filter was then applied to remove positions with extreme sequencing depths (median depth ± 3 M.A.D.). We also removed positions that were within 5 bp of GATK called indels."
