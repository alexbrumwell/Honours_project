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

Then to see how much coverage each genome has:
```bash
#new ref
cat /shared5/Alex/Y-chromosome_project/List_Y-chromosome_merged.txt | awk '{print $1}' | while read name; do
  paste <(echo $name) <(grep  "mean coverageData" /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}*stats/*txt) <(grep  "number of reads" /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}*stats/*txt)
done

#old ref
cat /shared5/Alex/Y-chromosome_project/List_Y-chromosome_merged.txt | awk '{print $1}' | while read name; do
  paste <(echo $name) <(grep  "mean coverageData" /shared5/Alex/Y-chromosome_project/BAMs/${name}*stats/*txt) <(grep  "number of reads" /shared5/Alex/Y-chromosome_project/BAMs/${name}*stats/*txt)
done
```



#### 9.0.1 Prelimnary consensus FASTA creation
Barbara says to try running it with the depth threshold set at 5x
```bash
#creating depth list
cat /shared5/Alex/Y-chromosome_project/List_Y-chromosome_merged.txt | awk '{print $1}' | while read name; do
  coverage=$(cat /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}*stats/genome_results.txt | grep "mean coverageData"  | awk '{print $4}' | sed 's/X//' | sed 's/,//' | cut -d "." -f 1);
  paste <(echo ${name}) <(echo ${coverage}) <(echo "${coverage} * 2" | bc);
done > /shared5/Alex/Y-chromosome_project/Y-chromosome_ID_depth_list_new.txt

#consensus FASTAs
cat /shared5/Alex/Y-chromosome_project/Y-chromosome_ID_depth_list_new.txt | while read line; do
  name=$(echo ${line} | awk '{print $1}');
  maxnumber=$(echo ${line} | awk '{print $(NF-1)}');
  angsd -i /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted_rmdp_MQ20.bam -minMapQ 30 -minQ 20 -setMinDepth 5 -setMaxDepth ${maxnumber} -remove_bads 1 -doFasta 2 -doCounts 1 -ref /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/GCA_002166905.2_LipY764_genomic.fna -out /shared5/Alex/Y-chromosome_project/FASTAs_LipY764/${name}_5x_LipY764_Y &
done

#first we unzip the compressed FASTAs
gunzip *fa.gz

#because there are many scaffolds, I remove the headings of the scaffolds
cat /shared5/Alex/Y-chromosome_project/Y-chromosome_ID_depth_list.txt | awk '{print $1}' | while read name; do
   sed -i '/>/d' /shared5/Alex/Y-chromosome_project/FASTAs_LipY764/${name}_5x_LipY764_Y.fa
done

#change headers
cat /shared5/Alex/Y-chromosome_project/Y-chromosome_ID_depth_list.txt | while read line; do
  name=$(echo $line | awk '{print $1}');
  breed=$(echo $line | awk '{print $NF}');
  echo ">${name}_${breed}_LipY764" | cat - /shared5/Alex/Y-chromosome_project/FASTAs_LipY764/${name}_5x_LipY764_Y.fa > tmp.fa && mv tmp.fa /shared5/Alex/Y-chromosome_project/FASTAs_LipY764/${name}_5x_LipY764_Y.fa
done
```

#### 9.0.3 Concatenating FASTAs
Then we compile all the FASTAs in one file and add the reference mitogenome as well:
```bash
#concatinating FASTAs
cat /shared5/Alex/Y-chromosome_project/Y-chromosome_ID_depth_list.txt | awk '{print $1}' | while read name; do
  cat /shared5/Alex/Y-chromosome_project/FASTAs_LipY764/${name}_5x_LipY764_Y.fa
done >> Exmoor_horse_5x_horseref_NEW_Y.fa


## adding reference sequence
cat /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/GCA_002166905.2_LipY764_genomic.fna Exmoor_horse_5x_horseref_NEW_Y.fa > tmp.fa && mv tmp.fa Exmoor_horse_5x_horseref_NEW_Y.fa
```

#### 9.0.4 Copy concatinated FASTA file to local disk
```bash
scp studentprojects@young.eng.gla.ac.uk:/shared5/Alex/Y-chromosome_project/FASTAs/Exmoor_horse_5x_horseref_NEW_Y.fa ./Desktop
```




## 9. Variant calling with Bcftools

#Anubhabs pipeline:
source activate popgen
bcftools mpileup -f REF BAM1 BAM2 BAM3  > name.pileup
bcftools call -c -O v --ploidy-file ploidy.txt -o name.vcf name.pileup
#ploIDY.TXT contains chromosomeY 1 ENDPOS M 1

```bash
#To compile list of BAM files:
cat /shared5/Alex/Y-chromosome_project/Y-chromosome_ID_depth_list_new.txt | cut -f 1 | while read name; do
  echo "/shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted_rmdp_MQ20.bam ";
done > tmp.txt
echo $(cat tmp.txt)

#bcftools mpileup:
bcftools mpileup -f /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/MH341179_Ychromosome.fasta /shared5/Alex/Y-chromosome_project/BAMs/E_102004_EKDN220034353-1A_H5GL5DSX5_L2_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/E_107013_EKDN220034352-1A_H5GL5DSX5_L2_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/E_21084_EKDN220034367-1A_H5J5MDSX5_L2_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/E_23279_merged_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/E_23416_merged_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/E_320005_merged_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/E_32023_merged_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/E_49031_EKDN220034370-1A_H5J5MDSX5_L2_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/E_519005_EKDN220034364-1A_H5J5MDSX5_L2_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/s2010_EDSW210003772-1a_H3WNKDSX2_L3_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/s2012_EDSW210003766-1a_H3WHWDSX2_L4_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/s479013_EDSW210003775-1a_H3WNKDSX2_L3_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/s49097_EDSW210003769-1a_H3FTWDSX2_L2_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/s49124_EDSW210003767-1a_H3WHWDSX2_L2_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/S276023_EDSW200015547-1a2a_HJFWVDSXY_L3L4_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR1527950_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR1527947_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR1305964_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR1545178_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR1545179_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR978597_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR978599_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR978603_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR2179545_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR2179555_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR1545180_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR2179549_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR2731060_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR1545181_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR2179547_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR1545190_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR863167_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR2731057_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR1527967_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR2179546_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR2179552_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR2179548_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR1545184_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR1545183_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR2179544_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR1527969_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR1527970_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR1545188_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR2179556_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR1545189_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR1735862_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR1545185_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR1545187_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/ERR1545186_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/SRR12719743_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/SRR12719745_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/SRR12719757_Y_sorted_markdup.bam /shared5/Alex/Y-chromosome_project/BAMs/SRR12719758_Y_sorted_markdup.bam > /shared5/Alex/Y-chromosome_project/variants_bcftools/Exmoors_2022_Y_oldref.pileup

#created ploidy file. I get the number from the third column by adding up the bases of the contigsfrom the fai index file
chromosomeY 1 6462487 M 1 > ploidy.txt

#bcftools call:
bcftools call -c -O v --ploidy-file /shared5/Alex/Y-chromosome_project/variants_bcftools/ploidy.txt -o /shared5/Alex/Y-chromosome_project/variants_bcftools/Exmoors_2022_Y_oldref.vcf /shared5/Alex/Y-chromosome_project/variants_bcftools/Exmoors_2022_Y_oldref.pileup
```

I do same thing but with new reference too:
```bash
#bcftools mpileup
bcftools mpileup -f /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/GCA_002166905.2_LipY764_genomic.fna /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_102004_EKDN220034353-1A_H5GL5DSX5_L2_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_107013_EKDN220034352-1A_H5GL5DSX5_L2_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_21084_EKDN220034367-1A_H5J5MDSX5_L2_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_23279_merged_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_23416_merged_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_320005_merged_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_32023_merged_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_49031_EKDN220034370-1A_H5J5MDSX5_L2_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_519005_EKDN220034364-1A_H5J5MDSX5_L2_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/s2010_EDSW210003772-1a_H3WNKDSX2_L3_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/s2012_EDSW210003766-1a_H3WHWDSX2_L4_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/s479013_EDSW210003775-1a_H3WNKDSX2_L3_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/s49097_EDSW210003769-1a_H3FTWDSX2_L2_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/s49124_EDSW210003767-1a_H3WHWDSX2_L2_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/S276023_EDSW200015547-1a2a_HJFWVDSXY_L3L4_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR1527950_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR1527947_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR1305964_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR1545178_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR1545179_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR978597_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR978599_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR978603_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR2179545_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR2179555_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR1545180_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR2179549_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR2731060_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR1545181_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR2179547_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR1545190_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR863167_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR2731057_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR1527967_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR2179546_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR2179552_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR2179548_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR1545184_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR1545183_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR2179544_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR1527969_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR1527970_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR1545188_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR2179556_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR1545189_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR1735862_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR1545185_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR1545187_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/ERR1545186_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/SRR12719743_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/SRR12719745_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/SRR12719757_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/SRR12719758_LipY764_sorted_rmdp_MQ20.bam > /shared5/Alex/Y-chromosome_project/variants_bcftools/Exmoors_2022_Y_newref.pileup

#bcftools call
bcftools call -c -O v --ploidy-file /shared5/Alex/Y-chromosome_project/variants_bcftools/ploidy.txt -o /shared5/Alex/Y-chromosome_project/variants_bcftools/Exmoors_2022_Y_newref.vcf /shared5/Alex/Y-chromosome_project/variants_bcftools/Exmoors_2022_Y_newref.pileup
```


# 10. VCF filtering
1.1 First step, quality filter and removing INDELS:
```bash
#oldref
vcftools --vcf /shared5/Alex/Y-chromosome_project/variants_bcftools/Exmoors_2022_Y_oldref.vcf --remove-indels --minQ 30 --minGQ 30 --out /shared5/Alex/Y-chromosome_project/variants_bcftools/Exmoors_2022_Y_oldref_rmvIndels_minQ30_minGQ30 --recode

#new reference
vcftools --vcf /shared5/Alex/Y-chromosome_project/variants_bcftools/Exmoors_2022_Y_newref.vcf --remove-indels --minQ 30 --minGQ 30 --out /shared5/Alex/Y-chromosome_project/variants_bcftools/Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30 --recode
```

1.2 To see how much missing information there is per individual:
`vcftools –gzvcf data.vcf.gz –missing-indv –out data`

2. Minimum allele count filter
vcftools --vcf Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30.recode.vcf --min-alleles 2 --mac 3 --out Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30_biallelic_mac3 --recode

3. Exclude sites where >10% of data is missing
vcftools --vcf Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30_biallelic_mac3.recode.vcf --max-missing 0.9 --out Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30_biallelic_mac3.recode.vcf

#To see missingness of individuals:
```bash vcftools --vcf Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30_biallelic_mac3.recode.vcf --missing-indv --out Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30_biallelic_mac3.recode.vcf
```

# 11. PLINK
## 11.1 Download PLINK
```bash
conda create -n plink -c bioconda plink
source activate plink
```

## 11.2 Create PLINK file
Create plink file:
```bash
vcftools --vcf Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30_biallelic_mac3.recode.vcf --plink --out Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30_biallelic_mac3
```

## 11.3 PCA
```bash
plink --file Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30_biallelic_mac3 --pca --out Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30_biallelic_mac3
```








## Other attempts to call variants
# Strelka
Anubhabs suggests I use Strelka.

I configure the configure strelka file first:
```bash
cat /shared5/Alex/Y-chromosome_project/Y-chromosome_ID_depth_list_new.txt | cut -f 1 | while read name; do
  echo "--bam /shared5/Alex/Y-chromosome_project/BAMs/${name}_Y_sorted_markdup.bam \\";
done

#output of configure_strelka_old:
configureStrelkaGermlineWorkflow.py \
--bam /shared5/Alex/Y-chromosome_project/BAMs/E_102004_EKDN220034353-1A_H5GL5DSX5_L2_Y_sorted_markdup.bam \
--bam /shared5/Alex/Y-chromosome_project/BAMs/E_107013_EKDN220034352-1A_H5GL5DSX5_L2_Y_sorted_markdup.bam \
--bam /shared5/Alex/Y-chromosome_project/BAMs/s479013_EDSW210003775-1a_H3WNKDSX2_L3_Y_sorted_markdup.bam \
(...)
--referenceFasta /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/MH341179_Ychromosome.fasta \
--runDir . \
--ploidy /shared5/Alex/Y-chromosome_project/variants/head_ploidy.vcf.gz


#for the ploidy file I created a header containing the sample IDs and other info and then:
bgzip -c head_ploidy.vcf > head_ploidy.vcf.gz
tabix -p vcf head_ploidy.vcf.gz

#run configure script:
bash configure_strelka.sh

#run runWorkflow
./runWorkflow.py -m local -j 10
```

Strelka won't run; getting error.


# GATK Haplotypecaller
```bash
#CODE FOR GATK
cat /shared5/Alex/Y-chromosome_project/List_Y-chromosome_merged.txt | awk '{print $1}' | while read name; do
  gatk /home/opt/miniconda2/pkgs/gatk-3.8-5/bin/GenomeAnalysisTK -R /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/GCA_002166905.2_LipY764_genomic.fna –I /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted_rmdp_MQ20.bam -O /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764.vcf
done

cat /shared5/Alex/Y-chromosome_project/List_Y-chromosome_merged.txt | awk '{print $1}' | while read name; do
  java -jar /shared3/Anubhab/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T HaplotypeCaller -R /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/GCA_002166905.2_LipY764_genomic.fna –I /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted_rmdp_MQ20.bam -O /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764.vcf
done

java -jar /shared3/Anubhab/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T CombineGVCFs -V sample1.g.vcf -V sample2.g.vcf
–V sampleX.g.vcf -o cohort.g.vcf
```




My original pipeline used ANGSD to call consensus FASTAs but most papers use GATK Haplotype caller.
Felkel et al., (2019) does:
```bash
java -jar GenomeAnalysisTK.jar -T HaplotypeCaller –I sample_mapped_sort_rmdup_MQ20.bam
-R Lip_Y_Assembly_reapr.fa > sample.g.vcf
java -jar GenomeAnalysisTK.jar -T CombineGVCFs -V sample1.g.vcf -V sample2.g.vcf
–V sampleX.g.vcf -o cohort.g.vcf
```


In Wallner et al. (2017) they use vcftools for variant calling:
"Variants for were called with Samtools.

    samtools-0.1.18/bcf-tools

    samtools-0.1.18/samtools view -h filename_trimmed_mapped_sorted_nodupSAMTminusS_MQ20.bam | samtools view -bS - | samtools mpileup -A -u -f reference.fa - | bcftools view -vcg - | vcfutils.pl varFilter -d 2 -D 200 > filename_trimmed_only_mapped_sorted_nodupSAMTminusS_MQ20.vcf"


Whilst the Canid Y-chromosome paper did: "We next recalled Y-chromosome genotypes in the male samples using the GATK haplotype caller using the EMIT_ALL_SITES flag. Following methods previously used for identifying repetitive sequence unsuitable for read mapping in primate Y-chromosomes, we used mapping and depth statistics from the VCF info field to identify callable regions on the canine Y-chromosome (Fig. 1) [26, 27]. Once a callable region was identified, we applied site level filtering to the remaining sequence: dropping maximum likelihood heterozygotes, missing sites, and positions with an MQ0/raw depth ratio > 0.10. A second depth filter was then applied to remove positions with extreme sequencing depths (median depth ± 3 M.A.D.). We also removed positions that were within 5 bp of GATK called indels."

From a Biostars post on transforming BAM to FASTA:
"# Get consensus fastq file
samtools mpileup -uf REFERENCE.fasta SAMPLE.bam | bcftools call -c | vcfutils.pl vcf2fq > SAMPLE_cns.fastq
# Convert .fastq to .fasta and set bases of quality lower than 20 to N
seqtk seq -aQ64 -q20 -n N SAMPLE_cns.fastq > SAMPLE_cns.fasta"

From a Biostars post on transforming a VCF to a FASTA file:
"You could call variants (using whatever variant calling software you like, GATK, freebayes etc.) from your .bam file and then use vcf-consensus (http://vcftools.sourceforge.net/perl_module.html#vcf-consensus) to build your consensus sequence. The code below should work:
cat ref.fa | vcf-consensus file.vcf.gz > out.fa"
