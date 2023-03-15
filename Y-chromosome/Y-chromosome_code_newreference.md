# Code for second Y-chromosomal analysis of Exmoor ponies

---
* Title: Code for updated Y-chromosome pipeline
* Author: Alexandra Brumwell
---

Code notebook for honours project on Exmoor pony populations genetics

--------------------------------------------------------------------------------------------
## Introduction
This is the updated pipeline for the studying the Y-chromosome in 15 male Exmoor ponies sequences.


## Downloading/obtaining a reference genome
The new reference was I used a Y-chromosome assembly (Genbank ID: GCA_002166905.2) from Felkel et al. (2019) that constituted 1.6 Mb of non-repetitive male specific regions. I downladed it into my directory using `wget -P /chosen_directory/ <copied_link>`.


## FASTQtoBAM pipeline
### 1. Mapping FASTQs
**Important to note**: The FASTQs I am working with are already trimmed.

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
  samtools view -h -b -F 4 -u /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764.sam | \
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
I use `samtools markdup`, but prior to this, the files need to be processed with `samtools fixmate` and then sorted by coordinate using `samtools sort`:

```bash
cat /shared5/Alex/Y-chromosome_project/List_ychromosome.txt | awk '{print $1}' | while read name ; do
  samtools fixmate -m /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted.bam - | \
  samtools sort - | \
  samtools markdup - /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted_markdup.bam &
done
```
> `samtools fixmate` options:
* `-m`: *Adds a mate score tags* These are used by markdup to select the best reads to keep.
* `-`: pipe output

> `samtools sort` options:
* *Note: you don't need to specify sorting by coordinate because that is the default option*
*`-`: pipe output

> `samtools markdup` options:
* `-`: pipe input


## 4. Merging BAM files
Since some of the sequences I have belong to the same individual, I have to merge the corresponding BAM files together.

`samtools merge` can be used for this:
```bash
samtools merge <output file> <input file1> <input file2>
```
Then, I create a new list (List_Y-chromosome_merged.txt) with the updated IDs of the sequences. E.g. The merged file IDs have now become "s479021_merged" instead of "s2479021_EDSW210003783-1a_H3WNKDSX2_L4".


## 5. Quality filtering
I filter reads based on options from Felkhel et al. (2019):

```bash
cat /shared5/Alex/Y-chromosome_project/List_Y-chromosome_merged.txt | awk '{print $1}' | while read name ; do
  samtools view -h -b -F 4 -q 20 /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted_markdup.bam -o /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted_markdup_MQ20.bam &
done
```
> `samtools view` options:
* `-h`: *Include the header in the output*
* `-b`: *Output in BAM format*
* `-F`: *Exclusion flag*. In this case the flag is set at "4" so it will exclude unmapped reads.
* `-q`: *Minimal mapping quality*. Basic quality filter that skips alignments with a mapping quality smaller than the set value (in this case 20)


## 6. Indexing
Now that the BAMs has been created, sorted and duplicates have been marked, they need to be indexed for future analyses.

```bash
cat /shared5/Alex/Y-chromosome_project/List_Y-chromosome_merged.txt | awk '{print $1}' | while read name ; do
  samtools index /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted_markdup_MQ20.bam &
done
```

## 7. Validating BAM files
It is good practice to make sure that the BAMs have been created correctly. I did this after step 1, 2, 3 and 5.

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
  qualimap bamqc -bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted_markdup_MQ20.bam --java-mem-size=200G 2>/shared5/Alex/Y-chromosome_project/BAMs_LipY764/log_files/${name}_qualimap.log & # Sets desired memory size
done
```

Then to see how the mean coverage and number of mapped reads from each sequence:
```bash
#new ref
cat /shared5/Alex/Y-chromosome_project/List_Y-chromosome_merged.txt | awk '{print $1}' | while read name; do
  paste <(echo $name) <(grep  "mean coverageData" /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}*stats/*txt) <(grep  "number of reads" /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}*stats/*txt)
done
```



## Prelimnary consensus FASTA creation
I initially tried to create consensus FASTAs for ech sequence, but these ended up being too large for analyses so they were not used in the end.

```bash
#creating depth list
cat /shared5/Alex/Y-chromosome_project/List_Y-chromosome_merged.txt | awk '{print $1}' | while read name; do
  coverage=$(cat /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}*stats/genome_results.txt | grep "mean coverageData"  | awk '{print $4}' | sed 's/X//' | sed 's/,//' | cut -d "." -f 1);
  paste <(echo ${name}) <(echo ${coverage}) <(echo "${coverage} * 2" | bc);
done > /shared5/Alex/Y-chromosome_project/Y-chromosome_ID_depth_list_new.txt

#creating consensus FASTAs
cat /shared5/Alex/Y-chromosome_project/Y-chromosome_ID_depth_list_new.txt | while read line; do

  #Defining variables:
  name=$(echo ${line} | awk '{print $1}');
  number=$(echo ${line} | awk '{print $(NF-1)}');
  maxnumber=$(echo ${line} | awk '{print $(NF)}');

  angsd -i /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted_markdup_MQ20.bam \
  -minMapQ 30 \ #minimum mapping quality
  -minQ 20 \ #minimum read quality
  -setMinDepth ${number} \ #minimum depth
  -setMaxDepth ${maxnumber} \ #maximum depth
  -remove_bads 1 \ #Discards 'bad' reads
  -doFasta 2 \ #creates FASTA using the most common (non N) base
  -doCounts 1 \ #Counts the number of A,C,G,T.
  -ref /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/GCA_002166905.2_LipY764_genomic.fna \ #reference
  -out /shared5/Alex/Y-chromosome_project/FASTAs_LipY764/${name}_LipY764_Y & #output file
done

#unzip the compressed FASTAs
gunzip *fa.gz

#change header IDs
cat /shared5/Alex/Y-chromosome_project/Y-chromosome_ID_depth_list.txt | while read line; do
  name=$(echo $line | awk '{print $1}');
  breed=$(echo $line | awk '{print $NF}');
  echo ">${name}_${breed}_LipY764" | cat - /shared5/Alex/Y-chromosome_project/FASTAs_LipY764/${name}_5x_LipY764_Y.fa > tmp.fa && mv tmp.fa /shared5/Alex/Y-chromosome_project/FASTAs_LipY764/${name}_LipY764_Y.fa
done

#concatenating FASTAs into one file
cat /shared5/Alex/Y-chromosome_project/Y-chromosome_ID_depth_list.txt | awk '{print $1}' | while read name; do
  cat /shared5/Alex/Y-chromosome_project/FASTAs_LipY764/${name}_LipY764_Y.fa
done >> Exmoor_horse_horseref_NEW_Y.fa


## adding reference sequence
cat /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/GCA_002166905.2_LipY764_genomic.fna Exmoor_horse_horseref_NEW_Y.fa > tmp.fa && mv tmp.fa Exmoor_horse_horseref_NEW_Y.fa
```



## 9. Variant calling with Bcftools
BCFtools mpileup and call (Li, 2011) was used to call variants from the created Y-chromosome BAM files.

#### 9.1 BCFtools mpileup
```bash
bcftools mpileup -f /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/GCA_002166905.2_LipY764_genomic.fna /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_102004_EKDN220034353-1A_H5GL5DSX5_L2_LipY764_sorted_markdup_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_107013_EKDN220034352-1A_H5GL5DSX5_L2_LipY764_sorted_markdup_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_21084_EKDN220034367-1A_H5J5MDSX5_L2_LipY764_sorted_markdup_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_23279_merged_LipY764_sorted_markdup_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_23416_merged_LipY764_sorted_markdup_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_320005_merged_LipY764_sorted_markdup_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_32023_merged_LipY764_sorted_markdup_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_49031_EKDN220034370-1A_H5J5MDSX5_L2_LipY764_sorted_markdup_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_519005_EKDN220034364-1A_H5J5MDSX5_L2_LipY764_sorted_markdup_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/s2010_EDSW210003772-1a_H3WNKDSX2_L3_LipY764_sorted_markdup_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/s2012_EDSW210003766-1a_H3WHWDSX2_L4_LipY764_sorted_markdup_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/s479013_EDSW210003775-1a_H3WNKDSX2_L3_LipY764_sorted_markdup_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/s49097_EDSW210003769-1a_H3FTWDSX2_L2_LipY764_sorted_markdup_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/s49124_EDSW210003767-1a_H3WHWDSX2_L2_LipY764_sorted_markdup_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/S276023_EDSW200015547-1a2a_HJFWVDSXY_L3L4_LipY764_sorted_markdup_MQ20.bam  > /shared5/Alex/Y-chromosome_project/variants_bcftools_onlyexmoors/Exmoors_2022_Y_newref.pileup
```

#### 9.2 BCFtools call
BCFtools call requieres a ploidy file that basically contains the number of bases in the sequences.
```bash
bcftools call -c -O v --ploidy-file /shared5/Alex/Y-chromosome_project/variants_bcftools_onlyexmoors/ploidy_newref.txt -o /shared5/Alex/Y-chromosome_project/variants_bcftools_onlyexmoors/Exmoors_2022_Y_newref.vcf /shared5/Alex/Y-chromosome_project/variants_bcftools_onlyexmoors/Exmoors_2022_Y_newref.pileup
```


# 10. VCF filtering
1.1 First step, quality filter and removing INDELS:
```bash
vcftools --vcf /shared5/Alex/Y-chromosome_project/variants_bcftools_onlyexmoors/Exmoors_2022_Y_newref.vcf --remove-indels --minQ 30 --minGQ 30 --out /shared5/Alex/Y-chromosome_project/variants_bcftools_onlyexmoors/Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30 --recode
```

1.2 To see how much missing information there is per individual:
`vcftools –gzvcf <data.vcf.gz> –missing-indv –out <data>`

2. Minimum allele count filter
vcftools --vcf Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30.recode.vcf --min-alleles 2 --mac 3 --out Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30_biallelic_mac3 --recode

3. Exclude sites where >10% of data is missing
vcftools --vcf Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30_biallelic_mac3.recode.vcf --max-missing 0.9

#To see missingness of individuals:
```bash
vcftools --vcf Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30_biallelic_mac3.recode.vcf --missing-indv --out Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30_biallelic_mac3.recode.vcf
```

OTHER VCFTOOLS CODE:
vcftools --gzvcf $SUBSET_VCF --depth --out $OUT #calculates depth per indivudal
vcftools --gzvcf $SUBSET_VCF --site-mean-depth --out $OUT #calculates depth per site
vcftools --gzvcf $SUBSET_VCF --missing-site --out $OUT #missing data per site
vcftools --gzvcf $SUBSET_VCF --missing-ind --out $OUT #missing data per individual


# 11. PLINK
First, create a PLINK file and then calculate eigen values and vectors:
```bash
#convert to plink format
vcftools --vcf Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30_biallelic_mac3.recode.vcf --plink --out Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30_biallelic_mac3

#PCA
plink --file Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30_biallelic_mac3 --pca --out Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30_biallelic_mac3
```
