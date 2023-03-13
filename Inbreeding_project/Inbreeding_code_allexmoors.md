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
cat /shared5/Alex/Inbreeding_project/List_allexmoors_FASTQ.txt | awk '{print $2}' | while read file ; do

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
cat /shared5/Alex/Inbreeding_project/List_allexmoors_FASTQ.txt | awk '{print $1}' | while read name ; do
  java -jar /shared5/Alex/picard-2.27.5/picard.jar ValidateSamFile -I /shared5/Alex/Inbreeding_project/BAMs/${name}*sam -O /shared5/Alex/Inbreeding_project/BAMs/${name}_sam_horseref_validate.txt --MODE SUMMARY -R /shared5/Alex/Inbreeding_project/horse_WG_reference/GCF_002863925.1_EquCab3.0_genomic.fna &
done
```



### 2. Converting to BAM and sorting
Now that the SAM files have been created, I convert them to BAM format so they take less space and sort them by reading name using `samtools view` and then `samtools sort`:

```bash
cat /shared5/Alex/Inbreeding_project/List_allexmoors_FASTQ.txt | awk '{print $1}' | while read name ; do
  samtools view -h -u -q 30 /shared5/Alex/Inbreeding_project/BAMs/${name}_horseref.sam |  samtools sort -n -o /shared5/Alex/Inbreeding_project/BAMs/${name}_sorted_horseref.bam &
done
```


### 3. Marking duplicates
```bash
# Using Samtools markduplicates:
cat /shared5/Alex/Inbreeding_project/List_allexmoors_FASTQ.txt | awk '{print $1}' | while read name ; do
  samtools fixmate -m /shared5/Alex/Inbreeding_project/BAMs/${name}_sorted_horseref.bam - | samtools sort - | samtools markdup - /shared5/Alex/Inbreeding_project/BAMs/${name}_markdup_sorted_horseref.bam &
done

# Originally tried using Picardtools Markduplicates but it didn't work; it didn't create any output file beacuse there were repeated reads or something like that
cat /shared5/Alex/Inbreeding_project/List_Sep_2022_FASTQ.txt | awk '{print $1}' | while read name ; do
  java -Xmx4g -jar /shared5/Alex/picard-2.27.5/picard.jar MarkDuplicates  \
    -I /shared5/Alex/Inbreeding_project/BAMs/${name}_sorted_horseref.bam \ # input file
    -O /shared5/Alex/Inbreeding_project/BAMs/${name}_rmdp_sorted_horseref.bam \ #output file
    -M /shared5/Alex/Inbreeding_project/BAMs/log_files/${name}_markdup_metrics.txt \ # text file with marked_dup_metrics
    --REMOVE_DUPLICATES true --ASSUME_SORT_ORDER coordinate --VALIDATION_STRINGENCY SILENT & \
done
```


### 4. Merging BAM files
Merge whichever individuals have several files:
```bash
#loop to print out code of files that need to be merged:
ll *markdup*bam | awk '{print $9}' | awk '{print substr($0,1,10)}' | uniq -c  | grep -v " 1 " | awk '{print $2}' | cut -d "_" -f 1,2 | while read file; do
  paste <(echo "samtools merge ${file}_merged_markdup_sorted_horseref.bam") <(ls ${file}*markdup*bam) ;
done
 #output:
samtools merge E_23279_merged_sorted_horseref.bam E_23279_EKDN220034363-1A_H5J5MDSX5_L2_markdup_sorted_horseref.bam E_23279_EKDN220034363-1A_H7VKMDSX5_L3_markdup_sorted_horseref.bam

#run code to merge files:
samtools merge -f E_23279_merged_markdup_sorted_horseref.bam E_23279_EKDN220034363-1A_H5J5MDSX5_L2_markdup_sorted_horseref.bam E_23279_EKDN220034363-1A_H7VKMDSX5_L3_markdup_sorted_horseref.bam & samtools merge -f E_23416_merged_markdup_sorted_horseref.bam E_23416_EKDN220034368-1A_H5J5MDSX5_L2_markdup_sorted_horseref.bam  E_23416_EKDN220034368-1A_H7VKMDSX5_L1_markdup_sorted_horseref.bam E_23416_EKDN220034368-1A_H7VL2DSX5_L2_markdup_sorted_horseref.bam & samtools merge -f E_23434_merged_markdup_sorted_horseref.bam E_23434_EKDN220034361-1A_H5J5MDSX5_L2_markdup_sorted_horseref.bam E_23434_EKDN220034361-1A_H7VKMDSX5_L3_markdup_sorted_horseref.bam & samtools merge -f E_320005_merged_markdup_sorted_horseref.bam E_320005_EKDN220034373-1A_H5J5MDSX5_L2_markdup_sorted_horseref.bam E_320005_EKDN220034373-1A_H7VKMDSX5_L1_markdup_sorted_horseref.bam E_320005_EKDN220034373-1A_H7VL2DSX5_L4_markdup_sorted_horseref.bam & samtools merge -f E_32023_merged_markdup_sorted_horseref.bam E_32023_EKDN220034372-1A_H5J5MDSX5_L2_markdup_sorted_horseref.bam E_32023_EKDN220034372-1A_H7VKMDSX5_L1_markdup_sorted_horseref.bam & samtools merge -f E_479023_merged_markdup_sorted_horseref.bam E_479023_EKDN220034374-1A_H5J5MDSX5_L2_markdup_sorted_horseref.bam  E_479023_EKDN220034374-1A_H7VKMDSX5_L1_markdup_sorted_horseref.bam  E_479023_EKDN220034374-1A_H7VL2DSX5_L2_markdup_sorted_horseref.bam & samtools merge -f E_512001_merged_markdup_sorted_horseref.bam E_512001_EKDN220034365-1A_H5J5MDSX5_L2_markdup_sorted_horseref.bam  E_512001_EKDN220034365-1A_H7VKMDSX5_L3_markdup_sorted_horseref.bam & samtools merge -f E_78170_merged_markdup_sorted_horseref.bam E_78170_EKDN220034369-1A_H5J5MDSX5_L2_markdup_sorted_horseref.bam E_78170_EKDN220034369-1A_H7VKMDSX5_L1_markdup_sorted_horseref.bam  E_78170_EKDN220034369-1A_H7VL2DSX5_L4_markdup_sorted_horseref.bam & samtools merge -f E_900588_merged_markdup_sorted_horseref.bam E_900588_EKDN220034362-1A_H5J5MDSX5_L2_markdup_sorted_horseref.bam E_900588_EKDN220034362-1A_H7VKVDSX5_L2_markdup_sorted_horseref.bam & samtools merge -f E_900717_merged_markdup_sorted_horseref.bam E_900717_EKDN220034357-1A_H5CYKDSX5_L4_markdup_sorted_horseref.bam  E_900717_EKDN220034357-1A_H5HJMDSX5_L1_markdup_sorted_horseref.bam  & samtools merge s479021_merged_markdup_sorted_horseref.bam s479021_EDSW210003765-1a_H3WHWDSX2_L4_sorted_horseref.bam s479021_EDSW210003765-2a_H2Y75DSX2_L2_sorted_horseref.bam s479021_EDSW210003765-2a_H2Y7CDSX2_L3_sorted_horseref.bam &  samtools merge s49127_merged_markdup_sorted_horseref.bam s49127_EDSW210003774-1a_H3FTWDSX2_L1_sorted_horseref.bam s49127_EDSW210003774-2a_H2Y75DSX2_L4_sorted_horseref.bam s49127_EDSW210003774-2a_H2Y7CDSX2_L3_sorted_horseref.bam
```


### 5. Indexing
Now that the BAMs has been created, sorted and duplicates have been marked, they need to be indexed for future analyses.
```bash
for file in /shared5/Alex/Inbreeding_project/BAMs*markdup*bam*; do
  samtools index $file
done
```

### 6. Qualimap
I don't run it in parallel because that was giving me errors.
```bash
cat /shared5/Alex/Inbreeding_project/List_Sep_2022_merged.txt | while read file; do
  qualimap bamqc --java-mem-size=200G -bam ${file}*markdup*bam
done
```

Checking output stats:
```bash

cat /shared5/Alex/Inbreeding_project/Allexmoor_merged.txt | awk '{print $1}' | while read name; do
  paste <(echo $name) <(grep  "mean coverageData" /shared5/Alex/Inbreeding_project/BAMs/${name}*stats/*txt) <(grep "number of reads" /shared5/Alex/Inbreeding_project/BAMs/${name}*stats/*txt)
done
```


Still missing:

grep: /shared5/Alex/Inbreeding_project/BAMs/s2012_EDSW210003766-1a_H3WHWDSX2_L4*stats/*txt: No such file or directory
grep: /shared5/Alex/Inbreeding_project/BAMs/s237009_EDSW210003788-1a_H3WNKDSX2_L2*stats/*txt: No such file or directory
grep: /shared5/Alex/Inbreeding_project/BAMs/s479021_merged*stats/*txt: No such file or directory
grep: /shared5/Alex/Inbreeding_project/BAMs/s49127_merged*stats/*txt: No such file or directory
grep: /shared5/Alex/Inbreeding_project/BAMs/S276023_EDSW200015547-1a2a_HJFWVDSXY_L3L4*stats/*txt: No such file or directory

### 7. Validating BAM Files
It is important to check the files:
```bash
cat /shared5/Alex/Inbreeding_project/List_Sep_2022_merged.txt | while read file; do
  mkdir -p /shared5/Alex/Inbreeding_project/BAMs/ValidateSamFile/
  java -jar /shared5/Alex/picard-2.27.5/picard.jar ValidateSamFile \
  -I /shared5/Alex/Inbreeding_project/BAMs/${file}_markdup_sorted_horseref.bam \
  -O /shared5/Alex/Inbreeding_project/BAMs/ValidateSamFile/${file}_markdup_sorted_horseref_errors.txt \
  -R /shared5/Alex/Inbreeding_project/horse_WG_reference/GCF_002863925.1_EquCab3.0_genomic.fna \
  --MODE SUMMARY &
done
```


## VARIANT CALLING using Strelka

### 1. Strelka
Preparing for strelka:
 ```bash
#Create a directory for each chromosome.
for num in {144..175}; do
  mkdir "/shared5/Alex/Inbreeding_project/variants_allexmoors/NC_009${num}.3"
done

#Make text files wich each line being a chromosome from the BED reference file:
cat /shared5/Alex/Inbreeding_project/horse_WG_reference/GCF_002863925.1_EquCab3.0_genomicScaff.bed | grep "NC_009"  | while read line; do
  scaffold=$(echo $line | awk '{print $1}');
  echo $line > /shared5/Alex/Inbreeding_project/variants_allexmoors/${scaffold}/${scaffold}.bed
done

#make sure bed files are tab-delimited:
for file in {144..175}; do
  cat NC_009${file}.3/NC_009${file}.3.bed | awk 'OFS="\t" {$1=$1; print}' > tmp.bed && mv tmp.bed NC_009${file}.3/NC_009${file}.3.bed;
done

#zip bed files for strelka
for file in {144..175}; do
  bgzip -c NC_009${file}.3/NC_009${file}.3.bed > NC_009${file}.3/NC_009${file}.3.bed.gz
  tabix -f -p bed NC_009${file}.3/NC_009${file}.3.bed.gz
done

#I copy the configure_strelka file form Run2 and edit it:
cat List_Sep_2022_merged.txt | while read name ; do
  echo "--bam /shared5/Alex/Inbreeding_project/BAMs/${name}_markdup_sorted_horseref.bam \\";
done
```

Running strelka:
```bash
for dir in /shared5/Alex/Inbreeding_project/variants/* ; do
  ${dir}/Out/runWorkflow.py -m local -j 2 &
done
```

### 2. Gather VCFs
Now the VCFs have been created for our files. We will merge them using the gathervcf script. Structure of the script is:
"picard GatherVcfs \
I=NC_052177/Out/results/variants/variants.vcf.gz \
I=NC_052178/Out/results/variants/variants.vcf.gz
O=merged_Exmoor2022horseref.vcf.gz"

In order to make this script:
```bash
for num in {144..175}; do
  echo "I=NC_009${num}.3/Out/results/variants/variants.vcf.gz \\"
done
```

Then, want to change the header of this new vcf file:
```bash
#select header and save in new file
less results/variants/variants.vcf.gz | head -n 50000 | grep -e "#" > head.txt

#Then I replace the line of "#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE1 SAMPLE2 SAMPLE3 SAMPLE4 SAMPLE5 SAMPLE6 SAMPLE7 SAMPLE8 SAMPLE9 SAMPLE10  SAMPLE11 SAMPLE12        SAMPLE13        SAMPLE14        SAMPLE15        SAMPLE16        SAMPLE17        SAMPLE18        SAMPLE19        SAMPLE20        SAMPLE21   SAMPLE22        SAMPLE23        SAMPLE24        SAMPLE25" with the sample names from the configure strelka file. So the final line is: "
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  E_102004        E_107013        E_21084 E_21149 E_23279 E_23416 E_23434 E_320005  E_32023  E_44009 E_458028        E_479023        E_49031 E_49057 E_49121 E_512001        E_519005        E_78170 E_900234        E_900585        E_900588  E_900600 E_900694        E_900717        E_900741"

#I make sure that the strings are separated by tab by using
tail head.txt | awk -v OFS="\t" '{$1=$1; print}'

#then need to merge the header with the VCF:
tabix -r head.txt merged_Exmoor2022horseref.vcf.gz > merged_Exmoor2022horseref_rename.vcf.gz

#Then we index this new vcf file:
tabix -p vcf merged_exmoor_8-2-2023.vcf.gz
```



Then we merge the 2022 and 2021 Exmoor VCFs:
```bash
export PATH=/shared3/Anubhab/miniconda3/bin:$PATH #path to activate vcftools

cd shared5/Alex/software/vcftools/src/perl
./vcf-merge /shared5/studentprojects/Qikai/Horse_database/fixmate_File/Exmoor_ponies/variants_rehead.vcf.gz /shared5/Alex/Inbreeding_project/variants/merged_Exmoor2022horseref_rename.vcf.gz > /shared5/Alex/Inbreeding_project/variants/merged_exmoor_8-2-2023.vcf

```

VCFTOOLS --> export PATH=/shared3/Anubhab/miniconda3/bin:$PATH



## VCF FILTERING
Now that the VCF file has been created it is important to filter the variants. I use `VCFtools` for this.

1.	The first filter removes insertions and deletions and low quality positions and variants:
```bash
vcftools --vcf /shared5/Alex/Inbreeding_project/variants/merged_exmoor_8-2-2023.vcf --remove-indels --minQ 30 --minGQ 30 --out /shared5/Alex/Inbreeding_project/variants/merged_exmoor_8-2-2023_rmvIndels_minQ30_minGQ30 --recode
```
* `--gzvcf`: *Input file name of compressed vcf*
* `--remove-indels`: *Include or exclude sites that contain an indel.* For these options "indel" means any variant that alters the length of the REF allele.
* `--minQ`: *Includes only sites with Quality value above this threshold*
* `-minGQ`: *Exclude all genotypes with a quality below the threshold specified*
* `--out`: *Specify output file name and location*
* `-recode`: *Generates a new file with applied filters*


2.	The second filter is a minimum allele count and also removes non-neutral alleles:
```bash
vcftools --vcf /shared5/Alex/Inbreeding_project/variants/merged_exmoor_8-2-2023_rmvIndels_minQ30_minGQ30.recode.vcf --mac 3 --hwe 0.05 --remove-filtered-all --out /shared5/Alex/Inbreeding_project/variants/merged_exmoor_8-2-2023_rmvIndels_minQ30_minGQ30_mac3_hwe0.05_PASSonly --recode
```
* `--vcf`: *Input file name*
* `--mac`: *Minimum allele count*. Minimum allele count set at 3, meaning that each allele must appear at least 3 times over all individuals for that site
* `--hwe`: *Assesses alleles for Hardy-Weinberg equilibrium.* Parameter is set 0.05 (i.e. p-value set at 0.05) so only retain alleles that are neutral.
* `--remove-filtered-all`: *Removes all sites with a FILTER flag other than PASS.*

3. Remove indiviudals with low depth
#To see Individual Missingness:
vcftools --vcf merged_exmoor_8-2-2023_rmvIndels_minQ30_minGQ30_mac3_hwe0.05_PASSonly.recode.vcf --missing-indv --out  merged_exmoor_8-2-2023_rmvIndels_minQ30_minGQ30_mac3_hwe0.05_PASSonly

Anubhab removed several of the indiviudals that had low coverage and the nwe aaply site-mean-depth

4. Remove the top and bottom 2.5% mean site depth
```bash
vcftools --vcf variants_rmvIndels_minQ30_minGQ30_mac3_hwe0.05_PASSonly_maxmissing0.8.recode.vcf --site-mean-depth --out variants_rmvIndels_minQ30_minGQ30_mac3_hwe0.05_PASSonly_maxmissing0.8
```
* `--site-mean-depth`: *Generates a file containing the mean depth per site averaged across all individuals*

Then we use R to calculate the top and bottom 2.5% quartile:
```R
R # This activates R in the command line

data = read.table("variants_rmvIndels_minQ30_minGQ30_mac3_hwe0.05_PASSonly_maxmissing0.8.ldepth.mean", header=T)
quantile(data$MEAN_DEPTH, probs=c(0.0275, 0.975)

# Output:
2.75%   97.5%
15.9286 24.1071

q() # quits R
```


Finally, to keep only the 95% quartile of mean loci depth, we use VCFtools again:
```bash
vcftools --vcf variants_rmvIndels_minQ30_minGQ30_mac3_hwe0.05_PASSonly_maxmissing0.8.recode.vcf --min-meanDP 9.7 --max-meanDP 14.4 --out variants_rmvIndels_minQ30_minGQ30_mac3_hwe0.05_PASSonly_maxmissing0.8_meanDPmid95percentile --recode
```
 We have created our final VCF and we can move onto calculating ROH.

 ## CALCULATING ROH
Then BCFtools ROH to calculate ROH in the indivudals
