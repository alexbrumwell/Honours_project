# Code for genetic diversity and inbreeding analysis of Exmoor ponies

---
* Title: Code for genetic diversity and inbreeding analysis of Exmoor ponies
* Author: Alexandra Brumwell
---

Code notebook for honours project on Exmoor pony populations genetics

--------------------------------------------------------------------------------------------
## Introduction
This aims to look at genetic diversity and inbreeding of the Exmoor pony sequences.

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
***Note on code formatting:*** *\"\\" at the end of a line means that the command continues in the following line of code. In other words, the five lines of code in the `bwa mem` command are usually written in one line but, in this case, it is more visual to see each argument on a different line*


### 2. Converting to BAM and sorting
Now that the SAM files have been created, I convert them to BAM format so they take less space and sort them by reading name using `samtools view` and then `samtools sort`:

```bash
cat /shared5/Alex/Inbreeding_project/List_allexmoors_FASTQ.txt | awk '{print $1}' | while read name ; do
  samtools view -h -u -q 30 /shared5/Alex/Inbreeding_project/BAMs/${name}_horseref.sam |  samtools sort -n -o /shared5/Alex/Inbreeding_project/BAMs/${name}_sorted_horseref.bam &
done
```
> `samtools view` options:
* `-q`: *Minimal mapping quality*. Basic quality filter that skips alignments with a mapping quality smaller than the set value (in this case 30)
* `-u`: *Output uncompressed data*. Preferred when piping to another samtools command.

> `samtools sort`options:
* `-n`: *Sort by read names rather than by chromosomal coordinates.*. We are sorting by read name for later steps
* `-o`: *Write the final sorted output to specified file name*



### 3. Marking duplicates
I use `samtools markdup`, but prior to this, the files need to be processed with `samtools fixmate` and then sorted by coordinate using `samtools sort`:

```bash
cat /shared5/Alex/Inbreeding_project/List_allexmoors_FASTQ.txt | awk '{print $1}' | while read name ; do
  samtools fixmate -m /shared5/Alex/Inbreeding_project/BAMs/${name}_sorted_horseref.bam - | samtools sort - | samtools markdup - /shared5/Alex/Inbreeding_project/BAMs/${name}_markdup_sorted_horseref.bam &
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


### 4. Merging BAM files
Since some of the BAM files I have come from the same indivdual, I merge them together. belong to the same individual, I have to merge the corresponding BAM files together.

I compile the IDs of the indivudals that are spearted into multiple files and use `samtools merge`:
```bash
samtools merge <output file> <input file1> input file2>
```


### 5. Indexing
Now that the BAMs has been created, sorted and duplicates have been marked, they need to be indexed for future analyses.
```bash
for file in /shared5/Alex/Inbreeding_project/BAMs*markdup*bam*; do
  samtools index $file
done
```

### 6. Calculating BAM coverage
I use `qualimap bamqc` to calculate the coverage of the BAM files [(qualimap bamqc manual)](http://qualimap.conesalab.org/doc_html/analysis.html#bamqc).

*Important to note:* There is no need to specify an output as qualimap creates an output folder containing the stats results.

```bash
cat /shared5/Alex/Inbreeding_project/List_allexmoors_merged.txt | awk '{print $1}' | while read name ; do
  qualimap bamqc -bam /shared5/Alex/Inbreeding_project/BAMs/${name}_*markdup*bam \
  --java-mem-size=200G \ #Sets desired memory size
done
```

To check the coverage of each file:

```bash
cat /shared5/Alex/Inbreeding_project/List_allexmoors_merged.txt | awk '{print $1}' | while read name; do
  paste <(echo $name) <(grep  "mean cov" /shared5/Alex/Inbreeding_project/BAMs/${name}*stats/*txt)
done
```

### 7. Validating BAM Files
It is important to check the files:
```bash
cat /shared5/Alex/Inbreeding_project/List_allexmoors_merged.txt | while read file; do
  mkdir -p /shared5/Alex/Inbreeding_project/BAMs/ValidateSamFile/
  java -jar /shared5/Alex/picard-2.27.5/picard.jar ValidateSamFile \
  -I /shared5/Alex/Inbreeding_project/BAMs/${file}_markdup_sorted_horseref.bam \
  -O /shared5/Alex/Inbreeding_project/BAMs/ValidateSamFile/${file}_markdup_sorted_horseref_errors.txt \
  -R /shared5/Alex/Inbreeding_project/horse_WG_reference/GCF_002863925.1_EquCab3.0_genomic.fna \
  --MODE SUMMARY &
done
```


### 8. Variant calling
Since a VCF had already been created using the 2020 and 2021 sequences, I only call variants from the 2022 sequences.

#### 8.1 Strelka
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

#I create a configure_strelka file with all the pathfiles to the bams and then run it :
bash /shared5/Alex/Inbreeding_project/variants_allexmoors/configure_strelka.sh
```

The configure strelka script has created a pyhton file in each variant directory so now I run strelka:
```bash
for dir in /shared5/Alex/Inbreeding_project/variants/* ; do
  ${dir}/Out/runWorkflow.py -m local -j 2 &
done
```
This creates a VCF file with the variants from the 2022 Exmoor pony sequences


#### 8.2 Gather VCFs
Then I want to compile the VCFs of each chromosome into one VCF file. I will merge them using the gathervcf script. Structure of the script is:
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
Then, I run the script:
```bash
bash gathervcf.sh
```

#### 8.3 Merge VCFs
Now I want to merge the VCF that contains the 2022 samples with the previously created VCF that has the 2020 and 2021 sequences:
```bash
./vcf-merge /shared5/studentprojects/Qikai/Horse_database/fixmate_File/Exmoor_ponies/variants_rehead.vcf.gz /shared5/Alex/Inbreeding_project/variants/merged_Exmoor2022horseref_rename.vcf.gz > /shared5/Alex/Inbreeding_project/variants/merged_exmoor_8-2-2023.vcf
```

### 9. VCF filtering
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

Then I use R to calculate the top and bottom 2.5% quartile:
```R
R # This activates R in the command line

data = read.table("variants_rmvIndels_minQ30_minGQ30_mac3_hwe0.05_PASSonly_maxmissing0.8.ldepth.mean", header=T)
quantile(data$MEAN_DEPTH, probs=c(0.0275, 0.975)

# Output:
2.75%   97.5%
15.9286 24.1071

q() # quits R
```

Finally, to keep only the 95% quartile of mean loci depth, I use VCFtools again:
```bash
vcftools --vcf variants_rmvIndels_minQ30_minGQ30_mac3_hwe0.05_PASSonly_maxmissing0.8.recode.vcf --min-meanDP 9.7 --max-meanDP 14.4 --out variants_rmvIndels_minQ30_minGQ30_mac3_hwe0.05_PASSonly_maxmissing0.8_meanDPmid95percentile --recode
```
 We have created our final VCF and we can move onto calculating ROH.


Unfortunately, this final, filtered file has very few sites and contains only 15 of the inital 55 individuals.

### 9. Homozygosity and inbreeding using VCFtools
Since the filtered merged file has veyr low quality data, I decide to use the inital VCF file I created that had only had sequences from 2022.

```bash
vcftools --gzvcf /shared5/Alex/Inbreeding_project/variants/merged_Exmoor2022horseref.vcf.gz --out /shared5/Alex/Inbreeding_project/variants/vcftools_summarystats/merged_Exmoor2022horseref --het
```
