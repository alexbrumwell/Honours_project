# Code for Run1 of phylogeny of Exmoor ponies and horses

---
* Title: Code for Run1 of Exmoor and horse phylogeny
* Author: Alexandra Brumwell
---

Code notebook for honours project on Exmoor pony populations genetics

--------------------------------------------------------------------------------------------
## Introduction
This sideproject aimed to look at the phylogeny of the Exmoor pony by comparing an Exmoor sequence with other horse breeds, Przewalskis and a donkey sequence (an outgroup).

### Downloading a reference genome (27/09/22)
I am mapping all the sequences to an outgroup (donkey), so I dow load the most complete donkey assembly that I found in NCBI. There is only [one assembly](https://www.ncbi.nlm.nih.gov/assembly/GCF_016077325.2 ) that is completeed at a chromosomal level so I download it.
To download it, I go to the ["FTP directory for RefSeq assembly"](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/077/325/GCF_016077325.2_ASM1607732v2/) link on the right-hand side of the page. Then I copy the link of the compressed FASTA file [(GCF_016077325.2_ASM1607732v2_genomic.fna.gz)](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/077/325/GCF_016077325.2_ASM1607732v2/GCF_016077325.2_ASM1607732v2_genomic.fna.gz).

In the command line I created a folder for the reference genome and then used `wget` to download the FASTA:
```bash
mkdir Donkey_ref # creates directory
wget  -P /Donkey_ref/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/077/325/GCF_016077325.2_ASM1607732v2/GCF_016077325.2_ASM1607732v2_genomic.fna.gz #downloads FASTA file and saves to "Donkey_ref" directory
```


## FASTQ to BAM pipeline

#### Donkey sequence (30/09/22)
**0. Downloading FASTQs**
I went to the  SRA repository on NCBI and searched for donkey genomes. I found one [(SRA ID: SRR7031483)](https://www.ncbi.nlm.nih.gov/sra/?term=SRR7031483) which had paired reads and was a reasonable size of 20 Gb. Then, I went to the European nucleotide archive and used the Run ID (SRX3963508) to search for the FASTQ files of the sample. I copied the FTP link and downladed the FASTQs using `wget`.

```bash
wget -P /shared5/Alex/Donkey/FASTQs/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR703/003/SRR7031483/SRR7031483_1.fastq.gz
wget -P /shared5/Alex/Donkey/FASTQs/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR703/003/SRR7031483/SRR7031483_2.fastq.gz
```

>* **Note:** When choosing the sample we made sure that the FASTQs are not from the same individual or project used for the donkey reference as that would bias the results
>* **Note:** In ENA there are two different FASTQ files that we can download: the generated FASTQ files and the submitted FASTQ files. I downloaded the  generated FASTQ files.


**1. Trimming FASTQs**
I use `trimgalore`
```bash
trim_galore --paired SRR7031483_1.fastq SRR7031483_2.fastq
```

**2. Mapping**
```bash
bwa mem /shared5/Alex/Donkey_ref/GCF_016077325.2_ASM1607732v2_genomic.fna.gz /shared5/Alex/Donkey/FASTQs/SRR7031483_1_val_1.fq /shared5/Alex/Donkey/FASTQs/SRR7031483_2_val_2.fq > /shared5/Alex/Donkey/BAMs/SRR7031483.sam
```

**3. Coverting to BAM and sorting by read names**
```bash
samtools view -q 30 -u /shared5/Alex/Donkey/BAMs/SRR7031483.sam | samtools sort -n -o /shared5/Alex/Donkey/BAMs/SRR7031483_sorted.bam
```

**4. Mark duplicates**
```bash
samtools fixmate -m /shared5/Alex/Donkey/BAMs/SRR7031483_sorted.bam - | samtools sort - | samtools markdup - /shared5/Alex/Donkey/BAMs/SRR7031483_markdup.bam
```

**5. Indexing**
```bash
samtools index /shared5/Alex/Donkey/BAMs/SRR7031483_markdup.bam
```

**6. Calculating BAM depth**
I used `qualimap` and this creates a folder that includes a file with the coverage data.
```bash
qualimap bamqc -bam /shared5/Alex/Donkey/BAMs/SRR7031483_markdup.bam --java-mem-size=200G
```

**7. Consensus FASTA file**
```bash
angsd -i /shared5/Alex/Donkey/BAMs/SRR7031483_markdup.bam -minMapQ 30 -minQ 20 -setMinDepth 9 -setMaxDepth 18 -remove_bads 1 -doFasta 2 -doCounts 1 -ref /shared5/Alex/Donkey_ref/GCF_016077325.2_ASM1607732v2_genomic.fna -out /shared5/Alex/Donkey/FASTAs/SRR7031483_DonkeyRef
```



#### Exmoor FASTQs
The Exmoor sequence I used was copied into my direcotry by ANubhab. The sample ID is s2010.

**1. Trimming FASTQs**
The FASTQ has already been trimmed

**2. Mapping**
```bash
bwa mem /shared5/Alex/Donkey_ref/GCF_016077325.2_ASM1607732v2_genomic.fna.gz /shared5/Alex/Exmoor/FASTQs/s2010_EDSW210003772-1a_H3WNKDSX2_L3_1_val_1.fq.gz /shared5/Alex/Exmoor/FASTQs/s2010_EDSW210003772-1a_H3WNKDSX2_L3_2_val_2.fq.gz > /shared5/Alex/Exmoor/BAMs/s2010_EDSW210003772.sam
```

**3. Coverting to BAM and sorting by read names**
```bash
samtools view -q 30 -u /shared5/Alex/Exmoor/BAMs/s2010_EDSW210003772.sam | samtools sort -n -o /shared5/Alex/Exmoor/BAMs/s2010_EDSW210003772_sorted.bam
```

**4. Mark duplicates**
```bash
samtools fixmate -m /shared5/Alex/Exmoor/BAMs/s2010_EDSW210003772_sorted.bam - | samtools sort - | samtools markdup - /shared5/Alex/Exmoor/BAMs/s2010_EDSW210003772_markdup.bam
```

**5. Indexing**
```bash
samtools index /shared5/Alex/Exmoor/BAMs/s2010_EDSW210003772_markdup.bamm
```

**6. Calculating BAM depth**
I used `qualimap` and this creates a folder that includes a file with the coverage data.
```bash
qualimap bamqc -bam /shared5/Alex/Exmoor/BAMs/s2010_EDSW210003772_markdup.bam --java-mem-size=200G
```

**7. Consensus FASTA file**
```bash
angsd -i /shared5/Alex/Exmoor/BAMs/s2010_EDSW210003772_markdup.bam -minMapQ 30 -minQ 20 -setMinDepth 9 -setMaxDepth 18 -remove_bads 1 -doFasta 2 -doCounts 1 -ref /shared5/Alex/Donkey_ref/GCF_016077325.2_ASM1607732v2_genomic.fna -out /shared5/Alex/Exmoor/FASTAs/s2010_EDSW210003772_DonkeyRef
```



#### Horse sequences (03/10/22-11/10/22)
A previous masters student studied horse breed diversity so I will be using the sequnces that they downloaded. The trimmed FASTQ files of these sequences are in this students directory (/shared5/studentprojects/Qikai/Horse_database/). I have a metadata list from this previous student that shows a total of 56 sequences representing 25 horse breeds. However, for this first phylogeny analysis I only need one sequence from each horse breed.

**0. Obtaining FASTQs**
1.	I review the metadata spreadsheet from the student and select the sequence with the highest coverage from each horse breed and then compile these sequences into the spreadsheet for this project called Run1. I now have a total of 25 horse sequnces, each representing a different horse breed.

2.	I create a text file containing all the IDs of the sequences I will be using and I save this textfile (horse_ID_list.txt) in my directory.

3.	Since the FASTQ files of these sequences are located in different directories, I use a loop to locate where each FASTQ file is and save the file paths in a text file called "horse_fastq_paths" :
```bash
cat /shared5/Alex/Horses/horse_ID_list.txt | awk '{print $1}' | while read name; do
  ls /shared5/studentprojects/Qikai/Horse_database/*/${name}*val* ;
done >> horse_fastq_paths.txt
\
#If there was file that wasn't located then I find it by using and then add it to "horse_fastq_paths.txt":
find shared5/studentprojects/Qikai/Horse_database/ -name “ERR1527947*”
```

4. I created a directory that would contain the soft links to all the FASTQs, so that I could visualize all the horse sequences in one  place:
```bash
mkdir /shared5/Alex/Horses/FASTQs/ #creates directory

cat /shared5/Alex/Horses/horse_fastq_paths.txt | while read line; do
  ln -s ${line} /shared5/Alex/Horses/FASTQs/ ;
done
```
**1. Trimming FASTQs**
The horse FASTQs are already trimmed so I skip this step.

**2. Mapping**
I used a loop to go over the horse sequence IDs so I can map them in parallel
```bash
cat /shared5/Alex/Horses/horse_ID_list.txt | awk '{print $1}' | while read name; do
bwa mem /shared5/Alex/Donkey_ref/GCF_016077325.2_ASM1607732v2_genomic.fna.gz /shared5/Alex/Horses/FASTQs/${name}_1_val_1.fq.gz /shared5/Alex/Horses/FASTQs/${name}_2_val_2.fq.gz > /shared5/Alex/Horses/BAMs/${name}.sam & # the "&" allows the loops to be run in parallel
done
```

* **Sidenote to self:** To check how large the SAM files created by the previous student were, I used:
```bash
cat /shared5/Alex/Horses/horse_ID_list.txt | cut -f 3 | cut -f 1 -d "_") | while read name;
do
file=$(echo "/shared5/studentprojects/Qikai/Horse_database/$name.sam");
ll -h $file;
done
```


**3. Coverting to BAM and sorting by read names**
```bash
cat /shared5/Alex/Horses/horse_ID_list.txt | awk '{print $1}' | while read name; do
samtools view -q 30 -u /shared5/Alex/Horses/BAMs/${name}.sam | samtools sort -n -o /shared5/Alex/Horses/BAMs/${name}_sorted.bam &
done
```

**4. Mark duplicates**
```bash
cat /shared5/Alex/Horses/horse_ID_list.txt | awk '{print $1}' | while read name; do
samtools fixmate -m /shared5/Alex/Horses/BAMs/${name}_sorted.bam - | samtools sort - | samtools markdup - /shared5/Alex/Horses/BAMs/${name}_markdup.bam &
done
```

**5. Indexing**
```bash
cat /shared5/Alex/Horses/horse_ID_list.txt | awk '{print $1}' | while read name; do
  samtools index /shared5/Alex/Horses/BAMs/${name}_markdup.bam &
  done
```


**6. Calculating BAM depth**
I used `qualimap` and this creates a folder that generates a file called "genome_results.txt" for each sequence with the coverage data:
```bash
cat /shared5/Alex/Horses/horse_ID_list.txt | awk '{print $1}' | while read name; do
  qualimap bamqc -bam /shared5/Alex/Horses/BAMs/${name}_markdup.bam --java-mem-size=200G &
done
```

**7. Consensus FASTA file**
Prior to obtain the consensus FASTA file, I need to get the mean and maximum depth of each BAM file:
```bash
cat /shared5/Alex/Horses/horse_ID_list.txt | awk '{print $1}' | while read name; do
  coverage=$(cat /shared5/Alex/Horses/BAMs/${name}_markdup_stats/genome_results.txt | grep "mean coverageData"  | awk '{print $4}' | sed 's/X//')
paste <(echo ${name}) <(echo ${coverage}) <(echo "${coverage} * 2" | bc);
done >> horse_ID_depth_list.txt

#horse_ID_depth_list.txt content:
ERR1527947      26.6281 53.2562
ERR1305963      19.9776 39.9552
ERR1512897      18.5663 37.1326
```

Now I can run the consensus FASTA command:
```bash
cat /shared5/Alex/Horses/horse_ID_depth_list.txt | while read line ; do
  #defining variables
  name=$(echo ${line} | awk '{print $1}');
  coverage=$(echo ${line} | awk '{print $2}');
  maxcoverage=$(echo ${line} | awk '{print $3}');

  angsd -i /shared5/Alex/Horses/BAMs/${name}_markdup.bam -minMapQ 30 -minQ 20 -setMinDepth ${coverage} -setMaxDepth ${maxcoverage} -remove_bads 1 -doFasta 2 -doCounts 1 -ref /shared5/Alex/Donkey_ref/GCF_016077325.2_ASM1607732v2_genomic.fna -out /shared5/Alex/Horses/FASTAs/${name}_DonkeyRef &
  done
```



#### Przewalski FASTQs
The previous master student also had 4 Przewalski genomes located in one of his directores (/shared5/studentprojects/Qikai/Horse_database/Prezwalski). I will be processing these 4 sequences but only using 1 sequence (the largest) for the variant calling.

**0. Obtaining FASTQs**
1.	I make soft links of the FASTQs into my directory
```bash
mkdir /shared5/Alex/Przewalski/FASTQs/

for name in $(ls /shared5/studentprojects/Qikai/Horse_database/Prezwaiski/*val*) ; do
  ln -s /shared5/studentprojects/Qikai/Horse_database/Prezwaiski/${name} /shared5/Alex/Przewalski/FASTQs/ ;
done
```

**1. Trimming FASTQs**
The Przewlaski FASTQs are already trimmed so I skip this step.

**2. Mapping**
```bash
cat /shared5/Alex/Przewalski/Przewalski_list.txt | awk '{print $1}' | while read name; do
  bwa mem /shared5/Alex/Donkey_ref/GCF_016077325.2_ASM1607732v2_genomic.fna.gz /shared5/Alex/Przewalski/FASTQs/${name}_1_val_1.fq.gz /shared5/Alex/Przewalski/FASTQs/${name}_2_val_2.fq.gz > /shared5/Alex/Przewalski/BAMs/${name}.sam &
done
```


**3. Coverting to BAM and sorting by read names**
```bash
cat /shared5/Alex/Przewalski/Przewalski_list.txt | awk '{print $1}' | while read name; do
  samtools view -q 30 -u /shared5/Alex/Przewalski/BAMs/${name}.sam | samtools sort -n -o /shared5/Alex/Przewalski/BAMs/${name}_sorted.bam &
done
```

**4. Mark duplicates**
```bash
cat /shared5/Alex/Przewalski/Przewalski_list.txt | awk '{print $1}' | while read name; do
  samtools fixmate -m /shared5/Alex/Przewalski/BAMs/${name}_sorted.bam - | samtools sort - | samtools markdup - /shared5/Alex/Przewalski/BAMs/${name}_markdup.bam &
done
```

**5. Indexing**
```bash
cat /shared5/Alex/Przewalski/Przewalski_list.txt | awk '{print $1}' | while read name; do
  samtools index /shared5/Alex/Przewalski/BAMs/${name}_markdup.bam &
done
```


**6. Calculating BAM depth**
I used `qualimap` and this creates a folder that generates a file called "genome_results.txt" for each sequence with the coverage data:
```bash
cat /shared5/Alex/Przewalski/Przewalski_list.txt | awk '{print $1}' | while read name; do
    qualimap bamqc -bam /shared5/Alex/Przewalski/BAMs/${name}_markdup.bam --java-mem-size=200G &
done
```

**7. Consensus FASTA file**
Prior to obtain the consensus FASTA file, I need to get the mean and maximum depth of each BAM file:
```bash
cat /shared5/Alex/Przewalski/Przewalski_list.txt | awk '{print $1}' | while read name; do
  coverage=$(cat /shared5/Alex/Przewalski/BAMs/${name}_markdup_stats/genome_results.txt | grep "mean coverageData"  | awk '{print $4}' | sed 's/X//')
paste <(echo ${name}) <(echo ${coverage}) <(echo "${coverage} * 2" | bc);
done >> przewalski_ID_depth_list.txt
```

Now I can run the consensus FASTA command:
```bash
cat /shared5/Alex/Przewalski/przewalski_ID_depth_list.txt | while read line ; do
  #defining variables
  name=$(echo ${line} | awk '{print $1}');
  coverage=$(echo ${line} | awk '{print $2}');
  maxcoverage=$(echo ${line} | awk '{print $3}');

  angsd -i /shared5/Alex/Przewalski/BAMs/${name}_markdup.bam -minMapQ 30 -minQ 20 -setMinDepth ${coverage} -setMaxDepth ${maxcoverage} -remove_bads 1 -doFasta 2 -doCounts 1 -ref /shared5/Alex/Donkey_ref/GCF_016077325.2_ASM1607732v2_genomic.fna -out /shared5/Alex/Przewalski/FASTAs/${name}_DonkeyRef &
  done
```



## Variant calling using Strelka
I used `strelka` for variant calling of the genomes because it is quicker than `VCFtools` and is optimized for Illumina data.

1.	First I created a directory for the first run of Strelka.
```bash
mkdir Run1
```

2.	Then I created a Strelka script that follows the follwoing template :

```bash
#content of strelka script:
configureStrelkaGermlineWorkflow.py \
--bam <BAM>\
--bam <BAM2> \
--referenceFasta <Reference FASTA> \
--runDir .
```
So to create the script I compiled a list of the BAMs and their files. The final script looked like this:
```bash
cat configure_strelka.sh

# Output
configureStrelkaGermlineWorkflow.py \
--bam /shared5/Alex/Horses/BAMs/ERR1527947_markdup.bam \
--bam /shared5/Alex/Horses/BAMs/ERR1305963_markdup.bam \
--bam /shared5/Alex/Horses/BAMs/ERR1512897_markdup.bam \
--bam /shared5/Alex/Horses/BAMs/ERR2179551_markdup.bam \
--bam /shared5/Alex/Horses/BAMs/ERR1545178_markdup.bam \
--bam /shared5/Alex/Horses/BAMs/ERR1545179_markdup.bam \
--bam /shared5/Alex/Horses/BAMs/ERR978603_markdup.bam \
--bam /shared5/Alex/Horses/BAMs/ERR2179542_markdup.bam \
--bam /shared5/Alex/Horses/BAMs/ERR2179553_markdup.bam \
--bam /shared5/Alex/Horses/BAMs/ERR2179549_markdup.bam \
--bam /shared5/Alex/Horses/BAMs/ERR2731058_markdup.bam \
--bam /shared5/Alex/Horses/BAMs/ERR2731057_markdup.bam \
--bam /shared5/Alex/Horses/BAMs/ERR1527967_markdup.bam \
--bam /shared5/Alex/Horses/BAMs/ERR2179546_markdup.bam \
--bam /shared5/Alex/Horses/BAMs/ERR2179552_markdup.bam \
--bam /shared5/Alex/Horses/BAMs/ERR1545183_markdup.bam \
--bam /shared5/Alex/Horses/BAMs/ERR2179544_markdup.bam \
--bam /shared5/Alex/Horses/BAMs/ERR1527968_markdup.bam \
--bam /shared5/Alex/Horses/BAMs/ERR868004_markdup.bam \
--bam /shared5/Alex/Horses/BAMs/ERR1545189_markdup.bam \
--bam /shared5/Alex/Horses/BAMs/ERR1735862_markdup.bam \
--bam /shared5/Alex/Horses/BAMs/ERR2203766_markdup.bam \
--bam /shared5/Alex/Horses/BAMs/ERR1305961_markdup.bam \
--bam /shared5/Alex/Horses/BAMs/ERR2179543_markdup.bam \
--bam /shared5/Alex/Horses/BAMs/ERR1545186_markdup.bam \
--bam /shared5/Alex/Przewalski/BAMs/SRR12719743_markdup.bam \
--bam /shared5/Alex/Exmoor/BAMs/s2010_EDSW210003772_markdup.bam \
--bam /shared5/Alex/Donkey/BAMs/SRR7031483_markdup.bam \
--referenceFasta /shared5/Alex/Donkey_ref/GCF_016077325.2_ASM1607732v2_genomic.fna \
--runDir .
```

3.	The I ran the script:
```bash
source activate strelka
bash configure_strelka.sh
```

4.	This created a python script for called “runWorkflow.py”, which I then executed:
```bash
./runWorkflow.py -m local -j 15
```
* `-m`: *Run mode*
* `-j`: *Number of jobs*. Default is the estimate total cores on this node for local mode.



## VCF filtering (31/10/22)
Now that the VCF file has been created it is important to filter the variants. I use `VCFtools` for this.

1.	The first filter removes insertions and deletions and low quality positions and variants:
```bash
vcftools --gzvcf variants.vcf.gz --remove-indels --minQ 30 --minGQ 30 --out ../../variants_rmvIndels_minQ30_minGQ30 –recode
```
* `--gzvcf`: *Input file name of compressed vcf*
* `--remove-indels`: *Include or exclude sites that contain an indel.* For these options "indel" means any variant that alters the length of the REF allele.
* `--minQ`: *Includes only sites with Quality value above this threshold*
* `-minGQ`: *Exclude all genotypes with a quality below the threshold specified*
* `--out`: *Specify output file name and location*
* `-recode`: *Generates a new file with applied filters*

2.	The second filter is a minimum allele count and also removes non-neutrak alleles:
```bash
vcftools --vcf variants_rmvIndels_minQ30_minGQ30.recode.vcf --mac 3 --hwe 0.05 --remove-filtered-all --out variants_rmvIndels_minQ30_minGQ30_mac3_hwe0.05_PASSonly –recode
```
* `--vcf`: *Input file name*
* `--mac`: *Minimum allele count*. Minimum allele count set at 3, meaning that each allele must appear at least 3 times over all individuals for that site
* `--hwe`: *Assesses alleles for Hardy-Weinberg equilibrium.* Parameter is set 0.05 (i.e. p-value set at 0.05) so only retain alleles that are neutral.
* `--remove-filtered-all`: *Removes all sites with a FILTER flag other than PASS.*

3.	The 3rd round of filtering keeps only sites that have at least 80% of data (i.e. the site is found in at least 80%of individuals)
```bash
vcftools --vcf variants_rmvIndels_minQ30_minGQ30_mac3_hwe0.05_PASSonly.recode.vcf --max-missing 0.8 --out variants_rmvIndels_minQ30_minGQ30_mac3_hwe0.05_PASSonly_maxmissing0.8 --recode
```
* `--max-missing`: *Exclude sites on the basis of the proportion of missing data*. Defined to be between 0 and 1, where 0 allows sites that are completely missing and 1 indicates no missing data allowed.


4.	The 4th round of filtering removes the top and bottom 2.5% mean loci depth

We will remove the top and bottom 2.5% site mean depth because they are not reliable; the top 2.5% might have an artificially large depth and the bottom 3.5% have very low coverage).



First, VCFtools calculates the mean depth for each site:
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
vcftools --vcf variants_rmvIndels_minQ30_minGQ30_mac3_hwe0.05_PASSonly_maxmissing0.8.recode.vcf --min-meanDP 15.9 --max-meanDP 24.1 --out variants_rmvIndels_minQ30_minGQ30_mac3_hwe0.05_PASSonly_maxmissing0.8_meanDPmid95percentile --recode
```
* `min-meanDP`: *Includes only sites with mean depth values (over all included individuals) greater than or equal to the value.*
* `max-meanDP`: *Includes only sites with mean depth values (over all included individuals) less than or equal to the "--max-meanDP" value.*
