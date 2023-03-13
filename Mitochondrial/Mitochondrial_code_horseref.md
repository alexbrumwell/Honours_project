# Code for mitochondrial analysis of Exmoor ponies using a horse mitogenome reference

---
* Title: Code for mitochodnrial analyses
* Author: Alexandra Brumwell
---

Code notebook for honours project on Exmoor pony populations genetics

--------------------------------------------------------------------------------------------
## Introduction
This is essentially the same pipeline as the previousmitochondrial analysis but using only the Exmoor sequences and mapping onto a horse mitogenome to get more detail.


## Downloading a reference genome
I dowdloaded the newest and highest quality horse mitogenome from NCBI (https://www.ncbi.nlm.nih.gov/nucleotide/NC_001640.1).



## FASTQtoBAM pipeline
### 0. FASTQC
Check the quality of the FASTQ files using FastQC.
```bash
cat /shared5/Alex/Exmoor_sequencing_data/all_exmoor_list.txt | awk '{print $2}' | while read file; do
fastqc ${file} -o /shared5/Alex/Exmoor_sequencing_data/FASTQC/
done
```




### 1. Mapping FASTQs
**Important to note**: The FASTQs I am working with are already trimmed.

I wrote a loop that goes over the list of samples we have and maps them to the donkey mitogenome using `bwa mem` :
```bash
cat /shared5/Alex/Mitochondrial_project/Mitochondrial_list.txt | grep "Exmoor" | awk '{print $2}' | while read file ; do

  # Defining variables:
  location=$(echo ${file} | awk -F "/" '{$NF=""}1' OFS="/") ; # variable of file path
  name=$(echo ${file} | awk -F "/" '{print $NF}' | sed 's/_1_val_1.fq.gz//'); # variable of sample name

  # Mapping
  bwa mem /shared5/Alex/Mitochondrial_project/horse_mtdna_reference/horse_mtdna_reference.fasta \ # reference genome in fasta format
  ${location}${name}_1_val_1.fq.gz \ # read 1
  ${location}${name}_2_val_2.fq.gz \ # read 2
  2> /shared5/Alex/Mitochondrial_project/BAMs_horseref/log_files/${name}.bwamem.log \ # creates a log file containing errors
  > /shared5/Alex/Mitochondrial_project/BAMs_horseref/${name}_mtdna_horseref.sam &  # output file location and name
done
```
***Note on code formatting:*** *\"\\" at the end of a line means that the command continues in the following line of code. In other words, the five lines of code in the `bwa mem` command are usually written in one line but, in this case, it is more visual to see each argument on a different line*



### 2. Converting to BAM and sorting
Now that the SAM files have been created, I convert them to BAM format so they take less space and sort them by reading name using `samtools view` and then `samtools sort`:

```bash
cat /shared5/Alex/Mitochondrial_project/Mitochondrial_list.txt | grep "Exmoor" | awk '{print $1}' | while read name ; do
  samtools view -q 30 -u /shared5/Alex/Mitochondrial_project/BAMs_horseref/${name}_mtdna_horseref.sam | \ # "-q" Skip alignments with MAPQ
  samtools sort -n -o /shared5/Alex/Mitochondrial_project/BAMs_horseref/${name}_sorted_mtdna_horseref.bam &
done
```
> `samtools view` options:
* `-q`: *Minimal mapping quality*. Basic quality filter that skips alignments with a mapping quality smaller than the set value (in this case 30)
* `-u`: *Output uncompressed data*. Preferred when piping to another samtools command.

> `samtools sort`options:
* `-n`: *Sort by read names rather than by chromosomal coordinates.*. We are sorting by read name for later steps
* `-o`: *Write the final sorted output to specified file name*



### 3. Merging BAM files
Since some of the sequences I have belong to the same individual, I have to merge the corresponding BAM files together.

I compile the IDs of the indivudals that are spearted into multiple files and use `samtools merge`:
```bash
samtools merge <output file> <input file1> input file2>

cat tmp.txt | while read name; do #tmp.txt contains the IDs of the files that need to be merged
  samtools merge ${name}_merged_sorted_mtdna_horseref.bam $(ls $name*sorted*);
done
```
Then I create a new list (List_mitochondrial_merged.txt) with the updated IDs of the sequences. E.g. The merged file IDs have now become "s479021_merged" instead of "s2479021_EDSW210003783-1a_H3WNKDSX2_L4"


### 4. Marking duplicates
I use `samtools markdup`, but prior to this, the files need to be processed with `samtools fixmate` and then sorted by coordinate using `samtools sort`:

```bash
cat /shared5/Alex/Mitochondrial_project/List_mitochondrial_merged.txt | grep "Exmoor" | awk '{print $1}' | while read name ; do
  samtools fixmate -m /shared5/Alex/Mitochondrial_project/BAMs_horseref/${name}_sorted_mtdna_horseref.bam - | \
  samtools sort - | \
  samtools markdup - /shared5/Alex/Mitochondrial_project/BAMs_horseref/${name}_markdup_sorted_mtdna_horseref.bam &
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


### 5. Indexing
Now that the BAMs has been created, sorted and duplicates have been marked, they need to be indexed for future analyses.

```bash
cat /shared5/Alex/Mitochondrial_project/List_mitochondrial_merged.txt | grep "Exmoor" | awk '{print $1}' | while read name ; do
  samtools index /shared5/Alex/Mitochondrial_project/BAMs_horseref/${name}_markdup_sorted_mtdna_horseref.bam &
done
```

### 6. Validating BAM files
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

cat /shared5/Alex/Mitochondrial_project/List_mitochondrial_merged.txt | grep "Exmoor" | awk '{print $1}' | while read name ; do
  mkdir -p /shared5/Alex/Mitochondrial_project/BAMs_horseref/ValidateSamFile/
  java -jar /shared5/Alex/picard-2.27.5/picard.jar ValidateSamFile -I /shared5/Alex/Mitochondrial_project/BAMs_horseref/${name}_mtdna_sorted_markdup.bam -O /shared5/Alex/Mitochondrial_project/BAMs_horseref/ValidateSamFile/${name}_mtdna_sorted_markdup_errors.txt -R /shared5/Alex/Donkey_ref/Donkey_mitogenome.fasta --MODE SUMMARY
done
```


### 7. Calculating BAM coverage
I use `qualimap bamqc` to calculate the coverage of the BAM files [(qualimap bamqc manual)](http://qualimap.conesalab.org/doc_html/analysis.html#bamqc).

*Important to note:* There is no need to specify an output as qualimap creates an output folder containing the stats results.

```bash
cat /shared5/Alex/Mitochondrial_project/List_mitochondrial_merged.txt | grep "Exmoor" | awk '{print $1}' | while read name ; do
  qualimap bamqc -bam /shared5/Alex/Mitochondrial_project/BAMs_horseref/${name}_markdup_sorted_mtdna_horseref.bam \
  --java-mem-size=200G \ #Sets desired memory size
  2>/shared5/Alex/Mitochondrial_project/BAMs_horseref/log_files/${name}_qualimap.log & #writes log file
done
```

To check the coverage of each file:
```bash
cat /shared5/Alex/Mitochondrial_project/List_mitochondrial_merged.txt | grep "Exmoor" | awk '{print $1}' | while read name; do
  paste <(echo ${name}) <(grep "mean coverageData" /shared5/Alex/Mitochondrial_project/BAMs_horseref/${name}*stats/genome_results.txt) 
done
```
Coverage data:
E_102004_EKDN220034353-1A_H5GL5DSX5_L2       mean coverageData = 1,098.599X
E_107013_EKDN220034352-1A_H5GL5DSX5_L2       mean coverageData = 737.6488X
E_21084_EKDN220034367-1A_H5J5MDSX5_L2        mean coverageData = 417.2181X
E_21149_EKDN220034359-1A_H5J5MDSX5_L2        mean coverageData = 738.6468X
E_23279_merged       mean coverageData = 1,982.7404X
E_23416_merged       mean coverageData = 1,429.9435X
E_23434_merged       mean coverageData = 2,204.9577X
E_320005_merged      mean coverageData = 1,714.2256X
E_32023_merged       mean coverageData = 2,064.4987X
E_44009_EKDN220034351-1A_H5GL5DSX5_L2        mean coverageData = 1,433.768X
E_458028_EKDN220034371-1A_H5GL5DSX5_L3       mean coverageData = 1,893.2809X
E_479023_merged      mean coverageData = 1,887.0448X
E_49031_EKDN220034370-1A_H5J5MDSX5_L2        mean coverageData = 1,246.0993X
E_49057_EKDN220034360-1A_H5GL5DSX5_L3        mean coverageData = 938.8903X
E_49121_EKDN220034366-1A_H5J5MDSX5_L2        mean coverageData = 1,602.552X
E_512001_merged      mean coverageData = 1,788.7761X
E_519005_EKDN220034364-1A_H5J5MDSX5_L2       mean coverageData = 1,706.5173X
E_78170_merged       mean coverageData = 1,965.0745X
E_900234_EKDN220034355-1A_H5J5MDSX5_L2       mean coverageData = 2,315.0523X
E_900585_EKDN220034354-1A_H5GL5DSX5_L2       mean coverageData = 2,153.8747X
E_900588_merged      mean coverageData = 1,796.6345X
E_900600_EKDN220034356-1A_H5GL5DSX5_L2       mean coverageData = 1,859.6243X
E_900694_EKDN220034350-1A_H5GL5DSX5_L2       mean coverageData = 2,318.5524X
E_900717_merged      mean coverageData = 2,752.2918X
E_900741_EKDN220034358-1A_H5J5MDSX5_L2       mean coverageData = 2,326.7842X
s12157_EDSW210003789-1a_H3WNKDSX2_L2         mean coverageData = 2,046.4336X
s2010_EDSW210003772-1a_H3WNKDSX2_L3          mean coverageData = 1,256.5758X
s2012_EDSW210003766-1a_H3WHWDSX2_L4          mean coverageData = 1,580.1834X
s21098_EDSW210003783-1a_H3WNKDSX2_L4         mean coverageData = 2,183.2298X
s21131_EDSW210003778-1a_H3WNKDSX2_L1         mean coverageData = 1,779.1316X
s21134_EDSW210003784-1a_H3WNKDSX2_L4         mean coverageData = 1,515.9487X
s23203_EDSW210003785-1a_H3WNKDSX2_L1         mean coverageData = 2,755.9809X
s235016_EDSW210003786-1a_H3WNKDSX2_L2        mean coverageData = 2,947.5282X
s237009_EDSW210003788-1a_H3WNKDSX2_L2        mean coverageData = 1,038.2477X
s276021_EDSW210003781-1a_H3WNKDSX2_L1        mean coverageData = 2,326.3255X
s335004_EDSW210003787-1a_H3WNKDSX2_L2        mean coverageData = 2,428.2085X
s479013_EDSW210003775-1a_H3WNKDSX2_L3        mean coverageData = 1,380.4795X
s479021_merged       mean coverageData = 1,844.7695X
s479022_EDSW210003782-1a_H3WNKDSX2_L4        mean coverageData = 1,739.3239X
s49011_EDSW210003773-1a_H3FTWDSX2_L1         mean coverageData = 1,804.4637X
s49022_EDSW210003776-1a_H3WNKDSX2_L3         mean coverageData = 407.9222X
s49052_EDSW210003777-1a_H3WNKDSX2_L3         mean coverageData = 1,315.4797X
s49097_EDSW210003769-1a_H3FTWDSX2_L2         mean coverageData = 1,573.3412X
s49124_EDSW210003767-1a_H3WHWDSX2_L2         mean coverageData = 1,700.6455X
s49125_EDSW210003779-1a_H3WNKDSX2_L1         mean coverageData = 2,435.7068X
s49127_merged        mean coverageData = 1,939.2587X
s49145_EDSW210003768-1a_H3WHWDSX2_L4         mean coverageData = 1,849.9449X
s49147_EDSW210003770-1a_H3FTWDSX2_L2         mean coverageData = 3,155.3542X
s49149_EDSW210003771-1a_H3FTWDSX2_L2         mean coverageData = 3,954.0727X
s78150_EDSW210003780-1a_H3WNKDSX2_L4         mean coverageData = 1,206.9502X
S21029_EDSW200015548-1a_HJFWVDSXY_L4         mean coverageData = 2,505.5127X
S21092_EDSW200015546-1a_HJFWVDSXY_L4         mean coverageData = 2,012.5363X
S21109_EDSW200015550-1a_HJFWVDSXY_L4         mean coverageData = 1,831.8222X
S21140_EDSW200015549-1a2a_HJFWVDSXY_L3L4             mean coverageData = 882.2166X
S276023_EDSW200015547-1a2a_HJFWVDSXY_L3L4            mean coverageData = 1,533.1698X
S49052_EDSW200015551-1a_HJFWVDSXY_L4         mean coverageData = 133.3785X

#To check size of the files:
```bash
cat /shared5/Alex/Mitochondrial_project/List_mitochondrial_merged.txt | grep "Exmoor" | awk '{print $1}' | while read name; do
  paste <(echo ${name}) <(ll -h ${name}*sorted_mtdna_horseref.bam | grep -v "markdup" | awk '{print $5}') <(ll -h ${name}*markdup_sorted_mtdna_horseref.bam | awk '{print $5}') <(grep "mean coverageData" /shared5/Alex/Mitochondrial_project/BAMs_horseref/${name}*stats/genome_results.txt | awk '{print $4}')
done
```

### 8. Creating consensus FASTA files
In order to call variants with VCFtools, we first need to create consensus FASTA files using `angsd`:

#### 8.1 Creating sample list with coverage data
ANGSD needs the minimum and maxium coverage values for each BAM, so we will create a new sample list that contains the names of each sample and also the mean and maxium coverage values for each BAM.

```bash
cat /shared5/Alex/Mitochondrial_project/List_mitochondrial_merged.txt | grep "Exmoor" | awk '{print $1}' | while read name; do
  coverage=$(cat /shared5/Alex/Mitochondrial_project/BAMs_horseref/${name}*stats/genome_results.txt | grep "mean coverageData"  | awk '{print $4}' | sed 's/X//' | sed 's/,//' | cut -d "." -f 1);
  paste <(echo ${name}) <(echo ${coverage}) <(echo "${coverage} * 2" | bc);
done > /shared5/Alex/Mitochondrial_project/Mitochondrial_ID_depth_list_horseref.txt
```

#### 8.2 Consensus FASTA command
IMPORTANT: I THINK I RERAN THIS BUT WITH MINIMUM DEPTH SET AT 100x
```bash
cat /shared5/Alex/Mitochondrial_project/Mitochondrial_ID_depth_list_horseref.txt | while read line; do

  #Defining variables:
  name=$(echo ${line} | awk '{print $1}');
  number=$(echo ${line} | awk '{print $(NF-1)}');
  maxnumber=$(echo ${line} | awk '{print $(NF)}');

  #consensus FASTA
  angsd -i /shared5/Alex/Mitochondrial_project/BAMs_horseref/${name}_markdup_sorted_mtdna_horseref.bam -minMapQ 30 -minQ 20 -setMinDepth ${number} -setMaxDepth ${maxnumber} -remove_bads 1 -doFasta 2 -doCounts 1 -ref /shared5/Alex/Mitochondrial_project/horse_mtdna_reference/horse_mtdna_reference.fasta -out /shared5/Alex/Mitochondrial_project/FASTAs_horseref/${name}_mtdna_horseref &
done
```

#### 8.3 Changing headers
We add the sequence ID to each FASTA header:
```bash
#first we unzip the compressed FASTAs
gunzip *fa.gz

#change headers
cat /shared5/Alex/Mitochondrial_project/Mitochondrial_ID_depth_list_horseref.txt | awk '{print $1}' | while read name; do
  sed -i "s/>NC_001640.1/>${name}_>NC_001640.1/g" /shared5/Alex/Mitochondrial_project/FASTAs_horseref/${name}_mtdna_horseref.fa
done
```

#### 8.4 Concatenating FASTAs
Then we compile all the FASTAs in one file and add the reference mitogenome as well:
```bash
#concatinating FASTAs
cat /shared5/Alex/Mitochondrial_project/Mitochondrial_ID_depth_list_horseref.txt | awk '{print $1}' | while read name; do
  cat /shared5/Alex/Mitochondrial_project/FASTAs_horseref/${name}_mtdna_horseref.fa
done >> Exmoor_mtdna_horseref.fa


## adding reference mitogenome
cat /shared5/Alex/Mitochondrial_project/horse_mtdna_reference/horse_mtdna_reference.fasta Exmoor_mtdna_horseref.fa > tmp.fa && mv tmp.fa Exmoor_mtdna_horseref.fa
```


## 9. Preparing NEXUS file
#### 9.1 Preparing (and aligning) sequences
I use Bioedit to view the sequences and check them.




#### 9.2 Determining haplotypes
I use DNAsp to determine how many different haplotypes there are. I do this by inputting the created FASTA file with only the Exmoor sequences and then I press "Generate > Haplotype data file". This creates a NEXUS file which I then have to add the TRAITS matrix to.

I then renumber the haplotypes using the numbering by Debbie.



#### 9.3 Create NEXUS file
I am going to create two NEXUS files depending on the TRAITS (i.e. population). First, I will create one based on the herd and the second on the maternal found of each sample. That way I can see how many haplotypes are found in each herd and maternal lineage.

I use EXCEL to add the traits part.
###### Herds
There are 20 herds.


###### Maternal lines
I create a text file which has all the possible maternal lines and then sort them and count how many there are using:
```bash
awk '{print $1}' all_maternal_lines.txt | sort | uniq -c
```
There are 23 different female founders.

## 10. Haplotype network
I use PopART for this.


## PCA
See R code
