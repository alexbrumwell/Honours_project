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
First, we need to have a reference genome to map the FASTQs on to.

I will be using a horse U-chromosome that was provide dto me by my post-doc.

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

```bash
cat /shared5/Alex/Y-chromosome_project/List_ychromosome.txt | awk '{print $2}' | while read file ; do

  # Defining variables:
  location=$(echo ${file} | awk -F "/" '{$NF=""}1' OFS="/") ; # variable of file path
  name=$(echo ${file} | awk -F "/" '{print $NF}' | sed 's/_1_val_1.fq.gz//'); # variable of sample name

  # Mapping
  bwa mem /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/horse_Y_chromosome.fasta \ # reference genome in fasta format
  ${location}${name}_1_val_1.fq.gz \ # read 1
  ${location}${name}_2_val_2.fq.gz \ # read 2
  2> shared5/Alex/Y-chromosome_project/BAMs/log_files/${name}.bwamem.log \ # creates a log file containing errors
  > /shared5/Alex/Y-chromosome_project/BAMs/${name}_Y.sam &  # output file location and name
done
```
***Note on code formatting:*** *\"\\" at the end of a line means that the command continues in the following line of code. In other words, the five lines of code in the `bwa mem` command are usually written in one line but, in this case, it is more visual to see each argument on a different line*



## 2. Converting to BAM and sorting
Now that the SAM files have been created, I convert them to BAM format so they take less space and sort them by reading name using `samtools view` and then `samtools sort`:

```bash
cat /shared5/Alex/Y-chromosome_project/List_ychromosome.txt | awk '{print $1}' | while read name ; do
  samtools view -q 30 -u /shared5/Alex/Y-chromosome_project/BAMs/${name}_Y.sam | \ # "-q" Skip alignments with MAPQ
  samtools sort -n -o /shared5/Alex/Y-chromosome_project/BAMs/${name}_Y_sorted.bam & # this creates a log file containing standard errors
done
```
> `samtools view` options:
* `-q`: *Minimal mapping quality*. Basic quality filter that skips alignments with a mapping quality smaller than the set value (in this case 30)
* `-u`: *Output uncompressed data*. Preferred when piping to another samtools command.

> `samtools sort`options:
* `-n`: *Sort by read names rather than by chromosomal coordinates.*. We are sorting by read name for later steps
* `-o`: *Write the final sorted output to specified file name*



## 3. Merging BAM files
Since some of the sequences I have belong to the same individual, I have to merge the corresponding BAM files together.

`samtools merge` can be used for this:
```bash
samtools merge <output file> <input file1> <input file2>
```
n I create a new list (List_Y-chromosome_merged.txt) with the updated IDs of the sequences. E.g. The merged file IDs have now become "s479021_merged" instead of "s2479021_EDSW210003783-1a_H3WNKDSX2_L4"


## 4. Marking duplicates
I use `samtools markdup`, but prior to this, the files need to be processed with `samtools fixmate` and then sorted by coordinate using `samtools sort`:

```bash
cat /shared5/Alex/Y-chromosome_project/List_Y-chromosome_merged.txt | awk '{print $1}' | while read name ; do
  samtools fixmate -m /shared5/Alex/Y-chromosome_project/BAMs/${name}_Y_sorted.bam - | \
  samtools sort - | \
  samtools markdup - /shared5/Alex/Y-chromosome_project/BAMs/${name}_Y_sorted_markdup.bam &
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


## 5. Indexing
Now that the BAMs has been created, sorted and duplicates have been marked, they need to be indexed for future analyses.

```bash
cat /shared5/Alex/Y-chromosome_project/List_Y-chromosome_merged.txt | awk '{print $1}' | while read name ; do
  samtools index /shared5/Alex/Y-chromosome_project/BAMs/${name}_Y_sorted_markdup.bam &
done
```

## 6. Validating BAM files
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

  /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/horse_Y_chromosome.fasta
```


## 7. Calculating BAM coverage
I use `qualimap bamqc` to calculate the coverage of the BAM files [(qualimap bamqc manual)](http://qualimap.conesalab.org/doc_html/analysis.html#bamqc).

*Important to note:* There is no need to specify an output as qualimap creates an output folder containing the stats results.

```bash
cat /shared5/Alex/Y-chromosome_project/List_Y-chromosome_merged.txt | awk '{print $1}' | while read name ; do
  qualimap bamqc -bam /shared5/Alex/Y-chromosome_project/BAMs/${name}_Y_sorted_markdup.bam \
  --java-mem-size=200G 2>/shared5/Alex/Y-chromosome_project/BAMs/log_files/${name}_qualimap.log & # Sets desired memory size
done
```


## 8. Creating consensus FASTA files
In order to call variants with VCFtools, we first need to create consensus FASTA files using `angsd`:

#### 8.1 Creating sample list with coverage data
ANGSD needs the minimum and maxium coverage values for each BAM, so we will create a new sample list that contains the names of each sample and also the minimum and maxium coverage values for each BAM.

```bash
cat /shared5/Alex/Y-chromosome_project/List_Y-chromosome_merged.txt | awk '{print $1}' | while read name; do
  coverage=$(cat /shared5/Alex/Y-chromosome_project/BAMs/${name}*stats/genome_results.txt | grep "mean coverageData"  | awk '{print $4}' | sed 's/X//' | sed 's/,//' | cut -d "." -f 1);
  paste <(echo ${name}) <(echo ${coverage}) <(echo "${coverage} * 2" | bc);
done > /shared5/Alex/Y-chromosome_project/Y-chromosome_ID_depth_list.txt
```


#### 8.2 Consensus FASTA command

```bash
cat /shared5/Alex/Y-chromosome_project/Y-chromosome_ID_depth_list.txt | while read line; do

  #Defining variables:
  name=$(echo ${line} | awk '{print $1}');
  number=$(echo ${line} | awk '{print $(NF-1)}');
  maxnumber=$(echo ${line} | awk '{print $(NF)}');

  #consensus FASTA
  angsd -i /shared5/Alex/Y-chromosome_project/BAMs/${name}_Y_sorted_markdup.bam -minMapQ 30 -minQ 20 -setMinDepth ${number} -setMaxDepth ${maxnumber} -remove_bads 1 -doFasta 2 -doCounts 1 -ref /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/horse_Y_chromosome.fasta -out /shared5/Alex/Y-chromosome_project/FASTAs/${name}_Y &
done
```

#### 8.3 Changing headers
We add the sequence ID to each FASTA header:
```bash
#first uncompress the FASTAs
gunzip *fa.gz

#add corresponding breed to chromosome_ID_depth_list.txt
paste <(cat Y-chromosome_ID_depth_list.txt) <(echo "Exmoor
Exmoor
...") > tmp.txt && mv tmp.txt Y-chromosome_ID_depth_list.txt

#change header
cat /shared5/Alex/Y-chromosome_project/Y-chromosome_ID_depth_list.txt | while read line; do
  name=$(echo $line | awk '{print $1}');
  breed=$(echo $line | awk '{print $4}');
  sed -i "s/>MH341179.1/>${name}_${breed}_MH341179.1/g" /shared5/Alex/Y-chromosome_project/FASTAs/${name}_Y.fa
done
```

#### 8.4 Concatenating FASTAs
Then we compile all the FASTAs in one file and add the reference mitogenome as well:
```bash
#concatinating FASTAs
cat /shared5/Alex/Y-chromosome_project/Y-chromosome_ID_depth_list.txt | awk '{print $1}' | while read name; do
  cat /shared5/Alex/Y-chromosome_project/FASTAs/${name}_Y.fa
done >> Exmoor_horse_Y.fa


## adding reference sequence
cat /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/horse_Y_chromosome.fasta Exmoor_horse_Y.fa > tmp.fa && mv tmp.fa Exmoor_horse_Y.fa
```

#### 8.5 Copy concatinated FASTA file to local disk
```bash
scp studentprojects@young.eng.gla.ac.uk:/shared5/Alex/Y-chromosome_project/FASTAs/Exmoor_horse_Y.fa ./Desktop
```



## 9. Preparing NEXUS file
#### 9.1 Aligning sequences
I use Bioedit to view the sequences and align them.

#### 9.2 Count haplotypes
Count how many haplotypes there are by using a tree or DNAsp

#### 9.3 Create NEXUS file 
