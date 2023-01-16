# Code for mitochondrial analysis of Exmoor ponies

---
* Title: Code for mitochodnrial analyses
* Author: Alexandra Brumwell
---

Code notebook for honours project on Exmoor pony populations genetics

--------------------------------------------------------------------------------------------
## Introduction
I will be studying mitochondrial DNA across all the Exmoor pony sequences, other horse breeds, Przewalskis and donkeys. The sequences I am have are: 56 Exmoor sequences, 54 horse sequences (maximum 5 genomes per breed), 4 Przewalski sequences and 1 donkey sequence (different to the reference donkey mitogenome).

Prior to starting any data processing or analyses I create a spreadsheet containing the metadata of the sequences like genome ID, date of sequencing, breed, file name, sex, etc. I save all this information in a [spreadsheet](ADD LINK TO SPREADSHEET IN GOOGLEDRIVE) in Excel or Googledrive.

##### Example of metadata table:

| Individual | File_name | Breed | Pony number | Pony name | Sequencing_date | Sex | Notes | Original_FASTQ_file_path |
|------------|-----------|-------|-------------|-----------|-----------------|-----|-------|--------------------------|
| ERR111 | ERR111_EYB1000 | Exmoor | 111/111 | Cheddar | Mar2020 | F | From unknown maternal line | /Path/to/FASTQ/file |



## Downloading a reference genome
First, we need to have a reference genome to map the FASTQs on to.

I will be using a donkey mitogenome as a reference because I am interested in looking at mtDNA diversity within Exmoors and compared to other horse breeds. Therefore I need to map the sequences to an outgroup (i.e. donkey).

I downloaded a high-quality donkey mitogenome from NCBI [(NCBI Reference Sequence: NC_001788.1)](https://www.ncbi.nlm.nih.gov/nucleotide/NC_001788.1). I downloaded the sequence mannually from NCBI and then uploaded them to my cluster but there are commands like `wget` or `ncbi-genome-download` that can download the sequence directly to a cluster via the command line. I have read that *NCBI Entrez Direct UNIX E-utilities* are also useful ot download sequences via the command line: (e.g. `esearch -db nucleotide -query "NC_001788.1" | efetch -format fasta > test.fasta` )


## Creation of sample list
The FASTQ to BAM pipeline I use loops through sample IDs, so first I create a text file with the all the samples I will use and save this text file where I will be working in. This text file contains two columns: the first is the name of the sample and the second is the path to the FASTQ file of the sample.

To create this text file I usually go to the directory where I have all the FASTQ files and then I use `paste` and `ls` to compile the file names and file paths.

##### Example of creating list of Exmoor samples from Sep2020:
```bash
paste <(ls -l /shared5/Alex/Exmoor_sequencing_data/Sep_2020/*val_1.fq.gz | awk '{print $9}' | awk -F "/" '{print $NF}' | sed 's/_1_val_1.fq.gz//' ) <(ls -l /shared5/Alex/Exmoor_sequencing_data/Sep_2020/*val_1.fq.gz | awk '{print $9}')

#Output:
S21029_EDSW200015548-1a_HJFWVDSXY_L4    /shared5/Alex/Exmoor_sequencing_data/Sep_2020/S21029_EDSW200015548-1a_HJFWVDSXY_L4_1_val_1.fq.gz
S21092_EDSW200015546-1a_HJFWVDSXY_L4    /shared5/Alex/Exmoor_sequencing_data/Sep_2020/S21092_EDSW200015546-1a_HJFWVDSXY_L4_1_val_1.fq.gz
S21109_EDSW200015550-1a_HJFWVDSXY_L4    /shared5/Alex/Exmoor_sequencing_data/Sep_2020/S21109_EDSW200015550-1a_HJFWVDSXY_L4_1_val_1.fq.gz
```



## FASTQtoBAM pipeline
### 1. Mapping FASTQs
**Important to note**: The FASTQs I am working with are already trimmed.

I wrote a loop that goes over the list of samples we have and maps them to the donkey mitogenome using `bwa mem` :
```bash
cat /shared5/Alex/Mitochondrial_project/Mitochondrial_list.txt | awk '{print $2}' | while read file ; do

  # Defining variables:
  location=$(echo ${file} | awk -F "/" '{$NF=""}1' OFS="/") ; # variable of file path
  name=$(echo ${file} | awk -F "/" '{print $NF}' | sed 's/_1_val_1.fq.gz//'); # variable of sample name

  # Mapping
  bwa mem /shared5/Alex/Donkey_ref/Donkey_mitogenome.fasta \ # reference genome in fasta format
  ${location}${name}_1_val_1.fq.gz \ # read 1
  ${location}${name}_2_val_2.fq.gz \ # read 2
  2> /shared5/Alex/Mitochondrial_project/BAMs/log_files/${name}.bwamem.log \ # creates a log file containing errors
  > /shared5/Alex/Mitochondrial_project/BAMs/${name}_mtdna.sam &  # output file location and name
done
```
***Note on code formatting:*** *\"\\" at the end of a line means that the command continues in the following line of code. In other words, the five lines of code in the `bwa mem` command are usually written in one line but, in this case, it is more visual to see each argument on a different line*

\
If running the loop for a subset of data (e.g. only the mar2021 and sep2020 exmoor samples), then `grep` can be used in the first line to select samples of interest:
```bash
cat /shared5/Alex/Mitochondrial_project/Mitochondrial_list.txt | grep 'Mar_2021\|Sep_2020' | awk '{print $2}' |  while read file ; do

  # Defining variables:
  location=$(echo ${file} | awk -F "/" '{$NF=""}1' OFS="/") ; # variable of file path
  name=$(echo ${file} | awk -F "/" '{print $NF}' | sed 's/_1_val_1.fq.gz//'); # variable of sample name

  # Mapping
  bwa mem /shared5/Alex/Donkey_ref/Donkey_mitogenome.fasta \ #reference genome in fasta format
  ${location}${name}_1_val_1.fq.gz \ # read 1
  ${location}${name}_2_val_2.fq.gz \ # read 2
  2> ${location}${name}.bwamem.log # creates a log file containing errors
  > /shared5/Alex/Mitochondrial_project/BAMs/${name}_mtdna.sam &  # Output file location and name
done
```


### 2. Converting to BAM and sorting
Now that the SAM files have been created, I convert them to BAM format so they take less space and sort them by reading name using `samtools view` and then `samtools sort`:

```bash
cat /shared5/Alex/Mitochondrial_project/Mitochondrial_list.txt | awk '{print $1}' | while read name ; do
  samtools view -q 30 -u /shared5/Alex/Mitochondrial_project/BAMs/${name}_mtdna.sam | \ # "-q" Skip alignments with MAPQ
  samtools sort -n -o /shared5/Alex/Mitochondrial_project/BAMs/${name}_mtdna_sorted.bam \
  2> /shared5/Alex/Mitochondrial_project/BAMs/log_files/${name}_mtdna_sorted.log & # this creates a log file containing standard errors
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
  samtools merge ${name}_mtdna_sorted_merged.bam $(ls $name*sorted*);
done
```
Then I create a new list (List_mitochondrial_merged.txt) with the updated IDs of the sequences. E.g. The merged file IDs have now become "s479021_merged" instead of "s2479021_EDSW210003783-1a_H3WNKDSX2_L4"


### 4. Marking duplicates
I use `samtools markdup`, but prior to this, the files need to be processed with `samtools fixmate` and then sorted by coordinate using `samtools sort`:

```bash
cat /shared5/Alex/Mitochondrial_project/List_mitochondrial_merged.txt | awk '{print $1}' | while read name ; do
  samtools fixmate -m /shared5/Alex/Mitochondrial_project/BAMs/${name}_mtdna_sorted.bam - | \
  samtools sort - | \
  samtools markdup - /shared5/Alex/Mitochondrial_project/BAMs/${name}_mtdna_sorted_markdup.bam &
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
cat /shared5/Alex/Mitochondrial_project/List_mitochondrial_merged.txt | awk '{print $1}' | while read name ; do
  samtools index /shared5/Alex/Mitochondrial_project/BAMs/${name}_mtdna_sorted_markdup.bam &
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

cat /shared5/Alex/Mitochondrial_project/List_mitochondrial_merged.txt | awk '{print $1}' | while read name ; do
  mkdir -p /shared5/Alex/Mitochondrial_project/BAMs/ValidateSamFile/
  java -jar /shared5/Alex/picard-2.27.5/picard.jar ValidateSamFile -I /shared5/Alex/Mitochondrial_project/BAMs/${name}_mtdna_sorted_markdup.bam -O /shared5/Alex/Mitochondrial_project/BAMs/ValidateSamFile/${name}_mtdna_sorted_markdup_errors.txt -R /shared5/Alex/Donkey_ref/Donkey_mitogenome.fasta --MODE SUMMARY
done
```


### 7. Calculating BAM coverage
I use `qualimap bamqc` to calculate the coverage of the BAM files [(qualimap bamqc manual)](http://qualimap.conesalab.org/doc_html/analysis.html#bamqc).

*Important to note:* There is no need to specify an output as qualimap creates an output folder containing the stats results.

```bash
cat /shared5/Alex/Mitochondrial_project/List_mitochondrial_merged.txt | awk '{print $1}' | while read name ; do
  qualimap bamqc -bam /shared5/Alex/Mitochondrial_project/BAMs/${name}_mtdna_sorted_markdup.bam \
  --java-mem-size=200G \ #Sets desired memory size
  2>/shared5/Alex/Mitochondrial_project/BAMs/log_files/${name}_qualimap.log & #writes log file
done
```


### 8. Creating consensus FASTA files
In order to call variants with VCFtools, we first need to create consensus FASTA files using `angsd`:

#### 8.1 Creating sample list with coverage data
ANGSD needs the minimum and maxium coverage values for each BAM, so we will create a new sample list that contains the names of each sample and also the mean and maxium coverage values for each BAM.

```bash
cat /shared5/Alex/Mitochondrial_project/List_mitochondrial_merged.txt | awk '{print $1}' | while read name; do
  coverage=$(cat /shared5/Alex/Mitochondrial_project/BAMs/${name}*stats/genome_results.txt | grep "mean coverageData"  | awk '{print $4}' | sed 's/X//' | sed 's/,//' | cut -d "." -f 1);
  paste <(echo ${name}) <(echo ${coverage}) <(echo "${coverage} * 2" | bc);
done > /shared5/Alex/Mitochondrial_project/Mitochondrial_ID_depth_list.txt
```

#### 8.2 Consensus FASTA command

```bash
cat /shared5/Alex/Mitochondrial_project/Mitochondrial_ID_depth_list.txt | while read line; do

  #Defining variables:
  name=$(echo ${line} | awk '{print $1}');
  number=$(echo ${line} | awk '{print $(NF-1)}');
  maxnumber=$(echo ${line} | awk '{print $(NF)}');

  #consensus FASTA
  angsd -i /shared5/Alex/Mitochondrial_project/BAMs/${name}_mtdna_sorted_markdup.bam -minMapQ 30 -minQ 20 -setMinDepth ${number} -setMaxDepth ${maxnumber} -remove_bads 1 -doFasta 2 -doCounts 1 -ref /shared5/Alex/Donkey_ref/Donkey_mitogenome.fasta -out /shared5/Alex/Mitochondrial_project/FASTAs/${name}_mtdna_DonkeyRef &
done
```

#### 8.3 Changing headers
We add the sequence ID to each FASTA header:
```bash
cat /shared5/Alex/Mitochondrial_project/List_mitochondrial_merged.txt | awk '{print $1}' | while read name; do
  sed -i "s/>NC_001788.1/>${name} NC_001788.1/g" /shared5/Alex/Mitochondrial_project/FASTAs/${name}_mtdna_DonkeyRef.fa
done
```

#### 8.4 Concatenating FASTAs
Then we compile all the FASTAs in one file and add the reference mitogenome as well:
```bash
#concatinating FASTAs
cat /shared5/Alex/Mitochondrial_project/List_mitochondrial_merged.txt | awk '{print $1}' | while read name; do
  cat /shared5/Alex/Mitochondrial_project/FASTAs/${name}_mtdna_DonkeyRef.fa
done >> horse_mtdna_DonkeyRef.fa


## adding reference mitogenome
cat /shared5/Alex/Donkey_ref/Donkey_mitogenome.fasta >> horse_mtdna_DonkeyRef.fa

```

## 9. Preparing NEXUS file
#### 9.1 Preparing (and aligning) sequences
I use Bioedit to view the sequences and align them.
Comments:
- Przewalski mitogenomes are all uncertain; they are all "N"s
- Exmoor sequence S49052 and Franches-Montagne sequence ERR978597 contain no information
- Certain regions like th D-loop contain no information and are unusable

I copy only donkey reference and Exmoor sequences in a new file called "Exmoor_mtdna_DonkeyRef.fa". I removed sequence S49052 from this file. Removed from position 16110 onwards because it is a hypervariable region and there were amny uncertain bases.



#### 9.2 Determining haplotypes
I use DNAsp to determine how many different haplotypes there are. I do this by inputting the created FASTA file with only the Exmoor sequences and then I press "Generate > Haplotype data file". This creates a NEXUS file which I then have to add the TRAITS matrix to.

#### 9.3 Create NEXUS file
I am going to create two NEXUS files depending on the TRAITS (i.e. population). First, I will create one based on the herd and the second on the maternal found of each sample. That way I can see how many haplotypes are found in each herd and maternal lineage.

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
