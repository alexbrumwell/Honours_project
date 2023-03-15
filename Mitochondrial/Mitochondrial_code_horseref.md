# Code for mitochondrial analysis of Exmoor ponies using a horse mitogenome reference

---
* Title: Code for mitochondrial analyses
* Author: Alexandra Brumwell
---

Code notebook for honours project on Exmoor pony populations genetics

--------------------------------------------------------------------------------------------
## Introduction
I will be studying mitochondrial DNA of 55 Exmoor pony sequences.

Prior to starting any data processing or analyses I create a spreadsheet containing the metadata of the sequences like genome ID, date of sequencing, breed, file name, sex, etc. I save all this information in a [spreadsheet](ADD LINK TO SPREADSHEET IN GOOGLEDRIVE) in Excel or Googledrive.

##### Example of metadata table:

| Sequence ID | File_name | Pony number | Pony name | Sequencing_date | Sex | Notes | Original_FASTQ_file_path |
|------------|-----------|-------------|-----------|-----------------|-----|-------|--------------------------|
| ERR111 | ERR111_EYB1000 | 111/111 | Cheddar | Mar2020 | F | From unknown maternal line | /Path/to/FASTQ/file |



## Downloading a reference genome
I dowdloaded the newest and highest quality horse reference mitogenome (NCBI reference: NC_001640.1).


## Creation of sample list
The FASTQ to BAM pipeline I use loops through sample IDs, so first I create a text file with the all the samples I will use and save this text file where I will be working in. This text file contains two columns: the first is the name of the sample and the second is the path to the FASTQ file of the sample.

To create this text file I usually go to the directory where I have all the FASTQ files and then I use `paste` and `ls` to compile the file names and file paths.

##### Example of creating list of Exmoor samples:
```bash
paste <(ls -l /shared5/Alex/Exmoor_sequencing_data/Sep_2020/*val_1.fq.gz | awk '{print $9}' | awk -F "/" '{print $NF}' | sed 's/_1_val_1.fq.gz//' ) <(ls -l /shared5/Alex/Exmoor_sequencing_data/Sep_2020/*val_1.fq.gz | awk '{print $9}')

#Output:
S21029_EDSW200015548-1a_HJFWVDSXY_L4    /shared5/Alex/Exmoor_sequencing_data/Sep_2020/S21029_EDSW200015548-1a_HJFWVDSXY_L4_1_val_1.fq.gz
S21092_EDSW200015546-1a_HJFWVDSXY_L4    /shared5/Alex/Exmoor_sequencing_data/Sep_2020/S21092_EDSW200015546-1a_HJFWVDSXY_L4_1_val_1.fq.gz
S21109_EDSW200015550-1a_HJFWVDSXY_L4    /shared5/Alex/Exmoor_sequencing_data/Sep_2020/S21109_EDSW200015550-1a_HJFWVDSXY_L4_1_val_1.fq.gz
```


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
Since some of the BAM files I have come from the same indivdual, I merge them together. belong to the same individual, I have to merge the corresponding BAM files together.

I compile the IDs of the indivudals that are spearted into multiple files and use `samtools merge`:
```bash
samtools merge <output file> <input file1> input file2>
```
Then I create a new list (List_mitochondrial_merged.txt) with the updated IDs of the sequences. E.g. The merged file IDs have now become "s479021_merged" instead of "s2479021_EDSW210003783-1a_H3WNKDSX2_L4" and "s2479021_EDSW210003783-1a_H3WNKDSX2_L3"


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
Now that the BAMs has been created, sorted and duplicates have been marked, they need to be indexed.

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
#Code for samtools quickcheck:
samtools quickcheck <input.bam>

#Code for ValidateSamFile:
java -jar /shared5/Alex/picard-2.27.5/picard.jar ValidateSamFile \
  -I <input.bam> \ # input file
  -O <output.textfile> \ # output textfile containing any errors
  -R <reference_sequence> \ # optional
  --MODE SUMMARY # outputs a summary table listing the numbers of all 'errors' and 'warnings'.

cat /shared5/Alex/Mitochondrial_project/List_mitochondrial_merged.txt | grep "Exmoor" | awk '{print $1}' | while read name ; do
  mkdir -p /shared5/Alex/Mitochondrial_project/BAMs_horseref/ValidateSamFile/
  java -jar /shared5/Alex/picard-2.27.5/picard.jar ValidateSamFile \
  -I /shared5/Alex/Mitochondrial_project/BAMs_horseref/${name}_mtdna_sorted_markdup.bam \
  -O /shared5/Alex/Mitochondrial_project/BAMs_horseref/ValidateSamFile/${name}_mtdna_sorted_markdup_errors.txt \
  -R /shared5/Alex/Mitochondrial_project/horse_mtdna_reference/horse_mtdna_reference.fasta \
  --MODE SUMMARY
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


### 8. Creating consensus FASTA files

#### 8.1 Creating sample list with coverage data
I used ANGSD for create consensus sequences. ANGSD needs the minimum and maxium coverage values for each BAM, so I created  a new sample list that contained the names of each sample and also the mean and maxium coverage values for each BAM.

```bash
cat /shared5/Alex/Mitochondrial_project/List_mitochondrial_merged.txt | grep "Exmoor" | awk '{print $1}' | while read name; do
  coverage=$(cat /shared5/Alex/Mitochondrial_project/BAMs_horseref/${name}*stats/genome_results.txt | grep "mean coverageData"  | awk '{print $4}' | sed 's/X//' | sed 's/,//' | cut -d "." -f 1);
  paste <(echo ${name}) <(echo ${coverage}) <(echo "${coverage} * 2" | bc);
done > /shared5/Alex/Mitochondrial_project/Mitochondrial_ID_depth_list_horseref.txt
```

#### 8.2 Consensus FASTA command
```bash
cat /shared5/Alex/Mitochondrial_project/Mitochondrial_ID_depth_list_horseref.txt | while read line; do

  #Defining variables:
  name=$(echo ${line} | awk '{print $1}');
  number=$(echo ${line} | awk '{print $(NF-1)}');
  maxnumber=$(echo ${line} | awk '{print $(NF)}');

  #consensus FASTA
  angsd -i /shared5/Alex/Mitochondrial_project/BAMs_horseref/${name}_markdup_sorted_mtdna_horseref.bam \ #input file
  -minMapQ 30 \ #minimum mapping quality
  -minQ 20 \ #minimum read quality
  -setMinDepth ${number} \ #minimum depth
  -setMaxDepth ${maxnumber} \ #maximum depth
  -remove_bads 1 \ #Discards 'bad' reads
  -doFasta 2 \ #creates FASTA using the most common (non N) base
  -doCounts 1 \ #Counts the number of A,C,G,T.
  -ref /shared5/Alex/Mitochondrial_project/horse_mtdna_reference/horse_mtdna_reference.fasta \ #reference file
  -out /shared5/Alex/Mitochondrial_project/FASTAs_horseref/${name}_mtdna_horseref & #output file
done
```

#### 8.3 Changing headers
I added the sequence ID to each FASTA file header:
```bash
#first unzip the compressed FASTAs
gunzip *fa.gz

#change headers
cat /shared5/Alex/Mitochondrial_project/Mitochondrial_ID_depth_list_horseref.txt | awk '{print $1}' | while read name; do
  sed -i "s/>NC_001640.1/>${name}_>NC_001640.1/g" /shared5/Alex/Mitochondrial_project/FASTAs_horseref/${name}_mtdna_horseref.fa
done
```

#### 8.4 Concatenating FASTAs
Then I compilee all the FASTAs in one file and added the reference mitogenome as well:
```bash
#concatenating FASTAs
cat /shared5/Alex/Mitochondrial_project/Mitochondrial_ID_depth_list_horseref.txt | awk '{print $1}' | while read name; do
  cat /shared5/Alex/Mitochondrial_project/FASTAs_horseref/${name}_mtdna_horseref.fa
done >> Exmoor_mtdna_horseref.fa


## adding reference mitogenome
cat /shared5/Alex/Mitochondrial_project/horse_mtdna_reference/horse_mtdna_reference.fasta Exmoor_mtdna_horseref.fa > tmp.fa && mv tmp.fa Exmoor_mtdna_horseref.fa
```

So my final file "Exmoor_mtdna_horseref.fa" is ready for analysis.
