# Code for first Y-chromosomal analysis of Exmoor ponies

---
* Title: Code for Y-chromosomal analyses
* Author: Alexandra Brumwell
---

Code notebook for honours project on Exmoor pony populations genetics

--------------------------------------------------------------------------------------------
## Introduction
I will be studying the Y-chromosome in 15 male Exmoor ponies sequences.

Prior to starting any data processing or analyses I create a spreadsheet containing the metadata of the sequences like genome ID, date of sequencing, breed, file name, sex, etc. I save all this information in a [spreadsheet](ADD LINK TO SPREADSHEET IN GOOGLEDRIVE) in Excel or Googledrive.

##### Example of metadata table:

| Individual | File_name | Breed | Pony number | Pony name | Sequencing_date | Sex | Notes | Original_FASTQ_file_path |
|------------|-----------|-------|-------------|-----------|-----------------|-----|-------|--------------------------|
| ERR111 | ERR111_EYB1000 | Exmoor | 111/111 | Cheddar | Mar2020 | F | From unknown maternal line | /Path/to/FASTQ/file |



## Downloading/obtaining a reference genome
We need to have a reference genome to map the FASTQs on to. I use the Y-chromosome assembly MH341179.1 (https://www.ncbi.nlm.nih.gov/nuccore/MH341179.1/) from the paper "Horse Y chromosome assembly displays unique evolutionary and putative stallion fertility genes" (Janecka et al., 2018). I downladed it into my directory using `wget -P /chosen_directory/ <copied_link>`.


## Creation of sample list
First I check the pedigree of the indivudals that each sequence is from to see which are male. Then I make a list of all the male Exmoors and save this in a text file called *List_ychromosome.txt*. I add the file path for read 1 of each FASTQ as a correspoding second column

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
  samtools view -q 30 -u /shared5/Alex/Y-chromosome_project/BAMs/${name}_Y.sam | \
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
Then, I create a new list (List_Y-chromosome_merged.txt) with the updated IDs of the sequences. E.g. The merged file IDs have now become "s479021_merged" instead of "s2479021_EDSW210003783-1a_H3WNKDSX2_L4"


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

I got low average coverage from the final BAMs so I made changes to the pipeline and used a different Y-chromosome pipeline, which are featured in the markdown file "Y-chromosome_code_newreference".
