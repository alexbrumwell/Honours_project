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
I use Bioedit to view the sequences and align them.
Comments:
- Sequence "S49052_EDSW200015551" contains no info so I remove it from the file

I get almost the same results as when mapped to the donkey reference. To do: compare the two haplotype networks.
