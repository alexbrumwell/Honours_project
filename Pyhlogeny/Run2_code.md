# Code for Run2 of horse

---
* Title: Code for Run2 of horse
* Author: Alexandra Brumwell
---

Code notebook for honours project on Exmoor pony populations genetics

--------------------------------------------------------------------------------------------
## Introduction
Run2 is the second part of the horse phylogeny sideproject. Run2 utilized the remaining horse and Przewalski genomes not utilized in  Run1.

## FASTQ to BAM pipeline

#### Remaining horse sequences (31/10/22)
There are a remaining 31 horse sequences that were not used in Run1 and will be used for Run2. These remaining sequences also need to be processed.

**0. Obtaining FASTQs**
1.	I review the metadata spreadsheet from the student and select the sequences that I haven't used in Run1 and compile these sequences into the spreadsheet for this project called Run2. I now have a total of 31 horse sequnces from various horse breeds.

2.	I create a text file containing all the IDs of the sequences I will be using and I save this textfile (remaining_horse_ID_list.txt) in my directory.

3.	Since the FASTQ files of these sequences are located in different directories, I use a loop to locate where each FASTQ file is and save the file paths in a text file called "horse_fastq_paths":
```bash
cat /shared5/Alex/Horses/Horses_remaining/remaining_horse_ID_list.txt | while read name; do
  ls /shared5/studentprojects/Qikai/Horse_database/*/${name}*val* ;
done >> remaining_horse_fastq_paths.txt
\
#If there was file that wasn't located then I find it by using and then add it to "remaining_horse_fastq_paths.txt":
find shared5/studentprojects/Qikai/Horse_database/ -name “ERR1527947*”
```

4. I created a directory that would contain the soft links to all the FASTQs, so that I could visualize all the horse sequences in one place:
```bash
mkdir /shared5/Alex/Horses/Horses_remaining/FASTQs/ #creates directory

cat /shared5/Alex/Horses/Horses_remaining/remaining_horse_fastq_paths.txt | while read line; do
  ln -s ${line} /shared5/Alex/Horses/Horses_remaining/FASTQs/ ;
done
```

**1. Trimming FASTQs**
The horse FASTQs are already trimmed so I skipped this step.

**2. Mapping**
I used a loop to go over the horse sequence IDs so I can map them in parallel
```bash
cat /shared5/Alex/Horses/Horses_remaining/remaining_horse_ID_list.txt | while read name; do
bwa mem /shared5/Alex/Donkey_ref/GCF_016077325.2_ASM1607732v2_genomic.fna.gz /shared5/Alex/Horses/Horses_remaining/FASTQs/${name}_1_val_1.fq.gz /shared5/Alex/Horses/Horses_remaining/FASTQs/${name}_2_val_2.fq.gz > /shared5/Alex/Horses/Horses_remaining/BAMs/${name}.sam & # the "&" allows the loops to be run in parallel
done
```

* **Sidenote to self:** To check how large the SAM files created by the previous student were, I used:
```bash
cat /shared5/Alex/Horses/Horses_remaining/remaining_horse_ID_list.txt | cut -f 3 | cut -f 1 -d "_") | while read name;
do
file=$(echo "/shared5/studentprojects/Qikai/Horse_database/$name.sam");
ll -h $file;
done
```


**3. Coverting to BAM and sorting by read names**
```bash
cat /shared5/Alex/Horses/Horses_remaining/remaining_horse_ID_list.txt | while read name; do
samtools view -q 30 -u /shared5/Alex/Horses/Horses_remaining/BAMs/${name}.sam | samtools sort -n -o /shared5/Alex/Horses/Horses_remaining/BAMs/${name}_sorted.bam &
done
```

**4. Mark duplicates**
```bash
cat /shared5/Alex/Horses/Horses_remaining/remaining_horse_ID_list.txt | while read name; do
samtools fixmate -m /shared5/Alex/Horses/Horses_remaining/BAMs/${name}_sorted.bam - | samtools sort - | samtools markdup - /shared5/Alex/Horses/Horses_remaining/BAMs/${name}_sorted_markdup.bam &
done
```

**5. Indexing**
```bash
cat /shared5/Alex/Horses/Horses_remaining/remaining_horse_ID_list.txt | while read name; do
  /shared5/Alex/Horses/Horses_remaining/BAMs/${name}_sorted_markdup.bam &
  done
```


**6. Calculating BAM depth**
I used `qualimap` and this creates a folder that generates a file called "genome_results.txt" for each sequence with the coverage data:
```bash
cat /shared5/Alex/Horses/Horses_remaining/remaining_horse_ID_list.txt | while read name; do
  qualimap bamqc -bam /shared5/Alex/Horses/BAMs/${name}_markdup.bam --java-mem-size=200G &
done
```

**7. Consensus FASTA file**
***I HAVE NOT DONE THIS STEP YET***
Prior to obtain the consensus FASTA file, I need to get the mean and maximum depth of each BAM file:
```bash
cat /shared5/Alex/Horses/Horses_remaining/remaining_horse_ID_list.txt | awk '{print $1}' | while read name; do
  coverage=$(cat /shared5/Alex/Horses/Horses_remaining/BAMs/${name}_markdup_stats/genome_results.txt | grep "mean coverageData"  | awk '{print $4}' | sed 's/X//')
paste <(echo ${name}) <(echo ${coverage}) <(echo "${coverage} * 2" | bc);
done >> remaining_horse_ID_depth_list.txt

#horse_ID_depth_list.txt content:
ERR1527947      26.6281 53.2562
ERR1305963      19.9776 39.9552
ERR1512897      18.5663 37.1326
```

Now I can run the consensus FASTA command:
```bash
cat /shared5/Alex/Horses/Horses_remaining/remaining_horse_ID_depth_list.txt | while read line ; do
  #defining variables
  name=$(echo ${line} | awk '{print $1}');
  coverage=$(echo ${line} | awk '{print $2}');
  maxcoverage=$(echo ${line} | awk '{print $3}');

  angsd -i /shared5/Alex/Horses/Horses_remaining/BAMs/${name}_markdup.bam -minMapQ 30 -minQ 20 -setMinDepth ${coverage} -setMaxDepth ${maxcoverage} -remove_bads 1 -doFasta 2 -doCounts 1 -ref /shared5/Alex/Donkey_ref/GCF_016077325.2_ASM1607732v2_genomic.fna -out /shared5/Alex/Horses/Horses_remaining/FASTAs/${name}_DonkeyRef &
  done
```



#### Przewalski FASTQs
These were processed when I was doing the FASTQtoBAM pipeline for Run1.



## Variant calling using Strelka
I used `strelka` for variant calling of the genomes because it is quicker than `VCFtools` and is optimized for Illumina data.

1. I convert the donkey reference FASTA file into a BED file using a script called *genome2bed.py*:
```bash
python2.7 genome2bed.py GCF_016077325.2_ASM1607732v2_genomic.fna GCF_016077325.2_ASM1607732v2_genomic.bed
```
I also created individual bed files for each scaffold and then indexed the BED files :
```bash
# Create individual bed files for scaffolds
cat GCF_016077325.2_ASM1607732v2_genomic.bed | grep "NC_0" | while read line; do
  name=$(echo ${line} | awk '{print $1}');
  echo $line | awk '{print $1"\t"$13"\t"$14}' > ${name}.txt;
done

#indexing bed files
for file in {177..208}; do
  bgzip NC_052${file}/NC_052${file}.1.bed > NC_052${file}/NC_052${file}.1.bed.gz;
  tabix -p bed  NC_052${file}/NC_052${file}.1.bed.gz;
done
```


2. I created a directory for the second run of Strelka and folders for each scaffold.
```bash
mkdir Run2

cd Run2

for file in {177..208}; do
  mkdir NC_052${file}
done
```



2.	Then I created a Strelka script that uses the remaining horse and Przewalski genomes and loops over the different scaffold folders :

```bash
#old and incorrect configure_strelka.sh content:
for dir in /shared5/Alex/Run2/*/ ; do
    cd $dir
    mkdir -p Out # the "-p" makes the command only create a directory if it doesn't already exit
    OUT_DIR=$(echo `pwd`/Out)
    CB=$(echo *bed.gz)
    CURR_BED=$(echo `pwd`/$CB)

    configureStrelkaGermlineWorkflow.py \
    --bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1305962_sorted_markdup.bam \
    --bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1305964_sorted_markdup.bam \
    --bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1306526_sorted_markdup.bam \
    --bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1527948_sorted_markdup.bam \
    --bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1527950_sorted_markdup.bam \
    --bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1527952_sorted_markdup.bam \
    --bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1527958_sorted_markdup.bam \
    --bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1527966_sorted_markdup.bam \
    --bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1527969_sorted_markdup.bam \
    --bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1527970_sorted_markdup.bam \
    --bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1527972_sorted_markdup.bam \
    --bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1545180_sorted_markdup.bam \
    --bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1545181_sorted_markdup.bam \
    --bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1545184_sorted_markdup.bam \
    --bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1545185_sorted_markdup.bam \
    --bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1545187_sorted_markdup.bam \
    --bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1545188_sorted_markdup.bam \
    --bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1545190_sorted_markdup.bam \
    --bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR2179545_sorted_markdup.bam \
    --bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR2179547_sorted_markdup.bam \
    --bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR2179548_sorted_markdup.bam \
    --bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR2179550_sorted_markdup.bam \
    --bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR2179554_sorted_markdup.bam \
    --bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR2179555_sorted_markdup.bam \
    --bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR2179556_sorted_markdup.bam \
    --bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR2731056_sorted_markdup.bam \
    --bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR2731060_sorted_markdup.bam \
    --bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR863167_sorted_markdup.bam \
    --bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR868003_sorted_markdup.bam \
    --bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR978597_sorted_markdup.bam \
    --bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR978599_sorted_markdup.bam \
    --bam /shared5/Alex/Przewalski/BAMs/SRR12719745_markdup.bam \
    --bam /shared5/Alex/Przewalski/BAMs/SRR12719757_markdup.bam \
    --bam /shared5/Alex/Przewalski/BAMs/SRR12719758_markdup.bam \
    --referenceFasta /shared5/Alex/Donkey_ref/GCF_016077325.2_ASM1607732v2_genomic.fna \
    --runDir .
    cd ../
    done

#correct version:
for dir in */
do
cd $dir
mkdir Out
OUT_DIR=$(echo `pwd`/Out)
CB=$(echo *bed.gz)
CURR_BED=$(echo `pwd`/$CB)

configureStrelkaGermlineWorkflow.py \
--bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1305962_sorted_markdup.bam \
--bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1305964_sorted_markdup.bam \
--bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1306526_sorted_markdup.bam \
--bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1527948_sorted_markdup.bam \
--bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1527950_sorted_markdup.bam \
--bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1527952_sorted_markdup.bam \
--bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1527958_sorted_markdup.bam \
--bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1527966_sorted_markdup.bam \
--bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1527969_sorted_markdup.bam \
--bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1527970_sorted_markdup.bam \
--bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1527972_sorted_markdup.bam \
--bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1545180_sorted_markdup.bam \
--bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1545181_sorted_markdup.bam \
--bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1545184_sorted_markdup.bam \
--bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1545185_sorted_markdup.bam \
--bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1545187_sorted_markdup.bam \
--bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1545188_sorted_markdup.bam \
--bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR1545190_sorted_markdup.bam \
--bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR2179545_sorted_markdup.bam \
--bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR2179547_sorted_markdup.bam \
--bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR2179548_sorted_markdup.bam \
--bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR2179550_sorted_markdup.bam \
--bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR2179554_sorted_markdup.bam \
--bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR2179555_sorted_markdup.bam \
--bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR2179556_sorted_markdup.bam \
--bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR2731056_sorted_markdup.bam \
--bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR2731060_sorted_markdup.bam \
--bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR863167_sorted_markdup.bam \
--bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR868003_sorted_markdup.bam \
--bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR978597_sorted_markdup.bam \
--bam /shared5/Alex/Horses/Horses_remaining/BAMs/ERR978599_sorted_markdup.bam \
--bam /shared5/Alex/Przewalski/BAMs/SRR12719745_markdup.bam \
--bam /shared5/Alex/Przewalski/BAMs/SRR12719757_markdup.bam \
--bam /shared5/Alex/Przewalski/BAMs/SRR12719758_markdup.bam \
--referenceFasta /shared5/Alex/Donkey_ref/GCF_016077325.2_ASM1607732v2_genomic.fna \
--callRegions $CURR_BED \
--runDir $OUT_DIR
cd ..
done
```


3.	The I ran the script:
```bash
source activate strelka
bash configure_strelka.sh
```

4.	This created python scripts for each scaffold called “runWorkflow.py”, which I then executed:
```bash
for dir in */ ; do
./${dir}runWorkflow.py -m local -j 2 &
done

#I rerun this but covering fewer chromosomes and more CPUs
for dir in {227..238}; do
  cd NC_12${dir}
  ./runWorkflow.py -m local -j 8 &
done

```
* `-m`: *Run mode*
* `-j`: *Number of jobs*. Default is the estimate total cores on this node for local mode.



## VCF filtering (31/10/22)
I HAVE NOT REACHED THIS STAGE YET
