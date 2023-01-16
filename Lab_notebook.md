# Honours Project Lab Notebook

---
* Title: Honours_project_labnotebook
* Author: Alexandra Brumwell
---

Code notebook for honours project on Exmoor pony populations genetics

----------------------------------------------------------------------------------------------
## 25/11/22
Submitted mapping of Exmoor sequences (March 2021 & Sep 2020) @poland (14:05h):
```bash
for file in $(cat /shared5/Alex/Exmoor_sequencing_data/old_exmoor_list.txt| awk '{print $1}'); do
  location=$(echo ${file} | awk -F "/" '{$NF=""}1' OFS="/") ;
  name=$(echo ${file} | awk -F "/" '{print $NF}' | sed 's/_1_val_1.fq.gz//');
  bwa mem /shared5/Alex/Donkey_ref/Donkey_mitogenome.fasta \
  ${location}${name}_1_val_1.fq.gz \
  ${location}${name}_2_val_2.fq.gz \
  2> /shared5/Alex/Mitochondrial_project/BAMs/log_files/${name}.bwamem.log \
  > /shared5/Alex/Mitochondrial_project/BAMs/${name}_mtdna.sam & \
done
```

## 28/11/22
Resubmitted sorting of horse mitogenomes on @young (submitted 14:22h):
```bash
for file in $(cat /shared5/Alex/Mitochondrial_project/Mitochondrial_list.txt | grep "ERR" | awk '{print $1}'); do
  name=$(echo ${file} | awk -F "/" '{print $NF}' | sed 's/_1_val_1.fq.gz//');
  samtools view -q 30 -u /shared5/Alex/Mitochondrial_project/BAMs/${name}_mtdna.sam | \
  samtools sort -n -o /shared5/Alex/Mitochondrial_project/BAMs/${name}_mtdna_sorted.bam 2> /shared5/Alex/Mitochondrial_project/BAMs/log_files/${name}_mtdna_sorted.log & \
done
```


## 29/11/22
Today I was mainly reading about Markdown files.

To check if the Exmoor (Mar21 and Sep2020) sequences have finished mapping by checking the logfiles:
```bash
cat /shared5/Alex/Mitochondrial_project/Mitochondrial_list.txt | grep 'Mar_2021\|Sep_2020' | while read line; do
name=$(echo ${line} | awk -F "/" '{print $NF}' | sed 's/_1_val_1.fq.gz//');
tail -n 1 /shared5/Alex/Mitochondrial_project/BAMs/log_files/$name*bwamem.log
done

# Output:
[main] Real time: 365278.453 sec; CPU: 10060.185 sec
[main] Real time: 343939.272 sec; CPU: 8433.604 sec
```


## 30/11/22
Submitted sorting of Exmoor (Mar21 and Sep2020) mitogenomes on @howe (submitted 12:22h):
```bash
cat /shared5/Alex/Mitochondrial_project/Mitochondrial_list.txt | grep 'Mar_2021\|Sep_2020' | awk '{print $2}' |  while read file ; do
  name=$(echo ${file} | awk -F "/" '{print $NF}' | sed 's/_1_val_1.fq.gz//');
  samtools view -q 30 -u /shared5/Alex/Mitochondrial_project/BAMs/${name}_mtdna.sam | \
  samtools sort -n -o /shared5/Alex/Mitochondrial_project/BAMs/${name}_mtdna_sorted.bam 1> /shared5/Alex/Mitochondrial_project/BAMs/log_files/${name}_mtdna_sorted.log & \
done
```

## 01/12/22
To check size of SAMs and BAMs created:
```bash
cat /shared5/Alex/Mitochondrial_project/Mitochondrial_list.txt | awk '{print $1}' | while read name; do paste <(ls -lh ${name}*sam | awk '{print $5}') <(ls -lh ${name}*bam | awk '{print $5}');
done

#Output:
69G     9.6M
70G     8.0M
```

## 02/12/22
To run ValidateSamFile on mitogenome sam files:
```bash
cat /shared5/Alex/Mitochondrial_project/Mitochondrial_list.txt | awk '{print $1}' | while read name; do
    java -jar /shared5/Alex/picard-2.27.5/picard.jar ValidateSamFile \
    -I /shared5/Alex/Mitochondrial_project/BAMs/${name}_mtdna.sam \
    -O /shared5/Alex/Mitochondrial_project/BAMs/ValidateSamFile_stats/${name}_mtdna_sam \
    -R /shared5/Alex/Donkey_ref/Donkey_mitogenome.fasta \
    --MODE SUMMARY
done
```

To check if sorted BAMs have finished creating:
```bash
cat /shared5/Alex/Mitochondrial_project/Mitochondrial_list.txt | awk '{print $1}' | while read name; do
  paste <(echo ${name}) <(ls -lh /shared5/Alex/Mitochondrial_project/BAMs/${name}_mtdna_sorted.bam | awk '{print $5}')
done
```

## 03/12/22
All the mitogenomes have finished mapping and sorting so now before marking duplicates, I have to merge files that are from the same indiviudals. The only seuqences that are separated into mulitple files are the Exmoor sequences from September 2022 and March 2021.
To check which individuals are separated into multiple files:
```bash
# To check if the Exmoor September 2022 individuals are repeated:
cat /shared5/Alex/Mitochondrial_project/Mitochondrial_list.txt | grep "Sep_2022" | awk '{print $1}' | cut -d "_" -f 2 | uniq -c | grep -wv "1" | awk '{print "E_$2"}'

# Output:
E_23279
E_23416
E_23434
E_320005
E_32023
E_479023
E_512001
E_78170
E_900588
E_900717


# To check if the Exmoor September 2020 and Mar 2021 individuals are repeated:
cat /shared5/Alex/Mitochondrial_project/Mitochondrial_list.txt | grep "Mar_2021" | awk '{print $1}' | cut -d "_" -f 1 | uniq -c |grep -wv "1" | awk '{print $2}'

# Output:
s479021
s49127

# Then I compile the output into a text file (tmp.txt) and use it to merge the files:
cat tmp.txt | while read name; do
  samtools merge ${name}_merged_mtdna_sorted.bam $(ls $name*sorted*);
done

#Then I create a new list (List_mitochondrial_merged.txt) with the updated IDs of the sequences. E.g. The merged file IDs have now become "s479021_merged" instead of "s2479021_EDSW210003783-1a_H3WNKDSX2_L4"
```

Submitted markdup of mitogenomes @ young (14:22h):
```bash
cat /shared5/Alex/Mitochondrial_project/List_mitochondrial_merged.txt | while read name ; do
  samtools fixmate -m /shared5/Alex/Mitochondrial_project/BAMs/${name}_mtdna_sorted.bam - | samtools sort - | samtools markdup - /shared5/Alex/Mitochondrial_project/BAMs/${name}_mtdna_sorted_markdup.bam 1> /shared5/Alex/Mitochondrial_project/BAMs/log_files/${name}_markdup.log &
 done
```

Submitted indexing of mitogenomes @young (14:27h):
```bash
cat /shared5/Alex/Mitochondrial_project/List_mitochondrial_merged.txt | awk '{print $1}' | while read name ; do
  samtools index /shared5/Alex/Mitochondrial_project/BAMs/${name}_mtdna_sorted_markdup.bam &
done
```

Submitted qualimap of mitogenomes @young (14:35h)
```bash
cat /shared5/Alex/Mitochondrial_project/List_mitochondrial_merged.txt | awk '{print $1}' | while read name ; do
  qualimap bamqc -bam /shared5/Alex/Mitochondrial_project/BAMs/${name}_mtdna_sorted_markdup.bam --java-mem-size=200G &
done
```

Submitted ValidateSamFile on markdup mitogenomes @howe (14:41h):
```bash
cat /shared5/Alex/Mitochondrial_project/List_mitochondrial_merged.txt | awk '{print $1}' | while read name ; do
  mkdir -p /shared5/Alex/Mitochondrial_project/BAMs/ValidateSamFile/
  java -jar /shared5/Alex/picard-2.27.5/picard.jar ValidateSamFile -I /shared5/Alex/Mitochondrial_project/BAMs/${name}_mtdna_sorted_markdup.bam -O /shared5/Alex/Mitochondrial_project/BAMs/ValidateSamFile/${name}_mtdna_sorted_markdup_errors.txt -R /shared5/Alex/Donkey_ref/Donkey_mitogenome.fasta --MODE SUMMARY
done
```
Creating the y-chromosome list. Finding the FASTQs from the genome ID:
```bash
cat /shared5/Alex/Y-chromosome_project/List_ychromosome.txt | head -n 22 | while read name; do
  find /shared5/Alex/Exmoor_sequencing_data/ -name "${name}*val_1.fq.gz";
done
```

Submitted mapping of Y-chromsomes (18:11h):
```bash
cat /shared5/Alex/Y-chromosome_project/List_ychromosome.txt | awk '{print $2}' | while read file ; do
  location=$(echo ${file} | awk -F "/" '{$NF=""}1' OFS="/") ;
  name=$(echo ${file} | awk -F "/" '{print $NF}' | sed 's/_1_val_1.fq.gz//');
  bwa mem /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/horse_Y_chromosome.fasta ${location}${name}_1_val_1.fq.gz ${location}${name}_2_val_2.fq.gz 2> /shared5/Alex/Y-chromosome_project/BAMs/log_files/${name}.bwamem.log > /shared5/Alex/Y-chromosome_project/BAMs/${name}_Y.sam &
done
```


## 05/12/22
To check if qualimap has finished creating files:
```bash
cat /shared5/Alex/Mitochondrial_project/List_mitochondrial_merged.txt | awk '{print $1}' | while read name; do
  paste <(echo ${name}) <(ls -lh /shared5/Alex/Mitochondrial_project/BAMs/${name}*stats/genome_results.txt)
done

# To check the log files:
cat /shared5/Alex/Mitochondrial_project/List_mitochondrial_merged.txt | awk '{print $1}' | while read name; do
  ls -lh /shared5/Alex/Mitochondrial_project/BAMs/log_files/${name}*qualimap.log
done
```

Apparently qualimap had only ran for 107 mitogenomes out of the total 117 mitogenomes so I killed the processes using `pkill java` and then reran qualimap:
```bash
cat /shared5/Alex/Mitochondrial_project/List_mitochondrial_merged.txt | awk '{print $1}' | while read name ; do
  qualimap bamqc -bam /shared5/Alex/Mitochondrial_project/BAMs/${name}_mtdna_sorted_markdup.bam --java-mem-size=200G 2>/shared5/Alex/Mitochondrial_project/BAMs/log_files/${name}_qualimap.log &
done
```
I seem to be getting errors for 10 mitogenomes:
>ERR1545188\
ERR2179547\
S49052_EDSW200015551-1a_HJFWVDSXY_L4\
S276023_EDSW200015547-1a2a_HJFWVDSXY_L3L4\
s78150_EDSW210003780-1a_H3WNKDSX2_L4\
ERR1527948\
s49149_EDSW210003771-1a_H3FTWDSX2_L2\
s21131_EDSW210003778-1a_H3WNKDSX2_L1\
E_900600_EKDN220034356-1A_H5GL5DSX5_L2\
E_900717_merged\

The log files of qualimap of these genomes state:
````java
OpenJDK 64-Bit Server VM warning: Ignoring option MaxPermSize; support was removed in 8.0
Failed to run bamqc

net.sf.samtools.SAMFormatException: SAM validation error: ERROR: Record 1, Read name A00627:132:HJV7LDSXY:3:1103:31602:23187, Second of pair flag should not be set for unpaired read.
        at net.sf.samtools.SAMUtils.processValidationErrors(SAMUtils.java:448)
        at net.sf.samtools.BAMFileReader$BAMFileIterator.advance(BAMFileReader.java:506)
        at net.sf.samtools.BAMFileReader$BAMFileIterator.<init>(BAMFileReader.java:465)
        at net.sf.samtools.BAMFileReader$BAMFileIterator.<init>(BAMFileReader.java:453)
        at net.sf.samtools.BAMFileReader.getIterator(BAMFileReader.java:255)
        at net.sf.samtools.SAMFileReader.iterator(SAMFileReader.java:315)
        at org.bioinfo.ngs.qc.qualimap.process.BamStatsAnalysis.run(BamStatsAnalysis.java:275)
        at org.bioinfo.ngs.qc.qualimap.main.BamQcTool.execute(BamQcTool.java:197)
        at org.bioinfo.ngs.qc.qualimap.main.NgsSmartTool.run(NgsSmartTool.java:187)
        at org.bioinfo.ngs.qc.qualimap.main.NgsSmartMain.main(NgsSmartMain.java:103)
````

I try to run qualimap on these samples individually:
```bash
 qualimap bamqc -bam /shared5/Alex/Mitochondrial_project/BAMs/ERR1545188_mtdna_sorted_markdup.bam --java-mem-size=200G 2> /shared5/Alex/Mitochondrial_project/BAMs/log_files/ERR1545188_qualimap.log
 ```
 But it still gives me the same error. So I checked the ValidateSamFile of this sample and saw these errors:
```bash
 ## HISTOGRAM    java.lang.String
Error Type      Count
ERROR:INVALID_FLAG_FIRST_OF_PAIR        2938
ERROR:INVALID_FLAG_MATE_UNMAPPED        5008
ERROR:INVALID_FLAG_SECOND_OF_PAIR       2681
ERROR:MATES_ARE_SAME_END        212
ERROR:MATE_CIGAR_STRING_INVALID_PRESENCE        170
ERROR:MISMATCH_FLAG_MATE_UNMAPPED       170
ERROR:MISSING_READ_GROUP        1
WARNING:RECORD_MISSING_READ_GROUP       120465
```
However these errors are also found in files that ran successfully in qualimap like *E_78170_merged_*:
```bash
## HISTOGRAM    java.lang.String
Error Type      Count
ERROR:INVALID_FLAG_FIRST_OF_PAIR        2305
ERROR:INVALID_FLAG_MATE_UNMAPPED        3612
ERROR:INVALID_FLAG_SECOND_OF_PAIR       2246
ERROR:MATES_ARE_SAME_END        455
ERROR:MATE_CIGAR_STRING_INVALID_PRESENCE        296
ERROR:MISMATCH_FLAG_MATE_NEG_STRAND     20
ERROR:MISMATCH_FLAG_MATE_UNMAPPED       296
ERROR:MISMATCH_MATE_ALIGNMENT_START     20
ERROR:MISMATCH_MATE_CIGAR_STRING        20
ERROR:MISSING_READ_GROUP        1
WARNING:RECORD_MISSING_READ_GROUP       222605
```

 I ran ValidateSamFile on this sample but the WGS rather than the mitogenome to see if this sample gave issues in other cases. But the WGS sample doesn't seem to have any unusual errors.



## 06/12/22

 Quick script to calculate mean coverage (found in an online forum):
```bash
# Script to calculate the average coverage of a genome sample (source: found online in a forum)

# compute the total length of the sample
bam="ERR2179547*markdup*bam"
tot=$(samtools view -H $bam | grep '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}')
echo $tot
#compute the coverage at each point and then calculate average coverage
sum=$(samtools depth $bam | awk '{sum+=$3} END {print sum}')
echo $sum
echo
avg=$(echo "${sum}/${tot}" | bc -l)
echo "The average coverage of the sample ${bam} is ${avg} x."
```


Resubmitted qualimap for mitogenomes (14:07) and this time it ran successfully on all samples even though I used the same code as last time:
```bash
cat /shared5/Alex/Mitochondrial_project/List_mitochondrial_merged.txt | awk '{print $1}' | while read name ; do
  qualimap bamqc -bam /shared5/Alex/Mitochondrial_project/BAMs/${name}_mtdna_sorted_markdup.bam --java-mem-size=200G 2>/shared5/Alex/Mitochondrial_project/BAMs/log_files/${name}_qualimap.log &
done
```

Submitted consensus FASTA for mitogenomes (14:33h):
```bash
#To calculate mean and max coverage:
cat /shared5/Alex/Mitochondrial_project/List_mitochondrial_merged.txt | awk '{print $1}' | while read name; do
  coverage=$(cat /shared5/Alex/Mitochondrial_project/BAMs/${name}*stats/genome_results.txt | grep "mean coverageData"  | awk '{print $4}' | sed 's/X//' | sed 's/,//' | cut -d "." -f 1);
  paste <(echo ${name}) <(echo ${coverage}) <(echo "${coverage} * 2" | bc);
done

# Create consensus FASTA:
cat /shared5/Alex/Mitochondrial_project/Mitochondrial_ID_depth_list.txt | while read line; do
  name=$(echo ${line} | awk '{print $1}');
  number=$(echo ${line} | awk '{print $(NF-1)}');
  maxnumber=$(echo ${line} | awk '{print $(NF)}');
  angsd -i /shared5/Alex/Mitochondrial_project/BAMs/${name}_mtdna_sorted_markdup.bam -minMapQ 30 -minQ 20 -setMinDepth ${number} -setMaxDepth ${maxnumber} -remove_bads 1 -doFasta 2 -doCounts 1 -ref /shared5/Alex/Donkey_ref/Donkey_mitogenome.fasta -out /shared5/Alex/Mitochondrial_project/FASTAs/${name}_mtdna_DonkeyRef &
done
```


## 07/12/22
ERROR - 1 mitogenome FASTA is empty --> tried to redo consensus fasta for one of the files. But it doesn't change anything:
```bash
angsd -i /shared5/Alex/Mitochondrial_project/BAMs/ERR1527968_mtdna_sorted_markdup.bam -minMapQ 30 -minQ 20 -setMinDepth 375 -setMaxDepth 750 -remove_bads 1 -doFasta 2 -doCounts 1 -ref /shared5/Alex/Donkey_ref/Donkey_mitogenome.fasta -out /shared5/Alex/Mitochondrial_project/FASTAs/ERR1527968_mtdna_DonkeyRef
```

Many of the Exmoor FASTAS are empty so I rerun ANGSD but with the minimum depth set at 50x instead and for all individuals:
```bash
cat /shared5/Alex/Mitochondrial_project/Mitochondrial_ID_depth_list.txt | while read line; do
  name=$(echo ${line} | awk '{print $1}');
  maxnumber=$(echo ${line} | awk '{print $(NF)}');
  angsd -i /shared5/Alex/Mitochondrial_project/BAMs/${name}_mtdna_sorted_markdup.bam -minMapQ 30 -minQ 20 -setMinDepth 50 -setMaxDepth ${maxnumber} -remove_bads 1 -doFasta 2 -doCounts 1 -ref /shared5/Alex/Donkey_ref/Donkey_mitogenome.fasta -out /shared5/Alex/Mitochondrial_project/FASTAs/${name}_mtdna_DonkeyRef &
done

#find which samples have mean coverage <50x
awk '$2<50' /shared5/Alex/Mitochondrial_project/Mitochondrial_ID_depth_list.txt
#output:
ERR1545180      22      44
ERR1545184      29      58
ERR1527970      40      80
ERR1545187      35      70

#I rerun angsd for these samples:
cat /shared5/Alex/Mitochondrial_project/Mitochondrial_ID_depth_list.txt | grep 'ERR1545180\|ERR1545184\|ERR1527970\|ERR1545187' | while read line; do
  name=$(echo ${line} | awk '{print $1}');
  number=$(echo ${line} | awk '{print $(NF-1)}');
  maxnumber=$(echo ${line} | awk '{print $(NF)}');
  angsd -i /shared5/Alex/Mitochondrial_project/BAMs/${name}_mtdna_sorted_markdup.bam -minMapQ 30 -minQ 20 -setMinDepth ${number} -setMaxDepth ${maxnumber} -remove_bads 1 -doFasta 2 -doCounts 1 -ref /shared5/Alex/Donkey_ref/Donkey_mitogenome.fasta -out /shared5/Alex/Mitochondrial_project/FASTAs/${name}_mtdna_DonkeyRef &
done
```

I have to add the sequence ID to the header of each FASTA mitogenome file:
```bash
cat /shared5/Alex/Mitochondrial_project/List_mitochondrial_merged.txt | while read line; do
  name=$(echo $line | awk '{print $1}')
  breed=$(echo $line | awk '{print $2}')
  sed -i "s/>NC_001788.1/>${name}_${breed} NC_001788.1/g" /shared5/Alex/Mitochondrial_project/FASTAs/${name}_mtdna_DonkeyRef.fa
done

#shorten header ID for Sep 2022 Exmoor sequences
cat /shared5/Alex/Mitochondrial_project/List_mitochondrial_merged.txt | grep "EKDN" | awk '{print $1}' | while read name; do
  ID=$(grep ">" /shared5/Alex/Mitochondrial_project/FASTAs/${name}_mtdna_DonkeyRef.fa | cut -d "_" -f 1,2)
  sed -i "s/>${name}/${ID}/g" /shared5/Alex/Mitochondrial_project/FASTAs/${name}_mtdna_DonkeyRef.fa
done

#shorten header ID for other Exmoor sequences
cat /shared5/Alex/Mitochondrial_project/List_mitochondrial_merged.txt | grep "EDSW" | awk '{print $1}' | while read name; do
  ID=$(grep ">" /shared5/Alex/Mitochondrial_project/FASTAs/${name}_mtdna_DonkeyRef.fa | cut -d "_" -f 1)
  sed -i "s/>${name}/${ID}/g" /shared5/Alex/Mitochondrial_project/FASTAs/${name}_mtdna_DonkeyRef.fa
done
```
Concatinating FASTAs :
```bash
#concatinating FASTAs:
cat /shared5/Alex/Mitochondrial_project/List_mitochondrial_merged.txt | awk '{print $1}' | while read name; do
  cat /shared5/Alex/Mitochondrial_project/FASTAs/${name}_mtdna_DonkeyRef.fa
  echo " "
done > horse_mtdna_DonkeyRef.fa

#add donkey reference fasta
cat /shared5/Alex/Donkey_ref/Donkey_mitogenome.fasta horse_mtdna_DonkeyRef.fa  >> tmp.fa && mv tmp.fa horse_mtdna_DonkeyRef.fa

#edit donkey reference header
sed -i 's/>NC_001788.1 Equus asinus mitochondrion, complete genome/>DonkeyRef NC_001788.1/g' horse_mtdna_DonkeyRef.fa
```

Submitted to Clustal Omega (https://www.ebi.ac.uk/Tools/services/web_clustalo/toolform.ebi)
#Clustal Omega: If you don't receive any email, please check the status of your job by following this link: toolresult.ebi?jobId=clustalo-E20221207-133418-0149-13690731-p1m

## 08/12/22
Submitted CollectAlignmentSummaryMetrics for mitogenomes (11:24):
```bash
textfile="/shared5/Alex/Mitochondrial_project/List_mitochondrial_merged.txt"
dir="/shared5/Alex/Mitochondrial_project/BAMs"

cat ${textfile} | awk '{print $1}' | while read name ; do
  mkdir -p ${dir}/CollectAlignmentSummaryMetrics/
  java -jar /shared5/Alex/picard-2.27.5/picard.jar CollectAlignmentSummaryMetrics \
  -I ${dir}/${name}_mtdna_sorted_markdup.bam \
  -O ${dir}/CollectAlignmentSummaryMetrics/${name}_mtdna_sorted_markdup_stats.txt \
  -R /shared5/Alex/Donkey_ref/Donkey_mitogenome.fasta \
  --VALIDATION_STRINGENCY SILENT
done
```

## 09/12/22
Selected only the reference and Exmoor FASTAs from the consensus:
```bash
grep -n ">" horse_mtdna_DonkeyRef.fa
#choose lines that only include reference and exmoor sequences:
sed -n 1,19057p horse_mtdna_DonkeyRef.fa > Exmoor_DonkeyRef.fa
```
Check: Job ID: clustalo-I20221209-134632-0861-65896462-p1m  >> ERROR NO INPUT

## 10/12/22
Call with Anubhab, he aligned the FASTAs for me using MUSCLE.
I downloaded PopART but the "minimum spanning network" option keeps crashing.

Submitted sorting of Y-chromosomes @young:
```bash
cat /shared5/Alex/Y-chromosome_project/List_ychromosome.txt | awk '{print $1}' | while read name ; do
  samtools view -q 30 -u /shared5/Alex/Y-chromosome_project/BAMs/${name}_Y.sam | samtools sort -n -o /shared5/Alex/Y-chromosome_project/BAMs/${name}_Y_sorted.bam &
done
```

## 13/12/22
Resubmitted sorting of Y-chromosomes @young (11:15h):
```bash
cat /shared5/Alex/Y-chromosome_project/List_ychromosome.txt | awk '{print $1}' | while read name ; do
  samtools view -q 30 -u /shared5/Alex/Y-chromosome_project/BAMs/${name}_Y.sam | samtools sort -n -o /shared5/Alex/Y-chromosome_project/BAMs/${name}_Y_sorted.bam &
done
```

## 18/12/22
Submitted samtools merge of Y-chromosome files:
```bash
samtools merge <output file> <input files>
```

Then I make a new ID list called "List_Y-chromosome_merged.txt"

Submitted markdup of Y-chromsome files (20:34h):
```bash
cat /shared5/Alex/Y-chromosome_project/List_Y-chromosome_merged.txt | awk '{print $1}' | while read name ; do
  samtools fixmate -m /shared5/Alex/Y-chromosome_project/BAMs/${name}_Y_sorted.bam - | \
  samtools sort - | \
  samtools markdup - /shared5/Alex/Y-chromosome_project/BAMs/${name}_Y_sorted_markdup.bam &
done
```



## 20/12/22
Submitted indexing of Y-chromosome:
```bash
cat /shared5/Alex/Y-chromosome_project/List_Y-chromosome_merged.txt | awk '{print $1}' | while read name ; do
  samtools index /shared5/Alex/Y-chromosome_project/BAMs/${name}_Y_sorted_markdup.bam &
done
```

Submitted Qualimap of Y-chromosomes (14:57h):
```bash
cat /shared5/Alex/Y-chromosome_project/List_Y-chromosome_merged.txt | awk '{print $1}' | while read name ; do
  qualimap bamqc -bam /shared5/Alex/Y-chromosome_project/BAMs/${name}_Y_sorted_markdup.bam --java-mem-size=200G 2>/shared5/Alex/Y-chromosome_project/BAMs/log_files/${name}_qualimap.log &
done
```

Check Qualimap of Y-chromosomes:
```bash
cat /shared5/Alex/Y-chromosome_project/List_Y-chromosome_merged.txt | awk '{print $1}' | while read name; do
  paste <(echo ${name}) <(ls -lh /shared5/Alex/Y-chromosome_project/BAMs/${name}*stats/genome_results.txt)
done

#the following indiivudals have not run correctly:
S276023_EDSW200015547-1a2a_HJFWVDSXY_L3L4_qualimap.log
ERR1735862_qualimap.log
ERR1545189_qualimap.log
ERR2179547_qualimap.log
ERR2731057_qualimap.log

I run qualimap individually on these sequences.
```

## 21/12/22
To check how many unmapped reads in BAM files:
```bash
cat /shared5/Alex/Y-chromosome_project/List_Y-chromosome_merged.txt | awk '{print $1}' | while read name; do
  samtools view -c -f 4 ${name}*markdup*bam;
done
```

Submitted consensus FASTA calling (11:28h):
```bash
#first make list with coverage data
cat /shared5/Alex/Y-chromosome_project/List_Y-chromosome_merged.txt | awk '{print $1}' | while read name; do
  coverage=$(cat /shared5/Alex/Y-chromosome_project/BAMs/${name}*stats/genome_results.txt | grep "mean coverageData"  | awk '{print $4}' | sed 's/X//' | sed 's/,//' | cut -d "." -f 1);
  paste <(echo ${name}) <(echo ${coverage}) <(echo "${coverage} * 2" | bc);
done > /shared5/Alex/Y-chromosome_project/Y-chromosome_ID_depth_list.txt

#then call consensus FASTAS
cat /shared5/Alex/Y-chromosome_project/Y-chromosome_ID_depth_list.txt | while read line; do
  name=$(echo ${line} | awk '{print $1}');
  number=$(echo ${line} | awk '{print $(NF-1)}');
  maxnumber=$(echo ${line} | awk '{print $(NF)}');
 angsd -i /shared5/Alex/Y-chromosome_project/BAMs/${name}_Y_sorted_markdup.bam -minMapQ 30 -minQ 20 -setMinDepth ${number} -setMaxDepth ${maxnumber} -remove_bads 1 -doFasta 2 -doCounts 1 -ref /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/horse_Y_chromosome.fasta -out /shared5/Alex/Y-chromosome_project/FASTAs/${name}_Y &
done

#uncompress the FASTAs
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

#concatented FASTAs
cat /shared5/Alex/Y-chromosome_project/Y-chromosome_ID_depth_list.txt | awk '{print $1}' | while read name; do
  cat /shared5/Alex/Y-chromosome_project/FASTAs/${name}_Y.fa
done >> Exmoor_horse_Y.fa

#added reference sequence
cat /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/horse_Y_chromosome.fasta Exmoor_horse_Y.fa > tmp.fa && mv tmp.fa Exmoor_horse_Y.fa

#copy to local disk
scp studentprojects@young.eng.gla.ac.uk:/shared5/Alex/Y-chromosome_project/FASTAs/Exmoor_horse_Y.fa ./Desktop
```

To check how many lines in each FASTA have clear base pairs:
```bash
cat /shared5/Alex/Y-chromosome_project/Y-chromosome_ID_depth_list.txt | awk '{print $1}' | while read name; do
grep -v "NNNNN" /shared5/Alex/Y-chromosome_project/FASTAs/${name}_Y.fa | wc -l ;
done
```

To get only the reference and Exmoor sequences from the Y-chromosome FASTA, I do:
```bash
grep -n ">" Exmoor_horse_Y.fa #to get the line number of each headers

sed -n 1,2978723p Exmoor_horse_Y.fa > Exmoor_Y.fa #prints only lines containing sequences that I'm interested in
```


## 24/12/22
Submitted mapping of Exmoor sequences to horse reference mitogenome (12:52h):
```bash
cat /shared5/Alex/Mitochondrial_project/Mitochondrial_list.txt | grep "Exmoor" | awk '{print $2}' | while read file ; do
  location=$(echo ${file} | awk -F "/" '{$NF=""}1' OFS="/") ;
  name=$(echo ${file} | awk -F "/" '{print $NF}' | sed 's/_1_val_1.fq.gz//');
  bwa mem /shared5/Alex/Mitochondrial_project/horse_mtdna_reference/horse_mtdna_reference.fasta
  ${location}${name}_1_val_1.fq.gz
  ${location}${name}_2_val_2.fq.gz
  2> /shared5/Alex/Mitochondrial_project/BAMs_horseref/log_files/${name}.bwamem.log
  > /shared5/Alex/Mitochondrial_project/BAMs_horseref/${name}_mtdna_horseref.sam &
done
```


## 27/12/22
Submitted sorting of mitochondrial Exmoor sequences mapped to horse ref:
```bash
cat /shared5/Alex/Mitochondrial_project/Mitochondrial_list.txt | grep "Exmoor" | awk '{print $1}' | while read name ; do
  samtools view -q 30 -u /shared5/Alex/Mitochondrial_project/BAMs_horseref/${name}_mtdna_horseref.sam | samtools sort -n -o /shared5/Alex/Mitochondrial_project/BAMs_horseref/${name}_sorted_mtdna_horseref.bam &
done
```

## 30/12/22
Merged Exmoor sequences horse_ref:
```bash
cat tmp.txt | while read name; do
  samtools merge ${name}_merged_sorted_mtdna_horseref.bam $(ls /shared5/Alex/Mitochondrial_project/BAMs_horseref/$name*sorted*bam);
done
```

Submitted markdup of Exmoor sequences horse_ref (10:30h):
```bash
cat /shared5/Alex/Mitochondrial_project/List_mitochondrial_merged.txt | grep "Exmoor" | awk '{print $1}' | while read name ; do
  samtools fixmate -m /shared5/Alex/Mitochondrial_project/BAMs_horseref/${name}_sorted_mtdna_horseref.bam - | samtools sort - | samtools markdup - /shared5/Alex/Mitochondrial_project/BAMs_horseref/${name}_markdup_sorted_mtdna_horseref.bam &
done
```

Submitted qualimap of Exmoor sequences horse_ref (10:38h):
```bash
cat /shared5/Alex/Mitochondrial_project/List_mitochondrial_merged.txt | grep "Exmoor" | awk '{print $1}' | while read name ; do
  qualimap bamqc -bam /shared5/Alex/Mitochondrial_project/BAMs_horseref/${name}_markdup_sorted_mtdna_horseref.bam --java-mem-size=200G 2>/shared5/Alex/Mitochondrial_project/BAMs_horseref/log_files/${name}_qualimap.log &
done
```

To access OneDrive files on computer:
`ll /drives/C/Users/asus/OneDrive\ -\ University\ of\ Glasgow/UNIVERSITY/5th\ YEAR/Honours\ project/Files/`


## 02/01/2023
I create consensus FASTAs but I get many undetermined bases so I rerun ANGSD but with a minimum depth threshold of 100x


Tried reruning Y-chromosome FASTAs but minimum depth set at 30x:
```bash
cat /shared5/Alex/Y-chromosome_project/Y-chromosome_ID_depth_list.txt | while read line; do
  name=$(echo ${line} | awk '{print $1}');
  number=$(echo ${line} | awk '{print $(NF-1)}');
  maxnumber=$(echo ${line} | awk '{print $(NF)}');
 angsd -i /shared5/Alex/Y-chromosome_project/BAMs/${name}_Y_sorted_markdup.bam -minMapQ 30 -minQ 20 -setMinDepth 30 -setMaxDepth ${maxnumber} -remove_bads 1 -doFasta 2 -doCounts 1 -ref /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/horse_Y_chromosome.fasta -out /shared5/Alex/Y-chromosome_project/FASTAs/${name}_30x_Y &
done
```

I rerun the Y-chromosome FASTAs but with 30x minimum depth and 250x max:
```bash
#consensus fasta
cat /shared5/Alex/Y-chromosome_project/Y-chromosome_ID_depth_list.txt | while read line; do
  name=$(echo ${line} | awk '{print $1}');
  angsd -i /shared5/Alex/Y-chromosome_project/BAMs/${name}_Y_sorted_markdup.bam -minMapQ 30 -minQ 20 -setMinDepth 30 -setMaxDepth 250 -remove_bads 1 -doFasta 2 -doCounts 1 -ref /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/horse_Y_chromosome.fasta -out /shared5/Alex/Y-chromosome_project/FASTAs/${name}_30x_Y &
done

#change headers
cat /shared5/Alex/Y-chromosome_project/Y-chromosome_ID_depth_list.txt | while read line; do
  name=$(echo $line | awk '{print $1}');
  breed=$(echo $line | awk '{print $4}');
  sed -i "s/>MH341179.1/>${name}_${breed}_MH341179.1/g" /shared5/Alex/Y-chromosome_project/FASTAs/${name}_30x_Y.fa
done

#concatente
cat /shared5/Alex/Y-chromosome_project/Y-chromosome_ID_depth_list.txt | awk '{print $1}' | while read name; do
  cat /shared5/Alex/Y-chromosome_project/FASTAs/${name}_30x_Y.fa
done >> Exmoor_horse_30x_Y.fa

#add reference
cat /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/horse_Y_chromosome.fasta Exmoor_horse_30x_Y.fa > tmp.fa && mv tmp.fa Exmoor_horse_30x_Y.fa

#copy
scp studentprojects@young.eng.gla.ac.uk:/shared5/Alex/Y-chromosome_project/FASTAs/Exmoor_horse_30x_Y.fa /drives/C/Users/asus/OneDrive\ -\ University\ of\ Glasgow/UNIVERSITY/5th\ YEAR/Honours\ project/Files/
```



## 16/01/23
Restarted Y-chromosome pipeline today.
Mapped to non-rep MSY assembly: (submitted 12:39h)
```bash
cat /shared5/Alex/Y-chromosome_project/List_ychromosome.txt | awk '{print $2}' | while read file ; do
  location=$(echo ${file} | awk -F "/" '{$NF=""}1' OFS="/") ;
  name=$(echo ${file} | awk -F "/" '{print $NF}' | sed 's/_1_val_1.fq.gz//');
  bwa mem /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/GCA_002166905.2_LipY764_genomic.fna ${location}${name}_1_val_1.fq.gz ${location}${name}_2_val_2.fq.gz > /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764.sam &
done
```

I also redo the Run2, but with a new configure_strelka file that Anubhab has written.
