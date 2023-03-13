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
  bwa mem /shared5/Alex/Y-chromosome_project/_chromosome_reference/_chromosome.fasta ${location}${name}_1_val_1.fq.gz ${location}${name}_2_val_2.fq.gz 2> /shared5/Alex/Y-chromosome_project/BAMs/log_files/${name}.bwamem.log > /shared5/Alex/Y-chromosome_project/BAMs/${name}_Y.sam &
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

## 18/01/23
I map the Sep_2022 Exmoor individuals to the EquCab3 reference: (submitted 12:15)
```bash
cat /shared5/Alex/Inbreeding_project/List_Sep_2022_FASTQ.txt | awk '{print $2}' | while read file ; do
  location=$(echo ${file} | awk -F "/" '{$NF=""}1' OFS="/") ;
  name=$(echo ${file} | awk -F "/" '{print $NF}' | sed 's/_1_val_1.fq.gz//');
  bwa mem /shared5/Alex/Inbreeding_project/horse_WG_reference/GCF_002863925.1_EquCab3.0_genomic.fna ${location}${name}_1_val_1.fq.gz ${location}${name}_2_val_2.fq.gz 2> /shared5/Alex/Inbreeding_project/BAMs/log_files/${name}.bwamem.log > /shared5/Alex/Inbreeding_project/BAMs/${name}_horseref.sam &
done
```

Remapped to non-rep MSY assembly but adding `-M` option: (submitted 12:50h)
```bash
cat /shared5/Alex/Y-chromosome_project/List_ychromosome.txt | awk '{print $2}' | while read file ; do
  location=$(echo ${file} | awk -F "/" '{$NF=""}1' OFS="/") ;
  name=$(echo ${file} | awk -F "/" '{print $NF}' | sed 's/_1_val_1.fq.gz//');
  bwa mem -M /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/GCA_002166905.2_LipY764_genomic.fna ${location}${name}_1_val_1.fq.gz ${location}${name}_2_val_2.fq.gz 2> /shared5/Alex/Y-chromosome_project/BAMs_LipY764/log_files/${name}.bwamem.log > /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764.sam &
done
```


# 19/01/23
Do pipeline to map 4 Exmoor individuals to the WG donkey reference for Run2 (E_102004, E_107013,E_21084, E_21149):
```bash
# Mapping
cat /shared5/Alex/Exmoor/List_Exmoor.txt | awk '{print $2}' | while read file ; do
  location=$(echo ${file} | awk -F "/" '{$NF=""}1' OFS="/") ;
  name=$(echo ${file} | awk -F "/" '{print $NF}' | sed 's/_1_val_1.fq.gz//');
  bwa mem /shared5/Alex/Donkey_ref/GCF_016077325.2_ASM1607732v2_genomic.fna.gz ${location}${name}_1_val_1.fq.gz ${location}${name}_2_val_2.fq.gz 2> /shared5/Alex/Exmoor/BAMs/log_files/${name}.bwamem.log > /shared5/Alex/Exmoor/BAMs/${name}_DonkeyRef.sam &
done

# Sorting
cat /shared5/Alex/Exmoor/List_Exmoor.txt | awk '{print $1}' | while read name ; do
  samtools view -q 30 -u /shared5/Alex/Exmoor/BAMs/${name}_DonkeyRef.sam | samtools sort -n -o /shared5/Alex/Exmoor/BAMs/${name}_sorted_DonkeyRef.bam &
done

# Markdup
cat /shared5/Alex/Exmoor/List_Exmoor.txt | awk '{print $1}' | while read name ; do
  samtools fixmate -m /shared5/Alex/Exmoor/BAMs/${name}_sorted_DonkeyRef.bam - | samtools sort - | samtools markdup - /shared5/Alex/Exmoor/BAMs/${name}_sorted_markdup_DonkeyRef.bam &
done

# Indexing
samtools index *markdup.bam

# Qualimap
cat /shared5/Alex/Exmoor/List_Exmoor.txt | awk '{print $1}' | while read name; do
  qualimap bamqc -bam /shared5/Alex/Exmoor/BAMs/${name}_sorted_markdup_DonkeyRef.bam --java-mem-size=200G 2>/shared5/Alex/Exmoor/BAMs/log_files/${name}_qualimap.log &
done
```


## 20/01/23
Continued pipieline of Y chromosome pipeline: (submitted 12:56h)
```bash
cat /shared5/Alex/Y-chromosome_project/List_ychromosome.txt | awk '{print $1}' | while read name ; do
  samtools view -h -b -F 4 -u /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764.sam |   samtools sort -o /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted.bam &
done
```



# 23/01/23
Continued pipieline of Y chromosome pipeline. Removed duplicates: (submitted 13:31h)
```bash
cat /shared5/Alex/Y-chromosome_project/List_ychromosome.txt | awk '{print $1}' | while read name ; do
  java -jar /shared5/Alex/picard-2.27.5/picard.jar MarkDuplicates -I /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted.bam -O /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted_rmdp.bam -M /shared5/Alex/Y-chromosome_project/BAMs_LipY764/log_files/${name}_markdup_metrics.txt --REMOVE_DUPLICATES true --ASSUME_SORT_ORDER coordinate --VALIDATION_STRINGENCY SILENT &
done
```

Continued pipeline of Inbreeding 2022 Exmoor genomes. Sorting (submitted 14:05h):
```bash
cat /shared5/Alex/Inbreeding_project/List_Sep_2022_FASTQ.txt | awk '{print $1}' | while read name ; do
  samtools view -h -b -F 4 -q 30 -u /shared5/Alex/Inbreeding_project/BAMs/${name}_horseref.sam |  samtools sort -o /shared5/Alex/Inbreeding_project/BAMs/${name}_sorted_horseref.bam &
done
```

# 24/01/23
Whilst running remove duplicates of the Ychromosome pipeline, some of the files failed because they ran out of memory. Error states:
`# There is insufficient memory for the Java Runtime Environment to continue.
# Native memory allocation (mmap) failed to map 8048869376 bytes for committing reserved memory.`

I will rerun the code on all the samples but this time reducing the allocated memory using `-Xmx4g` which sets memory to 4 Gb: (submitted 12:31h)
```bash
cat /shared5/Alex/Y-chromosome_project/List_ychromosome.txt | awk '{print $1}' | while read name ; do
  java -Xmx4g -jar /shared5/Alex/picard-2.27.5/picard.jar MarkDuplicates -I /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted.bam -O /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted_rmdp.bam -M /shared5/Alex/Y-chromosome_project/BAMs_LipY764/log_files/${name}_markdup_metrics.txt --REMOVE_DUPLICATES true --ASSUME_SORT_ORDER coordinate --VALIDATION_STRINGENCY SILENT --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000 &
done
```

Then I use samtools merge for Ychromsome pipeline. There are 4 indivudals that need to be merged:
```bash
samtools merge E_23279_merged_LipY764_sorted_rmdp.bam E_23279_EKDN220034363-1A_H5J5MDSX5_L2_LipY764_sorted_rmdp.bam E_23279_EKDN220034363-1A_H7VKMDSX5_L3_LipY764_sorted_rmdp.bam;

samtools merge E_23416_merged_LipY764_sorted_rmdp.bam E_23416_EKDN220034368-1A_H5J5MDSX5_L2_LipY764_sorted_rmdp.bam E_23416_EKDN220034368-1A_H7VKMDSX5_L1_LipY764_sorted_rmdp.bam E_23416_EKDN220034368-1A_H7VL2DSX5_L2_LipY764_sorted_rmdp.bam;

samtools merge E_320005_merged_LipY764_sorted_rmdp.bam E_320005_EKDN220034373-1A_H5J5MDSX5_L2_LipY764_sorted_rmdp.bam E_320005_EKDN220034373-1A_H7VKMDSX5_L1_LipY764_sorted_rmdp.bam E_320005_EKDN220034373-1A_H7VL2DSX5_L4_LipY764_sorted_rmdp.bam;

samtools merge E_32023_merged_LipY764_sorted_rmdp.bam E_32023_EKDN220034372-1A_H5J5MDSX5_L2_LipY764_sorted_rmdp.bam E_32023_EKDN220034372-1A_H7VKMDSX5_L1_LipY764_sorted_rmdp.bam
```
Then I submit quality filtering for Ychromosome pipeline: (submitted 15:39h)
```bash
cat /shared5/Alex/Y-chromosome_project/List_Y-chromosome_merged.txt | awk '{print $1}' | while read name ; do
  samtools view -h -b -F 4 -q 20 /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted_rmdp.bam -o /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted_rmdp_MQ20.bam &
done
```

Continued pipeline of Inbreeding 2022 Exmoor genomes. Markdup (submitted 12:45h):
```bash
cat /shared5/Alex/Inbreeding_project/List_Sep_2022_FASTQ.txt | awk '{print $1}' | while read name ; do
  java -Xmx4g -jar /shared5/Alex/picard-2.27.5/picard.jar MarkDuplicates -I /shared5/Alex/Inbreeding_project/BAMs/${name}_sorted_horseref.bam -O /shared5/Alex/Inbreeding_project/BAMs/${name}_rmdp_sorted_horseref.bam -M /shared5/Alex/Inbreeding_project/BAMs/log_files/${name}_markdup_metrics.txt --REMOVE_DUPLICATES true --ASSUME_SORT_ORDER coordinate --VALIDATION_STRINGENCY SILENT --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000 &
done
```
I try running ValidateSamFile to see what the issue is:
```bash
cat /shared5/Alex/Inbreeding_project/List_Sep_2022_FASTQ.txt | awk '{print $1}' | while read name ; do
  java -jar /shared5/Alex/picard-2.27.5/picard.jar ValidateSamFile -I /shared5/Alex/Inbreeding_project/BAMs/${name}*sam -O  /shared5/Alex/Inbreeding_project/BAMs/${name}_sam_horseref_validate.txt --MODE SUMMARY -R /shared5/Alex/Inbreeding_project/horse_WG_reference/GCF_002863925.1_EquCab3.0_genomic.fna &
done
```


I try rerunning it for just a single file.
java -Xmx4g -jar /shared5/Alex/picard-2.27.5/picard.jar MarkDuplicates -I /shared5/Alex/Inbreeding_project/BAMs/E_900588_EKDN220034362-1A_H7VKVDSX5_L2_sorted_horseref.bam -O /shared5/Alex/Inbreeding_project/BAMs/E_900588_EKDN220034362-1A_H7VKVDSX5_L2_rmdp_sorted_horseref.bam -M /shared5/Alex/Inbreeding_project/BAMs/log_files/E_900588_EKDN220034362-1A_H7VKVDSX5_L2_markdup_metrics.txt --REMOVE_DUPLICATES true --ASSUME_SORT_ORDER coordinate --VALIDATION_STRINGENCY SILENT --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000

Still get error for the Inbreeding genomes.



# 25/01/23
Continue pipeline for Ychromosome. Indexing and Qualimap. (submitted 12:30h):
```bash
#Indexing
cat /shared5/Alex/Y-chromosome_project/List_Y-chromosome_merged.txt | awk '{print $1}' | while read name ; do
  samtools index /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted_rmdp_MQ20.bam &
done

#Qualimap
cat /shared5/Alex/Y-chromosome_project/List_Y-chromosome_merged.txt | awk '{print $1}' | while read name ; do
  qualimap bamqc -bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted_rmdp_MQ20.bam --java-mem-size=200G 2>/shared5/Alex/Y-chromosome_project/BAMs_LipY764/log_files/${name}_qualimap.log &
done

#creating depth list
cat /shared5/Alex/Y-chromosome_project/List_Y-chromosome_merged.txt | awk '{print $1}' | while read name; do
  coverage=$(cat /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}*stats/genome_results.txt | grep "mean coverageData"  | awk '{print $4}' | sed 's/X//' | sed 's/,//' | cut -d "." -f 1);
  paste <(echo ${name}) <(echo ${coverage}) <(echo "${coverage} * 2" | bc);
done > /shared5/Alex/Y-chromosome_project/Y-chromosome_ID_depth_list_new.txt

#consensus FASTAs
cat /shared5/Alex/Y-chromosome_project/Y-chromosome_ID_depth_list_new.txt | while read line; do
  name=$(echo ${line} | awk '{print $1}');
  number=$(echo ${line} | awk '{print $(NF-1)}');
  maxnumber=$(echo ${line} | awk '{print $(NF)}');
  angsd -i /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted_rmdp_MQ20.bam -minMapQ 30 -minQ 20 -setMinDepth ${number} -setMaxDepth ${maxnumber} -remove_bads 1 -doFasta 2 -doCounts 1 -ref /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/GCA_002166905.2_LipY764_genomic.fna -out /shared5/Alex/Y-chromosome_project/FASTAs_LipY764/${name}_LipY764 &
done
```
I rerun the consensus FASTA but at 30x beacause I was getting no bases. Still get no bases.



# 26/01/23
Wrote up some of the mtDNA stuff for the report.

Now I want to work on the Y-chromosome BAMs and try to transform them to FASTAs.
```bash
# Attempt 1: using HaplotypeCaller (based on Felkel et al., 2019)
cat /shared5/Alex/Y-chromosome_project/List_Y-chromosome_merged.txt | awk '{print $1}' | while read name; do
  gatk /home/opt/miniconda2/pkgs/gatk-3.8-5/bin/GenomeAnalysisTK -R /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/GCA_002166905.2_LipY764_genomic.fna –I /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted_rmdp_MQ20.bam -O /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764.vcf
done

java -jar GenomeAnalysisTK.jar -T CombineGVCFs -V sample1.g.vcf -V sample2.g.vcf
–V sampleX.g.vcf -o cohort.g.vcf

#Gives ERROR --> bash: gatk: command not found

# Attempt 2: Using post from biostars post of BAM to FASTA
cat /shared5/Alex/Y-chromosome_project/List_Y-chromosome_merged.txt | awk '{print $1}' | while read name; do
  samtools mpileup -f /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/GCA_002166905.2_LipY764_genomic.fna /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted_rmdp_MQ20.bam | bcftools call -c | vcfutils.pl vcf2fq > /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_cns.fastq &&
  seqtk seq -aQ64 -q20 -n N SAMPLE_cns.fastq > SAMPLE_cns.fasta
done
#Gives error --> Use of uninitialized value $l in numeric lt (<) at /home/opt/miniconda2/bin/vcfutils.pl line 566.


# Attempt 3: Using mpileup (based on Wallner et al., 2017)
cat /shared5/Alex/Y-chromosome_project/List_Y-chromosome_merged.txt | awk '{print $1}' | while read name; do
  samtools view -h /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted_rmdp_MQ20.bam | samtools view -bS - | samtools mpileup -A -u -f /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/GCA_002166905.2_LipY764_genomic.fna - | bcftools view -vcg - | vcfutils.pl varFilter -d 2 -D 200 > /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted_rmdp_MQ20.vcf
done
#Gives error --> [E::cg] unknown type. Accepted types are snps, indels, mnps, other
```

Redid part of Inbreeding pipeline but using samtools fixmate and markduplicates because Picardtools wasn't working: (submitted 15:49h)
```bash
cat /shared5/Alex/Inbreeding_project/List_Sep_2022_FASTQ.txt | awk '{print $1}' | while read name ; do
  samtools view -h -u -q 30 /shared5/Alex/Inbreeding_project/BAMs/${name}_horseref.sam |  samtools sort -n -o /shared5/Alex/Inbreeding_project/BAMs/${name}_sorted_horseref.bam &
done
```



# 27/01/23
I try getting the consensus FASTAs for the Ychromosome again but this time with the minimum depth set at 5x:
```bash
#consensus FASTAs
cat /shared5/Alex/Y-chromosome_project/Y-chromosome_ID_depth_list_new.txt | while read line; do
  name=$(echo ${line} | awk '{print $1}');
  maxnumber=$(echo ${line} | awk '{print $(NF-1)}'); #CHANGE THIS
  angsd -i /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted_rmdp_MQ20.bam -minMapQ 30 -minQ 20 -setMinDepth 5 -setMaxDepth ${maxnumber} -remove_bads 1 -doFasta 2 -doCounts 1 -ref /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/GCA_002166905.2_LipY764_genomic.fna -out /shared5/Alex/Y-chromosome_project/FASTAs_LipY764/${name}_LipY764 &
done
```
Still not getting any coverage. I submitted the same consensus BAM but with the sorted_rmdp.bam files instead of the sorted_rmdp_MQ20.bam


Submitted markdup of inbreeding genomes: (submitted 18:31h)
```bash
cat /shared5/Alex/Inbreeding_project/List_Sep_2022_FASTQ.txt | awk '{print $1}' | while read name ; do
  samtools fixmate -m /shared5/Alex/Inbreeding_project/BAMs/${name}_sorted_horseref.bam - | samtools sort - | samtools markdup - /shared5/Alex/Inbreeding_project/BAMs/${name}_markdup_sorted_horseref.bam &
done
```


# 29/01/23
To compare read count of the old and new Ychromosome:
```bash
for file in /shared5/Alex/Y-chromosome_project/BAM*/*sorted_*bam; do
  paste <(samtools view -F 4 -c $file) <(echo $file)
done
```

Comparing read count examples:
- E_102004 Original file --> 4341253 // New file --> 4603163
  From the qualimap file, the original Ychromsome BAM shows:
  "mean coverageData = 35.2535X
    std coverageData = 252.3167X

    There is a 58.84% of reference with a coverageData >= 1X
    There is a 54.57% of reference with a coverageData >= 2X
    There is a 49.78% of reference with a coverageData >= 3X
    There is a 45.14% of reference with a coverageData >= 4X
    There is a 40.67% of reference with a coverageData >= 5X
    There is a 36.63% of reference with a coverageData >= 6X
    There is a 32.98% of reference with a coverageData >= 7X
    There is a 29.77% of reference with a coverageData >= 8X
    There is a 26.98% of reference with a coverageData >= 9X
    There is a 24.63% of reference with a coverageData >= 10X"

    From the new Ychromosome BAM qualimap file:
    " mean coverageData = 59.8409X
     std coverageData = 425.015X

     There is a 95.09% of reference with a coverageData >= 1X
     There is a 87.54% of reference with a coverageData >= 2X
     There is a 78.47% of reference with a coverageData >= 3X
     There is a 69.42% of reference with a coverageData >= 4X
     There is a 60.82% of reference with a coverageData >= 5X
     There is a 53.26% of reference with a coverageData >= 6X
     There is a 46.79% of reference with a coverageData >= 7X
     There is a 41.51% of reference with a coverageData >= 8X
     There is a 37.15% of reference with a coverageData >= 9X
     There is a 33.65% of reference with a coverageData >= 10X
"

It is confusing because the new Ychromosome files appear to have a higher coverage but the FASTAs do not show any bases. Whereas the old Ychromsome files have a lower coverage but actually produce usable FASTAs, I guess this might be an issue with ANGSD. In the meantime, I will run ANGSD conensus FASTA for the old Ychromosome files but at a lower coverage (5x):
```bash
# Consensus FASTAs
cat /shared5/Alex/Y-chromosome_project/Y-chromosome_ID_depth_list.txt | while read line; do
  name=$(echo ${line} | awk '{print $1}');
  maxnumber=$(echo ${line} | awk '{print $(NF)}');

  angsd -i /shared5/Alex/Y-chromosome_project/BAMs/${name}_Y_sorted_markdup.bam -minMapQ 30 -minQ 20 -setMinDepthInd 5 -setMaxDepth ${maxnumber} -remove_bads 1 -doFasta 2 -doCounts 1 -ref /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/MH341179_Ychromosome.fasta -out /shared5/Alex/Y-chromosome_project/FASTAs/${name}_5x_Y &
done

# checking reads from FASTA
for file in *5x*fa; do
  paste <(grep -v "NNN" $file | wc -l ) <(wc -l $file)
done
```
Weirdly this now produces no results in the FASTA. So I'm going to leave this aside for now.


Continuing pipeline for Inbreeding genomes. Merging files:
```bash
#loop to print out code of files that need to be merged:
ll *markdup*bam | awk '{print $9}' | awk '{print substr($0,1,10)}' | uniq -c  | grep -v " 1 " | awk '{print $2}' | cut -d "_" -f 1,2 | while read file; do
  paste <(echo "samtools merge ${file}_merged_markdup_sorted_horseref.bam") <(ls ${file}*markdup*bam) ;
done
 #output:
samtools merge E_23279_merged_sorted_horseref.bam E_23279_EKDN220034363-1A_H5J5MDSX5_L2_markdup_sorted_horseref.bam E_23279_EKDN220034363-1A_H7VKMDSX5_L3_markdup_sorted_horseref.bam

#run code to merge files:
samtools merge E_23279_merged_markdup_sorted_horseref.bam E_23279_EKDN220034363-1A_H5J5MDSX5_L2_markdup_sorted_horseref.bam E_23279_EKDN220034363-1A_H7VKMDSX5_L3_markdup_sorted_horseref.bam; samtools merge E_23416_merged_markdup_sorted_horseref.bam E_23416_EKDN220034368-1A_H5J5MDSX5_L2_markdup_sorted_horseref.bam  E_23416_EKDN220034368-1A_H7VKMDSX5_L1_markdup_sorted_horseref.bam E_23416_EKDN220034368-1A_H7VL2DSX5_L2_markdup_sorted_horseref.bam; samtools merge E_23434_merged_markdup_sorted_horseref.bam E_23434_EKDN220034361-1A_H5J5MDSX5_L2_markdup_sorted_horseref.bam E_23434_EKDN220034361-1A_H7VKMDSX5_L3_markdup_sorted_horseref.bam; samtools merge E_320005_merged_markdup_sorted_horseref.bam E_320005_EKDN220034373-1A_H5J5MDSX5_L2_markdup_sorted_horseref.bam E_320005_EKDN220034373-1A_H7VKMDSX5_L1_markdup_sorted_horseref.bam E_320005_EKDN220034373-1A_H7VL2DSX5_L4_markdup_sorted_horseref.bam; samtools merge E_32023_merged_markdup_sorted_horseref.bam E_32023_EKDN220034372-1A_H5J5MDSX5_L2_markdup_sorted_horseref.bam E_32023_EKDN220034372-1A_H7VKMDSX5_L1_markdup_sorted_horseref.bam; samtools merge E_479023_merged_markdup_sorted_horseref.bam E_479023_EKDN220034374-1A_H5J5MDSX5_L2_markdup_sorted_horseref.bam  E_479023_EKDN220034374-1A_H7VKMDSX5_L1_markdup_sorted_horseref.bam  E_479023_EKDN220034374-1A_H7VL2DSX5_L2_markdup_sorted_horseref.bam; samtools merge E_512001_merged_markdup_sorted_horseref.bam E_512001_EKDN220034365-1A_H5J5MDSX5_L2_markdup_sorted_horseref.bam  E_512001_EKDN220034365-1A_H7VKMDSX5_L3_markdup_sorted_horseref.bam; samtools merge E_78170_merged_markdup_sorted_horseref.bam E_78170_EKDN220034369-1A_H5J5MDSX5_L2_markdup_sorted_horseref.bam E_78170_EKDN220034369-1A_H7VKMDSX5_L1_markdup_sorted_horseref.bam  E_78170_EKDN220034369-1A_H7VL2DSX5_L4_markdup_sorted_horseref.bam; samtools merge E_900588_merged_markdup_sorted_horseref.bam E_900588_EKDN220034362-1A_H5J5MDSX5_L2_markdup_sorted_horseref.bam E_900588_EKDN220034362-1A_H7VKVDSX5_L2_markdup_sorted_horseref.bam; samtools merge E_900717_merged_markdup_sorted_horseref.bam E_900717_EKDN220034357-1A_H5CYKDSX5_L4_markdup_sorted_horseref.bam  E_900717_EKDN220034357-1A_H5HJMDSX5_L1_markdup_sorted_horseref.bam
```


# 30/01/23
Continuing pipeline for Inbreeding genomes. Indexing and Qualimap.
```bash
# Indexing
for file in /shared5/Alex/Inbreeding_project/BAMs*markdup*bam*; do
  samtools index $file
done

# Qualimap
for file in *markdup*bam; do
  qualimap bamqc -bam $file --java-mem-size=200G &
done
```


# 31/01/23
Qualimap of inbreeding genomes:
```bash
# Checking qualimap files
cat /shared5/Alex/Inbreeding_project/List_Sep_2022_merged.txt | while read name; do
  paste <(echo ${name}) <(ls -lh /shared5/Alex/Inbreeding_project/BAMs/${name}*stats/genome_results.txt)
  grep "mean coverageData" /shared5/Alex/Inbreeding_project/BAMs/${name}*stats/genome_results.txt
done

# Error --> Many haven't been created due to "Insufficient memory", so I rerun qualimap but without the java-mem-size and with number of threads set to 6. (submitted 13:44):
cat /shared5/Alex/Inbreeding_project/List_Sep_2022_merged.txt | while read file; do
  qualimap bamqc -nt 6 -bam ${file}_markdup_sorted_horseref.bam &
done

# `-nt 6` causes insufficient memory error. So I rerun qualimap again:
cat /shared5/Alex/Inbreeding_project/List_Sep_2022_merged.txt | while read file; do
  qualimap bamqc -bam ${file}_markdup_sorted_horseref.bam --java-mem-size=200G  &
done

#Qualimap was still giving memory errors so instead I ran the code again but isntead of simultaneously, I made it loop one at a time:
cat /shared5/Alex/Inbreeding_project/List_Sep_2022_merged.txt | while read file; do
  qualimap bamqc -bam ${file}_markdup_sorted_horseref.bam --java-mem-size=200G
done
```


# 01/02/23
Submitted ValidateSamFile for the Inbreeding genomes.
```bash
cat /shared5/Alex/Inbreeding_project/List_Sep_2022_merged.txt | while read file; do
  mkdir -p /shared5/Alex/Inbreeding_project/BAMs/ValidateSamFile/
  java -jar /shared5/Alex/picard-2.27.5/picard.jar ValidateSamFile -I /shared5/Alex/Inbreeding_project/BAMs/${file}_markdup_sorted_horseref.bam -O /shared5/Alex/Inbreeding_project/BAMs/ValidateSamFile/${file}_markdup_sorted_horseref_errors.txt -R /shared5/Alex/Inbreeding_project/horse_WG_reference/GCF_002863925.1_EquCab3.0_genomic.fna --MODE SUMMARY &
done
```
ValidateSamFile info text file not being created due to error:
"ERROR   2023-02-01 09:32:24     ValidateSamFile Value was put into PairInfoMap more than once.  1: A00551:538:H5J5MDSX5:2:1542:9118:22106"


Y-CHROMOSOME --> I talked with Anubhab about not being to create any consensus FASTAs and he say's that the variable "maxnumber" is not calling the values I want. This is because the final column in the "Y-chromosome_ID_depth_list.txt" was formed by breed names but some had spaces inbetween them so that was giving errors. I fixed this but substitituing spaces with underscores and now the correct values are being called so I will run the consensus FASTAs for both the new and old Y-chromosome files again.
```bash
#Old Ychromosome files
cat /shared5/Alex/Y-chromosome_project/Y-chromosome_ID_depth_list.txt | while read line; do
  name=$(echo ${line} | awk '{print $1}');
  maxnumber=$(echo ${line} | awk '{print $(NF-1)}');
  angsd -i /shared5/Alex/Y-chromosome_project/BAMs/${name}_Y_sorted_markdup.bam -minMapQ 30 -minQ 20 -setMinDepth 5 -setMaxDepth ${maxnumber} -remove_bads 1 -doFasta 2 -doCounts 1 -ref /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/MH341179_Ychromosome.fasta -out /shared5/Alex/Y-chromosome_project/FASTAs/${name}_5x_horseref_Y &
done

#New Ychromsome Files
cat /shared5/Alex/Y-chromosome_project/Y-chromosome_ID_depth_list_new.txt | while read line; do
  name=$(echo ${line} | awk '{print $1}');
  maxnumber=$(echo ${line} | awk '{print $(NF-1)}');
  angsd -i /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}_LipY764_sorted_rmdp_MQ20.bam -minMapQ 30 -minQ 20 -setMinDepth 5 -setMaxDepth ${maxnumber} -remove_bads 1 -doFasta 2 -doCounts 1 -ref /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/GCA_002166905.2_LipY764_genomic.fna -out /shared5/Alex/Y-chromosome_project/FASTAs_LipY764/${name}_5x_LipY764_Y &
done
```
Fortunately, the FASTAs created now are much better and have actual coverage.
```bash
for file in *5x*fa; do
  paste <(grep -v "NNN" $file | wc -l ) <(wc -l $file)
done

# Output of old Ychromosome files:
38938   189555 E_102004_EKDN220034353-1A_H5GL5DSX5_L2_5x_horseref_Y.fa
28685   189555 E_107013_EKDN220034352-1A_H5GL5DSX5_L2_5x_horseref_Y.fa
39416   189555 E_21084_EKDN220034367-1A_H5J5MDSX5_L2_5x_horseref_Y.fa
(...)

# Output of new Ychromosome files:
45191   130390 E_102004_EKDN220034353-1A_H5GL5DSX5_L2_5x_LipY764_Y.fa
26881   130390 E_107013_EKDN220034352-1A_H5GL5DSX5_L2_5x_LipY764_Y.fa
52996   130399 E_23279_merged_5x_LipY764_Y.fa
```

And checking qualimap we can see how much each Ychromosome file covers:
```bash
#check coverage of old Ychromosome files
cat /shared5/Alex/Y-chromosome_project/Y-chromosome_ID_depth_list_new.txt | awk '{print $1}' | while read name; do
  paste <(echo $name) <(grep  "coverageData >= 5X" /shared5/Alex/Y-chromosome_project/BAMs/${name}*stats/*txt)
done

#check coverage of new Ychromosome files
cat /shared5/Alex/Y-chromosome_project/Y-chromosome_ID_depth_list_new.txt | awk '{print $1}' | while read name; do
  paste <(echo $name) <(grep  "coverageData >= 5X" /shared5/Alex/Y-chromosome_project/BAMs_LipY764/${name}*stats/*txt)
done
```

I then concatenate the old Ychromosome files into respective compliated FASTAs so I can analyze them with DNAsp


# 05/02/23
Concatenating new Ychromosome files:
```bash

```

USEFUL separate words with TAB: tail head.txt | awk -v OFS="\t" '{$1=$1; print}'
USEFUL --> see available conda virtual environments: `conda info --envs`


# 13/02/23
In past few days, I have merged the vcfs of the Exmoors and the horses. And I am just started filtering the Exmoor vcfs:
```bash
#merged exmoor vcf
vcftools --vcf /shared5/Alex/Inbreeding_project/variants/merged_exmoor_8-2-2023.vcf --remove-indels --minQ 30 --minGQ 30 --out /shared5/Alex/Inbreeding_project/variants/merged_exmoor_8-2-2023_rmvIndels_minQ30_minGQ30 --recode

#merged horse vcf
vcftools --vcf /shared5/Alex/Run1/merged_donkey_8-2-2023.vcf --remove-indels --minQ 30 --minGQ 30 --out /shared5/Alex/Run1/merged_donkey_8-2-2023_rmvIndels_minQ30_minGQ30 --recode
```

Currently working on variant calling for the Ychromosome following reccomendations by Anubhab:
```bash
bcftools call -c -O v --ploidy-file /shared5/Alex/Y-chromosome_project/variants_bcftools/ploidy.txt -o /shared5/Alex/Y-chromosome_project/variants_bcftools/Exmoors_2022_Y_oldref.vcf /shared5/Alex/Y-chromosome_project/variants_bcftools/Exmoors_2022_Y_oldref.pileup
```



# 15/02/23
Did 2nd step of VCF filtering of Exmoor and horse VCFs:
```bash
#Exmoor ponies
vcftools --vcf /shared5/Alex/Inbreeding_project/variants/merged_exmoor_8-2-2023_rmvIndels_minQ30_minGQ30.recode.vcf --mac 3 --hwe 0.05 --remove-filtered-all --out /shared5/Alex/Inbreeding_project/variants/merged_exmoor_8-2-2023_rmvIndels_minQ30_minGQ30_mac3_hwe0.05_PASSonly --recode

#Horses
vcftools --vcf /shared5/Alex/Run1/merged_donkey_8-2-2023_rmvIndels_minQ30_minGQ30.recode.vcf --mac 3 --hwe 0.05 --remove-filtered-all --out /shared5/Alex/Run1/merged_donkey_8-2-2023_rmvIndels_minQ30_minGQ30_mac3_hwe0.05_PASSonly --recode
```

Ychromosome steps:
```bash
#Variant calling for newref Ychromsome:
bcftools call -c -O v --ploidy-file /shared5/Alex/Y-chromosome_project/variants_bcftools/ploidy.txt -o /shared5/Alex/Y-chromosome_project/variants_bcftools/Exmoors_2022_Y_oldref.vcf /shared5/Alex/Y-chromosome_project/variants_bcftools/Exmoors_2022_Y_oldref.pileup

#1st step of VCF filtering of newref Ychromosome:
vcftools --vcf /shared5/Alex/Y-chromosome_project/variants_bcftools/Exmoors_2022_Y_oldref.vcf --remove-indels --minQ 30 --minGQ 30 --out /shared5/Alex/Y-chromosome_project/variants_bcftools/Exmoors_2022_Y_oldref_rmvIndels_minQ30_minGQ30 --recode
#After filtering, kept 5677718 out of a possible 6265789 Sites

#1st step of VCF filtering of oldref Ychromosome:
vcftools --vcf /shared5/Alex/Y-chromosome_project/variants_bcftools/Exmoors_2022_Y_oldref.vcf --remove-indels --minQ 30 --minGQ 30 --out /shared5/Alex/Y-chromosome_project/variants_bcftools/Exmoors_2022_Y_oldref_rmvIndels_minQ30_minGQ30 --recode
#After filtering, kept 6291005 out of a possible 6584514 Sites
```

3rd step of VCF filtering - HAVEN'T RUN THIS YET BECAUSE THE SITES RETAINED WITH --MAXMISSING AT 0.8 ARE TOO LOW :
```bash
#Exmoors
vcftools --vcf /shared5/Alex/Inbreeding_project/variants/merged_exmoor_8-2-2023_rmvIndels_minQ30_minGQ30_mac3_hwe0.05_PASSonly.recode.vcf --max-missing 0.8 --out /shared5/Alex/Inbreeding_project/variants/merged_exmoor_8-2-2023_rmvIndels_minQ30_minGQ30_mac3_hwe0.05_PASSonly_maxmissing0.8 --recode
#after running this code without output and --max-missing 0.8 --> After filtering, kept 748 out of a possible 6710250 Sites
#after running this code without output and --max-missing 0.5 --> After filtering, kept 161057 out of a possible 6710250 Sites


#Horses
vcftools --vcf /shared5/Alex/Run1/merged_donkey_8-2-2023_rmvIndels_minQ30_minGQ30_mac3_hwe0.05_PASSonly.recode.vcf --max-missing 0.8 --out /shared5/Alex/Run1/merged_donkey_8-2-2023_rmvIndels_minQ30_minGQ30_mac3_hwe0.05_PASSonly_maxmissing0.8 --recode
#after running this code without output and --max-missing 0.8 --> After filtering, kept 4707998 out of a possible 12726242 Sites
#after running this code without output and --max-missing 0.5 -->
```



# 16/02/23
I copy the oldref semi-filtered Ychromosome VCF to see if I can get it to work in DNAsp:
```bash
scp studentprojects@poland.eng.gla.ac.uk:/shared5/Alex/Y-chromosome_project/variants_bcftools/Exmoors_2022_Y_oldref_rmvIndels_minQ30_minGQ30.recode.vcf /drives/C/Users/asus/OneDrive\ -\ University\ of\ Glasgow/UNIVERSITY/5th\ YEAR/Honours\ project/Files/Y-chromosome/
```
I don't know why but `scp` doesn't appear to be working and it doesn't copy anything, so instead I download a couple of example VCFs from the internet:
USEFUL --> download internet files --> `wget "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/pilot_data/release/2010_07/trio/snps/trio.2010_06.ychr.sites.vcf.gz"`
However DNAsp doesn't work with VCFs



# 17/02/23
Checking mtDNA haplotype networks.


# 23/02/23
Installing "multiqc". I cannot install it with conda or git so I manually install it:
```bash
curl -LOk https://github.com/ewels/MultiQC/archive/master.zip
unzip master.zip
cd MultiQC-master
pip install .
```
Trying to run MultiQC:
`/shared5/Alex/software/MultiQC-master/multiqc  /shared5/Alex/Exmoor_sequencing_data/Sep_2020`º
Doesn't work --> error --> "command not found"


REDO INBREEDING PIPELINE WITH ALL EXMOORS:
Mapping. Submitted (11:34h)
```bash
#mapping
cat /shared5/Alex/Inbreeding_project/List_allexmoors_FASTQ.txt | awk '{print $2}' | while read file ; do
  location=$(echo ${file} | awk -F "/" '{$NF=""}1' OFS="/") ;
  name=$(echo ${file} | awk -F "/" '{print $NF}' | sed 's/_1_val_1.fq.gz//');
  bwa mem /shared5/Alex/Inbreeding_project/horse_WG_reference/GCF_002863925.1_EquCab3.0_genomic.fna   ${location}${name}_1_val_1.fq.gz   ${location}${name}_2_val_2.fq.gz   2> /shared5/Alex/Inbreeding_project/BAMs/log_files/${name}.bwamem.log > /shared5/Alex/Inbreeding_project/BAMs/${name}_horseref.sam &
done

#sorting
cat /shared5/Alex/Inbreeding_project/List_allexmoors_FASTQ.txt | awk '{print $1}' | while read name ; do
  samtools view -h -u -q 30 /shared5/Alex/Inbreeding_project/BAMs/${name}_horseref.sam |  samtools sort -n -o /shared5/Alex/Inbreeding_project/BAMs/${name}_sorted_horseref.bam &
done

#markdup
cat /shared5/Alex/Inbreeding_project/List_Sep_2022_FASTQ.txt | awk '{print $1}' | while read name ; do
  samtools fixmate -m /shared5/Alex/Inbreeding_project/BAMs/${name}_sorted_horseref.bam - | samtools sort - | samtools markdup - /shared5/Alex/Inbreeding_project/BAMs/${name}_markdup_sorted_horseref.bam &
done
```

REDO Y-chromosome newref but with only Exmoor individuals:
```bash
bcftools mpileup -f /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/GCA_002166905.2_LipY764_genomic.fna /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_102004_EKDN220034353-1A_H5GL5DSX5_L2_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_107013_EKDN220034352-1A_H5GL5DSX5_L2_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_21084_EKDN220034367-1A_H5J5MDSX5_L2_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_23279_merged_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_23416_merged_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_320005_merged_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_32023_merged_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_49031_EKDN220034370-1A_H5J5MDSX5_L2_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_519005_EKDN220034364-1A_H5J5MDSX5_L2_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/s2010_EDSW210003772-1a_H3WNKDSX2_L3_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/s2012_EDSW210003766-1a_H3WHWDSX2_L4_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/s479013_EDSW210003775-1a_H3WNKDSX2_L3_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/s49097_EDSW210003769-1a_H3FTWDSX2_L2_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/s49124_EDSW210003767-1a_H3WHWDSX2_L2_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/S276023_EDSW200015547-1a2a_HJFWVDSXY_L3L4_LipY764_sorted_rmdp_MQ20.bam  > /shared5/Alex/Y-chromosome_project/variants_bcftools_onlyexmoors/Exmoors_2022_Y_newref.pileup

bcftools call -c -O v --ploidy-file /shared5/Alex/Y-chromosome_project/variants_bcftools_onlyexmoors/ploidy_newref.txt -o /shared5/Alex/Y-chromosome_project/variants_bcftools_onlyexmoors/Exmoors_2022_Y_newref.vcf /shared5/Alex/Y-chromosome_project/variants_bcftools_onlyexmoors/Exmoors_2022_Y_newref.pileup

#1st step of vcffiltering
vcftools --vcf /shared5/Alex/Y-chromosome_project/variants_bcftools_onlyexmoors/Exmoors_2022_Y_newref.vcf --remove-indels --minQ 30 --minGQ 30 --out /shared5/Alex/Y-chromosome_project/variants_bcftools_onlyexmoors/Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30 --recode

#2nd step of vcfftilering
vcftools --vcf Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30.recode.vcf --min-alleles 2 --mac 3 --out Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30_biallelic_mac3 --recode

#create plink file
vcftools --vcf Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30_biallelic_mac3.recode.vcf --plink --out Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30_biallelic_mac3

#create pca file
plink --file Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30_biallelic_mac3 --pca --out Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30_biallelic_mac3
```

Check:

- make PCA with Y chromosome data
- markduplicates of INbreeding pipeline




To access OneDrive files on computer:
cd /drives/C/Users/asus/OneDrive\ -\ University\ of\ Glasgow/UNIVERSITY/5th\ YEAR/Honours\ project/Files/
