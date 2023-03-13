## 9. Variant calling with Bcftools
Variant calling again but with only the Exmoor indiviudals:

#Anubhabs pipeline:
source activate popgen
bcftools mpileup -f REF BAM1 BAM2 BAM3  > name.pileup
bcftools call -c -O v --ploidy-file ploidy.txt -o name.vcf name.pileup
#ploIDY.TXT contains chromosomeY 1 ENDPOS M 1

```bash
#bcftools mpileup:
bcftools mpileup -f /shared5/Alex/Y-chromosome_project/horse_Y_chromosome_reference/GCA_002166905.2_LipY764_genomic.fna /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_102004_EKDN220034353-1A_H5GL5DSX5_L2_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_107013_EKDN220034352-1A_H5GL5DSX5_L2_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_21084_EKDN220034367-1A_H5J5MDSX5_L2_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_23279_merged_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_23416_merged_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_320005_merged_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_32023_merged_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_49031_EKDN220034370-1A_H5J5MDSX5_L2_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/E_519005_EKDN220034364-1A_H5J5MDSX5_L2_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/s2010_EDSW210003772-1a_H3WNKDSX2_L3_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/s2012_EDSW210003766-1a_H3WHWDSX2_L4_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/s479013_EDSW210003775-1a_H3WNKDSX2_L3_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/s49097_EDSW210003769-1a_H3FTWDSX2_L2_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/s49124_EDSW210003767-1a_H3WHWDSX2_L2_LipY764_sorted_rmdp_MQ20.bam /shared5/Alex/Y-chromosome_project/BAMs_LipY764/S276023_EDSW200015547-1a2a_HJFWVDSXY_L3L4_LipY764_sorted_rmdp_MQ20.bam  > /shared5/Alex/Y-chromosome_project/variants_bcftools_onlyexmoors/Exmoors_2022_Y_newref.pileup

#created ploidy file. I get the number from the third column by adding up the bases of the contigsfrom the fai index file
chromosomeY 1 6462487 M 1 > ploidy_newref.txt

#bcftools call:
bcftools call -c -O v --ploidy-file /shared5/Alex/Y-chromosome_project/variants_bcftools_onlyexmoors/ploidy_newref.txt -o /shared5/Alex/Y-chromosome_project/variants_bcftools_onlyexmoors/Exmoors_2022_Y_newref.vcf /shared5/Alex/Y-chromosome_project/variants_bcftools_onlyexmoors/Exmoors_2022_Y_newref.pileup
```


# 10. VCF filtering
1.1 First step, quality filter and removing INDELS:
```bash
vcftools --vcf /shared5/Alex/Y-chromosome_project/variants_bcftools_onlyexmoors/Exmoors_2022_Y_newref.vcf --remove-indels --minQ 30 --minGQ 30 --out /shared5/Alex/Y-chromosome_project/variants_bcftools_onlyexmoors/Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30 --recode
```

1.2 To see how much missing information there is per individual:
`vcftools –gzvcf <data.vcf.gz> –missing-indv –out <data>`

2. Minimum allele count filter
vcftools --vcf Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30.recode.vcf --min-alleles 2 --mac 3 --out Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30_biallelic_mac3 --recode

3. Exclude sites where >10% of data is missing
vcftools --vcf Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30_biallelic_mac3.recode.vcf --max-missing 0.9

#To see missingness of individuals:
```bash
vcftools --vcf Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30_biallelic_mac3.recode.vcf --missing-indv --out Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30_biallelic_mac3.recode.vcf
```

Check: https://evomics.org/learning/population-and-speciation-genomics/2020-population-and-speciation-genomics/first-steps-in-genomic-data-analysis/#ex3 for more info on vcf filtering
or https://speciationgenomics.github.io/filtering_vcfs/

OTHER VCFTOOLS CODE:
vcftools --gzvcf $SUBSET_VCF --depth --out $OUT #calculates depth per indivudal
vcftools --gzvcf $SUBSET_VCF --site-mean-depth --out $OUT #calculates depth per site
vcftools --gzvcf $SUBSET_VCF --missing-site --out $OUT #missing data per site
vcftools --gzvcf $SUBSET_VCF --missing-ind --out $OUT #missing data per individual

SET MAF AT 0.01?

# 11. PLINK
## 11.1 Download PLINK
```bash
conda create -n plink -c bioconda plink
source activate plink
```

## 11.2 Create PLINK file
Create plink file:
```bash
vcftools --vcf Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30_biallelic_mac3.recode.vcf --plink --out Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30_biallelic_mac3
```

## 11.3 PCA
```bash
plink --file Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30_biallelic_mac3 --pca --out Exmoors_2022_Y_newref_rmvIndels_minQ30_minGQ30_biallelic_mac3
```
