import sys,os

#1.download data (fastq.gz)
dataname=sys.argv[1]

os.system("/BiO/Install/sratoolkit.2.9.6-ubuntu64/bin/fastq-dump --gzip --split-3 {0}".format(dataname))

#2.unzip data (fastq.gz -> fastq)
os.system("gunzip {0}_1.fastq.gz".format(dataname))
os.system("gunzip {0}_2.fastq.gz".format(dataname))
#3.trimming data (fastq -> trim.fq)
os.system("java -jar /BiO/Install/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 -phred33 {0}_1.fastq {0}_2.fastq {0}.r1.trim.fq {0}.r1.unpair.fq {0}.r2.trim.fq {0}.r2.unpair.fq ILLUMINACLIP:/BiO/Install/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:151:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36".format(dataname))
#4.mapping  (trim.fq -> sam)
os.system("bwa mem -t 4 -R '@RG\\tPL:illumina\\tID:YUHL\\tSM:{0}\\tLB:Hiseq' /BiO/Education/WGS/REF/hg19.fa {0}.r1.trim.fq {0}.r2.trim.fq > {0}.sam".format(dataname))
#5.sam sorting -> bam
os.system("java -jar /BiO/Install/picard-tools-2.22.3/picard.jar AddOrReplaceReadGroups TMP_DIR=TEMP_PICARD VALIDATION_STRINGENCY=LENIENT SO=coordinate I={0}.sam O={0}_sorted.bam RGID={0} RGLB=Hiseq RGPL=Illumina RGPU=unit1 RGSM={0} CREATE_INDEX=true".format(dataname))

#6 remove duplicate (bam -> sam,metrics)
os.system("java -XX:ParallelGCThreads=4 -jar /BiO/Install/picard-tools-2.22.3/picard.jar MarkDuplicates TMP_DIR=TEMP_PICARD VALIDATION_STRINGENCY=LENIENT I={0}_sorted.bam O={0}_dedup.sam M={0}.duplicate_metrics REMOVE_DUPLICATES=true AS=true".format(dataname))
#7 sam->bam
os.system("java -XX:ParallelGCThreads=4 -jar /BiO/Install/picard-tools-2.22.3/picard.jar SortSam TMP_DIR=TEMP_PICARD VALIDATION_STRINGENCY=LENIENT SO=coordinate I={0}_dedup.sam O={0}_dedup.bam CREATE_INDEX=true".format(dataname))

#8~9 BQSR 1st pass
os.system("java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar BaseRecalibrator -R /BiO/Education/WGS/REF/hg19.fa -I {0}_dedup.bam --known-sites /BiO/Education/WGS/REF/dbsnp_138.hg19.vcf --known-sites /BiO/Education/WGS/REF/1000GENOMES-phase_3_indel.vcf -O {0}_recal_pass1.table".format(dataname))

os.system("java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar ApplyBQSR -R /BiO/Education/WGS/REF/hg19.fa -I {0}_dedup.bam --bqsr-recal-file {0}_recal_pass1.table -O {0}_recal_pass1.bam".format(dataname))

#10~11 BQSR 2nd pass
os.system("java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar BaseRecalibrator -R /BiO/Education/WGS/REF/hg19.fa -I {0}_recal_pass1.bam --known-sites /BiO/Education/WGS/REF/dbsnp_138.hg19.vcf --known-sites /BiO/Education/WGS/REF/1000GENOMES-phase_3_indel.vcf -O {0}_recal_pass2.table".format(dataname))

os.system("java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar ApplyBQSR -R /BiO/Education/WGS/REF/hg19.fa -I {0}_recal_pass1.bam -bqsr {0}_recal_pass2.table -O {0}_recal_pass2.bam".format(dataname))

#12 make vcf(snps+indels)
os.system("java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar HaplotypeCaller -R /BiO/Education/WGS/REF/hg19.fa -I {0}_recal_pass2.bam -O {0}.rawVariants.g.vcf -ERC GVCF --standard-min-confidence-threshold-for-calling 20".format(dataname))

#13 phasing
os.system("java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar GenotypeGVCFs -R /BiO/Education/WGS/REF/hg19.fa -V {0}.rawVarients.g.vcf -O {0}_genotype.vcf".format(dataname))

#14 extract snps,indels
os.system("java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar SelectVariants -R /BiO/Education/WGS/REF/hg19.fa -V {0}_genotype.vcf --select-type-to-include SNP -O {0}.rawSNPs.vcf".format(dataname))

os.system("java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar SelectVariants -R /BiO/Education/WGS/REF/hg19.fa -V {0}_genotype.vcf --select-type-to-include INDEL -O {0}.rawINDELs.vcf".format(dataname))

#15 filtered snps,indels
os.system("java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar VariantFiltration -R /BiO/Education/WGS/REF/hg19.fa -V {0}.rawSNPs.vcf -O {0}.rawSNPs.filtered.vcf --filter-name '.' --filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || HaplotypeScore > 13.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0'".format(dataname))

os.system("java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar VariantFiltration -R /BiO/Education/WGS/REF/hg19.fa -V {0}.rawINDELs.vcf -O {0}.rawINDELs.filtered.vcf --filter-name '.' --filter-expression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -8.0'".format(dataname))

#16 snps+indels = variant
os.system("java -Xmx8g -jar /BiO/Install/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar SortVcf -I {0}.rawSNPs.filtered.vcf -I {0}.rawINDELs.filtered.vcf -O {0}.filtered.variant.vcf".format(dataname))

#17 extract pass 
os.system("egrep "^#|PASS" {0}.filtered.variant.vcf > {0}.filtered.variant.pass.vcf".format(dataname))

#18 add annotation (annovar)
os.system("perl /BiO/Install/annovar/table_annovar.pl {0}.filtered.variant.pass.vcf /BiO/Education/WGS/humandb/ -buildver hg19 -out {0} -remove -protocol refGene,cytoBand,avsnp138,clinvar_20190305 -operation g,r,f,f -nastring . -vcfinput".format(dataname))



