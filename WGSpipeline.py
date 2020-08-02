import sys
import os

#1.download data (fastq.gz)
dataname=sys.argv[1]
os.system("/BiO/Install/sratoolkit.2.9.6-ubuntu64/bin/fastq-dump --gzip --split-3 {0}".format(dataname))

#2.unzip data (fastq.gz -> fastq)
os.system("gunzip *.gz")

#3.trimming data (fastq -> trim.fq)
os.system("java -jar /BiO/Install/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 -phred33 {0}_1.fastq {0}_2.fastq {0}.r1.trim.fq {0}.r1.unpair.fq {0}.r2.trim.fq {0}.r2.unpair.fq ILLUMINACLIP:/BiO/Install/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:151:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36".format(dataname))

#4.mapping  (trim.fq -> sam)
os.system("bwa mem -t 4 -R '@RG\tPL:illumina\tID:YUHL\tSM:SRR2069499\tLB:HiSeq' /BiO/Education/WGS/REF/hg19.fa {0}.r1.trim.fq {0}.r2.trim.fq > {0}.sam".format(dataname))

#5.sam sorting -> bam
os.system("java -jar /BiO/Install/picard-tools-2.22.3/picard.jar AddOrReplaceReadGroups TMP_DIR=TEMP_PICARD VALIDATION_STRINGENCY=LENIENT SO=coordinate I={0}.sam O=HSRR1002940_sorted.bam RGID=SRR1002940 RGLB=Hiseq RGPL=Illumina RGPU=unit1 RGSM=SRR1002940 CREATE_INDEX=true")


