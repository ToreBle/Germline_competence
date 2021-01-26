#!/bin/bash -l

#Alignment with bwa mem
bwa mem "/data/genomes/mm10.fa"  $ID".fastq.gz" > "bwa_"$ID".sam"
samtools view -b "bwa_"$ID".sam" > "bwa_"$ID".bam"
rm "bwa_"$ID".sam"
samtools sort "bwa_"$ID".bam" -o "bwa_"$ID"_complete.bam"
rm "bwa_"$ID".bam"
samtools index "bwa_"$ID"_complete.bam" "bwa_"$ID"_complete.bai"
samtools flagstat "bwa_"$ID"_complete.bam" > "bwa_"$ID"_flagstat_complete.txt"


#Remove duplicates with Picard tools
java -jar picard.jar MarkDuplicates INPUT="bwa_"$ID"_complete.bam" OUTPUT="bwa_"$ID"_deduplicat.bam" METRICS_FILE="picard_"$ID".txt" REMOVE_DUPLICATES=true
#rm "bwa_"$ID"_sorted.bam"
samtools sort "bwa_"$ID"_deduplicat.bam" -o "bwa_"$ID"_sorted.bam"
samtools index "bwa_"$ID"_sorted.bam" "bwa_"$ID"_sorted.bai"
samtools flagstat "bwa_"$ID"_sorted.bam" > "bwa_"$ID"_flagstat.txt"

#Convert in bigwig file with deeptools
bamCoverage --binSize 5 --normalizeUsing RPGC --effectiveGenomeSize 2150570000 --ignoreForNormalization chrX chrY -bl "mm10.blacklist_sorted.bed" -b "bwa_"$ID"_sorted.bam" -o "bwa_"$ID".bigWig"

#Peak calling with MAC2 (optional with specified INPUT.bam)
#macs2 callpeak -t "bwa_"$ID"_sorted.bam" -c INPUT.bam -f BAM -n "/scratch2/tbleckwe/macs/"$ID --broad -g mm -q 0.1 -m 5 50 --fix-bimodal --extsize 200

