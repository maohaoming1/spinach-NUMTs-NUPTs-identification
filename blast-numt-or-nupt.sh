#!/bin/bash



blastn -query mito.fasta -db hap1_db -evalue 1e-5 -out mito_hap1.txt -outfmt 6

blastn -query mito.fasta -db hap2_db -evalue 1e-5 -out mito_hap2.txt -outfmt 6

#The blast results were filtered identity > 80% and alignment length > 100  bp
awk '$3 > 80 && $4 > 100' mito_hap1.txt > mito_hap1.filter.txt
awk '$3 > 80 && $4 > 100' mito_hap2.txt > mito_hap2.filter.txt
awk '{print $1"\t"$7"\t"$8}' mito_hap1.filter.txt > mito_hap1.just-mito.txt
awk '{print $1"\t"$7"\t"$8}' mito_hap2.filter.txt > mito_hap2.just-mito.txt 
awk '{print $2"\t"$9"\t"$10}' mito_hap1.filter.txt > mito_hap1.just-nucl.txt
awk '{print $2"\t"$9"\t"$10}' mito_hap2.filter.txt > mito_hap2.just-nucl.txt
sort -k1,1 -k2,2n mito_hap1.just-mito.txt > mito_hap1.just-mito.sort.txt
sort -k1,1 -k2,2n mito_hap2.just-mito.txt > mito_hap2.just-mito.sort.txt
sort -k1,1 -k2,2n mito_hap2.just-nucl.txt > mito_hap2.just-nucl.sort.txt
sort -k1,1 -k2,2n mito_hap1.just-nucl.txt > mito_hap1.just-nucl.sort.txt
#Remove the overlap areas
python rm-overlap.new.py



#Extract the sequence after removing overlap area

seqkit subseq --bed mito_hap1.just-nucl.merge.txt -o mito_hap1.just-nucl.fa hap1.genome.fasta
seqkit subseq --bed mito_hap2.just-nucl.merge.txt -o mito_hap2.just-nucl.fa hap2.genome.fasta
seqkit subseq --bed mito_hap2.just-mito.merge.txt -o mito_hap2.just-mito.fa mito.fasta
seqkit subseq --bed mito_hap1.just-mito.merge.txt -o mito_hap1.just-mito.fa mito.fasta




##The fa file obtained in the fourth step has chloroplast sequence, so it needs to be compared with the chloroplast genome to remove the chloroplast sequence

blastn -query mito_hap1.just-nucl.fa -db pt_db -evalue 1e-5 -out mito_hap1.just-nucl-vs-pt.txt -outfmt 6
blastn -query mito_hap2.just-nucl.fa -db pt_db -evalue 1e-5 -out mito_hap2.just-nucl-vs-pt.txt -outfmt 6
blastn -query mito_hap2.just-mito.fa -db pt_db -evalue 1e-5 -out mito_hap2.just-mito-vs-pt.txt -outfmt 6
blastn -query mito_hap1.just-mito.fa -db pt_db -evalue 1e-5 -out mito_hap1.just-mito-vs-pt.txt -outfmt 6

awk '$3 > 80 && $4 > 100' mito_hap1.just-nucl-vs-pt.txt > mito_hap1.just-nucl-vs-pt.filter.txt
awk '$3 > 80 && $4 > 100' mito_hap2.just-nucl-vs-pt.txt > mito_hap2.just-nucl-vs-pt.filter.txt
awk '$3 > 80 && $4 > 100' mito_hap1.just-mito-vs-pt.txt > mito_hap1.just-mito-vs-pt.filter.txt
awk '$3 > 80 && $4 > 100' mito_hap2.just-mito-vs-pt.txt > mito_hap2.just-mito-vs-pt.filter.txt
awk '{print $1}' mito_hap1.just-nucl-vs-pt.filter.txt >  mito_hap1.just-nucl-vs-pt.remove.list
awk '{print $1}' mito_hap2.just-nucl-vs-pt.filter.txt >  mito_hap2.just-nucl-vs-pt.remove.list
awk '{print $1}' mito_hap1.just-mito-vs-pt.filter.txt > mito_hap1.just-mito-vs-pt.remove.list
awk '{print $1}' mito_hap2.just-mito-vs-pt.filter.txt > mito_hap2.just-mito-vs-pt.remove.list


seqkit grep -v -f mito_hap1.just-nucl-vs-pt.remove.list  mito_hap1.just-nucl.fa > mito_hap1.just-nucl-remove-pt.fa
seqkit grep -v -f mito_hap2.just-nucl-vs-pt.remove.list  mito_hap2.just-nucl.fa > mito_hap2.just-nucl-remove-pt.fa
seqkit grep -v -f mito_hap1.just-mito-vs-pt.remove.list mito_hap1.just-mito.fa > mito_hap1.just-mito-remove-pt.fa
seqkit grep -v -f mito_hap2.just-mito-vs-pt.remove.list mito_hap2.just-mito.fa > mito_hap2.just-mito-remove-pt.fa

seqkit rmdup -n -i mito_hap1.just-mito-remove-pt.fa -o mito_hap1.just-mito-remove-pt_dedup.fa 
seqkit rmdup -n -i mito_hap2.just-mito-remove-pt.fa -o mito_hap2.just-mito-remove-pt_dedup.fa 

seqkit rmdup -n -i mito_hap1.just-nucl-remove-pt.fa -o mito_hap1.just-nucl-remove-pt_dedup.fa
seqkit rmdup -n -i mito_hap2.just-nucl-remove-pt.fa -o mito_hap2.just-nucl-remove-pt_dedup.fa



#The initial NUMT sequence was obtained by removing the corresponding chloroplast sequence

#Calculate the demitoh of the NUMT sequence
minimap2 -ax map-hifi -t 10 mito_hap1.just-mito-remove-pt_dedup.fa hifi_reads.fastq.gz | samtools view -@ 10 -bS - | samtools sort -@ 10 -o mito_hap1.just-mito.sort.bam
samtools index mito_hap1.just-mito.sort.bam
samtools depth mito_hap1.just-mito.sort.bam > mito_hap1.just-mito.depth



minimap2 -ax map-hifi -t 10 mito_hap2.just-mito-remove-pt_dedup.fa hifi_reads.fastq.gz |  samtools view -@ 10 -bS - | samtools sort -@ 10 -o mito_hap2.just-mito.sort.bam
samtools index mito_hap2.just-mito.sort.bam
samtools depth mito_hap2.just-mito.sort.bam > mito_hap2.just-mito.depth
##the same to nucler fragment


awk '{sum[$1]+=$3; cnt[$1]++} END{for(i in sum) print i, sum[i]/cnt[i]}'  mito_hap1.just-mito.depth > mito_hap1.just-mito.mean.depth
awk '{sum[$1]+=$3; cnt[$1]++} END{for(i in sum) print i, sum[i]/cnt[i]}'  mito_hap2.just-mito.depth > mito_hap2.just-mito.mean.depth
awk '{sum[$1]+=$3; cnt[$1]++} END{for(i in sum) print i, sum[i]/cnt[i]}' mito_hap2.just-nucl.depth > mito_hap2.just-nucl.mean.depth
awk '{sum[$1]+=$3; cnt[$1]++} END{for(i in sum) print i, sum[i]/cnt[i]}' mito_hap1.just-nucl.depth > mito_hap1.just-nucl.mean.depth

#The NUMT fragments were finally determined based on nuclear and mitochondrial demitoh
