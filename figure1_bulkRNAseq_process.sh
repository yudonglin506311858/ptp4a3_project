mkdir fastqc_before
#质量控制
fastqc -t 50 *.gz -o ./fastqc_before
multiqc ./




mkdir fastqc_after
#质量控制
fastqc -t 30 -o ./fastqc_after/ *.gz


mkdir clean
mkdir hisat
mkdir count


for i in {mm_basophilic_1,mm_orthochromatic_1,mm_polychromatic_1,mm_proerythroblast_1,mm_basophilic_2,mm_orthochromatic_2,mm_polychromatic_2,mm_proerythroblast_2,mm_basophilic_3,mm_orthochromatic_3,mm_polychromatic_3,mm_proerythroblast_3};
do 
cutadapt --minimum-length 15 --max-n 8 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o ./clean/${i}_rmadp.fastq.gz  ${i}.fastq.gz >>filter.txt 2>&1
java -jar /data/yudonglin/software/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 30 -phred33 ./clean/${i}_rmadp.fastq.gz ./clean/${i}_fliter.fq.gz AVGQUAL:20 SLIDINGWINDOW:4:15 MINLEN:15 1>>filter.txt 2>&1;
hisat2 --threads 35 -x /data/yudonglin/reference/mm10/genome -U ./clean/${i}_fliter.fq.gz -S ./hisat/${i}.sam;
samtools view -S ./hisat/${i}.sam -b > ./hisat/${i}.bam;
samtools sort ./hisat/${i}.bam -o ./hisat/${i}_sorted.bam;
samtools index ./hisat/${i}_sorted.bam;
rm ./hisat/${i}.sam;
featureCounts -T 30 -t exon -g gene_id -a /data/yudonglin/reference/Mus_musculus.GRCm38.102.gtf -o ./count/${i}.count ./hisat/${i}_sorted.bam >>~/count.txt 2>&1
rm ./hisat/${i}.bam;
conda activate pyscenic
#将bam文件转换为bw文件
bamCoverage --bam ./hisat/${i}_sorted.bam -o ./hisat/${i}_sorted.bam.bw  --binSize 10 
conda deactivate
done




for i in {hs_early_basophilic_1,hs_late_basophilic_2,hs_orthochromatic_3,hs_proerythroblast_1,hs_early_basophilic_2,hs_late_basophilic_3,hs_polychromatic_1,hs_proerythroblast_2,hs_early_basophilic_3,hs_orthochromatic_1,hs_polychromatic_2,hs_proerythroblast_3,hs_late_basophilic_1,hs_orthochromatic_2,hs_polychromatic_3};
do 
cutadapt --minimum-length 15 --max-n 8 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o ./clean/${i}_rmadp.fastq.gz  ${i}.fastq.gz >>filter.txt 2>&1
java -jar /data/yudonglin/software/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 30 -phred33 ./clean/${i}_rmadp.fastq.gz ./clean/${i}_fliter.fq.gz AVGQUAL:20 SLIDINGWINDOW:4:15 MINLEN:15 1>>filter.txt 2>&1;
hisat2 --threads 35 -x /data/yudonglin/reference/hg19/genome -U ./clean/${i}_fliter.fq.gz -S ./hisat/${i}.sam;
samtools view -S ./hisat/${i}.sam -b > ./hisat/${i}.bam;
samtools sort ./hisat/${i}.bam -o ./hisat/${i}_sorted.bam;
samtools index ./hisat/${i}_sorted.bam;
rm ./hisat/${i}.sam;
featureCounts -T 30 -t exon -g gene_id -a /data/yudonglin/reference/Homo_sapiens.GRCh37.75.gtf -o ./count/${i}.count ./hisat/${i}_sorted.bam >>~/count.txt 2>&1
rm ./hisat/${i}.bam;
conda activate pyscenic
#将bam文件转换为bw文件
bamCoverage --bam ./hisat/${i}_sorted.bam -o ./hisat/${i}_sorted.bam.bw  --binSize 10 
conda deactivate
done
