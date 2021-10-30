# Plantago ovata scaffolding
Main stages for scaffolding:
1. Genomic DNA extraction and library preparation
2. Illumina sequencing
3. Quality control for sequencing with FASTQC
4. Quality control for library preparation with Phase Genomic pipeline
5. Alignment and generating Hi-C contact map
6. Automatic scaffolding
7. Manual correction


```
bgzip -c -l 9 Povata-HiC_combined_R1.fastq > Povata-HiC_combined_R1.fastq.gz
bgzip -c -l 9 Povata-HiC_combined_R2.fastq > Povata-HiC_combined_R2.fastq.gz
```
3. Quality control for sequencing with FASTQC
```
SAMPLES_HiC     = ["Povata-HiC_combined"]
READS           = ["R1", "R2"]
ENVS            = "envs/HiC.yaml"
LOGS            = "logs"
RAW             = "assembly/plantago_genome_sequences/HiC/raw"
QC              = "assembly/plantago_genome_sequences/HiC/qc"
MULTIQC         = "assembly/plantago_genome_sequences/HiC/multiqc"
TRIMMED         = "assembly/plantago_genome_sequences/HiC/trimmed"
QC_TRIMMED      = "assembly/plantago_genome_sequences/HiC/qc_trimmed"
MULTIQC_TRIMMED = "assembly/plantago_genome_sequences/HiC/multiqc_trimmed"
ADAPTERS        = "all_adapters.fa"
REFERENCES      = "references/Plantago.fasta.gz"

################
# Pseudo-rules #
################
rule all:
        input:
                expand(QC+"/{SAMPLE}_{READ}_fastqc.html", SAMPLE=SAMPLES_HiC, READ=READS),
                MULTIQC+"/multiqc_report.html",
                expand(TRIMMED+"/{SAMPLE}_{READ}.fastq.gz", SAMPLE=SAMPLES_HiC, READ=READS),
                expand(QC_TRIMMED+"/{SAMPLE}_{READ}_fastqc.html", SAMPLE=SAMPLES_HiC, READ=READS),
                MULTIQC_TRIMMED+"/multiqc_report.html",


################
# Rules Proper #
################
rule fastqc_raw:
        input:
                expand(RAW+"/{SAMPLE}_{READ}.fastq.gz", SAMPLE=SAMPLES_HiC, READ=READS)
        output:
                zip = QC+"/{SAMPLE}_{READ}_fastqc.zip",
                html = QC+"/{SAMPLE}_{READ}_fastqc.html"
        conda:
                ENVS
        log:
                LOGS+"/fastqc_raw/{SAMPLE}_{READ}"
        params:
                QC
        shell:
                """
                fastqc -o {params} --threads 2 {input}
                2> {log}
                """

rule multiqc_raw:
        input:
                expand(QC+"/{SAMPLE}_{READ}_fastqc.html", SAMPLE=SAMPLES_HiC, READ=READS)
        output:
                MULTIQC+"/multiqc_report.html"
        conda:
                ENVS
        log:
                LOGS+"/multiqc_raw/multiqc_report"
        params:
                indir = QC ,
                outdir = MULTIQC
        shell:
                """
                multiqc {params.indir} -o {params.outdir}
                2> {log}
                """

rule trimmomatic:
        input:
                r1 = RAW+"/{SAMPLES}_R1.fastq.gz",
                r2 = RAW+"/{SAMPLES}_R2.fastq.gz",
                adapters = ADAPTERS
        output:
                r1 = TRIMMED+"/{SAMPLES}_R1.fastq.gz",
                r1_unpaired = TRIMMED+"/{SAMPLES}_R1.unpaired.fastq.gz",
                r2 = TRIMMED+"/{SAMPLES}_R2.fastq.gz",
                r2_unpaired = TRIMMED+"/{SAMPLES}_R2.unpaired.fastq.gz"
        conda:
                ENVS
        log:
                r1 = LOGS+"/trimmomatic/{SAMPLES}_R1",
                r2 = LOGS+"/trimmomatic/{SAMPLES}_R1.unpaired",
                r1_unpaired = LOGS+"/trimmomatic/{SAMPLES}_R2",
                r2_unpaired = LOGS+"/trimmomatic/{SAMPLES}_R2.unpaired"
        shell:
                """
                trimmomatic PE -threads 2 {input.r1} {input.r2} {output.r1} \
                {output.r1_unpaired} {output.r2} {output.r2_unpaired} \
                ILLUMINACLIP:{input.adapters}:2:30:10:3:true LEADING:2 TRAILING:2 SLIDINGWINDOW:4:15 MINLEN:36
                """

rule fastqc_trimmed:
        input:
                expand(TRIMMED+"/{SAMPLE}_{READ}.fastq.gz", SAMPLE=SAMPLES_HiC, READ=READS)
        output:
                zip = QC_TRIMMED+"/{SAMPLE}_{READ}_fastqc.zip",
                html = QC_TRIMMED+"/{SAMPLE}_{READ}_fastqc.html"
        conda:
                ENVS
        log:
                LOGS+"/fastqc_TRIMMED/{SAMPLE}_{READ}"
        params:
                QC_TRIMMED
        shell:
                """
                fastqc -o {params} --threads 2 {input}
                2> {log}
                """

rule multiqc_trimmed:
        input:
                expand(QC_TRIMMED+"/{SAMPLE}_{READ}_fastqc.html", SAMPLE=SAMPLES_HiC, READ=READS)
        output:
                MULTIQC_TRIMMED+"/multiqc_report.html"
        conda:
                ENVS
        log:   
                LOGS+"/multiqc_trimmed/multiqc_report"
        params:
                indir = QC_TRIMMED ,
                outdir = MULTIQC_TRIMMED
        shell:
                """
                multiqc {params.indir} -o {params.outdir}
                2> {log}
                """
  ```
  
  ```
bwa index -a bwtsw -p clip_try1.fasta clip_try1.fasta

bwa mem -5SP -t 20 clip_try1.fasta \
/hpcfs/users/a1697274/bioinformatics/assembly/plantago_genome_sequences/HiC/trimmed/Povata-HiC_combined_R1.fastq.gz \
/hpcfs/users/a1697274/bioinformatics/assembly/plantago_genome_sequences/HiC/trimmed/Povata-HiC_combined_R2.fastq.gz \
| samblaster | samtools view -@ 20 -S -h -b -F 2316 > Povata-HiC_combined.bam

./hic_qc.py -b Povata-HiC_combined.bam -r --sample_type genome

matlock bamfilt -i Povata-HiC_combined.bam -o Povata-HiC_filter.bam

./hic_qc.py -b Povata-HiC_filter.bam -r --sample_type genome

matlock bam2 juicer Povata-HiC_filter.bam Povata-HiC_filter_matlock
 ```
 salsa.sh
 ```
 #'Sau3AI' : ['GATC', 'CTAG']
#https://github.com/marbl/SALSA

SALSA='/hpcfs/users/a1697274/genome/salsa/SALSA/run_pipeline.py'
PYTHON='/hpcfs/users/a1697274/.conda/envs/salsa2/bin/python2.7'
BAMtoBED='/hpcfs/users/a1697274/.conda/envs/salsa2/bin/bamToBed'
SAMTOOLS='/hpcfs/users/a1697274/.conda/envs/salsa2/bin/samtools'
ENZYME='GATC'
CONTIGS='/hpcfs/users/a1697274/genome/phase_genomic/hic_qc/old_files/clip_try1.fasta'
CONTIGS_LENGTH='/hpcfs/users/a1697274/genome/phase_genomic/hic_qc/old_files/clip_try1.fasta.fai'
CONTIGS_1='/hpcfs/users/a1697274/genome/phase_genomic/hic_qc/old_files/salsa_clip_try1.fasta'
CONTIGS_LENGTH_1='/hpcfs/users/a1697274/genome/phase_genomic/hic_qc/old_files/salsa_clip_try1.fasta.fai'
CONTIGS_BAM='/hpcfs/users/a1697274/genome/phase_genomic/hic_qc/old_files/Povata-HiC_filter.bam'
CONTIGS_BED='/hpcfs/users/a1697274/genome/phase_genomic/hic_qc/old_files/Povata-HiC_filter.bed'
CONTIGS_SORTED_BED='/hpcfs/users/a1697274/genome/phase_genomic/hic_qc/old_files/Povata-HiC_filter_sorted.bed'
CONTIGS_SORTED_BED_1='/hpcfs/users/a1697274/genome/phase_genomic/hic_qc/old_files/salsa_Povata-HiC_filter_sorted.bed'
CONTIGS_SCAFFOLDS='/hpcfs/users/a1697274/genome/salsa/old_file/scaffold'
TMP='/hpcfs/users/a1697274/genome/salsa/old_file/tmp'

### Creating and sorting bed file

$BAMtoBED -i $CONTIGS_BAM  > $CONTIGS_BED
sort -k 4 -T $TMP $CONTIGS_BED > $CONTIGS_SORTED_BED

### Getting contig length

$SAMTOOLS faidx $CONTIGS > $CONTIGS_LENGTH

#Because salsa cannot take input file with ":" in the header, we need to modify the header
sed '/^>/s/|arrow.*//' $CONTIGS > $CONTIGS_1
samtools faidx $CONTIGS_1 -o $CONTIGS_LENGTH_1
awk -v "OFS=\t" '{$1=$1;sub(/\|arrow.*/, "", $1); print}' $CONTIGS_SORTED_BED > $CONTIGS_SORTED_BED_1

### Option 2
#$PYTHON $SALSA -a $CONTIGS -l $CONTIGS_LENGTH -b $CONTIGS_SORTED_BED -e $ENZYME -o $CONTIGS_SCAFFOLDS -m yes

$PYTHON $SALSA -a $CONTIGS_1 -l $CONTIGS_LENGTH_1 -b $CONTIGS_SORTED_BED_1 -e $ENZYME -o $CONTIGS_SCAFFOLDS -m yes -i 20 -c 31

```

juicer.sh
```
module load arch/haswell
module load Java/13.0.2
module load LASTZ/1.04.00-foss-2016b
module load Python/3.6.1-foss-2016b
module load numpy/1.12.1-foss-2016b-Python-3.6.1
module load parallel/20180922-foss-2016b
module load SAMtools/1.10-foss-2016b
module load CUDA/7.5.18
./scripts/juicer.sh -g "salsa_clip_try1" -s "Sau3AI" -z references/salsa_clip_try1.fasta \
-y restriction_sites/salsa_clip_try1_Sau3AI.txt \
-p references/salsa_clip_try1.fasta.fai -D /hpcfs/users/a1697274/genome/juicer/Plantago \
-d /hpcfs/users/a1697274/genome/juicer/Plantago -e -t 20
```

3d_dna.sh
```
./run-asm-pipeline.sh -i 1000 -r 10 salsa_clip_try1.fasta merged_nodups.txt

run this after JBAT-ing

./run-asm-pipeline-post-review.sh -r salsa_clip_try1.final.review.assembly salsa_clip_try1.fasta merged_nodups.txt
```

