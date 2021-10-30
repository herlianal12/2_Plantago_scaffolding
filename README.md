# Plantago ovata scaffolding

Main stages in scaffolding:

1. Genomic DNA extraction, library preparation, and Illumina sequencing
2. Quality control sequencing using FASTQC
3. Quality control library preparation using Phase Genomic pipeline
4. Alignment and generating Hi-C contact map
5. Automatic scaffolding
6. Manual curation


**Step 1. Genomic DNA extraction and library preparation**

Description of this step can be found in this publication (link)

**Step 2. Quality control sequencing process***

Compressing files
```
bgzip -c -l 9 Povata-HiC_combined_R1.fastq > Povata-HiC_combined_R1.fastq.gz
bgzip -c -l 9 Povata-HiC_combined_R2.fastq > Povata-HiC_combined_R2.fastq.gz
```

Running snakemake pipeline for quality checking and trimming

The pipeline is written in Snakemake workflow management. Information about rules in snakemake can be obtained in this website https://snakemake.readthedocs.io/en/stable/. Snakefile templates and other files required for running snakemake in High Performance Computer (HPC) can be downloaded from https://github.com/UofABioinformaticsHub/snakemake_template.

All softwares or modules to execute this Snakefile are listed in envs/HiC.yaml

```
snakemake --profile profiles/slurm --use-singularity --use-conda --snakefile Snakefile_HiC --dry-run
snakemake --profile profiles/slurm --use-singularity --use-conda --snakefile Snakefile_HiC
```


**Step 3. Quality control library 
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

