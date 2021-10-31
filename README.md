# *Plantago ovata* scaffolding

Hi-C data was used to link contigs into a chromosome level assembly.

Main stages in scaffolding:

1. Genomic DNA extraction, library preparation, and Illumina sequencing
2. Installing bioinformatics tools
3. Quality control sequencing using FASTQC
4. Quality control library preparation using Phase Genomic pipeline
5. Scaffolding using SALSA2
6. Scaffolding using JUICER, 3D-DNA, and JBAT


**Step 1. Genomic DNA extraction and library preparation**

Description of this step can be found in this publication (link)

**Step 2. Installing software**

I used conda to install most of the bioinformatics tools and created several environments due to software incompatibility.
An example of how to create a conda environment:
```
conda create -n phase_genomics
conda activate phase_genomics
conda install -c bioconda samblaster=0.1.26
```
List of main tools for scaffolding:
- tabix 0.2.6 : https://github.com/samtools/tabix
- fastqc v0.11.9 : https://github.com/s-andrews/FastQC
- trimmomatic v0.39 : http://www.usadellab.org/cms/?page=trimmomatic
- bwa v0.7.17 : https://github.com/lh3/bwa
- samblaster v0.1.26 : https://github.com/GregoryFaust/samblaster
- samtools v1.9 : http://www.htslib.org/doc/1.9/samtools.html
- matlock v20181227 : https://github.com/phasegenomics/matlock 
- salsa2 : https://github.com/marbl/SALSA
- 3D-DNA : https://github.com/aidenlab/3d-dna
- Juicer : https://github.com/aidenlab/juicer
- Juicebox Assembly Tools (JBAT) : https://github.com/aidenlab/Juicebox

More tools can be found in Supplementary File Table 6 (link)

**Step 3. Quality control sequencing process**

Raw data has been deposited in SRA NCBI : SRR14643406

Compressing files
```
bgzip -c -l 9 Povata-HiC_combined_R1.fastq > Povata-HiC_combined_R1.fastq.gz
bgzip -c -l 9 Povata-HiC_combined_R2.fastq > Povata-HiC_combined_R2.fastq.gz
```

Running snakemake pipeline for quality checking and trimming

The pipeline is written in Snakemake workflow management. Information about rules in snakemake can be obtained on this website https://snakemake.readthedocs.io/en/stable/. Snakefile templates and other files required for running snakemake in High Performance Computer (HPC) can be downloaded from https://github.com/UofABioinformaticsHub/snakemake_template.

All tools or modules to execute this Snakefile are listed in envs/HiC.yaml

```
snakemake --profile profiles/slurm --use-singularity --use-conda --snakefile Snakefile_HiC --dry-run
snakemake --profile profiles/slurm --use-singularity --use-conda --snakefile Snakefile_HiC
```


**Step 4. Quality control library**

source: https://phasegenomics.github.io/2019/09/19/hic-alignment-and-qc.html
  ```
### Indexing genome sequences

bwa index -a bwtsw -p Plantago Plantago.fasta

### Aligning Hi-C read to the genome

bwa mem -5SP -t 20 Plantago.fasta Povata-HiC_combined_R1.fastq.gz \
Povata-HiC_combined_R2.fastq.gz | samblaster | samtools view -@ 20 -S -h -b -F 2316 > Povata-HiC_combined.bam


### QC before read filtering

./hic_qc.py -b Povata-HiC_combined.bam -r --sample_type genome

### Read filtering using matlock

matlock bamfilt -i Povata-HiC_combined.bam -o Povata-HiC_filter.bam

### QC after read filtering

./hic_qc.py -b Povata-HiC_filter.bam -r --sample_type genome
 ```

**Step 5: Scaffolding using SALSA2**

source: https://github.com/marbl/SALSA
 ```
 #'Sau3AI' : ['GATC', 'CTAG']
SALSA='run_pipeline.py'
PYTHON='python2.7'
BAMtoBED='bamToBed'
SAMTOOLS='samtools'
ENZYME='GATC'
CONTIGS='Plantago.fasta'
CONTIGS_LENGTH='Plantago.fasta.fai'
CONTIGS_BAM='Povata-HiC_filter.bam'
CONTIGS_BED='Povata-HiC_filter.bed'
CONTIGS_SORTED_BED='Povata-HiC_filter_sorted.bed'
CONTIGS_SCAFFOLDS='scaffold'
TMP='tmp'

### Creating and sorting bed file

$BAMtoBED -i $CONTIGS_BAM  > $CONTIGS_BED
sort -k 4 -T $TMP $CONTIGS_BED > $CONTIGS_SORTED_BED

### Getting contig length

$SAMTOOLS faidx $CONTIGS > $CONTIGS_LENGTH

### Scaffolding
$PYTHON $SALSA -a $CONTIGS -l $CONTIGS_LENGTH -b $CONTIGS_SORTED_BED -e $ENZYME -o $CONTIGS_SCAFFOLDS -m yes

```

**Step 6. Scaffolding using JUICER, 3D-DNA, and JBAT**

```
module load arch/haswell
module load Java/13.0.2
module load LASTZ/1.04.00-foss-2016b
module load Python/3.6.1-foss-2016b
module load numpy/1.12.1-foss-2016b-Python-3.6.1
module load parallel/20180922-foss-2016b
module load SAMtools/1.10-foss-2016b
module load CUDA/7.5.18
./scripts/juicer.sh -g "Plantago" -s "Sau3AI" -z Plantago.fasta \
-y restriction_sites/Plantago_Sau3AI.txt \
-p references/Plantago.fasta.fai -D juicer/Plantago \
-d juicer/Plantago -e -t 20
```

source: https://github.com/aidenlab/3d-dna
```
./run-asm-pipeline.sh -i 1000 -r 10 Plantago.fasta merged_nodups.txt

run this after JBAT-ing

./run-asm-pipeline-post-review.sh -r Plantago.final.review.assembly Plantago.fasta merged_nodups.txt
```

