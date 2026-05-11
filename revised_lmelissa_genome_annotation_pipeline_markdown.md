---
output:
  word_document: default
  html_document: default
---

# Genome Annotation Pipeline for *Lycaeides melissa*

This document describes a reproducible genome annotation workflow for the chromosome-level genome assembly of the butterfly *Lycaeides melissa*. The pipeline integrates repeat identification, repeat masking, RNA-seq alignment, gene prediction using BRAKER3, and functional annotation using InterProScan.

The workflow is designed as a general reference for annotating high-quality eukaryotic genomes using modern annotation tools and high-performance computing (HPC) environments.

---

# Genome Annotation using BRAKER3

[BRAKER3](https://pubmed.ncbi.nlm.nih.gov/37398387/) is an automated genome annotation pipeline that integrates RNA-seq evidence and protein homology for gene prediction in high-quality genome assemblies.

This workflow includes:

1. *De novo* repeat identification
2. Repeat masking
3. RNA-seq alignment
4. Gene prediction with BRAKER3
5. Functional annotation with InterProScan

---

# Software Requirements

The following software packages were used:

1. [BRAKER3](https://github.com/Gaius-Augustus/BRAKER)
2. [RepeatModeler](http://www.repeatmasker.org/RepeatModeler/)
3. [RepeatMasker](http://www.repeatmasker.org/RepeatMasker/)
4. [RepeatScout](http://bix.ucsd.edu/repeatscout/)
5. [AUGUSTUS](https://bioinf.uni-greifswald.de/augustus/)
6. [BUSCO](https://busco.ezlab.org/)
7. [BEDTools](https://bedtools.readthedocs.io/en/latest/)
8. [STAR](https://github.com/alexdobin/STAR)
9. [InterProScan](https://github.com/ebi-pf-team/interproscan)
10. Perl
11. Python3
12. Singularity/Apptainer

---

# Input Data

## Reference Genome

Chromosome-level *Lycaeides melissa* genome assembly in FASTA format.

## RNA-seq Data

Paired-end RNA-seq libraries representing multiple tissues and developmental stages.

## Protein Databases

Protein FASTA files from closely related Lepidopteran species for evidence-based gene prediction.

---

# 1. Repeat Identification

Repeat identification was performed using both RepeatModeler and RepeatScout to generate *de novo* repeat libraries.

---

## A. RepeatModeler

### Build RepeatModeler Database

```bash
BuildDatabase \
  -name lmelissa \
  -engine ncbi \
  reference_genome.fasta
```

### Run RepeatModeler

```bash
RepeatModeler \
  -database lmelissa \
  -engine ncbi \
  -LTRStruct \
  -pa 16
```

### Example SLURM Script

```bash
#!/bin/bash
#SBATCH --job-name=repeatmodeler
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=180:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL

module load tetools

BuildDatabase -name lmelissa -engine ncbi reference_genome.fasta

RepeatModeler \
  -database lmelissa \
  -engine ncbi \
  -LTRStruct \
  -pa 16
```

---

## B. RepeatScout

RepeatScout was additionally used to identify repetitive elements not captured by RepeatModeler.

### Build Repeat Frequency Table

```bash
build_lmer_table \
  -sequence reference_genome.fasta \
  -freq genome.freq
```

### Generate Repeat Library

```bash
RepeatScout \
  -sequence reference_genome.fasta \
  -output genome.repseq.fa \
  -freq genome.freq
```

---

## RepeatScout Filtering

### Filter Stage 1

```bash
cat genome.repseq.fa | filter-stage-1.prl > genome.repseq.f1.fa
```

### Mask Genome Using Preliminary Library

```bash
RepeatMasker \
  -e ncbi \
  -lib genome.repseq.f1.fa \
  -dir masking_round1 \
  reference_genome.fasta
```

### Filter Stage 2

```bash
cat genome.repseq.f1.fa | \
  filter-stage-2.prl --cat=masking_round1/reference_genome.fasta.out \
  > genome.repseq.f2.fa
```

---

# 2. Repeat Masking

RepeatMasker was used to generate both hard-masked and soft-masked genome assemblies.

- Hard-masked genome: used for alignment tools such as STAR
- Soft-masked genome: used for gene prediction with BRAKER3 and AUGUSTUS

---

## RepeatMasker using RepeatModeler Library

```bash
RepeatMasker \
  -pa 8 \
  -lib consensi.fa.classified \
  -dir repeatmasker_output \
  reference_genome.fasta
```

---

## Generate Hard-Masked Genome

```bash
RepeatMasker \
  -pa 8 \
  -lib genome.repseq.f2.fa \
  -dir hardmask_output \
  reference_genome.fasta
```

---

## Generate Soft-Masked Genome

```bash
RepeatMasker \
  -pa 8 \
  -xsmall \
  -lib genome.repseq.f2.fa \
  -dir softmask_output \
  reference_genome.fasta
```

The `-xsmall` option converts repetitive regions to lowercase bases, producing a soft-masked genome suitable for BRAKER3.

---

# 3. RNA-seq Alignment

RNA-seq reads were aligned to the repeat-masked genome assembly using STAR.

---

## Build STAR Genome Index

```bash
STAR \
  --runMode genomeGenerate \
  --runThreadN 12 \
  --genomeDir lmelissa_starindex \
  --genomeFastaFiles softmasked_genome.fasta
```

---

## Concatenate RNA-seq Reads

```bash
cat *R1*.fq.gz > merged_R1.fastq.gz
cat *R2*.fq.gz > merged_R2.fastq.gz
```

---

## Align RNA-seq Reads

```bash
STAR \
  --twopassMode Basic \
  --runThreadN 12 \
  --genomeDir lmelissa_starindex \
  --readFilesIn merged_R1.fastq.gz merged_R2.fastq.gz \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix lmelissa_
```

The resulting BAM alignment file is used as transcript evidence for BRAKER3.

---

# 4. Gene Prediction with BRAKER3

BRAKER3 was run using Singularity/Apptainer containers together with RNA-seq alignments.

---

## Build BRAKER3 Container

```bash
singularity build braker3.sif docker://teambraker/braker3:latest
```

---

## Configure Environment Variables

```bash
export BRAKER_SIF=/path/to/braker3.sif
export AUGUSTUS_CONFIG_PATH=/path/to/augustus/config
```

---

## Run BRAKER3

```bash
singularity exec -B ${PWD}:${PWD} ${BRAKER_SIF} braker.pl \
  --genome=softmasked_genome.fasta \
  --bam=lmelissa_Aligned.sortedByCoord.out.bam \
  --softmasking \
  --threads=16 \
  --workingdir=braker_lmelissa
```

---

## Example SLURM Submission Script

```bash
#!/bin/bash
#SBATCH --job-name=braker3
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=240:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL

module load singularity

export BRAKER_SIF=/path/to/braker3.sif
export AUGUSTUS_CONFIG_PATH=/path/to/augustus/config

singularity exec -B ${PWD}:${PWD} ${BRAKER_SIF} braker.pl \
  --genome=softmasked_genome.fasta \
  --bam=lmelissa_Aligned.sortedByCoord.out.bam \
  --softmasking \
  --threads=16 \
  --workingdir=braker_lmelissa
```

---

# 5. Functional Annotation with InterProScan

Predicted coding sequences were functionally annotated using InterProScan.

---

## Extract CDS Features from GTF

```bash
awk '$3 == "CDS"' braker.gtf > CDS_features.gtf
```

---

## Convert CDS Coordinates to FASTA

```bash
bedtools getfasta \
  -fi softmasked_genome.fasta \
  -bed CDS_features.gtf \
  -fo CDS_sequences.fasta
```

---

## Split FASTA Files for Parallelization

```bash
awk 'BEGIN {n_seq=0;count=1} /^>/ {
    if(n_seq%10000==0){
        file=sprintf("CDS_%d.fasta",count);
        count++;
    }
    print >> file;
    n_seq++;
    next;
} { print >> file; }' < CDS_sequences.fasta
```

---

## Run InterProScan

```bash
interproscan.sh \
  -t n \
  -i CDS_1.fasta \
  -f tsv \
  -goterms \
  -iprlookup \
  -pathways
```

---

## Parallel InterProScan Execution Script

```perl
use Parallel::ForkManager;

my $max = 16;
my $pm = Parallel::ForkManager->new($max);

FILES:
foreach my $fa (@ARGV){
    $pm->start and next FILES;

    system "interproscan.sh -t n -i $fa -f tsv -goterms -iprlookup -pathways";

    $pm->finish;
}

$pm->wait_all_children;
```

---

# 6. Processing Functional Annotation Output

Custom Python scripts were used to parse and reformat InterProScan output into tabular annotation files for downstream analyses.

Example command:

```bash
python3 extract_annot_info_from_iprfile.py \
  finalCDS_ipr.tsv \
  final_func_annot_info.tsv
```

---

# 7. SNP Annotation Workflow

Structural and functional annotations for SNPs were generated using custom Python scripts.

---

## Structural Annotation

Structural annotations were generated using the final BRAKER3 GTF file.

Example command:

```bash
python3 create_snp_annotations.py \
  --gtf braker.gtf \
  --snp_file snps.txt \
  --output structural_annotations.tsv
```

---

## Functional Annotation

Functional annotations were generated using formatted InterProScan outputs.

Example command:

```bash
python3 snp_func_annot_2.py \
  --snp_file snps.txt \
  --annotation_file final_func_annot_info.tsv \
  --output_file annotated_snps.tsv
```

---

# Recommended Output Files

The final annotation workflow generates the following primary outputs:

| Output | Description |
|---|---|
| `softmasked_genome.fasta` | Soft-masked genome assembly |
| `hardmasked_genome.fasta` | Hard-masked genome assembly |
| `braker.gtf` | Final BRAKER3 gene models |
| `braker.aa` | Predicted protein sequences |
| `braker.codingseq` | Predicted coding sequences |
| `InterProScan.tsv` | Functional annotation output |
| `final_func_annot_info.tsv` | Parsed functional annotation table |
| `annotated_snps.tsv` | Final SNP annotation table |

---

# Notes and Recommendations

1. Always assess annotation completeness using BUSCO.
2. Soft-masked genomes are strongly recommended for BRAKER3.
3. RNA-seq data from multiple tissues and developmental stages substantially improves annotation quality.
4. Parallelization is highly recommended for InterProScan due to runtime requirements.
5. Maintain reproducible software environments using containers or Conda environments.

---

# References

- Hoff KJ et al. 2023. BRAKER3: Fully automated genome annotation.  
- RepeatModeler and RepeatMasker documentation  
- InterProScan documentation  
- STAR aligner documentation  
- BEDTools documentation

