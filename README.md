# Peanut-PanMAGIC-Platform

This repository hosts pipelines for peanut genome assembly and pangenome analysis.

The peanut pangenome was constructed using 18 peanut (*Arachis hypogaea*) lines
representing diverse cultivated germplasm. Graph-based integration of whole-genome
assemblies captures SNPs, structural variants, and presence–absence variation,
enabling comparative analyses, introgression detection, and variant density
profiling across chromosomes.

---

## Genome Assembly Pipeline

This section documents the genome assembly workflow used to generate assemblies
for the 18 parental peanut genomes used in pangenome construction.

---

## Biological Material and Sequencing

For the 18 parental genomes, fresh leaf tissues were collected from each parent
to obtain at least 5 g of fresh tissue. DNA was extracted from the collected tissue
and sequencing libraries were prepared for PacBio Sequel II HiFi sequencing.
Each genotype was sequenced to approximately 20× coverage.

Genome assembly was conducted using **hifiasm**. The assembled genomes were
scaffolded using the reference Tifrunner genome (version 2) with **RagTag**.

---

## Steps Involved

---

### 1. Adapter Filtration

Adapters were filtered using **HiFiAdapterFilt**
(https://github.com/sheinasim/HiFiAdapterFilt).

```bash
ml HiFiAdapterFilt
bash hifiadapterfilt.sh -t 30 -p Prefix -o Outputdirectory
```

---

### 2. Genome Assembly with hifiasm

#### A. Homozygous assemblies

i. A total of 16 out of 18 genomes were assembled using the following pipeline:

```bash
module load Hifiasm
hifiasm -o Name.hifiasm -t 64 -lo Name.filt.fastq.gz
```

The `-lo` flag was specified for homozygous genomes.

Primary contigs were converted to FASTA format using:

```bash
awk '/^S/{print ">"$2;print $3}' C99R.HiFi.bp.p_ctg.gfa > C99R.HiFi.bp.p_ctg.fa
```

ii. Due to tandem duplication in some regions of the genomes GA-12Y and MarcI,
contigs smaller than 100 kb were filtered from the primary assemblies, resulting
in improved dot plot interpretation.

```bash
ml SeqKit
seqkit seq -m 100000 contig.fa > filtered_contig.fa
```

---

#### B. Heterozygous assembly with Hi-C integration (NMValA)

Dot plots generated using the `-lo` flag showed potential tandem duplication or
heterozygosity. Raw reads of New Mexico Valencia (NMValA) were analyzed using
GenomeScope2, which indicated that the genome was heterozygous.

The genome was therefore assembled by phasing into two haplotypes. Hi-C reads
were incorporated, and the coverage peak was adjusted to 27 based on
GenomeScope2 analysis.

```bash
ml hifiasm

hifiasm -o NMValA.HIFIasm -t64 --hom-cov 27 \
  --h1 NMValA_OmniC_R1.fastq.gz \
  --h2 NMValA_OmniC_R2.fastq.gz \
  NMValA_HIFI.filt.fastq.gz
```

---

#### C. Chromosome-specific heterozygosity (CC477)

CC477 showed homozygosity across most chromosomes; however, chromosome 9 exhibited
heterozygosity or unexpected structural variation. Reads corresponding to
chromosome 9 were extracted and assembled as heterozygous (default hifiasm),
while the remaining chromosomes were assembled as homozygous using the `-l0`
flag.

Hap1 and Hap2 assemblies were generated, differing only in chromosome 9.

---

##### Mapping with reference genome

```bash
minimap2 -d $workdir/Tifrunner.mmi $Tifref
minimap2 -ax asm10 Tifrunner.mmi $workdir/CC477_HIFI.filt.fastq.gz -t 64 > CC477_HIFI.sam
samtools view -Sb CC477_HIFI.sam -@ 64 > CC477_HIFI.bam
samtools sort CC477_HIFI.bam -@ 64 -o CC477_HIFI_sorted.bam
samtools index CC477_HIFI_sorted.bam
```

---

##### Extracting chromosome 9

```bash
samtools view -b -h CC477_HIFI_sorted.bam "arahy.Tifrunner.gnm2.chr09" > CC477_HIFI_chr09.bam
samtools sort CC477_HIFI_chr09.bam -@ 64 -o CC477_HIFI_chr09_sorted.bam
samtools index CC477_HIFI_chr09_sorted.bam
```

---

##### Testing the BAM file

```bash
samtools view -c -F 4 CC477_HIFI_chr09.bam
samtools quickcheck CC477_HIFI_chr09.bam
```

---

##### BAM to FASTQ conversion

```bash
samtools view -b CC477_HIFI_chr09.bam | samtools fastq - > CC477_HIFI_chr09.fastq
```

---

##### Genome assembly of chromosome 9

```bash
ml hifiasm
hifiasm -o CC477.chr09.HIFIasm -t64 --hom-cov 26 CC477_HIFI_Chr09.fastq.gz
```

---

##### Extracting all but chromosome 9

```bash
samtools view -H CC477_HIFI_sorted.bam | grep -v 'arahy.Tifrunner.gnm2.chr09' > header.txt
samtools view CC477_HIFI_sorted.bam | grep -v 'arahy.Tifrunner.gnm2.chr09' | \
  cat header.txt - | samtools view -b - > CC477_HIFI_No_Chr09.bam

samtools sort CC477_HIFI_No_Chr09.bam -@ 64 -o CC477_HIFI_No_Chr09.sorted.bam
samtools index CC477_HIFI_No_Chr09.sorted.bam
mv CC477_HIFI_No_Chr09.sorted.bam NO.Chr09.CC477.HIFI.sorted.bam
```

---

##### Testing BAM file

```bash
ml SAMtools
samtools quickcheck NO.Chr09.CC477.HIFI.sorted.bam
samtools index NO.Chr09.CC477.HIFI.sorted.bam
samtools idxstats NO.Chr09.CC477.HIFI.sorted.bam
```

---

##### BAM to FASTQ conversion

```bash
samtools view -b NO.Chr09.CC477.HIFI.sorted.bam | samtools fastq - > NO.Chr09.CC477.HIFI.sorted.fastq
```

---

##### Genome assembly of remaining chromosomes

```bash
ml hifiasm/0.19.6-GCCcore-11.3.0
hifiasm -o NO.Chr09.CC477.HIFIasm -t 64 -l0 NO.Chr09.CC477.HIFI.sorted.fastq
```

---

##### Concatenating haplotypes

```bash
cat CC477.Chr09.HIFIasm.bp.hap1.p_ctg.fa NO.Chr09.CC477.HIFIasm.bp.p_ctg.fa > CC477.complete.Hap1.fasta
cat CC477.Chr09.HIFIasm.bp.hap2.p_ctg.fa NO.Chr09.CC477.HIFIasm.bp.p_ctg.fa > CC477.complete.Hap2.fasta
```

---

## Scaffolding to the Progenitor Reference Genome

The assembled genomes were scaffolded using the progenitor genome (version 1),
where *A. duranensis* and *A. ipaensis* were concatenated into a single reference.

```bash
ml RagTag
ragtag.py scaffold -t 64 -r -o Folder_name WildGenome.fa contigs.fa
```

---

## Genome Alignment and Dot Plot Analysis

Primary contigs or Hap1/Hap2 contigs were plotted against the Tifrunner genome
(version 2). Scaffolded genomes were also compared against Tifrunner V2.

```bash
module load MUMmer
nucmer -t 64 Tifref.fasta ragtag.scaffold.fasta -p mummerplot
```

```bash
ml R
show-coords -c mummerplot.delta > mummerplot.coords
./mummerCoordsDotPlotly.R -i mummerplot.coords -o mummerpicture -s -t -m 15000 -q 200000 -k 20 -l
```

---

## Pangenome Graph Construction

The peanut pangenome was constructed using **Minigraph-Cactus (v2.9.0)**, a
graph-based pangenome framework that integrates multiple genome assemblies into
a single reference graph.

Genome assemblies from the parental lines were provided as input, and **Tifseq**
was used as the reference genome for graph construction. The resulting pangenome
graph captures sequence variation and structural diversity across genomes.

The pipeline generates multiple graph and variant representations to support
downstream analyses, including graph-based read mapping and variant projection.

---

### Pangenome construction command (example)

```bash
cactus-pangenome \
  --workDir cactus_wd \
  ./js \
  ./sample_names \
  --outDir output_directory \
  --outName pangenome_image \
  --reference Tifseq \
  --vcf \
  --giraffe \
  --gfa \
  --gbz \
  --maxCores 80
```

---

## Output Formats

The pangenome construction produced the following outputs:

- **VCF**: projected variants relative to the reference genome  
- **GFA**: graph representation of the pangenome  
- **GBZ**: compressed graph index for efficient storage and traversal  
- **Giraffe indexes**: for graph-based read mapping  

---
------------------------------------------------------------
Pangenome Visualization
------------------------------------------------------------

Pangenome alignments can be visualized using multiple alignment viewers such as
Mauve or Geneious Prime. These tools enable interactive inspection of synteny,
structural variation, rearrangements, and introgressed regions across genomes
included in the pangenome.

Visualization is performed by exporting pangenome alignments into standard
multiple-alignment formats (MAF / XMFA).

------------------------------------------------------------
Visualization workflow
------------------------------------------------------------

```

Pangenome graph / HAL alignment
  ↓
Export multiple alignment (MAF)
  ↓
Convert MAF to XMFA
  ↓
Visualize in Mauve or Geneious Prime

```
------------------------------------------------------------
Exporting alignments from HAL (example)
------------------------------------------------------------

halStats --genomes chr09.hal

hal2maf \
  --targetGenomes Tifseq,TifNV,C431 \
  --maxBlockLen 10000000 \
  chr09.hal \
  chr09.maf

------------------------------------------------------------
Converting MAF to XMFA (Python)
------------------------------------------------------------

XMFA format is supported by both Mauve and Geneious Prime and allows
interactive visualization of multiple genome alignments.

python maf2xmfa.py -i chr09.maf -o chr09.xmfa

------------------------------------------------------------
Visualization
------------------------------------------------------------

The resulting XMFA file can be opened directly in:

- Mauve: for visualization of locally collinear blocks, rearrangements,
  and large-scale structural variation.

- Geneious Prime: for integrated visualization alongside genome annotations
  and sequence features.

## Notes

- This repository documents workflows and commands only  
- Raw sequencing data, assemblies, and intermediate files are not included  
- File paths and HPC-specific settings are intentionally generalized
