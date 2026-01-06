# Peanut-PanMAGIC-Platform
This pangenome was constructed using 18 peanut (Arachis hypogaea) lines representing diverse cultivated germplasm. Graph-based integration of whole-genome assemblies captures SNPs, structural variants, and presence–absence variation, enabling comparative analyses, introgression detection, and variant density profiling across chromosomes.

Peanut_PanMAGIC_Platform

This repository hosts pipelines for peanut genome assembly and pangenome analysis.



# Genome Assembly Pipeline

This branch documents the genome assembly workflow used for generating assemblies
of the 18 parental peanut genomes used in the pangenome construction.

---

## Biological Material and Sequencing

For the 18 parental genomes, fresh leaf tissues were collected from each parent to
obtain at least 5 g of fresh tissue. DNA was extracted from the collected tissue
and sequencing libraries were prepared for PacBio Sequel II HiFi sequencing.
Each genotype was sequenced to approximately 20× coverage.

The genome assembly was conducted using **hifiasm**. The assembled genomes were
scaffolded with the progenitor peanut genomes (*Arachis ipaensis* and
*Arachis duranensis*) using **RagTag**.

---

## Steps Involved

---

## 1. Adapter Filtration

Adapters were filtered using **HiFiAdapterFilt**
(https://github.com/sheinasim/HiFiAdapterFilt).

```bash
ml HiFiAdapterFilt/2.0.1-gompi-2020b
bash hifiadapterfilt.sh -t 30 -p Prefix -o Outputdirectory
```

---

## 2. Genome Assembly with hifiasm

### A. Homozygous assemblies

i. A total of 16 out of 18 genomes were assembled using the following pipeline:

```bash
module load Hifiasm/0.13-foss-2019b
hifiasm -o Name.hifiasm -t 64 -lo Name.filt.fastq.gz
```

The `-lo` flag was specified for homozygous genomes.

Primary contigs were converted to FASTA format using:

```bash
awk '/^S/{print ">"$2;print $3}' C99R.HiFi.bp.p_ctg.gfa > C99R.HiFi.bp.p_ctg.fa
```

---

ii. Due to tandem duplication in some regions of the genomes GA-12Y and MarcI,
the primary contigs had contigs smaller than 100 kb filtered out. This filtering
resulted in improved dot plots.

```bash
ml SeqKit/2.5.1
seqkit seq -m 100000 contig.fa > filtered_contig.fa
```

---

### B. Heterozygous assembly with Hi-C integration (NMValA)

Dot plots generated using the `-lo` flag showed potential tandem duplication or
heterozygosity. The raw reads of New Mexico Valencia (NMValA) were analyzed using
GenomeScope2, which indicated that the genome was heterozygous.

The genome was therefore assembled by phasing into two haplotypes. Hi-C reads
were used, and the coverage peak was adjusted to 27 after observing the coverage
distribution using GenomeScope2.

```bash
ml hifiasm/0.19.6-GCCcore-11.3.0

hifiasm -o NMValA.HIFIasm -t64 --hom-cov 27 \
  --h1 NMValA_OmniC_R1.fastq.gz \
  --h2 NMValA_OmniC_R2.fastq.gz \
  NMValA_HIFI.filt.fastq.gz
```

---

### C. Chromosome-specific heterozygosity (CC477)

CC477 showed homozygosity across most chromosomes; however, chromosome 9 appeared
heterozygous or contained unexpected structural variants. Reads corresponding to
chromosome 9 were extracted and assembled as heterozygous (default hifiasm),
while the remaining chromosomes were assembled as homozygous using the `-l0`
flag.

Hap1 and Hap2 assemblies were generated, differing only in chromosome 9.

---

### Mapping with reference genome

```bash
minimap2 -d $workdir/Tifrunner.mmi $Tifref
minimap2 -ax asm10 Tifrunner.mmi $workdir/CC477_HIFI.filt.fastq.gz map-hifi -t 64 > CC477_HIFI.sam
samtools view -Sb CC477_HIFI.sam -@ 64 > CC477_HIFI.bam
samtools sort CC477_HIFI.bam -@ 64 -o CC477_HIFI_sorted.bam
samtools index CC477_HIFI_sorted.bam
```

---

### Extracting chromosome 9

```bash
samtools view -b -h CC477_HIFI_sorted.bam "arahy.Tifrunner.gnm2.chr09" > CC477_HIFI_chr09.bam
samtools sort CC477_HIFI_chr09.bam -@ 64 -o CC477_HIFI_chr09_sorted.bam
samtools index CC477_HIFI_chr09_sorted.bam
```

---

### Testing the BAM file

```bash
samtools view -c -F 4 CC477_HIFI_chr09.bam
samtools quickcheck CC477_HIFI_chr09.bam
```

---

### BAM to FASTQ conversion

```bash
samtools view -b CC477_HIFI_chr09.bam | samtools fastq - > CC477_HIFI_chr09.fastq
```

---

### Genome assembly of chromosome 9

```bash
ml hifiasm/0.19.6-GCCcore-11.3.0
hifiasm -o CC477.chr09.HIFIasm -t64 --hom-cov 26 CC477_HIFI_Chr09.fastq.gz
```

---

### Extracting all but chromosome 9

```bash
samtools view -H CC477_HIFI_sorted.bam | grep -v 'arahy.Tifrunner.gnm2.chr09' > header.txt
samtools view CC477_HIFI_sorted.bam | grep -v 'arahy.Tifrunner.gnm2.chr09' | \
  cat header.txt - | samtools view -b - > CC477_HIFI_No_Chr09.bam

samtools sort CC477_HIFI_No_Chr09.bam -@ 64 -o CC477_HIFI_No_Chr09.sorted.bam
samtools index CC477_HIFI_No_Chr09.sorted.bam
```

```bash
mv CC477_HIFI_No_Chr09.sorted.bam NO.Chr09.CC477.HIFI.sorted.bam
```

---

### Testing BAM file

```bash
ml SAMtools/1.18-GCC-12.3.0
samtools quickcheck NO.Chr09.CC477.HIFI.sorted.bam
samtools index NO.Chr09.CC477.HIFI.sorted.bam
samtools idxstats NO.Chr09.CC477.HIFI.sorted.bam
```

---

### BAM to FASTQ conversion

```bash
samtools view -b NO.Chr09.CC477.HIFI.sorted.bam | samtools fastq - > NO.Chr09.CC477.HIFI.sorted.fastq
```

---

### Genome assembly of remaining chromosomes

```bash
ml hifiasm/0.19.6-GCCcore-11.3.0
hifiasm -o NO.Chr09.CC477.HIFIasm -t 64 -l0 NO.Chr09.CC477.HIFI.sorted.fastq
```

---

### Concatenating haplotypes

```bash
cat CC477.Chr09.HIFIasm.bp.hap1.p_ctg.fa NO.Chr09.CC477.HIFIasm.bp.p_ctg.fa > CC477.complete.Hap1.fasta
cat CC477.Chr09.HIFIasm.bp.hap2.p_ctg.fa NO.Chr09.CC477.HIFIasm.bp.p_ctg.fa > CC477.complete.Hap2.fasta
```

---

## Scaffolding to the progenitor reference genome

The assembled genomes were scaffolded using the progenitor genome V1, where
*A. duranensis* and *A. ipaensis* were concatenated into a single reference file.

```bash
ml RagTag/2.1.0
ragtag.py scaffold -t 64 -r -o Folder_name WildGenome.fa contigs.fa
```

---

## Genome alignment observation with dot plots

Primary contigs or Hap1/Hap2 contigs were plotted against the Tifrunner genome
(V2). Scaffolded genomes were also compared to the Tifrunner V2 genome.

```bash
module load MUMmer/4.0.0rc1-GCCcore-11.3.0
nucmer -t 64 Tifref.fasta ragtag.scaffold.fasta -p mummerplot
```

```bash
ml R/4.1.2-foss-2021b
show-coords -c mummerplot.delta > mummerplot.coords
./mummerCoordsDotPlotly.R -i mummerplot.coords -o mummerpicture -s -t -m 15000 -q 200000 -k 20 -l
```

---
