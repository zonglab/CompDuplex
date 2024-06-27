
# CompDuplex Pipeline README

## Data Acquisition

Data acquisition follows the protocol described in the [CompDuplex: Accurate Detection of Somatic Mutations](https://www.protocols.io/view/compduplex-accurate-detection-of-somatic-mutations-kxygx3x4og8j/v1).

## Running the Pipeline

To run this pipeline, execute the following command:
```bash
bash SCRIPTs/s1_CompDup_00_submit_job.sh <working_dir> <fq___prefix> <normal_bulk_path>
```

### Parameters:
- `<working_dir>`: Ideally, the absolute path to where you want to store the analysis results.
- `<fq___prefix>`: Path and prefix to the sequenced fastq files, e.g., `/PATH/TO/PREFIX` corresponding to `/PATH/TO/PREFIX_R1[_001].fastq.gz`.
- `<normal_bulk>`: Path to the matched bulk 30X WGS results. Must include:
  - `05_add_RG.bam`
  - `15_fin_hetero.vcf`
  - `15_fin___homo.vcf`

## Test Dataset

A test dataset is included surrounding hg19 chr14:100,000,000-108,000,000. This dataset was generated in silico by subsampling the reads from a CompDuplex library sequenced on Illumina NovaSeq 6000 S4.

## Analysis Rationale

### Germline Heterozygous Variant Calling from Traditional 30X Bulk WGS

- The 3’ Nextera adapter sequences were trimmed from paired-end sequencing reads using `cutadapt v3.4` with the arguments `-a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -m 50`.
- The reads were mapped to the hg19 reference genome using `bwa-mem v0.7.13-r1126`.
- PCR duplicates were marked using `samtools v1.3.1` `markdup` command.
- Germline variants were called using `samtools mpileup -q 50 -Q 20` followed by `bcftools call`.
- Only heterozygous single nucleotide variants with a calling quality score higher than 50 were kept.
- Tandem regions, centromere regions, and homopolymer regions were filtered out.

### Somatic Mutation Calling from CompDuplex-seq and NanoSeq

#### Preprocessing:
- For CompDuplex-seq, the 3’ Nextera adapter sequences were trimmed using `cutadapt v3.4` with the arguments `-a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -m 25`.
- For NanoSeq, the 3 bases UMI on Read 1 and Read 2 were extracted first,
  then 3’ Truseq adapter sequences were trimmed using `cutadapt v3.4` with the arguments `-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -m 25`.
- The reads were mapped to the hg19 reference genome using `bwa-mem v0.7.13-r1126`.
- Appropriately mapped read pairs mapped to chromosome 1-22, X, and Y (total size: 3,095,677,412 bp) with mapping quality no less than 60 were kept.
- Optical duplicates were removed with `v3.0.0 MarkDuplicates` with the arguments `–OPTICAL_DUPLICATE_PIXEL_DISTANCE=12000` and `--REMOVE_SEQUENCING_DUPLICATES`.
- Raw variant calling was performed by `samtools mpileup` with the argument `--ff QCFAIL`, followed by `bcftools call` with the arguments `-mv -p 0.99 -P 0.99`. All called variants were retained regardless of the variant-calling quality score.
- Tandem regions, centromere regions, and homopolymer regions were filtered out (a total genome mask of 132,519,574 bp), resulting in a list of raw variants.

#### Unique molecule identifier:
- Reads originating from the same gDNA fragment molecules were identified by the same mapping coordinates.
- Reads originating from the same gDNA fragment molecules were identified by the same mapping coordinates and UMIs.

A custom Python script was used to traverse each raw variant. A variant would be called a mutation if it fulfilled the following criteria:
1. Covered in both paired-read directions.
2. At least two reads for each direction with sequencing quality no less than 30.
3. All q>=30 bases support the variant base.
4. More than 60% of all mapped reads at this locus have a q>=30 score.
5. The variant base is localized at least 8 bases away from the 3’ end of the insert gDNA fragment to filter out mapping errors.
6. The variant base is localized more than 13 bases away from the 5’ end of the insert gDNA fragment to filter out mapping errors and newly synthesized bases by Q5 DNA polymerase gap filling.
7. The 5’ soft clip, if present, should not exceed 20 bases in length.

A full list of called somatic mutations is provided in Supplementary Table 6.

### Allele Coverage Depth

On the basis of single molecule detection, by dividing the number of called germline heterozygous mutations (without deduplication of the same mutations called from different original DNA fragments) by the total number of germline heterozygous mutations, we can determine the allele coverage depth for mutation calling (Supplementary Table 1, 3, 4). 

### Somatic Mutation Calling Criteria

A mutation would be called a somatic mutation if:
1. The locus was covered by at least 10 reads in the traditional 30X bulk WGS with sequencing quality >=30 (the fraction of genomic loci fulfilling this requirement will be used to correct the genome-wide somatic mutation load estimation), but zero read supports the variant call.
2. The variant is not annotated in the dbSNP database.

### Somatic Mutation Burden

Somatic mutation burden (in the unit of mutation counts per cell, equivalent to 5,720,696,467 bp for the male or 5,812,746,744 for the female after the genome mask) was calculated by dividing the number of somatic mutations detected by the allele coverage depth.

### Note on Embryonic Mutations

For equal footing comparison of cord blood mutation burdens and mutational signatures, this pipeline did not filter out the low-frequency (variant allele frequency > 0.01) embryonic mutations as described in Abascal et al.

## Package Requirements

- `cutadapt 4.8`
- `bwa 0.7.13-r1126`
- `picard 3.1.1`
- `samtools 1.12`
- `bcftools 1.3`
- `seqtk 1.4-r122`
- `SnpSift 4.1k`

### Python Packages

- `pysam 0.22.1`
- `numpy 1.26.4`
- `pandas 2.2.2`

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

Special thanks to the contributors and the scientific community for their invaluable support and guidance.
