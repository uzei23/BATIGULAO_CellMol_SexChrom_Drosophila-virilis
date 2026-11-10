# Sex-Chromosome Dosage Analysis in *Drosophila virilis*

## Project Overview

This project estimates sex-chromosome dosage in *Drosophila virilis* using publicly available genome assembly and whole-genome sequencing (WGS) data. The analysis uses NCBI, Galaxy, and R to compare the sequencing read depth of the X chromosome with the mean read depth of the autosomes.

## Study Information

- **Species:** *Drosophila virilis*
- **Reference assembly:** GCF_030788295.1
- **Assembly name:** Dvir_AGI_RSII-ME
- **Assembly BioSample:** SAMN35709873
- **WGS Run:** SRR6399436
- **WGS Experiment:** SRX3492610
- **WGS BioSample:** SAMN08216748
- **Recorded WGS sample sex:** Male
- **Strain:** 15010-1051.87

The assembly and WGS data have different BioSample accessions. They are from the same recorded strain, but they are not confirmed to represent the same individual.

## Data Sources

- [NCBI Reference Assembly](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_030788295.1/)
- [NCBI SRA Experiment](https://www.ncbi.nlm.nih.gov/sra/SRX3492610)
- [NCBI SRA Run](https://trace.ncbi.nlm.nih.gov/Traces/?run=SRR6399436)
- [NCBI WGS BioSample](https://www.ncbi.nlm.nih.gov/biosample/SAMN08216748)

## Methods

1. Downloaded the chromosome-level reference genome from NCBI.
2. Retrieved paired-end WGS reads for SRR6399436.
3. Evaluated read quality using FastQC.
4. Aligned paired-end reads to the reference genome using BWA-MEM2 in Galaxy.
5. Generated chromosome coverage statistics using Samtools coverage.
6. Imported the coverage data and chromosome classifications into R.
7. Calculated the mean autosomal depth and X-to-autosome (X:A) depth ratio.
8. Generated a chromosome relative-depth table and graph.

Mapping and base quality thresholds were set to 20.

## Chromosome Classification

| Reference sequence | Chromosome | Group |
|---|---|---|
| NC_091543.1 | X | X |
| NC_091544.1 | 2 | A |
| NC_091545.1 | 3 | A |
| NC_091546.1 | 4 | A |
| NC_091547.1 | 5 | A |
| NC_091548.1 | 6 | A |

Only the six primary chromosomes were included in the main analysis. Unplaced scaffolds were excluded.

## Results

- **X chromosome mean depth:** 3.59703
- **Mean autosomal depth:** 6.57621
- **X:A ratio:** 0.547
- **X:A ratio excluding chromosome 6:** 0.510

The observed X:A ratio is close to 0.5, consistent with one X chromosome relative to two autosomal copies. This result agrees with the male sex recorded for the WGS sample.

## Biological Interpretation

In *Drosophila*, sex determination depends primarily on the ratio of X chromosomes to sets of autosomes. An X:A ratio near 0.5 is consistent with a single X chromosome relative to two autosomal sets. The Y chromosome is primarily required for male fertility rather than being the primary determinant of sex.

## Limitations

1. The reference assembly and WGS data have different BioSample accessions, so they are not confirmed to originate from the same individual.
2. Mapping bias and repetitive sequences may influence read-depth estimates.
3. Chromosome 6 showed lower coverage than the other autosomes.
4. Relative read depth alone cannot independently confirm biological sex.
5. The interpretation thresholds are simplified classroom guidelines.

## Repository Contents

- `sex_depth_analysis.R` - R script for chromosome-depth analysis.
- `data/` - Chromosome classification and Galaxy coverage data.
- `results/chromosome_relative_depth.csv` - Relative chromosome-depth results.
- `results/chromosome_depth.png` - Chromosome-depth graph.
- `results/sex_inference.txt` - Summary of sex-chromosome dosage interpretation.

## Galaxy Workflow Documentation

Figure 1. Retrieval of the chromosome-level reference genome of Drosophila virilis from NCBI.
<img width="1702" height="921" alt="image" src="https://github.com/user-attachments/assets/ace487bd-eaea-4740-8eaa-06d08eecdf25" />

Figure 2. Selection of paired-end whole-genome sequencing reads from the NCBI Sequence Read Archive.
<img width="987" height="652" alt="image" src="https://github.com/user-attachments/assets/40ece2de-066d-4e34-abc3-922c482632d5" />

Figure 3. Importation of the reference genome and paired-end sequencing reads into Galaxy.
<img width="1141" height="875" alt="image" src="https://github.com/user-attachments/assets/88bd07ef-9407-4b50-888a-5e0be0e497d9" />

Figure 4. Assessment of paired-end read quality using FastQC. The per-base sequence quality profiles were within the acceptable range.
<img width="1786" height="872" alt="image" src="https://github.com/user-attachments/assets/9e55244e-48c0-4462-8a3f-e331403155fd" />

<img width="1786" height="875" alt="image" src="https://github.com/user-attachments/assets/99bbec04-411e-47f8-b975-eb108bf87d70" />

Figure 5. Alignment of paired-end WGS reads against the D. virilis reference genome using BWA-MEM2.
<img width="1787" height="877" alt="image" src="https://github.com/user-attachments/assets/2d90ec3a-563f-43e5-81de-b844ecabf76b" />

<img width="1782" height="877" alt="image" src="https://github.com/user-attachments/assets/79f6c8e1-4d3c-43ee-9147-89ba2fa92aca" />


Figure 6. Generation of a coordinate-sorted BAM file containing the mapped sequencing reads.
<img width="590" height="366" alt="image" src="https://github.com/user-attachments/assets/28d75774-fb43-43fb-8c29-b5ea1e925371" />

Figure 7. Calculation of chromosome-level sequencing coverage using Samtools coverage with mapping and base quality thresholds of 20.
<img width="1792" height="878" alt="image" src="https://github.com/user-attachments/assets/a2d5c835-0bd7-4537-9af7-f76302745f9e" />

Figure 8. Chromosome coverage output used for downstream relative-depth calculations in R.
<img width="1793" height="878" alt="image" src="https://github.com/user-attachments/assets/fe45ede7-ee02-4584-8b54-a953fbbb55d7" />

## Questions for the Report

**1. What species did you choose, and what evidence supports a diploid XX/XY or ZZ/ZW system?**

The selected species was *Drosophila virilis* (fruit fly), which has an XX/XY sex chromosome system. Females typically have two X chromosomes (XX), while males have one X and one Y chromosome (XY). Its chromosome-level reference assembly identifies chromosome X separately from the autosomes, supporting its classification as an XX/XY species.

**2. What are your assembly accession, WGS run, and both BioSample accessions? Are the donors the same?**

- Assembly accession: GCF_030788295.1
- WGS run accession: SRR6399436
- Assembly BioSample accession: SAMN35709873
- WGS BioSample accession: SAMN08216748

The donors are not confirmed to be the same because the assembly and WGS data have different BioSample accessions. However, both are associated with the same strain, 15010-1051.87.

**3. Which sequences were used as autosomes and as X or Z? How did you verify their identities?**

Chromosome X (NC_091543.1) was used as the sex chromosome, while chromosomes 2, 3, 4, 5, and 6 (NC_091544.1 to NC_091548.1) were used as autosomes. Their identities were verified using the NCBI assembly report and chromosome annotations of *Drosophila virilis*. Unplaced scaffolds and mitochondrial sequences were excluded from the primary analysis.

**4. What are the median autosomal depth, X:A or Z:A ratio, and inferred sex-chromosome complement?**

The median autosomal mean depth was approximately 7.04858, based on chromosomes 2, 3, 4, 5, and 6. The mean autosomal depth was approximately 6.57621, while the X chromosome mean depth was 3.59703. This resulted in an X:A ratio of approximately 0.547.

The ratio is close to 0.5, suggesting one X chromosome relative to two copies of the autosomes. Therefore, the inferred sex-chromosome complement is XY, consistent with a male individual.

**5. Does your result agree with the recorded sex of the WGS donor, if reported? Explain any discrepancy or uncertainty.**

Yes. The inferred XY complement agrees with the recorded male sex of the WGS donor in NCBI BioSample SAMN08216748. However, this inference is based on relative sequencing depth and does not directly confirm the presence of a Y chromosome. Differences in coverage, mapping bias, and the use of an assembly from a different BioSample may introduce uncertainty.

**6. Why is a chromosome-level assembly alone insufficient for a reliable read-depth sex inference?**

A chromosome-level assembly provides the reference sequences and chromosome identities, but it does not provide information about the actual sequencing depth of an individual sample. Reliable sex inference requires sequencing reads from a known biological sample, alignment of those reads to the reference genome, and comparison of coverage between the sex chromosome and autosomes. Factors such as sequencing depth, mapping bias, repetitive regions, and differences between the reference assembly donor and WGS donor can also affect the results.


## References

- NCBI Genome Data Viewer and Datasets: *Drosophila virilis*, GCF_030788295.1.
- NCBI Sequence Read Archive: SRR6399436.
- Galaxy Training Materials and Galaxy tools used in this analysis.
- Samtools documentation.
- BWA-MEM2 documentation.
