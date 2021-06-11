# SCCS Vitamin D Binding Protein Study

This repository contains the code utilized for the manuscript: _"Genome-wide association study (GWAS) and functional analysis of genetic variants associated with circulating vitamin D binding protein in African Americans"_ (to be submitted)

The directory structure is:

- `01_quality_control/` - QC manipulations to get a first dataset
- `02_merge_imputation/` - All manipulations to impute the different cohorts, and then join the imputed results to build a full genotype dataset.
- `03_gwas_analysis/` - Data analysis presented in the manuscript
- `04_functional_analysis/` - Manual figures of FUMA analysis.

---

The requirements to use this are:

- [HTCondor v8.6.13](https://research.cs.wisc.edu/htcondor/)
- [PLINK v1.9](https://www.cog-genomics.org/plink/1.9/)
- [PLINK2 v2.0](https://www.cog-genomics.org/plink/2.0/)
- [Minimac4 v1.0.1](https://genome.sph.umich.edu/wiki/Minimac4)
- [R v4.0.2](https://www.r-project.org/)

---
