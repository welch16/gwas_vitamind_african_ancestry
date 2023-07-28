# Genome-wide association study (GWAS) of circulating vitamin D outcomes among individuals of African ancestry

[![DOI](https://zenodo.org/badge/373951249.svg)](https://zenodo.org/badge/latestdoi/373951249)


L Parlato, R Welch, IM Ong, J Long, L Cai, MD Steinwandel, WJ Blot, W Zheng, S Warren Andersen, _"Genome-wide association study (GWAS) of circulating vitamin D outcomes among individuals of African ancestry"_ (2023), The American Journal of Clinical Nutrition, 117-2, 308-316

The directory structure is:

- `01_quality_control/` - QC manipulations to get a first dataset
- `02_merge_imputation/` - All manipulations to impute the different cohorts, and then join the imputed results to build a full genotype dataset.
- `03_gwas_analysis/` - Data analysis presented in the manuscript

---

The requirements to use this are:

- [HTCondor v8.6.13](https://research.cs.wisc.edu/htcondor/)
- [PLINK v1.9](https://www.cog-genomics.org/plink/1.9/)
- [PLINK2 v2.0](https://www.cog-genomics.org/plink/2.0/)
- [Minimac4 v1.0.1](https://genome.sph.umich.edu/wiki/Minimac4)
- [R v4.0.2](https://www.r-project.org/)

---
