# References and data terms

## Data and annotations

- [TCGA-BRCA PanCancer Atlas study](https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018).
- [Pinned DataHub study metadata](https://github.com/cBioPortal/datahub/blob/dda832891555904730395d476c5c888b2e1553e1/public/brca_tcga_pan_can_atlas_2018/meta_study.txt) and [expression metadata](https://github.com/cBioPortal/datahub/blob/dda832891555904730395d476c5c888b2e1553e1/public/brca_tcga_pan_can_atlas_2018/meta_mrna_seq_v2_rsem.txt). The latter identifies the RSEM expression profile and source filename.
- [NCI PanCancer Atlas publication resources](https://gdc.cancer.gov/about-data/publications/pancanatlas).
- [cBioPortal data processing FAQ](https://docs.cbioportal.org/user-guide/faq/) and [DataHub repository](https://github.com/cBioPortal/datahub).
- Hoadley et al. (2018), *Cell of Origin Patterns Dominate the Molecular Classification of 10,000 Tumors from 33 Types of Cancer*, Cell. [DOI: 10.1016/j.cell.2018.03.022](https://doi.org/10.1016/j.cell.2018.03.022). PanCancer molecular classification context.

## PAM50 panel

- Parker et al. (2009), *Supervised Risk Predictor of Breast Cancer Based on Intrinsic Subtypes*, Journal of Clinical Oncology. [DOI: 10.1200/JCO.2008.18.1370](https://doi.org/10.1200/JCO.2008.18.1370).
- [genefu PAM50 reference object](https://github.com/bhklab/genefu/blob/76fc3ca34eee4b55a01219d37da95b7faddaf28f/data/pam50.rda), source of the 50 gene names/Entrez IDs only. No centroid coefficients or genefu code are distributed here.

## Implementation

- [scikit-learn guidance on leakage](https://scikit-learn.org/stable/common_pitfalls.html#data-leakage).
- [RepeatedStratifiedKFold API](https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.RepeatedStratifiedKFold.html).
- [Pipeline API](https://scikit-learn.org/stable/modules/generated/sklearn.pipeline.Pipeline.html).

These public source files are downloaded directly from the provider and remain subject to their source terms. Consult [cBioPortal's permissions/data FAQ](https://docs.cbioportal.org/user-guide/faq/) and [GDC data access policies](https://gdc.cancer.gov/access-data/data-access-processes-and-tools) before downstream reuse. No blanket new license for TCGA/cBioPortal data or a clinical assay is granted by this repository. A project source-code license has not been selected on behalf of all team members.
