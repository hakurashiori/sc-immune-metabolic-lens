# Single-cell immune-metabolic lens on vaccine responsiveness

Analysis code accompanying the manuscript:

> **A pre-existing MT-ND3-linked OXPHOS-dominant monocyte state constrains vaccine responsiveness in immunocompromised patients**  
> Zi-Yun Peng, Cheng-Hsu Chen, Chih-Le Lin, Hung-Ting Liao, Jing-Rong Wang, and Tai-Ming Ko.

This repository contains R and Python scripts for single-cell transcriptomic, immunometabolic, and immune-repertoire analyses of vaccine responsiveness. The study centers on the **pre-vaccination monocyte state** associated with subsequent serological non-response in kidney transplant recipients (KTRs), with analyses of longitudinal samples, independent vaccine-response cohorts, and other immunocompromised settings.

## Study overview

The study profiled **438,088 peripheral blood mononuclear cell (PBMC) transcriptomes** after quality control, including longitudinal samples from **43 KTRs**, together with healthy-reference datasets. The primary KTR cohort received **mRNA-1273**; additional analyses included booster-associated samples and an independent **BNT162b2** KTR cohort.

The pre-vaccination discovery comparison comprised **9 non-responders (NRs)** and **8 responders (Rs)**, classified using post-vaccination antibody measurements. Monocytes were **computationally subset from PBMCs**, not physically purified before sequencing. An initial PBMC-level classical-monocyte annotation was followed by higher-resolution analyses of classical, intermediate, dendritic-cell-like, and macrophage-like monocyte states.

### Main findings

- **Pre-existing immune-metabolic state:** Monocytes from subsequent non-responders showed an MT-ND3-linked, oxidative-phosphorylation (OXPHOS)-biased transcriptional program before vaccination.
- **Immune-metabolic coordination:** Subset-resolved analyses identified reduced alignment between mitochondrial/metabolic activity and antigen-presentation or immune-activation programs in non-responders. Higher metabolic gene activity should not be interpreted as evidence of greater antigen-presenting function.
- **Serological association:** Baseline monocyte NR-signature scores were inversely associated with antibody titers after the first dose (**Spearman ρ = −0.58; P = 0.014; n = 17 donors**). This is a within-discovery-cohort association, not an independent prospective validation of a clinical prediction test.
- **Longitudinal and adaptive context:** Trajectory and cell-communication analyses suggested delayed or altered monocyte-state remodeling and an inferred immunoregulatory milieu. Antigen-specific CD8⁺ T-cell repertoire findings are exploratory and do not establish monocyte-to-T-cell causality.
- **Cross-cohort assessment:** NR-associated signatures were examined in independent vaccine-response cohorts and projected across immunocompromised settings, including autoimmune disease, cancer, HIV, diabetes, and kidney transplantation. These observations support broader relevance but do not establish uniform predictive performance across diseases.
- **Candidate regulatory mechanisms:** Network modeling and *in silico* perturbation nominated **FOXO1** as a possible regulatory node. **MT-ND3** and potential metabolic modulators, including **metformin**, are research hypotheses for future functional or controlled clinical evaluation—not demonstrated causal drivers or treatments.

## Repository structure

The analysis scripts are organized under [`sc-immtblc/`](sc-immtblc/):

| Directory | Contents |
| --- | --- |
| [`clustering/`](sc-immtblc/clustering/) | Seurat/Harmony data integration, clustering, and immune-cell annotation |
| [`scVI/`](sc-immtblc/scVI/) | scVI latent representation and Leiden reclustering |
| [`Cohernce_score/`](sc-immtblc/Cohernce_score/) | Coherence analysis of mitochondrial and antigen-presentation expression programs |
| [`Pathway__enrichment/`](sc-immtblc/Pathway__enrichment/) | Differential-expression and pathway-enrichment analyses |
| [`WGCNA/`](sc-immtblc/WGCNA/) | hdWGCNA co-expression modules and candidate hub genes |
| [`Feature_selection/`](sc-immtblc/Feature_selection/) | Exploratory feature selection and classification |
| [`trajectory/`](sc-immtblc/trajectory/) | Slingshot and Palantir cell-state trajectory analyses |
| [`cellchat/`](sc-immtblc/cellchat/) | Inferred cell–cell communication using CellChat |
| [`scRepertoire/`](sc-immtblc/scRepertoire/) | Antigen-specific T-cell receptor repertoire analysis |
| [`Broad_IC/`](sc-immtblc/Broad_IC/) | NR-signature projection and metabolic associations across immunocompromised datasets |
| [`Celloracle/`](sc-immtblc/Celloracle/) | Gene regulatory network analysis and *in silico* perturbation |

Additional information on available cell-level metadata is provided in [`sc-immtblc/README.md`](sc-immtblc/README.md), and cross-condition metadata are described in [`sc-immtblc/Broad_IC/README.md`](sc-immtblc/Broad_IC/README.md). Directory names are retained as they appear in the repository.

## Getting started

Clone the repository:

```bash
git clone https://github.com/hakurashiori/sc-immune-metabolic-lens.git
cd sc-immune-metabolic-lens
```

The scripts use both **R** (including Seurat, Harmony, CellChat, hdWGCNA, fgsea, Slingshot, and related packages) and **Python** (including Scanpy, scvi-tools, Palantir, CellOracle, scikit-learn, and XGBoost). Consult the imports at the top of the relevant scripts for task-specific dependencies.

**Before running a script:**

1. Obtain the appropriate processed data and check which RNA/ADT/HTO assays, annotation columns, and embeddings are available in the selected object.
2. Change example or local paths (e.g., `./data/sample1.rds`, `raw_count.h5ad`, or `PATH/TO/*.gmt`) to the corresponding local inputs and reference gene sets.
3. Check group labels, batch variables, cluster identifiers, and required metadata against the chosen dataset. Some scripts assume objects have already been loaded into the R environment.
4. Run only the modules appropriate for the available data. Certain Python scripts contain notebook-style commands and may require adaptation to run as standalone `.py` files.

These files document analysis components rather than a fully automated, version-locked, end-to-end pipeline. Paths and parameter settings may need adjustment. Analyses intended to support donor-level inference must account for donor/sample identity to avoid treating cells from the same donor as independent biological replicates.

## Data availability and privacy

- **NCBI Gene Expression Omnibus (GEO):** [GSE347367](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE347367) — study-associated, de-identified single-cell data deposited as described in the manuscript. Public access is subject to the GEO release status.
- **figshare:** [Supporting processed data](https://figshare.com/s/4a855295d192b3aed5ff) — the data link specified in the current manuscript.
- **Supporting Data Values:** Values underlying the manuscript's graphs are supplied in the study's supporting materials.
- **Raw sequencing reads and identifiable clinical data:** Not distributed through this GitHub repository. Human-participant data access and reuse remain subject to privacy, consent, and institutional/ethical restrictions.

This repository contains **analysis code, not patient-level datasets**. The availability of individual data modalities varies across study components; the presence of an analysis script does not imply that all of its inputs are publicly downloadable.

## Interpretation and limitations

The reported NR-associated monocyte signature is an **observational, pre-vaccination host-state indicator**. It does not establish that MT-ND3 expression, OXPHOS activity, FOXO1, or inferred cell–cell signaling directly causes failure of seroconversion. Medication exposures and other clinical factors may modify these associations. Functional perturbation, independent donor-level validation, and prospective studies are required before clinical risk prediction or interventions can be recommended.

## Citation

If you use this repository, please cite the accompanying manuscript:

**Peng Z-Y, Chen C-H, Lin C-L, Liao H-T, Wang J-R, Ko T-M.** *A pre-existing MT-ND3-linked OXPHOS-dominant monocyte state constrains vaccine responsiveness in immunocompromised patients.* Manuscript under editorial consideration. Please use the final journal citation when available.

## License

See [LICENSE](LICENSE) for the repository license.
