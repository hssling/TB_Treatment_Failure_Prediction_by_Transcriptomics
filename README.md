# TB Treatment Failure Prediction

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

**Predicting Tuberculosis Treatment Failure from Baseline Blood Transcriptomics using Machine Learning**

---

## 📄 Publication

**Title**: Baseline Blood Transcriptomic Signatures Predict Treatment Failure in Tuberculosis: A Machine Learning Study

**Authors**: [Author List]

**Journal**: EBioMedicine (submitted January 2026)

**Preprint**: [bioRxiv/medRxiv DOI] (if applicable)

---

## 🎯 Overview

This repository contains all code and analysis for predicting tuberculosis (TB) treatment failure from baseline (pre-treatment) blood gene expression using machine learning.

### Key Findings

- **Performance**: XGBoost achieved AUC 0.794 in nested cross-validation
- **Top Biomarker**: USP30 (ubiquitin-specific protease)
- **Clinical Utility**: At 90% sensitivity → 60% specificity
- **Novelty**: First study to predict TB treatment outcomes from baseline transcriptomics

---

## 📊 Data

### Source
- **Dataset**: GSE89403 (Berry et al., South Africa)
- **Platform**: RNA-seq (Illumina)
- **Samples**: 254 baseline samples with treatment outcomes
- **Outcomes**: 247 cures, 7 failures (2.8% failure rate)
- **Features**: 16,147 genes

### Access
Data are publicly available from GEO:
```
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89403
```

---

## 🚀 Quick Start

### Prerequisites
```bash
Python 3.9+
pandas >= 1.3.0
numpy >= 1.21.0
scikit-learn >= 1.0.0
xgboost >= 1.5.0
shap >= 0.41.0
matplotlib >= 3.4.0
```

### Installation
```bash
# Clone repository
git clone https://github.com/[username]/TB-Treatment-Failure-Prediction.git
cd TB-Treatment-Failure-Prediction

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### Run Analysis
```bash
# Complete pipeline (from data download to results)
python run_pipeline.py

# Individual steps
python workflow/scripts/fetch_geo_metadata.py
python workflow/scripts/ingest_expression.py
python workflow/scripts/build_dataset.py
python workflow/scripts/train_models.py
python workflow/scripts/external_validation.py

# Generate figures
python generate_shap_plots.py
python generate_roc_curves.py
python pathway_enrichment.py
```

---

## 📁 Repository Structure

```
TB-Treatment-Failure-Prediction/
│
├── README.md                          # This file
├── LICENSE                            # MIT License
├── requirements.txt                   # Python dependencies
├── environment.yml                    # Conda environment
├── .gitignore                         # Git ignore rules
│
├── config/                            # Configuration files
│   ├── config.yaml                    # Main configuration
│   └── label_rules.yaml               # Outcome label mapping
│
├── workflow/                          # Analysis pipeline
│   └── scripts/
│       ├── fetch_geo_metadata.py      # Download GEO metadata
│       ├── ingest_expression.py       # Load expression data
│       ├── build_cohorts_table.py     # Build cohort table
│       ├── build_dataset.py           # Create ML dataset
│       ├── train_models.py            # Train ML models
│       └── external_validation.py     # External validation
│
├── analysis/                          # Analysis scripts
│   ├── generate_shap_plots.py         # SHAP analysis
│   ├── generate_roc_curves.py         # ROC curve plots
│   ├── pathway_enrichment.py          # Pathway analysis
│   └── map_gene_symbols.py            # Gene ID mapping
│
├── outputs/                           # Generated outputs
│   ├── metadata/                      # GEO metadata
│   ├── expression/                    # Expression matrices
│   ├── dataset/                       # ML-ready datasets
│   └── models/                        # Trained models
│
├── reports/                           # Results and figures
│   ├── figures/                       # Publication figures
│   └── tables/                        # Publication tables
│
├── docs/                              # Documentation
│   ├── METHODS.md                     # Detailed methods
│   ├── RESULTS.md                     # Complete results
│   └── TUTORIAL.md                    # Step-by-step tutorial
│
└── manuscript/                        # Manuscript files
    ├── MANUSCRIPT.md                  # Main manuscript
    ├── REFERENCES.md                  # References with DOIs
    ├── AUTHORS.md                     # Author contributions
    └── COVER_LETTER.md                # Submission cover letter
```

---

## 🔬 Methods

### Machine Learning Approach

**Models Tested**:
1. Logistic Regression (L2-regularized)
2. Random Forest (100 trees)
3. XGBoost (gradient boosting)

**Cross-Validation**:
- Nested stratified 3-fold CV
- Outer loop: Performance estimation
- Inner loop: Hyperparameter tuning

**Class Imbalance Handling**:
- `class_weight='balanced'` (Logistic Regression, Random Forest)
- `scale_pos_weight=35.3` (XGBoost)

**Feature Interpretation**:
- SHAP (SHapley Additive exPlanations)
- Pathway enrichment (Enrichr API)

### Data Processing

1. **Quality Control**: Remove low-variance genes
2. **Filtering**: Baseline samples only (prevent data leakage)
3. **Normalization**: Log2-transformed expression values
4. **Feature Selection**: All 16,147 genes (no pre-filtering)

---

## 📈 Results

### Model Performance

| Model | Mean AUC | 95% CI | PR-AUC |
|-------|----------|--------|--------|
| **XGBoost** | **0.794** | 0.699-0.854 | 0.112 |
| Random Forest | 0.774 | 0.605-0.878 | 0.138 |
| Logistic Regression | 0.622 | 0.548-0.683 | 0.059 |

### Top Predictive Genes

| Rank | Gene | Mean |SHAP| | Function |
|------|------|-------------|----------|
| 1 | USP30 | 1.369 | Ubiquitin-specific protease |
| 2 | TMEM132D | 0.655 | Transmembrane protein |
| 3 | CRIP2 | 0.500 | Cysteine-rich protein |
| 4 | BRF1 | 0.456 | RNA polymerase III transcription |
| 5 | TYW1 | 0.391 | tRNA modification |

### Clinical Utility

At 90% sensitivity:
- **Specificity**: 60%
- **PPV**: 6% (given 2.8% prevalence)
- **NPV**: 99.7%

---

## 📊 Figures

All figures are available in `reports/figures/`:

1. **ROC Curves** (`roc_curves_combined.png`)
2. **SHAP Importance** (`shap_importance_bar.png`)
3. **SHAP Summary** (`shap_summary_plot.png`)
4. **SHAP Dependence** (`shap_dependence_*.png`)

---

## 🔄 Reproducibility

### Docker Container
```bash
# Pull Docker image
docker pull [dockerhub/tb-treatment-failure]:latest

# Run analysis
docker run -v $(pwd)/outputs:/outputs [dockerhub/tb-treatment-failure]
```

### Conda Environment
```bash
# Create environment
conda env create -f environment.yml
conda activate tb-treatment-failure

# Run analysis
python run_pipeline.py
```

### Manual Setup
```bash
# Install dependencies
pip install -r requirements.txt

# Download data (automatic)
python workflow/scripts/fetch_geo_metadata.py

# Run pipeline
python run_pipeline.py
```

---

## 📝 Citation

If you use this code or data, please cite:

```bibtex
@article{author2026treatment,
  title={Baseline Blood Transcriptomic Signatures Predict Treatment Failure in Tuberculosis: A Machine Learning Study},
  author={[Author List]},
  journal={EBioMedicine},
  year={2026},
  doi={10.xxxx/xxxxx}
}
```

---

## 🤝 Contributing

We welcome contributions! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

---

## 📧 Contact

**Corresponding Author**: [Name]  
**Email**: [email@institution.edu]  
**Institution**: [Institution]

**Issues**: Please use GitHub Issues for bug reports and feature requests.

---

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 🙏 Acknowledgments

- Berry laboratory for making GSE89403 data publicly available
- Original study participants
- [Institution] for computational resources
- Reviewers and collaborators for helpful feedback

---

## 🔗 Links

- **Publication**: [EBioMedicine DOI]
- **Preprint**: [bioRxiv/medRxiv DOI]
- **Data**: [GEO GSE89403](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89403)
- **Models**: [Zenodo DOI]
- **Docker**: [Docker Hub]

---

## 📊 Project Status

- ✅ Data analysis complete
- ✅ Manuscript submitted
- ⏳ Under peer review
- ⏳ Awaiting publication

**Last Updated**: December 25, 2025

---

## ⚠️ Disclaimer

This is a research tool and should not be used for clinical decision-making without proper validation. The model was developed on a small cohort (n=254, 7 failures) and requires prospective validation before clinical use.

---

## 🔮 Future Work

1. **External validation** in independent cohorts
2. **Prospective validation** in clinical trials
3. **Integration with clinical factors** (HIV, drug resistance)
4. **Mechanistic studies** of top genes (USP30, TMEM132D, CRIP2)
5. **Multi-site validation** across geographies

---

**Star ⭐ this repository if you find it useful!**
