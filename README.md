# ABIDE Comprehensive Analysis Pipeline

This repository contains a unified R script that combines all analysis steps for autism research using the ABIDE (Autism Brain Imaging Data Exchange) dataset. The script consolidates multiple individual R files into a single, efficient, and executable workflow.

## 🎯 Overview

The comprehensive analysis pipeline performs the following steps:

1. **Data Extraction & Processing**: ABIDE neuroimaging data extraction and combination
2. **Phenotype Data Preprocessing**: Clinical/behavioral data normalization and feature creation  
3. **Correlation Analysis**: Brain connectivity analysis
4. **Data Merging**: Combining neuroimaging and phenotypic data
5. **Principal Component Analysis (PCA)**: With varimax and promax rotations
6. **Exploratory Data Analysis (EDA)**: Statistical summaries and visualizations
7. **MANCOVA**: Multivariate analysis of covariance
8. **Support Vector Machine (SVM)**: Classification analysis
9. **Comprehensive Reporting**: Summary statistics and visualizations

## 📋 Requirements

### Required R Packages
The script will automatically install missing packages:
- `readxl`, `readr`, `stringr`
- `dplyr`, `data.table`, `tidyr`, `reshape2`
- `ggplot2`, `corrplot`
- `GPArotation`, `psych`
- `car`, `e1071`, `caret`, `pROC`
- `foreach`, `doParallel`

### Required Data Files
Place these files in your working directory:
- `Phenotypic_V1_0b_v1.csv` - ABIDE phenotypic data
- `CC400_combined.csv` - Combined neuroimaging data **OR**
- `CC400_CPAC/` folder with individual CSV files

## 🚀 Usage

### Method 1: Direct Execution
```bash
Rscript ABIDE_Comprehensive_Analysis.R
```

### Method 2: R Console
```r
source("ABIDE_Comprehensive_Analysis.R")
```

### Method 3: RStudio
1. Open `ABIDE_Comprehensive_Analysis.R` in RStudio
2. Click "Source" or press Ctrl+Shift+S (Cmd+Shift+S on Mac)

## 📁 Output Structure

The script creates organized output directories:

```
├── 01_phenotype_data/
│   └── phenotype_data.csv
├── 02_correlation_analysis/
│   └── correlation_results_wide.csv
├── 03_merged_data/
│   └── merged_data.csv
├── 04_pca/
│   ├── explained_variance_pca.csv
│   ├── varimax_scores.csv
│   └── promax_scores.csv
├── 05_eda/
│   ├── dx_group_summary.csv
│   ├── age_distribution.png
│   ├── dx_group_distribution.png
│   └── social_comm_comparison.png
├── 06_mancova/
│   ├── mancova_communication_results.csv
│   └── mancova_social_results.csv
├── 07_svm/
│   ├── svm_performance_results.csv
│   └── svm_confusion_matrix.csv
├── varimax_alldata.csv
├── promax_alldata.csv
└── ANALYSIS_SUMMARY.txt
```

## 🔧 Configuration

You can modify these parameters at the top of the script:

```r
set.seed(123)  # For reproducibility
top_pcs <- 20  # Number of top principal components to use
test_proportion <- 0.1  # Proportion for test set in SVM
```

## 📊 Key Features

### Efficiency Improvements
- **Parallel Processing**: Utilizes multiple CPU cores
- **Error Handling**: Robust error handling with informative messages
- **Memory Management**: Efficient data handling for large datasets
- **Modular Functions**: Reusable functions reduce code duplication

### Robustness
- **Data Quality Checks**: Automated validation at each step
- **Flexible Input**: Works with existing processed files or raw data
- **Graceful Degradation**: Continues analysis even if some steps fail
- **Comprehensive Logging**: Detailed progress reporting

### Output Quality
- **Organized Results**: Structured output directories
- **Comprehensive Summary**: Detailed analysis report
- **Visualizations**: Automatic generation of key plots
- **Reproducible**: Fixed random seed for consistent results

## 📈 Analysis Steps Details

### 1. Phenotype Processing
- Loads and cleans clinical/behavioral data
- Normalizes assessment scores (0-1 scale)
- Creates composite social and communication scores
- Handles missing data appropriately

### 2. Neuroimaging Processing
- Combines individual CPAC connectivity files
- Performs correlation analysis on brain connectivity
- Creates wide-format correlation matrices

### 3. Data Integration
- Merges phenotypic and neuroimaging data
- Ensures proper subject matching
- Creates analysis-ready datasets

### 4. Principal Component Analysis
- Performs PCA on correlation features
- Applies Varimax and Promax rotations
- Calculates component scores for further analysis

### 5. Statistical Analysis
- MANCOVA with interaction terms
- Tests for group differences in PC space
- Separate analyses for social and communication domains

### 6. Machine Learning
- SVM classification (autism vs. control)
- Cross-validation with hyperparameter tuning
- Performance evaluation with ROC/AUC

## 🛠 Troubleshooting

### Common Issues

1. **Missing Data Files**
   - Ensure `Phenotypic_V1_0b_v1.csv` is in working directory
   - Check that neuroimaging data is available

2. **Memory Issues**
   - Reduce `top_pcs` parameter for large datasets
   - Close other R sessions

3. **Package Installation**
   - Update R to latest version
   - Install packages manually if automatic installation fails

4. **Permission Errors**
   - Ensure write permissions in working directory
   - Run R as administrator if needed

### Performance Tips
- Use SSD storage for faster I/O
- Allocate sufficient RAM (8GB+ recommended)
- Close unnecessary applications during analysis

## 📝 Citation

If you use this analysis pipeline, please cite the original ABIDE project:
- Di Martino, A., et al. (2014). The autism brain imaging data exchange: towards a large-scale evaluation of the intrinsic brain architecture in autism. Molecular psychiatry, 19(6), 659-667.

## 📄 License

This project is provided under the MIT License. See LICENSE file for details.

## 🤝 Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request with detailed description

## 📧 Support

For questions or issues:
1. Check the troubleshooting section
2. Review the ANALYSIS_SUMMARY.txt output
3. Open an issue with detailed error messages

---

**Note**: This script consolidates and optimizes the functionality from multiple individual R files into a single, efficient workflow while maintaining all original analysis capabilities.