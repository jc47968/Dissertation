# ABIDE Comprehensive Analysis - Combined Script

This repository contains a comprehensive R script that combines all the individual analysis workflows from the original ABIDE dataset analysis into one efficient and executable code.

## Overview

The script integrates the following analysis components:
1. **Data Loading and Preprocessing** - Phenotype data processing, normalization, and correlation analysis
2. **Principal Component Analysis** - PCA with Varimax and Promax rotations
3. **Exploratory Data Analysis** - Descriptive statistics and visualizations
4. **Statistical Analysis** - MANCOVA for communication and social variables
5. **Machine Learning** - Support Vector Machine classification
6. **Permutation Testing** - Robust statistical inference
7. **Reporting** - Automated result compilation and visualization

## Requirements

### R Version
- R 4.0.0 or higher

### Required R Packages
The script will automatically install missing packages, but you can install them manually:

```r
install.packages(c(
  "dplyr", "data.table", "readr", "tidyr", "reshape2",
  "ggplot2", "corrplot", "GPArotation", "psych", 
  "e1071", "caret", "pROC", "MVN", "car", "biotools",
  "vegan", "igraph", "ggraph", "tidygraph", "ggimage",
  "png", "grid"
))
```

### Required Data Files
Place these files in your working directory:
- `Phenotypic_V1_0b_v1.csv` - ABIDE phenotypic data (required)
- `varimax_alldata.csv` - Processed Varimax data (optional, will be generated if missing)
- `promax_alldata.csv` - Processed Promax data (optional, will be generated if missing)
- `correlation_results_wide.csv` - Correlation analysis results (optional)

## Usage

### Method 1: Direct Execution
```bash
# Make the script executable (Linux/Mac)
chmod +x ABIDE_Comprehensive_Analysis_Combined.R

# Run the script
./ABIDE_Comprehensive_Analysis_Combined.R
```

### Method 2: From R Console
```r
# Set working directory to where your data files are located
setwd("/path/to/your/data")

# Source the script
source("ABIDE_Comprehensive_Analysis_Combined.R")
```

### Method 3: From Command Line
```bash
# Run with Rscript
Rscript ABIDE_Comprehensive_Analysis_Combined.R
```

## Output Structure

The script creates an organized output directory structure:

```
output/
├── data/
│   ├── phenotype_data.csv
│   ├── varimax_scores.csv
│   ├── promax_scores.csv
│   ├── varimax_alldata.csv
│   └── promax_alldata.csv
├── results/
│   ├── explained_variance_pca.csv
│   ├── demographics_summary.csv
│   ├── pc_summary_by_group.csv
│   ├── mancova_varimax_comm.csv
│   ├── mancova_promax_comm.csv
│   ├── mancova_varimax_social.csv
│   ├── mancova_promax_social.csv
│   ├── svm_performance.csv
│   ├── perm_mancova_varimax_comm.csv
│   └── perm_mancova_promax_comm.csv
├── plots/
│   ├── scree_plot.png
│   ├── age_distribution.png
│   ├── dx_group_distribution.png
│   ├── pc_correlation_matrix.png
│   └── roc_curves.png
└── analysis_report.txt
```

## Script Features

### Intelligent Data Handling
- Automatically detects available data files
- Gracefully handles missing optional files
- Provides informative progress messages
- Implements error handling for robust execution

### Comprehensive Analysis Pipeline
- **Data Preprocessing**: Normalization, cleaning, composite variable creation
- **PCA Analysis**: Full PCA with both rotation methods and variance analysis
- **Statistical Testing**: MANCOVA with assumption testing and permutation validation
- **Machine Learning**: SVM with cross-validation and performance metrics
- **Visualization**: Automated plot generation for key findings

### Optimized Performance
- Modular function design for clarity and maintainability
- Efficient data processing with data.table
- Parallel-ready code structure
- Memory-efficient operations

## Analysis Workflow

1. **Setup**: Loads libraries and creates output directories
2. **Data Processing**: Loads and preprocesses phenotype data, handles correlation data
3. **PCA**: Performs principal component analysis with rotations
4. **EDA**: Generates descriptive statistics and basic visualizations
5. **MANCOVA**: Conducts multivariate analysis of covariance
6. **SVM**: Trains and evaluates support vector machine models
7. **Permutation**: Performs permutation testing for robust inference
8. **Reporting**: Generates comprehensive analysis report

## Troubleshooting

### Common Issues

1. **Missing Data Files**
   - Ensure `Phenotypic_V1_0b_v1.csv` is in your working directory
   - Check file permissions and path accessibility

2. **Package Installation Errors**
   - Run R as administrator/sudo if needed
   - Update R to the latest version
   - Install packages manually if automatic installation fails

3. **Memory Issues**
   - Increase R memory limit: `memory.limit(size = 8000)` (Windows)
   - Close other applications to free up RAM
   - Consider running on a machine with more memory

4. **Analysis Skipped Messages**
   - These are normal when optional data files are missing
   - The script will run available analyses and skip missing components

### Performance Tips
- Use SSD storage for faster I/O operations
- Increase available RAM for large datasets
- Run in a clean R session for optimal performance

## Customization

You can modify the script parameters:
- Change the number of principal components (default: 20)
- Adjust SVM cross-validation folds (default: 5)
- Modify permutation test iterations (default: 1000)
- Update output file formats and locations

## Support

For issues or questions:
1. Check the console output for specific error messages
2. Verify all required data files are present
3. Ensure R packages are properly installed
4. Review the generated `analysis_report.txt` for completion status

## Citation

If you use this script in your research, please cite the original ABIDE dataset and relevant methodological papers for the statistical and machine learning methods employed.