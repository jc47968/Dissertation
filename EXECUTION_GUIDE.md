# ABIDE Comprehensive Analysis - Execution Guide

## 🚀 Quick Start

This guide provides step-by-step instructions for executing the comprehensive ABIDE analysis pipeline.

## 📋 System Requirements

### 1. R Installation
Ensure R (version 4.0+) is installed on your system:

#### Windows:
- Download from: https://cran.r-project.org/bin/windows/base/
- Install R and optionally RStudio

#### macOS:
```bash
# Using Homebrew
brew install r

# Or download from: https://cran.r-project.org/bin/macosx/
```

#### Linux (Ubuntu/Debian):
```bash
sudo apt update
sudo apt install r-base r-base-dev
```

#### Linux (CentOS/RHEL):
```bash
sudo yum install R
```

### 2. Verify R Installation
```bash
R --version
# Should show R version 4.0.0 or higher
```

## 📁 Data Preparation

### Required Files
Ensure these files are in your working directory:

1. **Phenotypic_V1_0b_v1.csv** - ABIDE phenotypic data
   - Download from: https://fcp-indi.s3.amazonaws.com/data/Projects/ABIDE/Phenotypic_V1_0b_v1.csv

2. **Neuroimaging Data** (one of the following):
   - **Option A**: `CC400_combined.csv` - Pre-combined neuroimaging data
   - **Option B**: `CC400_CPAC/` folder with individual CSV files

### Data Structure Check
```bash
# Verify your directory structure
ls -la
# Should show:
# - ABIDE_Comprehensive_Analysis.R
# - Phenotypic_V1_0b_v1.csv
# - CC400_combined.csv (or CC400_CPAC/ folder)
```

## 🎯 Execution Methods

### Method 1: Command Line (Recommended)
```bash
# Navigate to your project directory
cd /path/to/your/project

# Execute the script
Rscript ABIDE_Comprehensive_Analysis.R
```

### Method 2: R Console
```r
# Set working directory
setwd("/path/to/your/project")

# Source the script
source("ABIDE_Comprehensive_Analysis.R")
```

### Method 3: RStudio
1. Open RStudio
2. Open `ABIDE_Comprehensive_Analysis.R`
3. Set working directory: `Session > Set Working Directory > To Source File Location`
4. Click `Source` button or press `Ctrl+Shift+S` (Windows/Linux) or `Cmd+Shift+S` (Mac)

## ⏱️ Expected Runtime

| Dataset Size | Estimated Time | Memory Usage |
|-------------|---------------|--------------|
| Small (< 200 subjects) | 5-15 minutes | 2-4 GB |
| Medium (200-500 subjects) | 15-45 minutes | 4-8 GB |
| Large (500+ subjects) | 45-120 minutes | 8-16 GB |

## 📊 Progress Monitoring

The script provides detailed progress information:

```
Loading required libraries...
Setting up configuration...
Defining utility functions...

=== Step 1: Processing Phenotype Data ===
=== Data Quality Check for Raw Phenotype Data ===
Dimensions: 1112 x 74
Missing values: 15234
Complete cases: 568

=== Step 2: Processing Neuroimaging Data and Correlation Analysis ===
Performing correlation analysis...

=== Step 3: Merging Phenotype and Neuroimaging Data ===
=== Data Quality Check for Merged Data ===
Dimensions: 486 x 76652

... and so on for all 9 steps
```

## 📁 Output Verification

After successful execution, verify these outputs:

### Directory Structure
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
│   └── *.png files
├── 06_mancova/
│   └── *.csv files
├── 07_svm/
│   └── *.csv files
├── varimax_alldata.csv
├── promax_alldata.csv
└── ANALYSIS_SUMMARY.txt
```

### Key Output Files
- **ANALYSIS_SUMMARY.txt**: Comprehensive analysis report
- **varimax_alldata.csv**: Final dataset with Varimax PC scores
- **promax_alldata.csv**: Final dataset with Promax PC scores

## 🛠️ Troubleshooting

### Common Issues and Solutions

#### Issue 1: Package Installation Errors
```r
# If automatic installation fails, install manually:
install.packages(c("dplyr", "data.table", "ggplot2", "caret", "e1071"))
```

#### Issue 2: Memory Issues
```r
# For large datasets, increase memory limit (Windows):
memory.limit(size = 8000)  # 8GB

# Or reduce the number of PCs:
# Edit line in script: top_pcs <- 10  # instead of 20
```

#### Issue 3: Permission Errors
```bash
# Ensure write permissions
chmod 755 .
# Run R as administrator (Windows) or with sudo (Linux)
```

#### Issue 4: Missing Data Files
```
Error: No phenotype data available. Please ensure Phenotypic_V1_0b_v1.csv is in the working directory.
```
**Solution**: Download and place the required data files in your working directory.

#### Issue 5: Correlation Analysis Fails
```
Warning: No neuroimaging data found. Skipping neuroimaging analysis.
```
**Solution**: This is normal if you only have phenotypic data. The script will continue with available data.

## 🔍 Validation Checks

### 1. Data Integrity
```r
# Check phenotype data
phenotype_data <- read.csv("01_phenotype_data/phenotype_data.csv")
summary(phenotype_data)

# Check final datasets
varimax_data <- read.csv("varimax_alldata.csv")
dim(varimax_data)  # Should have subjects x (demographic + PC columns)
```

### 2. Results Validation
```r
# Check analysis summary
cat(readLines("ANALYSIS_SUMMARY.txt"), sep = "\n")

# Check SVM performance
svm_results <- read.csv("07_svm/svm_performance_results.csv")
print(svm_results)
```

## 📈 Performance Optimization

### For Large Datasets:
1. **Use SSD storage** for faster I/O
2. **Increase RAM** (16GB+ recommended)
3. **Close other applications** during analysis
4. **Use parallel processing** (already enabled in script)

### For Limited Resources:
1. **Reduce top_pcs** parameter (e.g., from 20 to 10)
2. **Process in batches** if memory issues persist
3. **Use server/cloud computing** for very large datasets

## 📧 Support

### Getting Help:
1. **Check this guide** for common solutions
2. **Review ANALYSIS_SUMMARY.txt** for execution details
3. **Check R console output** for specific error messages
4. **Verify data file formats** match expected structure

### Error Reporting:
When reporting issues, include:
- R version: `R.version.string`
- Operating system
- Error messages (full text)
- Data file sizes and formats
- Available system memory

## ✅ Success Indicators

Your analysis completed successfully if you see:

```
=== Analysis Complete ===
Results have been saved to respective output directories.
Summary report available in: ANALYSIS_SUMMARY.txt
🎉 ABIDE Comprehensive Analysis Pipeline completed successfully! 🎉
```

And the following files exist:
- `ANALYSIS_SUMMARY.txt`
- `varimax_alldata.csv`
- `promax_alldata.csv`
- Output directories (01-07) with analysis results

---

**Congratulations!** Your ABIDE autism research analysis is complete. Review the `ANALYSIS_SUMMARY.txt` file for detailed results and proceed with your research interpretation.