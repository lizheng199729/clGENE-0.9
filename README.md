# clGENE a framework designed to integrate multi-omics data based on clinical phenotypic information
![幻灯片1](https://github.com/lizheng199729/clGENE/assets/138444869/283902d1-f83c-46e3-a725-d13fc8b20639)


## Introduction

`clGENE` is an innovative R package designed for the comprehensive analysis of various omics data. It facilitates the integrated analysis of multi-omics datasets, aiming to assist researchers in exploring the potential molecular mechanisms that underpin similar phenotypes across different diseases.

## Features

- **Standardized Input**: Accepts a standardized expression matrix for genes, proteins, or metabolites.
- **Phenotype Analysis**: Leverages clinical phenotype data to pinpoint gene expressions associated with specific clinical conditions.
- **Advanced Algorithms**: Implements PCA algorithms and Euclidean distance measures to filter gene clusters.
- **Visualization Tools**: Includes three distinct methods for visualizing data and results.
- **Open Source**: Provides full access to its source code for extended flexibility and customization.

## Getting Started

### Installation

To install the latest version of `clGENE`, you can use the following command in R:

```R
# install.packages("devtools")
devtools::install_github("lizheng199729/clGENE")
# Load clGENE
library(clGENE)

# library clGENE needs to use the package
library(dplyr)
library(Matrix)
library(NMF)
library(pheatmap)
# Reading sample data
data(package = "clGENE")
# Using the orthogonal algorithm integrating data
con_data<- Orthogonal(con_Raw_data)
stage1_data<- Orthogonal(stage1_Raw_data)
common_rows <- intersect(rownames(con_data), rownames(stage1_data))
# Merge data frames using cbind
exp <- cbind(
  con_data[common_rows, ],
  stage1_data[common_rows, ]
)
# Releasing the group file
group1 <- group$group
new_matrix <- rbind(group$group, exp)
rownames(new_matrix)[1] <- "group"
# Using the second step, the genes most relevant to the clinical data were calculated
hubGENE2 <-TOPGENE(new_matrix, threshold = 0.5, genenumber = 20)
# Use the third step, screening key molecular processes associated with clinical information
Program2 <- Program(hubGENE2[[1]],Clustering = 10)

