🧬 CRC Transcriptomic Biomarker Discovery
From RNA-seq Data to Predictive Models for Precision Oncology

🚀 Overview
This project presents an end-to-end transcriptomic analysis pipeline designed to identify robust biomarkers that distinguish tumor vs. non-tumor tissues in colorectal cancer (CRC).
It integrates:
📊 Statistical modeling (DESeq2)
🧠 Network biology (WGCNA)
🤖 Machine learning (Random Forest)
👉 The final outcome is a biomarker panel with strong predictive power (AUC ~0.97), demonstrating its potential for clinical stratification and drug development.

🧠 Business Case
🎯 Problem
In oncology and drug development, a key challenge is:
How to reliably distinguish tumor from non-tumor tissue using molecular data?
This is critical for:
Patient stratification
Treatment selection
Biomarker-driven clinical trials
Target discovery

💡 Solution
This project builds a robust transcriptomic signature by combining:
Differential expression (DEGs)
Co-expression networks (WGCNA)
Hub gene prioritization
Machine learning validation

🧪 Value for Industry
✔ Identify clinically relevant biomarkers (FDR < 0.05)
✔ Discover biologically coherent gene modules
✔ Prioritize high-impact targets for drug development
✔ Build predictive models for precision medicine

📂 Dataset
Source: GEO (GSE50760 – synthetic subset)
Samples: 54 total
18 Normal
18 Tumor
18 Metastasis
Features: ~10,000 genes
Design: paired samples (patients) → higher statistical power

⚙️ Pipeline Overview
1️⃣ Data Integration & QC
Count matrix + metadata alignment
Sample consistency checks
Filtering low-expression genes
2️⃣ Normalization
DESeq2 normalization
Variance stabilizing transformation (VST)
3️⃣ Exploratory Data Analysis (EDA)
PCA → clear separation of conditions
Distribution analysis
Sample clustering
4️⃣ Differential Expression Analysis
Tumor vs Normal comparison
Filtering:
FDR < 0.05
|log2FC| > 1
📈 Output:
Volcano plot
DEG tables
5️⃣ Network Analysis (WGCNA)
Gene co-expression modules
Module–trait correlation heatmap
🔥 Key insight:
Modules strongly associated with:
Tumor
Normal
Metastasis
6️⃣ Hub Gene Identification
Custom metric:
Hub Score = |kME| × |log2FC| × (−log10(padj))
Where:
kME → module connectivity (network importance)
log2FC → effect size
padj → statistical significance
👉 Result: biologically meaningful "leader genes"
7️⃣ Biomarker Panel Construction
Filtering criteria:
padj < 0.01
|log2FC| > 1
|kME| > 0.7
📌 Output:
Top 20–30 biomarkers
8️⃣ Machine Learning Validation
Model: Random Forest
Feature set: selected biomarkers
Evaluation:
OOB error
ROC curve

📊 Key Results
🧬 Biological Insights
Clear transcriptomic separation between:
Normal vs Tumor vs Metastasis
Identification of coordinated gene programs (modules)
Functional relevance of biomarkers (e.g. ECM remodeling, ion transport, immune response)

🤖 Predictive Performance
AUC: ~0.975
Strong classification of tumor vs normal
Stable model (OOB error converges early)

🧪 Example Biomarkers
COL10A1
MMP1
BEST4
FOXQ1
SPIB
OTOP family

👉 These genes are linked to:
Tumor progression
Microenvironment remodeling
Cellular differentiation

📈 Visualizations Included
PCA plot
Volcano plot
Violin plots (top biomarkers)
DEG heatmap
WGCNA module heatmap
Hub gene ranking
Random Forest importance plot
ROC curve

🧬 Why This Matters
This project goes beyond standard DEG analysis:
Approach	Insight Level
DEGs only	Basic
DEGs + WGCNA + ML	🚀 Industry-grade / publication-level

🧠 Key Takeaway
This pipeline enables the discovery of robust, biologically meaningful and clinically actionable transcriptomic biomarkers, ready to be translated into:
Diagnostic tools
Predictive models
Drug targets

🛠 Tech Stack
R
DESeq2
WGCNA
ComplexHeatmap
ggplot2
Python
scikit-learn
pandas
matplotlib / seaborn
Environment
Google Colab

📌 Future Applications
Integration with:
Genomics (mutations)
Clinical variables
Drug response data
Expansion to:
Multi-omics
Survival prediction
Therapy response modeling

👨‍🔬 Author
David Villafañe
Biotechnologist | PhD in Biological Sciences
Specializing in Genomics, Bioinformatics & Translational Medicine
Linkedin: https://www.linkedin.com/in/davidvillafanie/

💬 Final Statement
“Transforming transcriptomic data into actionable precision medicine insights.”

Notion :https://www.notion.so/Transcriptomic-Biomarker-Discovery-Predictive-Modeling-in-Colorectal-Cancer-351e81ae11b7802db598cf2b23b7c5ab
