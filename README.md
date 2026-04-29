🧬 CRC Transcriptomic Biomarker Discovery Pipeline
This project implements an end-to-end RNA-seq analysis pipeline to identify transcriptional biomarkers in colorectal cancer (CRC), integrating:
Differential Expression (DESeq2)
Co-expression Network Analysis (WGCNA)
Biomarker prioritization
Machine Learning modeling

🎯 Objective
To identify robust transcriptomic biomarkers from:
Normal tissue
Primary tumor
Metastatic samples
and integrate them with clinical/genomic features to improve:
Drug response prediction
Treatment stratification

🧪 Pipeline Overview
Data loading (RNA-seq counts + metadata)
Gene annotation (Entrez → Symbol)
Quality control & normalization
Differential expression analysis
Biomarker prioritization
Network analysis (WGCNA)
Hub gene identification
Final biomarker panel selection
Heatmap visualization
ML-ready dataset generation

🧠 Key Features
Multi-condition analysis (Normal vs Tumor vs Metastasis)
Integrated scoring system for biomarker prioritization
Network-aware biomarker selection (WGCNA + DEGs)
ML-ready outputs
Fully reproducible (Colab-compatible)

📊 Outputs
Differential expression tables
Volcano plots
PCA plots
Heatmaps
WGCNA module-trait relationships
Final biomarker panel
ML-ready feature matrix

🚀 Technologies
Python (pandas, mygene)
R (DESeq2, WGCNA, ComplexHeatmap)
Google Colab

🧬 Example Results

8000 DEGs identified
Network modules strongly associated with tumor state
High-confidence biomarker panel (n=30)
Clear separation of biological conditions (PCA, heatmaps)

📈 Future Work
Integration with clinical variables
Drug response modeling
Deep learning approaches
External dataset validation (TCGA)

👨‍🔬 Author
David Villafañe
Biotechnologist & PhD in Biological Sciences
Focused on Genomics, Bioinformatics & Precision Medicine
