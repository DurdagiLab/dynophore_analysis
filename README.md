# Dynamic e-Pharmacophore Analysis Tool

# Description
This Python-based tool offers an automated and comprehensive pipeline for the analysis of dynamic pharmacophore hypotheses derived from molecular dynamics (MD) simulations performed using Desmond (Schrödinger Suite). Specifically developed for structure-based drug design (SBDD) workflows, the tool enables the temporal and structural assessment of pharmacophoric feature evolution throughout an MD trajectory.

It is fully compatible with pharmacophore feature tables generated via Phase and integrates RMSD-based structural correlation through `trajrmsd.dat` files obtained from VMD-aligned trajectories.

Note: Prior to running this analysis, the MD trajectory must be structurally aligned to the average structure using protein  `backbone` or  `Cα-atoms` in VMD. The corresponding RMSD profile (`trajrmsd.dat `) must be generated and placed in the working directory.

# Key Functionalities
- `Pharmacophore Feature Parsing` - Automated extraction and standardization of features from frame-specific tables;
- `Abstraction of Hypotheses` - Removal of numerical labels to generate sequence-invariant feature representations;
- `Structural Mapping` - Integration of frame-wise pharmacophore patterns with structural deviation (RMSD) metrics;
- `Quantitative Profiling` - Calculation of frequency, percentage occurrence, and relative prominence of pharmacophore hypotheses across the trajectory;
- `Structural Representativeness` - Identification of the most structurally representative frame (lowest RMSD) for each hypothesis;
- `Visual Analytics` - Generation of high-resolution pie charts, bar plots, and stacked bar plots for dominant patterns;
- `Reporting` - Export of an interactive HTML summary and a printable PDF report;
- `Automated Output Management` - Aggregation of all outputs in a dedicated results folder (`DYNOPHORE_RESULTS`);
- `Lead Hypotheses Extraction`: Retrieval of the top three most frequently observed hypotheses (length ≥ 4) with export of their `.phypo` representations to the `BEST_HYPOTHESES` directory.

# System Requirements
Python Version: 3.7 or higher


# Required Python Packages
All dependencies are listed in the requirements.txt file:
- `pandas>=1.1` 
- `numpy>=1.18` 
- `matplotlib>=3.3`
- `seaborn>=0.11`
- `weasyprint>=53.0`

# Installation Instructions:
We strongly recommend using a virtual environment (e.g., venv or conda) for package isolation. Once activated, install the dependencies via:

> pip install -r requirements.txt

***Note: On some systems, weasyprint may require additional system-level libraries such as cairo, pango, or gdk-pixbuf. For installation help, consult the WeasyPrint installation guide.

# Input Requirements
- trajrmsd.dat: RMSD profile of the aligned MD trajectory (generated via VMD
- <working_directory>/DYNOPHORE_ANALYSIS/PROCESSED_FILES/: Directory containing pharmacophore feature tables (CSV format) named as <frame>_hypo_features_table.csv
- <working_directory>/DYNOPHORE_ANALYSIS/saved_HYPOTHESIS/: Directory containing .phypo pharmacophore files for corresponding frames

# Output
- DYNOPHORE_RESULTS/ folder containing:
  - PDF and HTML reports
  - All generated figures
  - BEST_HYPOTHESES/ folder with .phypo files of the top 3 most frequent hypotheses (length ≥ 4)

# Citation
If you utilize this tool in your research or publication, please cite the following institution:

Computational Drug Design Center (HITMER), Faculty of Pharmacy, Bahçeşehir University, Istanbul, Turkey
