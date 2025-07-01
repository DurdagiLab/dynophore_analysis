# Dynamic e-Pharmacophore Analysis Tool


# Description

This Python script provides an automated and robust pipeline for the analysis of dynamic pharmacophore hypotheses derived from molecular dynamics (MD) simulations performed using Desmond (Schrödinger Suite). Tailored for structure-based drug design workflows, this tool enables temporal and structural profiling of pharmacophoric feature evolution throughout the simulation. It is designed to be fully compatible with the output formats of Schrödinger’s Desmond engine and integrates seamlessly with pharmacophore hypotheses generated via Phase.

**Important: Prior to hypothesis analysis, the trajectory must be structurally aligned in VMD using either the `protein backbone` or `Cα-atoms` with respect to the average structure. Following alignment, the RMSD profile is computed, and the resulting `trajrmsd.dat` file is used for structural correlation within this script.


# Key Features

- Parsing and standardization of pharmacophoric annotations from frame-specific feature tables;
- Sequence-preserving abstraction of features by removing numeric labels;
- Mapping of trajectory frames to pharmacophore hypotheses in conjunction with structural deviation metrics (RMSD);
- Quantitative assessment of hypothesis frequency and relative prevalence across the trajectory;
- Identification of the structurally most representative frame (i.e., with the lowest RMSD) for each unique hypothesis;
- Generation of high-resolution, publication-quality visualizations (pie chart, bar plot, and stacked bar plot)
  illustrating dominant pharmacophore patterns;
- Export of an interactive HTML summary and a printable PDF report of all computed statistics;
- Aggregation of all output files in a dedicated results directory (`DYNOPHORE_RESULTS`);
- Extraction of the top three most frequently observed hypotheses (with length ≥ 4) and export of their corresponding
  `.phypo` representations into a `BEST_HYPOTHESES` folder.

# Requirements

- Python 3.7+




- External dependencies (install via pip):

pip install -r requirements.txt
