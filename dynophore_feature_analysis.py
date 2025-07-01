"""
###############################################################################################################################
Title: Dynamic e-Pharmacophore Analysis Tool

Author: Dr. Mine Isaoglu  
Principal Investigator: Prof. Dr. Serdar Durdagi
Affiliation: Computational Drug Design Center (HITMER), Faculty of Pharmacy, Bahçeşehir University, Istanbul, Turkey  
             
Version: May 2025
###############################################################################################################################
"""

import os
import re
import shutil
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from weasyprint import HTML

# === CONFIGURATION PARAMETERS === #
top_n = 10
exclude_lengths = [1, 2, 3]  # Feature lengths to exclude from summary and plots

# === DIRECTORY PATHS === #
base_dir = os.getcwd()
rmsd_path = os.path.join(base_dir, 'trajrmsd.dat')
csv_dir = os.path.join(base_dir, 'DYNOPHORE_ANALYSIS', 'PROCESSED_FILES')

# === COLOR PALETTE SELECTION === #
color_palettes = [
    "Set1", "Set2", "Set3", "Pastel1", "Pastel2", "tab10", "tab20",
    "Paired", "Dark2", "Accent", "coolwarm", "Spectral", "YlGnBu"
]
print("Available color palettes:")
for i, cmap in enumerate(color_palettes, 1):
    print(f"{i}. {cmap}")
try:
    cmap_index = int(input("\nEnter the number of the colormap you want to use (default: Set2): "))
    palette_name = color_palettes[cmap_index - 1] if 1 <= cmap_index <= len(color_palettes) else "Set2"
except:
    palette_name = "Set2"
print(f"Selected colormap: {palette_name}")

# === FEATURE TABLE PARSING === #
csv_files = sorted([
    os.path.join(csv_dir, f) for f in os.listdir(csv_dir)
    if re.match(r'\d+_hypo_features_table\.csv$', f)
])

frame_feature_map = {}
for csv_path in csv_files:
    match = re.search(r'(\d+)_hypo_features_table\.csv', os.path.basename(csv_path))
    if match:
        frame_id = int(match.group(1))
        try:
            with open(csv_path, 'r') as f:
                header = f.readline()
            delimiter = ';' if header.count(';') > header.count(',') else ','
            df = pd.read_csv(csv_path, delimiter=delimiter)
            feature_col = 'Feature_label' if 'Feature_label' in df.columns else df.columns[1]
            features_raw = df[feature_col].astype(str).tolist()
            features_clean = ''.join([re.sub(r'\d+', '', feat) for feat in features_raw])
            frame_feature_map[frame_id] = features_clean
        except Exception as e:
            print(f"Error reading {csv_path}: {e}")

# === RMSD FILE PARSING === #
rmsd_df = pd.read_csv(rmsd_path, sep=r'\s+', comment='#', header=None, names=["Frame", "RMSD"])
rmsd_df = rmsd_df.iloc[1:].copy()
rmsd_df["Frame"] = rmsd_df["Frame"].astype(int)
rmsd_df["RMSD"] = rmsd_df["RMSD"].astype(float)
rmsd_df["Frame_ID"] = rmsd_df["Frame"] + 1
rmsd_df = rmsd_df[rmsd_df["Frame_ID"] != 5002]
rmsd_df['Features'] = rmsd_df['Frame_ID'].map(frame_feature_map).fillna('-')

# === HYPOTHESIS STATISTICS === #
total_frames = len(rmsd_df)
feature_counts = rmsd_df[rmsd_df['Features'] != '-']['Features'].value_counts().reset_index()
feature_counts.columns = ['Features', 'Count']
feature_counts['Percent'] = 100 * feature_counts['Count'] / total_frames
feature_counts['Length'] = feature_counts['Features'].apply(len)

if feature_counts.empty:
    print("No matched pharmacophore hypotheses found in frames.")
    exit()

summary_records = []
for _, row in feature_counts.iterrows():
    feat = row['Features']
    subset = rmsd_df[rmsd_df['Features'] == feat]
    min_row = subset.loc[subset['RMSD'].idxmin()]
    summary_records.append({
        'Features': feat,
        'Lowest_RMSD': round(min_row['RMSD'], 3),
        'Frame': int(min_row['Frame_ID']),
        'Count': int(row['Count']),
        'Length': len(feat),
        'Percent': round(row['Percent'], 1)
    })
summary_df = pd.DataFrame(summary_records).sort_values(by='Count', ascending=False)

# === HTML REPORT GENERATION === #
summary_lines = []
for length in sorted(summary_df['Length'].unique()):
    if length in exclude_lengths:
        continue
    group = summary_df[summary_df['Length'] == length]
    summary_lines.append(f"<h3>{length}-Feature Combinations</h3><ul>")
    for _, row in group.iterrows():
        summary_lines.append(
            f"<li><b>{row['Features']}</b>: {row['Count']} frames ({row['Percent']}%), "
            f"Lowest RMSD = {row['Lowest_RMSD']} at frame {row['Frame']}.</li>"
        )
    summary_lines.append("</ul>")

mapping_table_html = rmsd_df[['Frame_ID', 'RMSD', 'Features']].to_html(index=False)
html_summary = (
    "<h2>Generated Pharmacophore Hypotheses</h2>" +
    ''.join(summary_lines) +
    "<h2>Frame-wise Hypothesis Mapping</h2>" +
    mapping_table_html
)

html_path = os.path.join(base_dir, 'Feature_Summary_Report.html')
final_report = os.path.join(base_dir, 'Feature_Analysis_Report.pdf')
with open(html_path, 'w') as f:
    f.write(html_summary)
HTML(html_path).write_pdf(final_report)

# === VISUALIZATIONS (EXCLUDING LENGTHS) === #
plot_df = summary_df[~summary_df['Length'].isin(exclude_lengths)].nlargest(top_n, 'Percent').copy()
labels = [f"{f} ({p:.1f}%)" for f, p in zip(plot_df['Features'], plot_df['Percent'])]
colors = sns.color_palette(palette_name, len(plot_df))

# Pie chart
plt.figure(figsize=(9, 9))
plt.pie(plot_df['Percent'], labels=labels, startangle=140, colors=colors,
        wedgeprops=dict(width=0.5), textprops={'fontsize': 11})
plt.gca().add_artist(plt.Circle((0, 0), 0.3, color='white'))
plt.title(f"Top {top_n} Pharmacophore Hypotheses", fontsize=15, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(base_dir, 'Feature_PieChart.png'), dpi=300)
plt.close()

# Bar chart
plt.figure(figsize=(10, 5))
ax = sns.barplot(data=plot_df, x='Features', y='Percent', hue='Features', palette=palette_name)
if ax.legend_ is not None:
    ax.legend_.remove()
plt.title(f"Top {top_n} Hypotheses by Frequency", fontsize=14, fontweight='bold')
plt.xlabel("Hypothesis (Feature Sequence)", fontsize=12, fontweight='bold')
plt.ylabel("Percentage (%)", fontsize=12, fontweight='bold')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig(os.path.join(base_dir, 'Feature_BarChart.png'), dpi=300)
plt.close()

# Stacked bar chart
plot_df['Length'] = plot_df['Features'].apply(len)
length_grouped = plot_df.groupby(['Length', 'Features'])['Percent'].sum().unstack().fillna(0)
length_grouped.T.plot(kind='bar', stacked=True, colormap=palette_name, figsize=(10, 6))
plt.title(f"Stacked Bar Chart of Top {top_n} Hypotheses by Length", fontsize=14, fontweight='bold')
plt.xlabel("Feature Combination", fontsize=12, fontweight='bold')
plt.ylabel("Percentage (%)", fontsize=12, fontweight='bold')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig(os.path.join(base_dir, 'Feature_StackedBarChart.png'), dpi=300)
plt.close()

# === ORGANIZE OUTPUTS === #
results_dir = os.path.join(base_dir, 'DYNOPHORE_RESULTS')
os.makedirs(results_dir, exist_ok=True)

result_files = [
    'Feature_Analysis_Report.pdf',
    'Feature_PieChart.png',
    'Feature_BarChart.png',
    'Feature_StackedBarChart.png',
    'Feature_Summary_Report.html'
]
for file in result_files:
    src = os.path.join(base_dir, file)
    dst = os.path.join(results_dir, file)
    if os.path.exists(src):
        shutil.move(src, dst)

# === COPY BEST HYPOTHESES WITH LENGTH >= 4 === #
best_hypo_dir = os.path.join(base_dir, 'BEST_HYPOTHESES')
os.makedirs(best_hypo_dir, exist_ok=True)
saved_hypo_dir = os.path.join(base_dir, 'DYNOPHORE_ANALYSIS', 'saved_HYPOTHESIS')

filtered_summary = summary_df[summary_df['Length'] >= 4].head(3)
for _, row in filtered_summary.iterrows():
    frame_id = row['Frame']
    hypo_file = f"{frame_id}_hypo.phypo"
    src_path = os.path.join(saved_hypo_dir, hypo_file)
    dst_path = os.path.join(best_hypo_dir, hypo_file)
    if os.path.exists(src_path):
        shutil.copy2(src_path, dst_path)
    else:
        print(f"Warning: {hypo_file} not found in saved_HYPOTHESIS.")

print("Analysis complete.")
print(f"Final PDF report saved as: {os.path.join(results_dir, 'Feature_Analysis_Report.pdf')}")
