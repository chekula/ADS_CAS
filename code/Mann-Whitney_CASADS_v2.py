#pip3 install pandas matplotlib seaborn numpy scipy

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import mannwhitneyu, ttest_ind, probplot
from matplotlib.gridspec import GridSpec

# print current directory
import os
print(os.getcwd())

# Read and process the data
hl = pd.read_csv('../data/HL.csv')
hl = hl[hl['group'].isin(['Neurite', 'Soma'])]
hl_pivoted = hl.pivot(index='gene_id', columns='group', values='Half_life')
hl_pivoted.columns = ['hl.Neurite', 'hl.Soma']

# Set up the plotting area
fig = plt.figure(figsize=(8, 10))
gs = GridSpec(5, 1, height_ratios=[1, 1, 0.1, 1, 1])

# Plot density for hl.Soma with mean and median lines and labels
ax1 = fig.add_subplot(gs[0])
sns.kdeplot(hl_pivoted['hl.Soma'].dropna(), fill=True, color='blue', alpha=0.5, ax=ax1)
ax1.axvline(hl_pivoted['hl.Soma'].mean(), color='red', linestyle='--', linewidth=1)
ax1.axvline(hl_pivoted['hl.Soma'].median(), color='green', linestyle='--', linewidth=1)
ax1.text(hl_pivoted['hl.Soma'].mean() + 0.5, 0.05, 'Mean', color='red', rotation=90, va='center')
ax1.text(hl_pivoted['hl.Soma'].median() + 0.5, 0.05, 'Median', color='green', rotation=90, va='center')
ax1.set_title('Density Plot of hl.Soma')
ax1.set_xlabel('hl.Soma')
ax1.set_ylabel('Density')

# Plot density for hl.Neurite with mean and median lines and labels
ax2 = fig.add_subplot(gs[1])
sns.kdeplot(hl_pivoted['hl.Neurite'].dropna(), fill=True, color='purple', alpha=0.5, ax=ax2)
ax2.axvline(hl_pivoted['hl.Neurite'].mean(), color='red', linestyle='--', linewidth=1)
ax2.axvline(hl_pivoted['hl.Neurite'].median(), color='green', linestyle='--', linewidth=1)
ax2.text(hl_pivoted['hl.Neurite'].mean() + 0.5, 0.05, 'Mean', color='red', rotation=90, va='center')
ax2.text(hl_pivoted['hl.Neurite'].median() + 0.5, 0.05, 'Median', color='green', rotation=90, va='center')
ax2.set_title('Density Plot of hl.Neurite')
ax2.set_xlabel('hl.Neurite')
ax2.set_ylabel('Density')

# Perform Mann-Whitney U Test
mann_whitney_result = mannwhitneyu(hl_pivoted['hl.Soma'].dropna(), hl_pivoted['hl.Neurite'].dropna(), alternative='two-sided')

# Perform Welch's t-test
welch_test_result = ttest_ind(hl_pivoted['hl.Soma'].dropna(), hl_pivoted['hl.Neurite'].dropna(), equal_var=False)

# Extract p-value from the Mann-Whitney U test result
p_value = mann_whitney_result.pvalue
p_value_text = f"Mann-Whitney p-value: {p_value:.10f}"

# Create a text annotation for the p-value
ax3 = fig.add_subplot(gs[2])
ax3.text(0.5, 0.5, p_value_text, fontsize=12, ha='center', va='center')
ax3.axis('off')

# Q-Q plot for hl.Soma
ax4 = fig.add_subplot(gs[3])
probplot(hl_pivoted['hl.Soma'].dropna(), dist="norm", plot=ax4)
ax4.set_title("Q-Q Plot of hl.Soma")

# Q-Q plot for hl.Neurite
ax5 = fig.add_subplot(gs[4])
probplot(hl_pivoted['hl.Neurite'].dropna(), dist="norm", plot=ax5)
ax5.set_title("Q-Q Plot of hl.Neurite")

# Save the figure
plt.savefig("1C_delta.hl_histo_forCAS.ADS.pdf", bbox_inches='tight')

# Show the plot
plt.show()

# Print the results of the tests
print("Mann-Whitney U Test result:", mann_whitney_result)
print("Welch's t-test result:", welch_test_result)