import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

rna = pd.read_csv("rnaseq.csv")

plt.figure(figsize=(10, 6))

sns.histplot(data=rna, x='expression', bins=30, kde=True, color='skyblue')
plt.xlabel('Expression')
plt.title('Expression Histogram')
plt.tight_layout()
plt.show()

plt.figure(figsize=(10, 6))

# Mutate expression_log
rna['expression_log'] = np.log2(rna['expression'] + 1)

sns.histplot(data=rna, x='expression_log', bins=30, kde=True, color='lightcoral')
plt.xlabel('Expression (log2)')
plt.title('Expression Log Histogram')
plt.tight_layout()
plt.show()

# Perform gene expression analysis
rna_fc = rna.groupby(['gene', 'time', 'gene_biotype']).agg(mean_exp=('expression_log', 'mean')).reset_index()
rna_fc = rna_fc.pivot(index=['gene', 'gene_biotype'], columns='time', values='mean_exp').reset_index()
rna_fc['time_8_vs_0'] = rna_fc[8] - rna_fc[0]
rna_fc['time_4_vs_0'] = rna_fc[4] - rna_fc[0]

plt.figure(figsize=(12, 8))

sns.scatterplot(data=rna_fc, x='time_4_vs_0', y='time_8_vs_0', alpha=0.6, hue='gene_biotype', palette='Set2')
plt.xlabel('Change in Expression from Time 0 to Time 4')
plt.ylabel('Change in Expression from Time 0 to Time 8')
plt.title('Gene Expression Change Comparison')
plt.tight_layout()
plt.show()

plt.figure(figsize=(12, 6))

sns.boxplot(data=rna, x='sample', y='expression_log', color='lightgreen', showmeans=True)
sns.stripplot(data=rna, x='sample', y='expression_log', color='black', alpha=0.2)
plt.xlabel('Sample')
plt.ylabel('Expression (log2)')
plt.title('Gene Expression by Sample')
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()

# Select top 10 genes
genes_selected = rna_fc.sort_values(by='time_8_vs_0', ascending=False).head(10)['gene']

sub_rna = rna[rna['gene'].isin(genes_selected)]

mean_exp_by_time = sub_rna.groupby(['gene', 'time'])['expression_log'].mean().reset_index()

plt.figure(figsize=(8, 6))

sns.lineplot(data=mean_exp_by_time, x='time', y='expression_log', hue='gene', palette='tab10')
plt.xlabel('Time (Duration of infection)')
plt.ylabel('Mean gene Expression (log2)')
plt.title('Mean Expression by Time')
plt.tight_layout()
plt.show()

mean_exp_by_time_sex = sub_rna.groupby(['gene', 'time', 'sex'])['expression_log'].mean().reset_index()


plt.clf()

plt.figure(figsize=(12, 8))

genes_selected = rna.groupby('gene')['expression_log'].mean().nlargest(10).index
sub_rna = rna[rna['gene'].isin(genes_selected)]

# Group by gene, time, and sex and calculate mean expression
mean_exp_by_time_sex = sub_rna.groupby(['gene', 'time', 'sex'])['expression_log'].mean().reset_index()

for i, gene in enumerate(genes_selected, 1):
    plt.subplot(3, 4, i)
    gene_data = mean_exp_by_time_sex[mean_exp_by_time_sex['gene'] == gene]
    sns.lineplot(data=gene_data, x='time', y='expression_log', hue='sex', style='sex', markers=True)
    plt.title(gene)
    plt.xlabel('Time (Duration of infection)')
    plt.ylabel('Mean gene Expression (log2)')
    plt.tight_layout()

plt.show()

