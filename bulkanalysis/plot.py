#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
# SPDX-FileCopyrightText: Copyright © 2026 University of Bern
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@unibe.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

def plot_gene_expression_distribution(df_genes, title='Gene expression', save_file=None):
    """Histogram of log2(counts + 1) of every sample (one step line per column).

    Args:
        df_genes: gene count matrix, samples in columns.
        title: plot title. Defaults to "Gene expression".
        save_file: if given, path where the figure is saved.
    """
    df_genes_log = np.log2(df_genes+1)
    plt.figure()
    sns.histplot(df_genes_log, stat='proportion', alpha=0.05, element='step')
    plt.legend('',frameon=False)
    plt.title(title, weight='bold')
    plt.xlabel("log2(read counts + 1)")
    if save_file is not None:
        plt.savefig(save_file, bbox_inches='tight')

def boxplot_pca(df_pca, treatment_col=None, treatment_order=None, PCs=None, palette=None):
    """Box plots of the principal component values per treatment, one panel per component.

    Args:
        df_pca: one row per sample, with `PC1`, `PC2`, ... columns and a treatment column.
        treatment_col: treatment column.
        treatment_order: order of the treatments on the x-axis (default: order of appearance).
        PCs: 0-based indices of the components to plot (default: every PC column).
        palette: colours per treatment.
    """
    if PCs is None:
        PCs = range(len([col for col in df_pca.columns if col.startswith('PC')]))
    
    nb_pcs = len(PCs)

    fig, ax = plt.subplots(1, nb_pcs, figsize=(2*nb_pcs, 5))
    dfs=[]
    if treatment_order is not None:
        for treatment in treatment_order:
            df = df_pca[df_pca[treatment_col] == treatment]
            dfs.append(df)
        df_pca = pd.concat(dfs)
    for ax_nb, i in enumerate(PCs):
        sns.boxplot(data=df_pca, x=treatment_col, y='PC'+str(i+1), ax=ax[ax_nb], palette=palette)
        ax[ax_nb].set_xticks(ticks=list(range(len(df_pca[treatment_col].value_counts(sort=False).index))), labels=df_pca[treatment_col].value_counts(sort=False).index, rotation=45)
        ax[ax_nb].set_title('PC'+str(i+1), weight='bold')
    plt.suptitle(f"PC values distribution across samples depending on {treatment_col}", weight='bold')
    plt.tight_layout()