#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
# SPDX-FileCopyrightText: Copyright © 2026 University of Bern
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@unibe.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import combinations
from statannotations.Annotator import Annotator
from scipy.stats import ttest_ind
import csv



def load_DaPars_outputs(path: str):
    """Concatenate the DaPars result tables (.txt) of every chromosome folder (name containing "chr") of `path`; column
    names are cut to their last path component.

    Args:
        path: DaPars output directory.

    Returns:
        pd.DataFrame: the DaPars results of all chromosomes.
    """
    # Load DaPars output
    dfs = []
    for dir in os.listdir(path):
        if 'chr' in dir:
            for file in os.listdir(os.path.join(path, dir)):
                # df = pd.read_csv(os.path.join(dir))
                if file.endswith('.txt'):
                    dfs.append(pd.read_csv(os.path.join(path, dir, file), sep='\t'))

    df = pd.concat(dfs)
    df.columns = [col.split('/')[-1] for col in df.columns]
    return df

def PDUI_boxplots_distributions(df: pd.DataFrame, treatments_dict: dict, fs: int = 17, **kwargs):
    """Box plots of the PDUI (DaPars) per treatment, with the Welch's t-test p-value of every pair of treatments.

    A `<treatment>_PDUI` column (mean of the `<sample>_PDUI` columns of the treatment) is added to `df`.

    Args:
        df: DaPars results with one `<sample>_PDUI` column per sample.
        treatments_dict: {treatment: sample names}.
        fs: font size of the p-values. Defaults to 17.
        **kwargs: passed to seaborn.boxplot.
    """
    ## Mean of the DPUI for the replicates
    for treatment, ids in treatments_dict.items():    
        df[f"{treatment}_PDUI"] = df[[f"{id}_PDUI" for id in ids]].mean(axis=1)

    fig, ax = plt.subplots(1, 1, figsize=(2, 2.5*len(treatments_dict)))
    # create stats annotations
    pairs = list(combinations([f"{treatment}_PDUI" for treatment in treatments_dict], 2))

    # compute p-values
    p_values = []
    for p1, p2 in pairs:
        p_values.append(ttest_ind(df[p1], 
                        df[p2], 
                        equal_var=False, 
                        alternative='two-sided',
                        nan_policy='omit')[1])
    formatted_pvalues = [f'P={pvalue:.2e}' for pvalue in p_values]

    plotting_parameters = {'data': df[[f"{treatment}_PDUI" for treatment in treatments_dict.keys()]],
                        'linewidth': 2.3, **kwargs}

    sns.boxplot(**plotting_parameters, ax=ax,boxprops={'edgecolor': 'k'},
                medianprops={'color':'black'},
                whiskerprops={'color':'k'},
                capprops={'color':'k'})
    annotator = Annotator(ax, pairs, **plotting_parameters)
    annotator.configure(fontsize=fs)

    annotator.set_custom_annotations(formatted_pvalues)
    annotator.annotate()
    plt.xticks(labels=treatments_dict.keys(), ticks=list(range(len(treatments_dict))), rotation=90);

