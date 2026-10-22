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

def create_design_matrix_for_DGE(df_counts_genes, treatments, path_to_results):
    """Write the inputs of a differential expression analysis (e.g. edgeR) to `<path_to_results>/DGE/`:
    `df_bulk_for_dge.csv` (raw counts of the samples of `treatments` present in `df_counts_genes`) and
    `design_for_dge.csv` (Sample, Treatment, Replicate).

    Args:
        df_counts_genes: gene count matrix, samples in columns.
        treatments: {treatment: sample names}.
        path_to_results: output directory.
    """
    # Take the raw counts 
    if not os.path.exists(os.path.join(path_to_results, "DGE")):
        os.makedirs(os.path.join(path_to_results, "DGE"))



    treatment_col = []
    replicate_col = []
    sample_col = []
    for name, ids in treatments.items():
        ids = [id for id in ids if id in df_counts_genes.columns]
        sample_col = sample_col + ids
        treatment_col = treatment_col + [name]*len(ids)
        replicate_col = replicate_col + list(range(1 ,len(ids)+1))


    df_counts_genes[sample_col].to_csv(os.path.join(path_to_results, "DGE/df_bulk_for_dge.csv"))
    print(len(sample_col))
    print(len(treatment_col))
    print(len(replicate_col))
    design = pd.DataFrame({'Sample': sample_col, 
                           'Treatment': treatment_col, 
                           'Replicate': replicate_col})

    design.to_csv(os.path.join(path_to_results, "DGE/design_for_dge.csv"))




