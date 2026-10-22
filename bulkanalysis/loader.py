#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
# SPDX-FileCopyrightText: Copyright © 2026 University of Bern
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@unibe.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import pandas as pd


def bulk_aggregate_featureCounts_outputs(files: list, names: list, save_path: None):
    """Merge featureCounts outputs (one file per sample) into one count matrix.

    The sample column of each file (named after its BAM file, read from the file's first line) is renamed with the
    corresponding name of `names`.

    Args:
        files: featureCounts output files.
        names: sample names, in the order of `files`.
        save_path: if not None, the matrix is also written there (CSV).

    Returns:
        pd.DataFrame: counts, features (`Geneid`) in rows, samples in columns.
    """

    bamfile_names = []
    dfs = []
    for i, file in enumerate(files):
        with open(file) as f:
            first_line = f.readline()
            bamfile_name = [bamfile for bamfile in first_line.split('"') if bamfile.endswith("bam")][0]
        
        df = pd.read_csv(file, 
                            sep='\t', 
                            skiprows=1)
        df.rename(columns={bamfile_name: names[i]}, inplace=True)
        if i > 0:
            df.drop('Geneid', axis=1, inplace=True)
        dfs.append(df)
        

    
    df_final = pd.concat(dfs, axis=1)
    df_final.set_index('Geneid', inplace=True)

    if save_path is not None:
        df_final[names].to_csv(save_path)
    return df_final[names]