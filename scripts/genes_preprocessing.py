#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
# SPDX-FileCopyrightText: Copyright © 2026 University of Bern
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@unibe.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import argparse
import logging
import os
import sys
import yaml
import numpy as np
import pandas as pd
from bulkanalysis.plot import plot_gene_expression_distribution
from bulkanalysis.filtering import gmm_threshold, filter_samples_based_on_thresholds, remove_outliers_samples
from kallistomanager.loader import merge_kallisto_counts, from_whole_transcriptome_quantification_to_coding_genes
from bulkanalysis.normalization import quantile_normalize
from bulkanalysis.dge import create_design_matrix_for_DGE

import matplotlib
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42


def format_genes_from_kallisto_whole_transcriptome(df_counts_path,
                                             gene_names_path,
                                             gene_info_path, data_list=None):
    """Merge the kallisto outputs of `df_counts_path` (or of `data_list`) into `<df_counts_path>/counts_matrix.csv`,
    then sum the transcripts per protein-coding gene (kallistomanager).

    Args:
        df_counts_path: folder of the kallisto outputs, where the merged matrix is also written.
        gene_names_path: GENCODE gene names (e.g. Gencode_geneNames_hg38V44.txt).
        gene_info_path: NCBI gene information (e.g. Homo_sapiens.gene_info), to keep the protein-coding genes.
        data_list: optional list of kallisto `abundance.tsv` files to merge instead of the folder.

    Returns:
        pd.DataFrame: gene count matrix, genes in rows, samples in columns.
    """
    merge_kallisto_counts(data_path=df_counts_path, data_list=data_list, save_dir=df_counts_path)
    df_genes = from_whole_transcriptome_quantification_to_coding_genes(df_counts_path=os.path.join(df_counts_path, 
                                                                                                   "counts_matrix.csv"),
                                                                                            gene_names_path=gene_names_path,
                                                                                            gene_info_path=gene_info_path)
    return df_genes

def readjust_treatments(actual_samples, treatments):
    """Remove from each treatment of `treatments` ({treatment: sample names}) the samples absent from `actual_samples`,
    printing each removal.
    """
    # build new lists: removing items from the list being iterated skipped the sample following each removed one
    for treatment, samples in treatments.items():
        for sample in samples:
            if sample not in actual_samples:
                print(f"Removing sample {sample} from treatment {treatment} as it is not in the actual samples.")
        treatments[treatment] = [sample for sample in samples if sample in actual_samples]
    return treatments


def main():
    """Gene-level preprocessing of bulk RNA-seq data, configured by a YAML file (see README and
    config/genes_preprocessing_template.yaml): merge the quantification, filter the reliably expressed genes per sample
    (Gaussian mixture threshold), quantile-normalise, and write the filtered matrix, the inputs of a differential
    expression analysis (DGE/) and the QC figures. The stdout of the run is written to `<path_to_results>/message.log`.
    """
    FLAGS = argparse.ArgumentParser(description='Length analysis for bulk RNAseq data.')
    FLAGS.add_argument('--config_file', help='Location of config file')
    args = FLAGS.parse_args()
    with open(args.config_file, 'r') as stream:
        config = yaml.safe_load(stream)


    if not os.path.exists(config['path_to_results']):
        os.makedirs(config['path_to_results'])

    old_stdout = sys.stdout

    log_file = open(os.path.join(config['path_to_results'], "message.log"),"w")

    sys.stdout = log_file


    if not os.path.exists(os.path.join(config['path_to_results'], 'Figures')):
        os.makedirs(os.path.join(config['path_to_results'], 'Figures'))

    path_to_figures = os.path.join(os.path.join(config['path_to_results'], 'Figures'))

    if 'data_list' in config.keys():
        data_list = config['data_list']
    else: 
        data_list = None

    if config['data_origin'] == 'kallisto_whole_transcriptome':
        # check if df_counts_path is a path or a file
        if os.path.isdir(config['df_counts_path']):
            df_genes = format_genes_from_kallisto_whole_transcriptome(df_counts_path=config['df_counts_path'],
                                                                       gene_names_path=config['gene_names_path'],
                                                                       gene_info_path=config['gene_info_path'], 
                                                                       data_list=data_list)
        elif os.path.isfile(config['df_counts_path']):
            df_genes = pd.read_csv(config['df_counts_path'], index_col=0)
        else:
            raise ValueError(f"df_counts_path {config['df_counts_path']} is neither a directory nor a file.")
    else:
        raise NotImplementedError("This script does not fit with the provided data origin")
    
    print(f"Samples at this step 1: {df_genes.columns}")

    
    
    if "remove_outliers" not in config.keys():
        remove_outliers = True
    else:
        remove_outliers = config['remove_outliers']
        
    if "quantile_threshold" not in config.keys():
        quantile_threshold = 0.95
    else: 
        quantile_threshold = config['quantile_threshold']
    print(f"Using threshold: {quantile_threshold}")
    
    # Here add possibility to aggregate SRR ids from similar samples
    if 'aggregation' in config.keys():
        for key, values in config['aggregation'].items():
            df_genes[key] = df_genes[values].sum(axis=1)
            df_genes.drop(columns=values, inplace=True)
            
    # Keep only the samples defined in the treatments dictionary
    if 'treatments' in config.keys():
        treatments = config['treatments']
        all_samples = []
        for treatment, samples in treatments.items():
            all_samples.extend(samples)
        df_genes = df_genes[all_samples]

    print(f"Samples at this step 2: {df_genes.columns}")
    
    print(f"Before filtering, we have {len(df_genes)} genes.")
    print(f"Before filtering, we have treatments: {config['treatments']}")

    df_genes = df_genes[df_genes.sum(axis=1) != 0]
    plot_gene_expression_distribution(df_genes=df_genes, 
                                      title="Gene expression before filtering", 
                                      save_file=os.path.join(path_to_figures,
                                                             f"gene_exp_before_filtering.{config['figures_extension']}"))
    

    df_genes_log = np.log2(df_genes+1)
    ## Compute the thresholds for reliably expressed genes for each sample
    thresholds = {}
    for col in df_genes_log.columns:

        threshold = gmm_threshold(data= np.array(df_genes_log[df_genes_log[col] != 0][col]), 
                                        quantile_threshold=quantile_threshold,
                                        background_quantile=True, 
                                        density=True, show_plot=False, data_name=col, 
                                        saving_folder=path_to_figures, 
                                        extension=config['figures_extension'],
                                        return_int=False)
        thresholds[col] = threshold
    print(f"Samples at this step 3: {df_genes.columns}")
    df_genes = remove_outliers_samples(df_genes=df_genes,
                                        thresholds=thresholds,
                                        save_path=path_to_figures, 
                                        extension=config['figures_extension'],
                                        remove=remove_outliers)
    print(f"Samples at this step 4: {df_genes.columns}")
    
    treatments = readjust_treatments(actual_samples=df_genes.columns, treatments=config['treatments'])
    print(f"All treatments: {config['treatments']}")
    print(f"All thresholds: {thresholds}")
    print(f"Actual samples: {df_genes.columns}")
    print(f"Actual treatments: {treatments}")

    # Remove the items in thresholds that are not in the actual samples
    thresholds = {key: value for key, value in thresholds.items() if key in df_genes.columns}
    print(f"Actual thresholds: {thresholds}")

    df_genes_filtered = filter_samples_based_on_thresholds(df_genes=df_genes, 
                                                           treatments=treatments, 
                                                           thresholds=thresholds, 
                                                           pct_in_treatment=config['pct_in_treatment'], 
                                                           save_file=os.path.join(path_to_figures, 
                                                                                  f"gene_exp_after_filtering.{config['figures_extension']}"))
   

    

    create_design_matrix_for_DGE(df_counts_genes=df_genes_filtered, treatments=treatments, path_to_results=config['path_to_results'])
    df_filtered_normalized = quantile_normalize(df_genes_filtered)
    plot_gene_expression_distribution(df_genes=df_filtered_normalized,
                                      title="Gene expression after filtering and quantile normalization",
                                      save_file=os.path.join(path_to_figures,
                                                             f"gene_exp_after_filtering_and_normalization.{config['figures_extension']}"))
    
    print(f"After filtering and normalization, we have {len(df_filtered_normalized)} remaining genes.")

    df_filtered_normalized.to_csv(os.path.join(config['path_to_results'], 'df_genes_filtered_normalized.csv'))


    sys.stdout = old_stdout

    log_file.close()

if __name__ == "__main__":
    main()
