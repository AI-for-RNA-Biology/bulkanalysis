<!--
#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
# SPDX-FileCopyrightText: Copyright © 2026 University of Bern
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@unibe.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
-->

# Project

This is a package to help with the analysis of bulk RNA sequencing data. This package is built on scanpy.

## Modules
- `loader`: merge featureCounts outputs into a count matrix (`bulk_aggregate_featureCounts_outputs`).
- `filtering`: gene and sample filtering (Gaussian mixture expression thresholds, expression in all replicates,
  genes overlapping a gene on the opposite strand, outlier samples).
- `normalization`: quantile and CPM normalisation (`quantile_normalize` uses the R package `limma` through `rpy2`,
  `CPM_normalize` uses `bioinfokit`: install the `normalization` extra, `pip install ".[normalization]"`; the
  `genes_preprocessing` script needs it).
- `dge`: inputs of a differential expression analysis (count matrix and design).
- `dapars`: load DaPars outputs and compare PDUI between treatments.
- `plot`: expression distributions and PCA box plots.

`bulkanalysis` depends on the `genomemanager` and `kallistomanager` packages of
[AI-for-RNA-Biology](https://github.com/AI-for-RNA-Biology), installed with it.

## Origin

bulkanalysis was first developed at the [Idiap Research Institute](https://www.idiap.ch), where [version 0.1](https://github.com/idiap/bulkanalysis) was released. It has been developed since at the University of Bern: this repository (version 0.2.0 and later) continues it.

## Installation of this package in another project

To install this library in another Python project, execute simply:

```bash
pip install git+https://github.com/AI-for-RNA-Biology/bulkanalysis.git
```

## Scripts
### `genes_preprocessing`
`genes_preprocessing.py` aims at performing basic filtering and normalization on bulk RNA sequencing data. To run it, udpate the file `config/genes_preprocessing_template.yaml` with the right paths. You can then run the script as follows:
```bash
python3 scripts/genes_preprocessing.py --config_file config/genes_preprocessing_template.yaml
```

The config file should contain:
- `data_origin`: Origin of the transcripts matrix. For now, the only supported option is `kallisto_whole_transcriptome`, meaning that the transcripts matrix must come from a kallisto quantification on a whole transcriptome. Other options might be supported in the future.
- `df_counts_path`: Path to the matrix of transcripts.
- `gtf_file_no_focus`: Original GTF used for the whole transcriptome quantification.
- `gene_names_path`: GENCODE genes names with gene symbols. Example: "Gencode_geneNames_hg38V44.txt"
- `gene_info_path`: Genes symbol with their information, in particular whether they are protein-coding. Example: "Homo_sapiens.gene_info"
- `treatments`: dictionnary with name of the treatments in keys and list of corresponding sample names in keys.
- `path_to_results`: Directory where to save the results.
- `figures_extension`: Extension you want to save your figures with, e.g "pdf", "png",...
- `pct_in_treatment`: Percentage of samples within a treatment group in which a gene should be reliably expressed to be kept.

### `aggregate_featureCounts_output`
`aggregate_featureCounts_output.py` aims at merging the outputs from `featureCounts` for multiple samples.
 To run it, run the script as follows:
```bash
python3 scripts/aggregate_featureCounts_output.py -f sample1.txt sample2.txt sample3.txt -n sample_name1 sample_name2 sample_name3 -s df_counts.csv

```
with:
- `sample1.txt sample2.txt sample3.txt` being the output files of the function [featureCounts](#https://subread.sourceforge.net/featureCounts.html)
- `sample_name1 sample_name2 sample_name3` being the names of the samples you want to appear in the final matrix
- `df_counts.csv`: name of the file where to save the final matrix.