import matplotlib.pyplot as plt
import os
import sys
import getopt
import numpy as np
import pandas as pd
import seaborn as sns
from src.readfile import save_plot_pdf, GenomeIndex, Histmap


def _read_genes_n_index(gene_file, input_index):
    genes = pd.read_table(gene_file, names=['chr', 'start', 'end', 'geneid'])

    genome_index = GenomeIndex(input_index)

    genes.to_csv('gene_list.tmp', columns=['chr', 'start', 'end', 'geneid'],
                 sep="\t", header=False, index=False)
    genome_index.indexes.to_csv('genome_window.tmp', columns=['chr', 'start', 'end', 'index'],
                                sep="\t", header=False, index=False)

    os.system(f"bedtools intersect -wa -wb -a genome_window.tmp -b gene_list.tmp > gene_inter_bins.bed")

    inter_genes = pd.read_table("gene_inter_bins.bed",
                                names=['chr', 'start', 'end', 'index', 'chr1', 'start1', 'end1', 'gene_name'])
    inter_genes.drop(columns=['chr1', 'start1', 'end1'], inplace=True)

    inter_genes = inter_genes.groupby(by=['chr', 'start', 'end', 'index']).aggregate(lambda x: ",".join(x))
    inter_genes.reset_index(inplace=True)

    inter_genes.index = inter_genes['index']

    os.system("rm gene_inter_bins.bed genome_window.tmp gene_list.tmp")

    return inter_genes.index


def gene_fold_changes(argv):
    title = ""
    try:
        opts, args = getopt.getopt(argv[1:], "t:")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: D2 gene <hist_file> <mark_idx_file> <gene_file> <out_file>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -t Str         Title. default: empty \n")
        # sys.stderr.write("  -v Float       Vmax and Vmin. default: 1.5 \n")
        return 1
    for o, a in opts:
        if o == "-t":
            title = a
    hist_file, index_file, gene_file, out_file = args[:4]

    histmaps = Histmap(hist_file)
    histmaps.average_hist()

    gene_index = _read_genes_n_index(gene_file, index_file)

    inter_index = histmaps.histmap.index.intersection(gene_index)
    gene_histmap = histmaps.histmap.loc[inter_index, :].copy()

    hist_sum_arr = np.array(gene_histmap).mean(axis=0)
    hist_sum_arr = hist_sum_arr.reshape(histmaps.shape)

    all_hist_sum_arr = np.array(histmaps.histmap).mean(axis=0)
    all_hist_sum_arr = all_hist_sum_arr.reshape(histmaps.shape)

    nonzero_idx = (hist_sum_arr > 0) & (all_hist_sum_arr > 0)
    hist_sum_arr[nonzero_idx] = np.log(hist_sum_arr[nonzero_idx]) - np.log(all_hist_sum_arr[nonzero_idx])

    fig, ax = plt.subplots(1, 1, figsize=(5, 5))

    g = sns.heatmap(data=hist_sum_arr, cmap='vlag', ax=ax,
                    vmin=-1, vmax=1, center=0, cbar=False,
                    mask=histmaps.mask.reshape(histmaps.shape))

    ax.invert_yaxis()
    ax.set_aspect('equal')
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.set_yticks([])
    ax.set_xticks([])
    plt.ylabel('Distance to periphery')
    plt.xlabel('DNA Density')
    plt.title(title)

    plt.tight_layout()

    if out_file is None:
        plt.show()
    else:
        save_plot_pdf(out_file)
