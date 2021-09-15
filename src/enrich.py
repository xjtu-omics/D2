import matplotlib.pyplot as plt
import sys
import glob
import getopt
import numpy as np
import pandas as pd
import seaborn as sns
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import pdist
from src.readfile import Values, save_plot_pdf, Histmap


def to_value_hist(hist_maps, value_file):
    hist_maps.average_hist()
    value_hists = {}

    i = 0
    values = Values(value_file)
    for value_name in values.value_names:
        value_col = values[value_name]
        value_hist_df = pd.concat([hist_maps.histmap, value_col], axis=1, join='inner')

        # Multiply bin histmap with experimental value
        value_hist = np.multiply(np.array(value_hist_df.iloc[:, :-1]).T,
                                 np.array(value_hist_df.iloc[:, -1]))
        # Average by histmap
        value_hist = np.sum(value_hist.T, axis=0)
        mean_hist = np.sum(value_hist_df.iloc[:, :-1], axis=0)
        value_hist[mean_hist > 0] = value_hist[mean_hist > 0] / mean_hist[mean_hist > 0]

        value_hists[value_name] = (value_hist - value_hist.mean()) / value_hist.std()

        i += 1
        print(f'Value {value_name}, {i}/{len(values.value_names)} DONE.')

    value_hists = pd.DataFrame(value_hists).T
    value_hists = Histmap(input_hist=value_hists, input_bins=hist_maps.bins)

    return value_hists


def plot_value_hists(argv):
    vmax = 1.5
    title = ""
    try:
        opts, args = getopt.getopt(argv[1:], "t:v:d:")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: D3 enrich <hist_file> <mark_idx_file> <output>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -t Str         Tissue name. default: empty \n")
        sys.stderr.write("  -v Float       Vmax and Vmin. default: 1.5 \n")
        return 1
    for o, a in opts:
        if o == "-t":
            title = a
        if o == "-v":
            vmax = a
    hist_file, value_file, output = args[:3]

    hist_maps = Histmap(hist_file)
    mask = hist_maps.mask

    value_hists = to_value_hist(hist_maps, value_file)

    for value_name in value_hists.histmap.index:
        value_hist = np.array(value_hists.histmap.loc[value_name, :])
        value_hist = value_hist.reshape(value_hists.shape)

        ax = sns.heatmap(data=value_hist, cmap='vlag', vmin=-vmax, vmax=vmax, center=0,
                         xticklabels=value_hists.x_draw_bins, yticklabels=value_hists.y_draw_bins,
                         mask=mask.reshape(value_hists.shape))
        ax.invert_yaxis()
        plt.ylabel('Distance to periphery')
        plt.xlabel('DNA Density')
        plt.title(f'{title} {value_name}')
        save_plot_pdf(f'{output}_{value_name}_histplot.pdf')
        # plt.show()
        plt.close()


def hierarchy_hist(argv):
    vmax = 2
    flip_cluster = False
    try:
        opts, args = getopt.getopt(argv[1:], "v:f:")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: D3 enrich <hist_file> <mark_idx_file> <output>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -v Float       Vmax and Vmin. default: 2 \n")
        sys.stderr.write("  -f Flip        Either flip the ranking result. default: False \n")
        return 1
    for o, a in opts:
        if o == "-v":
            vmax = a
        if o == "-f":
            flip_cluster = int(a)
    hist_file, value_file, output = args[:3]

    hist_maps = Histmap(hist_file)
    mask = hist_maps.mask
    value_hists = to_value_hist(hist_maps, value_file)

    ###########################
    # value cluster plot
    value_histmap = np.array(value_hists.histmap)[:, ~mask]
    d = pdist(value_histmap.T)
    L = sch.linkage(d, method='ward', optimal_ordering=True)

    g = sns.clustermap(value_histmap, cmap='vlag', vmin=-vmax, vmax=vmax,
                       col_linkage=L, xticklabels=False,
                       yticklabels=value_hists.histmap.index, row_cluster=False)
    if output is None:
        plt.show()
    else:
        save_plot_pdf(f'{output}_value_hierarchy.pdf')
    plt.close()

    ###########################
    # rank plot
    idx, mask_idx2hist_idx = 0, {}
    for i in range(len(mask)):
        if not mask[i]:
            mask_idx2hist_idx[idx] = i
            idx += 1
    reorder_ind = np.array(g.dendrogram_col.reordered_ind)
    reorder_map = np.zeros(value_hists.histmap.shape[1], dtype=int)
    for i in range(len(reorder_ind)):
        if flip_cluster != 0:
            reorder_map[mask_idx2hist_idx[reorder_ind[i]]] = value_hists.histmap.shape[1] - i
        else:
            reorder_map[mask_idx2hist_idx[reorder_ind[i]]] = i
    reorder_map = reorder_map.reshape(value_hists.shape)
    ax = sns.heatmap(data=reorder_map, cmap='Spectral', mask=mask.reshape(value_hists.shape),
                     xticklabels=value_hists.x_draw_bins, yticklabels=value_hists.y_draw_bins)
    ax.invert_yaxis()
    plt.ylabel('Distance to periphery')
    plt.xlabel('DNA Density')

    if output is None:
        plt.show()
    else:
        save_plot_pdf(f'{output}_hierarchy_hist.pdf')
    plt.close()

    ###########################
    # cluster plot
    fcluster = sch.fcluster(L, t=4, criterion="maxclust")

    fcluster_map = np.zeros(value_hists.histmap.shape[1], dtype=int)
    fcluster_map[~mask] = fcluster
    fcluster_map = fcluster_map.reshape(value_hists.shape)

    ax = sns.heatmap(data=fcluster_map, cmap='Spectral', mask=mask.reshape(value_hists.shape),
                     xticklabels=value_hists.x_draw_bins, yticklabels=value_hists.y_draw_bins)
    ax.invert_yaxis()
    plt.ylabel('Distance to periphery')
    plt.xlabel('DNA Density')

    if output is None:
        plt.show()
    else:
        save_plot_pdf(f'{output}_hierarchy_cluster_hist.pdf')
    plt.close()

