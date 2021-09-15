from src.readfile import Histmap, save_plot_pdf
import sys
import getopt
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def get_one_weighted_std(x, map_idx=None):
    x = np.array(x)

    if x.sum() == 0:
        return 0, 0

    average = np.average(map_idx, weights=x)
    variance = np.average((map_idx - average) ** 2, weights=x)
    return average, np.sqrt(variance)


def get_spread(histmaps, size=None):
    map_idx_df = histmaps.get_map_idx_df()

    if size is not None:
        histmaps.histmap = histmaps.histmap.iloc[:size, :]

    x_spread = histmaps.histmap.apply(get_one_weighted_std, axis=1, result_type='expand',
                                      map_idx=np.array(map_idx_df[:, 0]))
    x_spread.rename(columns={0: 'mean', 1: 'std'}, inplace=True)

    y_spread = histmaps.histmap.apply(get_one_weighted_std, axis=1, result_type='expand',
                                      map_idx=np.array(map_idx_df[:, 1]))
    y_spread.rename(columns={0: 'mean', 1: 'std'}, inplace=True)

    spread = {}
    x_draw_bins, y_draw_bins = histmaps.x_draw_bins, histmaps.y_draw_bins
    spread['den_mean'] = x_spread['mean'] * (x_draw_bins[1] - x_draw_bins[0]) + x_draw_bins[0]
    spread['den_std'] = x_spread['std'] * (x_draw_bins[1] - x_draw_bins[0]) + x_draw_bins[0]
    spread['dtp_mean'] = y_spread['mean'] * (y_draw_bins[1] - y_draw_bins[0]) + y_draw_bins[0]
    spread['dtp_std'] = y_spread['std'] * (y_draw_bins[1] - y_draw_bins[0]) + y_draw_bins[0]

    spread = pd.DataFrame(spread)
    spread.index = histmaps.histmap.index
    return spread


def ave_std(argv):
    fig_output = None
    try:
        opts, args = getopt.getopt(argv[1:], "f:")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: D3 ave [options] <hist_file> <out_file>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -f STR         Figure output. default: None.\n")
        return 1
    for o, a in opts:
        if o == "-f":
            fig_output = a
            if fig_output == 'None':
                fig_output = None
    hist_file, out_file = args[:2]

    histmaps = Histmap(hist_file)
    spread = get_spread(histmaps)
    spread.to_csv(out_file, sep="\t")

    if fig_output is not None:
        ran_size = 10000
        sub_spread = spread.iloc[np.random.choice(np.arange(spread.shape[0]), ran_size), :]

        sns.kdeplot(data=sub_spread, x='den_mean', y='den_std', color="#cbd5e8",
                    levels=5, lw=.1, thresh=.25)
        sns.histplot(data=spread, x='den_mean', y='den_std',
                     cbar=True, cbar_kws=dict(shrink=.5), thresh=5)
        plt.xlabel('Density Average')
        plt.ylabel('Density Std')
        sns.despine(offset=2, trim=True)
        # plt.show()
        save_plot_pdf(f'{fig_output}_den_mean_std.pdf')
        plt.close()

        sns.kdeplot(data=sub_spread, x='dtp_mean', y='dtp_std', color="#cbd5e8",
                    levels=5, lw=.1, thresh=.25)
        sns.histplot(data=spread, x='dtp_mean', y='dtp_std',
                     cbar=True, cbar_kws=dict(shrink=.5), thresh=5)
        plt.xlabel('Density Average')
        plt.ylabel('Density Std')
        sns.despine(offset=2, trim=True)
        # plt.show()
        save_plot_pdf(f'{fig_output}_dtp_mean_std.pdf')
        plt.close()
