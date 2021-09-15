import sys
import getopt
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from src.readfile import DenDtps


def define_hist_bins(den_dtps):
    borders = []

    dtp_array = np.array(den_dtps[[f'dtp_{cell}' for cell in den_dtps.cell_names]])
    dtp_array = dtp_array.reshape((-1,))
    dtp_array = dtp_array[~np.isnan(dtp_array)]
    dtp_array = dtp_array[dtp_array > 0]
    borders.append(np.min(dtp_array))
    borders.append(np.percentile(dtp_array, 95))

    print(f'DisTP border: {borders[0]} ~ {borders[1]}')

    den_array = np.array(den_dtps[[f'den_{cell}' for cell in den_dtps.cell_names]])
    den_array = den_array.reshape((-1,))
    den_array = den_array[~np.isnan(den_array)]
    borders.append(den_array.mean() - den_array.std() * 1.96)
    borders.append(den_array.mean() + den_array.std() * 1.96)

    print(f'Density border: {borders[2]} ~ {borders[3]}')

    return borders


def draw_2d_hist(den_dtps, out_file=None, ran_size=50000, border=None):
    dtp_array = np.array(den_dtps[[f'dtp_{cell}' for cell in den_dtps.cell_names]])
    dtp_array = dtp_array.reshape((-1,))

    den_array = np.array(den_dtps[[f'den_{cell}' for cell in den_dtps.cell_names]])
    den_array = den_array.reshape((-1,))

    dtp_array, den_array = dtp_array[~np.isnan(dtp_array)], den_array[~np.isnan(dtp_array)]
    dtp_array, den_array = dtp_array[~np.isnan(den_array)], den_array[~np.isnan(den_array)]
    dtp_array, den_array = dtp_array[dtp_array > 0], den_array[dtp_array > 0]

    index = np.random.choice(np.arange(len(dtp_array)), ran_size, replace=False)
    dtp_array, den_array = dtp_array[index], den_array[index]

    g = sns.JointGrid(x=den_array, y=dtp_array, marginal_ticks=True)
    g.plot_joint(sns.kdeplot, color="#1f78b4")
    g.plot_joint(plt.scatter, c='lightgray', s=5)
    g.plot_marginals(sns.histplot, bins=40, kde=True, color="#1f78b4",
                     edgecolor='white', stat='probability')

    g.ax_joint.set_xlim([0, 4.5])
    g.ax_marg_x.set_xlim([0, 4.5])
    g.ax_joint.set_ylim([-1, 25])
    g.ax_marg_y.set_ylim([-1, 25])

    if border is not None:
        g.ax_joint.axhline(border[0], c='#fdc086')
        g.ax_joint.axhline(border[1], c='#fdc086')
        g.ax_joint.axvline(border[2], c='#fdc086')
        g.ax_joint.axvline(border[3], c='#fdc086')

        g.ax_marg_y.axhline(border[0], c='#fdc086')
        g.ax_marg_y.axhline(border[1], c='#fdc086')

        g.ax_marg_x.axvline(border[2], c='#fdc086')
        g.ax_marg_x.axvline(border[3], c='#fdc086')

        g.ax_marg_x.set_title(f'Den Min: {np.round(border[2], 2)}, '
                              f'Den Max: {np.round(border[3], 2)}\n'
                              f'DTP Min: {np.round(border[0], 2)}, '
                              f'DTP Max: {np.round(border[1], 2)}')

    g.ax_joint.set_xlabel('Density')
    g.ax_joint.set_ylabel('DisTP')

    plt.tight_layout()

    sns.despine(offset=2, trim=True, ax=g.ax_marg_y)
    sns.despine(offset=2, trim=True, ax=g.ax_marg_x)

    if out_file is None:
        plt.show()
    else:
        plt.savefig(f"{out_file}.png")
    plt.close()


def draw_hist(argv):
    try:
        opts, args = getopt.getopt(argv[1:], "i:")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: D3 sta <den_dtp_dir> <index_file> <out_file>\n")
        return 1
    # for o, a in opts:
    #     if o == "-d":
    #         debug = int(a)
    den_dir, index_file, out_file = args[:3]

    den_dtps = DenDtps(den_dir, index_file, join='outer')

    border = define_hist_bins(den_dtps)
    draw_2d_hist(den_dtps, out_file=out_file, border=border)


def to_histmap(argv):
    n_bins = 15 + 1
    den_min, den_max = 1, 3.2
    dtp_min, dtp_max = 1.21, 16
    try:
        opts, args = getopt.getopt(argv[1:], "i:")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: D3 map [options] <den_dtp_dir> <index_file> <out_file>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -n Int         Bin number. default: 15.\n")
        sys.stderr.write("  -ei Float      Density min. default: 1.\n")
        sys.stderr.write("  -ea Float      Density max. default: 3.2.\n")
        sys.stderr.write("  -ti Float      DisTP min. default: 1.21.\n")
        sys.stderr.write("  -ta Float      DisTP max. default: 16.\n")
        return 1
    for o, a in opts:
        if o == "-n":
            n_bins = int(a) + 1
        if o == "-ei":
            den_min = float(a)
        if o == "-ea":
            den_max = float(a)
        if o == "-ti":
            dtp_min = float(a)
        if o == "-ta":
            dtp_max = float(a)
    den_dir, index_file, out_file = args[:3]

    hist_bins = (np.linspace(den_min, den_max, n_bins),
                 np.linspace(dtp_min, dtp_max, n_bins))

    den_dtps = DenDtps(den_dir, index_file, join='outer')

    histmaps = den_dtps.get_den_dtp_histmap(hist_bins)
    histmaps.output_histmap(out_file, format="%d", )
