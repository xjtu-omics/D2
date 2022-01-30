import matplotlib.pyplot as plt
import os
import sys
import getopt
import numpy as np
import pandas as pd
import seaborn as sns
from src.readfile import save_plot_pdf, GenomeIndex, Histmap

mouse_active_regions = np.array(list(range(120, 126)) + list(range(135, 142)) + \
                                list(range(150, 158)) + list(range(165, 174)) + \
                                list(range(180, 189)) + list(range(195, 205)) + \
                                list(range(210, 220)))
human_active_regions = np.array(list(range(47, 52)) + list(range(60, 67)) + \
                                list(range(75, 82)) + list(range(90, 98)) + \
                                list(range(105, 113)) + list(range(120, 128)) + \
                                list(range(135, 143)) + list(range(150, 158)) + \
                                list(range(165, 173)) + list(range(180, 188)) + \
                                list(range(195, 203)) + list(range(210, 218)))

def compute_activation_index(argv):
    active_regions = 'mouse'
    try:
        opts, args = getopt.getopt(argv[1:], "a:")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: D2 act <hist_file> <out_file>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -a Str       'mouse' or 'human' for default active regions, or input hierarchy cluster file. default: 'mouse'\n")
        # sys.stderr.write("  -f Flip        Either flip the ranking result. default: False \n")
        return 1
    for o, a in opts:
        if o == "-a":
            active_regions = a
        # if o == "-f":
        #     flip_cluster = int(a)
    hist_file, out_file = args[:2]

    if active_regions == 'mouse':
        active_regions = mouse_active_regions
    elif active_regions == 'human':
        active_regions = human_active_regions
    else:
        hiera_map = pd.read_table(active_regions, names=['cluster'])
        active_regions = np.array(hiera_map[(hiera_map['cluster'] == 1) |
                                            (hiera_map['cluster'] == 2)].index)

    histmap = Histmap(hist_file)

    all_hist_sum_arr = np.array(histmap.histmap).mean(axis=0)
    all_hist_sum_arr /= np.sum(all_hist_sum_arr)

    histmap.average_hist()
    a = np.log(histmap.histmap) - np.log(all_hist_sum_arr)
    a[histmap.histmap == 0] = 0
    histmap.histmap = a

    active_index = histmap.histmap.iloc[:, active_regions]
    active_index = np.sum(active_index, axis=1)

    active_index = (active_index - np.mean(active_index)) / np.std(active_index)

    active_index = pd.DataFrame(active_index, index=histmap.histmap.index, columns=['active_index'])
    active_index.dropna(inplace=True)
    active_index.to_csv(out_file, sep="\t")
