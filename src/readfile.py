import pandas as pd
import numpy as np
import pickle
import glob
import os
import warnings
import seaborn as sns
import matplotlib.pyplot as plt

###################################
# paras
index_names = ['chr', 'start', 'end', 'index']
dg_names = ['chr', 'start', 'end', 'x', 'y', 'z', 'filter']

den_dtp_suffix = ".den_dtp.txt"
den_names = ['filter', 'x_loc', 'y_loc', 'z_loc', 'den', 'dtp', 'min']

min_bin_count = 200


###################################
# Utilities
def save_plot_pdf(out_file):
    import matplotlib as mpl

    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42

    import matplotlib.pyplot as plt

    if out_file[-4:] != '.pdf':
        out_file = out_file + '.pdf'

    plt.tight_layout()
    plt.savefig(out_file)
    plt.close()


###################################
# Basic Class
class GenomeIndex:
    def __init__(self, input_index, n_chro=None):
        if isinstance(input_index, str):
            self.indexes = pd.read_table(input_index, names=index_names,
                                         dtype={'chr': str, 'index': int})

        elif isinstance(input_index, pd.DataFrame):
            self.indexes = input_index

        else:
            raise ValueError('Unrecognized input index type.')

        self.n_chro = n_chro

        self.chros = None
        self.chr_seps = None
        self.chr_lens = None

        self.all_len = None
        self.resolution = None

        self.get_index_features()

    def __getitem__(self, item):
        return self.indexes[item]

    def __setitem__(self, key, value):
        self.indexes[key] = value

    def get_index_features(self):
        self.indexes.index = self.indexes['index']

        self.all_len = np.max(self.indexes['index']) + 1

        self.chros = list(set(self.indexes['chr']))
        self.chros.sort()
        if self.n_chro is not None:
            self.chros = self.chros[:self.n_chro]
        else:
            self.n_chro = len(self.chros)

        self.chr_seps, self.chr_lens = {}, {}
        for chro in self.chros:
            chr_indexs = self.indexes[self.indexes['chr'] == chro]
            self.chr_seps[chro] = [np.min(chr_indexs['index']), np.max(chr_indexs['index'])]
            self.chr_lens[chro] = self.chr_seps[chro][1] - self.chr_seps[chro][0] + 1

        self.resolution = self.indexes.iloc[0, 2] - self.indexes.iloc[0, 1]

        return self.indexes

    def get_chr_index(self):

        for chro in self.chros:
            self.indexes.loc[self.indexes['chr'] == chro, 'chr_index'] = \
                self.indexes.loc[self.indexes['chr'] == chro, 'index'] - self.chr_seps[chro][0]

    def get_new_index(self, new_bed):
        new_chros = list(set(new_bed['chr']))

        if new_chros == self.chros:
            new_gen_index = self

        else:
            chr_gen_index = np.zeros(self.indexes.shape[0], dtype=bool)
            for chro in new_chros:
                chr_gen_index = chr_gen_index | (self.indexes['chr'] == chro)

            new_gen_index = self.indexes[chr_gen_index]
            new_gen_index = GenomeIndex(input_index=new_gen_index)

        return new_gen_index

    def get_index(self, chro, start, end):
        chr_index = self.indexes[self.indexes['chr'] == chro]

        start_idx = chr_index[chr_index['start'] >= start]
        start_idx = start_idx.index.array[0]

        end_idx = chr_index[chr_index['end'] <= end]
        end_idx = end_idx.index.array[-1]

        return start_idx, end_idx

    def to_output(self, output):
        self.indexes.to_csv(output, sep="\t", index=False, header=False,
                            columns=['chr', 'start', 'end', 'index'])


class Values:
    def __init__(self, input_value):
        if isinstance(input_value, str):
            self.values = pd.read_table(input_value, comment="#")
            self.values.index = self.values['index']

            self.value_names = list(self.values.columns.values)

            for col in index_names:
                if col in self.value_names:
                    self.value_names.remove(col)

        elif isinstance(input_value, pd.DataFrame):
            self.values = input_value
            self.value_names = list(self.values.columns.values)

            for col in index_names:
                if col in self.value_names:
                    self.value_names.remove(col)

        else:
            raise ValueError('Input value type unrecognized.')

    def __getitem__(self, item):
        return self.values[item]

    def __setitem__(self, key, value):
        self.values[key] = value

    def select_values(self, include_pattern):
        value_names = [value_name for value_name in self.value_names if include_pattern in value_name]
        values = self.values[value_names]

        return Values(input_value=values)

    def normalize_values_zscore(self):
        for col in self.value_names:
            da_ar = self.values[col]
            self.values[col] = (da_ar - np.mean(da_ar)) / np.std(da_ar)


class Dg:
    def __init__(self, input_dg, input_window=None):
        if isinstance(input_dg, str):
            self.dg_df = pd.read_table(input_dg, names=dg_names)

        elif isinstance(input_dg, pd.DataFrame):
            self.dg_df = input_dg

        else:
            raise ValueError(f'Unrecognized input format: {input_dg.type}.')

        if 'index' not in self.dg_df.columns and input_window is None:
            raise ValueError(f'Please input window for dg.')

        elif input_window is not None:
            if isinstance(input_window, GenomeIndex):
                self.window = input_window
            else:
                self.window = GenomeIndex(input_window)

            self.concat_with_window()

        self.shape = self.dg_df.shape

    def __getitem__(self, item):
        return self.dg_df[item]

    def __setitem__(self, key, value):
        self.dg_df[key] = value

    def concat_with_window(self):
        dg_df = self.dg_df.copy()

        dg_df.index = pd.MultiIndex.from_arrays([dg_df['chr'], dg_df['start']])
        dg_df.drop(columns=['chr', 'start', 'end'], inplace=True)

        window = self.window.indexes.copy()
        window.index = pd.MultiIndex.from_arrays([window['chr'], window['start']])

        dg_df = pd.concat([window, dg_df], join='inner', axis=1)
        dg_df.index = dg_df['index']

        self.dg_df = dg_df

    def concat_with_den_dtp(self, den_dtp):
        den_dtp_df = den_dtp.den_dtp.copy()
        den_dtp_df.drop(columns=index_names, inplace=True)
        den_dtp_df.drop(columns=['filter'], inplace=True)

        self.dg_df = pd.concat([self.dg_df, den_dtp_df], axis=1, join='inner')


class DenDtp:
    def __init__(self, input_den, input_index=None):
        # Get gen_index
        if input_index is None:
            gen_index = None
        elif isinstance(input_index, GenomeIndex):
            gen_index = input_index
        else:
            gen_index = GenomeIndex(input_index)

        self.gen_index = gen_index

        # Get den_dtp
        if isinstance(input_den, str):
            den_dtp = pd.read_table(input_den,
                                    names=index_names + den_names,
                                    dtype={'filter': float})
            den_dtp.index = den_dtp['index']
            den_dtp = den_dtp
        elif isinstance(input_den, pd.DataFrame):
            den_dtp = input_den.copy()
        else:
            raise ValueError('No cluster bed input.')

        self.den_dtp = den_dtp

    def __getitem__(self, item):
        return self.den_dtp[item]

    def __setitem__(self, key, value):
        self.den_dtp[key] = value

    def out_to_file(self, output):
        if isinstance(self.den_dtp, pd.DataFrame):
            self.den_dtp.to_csv(output, sep="\t", index=False, header=False,
                                columns=index_names + den_names)
        else:
            raise ValueError(f'self.den_dtp is not pd.Dataframe type, but {self.den_dtp.type}')

    def concat_with_dg(self, dg):
        dg_df = dg.dg_df.copy()
        dg_df.drop(columns=['chr', 'start', 'end', 'index', 'filter'], inplace=True)

        self.den_dtp = pd.concat([self.den_dtp, dg_df], axis=1, join='inner')

    def yield_one_chro(self):
        for chro in self.gen_index.chros:
            new_bed = self.den_dtp[self.den_dtp['chr'] == chro].copy()
            new_gen_index = GenomeIndex(self.gen_index.indexes[self.gen_index.indexes['chr'] == chro].copy())

            new_bed = DenDtp(new_bed, new_gen_index)

            yield chro, new_bed

    def normalize_minmax(self):
        da_ar = self.den_dtp['den']
        self.den_dtp['den'] = (da_ar - np.min(da_ar)) / (np.max(da_ar) - np.min(da_ar))

        da_ar = self.den_dtp['dtp']
        self.den_dtp['dtp'] = (da_ar - np.min(da_ar)) / (np.max(da_ar) - np.min(da_ar))

    def tid_off_low_map(self):
        self.den_dtp = self.den_dtp[self.den_dtp['filter'] == 1]


class DenDtps:
    def __init__(self, input_den, input_index=None, cell_names=None, read_den=True,
                 join='inner', include_pat=None, exclude_pat=None, n_cell=None,
                 tid_off_low_map=True):
        ##################
        # Read index
        if input_index is None:
            gen_index = None
        elif isinstance(input_index, GenomeIndex):
            gen_index = input_index
        else:
            gen_index = GenomeIndex(input_index)

        self.gen_index = gen_index

        ##################
        # Read cluster beds
        if isinstance(input_den, str):
            # Get cell names and bed files
            cell_names = glob.glob(input_den + f'/*{den_dtp_suffix}')

            if include_pat is not None:
                cell_names = [cell_name for cell_name in cell_names if include_pat in cell_name]

            if exclude_pat is not None:
                cell_names = [cell_name for cell_name in cell_names if exclude_pat not in cell_name]

            if n_cell is not None:
                cell_names = cell_names[:n_cell]

            cell_names = [cell_name.split("/")[-1].split(den_dtp_suffix)[0] for cell_name in cell_names]

            self.bed_files = [os.path.join(input_den, cell_name + den_dtp_suffix) for cell_name in cell_names]

            self.cell_names = cell_names
            self.join = join

            if read_den:
                den_dtp = self.get_den_dtps()
            else:
                den_dtp = None

        elif isinstance(input_den, pd.DataFrame):
            if cell_names is not None:
                den_dtp = input_den
                self.cell_names = cell_names
            else:
                raise ValueError('Create ClusterBeds by pd.Dataframe must assign cell_names.')

        else:
            raise ValueError('Cell names should be inputted with cluster beds.')

        self.den_dtp = den_dtp

        if den_dtp is not None:
            self.den_dtp['n_cells'] = np.nansum(den_dtp.loc[:, [f'filter_{cell}' for cell in cell_names]], axis=1)

        if tid_off_low_map and read_den:
            self.tid_off_low_map()

    def __getitem__(self, item):
        return self.den_dtp[item]

    def __setitem__(self, key, value):
        self.den_dtp[key] = value

    def get_den_dtps(self):
        # Read index file as original cluster bed
        den_dtp = self.gen_index.indexes.copy()

        # Read cluster bed file
        i = 0
        for cell_name, bed_file in zip(self.cell_names, self.bed_files):
            one_den_stp = DenDtp(bed_file, input_index=self.gen_index).den_dtp

            for col in den_names:
                den_dtp.loc[one_den_stp.index, f'{col}_{cell_name}'] = one_den_stp[col]

            i += 1
            print(f"Read cluster bed: {i}/{len(self.bed_files)}. {bed_file.split('/')[-1]}")

        if self.join == 'inner':
            den_dtp.dropna(how='any', inplace=True)

        # den_dtp.fillna(value=nan_value, inplace=True)
        print('Read cluster bed DONE.')

        return den_dtp

    def yield_one_cell(self, yield_from_file=True):

        if yield_from_file:
            for cell_name, bed_file in zip(self.cell_names, self.bed_files):
                new_bed = DenDtp(bed_file, self.gen_index)
                yield cell_name, new_bed

        else:
            for cell_name in self.cell_names:
                chang_name_dic = {f'{col}_{cell_name}': col for col in den_names}
                new_bed = self.den_dtp[index_names + list(chang_name_dic.keys())].copy()
                new_bed.rename(columns=chang_name_dic, inplace=True)

                new_bed = DenDtp(new_bed, self.gen_index)
                yield cell_name, new_bed

    def yield_one_chro(self):
        for chro in self.gen_index.chros:
            new_bed = self.den_dtp[self.den_dtp['chr'] == chro].copy()
            new_gen_index = GenomeIndex(self.gen_index.indexes[self.gen_index.indexes['chr'] == chro].copy())

            new_bed = DenDtps(new_bed, new_gen_index, cell_names=self.cell_names)

            yield chro, new_bed

    def mean_col(self, col_names=None, include_nan_row=True):
        col_names = ['den', 'dtp'] if col_names is None else col_names

        for col in col_names:
            value_array = self.den_dtp[[f'{col}_{cell_name}' for cell_name in self.cell_names]]

            if not include_nan_row:
                value_array.dropna(how='any', inplace=True)

                col_mean = np.mean(value_array, axis=1)
                self.den_dtp.loc[value_array.index, f'{col}_mean'] = col_mean

            else:
                value_array = self.den_dtp[[f'{col}_{cell_name}' for cell_name in self.cell_names]]
                value_array = np.array(value_array)

                col_num = np.count_nonzero(~np.isnan(value_array), axis=1)
                col_mean = np.sum(np.nan_to_num(value_array), axis=1)

                non_zero_idx = col_num > 0
                col_mean[non_zero_idx] = col_mean[non_zero_idx] / col_num[non_zero_idx]

                self.den_dtp[f'{col}_mean'] = np.nan_to_num(col_mean)

    def get_den_dtp_histmap(self, bins, size=None):

        def _get_histmap_one_row(x, bins=None):
            a = x[[f'den_{cell}' for cell in self.cell_names]]
            b = x[[f'dtp_{cell}' for cell in self.cell_names]]
            map, _, _ = np.histogram2d(a, b, bins=bins)
            return map.T.reshape((-1,))

        if size is not None:
            den_dtp = self.den_dtp.iloc[np.random.choice(self.den_dtp.shape[0], size), :]
        else:
            den_dtp = self.den_dtp

        histmaps = den_dtp.apply(_get_histmap_one_row, bins=bins,
                                 axis=1, result_type='expand')

        return Histmap(histmaps, bins)

    def normalize_minmax(self):
        for cell_name in self.cell_names:
            da_ar = self.den_dtp[f'den_{cell_name}']
            self.den_dtp[f'den_{cell_name}'] = (da_ar - np.min(da_ar)) / (np.max(da_ar) - np.min(da_ar))

            da_ar = self.den_dtp[f'dtp_{cell_name}']
            self.den_dtp[f'dtp_{cell_name}'] = (da_ar - np.min(da_ar)) / (np.max(da_ar) - np.min(da_ar))

    def tid_off_low_map(self, ratio=0.8):
        filter_sum = np.nansum(self.den_dtp[[f'filter_{cell}' for cell in self.cell_names]], axis=1)

        self.den_dtp = self.den_dtp[filter_sum >= ratio * len(self.cell_names)]


class Histmap:
    def __init__(self, input_hist, input_bins=None):
        self.histmap, self.bins = None, None

        if isinstance(input_hist, pd.DataFrame):
            self.histmap = input_hist
            self.bins = input_bins

            if input_bins is None:
                raise ValueError('Please input bins for Histmap.')

        elif isinstance(input_hist, str):
            self._read_histmap_from_file(input_hist)

        else:
            raise ValueError('Unrecognized input histmap.')

        self.cell_nums = np.sum(self.histmap, axis=1)

        self.x_draw_bins = np.round((self.bins[0][1:] + self.bins[0][:-1]) / 2, 2)
        self.y_draw_bins = np.round((self.bins[1][1:] + self.bins[1][:-1]) / 2, 2)

        self.has_average = False
        self.shape = (len(self.x_draw_bins), len(self.y_draw_bins))

        self.map_idx = np.arange(self.shape[0] * self.shape[1]).reshape(self.shape)
        self.mask = np.array(self.histmap).sum(axis=0) <= min_bin_count

    def __getitem__(self, item):
        return self.histmap[item]

    def __setitem__(self, key, value):
        self.histmap[key] = value

    def _read_histmap_from_file(self, hist_file):
        xbins, ybins = None, None

        file_handle = open(hist_file, 'r')
        for line in file_handle.readlines():
            if line[0] != "#":
                break

            if line[:6] == '#xbins':
                line = line.split('#xbins: ')[1]
                xbins = np.fromstring(line, sep="\t")

            if line[:6] == '#ybins':
                line = line.split('#ybins: ')[1]
                ybins = np.fromstring(line, sep="\t")
        file_handle.close()

        if xbins is None or ybins is None:
            raise ValueError(f'Cannot find xbins or ybins in {hist_file}.')

        self.bins = (xbins, ybins)

        self.histmap = pd.read_table(hist_file, comment="#", index_col=0).astype(float)

    def get_map_idx_df(self):
        # x in map_idx_df denotes x-axis position, which is col index.
        # y in map_idx_df denotes y-axis position, which is row index.
        map_idx_df = {'x': [], 'y': []}

        for i in range(self.map_idx.min(), self.map_idx.max() + 1):
            x, y = np.where(self.map_idx == i)

            map_idx_df['x'].append(y[0])
            map_idx_df['y'].append(x[0])

        map_idx_df = np.array(pd.DataFrame(map_idx_df))

        return map_idx_df

    def average_hist(self):
        if self.has_average:
            return 0

        hist_sum = np.sum(np.array(self.histmap), axis=1)
        non_zero_idx = hist_sum > 0

        data = self.histmap.values[non_zero_idx, :].T / hist_sum[non_zero_idx]
        self.histmap.values[non_zero_idx, :] = data.T

        self.has_average = True

    def output_histmap(self, out_file, format=None):
        if os.path.exists(out_file):
            os.remove(out_file)

        out_file_handle = open(out_file, 'a')

        x_bins, y_bins = self.bins
        out_file_handle.write('#xbins: ')
        np.savetxt(out_file_handle, x_bins, delimiter='\t', newline="\t")
        out_file_handle.write('\n#ybins: ')
        np.savetxt(out_file_handle, y_bins, delimiter='\t', newline="\t")
        out_file_handle.write('\n')

        self.histmap.to_csv(out_file_handle, sep="\t", float_format=format)

        out_file_handle.close()

    def draw_hist(self, site=None, log=False):
        if site is None:
            histmap = self.histmap[self.cell_nums == np.max(self.cell_nums)]
            histmap = np.sum(histmap, axis=0)
        else:
            histmap = self.histmap.loc[site, :]

        histmap = np.array(histmap).reshape(self.shape)

        from matplotlib.colors import LogNorm
        norm = LogNorm() if log else None

        sns.heatmap(data=histmap, norm=norm,
                    xticklabels=self.x_draw_bins, yticklabels=self.y_draw_bins)
        plt.ylabel('Distance to periphery')
        plt.xlabel('DNA Density')

        if site is not None:
            plt.title(site)

        plt.show()
        plt.close()
