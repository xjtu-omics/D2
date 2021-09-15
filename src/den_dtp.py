from src.readfile import Dg, save_plot_pdf
import os
import sys
import glob
import getopt
import pickle
import warnings
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import ndimage
from numpy.linalg import norm
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
from scipy.stats import entropy


###################
# Assign points
def _assign_points(value_map, location_dict, point_num, value=None):
    out_array = np.zeros(point_num)

    points_pos = np.array(np.where(value_map)).T

    for i in range(points_pos.shape[0]):
        location = tuple(points_pos[i, :])

        if location in location_dict:
            if value is None:
                out_array[np.array(location_dict[location])] = value_map[location]
            else:
                out_array[np.array(location_dict[location])] = value

    return out_array


###################
# Get Density Map
def _define_bin_size(points, dis_number=7):
    ran_size = 1000

    point_matrix = np.array(points[['x', 'y', 'z']])

    dis_percent = dis_number * 100 / points.shape[0]

    # If the matrix size is too large, randomly sample part of it for bin_size calculation.
    if point_matrix.shape[0] > ran_size:
        ran_index = np.random.choice(point_matrix.shape[0] - 1, size=ran_size,
                                     replace=False)
        sub_matrix = point_matrix[ran_index, :]
    else:
        sub_matrix = point_matrix

    distance_matrix = cdist(sub_matrix, point_matrix, metric='euclidean')

    # Substitute zero distance (point to itself) to inf.
    # distance_matrix[distance_matrix == 0] = np.inf

    bin_size = np.percentile(distance_matrix, dis_percent, axis=1)
    bin_size = np.mean(bin_size)

    return bin_size


def _cal_density_map(points, bin_size):
    # Divide 3D space into evenly distributed bins.
    x_range = [min(points['x']), max(points['x'])]
    x_num = int((x_range[1] - x_range[0]) / bin_size) + 1

    y_range = [min(points['y']), max(points['y'])]
    y_num = int((y_range[1] - y_range[0]) / bin_size) + 1

    z_range = [min(points['z']), max(points['z'])]
    z_num = int((z_range[1] - z_range[0]) / bin_size) + 1

    # extended_bin_num = int(np.min([x_num, y_num, z_num]) / 2) + 1
    extended_bin_num = 0

    x_loc = (np.array(points['x']) - x_range[0]) / bin_size
    y_loc = (np.array(points['y']) - y_range[0]) / bin_size
    z_loc = (np.array(points['z']) - z_range[0]) / bin_size

    x_loc = x_loc.astype(int) + extended_bin_num
    y_loc = y_loc.astype(int) + extended_bin_num
    z_loc = z_loc.astype(int) + extended_bin_num
    points['x_loc'] = x_loc
    points['y_loc'] = y_loc
    points['z_loc'] = z_loc

    locations = list(zip(*(x_loc, y_loc, z_loc)))

    density_map = np.zeros((x_num + extended_bin_num * 2,
                            y_num + extended_bin_num * 2,
                            z_num + extended_bin_num * 2),
                           dtype=int)
    location_dict = {}
    for i in range(len(locations)):
        location = locations[i]

        if location in location_dict:
            location_dict[location].append(i)

        else:
            location_dict[location] = [i]

        density_map[location] += 1

    return density_map, location_dict


###################
# Get Outside Periphery
def _find_out_one_row(row, k_length=5, k_num=3):
    ###################
    # Left
    # Make sure there are at least k_num points within k_length from peri
    k = np.zeros(k_length * 2 + 1)
    k[:k_length + 1] = 1

    edge = ndimage.convolve(row, k, mode='constant', cval=0)
    left_edges = edge >= k_num

    # Make sure the point next to peri must not be empty
    is_filled = np.array(list(row[1:] > 0) + [0])
    left_edges = np.logical_and(left_edges, is_filled)

    if not np.any(np.sum(left_edges)):
        return np.ones_like(row), np.zeros_like(row)

    left_edge = np.min(np.where(left_edges)[0])

    ###################
    # Right
    k = np.zeros(k_length * 2 + 1)
    k[k_length:] = 1

    edge = ndimage.convolve(row, k, mode='constant', cval=0)
    right_edges = edge >= k_num

    # Make sure the point next to peri must not be empty
    is_filled = np.array([0] + list(row[:-1] > 0))
    right_edges = np.logical_and(right_edges, is_filled)

    if not np.any(np.sum(right_edges)):
        return np.ones_like(row), np.zeros_like(row)

    right_edge = np.max(np.where(right_edges)[0])

    ###################
    # out_points
    new_row = np.zeros_like(row)
    new_row[:left_edge + 1] = 1
    new_row[right_edge:] = 1

    # edges
    edge_row = np.zeros_like(row)
    edge_row[left_edge] = 1
    edge_row[right_edge] = 1

    return new_row, edge_row


def _find_out(density_map):
    new_density_map = density_map.copy()

    # y-z
    density_map_yz = np.zeros_like(new_density_map)
    edge_map_yz = np.zeros_like(new_density_map)
    for z in range(new_density_map.shape[2]):
        for y in range(new_density_map.shape[1]):
            row = new_density_map[:, y, z]
            # row[row != 0] = 1
            row, edge_row = _find_out_one_row(row)
            density_map_yz[:, y, z] = row
            edge_map_yz[:, y, z] = edge_row

    # x-y
    density_map_xy = np.zeros_like(new_density_map)
    edge_map_xy = np.zeros_like(new_density_map)
    for x in range(new_density_map.shape[0]):
        for y in range(new_density_map.shape[1]):
            row = new_density_map[x, y, :]
            # row[row != 0] = 1
            row, edge_row = _find_out_one_row(row)
            density_map_xy[x, y, :] = row
            edge_map_xy[x, y, :] = edge_row

    # x-z
    density_map_xz = np.zeros_like(new_density_map)
    edge_map_xz = np.zeros_like(new_density_map)
    for z in range(new_density_map.shape[2]):
        for x in range(new_density_map.shape[0]):
            row = new_density_map[x, :, z]
            # row[row != 0] = 1
            row, edge_row = _find_out_one_row(row)
            density_map_xz[x, :, z] = row
            edge_map_xz[x, :, z] = edge_row

    out_map = np.logical_or(np.logical_or(density_map_xy, density_map_yz), density_map_xz)
    edge_map = np.logical_or(np.logical_or(edge_map_xy, edge_map_yz), edge_map_xz)

    return out_map, edge_map


###################
# Smooth
def _smooth_density_map(density_map, out_map, k_dimension=3):
    norm_func = lambda x: 1 / x if x > 0 else 2

    ############
    # Get smoothed density density_map
    k_radius = int(k_dimension / 2)
    k_mtx = np.zeros((k_dimension, k_dimension, k_dimension))

    for i in range(k_dimension):
        for j in range(k_dimension):
            for k in range(k_dimension):
                dis = norm((i - k_radius, j - k_radius, k - k_radius))
                k_mtx[i, j, k] = norm_func(dis)

    smooth_map = ndimage.convolve(density_map, k_mtx, mode='constant', cval=0)
    smooth_map = np.array(smooth_map).astype(float)

    ############
    # Get valid points number
    inside_points = np.logical_not(out_map).astype(int)
    point_counts = ndimage.convolve(inside_points, k_mtx, mode='constant', cval=0)
    point_counts = np.array(point_counts).astype(float)

    ############
    # Get average density
    non_zero_index = point_counts > 0
    smooth_map[non_zero_index] = np.nan_to_num(smooth_map[non_zero_index] / point_counts[non_zero_index])
    smooth_map[out_map] = 0

    return smooth_map


###################
# Get distance to periphery
def _get_dis_to_out(out_map, edge_map, dis_num=10):
    points_index = np.where(np.logical_not(out_map))
    points_pos = np.array(points_index).T
    edge_pos = np.array(np.where(edge_map)).T

    dis_to_out = cdist(points_pos, edge_pos)
    new_dis_to_out = np.zeros(dis_to_out.shape[0])
    dis_per = dis_num / dis_to_out.shape[1]
    for i in range(dis_to_out.shape[0]):
        da_ar = dis_to_out[i, :]
        new_dis_to_out[i] = (np.min(da_ar) + np.quantile(da_ar, dis_per)) / 2
        # new_dis_to_out[i] = np.mean(da_ar[da_ar <= out_dis_range + min_dis])
    del dis_to_out

    dis_to_out_map = np.zeros_like(out_map, dtype=float)
    dis_to_out_map[points_index] = new_dis_to_out

    return dis_to_out_map


###################
# Get intermingle
def _get_intermingle(points, map_shape):
    inter_range = 3  # Must be odd
    half_range = (inter_range - 1) / 2

    mingle_map = np.zeros(map_shape)

    inter_counts = points.dg_df.groupby(['x_loc', 'y_loc', 'z_loc', 'chr']).size()
    inter_counts = inter_counts.reset_index()
    for x_idx in range(inter_range):
        sub_x_idx = x_idx - half_range
        x_max = int((mingle_map.shape[0] - 1 - x_idx) / inter_range) * inter_range + half_range

        inter_counts_x = inter_counts[(inter_counts['x_loc'] >= sub_x_idx) &
                                      (inter_counts['x_loc'] <= x_max)].copy()
        inter_counts_x['tran_x'] = ((inter_counts_x['x_loc'] - sub_x_idx) / inter_range).astype(int)

        for y_idx in range(inter_range):
            sub_y_idx = y_idx - half_range
            y_max = int((mingle_map.shape[1] - 1 - y_idx) / inter_range) * inter_range + half_range

            inter_counts_xy = inter_counts_x[(inter_counts_x['y_loc'] >= sub_y_idx) &
                                             (inter_counts_x['y_loc'] <= y_max)].copy()
            inter_counts_xy['tran_y'] = ((inter_counts_xy['y_loc'] - sub_y_idx) / inter_range).astype(int)

            for z_idx in range(inter_range):
                sub_z_idx = z_idx - half_range
                z_max = int((mingle_map.shape[2] - 1 - z_idx) / inter_range) * inter_range + half_range

                inter_counts_xyz = inter_counts_xy[(inter_counts_xy['z_loc'] >= sub_z_idx) &
                                                   (inter_counts_xy['z_loc'] <= z_max)].copy()
                inter_counts_xyz['tran_z'] = ((inter_counts_xyz['z_loc'] - sub_z_idx) / inter_range).astype(int)

                inter_counts_xyz = inter_counts_xyz.groupby(['tran_x', 'tran_y', 'tran_z', 'chr']).agg({0: 'sum'})
                inter_counts_xyz.reset_index(inplace=True)

                inter_counts_xyz = inter_counts_xyz.groupby(['tran_x', 'tran_y', 'tran_z']).agg({0: entropy})
                inter_counts_xyz.reset_index(inplace=True)

                map_x_idx = np.array(inter_counts_xyz['tran_x']) * inter_range + x_idx
                map_y_idx = np.array(inter_counts_xyz['tran_y']) * inter_range + y_idx
                map_z_idx = np.array(inter_counts_xyz['tran_z']) * inter_range + z_idx

                mingle_map[map_x_idx, map_y_idx, map_z_idx] = np.array(inter_counts_xyz.loc[:, 0])

    return mingle_map


###################
# Draw maps
def _2d_scatter(maps, map_types, map_names, nrows=2, ncols=2, fig_size=(5, 5),
                axis=0, axis_range=None, output=None):
    axis_range = [0, maps[0].shape[axis]] if axis_range is None else axis_range

    if len(maps) > nrows * ncols:
        raise ValueError(f'length of maps {len(maps)} > numer of subplots {nrows * ncols}.')

    #################
    # normalize and transfer
    new_maps = []
    for amap, map_type in zip(maps, map_types):

        if map_type == 'map':
            amap = amap.astype(float)
            amap = (amap - np.min(amap)) / (np.max(amap) - np.min(amap))
            new_maps.append(amap)

        elif map_type == 'point':
            point_map = np.ones_like(maps[0]) * -1
            if amap.shape[1] == 4:
                point_map[amap[:, 0], amap[:, 1], amap[:, 2]] = amap[:, 3]
            else:
                point_map[amap[:, 0], amap[:, 1], amap[:, 2]] = 0
            new_maps.append(point_map)

        elif map_type == 'grad':
            amap = amap.astype(float)
            amap[amap > 0] = amap[amap > 0] / np.max(amap[amap > 0])
            amap[amap < 0] = amap[amap < 0] / -np.min(amap[amap < 0])
            new_maps.append(amap)

        else:
            raise ValueError(f'Wrong map_type: {map_type}')

    #################
    # draw
    vmaxes = [np.max(amap) for amap in new_maps]

    for axis_num in range(axis_range[0], axis_range[1]):

        plt.figure(figsize=fig_size)

        n_map = 0
        for amap, map_type, map_name in zip(new_maps, map_types, map_names):
            n_map += 1
            plt.subplot(nrows, ncols, n_map)

            if axis == 0:
                plot_map = amap[axis_num, :, :]
            elif axis == 1:
                plot_map = amap[:, axis_num, :]
            elif axis == 2:
                plot_map = amap[:, :, axis_num]
            else:
                raise ValueError(f'Only Three dimensions, with input dimension:{axis}')

            if map_type == 'map':
                sns.heatmap(plot_map, vmax=1, vmin=0, cbar=False)

            elif map_type == 'grad':
                sns.heatmap(plot_map, vmax=1, vmin=-1, center=0, cmap='bwr', cbar=False)

            else:
                cluster_num = vmaxes[n_map - 1]
                if cluster_num > 0:
                    colors = sns.color_palette("Set2", cluster_num + 1)
                    colors[0] = (0, 0, 0)
                    sns.heatmap(plot_map, vmax=cluster_num, cmap=colors)
                else:
                    sns.heatmap(plot_map)

            plt.axis('off')
            plt.title(map_name)

        if output is None:
            plt.show()
        else:
            plt.savefig(output + f'_axis{axis}_{axis_num}th.png')
        plt.close()


def cal_den_dtp(argv):
    debug = 0
    write_log = 0
    try:
        opts, args = getopt.getopt(argv[1:], "w:d:")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: D3 D3 [options] <3dg_file> <index_file> <output>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -d Bool         If is debug model. default: 0.\n")
        sys.stderr.write("  -w Bool         If write to log. default: 0.\n")
        return 1
    for o, a in opts:
        if o == "-d":
            debug = int(a)
        if o == "-w":
            write_log = int(a)
    dg_file, index_file, output = args[:3]

    if write_log != 0:
        log_file = open(f'{output}.log', 'w')
        sys.stdout = log_file
        sys.stderr = log_file

    points = Dg(dg_file, input_window=index_file)
    print(f'Read file DONE: {dg_file}')

    bin_size = _define_bin_size(points)
    print('Bin size is ' + str(round(bin_size, 4)))

    density_map, location_dict = _cal_density_map(points, bin_size)

    out_map, edge_map = _find_out(density_map)
    is_peri = _assign_points(out_map, location_dict, points.shape[0], value=1)
    points['dc'] = is_peri
    density_map[out_map] = 0  # clear out points
    print(f'Get out points DONE.')

    smooth_map = _smooth_density_map(density_map, out_map)
    den_array = _assign_points(smooth_map, location_dict, points.shape[0])
    points['density'] = den_array
    print(f'Smooth density DONE.')

    dis_to_out_map = _get_dis_to_out(out_map, edge_map)
    dis_to_out_array = _assign_points(dis_to_out_map, location_dict, points.shape[0])
    points['dis_to_out'] = dis_to_out_array
    print(f'Get dis_to_out DONE.')

    points.dg_df.to_csv(f'{output}.den_dtp.txt', sep="\t", index=False, header=False,
                        columns=['chr', 'start', 'end', 'index', 'filter',
                                 'x_loc', 'y_loc', 'z_loc', 'density', 'dis_to_out'])
    print('Save result DONE')

    os.system(f'mkdir -p {output}_2d_scatter')
    _2d_scatter([density_map, out_map, smooth_map, dis_to_out_map],
                ['map', 'map', 'map', 'map'],
                ['Raw Density', 'Out Cubes', 'Smoothed Density', 'DisTP'],
                output=f"{output}_2d_scatter/{output.split('/')[-1]}")
    print('Plot result DONE')

    if debug:
        debug_variables = {'bin_size': bin_size, 'density_map': density_map,
                           'location_dic': location_dict, 'out_map': out_map, 'edge_map': edge_map,
                           'smooth_map': smooth_map, 'dis_to_out_map': dis_to_out_map,
                           'points': points}
        pickle.dump(debug_variables, open(f'{output}_debug.pkl', 'wb'))
    print('All DONE')

    if write_log:
        log_file.close()


def cal_den_dtps(argv):
    debug = 0
    write_log = 0
    try:
        opts, args = getopt.getopt(argv[1:], "d:w:")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: D3 D3s [options] <3dg_dir> <index_file> <out_dir>\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -d INT         Debug model: 1 for yes, 0 for no. default: 0.\n")
        sys.stderr.write("  -w Bool         If write to log. default: 0.\n")
        return 1
    for o, a in opts:
        if o == "-d":
            debug = int(a)
        if o == "-w":
            write_log = int(a)
    dg_dir, index_file, out_dir = args[:3]

    dg_files = glob.glob(f'{dg_dir}/*')

    os.system(f'mkdir -p {out_dir}/den_dtp')

    for dg_file in dg_files:
        dg_name = dg_file.split('/')[-1].split('.')[0]
        output = f'{out_dir}/{dg_name}'
        os.system(f'mkdir -p {output}')
        cal_den_dtp(['D3', '-d', str(debug), '-w', str(write_log),
                     dg_file, index_file, f'{output}/{dg_name}'])

        os.system(f'cp {output}/{dg_name}.den_dtp.txt {out_dir}/den_dtp/')
