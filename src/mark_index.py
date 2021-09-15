import pandas as pd
import numpy as np
import getopt
import glob
import sys
import os


def index_mark(argv):
    try:
        opts, args = getopt.getopt(argv[1:], "i:")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: D3 mark <mark_file> <index_file> <out_file>\n")
        # sys.stderr.write("Options:\n")
        # sys.stderr.write("  -f STR         Figure output. default: None.\n")
        return 1
    # for o, a in opts:
    #     if o == "-f":
    #         fig_output = a
    mark_file, index_file, out_file = args[:3]

    value_name = mark_file.split("/")[-1].split('.')[0]
    print(f'Index Mark: {value_name}.')

    os.system(f"bedtools intersect -a {index_file} -b {mark_file} -wa -wb"
              f" > {mark_file}.intersect.tmp")

    value_bed = pd.read_table(f"{mark_file}.intersect.tmp",
                              names=['chr', 'start', 'end', 'index', 'chr1', 'start1', 'end1', value_name])
    value_bed.drop(columns=['chr', 'start', 'end', 'chr1', 'start1', 'end1'], inplace=True)

    os.system(f"rm {mark_file}.intersect.tmp")

    value_bed = value_bed.groupby(['index']).mean()

    result_bed = pd.read_table(index_file, names=['chr', 'start', 'end', 'index'])
    result_bed = pd.concat([result_bed, value_bed], axis=1, join='outer')
    result_bed.fillna(value=0, inplace=True)

    result_bed.to_csv(out_file, sep="\t", index=False,
                      columns=['chr', 'start', 'end', 'index', value_name])


def index_marks(argv):
    try:
        opts, args = getopt.getopt(argv[1:], "i:")
    except getopt.GetoptError as err:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1
    if len(args) == 0:
        sys.stderr.write("Usage: D3 marks <mark_dir> <index_file> <out_file>\n")
        # sys.stderr.write("Options:\n")
        # sys.stderr.write("  -f STR         Figure output. default: None.\n")
        return 1
    # for o, a in opts:
    #     if o == "-f":
    #         fig_output = a
    mark_dir, index_file, out_file = args[:3]

    value_files = glob.glob(f'{mark_dir}/*')

    value_names = [value_file.split("/")[-1].split('.')[0] for value_file in value_files]
    result_bed = pd.read_table(index_file, names=['chr', 'start', 'end', 'index'])
    i = 0
    for value_file, value_name in zip(value_files, value_names):
        i += 1
        print(f'Read value: {i}/{len(value_files)}. {value_name}')
        os.system(f"bedtools intersect -a {index_file} -b {value_file} -wa -wb > {value_file}.intersect.tmp")

        value_bed = pd.read_table(f"{value_file}.intersect.tmp",
                                  names=['chr', 'start', 'end', 'index', 'chr1', 'start1', 'end1', value_name])
        value_bed.drop(columns=['chr', 'start', 'end', 'chr1', 'start1', 'end1'], inplace=True)

        os.system(f"rm {value_file}.intersect.tmp")

        value_bed = value_bed.groupby(['index']).mean()

        result_bed = pd.concat([result_bed, value_bed], axis=1, join='outer')
        result_bed.fillna(value=0, inplace=True)

    result_bed.to_csv(out_file, sep="\t", index=False)

    return result_bed
