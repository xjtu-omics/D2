# !/usr/bin/env python

import os
import sys

cwd = os.getcwd()
sys.path.append(cwd)


def main():
    if len(sys.argv) == 1:
        sys.stderr.write("Usage: D2 <command> <arguments>\n")
        sys.stderr.write("Commands:\n")
        sys.stderr.write("  D2       Compute DNA density and DisTP of one cell\n")
        sys.stderr.write("  D2s      Compute DNA density and DisTP of multiple cells\n")
        sys.stderr.write("  sta      Histogram density and DisTP to determine the boundaries for D2 plot\n")
        sys.stderr.write("  map      Put genomic bins into D2 plot\n")
        sys.stderr.write("  ave      Obtain the mean and SD of density and DisTP\n")
        sys.stderr.write("\n")
        sys.stderr.write("  mark     Index genomic marker\n")
        sys.stderr.write("  marks    Index genomic markers\n")
        sys.stderr.write("  enrich   Compute marker fold change on D2 plot\n")
        sys.stderr.write("  hiera    Hierarchy cluster of physical states\n")
        sys.stderr.write("\n")
        sys.stderr.write("  gene     Compute gene fold changes on D2 plot\n")
        sys.stderr.write("  act      Compute activation index\n")
        return 1

    if sys.argv[1] == "D2":
        from src.den_dtp import cal_den_dtp
        cal_den_dtp(sys.argv[1:])
    elif sys.argv[1] == "D2s":
        from src.den_dtp import cal_den_dtps
        cal_den_dtps(sys.argv[1:])
    elif sys.argv[1] == "sta":
        from src.histmap import draw_hist
        draw_hist(sys.argv[1:])
    elif sys.argv[1] == "map":
        from src.histmap import to_histmap
        to_histmap(sys.argv[1:])
    elif sys.argv[1] == "ave":
        from src.ave_std import ave_std
        ave_std(sys.argv[1:])
    elif sys.argv[1] == "mark":
        from src.mark_index import index_mark
        index_mark(sys.argv[1:])
    elif sys.argv[1] == "marks":
        from src.mark_index import index_marks
        index_marks(sys.argv[1:])
    elif sys.argv[1] == "enrich":
        from src.enrich import plot_value_hists
        plot_value_hists(sys.argv[1:])
    elif sys.argv[1] == "hiera":
        from src.enrich import hierarchy_hist
        hierarchy_hist(sys.argv[1:])
    elif sys.argv[1] == "gene":
        from src.gene_fold import gene_fold_changes
        gene_fold_changes(sys.argv[1:])
    elif sys.argv[1] == "act":
        from src.act_index import compute_activation_index
        compute_activation_index(sys.argv[1:])
    else:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1


if __name__ == "__main__":
    main()
