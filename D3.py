# !/usr/bin/env python

import sys


def main():
    if len(sys.argv) == 1:
        sys.stderr.write("Usage: D3 <command> <arguments>\n")
        sys.stderr.write("Commands:\n")
        sys.stderr.write("  D3       Compute DNA density and DisTP\n")
        sys.stderr.write("  D3s       Compute DNA density and DisTP\n")
        sys.stderr.write("  sta      Histogram density and DisTP, show the border of outliers \n")
        sys.stderr.write("  map      Put genomic bins into density-DisTP matrix\n")
        sys.stderr.write("  ave      Get the mean and sd of density and DisTP\n")
        sys.stderr.write("\n")
        sys.stderr.write("  mark     Index genomic marker\n")
        sys.stderr.write("  marks     Index genomic marker\n")
        sys.stderr.write("  enrich   Enrich marker at density-DisTP matrix\n")
        sys.stderr.write("  hiera    Hierarchy cluster density states by markers\n")
        # sys.stderr.write("\n")
        # sys.stderr.write("  gene     Genes enrichment at matrix\n")
        # sys.stderr.write("  move     Genes movement from type1 to type2\n")
        return 1

    if sys.argv[1] == "D3":
        from src.den_dtp import cal_den_dtp
        cal_den_dtp(sys.argv[1:])
    elif sys.argv[1] == "D3s":
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
    else:
        sys.stderr.write("[E::" + __name__ + "] unknown command\n")
        return 1


if __name__ == "__main__":
    main()
