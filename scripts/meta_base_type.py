#!/usr/bin/env python
import numpy as np
import sys
from footprint_hist_parser import *
from global_params import *
from meta_profile import *

def create_start_base(tlist, cds_range, tid_included, ibegin, iend, tseq):
    """
    create count accumulation for each base location with different read length
    exclude transcript if id not in tid_included (a set)
    ibegin: index bound to include before START (negative)
    iend: number of bases to include after START (positive)
    """
    print "\ncreate meta hist"
    meta_hist = {}
    for rid in tlist:
        tid = tlist[rid]['tid']
        if tid not in tid_included: continue
        start, end = cds_range[tid]
        for rlen in tlist[rid]['prof']:
            meta_hist.setdefault(rlen, {})
            for pos, cnt in tlist[rid]['prof'][rlen]:
                i = pos - start
                if i < ibegin or i > iend : continue
                base = tseq[tid][pos]
                meta_hist[rlen].setdefault(base, {})
                meta_hist[rlen][base].setdefault(i,0)
                meta_hist[rlen][base][i] += cnt
        sys.stdout.write("processed transcript {0}.\t\r".format(rid))
        sys.stdout.flush()
    sys.stdout.write("\n")
    return meta_hist

def create_stop_base(tlist, cds_range, tid_included, ibegin, iend, tseq):
    """
    create count accumulation for each base location with different read length
    exclude transcript if id not in tid_included (a set)
    ibegin: index bound to include before START (negative)
    iend: number of bases to include after START (positive)
    """
    print "\ncreate meta hist"
    meta_hist = {}
    for rid in tlist:
        tid = tlist[rid]['tid']
        if tid not in tid_included: continue
        start, end = cds_range[tid]
        for rlen in tlist[rid]['prof']:
            meta_hist.setdefault(rlen, {})
            for pos, cnt in tlist[rid]['prof'][rlen]:
                i = pos + rlen - start
                if i < ibegin or i > iend : continue
                base = tseq[tid][pos+rlen]
                meta_hist[rlen].setdefault(base, {})
                meta_hist[rlen][base].setdefault(i,0)
                meta_hist[rlen][base][i] += cnt
        sys.stdout.write("processed transcript {0}.\t\r".format(rid))
        sys.stdout.flush()
    sys.stdout.write("\n")
    return meta_hist

def plot_base_pos_hist(meta_hist, ibegin, iend, fn_prefix):
    print "plotting meta pos hist..."
    x_pos = np.arange(ibegin, iend+1)
    figwidth = len(x_pos)/10.0*1.5
    plt.figure(figsize=(figwidth,6))
    plt.rc('axes', color_cycle=['b', 'c', 'r', 'm'])
    for rlen in meta_hist:
        plt.cla()
        rcParams['figure.figsize'] = figwidth, 6
        pbase = {}
        y_cnt = np.zeros(len(x_pos))
        for base in meta_hist[rlen]:
            pbase[base] = get_pos_hist(meta_hist[rlen][base], ibegin, iend)
            y_cnt += pbase[base]
        if np.any(y_cnt == 0): continue
        print rlen,
        for base in pbase:
            print base,
            ycnt_base = pbase[base]/y_cnt
            plt.plot(x_pos, ycnt_base, '-o', markeredgecolor='none', alpha=0.5)
        print "\n",
        plt.xlabel("offset from start")
        plt.ylabel("nucleotide portion")
        plt.xticks(range(-100, iend+1, 5))
        plt.xlim((ibegin-1, iend+1))
        plt.savefig("{0}_{1}.pdf".format(fn_prefix, rlen), bbox_inches="tight")
        plt.close()
    return 

#=============================
# main
#=============================
def plot_rlen_hist_pipe():
    from deblur_result_io import ensure_dir, get_file_core
    if len(sys.argv) != 5:
        print "Usage: python meta_base_type.py input_rlen.hist ref.fa cds_range output_dir"
        exit(1)
    hist_fn = sys.argv[1]
    ref_fn = sys.argv[2]
    cds_txt = sys.argv[3]
    odir = sys.argv[4]
    min_sample_cnt = 100
    imax = 350
    ensure_dir(odir)
    cds_range = get_cds_range(cds_txt)
    tid_select = filter_transcript_by_length(cds_range, imax)
    tlist = parse_rlen_hist(hist_fn)
    rlen2cnt = get_length_count(tlist)
    rlen2portion = get_length_distribution(rlen2cnt, rlen_min, rlen_max)
    tot_portion = 0
    for rlen in sorted(rlen2portion.keys()):
        print "{0}-mers: {1:.2%}".format(rlen, rlen2portion[rlen])
        tot_portion += rlen2portion[rlen]
    print "total: {0:.2%}".format(tot_portion)
    sframe = get_frame_str_from_tlist(tlist, cds_range)
    tseq = get_tseq(ref_fn, None)
    meta_hist = create_start_base(tlist, cds_range, tid_select, utr5_offset, imax, tseq)
    #meta_hist = create_stop_base(tlist, cds_range, tid_select, utr5_offset, imax, tseq)
    fn_prefix = odir+"/"+get_file_core(hist_fn)+'_start_base'
    plot_base_pos_hist(meta_hist, utr5_offset, imax, fn_prefix)

if __name__ == "__main__": plot_rlen_hist_pipe()
