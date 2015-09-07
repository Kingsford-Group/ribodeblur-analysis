#!/usr/bin/env python
import numpy as np
import sys
from multiprocessing import Pool
from meta_profile import *
from deblur_utils import *
from deblur_result_io import *
from global_params import *

def plot_vblur(vblur, ofname=None):
    plot_cnt = len(vblur)
    j = 1
    fig = plt.figure(figsize=(10,13))
    rlen_list = sorted(vblur.keys())
    for rlen in rlen_list:
        ax = fig.add_subplot(6,3,j)
        plt.plot(vblur[rlen], 'm-o', markeredgecolor='none')
        plt.title("{0}".format(rlen))
        j += 1
    plt.tight_layout()
    if ofname !=None:
        plt.savefig(ofname, bbox_inches='tight')
        plt.close()
    else:
        plt.show()

def get_max_frame_percent(vec, percentile=95):
    fsum = np.zeros(3)
    for i in xrange(3):
        vsub = vec[i::3]
        threshold = np.percentile(vsub, percentile)
        fsum[i] = np.sum(vsub[vsub<=threshold])
    return max(fsum)/float(sum(fsum))

def get_min_mean_frame_cnt(vec, percentile=95):
    fmean = np.zeros(3)
    for i in xrange(3):
        vsub = vec[i::3]
        threshold = np.percentile(vsub, percentile)
        fmean[i] = np.mean(vsub[vsub<=threshold])
    return min(fmean)

def get_vblur_rlen_range(mobs):
    vrlen_min = rlen_min
    vrlen_max = rlen_max
    for rlen in xrange(rlen_min, rlen_max+1):
        if rlen not in mobs: continue
        vrlen_min = rlen
        min_frame_cnt = get_min_mean_frame_cnt(mobs[rlen],100)
        max_frame_portion = get_max_frame_percent(mobs[rlen], 100)
        if min_frame_cnt >= lowest_frame_cnt and max_frame_portion >= lowest_frame_percent:
            break
    for rlen in xrange(rlen_max, rlen_min, -1):
        if rlen not in mobs: continue
        vrlen_max = rlen
        min_frame_cnt = get_min_mean_frame_cnt(mobs[rlen],100)
        max_frame_portion = get_max_frame_percent(mobs[rlen], 100)
        # print "{0} {1:.2%} {2:.0f}".format(rlen, max_frame_portion, min_frame_cnt)
        if min_frame_cnt >= lowest_frame_cnt and max_frame_portion >= lowest_frame_percent:
            break
    # for rlen in xrange(rlen_min, rlen_max+1):
    #     if rlen not in mobs: continue
    #     min_frame_cnt = get_min_mean_frame_cnt(mobs[rlen],100)
    #     max_frame_portion = get_max_frame_percent(mobs[rlen], 100)
    #     print "read length: {0} min frame cnt: {1:.0f} frame skew: {2:.2%}".format(rlen, min_frame_cnt, max_frame_portion)
    return vrlen_min, vrlen_max

def meta_pipeline(tlist, cds_range, istart, istop, rlen_min, rlen_max, converge_cutoff, obj_pdf):
    # train vblur on meta profile
    tid_select = filter_transcript_by_length(cds_range, istop)
    pos_hist = create_rlen_meta_profile(tlist, cds_range, tid_select, istart, istop)
    mobs = get_cobs(pos_hist, rlen_min, rlen_max, 0, istop)    
    vrlen_min, vrlen_max = get_vblur_rlen_range(mobs)
    mobs_hc = { rlen:mobs[rlen] for rlen in xrange(vrlen_min, vrlen_max+1) } 
    estep = True
    b, ptrue, eps = train_vblur_from_meta_profiles(mobs_hc, klist, low, percentile, converge_cutoff, estep, obj_pdf)
    return b, ptrue, eps
    
def main():
    if len(sys.argv) != 4:
        print "Usage: python train_vblur_from_meta.py input_rlen.hist cds_range.txt output_dir"
        exit(1)
    hist_fn = sys.argv[1]
    cds_txt = sys.argv[2]
    odir = sys.argv[3]
    ensure_dir(odir)
    cds_range = get_cds_range(cds_txt)
    tlist = parse_rlen_hist(hist_fn)
    b, ptrue, eps = meta_pipeline(tlist, cds_range, utr5_offset, imax, rlen_min, rlen_max, converge_cutoff, odir+get_file_core(hist_fn)+"_train_vblur_trial.pdf")
    plot_vblur(b, odir+get_file_core(hist_fn)+"_vblur_single.pdf")
    write_vblur(b, odir+get_file_core(hist_fn)+".vblur")
    # # true signal deblur and plot
    # mtrue = estimate_ctrue(ptrue, eps, mobs)
    # for rlen in mobs:
    #     mobs_estimate = estimate_pobs_single(b[rlen], rlen-28, mtrue[rlen])
    #     fn = odir+"deblur_meta_rlen_{0}.pdf".format(rlen)
    #     plot_profile(mobs[rlen], ptrue * sum(mobs[rlen]), mobs_estimate, mtrue[rlen], rlen)

if __name__ == "__main__": main()

