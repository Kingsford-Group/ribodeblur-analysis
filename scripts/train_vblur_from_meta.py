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

def meta_pipeline(tlist, cds_range, istart, istop, rlen_min, rlen_max, converge_cutoff, obj_pdf):
    # train vblur on meta profile
    tid_select = filter_transcript_by_length(cds_range, istop)
    pos_hist = create_rlen_meta_profile(tlist, cds_range, tid_select, istart, istop)
    mobs = get_cobs(pos_hist, rlen_min, rlen_max, 0, istop)    
    estep = True
    b, ptrue, eps = train_vblur_from_meta_profiles(mobs, low, percentile, converge_cutoff, estep, obj_pdf)
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
    b, ptrue, eps = meta_pipeline(tlist, cds_range, utr5_offset, imax, rlen_min, rlen_max, converge_cutoff, odir+"train_vblur_trial.pdf")
    plot_vblur(b, odir+"vblur_single.pdf")
    write_vblur(b, odir+get_file_core(hist_fn)+".vblur")
    # # true signal deblur and plot
    # mtrue = estimate_ctrue(ptrue, eps, mobs)
    # for rlen in mobs:
    #     mobs_estimate = estimate_pobs_single(b[rlen], rlen-28, mtrue[rlen])
    #     fn = odir+"deblur_meta_rlen_{0}.pdf".format(rlen)
    #     plot_profile(mobs[rlen], ptrue * sum(mobs[rlen]), mobs_estimate, mtrue[rlen], rlen)

if __name__ == "__main__": main()


        

