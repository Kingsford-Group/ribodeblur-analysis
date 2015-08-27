#!/usr/bin/env python
import sys
from global_params import *
from deblur_utils import *
from deblur_result_io import *
from synthetic_frameshift import read_shift_points
from footprint_hist_parser import get_transcript_profiles
from profile_evaluation import construct_all_cobs, get_frame_portion

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.size'] = 16
rcParams['xtick.major.size'] = 5
rcParams['ytick.major.size'] = 5
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42

if __name__ == "__main__":
    if len(sys.argv) != 7:
        print "Usage: python frameshift_evaluation.py rlen.hist cds_range.txt rlen.vblur rlen.eps shiftpoints.fsp output_dir"
        exit(1)
    hist_fn = sys.argv[1]
    cds_txt = sys.argv[2]
    vblur_fname = sys.argv[3]
    deblur_fname = sys.argv[4]
    shiftpoint_fname = sys.argv[5]
    odir = sys.argv[6]
    ensure_dir(odir)
    fname = get_file_core(hist_fn)
    print "get cds range"
    cds_range = get_cds_range(cds_txt)
    print "parse read len hist file"
    tlist = parse_rlen_hist(hist_fn)
    print "get pre-computed blur vector"
    b = read_vblur(vblur_fname)
    print "get pre-computed deblur results"
    ptrue, eps = read_essentials(deblur_fname)
    shift_points = read_shift_points(shiftpoint_fname)
    print "construct cobs all at once"
    tprofile = get_transcript_profiles(tlist, cds_range, utr5_offset, utr3_offset)
    cobs = construct_all_cobs(tprofile, cds_range, utr5_offset, utr3_offset, rlen_min, rlen_max)

    tid_list = np.array(cobs.keys())
    tlen = np.array([ cds_range[tid][1]-cds_range[tid][0] for tid in tid_list])
    order = np.argsort(tlen)
    tid_list = tid_list[order]
    i = 0
    f0_before = []
    f0_after = []
    f1_before = []
    f1_after = []
    start = utr5_offset
    for tid in tid_list:
        cobs_merge = merge_profiles(cobs[tid])
        ctrue = estimate_ctrue(ptrue[tid], eps[tid], cobs[tid])
        ctrue_merge = merge_profiles(ctrue)
        cds_begin, cds_end  = cds_range[tid]
        ifs = shift_points[tid] - start
        if np.any(cobs_merge[:ifs]!=0):
            f0_before.append(get_frame_portion(cobs_merge[:ifs],0))
            f0_after.append(get_frame_portion(ctrue_merge[:ifs],0))
        if np.any(cobs_merge[ifs:]!=0):
            f1_before.append(get_frame_portion(cobs_merge[ifs:],1))
            f1_after.append(get_frame_portion(ctrue_merge[ifs:],1))

        i += 1
        if i<20:
            stop = (cds_end-cds_begin)+utr3_offset
            x_pos = np.arange(start, stop+1)
            start_in_frame = start - start%3
            xgrid = np.arange(start_in_frame - 0.4, stop+1, 3)
            figwidth = len(ctrue_merge)/10.0*1.5
            fig = plt.figure(figsize=(figwidth, 13))
            ax = fig.add_subplot(2,1,1)
            plt.bar(x_pos-0.4, cobs_merge, width=0.8, color='b', edgecolor='white', alpha=0.3)
            yheight = max(cobs_merge)+1 
            plt.bar(xgrid, [yheight]*len(xgrid), width=0.8, color='r', edgecolor='white', alpha=0.1)
            plt.bar(shift_points[tid]-0.4, yheight, width=0.8, color='r', edgecolor='white', alpha=0.6)
            plt.ylim((0,yheight))
            plt.xlim((x_pos[0]-1, x_pos[-1]+1))
            plt.title("Before deblur".format(tid))
            ax = fig.add_subplot(2,1,2)
            plt.bar(x_pos-0.4, ctrue_merge, width=0.8, color='b', edgecolor='white', alpha=0.3)
            yheight = max(ctrue_merge)+1 
            plt.bar(xgrid, [yheight]*len(xgrid), width=0.8, color='r', edgecolor='white', alpha=0.1)
            plt.bar(shift_points[tid]-0.4, yheight, width=0.8, color='r', edgecolor='white', alpha=0.6)
            plt.ylim((0,yheight))
            plt.xlim((x_pos[0]-1, x_pos[-1]+1))
            plt.title("After deblur")
            plt.tight_layout()
            plt.savefig("{0}{1}_fs.pdf".format(odir, tid), bbox_inches='tight')
    fig = plt.figure(figsize=(12,6))
    ax = fig.add_subplot(1,2,1)
    ns, bins, patches = plt.hist(f0_before, 50, histtype='stepfilled', color='c', edgecolor='c', alpha=0.4)
    ns, bins, patches = plt.hist(f0_after, 50, histtype = 'stepfilled', color='b', edgecolor='b', hatch='/', alpha=0.3)
    plt.xlabel("frame 0 percentage before frameshift point")
    plt.ylabel("number of transcripts")
    plt.legend(["before deblur", "after deblur"], loc=2, frameon=False, fontsize=16)
    print "improve rate: {0:.2%}".format(np.mean(np.array(f0_before)<np.array(f0_after)))
    ax = fig.add_subplot(1,2,2)
    ns, bins, patches = plt.hist(f1_before, 50, histtype='stepfilled', color='c', edgecolor='c', alpha=0.4)
    ns, bins, patches = plt.hist(f1_after, 50, histtype = 'stepfilled', color='b', edgecolor='b', hatch='/', alpha=0.3)
    plt.xlabel("frame 1 percentage after frameshift point")
    print "improve rate: {0:.2%}".format(np.mean(np.array(f1_before)<np.array(f1_after)))
    plt.tight_layout()
    plt.savefig("{0}fs_portion.png".format(odir), bbox_inches='tight')
                                         

    
