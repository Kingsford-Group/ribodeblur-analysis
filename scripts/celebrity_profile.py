#!/usr/bin/env python
import sys
import numpy as np
from global_params import *
from deblur_utils import *
from deblur_result_io import *
from meta_profile import get_frame_str_vec
from footprint_hist_parser import get_cds_range, parse_rlen_hist, get_transcript_profiles
from profile_evaluation import get_frame_portion, construct_all_cobs

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.size'] = 20
rcParams['xtick.major.size'] = 5
rcParams['ytick.major.size'] = 5
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42

def batch_distance(tid_list, cobs, ctrue):
    dist = []
    for tid in tid_list:
        obs_merge = merge_profiles(cobs[tid])
        true_merge = merge_profiles(ctrue[tid])
        diff = obs_merge - true_merge
        dist.append( np.linalg.norm(diff)/sum(obs_merge) )
    return np.array(dist)

def get_rlen_portion(read_list, rlen=28):
    rlen_cnt = {}
    for rl, pos_list in read_list.iteritems():
        crlen = 0
        for pos, cnt in pos_list:
            crlen += cnt
        rlen_cnt[rl] = crlen
    ctot = sum(rlen_cnt.values())
    return float(rlen_cnt[rlen])/ctot

def batch_rlen_portion(tid_list, tprofile, rlen=28):
    return np.array([ get_rlen_portion(tprofile[tid], rlen) for tid in tid_list])

def plot_transcript(plist, start, end, tid, ofname=None):
    rlen_list = sorted(plist.keys())
    plot_cnt = len(rlen_list)+1
    x_pos = np.arange(start, end+1)
    start_in_frame = start - start%3
    xgrid = np.arange(start_in_frame-0.4, end+1, 3)
    y_tot = np.zeros(len(x_pos))
    figwidth = len(x_pos)/10.0*1.5
    fig = plt.figure(figsize=(figwidth, 13))
    for i in xrange(plot_cnt-1):
        y_cnt = plist[rlen_list[i]]
        y_tot += y_cnt
        ax = fig.add_subplot(plot_cnt, 1, i+1)
        plt.bar(x_pos-0.4, y_cnt, width=0.8, color='b', edgecolor='white', alpha=0.3)
        frame_str = get_frame_str_vec(y_cnt)
        plt.title("{0}\n{1}".format(rlen_list[i], frame_str))
        yheight = [max(y_cnt)+1]*len(xgrid)
        plt.bar(xgrid, yheight, width=0.8, color='r', edgecolor='white', alpha=0.1)
        plt.xticks(np.arange(start-1,end+1, 10))
        plt.xlim((start-1, end+1))
        plt.ylim((0,max(y_cnt)))
    ax = fig.add_subplot(plot_cnt, 1, plot_cnt)
    plt.bar(x_pos-0.4, y_tot, width=0.8, color='b', edgecolor='white', alpha=0.3)
    frame_str = get_frame_str_vec(y_tot)
    plt.title("merged\n{0}".format(frame_str))
    yheight = [max(y_tot)+1]*len(xgrid)
    plt.bar(xgrid, yheight, width=0.8, color='r', edgecolor='white', alpha=0.1)
    plt.xticks(np.arange(start-1,end+1, 10))
    plt.xlim((start-1, end+1))
    plt.ylim((0,max(y_tot)))
    plt.tight_layout()
    if ofname!=None:
        plt.savefig(ofname, bbox_inches='tight')
        plt.close()
    else:
        plt.show()

def floppy_plot_transcript(plist, start, end, tid, ofname=None):
    rlen_list = sorted(plist.keys())
    plot_cnt = len(rlen_list)+1
    x_pos = np.arange(start, end+1)
    start_in_frame = start - start%3
    xgrid = np.arange(start_in_frame-0.4, end+1, 3)
    y_tot = np.zeros(len(x_pos))
    fig = plt.figure()
    for i in xrange(plot_cnt-1):
        y_cnt = plist[rlen_list[i]]
        y_tot += y_cnt
        ax = fig.add_subplot(plot_cnt, 1, i+1)
        plt.bar(x_pos-0.4, y_cnt, width=0.8, color='b', edgecolor='white', alpha=0.6)
        frame_str = get_frame_str_vec(y_cnt)
        if i == 0:
            plt.title("{0}\n{1}\n{2}".format(tid, rlen_list[i], frame_str))
        else:
            plt.title("{0}\n{1}".format(rlen_list[i], frame_str))
        yheight = [max(y_cnt)+1]*len(xgrid)
        #plt.bar(xgrid, yheight, width=0.8, color='r', edgecolor='white', alpha=0.1)
        # plt.xticks(np.arange(start-1,end+1, 10))
        plt.xlim((start-1, end+1))
        plt.ylim((0,max(y_cnt)))
    ax = fig.add_subplot(plot_cnt, 1, plot_cnt)
    plt.bar(x_pos-0.4, y_tot, width=0.8, color='b', edgecolor='white', alpha=0.6)
    frame_str = get_frame_str_vec(y_tot)
    plt.title("merged\n{0}".format(frame_str))
    yheight = [max(y_tot)+1]*len(xgrid)
    #plt.bar(xgrid, yheight, width=0.8, color='r', edgecolor='white', alpha=0.1)
    # plt.xticks(np.arange(start-1,end+1, 10))
    plt.xlim((start-1, end+1))
    plt.ylim((0,max(y_tot)))
    #plt.tight_layout()
    if ofname!=None:
        plt.savefig(ofname, bbox_inches='tight')
        plt.close()
    else:
        plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print "Usage: python frameshift_celebrity.py rlen.hist cds_range.txt rlen.vblur rlen.eps output_dir"
        exit(1)
    hist_fn = sys.argv[1]
    cds_txt = sys.argv[2]
    vblur_fname = sys.argv[3]
    deblur_fname = sys.argv[4]
    odir = sys.argv[5]
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
    print "construct cobs all at once"
    tprofile = get_transcript_profiles(tlist, cds_range, utr5_offset, utr3_offset)
    #cobs = build_cobs_with_shifts(tprofile, cds_range, utr5_offset, utr3_offset, rlen_min, rlen_max, klist)
    cobs = construct_all_cobs(tprofile, cds_range, utr5_offset, utr3_offset, rlen_min, rlen_max)
    print "construct ctrue all at once"
    ctrue = batch_build_ctrue(ptrue, eps, cobs)
    tid_list = np.array(cobs.keys())

    # # different orders of transcripts
    # # frame porition
    # frame_list = np.array([ get_frame_portion(merge_profiles(ctrue[tid])) for tid in tid_list])
    # order = np.argsort(frame_list)
    # # transcript length
    # tlen_list = np.array([ len(ptrue[tid]) for tid in tid_list ])
    # order = np.argsort(tlen_list)
    # # euclidean distance between cobs and ctrue divided by transcript length
    # dis_list = batch_distance(tid_list, cobs, ctrue)
    # dis_list = dis_list/tlen_list
    # order = np.argsort(dis_list)[::-1]
    # 28-mer portion
    # good_rlen_portion = batch_rlen_portion(tid_list, tprofile, 28)
    # order = np.argsort(good_rlen_portion)

    print "plotting"
    for tid in ['YOR302W']:#tid_list[order][:10]:
        start, end = cds_range[tid]
        plot_transcript(cobs[tid], utr5_offset, (end-start)+utr3_offset, tid, "{0}{1}_before.pdf".format(odir,tid))        
        plot_transcript(ctrue[tid], utr5_offset, (end-start)+utr3_offset, tid, "{0}{1}_after.pdf".format(odir,tid))                
    
