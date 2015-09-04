#!/usr/bin/env python
import numpy as np
import sys
from footprint_hist_parser import *
from global_params import *
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.size'] = 12
rcParams['xtick.major.size'] = 5
rcParams['ytick.major.size'] = 5
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42

#=============================
# utils
#=============================
def unzip_dic_tuples(dic):
    keys, vals = zip(*dic.items())
    order = np.argsort(keys)
    return np.array(keys)[order], np.array(vals)[order]

def get_max_len(cds_range):
    max_len = 0
    for tid, (start, end) in cds_range.iteritems():
        if end > max_len:
            max_len = end
    return max_len
    
def is_sparse(vec, threshold):
    order = np.argsort(vec)
    # largest peak not included when determining sparsity
    if np.sum(vec[order][:-1])>= threshold:
        return False
    else:
        return True

def get_ylim(vec, pcnt=3):
    """ 
    pcnt: which peak to use as the ylim
    """
    order = np.argsort(vec)
    idx = order[-pcnt]
    return vec[idx]

def get_pos_hist(pos_hist, start, end):
    y_cnt = np.zeros(end-start+1)
    for pos, cnt in pos_hist.iteritems():
        if pos < start or pos > end: continue
        idx = pos - start
        y_cnt[idx] = cnt
    return y_cnt

def get_cobs(pos_hist, rlstart, rlend, istart, iend):
    pobs = {}
    for rlen in xrange(rlstart, rlend+1):
        pobs[rlen] = get_pos_hist(pos_hist[rlen], istart, iend)
    return pobs

def get_frames(vec):
    fsum = np.array([ np.sum(vec[i::3]) for i in xrange(3) ])
    tsum = float(sum(fsum))
    f0 = np.argmax(fsum)
    fdic = { f0:0, (f0+1)%3:1, (f0+2)%3:2 }
    print "frames: {0:.2%} {1:.2%} {2:.2%}".format(fsum[0]/tsum, fsum[1]/tsum, fsum[2]/tsum)
    return fdic

def get_frame_str_vec(vec):
    fsum = np.array([ np.sum(vec[i::3]) for i in xrange(3) ])
    tsum = float(sum(fsum))
    f0 = np.argmax(fsum)
    fdic = { f0:0, (f0+1)%3:1, (f0+2)%3:2 }
    sframe =  "{0:.2%} {1:.2%} {2:.2%}".format(fsum[0]/tsum, fsum[1]/tsum, fsum[2]/tsum)
    return sframe  

def get_frame_str(pobs):
    return { rlen: get_frame_str_vec(vec) for rlen, vec in pobs.iteritems() }

def get_frame_str_from_tlist(tlist, cds_range):
    sframe = {}
    for rid in tlist:
        for rlen in tlist[rid]['prof']:
            tid = tlist[rid]['tid']
            start, end = cds_range[tid]
            for pos, cnt in tlist[rid]['prof'][rlen]:
                sframe.setdefault(rlen, np.zeros(3))
                sframe[rlen][(pos-start)%3] += cnt
    return { rlen: "{0:.2%} {1:.2%} {2:.2%}".format(f[0]/sum(f), f[1]/sum(f), f[2]/sum(f)) \
             for rlen,f in sframe.iteritems() }

#=============================
# transcript pre-selection
#=============================
def filter_transcript_by_length(cds_range, length_cutoff):
    tid_list = np.array(cds_range.keys())
    tlen = np.array([ v[1]-v[0] for v in cds_range.values()])
    return set(tid_list[tlen > length_cutoff])

#=============================
# generate meta profile
#=============================
def create_rlen_meta_profile(tlist, cds_range, tid_included, ibegin, iend):
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
                meta_hist[rlen].setdefault(i,0)
                meta_hist[rlen][i] += cnt
        sys.stdout.write("processed transcript {0}.\t\r".format(rid))
        sys.stdout.flush()
    sys.stdout.write("\n")
    return meta_hist

#=============================
# read length distribution
#=============================
def get_length_count(tlist):
    rlen2cnt = {}
    for rid in tlist:
        for rlen in tlist[rid]['prof']:
            tid = tlist[rid]['tid']
            for pos, cnt in tlist[rid]['prof'][rlen]:
                rlen2cnt.setdefault(rlen, 0)
                rlen2cnt[rlen] += cnt
    return rlen2cnt

def get_length_distribution(rlen2cnt, rlen_min, rlen_max):
    rlen2cnt_filter = { rlen:cnt for rlen,cnt in rlen2cnt.iteritems() \
                        if rlen >= rlen_min and rlen <= rlen_max }
    tot_cnt = float(sum(rlen2cnt_filter.values()))
    return { rlen: cnt/tot_cnt for rlen,cnt in rlen2cnt_filter.iteritems() }

#=============================
# plot meta profile
#=============================
def plot_pos_hist(meta_hist, sframe, ibegin, iend, min_cnt, fn_prefix):
    print "plotting pos hist..."
    rlen_list = sorted(meta_hist.keys())
    x_pos = np.arange(ibegin, iend+1)
    figwidth = len(x_pos)/10.0*1.5
    for rlen in rlen_list:
        y_cnt = get_pos_hist(meta_hist[rlen], ibegin, iend)
        if is_sparse(y_cnt, min_cnt): continue
        fig = plt.figure(figsize=(figwidth,6))
        ax = fig.add_subplot(1,1,1)
        max_cnt = get_ylim(y_cnt)
        ymax = max_cnt + 1
        plt.bar(-rlen+16-0.4, ymax, width=0.8, color='k', edgecolor='white', alpha=0.4)
        plt.bar(x_pos-0.4, y_cnt, width=0.8, color='k', edgecolor='white', alpha=0.4)
        plt.xlabel("offset from start")
        plt.ylabel("footprint count")
        plt.xticks(range(-100, iend+1, 5))
        plt.xlim((ibegin-1, iend+1))
        plt.ylim((0.1, ymax))
        plt.title("{0}\n{1}".format(rlen, sframe[rlen]), size=20)
        plt.savefig("{0}_{1}_pos_hist.pdf".format(fn_prefix, rlen), bbox_inches="tight")
        plt.close()

def plot_meta_pos_hist(meta_hist, ibegin, iend, fn_prefix):
    print "plotting meta pos hist..."
    x_pos = np.arange(ibegin, iend+1)
    y_cnt = np.array([0]*len(x_pos))
    figwidth = len(x_pos)/10.0*1.5
    for rlen in meta_hist:
        y_cnt_rlen = get_pos_hist(meta_hist[rlen], ibegin, iend)
        y_cnt += y_cnt_rlen
    fig = plt.figure(figsize=(figwidth,6))        
    ax = fig.add_subplot(1,1,1)
    max_cnt = get_ylim(y_cnt)
    ymax = max_cnt + 1
    plt.bar(x_pos-0.4, y_cnt, width=0.8, color='b', edgecolor='white', alpha=0.5)
    imax = np.argmax(y_cnt)
    plt.xlabel("offset from start")
    plt.ylabel("footprint count")
    plt.xticks(range(-100, iend+1, 5))
    plt.xlim((ibegin-1, iend+1))
    plt.ylim((0.1, ymax))
    plt.savefig("{0}_{1}_pos_hist.pdf".format(fn_prefix, "meta"), bbox_inches="tight")
    plt.close()

#=============================
# main
#=============================
def plot_rlen_hist_pipe():
    from deblur_result_io import ensure_dir, get_file_core
    if len(sys.argv) != 4:
        print "Usage: python meta_profile.py input_rlen.hist cds_range output_dir"
        exit(1)
    hist_fn = sys.argv[1]
    cds_txt = sys.argv[2]
    odir = sys.argv[3]
    # imax = 50
    min_sample_cnt = 100
    ensure_dir(odir)
    cds_range = get_cds_range(cds_txt)
    tid_select = filter_transcript_by_length(cds_range, imax)
    tlist = parse_rlen_hist(hist_fn)
    rlen2cnt = get_length_count(tlist)
    rlen2portion = get_length_distribution(rlen2cnt, rlen_min, rlen_max)
    tot_portion = 0
    for rlen in xrange(25, 36):
        print "{0}-mers: {1:.2%}".format(rlen, rlen2portion[rlen])
        tot_portion += rlen2portion[rlen]
    print "total: {0:.2%}".format(tot_portion)
    exit(1)
    sframe = get_frame_str_from_tlist(tlist, cds_range)
    meta_hist = create_rlen_meta_profile(tlist, cds_range, tid_select, utr5_offset, imax)
    fn_prefix = odir+"/"+get_file_core(hist_fn)+"_start_{0}".format(imax)
    plot_pos_hist(meta_hist, sframe, utr5_offset, imax, min_sample_cnt, fn_prefix)    
    plot_meta_pos_hist(meta_hist, utr5_offset, imax, fn_prefix)

if __name__ == "__main__": plot_rlen_hist_pipe()
