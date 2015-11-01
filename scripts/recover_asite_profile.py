#!/usr/bin/env python
import numpy as np
import sys
from multiprocessing import Pool
from meta_profile import *
from deblur_utils import *
from deblur_result_io import *
from global_params import *

def shift_rlen_profile(vec, offset):
    vlen = len(vec)
    new_vec = np.zeros(vlen)
    for i in xrange(vlen):
        idx_adj = i-offset
        if idx_adj < 0 or idx_adj >= vlen: continue
        new_vec[idx_adj] = vec[i]
    return new_vec

def shift_transcript_profile(cobs, klist):
    return { rlen: shift_rlen_profile(prof, klist[rlen]) for rlen, prof in cobs.iteritems() }

def get_dominant_frame(vec):
    return np.argmax([ np.sum(vec[i::3]) for i in range(3) ])

def separate_high_cover_rlen_profile(cobs, cover_ratio=0.5, cnt_threshold=0):
    tlen = len(cobs.values()[0])
    pos_covered = np.ceil(cover_ratio*tlen/3)
    cobs_hc = {}
    cobs_lc = {}
    for rlen, prof in cobs.iteritems():
        # all positions with high enough count
        hcnt_vec = prof > cnt_threshold
        # dominant frame with the highest read count sum
        f0 = get_dominant_frame(prof)
        # in-frame codon position coverage > cover_ratio
        if np.sum(hcnt_vec[f0::3]) > pos_covered:
            cobs_hc[rlen] = prof
        else:
            cobs_lc[rlen] = prof
    return cobs_hc, cobs_lc

class single_transcript_asite_shift(object):
    def __init__(self, blur_vec, klist, converge_cutoff, cover_ratio, cnt_threshold):
        self.b = blur_vec
        self.k = klist
        self.c = converge_cutoff
        self.cover_ratio = cover_ratio
        self.cnt_threshold = cnt_threshold

    def __call__(self, params):
        tid = params[0]
        cobs = params[1]
        cobs_hc, cobs_lc = separate_high_cover_rlen_profile(cobs, self.cover_ratio, self.cnt_threshold)
        ctrue = {}
        if len(cobs_lc) != 0:
            ctrue_lc = shift_transcript_profile(cobs_lc, self.k)
            ctrue.update(ctrue_lc)
        if len(cobs_hc)!=0:
            ptrue, eps = recover_true_profile(cobs_hc, self.k, self.b, 0, 100, self.c)
            if ptrue is None: 
                ctrue_hc = shift_transcript_profile(cobs_hc, self.k)
            else:
                ctrue_hc = estimate_ctrue(ptrue, eps, cobs_hc)
            ctrue.update(ctrue_hc)
        else:
            ptrue = None
            eps = None
        ctrue_merge = merge_profiles(ctrue)
        return tid, ctrue_merge, ptrue, eps

class single_transcript_asite_deblur(object):
    def __init__(self, blur_vec, klist, converge_cutoff, cover_ratio, cnt_threshold):
        self.b = blur_vec
        self.k = klist
        self.c = converge_cutoff
        self.cover_ratio = cover_ratio
        self.cnt_threshold = cnt_threshold

    def __call__(self, params):
        tid = params[0]
        cobs = params[1]
        cobs_hc, cobs_lc = separate_high_cover_rlen_profile(cobs, self.cover_ratio, self.cnt_threshold)
        ctrue = {}
        if len(cobs_hc)!=0:
            ptrue, eps = recover_true_profile(cobs_hc, self.k, self.b, 0, 100, self.c)
            # deblur failed, merge high-coverage profiles to low-coverage profiles
            if ptrue is None: 
                cobs_lc.update(cobs_hc)
            else:
                ctrue_hc = estimate_ctrue(ptrue, eps, cobs_hc)
                ctrue.update(ctrue_hc)
                # deblur partially failed, merge failed profiles to low-coverage profiles
                if len(eps) != len(cobs_hc):
                    cobs_failed = {rlen:prof for rlen,prof in cobs_hc.iteritems() if rlen not in eps}
                    cobs_lc.update(cobs_failed)
        else:
            ptrue = None
            eps = None
        if len(cobs_lc) != 0:
            # deblur the rest all together
            # use rlen=0 to store this one
            ctrue_lc = recover_sparse_true_profile(cobs_lc, self.k, self.b)
            ctrue[0] = ctrue_lc
        ctrue_merge = merge_profiles(ctrue)
        return tid, ctrue_merge, ptrue, eps

def batch_Asite_recovery_parallel(tprofile, cds_range, utr5_offset, utr3_offset, rlen_min, rlen_max, blur_vec, klist, converge_cutoff, cover_ratio, cnt_threshold, nproc):
    cobs_all = np.array([ build_cobs_for_deblur(prof, utr5_offset, (cds_range[tid][1]-cds_range[tid][0])+utr3_offset, rlen_min, rlen_max) for tid,prof in tprofile.iteritems() ])
    cobs_len_all = np.array(map(len, cobs_all))
    tid_all = np.array(tprofile.keys())
    cobs_in = cobs_all[cobs_len_all!=0]
    tid_in = tid_all[cobs_len_all!=0]
    print "total transcripts with legit cobs {0}".format(len(tid_in))
    print "batch A-site recovery"
    pool = Pool(processes=nproc)
    results = [ r for r in pool.imap_unordered(single_transcript_asite_deblur(blur_vec, klist, converge_cutoff, cover_ratio, cnt_threshold), zip(tid_in, cobs_in), 100) ]
    pool.close()
    pool.join()
    tid_list, ctrue_list, ptrue_list, eps_list = zip(*results)
    tid_list = np.array(tid_list)
    ptrue_list = np.array(ptrue_list)
    eps_list = np.array(eps_list)
    valid = np.array(map(lambda x: x != None, ptrue_list))
    ctrue = dict(zip(tid_list, ctrue_list))
    ptrue = dict(zip(tid_list[valid], ptrue_list[valid]))
    eps = dict(zip(tid_list[valid], eps_list[valid]))
    print "total deblurred transcripts: {0}".format(len(ptrue))
    print "total processed transcripts: {0}".format(len(ctrue))
    return ctrue, ptrue, eps

def batch_Asite_recovery(tprofile, cds_range, utr5_offset, utr3_offset, rlen_min, rlen_max, blur_vec, klist, converge_cutoff, cover_ratio, cnt_threshold, nproc):
    ctrue = {}
    ptrue = {}
    eps = {}
    i = 0
    tlen = np.array([ cds_range[tid][1]-cds_range[tid][0] for tid in cds_range ])
    tid_list = np.array(cds_range.keys())
    order = np.argsort(tlen)
    tid_list = tid_list[order]
    for tid in tid_list:
        if tid not in tprofile: continue
        i += 1
        sys.stdout.write("deblur {0}th transcript: {1}\t\t\r".format(i, tid))
        sys.stdout.flush()        
        start, end = cds_range[tid]
        cobs = build_cobs_for_deblur(tprofile[tid], utr5_offset, (end-start)+utr3_offset, rlen_min, rlen_max)
        if len(cobs) == 0: continue
        deblurer = single_transcript_asite_deblur(blur_vec, klist, converge_cutoff, cover_ratio, cnt_threshold)
        tid, ctrue_tid, ptrue_tid, eps_tid = deblurer((tid,cobs))
        ctrue[tid] = ctrue_tid
        ptrue[tid] = ptrue_tid
        eps[tid] = eps_tid
        cobs_merge = merge_profiles(cobs)
        x_pos = np.arange(len(ctrue_tid))
        figwidth = len(x_pos)/10.0*1.5
        fig = plt.figure(figsize=(figwidth, 13))
        ax = fig.add_subplot(2,1,1)
        plt.bar(x_pos-0.4, cobs_merge, width=0.8, color='b', edgecolor='white', alpha=0.3)
        yheight = max(cobs_merge)+1
        plt.ylim((0,yheight))
        xgrid = x_pos[::3]
        plt.bar(xgrid-0.4, [yheight]*len(xgrid), width=0.8, color='r', edgecolor='white', alpha=0.1)
        plt.xlim((x_pos[0]-1, x_pos[-1]+1))
        plt.title("Before deblur".format(tid))
        ax = fig.add_subplot(2,1,2)
        plt.bar(x_pos-0.4, ctrue_tid, width=0.8, color='b', edgecolor='white', alpha=0.3)
        yheight = max(ctrue_tid)+1
        plt.bar(xgrid-0.4, [yheight]*len(xgrid), width=0.8, color='r', edgecolor='white', alpha=0.1)
        plt.ylim((0,yheight))
        plt.xlim((x_pos[0]-1, x_pos[-1]+1))
        plt.title("After deblur")
        plt.tight_layout()
        plt.savefig("{0}.pdf".format(tid), bbox_inches='tight')
        if i>10: break
    print "\ntotal deblurred transcripts: {0}".format(len(ptrue))
    return ctrue, ptrue, eps

def batch_build_Aprof(prof_dic, cds_range, utr5_offset, asite_offset):
    aprof = {}
    for tid, prof in prof_dic.iteritems():
        cds_start, cds_end = cds_range[tid]
        istart = utr5_offset - asite_offset
        iend = istart + ( cds_end - cds_start )
        if np.any(prof[istart: iend]!=0):
            aprof[tid] = prof[istart: iend]
    return aprof

def main():
    if len(sys.argv) != 5:
        print "Usage: python recover_asite_profile.py input_rlen.hist input_rlen.vblur cds_range.txt output_dir"
        exit(1)
    hist_fn = sys.argv[1]
    vblur_txt = sys.argv[2]
    cds_txt = sys.argv[3]
    odir = sys.argv[4]
    ensure_dir(odir)
    print "get pre-computed blur vector"
    b = read_vblur(vblur_txt)
    vrlen_min, vrlen_max = get_rlen_range_from_vblur(b)
    cds_range = get_cds_range(cds_txt)
    tlist = parse_rlen_hist(hist_fn)
    # build profile for each transcript per read length
    tprofile = get_transcript_profiles(tlist, cds_range, utr5_offset, utr3_offset)
    ctrue, ptrue, eps = batch_Asite_recovery_parallel(tprofile, cds_range, utr5_offset, utr3_offset, vrlen_min, vrlen_max, b, klist, converge_cutoff, cover_ratio, cnt_threshold, nproc)
    base_prof = batch_build_Aprof(ctrue, cds_range, -utr5_offset, asite_offset)
    ofname = odir + get_file_core(hist_fn) + "_cds.base"
    write_cds_profile(base_prof, ofname)
    ofname = odir + get_file_core(hist_fn) + ".eps"
    write_essentials(ptrue, eps, ofname)

if __name__ == "__main__": main()


        

