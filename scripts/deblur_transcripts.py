#!/usr/bin/env python
import numpy as np
import sys
from multiprocessing import Pool
from meta_profile import *
from deblur_utils import *
from deblur_result_io import *
from global_params import *

def build_profile_from_list(pos_list, start, end):
    profile = np.zeros(end-start+1)
    for pos, cnt in pos_list:
        profile[pos-start] = cnt
    return profile

def build_cobs_for_deblur(clist, start, end, rlen_min, rlen_max):
    cobs = {}
    for rlen, pos_list in clist.iteritems():
        if rlen < rlen_min or rlen > rlen_max: continue
        cobs[rlen] = build_profile_from_list(pos_list, start, end)
    return cobs
    
def batch_Asite_recovery(tprofile, cds_range, utr5_offset, utr3_offset, rlen_min, rlen_max, blur_vec, converge_cutoff):
    ptrue = {}
    eps = {}
    i = 0
    for tid in tprofile:
        i += 1
        sys.stdout.write("deblur {0}th transcript: {1}\t\t\r".format(i, tid))
        sys.stdout.flush()        
        start, end = cds_range[tid]
        cobs = build_cobs_for_deblur(tprofile[tid], utr5_offset, (end-start)+utr3_offset, rlen_min, rlen_max)
        if len(cobs) == 0: continue
        ptrue_tid, eps_tid = recover_true_profile(cobs, blur_vec, 0, 100, converge_cutoff, odir+"{0}_obj.pdf".format(tid))
        if np.all(ptrue_tid==0): continue
        ptrue[tid] = ptrue_tid
        eps[tid] = eps_tid
    print "\ntotal deblurred transcripts: {0}".format(len(ptrue))
    return ptrue, eps

class single_transcript_asite(object):
    def __init__(self, blur_vec, converge_cutoff):
        self.b = blur_vec
        self.c = converge_cutoff
    def __call__(self, params):
        tid = params[0]
        cobs = params[1]
        ptrue, eps = recover_true_profile(cobs, self.b, 0, 100, self.c)
        return tid, ptrue, eps

def batch_Asite_recovery_parallel(tprofile, cds_range, utr5_offset, utr3_offset, rlen_min, rlen_max, blur_vec, converge_cutoff, nproc):
    cobs_all = np.array([ build_cobs_for_deblur(prof, utr5_offset, (cds_range[tid][1]-cds_range[tid][0])+utr3_offset, rlen_min, rlen_max) for tid,prof in tprofile.iteritems() ])
    cobs_len_all = np.array(map(len, cobs_all))
    tid_all = np.array(tprofile.keys())
    cobs_in = cobs_all[cobs_len_all!=0]
    tid_in = tid_all[cobs_len_all!=0]
    pool = Pool(processes=nproc)
    results = [ r for r in pool.imap_unordered(single_transcript_asite(blur_vec, converge_cutoff), zip(tid_in, cobs_in), 10) ]
    pool.close()
    pool.join()
    tid_list, ptrue_list, eps_list = zip(*results)
    tid_list = np.array(tid_list)
    ptrue_list = np.array(ptrue_list)
    eps_list = np.array(eps_list)
    valid = np.array(map(lambda x: x != None, ptrue_list))
    ptrue = dict(zip(tid_list[valid], ptrue_list[valid]))
    eps = dict(zip(tid_list[valid], eps_list[valid]))
    print "\ntotal deblurred transcripts: {0}".format(len(ptrue))
    return ptrue, eps

def main():
    if len(sys.argv) != 5:
        print "Usage: python deblur_transcripts.py input_rlen.hist input_rlen.vblur cds_range.txt output_dir"
        exit(1)
    hist_fn = sys.argv[1]
    vblur_txt = sys.argv[2]
    cds_txt = sys.argv[3]
    odir = sys.argv[4]
    ensure_dir(odir)
    print "get pre-computed blur vector"
    b = read_vblur(vblur_txt)
    rlen_list = sorted(b.keys())
    vrlen_min = rlen_list[0]
    vrlen_max = rlen_list[-1]
    cds_range = get_cds_range(cds_txt)
    tlist = parse_rlen_hist(hist_fn)
    # build profile for each transcript per read length
    tprofile = get_transcript_profiles(tlist, cds_range, utr5_offset, utr3_offset)
    print "batch A-site recovery"
    ptrue, eps = batch_Asite_recovery_parallel(tprofile, cds_range, utr5_offset, utr3_offset, vrlen_min, vrlen_max, b, converge_cutoff, nproc)
    ofname = odir+get_file_core(hist_fn)+".eps"
    write_essentials(ptrue, eps, ofname)

if __name__ == "__main__": main()


        

