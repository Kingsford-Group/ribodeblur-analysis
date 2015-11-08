#!/usr/bin/env python
import sys
import numpy as np
from footprint_hist_parser import parse_rlen_hist, get_cds_range
from deblur_result_io import read_cds_profile

def tot_frame_from_tlist(tlist, tid_list, cds_range):
    frame = np.zeros(3)
    for rid in tlist:
        for rlen in tlist[rid]['prof']:
            tid = tlist[rid]['tid']
            if tid not in tid_list: continue
            start, end = cds_range[tid]
            for pos, cnt in tlist[rid]['prof'][rlen]:
                frame[(pos-start)%3] += cnt
    tot_reads = np.sum(frame)
    print "total reads: {0:.0f}".format(tot_reads)
    frame /= tot_reads
    return frame

def tot_frame_from_prof(prof):
    frame = np.zeros(3)
    for tid in prof:
        for i in range(3):
            frame[i] += np.sum(prof[tid][i::3])
    tot_reads = np.sum(frame)
    print "total reads: {0:.0f}".format(tot_reads)
    frame /= tot_reads
    return frame
    
def main():
    if len(sys.argv) != 4:
        print "Usage: python frame_improvement.py cds_range rlen.hist deblur.out"
        exit(1)
    
    cds_fn = sys.argv[1]
    rlen_hist = sys.argv[2]
    deblur_fn = sys.argv[3]
    
    cds_range = get_cds_range(cds_fn)
    print "frame portion"
    print "after deblur:"
    prof = read_cds_profile(deblur_fn)
    f = tot_frame_from_prof(prof)
    print "{0:.2%} {1:.2%} {2:.2%}\n".format(f[0], f[1], f[2])    

    print "before deblur:"
    tid_list = prof.keys()
    tlist = parse_rlen_hist(rlen_hist)
    f = tot_frame_from_tlist(tlist, tid_list, cds_range)
    print "{0:.2%} {1:.2%} {2:.2%}\n".format(f[0], f[1], f[2])
    
if __name__ == "__main__": main()
