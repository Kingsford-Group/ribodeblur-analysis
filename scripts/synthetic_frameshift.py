#!/usr/bin/env python
import sys
import os
import random
from footprint_hist_parser import *
from global_params import *

def frameshift(plist, pos):
    plist_shift = {}
    for rlen, pos_list in plist.iteritems():
        plist_shift[rlen] = [ (i, cnt) if i<pos else (i+1, cnt) for i,cnt in pos_list ]
    return plist_shift

def batch_frameshifts(tprofile, cds_range):
    random.seed(4238230948)
    tshifts = {}
    shift_points = {}
    for tid, plist in tprofile.iteritems():
        start, stop = cds_range[tid]
        pos = random.randrange(start,stop,3) - start
        tshifts[tid] = frameshift(plist, pos)
        shift_points[tid] = pos
    return tshifts, shift_points

def write_shift_points(shift_points, ofname):
    tf = open(ofname, 'wb')
    text = [ "{0}\t{1}\n".format(tid, pos) for tid, pos in shift_points.iteritems() ]
    tf.writelines(text)
    tf.close()

def read_shift_points(fname):
    shift_points = {}
    tf = open(fname)
    for line in tf:
        tid, pos = line.rstrip().split()
        shift_points[tid] = int(pos)
    tf.close()
    return shift_points

def main():
    if len(sys.argv) != 4:
        print "Usage: python synthetic_frameshift.py high_coverage.hist cds_range.txt ofname.hist"
        exit(1)
    hist_fn = sys.argv[1]
    cds_txt = sys.argv[2]
    ofname = sys.argv[3]
    cds_range = get_cds_range(cds_txt)
    tlist = parse_rlen_hist(hist_fn)
    tprofile = get_transcript_profiles(tlist, cds_range, utr5_offset, utr3_offset)
    pframeshifts, fspts = batch_frameshifts(tprofile, cds_range)
    tid2rid = { t['tid']: rid for rid, t in tlist.iteritems() }
    write_rlen_hist(pframeshifts, cds_range, tid2rid, ofname)
    ofname = os.path.splitext(ofname)[0]+".fsp"
    write_shift_points(fspts, ofname)

if __name__ == "__main__": main()
