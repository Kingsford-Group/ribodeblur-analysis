#!/usr/bin/env python
import sys
import math
from global_params import *
from footprint_hist_parser import *

def get_high_cover_rlen_profile(plist, tlen, cover_ratio, cnt_threshold):
    transcript = {}
    # in-frame codon position coverage > cover_ratio
    pos_covered = math.ceil(tlen/3*cover_ratio)
    for rlen, profile in plist.iteritems():
        pos_cnt = [0]*3
        for loc, cnt in profile:
            if cnt > cnt_threshold:
                pos_cnt[loc%3] += 1
        if max(pos_cnt) > pos_covered:
            transcript[rlen] = profile
    return transcript

def filter_transcript_profiles(tprofile, cds_range, cnt_threshold, cover_ratio):
    """
    filter high covered interesting profiles with different rlen for visualization
    """
    print "filter transcript profiles"
    pcelebrity = {}
    i = 0
    for tid, plist in tprofile.iteritems():
        start, end = cds_range[tid]
        tlen = end-start
        transcript = get_high_cover_rlen_profile(plist, tlen, cover_ratio, cnt_threshold)
        if len(transcript)>1 and 28 in transcript:
            pcelebrity[tid] = transcript.copy()
        if 28 not in transcript and len(transcript)>0:
            print tid, "no 28-mers", transcript.keys()
        i += 1
        sys.stdout.write("processed transcript {0}.\t\r".format(i))
        sys.stdout.flush()
    sys.stdout.write("\n")
    print "total celebrity genes: {0}".format(len(pcelebrity))
    return pcelebrity

def main():
    if len(sys.argv) != 6:
        print "Usage: python high_cover_profile.py input_rlen.hist cds_range.txt cover_ratio cnt_threshold ofname.hist"
        exit(1)
    hist_fn = sys.argv[1]
    cds_txt = sys.argv[2]
    cover_ratio = float(sys.argv[3])
    cnt_threshold = float(sys.argv[4])
    ofname = sys.argv[5]
    cds_range = get_cds_range(cds_txt)
    tlist = parse_rlen_hist(hist_fn)
    tprofile = get_transcript_profiles(tlist, cds_range, utr5_offset, utr3_offset)
    pcelebrity = filter_transcript_profiles(tprofile, cds_range, cnt_threshold, cover_ratio)
    tid2rid = { t['tid']: rid for rid, t in tlist.iteritems() }
    write_rlen_hist(pcelebrity, cds_range, tid2rid, ofname)

if __name__ == "__main__": main()
