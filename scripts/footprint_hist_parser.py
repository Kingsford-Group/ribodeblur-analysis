#!/usr/bin/env python
import sys
from Bio import SeqIO

def parse_rlen_hist(fname):
    print "parsing footprint counts grouped by read length"
    tlist = {}
    transcript = {}
    tf = open(fname)
    i = 0
    line = tf.readline()
    while line:
        if line.startswith("refID: "):
            if transcript:
                rid = transcript["id"]
                i += 1
                sys.stdout.write("processed transcript {0}.\t\r".format(i))
                sys.stdout.flush()
                tlist[rid] = transcript.copy()
                transcript.clear()
            rid = int(line.lstrip("refID: ").rstrip("\n"))
            transcript["id"] = rid
        elif line.startswith("tid: "):
            tid = line.lstrip("tid: ").rstrip("\n")
            transcript["tid"] = tid
            transcript["prof"] = {}
        # profile len grouped by read length
        # format: read len: pos count, ...
        # [pos count] is the adjacency list way to save space
        else:
            read_len, pc_str = line.rstrip("\n").split(": ")
            read_len = int(read_len)
            prof = []
            pc_pairs = pc_str.rstrip(", ").split(", ")
            for pc in pc_pairs:
                pos, cnt = pc.split(" ")
                prof.append([int(pos), float(cnt)])
            transcript["prof"][read_len] = prof
        line = tf.readline()
    if transcript:
        rid = transcript["id"]
        tlist[rid] = transcript
    tf.close()
    sys.stdout.write("\n")
    return tlist

def get_transcript_profiles(tlist, cds_range, ibegin, iend):
    """ 
    create read length specific count profiles for the read start loci for each transcript
    ibegin: index bound to include before START (negative)
    iend: bases to include after STOP (positive)
    """
    print "create transcript profiles"
    tprofile = {}
    for rid in tlist:
        tid = tlist[rid]['tid']
        start, end = cds_range[tid]
        tlen = end-start
        for rlen in tlist[rid]['prof']:
            for pos, cnt in tlist[rid]['prof'][rlen]:
                i = pos - start
                if i < ibegin or i > tlen+iend : continue
                tprofile.setdefault(tid,{})
                tprofile[tid].setdefault(rlen, [])
                tprofile[tid][rlen].append((i,cnt))
        sys.stdout.write("processed transcript {0}.\t\r".format(rid))
        sys.stdout.flush()
    sys.stdout.write("\n")
    return tprofile

def write_rlen_hist(tprofile, cds_range, tid2rid, ofname):
    tf = open(ofname, 'wb')
    i = 0
    for tid, plist in tprofile.iteritems():
        start, stop = cds_range[tid]
        text = [ "refID: {0}\n".format(tid2rid[tid]), 
                 "tid: {0}\n".format(tid) ] + \
            [ "{0}: {1}\n".format(rlen, ", ".join(map(lambda x: "{0} {1:.0f}".format(x[0]+start, x[1]), pos_list))) for rlen, pos_list in plist.iteritems() ]
        tf.writelines(text)
        i += 1
        sys.stderr.write("processed transcript {0}.\t\t\r".format(i))
        sys.stderr.flush()
    sys.stderr.write("\n")
    tf.close()

#=============================
# stats on read length
# and frame distribution
#=============================
def print_stats(tlist, cds_range):
    print "get read frame stats"
    fcnt = [0]*3
    rlmin = 100
    rlmax = 0
    for rid in tlist:
        for rlen in tlist[rid]['prof']:
            if rlen < rlmin:
                rlmin = rlen
            if rlen > rlmax:
                rlmax = rlen
            tid = tlist[rid]['tid']
            start, stop = cds_range[tid]
            for pos, cnt in tlist[rid]['prof'][rlen]:
                fcnt[(pos-start)%3] += cnt
    tcnt = float(sum(fcnt))
    print "total footprints: {0:.0f}".format(tcnt)
    print "frame distribution: {0:.2%} {1:.2%} {2:.2%}".format(fcnt[0]/tcnt, fcnt[1]/tcnt, fcnt[2]/tcnt)
    print "read length range: {0}-{1}".format(rlmin, rlmax)
    return rlmin, rlmax

def get_pep_seq(fn):
    print "getting peptide sequence..."
    ifile = open(fn, 'r')
    pep_seq = {}
    for rec in SeqIO.parse(ifile, "fasta"):
        tid = rec.id.split()[0]
        pseq = str(rec.seq)
        pep_seq[tid] = pseq
    ifile.close()
    return pep_seq

def get_cds_range(fn):
    print "getting cds range.."
    ifile = open(fn, 'r')
    cds_range = {}
    for line in ifile:
        tid, start, stop = line.strip().split()
        cds_range[tid] = (int(start), int(stop))
    ifile.close()
    return cds_range

def get_tseq(fn, cds_range):
    print "getting transcript sequences..."
    ifile = open(fn, 'r')
    tseq = {}
    for rec in SeqIO.parse(ifile, "fasta"):
        tid = rec.id.split()[0]
        if cds_range == None:
            tseq[tid] = str(rec.seq)
        else:
            start, stop = cds_range[tid]
            seq = str(rec.seq[start:stop])
            tseq[tid] = seq
    ifile.close()
    return tseq

def get_rates(fn, line_sep="\n", sep="\t"):
    rate = {}
    ifile = open(fn, 'r')
    for line in ifile:
        c, r = line.strip(line_sep).split(sep)
        rate[c] = float(r)
    ifile.close()
    return rate

def write_rates(fn, rate):
    ofile = open(fn, 'w')
    text = ["{0}\t{1}\n".format(c,rate[c]) for c in sorted(rate.keys()) ]
    ofile.writelines(text)
    ofile.close()

if __name__ == "__main__":
    hist_fn = "../data/BY_FP_rlen_5p.hist"
    tlist = parse_rlen_hist(hist_fn)
