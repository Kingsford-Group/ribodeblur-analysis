#!/usr/bin/env python
import numpy as np
import scipy.stats
from multiprocessing import Pool
from footprint_hist_parser import parse_rlen_hist, get_transcript_profiles
from deblur_utils import *
from deblur_result_io import *
from global_params import *

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.size'] = 12
rcParams['xtick.major.size'] = 5
rcParams['ytick.major.size'] = 5
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42

def construct_all_cobs(tprofile, cds_range, utr5_offset, utr3_offset, rlen_min, rlen_max):
    cobs = {}
    i = 0
    for tid, prof in tprofile.iteritems():
        start, end = cds_range[tid]
        cobs[tid] = build_cobs_for_deblur(prof, utr5_offset, (end-start)+utr3_offset, rlen_min, rlen_max)
        i += 1
        sys.stdout.write("processed transcript {0}.\t\r".format(i))
        sys.stdout.flush()
    sys.stdout.write("\n")
    return cobs

class single_transcript_obj(object):
    def __init__(self, blur_vec, klist):
        self.b = blur_vec
        self.klist = klist
    def __call__(self, params):
        tid, cobs, ptrue, eps = params
        abd = { rlen: sum(p) for rlen, p in cobs.iteritems() }
        pobs = { rlen: p/float(abd[rlen]) for rlen, p in cobs.iteritems() \
                 if rlen in eps }
        select = { rlen: p>0 for rlen, p in cobs.iteritems() }
        ptrue_init = initiate_ptrue(cobs)
        before_tot, before_obs, before_true = compute_tot_obj(ptrue_init, abd, pobs, self.klist, select, self.b, None, False)
        after_tot, after_obs, after_true = compute_tot_obj(ptrue, abd, pobs, self.klist, select, self.b, eps, False)
        return tid, before_tot, before_obs, before_true, after_tot, after_obs, after_true

def compute_objs(b, klist, cobs, ptrue, eps):
    print "in compute_objs"
    params_list = [ (tid, cobs[tid], ptrue[tid], eps[tid]) for tid in eps if len(eps[tid])>1 ]
    pool = Pool(processes=nproc)
    # 0-tid 1-before_tot 2-before_obs 3-before_true 4-after_tot 5-after_obs 6-after_true
    result = [ r for r in pool.imap_unordered(single_transcript_obj(b,klist), params_list, 10) ]
    pool.close()
    pool.join()
    print "done compute objs"
    obj_before = { result[i][0] : result[i][2] for i in xrange(len(result)) }
    obj_after = { result[i][0] : result[i][5] for i in xrange(len(result)) }
    return obj_before, obj_after

def compare_pobs_fit(b, klist, cobs, ptrue, eps, ofname):
    abd_list = []
    improve_list = []
    obj_before, obj_after = compute_objs(b, klist, cobs, ptrue, eps)
    i = 0
    for tid in obj_before:
        abd = sum(map(sum, cobs[tid].values()))/float(len(cobs[tid].values()[0]))
        relative_improvement = (obj_before[tid]-obj_after[tid])/float(obj_before[tid])
        abd_list.append(abd)
        improve_list.append(relative_improvement)
        i += 1
        sys.stderr.write("processed transcript {0}.\t\r".format(i))
        sys.stderr.flush()
    sys.stderr.write("\n")
    plt.figure()
    abd_list = np.array(abd_list)
    improve_list = np.array(improve_list)
    abd_filter = abd_list[abd_list>1]
    improve_filter = improve_list[abd_list>1]
    ns, bins, patches = plt.hist(improve_filter, 25, color='b', alpha=0.3)
    plt.xlabel("relative improvement of least square\nbetween model profiles and observed profiles", fontsize=16)
    plt.ylabel("number of transcripts", fontsize=16)
    plt.savefig(ofname, bbox_inches = 'tight')
    plt.close()
    print "improving rate: {0:.2%} total: {1}".format(np.mean(improve_list>0), len(improve_list))
    print "average improving level: {0:.2%} total: {1}".format(np.mean(improve_filter), len(improve_filter))

def get_frame_portion(vec, frame=0):
    return sum(vec[frame::3])/float(sum(vec)) if np.any(vec!=0) else 0

def compare_frame(cobs, ptrue, eps, ofname):
    frame_before = []
    frame_after = []
    i = 0
    for tid in ptrue:
        if len(eps[tid])<2: continue
        cobs_merge = merge_profiles({rlen:prof for rlen, prof in cobs[tid].iteritems() if rlen in eps[tid] })
        ctrue = estimate_ctrue(ptrue[tid], eps[tid], cobs[tid])
        ctrue_merge = merge_profiles(ctrue)
        frame_before.append(get_frame_portion(cobs_merge))
        frame_after.append(get_frame_portion(ctrue_merge))
        i += 1
        sys.stderr.write("processed transcript {0}.\t\r".format(i))
        sys.stderr.flush()
    sys.stderr.write("\n")
    plt.figure()
    ns,bins,patches = plt.hist(frame_before, 25, histtype = 'stepfilled', color='c', edgecolor='c', alpha = 0.4)
    ns,bins,patches = plt.hist(frame_after, 25, histtype = 'stepfilled', color='b', edgecolor='b', hatch='/', alpha = 0.3)
    plt.xlabel("in-frame frequency", fontsize=16)
    plt.ylabel("number of transcripts", fontsize=16)
    plt.legend(["before deblur", "after deblur"], loc="upper left", frameon=False, fontsize=16)
    plt.savefig(ofname, bbox_inches = 'tight')
    plt.close()
    print "improving rate: {0:.2%} total: {1}".format(np.mean(np.array(frame_before)<=np.array(frame_after)), len(frame_before))
    print "average frame portion: before: {0:.2f} after: {1:.2f} mwu: {2:.2e}".format(np.mean(frame_before), np.mean(frame_after), scipy.stats.mannwhitneyu(frame_before, frame_after)[1])

def correlate_true_before_after(cobs, ptrue, eps, ofname):
    corr_list = []
    for tid in ptrue:
        if len(eps[tid])<2: continue
        if 28 not in cobs[tid]: continue
        ptrue_init = cobs[tid][28]
        ctrue = estimate_ctrue(ptrue[tid], eps[tid], cobs[tid])
        ctrue_merge = merge_profiles(ctrue)
        corr_list.append(scipy.stats.pearsonr(ptrue_init, ctrue_merge)[0])
    plt.figure()
    ns,bins, patches = plt.hist(corr_list, bins=25, color='b', alpha=0.3)
    plt.xlabel("correlation between 28-mer profile and deblured profile")
    plt.ylabel("number of transcripts")
    plt.savefig(ofname, bbox_inches='tight')
    plt.close()
    print "average correlation: {0}".format(np.mean(corr_list))

def evaluate_correlation(pcenter, plist):
    return { rlen: scipy.stats.pearsonr(prof, pcenter)[0] for rlen,prof in plist.iteritems() }

def include_corr_to_hist(chist, corr_list):
    for rlen, corr in corr_list.iteritems():
        chist.setdefault(rlen, []).append(corr)

def compare_ptrue_correlation(cobs, cobs_shift, ptrue, eps, ofname):
    chist_before = {}
    chist_middle = {}
    chist_after = {}
    chist_final = {}
    i = 0
    for tid in ptrue:
        if len(eps[tid])<2: continue
        ptrue_init = initiate_ptrue(cobs[tid])
        rlen_list = list(set(eps[tid].keys()) & set(cobs[tid].keys()))
        cobs_tmp = { rlen : cobs[tid][rlen] for rlen in rlen_list }
        cobs_mtmp = { rlen : cobs_shift[tid][rlen] for rlen in rlen_list }
        eps_tmp = { rlen : eps[tid][rlen] for rlen in rlen_list }
        corr_before = evaluate_correlation(ptrue_init, cobs_tmp)
        include_corr_to_hist(chist_before, corr_before)
        corr_mid = evaluate_correlation(ptrue_init, cobs_mtmp)
        include_corr_to_hist(chist_middle, corr_mid)
        ctrue = estimate_ctrue(ptrue[tid], eps_tmp, cobs_tmp)
        corr_after = evaluate_correlation(ptrue_init, ctrue)
        include_corr_to_hist(chist_after, corr_after)
        corr_final = evaluate_correlation(ptrue[tid], ctrue)
        include_corr_to_hist(chist_final, corr_final)
        i += 1
        sys.stdout.write("processed transcript {0}.\t\r".format(i))
        sys.stdout.flush()
    sys.stdout.write("\n")
    data_before = []
    data_mid = []
    data_after = []
    data_final = []
    for rlen in chist_before:
        data_before.append(chist_before[rlen])
        data_mid.append(chist_middle[rlen])
        data_after.append(chist_after[rlen])
        data_final.append(chist_final[rlen])
        print "rlen: {0} total: {1} improving correlation: {2:.2%}".format(rlen, len(chist_before[rlen]), np.mean(np.array(chist_before[rlen])<=np.array(chist_final[rlen])))

    fig = plt.figure(figsize=(20,5))
    ax = fig.add_subplot(1,4,1)
    ax.set_title("observed", fontsize=16)
    plt.boxplot(data_before, labels=chist_before.keys())
    ax.set_ylim((-0.2, 1.02))
    ax.set_ylabel("Pearson correlation", fontsize=16)
    ax = fig.add_subplot(1,4,2)
    ax.set_title("shifted", fontsize=16)
    ax.set_xlabel("read length", fontsize=16)
    plt.boxplot(data_mid, labels=chist_before.keys())
    ax.set_ylim((-0.2, 1.02))
    ax = fig.add_subplot(1,4,3)
    ax.set_title("deblurred", fontsize=16)
    plt.boxplot(data_after, labels=chist_before.keys())
    ax.set_ylim((-0.2, 1.02))
    ax = fig.add_subplot(1,4,4)
    ax.set_title("final", fontsize=16)
    plt.boxplot(data_final, labels=chist_before.keys())
    ax.set_ylim((-0.2, 1.02))
    plt.tight_layout()
    plt.savefig(ofname, bbox_inches='tight')

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print "Usage: python profile_evaluation.py input_rlen.hist cds_range.txt rlen.vblur rlen.eps output_dir"
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
    vrlen_min, vrlen_max = get_rlen_range_from_vblur(b)
    print "get pre-computed deblur results"
    ptrue, eps = read_essentials(deblur_fname)
    print "construct cobs all at once"
    tprofile = get_transcript_profiles(tlist, cds_range, utr5_offset, utr3_offset)
    cobs = construct_all_cobs(tprofile, cds_range, utr5_offset, utr3_offset, vrlen_min, vrlen_max)
    cobs_shift = build_cobs_with_shifts(tprofile, cds_range, utr5_offset, utr3_offset, vrlen_min, vrlen_max, klist)
    print "compare least square of pobs fitting"
    compare_pobs_fit(b, klist, cobs, ptrue, eps, odir+fname+"_cmp_pobs.pdf")
    print "compare frame distribution"
    compare_frame(cobs, ptrue, eps, odir+fname+"_cmp_frame.png")
    print "compare correlation of ptrue"
    compare_ptrue_correlation(cobs, cobs_shift, ptrue, eps, odir+fname+"_cmp_corr.pdf")
    print "compare ptrue with 28-mers"
    correlate_true_before_after(cobs, ptrue, eps, odir+fname+"_cmp_true.pdf")
