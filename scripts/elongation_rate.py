#!/usr/bin/env python
import sys
import numpy as np
import scipy.stats
from global_params import *
from deblur_utils import *
from deblur_result_io import *
from footprint_hist_parser import get_cds_range, parse_rlen_hist, get_transcript_profiles, get_tseq
from profile_evaluation import construct_all_cobs

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.size'] = 12
rcParams['xtick.major.size'] = 3
rcParams['ytick.major.size'] = 3
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42

def batch_build_Aprof(prof_dic, cds_range, utr5_offset, asite_offset):
    aprof = {}
    for tid, prof in prof_dic.iteritems():
        cds_start, cds_end = cds_range[tid]
        istart = utr5_offset - asite_offset
        iend = istart + ( cds_end - cds_start )
        aprof[tid] = prof[istart: iend]
    return aprof

def merge_frames(p):
    """ merge off-frames to the closest frame0"""
    if len(p)<3: return p
    p_sum = np.array(p[::3])
    pplus = np.array(p[1::3])
    p_sum += pplus
    pplus = np.zeros(len(p_sum))
    pplus[1:] = p[2::3][:-1]
    p_sum += pplus
    return p_sum

def sum_frames(p):
    """ sum frame 1 and 2 to frame0 """
    return [ sum(p[i:i+3]) for i in xrange(0,len(p),3) ]

def in_frames(p):
    return p[::3]

def codonp_from_basep(prof_dic, merge_func):
    """convert base counts to codon counts"""
    return { tid: merge_func(prof) for tid,prof in prof_dic.iteritems() }

def get_mean_rc(vec, threshold=1):
    return np.mean(vec[vec>threshold])

def append_cnt_to_codon(dic, codon, cnt):
    dic.setdefault(codon, []).append(cnt)

def include_loc(p, seq, loc, dic):
    codon = seq[loc*3 : (loc+1)*3]
    cnt = p[loc]
    append_cnt_to_codon(dic, codon, cnt)

def group_cnts_by_codon(prof_dic, tseq, rc_threshold=1, start=20, stop=-20):
    codon_dic = {}
    for tid, rprof in prof_dic.iteritems():
        assert len(rprof)*3 == len(tseq[tid])
        # to be consistent with Tuller, skip first and last 20 codons
        if stop == 0: stop = len(rprof)
        rprof = np.array(rprof)
        plen = len(rprof[start:stop])
        if plen==0: continue
        nrc = get_mean_rc(rprof[start:stop], rc_threshold)
        prof = rprof/nrc
        # group cnts by codons
        for c in xrange(start, start+plen):
            if rprof[c] > rc_threshold :
                include_loc(prof, tseq[tid], c, codon_dic)
    return codon_dic

def fit_lognormal(rc_list):
    params = {}
    for codon in rc_list:
        shape, loc, scale = scipy.stats.lognorm.fit(rc_list[codon])
        params[codon] = (shape, loc, scale)
    return params

def lognormal_pdf(x, params):
    """ wrapper function for plotting"""
    shape, loc, scale = params
    return scipy.stats.lognorm.pdf(x,shape,loc, scale)

def generate_codon_list_alphabetical():
    codon_list = []
    base_list = ['T', 'C', 'G', 'A']
    for b0 in base_list:
        for b1 in base_list:
            for b2 in base_list:
                codon = b0+b1+b2
                codon_list.append(codon)
    return codon_list

def exclude_keys(dic, key_list):
    return { k:v for k,v in dic.iteritems() if k not in key_list }

def fit_cnt_to_lognormal(prof_dic, cds_range, utr5_offset, asite_offset, tseq):
    mprof = { tid: merge_profiles(plist) for tid, plist in prof_dic.iteritems() }
    base_prof = batch_build_Aprof(mprof, cds_range, -utr5_offset, asite_offset)    
    print "covert base profiles to codon profiles"
    codon_prof = codonp_from_basep(base_prof, in_frames)
    print "group counts by codons"
    codon_dic = group_cnts_by_codon(codon_prof, tseq, rc_threshold=1, start=20, stop=-20)
    print "fit lognormal model"
    params = fit_lognormal(codon_dic)
    return params

def get_lognormal_stats(params):
    means = {}
    varrs = {}
    skews = {}
    for c, param in params.iteritems():
        shape, loc, scale = param
        mean, var, skew, kurt = scipy.stats.lognorm.stats(shape, loc, scale, moments='mvsk')
        means[c] = mean
        varrs[c] = var
        skews[c] = skew
    return means, varrs, skews

if __name__ == "__main__":
    if len(sys.argv) != 7:
        print "Usage: python elongation_rate.py transcript.fasta cds_range.txt rlen.hist rlen.vblur rlen.eps output_dir"
        exit(1)
    tfasta = sys.argv[1]
    cds_txt = sys.argv[2]
    hist_fn = sys.argv[3]
    vblur_fname = sys.argv[4]
    deblur_fname = sys.argv[5]
    odir = sys.argv[6]
    ensure_dir(odir)
    fname = get_file_core(hist_fn)
    print "get cds range"
    cds_range = get_cds_range(cds_txt)
    tseq = get_tseq(tfasta, cds_range)
    print "parse read len hist file"
    tlist = parse_rlen_hist(hist_fn)
    print "get pre-computed blur vector"
    b = read_vblur(vblur_fname)
    vrlen_min, vrlen_max = get_rlen_range_from_vblur(b)
    print "get pre-computed deblur results"
    ptrue, eps = read_essentials(deblur_fname)
    print "construct cobs all at once"
    tprofile = get_transcript_profiles(tlist, cds_range, utr5_offset, utr3_offset)
    # cobs = build_cobs_with_shifts(tprofile, cds_range, utr5_offset, utr3_offset, rlen_min, rlen_max, klist)
    cobs = construct_all_cobs(tprofile, cds_range, utr5_offset, utr3_offset, vrlen_min, vrlen_max)
    print "construct ctrue all at once"
    ctrue_rlen = batch_build_ctrue(ptrue, eps, cobs)
    mprof = { tid: merge_profiles(plist) for tid, plist in ctrue_rlen.iteritems() }
    base_prof = batch_build_Aprof(mprof, cds_range, -utr5_offset, asite_offset)    
    print "covert base profiles to codon profiles"
    codon_prof = codonp_from_basep(base_prof, in_frames)
    print "group counts by codons"
    codon_deblur = group_cnts_by_codon(codon_prof, tseq, rc_threshold=1, start=20, stop=-20)
    print "fit lognormal model"
    params_deblur = fit_lognormal(codon_deblur)
    mean_deblur, var_deblur, skew_deblur = get_lognormal_stats(params_deblur)

    mprof = { tid: merge_profiles(plist) for tid, plist in cobs.iteritems() }
    base_prof = batch_build_Aprof(mprof, cds_range, -utr5_offset, asite_offset)    
    print "covert base profiles to codon profiles"
    codon_prof = codonp_from_basep(base_prof, in_frames)
    print "group counts by codons"
    codon_blur = group_cnts_by_codon(codon_prof, tseq, rc_threshold=1, start=20, stop=-20)
    print "fit lognormal model"
    params_blur = fit_lognormal(codon_blur)
    mean_blur, var_blur, skew_blur = get_lognormal_stats(params_blur)

    # plot_codon_hist_fitting(codon_dic, params, lognormal_pdf, fn_prefix=odir+"fp")
    tai_raw = get_rates("Tuller_tai.txt")#"codon_usage.txt")
    stop_codons = ["TAA", "TGA", "TAG"]
    tai_raw = exclude_keys(tai_raw, stop_codons)
    clist = np.array(sorted(tai_raw.keys()))

    # deblur_list = []
    # for c in clist:
    #     data_c = codon_deblur[c]
    #     param_c = params_deblur[c]
    #     d, p = scipy.stats.kstest(data_c, scipy.stats.lognorm.cdf, param_c)
    #     print c, p
    #     deblur_list.append(d)

    # blur_list = []
    # for c in clist:
    #     data_c = codon_blur[c]
    #     param_c = params_blur[c]
    #     d, p = scipy.stats.kstest(data_c, scipy.stats.lognorm.cdf, param_c)
    #     print c, p
    #     blur_list.append(d)

    #tai_list = [ 1.0/tai_inverse[c] for c in clist ]
    tai_list = [ tai_raw[c] for c in clist ]
    deblur_list = [ np.var(codon_deblur[c]) for c in clist ]
    blur_list = [ np.var(codon_blur[c]) for c in clist ]
    diff_list = np.array([ deblur_list[i] - blur_list[i] for i in xrange(len(deblur_list)) ])
    print "larger variance in deblur:", clist[diff_list>0], "{0:.2%}".format(sum(diff_list>0)/float(len(diff_list)))

    deblur_list = [ skew_deblur[c] for c in clist ]
    blur_list = [ skew_blur[c] for c in clist ]
    print "deblur"
    print "pearson", scipy.stats.pearsonr(tai_list, deblur_list)
    print "spearman", scipy.stats.spearmanr(tai_list, deblur_list)
    print "blur"
    print "pearson", scipy.stats.pearsonr(tai_list, blur_list)
    print "spearman", scipy.stats.spearmanr(tai_list, blur_list)
    
    # # visualizing correlation
    # plt.figure()
    # plt.plot(tai_list, blur_list, 'ro', markeredgecolor='None', alpha=0.2, markersize=10, label='before deblur')
    # plt.plot(tai_list, deblur_list, 'bo', markeredgecolor='None', alpha=0.2, markersize=10, label='after deblur')
    # linear_fit = np.polyfit(tai_list, deblur_list,1)
    # fit_fn = np.poly1d(linear_fit)
    # err = np.array(deblur_list)-fit_fn(tai_list)
    # high_cut = np.percentile(err,90)
    # low_cut = np.percentile(err,10)
    # for i in xrange(len(clist)):
    #     if err[i] > high_cut or err[i] < low_cut :
    #         codon = clist[i]
    #         x = tai_list[i]
    #         y = deblur_list[i]
    #         plt.annotate(codon, xy=(x,y))
    # plt.xlabel('tAI')
    # plt.ylabel('estimated codon decoding time')
    # plt.legend()
    # plt.show()

    # plt.plot(blur_list, deblur_list, 'bo', alpha=0.2, markeredgecolor='None', markersize=10)
    # linear_fit = np.polyfit(blur_list, deblur_list,1)
    # fit_fn = np.poly1d(linear_fit)
    # err = np.array(deblur_list)-fit_fn(blur_list)
    # high_cut = np.percentile(err,90)
    # low_cut = np.percentile(err,10)
    # for i in xrange(len(clist)):
    #     if err[i] > high_cut or err[i] < low_cut :
    #         codon = clist[i]
    #         x = blur_list[i]
    #         y = deblur_list[i]
    #         plt.annotate(codon, xy=(x,y))
    # xmin=min(min(blur_list),min(deblur_list))
    # xmax=max(max(blur_list),max(deblur_list))
    # plt.plot(range(xmin,xmax),range(xmin,xmax), 'r-', linewidth=3,alpha=0.5)
    # plt.plot(blur_list, fit_fn(blur_list), 'r-', linewidth=3, alpha=0.5)
    # plt.xlim((min(blur_list), max(blur_list)))
    # plt.xlabel('estimated codon decoding time before deblur')
    # plt.ylabel('estimated codon decoding time after deblur')
    # plt.show()
    
    from wobble_pairing import *
    get_diff = lambda codon_dic, wobble_pairs: [ (codon_dic[cwb]-codon_dic[cwc])/float(codon_dic[cwc]) for cwc, cwb in wobble_pairs ]
    get_2base = lambda wobble_pairs: [ cwc[:2] for cwc, cwb in wobble_pairs ]
    ymin = -0.1
    ymax = 0.3
    print "GCU"
    blur_list = np.array(get_diff(skew_blur, GCU_pairs))
    deblur_list = np.array(get_diff(skew_deblur, GCU_pairs))
    names = np.array(get_2base(GCU_pairs))
    order = np.argsort(deblur_list)[::-1]
    blur_list = blur_list[order]
    deblur_list = deblur_list[order]
    names = names[order]
    plt.figure(figsize=(len(names)/2.0, 5))
    x = range(0,len(names)*2, 2)
    plt.bar(x,blur_list, color='b', edgecolor='white', alpha=0.2)
    x = range(1,len(names)*2+1, 2)
    plt.bar(x,deblur_list, color='b', edgecolor='white', alpha=0.6)
    plt.xticks(x, names, rotation='vertical', fontsize=12)
    plt.ylim((ymin,ymax))
    plt.savefig(odir+"skew_GCU_diff.pdf", bbox_inches='tight')
    plt.close()
    print "mean diff: deblur: {0} blur: {1}".format(np.mean(deblur_list), np.mean(blur_list))
    print "ICU"    
    blur_list = np.array(get_diff(skew_blur, ICU_pairs))
    deblur_list = np.array(get_diff(skew_deblur, ICU_pairs))
    names = np.array(get_2base(ICU_pairs))
    order = np.argsort(deblur_list)[::-1]
    blur_list = blur_list[order]
    deblur_list = deblur_list[order]
    names = names[order]
    plt.figure(figsize=(len(names)/2.0, 5))
    x = range(0,len(names)*2, 2)
    plt.bar(x,blur_list, color='r', edgecolor='white', alpha=0.2)
    x = range(1,len(names)*2+1, 2)
    plt.bar(x,deblur_list, color='r', edgecolor='white', alpha=0.6)
    plt.xticks(x, names, rotation='vertical', fontsize=12)
    plt.ylim((ymin,ymax))
    plt.savefig(odir+"skew_ICU_diff.pdf", bbox_inches='tight')
    plt.close()
    print "mean diff: deblur: {0} blur: {1}".format(np.mean(deblur_list), np.mean(blur_list))
    print "UAG"    
    blur_list = np.array(get_diff(skew_blur, UAG_pairs))
    deblur_list = np.array(get_diff(skew_deblur, UAG_pairs))
    names = np.array(get_2base(UAG_pairs))
    order = np.argsort(deblur_list)[::-1]
    blur_list = blur_list[order]
    deblur_list = deblur_list[order]
    names = names[order]
    plt.figure(figsize=(len(names)/2.0, 5))
    x = range(0,len(names)*2, 2)
    plt.bar(x,blur_list, color='g', edgecolor='white', alpha=0.2)
    x = range(1,len(names)*2+1, 2)
    plt.bar(x,deblur_list, color='g', edgecolor='white', alpha=0.6)
    plt.xticks(x, names, rotation='vertical', fontsize=12)
    plt.ylim((ymin,ymax))
    plt.savefig(odir+"skew_UAG_diff.pdf", bbox_inches='tight')
    plt.close()
    print "mean diff: deblur: {0} blur: {1}".format(np.mean(deblur_list), np.mean(blur_list))
    blur_gcu = np.array(get_diff(skew_blur, GCU_pairs))
    blur_icu = np.array(get_diff(skew_blur, ICU_pairs))
    print "mwu: blur", scipy.stats.mannwhitneyu(blur_gcu, blur_icu)
    deblur_gcu = np.array(get_diff(skew_deblur, GCU_pairs))
    deblur_icu = np.array(get_diff(skew_deblur, ICU_pairs))
    print "mwu: deblur", scipy.stats.mannwhitneyu(deblur_gcu, deblur_icu)

    gcu_time_blur = []
    gcu_time_deblur = []
    for cwc, cwb in GCU_pairs:
        gcu_time_blur.extend([skew_blur[cwc], skew_blur[cwb]])
        gcu_time_deblur.extend([skew_deblur[cwc], skew_deblur[cwb]])
    icu_time_blur = []
    icu_time_deblur = []
    for cwc, cwb in ICU_pairs:
        icu_time_blur.extend([skew_blur[cwc], skew_blur[cwb]])
        icu_time_deblur.extend([skew_deblur[cwc], skew_deblur[cwb]])
    uag_time_blur = []
    uag_time_deblur = []
    for cwc, cwb in UAG_pairs:
        uag_time_blur.extend([skew_blur[cwc], skew_blur[cwb]])
        uag_time_deblur.extend([skew_deblur[cwc], skew_deblur[cwb]])
    print "gcu, icu, uag"
    print np.mean(gcu_time_blur), np.mean(icu_time_blur), np.mean(uag_time_blur)
    print np.mean(gcu_time_deblur), np.mean(icu_time_deblur), np.mean(uag_time_blur)    
