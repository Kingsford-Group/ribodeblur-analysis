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

def plot_codon_hist_fitting(rc_list, params, pdf_func, xmin=1e-7, xmax=3, fn_prefix="fp"):
    print "plot histograms"
    fig = plt.figure(facecolor='white', figsize=(16, 16))
    clist = generate_codon_list_alphabetical()
    for i in xrange(len(clist)):
        codon = clist[i]
        if codon not in rc_list: continue
        nrc_list = np.array(rc_list[codon])
        ax = fig.add_subplot(8, 8, i)
        ns,bins,patches = plt.hist(nrc_list, np.linspace(xmin,xmax,50), normed = True, histtype='stepfilled', log=False, color = 'c', edgecolor='c', alpha=0.4)
        ymax = max(ns)
        x = np.linspace(xmin, xmax, 10000)
        y = pdf_func(x, params[codon])
        plt.plot(x,y,'r-',lw=2, alpha=0.5)
        ax.set_title(codon,fontsize=10)
        ax.set_xlim((0,xmax))
        ax.set_ylim((0,ymax+0.1))
        plt.xticks(range(0,xmax+1))
    plt.tight_layout()
    plt.savefig(fn_prefix+"_codon_hist.png", bbox_inches='tight')
    plt.close()

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
    
    from codon_table import codon2aa, aa2shortname
    base_list = ['T', 'C', 'A', 'G']
    data_tab = []
    name_tab = []
    for b1 in base_list:
        for b3 in base_list:
            data_row = []
            name_row = []
            for b2 in base_list:
                codon = b1+b2+b3
                if codon in skew_deblur:
                    data_row.append(skew_deblur[codon])
                else:
                    data_row.append(0)
                name_row.append('{0} ({1})'.format(codon, aa2shortname[codon2aa[codon]]))
            data_tab.append(data_row)
            name_tab.append(name_row)
            
    data_tab = np.array(data_tab)
    print data_tab.shape
    vmin=np.min(data_tab[data_tab!=0])
    vmax=np.max(data_tab)
    print vmin, vmax 
    fig = plt.figure()
    ax = fig.add_subplot(111)
    heatmap = plt.pcolor(data_tab, cmap=matplotlib.cm.Blues, vmin=vmin, vmax=vmax)
    for y in range(data_tab.shape[0]):
        for x in range(data_tab.shape[1]):
            if data_tab[y][x] < (2*vmax-vmin)/2.0 :
                c = (0,0,0)
            else:
                c = (1,1,1)
            plt.text(x+0.5, y+0.5, name_tab[y][x], horizontalalignment='center', verticalalignment='center', color=c)
    plt.colorbar(heatmap)
    plt.xticks([])
    plt.yticks([])
    plt.savefig('test.pdf', bbox_inches='tight')
    

