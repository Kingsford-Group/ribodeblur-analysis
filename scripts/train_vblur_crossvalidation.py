#!/usr/bin/env python
import scipy.stats
from meta_profile import *
from deblur_utils import *
from global_params import *
from operator import itemgetter
colormap = plt.cm.gist_rainbow

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.size'] = 20
rcParams['xtick.major.size'] = 5
rcParams['ytick.major.size'] = 5
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42

#=============================
# train on subset
#=============================
def train_vblur_on_subsample_loc(cobs,estep):
    # shuffle pobs order and randomly select 1/5 loci for model training
    np.random.seed(970743242)
    order = range(len(cobs[28]))
    np.random.shuffle(order)
    interval = len(order)/5
    abd = get_abundance(cobs)
    blist = []
    for i in xrange(6):
        # reset ptrue
        ptrue = initiate_ptrue(cobs)
        if i == 5:
            pos_list = None
        else:
            pos_list = order[i*interval:(i+1)*interval]
        # blur kernel initial estimation
        b = {}
        for rlen in cobs:
            threshold = np.percentile(cobs[rlen], percentile)
            k = rlen-28
            ctrue = ptrue * abd[rlen]
            vblur = single_kernel_width(ctrue, cobs[rlen], k, threshold, pos_list, False)
            b[rlen] = vblur
        # train vblur on a subset of loci    
        b,ptrue,eps = train_blur_vec(cobs, ptrue, abd, b, low, percentile, converge_cutoff, pos_list, estep, None)
        blist.append(b)
    return blist

def train_vblur_on_subsample_transcript(tid_list, tlist, cds_range, ibegin, iend, estep):
    tid_list = np.array(list(tid_list))
    # shuffle transcript and randomly select 1/5 transcripts for model training
    np.random.seed(654261789)
    order = range(len(tid_list))
    np.random.shuffle(order)
    interval = len(order)/5
    pos_list = None
    blist = []
    for i in xrange(6):
        if i == 5:
            tid_select = set(tid_list)
        else:
            select = order[i*interval:(i+1)*interval]
            tid_select = set(tid_list[select])                  
        pos_hist = create_rlen_meta_profile(tlist, cds_range, tid_select, ibegin, iend)
        cobs = get_cobs(pos_hist, rlen_min, rlen_max, 0, imax)
        abd = get_abundance(cobs)
        ptrue = initiate_ptrue(cobs)
        # blur kernel initial estimation
        b = {}
        for rlen in cobs:
            threshold = np.percentile(cobs[rlen], percentile)
            k = rlen-28
            ctrue = ptrue * abd[rlen]
            vblur = single_kernel_width(ctrue, cobs[rlen], k, threshold, pos_list, False)
            b[rlen] = vblur
        # train vblur on a subset of transcripts
        b,ptrue,eps = train_blur_vec(cobs, ptrue, abd, b, low, percentile, converge_cutoff, pos_list, estep, None)
        blist.append(b)
    return blist

#=============================
# plotting
#=============================
def compute_square_error(vblur,k,pobs,ptrue):
    """
    sum_i | pobs[i] - ptrue[i-k]*vblur |_2
    only used in the training step
    """
    pobs_estimate = np.zeros(len(pobs))
    w = len(vblur)
    blur_start, blur_end = get_blur_range(w, k, len(ptrue))
    true_start, true_end = get_true_range(w, k, len(ptrue))
    b = w - 1
    if w%2 == 0:
        start = b+1
    else:
        start = b
    pobs_estimate[blur_start: blur_end] = np.convolve(vblur, ptrue[true_start:true_end])[start:-b]
    res = pobs[blur_start: blur_end]-pobs_estimate[blur_start: blur_end]
    lse = np.sqrt(np.dot(res,res))
    return lse

def square_error_list(pobs_list, ptrue, blur_list):
    obj_list = []
    for rlen in pobs_list:
        k = rlen - 28
        obj = compute_square_error(blur_list[rlen], k, pobs_list[rlen], ptrue)
        obj_list.append(obj)
    return obj_list

def evaluate_vblur_obj(pobs, ptrue, blist, ofname=None):
    obj_list = []
    for b in blist:
        objs = square_error_list(pobs, ptrue, b)
        obj_list.append(objs)
    plt_cnt = len(blist)
    rlen_list = blist[0].keys()
    plt.figure()
    c = [ colormap(i) for i in np.linspace(0,1,plt_cnt) ]
    for i in xrange(len(obj_list)):
        if i != len(obj_list)-1:
            plt.scatter(rlen_list, obj_list[i], marker='+', alpha=0.5, color=c[i], s=100)
        else:
            plt.scatter(rlen_list, obj_list[i], marker='o', alpha=0.5, color=c[i], edgecolor='None', s=100)
    plt.xlabel("read length")
    plt.ylabel("least square")
    plt.xlim((rlen_list[0]-0.1, rlen_list[-1]+0.1))
    if ofname!=None:
        plt.savefig(ofname, bbox_inches='tight')
        plt.close()
    else:
        plt.show()    

def evaluate_vblur_pearsonr(pobs, ptrue, blist, ofname=None):
    pr_list = []
    for b in blist:
        prs = []
        for rlen in pobs:
            k = rlen-28
            pobs_est = estimate_pobs_single(b[rlen],k, ptrue)
            blur_start, blur_end = get_blur_range(len(b[rlen]), k, len(pobs[rlen]), True)
            pr = scipy.stats.pearsonr(pobs_est[blur_start:blur_end], pobs[rlen][blur_start:blur_end])[0]
            prs.append(pr)
        pr_list.append(prs)

    rlen_list = pobs.keys()
    plt_cnt = len(blist)
    c = [ colormap(i) for i in np.linspace(0,1,plt_cnt) ]
    fig = plt.figure()
    for i in xrange(len(pr_list)):
        if i!=len(blist)-1:
            plt.scatter(rlen_list, pr_list[i], marker='+', alpha=0.5, c=c[i], s=100)
        else:
            plt.scatter(rlen_list, pr_list[i], marker='o', alpha=0.5, c=c[i], s=100, edgecolor='None')
    plt.ylim((0, 1.1))
    plt.xlim((rlen_list[0]-0.1,rlen_list[-1]+0.1))
    plt.xlabel('read length')
    plt.ylabel('Pearson correlation')
    if ofname!=None:
        plt.savefig(ofname, bbox_inches='tight')
        plt.close()
    else:
        plt.show()    

def evaluate_vblur_pearsonr_frame(pobs, ptrue, blist, ofname=None):
    pr_list = []
    for b in blist:
        prs = []
        for rlen in pobs:
            k = rlen-28
            pobs_est = estimate_pobs_single(b[rlen], k, ptrue)
            pr = evaluate_correlation(pobs[rlen], pobs_est)
            prs.append(pr)
        frame_list = [ map(itemgetter(i), prs) for i in xrange(4) ]
        pr_list.append(frame_list)

    rlen_list = pobs.keys()
    plt_cnt = len(blist)
    c = [ colormap(i) for i in np.linspace(0,1,plt_cnt) ]
    t = ['f0', 'f1', 'f2', 'all']
    fig = plt.figure()
    for i in xrange(4):
        ax = fig.add_subplot(2,2,i+1)
        for j in xrange(len(pr_list)):
            f = pr_list[j][i]
            plt.plot(rlen_list, f, marker='o', alpha=0.5, c=c[j])
        plt.title(t[i])
        plt.ylim((-0.5, 1.1))
        plt.xlim((13,32))
    if ofname!=None:
        plt.savefig(ofname, bbox_inches='tight')
        plt.close()
    else:
        plt.show()    

def plot_vblur_batch_align(blist, ofname=None):
    plot_cnt = len(blist)
    c = [ colormap(i) for i in np.linspace(0,1,plot_cnt) ]
    j = 1
    fig = plt.figure(figsize=(10,13))
    rlen_list = blist[0].keys()
    for rlen in rlen_list:
        ax = fig.add_subplot(6,3,j)
        for i in xrange(plot_cnt):
            vblur = blist[i][rlen]
            if i == 0:
                x = range(len(vblur))
                center = vblur
            else:
                k = get_shift(center, vblur)
                x = range(k, k+len(vblur))
            plt.plot(x, vblur, '-o', color=c[i], alpha=0.5)
        plt.title("{0}".format(rlen))
        j += 1
    plt.tight_layout()
    if ofname!=None:
        plt.savefig(ofname, bbox_inches='tight')
        plt.close()
    else:
        plt.show()    
    
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print "Usage: python train_vblur_crossvalidation rlen.hist cds.txt odir"
        exit(1)
    hist_fn = sys.argv[1]
    cds_txt = sys.argv[2]
    odir = sys.argv[3]
    cds_range = get_cds_range(cds_txt)
    tid_select = filter_transcript_by_length(cds_range, imax)
    tlist = parse_rlen_hist(hist_fn)
    pos_hist = create_rlen_meta_profile(tlist, cds_range, tid_select, utr5_offset, imax)
    mobs = get_cobs(pos_hist, rlen_min, rlen_max, 0, imax)
    abd = get_abundance(mobs)
    pobs = { rlen : prof/float(abd[rlen]) for rlen,prof in mobs.iteritems() }
    ptrue = initiate_ptrue(mobs)

    print "train vblur on subloci with e-step"
    blist = train_vblur_on_subsample_loc(mobs, estep=True)
    evaluate_vblur_obj(pobs, ptrue, blist, odir+"subloc_pobs_obj.pdf")
    evaluate_vblur_pearsonr(pobs, ptrue, blist, odir+"subloc_pobs_pr.pdf")

    print "train vblur on sub transcripts with e-step"    
    blist = train_vblur_on_subsample_transcript(tid_select, tlist, cds_range, utr5_offset, imax, estep=True)
    evaluate_vblur_obj(pobs, ptrue, blist, odir+"subtrans_pobs_obj.pdf")
    evaluate_vblur_pearsonr(pobs, ptrue, blist, odir+"subtrans_pobs_pr.pdf")
    
    for rlen in mobs:
        mobs_estimate = estimate_pobs_single(blist[-1][rlen], klist[rlen], ptrue*abd[rlen])
        plot_profile(mobs[rlen], ptrue*abd[rlen], mobs_estimate, ptrue*abd[rlen], rlen, "{0}meta_cmp_{1}.pdf".format(odir, rlen))
            



