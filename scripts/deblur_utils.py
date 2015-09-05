#!/usr/bin/env python
import scipy.optimize
import scipy.sparse.linalg
import numpy as np
from meta_profile import *

#=============================
# utils
#=============================
def get_abundance(pobs):
    return { rlen: sum(p) for rlen, p in pobs.iteritems() }

def merge_profiles(plist):
    return np.sum(plist.values(), axis=0)

def initiate_ptrue(pobs):
    if 28 in pobs:
        ptrue_merge = np.zeros(len(pobs[28]))
        ptrue_merge[0::3] = pobs[28][0::3]
    else:
        ptrue_merge = np.zeros(len(pobs.values()[0]))
        psum = merge_profiles(pobs)
        ptrue_merge[0::3] = psum[0::3]
    if np.all(ptrue_merge==0):
        return ptrue_merge
    else:
        return ptrue_merge/sum(ptrue_merge)

def estimate_pobs_single(vblur, k, ptrue):
    A = build_blur_A(len(vblur), k, ptrue)
    pobs_estimate = np.dot(A, vblur)
    return pobs_estimate

def estimate_ctrue(ptrue, eps, cobs):
    ctrue = {}
    rlen_list = list(set(eps.keys()) & set(cobs.keys()))
    for rlen in rlen_list:
        ctruer = sum(cobs[rlen]) * ( ptrue - eps[rlen])
        # adjust negatives        
        ctruer[ctruer<0] = 0
        ctrue[rlen] = ctruer
    return ctrue

def batch_build_ctrue(ptrue, eps, cobs):
    return { tid: estimate_ctrue(ptrue[tid], eps[tid], cobs[tid]) for tid in ptrue }

def batch_distance(tid_list, cobs, ctrue):
    dist = []
    for tid in tid_list:
        obs_merge = merge_profiles(cobs[tid])
        true_merge = merge_profiles(ctrue[tid])
        diff = obs_merge - true_merge
        dist.append( np.linalg.norm(diff)/sum(obs_merge) )
    return np.array(dist)

def compute_total_least_square(vblur, k, pobs, ptrue, eps, selected, train=True):
    """
    sum_i (pobs[i] - ptrue[i-k]*vblur)^2 + sum_j (ptrue[j]-ptrue_merge[j])^2 
    """
    ptrue_estimate = ptrue - eps
    pobs_estimate = estimate_pobs_single(vblur, k, ptrue_estimate)
    true_start, true_end = get_true_range(len(vblur), k, len(ptrue), train)
    eps_obs = (pobs-pobs_estimate)[selected]
    lsq_obs = np.dot(eps_obs,eps_obs)
    eps_true = eps[true_start: true_end]
    lsq_true = np.dot(eps_true, eps_true)    
    lsq_tot = lsq_obs + lsq_true
    # print "in tls: lsq_obs: {0:.2e} lsq_true: {1:.2e} lsq_tot: {2:.2e}".format(lsq_obs, lsq_true, lsq_tot)
    return lsq_tot, lsq_obs, lsq_true

def euclidean_distance_after_extend_short_vec(vs, vl):
    ls = len(vs)
    ll = len(vl)
    vnew = np.zeros(ll)
    start = (ll-ls)/2
    if (ll-ls)%2 == 0:
        vnew[start:-start] = vs
        res = vnew - vl
        return np.sqrt(np.dot(res, res))
    else:
        vnew[start:start+ls] = vs
        res_pre = vnew - vl
        d_pre = np.dot(res_pre, res_pre)
        vpost = np.zeros(ll)
        vpost[start+1:start+1+ls] = vs
        res_post = vpost - vl
        d_post = np.dot(res_post, res_post)
        return np.sqrt(d_pre) if d_pre < d_post else np.sqrt(d_post)
        
def check_convergence(vec, vec_pre, eps):
    l = len(vec)
    l_pre = len(vec_pre)
    if l < l_pre:
        lse = euclidean_distance_after_extend_short_vec(vec, vec_pre)
    elif l > l_pre:
        lse = euclidean_distance_after_extend_short_vec(vec_pre, vec)
    else:
        res = vec - vec_pre
        lse = np.sqrt(np.dot(res, res))
    return lse < eps

#=============================
# Methods evaluation
#=============================
def evaluate_correlation(true_vec, est_vec):
    import scipy.stats
    fdic = get_frames(true_vec)
    pr_list = np.zeros(4)
    for i in xrange(3):
        pof = true_vec[i::3]
        poef = est_vec[i::3]
        # poef!=0 filter out non-recovered zero entries at the beginiing and end
        pearsonr = scipy.stats.pearsonr(poef[poef!=0], pof[poef!=0])[0]
        spearr = scipy.stats.spearmanr(poef[poef!=0], pof[poef!=0])[0]
        print "{0}(f{1}): pr {2:.2} sr {3:.2}".format(i, fdic[i], pearsonr, spearr)
        pr_list[fdic[i]] = pearsonr
    pearsonr = scipy.stats.pearsonr(true_vec[est_vec!=0], est_vec[est_vec!=0])[0]
    spearr = scipy.stats.spearmanr(true_vec[est_vec!=0], est_vec[est_vec!=0])[0]
    pr_list[3] = pearsonr
    print "all: pr {1:.2} sr {2:.2}".format(i, pearsonr, spearr)
    return pr_list

def plot_profile(pobs, ptrue, pobs_estimate, ptrue_estimate, rlen=28, fn=None):
    fig_width = len(ptrue)/10.0*1.5
    fig = plt.figure(figsize=(10,12))
    ax = fig.add_subplot(4,1,1)
    in_range = lambda i, start, end: (i>=start) and (i<end)
    color = [ 'b' if pobs_estimate[i]!=0  else 'k' for i in xrange(len(pobs_estimate)) ]
    plt.bar(np.arange(len(pobs))-0.4, pobs, width=0.8, color=color, edgecolor='white', alpha=0.3)
    plt.title("read length: {0} \nobserved counts".format(rlen), fontsize=12)
    plt.xlim((-1, len(pobs)+1))
    max_cnt = get_ylim(pobs) + 1
    plt.ylim((0, max_cnt))
    ax = fig.add_subplot(4,1,2)
    plt.bar(np.arange(len(pobs_estimate))-0.4, pobs_estimate, width=0.8, color='b', edgecolor='white', alpha=0.3)
    plt.title("estimated observed counts", fontsize=12)
    plt.xlim((-1, len(pobs)+1))
    ax = fig.add_subplot(4,1,3)
    color = [ 'b' if ptrue[i]!=ptrue_estimate[i] else 'k' for i in xrange(len(ptrue)) ]
    plt.bar(np.arange(len(ptrue))-0.4, ptrue, width=0.8, color=color, edgecolor='white', alpha=0.3)
    plt.title("true counts", fontsize=12)
    plt.xlim((-1,len(ptrue)+1))
    ax = fig.add_subplot(4,1,4)
    plt.bar(np.arange(len(ptrue))-0.4, ptrue_estimate, width=0.8, color=color, edgecolor='white', alpha=0.3)
    plt.title("estimated true counts", fontsize=12)
    plt.xlim((-1,len(ptrue)+1))
    plt.tight_layout()
    if fn == None:
        plt.show()
    else:
        plt.savefig(fn, bbox_inches='tight')
        plt.close()
    
#=============================
# blur kernel estimation
#=============================
def get_true_range(blur_width, k, ptrue_len, train=True):
    if train == False:
        return 0, ptrue_len
    half_width = blur_width/2
    blur_start, blur_end = get_blur_range(blur_width, k, ptrue_len)
    true_start = blur_start - half_width - k
    true_end = blur_end + half_width - k
    return true_start, true_end

def get_blur_range(blur_width, k, ptrue_len, train=True):
    if train == False:
        return 0, ptrue_len
    half_width = blur_width/2
    idx_start = max(0, half_width + k)
    idx_end = min(ptrue_len - half_width + k, ptrue_len)
    return idx_start, idx_end

def build_blur_A(blur_width, k, ptrue):
    """
    build A matrix based on ptrue for nnls
    so that vblur is shifted to be around the center of the location
    pblur[i] = sum_j(-w,w) vblur[j]ptrue[i-j]
    adjust j to be (0,2w+1):
    pblur[i] = sum_j(0,2w+1) vblur[j]ptrue[i+w-j]
    where ptrue is shifted right by k to align with pobs
    pblur is built to have the same length as ptrue
    """
    half_width = blur_width/2
    nrow = len(ptrue)
    A = np.zeros((nrow, blur_width))
    for i in xrange(nrow):
        for j in xrange(blur_width):
            itrue = i-j+half_width-k
            if itrue < 0 or itrue >= nrow: continue
            A[i][j] = ptrue[itrue]
    return A

def filter_cobs(cobs, low, high, blur_start, blur_end, pos_list=None):
    # all loci are included inintially
    selected = np.array([True]*len(cobs))
    # include loci in pos_list
    if pos_list != None and len(pos_list)!=0:
        selected[pos_list] = False
        selected = ~selected
    # exclude loci outside of blur_start and blur_end
    selected[:blur_start] = False
    selected[blur_end:] = False
    # exculde outliers
    selected[cobs>high] = False
    selected[cobs<=low] = False
    return selected

def get_normalized_blur_kernel(A, b):
    x, res = scipy.optimize.nnls(A,b)
    x /= np.sum(x)
    res = np.dot(A,x) - b
    obj = np.sqrt(np.sum(res*res))
    # vec_str = " ".join(map(lambda xi: "{:.2f}".format(xi), x))
    # print "obj: {0:.2f}, # data: {1}, vblur: {2}".format(obj,len(b),vec_str)
    return x,obj

def get_normalized_blur_kernal_rerr(A, b):
    A = A[b!=0]
    b = b[b!=0]
    A = (A.T/b).T
    b = np.array([1]*len(b))
    x, res = scipy.optimize.nnls(A,b)
    x /= np.sum(x)
    res = np.dot(A,x) - b
    obj = np.sqrt(np.sum(res*res))
    return x,obj

def single_kernel_width(ptrue, pobs, k, threshold, pos_list=None, re=False):
    """ 
    cannot change interface to passing in selected
    because selected depends on vblur width
    """
    L_opt = np.inf
    vblur_opt = [1]
    # for w in xrange(2,31):
    # fix w to be 31
    for w in [31]:
        A = build_blur_A(w, k, ptrue)
        blur_start, blur_end = get_blur_range(w, k, len(ptrue))
        selected = filter_cobs(pobs, -1, threshold, blur_start, blur_end, pos_list)
        A_filter = A[selected]
        pobs_filter = pobs[selected]
        if re:
            x, L = get_normalized_blur_kernal_rerr(A_filter, pobs_filter)
        else:
            x,L = get_normalized_blur_kernel(A_filter, pobs_filter)
        if L < L_opt:
            L_opt = L
            vblur_opt = x
    # print "in nnls: eps_obs: {0:.2e} data points: {1}".format(L_opt*L_opt, sum(selected))
    # vec_str = " ".join(map(lambda xi: "{:.2f}".format(xi), vblur_opt))
    # print "vblur: [{0}]".format(vec_str)
    return vblur_opt

#=============================
# read count deblur
#=============================
def build_true_A(vblur, k, ptrue_len, pobs_len, train=True):
    """
    build A matrix based on vblur for leastsq
    so that vblur is shifted to be around the center of the location
    pblur[i] = sum_j(-w,w) vblur[j]ptrue[i-j]
    adjust j to be (0,2w+1):
    pblur[i] = sum_j(0,2w+1) vblur[j]ptrue[i+w-j]
    k is used to shift left this center location i to an offset
    """
    blur_width = len(vblur)
    half_width = blur_width/2
    true_start, true_end = get_true_range(blur_width, k, ptrue_len, train)
    A = np.zeros((pobs_len, ptrue_len))
    for i in xrange(pobs_len):
        for j in xrange(blur_width):
            col = i-j+half_width-k
            if col < 0 or col >= pobs_len: continue
            A[i][col] = vblur[j]
    return A[:, true_start:true_end]

def deblur_eps(vblur, k, pobs, ptrue, selected, train=True):
    """
    find the right slack variables (eps) such that
    min |pobs - (ptrue - eps) * vblur |_2^2 + |eps|_2^2
    both pobs and ptrue are normalized
    """
    true_start, true_end = get_true_range(len(vblur), k, len(ptrue), train)
    pobs_estimate = estimate_pobs_single(vblur, k, ptrue)
    # re-normalize pobs_estimate such that the selected ones sum to 1
    if train == False:
        s = sum(pobs_estimate[selected])
        if s == 0:
            return None
        else:
            pobs_estimate = pobs_estimate/float(s)
    A = build_true_A(vblur, k, len(ptrue), len(pobs), train)
    eobs = (pobs_estimate - pobs)
    # only include selected points
    A_filter = A[selected]
    eobs_filter = eobs[selected]
    # solve the least square problem with ridge regression
    sol = scipy.sparse.linalg.lsmr(A_filter, eobs_filter, damp=1, show=False)
    # filling up slack variables
    eps = np.zeros(len(ptrue))
    eps[true_start: true_end] = sol[0]
    # obj_before = np.dot(eobs_filter, eobs_filter)
    # res = np.dot(A_filter, sol[0]) - eobs_filter
    # eps_obs = np.dot(res,res)
    # eps_true = np.dot(eps,eps)
    # obj_after = eps_obs + eps_true
    # print "in deblur_eps: before: {0:.2e}, after: {1:.2e}".format(obj_before, obj_after)
    # print "eps_obs: {0:.2e} eps_true: {1:.2e}".format(eps_obs, eps_true)
    # # re-adjust invalid slacks to be valid
    # negatives = ptrue - eps < 0
    # eps[negatives] = ptrue[negatives]
    # print "data points: total: {0} selected: {1} invalid: {2}".format(len(pobs), sum(selected), sum(negatives))
    # res = np.dot(A_filter, eps[true_start: true_end]) - eobs_filter
    # eps_obs = np.dot(res, res)
    # eps_true = np.dot(eps, eps)
    # obj_after = eps_obs + eps_true
    # print "final obj: eps_obs: {0:.2e} eps_true: {1:.2e} total: {2:.2e}".format(eps_obs, eps_true, obj_after)
    return eps

#=============================
# true signal expectation
#=============================
def expected_ptrue(ptrue_pre, eps, weight):
    # compute average epsilon to mimimize epsilon variance
    eps_expected = sum([ weight[rlen] * e for rlen, e in eps.iteritems() ] )
    eps_expected /= float(sum(weight.values()))
    ptrue = ptrue_pre - eps_expected
    # reset negatives
    ptrue[ptrue<0] = 0
    # re-normalization
    ptrue /= sum(ptrue)
    return ptrue

def plot_obj(tot_list, obs_list, true_list, ofname):
    if ofname==None:
        return False
    plt.figure()
    plt.plot(tot_list, 'r-*', markersize = 20, markeredgecolor = 'none')
    plt.plot(obs_list,'b-o', markersize = 12, markeredgecolor='none')
    plt.plot(true_list, 'c-D', markersize = 10, markeredgecolor='none')
    plt.legend(['eps total', 'eps observed', 'eps true'], loc='upper right', frameon=False, ncol=3)
    plt.xlim((-0.2, len(tot_list)-1))
    plt.xlabel('iteration')
    plt.ylabel('least square')
    plt.savefig(ofname, bbox_inches='tight')
    plt.close()

def select_single_rlen_loci(cobs_rlen, vblur, k, percentile, low, pos_list, train):
    threshold = np.percentile(cobs_rlen, percentile) if percentile<100 else np.inf
    blur_start, blur_end = get_blur_range(len(vblur), k, len(cobs_rlen), train)
    return filter_cobs(cobs_rlen, low, threshold, blur_start, blur_end, pos_list)
    
def select_loci(cobs, b, klist, percentile=98.35, low=-1, train=True, pos_list=None):
    return { rlen: select_single_rlen_loci(cobs_rlen, b[rlen], klist[rlen], percentile, low, pos_list, train)
             for rlen,cobs_rlen in cobs.iteritems() }

def compute_tot_obj(ptrue, abd, pobs, klist, select, b, eps_rlen, train):
    obj_obs = 0
    obj_true = 0
    obj_tot = 0
    for rlen in pobs:
        if eps_rlen == None or rlen not in eps_rlen:
            eps = np.zeros(len(ptrue))
        else:
            eps = eps_rlen[rlen]
        vblur = [1] if b==None else b[rlen]
        eps_tot, eps_obs, eps_true = compute_total_least_square(vblur, klist[rlen], pobs[rlen], ptrue, eps, select[rlen], train)
        obj_obs += abd[rlen] * eps_obs
        obj_true += abd[rlen] * eps_true
        obj_tot += abd[rlen] * eps_tot
    return obj_tot, obj_obs, obj_true

def train_blur_vec(cobs, ptrue, abd, b, klist, low, percentile, converge_cutoff=1e-1, pos_list=None, estep=True, ofname=None):
    """
    EM frame work on training vblur
    Max params: vblur, ptrue_eps 
    Expected latent: ptrue_merge
    """
    obs_list = []
    true_list = []
    tot_list = []
    pobs = { rlen : prof/float(abd[rlen]) for rlen,prof in cobs.iteritems() }
    select = select_loci(cobs, b, klist, percentile, low)
    # initial evaluation of objective function
    obj_tot, obj_obs, obj_true = compute_tot_obj(ptrue, abd, pobs, klist, select, b, None, True)
    obs_list.append(obj_obs)
    true_list.append(obj_true)
    tot_list.append(obj_tot)
    converge_cnt = 0
    i = 0
    while converge_cnt < len(cobs):
        i += 1
        converge_cnt = 0
        sys.stdout.write("vblur training EM iteration: {0}.\t\r".format(i))
        sys.stdout.flush()
        eps_rlen = {}
        # M-step: paramter estimation
        for rlen in pobs:
            vblur = b[rlen]
            k = klist[rlen]
            selected = select[rlen]
            # print "rlen: {0} default k: {1} data points: {2}".format(rlen, klist[rlen], sum(selected))
            # parameter 1: slack vars in ptrue
            # true signal deblur
            eps = deblur_eps(vblur, k, pobs[rlen], ptrue, selected, True)
            # print "original: ",
            # eps_tot, eps_obs, eps_true = compute_total_least_square(vblur, k, pobs[rlen], ptrue, np.zeros(len(ptrue)), selected, True)            
            # print "ptrue step: ",
            # eps_tot,eps_obs,eps_true = compute_total_least_square(vblur, k, pobs[rlen], ptrue, eps, selected, True)
            ctrue_estimated = abd[rlen]*(ptrue - eps)
            # parameter 2: blur kernel re-estimation
            # print "vblur step: ",
            # print "before:",
            # eps_tot,eps_obs,eps_true = compute_total_least_square(vblur, k, cobs[rlen], ctrue_estimated, np.zeros(len(ctrue_estimated)), selected, True)
            threshold = np.percentile(pobs[rlen]*abd[rlen], percentile)
            vblur = single_kernel_width(ctrue_estimated, pobs[rlen]*abd[rlen], k, threshold, pos_list)
            # print "after:",
            # eps_tot,eps_obs,eps_true = compute_total_least_square(vblur,k,cobs[rlen], ctrue_estimated, np.zeros(len(ctrue_estimated)), selected, True)
            # print "final m-step:",
            # eps_tot,eps_obs,eps_true = compute_total_least_square(vblur, k, pobs[rlen], ptrue, eps, selected, True)
            converge_cnt += check_convergence(b[rlen], vblur, converge_cutoff)
            if len(vblur) != len(b[rlen]):
                select[rlen] = select_single_rlen_loci(cobs[rlen], vblur, k, percentile, low, pos_list, True)
            # update blur kernel
            b[rlen] = vblur
            # update slacks for true estimates
            eps_rlen[rlen] = eps
        # keep a record of m-step obj
        #obj_mstep,_,_ = compute_tot_obj(ptrue, abd, pobs, klist, select, b, eps_rlen, True)
        # E-step: ptrue re-estimation
        if estep == True:
            # print "E-step:"
            ptrue_curr = expected_ptrue(ptrue, eps_rlen, abd)
            for rlen in cobs:
                # print "rlen: {0}".format(rlen),
                ptrue_estimate = ptrue - eps_rlen[rlen]
                eps = ptrue_curr - ptrue_estimate
                # eps_tot, eps_obs, eps_true = compute_total_least_square(b[rlen], klist[rlen], pobs[rlen], ptrue_curr, eps, select[rlen], True)
                # update slacks for true estimates (only change in E-step)
                eps_rlen[rlen] = eps
            ptrue = ptrue_curr
            obj_tot, obj_obs, obj_true = compute_tot_obj(ptrue, abd, pobs, klist, select, b, eps_rlen, True)
            # print "{0:.2e} {1:.2e}".format(obj_mstep, obj_tot)
        # keep a record of the objective function from the M-step
        obs_list.append(obj_obs)
        true_list.append(obj_true)
        tot_list.append(obj_tot)
    sys.stdout.write("\n")
    plot_obj(tot_list, obs_list, true_list, ofname)    
    print "final E-step obj: ", tot_list[-1]
    return b, ptrue, eps_rlen

def train_vblur_from_meta_profiles(cobs, klist, low, percentile, converge_cutoff = 1e-1, estep=True, ofname=None):
    """wrapper function to train vblur on meta profiles"""
    abd = get_abundance(cobs)
    ptrue = initiate_ptrue(cobs)
    b = {}
    # HARD CODED include all loci within range [0,iend]
    pos_list = None
    # initialize vblur for different read lengths
    for rlen in cobs:
        # HARD CODED outlier percentile
        threshold = np.percentile(cobs[rlen], percentile)
        # HARD CODED profile shift offset compare to read length 28
        k = klist[rlen]
        ctrue = ptrue * abd[rlen]
        vblur = single_kernel_width(ctrue, cobs[rlen], k, threshold, pos_list, False)
        b[rlen] = vblur
    b, ptrue, eps = train_blur_vec(cobs, ptrue, abd, b, klist, low, percentile, converge_cutoff, pos_list, estep, ofname)
    return b, ptrue, eps

def recover_true_profile(cobs, b, low, percentile, converge_cutoff, ofname=None):
    abd = get_abundance(cobs)
    ptrue = initiate_ptrue(cobs)
    if np.all(ptrue==0):
        print "transcript too sparse"        
        return None, None
    pobs = { rlen : prof/float(abd[rlen]) for rlen,prof in cobs.iteritems() }
    klist = { rlen: rlen - 28 for rlen in cobs }
    select = select_loci(cobs, b, klist, percentile, low, False)
    obs_list = []
    true_list = []
    tot_list = []
    # initial evaluation of objective function
    obj_tot, obj_obs, obj_true = compute_tot_obj(ptrue, abd, pobs, klist, select, b, None, False)
    obs_list.append(obj_obs)
    true_list.append(obj_true)
    tot_list.append(obj_tot)

    # EM framework
    converge = False
    i = 0
    while converge == False:
        i += 1
        # sys.stdout.write("deblur testing EM iteration: {0}.\t\r".format(i))
        # sys.stdout.flush()        
        eps_rlen = {}
        # M-step: paramter estimation
        for rlen in pobs:
            vblur = b[rlen]
            k = rlen - 28
            selected = select[rlen]
            # true signal deblur
            eps = deblur_eps(vblur, k, pobs[rlen], ptrue, selected, False)
            if eps is None: continue
            # print "original: ",
            # eps_tot, eps_obs, eps_true = compute_total_least_square(vblur, k, pobs[rlen], ptrue, np.zeros(len(ptrue)), selected, False)            
            # print "ptrue step: ",
            # eps_tot,eps_obs,eps_true = compute_total_least_square(vblur, k, pobs[rlen], ptrue, eps, selected, False)
            eps_rlen[rlen] = eps
        if len(eps_rlen) == 0:
            print "transcript too sparse"
            return None, None
        obj_mstep,_,_ = compute_tot_obj(ptrue, abd, pobs, klist, select, b, eps_rlen, True)
        # E-step: ptrue re-estimation
        ptrue_curr = expected_ptrue(ptrue, eps_rlen, abd)
        for rlen in eps_rlen:
            ptrue_estimate = ptrue - eps_rlen[rlen]
            eps = ptrue_curr - ptrue_estimate
            # eps_tot, eps_obs, eps_true = compute_total_least_square(b[rlen], klist[rlen], pobs[rlen], ptrue_curr, eps, select[rlen], False)
            # update slacks for true estimates (only change in E-step)
            eps_rlen[rlen] = eps
        ptrue = ptrue_curr
        obj_tot, obj_obs, obj_true = compute_tot_obj(ptrue, abd, pobs, klist, select, b, eps_rlen, False)
        if obs_list[0] < obj_obs : 
            print "deblur cannot model this profile!"
            return None, None
        if (obs_list[-1] - obj_obs)/float(obs_list[-1]) < converge_cutoff:
            converge = True
        # print "{0:.2e} {1:.2e}".format(obj_mstep, obj_tot)        
        # keep a record of the objective function from the M-step
        obs_list.append(obj_obs)
        true_list.append(obj_true)
        tot_list.append(obj_tot)
        if i > 50:
            print "fail to converge"
            return None, None
    # sys.stdout.write("\n")
    plot_obj(tot_list, obs_list, true_list, ofname)    
    # print "final M-step obj: ", tot_list[-1]
    return ptrue, eps_rlen
