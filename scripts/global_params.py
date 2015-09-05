#!/usr/bin/env python
utr5_offset = -24
utr3_offset = 0
asite_offset = 15
imax = 350  # right most position form the 5' end to sample for historgram
converge_cutoff = 1e-2 # minimum threshold for the euclidean distance of vblur between iterations
rlen_min = 14
rlen_max = 35
low = -1 # cobs filter lower bound
percentile = 98.35 # cobs outlier filter percentile
nproc = 30 # number of processes to use
# for filter highly covered profiles
cover_ratio = 0.5
cnt_threshould = 0
klist = {rlen: 28-rlen for rlen in xrange(rlen_min, rlen_max+1) }
lowest_frame_percent = 0.4 # threshold to accept a read length
lowest_frame_cnt = 10 # meta profile has to have enough coverage to train blur vector
