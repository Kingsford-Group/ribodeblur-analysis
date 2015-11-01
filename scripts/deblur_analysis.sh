hist_fn=../outputs/BY_FP_rlen_5p.hist
vblur_fn=../figures/Albert14/BY_FP_rlen_5p.vblur
cds_range=~/scratch/scratch2/ribomap-playground/yeast/ref/cds_range.txt
output_dir=../figures/test/
python recover_asite_profile.py ${hist_fn} ${vblur_fn} ${cds_range} ${output_dir}
exit
hist_fn=../outputs/BY_FP_hc.hist
vblur_fn=../figures/Albert14/BY_FP_rlen_5p.vblur
cds_range=~/scratch/scratch2/ribomap-playground/yeast/ref/cds_range.txt
output_dir=../figures/test/


python elongation_rate.py ~/scratch/scratch2/ribomap-playground/yeast/ref/protein_coding_100_filtered.fasta ~/scratch/scratch2/ribomap-playground/yeast/ref/cds_range.txt ../outputs/BY_FP_hc.hist ../figures/Albert14/BY_FP_rlen_5p.vblur ../figures/Albert14/BY_FP_hc.eps  ../figures/test/
