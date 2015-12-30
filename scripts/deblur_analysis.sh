hist_fn=../outputs/BY_FP_rlen_5p.hist
vblur_fn=../figures/Albert14/BY_FP_rlen_5p.vblur
eps_fn=../figures/Albert14/BY_FP_hc.eps
cds_range=~/scratch/scratch2/ribomap-playground/yeast/ref/cds_range.txt
ref_fa=~/scratch/scratch2/ribomap-playground/yeast/ref/protein_coding_100_filtered.fasta
output_dir=../figures/test/

python meta_base_type.py ${hist_fn} ${ref_fa} ${cds_range} ${output_dir}

exit
python recover_asite_profile.py ${hist_fn} ${vblur_fn} ${cds_range} ${output_dir}
python profile_evaluation.py ${hist_fn} ${cds_range} ${vblur_fn} ${eps_fn} ${output_dir}

hist_fn=../outputs/BY_FP_hc.hist
python elongation_rate.py ${ref_fa} ${cds_range} ${hist_fn} ${vblur_fn} ${eps_fn} ${output_dir}
