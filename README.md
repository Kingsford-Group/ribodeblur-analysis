Overview
------
Ribosome profiling quantitatively captures ribosome footprints during translation. The resulting profiles of ribosome locations are widely used to study translational speed. However, an accurate estimation of the ribosome location depends on identifying the actively translated A-site from ribosome profiling reads, a problem that was previously unsolved. Here, we propose a novel method to estimate the ribosome A-site position signals from ribosome profiling data. Our model, allows more footprint data to be used, accurately explains the 3-nt periodicity of ribosome profiling reads of various lengths, and results in profiles that are correctly highly skewed towards a single frame within a codon, while being consistent across different read lengths. The method retains sub-codon resolution in the recovered A-site position signals, making it possible to detect off-frame translational events, such as frameshifting. Using these refined profiles, we show that wobble pairing codons are translated slower than their synonymous codons with Watson-Crick pairing. Such results provide evidence that protein synthetic rate can be tuned by synonymous codon usage bias.

Citation
------
This page provides access to the source code and the instructions for running the pipelines to produce the results in the following manuscript:

__Hao Wang, Joel McManus and Carl Kingsford__. *Accurate recovery of ribosome position signals reveals slow translation of wobble-pairing codons in yeast*. Under review (2015).

Deblur pipeline
------
### Data preprocessing
#### Prerequisites
* [SRA Toolkit](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software)
* [STAR (v2.4.0j)](https://github.com/alexdobin/STAR/releases/tag/STAR_2.4.0j)
* python 2.7 
* python package [Biopython](http://biopython.org/wiki/Main_Page)
* C++ compiler that support c++11 features (for instance g++ >= 4.7)
* C++ library [seqan (v1.4.1)](http://www.seqan.de/)

#### Compile C++ source code

     cd src
     make read_len_hist INC="-I/opt/local/include"

#### Processing Pipeline
1. Download data, align reads to transcriptome by running `scripts/BY_pipeline.sh`.
Please change in `BY_pipeline.sh` the line of `$work_dir` to be the working directory of your choice, and the line of `$bin_dir` to be where STAR is installed. The reference fasta is saved as `protein_coding_100_filtered.fasta` under the directory `$work_dir/ref/`. This reference fasta is needed for estimating the decoding time from the ribosome profilng data. A CDS range file `cds_range.txt` is also automatically generated under the directory `$work_dir/ref/`. This is a plain text file that defines the start and stop of the coding regions (CDS) for each transcript. This file is needed for all of the analysis below.
2. Prepare alignments to be a list of position, read count list by running `src/read_len_hist align.bam input_rlen.hist`. This produces file `input_rlen.hist`. The format of `input_rlen.hist` can be found below.
3. Plot meta profiles: `python scripts/meta_profile.py input_rlen.hist cds_range.txt output_dir`. This produces all subfigures for Figure S1. 

#### Convertion from bam alignment to `input_rlen.hist`
Since ribo-seq data is usually very sparse, the alignment bam file is converted into a `*.hist` file with the following entries:
~~~~~
	refID: 0
	tid: YAL001C
	9: 1 1,
	16: 2500 1,
	18: 97 1, 143 1,
	19: 97 16, 151 2, 157 1, ...
	20: 105 1, 162 1, 672 1, ...
	21: 149 1, 311 3, 338 1, ...
	...
~~~~~
`refID` is the reference id of a transcript in the fasta file, `tid` is the transcript name appeared in the reference fasta, numbers before the colon is the read length, then pairs of read start position and read count are listed behind the read length, each position-count pair is separated by comma.

### Deblur testing
#### Prerequisites
* python package: matplotlib, numpy, scipy, and Biopython

#### Processing Pipeline
1. Create high coverage profiles with at least 50% loci have nonzero counts: `python scripts/high_cover_profile.py input_rlen.hist cds_range.txt coverage_ratio count_threshold high_coverage.hist`, where `coverage_ratio = 0.5`, `count_threshold = 0`. The high coverage profiles are output to `high_coverage.hist`. 
2. Train blur vector on meta profiles: `python scripts/train_vblur_from_meta.py input_rlen.hist cds_range.txt output_dir`. This outputs the learned blur vector for each read length to file `input_rlen.vblur` under folder `output_dir`.
3. Deblur individual transcript profiles: `python scripts/deblur_transcripts.py high_coverage.hist input_rlen.vblur cds_range.txt output_dir`. This outputs the consensus profile and the deviations to `high_coverage.eps` under folder `output_dir`.
4. Evaluating deblur results: `python scripts/profile_evaluation.py high_coverage.hist cds_range.txt input_rlen.vblur high_coverage.eps output_dir`. This loads in pre-computed deblur results and plots evaluation results to `output_dir`. This generates Figure 2, 3, S5 in the paper.

#### blur vectors in `*.vblur`
Each line starts with the read length, then a colon, then a vector of numbers representing the probability of blurring effect learned from the meta profiles for that read length.

#### consensus and deviations in `*.eps`
Each transcript entry starts with a line like `tid: YBR177C`, then a list of deviations for each read length for this transcript is provided with `read length : deviation vector`, where the consensus is given when `read length = 0`.

Synthetic frameshift pipeline
------
1. create synthetic frameshifts by running `python scripts/synthetic_frameshift.py high_coverage.hist cds_range.txt high_coverage_fs.hist`. This produces a set of frameshifted profiles on profiles with high coverage. The outputs are: `high_coverage_fs.hist` for the frameshifted profiles, and `high_coverage_fs.fsp` for the locations of the frameshift point.
2. Deblur these profiles: `python scripts/deblur_transcripts.py high_coverage_fs.hist input_rlen.vblur cds_range.txt output_dir`. This outputs the consensus profile and the deviations to `high_coverage_fs.eps` under folder `output_dir`.
1. Evaluate frame distribution: `python scripts/frameshift_evaluation.py high_coverage_fs.hist cds_range.txt input_rlen.vblur high_coverage_fs.eps high_coverage_fs.fsp output_dir`. This produces Figure 4 in the paper.

Codon decoding rate estimation pipeline
------
`python scripts/elongation_rate.py transcript.fasta cds_range.txt high_coverage.hist input_rlen.vblur high_coverage.eps output_dir`. This produces Figure 5 in the paper.
