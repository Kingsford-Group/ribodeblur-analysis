#!/bin/bash
# pre-compiled STAR indices should be provided otherwise STAR won't run
riboseq_fq=$1
align_dir=$2
adapter=$3
if [ -z "${riboseq_fq}" ] || [ -z "${align_dir}" ]; then
    echo "Usage: ./read_align.sh ribo-seq.fq.gz output_dir adapter(optional)"
    exit
fi

if [ -z "${adapter}" ]; then
    adapter=N
fi
#=============================
# default parameters
#=============================
cur_dir=`dirname $0`
bin_dir=/home/hw1/scratch/software_testing/ribomap/bin/
export PATH=${bin_dir}:$PATH
mkdir -p ${align_dir}
contaminant_idx=/home/hw1/scratch/scratch2/ribomap-playground/yeast/StarIndex/contaminant/
transcript_idx=/home/hw1/scratch/scratch2/ribomap-playground/yeast/StarIndex/by_transcript/
nproc=30
#============================================
# step 1: filter rrna
#============================================
echo "filtering contaminated sequences in riboseq"
ribo_core=`basename ${riboseq_fq}`
ribo_core=${ribo_core%%.*}
common_params="--runThreadN ${nproc} --clip3pAdapterSeq ${adapter} --seedSearchLmax 10 --outFilterMultimapScoreRange 0 --outFilterMultimapNmax 255 --outFilterMismatchNmax 1 --outFilterIntronMotifs RemoveNoncanonical"
oriboprefix=${align_dir}${ribo_core}_rrna_
ribo_nrrna_fa=${oriboprefix}Unmapped.out.mate1
STAR --genomeDir ${contaminant_idx} --readFilesIn ${riboseq_fq} --outFileNamePrefix ${oriboprefix} --outStd SAM --outReadsUnmapped Fastx --outSAMmode NoQS "--readFilesCommand zcat <" ${common_params} > /dev/null 
#============================================
# step 2: map to transcriptome
#============================================
echo "aligning riboseq to the transcriptome"
oriboprefix=${align_dir}${ribo_core}_transcript_
STAR --genomeDir ${transcript_idx} --readFilesIn ${ribo_nrrna_fa} --outFileNamePrefix ${oriboprefix} --outSAMtype BAM Unsorted --outSAMmode NoQS --outSAMattributes NH NM ${common_params}
