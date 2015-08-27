#!/bin/bash
work_dir=/home/hw1/scratch/scratch2/ribomap-playground/yeast/
sra_dir=${work_dir}data/sra/
fasta_dir=${work_dir}data/fasta/
ref_dir=${work_dir}ref/
bin_dir=/home/hw1/scratch/software_testing/ribomap/bin/
export PATH=${bin_dir}:$PATH
contaminant_idx=${work_dir}StarIndex/contaminant/
transcript_idx=${work_dir}StarIndex/by_transcript/
align_dir=${work_dir}alignment/
sm_odir=${work_dir}sm_quant/BY/
riboseq_fq=${fasta_dir}BY_FP.fastq.gz
rnaseq_fq=${fasta_dir}BY_mRNA.fastq.gz
nproc=30
adapter=CTGTAGGCACCATCAAT
#============================================================
# step 1 prepare footprint data
# step 1.1: download data
#============================================================
# BY_mRNA
rnaseq_url=ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX476%2FSRX476344/SRR1177156/SRR1177156.sra
# BY_FP
riboseq_url=ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX476%2FSRX476345/SRR1177157/SRR1177157.sra
wget -P ${sra_dir} -N ${rnaseq_url}
wget -P ${sra_dir} -N ${riboseq_url}
#============================================================
# step 1.2: convert sra to fastq
#============================================================
echo "converting sra to fastq..."
fastq-dump ${sra_dir}SRR1177156.sra -O ${fasta_dir} --gzip
fastq-dump ${sra_dir}SRR1177157.sra -O ${fasta_dir} --gzip
#============================================================
# step 1.3: rename files
#============================================================
echo "renaming fastq files to be more informative..."
mv ${fasta_dir}SRR1177156.fastq.gz ${rnaseq_fq}
mv ${fasta_dir}SRR1177157.fastq.gz ${riboseq_fq}
#============================================================
# step 2 prepare reference
# step 2.1: download data
#============================================================
# contaminants
contaminant_url=http://downloads.yeastgenome.org/sequence/S288C_reference/rna/rna_coding.fasta.gz
wget -P ${ref_dir} -N ${contaminant_url}
# transcriptome
wget -P ${ref_dir} -N http://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding.fasta.gz
wget -P ${ref_dir} -N http://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_genomic_1000.fasta.gz
wget -P ${ref_dir} -N http://downloads.yeastgenome.org/sequence/S288C_reference/orf_protein/orf_trans.fasta.gz
gunzip -f ${ref_dir}*.gz
python filter_yeast_transcript.py ${ref_dir} 100 #according to Joel
contaminant_fa=${ref_dir}`basename ${contaminant_url} .gz`
transcript_fa=${ref_dir}protein_coding_100_filtered.fasta
#============================================================
# step 2.2: build star index
#============================================================
mkdir -p ${contaminant_idx}
mkdir -p ${transcript_idx}
STAR --runThreadN $nproc --runMode genomeGenerate --genomeDir ${contaminant_idx} --genomeFastaFiles ${contaminant_fa} --genomeSAindexNbases 8 --genomeChrBinNbits 11
STAR --runThreadN $nproc --runMode genomeGenerate --genomeDir ${transcript_idx} --genomeFastaFiles ${transcript_fa} --genomeSAindexNbases 11 --genomeChrBinNbits 12
#============================================================
# step 3 align reads
# step 3.1: filter rrna
#============================================================
echo "filtering contaminated sequences in riboseq"
common_params="--runThreadN ${nproc} --clip3pAdapterSeq ${adapter} --seedSearchLmax 10 --outFilterMultimapScoreRange 0 --outFilterMultimapNmax 255 --outFilterMismatchNmax 1 --outFilterIntronMotifs RemoveNoncanonical"
rna_core=`basename ${rnaseq_fq}`
rna_core=${rna_core%%.*}
ornaprefix=${align_dir}${rna_core}_rrna_
rna_nrrna_fa=${ornaprefix}Unmapped.out.mate1
STAR --genomeDir ${contaminant_idx} --readFilesIn ${rnaseq_fq} --outFileNamePrefix ${ornaprefix} --outStd SAM --outReadsUnmapped Fastx --outSAMmode NoQS "--readFilesCommand zcat <" ${common_params} > /dev/null
ribo_core=`basename ${riboseq_fq}`
ribo_core=${ribo_core%%.*}
oriboprefix=${align_dir}${ribo_core}_rrna_
ribo_nrrna_fa=${oriboprefix}Unmapped.out.mate1
STAR --genomeDir ${contaminant_idx} --readFilesIn ${riboseq_fq} --outFileNamePrefix ${oriboprefix} --outStd SAM --outReadsUnmapped Fastx --outSAMmode NoQS "--readFilesCommand zcat <" ${common_params} > /dev/null
#============================================================
# step 3.2: map to transcriptome
#============================================================
echo "aligning riboseq to the transcriptome"
ornaprefix=${align_dir}${rna_core}_transcript_
STAR --genomeDir ${transcript_idx} --readFilesIn ${rna_nrrna_fa} --outFileNamePrefix ${ornaprefix} --outSAMtype BAM Unsorted --outSAMmode NoQS --outSAMattributes NH NM ${common_params}
rna_bam=${ornaprefix}Aligned.out.bam
oriboprefix=${align_dir}${ribo_core}_transcript_
STAR --genomeDir ${transcript_idx} --readFilesIn ${ribo_nrrna_fa} --outFileNamePrefix ${oriboprefix} --outSAMtype BAM Unsorted --outSAMmode NoQS --outSAMattributes NH NM ${common_params}
