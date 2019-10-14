import os
#from Caleb 09.10.2019
# project-specific parameters
threads = 8
project_dir = '/project/klynclab/ranzani_2015/'
bamdir = os.path.join(project_dir, 'bams')
ref_genome = "/project/klynclab/GENOMES/human/hg38_GRCh38_Reference_Genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
star_genome = '/project/klynclab/GENOMES/human/hg38_GRCh38_Reference_Genome/star_genomes/hg38_len100_fromGTF/'
salmon_genome_dir = '/project/klynclab/GENOMES/human/hg38_GRCh38_Reference_Genome/salmon_genomes/hg38_salmon_transcriptome_index'
salmon_library_type = 'A'  # automatic detection of strandness for Salmon


# Software (in lab software folder)
moduleloadline = 'module load java/openjdk-1.8.0 python/3.6.3 STAR/2.5.2a FastQC-0.11.2'
# includes cutadapt...
trimgalore = '/project/klynclab/software/TrimGalore-0.5.0/trim_galore'
samtools = '/project/klynclab/software/samtools-1.9/samtools'
fastqdump = '/project/klynclab/software/sratoolkit.2.9.2-centos_linux64/bin/fastq-dump'
# Locally=installed software
salmon = '/home/cradens/bin/salmon-0.11.3-linux_x86_64/bin/salmon'


def check_file_exists(filepath):
    try:
        os.path.isfile(filepath)
    except:
        raise RuntimeError("%s doesn't exist..." % filepath)


samples = list()
with open(os.path.join(project_dir, 'sras.txt'), 'r') as openfile:
    for line in openfile:
        sample = line.strip('\r\n')
        if len(sample) > 0:
            samples.append(sample)


# (1) download sras on mercury
# ssh username@mercury.pmacs.upenn.edu
# bash /project/klynclab/software/rna_seq_processing/download_sras.sh
# sras.txt 6 /project/klynclab/ranzani_2015/download_paired.txt


# (2) dump sra to fastq
skip = True
if not skip:
    with open(os.path.join(project_dir, 'submit_sra_to_fastq_paired.sh'), 'w') as sra_to_fastqouthandle:
        for sample in samples:
            sra_to_fastq_fp = os.path.join(
                project_dir, sample, 'sra_to_fastq_%s.sh' % (sample))
            with open(sra_to_fastq_fp, 'w') as outhandle:
                outhandle.write('cd %s\n' %
                                (os.path.join(project_dir, sample)))
                outhandle.write('%s %s.sra --split-3 --gzip\n' %
                                (fastqdump, sample))
                outhandle.write('rm %s.sra\n' % (sample))
            errfp = os.path.join(project_dir, sample,
                                 'error.sra_to_fastq.' + sample)
            outfp = os.path.join(project_dir, sample,
                                 'out.sra_to_fastq.' + sample)
            sra_to_fastqouthandle.write('cd %s\n' % (
                os.path.join(project_dir, sample)))
            sra_to_fastqouthandle.write(
                'bsub -J %s -e %s -o %s bash %s\n' % ('sra_' + sample, errfp, outfp, sra_to_fastq_fp))
# bash /project_dir/submit_sra_to_fastq.sh


# (3) trim and fastqc
skip = True
if not skip:
    with open(os.path.join(project_dir, 'submit_trim_QC_fastq_paired.sh'), 'w') as trim_QC_fastqouthandle:
        for sample in samples:
            trim_QC_fastq_fp = os.path.join(
                project_dir, sample, 'trim_QC_fastq_%s.sh' % (sample))
            with open(trim_QC_fastq_fp, 'w') as outhandle:
                outhandle.write('%s\n' % moduleloadline)
                outhandle.write('cd %s\n' %
                                (os.path.join(project_dir, sample)))
                fq_1 = '%s_1.fastq.gz' % (sample)
                fq_2 = '%s_2.fastq.gz' % (sample)
                check_file_exists(fq_1)
                check_file_exists(fq_2)
                outhandle.write(
                    '%s --paired --stringency 5 --length 35 -q 20 %s %s\n' % (trimgalore, fq_1, fq_2))
                outhandle.write('mkdir ./fastqc_trimmed\n')
                trimmed_fq_1 = '%s_1_val_1.fq.gz' % (sample)
                trimmed_fq_2 = '%s_2_val_2.fq.gz' % (sample)
                outhandle.write(
                    'fastqc %s -o ./fastqc_trimmed &\n' % (trimmed_fq_1))
                outhandle.write(
                    'fastqc %s -o ./fastqc_trimmed &\n' % (trimmed_fq_2))
                outhandle.write('wait\n')
            errfp = os.path.join(project_dir, sample,
                                 'error.trim_QC_fastq.' + sample)
            outfp = os.path.join(project_dir, sample,
                                 'out.trim_QC_fastq.' + sample)
            trim_QC_fastqouthandle.write(
                'cd %s\n' % (os.path.join(project_dir, sample)))
            trim_QC_fastqouthandle.write(
                'bsub -J %s -e %s -o %s bash %s\n' % ('trim_' + sample, errfp, outfp, trim_QC_fastq_fp))
# bash /project_dir/submit_trim_QC_fastq.sh


# (4) align
skip = True
if not skip:
    with open(os.path.join(project_dir, 'submit_star_paired.sh'), 'w') as starouthandle:
        starouthandle.write('mkdir %s\n' % (bamdir))
        for sample in samples:
            star_fp = os.path.join(project_dir, sample,
                                   'star_%s.sh' % (sample))
            with open(star_fp, 'w') as outhandle:
                outhandle.write('%s\n' % moduleloadline)
                outhandle.write('cd %s\n' %
                                (os.path.join(project_dir, sample)))
                trimmed_fq_1 = '%s_1_val_1.fq.gz' % (sample)
                trimmed_fq_2 = '%s_2_val_2.fq.gz' % (sample)
                check_file_exists(trimmed_fq_1)
                check_file_exists(trimmed_fq_2)
                outhandle.write('STAR --genomeDir %s  --readFilesIn %s %s --runThreadN %s --outSAMtype BAM Unsorted --outFileNamePrefix ./%s. --outSAMattributes All --alignSJoverhangMin 8 --readFilesCommand zcat --outSAMunmapped Within\n' %
                                (star_genome, trimmed_fq_1, trimmed_fq_2, threads, sample))
                # sort bam, writng output to bam dirs
                final_bam_fp = os.path.join(bamdir, sample + ".bam")
                outhandle.write('%s sort -@ %s -o %s ./%s.Aligned.out.bam\n' %
                                (samtools, threads, final_bam_fp, sample))
                outhandle.write('%s index -@ %s %s\n' %
                                (samtools, threads, final_bam_fp))
            errfp = os.path.join(project_dir, sample, 'error.star.' + sample)
            outfp = os.path.join(project_dir, sample, 'out.star.' + sample)
            starouthandle.write('cd %s\n' %
                                (os.path.join(project_dir, sample)))
            starouthandle.write('bsub -J %s -e %s -o %s -n %s -R "span[ptile=%s]" -M 50000 bash %s\n' % (
                'star_' + sample, errfp, outfp, threads, threads, star_fp))
# bash /project_dir/submit_star.sh


# (5) salmon
skip = True
if not skip:
    with open(os.path.join(project_dir, 'submit_salmon_paired.sh'), 'w') as salmon_outhandle:
        for sample in samples:
            salmon_fp = os.path.join(
                project_dir, sample, 'salmon_%s.sh' % (sample))
            trimmed_fq_1 = '%s_1_val_1.fq.gz' % (sample)
            trimmed_fq_2 = '%s_2_val_2.fq.gz' % (sample)
            check_file_exists(trimmed_fq_1)
            check_file_exists(trimmed_fq_2)
            with open(salmon_fp, 'w') as outhandle:
                outhandle.write('cd %s\n' %
                                (os.path.join(project_dir, sample)))
                outhandle.write('%s quant -i %s -l %s -p %s -1 %s -2 %s -o salmon_quant_%s\n' % (
                    salmon, salmon_genome_dir, salmon_library_type, threads, trimmed_fq_1, trimmed_fq_2, sample))
            errfp = os.path.join(project_dir, sample, 'error.salmon.' + sample)
            outfp = os.path.join(project_dir, sample, 'out.salmon.' + sample)
            salmon_outhandle.write('cd %s\n' % (
                os.path.join(project_dir, sample)))
            salmon_outhandle.write('bsub -J %s -e %s -o %s -n %s -R "span[ptile=%s]" -M 50000 bash %s\n' % (
                'salmon_' + sample, errfp, outfp, threads, threads, salmon_fp))
    with open(os.path.join(project_dir, 'copy_salmon_lines.txt'), 'w') as copy_outhandle:
        for sample in samples:
            quant_sf_fp = os.path.join(
                project_dir, '%s/salmon_quant_%s/quant.sf' % (sample, sample))
            copy_outhandle.write(
                'scp cradens@mercury.pmacs.upenn.edu:%s ./%s.quant.sf\n' % (quant_sf_fp, sample))

