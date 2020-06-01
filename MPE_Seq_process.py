'''
Created on Sep 16, 2019

@author: kyang
'''

import os
import re
from typing import List
import ast
from glob import glob
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import numpy as np
import pandas as pd
import seaborn as sns

import multiprocessing

sns.set()
import matplotlib
import matplotlib.pyplot as plt
from subprocess import check_output


# Above the divider is the collated pipeline I use for KY002, below is for KY001
# For KY002, we first want some sanity checks. Here are some initial objectives:
# Overall sequence composition - how does the library look overall?
# Do this analysis per library, and then comparing the library prep methods stim/unstim to each other
# %Mapped and deduplicated, final number of reads left
# Percent nucleotide composition of mapped sequences (should this be all mapped or just targeted events?)
# PCA analysis of 1) all genes on target/off target read count, 2) psi values of targeted genes

# On vs off-target
# Extract the targeted genes and LSVs out of the mapped file
# Per LSV, make scatter plots of number of "on target" and "off target" primed reads versus the following

# Per targeted LSV, how many

# initial summary stats
def picard_summarize_read_length(
        liblist="Biotin-RepA-Stim Biotin-RepA-Unstim Biotin-RepC-Stim Biotin-RepD-Stim Biotin-RepD-Unstim TSO-RepA-Stim TSO-RepA-Unstim TSO-RepC-Stim TSO-RepC-Unstim TSO-RepD-Stim TSO-RepD-Unstim"):
    # this makes a script that run picard per library on a few different things:
    # 1)the mapped deduplicated reads
    # 2) the mapped deduplicated reads for just read 1
    # 3) the mapped overall reads before deduplication
    # 4) the unmapped reads (need to use fastqtosam first to convert the unmapped fastq's to a bam)
    pass


def summarize_overall_stats(root_path="/Users/kyang/Box Sync/Rotation_2/MPE-Seq/KY002/STAR_statistics",
                            liblist="Biotin-RepA-Stim Biotin-RepA-Unstim Biotin-RepC-Stim Biotin-RepD-Stim Biotin-RepD-Unstim TSO-RepA-Stim TSO-RepA-Unstim TSO-RepC-Stim TSO-RepC-Unstim TSO-RepD-Stim TSO-RepD-Unstim"):
    # Iterate through files in root_path and 1) pull out relevant data from the STAR final log
    # file into a summary tsv file with number input reads, %uniquely mapped, %unmapped, and %multi-mapping reads
    # ADDITIONALLY parse output of log file for TSO to calculate percent of reads deduplicated
    # root_path is the directory where all files are (STAR and deduplication logs)
    # liblist is a list of space-separated library names
    lib_types = []
    for x in ["Biotin", "TSO"]:
        for y in ["Stim", "Unstim"]:
            lib_types.append((x, y))
    mapped_count_list = []
    total_reads_list = []
    write_out = []
    lib_df = pd.DataFrame(columns=['Library', 'Type', 'Total Reads', 'Mapped %', 'Unmapped %',
                                   '% of Total Mapped Remaining after Deduplication'])
    for lib in liblist.split(" ")[:]:
        os.chdir(root_path)
        mapped_count = 0
        unmapped_count = 0
        # First, open the STAR log
        with open(lib + ".Log.final.out") as STAR_log:
            for line in STAR_log:
                if "Uniquely mapped reads %" in line or "% of reads mapped to multiple loci" in line or "% of reads mapped to too many loci" in line:
                    mapped_count += float(line.split("|\t")[-1][:-2])
                if "% of reads unmapped" in line:
                    unmapped_count += float(line.split("|\t")[-1][:-2])
                if "Number of input reads" in line:
                    total_reads = int(line.split("|\t")[-1][:-1])
        # Next, open the deduplicated counts file
        temp_denom = 0
        temp_num = 0
        with open(lib + "_UMI_deduplicated.log") as dedup_log:
            for line in dedup_log:
                if "INFO Reads: Input Reads:" in line:
                    temp_denom = int(line.split(" ")[6][:-1])
                if "INFO Number of reads out:" in line:
                    temp_num = int(line.split(" ")[-1][:-1])
        temp_frac = temp_num / temp_denom * 100
        write_out.append("\t".join(
            [lib, str(total_reads / 1000000), str(round(mapped_count, 3)), str(round(unmapped_count, 3)),
             str(round(temp_frac, 3)) + "\n"]))
        for type in lib_types:
            if type[0] in lib and type[1] in lib:
                break
        lib_df = lib_df.append(
            {'Library': lib, 'Type': type[0] + ",\n" + type[1], 'Total Reads': total_reads, 'Mapped %': mapped_count,
             'Unmapped %': unmapped_count, '% of Total Mapped Remaining after Deduplication': temp_frac},
            ignore_index=True)
    # For a few of the columns in lib_df, need to make boxplot/stripplots for these
    fig, axs = plt.subplots(1, 4)
    sns.boxplot(x='Type', y='Total Reads', ax=axs[0], data=lib_df)
    sns.stripplot(x='Type', y='Total Reads', ax=axs[0], data=lib_df, color=".3")
    sns.boxplot(x='Type', y='Mapped %', ax=axs[1], data=lib_df)
    sns.stripplot(x='Type', y='Mapped %', ax=axs[1], data=lib_df, color=".3")
    axs[1].set_ylim((0, 100))
    sns.boxplot(x='Type', y='Unmapped %', ax=axs[2], data=lib_df)
    sns.stripplot(x='Type', y='Unmapped %', ax=axs[2], data=lib_df, color=".3")
    axs[2].set_ylim((0, 100))
    sns.boxplot(x='Type', y='% of Total Mapped Remaining after Deduplication', ax=axs[3], data=lib_df)
    sns.stripplot(x='Type', y='% of Total Mapped Remaining after Deduplication', ax=axs[3], data=lib_df, color=".3")
    axs[3].set_ylim((0, 100))
    # plt.tight_layout()
    plt.show()
    os.chdir(root_path)
    with open("summary_stats", 'w+') as inSum:
        inSum.write("\t".join(
            ['Library', 'Total Reads', 'Mapped %', 'Unmapped %', '% of Total Mapped Remaining after Deduplication\n']))
        inSum.write("".join(write_out))


def gen_lsv_bed_1(gtf_path="/Users/kyang/Box Sync/Rotation_2/ref/gencode.v31.annotation.gtf.bed",
                  all_chr_bed_path="/Users/kyang/Box Sync/Rotation_2/ref/chr_lengths.bed.txt",
                  tsv_path="/Users/kyang/Box Sync/Rotation_2/Primers/KY002/",
                  out_path="/Users/kyang/Box Sync/Rotation_2/Primers/KY002/KY002_on_target_",
                  out_short="/Users/kyang/Box Sync/Rotation_2/Primers/KY002/"):
    # given the folder to the .tsvs of KY002, return separately 2 files:
    # 1) a list of the gene ENSG ID's and 2) a bed file of targeted junctions
    # currently this method is limited to target exons specifically
    # gtf_path = path to the folder containing the gencode gtf file converted to bed via bedops (see Analysis/11_29_19.ipynb)
    # tsv_path = path to the folder of all tsvs from previous classify_process method in OligoMiner_wrapper
    # out_path = the prefix for generating the 2 files, one of which will be a simple text file (list of ENSGs) and
    # the other of which will be a properly formatted .bed file

    # Additionally, it would be helpful to have a third file which outputs coordinates to determine
    # reads which are extended at least from 1 base pair upstream of the primer
    # up to 1 base pair downstream of the exon junction coordinate

    # The last useful measures would be coordinates to compute psi values:
    # Namely this would be coordinates of 1 splicing event over the other... think about this more
    # it will be a fraction of anything that overlaps with 1 basepair before the junction
    # bed.write("chr\tlower_bound\tupper_bound\n")
    gtf_ref = open(gtf_path).readlines()
    os.chdir(tsv_path)
    file_list = [x for x in os.listdir() if ".tsv" in x]
    gtf_gene_ref = []
    gtf_exon_ref = []
    with open(gtf_path) as inG:
        for line in inG:
            if line.split("\t")[7] == "gene":
                gtf_gene_ref.append(line)
            if line.split("\t")[7] == "exon":
                gtf_exon_ref.append(line)
    print(str(len(gtf_gene_ref)))
    with open(out_short + "all_genes.bed", 'w+') as inG:
        inG.write("\n".join(gtf_gene_ref))
    with open(out_short + "all_exons.bed", 'w+') as inG:
        inG.write("\n".join(gtf_exon_ref))
    if True:
        count = 0
        for f in file_list:
            bed = open(out_path + f.split(".")[0] + "_LSVs.bed", 'w+')
            genes = open(out_path + f.split(".")[0] + "_genes.bed", 'w+')
            extended = open(out_path + f.split(".")[0] + "_extended.bed", 'w+')
            bedlist = []
            extendlist = []
            genelist = []
            with open(f) as inF:
                for line in inF:
                    if "ENSG" in line:
                        subline_list: List[str] = line.split("\t")[1].split(":")
                        if "constitutive" in line:
                            subline_list[2] = "".join([subline_list[2].split("_")[0] + subline_list[2][-3:]])
                        # First, find the corresponding gene so we can extract info about chromosome
                        # and also write out the bed file for the gene
                        gene = subline_list[0].split(".")[0]
                        linelist = []
                        for line2 in gtf_ref:
                            if gene in line2 and line2.split("\t")[7] == "gene":
                                linelist = line2.split("\t")
                                chrom = linelist[0]
                                genelist.append("\t".join(linelist[:3] + [gene]))
                                if line2 in gtf_gene_ref:
                                    count += 1
                                    gtf_gene_ref.remove(line2)
                                break
                        else:
                            print("Error")
                        # if positive stranded, junction corresponds to upstream exon boundary
                        if "+" in subline_list[-1]:
                            bedlist.append(chrom + "\t" + subline_list[2].split("-")[0] + "\t" + str(
                                int(subline_list[2].split("-")[0]) + 1) + "\t" + line.split("\t")[1])
                            extendlist.append(chrom + "\t" + subline_list[2].split("-")[0] + "\t" + str(
                                int(subline_list[2].split("-")[0]) + int(line.split("\t")[2])) + "\t" +
                                              line.split("\t")[1])
                        # if negative stranded, junction corresponds to downstream exon boundary
                        else:
                            bedlist.append(chrom + "\t" + str(int(subline_list[2].split("-")[1][:-1]) - 1) + "\t" +
                                           subline_list[2].split("-")[1][:-1] + "\t" + line.split("\t")[1])
                            extendlist.append(chrom + "\t" + str(
                                int(subline_list[2].split("-")[1][:-1]) - int(line.split("\t")[2])) + "\t" +
                                              subline_list[2].split("-")[1][:-1] + "\t" + line.split("\t")[1])
            bed.write("\n".join(bedlist))
            genes.write("\n".join(genelist))
            extended.write("\n".join(extendlist))
            bed.close()
            genes.close()
            extended.close()
        with open(out_short + "off_target_genes.bed", 'w+') as off_target:
            gtf_gene_ref_mod = ["\t".join(x.split("\t")[:4]) for x in gtf_gene_ref]
            print(len(gtf_gene_ref_mod))
            off_target.write("\n".join(
                gtf_gene_ref_mod))  # NEW: We want a couple other beds as well for extracting off-target reads  # Let
            # 's make two beds: one is a bed of all possible off-target non-genic regions (, another is a bed of all
            # possible  # off-target exonic regions, the last is bed of all possible off-target intronic regions


def intersect_script(out_dir="/Users/kyang/Box Sync/Rotation_2/MPE-Seq/KY002/",
                     scratch_dir=" /scr1/users/yangk4/KY002/", n_threads=3,
                     bed_string="KY002_on_target_hi_confidence_large_change_primers_extended.bed KY002_on_target_hi_confidence_large_change_primers_genes.bed KY002_on_target_hi_confidence_large_change_primers_LSVs.bed KY002_on_target_hi_confidence_small_change_primers_extended.bed KY002_on_target_hi_confidence_small_change_primers_genes.bed KY002_on_target_hi_confidence_small_change_primers_LSVs.bed KY002_on_target_low_confidence_large_change_primers_extended.bed KY002_on_target_low_confidence_large_change_primers_genes.bed KY002_on_target_low_confidence_large_change_primers_LSVs.bed KY002_on_target_low_confidence_small_change_primers_extended.bed KY002_on_target_low_confidence_small_change_primers_genes.bed KY002_on_target_low_confidence_small_change_primers_LSVs.bed",
                     off_target_bed="off_target_genes.bed",
                     folder_string="Biotin-RepA-Stim/ Biotin-RepA-Unstim/ Biotin-RepC-Stim/ Biotin-RepD-Stim/ Biotin-RepD-Unstim/ TSO-RepA-Stim/ TSO-RepA-Unstim/ TSO-RepC-Stim/ TSO-RepC-Unstim/ TSO-RepD-Stim/ TSO-RepD-Unstim/"):
    # scratch_dir is the directory containing .bed in the top level directory and subdirectories of bams from MPE_Seq_pipeline
    # NOTE: you CANNOT use symlinks for scratch dir otherwise qsub has an Eqw error I think because of permissions
    # use the hard coded path
    # bed_string is the output from echo *.bed in bash EXCLUDING the off-target bed
    # off_target_bed is the name of the off-target bed
    # folder_string is the output from echo */ in bash
    # n_threads is number of processes to allow for running
    # Now that appropriate beds have been generated, generate a script which can run bedtools across all these beds in parallel
    fw = open(out_dir + "bed_process", 'w+')
    fw.write('''#! /bin/bash
#$ -V
#$ -wd %s
#$ -l h_vmem=%sG
#$ -l m_mem_free=%sG
#$ -pe smp %s -binding linear:%s
module load gcc8
module load BEDTools
cd %s
''' % (scratch_dir, str(30), str(30), str(n_threads + 1), str(n_threads + 1), scratch_dir))
    srr_list = folder_string.split(" ")
    for folder in srr_list:
        templist = []
        folder = folder.split("/")[0]
        with open(out_dir + 'intersect_' + folder, 'w+') as inF:
            for bed in bed_string.split(" "):
                templist.append("bedtools intersect -b %s/%s_deduplicated.bam -a %s -c -bed >> %s/%s" % (
                folder, folder, bed, folder, "intersect_" + bed))
            inF.write("\n".join(templist))
    chunk_list = [srr_list[i:i + n_threads] for i in range(0, len(srr_list), n_threads)]
    for n in range(len(chunk_list)):
        fw.write('#Processing chunk %s\n' % n)
        chunk = chunk_list[n]  # this is not right
        for samp in chunk:
            fw.write('bash intersect_%s &\n' % samp.split("/")[0])
        fw.write('wait\n')
    fw.close()
    fw2 = open(out_dir + "bed_process_off_target", 'w+')
    fw2.write('''#! /bin/bash
#$ -V
#$ -wd %s
#$ -l h_vmem=%sG
#$ -l m_mem_free=%sG
#$ -pe smp %s -binding linear:%s
module load gcc8
module load BEDTools
cd %s
''' % (scratch_dir, str(30), str(30), str(n_threads + 1), str(n_threads + 1), scratch_dir))
    srr_list = folder_string.split(" ")
    chunk_list = [srr_list[i:i + n_threads] for i in range(0, len(srr_list), n_threads)]
    for n in range(len(chunk_list)):
        fw2.write('#Processing chunk %s\n' % n)
        chunk = chunk_list[n]  # this is not right
        for samp in chunk:
            folder = samp.split("/")[0]
            tempstring = "bedtools intersect -b %s/%s_deduplicated.bam -a %s -c -bed >> %s/%s" % (
                folder, folder, off_target_bed, folder, "intersect_" + off_target_bed)
            fw2.write('%s &\n' % tempstring)
        fw2.write('wait\n')
    fw2.close()


def summarize_intersect_pca(intersect_path="/Users/kyang/Box Sync/Rotation_2/MPE-Seq/KY002/intersect_beds",
                            read_counts_path="/Users/kyang/Box Sync/Rotation_2/MPE-Seq/KY002/library_read_counts",
                            kallisto_path="/Users/kyang/Box Sync/Rotation_2/RNA-Seq/JSL1/",
                            cdna_fa_path="/Users/kyang/Box Sync/Rotation_2/ref/Homo_sapiens.GRCh38.cdna.all.fa"):
    # intersect_path is the path to a folder containing subdirectories, each of which
    # contains the bed files for a single library prep
    # read_counts_path is the path to the number of read counts in every library (see block comment below)
    '''
    #to get readcounts run the following bash script:
    dirs=$(echo */)
    for i in $dirs
    do
        l1=$(samtools idxstats $i*_deduplicated.bam | awk '{total += $3} END{print total}')
        echo -e ${i::-1}"\t"$l1 >> library_read_counts
    done
    '''
    # Kallisto_path is the path to a folder containing Stim and Unstim average TPM output from Kallisto.
    # gtf_path is the path to the gencode gtf converted to bed (could have implemented with actual gtf but easier since I
    # already have the bed from previous processing step)
    # after the above intersect script has been run, output some summary statistics
    # and basic graphs summarizing number of unextended, extended, and per library
    # consider combining TSO libraries into 1 average and Pleiss libraries into a 2nd average
    # if they look broadly similar

    # To determine if they look broadly similar, do PCA for the libraries for number of reads corresponding to
    # each gene out of total number of reads
    # See https://towardsdatascience.com/pca-using-python-scikit-learn-e653f8989e60
    # Also do for number of reads corresponding to LSVs and number of reads corresponding to unextended primers.
    read_count_dict = dict()
    with open(read_counts_path) as inR:
        for line in inR:
            linelist = line.split("\t")
            read_count_dict[linelist[0]] = int(linelist[1])
    # X is a M*N-dimensional pandas matrix df where M is the number of libraries and N is the different genes
    # with Y being the libraries
    # Parse the data to get this desired matrix (first testing for just gene expression to make sure it works)
    os.chdir(intersect_path)
    intersect_dirlist = [x for x in os.listdir() if os.path.isdir(x)]  # list directories to iterate here
    intersect_dirlist.sort(key=lambda x: (x.split("-")[0], x.split("-")[2], x.split("-")[1]))
    # Do each of genes, extended, LSV
    X_all = []  # X_all is a list of arrays with each array being either  genes, extended, LSVs, or off-target
    for k in ["genes", "extended", "LSVs"]:
        lib_names = []
        c2 = 0
        for f in intersect_dirlist:
            c2 += 1
            lib_names.append(f.split("/")[0])  # lib names is the list tracking names of things
            os.chdir(intersect_path)
            os.chdir(f)
            # Need to concatenate on-target-beds into same row
            on_target_beds = [g for g in os.listdir()[:] if "on_target" in g and k in g]
            c1 = 0
            for g in on_target_beds:
                c1 += 1
                temp = pd.read_table(g, header=None)
                headers = temp[temp.columns[3:4]].T
                temp = pd.DataFrame(temp[temp.columns[4:5]].T.values, columns=headers.iloc[0].values)
                temp = temp.loc[:, ~temp.columns.duplicated()]  # remove dup columns
                if c1 == 1:
                    final1 = temp
                else:
                    # this concatenates rows while removing dupes
                    final1 = pd.concat([final1, temp], axis=1)
                    final1 = final1.loc[:,
                             ~final1.columns.duplicated()]  # After one row is finished, move on to the next one
            final1 = final1 / read_count_dict[
                f.split("/")[0]]  # divide by total read counts to normalize between libraries
            if c2 == 1:
                final_on = final1
            else:
                final_on = final_on.append(final1)  # divide by read counts to normalize between libraries
        if k == "genes":
            on_target_pd = final_on
        if k == "extended":
            extended_pd = final_on
        X_all.append(final_on.rename_axis('ID').values)

    # Now, just 1 off-target-bed (but do same thing for it)
    c2 = 0
    for f in intersect_dirlist:
        c2 += 1
        os.chdir(intersect_path)
        os.chdir(f)
        off_target_bed = [g for g in os.listdir()[:] if "off_target" in g]
        for g in off_target_bed:
            temp = pd.read_table(g, header=None)
            temp = temp.loc[:, ~temp.columns.duplicated()]
            headers = temp[temp.columns[3:4]].T
            temp = pd.DataFrame(temp[temp.columns[4:5]].T.values, columns=headers.iloc[0].values)
            final1 = temp
            final1 = final1 / read_count_dict[f.split("/")[
                0]]  # add a pseudocount of 1 and divide by total read counts to normalize between libraries
        if c2 == 1:
            final_off = final1
        else:
            final_off = final_off.append(final1)
    off_target_pd = final_off
    off_target_pd.columns = [x.split(".")[0] for x in off_target_pd.columns]
    # print(str(off_target_pd))
    X_all.append(final_off.rename_axis('ID').values)

    # Code here is a bit messy because I had to add a new plot for nongenic regions

    library_prep = []
    library_prep2 = []
    for l in lib_names:
        if "TSO" in l:
            l1 = "TSO"
        else:
            l1 = "Biotin"
        if "Unstim" in l:
            l2 = "Unstim"
        else:
            l2 = "Stim"
        library_prep.append(l1 + ", " + l2)
        library_prep2.append(l1 + ",\n" + l2)
    if False:
        # this does PCA
        for index, X in enumerate(X_all):
            # Standardscaler to standardize
            X = StandardScaler().fit_transform(X)
            pca = PCA(n_components=2)
            principalComponents = pca.fit_transform(X)
            print(str(principalComponents))
            principalDf = pd.DataFrame(data=principalComponents,
                                       columns=['principal component 1', 'principal component 2'])
            principalDf.insert(0, "library_type", library_prep, True)
            # just plotting stuff
            # We color different stim and unstim conditions differently
            # Label all points with what they are called
            fig, axe = plt.subplots()
            d = sns.scatterplot(x="principal component 1", y="principal component 2", hue="library_type",
                                data=principalDf, ax=axe)
            # add labels, see https://stackoverflow.com/questions/14432557/matplotlib-scatter-plot-with-different-text-at-each-data-point
            for i, txt in enumerate(lib_names):
                d.annotate(txt, (np.asarray(principalDf["principal component 1"])[i],
                                 np.asarray(principalDf["principal component 2"])[i]))

            axe.set(xlabel="principal component 1 (% Variance Explained = " + str(
                np.round(pca.explained_variance_ratio_[0], 4)) + ")",
                    ylabel="principal component 2 (% Variance Explained = " + str(
                        np.round(pca.explained_variance_ratio_[1], 4)) + ")")
            os.chdir(intersect_path)
            plt.savefig(["genes_PCA", "extended_PCA", "LSVs_PCA", "off_target_PCA"][index] + ".png", dpi=200)
            plt.show()
    # We also want to make a simple boxplot showing % of on target and off target reads overall for each of
    # 4 conditions
    # sum across all columns of extended_pd and off_target_pd
    off_target_series = off_target_pd.sum(axis=1).reset_index(drop=True)
    on_target_series = extended_pd.sum(axis=1).reset_index(drop=True)
    off_target_percent = (off_target_series / (off_target_series + on_target_series) * 100).rename("Off Target %")
    on_target_percent = (on_target_series / (off_target_series + on_target_series) * 100).rename("On Target %")
    combined_pd = on_target_percent.to_frame().join(off_target_percent)
    combined_pd.insert(0, "Library Type", library_prep2, True)
    fig, axs = plt.subplots(1, 2)
    sns.boxplot(x='Library Type', y='On Target %', ax=axs[0], data=combined_pd)
    sns.stripplot(x='Library Type', y='On Target %', ax=axs[0], data=combined_pd, color=".3")
    axs[0].set_ylim((0, None))
    sns.boxplot(x='Library Type', y='Off Target %', ax=axs[1], data=combined_pd)
    sns.stripplot(x='Library Type', y='Off Target %', ax=axs[1], data=combined_pd, color=".3")
    axs[1].set_ylim((None, 100))
    os.chdir(intersect_path)
    plt.tight_layout()
    plt.savefig("On_off_target_percent_summary.png", dpi=200)
    plt.show()

    # We can try to correlate the number on-target reads with Kallisto gene expression and do the same thing with off-target reads
    # Depending on how this looks, we can also normalize by Kallisto gene expression to calculate amount enriched
    if False:
        plt.rcParams.update({'font.size': 7})
        plt.rc('axes', titlesize=10)  # fontsize of the axes title
        plt.rc('axes', labelsize=7)  # fontsize of the x and y labels
        plt.rc('xtick', labelsize=7)  # fontsize of the tick labels
        plt.rc('ytick', labelsize=7)  # fontsize of the tick labels

        GTF_lines = dict()
        # with open(gtf_path) as inG:
        #    for line in inG:
        #        if line[7] == "transcript":
        #            GTF_lines.append(line)
        with open(cdna_fa_path) as inG:
            for line in inG:
                if line[0] == ">":
                    linelist = line.split(" ")
                    GTF_lines[linelist[0][1:]] = linelist[3].split(":")[1].split(".")[0]
        # for stim and unstim Kallisto files,
        # we want to build pandas series containing the sum of TPM for genes we are interested in
        # We want ENSGs as columns and a single row corresponding to the sum
        unstim_tpm = pd.Series()
        with open(kallisto_path + "Unstim_TPM.txt") as inK:
            for i, line in enumerate(inK):
                if i == 0:
                    continue
                linelist = line.split("\t")
                gene = GTF_lines[linelist[0]]
                if gene in unstim_tpm.index:
                    unstim_tpm[gene] += float(linelist[1])
                else:
                    unstim_tpm[gene] = float(linelist[1])
        stim_tpm = pd.Series()
        with open(kallisto_path + "Stim_TPM.txt") as inK:
            for i, line in enumerate(inK):
                if i == 0:
                    continue
                linelist = line.split("\t")
                # for GTF_line in GTF_lines:
                #     if linelist[0] in GTF_line:
                #         gene = GTF_line.split("\t")[3]
                gene = GTF_lines[linelist[0]]
                if gene in stim_tpm.index:
                    stim_tpm[gene] += float(linelist[1])
                else:
                    stim_tpm[gene] = float(linelist[
                                               1])  # else:  #     if float(linelist[1]) > 0:  #         raise ValueError('A ENST ID was not found in the GTF AKA '+linelist[0])
        # iterate over each row of the gene matrix, making 3*3+2*1 subplots
        # we want to create a new plot for each row
        c = 0
        for df in [on_target_pd, off_target_pd]:
            for tpm_vals in [unstim_tpm, stim_tpm]:
                fig, axs = plt.subplots(3, 4, sharey=True)
                fig.suptitle(["On Target vs Bulk Unstim", "On Target vs Bulk Stim", "Off Target vs Bulk Unstim",
                              "Off Target vs Bulk Stim"][c], fontsize='large')
                c += 1
                for (i0, (i2, row)) in enumerate(df.iterrows()):
                    i = i0
                    if i0 >= 4:
                        i += 1
                    # return[tpm_vals,row]
                    cmap = sns.color_palette(["Blues", "Greens", "Reds", "Oranges"][int(i / 3)])
                    # sns.set_palette(sns.color_palette(["Blues","Greens","Reds","Oranges"][int(i / 3)]))
                    # make a new dataframe which has the tpm values as row 1
                    # and the row of the original genes np array as row 2
                    # this dataframe should ONLY have values included in
                    # then plot as a scatterplot in ax.subplot
                    temp_frame = pd.DataFrame(zip(tpm_vals[row.index], row),
                                              columns=["Bulk RNASeq TPM", "MPE-Seq Read Count"])
                    # calculate r-value and add to upper left of graph
                    corr_val = round(temp_frame["Bulk RNASeq TPM"].corr(temp_frame["MPE-Seq Read Count"]), 4)
                    axs[i % 3, int(i / 3)].set_xlim((0, row.max() * 1.01))
                    sns.scatterplot(x="MPE-Seq Read Count", y="Bulk RNASeq TPM", ax=axs[i % 3, int(i / 3)],
                                    data=temp_frame, palette=cmap)
                    axs[i % 3, int(i / 3)].set_title(lib_names[i0])
                    axs[i % 3, int(i / 3)].text(0.1, 0.9, "R= " + str(corr_val), ha='center', va='center',
                                                transform=axs[i % 3, int(
                                                    i / 3)].transAxes)  # return temp_frame  # sns.scatterplot(x="MPE-Seq Read Count", y="Bulk RNASeq TPM",data=temp_frame)
                plt.tight_layout()
                plt.show()  # return [tpm_vals, off_target_pd]  # return temp_frame


def calc_psi():
    # PSI calculations
    # -Calculate psi for all junctions as a simple measure independent of MAJIQ and see how well it agrees
    # Per LSV, we want the intersect of a given splice junction immediately be
    pass

    # Can do the same thing except with fold enrichment vs Kallisto-estimated RNASeq-expression


def blast_primed_regions(db=False, primerpath="/Users/kyang/Box Sync/Rotation_2/Primers/KY002/",
                         rootpath="/Users/kyang/Box Sync/Rotation_2/MPE-Seq/KY002/", prefix="/scr1/users/yangk4/KY002/",
                         # folder_string="Biotin-RepA-Stim",
                         folder_string="Biotin-RepA-Stim Biotin-RepA-Unstim Biotin-RepC-Stim Biotin-RepD-Stim Biotin-RepD-Unstim TSO-RepA-Stim TSO-RepA-Unstim TSO-RepC-Stim TSO-RepC-Unstim TSO-RepD-Stim TSO-RepD-Unstim",
                         nthreads="5", group_size=2):
    # Given an input file of primers, generate a bash file for HPC which will for each library folder,
    # a) generate a blast db for each bam file
    # 1) blast each primer against the entire bam file and extract "hit" reads with
    # samtools view | grep
    # 2) divide each primer .fa file with reads into 3 different files:
    # a) on target mRNA mapped, b) off target mRNA, mapped
    # using samtools view chr :
    # 3) run analyses for each
    # db = False if bam not generated yet into blastdb, otherwise the path provided will be the path to db
    # primerpath = FASTA file of primers to BLASTN-short for
    # rootpath = file on LOCAL computer where script is generated
    # prefix = path prefix on HOST computer where script is run, it should contain directories containing
    # the fastq files for each library
    # folder_string = the list of subdirectories to cd into from rootpath

    # First, generate a primer FASTA
    os.chdir(primerpath)
    tsvs = [f for f in os.listdir()[:] if ".tsv" in f]
    primer_list = []
    for f in tsvs:
        with open(f) as inF:
            for line in inF:
                if "ENSG" in line:
                    primer_list.append(">" + line.split("\t")[1] + "\n" + line.split("\t")[4])
    with open(primerpath + "primers.fa", 'w+') as inP:
        inP.write("\n".join(primer_list))
    with open(rootpath + "extract_primer", 'w+') as inBlast:
        if db == False:
            inBlast.write('''#! /bin/bash
#$ -V
#$ -wd %s
#$ -l h_vmem=8G
#$ -l m_mem_free=8G
#$ -pe smp %s -binding linear:%s
module load SAMtools
folders="%s"
cd %s
for c in $folders;
do
    echo "cd $c/" >> blast_$c
    #convert bam to fastq (extract just read 1), need this for blast database compilation
    echo "samtools fasta -1 ${c}_deduplicated_read1.fa -n ${c}_deduplicated.bam >/dev/null" >> blast_$c
    #extract uniquely mapping Read 1 to new sam file
    echo "samtools view -f 64 -F 256 -q 255 ${c}_deduplicated.bam > ${c}_deduplicated.sam" >> blast_$c
    #Run makeblastdb for the bam file
    echo "mkdir blast" >> blast_$c
    echo "makeblastdb -in ${c}_deduplicated_read1.fa -dbtype nucl -hash_index -out blast/$c -max_file_sz 4GB -logfile blast/index.log" >> blast_$c
    #BLAST primer file against generated blastdbs
    echo "blastn -query %s -task blastn-short -db blast/$c -num_threads %s -max_target_seqs 10000000 -outfmt \\\"6 qseqid sseqid btop pident length mismatch gapopen qstart qend sstart send evalue\\\" -out blast/primer_query.tsv" >> blast_$c
done
''' % (prefix, str(group_size * int(nthreads)), str(group_size * int(nthreads)), folder_string, prefix,  # intial
       prefix + "primers.fa", nthreads  # blast params
       ))
        srr_list = folder_string.split(" ")
        chunk_list = [srr_list[i:i + group_size] for i in range(0, len(srr_list), group_size)]
        for n in range(len(chunk_list)):
            inBlast.write('#Processing chunk %s\n' % n)
            arr_string = ''
            chunk = chunk_list[n]  # this is not right
            for samp in chunk:
                inBlast.write('bash blast_%s &\n' % samp)
                arr_string += '"%s" ' % samp
            inBlast.write('wait\n')


def extract_primer_origin(prefix="/scr1/users/yangk4/KY002/",
                          fold_string="Biotin-RepA-Stim Biotin-RepA-Unstim Biotin-RepC-Stim Biotin-RepD-Stim Biotin-RepD-Unstim TSO-RepA-Stim TSO-RepA-Unstim TSO-RepC-Stim TSO-RepC-Unstim TSO-RepD-Stim TSO-RepD-Unstim", ):
    # given that primers have been blasted against the deduplicated fastq, write out a new file with the following
    # columns per read: 1) the name of primer, 2) the name of the read, 2) the sequence of the read,
    # 3) the aligned genomic coordinates from the original bam, 4) positive vs negative strand, 5) ON vs
    # OFF target gene-based alignment (if not on, infer it is off)
    # number 3 and 4 follow bed format
    # (can tell from the binary-decoded FLAG) '{0:08b}'.format(num)
    # 5) the cigar representation from original bam, 6) the btop representation from the new blast
    # Also, we want to report how many total unique mapped reads we are able to account for from the entire
    # deduplicated bam file via this BLAST approach
    os.chdir(prefix)
    l = os.listdir()
    on_target = dict()
    for file in l:
        if "on_target" in file and "extended.bed" in file:
            with open(file) as inF:
                for line in inF:
                    linelist = line.split("\t")
                    if "\n" == linelist[3][-1]:
                        linelist[3] = linelist[3][:-1]
                    on_target[linelist[3]] = {'chr': linelist[0], 'start': linelist[1], 'end': linelist[2]}
    primer_dict = dict()
    with open("primers.fa") as inFasta:
        odd = -1
        for line in inFasta:
            if "\n" == line[-1]:
                line = line[:-1]
            odd *= -1
            if odd == 1:
                temp = line[1:]
            else:
                primer_dict[temp] = str(len(line))
    for fold in fold_string.split(" "):
        os.chdir(prefix)
        os.chdir(fold)
        # sam_pd = pd.DataFrame(columns=['flag','chr','start','cigar','Read Sequence'])
        sam_pd = dict()
        fast_dict = dict()
        with open(fold + "_deduplicated_read1.fa") as inFasta:
            odd = -1
            for line in inFasta:
                odd *= -1
                if odd == 1:
                    temp = line[1:-1]
                else:
                    fast_dict[temp] = line[:-1]
        inSam = open(fold + "_deduplicated.sam")
        k = 0
        print("Opening sam for " + fold)
        for line in inSam.readlines():
            k += 1
            linelist = line.split("\t")
            # if sam_pd.loc[linelist[0]] in fast_dict[linelist[0]]}
            if k % 500000 == 0:
                print("Processed " + str(k) + " lines of sam")
            sam_pd[linelist[0]] = {'flag': linelist[1], 'chr': linelist[2], 'start': linelist[3], 'end': linelist[8],
                                   'cigar': linelist[5], 'Read Sequence': fast_dict[linelist[0]]}
        inSam.close()
        write_pd = list()
        inP = open("blast/primer_query.tsv")
        k = 0
        read_set = set()
        # print(on_target.to_string)
        reject = 0
        reject2 = 0
        for line in inP.readlines():
            k += 1
            read_set.add(linelist[1])
            if k % 500000 == 0:
                print("Processed " + str(k) + " lines of primer_query")
            linelist = line.split("\t")
            if int(linelist[-3]) > 40:  # reject an alignment if it aligns to more than 40 bp out
                reject += 1
                continue
            if linelist[
                1] in sam_pd.keys():  # skip over if it is not a valid, uniquely mapping read and is a multi mapper
                read_series = sam_pd[linelist[1]]
            else:
                continue
            if int(primer_dict[linelist[0]]) >= len(
                    sam_pd[linelist[1]]['Read Sequence']):  # reject a read if it is just unextended primer
                reject2 += 1
                continue
            on_target_series = on_target[linelist[0]]
            strandedness = "+" if int('{0:08b}'.format(int(read_series['flag']))[-5]) == 0 else "-"
            # defined by coverage
            read_bound_1 = min(int(read_series['start']), int(read_series['start']) + int(read_series['end']))
            read_bound_2 = max(int(read_series['start']), int(read_series['start']) + int(read_series['end']))
            if linelist[0][-2] == "+":
                # is_on_target = True if (read_series['chr'] == on_target_series['chr'] and
                #                        int(on_target_series['start']) - 50 - int(primer_dict[linelist[0]]) <= int(read_series['start']) <= int(
                #            on_target_series['end'])) else False
                is_on_target = True if read_series['chr'] == on_target_series['chr'] and (
                    check_overlap(read_bound_1, read_bound_2,
                                  int(on_target_series['start']) - 50 - int(primer_dict[linelist[0]]),
                                  int(on_target_series['end']))) else False
            else:
                # is_on_target = True if (read_series['chr'] == on_target_series['chr'] and
                #                        int(on_target_series['end']) + 50 + int(primer_dict[linelist[0]]) >= int(read_series['start']) >= int(
                #            on_target_series['start'])) else False
                is_on_target = True if read_series['chr'] == on_target_series['chr'] and (
                    check_overlap(read_bound_1, read_bound_2, int(on_target_series['start']),
                                  int(on_target_series['end']) + 50 + int(primer_dict[linelist[0]]))) else False
            # the "end" is
            # if read_series['chr'] == on_target_series['chr']:
            #    print(str(on_target_series))
            #    print(str(read_series))
            #    print(str(is_on_target))
            # write_pd.append({'Primer':linelist[0],'Primer Length':primer_dict[linelist[0]],'Read':linelist[1],
            #                 'Read Sequence':sam_pd[linelist[1]]['Read Sequence'],'Read Mapped Location',
            #                 'On/Off Target':is_on_target,'Read Strandedness':strandedness,
            #                 'Read CIGAR':read_series['cigar'],'Primer BTOP':linelist[2]})
            write_pd.append("\t".join([str(linelist[0]), str(primer_dict[linelist[0]]), str(linelist[1]),
                                       str(sam_pd[linelist[1]]['Read Sequence']),
                                       str(read_series['chr']) + ":" + str(read_bound_1) + "-" + str(read_bound_2),
                                       str(is_on_target), str(strandedness), str(read_series['cigar']),
                                       str(linelist[2])]))
        with open("primer_summary.tsv", 'w+') as inP:
            inP.write("\n".join(write_pd))
        print("Rejected " + str(reject) + " Alignments")
        print("Rejected " + str(reject2) + " Unextended Primers")
        print(str(len(read_set)) + " Unique Reads")


def check_overlap(x, y, x_star, y_star):
    # helper function for extract_primer_origin, this simply checks that a read overlaps with a given interval
    # X and y are the bounds of the read to be checked and
    # x_star and y_star are the bounds of the interval to check for overlap with
    if (x <= x_star + 2 and y >= x_star - 2) or (y >= y_star - 2 and x <= y_star + 2) or (
            x_star - 2 <= x <= y_star + 2) or (y_star + 2 >= y >= x_star - 2):
        return True
    else:
        return False


def corr_primers(prefix="/scr1/users/yangk4/KY002/", # fold_string="Biotin-RepA-Stim Biotin-RepC-Stim"
                 fold_string="Biotin-RepA-Stim Biotin-RepC-Stim Biotin-RepD-Stim Biotin-RepA-Unstim Biotin-RepD-Unstim TSO-RepA-Stim TSO-RepC-Stim TSO-RepD-Stim TSO-RepA-Unstim TSO-RepC-Unstim TSO-RepD-Unstim"
                 # fold_string="Biotin-RepA-Stim Biotin-RepA-Unstim Biotin-RepC-Stim Biotin-RepD-Stim Biotin-RepD-Unstim TSO-RepA-Stim TSO-RepA-Unstim TSO-RepC-Stim TSO-RepC-Unstim TSO-RepD-Stim TSO-RepD-Unstim",
                 ):
    # given the directory of tab-delimited files containing reads organized by sequence, we want a few basic measures
    # number of reads accountable from primers
    # 1) the number of ON target and OFF target alignments per primer, plotted as a stacked bar graph
    # 2) correlation of ON target vs off target alignments per primer
    # 3) the number of matching nucleotides per primer by distribution
    # 4) the number of mismatched nucleotides per primer
    # 5) the number of consecutive non-mismatched nucleotides from the 3' end per primer
    # 1a) Number of ON target and OFF target alignments per primer, side by side bar graphs in ascending order
    # 1b) Number of on target vs number of off target alignments per primer, side by side correlation plots
    os.chdir(prefix)
    dirlist = os.listdir()
    p_duplex = dict()
    tm = dict()
    gc = dict()
    tss_dist = dict()
    consec_df_list, consec_on_df_list, consec_off_df_list, target_df_list = [], [], [], []

    for f in os.listdir():
        if "_primers.tsv" in f:
            with open(f) as inf:
                for i, line in enumerate(inf):
                    if "ENSG" not in line:
                        continue
                    if line[-1] == "\n":
                        line = line[:-1]
                    linelist = line.split("\t")
                    if linelist[7] == "unique":
                        p_duplex[linelist[1]] = float(0)
                    else:
                        p_duplex[linelist[1]] = float(linelist[7])
                    tm[linelist[1]] = float(linelist[6])
                    g_or_c = [x for x in linelist[4] if x == 'C' or x == 'G']
                    gc[linelist[1]] = float(len(g_or_c) / len(linelist[4]))
                    tss_dist[linelist[1]] = float(linelist[2])

    for fold in fold_string.split(" "):
        os.chdir(prefix)
        os.chdir(fold)
        inP = open("primer_summary.tsv")
        current = "blah"
        target_pd = pd.DataFrame(
            columns=['On Target Read Count', 'Off Target Read Count', 'Average Consecutive Nucleotides', 'P(duplex)',
                     'Tm'])
        consec_list = dict()
        consec_off_list = dict()
        consec_on_list = dict()
        consec_nucleotides_all = list()
        on_target_count = 0
        off_target_count = 0
        k = 0
        for line in inP.readlines():
            linelist = line.split("\t")
            if linelist[0] != current:
                if current != "blah":
                    # print(str([x for x in consec_list[current] if isinstance(x, str)]))
                    consec_list[current] = {'Primer Length': primer_length,
                                            'Average Consecutive Nucleotides': sum(consec_list[current]) / len(
                                                consec_list[current]), 'Tm': tm[current], 'GC Content': gc[current],
                                            'Distance from Junction': tss_dist[current]}
                    if len(consec_on_list[current]) > 0:
                        consec_on_list[current] = {'Primer Length': primer_length,
                                                   'Average Consecutive Nucleotides': sum(
                                                       consec_on_list[current]) / len(consec_on_list[current]),
                                                   'Tm': tm[current], 'GC Content': gc[current],
                                                   'Distance from Junction': tss_dist[current]}
                    else:
                        temp = consec_on_list.pop(current)
                    consec_off_list[current] = {'Primer Length': primer_length,
                                                'Average Consecutive Nucleotides': sum(consec_off_list[current]) / len(
                                                    consec_off_list[current]), 'Tm': tm[current],
                                                'GC Content': gc[current], 'Distance from Junction': tss_dist[current]}
                    target_pd.loc[current] = {'On Target Read Count': on_target_count,
                                              'Off Target Read Count': off_target_count,
                                              'Average Consecutive Nucleotides': consec_list[current][
                                                  'Average Consecutive Nucleotides'], 'P(duplex)': p_duplex[current],
                                              'Tm': tm[current]}

                on_target_count = 0
                off_target_count = 0
                consec_list[linelist[0]] = list()
                consec_on_list[linelist[0]] = list()
                consec_off_list[linelist[0]] = list()
                primer_length = float(linelist[1])
                current = linelist[0]
            if linelist[5] == "True":
                on_target_count += 1
            else:
                off_target_count += 1
            cigar = linelist[7]
            btoplist = re.findall('\d*', linelist[8][:-1])[:-1]
            btop = 0
            for b in btoplist:
                if b == '':
                    btop += 0.5
                else:
                    btop += int(b)
            if linelist[6] == "+":  # positive stranded, expect beginning match
                c = cigar
                for char in c:
                    if char.isalpha():
                        if char == 'M':  # very beginning does match without soft-clip so just return btop
                            consec = btop
                        elif char == 'S' or char == 'D' or char == 'I':
                            if char == 'D' or char == 'I':
                                print(cigar)
                            consec = max(0, btop - int(cigar.split(char)[0]))
                            # if btop - int(cigar.split(char)[0]) < 0: #formula is btop minus softclipped to get
                            #    print(str([btop,cigar]))
                            break  # part of primer which actually matches
                        else:
                            raise ValueError('Char ' + char + " was observed")
                        break
            else:
                for char in reversed(cigar):
                    if char.isalpha():
                        if char == 'M':
                            consec = btop
                        elif char == 'S' or char == 'D' or char == 'I':
                            if char == 'D' or char == 'I':
                                print(cigar)
                            consec = max(0, btop - int([x for x in re.findall('\d*', cigar) if x != ''][-1]))
                        else:
                            raise ValueError('Char ' + char + " was observed")
                        break
            # consec_dict[linelist[2]] = {'Primer':linelist[0],'Primer Length':linelist[1],'Consecutive Nucleotides':int(char)}
            # Append to both a primer-specific list and a global list
            consec_list[linelist[0]].append(int(consec))
            consec_nucleotides_all.append(int(consec))
            # Append to a on-target specific list and a off-target specific list depending
            if linelist[5] == "True":
                consec_on_list[linelist[0]].append(int(consec))
            else:
                consec_off_list[linelist[0]].append(int(consec))

            # above lines writes out number of consecutive nucleotides in primer at 3' end
        consec_list[current] = {'Primer Length': primer_length,
                                'Average Consecutive Nucleotides': sum(consec_list[current]) / len(
                                    consec_list[current]), 'Tm': tm[current], 'GC Content': gc[current],
                                'Distance from Junction': tss_dist[current]}
        consec_on_list[current] = {'Primer Length': primer_length,
                                   'Average Consecutive Nucleotides': sum(consec_on_list[current]) / len(
                                       consec_on_list[current]), 'Tm': tm[current], 'GC Content': gc[current],
                                   'Distance from Junction': tss_dist[current]}
        consec_off_list[current] = {'Primer Length': primer_length,
                                    'Average Consecutive Nucleotides': sum(consec_off_list[current]) / len(
                                        consec_off_list[current]), 'Tm': tm[current], 'GC Content': gc[current],
                                    'Distance from Junction': tss_dist[current]}
        target_pd.loc[current] = {'On Target Read Count': on_target_count, 'Off Target Read Count': off_target_count,
                                  'Average Consecutive Nucleotides': consec_list[current][
                                      'Average Consecutive Nucleotides'], 'P(duplex)': p_duplex[current],
                                  'Tm': tm[current]}
        # print("Number of soft-clipped reads that couldn't be resolved: "+str(k))
        # print(str(consec_list))
        consec_df = pd.DataFrame(consec_list).T
        consec_off_df = pd.DataFrame(consec_off_list).T
        consec_on_df = pd.DataFrame(consec_on_list).T
        consec_df_list.append(consec_df)
        consec_off_df_list.append(consec_off_df)
        consec_on_df_list.append(consec_on_df)
        target_df_list.append(target_pd)
    # Now plot everything
    plt.rcParams.update({'font.size': 7})
    plt.rc('axes', titlesize=8)  # fontsize of the axes title
    plt.rc('axes', labelsize=7)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=7)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=7)  # fontsize of the tick labels
    fig, axe = plt.subplots(3, 4)
    fig2, axe2 = plt.subplots(3, 4)
    fig3, axe3 = plt.subplots(3, 4)
    fig4, axe4 = plt.subplots(3, 4)
    fig5, axe5 = plt.subplots(3, 4)
    fig6, axe6 = plt.subplots(3, 4)
    fig7, axe7 = plt.subplots(3, 4)
    fig8, axe8 = plt.subplots(3, 4)
    fig9, axe9 = plt.subplots(3, 4)
    fig10, axe10 = plt.subplots(3, 4)
    fig11, axe11 = plt.subplots(3, 4)
    fig12, axe12 = plt.subplots(3, 4)
    fig13, axe13 = plt.subplots(3, 4)
    fig14, axe14 = plt.subplots(3, 4)
    fig15, axe15 = plt.subplots(3, 4)
    # fig16, axe16 = plt.subplots(3, 4)
    # fig17, axe17 = plt.subplots(3, 4)
    fig01, ax1 = plt.subplots(3, 4)
    fig02, ax2 = plt.subplots(3, 4)
    fig03, ax3 = plt.subplots(3, 4)
    fig04, ax4 = plt.subplots(3, 4)
    fig05, ax5 = plt.subplots(3, 4)
    fig06, ax6 = plt.subplots(3, 4)
    fig07, ax7 = plt.subplots(3, 4)
    fig08, ax8 = plt.subplots(3, 4)
    for (i_temp, (consec_df, consec_off_df, consec_on_df, target_pd)) in enumerate(
            zip(consec_df_list, consec_off_df_list, consec_on_df_list, target_df_list)):
        if i_temp >= 4:
            i = i_temp + 1
        else:
            i = i_temp
        # Consecutive nucleotides without distinguishing on vs off
        sns.distplot(consec_df['Average Consecutive Nucleotides'], kde=False, ax=axe[i % 3, int(i / 3)])

        # consec nucleotides vs off, not distinguishing on vs off
        sns.scatterplot(x='Average Consecutive Nucleotides', y='Off Target Read Count', data=target_pd,
                        ax=axe2[i % 3, int(i / 3)])
        corr_val = round(target_pd['Average Consecutive Nucleotides'].corr(target_pd['Off Target Read Count']), 4)
        axe2[i % 3, int(i / 3)].text(0.1, 0.9, "R= " + str(corr_val), ha='center', va='center',
                                     transform=axe2[i % 3, int(i / 3)].transAxes)

        # consec nucleotides vs on, not distinguishing on vs off
        sns.scatterplot(x='Average Consecutive Nucleotides', y='On Target Read Count', data=target_pd,
                        ax=axe3[i % 3, int(i / 3)])
        corr_val = round(target_pd['Average Consecutive Nucleotides'].corr(target_pd['On Target Read Count']), 4)
        axe3[i % 3, int(i / 3)].text(0.1, 0.9, "R= " + str(corr_val), ha='center', va='center',
                                     transform=axe3[i % 3, int(i / 3)].transAxes)

        # Consecutive nucleotides for on target only
        sns.distplot(consec_on_df['Average Consecutive Nucleotides'], kde=False, ax=axe4[i % 3, int(i / 3)])

        # consec nucleotides vs on for on target only
        sns.scatterplot(x=consec_on_df['Average Consecutive Nucleotides'], y=target_pd['On Target Read Count'],
                        ax=axe5[i % 3, int(i / 3)])
        corr_val = round(consec_on_df['Average Consecutive Nucleotides'].corr(target_pd['On Target Read Count']), 4)
        axe5[i % 3, int(i / 3)].text(0.1, 0.9, "R= " + str(corr_val), ha='center', va='center',
                                     transform=axe4[i % 3, int(i / 3)].transAxes)
        axe5[i % 3, int(i / 3)].set(xlabel='Average Consecutive Nucleotides')

        # Consecutive nucleotides for off target only
        sns.distplot(consec_off_df['Average Consecutive Nucleotides'], kde=False, ax=axe6[i % 3, int(i / 3)])

        # consec nucleotides vs on for off target only
        sns.scatterplot(x=consec_off_df['Average Consecutive Nucleotides'], y=target_pd['Off Target Read Count'],
                        ax=axe7[i % 3, int(i / 3)])
        corr_val = round(consec_off_df['Average Consecutive Nucleotides'].corr(target_pd['Off Target Read Count']), 4)
        axe7[i % 3, int(i / 3)].text(0.1, 0.9, "R= " + str(corr_val), ha='center', va='center',
                                     transform=axe7[i % 3, int(i / 3)].transAxes)

        # consec nucleotides normalized by primer length
        sns.scatterplot(x=consec_on_df['Average Consecutive Nucleotides'] / consec_on_df['Primer Length'],
                        y=target_pd['On Target Read Count'], ax=axe8[i % 3, int(i / 3)])
        axe8[i % 3, int(i / 3)].set(xlabel='Normalized Avg Consec Nucleotides')
        temp = consec_on_df['Average Consecutive Nucleotides'] / consec_on_df['Primer Length']
        corr_val = round(temp.corr(target_pd['On Target Read Count']), 4)
        axe8[i % 3, int(i / 3)].text(0.1, 0.9, "R= " + str(corr_val), ha='center', va='center',
                                     transform=axe8[i % 3, int(i / 3)].transAxes)

        # consec nucleotides normalized by primer length

        sns.scatterplot(x=consec_off_df['Average Consecutive Nucleotides'] / consec_off_df['Primer Length'],
                        y=target_pd['Off Target Read Count'], ax=axe9[i % 3, int(i / 3)])
        axe9[i % 3, int(i / 3)].set(xlabel='Normalized Avg Consec Nucleotides')
        temp = consec_off_df['Average Consecutive Nucleotides'] / consec_off_df['Primer Length']
        corr_val = round(temp.corr(target_pd['Off Target Read Count']), 4)
        axe9[i % 3, int(i / 3)].text(0.1, 0.9, "R= " + str(corr_val), ha='center', va='center',
                                     transform=axe9[i % 3, int(i / 3)].transAxes)

        # primer GC content vs off and on target

        sns.scatterplot(x=consec_off_df['GC Content'], y=target_pd['Off Target Read Count'],
                        ax=axe10[i % 3, int(i / 3)])
        corr_val = round(consec_off_df['GC Content'].corr(target_pd['Off Target Read Count']), 4)
        axe10[i % 3, int(i / 3)].text(0.1, 0.9, "R= " + str(corr_val), ha='center', va='center',
                                      transform=axe10[i % 3, int(i / 3)].transAxes)

        sns.scatterplot(x=consec_on_df['GC Content'], y=target_pd['On Target Read Count'], ax=axe11[i % 3, int(i / 3)])
        corr_val = round(consec_on_df['GC Content'].corr(target_pd['On Target Read Count']), 4)
        axe11[i % 3, int(i / 3)].text(0.1, 0.9, "R= " + str(corr_val), ha='center', va='center',
                                      transform=axe11[i % 3, int(i / 3)].transAxes)

        # primer Tm vs off and on target

        sns.scatterplot(x=consec_off_df['Tm'], y=target_pd['Off Target Read Count'], ax=axe12[i % 3, int(i / 3)])
        corr_val = round(consec_off_df['Tm'].corr(target_pd['Off Target Read Count']), 4)
        axe12[i % 3, int(i / 3)].text(0.1, 0.9, "R= " + str(corr_val), ha='center', va='center',
                                      transform=axe12[i % 3, int(i / 3)].transAxes)

        sns.scatterplot(x=consec_on_df['Tm'], y=target_pd['On Target Read Count'], ax=axe13[i % 3, int(i / 3)])
        corr_val = round(consec_on_df['Tm'].corr(target_pd['On Target Read Count']), 4)
        axe13[i % 3, int(i / 3)].text(0.1, 0.9, "R= " + str(corr_val), ha='center', va='center',
                                      transform=axe13[i % 3, int(i / 3)].transAxes)

        # Distance from junction vs off and on target

        sns.scatterplot(x=consec_off_df['Distance from Junction'], y=target_pd['Off Target Read Count'],
                        ax=axe14[i % 3, int(i / 3)])
        corr_val = round(consec_off_df['Distance from Junction'].corr(target_pd['Off Target Read Count']), 4)
        axe14[i % 3, int(i / 3)].text(0.1, 0.9, "R= " + str(corr_val), ha='center', va='center',
                                      transform=axe14[i % 3, int(i / 3)].transAxes)

        sns.scatterplot(x=consec_on_df['Distance from Junction'], y=target_pd['On Target Read Count'],
                        ax=axe15[i % 3, int(i / 3)])
        corr_val = round(consec_on_df['Distance from Junction'].corr(target_pd['On Target Read Count']), 4)
        axe15[i % 3, int(i / 3)].text(0.1, 0.9, "R= " + str(corr_val), ha='center', va='center',
                                      transform=axe15[i % 3, int(i / 3)].transAxes)

        # most common motifs mono- di- and tri- nucleotides per primer in top 25th percentile

        # P(duplex) vs off and on-target, filter unique
        # off target
        sns.scatterplot(x='P(duplex)', y='Off Target Read Count', data=target_pd[target_pd['P(duplex)'] > 0],
                        ax=ax5[i % 3, int(i / 3)])  #
        ax5[i % 3, int(i / 3)].set_ylim((0, None))
        ax5[i % 3, int(i / 3)].set_xlim((0, None))
        corr_val = round(target_pd[target_pd['P(duplex)'] > 0]['P(duplex)'].corr(
            target_pd[target_pd['P(duplex)'] > 0]['On Target Read Count']), 4)
        ax5[i % 3, int(i / 3)].text(0.1, 0.9, "R= " + str(corr_val), ha='center', va='center',
                                    transform=ax5[i % 3, int(i / 3)].transAxes)
        # on target

        # print(str(target_pd))
        sns.scatterplot(x='P(duplex)', y='On Target Read Count', data=target_pd[target_pd['P(duplex)'] > 0],
                        ax=ax6[i % 3, int(i / 3)])
        ax6[i % 3, int(i / 3)].set_ylim((0, None))
        ax6[i % 3, int(i / 3)].set_xlim((0, None))
        corr_val = round(target_pd[target_pd['P(duplex)'] > 0]['P(duplex)'].corr(
            target_pd[target_pd['P(duplex)'] > 0]['On Target Read Count']), 4)
        ax6[i % 3, int(i / 3)].text(0.1, 0.9, "R= " + str(corr_val), ha='center', va='center',
                                    transform=ax6[i % 3, int(i / 3)].transAxes)

        # Distribution of unique vs all other P(duplex), off target read count
        target_pd = target_pd.assign(
            is_unique=lambda dataframe: dataframe['P(duplex)'].map(lambda p: "unique" if p == 0 else "non-unique"))
        sns.boxplot(x='is_unique', y='Off Target Read Count', ax=ax7[i % 3, int(i / 3)], data=target_pd)
        # medians = target_pd.groupby(['is_unique'])['Off Target Read Count'].median().values()
        # median_labels = [str(np.round(s,2)) for s in medians]
        # pos = range(len(medians))
        # for tick, label in zip(pos,ax7[i % 3, int(i / 3)].get_xticklabels()):
        #    ax7[i % 3, int(i / 3)].text(pos[tick],medians[tick] + 0.5, median_labels[tick],
        #    horizontalalignment='center', size='x-small', color='w', weight='semibold')

        # Distribution of unique vs all other P(duplex), on target read count
        # target_pd = target_pd.assign(
        #    is_unique=lambda dataframe: dataframe['P(duplex)'].map(lambda p: "unique" if p == 0 else "non-unique"))
        sns.boxplot(x='is_unique', y='On Target Read Count', ax=ax8[i % 3, int(i / 3)], data=target_pd)
        # medians = target_pd.groupby(['is_unique'])['On Target Read Count'].median().values()
        # median_labels = [str(np.round(s, 2)) for s in medians]
        # pos = range(len(medians))
        # for tick, label in zip(pos, ax8[i % 3, int(i / 3)].get_xticklabels()):
        #    ax8[i % 3, int(i / 3)].text(pos[tick], medians[tick] + 0.5, median_labels[tick],
        #                                horizontalalignment='center', size='x-small', color='w', weight='semibold')
        # sns.stripplot(x='is_unique', y='Off Target Read Count', ax=ax7[i % 3, int(i / 3)], data=target_pd, color=".3")

        # test correlation of on target vs off target: does correlation exist?

        # print(str(target_pd))
        sns.scatterplot(x='On Target Read Count', y='Off Target Read Count', data=target_pd, ax=ax1[i % 3, int(i / 3)])
        ax1[i % 3, int(i / 3)].set_ylim((0, None))
        ax1[i % 3, int(i / 3)].set_xlim((0, None))
        corr_val = round(target_pd['Average Consecutive Nucleotides'].corr(target_pd['Off Target Read Count']), 4)
        ax1[i % 3, int(i / 3)].text(0.1, 0.9, "R= " + str(corr_val), ha='center', va='center',
                                    transform=ax1[i % 3, int(i / 3)].transAxes)

        sns.distplot(target_pd['On Target Read Count'], kde=False, ax=ax2[i % 3, int(i / 3)])
        sns.distplot(target_pd['Off Target Read Count'], kde=False, ax=ax3[i % 3, int(i / 3)])
        sns.distplot(target_pd['On Target Read Count'] / (target_pd['Off Target Read Count'] + 1), kde=False,
                     ax=ax4[i % 3, int(i / 3)])
        ax4[i % 3, int(i / 3)].set(xlabel='On Target/Off Target Ratio')
        # ax4[i % 3, int(i / 3)].set_ylim((0, 3))
        ax4[i % 3, int(i / 3)].set_xlim((0, None))

        # print(str((target_pd['On Target Read Count']/(target_pd['Off Target Read Count']+1)).max()))
    # name titles
    for f in [axe, axe2, axe3, axe4, axe5, axe6, axe7, axe8, axe9, axe10, axe11, axe12, axe13, axe14, axe15, ax6, ax7,
              ax8, ax1, ax2, ax3, ax4, ax5]:
        for i, sub_title in enumerate(fold_string.split(" ")):
            f[i % 3, int(i / 3)].set_title(sub_title, fontsize=8)

    for i, f in enumerate(
            [fig, fig2, fig3, fig4, fig5, fig6, fig7, fig8, fig9, fig10, fig11, fig12, fig13, fig14, fig15, fig05,
             fig06, fig07, fig08, fig01, fig02, fig03, fig04]):
        f.set_size_inches(11, 6)
        f.tight_layout()
        f.savefig(prefix + "Fig_combined_" + str(i) + ".png", dpi=200)
        f.show()


def extract_suitable_off_targets(off_target_exon_fold="/Users/kyang/Box Sync/Rotation_2/MPE-Seq/KY002/off_target_exons",
                                 libsizes="/Users/kyang/Box "
                                          "Sync/Rotation_2/MPE-Seq/KY002/off_target_exons/lib_sizes"):
    # after confirming what non-genic regions are and that most things are genic, we will look at genic regions
    # and try to figure out suitable genic regions for primers to be designed to
    # Edit 1/6/20: it looks like only 1-2% map non-genic (need to make plot proving this), and most of the "non-genic"
    # contaminant is really just
    # unannotated genic stuff.
    # off_target_exon_fold is the folder containing downloaded off target exon bed count files output from bedtools
    # intersect.
    # libsizes is the path to the file containing library sizes
    # Looking at genic regions, what we want is an exon-based analysis of off-target genes so we need to
    # run intersect script after
    # first plotting a histogram of each exon normalized based on exon length
    os.chdir(off_target_exon_fold)
    filelist = [x for x in os.listdir() if "bed" in x.split(".") and "off_target_exons" in x]
    c = 0
    sizedict = dict()
    with open(libsizes) as inF:
        for i, line in enumerate(inF):
            if i % 2 == 1:
                sizedict[temp] = int(line[:-1])
            temp = line[:-1]
    exonDict = dict()
    for j, file in enumerate(filelist):
        print("On " + file)
        size = sizedict[file.split("_")[0]]
        with open(file) as inF:
            F_read = inF.readlines()
        for i, line in enumerate(F_read):
            if i % 100000 == 0:
                print("on line " + str(i))
            linelist = line.split("\t")
            exonID = (linelist[4], int(linelist[5]), int(linelist[6]))
            # even though some exon IDs are redundant doesn't matter if we just join dictionaries at the end
            # since the redundant exon IDs will all have the same count anyways
            if j == 0:
                exonDict[exonID] = [int(linelist[-1][:-1]) / size]
            elif len(exonDict[exonID]) < j + 1:
                exonDict[exonID].append(int(linelist[-1][
                                            :-1]) / size)  # if exonID not in exonDict.keys():  #    templist.append(exonID)  #    exonDict[exonID] = [int(linelist[-1][:-1])/size]  # elif exonID not in templist:  #    templist.append(exonID)  #    exonDict[exonID].append(int(linelist[-1][:-1])/size)
    # pool = multiprocessing.Pool(processes=nthreads)
    # resulting_dicts = [pool.apply()]
    # so far, this adds the count normalized by library size per library
    # what we want to output is the mean count and the variance of count (i.e. measured by SD)
    # we will write out these mean and variance counts to a new file (after being normalized by exon size)
    # and also create a histogram of mean size-normalized counts
    avglist = []
    cvlist = []
    outlist = []
    with open("exon_summary", "w+") as exonSum:
        outlist.append("exon\taverage\tstd\tlength\tlength normalized average\tlength normalized std\tCV")
        for k, v in exonDict.items():
            if len(v) != len(filelist):
                c += 1
                continue
            average = sum(v) / len(v)
            std = np.std(v)
            length = k[2] - k[1]
            norm_average = average / length
            norm_std = std / length
            cv = norm_std / norm_average * 100
            avglist.append(norm_average)
            cvlist.append(norm_std)
            outlist.append("\t".join(
                [k[0] + ":" + str(k[1]) + "-" + str(k[2]), str(average), str(std), str(length), str(norm_average),
                 str(norm_std), str(cv)]))
        exonSum.write("\n".join(outlist))


def plot_suitable_off_targets_stats(off_target_sum_file="/Users/kyang/Box "
                                                        "Sync/Rotation_2/MPE-Seq/KY002/off_target_exons/exon_summary",
                                    prefix="/Users/kyang/Box Sync/Rotation_2/MPE-Seq/KY002/off_target_exons/",
                                    libsizes="/Users/kyang/Box "
                                              "Sync/Rotation_2/MPE-Seq/KY002/off_target_exons/lib_sizes"):
    # based on results from the above method, we will try to figure out some rules for selection of the exons to do
    # RT-PCR against. I'll try to select ~5-10 top candidate exons.
    # off_target_sum_file is the path to the file output from the preivous function
    #prefix is the folder to which files should be saved, and in which the non genic beds from bedtools
    # intersect are located

    # Filter based on some percentile cutoff first for 1) low coefficient of variation in expression, and 2) high
    # expression we can try for instance >= 95th percentile in both categories only (hopefully there will be some
    # left still) then, after obtaining a smaller filter list, rank by homogeneity in terms of proportion of the most
    # represented primer per event, we will hopefully obtain events that are reproducible (low variance), detectable (
    # high expression), and result from a single primer (homogeneous)
    df = pd.read_csv(off_target_sum_file, sep="\t")
    df.dropna()
    sizedict = dict()
    with open(libsizes) as inF:
        for i, line in enumerate(inF):
            if i % 2 == 1:
                sizedict[temp] = int(line[:-1])
            temp = line[:-1]
    os.chdir(prefix)
    filelist = [x for x in os.listdir() if "bed" in x.split(".") and "non_genic_regions" in x]
    non_genic_count = dict()
    for f in filelist:
        tempdf = pd.read_csv(f,sep="\t",header=None,usecols=[3])
        non_genic_count[f.split("_")[0]] = tempdf.sum(axis=0)[3]
    non_genic_count_percents = list()
    for k in sizedict.keys():
        non_genic_count_percents.append((sizedict[k]-non_genic_count[k])/sizedict[k])
    fig0,ax0 = plt.subplots()
    fig01, ax1 = plt.subplots(1, 3)
    fig01_2, ax1_2 = plt.subplots(1, 3)
    # fig01_3, ax1_3 = plt.subplots()
    fig02, ax2 = plt.subplots()
    fig04, ax4 = plt.subplots()
    fig03, ax3 = plt.subplots()
    fig05, ax5 = plt.subplots()
    #Plot a boxplot showing distribution of number of genic vs non-genic reads across all libraries
    #Just show genic % since non-genic is 100% minus everything else
    percent_df = pd.DataFrame(non_genic_count_percents,columns=["genic percentage"])
    sns.boxplot(y="genic percentage",data=percent_df,ax=ax0)

    # Plot a boxplot for the normalized means, coefficients of variance, lengths
    sns.boxplot(y="length normalized average", data=df, ax=ax1[0])
    sns.boxplot(y="length", data=df, ax=ax1[1])
    sns.boxplot(y="CV", data=df, ax=ax1[2])

    # Same plot but with outliers removed
    sns.boxplot(y="length normalized average", data=df, ax=ax1_2[0], showfliers=False)
    sns.boxplot(y="length", data=df, ax=ax1_2[1], showfliers=False)
    sns.boxplot(y="CV", data=df, ax=ax1_2[2], showfliers=False)

    ##one more plot to look at general distribution -- plot a kde of the average normalized count with the CDF
    ## on top of it
    # sns.kdeplot(df["average"], ax=ax1_3, cumulative=True, shade=True)
    ##plot SD versus mean with R coefficient to prove that looking at the coefficient of variance is a bit better
    ##than just looking at the raw SD here. CV also removes the effects of length as well.

    sns.scatterplot(x="std", y="average", data=df, ax=ax2)
    # also add R corr
    corr_val = round(df["std"].corr(df["average"]), 4)
    ax2.text(0.1, 0.9, "R= " + str(corr_val), ha='center', va='center', transform=ax2.transAxes)
    ax2.set_ylim((0, 0.01))
    ax2.set_xlim((0, 0.005))

    # Now do same thing but CV vs average instead
    sns.scatterplot(x="CV", y="average", data=df, ax=ax4)
    corr_val = round(df["CV"].corr(df["average"]), 4)
    ax4.text(0.1, 0.9, "R= " + str(corr_val), ha='center', va='center', transform=ax4.transAxes)
    ax4.set_ylim((0, 0.01))
    ax4.set_xlim((0, None))

    # plot mean versus length to see if it is even necessary to normalize by length -- appears not
    sns.scatterplot(x="average", y="length", data=df, ax=ax3)
    ax3.set_xlim((0, 0.01))
    ax3.set_ylim((0, None))
    # also add R corr
    corr_val = round(df["average"].corr(df["length"]), 4)
    ax3.text(0.1, 0.9, "R= " + str(corr_val), ha='center', va='center', transform=ax3.transAxes)

    #repeat with filtered set
    df_temp = df[df["CV"] <= df["CV"].quantile(0.05)]
    sns.scatterplot(x="average", y="length", data=df_temp, ax=ax5)
    ax5.set_xlim((0, 0.01))
    ax5.set_ylim((0, None))
    # also add R corr
    corr_val = round(df_temp["average"].corr(df_temp["length"]), 4)
    ax5.text(0.1, 0.9, "R= " + str(corr_val), ha='center', va='center', transform=ax5.transAxes)

    for i, f in enumerate([fig0, fig01, fig01_2, fig02, fig04, fig03, fig05]):
        f.set_size_inches(11, 6)
        f.tight_layout()
        f.savefig(prefix + "Fig_exon_off_target_" + str(i) + ".png", dpi=200)  # f.show()

    # Now, we want to create some rules for filtering

    df_filtered = df[(df["length normalized average"] >= df["length normalized average"].quantile(0.99)) & (
                df["CV"] <= df["CV"].quantile(0.05)) & (df['length'] >= 100)]
    print(str(df_filtered))
    df_filtered.to_csv(prefix + "filtered_exons_out.tsv", sep="\t", index=False)


def calc_suitable_off_targets_homogeneity(filtered_file="/home/yangk4/scratch/KY002/off_target_exons/filtered_exons_out.tsv",
                                          rootdir="/home/yangk4/scratch/KY002/",
                                          fold_string="Biotin-RepA-Stim Biotin-RepC-Stim Biotin-RepD-Stim "
                                                      "Biotin-RepA-Unstim Biotin-RepD-Unstim TSO-RepA-Stim "
                                                      "TSO-RepC-Stim TSO-RepD-Stim TSO-RepA-Unstim TSO-RepC-Unstim "
                                                      "TSO-RepD-Unstim"):
    # This takes the filtered output from the previous function
    # it adds a new column to calculate how homogeneous each exon is in terms of primer coverage
    # then orders the exons from greatest to least homogeneity
    # filtered_file is the path to the filtered output from the previous function
    # rootdir is the main directory all the other folders are in
    # fold_string is the path to the folders to iterate over
    # Here, we choose to ignore which library each read came from and instead look at the overall percent across all
    # libraries collectively
    dfList = []
    # process tsvs for all libraries
    for fold in fold_string.split(" "):
        print("processing primer summary tsv for " + fold)
        os.chdir(rootdir)
        os.chdir(fold)
        temp_frame = pd.read_csv("primer_summary.tsv", sep="\t", header=None, usecols=[0, 4])
        temp_list = [x for x in temp_frame[4].str.replace(":", "-").str.split("-").tolist() if len(x) == 3]  # parse
        # chr, start, and stop, there are some empty strings in some for some reason
        #update: empty strings because of mapping to strange things like circular chromosome (mitochondrial) or chimeras
        #not very many like 10-20 per library however so we can just discard them
        new_frame = pd.DataFrame(temp_list, columns=['chr', 'start', 'stop'])
        new_frame['start'] = new_frame['start'].astype(int)
        new_frame['stop'] = new_frame['stop'].astype(int)
        new_frame.insert(0, "primer", temp_frame[0], True)
        dfList.append(new_frame)
    # now make a single giant dataframe out of dfList
    df_final = dfList[0].append([dfList[x] for x in range(1, len(dfList))])
    # process the filtered file also into dataframe
    exon_frame = pd.read_csv(filtered_file, sep="\t")
    # add a few new columns to eventually write out in exon frame
    # exon_frame["primer percentages"] = np.nan
    # exon_frame["primer names"] = np.nan
    percent_list = []
    name_list = []
    # now we want to filter the primer summary tsv for overlap per exon in exon_frame
    for exon in exon_frame['exon']:
        exonlist = exon.replace(":", "-").split("-")
        exonchr = exonlist[0]
        exonstart = int(exonlist[1])
        exonstop = int(exonlist[2])
        # now filter the df_final for reads overlapping exons
        # see https://stackoverflow.com/questions/6821156/how-to-find-range-overlap-in-python
        summarySeries = df_final[(((df_final['start'] < exonstop) & (df_final['stop'] > exonstart)) | ((
                    df_final['stop'] > exonstart) & (exonstop > df_final['start']))) & (
                                             df_final['chr'] == exonchr)]['primer'].value_counts(normalize=True)
        percent_list.append(summarySeries.tolist())
        name_list.append(summarySeries.index.tolist())
    exon_frame["primer percentages"] = percent_list
    exon_frame["primer names"] = name_list
    # make a temp column to sort
    exon_frame['max'] = exon_frame["primer percentages"].apply(max)
    exon_frame = exon_frame.sort_values('max', ascending=False).drop('max', 1)
    # now write out the results to a new csv
    print(str(exon_frame))
    exon_frame.to_csv(rootdir + "off_target_exons/filtered_exons_out_with_homogeneity_calc.tsv", sep="\t", index=False)

def count_on_target_primers(rootdir="/home/yangk4/scratch/KY002/",
                                          fold_string="Biotin-RepA-Stim Biotin-RepC-Stim Biotin-RepD-Stim "
                                                      "Biotin-RepA-Unstim Biotin-RepD-Unstim TSO-RepA-Stim "
                                                      "TSO-RepC-Stim TSO-RepD-Stim TSO-RepA-Unstim TSO-RepC-Unstim "
                                                      "TSO-RepD-Unstim"):
    #while corr_primer is used to make graphs, it doesn't output anything that gives a summary of output stats
    #so this is a method which counts on target primers and outputs a list of library depth-normalized counts,
    #average, std, and coefficient of variation for number of on target reads
    #This could also be easily adapted to count off target as well although it doesn't right now

    #use the primers.fa in the folder previously used to blast primers against the bams as a list of all primers to
    # iterate through
    primerlist = []
    os.chdir(rootdir)
    with open("primers.fa") as inP:
        for line in inP:
            if ">" in line:
                primerlist.append(line[1:-1])
    sizedict = dict()
    with open("lib_sizes") as inF:
        for i, line in enumerate(inF):
            if i % 2 == 1:
                sizedict[temp] = int(line[:-1])
            temp = line[:-1]
    on_target_counts = []
    lib_lengths = []
    for fold in fold_string.split(" "):
        os.chdir(rootdir)
        os.chdir(fold)
        temp_frame = pd.read_csv("primer_summary.tsv", sep="\t", header=None, usecols=[0, 5])
        temp_frame[5] = temp_frame[5].astype(bool)
        temp_frame = temp_frame[temp_frame[5] == True]
        on_target_counts.append(temp_frame[0].value_counts())
        lib_lengths.append(sizedict[fold])
    outlist = []
    outlist.append("Primer\tCounts per Library\tNormalized Average\tNormalized Std\tCV")
    for primer in primerlist:
        temp_list = []
        for i in range(len(on_target_counts)):
            if primer in on_target_counts[i].index:
                temp_list.append(on_target_counts[i][primer])
            else:
                temp_list.append(0)
        avg = np.average(temp_list)/lib_lengths[i]
        std = np.std(temp_list)/lib_lengths[i]
        CV = std/avg*100
        outlist.append("\t".join([primer,str(temp_list),str(avg),str(std),str(CV)]))
    os.chdir(rootdir)
    with open("primer_on_target_count_summary.tsv",'w+') as inP:
        inP.write("\n".join(outlist))
def combine_primer_tsvs(filt_file = "/Users/kyang/Box "
                                    "Sync/Rotation_2/MPE-Seq/KY002/off_target_exons/filtered_exons_out_with_homogeneity_calc.tsv",
                        primer_file = "/Users/kyang/Box "
                                      "Sync/Rotation_2/MPE-Seq/KY002/off_target_exons/primer_on_target_count_summary"
                                      ".tsv",
                        write_dir = "/Users/kyang/Box "
                                      "Sync/Rotation_2/MPE-Seq/KY002/off_target_exons/"):
    #just a helper function to help myself interpret results
    #it just adds
    primer_df = pd.read_csv(primer_file, sep="\t", usecols=["Primer","Normalized Average"])
    primer_df['pct_rank'] = primer_df["Normalized Average"].rank(pct = True)
    filt_df = pd.read_csv(filt_file, sep="\t")
    filt_df['primer names'] = filt_df['primer names'].apply(ast.literal_eval)
    pct_list = []
    norm_list = []
    for k in filt_df['primer names']:
        temp = primer_df[primer_df['Primer'] == k[0]]
        pct_list.append(temp['pct_rank'].iloc[0])
        norm_list.append(temp["Normalized Average"].iloc[0])
    filt_df['primer on target normalized average'] = norm_list
    filt_df['primer on target pct rank'] = pct_list
    filt_df.to_csv(write_dir+"primer_on_target_count_combined_summary.tsv",sep="\t",index=False)

###DIVIDER

# The pipeline is currently as follows for 50 primer experiment:
# gen_lsv_bed -> #
def gen_lsv_bed(tsv_path="/Users/kyang/Dropbox/Rotation_2/Primers/caleb_primers/first_50_oligos_moreinfo.tsv",
                bed_path="/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/50_lsvs_n_plus_1.bed"):
    # given the .tsv of 50 oligos info, return an LSV bed
    # bed
    # currently this method is limited to target exons specifically
    # update 12_3_19: This was done wrong!! We reversed the order of coordinates for the upstream positive stranded
    # boundary
    bed = open(bed_path, 'w+')
    # bed.write("chr\tlower_bound\tupper_bound\n")
    with open(tsv_path) as inF:
        for line in inF:
            if "LSV_ID" not in line:
                linelist = line.split("\t")
                if linelist[17] != "":
                    # print(linelist[17])
                    # parse this into a bed file appropriate format and write out
                    # if positive stranded, junction corresponds to upstream exon boundary
                    subline = linelist[17].split(":")
                    # print(str(subline))
                    if subline[-1] == "+":
                        print(subline[0] + "\t" + str(int(subline[1].split("-")[0]) - 1) + "\t" + subline[1].split("-")[
                            0] + "\n")
                        bed.write(
                            subline[0] + "\t" + str(int(subline[1].split("-")[0]) - 1) + "\t" + subline[1].split("-")[
                                0] + "\n")
                    # if negative stranded, junction corresponds to downstream exon boundary
                    else:
                        print("pass")
                        print(subline[0] + "\t" + subline[1].split("-")[1] + "\t" + str(
                            int(subline[1].split("-")[1]) + 1) + "\n")
                        bed.write(subline[0] + "\t" + subline[1].split("-")[1] + "\t" + str(
                            int(subline[1].split("-")[1]) + 1) + "\n")

    bed.close()


def rewrite_lsv_bed(bed_path="/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/50_lsvs.bed",
                    rootdir="/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/"):
    # for some inane reason the bed file from gen_lsv_bed -> Liftover doesn't work so traverse it line by line writing a new file
    # then call bedtool afterwards
    os.chdir(rootdir)
    with open("temp.bed", 'w+') as temp:
        with open(bed_path) as inF:
            for line in inF:
                temp.write(line)


def traverse_lsv_bed(bed_path="/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/50_lsvs.bed",
                     rootdir="/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/"):
    # this is he funciton call to bedtools to see the interseciton between the bam and he temp file, to see
    # how many "ON target" reads there are
    os.chdir(rootdir)
    call("bedtools intersect -b KY001_deduplicated.bam -a temp.bed -c -bed >> intersect_lsvs.bed", shell=True)


def traverse_lsv_bed2(bed_path="/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/50_lsvs.bed",
                      rootdir="/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/"):
    # for some inane reason the bed file from gen_lsv_bed -> Liftover doesn't work so traverse it line by line in calling bedtools
    # otherwise, same funciton as above
    os.chdir(rootdir)
    with open(bed_path) as inF:
        for line in inF:
            with open("temp.bed", 'w+') as temp:
                temp.write(line)
            call("bedtools intersect -b KY001_deduplicated.bam -a temp.bed -c -bed >> intersect_lsvs.bed", shell=True)


def gen_primerable_seq_fasta(
        tsv_path="/Users/kyang/Dropbox/Rotation_2/Primers/caleb_primers/first_50_oligos_moreinfo_alphabet.tsv",
        fastq_path="/Users/kyang/Dropbox/Rotation_2/Primers/caleb_primers/50_lsvs_alphabet.fastq"):
    # given the .tsv of 50 oligos info, return an LSV fasta of "primerable region"
    fastq = open(fastq_path, 'w+')
    # bed.write("chr\tlower_bound\tupper_bound\n")
    with open(tsv_path) as inF:
        for line in inF:
            if "LSV_ID" not in line:
                linelist = line.split("\t")
                if linelist[13] != "":
                    print(linelist[13])
                    fastq.write(">" + linelist[2] + "\n" + linelist[13] + "\n")

    fastq.close()


def lsv_histo_data(bed_path="/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/50_lsvs_gene_labeled.txt",
                   rootdir="/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/lsv_histograms/"):
    # given the bed from gen_lsv_bed, iterate through this bed per LSV.
    # Per LSV, note and write out sum of bedtools intersect output with .bam for every 100 nucleotide range
    # from -2000 to +2000 bp away from the LSV.
    # Then, create a line histogram of reads from -2000 to +2000 bp per LSV
    os.chdir(rootdir)
    with open(bed_path) as inbed:
        for line in inbed:
            linelist = line.split("\t")
            zero = int(linelist[2][:-1])
            with open(rootdir + linelist[0] + "_data.txt", 'w+'):
                with open(rootdir + "temp.bed", 'w+') as temp_bed:
                    for x in range(-2000, 2001, 100):
                        temp = zero + x
                        temp_bed.write(linelist[1] + "\t" + str(temp) + "\t" + str(temp + 100) + "\n")
            call("bedtools intersect -b ../KY001_deduplicated.bam -a temp.bed -c -bed >> " + linelist[0] + "_data.txt",
                 shell=True)


def lsv_histogram(rootdir="/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/lsv_histograms/"):
    # given the output from lsv_histo_data,
    x_axis = [x for x in range(-2000, 2001, 100)]
    os.chdir(rootdir)
    for file in os.listdir(rootdir):
        if ".txt" in file:
            templist = []
            with open(file) as inF:
                for line in inF:
                    if line != "":
                        linelist = line.split("\t")
                        templist.append(int(linelist[-1][:-1]))
            print(len(x_axis))
            print(str(x_axis))
            print(len(templist))
            plt.plot(x_axis, templist)
    plt.savefig("./histo.png")
    plt.show()


def summarize_unextended(rootdir="/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/unextended_primer/"):
    os.chdir(rootdir)
    out = open(rootdir + "summary.txt", 'w+')
    for file in os.listdir(rootdir):
        if ".txt" in file and "_1_" in file:
            with open(file) as inF:
                for line in inF:
                    if "Matched" in line:
                        out.write(line.split("\t")[1] + "\n")


def off_target_coverage_bed(root="/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/off-target/",
                            chromosome_index="/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/off-target/chr_index"):
    # Generate bed with counts to plot histograms of 1000 KB and 10000 KB bins for each chromosome
    # Compare bulk RNA Seq distribution vs targeted Seq

    # First, get the max genomic coordinates of each chromosome
    os.chdir(root)
    index = dict()
    with open(chromosome_index) as inD:
        for line in inD:
            index[line.split("\t")[0]] = line.split("\t")[1][:-1]
    # write out a bed with all of the chosen increment for the bins
    print(str(index))
    bin_bed = open(root + "bins.bed", 'w+')
    n = 1000  # define bin size here
    for (k, v) in index.items():
        for x in range(1, int(v), n):
            bin_bed.write("chr" + k + "\t" + str(x) + "\t" + str(x + 1000) + "\n")
        bin_bed.write(k + "\t" + str(x) + "\t" + str(v) + "\n")
    bin_bed.close()
    # run bed intersect on this bed with the bam file
    bin_bed_counts = open(root + "bin_counts.bed", 'w+')
    bin_bed_counts.close()
    test = open(root + "test.bam", 'w+')
    test.close()
    call("bedtools intersect -a  ../KY001_deduplicated.bam -b ../50genes.bed -v >> KY001_off_target.bam", shell=True)
    call("bedtools intersect -b KY001_off_target.bam -a bins.bed -c -bed >> bin_counts.bed", shell=True)
    # now plot histogram of the graphs to get a sense of distribution
    bin_reads = []
    with open("bin_counts.bed") as countfile:
        for line in countfile:
            linelist = line.split("\t")
            if "\n" in linelist[-1]:
                linelist[-1] = linelist[-1][:-1]
            count = int(linelist[-1])
            if count > 0 and count < 100:
                bin_reads.append(count)
    plt.hist(bin_reads, cumulative=True, histtype='step', alpha=0.8, color='k')
    print(str(sum(bin_reads)))
    plt.savefig("./histo.png")  # plt.show()


def off_target_coverage_bed_process(root="/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/off-target/",
                                    chromosome_index="/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/off-target/chr_index"):
    # Generate bed with counts to plot histograms of 1000 KB and 10000 KB bins for each chromosome
    # Compare bulk RNA Seq distribution vs targeted Seq

    # First, get the max genomic coordinates of each chromosome
    os.chdir(root)
    # now plot histogram of the graphs to get a sense of distribution
    bin_reads = []
    total_bin_reads = []
    high_bins = open("high_bins.txt", 'w+')
    with open("bin_counts.bed") as countfile:
        for line in countfile:
            linelist = line.split("\t")
            if "\n" in linelist[-1]:
                linelist[-1] = linelist[-1][:-1]
            count = int(linelist[-1])
            if count > 1000:
                bin_reads.append(count)
                high_bins.write(line)
            if count > 0:
                total_bin_reads.append(count)
    plt.hist(bin_reads, alpha=0.8, color='k')
    print(str(sum(bin_reads)))
    print(str(len(bin_reads)))
    print(str(sum(total_bin_reads)))
    print(str(len(total_bin_reads)))
    plt.savefig("./histo.png")
    high_bins.close()  # plt.show()


def gen_primer_fasta(primerpath="/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/blast/50_primers",
                     rootpath="/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/blast/"):
    # from caleb file, generate a primer FASTA file
    os.chdir(rootpath)
    with open(primerpath) as inF:
        with open("50primers.fasta", 'w+') as inG:
            for line in inF:
                linelist = line.split("\t")
                inG.write(">" + linelist[0] + "\n" + linelist[1][:-1] + "\n")


def gen_primer_list2(primerpath="/Users/kyang/Dropbox/Rotation_2/MPE-Seq/s_pombe/s_pombe_primers.txt",
                     rootpath="/Users/kyang/Dropbox/Rotation_2/MPE-Seq/s_pombe/"):
    # from caleb file, generate a primer FASTA file
    os.chdir(rootpath)
    with open(primerpath) as inF:
        with open("primers.txt", 'w+') as inG:
            for line in inF:
                inG.write(line[line.rfind("N") + 1:-1] + ",")


def blast_off_target_unmapped(db=False, primerpath="/scr1/users/yangk4/KY001/50primers.fasta",
                              rootpath="/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/blast/",
                              prefix="/scr1/users/yangk4/KY001/", dbpath="/scr1/users/yangk4/KY001/blastdb",
                              nthreads="5"):
    # Given an input file of primers, generate a bash file for HPC which will
    # blast each primer against unmapped and targeted separately/in parallel
    # db = False if fasta not generated yet into blastdb, otherwise the path provided will be the path to db
    # primerpath = FASTA file of primers to BLASTN-short for
    # rootpath = file on LOCAL computer where script is generated
    # prefix = path prefix on HOST computer where script is run
    # dbpath = path to the dbs
    # numthreads = number of threads to use
    with open(rootpath + "blast_off_target_unmapped", 'w+') as inBlast:
        if db == False:
            inBlast.write('''# !/bin/bash
#$ -l h_vmem=8G
#$ -l m_mem_free=8G
#$ -pe smp 10 -binding linear:10
cd %s
#Run makeblastdb for both the unmapped reads as well as off_target reads
mkdir blastdb
/home/yangk4/ncbi-blast-2.9.0+/bin/makeblastdb -in %sKY001_off_target.fa -dbtype nucl -hash_index -out %sblastdb/KY001_off_target -max_file_sz 4GB -logfile %sblastdb/off_target.log &
/home/yangk4/ncbi-blast-2.9.0+/bin/makeblastdb -in %sKY001_unmapped.fa -dbtype nucl -hash_index -out %sblastdb/KY001_unmapped -max_file_sz 4GB -logfile %sblastdb/unmapped.log &
wait
#Now BLAST the primer file against the generated blastdbs
mkdir blastn
/home/yangk4/ncbi-blast-2.9.0+/bin/blastn -query %s -task blastn-short -db %sblastdb/KY001_off_target -num_threads %s -outfmt 1 -num_descriptions 100000 -num_alignments 100000 -out %sblastn/off_target_blast &
/home/yangk4/ncbi-blast-2.9.0+/bin/blastn -query %s -task blastn-short -db %sblastdb/KY001_unmapped -num_threads %s -outfmt 1 -num_descriptions 100000 -num_alignments 100000 -out %sblastn/unmapped_blast &
wait
''' % (prefix, prefix, prefix, prefix, prefix, prefix, prefix, primerpath, prefix, nthreads, prefix, primerpath, prefix,
       nthreads, prefix))


def blast_summarize(rootpath="/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/blast/"):
    # from the BLAST text file, summarize number of counts for off and on target, and also summarize the total
    # Also make graphs to compare number of on target reads, off target, and unmapped reads per gene
    # rootpath = path to both off and on target BLAST output text files
    os.chdir(rootpath)
    totaldict = dict()
    index = False
    inQuery = False
    output = open(rootpath + "summary.txt", 'w+')
    for i in ['off_target_blast', 'unmapped_blast']:
        with open(rootpath + i) as inF:
            for line in inF:
                if "Query=" in line:
                    if index == True:
                        totaldict[Query] = count
                    inQuery = True
                    index = True
                    Query = line.split("=")[1][1:-1]
                    count = 0
                if inQuery == True and ":" in line:
                    count += 1
        output.write(i + "\n")
        for k, v in totaldict.items():
            output.write(k + "\t" + str(v) + "\n")
        plt.hist([int(v) for v in totaldict.values()], label=i.replace("_", " "))
    plt.legend()
    plt.xlabel("Number of Reads mapping by BLAST")
    plt.savefig(i + "_histogram.png")
    output.close()


if __name__ == '__main__':
    extract_primer_origin()
