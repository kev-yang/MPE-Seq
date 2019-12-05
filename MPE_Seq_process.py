'''
Created on Sep 16, 2019

@author: kyang
'''

import os
from typing import List
from glob import glob
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import numpy as np
import pandas as pd
import seaborn as sns
sns.set()
import matplotlib.pyplot as plt
from subprocess import check_output
#Above the divider is the collated pipeline I use for KY002, below is for KY001
#For KY002, we first want some sanity checks. Here are some initial objectives:
#Overall sequence composition - how does the library look overall?
#Do this analysis per library, and then comparing the library prep methods stim/unstim to each other
#%Mapped and deduplicated, final number of reads left
#Percent nucleotide composition of mapped sequences (should this be all mapped or just targeted events?)
#PCA analysis of 1) all genes on target/off target read count, 2) psi values of targeted genes

#On vs off-target
#Extract the targeted genes and LSVs out of the mapped file
#Per LSV, make scatter plots of number of "on target" and "off target" primed reads versus the following

#Per targeted LSV, how many
def gen_lsv_bed_1(gtf_path="/Users/kyang/Box Sync/Rotation_2/ref/gencode.v31.annotation.gtf.bed",
                tsv_path="/Users/kyang/Box Sync/Rotation_2/Primers/KY002/",
                out_path="/Users/kyang/Box Sync/Rotation_2/Primers/KY002/KY002_on_target_"):
    # given the folder to the .tsvs of KY002, return separately 2 files:
    #1) a list of the gene ENSG ID's and 2) a bed file of targeted junctions
    # currently this method is limited to target exons specifically
    #gtf_path = path to the folder containing the gencode gtf file converted to bed via bedops (see Analysis/11_29_19.ipynb)
    #tsv_path = path to the folder of all tsvs from previous classify_process method in OligoMiner_wrapper
    #out_path = the prefix for generating the 2 files, one of which will be a simple text file (list of ENSGs) and
    #the other of which will be a properly formatted .bed file

    #Additionally, it would be helpful to have a third file which outputs coordinates to determine
    #reads which are extended at least from 1 base pair upstream of the primer
    #up to 1 base pair downstream of the exon junction coordinate

    #The last useful measures would be coordinates to compute psi values:
    #Namely this would be coordinates of 1 splicing event over the other... think about this more
    #it will be a fraction of anything that overlaps with 1 basepair before the junction
    #bed.write("chr\tlower_bound\tupper_bound\n")
    gtf_ref = open(gtf_path).readlines()
    os.chdir(tsv_path)
    file_list = [x for x in os.listdir() if ".tsv" in x]
    gtf_gene_ref = []
    with open(gtf_path) as inG:
        for line in inG:
            if line.split("\t")[7] == "gene" and "protein_coding" in line:
                gtf_gene_ref.append(line)
    print(str(len(gtf_gene_ref)))
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
                        subline_list[2] = "".join([subline_list[2].split("_")[0]+subline_list[2][-3:]])
                    #First, find the corresponding gene so we can extract info about chromosome
                    #and also write out the bed file for the gene
                    gene = subline_list[0].split(".")[0]
                    linelist = []
                    for line2 in gtf_ref:
                        if gene in line2 and line2.split("\t")[7] == "gene":
                            linelist = line2.split("\t")
                            chrom = linelist[0]
                            genelist.append("\t".join(linelist[:3]+[gene]))
                            if line2 in gtf_gene_ref:
                                count += 1
                                gtf_gene_ref.remove(line2)
                            break
                    else:
                        print("Error")
                    # if positive stranded, junction corresponds to upstream exon boundary
                    if "+" in subline_list[-1]:
                        bedlist.append(chrom + "\t" + subline_list[2].split("-")[0] + "\t" +
                                       str(int(subline_list[2].split("-")[0]) + 1) + "\t" + line.split("\t")[1])
                        extendlist.append(chrom + "\t" + subline_list[2].split("-")[0] + "\t" +
                                  str(int(subline_list[2].split("-")[0]) + int(line.split("\t")[2])) + "\t" + line.split("\t")[1])
                    # if negative stranded, junction corresponds to downstream exon boundary
                    else:
                        bedlist.append(chrom + "\t" + str(int(subline_list[2].split("-")[1][:-1]) - 1) + "\t" +
                                       subline_list[2].split("-")[1][:-1] + "\t" + line.split("\t")[1])
                        extendlist.append(chrom + "\t" + str(int(subline_list[2].split("-")[1][:-1]) - int(line.split("\t")[2])) + "\t" +
                                          subline_list[2].split("-")[1][:-1] + "\t" + line.split("\t")[1])
        bed.write("\n".join(bedlist))
        genes.write("\n".join(genelist))
        extended.write("\n".join(extendlist))
        bed.close()
        genes.close()
        extended.close()
    with open(out_path + "off_target_genes.bed", 'w+') as off_target:
        gtf_gene_ref_mod = ["\t".join(x.split("\t")[:4]) for x in gtf_gene_ref]
        off_target.write("\n".join(gtf_gene_ref_mod))
        print(str(count))
def intersect_script(out_dir="/Users/kyang/Box Sync/Rotation_2/MPE-Seq/KY002/",
                     scratch_dir=" /scr1/users/yangk4/scratch/KY002/",
                     n_threads=3,
                     bed_string="KY002_on_target_hi_confidence_large_change_primers_extended.bed KY002_on_target_hi_confidence_large_change_primers_genes.bed KY002_on_target_hi_confidence_large_change_primers_LSVs.bed KY002_on_target_hi_confidence_small_change_primers_extended.bed KY002_on_target_hi_confidence_small_change_primers_genes.bed KY002_on_target_hi_confidence_small_change_primers_LSVs.bed KY002_on_target_low_confidence_large_change_primers_extended.bed KY002_on_target_low_confidence_large_change_primers_genes.bed KY002_on_target_low_confidence_large_change_primers_LSVs.bed KY002_on_target_low_confidence_small_change_primers_extended.bed KY002_on_target_low_confidence_small_change_primers_genes.bed KY002_on_target_low_confidence_small_change_primers_LSVs.bed",
                     off_target_bed="off_target_genes.bed",
                     folder_string="Biotin-RepA-Stim/ Biotin-RepA-Unstim/ Biotin-RepC-Stim/ Biotin-RepD-Stim/ Biotin-RepD-Unstim/ TSO-RepA-Stim/ TSO-RepA-Unstim/ TSO-RepC-Stim/ TSO-RepC-Unstim/ TSO-RepD-Stim/ TSO-RepD-Unstim/"):
    #scratch_dir is the directory containing .bed in the top level directory and subdirectories of bams from MPE_Seq_pipeline
    #NOTE: you CANNOT use symlinks for scratch dir otherwise qsub has an Eqw error I think because of permissions
    #use the hard coded path
    #bed_string is the output from echo *.bed in bash EXCLUDING the off-target bed
    #off_target_bed is the name of the off-target bed
    #folder_string is the output from echo */ in bash
    #n_threads is number of processes to allow for running
    #Now that appropriate beds have been generated, generate a script which can run bedtools across all these beds in parallel
    fw = open(out_dir+ "bed_process", 'w+')
    fw.write(
'''#! /bin/bash
#$ -V
#$ -wd %s
#$ -l h_vmem=%sG
#$ -l m_mem_free=%sG
#$ -pe smp %s -binding linear:%s
module load gcc8
module load BEDTools
cd %s
''' %(scratch_dir,str(30),str(30),str(n_threads+1),str(n_threads+1),scratch_dir))
    srr_list = folder_string.split(" ")
    for folder in srr_list:
        templist = []
        folder = folder.split("/")[0]
        with open(out_dir+'intersect_'+folder,'w+') as inF:
            for bed in bed_string.split(" "):
                templist.append("bedtools intersect -b %s/%s_deduplicated.bam -a %s -c -bed >> %s/%s"
                                %(folder,folder,bed,folder,"intersect_"+bed))
            inF.write("\n".join(templist))
    chunk_list = [srr_list[i:i + n_threads] for i in range(0, len(srr_list), n_threads)]
    for n in range(len(chunk_list)):
        fw.write('#Processing chunk %s\n'% n)
        chunk = chunk_list[n] #this is not right
        for samp in chunk:
            fw.write('bash intersect_%s &\n'% samp.split("/")[0])
        fw.write('wait\n')
    fw.close()
    fw2 = open(out_dir+ "bed_process_off_target", 'w+')
    fw2.write(
'''#! /bin/bash
#$ -V
#$ -wd %s
#$ -l h_vmem=%sG
#$ -l m_mem_free=%sG
#$ -pe smp %s -binding linear:%s
module load gcc8
module load BEDTools
cd %s
''' %(scratch_dir,str(30),str(30),str(n_threads+1),str(n_threads+1),scratch_dir))
    srr_list = folder_string.split(" ")
    chunk_list = [srr_list[i:i + n_threads] for i in range(0, len(srr_list), n_threads)]
    for n in range(len(chunk_list)):
        fw2.write('#Processing chunk %s\n'% n)
        chunk = chunk_list[n] #this is not right
        for samp in chunk:
            folder = samp.split("/")[0]
            tempstring = "bedtools intersect -b %s/%s_deduplicated.bam -a %s -c -bed >> %s/%s" % (
            folder, folder, off_target_bed, folder, "intersect_" + off_target_bed)
            fw2.write('%s &\n'%tempstring)
        fw2.write('wait\n')
    fw2.close()
def summarize_intersect_pca(intersect_path="/Users/kyang/Box Sync/Rotation_2/MPE-Seq/KY002/intersect_beds",
                        read_counts_path="/Users/kyang/Box Sync/Rotation_2/MPE-Seq/KY002/library_read_counts"):
    #intersect_path is the path to a folder containing subdirectories, each of which
    #read_counts_path is the path to the number of read counts in every library (see block comment below)
    '''
    #to get readcounts run the following bash script:
    dirs=$(echo */)
    for i in $dirs
    do
        l1=$(samtools idxstats $i*_deduplicated.bam | awk '{total += $3} END{print total}')
        echo -e ${i::-1}"\t"$l1 >> library_read_counts
    done
    '''
    #
    #contains the bed files for a single library prep
    #after the above intersect script has been run, output some summary statistics
    #and basic graphs summarizing number of unextended, extended, and per library
    #consider combining TSO libraries into 1 average and Pleiss libraries into a 2nd average
    #if they look broadly similar


    #To determine if they look broadly similar, do PCA for the libraries for number of reads corresponding to
    #each gene out of total number of reads
    #See https://towardsdatascience.com/pca-using-python-scikit-learn-e653f8989e60
    #Also do for number of reads corresponding to LSVs and number of reads corresponding to unextended primers.
    read_count_dict = dict()
    with open(read_counts_path) as inR:
        for line in inR:
            linelist=line.split("\t")
            read_count_dict[linelist[0]] = int(linelist[1])
    #X is a M*N-dimensional matrix where M is the number of libraries and N is the different genes
    #with Y being the libraries
    #Parse the data to get this desired matrix (first testing for just gene expression to make sure it works)
    os.chdir(intersect_path)
    intersect_dirlist = [x for x in os.listdir() if os.path.isdir(x)] #list directories to iterate here
    #Do each of genes, extended, LSV
    X_all = [] # X_all is a list of arrays with each array being either  genes, extended, LSVs, or off-target
    for k in ["genes", "extended", "LSVs"]:
        lib_names = []
        c2 = 0
        for f in intersect_dirlist:
            c2 += 1
            lib_names.append(f.split("/")[0])  # lib names is the list tracking names of things
            os.chdir(intersect_path)
            os.chdir(f)
            #Need to concatenate on-target-beds into same row
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
                    final1 = pd.concat([final1,temp],axis=1)
                    final1 = final1.loc[:,~final1.columns.duplicated()]
                # After one row is finished, move on to the next one
            final1 = final1/read_count_dict[f.split("/")[0]] #divide by total read counts to normalize between libraries
            if c2 == 1:
                final_on = final1
            else:
                final_on = final_on.append(final1) #divide by read counts to normalize between libraries
        X_all.append(final_on.rename_axis('ID').values)
    #Now, just 1 off-target-bed (but do same thing for it)
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
            final1 = final1/read_count_dict[f.split("/")[0]] # add a pseudocount of 1 and divide by total read counts to normalize between libraries
        if c2 == 1:
            final_off = final1
        else:
            final_off = final_off.append(final1)
    X_all.append(final_off.rename_axis('ID').values)

    #this does PCA
    for X in X_all:
        #Standardscaler to standardize
        X = StandardScaler().fit_transform(X)
        pca = PCA(n_components = 2)
        principalComponents = pca.fit_transform(X)
        principalDf = pd.DataFrame(data=principalComponents
                                   ,columns=['principal component 1', 'principal component 2'])
        library_prep = []
        for l in lib_names:
            if "TSO" in l:
                l1 = "TSO"
            else:
                l1 = "Biotin"
            if "Unstim" in l:
                l2 = "Unstim"
            else:
                l2 = "Stim"
            library_prep.append(l1+", "+l2)
        principalDf.insert(0,"library_type",library_prep,True)
        #just plotting stuff
        #We color different stim and unstim conditions differently
        #Label all points with what they are called
        d = sns.scatterplot(x="principal component 1",y="principal component 2",hue="library_type",data=principalDf)
        #add labels, see https://stackoverflow.com/questions/14432557/matplotlib-scatter-plot-with-different-text-at-each-data-point
        for i, txt in enumerate(lib_names):
            d.annotate(txt, (np.asarray(principalDf["principal component 1"])[i], np.asarray(principalDf["principal component 2"])[i]))
        plt.show()
    #We can try to correlate the on-target reads with various measures
    #and do the same thing with off-target reads
def intersect_
    #Further process X_all make ensembl averages per


    #Now that we have these


    #Can do the same thing except with fold enrichment vs Kallisto-estimated RNASeq-expression
    
    #PSI calculations
    #-Calculate psi for all junctions as a simple measure independent of MAJIQ and see how well it agrees
    #Per LSV, we want the intersect of a given splice junction immediately be
    
    
    
    
    ###DIVIDER
    
    #The pipeline is currently as follows for 50 primer experiment:
    #gen_lsv_bed ->
    def gen_lsv_bed(tsv_path="/Users/kyang/Dropbox/Rotation_2/Primers/caleb_primers/first_50_oligos_moreinfo.tsv",
                    bed_path="/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/50_lsvs_n_plus_1.bed"):
        #given the .tsv of 50 oligos info, return an LSV bed
        #bed 
        #currently this method is limited to target exons specifically
        #update 12_3_19: This was done wrong!! We reversed the order of coordinates for the upstream positive stranded
        #boundary
        bed = open(bed_path,'w+')
        #bed.write("chr\tlower_bound\tupper_bound\n")
        with open(tsv_path) as inF:
            for line in inF:
                if "LSV_ID" not in line:
                    linelist = line.split("\t")
                    if linelist[17] != "":
                        #print(linelist[17])
                        #parse this into a bed file appropriate format and write out
                        #if positive stranded, junction corresponds to upstream exon boundary
                        subline = linelist[17].split(":")
                        #print(str(subline))
                        if subline[-1] == "+":
                            print(subline[0]+"\t"+str(int(subline[1].split("-")[0])-1)+"\t"+subline[1].split("-")[0]+"\n")
                            bed.write(subline[0]+"\t"+str(int(subline[1].split("-")[0])-1)+"\t"+subline[1].split("-")[0]+"\n")
                        #if negative stranded, junction corresponds to downstream exon boundary
                        else:
                            print("pass")
                            print(subline[0]+"\t"+subline[1].split("-")[1]+"\t"+str(int(subline[1].split("-")[1])+1)+"\n") 
                            bed.write(subline[0]+"\t"+subline[1].split("-")[1]+"\t"+str(int(subline[1].split("-")[1])+1)+"\n")    
                                           
        bed.close()
    def rewrite_lsv_bed(bed_path="/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/50_lsvs.bed",
                         rootdir="/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/"):
            #for some inane reason the bed file from gen_lsv_bed -> Liftover doesn't work so traverse it line by line writing a new file
            #then call bedtool afterwards
            os.chdir(rootdir)
            with open("temp.bed",'w+') as temp:
                with open(bed_path) as inF:
                    for line in inF:
                            temp.write(line)
    def traverse_lsv_bed(bed_path="/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/50_lsvs.bed",
                         rootdir="/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/"):
        #this is he funciton call to bedtools to see the interseciton between the bam and he temp file, to see
        #how many "ON target" reads there are
        os.chdir(rootdir)
        call("bedtools intersect -b KY001_deduplicated.bam -a temp.bed -c -bed >> intersect_lsvs.bed", shell=True)
    
    
    def traverse_lsv_bed2(bed_path="/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/50_lsvs.bed",
                         rootdir="/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/"):
        # for some inane reason the bed file from gen_lsv_bed -> Liftover doesn't work so traverse it line by line in calling bedtools
        #otherwise, same funciton as above
        os.chdir(rootdir)
        with open(bed_path) as inF:
            for line in inF:
                with open("temp.bed", 'w+') as temp:
                    temp.write(line)
                call("bedtools intersect -b KY001_deduplicated.bam -a temp.bed -c -bed >> intersect_lsvs.bed", shell=True)
    
    def gen_primerable_seq_fasta(tsv_path="/Users/kyang/Dropbox/Rotation_2/Primers/caleb_primers/first_50_oligos_moreinfo_alphabet.tsv",
            fastq_path="/Users/kyang/Dropbox/Rotation_2/Primers/caleb_primers/50_lsvs_alphabet.fastq"):
        #given the .tsv of 50 oligos info, return an LSV fasta of "primerable region"
        fastq = open(fastq_path,'w+')
        #bed.write("chr\tlower_bound\tupper_bound\n")
        with open(tsv_path) as inF:
            for line in inF:
                if "LSV_ID" not in line:
                    linelist = line.split("\t")
                    if linelist[13] != "":
                        print(linelist[13])
                        fastq.write(">"+linelist[2]+"\n"+linelist[13]+"\n")  
                                           
        fastq.close()
        
        
    def lsv_histo_data(bed_path="/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/50_lsvs_gene_labeled.txt",
                      rootdir="/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/lsv_histograms/"):
        #given the bed from gen_lsv_bed, iterate through this bed per LSV.
        #Per LSV, note and write out sum of bedtools intersect output with .bam for every 100 nucleotide range
        #from -2000 to +2000 bp away from the LSV.
        #Then, create a line histogram of reads from -2000 to +2000 bp per LSV
        os.chdir(rootdir)
        with open(bed_path) as inbed:
            for line in inbed:
                linelist = line.split("\t")
                zero = int(linelist[2][:-1])
                with open(rootdir+linelist[0]+"_data.txt",'w+'):
                    with open(rootdir+"temp.bed",'w+') as temp_bed:
                        for x in range(-2000,2001,100):
                            temp = zero+x
                            temp_bed.write(linelist[1]+"\t"+str(temp)+"\t"+str(temp+100)+"\n")
                call("bedtools intersect -b ../KY001_deduplicated.bam -a temp.bed -c -bed >> "+linelist[0]+"_data.txt",shell=True)
    def lsv_histogram(rootdir="/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/lsv_histograms/"):
        #given the output from lsv_histo_data, 
        x_axis = [x for x in range(-2000,2001,100)]
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
                plt.plot(x_axis,templist)
        plt.savefig("./histo.png")
        plt.show()
                   
    def summarize_unextended(rootdir="/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/unextended_primer/"):
        os.chdir(rootdir)
        out = open(rootdir+"summary.txt",'w+')
        for file in os.listdir(rootdir):
            if ".txt" in file and "_1_" in file:
                with open(file) as inF:
                    for line in inF:
                        if "Matched" in line:
                            out.write(line.split("\t")[1]+"\n")
    def off_target_coverage_bed(root = "/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/off-target/",
                                chromosome_index = "/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/off-target/chr_index"):
        #Generate bed with counts to plot histograms of 1000 KB and 10000 KB bins for each chromosome
        #Compare bulk RNA Seq distribution vs targeted Seq
        
        #First, get the max genomic coordinates of each chromosome
        os.chdir(root)
        index = dict()
        with open(chromosome_index) as inD:
            for line in inD:
                index[line.split("\t")[0]] = line.split("\t")[1][:-1]
        #write out a bed with all of the chosen increment for the bins
        print(str(index))
        bin_bed = open(root+"bins.bed",'w+')
        n = 1000 #define bin size here
        for (k,v) in index.items():
            for x in range(1,int(v),n):
                bin_bed.write("chr"+k+"\t"+str(x)+"\t"+str(x+1000)+"\n")
            bin_bed.write(k+"\t"+str(x)+"\t"+str(v)+"\n")
        bin_bed.close()
        #run bed intersect on this bed with the bam file
        bin_bed_counts = open(root+"bin_counts.bed",'w+')
        bin_bed_counts.close()
        test = open(root+"test.bam",'w+')
        test.close()
        call("bedtools intersect -a  ../KY001_deduplicated.bam -b ../50genes.bed -v >> KY001_off_target.bam",shell=True)
        call("bedtools intersect -b KY001_off_target.bam -a bins.bed -c -bed >> bin_counts.bed",shell=True)
        #now plot histogram of the graphs to get a sense of distribution
        bin_reads = []
        with open("bin_counts.bed") as countfile:
            for line in countfile:
                linelist = line.split("\t")
                if "\n" in linelist[-1]:
                    linelist[-1] = linelist[-1][:-1]
                count = int(linelist[-1])
                if count >0 and count < 100:
                    bin_reads.append(count)
        plt.hist(bin_reads, cumulative=True,histtype='step', alpha=0.8, color='k')
        print(str(sum(bin_reads)))
        plt.savefig("./histo.png")
        #plt.show()
                        
    def off_target_coverage_bed_process(root = "/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/off-target/",
                                chromosome_index = "/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/off-target/chr_index"):
        #Generate bed with counts to plot histograms of 1000 KB and 10000 KB bins for each chromosome
        #Compare bulk RNA Seq distribution vs targeted Seq
        
        #First, get the max genomic coordinates of each chromosome
        os.chdir(root)
        #now plot histogram of the graphs to get a sense of distribution
        bin_reads = []
        total_bin_reads = []
        high_bins = open("high_bins.txt",'w+')
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
        high_bins.close()
        #plt.show()
    def gen_primer_fasta(primerpath = "/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/blast/50_primers",
                         rootpath = "/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/blast/"):
        #from caleb file, generate a primer FASTA file
        os.chdir(rootpath)
        with open(primerpath) as inF:
            with open("50primers.fasta",'w+') as inG:
                for line in inF:
                    linelist = line.split("\t")
                    inG.write(">"+linelist[0]+"\n"+linelist[1][:-1]+"\n")
    
    def gen_primer_list2(primerpath = "/Users/kyang/Dropbox/Rotation_2/MPE-Seq/s_pombe/s_pombe_primers.txt",
                         rootpath = "/Users/kyang/Dropbox/Rotation_2/MPE-Seq/s_pombe/"):
        #from caleb file, generate a primer FASTA file
        os.chdir(rootpath)
        with open(primerpath) as inF:
            with open("primers.txt",'w+') as inG:
                for line in inF:
                    inG.write(line[line.rfind("N")+1:-1]+",")
    
    def blast_off_target_unmapped(db = False,
                                  primerpath = "/scr1/users/yangk4/KY001/50primers.fasta",
                                  rootpath = "/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/blast/",
                                  prefix = "/scr1/users/yangk4/KY001/",
                                  dbpath = "/scr1/users/yangk4/KY001/blastdb",
                                  nthreads = "5"):
        #Given an input file of primers, generate a bash file for HPC which will
        #blast each primer against unmapped and targeted separately/in parallel
        #db = False if fasta not generated yet into blastdb, otherwise the path provided will be the path to db
        #primerpath = FASTA file of primers to BLASTN-short for
        #rootpath = file on LOCAL computer where script is generated
        #prefix = path prefix on HOST computer where script is run
        #dbpath = path to the dbs
        #numthreads = number of threads to use
        with open(rootpath+"blast_off_target_unmapped",'w+') as inBlast:
            if db == False:
                inBlast.write(
    '''# !/bin/bash
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
'''% (prefix,prefix,prefix,prefix,prefix,prefix,prefix,primerpath,prefix,nthreads,prefix,primerpath,prefix,nthreads,prefix)
)

def blast_summarize(rootpath = "/Users/kyang/Dropbox/Rotation_2/MPE-Seq/09_10_19/blast/"):
    #from the BLAST text file, summarize number of counts for off and on target, and also summarize the total
    #Also make graphs to compare number of on target reads, off target, and unmapped reads per gene
    #rootpath = path to both off and on target BLAST output text files
    os.chdir(rootpath)
    totaldict = dict()
    index = False
    inQuery = False
    output = open(rootpath+"summary.txt",'w+')
    for i in ['off_target_blast','unmapped_blast']:
        with open(rootpath+i) as inF:
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
        output.write(i+"\n")
        for k,v in totaldict.items():
            output.write(k+"\t"+str(v)+"\n")
        plt.hist([int(v) for v in totaldict.values()],label=i.replace("_"," "))
    plt.legend()
    plt.xlabel("Number of Reads mapping by BLAST")
    plt.savefig(i+"_histogram.png")
    output.close()