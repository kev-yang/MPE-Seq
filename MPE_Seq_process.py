'''
Created on Sep 16, 2019

@author: kyang
'''
from subprocess import call
import os
import matplotlib.pyplot as plt
import numpy as np
#Above the divider is the collated pipeline I use for KY002, below is for KY001
#For KY002, we first want some sanity checks. Here are some initial objectives:
#Overall sequence composition - how does the library look overall?
#Do this analysis per library, and then comparing the library prep methods stim/unstim to each other
#%Mapped and deduplicated, final number of reads left
#Percent nucleotide composition of mapped sequences (should this be all mapped or just targeted events?)
#PCA analysis of 1) all genes, 2) psi values of targeted genes

#On vs off-target
#Extract the targeted genes and LSVs out of the mapped file
#Per LSV, make scatter plots of number of "on target" and "off target" primed reads versus the following

#Per targeted LSV, how many
def gen_lsv_bed(gtf_path="/Users/kyang/Box Sync/Rotation_2/ref/gencode.v31.annotation.gtf.bed",
                tsv_path="/Users/kyang/Box Sync/Rotation_2/Primers/KY002",
                out_path="/Users/kyang/Box Sync/Rotation_2/Primers/KY002/KY002_on_target"):
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
    bed = open(out_path+"_LSVs.bed", 'w+')
    genes = open(out_path+"_genes.bed", 'w+')
    extended = open(out_path+"_extended.bed", 'w+')
    #bed.write("chr\tlower_bound\tupper_bound\n")
    gtf_ref = open(gtf_path).readlines()
    genelist,bedlist,extendlist = [],[],[]
    with open(tsv_path) as inF:
        file_list = [x for x in os.listdir() if ".tsv" in x]
        for f in file_list:
            with open(f) as inF:
                for line in inF:
                    if "ENSG" in line:
                        subline_list = line.split("\t")[1].split(":")
                        #First, find the corresponding gene so we can extract info about chromosome
                        #and also write out the bed file for the gene
                        gene = subline_list[0].split(".")[0]
                        linelist = []
                        for line in gtf_ref:
                            if gene in line and "gene" in line:
                                linelist = line.split("\t")
                                chrom = linelist[0]
                                genelist.append("\t".join(linelist[:3]+[gene]))
                                break
                        else:
                            print("ERROR")
                        # if positive stranded, junction corresponds to upstream exon boundary
                        if "+" in subline_list[-1]:
                            bedlist.append(chrom + "\t" + str(int(subline_list[2].split("-")[0])-1) + "\t" +
                                      subline_list[2].split("-")[0] + "\t"+ subline_list.join(":"))
                            extendlist.append(chrom + "\t" + str(int(subline_list[2].split("-")[0])-(50-int(linelist[3]))) + "\t" +
                                      subline_list[2].split("-")[0] + "\t"+ subline_list.join(":"))
                        # if negative stranded, junction corresponds to downstream exon boundary
                        else:
                            bedlist.append(chrom + "\t" + str(int(subline_list[2].split("-")[1])+1) + "\t" +
                                      subline_list[2].split("-")[1] + "\t"+ subline_list.join(":"))
                            extendlist.append(chrom + "\t" + str(int(subline_list[2].split("-")[1]) + (50 - int(linelist[3]))) + "\t" +
                                              subline_list[2].split("-")[1] + "\t" + subline_list.join(":"))
    bed.write(bedlist.join("\n"))
    genes.write(genelist.join("\n"))
    extended.write(extendlist.join("\n"))
#PSI calculations
#-Calculate psi for all junctions as a simple measure independent of MAJIQ and see how well it agrees




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