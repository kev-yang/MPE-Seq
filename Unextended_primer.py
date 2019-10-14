'''
Created on Sep 14, 2019

@author: kyang
'''

import shutil
from subprocess import *
#This script will extract unextended primer via several calls to bbduk, allowing for n mismatches
#First, call 

def pretrim(output_file = "mpe_seq",
                           files_to_process="process_list",
                           root_dir="/scr1/users/yangk4/mpe_seq/",
                           nprocs = "8",
                           memory = "40"):
    #pre-trim to generate the required fasta for remove_primer
    fw = open(output_file+"_unextended_process", 'w+')
    fw.write(
'''# !/bin/bash
#$ -l h_vmem=5G
#$ -l m_mem_free=5G
#$ -pe smp 8 -binding linear:8
#Now do real things
cd %s
mkdir STAR
for i in `cat %s`;
do
    echo "#Trim with BBDUK"
    echo "mkdir ./trimmed" >> DL_and_process_$i
    echo "/home/yangk4/bbmap/bbduk.sh ref=/home/yangk4/bbmap/resources/adapters.fa in1=$i\_1.fastq in2=$i\_2.fastq out1=./trimmed/$i\_1_unextended_trim.fastq out2=./trimmed/$i\_2_unextended_trim.fastq ktrim=r k=23 mink=11 hdist=1 tpe tbo -Xmx%sg threads=%s" >> DL_and_process_$i
    echo "#Run UMI Tools Extract"
    echo "source /home/yangk4/majiq_2_install/env/bin/activate" >> DL_and_process_$i
    echo "umi_tools extract -I ./trimmed/{$i}_1_unextended_trim.fastq --bc-pattern=NNNNNNN --read2-in=./trimmed/${i}_2_unextended_trim.fastq --stdout=${i}_extract_u_1.fastq --read2-out=${i}_extract_u_2.fastq --log=$i\_UMI_extracted_unextended.log" >> DL_and_process_$i
    echo "deactivate" >> DL_and_process_$i
done

''' % (root_dir,files_to_process, #general params
memory,nprocs, #bbduck params
))
    srr_list = open(files_to_process, 'r').read().splitlines()
    chunk_list = [srr_list[i:i + group_size] for i in range(0, len(srr_list), group_size)]
    for n in range(len(chunk_list)):
        fw.write('#Processing chunk %s\n'% n)
        arr_string = ''
        chunk = chunk_list[n] #this is not right
        for samp in chunk:
            fw.write('bash DL_and_process_%s &\n'% samp)
            arr_string += '"%s" ' % samp
        fw.write('wait\n')
    fw.close()

def remove_primer(fastqpath = "/scr1/users/yangk4/KY001/KY001_extract_u",
                  primerlist = "/scr1/users/yangk4/KY001/50_primers",mismatch="7"):
    
    #Currently this script assumes paired end reads
    #fastqpath is the path to the .fastq prior to "_1.fastq" and "_2.fastq"
    #primerlist is the path to the list of primers to query for
    
    #Do some pre-processing on primer list to "trim" them to just the primed region and arrange them into sublists by length
    primerl = []
    with open(primerlist) as inF:
        for line in inF:
            primerl.append(line.split("\n")[0])
    sorted_primers = []
    for i in range(min([len(x) for x in primerl]),max([len(x) for x in primerl])+1):
        temp = []
        for j in primerl:
            if len(j) == i:
                temp.append(j)
        if len(temp) > 0:
            sorted_primers.append(temp)
    #Now, iterate over each sorted sublist arranged by length
    fileprefix = "/".join(fastqpath.split("/")[:-1])+"/"
    for temp in sorted_primers:
        #First, copy the fastq to some temp files and call bbduk to write out filtered reads to separate file
        shutil.copyfile(fastqpath+"_1.fastq","/"+fileprefix+"temp_1.fastq")
        shutil.copyfile(fastqpath+"_2.fastq","/"+fileprefix+"temp_2.fastq")
        with open(fileprefix+"temp_primers.fa",'w+') as primerF:
            k = 0
            for i in temp:
                k+=1
                primerF.write(">Primer_"+str(k)+"\n"+i+"\n")
            primerF.write("\n")
        print("/home/yangk4/bbmap/bbduk.sh "+fileprefix+"temp_1.fastq"+" out=clean.fq outm="+
             fastqpath+"_minlen_1.fastq"+" minlength="+str(len(temp[0])))
        call("/home/yangk4/bbmap/bbduk.sh "+fileprefix+"temp_1.fastq"+" out=clean.fq outm="+
             fastqpath+"_minlen_1.fastq"+" minlength="+str(len(temp[0])),shell=True)
        call("/home/yangk4/bbmap/bbduk.sh "+fileprefix+"temp_2.fastq"+" out=clean.fq outm="+
             fastqpath+"_minlen_2.fastq"+" minlength="+str(len(temp[0])),shell=True)
        #After filtering reads below a minimum length, look for matches to the templist
        call("/home/yangk4/bbmap/bbduk.sh "+fastqpath+"_minlen_1.fastq"+" out=clean.fq ref="+
             fileprefix+"temp_primers.fa"+" hdist=3"+" k="+str(len(temp[0])-10)+" stats="+str(len(temp[0]))+"_1_stat.txt",shell=True)
        call("/home/yangk4/bbmap/bbduk.sh "+fastqpath+"_minlen_2.fastq"+" out=clean.fq ref="+
             fileprefix+"temp_primers.fa"+" hdist=3"+" k="+str(len(temp[0])-10)+" stats="+str(len(temp[0]))+"_2_stat.txt",shell=True)

def remove_primer(fastqpath = "/scr1/users/yangk4/KY001/KY001_extract_u",
                  primerlist = "/scr1/users/yangk4/KY001/50_primers",mismatch="7"):
    
    #Currently this script assumes paired end reads
    #fastqpath is the path to the .fastq prior to "_1.fastq" and "_2.fastq"
    #primerlist is the path to the list of primers to query for
    
    #Do some pre-processing on primer list to "trim" them to just the primed region and arrange them into sublists by length
    primerl = []
    with open(primerlist) as inF:
        for line in inF:
            primerl.append(line.split("\n")[0])
    sorted_primers = []
    for i in range(min([len(x) for x in primerl]),max([len(x) for x in primerl])+1):
        temp = []
        for j in primerl:
            if len(j) == i:
                temp.append(j)
        if len(temp) > 0:
            sorted_primers.append(temp)
    #Now, iterate over each sorted sublist arranged by length
    fileprefix = "/".join(fastqpath.split("/")[:-1])+"/"
    for temp in sorted_primers:
        #First, copy the fastq to some temp files and call bbduk to write out filtered reads to separate file
        shutil.copyfile(fastqpath+"_1.fastq","/"+fileprefix+"temp_1.fastq")
        shutil.copyfile(fastqpath+"_2.fastq","/"+fileprefix+"temp_2.fastq")
        with open(fileprefix+"temp_primers.fa",'w+') as primerF:
            k = 0
            for i in temp:
                k+=1
                primerF.write(">Primer_"+str(k)+"\n"+i+"\n")
            primerF.write("\n")
        print("/home/yangk4/bbmap/bbduk.sh "+fileprefix+"temp_1.fastq"+" out=clean.fq outm="+
             fastqpath+"_minlen_1.fastq"+" minlength="+str(len(temp[0])))
        call("/home/yangk4/bbmap/bbduk.sh "+fileprefix+"temp_1.fastq"+" out=clean.fq outm="+
             fastqpath+"_minlen_1.fastq"+" minlength="+str(len(temp[0])),shell=True)
        call("/home/yangk4/bbmap/bbduk.sh "+fileprefix+"temp_2.fastq"+" out=clean.fq outm="+
             fastqpath+"_minlen_2.fastq"+" minlength="+str(len(temp[0])),shell=True)
        #After filtering reads below a minimum length, look for matches to the templist
        call("/home/yangk4/bbmap/bbduk.sh "+fastqpath+"_minlen_1.fastq"+" out=clean.fq ref="+
             fileprefix+"temp_primers.fa"+" hdist=3"+" k="+str(len(temp[0])-10)+" stats="+str(len(temp[0]))+"_1_stat.txt",shell=True)
        call("/home/yangk4/bbmap/bbduk.sh "+fastqpath+"_minlen_2.fastq"+" out=clean.fq ref="+
             fileprefix+"temp_primers.fa"+" hdist=3"+" k="+str(len(temp[0])-10)+" stats="+str(len(temp[0]))+"_2_stat.txt",shell=True)