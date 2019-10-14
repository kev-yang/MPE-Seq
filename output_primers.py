def MPE_Seq_paired_process(output_file = "KY001",
                           files_to_process="process_list",
                           root_dir="/scr1/users/yangk4/KY001/",
                           group_size = 100, 
                           genome = 'hg38',
                           genomepath = "/scr1/users/yangk4/ref/genomes/hg38",
                           gffpath = "/scr1/users/yangk4/ref/db/s288c.gff3",
                           gtfpath = "/scr1/users/yangk4/ref/db/s288c.gtf",
                           readlen =  "75",
                           nprocs = "5",
                           memory = "40"):


    #this method is for CHOP Respublica, it processes already downloaded data
    #output_file = the output script to download onto the server
    #files_to_process = text file list of .fastq file prefixes to run this on in parallel (ONLY prefix before .fastq)
    #For paired files, name both reads the same prefix but end in _1.fastq for Read 1 and _2.fastq for Read 2
    #root_dir = root directory containing the .fastq files to process
    #group_size = For large-scale parallel processing, how mnay to run at a time (not currently implemented), NOT a string
    #genome = name of genome used for MAJIQ (e.g. mm10, hg38)
    #genomepath = path of genome for STAR and MAJIQ
    #gffpath = path of gff file for MAJIQ
    #gtfpath = path of gtf file for RMATS (use cufflinks to convert gff3 to gtf, can also use gffread to convert gtf to gff3)
    #readlen = read length for MAJIQ
    #nprocs = number of threads to use for STAR, MAJIQ, BBDuk PER process
    #memory = total amount of memory requested for BBDuk, in gigabytes
    #BBDuk is hard-coded to use 7 threads so request at least that many
    #I didn't code the directories for installs as options, so you made need to go back and change those.
    #I have the following environments: Python 3 venv under ~/majiq_2_install to run UMI tools and MAJIQ, 
    #Python 2 venv under ~/rmats_install to run rMATS.
    fw = open(output_file+"_process", 'w+')
    
#     fw1.write(
# '''cd %s
# for i in `cat %s`;
# do
#     echo "wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${i:0: 6}/$i/$i.sra" >>DL_only_$i
# done\n''' % (out_dir, SRR_list_server_fileName))
#     
#     srr_list = open(SRR_list_file, 'r').read().splitlines()
#     chunk_list = [srr_list[i:i + group_size] for i in range(0, len(srr_list), group_size)]
#     
#     for n in range(len(chunk_list)):
#         fw1.write('#Processing chunk %s\n'% n)
#         arr_string = ''
#         chunk = chunk_list[n]
#         for samp in chunk:
#             fw1.write('bash DL_only_%s &\n'% samp)
#             arr_string += '"%s" ' % samp
#         fw1.write('wait\n')
#NOTE: Caleb normally has minlength=35 trimmed out, I have removed this line for now
#Apparently UMI tools also requires samtools indexing.

    fw.write(
'''# !/bin/bash
#$ -l h_vmem=8G
#$ -l m_mem_free=8G
#$ -pe smp 10
module unload binutils
module load STAR
module load SRA-Toolkit
module load jdk
module load python/3.6
module load gcc8
module load SAMtools
#Now do real things
cd %s
mkdir STAR
for i in `cat %s`;
do
    echo "#Run FASTQC"
    echo "mkdir ./fastqc_1" >> DL_and_process_$i
    echo "mkdir ./fastqc_2" >> DL_and_process_$i
    echo "/home/yangk4/FastQC/fastqc $i\_1.fastq -o ./fastqc_1 &" >> DL_and_process_$i
    echo "/home/yangk4/FastQC/fastqc $i\_2.fastq -o ./fastqc_2 &" >> DL_and_process_$i
    echo "wait" >> DL_and_process_$i
    echo "#Trim with BBDUK"
    echo "mkdir ./trimmed" >> DL_and_process_$i
    echo "/home/yangk4/bbmap/bbduk.sh ref=/home/yangk4/bbmap/resources/adapters.fa in1=$i\_1.fastq in2=$i\_2.fastq out1=./trimmed/$i\_1_trimmed.fastq out2=./trimmed/$i\_2_trimmed.fastq ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=20 qin=33 -Xmx%sg threads=%s" >> DL_and_process_$i
    echo "#Run FASTQC again"
    echo "mkdir ./trimmed/fastqc_1" >> DL_and_process_$i
    echo "mkdir ./trimmed/fastqc_2" >> DL_and_process_$i
    echo "/home/yangk4/FastQC/fastqc ./trimmed/$i\_1_trimmed.fastq -o ./trimmed/fastqc_1 &" >> DL_and_process_$i
    echo "/home/yangk4/FastQC/fastqc ./trimmed/$i\_2_trimmed.fastq -o ./trimmed/fastqc_2 &" >> DL_and_process_$i
    echo "#Run UMI Tools Extract"
    echo "source /home/yangk4/majiq_2_install/env/bin/activate" >> DL_and_process_$i
    echo "umi_tools extract -I ./trimmed/$i\_1_trimmed.fastq --bc-pattern=NNNNNNN --read2-in=./trimmed/$i\_2_trimmed.fastq --stdout=$i\_1_extract.fastq --read2-out=$i\_2_extract.fastq --log=$i\_UMI_extracted.log" >> DL_and_process_$i
    echo "deactivate" >> DL_and_process_$i
    echo "#Run STAR"
    echo "STAR --genomeDir %s  --readFilesIn ./$i\_1_extract.fastq ./$i\_2_extract.fastq --runThreadN %s --outSAMtype BAM Unsorted --outFileNamePrefix ./$i. --outSAMattributes All --alignSJoverhangMin 8 --outSAMunmapped Within" >> DL_and_process_$i
    echo "#Remove PCR Duplicates"
    echo "samtools sort -@ %s -o $i.Aligned.sorted.bam $i.Aligned.out.bam" >> DL_and_process_$i
    echo "samtools index ./$i.Aligned.sorted.bam" >> DL_and_process_$i
    echo "source /home/yangk4/majiq_2_install/env/bin/activate" >> DL_and_process_$i
    echo "umi_tools dedup -I $i.Aligned.sorted.bam --paired -S $i\_deduplicated.bam --log=$i\_UMI_deduplicated.log" >> DL_and_process_$i
    echo "deactivate" >> DL_and_process_$i
    echo "#Index with samtools for MAJIQ"
    echo "samtools index $i\_deduplicated.bam" >> DL_and_process_$i
    echo "#Copy duplicate for rMATS"
    echo "cp ${i}_deduplicated.bam ${i}_deepfake.bam" >> DL_and_process_$i
    echo "#Run bedtools" >> DL_and_process_$i
    echo "#Run rMATS" >> DL_and_process_$i
    echo "source /home/yangk4/rmats_install/rmatsenv/bin/activate" >> DL_and_process_$i
    echo "${i}_deduplicated.bam" >> DL_and_process_${i}_b1
    echo "${i}_deepfake.bam" >> DL_and_process_${i}_b2
    mkdir ./rMATS_$i
    echo "python2 /home/yangk4/rmats_install/rMATS.4.0.2/rMATS-turbo-Linux-UCS2/rmats.py --b1 DL_and_process_${i}_b1 --b2 DL_and_process_${i}_b2 --gtf %s --od ./rMATS_$i --nthread %s" >> DL_and_process_$i
    echo "deactivate" >> DL_and_process_$i
    echo "#Run MAJIQ" >> DL_and_process_$i
    echo "source /home/yangk4/majiq_2_install/env/bin/activate" >> DL_and_process_$i
    mkdir ./majiq_$i
    echo -e "[info]\\nreadlen=%s\\nbamdirs=./\\ngenome=%s\\ngenome_path=%s\\n[experiments]\\nWT=${i}_deduplicated" >> ./majiq_$i/settings.txt
    echo "majiq build %s -c ./majiq_$i/settings.txt --output ./majiq_$i/build --nproc %s" >> DL_and_process_$i
    echo "majiq psi ./majiq_$i/build/$i\_deduplicated.majiq --nproc %s --output ./majiq_$i --name $i" >> DL_and_process_$i
    echo "voila psi ./majiq_$i/$i.psi.voila --splicegraph ./majiq_$i/build/splicegraph.sql -o ./majiq_$i/voila/$i" >> DL_and_process_$i
    echo "deactivate" >> DL_and_process_$i
done

''' % (root_dir,files_to_process, #general params
memory,nprocs, #bbduck params
genomepath,nprocs, #STAR params
str(int(nprocs)-1), #samtools sort params
gtfpath, nprocs, #rMATS params
readlen,genome,genomepath, gffpath, nprocs, nprocs #MAJIQ params
))
    
    
#echo "/home/kyang1/TrimGalore-0.4.5/trim_galore --paired -stringency 5 -length 35 -q 20 -o ./ $i\_1.fastq $i\_2.fastq" >> DL_and_process_$i
#echo "mkdir ./fastqc_1_trimmed" >> DL_and_process_$i
#echo "mkdir ./fastqc_2_trimmed" >> DL_and_process_$i
#echo "fastqc $i\_1_val_1.fq -o ./fastqc_1_trimmed &" >> DL_and_process_$i
#echo "fastqc $i\_2_val_2.fq -o ./fastqc_2_trimmed" >> DL_and_process_$i
#echo "wait" >> DL_and_process_$i
    
#echo "rm ./*.fastq" >> DL_and_process_$i
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
#         fw.write('declare -a arr=(%s)\n' % arr_string.strip())
#         fw.write(
# '''
# for i in "${arr[@]}";
# do
#     mkdir ./STAR/$i
#     cd ./STAR/$i
#     STAR --genomeDir /project/barash_hdr1/STAR_genomes/%s/ --readFilesIn ../../$i\_trimtest/$i\_1_val_1.fq ../../$i\_trimtest/$i\_2_val_2.fq --runThreadN 7 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./$i. --outSAMattributes All --alignSJoverhangMin 8
#     rm ../../$i\_trimtest/*.fq
#     cd ../../
# done\n''' % genome)
    fw.close()