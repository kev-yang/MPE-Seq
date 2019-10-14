import argparse
import os
import sys
import subprocess as sp


def check_file_exists(filepath):
    if not os.path.isfile(filepath):
        raise RuntimeError("%s doesn't exist..." % filepath)


def is_exe(filepath):
    check_file_exists(filepath)
    if not os.access(filepath, os.X_OK):
        raise RuntimeError("%s isn't executable..." % filepath)


def check_folder_exists(folderpath):
    if not os.path.isdir(folderpath):
        raise RuntimeError("%s doesn't exist..." % folderpath)


def gen_prefetch_lines(sras,
                       logouts,
                       sratools_bin_dir,
                       sra_project_working_directory):
    """
    sras is a list of sra ids
    logouts: a list of logout prefixes
    """
    prefetch = os.path.join(sratools_bin_dir, "prefetch")
    is_exe(prefetch)
    prefetch_lines = "\n\necho '### prefetch ###'\n"
    # need to be in sra working directory for prefetch to work
    prefetch_lines += "cd %s\n" % sra_project_working_directory
    for sra, logout in zip(sras, logouts):
        prefetch_lines += "%s " % prefetch
        prefetch_lines += "%s " % sra
        prefetch_lines += "--max-size 2000000000 "
        prefetch_lines += "--verbose "
        prefetch_lines += "> %s " % (logout + ".stdout")
        prefetch_lines += "2> %s\n" % (logout + ".error")
    prefetch_lines += "cd -\n"
    return(prefetch_lines)


def gen_vdb_validate_lines(sras,
                           sratools_bin_dir,
                           sra_project_working_directory,
                           logouts):
    """
    sras is a list of sra ids
    logouts: a list of logout prefixes

    """
    vdb_validate = os.path.join(sratools_bin_dir, "vdb-validate")
    is_exe(vdb_validate)
    vdb_validate_lines = "\n\necho '### vdb_validate ###'\n"
    # need to be in sra working directory for prefetch to work
    vdb_validate_lines += "cd %s\n" % sra_project_working_directory
    for sra, logout in zip(sras, logouts):
        vdb_validate_lines += "%s " % vdb_validate
        vdb_validate_lines += "-I no "
        vdb_validate_lines += "%s " % sra
        # 2> redirects std output to file
        vdb_validate_lines += "> %s " % (logout + ".stdout")
        # 2> redirects error output to file
        vdb_validate_lines += "2> %s\n" % (logout + ".error")
    vdb_validate_lines += "cd -\n"
    return(vdb_validate_lines)


def gen_zip_lines(files,
                  pigz_executable,
                  pigz_threads,
                  compression_level=6):
    """
    files is a list of files, or a nested list of files to be pigz compressed
    """
    is_exe(pigz_executable)
    pigz_lines = ""
    for file in files:
        if isinstance(file, list):
            pigz_lines += gen_zip_lines(files=file,
                                        pigz_executable=pigz_executable,
                                        pigz_threads=pigz_threads,
                                        compression_level=compression_level)
        else:
            pigz_lines += "%s " % pigz_executable
            pigz_lines += "-p %s " % pigz_threads
            # compression level (6 is default of gzip)
            pigz_lines += "-%s " % compression_level
            pigz_lines += "%s\n" % file
    return pigz_lines


def gen_fastqdump_lines(sras,
                        outdirectories,
                        logouts,
                        sratools_bin_dir,
                        sra_project_working_directory,
                        fasterqdump_threads,
                        pigz_executable,
                        pigz_threads,
                        is_paired):
    """
    sras is a list of sra ids
    outdirectories, logouts are lists of sample length as sras
    logouts: a list of logout prefixes
    """
    # fastqdump
    fasterqdump = os.path.join(sratools_bin_dir, "fasterq-dump")
    is_exe(fasterqdump)
    fastqdump_lines = "\n\necho '### fastq dump ###'\n"
    # need to be in sra working directory for fasterq-dump to work, I think
    fastqdump_lines += "cd %s\n" % sra_project_working_directory
    for this_sra, outd, logout in zip(sras, outdirectories, logouts):
        fastqdump_lines += "%s --verbose " % fasterqdump
        fastqdump_lines += "--threads %s " % fasterqdump_threads
        fastqdump_lines += "--split-3 "
        fastqdump_lines += "--outdir %s " % outd
        fastqdump_lines += "%s " % this_sra
        fastqdump_lines += "> %s " % (logout + ".stdout")
        fastqdump_lines += "2> %s\n" % (logout + ".error")
        # gzip the fastqs (deletes input file after zipping)
        if is_paired:
            the_fastqs = [os.path.join(outd, this_sra + r + ".fastq")
                          for r in ["_1", "_2"]]
        else:
            the_fastqs = os.path.join(outd, this_sra + ".fastq")
        # if pigz specified, use it! much faster...
        if os.path.exists(pigz_executable):
            fastqdump_lines += gen_zip_lines(files=[the_fastqs],
                                             pigz_executable=pigz_executable,
                                             pigz_threads=pigz_threads,
                                             compression_level=6)
        else:
            fastqdump_lines += "gzip %s\n" % (the_fastqs)
    fastqdump_lines += "cd -\n"
    return fastqdump_lines


def gen_trim_lines(fastq_files,
                   log_outs,
                   bbduk_executable,
                   bbduk_adaptor,
                   bbduk_threads,
                   bbduk_memory_in_gb,
                   cleanup_untrimmed=True):
    """
    fastq_files is a list of files, or pair of files if paired-end, or both!
        e.g. [mysample1.fastq.gz, 
                mysample2.fastq.gz, 
                [mysample3_R1.fastq.gz, mysample3_R2.fastq.gz],
                mysample4.fastq.gz,
                [mysample5_R1.fastq.gz,mysample5_R2.fastq.gz]]
    logouts: a list of logout prefixes
    """
    is_exe(bbduk_executable)
    check_file_exists(bbduk_adaptor)
    trim_lines = "\n\necho '### trim ###'\n"
    for fastqs, log_outs in zip(fastq_files, log_outs):
        is_pair = True if isinstance(
            fastqs, list) and len(fastqs) == 2 else False
        if is_pair:
            out1 = fastqs[0].replace(".fastq", ".trimmed.fastq")
            out2 = fastqs[1].replace(".fastq", ".trimmed.fastq")
        else:
            out1 = fastqs[0].replace(".fastq", ".trimmed.fastq")
        # trim
        trim_lines += "%s " % bbduk_executable
        trim_lines += "in=%s " % fastqs[0]
        if is_pair:
            trim_lines += "in2=%s " % fastqs[1]
        trim_lines += "out=%s " % out1
        if is_pair:
            trim_lines += "out2=%s " % out2
        trim_lines += "ref=%s " % bbduk_adaptor
        # trimq is trimming based on quality score, standard is 20
        # qin=33 is based on the assumption that samples are Sanger +33 quality
        # scores...
        trim_lines += "ktrim=r k=23 mink=11 hdist=1 minlength=35 tpe tbo qtrim=r trimq=20 qin=33 "
        trim_lines += "-Xmx%sg " % bbduk_memory_in_gb
        trim_lines += "threads=%s " % bbduk_threads
        # save error log (has trim statistics)
        logout = log_outs[0]
        trim_lines += "> %s " % (logout + ".stdout")
        trim_lines += "2> %s\n" % (logout + ".error")
        # cleanup
        if cleanup_untrimmed:
            trim_lines += "rm %s\n" % fastqs[0]
            if is_pair:
                trim_lines += "rm %s\n" % fastqs[1]
    return trim_lines


def gen_fastqc_lines(fastq_files,
                     sample_ids,
                     outdirectories,
                     fastqc_executable):
    """
    fastq_files is a list of files, or pair of files if paired-end, or both!
        e.g. [mysample1.fastq.gz, 
                mysample2.fastq.gz, 
                [mysample3_R1.fastq.gz, mysample3_R2.fastq.gz],
                mysample4.fastq.gz,
                [mysample5_R1.fastq.gz,mysample5_R2.fastq.gz]]
    sample_ids: list of names same length as fastq_files
    outdirectories: a list of output directories
    """
    is_exe(fastqc_executable)
    if not len(fastq_files) == len(sample_ids) == len(outdirectories):
        raise RuntimeError("Expected same # of fastqs, sample_ids, outdirs...")
    fastqc_lines = "\n\necho '### fastqc ###'\n"
    for fastq, s, out in zip(fastq_files, sample_ids, outdirectories):
        # if paired...
        if isinstance(fastq, list):
            for f in fastq:
                fastqc_lines += "%s " % fastqc_executable
                fastqc_lines += "%s " % f
                fastqc_lines += "-o %s\n" % out
        else:  # else not paired
            fastqc_lines += "%s " % fastqc_executable
            fastqc_lines += "%s " % f
            fastqc_lines += "-o %s\n" % out
    return fastqc_lines


def gen_align_lines(fastq_files,
                    sample_ids,
                    outdirectories,
                    logouts,
                    star_executable,
                    star_genome_dir,
                    star_threads,
                    samtools_executable,
                    cleanup_unsorted=True):
    """
    fastq_files is a list of files, or pair of files if paired-end, or both!
        e.g. [mysample1.fastq.gz, 
                mysample2.fastq.gz, 
                [mysample3_R1.fastq.gz, mysample3_R2.fastq.gz],
                mysample4.fastq.gz,
                [mysample5_R1.fastq.gz,mysample5_R2.fastq.gz]]
    sample_ids: list of names same length as fastq_files, used for STAR prefix
    outdirectories: a list of output directories, used for STAR prefix
    logouts: a list of logout prefixes
    cleanup_unsorted: rm unsorted bam?
    """
    is_exe(star_executable)
    is_exe(samtools_executable)
    if not len(fastq_files) == len(sample_ids) == len(outdirectories):
        raise RuntimeError("Expected same # of fastqs, sample_ids, outdirs...")
    align_lines = "\n\necho '### STAR alignment ###'\n"
    for fastq, sample_id, outdir, logout in zip(fastq_files, sample_ids, outdirectories, logouts):
        staroutfiles = (
            os.path.join(
                outdir,
                sample_id + ".")
        )
        logout = logout[0]
        alignedbam = staroutfiles + "Aligned.out.bam"
        final_bam_fp = os.path.join(outdir, sample_id + ".bam")
        # star
        align_lines += "%s --genomeDir " % star_executable
        align_lines += "%s " % star_genome_dir
        align_lines += "--readFilesIn %s " % (" ".join(fastq))
        align_lines += "--runThreadN %s " % (star_threads)
        align_lines += "--outSAMtype BAM Unsorted "
        align_lines += "--outFileNamePrefix %s " % staroutfiles
        align_lines += "--outSAMattributes All --alignSJoverhangMin 8 "
        align_lines += "--readFilesCommand zcat --outSAMunmapped Within "
        align_lines += "> %s " % (logout + ".stdout")
        align_lines += "2> %s\n" % (logout + ".error")

        # samtools sort
        align_lines += "%s sort " % samtools_executable
        align_lines += "-@ %s " % star_threads
        align_lines += "-o %s " % final_bam_fp
        align_lines += "%s\n" % alignedbam
        # samtools index
        align_lines += "%s index " % samtools_executable
        align_lines += "-@ %s " % star_threads
        align_lines += "%s\n" % final_bam_fp
        # cleanup non-sorted bam
        if cleanup_unsorted:
            align_lines += "rm %s\n" % alignedbam
    return(align_lines)


def gen_salmon_lines(fastq_files,
                     sample_ids,
                     outdirectories,
                     logouts,
                     salmon_executable,
                     salmon_genome,
                     salmon_libtype,
                     salmon_threads):
    """
    fastq_files is a list of files, or pair of files if paired-end, or both!
        e.g. [mysample1.fastq.gz, 
                mysample2.fastq.gz, 
                [mysample3_R1.fastq.gz, mysample3_R2.fastq.gz],
                mysample4.fastq.gz,
                [mysample5_R1.fastq.gz,mysample5_R2.fastq.gz]]
    sample_ids: list of names same length as fastq_files
    outdirectories: a list of output directories
    logouts: a list of logout prefixes

    """
    is_exe(salmon_executable)
    check_folder_exists(salmon_genome)
    if not len(fastq_files) == len(sample_ids) == len(outdirectories):
        raise RuntimeError("Expected same # of fastqs, sample_ids, outdirs...")
    salmon_lines = "\n\necho '### Salmon counts ###'\n"
    for fastq, sample_id, outdir, logout in zip(fastq_files, sample_ids, outdirectories, logouts):
        # check if pair of fastqs, or single fastq
        if isinstance(fastq, list) and len(fastq) == 2:
            is_pair = True
        else:
            is_pair = False
        logout = logout[0]
        salmon_lines += "%s quant " % salmon_executable
        salmon_lines += "-i %s " % salmon_genome
        salmon_lines += "-l %s " % salmon_libtype
        salmon_lines += "-p %s " % salmon_threads
        # devs say it improves sensitivity+specificity
        salmon_lines += "--validateMappings "
        if is_pair:
            salmon_lines += "-1 %s " % fastq[0]
            salmon_lines += "-2 %s " % fastq[1]
        else:
            salmon_lines += "-r %s " % fastq[0]
        salmon_res_dir = os.path.join(outdir, "salmon_quant_" + sample_id)
        salmon_lines += "-o %s " % salmon_res_dir
        salmon_lines += "> %s " % (logout + ".stdout")
        salmon_lines += "2> %s\n" % (logout + ".error")
        # rename quant.sf file to sample_id.quant.sf
        quant_sf = os.path.join(salmon_res_dir, "quant.sf")
        renamed_quant_sf = os.path.join(
            salmon_res_dir, "%s.quant.sf" % sample_id)
        salmon_lines += "mv %s %s\n" % (quant_sf, renamed_quant_sf)
    return(salmon_lines)


def gen_majiq_config(samplenames,
                     bamdirs,
                     strandnesses,
                     outdirectories,
                     readlength,
                     genome_name):
    inputs = [samplenames, bamdirs, strandnesses, outdirectories]
    it = iter(inputs)
    if not all(len(l) == len(inputs[0]) for l in it):
        raise ValueError('not all inputs have same length!')
    if readlength <= 0:
        raise ValueError("readlength needs to be greater than 0...")
    settings_fps = list()
    for name, bdir, strand, outd in zip(samplenames,
                                        bamdirs,
                                        strandnesses,
                                        outdirectories):
        settings_fp = os.path.join(outd, name + ".majiq.settings.txt")
        with open(settings_fp, "w") as handle:
            ll = handle.write("[info]\n")
            ll = handle.write("readlen=%s\n" % readlength)
            ll = handle.write("bamdirs=%s\n" % bdir)
            ll = handle.write("genome=%s\n" % genome_name)
            ll = handle.write("strandness=%s\n\n" % strand)
            ll = handle.write("[experiments]\n")
            ll = handle.write("%s=%s\n" % (name, name))
        settings_fps.append(settings_fp)
    return(settings_fps)


def gen_majiq_lines(bamnames,
                    bamdirs,
                    outdirectories,
                    strandnesses,
                    config_dirs,
                    logouts,
                    majiq_threads,
                    transcriptome_annotation,
                    max_readlen,
                    genome_name):
    """
    BTW, now that --integrator works, generate runlines for majiq to run on subsets of samples
    Later, you can --intergrator all your separate majiq builds!

    THESE ARE ALL LISTS OF SAME LENGTH:
    bamnames: bam basenames, without .bam
    bamdirs
    outdirectories
    strandnesses
    config_dirs: where to save config files?
    logouts: list of logout prefixes

    NOT LISTS:
    majiq_threads
    transcriptome_annotation
    max_readlen
    genome

    For each, bam, generate runlines and config file to run majiq build on that one bam.
    """
    inputs = [bamnames, bamdirs, outdirectories, strandnesses, config_dirs]
    it = iter(inputs)
    if not all(len(l) == len(inputs[0]) for l in it):
        raise ValueError('not all inputs have same length!')
    majiq_lines = "\n\necho '### majiq ###'\n"
    config_fps = gen_majiq_config(samplenames=bamnames,
                                  bamdirs=bamdirs,
                                  strandnesses=strandnesses,
                                  outdirectories=config_dirs,
                                  readlength=max_readlen,
                                  genome_name=genome_name)
    for name, bdir, odir, strand, config_fp, logout in zip(bamnames,
                                                           bamdirs,
                                                           outdirectories,
                                                           strandnesses,
                                                           config_fps,
                                                           logouts):
        majiq_lines += "majiq build --conf %s " % config_fp
        majiq_lines += "--nproc %s " % majiq_threads
        # only write junc files (increment join later!)
        majiq_lines += "--junc-files-only "
        majiq_lines += "-o %s " % odir
        majiq_lines += "%s " % transcriptome_annotation
        majiq_lines += "> %s " % (logout + ".stdout")
        majiq_lines += "2> %s\n" % (logout + ".error")
    return(majiq_lines)


def make_folder(dirpath):
    if not os.path.isdir(dirpath):
        os.mkdir(dirpath)


def initialize_folders(outdir):
    if not os.path.isdir(outdir):
        try:
            os.mkdir(outdir)
        except:
            raise RuntimeError("%s cannot be created..." % outdir)
    logs = os.path.join(outdir, "logs")
    make_folder(logs)
    fastqcs = os.path.join(outdir, "fastqc")
    make_folder(fastqcs)
    build = os.path.join(outdir, "majiq_build")
    make_folder(build)
    config = os.path.join(outdir, "majiq_configs")
    make_folder(config)
    star = os.path.join(outdir, "star")
    make_folder(star)
    return logs, fastqcs, build, config, star


def isPairedSRA(sra_id,
                sratools_bin_dir,
                sra_project_working_directory):
    """
    uses fastq-dump to check first read, if it is paired, it will have 8 lines (9 newlines)...
    based off of: https://www.biostars.org/p/139422/#139495
    """
    fastqdump = os.path.join(sratools_bin_dir, "fastq-dump")
    is_exe(fastqdump)
    cwd = os.getcwd()
    try:
        # need to be in sra_tools working dir
        os.chdir(sra_project_working_directory)
        contents = sp.getoutput(
            "%s -X 1 -Z --split-spot %s" % (fastqdump, sra_id))
    except:
        os.chdir(cwd)
        raise Exception("Error running fastq-dump on", sra_id)
    os.chdir(cwd)
    if(contents.count("\n") == 5):
        return False
    elif(contents.count("\n") == 9):
        return True
    else:
        raise Exception("Unexpected output from fast-dump on ", sra_id)


def in_directory(potential_sub_dir, directory):
    # make both absolute
    directory = os.path.join(os.path.realpath(directory), '')
    potential_sub_dir = os.path.join(os.path.realpath(potential_sub_dir), '')
    check_folder_exists(directory)
    # return true, if the common prefix of both is equal to directory
    # e.g. potential_sub_dir /a/b/c/d/ and directory is /a/b/, the common
    # prefix is /a/b/
    return os.path.commonprefix([potential_sub_dir, directory]) == directory


def gen_cleanup_lines(files_to_remove):
    clean_lines = ""
    for file in files_to_remove:
        if file.count("*") > 0:
            raise ValueError("'Gonna just stop you right there bro' \
                - friendly neighborhood spiderbro")
        clean_lines += "rm %s\n" % file
    return clean_lines


def proc_lines_for_sra(sra_id,
                       out_directory,
                       sra_project_working_directory,
                       paired_or_single,
                       ### sra tools info ###
                       sratools_bin_dir,
                       ### gzipping info ###
                       pigz_executable,
                       pigz_threads,
                       ### trimming info ###
                       bbduk_executable,
                       bbduk_adaptor,
                       bbduk_threads,
                       bbduk_memory_in_gb,
                       ### fastqc info ###
                       fastqc_executable,
                       fasterqdump_threads,
                       ### aligner info ###
                       star_executable,
                       star_genome_dir,
                       star_threads,
                       ### samtools info ###
                       samtools_executable,
                       ### gene expression quantification info ###
                       salmon_executable,
                       salmon_genome,
                       salmon_libtype,
                       salmon_threads,
                       ### majiq info ###
                       max_readlen,
                       transcriptome_annotation,
                       strandness,
                       genome_name,
                       majiq_threads
                       ):
    # check that arguments.sra_project_working_directory exists
    check_folder_exists(sra_project_working_directory)
    # check if user set outdir to be in the sra-tools working directory
    if in_directory(out_directory, sratools_bin_dir):
        raise RuntimeError("Please change output to a subdirectory *not* in\
         your sra tools project directory...")
    # create /logs /etc folders inside of out_directory
    logdir, fastqcdir, builddir, configdir, stardir = initialize_folders(
        out_directory)
    # If user unsure whether sample is paired or didn't bother to set this arg,
    #   try and guess whether sample is paird or not from SRA ID
    if paired_or_single == "unknown":
        try:
            is_paired = isPairedSRA(
                sra_id, sratools_bin_dir, sra_project_working_directory)
            if is_paired:
                thisis = "paired"
            else:
                thisis = "single"
            print("Determined %s is %s-end ..." % (sra_id, thisis))
        except:
            raise RuntimeError("Unable to determine if %s is paired, \
                please specify with --paired_or_single" % sra_id)
    elif paired_or_single == "paired":
        is_paired = True
    else:
        is_paired = False
    # filepaths of intermediate files
    sra_fp = os.path.join(
        sra_project_working_directory, "sra", sra_id + ".sra")
    sra_cache_fp = sra_fp + ".vdbcache.cache"
    prefetch_log = os.path.join(logdir, sra_id + ".prefetch")
    vdb_val_log = os.path.join(logdir, sra_id + ".vdb_validate")
    fastqdump_log = os.path.join(logdir, sra_id + ".fastq_dump")
    trim_log = os.path.join(logdir, sra_id + ".bbduk_trim")
    align_log = os.path.join(logdir, sra_id + ".star_align")
    salmon_log = os.path.join(logdir, sra_id + ".salmon")
    majiq_log = os.path.join(logdir, sra_id + ".majiq")

    if is_paired:
        f1 = os.path.join(out_directory,
                          sra_id + "_1.fastq.gz")
        f2 = os.path.join(out_directory,
                          sra_id + "_2.fastq.gz")
        fastq_files = [f1, f2]
        t1 = f1.replace(".fastq", ".trimmed.fastq")
        t2 = f2.replace(".fastq", ".trimmed.fastq")
        trimmed_fastqs = [t1, t2]
        trim_logs = [trim_log, trim_log]
        salmon_logs = [salmon_log, salmon_log]
        align_logs = [align_log, align_log]
    else:
        f1 = os.path.join(out_directory, sra_id + ".fastq.gz")
        fastq_files = [f1]
        t1 = f1.replace(".fastq", ".trimmed.fastq")
        trimmed_fastqs = [t1]
        trim_logs = [trim_log]
        salmon_logs = [salmon_log]
        align_logs = [align_log]
    bam_fp = os.path.join(stardir, sra_id + ".bam")
    prefetch_lines = gen_prefetch_lines(sras=[sra_id],
                                        logouts=[prefetch_log],
                                        sratools_bin_dir=sratools_bin_dir,
                                        sra_project_working_directory=sra_project_working_directory)
    validate_lines = gen_vdb_validate_lines(sras=[sra_id],
                                            sratools_bin_dir=sratools_bin_dir,
                                            sra_project_working_directory=sra_project_working_directory,
                                            logouts=[vdb_val_log])
    fastqdump_lines = gen_fastqdump_lines(sras=[sra_id],  # just give SRA ID. fasterqdump knows where local sra is
                                          outdirectories=[out_directory],
                                          logouts=[fastqdump_log],
                                          sratools_bin_dir=sratools_bin_dir,
                                          sra_project_working_directory=sra_project_working_directory,
                                          fasterqdump_threads=fasterqdump_threads,
                                          pigz_executable=pigz_executable,
                                          pigz_threads=pigz_threads,
                                          is_paired=is_paired)
    cleanup_sras = gen_cleanup_lines(files_to_remove=[sra_fp, sra_cache_fp])
    trim_lines = gen_trim_lines(fastq_files=[fastq_files],
                                log_outs=[trim_logs],
                                bbduk_executable=bbduk_executable,
                                bbduk_adaptor=bbduk_adaptor,
                                bbduk_threads=bbduk_threads,
                                bbduk_memory_in_gb=bbduk_memory_in_gb,
                                cleanup_untrimmed=True)
    fastqc_lines = gen_fastqc_lines(fastq_files=[trimmed_fastqs],
                                    sample_ids=[sra_id],
                                    outdirectories=[fastqcdir],
                                    fastqc_executable=fastqc_executable)
    salmon_lines = gen_salmon_lines(fastq_files=[trimmed_fastqs],
                                    sample_ids=[sra_id],
                                    outdirectories=[out_directory],
                                    logouts=[salmon_logs],
                                    salmon_executable=salmon_executable,
                                    salmon_genome=salmon_genome,
                                    salmon_libtype=salmon_libtype,
                                    salmon_threads=salmon_threads)
    align_lines = gen_align_lines(fastq_files=[trimmed_fastqs],
                                  sample_ids=[sra_id],
                                  outdirectories=[stardir],
                                  logouts=[align_logs],
                                  star_executable=star_executable,
                                  star_genome_dir=star_genome_dir,
                                  star_threads=star_threads,
                                  samtools_executable=samtools_executable,
                                  cleanup_unsorted=True)
    majiq_lines = gen_majiq_lines(bamnames=[sra_id],
                                  bamdirs=[stardir],
                                  outdirectories=[builddir],
                                  logouts=[majiq_log],
                                  strandnesses=[strandness],
                                  config_dirs=[configdir],
                                  majiq_threads=majiq_threads,
                                  transcriptome_annotation=transcriptome_annotation,
                                  max_readlen=max_readlen,
                                  genome_name=genome_name)
    cleanup_trimmed = gen_cleanup_lines(files_to_remove=trimmed_fastqs)
    cleanup_bams = gen_cleanup_lines(files_to_remove=[bam_fp])

    res = dict()
    res["suggested_run_order"] = ["prefetch_lines",
                                  "validate_lines",
                                  "fastqdump_lines",
                                  "cleanup_sras",
                                  "trim_lines",
                                  "fastqc_lines",
                                  "salmon_lines",
                                  "align_lines",
                                  "majiq_lines",
                                  "cleanup_trimmed",
                                  "cleanup_bams"]
    res["runlines"] = dict()
    res["runlines"]["prefetch_lines"] = prefetch_lines
    res["runlines"]["validate_lines"] = validate_lines
    res["runlines"]["fastqdump_lines"] = fastqdump_lines
    res["runlines"]["cleanup_sras"] = cleanup_sras
    res["runlines"]["trim_lines"] = trim_lines
    res["runlines"]["fastqc_lines"] = fastqc_lines
    res["runlines"]["majiq_lines"] = majiq_lines
    res["runlines"]["align_lines"] = align_lines
    res["runlines"]["salmon_lines"] = salmon_lines
    res["runlines"]["cleanup_trimmed"] = cleanup_trimmed
    res["runlines"]["cleanup_bams"] = cleanup_bams
    return res

