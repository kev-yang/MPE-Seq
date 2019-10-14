import sys
sys.path.insert(
    0, '/mnt/isilon/thomas-tikhonenko_lab/user/radensc/rsynced_scripts/scripts')
from rna_seq_processing_steps import *
from check_processed import *
from qsub_help import *
import math
import pandas as pd

email = "cradens@biociphers.org"
working_directory = "/scr1/users/radensc/jsl1/"
majiq_env = "/home/radensc/bin/majiq2/majiq2_dev/bin/activate"
number_of_jobs_downloading = 8
max_number_of_jobs_processing = 8

n_parallel_cores = 4
memory_per_core_Gb = 10
total_memory_limit_Gb = memory_per_core_Gb * n_parallel_cores
total_memory_limit_Mb = (total_memory_limit_Gb) * 2**10

runinfo = os.path.join(working_directory, "run_info.txt")
run_info = pd.read_csv(runinfo, delimiter=",")


processed_reads_dir = os.path.join(working_directory, "processed_reads")
make_folder(processed_reads_dir)
all_script_dir = os.path.join(working_directory, "scripts")
make_folder(all_script_dir)
dl_scripts_dir = os.path.join(all_script_dir, "download_reads_scripts")
make_folder(dl_scripts_dir)
dl_scripts_logs_dir = os.path.join(dl_scripts_dir, "logs")
make_folder(dl_scripts_logs_dir)
proc_scripts_dir = os.path.join(all_script_dir, "process_reads_scripts")
make_folder(proc_scripts_dir)
proc_scripts_logs_dir = os.path.join(proc_scripts_dir, "logs")
make_folder(proc_scripts_logs_dir)


total_number_to_process = len(run_info["Run"])  # len(failed)
max_number_of_jobs_processing= total_number_to_process # temp
print("%s Total SRAs" % total_number_to_process)
number_downloads_per_job = int(
    math.ceil(total_number_to_process / number_of_jobs_downloading))
print("Download %s at a time (up to %s per job)" %
      (number_of_jobs_downloading, number_downloads_per_job))
number_processed_per_job = int(
    math.ceil(total_number_to_process / max_number_of_jobs_processing))
actual_number_of_jobs_processing = int(
    math.ceil(total_number_to_process / number_processed_per_job))
print("Process %s at a time (up to %s per job)" %
      (actual_number_of_jobs_processing, number_processed_per_job))

download_jobs = list()
process_jobs = list()


i = 0
this_download_job = list()
process_job = list()
n_dl_seen = 0
n_proc_seen = 0

for sra_id in run_info["Run"]:
    # if sra_id not in failed:
    #     continue
    i += 1
    outdir = os.path.join(processed_reads_dir, sra_id)
    logdir = os.path.join(outdir, "logs")
    scriptdir = os.path.join(outdir, "scripts")
    runlines = proc_lines_for_sra(sra_id=sra_id,
                                  out_directory=outdir,
                                  sra_project_working_directory="/scr1/userss/radensc/ncbi/",
                                  paired_or_single="paired",
                                  ### sra tools info ###
                                  sratools_bin_dir="/home/radensc/bin/sratoolkit.2.9.6-centos_linux64/bin/",

                                  ### gzipping info ###
                                  pigz_executable="/home/radensc/bin/pigz/pigz-2.4/pigz",
                                  pigz_threads=4,

                                  ### trimming info ###
                                  bbduk_executable="/home/radensc/bin/bbduk/bbmap/bbduk.sh",
                                  bbduk_adaptor="/home/radensc/bin/bbduk/bbmap/resources/adapters.fa",
                                  bbduk_threads=4,
                                  bbduk_memory_in_gb=10,

                                  ### fastqc info ###
                                  fastqc_executable="/home/radensc/bin/fastqc/FastQC/fastqc",
                                  fasterqdump_threads=4,

                                  ### aligner info ###
                                  star_executable="/home/radensc/bin/STAR/bin/Linux_x86_64/STAR",
                                  star_genome_dir="/scr1/users/radensc/genomic_info/human/hg38_GRCh38_Reference_Genome/star_genomes/hg38_len100_fromGTF",
                                  star_threads=4,

                                  ### samtools info ###
                                  samtools_executable="/home/radensc/bin/samtools-1.9/samtools",

                                  # gene expression quantification
                                  # info ###
                                  salmon_executable="/home/radensc/bin/salmon/salmon-latest_linux_x86_64/bin/salmon",
                                  salmon_genome="/scr1/users/radensc/genomic_info/human/hg38_GRCh38_Reference_Genome/salmon_info/salmon_hg38_index",
                                  salmon_libtype="A",
                                  salmon_threads=4,

                                  ### majiq info ###
                                  max_readlen=151,
                                  transcriptome_annotation="/scr1/users/radensc/genomic_info/human/hg38_GRCh38_Reference_Genome/transcriptome_annotation/Homo_sapiens.GRCh38.94.chr.gff3",
                                  strandness="None",
                                  genome_name="hg38",
                                  majiq_threads=4
                                  )
    this_download_job.append(runlines["runlines"]["prefetch_lines"])
    runline_name = ["validate_lines",
                    "fastqdump_lines",
                    "cleanup_sras",
                    "trim_lines",
                    "fastqc_lines",
                    "salmon_lines",
                    "align_lines",
                    "majiq_lines"]#,
                    #"cleanup_trimmed"]
    runline = [runlines["runlines"][x] for x in runline_name]
    runline = "".join(runline)
    if len(process_job) == 0:
        prepend = "\n\n### Activate MAJIQ environment ###\n"
        prepend += "source " + majiq_env + "\n\n"
        runline = prepend + runline
    process_job.append(runline)

    # every number_downloads_per_job'th, save download lines, re-initialize new
    # job
    if i % number_downloads_per_job == 0:
        download_jobs.append("".join(this_download_job))
        this_download_job = list()
    if i % number_processed_per_job == 0:
        process_jobs.append("".join(process_job))
        process_job = list()
    if i == total_number_to_process:
        if len(this_download_job) > 0:
            download_jobs.append("".join(this_download_job))
        if len(process_job) > 0:
            process_jobs.append("".join(process_job))
        break


print("%s download jobs..." % len(download_jobs))
print("%s process_jobs..." % len(process_jobs))

all_download_job_scripts = list()
for job in range(number_of_jobs_downloading):
    download_job_script = os.path.join(
        dl_scripts_dir, "download_%s.sh" % job)
    all_download_job_scripts.append(download_job_script)
    with open(download_job_script, "w") as handle:
        handle.write(download_jobs[job])


submit_dl_script = os.path.join(dl_scripts_dir, "submit_download_jobs.sh")
with open(submit_dl_script, "w") as handle:
    handle.write("cd %s\n" % dl_scripts_dir)
    for job in range(number_of_jobs_downloading):
        dl_script = all_download_job_scripts[job]
        dl_script_base = os.path.basename(dl_script)
        elog = os.path.join(dl_scripts_logs_dir, dl_script_base + ".error")
        olog = os.path.join(dl_scripts_logs_dir, dl_script_base + ".out")
        wait = "n"
        # wait for last submission to complete, so that master subission waits,
        # too
        if job == (number_of_jobs_downloading - 1):
            wait = "-sync y "
        else:
            wait = "-sync n "
        submissionline = qsub_job_submission_line(jobname="DL_%s" % job,
                                                  script=dl_script,
                                                  elog=elog,
                                                  olog=olog,
                                                  email="-M cradens@biociphers.org ",
                                                  email_if="-m as ",
                                                  wait_to_finish=wait,
                                                  pe="-pe smp 1 ",
                                                  memory_per_core="1G")
        handle.write(submissionline)


all_process_job_scripts = list()
for job in range(actual_number_of_jobs_processing):
    process_job_script = os.path.join(proc_scripts_dir, "process_%s.sh" % job)
    all_process_job_scripts.append(process_job_script)
    with open(process_job_script, "w") as handle:
        handle.write(process_jobs[job])


submit_process_script = os.path.join(
    proc_scripts_dir, "submit_process_jobs.sh")
with open(submit_process_script, "w") as handle:
    handle.write("cd %s\n" % proc_scripts_dir)
    for job in range(actual_number_of_jobs_processing):
        proc_script = all_process_job_scripts[job]
        proc_script_base = os.path.basename(proc_script)
        elog = os.path.join(proc_scripts_logs_dir, proc_script_base + ".error")
        olog = os.path.join(proc_scripts_logs_dir, proc_script_base + ".out")
        submissionline = qsub_job_submission_line(jobname="PROC_%s" % job,
                                                  script=proc_script,
                                                  elog=elog,
                                                  olog=olog,
                                                  email="-M cradens@biociphers.org ",
                                                  email_if="-m as ",
                                                  wait_to_finish="-sync n ",
                                                  pe="-pe smp %s " % n_parallel_cores,
                                                  memory_per_core="%sG" % memory_per_core_Gb)
        handle.write(submissionline)


submit_master = os.path.join(
    all_script_dir, "submit_download_read_proccesing_scripts.sh")
with open(submit_master, "w") as handle:
    submit_master_base = os.path.basename(submit_master)
    errlog = os.path.join(dl_scripts_logs_dir,
                          submit_master_base + ".downloads.error")
    outlog = os.path.join(dl_scripts_logs_dir,
                          submit_master_base + ".downloads.out")
    handle.write("cd %s\n" % all_script_dir)
    submissionline = qsub_job_submission_line(jobname="DL_BOSS",
                                              script=submit_dl_script,
                                              elog=errlog,
                                              olog=outlog,
                                              email="-M cradens@biociphers.org ",
                                              email_if="-m aes ",
                                              wait_to_finish="-sync y ",
                                              pe="-pe smp 1 ",
                                              memory_per_core="500M")
    handle.write(submissionline)
    errlog = os.path.join(proc_scripts_logs_dir,
                          submit_master_base + ".readprocessing.error")
    outlog = os.path.join(proc_scripts_logs_dir,
                          submit_master_base + ".readprocessing.out")
    submissionline = qsub_job_submission_line(jobname="PROC_BOSS",
                                              script=submit_process_script,
                                              elog=errlog,
                                              olog=outlog,
                                              email="-M cradens@biociphers.org ",
                                              email_if="-m aes ",
                                              wait_to_finish="-sync y ",
                                              pe="-pe smp 1 ",
                                              memory_per_core="500M")
    handle.write(submissionline)
