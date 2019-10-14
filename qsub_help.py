def qsub_job_submission_line(jobname,
                             script,
                             elog,
                             olog,
                             email="",
                             email_if="-m as ",
                             wait_to_finish="-sync n ",
                             pe="-pe smp 1 ",
                             memory_per_core="500M"):
    """
    jobname
    script
    elog
    olog
    n_cores
    memory_per_core_Gb
    # THESE REQUIRE TRAILING SPACES:
    pe "-pe smp 1 "
    email "-M username@domain.com "
    email_if "-m aes "
        a: if aborted
        s: if suspended
        e: when job ends
    wait_to_finish "-sync n " or "-sync y "
    """
    jobline = "qsub "
    jobline += "-N %s " % jobname
    jobline += "-e %s " % elog
    jobline += "-o %s " % olog
    jobline += email
    jobline += email_if
    jobline += wait_to_finish
    jobline += pe
    jobline += "-l h_vmem=%s " % memory_per_core
    jobline += "-l m_mem_free=%s " % memory_per_core
    jobline += "-cwd "
    jobline += "%s\n" % script
    return(jobline)