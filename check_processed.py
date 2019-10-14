import os
import sys


def check_vdb_validate(vdb_validate_results):
    """
    vdb_validate_results is a list of files, each having output generated
        by vdb-validate 2> result
    """
    res = list()
    for vdb_res in vdb_validate_results:
        if not os.path.exists(vdb_res):
            print("%s doesn't exist ..." % vdb_res)
            res.append(False)
            continue
        md5_oks = 0  # 5
        appears_valid = 0  # 1
        is_consistent = 0  # 1
        with open(vdb_res, "r") as handle:
            for line in handle:
        #         if "info: Column" in line:
        #             if not line.strip("\n\r").endswith("ok"):
        #                 print("%s failed vdb-validate Column..." % vdb_res)
        #                 res.append(False)
        #                 break
        #         if "info: Table" in line:
        #             if not line.strip("\n\r").endswith("ok") and not line.strip("\n\r").endswith("is consistent"):
        #                 print("%s failed vdb-validate Table..." % vdb_res)
        #                 res.append(False)
        #                 break
        #         if "info: Database" in line:
        #             if not "metadata: " in line:
        #                 if line.strip("\n\r").endswith("is consistent"):
        #                     print(
        #                         "%s failed vdb-validate not consistent..." % vdb_res)
        #                     res.append(False)
        #                     break
        # res.append(True)
                if "md5 ok" in line.strip("\n\r"):
                    md5_oks += 1
                elif "appears valid" in line.strip("\n\r"):
                    appears_valid += 1
                elif "is consistent" in line.strip("\n\r"):
                    is_consistent += 1
        if md5_oks < 1 or appears_valid < 1 or is_consistent < 1:
            print("%s didn't pass md5=%s, appears_valid=%s, is_consistent=%s..." % (vdb_res,
                                                                                md5_oks,
                                                                                appears_valid,
                                                                                is_consistent))
            res.append(False)
        else:
            res.append(True)
    return(res)


def check_fastq_dump(fastqdump_errorlogs, fastqdump_outlogs):
    error_lines = 0
    res = []
    spots = "spots read"
    readsread = "reads read"
    readwritten = "reads written"
    error_lines = 0
    for fastqdump_errorlog, fastqdump_outlog in zip(fastqdump_errorlogs, fastqdump_outlogs):
        if not os.path.exists(fastqdump_errorlog):
            print("%s doesn't exist ..." % fastqdump_errorlog)
            res.append(False)
            continue
        if not os.path.exists(fastqdump_outlog):
            print("%s doesn't exist ..." % fastqdump_outlog)
            res.append(False)
            continue
        with open(fastqdump_errorlog, "r") as handle:
            for line in handle:
                error_lines += 1
        if error_lines > 0:
            print("%s lines in %s" % (error_lines, fastqdump_errorlog))
            res.append(False)
            continue
        else:
            count = 0
            with open(fastqdump_outlog, "r") as handle:
                for line in handle:
                    if spots in line:
                        count += 1
                    if readsread in line:
                        count += 1
                    if readwritten in line:
                        count += 1
            if count != 3:
                print("%s/3 lines passed for %s" % (count, fastqdump_outlog))
                res.append(False)
                continue
        res.append(True)
    return(res)


def check_trim(trim_results, expected_number_lines=27):
    """
    trim_results is a list of files, each having output generated
        by bbduk on one sra. If fewer than 27, bbduk probably failed.
    """
    res = list()
    for trim_result in trim_results:
        if not os.path.exists(trim_result):
            print("%s doesn't exist!" % trim_result)
            res.append(False)
            continue
        lines = 0
        with open(trim_result, "r") as handle:
            for line in handle:
                lines += 1
        if lines < 27:
            print("%s didn't finish..." % trim_result)
            res.append(False)
            continue
        elif lines > 27:
            print("%s had more lines than expected (%s) ... " %
                  (trim_result, lines))
            res.append(False)
            continue
        else:
            res.append(True)
    return(res)


def check_majiq_build(majiq_logs, date):
    success = "MAJIQ Builder is ended successfully!"
    res = list()
    for log in majiq_logs:
        if not os.path.exists(log):
            print("%s doesn't exist ..." % log)
            res.append(False)
            continue
        passed = False
        with open(log, "r") as handle:
            for line in handle:
                if success in line:
                    if line.startswith(date):
                        passed = True
                        res.append(True)
                        break
        if not passed:
            print("%s failed ... " % log)
            res.append(False)
    return(res)
