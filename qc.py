import os
import logging
import subprocess
import pipelineUtil

def fastqc(fastqc_path, reads_1, reads_2, rg_id_dir, analysis_id, logger):

    if not os.path.isdir(rg_id_dir):
        raise Exception("Invalid directory: %s")

    fastqc_results = "%s" %(os.path.join(rg_id_dir, "fastqc_results"))
    if not os.path.isdir(fastqc_results):
        os.mkdir(fastqc_results)
    cmd = [fastqc_path, reads_1, reads_2, '--outdir', fastqc_results, '--extract']
    print cmd
    pipelineUtil.log_function_time("FastQC", analysis_id, cmd, logger)
    #for dirname in os.listdir(fastqc_results):
        #dirname = os.path.join(fastqc_results, dirname)
        #if os.path.isdir(dirname):
           # for filename in os.listdir(dirname):
             #   if filename == "summary.txt":
                    #filename = "%s_%s_%s" %(analysis_id, os.path.basename(rg_id_dir), filename)
                    #filename = os.path.join(dirname, filename)


