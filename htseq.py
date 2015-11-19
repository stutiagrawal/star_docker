import os
import subprocess
import signal

def htseq_count(bam, analysis_id, annotation, outdir, logger):
    """ Get raw counts using HTseq """

    if os.path.isfile(bam) and os.path.isfile(annotation):

        mapped_cmd = ["samtools", "view", "-F", "4", bam]
        count_cmd = ["htseq-count", "-m", "intersection-nonempty", "--idattr",
                    "gene_id", "-r", "pos", "-", annotation]


        #Python ignores SIGPIPE when started, therefore SIGPIPE needs to be restored.
        #More at: https://blog.nelhage.com/2010/02/a-very-subtle-bug/
        mapped = subprocess.Popen(mapped_cmd, stdout=subprocess.PIPE,
                                  preexec_fn=lambda:signal.signal(signal.SIGPIPE, signal.SIG_DFL))

        logger.info("Starting HTseq")
        out_file_name = os.path.join(outdir, "%s.htseq.counts" %analysis_id)
        with open(out_file_name, "w") as outfile:
            child = subprocess.Popen(count_cmd, stdin=mapped.stdout, stdout=outfile,
                                              stderr=subprocess.PIPE)
            stdout, stderr = child.communicate()
            exit_code = child.returncode

        outfile.close()

        if logger != None:
            stderr = stderr.split("\n")
            for line in stderr:
                logger.info(line)

        return (exit_code, out_file_name)

    else:
        raise Exception("Invalid BAM file %s or annotation file %s. Please check the file exists and the path is correct." %(bam, annotation))

