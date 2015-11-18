import logging
import os
import setupLog
import htseq
import argparse

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Raw counts using HT-seq")
    required = parser.add_argument_group("Required input paramters")
    required.add_argument("--bam", default=None, help="path to BAM file", required=True)
    required.add_argument("--genome_annotation", default=None, help="path to annotation file", required=True)
    required.add_argument("--outdir", default="./", help="path to output directory")

    optional = parser.add_argument_group("optional input parameters")
    optional.add_argument("--id", default="unknown", help="unique identifer")
    optional.add_argument("--tobucket", default="s3://bioinformatics_scratch/")

    args = parser.parse_args()

    log_file = "%s.htseq.log" % os.path.join(args.outdir, args.id)
    logger = setupLog.setup_logging(logging.INFO, args.id, log_file)

    if not os.path.isfile(args.bam):
        logger.info("Downloading %s" %args.id)
        pipelineUtil.download_from_cleversafe(logger, args.bam, args.outdir)

    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)

    exit_code, out_file_name = htseq.htseq_count(args.bam, args.id, args.genome_annotation, args.outdir, logger)

    if not exit_code:
        pipelineUtil.upload_to_cleversafe(logger, args.tobucket, out_file_name)
