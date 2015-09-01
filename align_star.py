import pipelineUtil
import qc
import argparse
import glob
import os
import re
import post_alignment_qc
import setupLog
import logging

def decompress(filename, workdir):
    """ Unpack fastq files """

    if filename.endswith(".tar"):
        cmd = ['tar', 'xvf', filename, '-C', workdir]
    elif filename.endswith(".gz"):
        cmd = ['tar', 'xzvf', filename, '-C', workdir]
    elif filename.endswith(".bz"):
        cmd = ['tar', 'xvjf', filename, '-C', workdir]
    else:
        raise Exception('Unknown input file extension for file %s' % filename)
    pipelineUtil.log_function_time("tar", filename, cmd)


def scan_workdir_helper(dirname, extension):
    fastq_files = glob.glob(os.path.join(dirname, "*_[12].%s"%(extension)))
    all_read_groups = list()
    read_group_set = dict()
    if not (fastq_files == []):
        for filename in fastq_files:
            rg_id = re.sub(r'_[12].%s$'%(extension), '', filename)
            read_group_set[rg_id] = read_group_set.get(rg_id, 0) + 1
        if not all(i == 2 for i in read_group_set.values()):
            raise Exception("Missing Pair")
        print read_group_set
        for rg_id in read_group_set.keys():
            reads_1 = "%s_1.%s" %(rg_id, extension)
            reads_2 = "%s_2.%s" %(rg_id, extension)
            read_pair = (os.path.basename(rg_id), reads_1, reads_2)
            all_read_groups.append(read_pair)

    return all_read_groups

def scan_workdir(dirname):
    """ Select the unpacked fastq files """

    print dirname
    fastq_files = scan_workdir_helper(dirname, "fastq")
    if fastq_files == []:
        fastq_files = scan_workdir_helper(dirname, "fastq.gz")
        if fastq_files == []:
            fastq_files = scan_workdir_helper(dirname, "fastq.bz")
    return fastq_files

def post_aln_qc(args, bam_file, logger=None):
    """ perform post alignment quality check """

    #validate the post-alignment BAM file
    post_alignment_qc.validate_bam_file(args.picard, bam_file, args.id, args.workDir, logger)

    #collect RNA-seq metrics
    post_alignment_qc.collect_rna_seq_metrics(args.picard, bam_file, args.id,
                                              args.workDir, args.ref_flat, logger)
    #run rna_seq_qc from broad institute
    post_alignment_qc.bam_index(bam_file, args.id, logger)

    exit_code = post_alignment_qc.rna_seq_qc(args.rna_seq_qc_path, bam_file, args.id, args.workDir,
                                args.ref_genome,args.rna_seq_qc_annotation, logger)

    if not(exit_code == 0):

        reordered_bam = post_alignment_qc.reorder_bam(args.picard, bam_file, args.id, args.workDir,
                                                args.ref_genome, logger)
        post_alignment_qc.bam_index(reordered_bam, args.id, logger)
        post_alignment_qc.rna_seq_qc(args.rna_seq_qc_path, reordered_bam, args.id, args.workDir,
                                args.ref_genome,args.rna_seq_qc_annotation, logger)

def bam_to_fastq(fastq_dir, bam_file, analysis_id, logger=None):
    """ Convert input BAM to Fastq files """

    tmp_fastq = os.path.join(fastq_dir, 'tmp')

    cmd = ['bamtofastq', 'filename=%s' %bam_file, 'outputdir=%s' %fastq_dir,
            'tryoq=1', 'collate=1', 'outputperreadgroup=1', 'T=%s' %tmp_fastq]

    exit_code = pipelineUtil.log_function_time('Biobambam', analysis_id, cmd, logger)

    if exit_code == 0:
        for filename in os.listdir(fastq_dir):
            if filename.endswith(".fq"):
                new_filename = filename.replace(".fq", ".fastq")
                os.rename(os.path.join(fastq_dir, filename), os.path.join(fastq_dir, new_filename))
    else:
        logger.error("Biobambam BamToFastq conversion of %s returned a non-zero exit code %s"
                    %(analysis_id, exit_code))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extension of RNA-seq alignment with QC.")
    required = parser.add_argument_group("Required input parameters")
    required.add_argument("--genomeDir", default=None, help="Directory containing the reference genome index", required=True)
    required.add_argument("--input_file", default=None, help="Input file containing the sequence information", required=True)
    required.add_argument("--ref_genome", default=None, help="path to reference genome", required=True)
    required.add_argument("--genome_annotation", default=None, help="path to genome annotation file", required=True)
    required.add_argument("--rna_seq_qc_annotation", default=None, help="annotation file for RNAseq-QC", required=True)

    optional = parser.add_argument_group("optional input parameters")
    optional.add_argument('--fastqc_path', help='path to fastqc', default='/home/ubuntu/bin/FastQC/fastqc')
    optional.add_argument("--out", default="out.bam", help="Name of the output BAM file")
    optional.add_argument("--workDir", default="./", help="Work directory")
    optional.add_argument("--metaDataTab", default=None, help="File containing metadata for the alignment header")
    optional.add_argument("--id", default=None, help="Analysis ID to be considered in the metadata file")
    optional.add_argument("--keepJunctions", default=False, action='store_true', help="keeps the junction file as {--out}.junctions")
    optional.add_argument("--useTMP", default=None, help="environment variable that is used as prefix for temprary data")
    #optional.add_argument("-h", "--help", action='store_true', help="show this help message and exit")
    optional.add_argument("--icgc_pipeline", default="/home/ubuntu/icgc_rnaseq_align/star_align.py", help="path to icgc star alignment pipeline")
    optional.add_argument("--picard", default="/home/ubuntu/bin/picard-tools-1.136/picard.jar", help="path to picard")
    optional.add_argument("--rna_seq_qc_path", default="/home/ubuntu/bin/RNA-SeQC_v1.1.8.jar", help="path to RNASeq-QC")
    optional.add_argument("--ref_flat", default="/home/ubuntu/SCRATCH/grch38/gencode.v21.annotation.ref_flat_final", help="path to ref flat file")

    star = parser.add_argument_group("STAR input parameters")
    star.add_argument("--runThreadN", type=int, default=4, help="Number of threads")
    star.add_argument("--outFilterMultimapScoreRange", type=int, default=1, help="outFilterMultimapScoreRange")
    star.add_argument("--outFilterMultimapNmax", type=int, default=20, help="outFilterMultimapNmax")
    star.add_argument("--outFilterMismatchNmax", type=int, default=10, help="outFilterMismatchNmax")
    star.add_argument("--alignIntronMax", type=int, default=500000, help="alignIntronMax")
    star.add_argument("--alignMatesGapMax", type=int, default=1000000, help="alignMatesGapMax")
    star.add_argument("--sjdbScore", type=int, default=2, help="sjdbScore")
    star.add_argument("--limitBAMsortRAM", type=int, default=0, help="limitBAMsortRAM")
    star.add_argument("--alignSJDBoverhangMin", type=int, default=1, help="alignSJDBoverhangMin")
    star.add_argument("--genomeLoad", default="NoSharedMemory", help="genomeLoad")
    star.add_argument("--genomeFastaFiles", default=None, help="genome sequence in fasta format to rebuild index")
    star.add_argument("--outFilterMatchNminOverLread", type=float, default=0.33, help="outFilterMatchNminOverLread")
    star.add_argument("--outFilterScoreMinOverLread", type=float, default=0.33, help="outFilterScoreMinOverLread")
    star.add_argument("--twopass1readsN", type=int, default=-1, help="twopass1readsN (-1 means all reads are used for remapping)")
    star.add_argument("--sjdbOverhang", type=int, default=100, help="sjdbOverhang (only necessary for two-pass mode)")
    star.add_argument("--outSAMstrandField", default="intronMotif", help="outSAMstrandField")
    star.add_argument("--outSAMattributes", default=["NH", "HI", "NM", "MD", "AS", "XS"], help="outSAMattributes")
    star.add_argument("--outSAMunmapped", default="Within", help="outSAMunmapped")
    star.add_argument("--outSAMtype", default=["BAM", "SortedByCoordinate"], help="outSAMtype")
    star.add_argument("--outSAMheaderHD", default=["@HD", "VN:1.4"], help="outSAMheaderHD")
    star.add_argument("--outSAMattrRGline", default=None, help="RG attribute line submitted to outSAMattrRGline")
    star.add_argument("--outSAMattrRGfile", default=None, help="File containing the RG attribute line submitted to outSAMattrRGline")

    args = parser.parse_args()

    log_file = "%s.log" % os.path.join(args.workDir, "%s" %args.id)
    logger = setupLog.setup_logging(logging.INFO, "%s" %args.id, log_file)

    if not args.workDir:
        os.mkdir(args.workDir)

    fastq_dir = os.path.join(args.workDir, '%s_fastq_files' %args.id)

    if not os.path.isdir(fastq_dir):
        os.mkdir(fastq_dir)

        if(args.input_file.endswith('bam')):
            bam_to_fastq(fastq_dir,  args.input_file, args.id,  logger)
        else:
            decompress(args.input_file, fastq_dir)

    read_group_pairs = scan_workdir(fastq_dir)
    for (rg_id, reads_1, reads_2) in read_group_pairs:
        rg_id_dir = os.path.join(args.workDir, rg_id)
        if not os.path.isdir(rg_id_dir):
            os.mkdir(rg_id_dir)
            qc.fastqc(args.fastqc_path, reads_1, reads_2, rg_id_dir, rg_id, logger)

    cmd = ["python", args.icgc_pipeline,
            "--genomeDir", str(args.genomeDir),
            "--FastqFileIn", str(fastq_dir),
            "--workDir", str(args.workDir),
            "--out", str(args.out),
            "--genomeFastaFiles", str(args.genomeFastaFiles),
            "--runThreadN", str(args.runThreadN),
            "--outFilterMultimapScoreRange", str(args.outFilterMultimapScoreRange),
            "--outFilterMultimapNmax", str(args.outFilterMultimapNmax),
            "--outFilterMismatchNmax", str(args.outFilterMismatchNmax),
            "--alignIntronMax", str(args.alignIntronMax),
            "--alignMatesGapMax", str(args.alignMatesGapMax),
            "--sjdbScore", str(args.sjdbScore),
            "--limitBAMsortRAM", str(args.limitBAMsortRAM),
            "--alignSJDBoverhangMin", str(args.alignSJDBoverhangMin),
            "--genomeLoad", str(args.genomeLoad),
            "--outFilterMatchNminOverLread", str(args.outFilterMatchNminOverLread),
            "--outFilterScoreMinOverLread", str(args.outFilterScoreMinOverLread),
            "--twopass1readsN", str(args.twopass1readsN),
            "--sjdbOverhang", str(args.sjdbOverhang),
            "--outSAMstrandField", str(args.outSAMstrandField),
            "--outSAMunmapped", str(args.outSAMunmapped)
            ]

    if args.keepJunctions:
        cmd = cmd.append("--keepJunctions")
        cmd = cmd.append(str(args.keepJunctions))

    if not args.metaDataTab == None:
        cmd = cmd.append("--metaDataTab")
        cmd = cmd.append(str(args.metaDataTab))

    if not args.outSAMattrRGline == None:
        cmd = cmd.append("--outSAMattrRGline")
        cmd = cmd.append(str(args.outSAMattrRGline))

    if not args.outSAMattrRGfile == None:
        cmd = cmd.append("--outSAMattrRGfile")
        cmd = cmd.append(str(args.outSAMattrRGfile))

    logger.info('Starting Alignment with STAR')

    exit_code = pipelineUtil.log_function_time("STAR_ALIGN", args.id, cmd, logger)
    if exit_code == 0:
        logger.info('Starting post alignment QC')
        post_aln_qc(args, args.out, logger)
    else:
        logger.error('STAR returned a non-zero exit code %s' %exit_code)


