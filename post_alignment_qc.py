import pipelineUtil
import os

def validate_bam_file(picard_path, bam_file, uuid, outdir, logger=None):
   """ Validate resulting post-alignment BAM file """

   if os.path.isfile(picard_path) and os.path.isfile(bam_file):
        tmp_dir = os.path.join(outdir, 'tmp')
        if not os.path.isdir(tmp_dir):
            os.mkdir(tmp_dir)
        cmd = ['java', '-jar', picard_path, "ValidateSamFile", "I=%s" %bam_file,
               "O=%s" %os.path.join(outdir, "%s.validate" %uuid), "VALIDATION_STRINGENCY=LENIENT",
               "TMP_DIR=%s" %tmp_dir]
        pipelineUtil.log_function_time("ValidateSAM", uuid, cmd, logger)
   else:
       raise Exception("Invalid path to picard or BAM")

def collect_rna_seq_metrics(picard_path, bam_file, uuid, outdir, ref_flat, logger=None):
    """ Collect RNA-seq metrics using Picard """

    if os.path.isfile(picard_path) and os.path.isfile(bam_file):
        tmp_dir = os.path.join(outdir, 'tmp')
        if not os.path.isdir(tmp_dir):
            os.mkdir(tmp_dir)
        cmd = ['java', '-jar', picard_path, "CollectRnaSeqMetrics", "METRIC_ACCUMULATION_LEVEL=READ_GROUP",
                "I=%s" %bam_file, "O=%s" %os.path.join(outdir, "%s.rna_seq_metrics.txt" %uuid), "STRAND=NONE",
                "REF_FLAT=%s" %ref_flat, "VALIDATION_STRINGENCY=LENIENT", "TMP_DIR=%s" %tmp_dir]
        pipelineUtil.log_function_time("RNAseq_metrics", uuid, cmd, logger)
    else:
        raise Exception("Invalid path to picard or bam")

def bam_index(bam_file, uuid, logger=None):
    """ Index the resultant post alignment BAM file """

    if os.path.isfile(bam_file):
        cmd = ['samtools', 'index', '-b', bam_file]
        pipelineUtil.log_function_time("BamIndex", uuid, cmd, logger)
    else:
        raise Exception("invalid bam file")
    assert(os.path.isfile('%s.bai' %bam_file))

def reorder_bam(picard_path, bam_file, uuid, outdir, ref_genome, logger=None):
    """ Reorder the BAM file according to the reference genome """

    if os.path.isfile(bam_file) and os.path.isfile(picard_path) and os.path.isfile(ref_genome):
        outbam = os.path.join(outdir, '%s.reorder.bam' %uuid)
        tmp_dir = os.path.join(outdir, 'tmp')
        if not os.path.isdir(tmp_dir):
            os.mkdir(tmp_dir)
        cmd = ['java', '-jar', picard_path, 'ReorderSam', 'I=%s' %bam_file, 'O=%s' %outbam, 'R=%s' %ref_genome,
                'VALIDATION_STRINGENCY=LENIENT', 'TMP_DIR=%s' %tmp_dir]
        pipelineUtil.log_function_time("picard_reorder_sam", uuid, cmd, logger)
    else:
        raise Exception("invalid bam, picard path or reference genome")
    return outbam

def rna_seq_qc(rna_seq_qc_path, bam_file, uuid, outdir, ref_genome, gtf, logger=None):
    """ Perform RNA-seqQC on post alignment BAM file """

    if os.path.isfile(bam_file) and os.path.isfile(rna_seq_qc_path) and os.path.isfile(gtf):
        cmd = ['java', '-jar', rna_seq_qc_path, '-o', outdir, '-r', ref_genome, '-s',
                '%s|%s|%s' %(uuid, bam_file, uuid), '-t', gtf]
        exit_code = pipelineUtil.log_function_time('RNAseq_qc', uuid, cmd, logger)
    else:
        raise Exception("Invalid path to rnaseq-qc or bam")
    print "exit_code received by rna_seq_qc = %s" %exit_code
    return exit_code

def add_or_replace_read_group(picard_path, bam_file,  outdir, uuid, rg_id, rg_lb="Unknown", rg_pl="Unknown", rg_pu="Unknown",rg_sm="Unknown", logger=None):
    outbam = '%s.addRG.bam' %os.path.join(outdir, uuid)
    if os.path.isfile(bam_file) and os.path.isfile(picard_path):
        tmp_dir = os.path.join(outdir, 'tmp')
        if not os.path.isdir(tmp_dir):
            os.mkdir(tmp_dir)
        cmd = ['java', '-jar', picard_path, 'AddOrReplaceReadGroups', 'I=%s' %bam_file, 'O=%s' %outbam,
                'RGID=%s'%rg_id, 'RGLB=%s' %rg_lb, 'RGPL=%s' %rg_pl, 'RGPU=%s' %rg_pu, 'RGSM=%s' %rg_sm,
                'VALIDATION_STRINGENCY=LENIENT','TMP_DIR=%s' %tmp_dir]
        pipelineUtil.log_function_time('AddOrReplaceReadGroups', uuid, cmd, logger)
    print outbam
    if os.path.isfile(outbam):
        print "returning file now %s" %outbam
        return outbam
    else:
        raise Exception('Could not add or replace read groups. Check log file for errors')

